#!/usr/bin/env Rscript
# Joint baseline -> measured-program displacement -> function rules. Pair inputs
# remain separate; one symmetric interaction rule is shared across every pair.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(glmnet)
  library(digest)
  library(ggplot2)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
root <- '/work/results/program_pilot/83_lincs_and_tahoe_joint_program_bridge'
previous <- '/work/results/program_pilot/78_native_dose_prism_functional_interpolation'
reference <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query'
original <- file.path(reference, 'inputs/original')
script_path <- '/work/src/program_pilot/83_lincs_and_tahoe_joint_program_bridge.R'
prism_dir <- '/work/cache/tahoe_prism_full_chain/prism'
native_doses <- c(.05, .5, 5)
anchor <- 'Lapatinib ditosylate'
baseline_rank <- 8L
molecular_ridge <- 1
decoder_rank <- 4L
latent_pairs <- which(upper.tri(matrix(0, decoder_rank, decoder_rank), diag = TRUE), arr.ind = TRUE)
nonlinear_features <- 80L + nrow(latent_pairs)
molecular_weight <- 10
observed_function_weight <- 1
within_function_weight <- .8
decoder_penalty <- .003
interaction_penalty <- .01
interaction_features <- 40L + nrow(latent_pairs)
maximum_iterations <- 500L
seed <- 77012L
families <- c('Tahoe only', 'LINCS only', 'Tahoe + LINCS')
assays <- c('secondary', 'primary', 'nci_FF', 'nci_FG')
source_totals <- c(prism_secondary = .50, prism_primary = .10, nci_single = .15,
  combination_total = .125, combination_departure = .125)
set.seed(seed)
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
options(error = function() {
  status('FAILED; see run log')
  traceback(5)
  quit(save = 'no', status = 1L)
})
pearson <- function(x, y) if (length(x) >= 4L && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
key <- function(cell, drug, dose) paste(cell, drug, sprintf('%.10g', dose), sep = '|')
dose_weights <- function(dose) {
  stopifnot(all(dose >= .05 - 1e-10), all(dose <= 5 + 1e-10))
  w <- matrix(0, length(dose), 3L)
  for (i in seq_along(dose)) {
    lo <- max(which(native_doses <= dose[i] + 1e-10))
    hi <- min(which(native_doses >= dose[i] - 1e-10))
    if (lo == hi) w[i, lo] <- 1 else {
      fraction <- log(dose[i] / native_doses[lo]) / log(native_doses[hi] / native_doses[lo])
      w[i, lo] <- 1 - fraction
      w[i, hi] <- fraction
    }
  }
  stopifnot(max(abs(rowSums(w) - 1)) < 1e-12)
  w
}
status('Preparing matched cells for the common nine drugs and three molecular sources')
setwd(original)
source('modular_pair_pilot/bank_operator.R')
archive <- readRDS('modular_pair_pilot/inputs.rds')
PREP <- archive$prep
settings78 <- qread(file.path(previous, 'models/settings.qs'))
lincs <- qread(file.path(root, 'data/lincs_program_observations.qs'))
selected_drugs <- settings78$selected_drugs[settings78$selected_drugs %in% lincs$panel]
functional_drugs <- selected_drugs
DRUGS <- archive$drugs[archive$drugs %in% selected_drugs]
PREP$m <- PREP$m[PREP$m$drug %in% selected_drugs, ]
SINGLE <- fread(file.path(previous, 'data/individual_training_outcomes.tsv'))
d <- fread(file.path(previous, 'data/combination_training_outcomes.tsv'))
SINGLE <- SINGLE[drug %in% selected_drugs]
d <- d[drugA %in% selected_drugs & drugB %in% selected_drugs]
evaluation <- fread(file.path(previous, 'data/held_single_predictions_frozen.tsv'))[
  family == 'Joint singles and combinations', .(evaluation_id = single_id, depmap_id, drug, dose,
    y = observed_log2, assay = 'secondary')]
evaluation <- evaluation[drug %in% selected_drugs]
old_fold <- fread(file.path(previous, 'data/recipient_folds.tsv'))
recipient_fold <- setNames(old_fold$fold, old_fold$depmap_id)
stopifnot(length(selected_drugs) == 9L, !AU %in% c(SINGLE$depmap_id, d$depmap_id, PREP$m$depmap_id),
  sum(evaluation$drug == anchor & evaluation$dose == 2.5) == 441L)
write_table(old_fold, 'recipient_folds')
write_table(SINGLE, 'individual_training_outcomes')
write_table(d, 'combination_training_outcomes')
write_table(SINGLE[, .(rows = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug)), by = source], 'functional_data_coverage')
rm(archive)
invisible(gc())

# Eight baseline directions are derived only from training RNA recipients. The
# same drug/dose slope coefficients receive RNA and functional supervision.
# Measured RNA means and the 40 gene-program coordinates are never optimized.
prepare <- function(si, pi, excluded, molecular_source) {
  tahoe <- as.data.table(PREP$m)[!depmap_id %in% excluded, c('depmap_id', 'drug', 'dose', P), with = FALSE]
  tahoe[, molecular_source := 'Tahoe']
  lin <- copy(lincs$metadata)[!depmap_id %in% excluded, c('depmap_id', 'drug', 'dose', P), with = FALSE]
  # These are full-coordinate estimates from partial observations, not measured
  # full-program labels. The observation fractions remain in RNA supervision.
  full_estimate <- sweep(as.matrix(lin[, P, with = FALSE]), 2L, diag(lincs$observer), '/')
  lin[, (P) := as.data.table(full_estimate)]
  lin[, molecular_source := 'LINCS']
  molecular <- if (molecular_source == 'Tahoe only') tahoe else
    if (molecular_source == 'LINCS only') lin else rbindlist(list(tahoe, lin), use.names = TRUE)
  ref_ids <- sort(unique(molecular$depmap_id))
  baseline <- getmatrix(PREP$raw, ref_ids, P)
  baseline_mean <- colMeans(baseline)
  baseline_sd <- apply(baseline, 2L, sd)
  stopifnot(all(baseline_sd > 0), length(ref_ids) >= 6L)
  centered <- molecular[, lapply(.SD, function(z) z - mean(z)), by = .(drug, dose, molecular_source), .SDcols = P]
  response_sd <- vapply(centered[, P, with = FALSE], sd, numeric(1))
  stopifnot(all(response_sd > 0))
  training_baseline <- sweep(sweep(baseline, 2L, baseline_mean), 2L, baseline_sd, '/')
  baseline_projection <- svd(training_baseline, nu = 0, nv = baseline_rank)$v
  baseline_projection_scale <- sqrt(colMeans((training_baseline %*% baseline_projection)^2))
  stopifnot(all(baseline_projection_scale > 1e-8))
  ss <- SINGLE[si]
  pp <- d[pi]
  entities <- unique(rbind(ss[, .(depmap_id, drug)], pp[, .(depmap_id, drug = drugA)],
    pp[, .(depmap_id, drug = drugB)], molecular[, .(depmap_id, drug)]))
  query <- entities[, .(dose = native_doses), by = .(depmap_id, drug)]
  query_key <- key(query$depmap_id, query$drug, query$dose)
  nquery <- nrow(query)
  query_baseline <- getmatrix(PREP$raw, query$depmap_id, P)
  s <- sweep(sweep(sweep(query_baseline, 2L, baseline_mean), 2L, baseline_sd, '/') %*%
    baseline_projection, 2L, baseline_projection_scale, '/')
  molecular_response <- sweep(as.matrix(molecular[, P, with = FALSE]), 2L, response_sd, '/')
  molecular_index <- match(key(molecular$depmap_id, molecular$drug, molecular$dose), query_key)
  stopifnot(!anyNA(molecular_index), !any(query$depmap_id %in% excluded))
  group <- (match(query$drug, selected_drugs) - 1L) * 3L + match(query$dose, native_doses)
  group_rows <- split(seq_len(nquery), factor(group, levels = seq_len(length(selected_drugs) * 3L)))
  molecular_group <- (match(molecular$drug, selected_drugs) - 1L) * 3L + match(molecular$dose, native_doses)
  group_molecular_rows <- split(seq_len(nrow(molecular)), factor(molecular_group, levels = seq_along(group_rows)))
  group_baseline_mean <- matrix(0, length(group_rows), baseline_rank)
  group_response_mean <- matrix(0, length(group_rows), 40L)
  coefficient <- array(0, c(baseline_rank, 40L, length(group_rows)))
  molecular_weights <- matrix(0, nrow(molecular), 40L)
  initial_response <- matrix(0, nrow(query), 40L)
  for (g in seq_along(group_rows)) {
    ii <- group_rows[[g]]
    mi <- group_molecular_rows[[g]]
    stopifnot(length(mi) >= 3L)
    anchor_rows <- mi[molecular$molecular_source[mi] == 'Tahoe']
    if (!length(anchor_rows)) anchor_rows <- mi
    group_baseline_mean[g, ] <- colMeans(s[molecular_index[anchor_rows], , drop = FALSE])
    group_response_mean[g, ] <- colMeans(molecular_response[anchor_rows, , drop = FALSE])
    s[ii, ] <- sweep(s[ii, , drop = FALSE], 2L, group_baseline_mean[g, ])
    sources <- unique(molecular$molecular_source[mi])
    for (src in sources) {
      mr <- mi[molecular$molecular_source[mi] == src]
      confidence <- if (src == 'Tahoe') rep(1, 40L) else diag(lincs$observer) / mean(diag(lincs$observer))
      molecular_weights[mr, ] <- matrix(confidence / (length(sources) * length(mr)),
        length(mr), 40L, byrow = TRUE)
    }
    for (p in 1:40) {
      w <- molecular_weights[mi, p]
      fit <- glmnet(s[molecular_index[mi], , drop = FALSE],
        molecular_response[mi, p] - group_response_mean[g, p], weights = w,
        alpha = 0, lambda = molecular_ridge / sum(w), intercept = FALSE,
        standardize = FALSE, thresh = 1e-10)
      stopifnot(fit$jerr == 0L)
      coefficient[, p, g] <- as.numeric(coef(fit))[-1L]
    }
    initial_response[ii, ] <- sweep(s[ii, , drop = FALSE] %*% coefficient[, , g],
      2L, group_response_mean[g, ], '+')
  }
  dec_eigen <- eigen(crossprod(molecular_response), symmetric = TRUE)
  u <- dec_eigen$vectors[, seq_len(decoder_rank), drop = FALSE]
  u <- sweep(u, 2L, sqrt(colMeans((molecular_response %*% u)^2)), '/')
  query_indices <- function(cell, drug) {
    ix <- vapply(native_doses, function(dose) match(key(cell, drug, dose), query_key), integer(length(cell)))
    stopifnot(!anyNA(ix))
    ix
  }
  sparse_map <- function(cell, drug, dose, assay) {
    ix <- query_indices(cell, drug)
    columns <- ix + nquery * (match(assay, assays) - 1L)
    w <- dose_weights(dose)
    sparseMatrix(i = rep(seq_along(cell), 3L), j = as.vector(columns), x = as.vector(w),
      dims = c(length(cell), nquery * length(assays)))
  }
  ms <- sparse_map(ss$depmap_id, ss$drug, ss$dose, ss$assay)
  ma <- sparse_map(pp$depmap_id, pp$drugA, pp$doseA, paste0('nci_', pp$SCREENER))
  mb <- sparse_map(pp$depmap_id, pp$drugB, pp$doseB, paste0('nci_', pp$SCREENER))
  ia <- query_indices(pp$depmap_id, pp$drugA)
  ib <- query_indices(pp$depmap_id, pp$drugB)
  wa <- dose_weights(pp$doseA)
  wb <- dose_weights(pp$doseB)
  grid <- rbindlist(lapply(1:3, function(i) rbindlist(lapply(1:3, function(j)
    data.table(row = seq_len(nrow(pp)), a = ia[, i], b = ib[, j],
      assay = match(paste0('nci_', pp$SCREENER), assays), weight = wa[, i] * wb[, j])))))
  grid <- grid[weight > 0]
  pair_points <- unique(grid[, .(a, b, assay)])
  point_key <- paste(pair_points$a, pair_points$b, pair_points$assay)
  grid[, point := match(paste(a, b, assay), point_key)]
  mp <- sparseMatrix(i = grid$row, j = grid$point, x = grid$weight, dims = c(nrow(pp), nrow(pair_points)))
  metadata <- rbindlist(list(ss[, .(depmap_id, source, group, control_key, y)],
    pp[, .(depmap_id, source = 'combination_total',
      group = paste(SCREENER, pair, doseA, doseB, sep = '|'), control_key = paste0('total:', row_id), y)],
    pp[, .(depmap_id, source = 'combination_departure',
      group = paste(SCREENER, pair, doseA, doseB, sep = '|'), control_key = paste0('departure:', row_id), y = I)]))
  metadata[, center_group := paste(source, group, sep = '|')]
  w <- numeric(nrow(metadata))
  for (src in unique(metadata$source)) {
    ii <- which(metadata$source == src)
    z <- metadata[ii]
    unit <- 1 / as.numeric(table(z$control_key)[z$control_key])
    totals <- tapply(unit, z$group, sum)
    w[ii] <- source_totals[src] * unit / as.numeric(totals[z$group]) / uniqueN(z$group)
  }
  w <- w / sum(w)
  obs <- merge(ss[assay == 'secondary'], molecular[molecular_source == 'Tahoe',
    .(depmap_id, drug, dose, molecular_row = which(molecular$molecular_source == 'Tahoe'))],
    by = c('depmap_id', 'drug', 'dose'))
  obs_weight <- if (nrow(obs)) macro_weights(paste(obs$drug, obs$dose)) else numeric()
  make_loss_geometry <- function(y, weight, group) {
    if (!length(y)) return(list(y = numeric(), weight = numeric(), group = integer(),
      group_weight = numeric(), within_scale = 1, absolute_scale = 1))
    g <- match(group, unique(group))
    gw <- as.numeric(rowsum(weight, g, reorder = FALSE))
    ym <- as.numeric(rowsum(y * weight, g, reorder = FALSE)) / gw
    list(y = y, weight = weight, group = g, group_weight = gw,
      within_scale = max(sum(weight * (y - ym[g])^2), 1e-8), absolute_scale = max(sum(weight * y^2), 1e-8))
  }
  pair_initial_a <- initial_response[pair_points$a, , drop = FALSE]
  pair_initial_b <- initial_response[pair_points$b, , drop = FALSE]
  list(s = s, initial_response = initial_response, initial_coefficient = coefficient,
    group_response_mean = group_response_mean, group_molecular_rows = group_molecular_rows,
    molecular_weights = molecular_weights, molecular_source = molecular_source,
    query = query, nquery = nquery, group_rows = group_rows, ngroup = length(group_rows), u = u,
    molecular_response = molecular_response, molecular_index = molecular_index,
    observed_response = molecular_response[obs$molecular_row, , drop = FALSE],
    ms = ms, ma = ma, mb = mb, mp = mp, pair_points = pair_points,
    ns = nrow(ss), np = nrow(pp), metadata = metadata,
    function_geometry = make_loss_geometry(metadata$y, w, metadata$center_group),
    observed_geometry = make_loss_geometry(obs$y, obs_weight, paste(obs$drug, obs$dose)),
    observed_cells = unique(obs$depmap_id),
    context = list(response_sd = response_sd, baseline_mean = baseline_mean,
      baseline_sd = baseline_sd, baseline_projection = baseline_projection,
      baseline_projection_scale = baseline_projection_scale, group_baseline_mean = group_baseline_mean,
      excluded_cells = excluded, reference_cells = ref_ids),
    decoder_rms = sqrt(colMeans(initial_response^2)),
    nonlinear_rms = pmax(sqrt(colMeans(phi(initial_response, u)^2)), 1e-6),
    interaction_rms = pmax(sqrt(colMeans(pair_features(pair_initial_a, pair_initial_b, u)^2)), 1e-6))
}

phi <- function(x, u) {
  latent <- tanh(x %*% u)
  cbind(tanh(x), x^2 / (1 + x^2),
    latent[, latent_pairs[, 1], drop = FALSE] * latent[, latent_pairs[, 2], drop = FALSE])
}
decoder <- function(x, par, assay = NULL) {
  z <- x %*% par$u
  h <- phi(x, par$u)
  common <- as.numeric(x %*% par$w + h %*% par$a)
  if (is.null(assay)) {
    output <- matrix(common, nrow(x), length(assays))
    output[, -1L] <- output[, -1L] + x %*% par$deviation
  } else {
    output <- common
    for (j in 2:length(assays)) {
      ii <- which(assay == j)
      if (length(ii)) output[ii] <- output[ii] + as.numeric(x[ii, , drop = FALSE] %*% par$deviation[, j - 1L])
    }
  }
  list(output = output, z = z, h = h, x = x)
}
decoder_back <- function(node, gradient, par, assay = NULL) {
  if (is.null(assay)) {
    g <- rowSums(gradient)
    gd <- gradient[, -1L, drop = FALSE]
  } else {
    g <- as.numeric(gradient)
    gd <- vapply(2:length(assays), function(j) g * as.numeric(assay == j), numeric(length(g)))
  }
  dz <- matrix(0, nrow(node$x), decoder_rank)
  latent <- tanh(node$z)
  for (j in seq_len(nrow(latent_pairs))) {
    a <- latent_pairs[j, 1]
    b <- latent_pairs[j, 2]
    dz[, a] <- dz[, a] + g * par$a[80L + j] * latent[, b] * (1 - latent[, a]^2)
    dz[, b] <- dz[, b] + g * par$a[80L + j] * latent[, a] * (1 - latent[, b]^2)
  }
  dx <- g * (matrix(par$w, nrow(node$x), 40L, byrow = TRUE) +
    matrix(par$a[1:40], nrow(node$x), 40L, byrow = TRUE) * (1 - tanh(node$x)^2) +
    matrix(par$a[40L + 1:40], nrow(node$x), 40L, byrow = TRUE) * (2 * node$x / (1 + node$x^2)^2))
  list(x = dx + dz %*% t(par$u) + gd %*% t(par$deviation),
    w = as.numeric(crossprod(node$x, g)), a = as.numeric(crossprod(node$h, g)),
    u = crossprod(node$x, dz), deviation = crossprod(node$x, gd))
}
function_loss <- function(predicted, geometry) {
  if (!length(geometry$y)) return(list(value = 0, gradient = numeric()))
  error <- predicted - geometry$y
  mean_error <- as.numeric(rowsum(error * geometry$weight, geometry$group, reorder = FALSE)) / geometry$group_weight
  centered <- error - mean_error[geometry$group]
  value <- within_function_weight * sum(geometry$weight * centered^2) / geometry$within_scale +
    (1 - within_function_weight) * sum(geometry$weight * error^2) / geometry$absolute_scale
  gradient <- 2 * geometry$weight * (within_function_weight * centered / geometry$within_scale +
    (1 - within_function_weight) * error / geometry$absolute_scale)
  list(value = value, gradient = gradient)
}

# Both inputs are retained. Products are symmetric and vanish at a zero input;
# no drug, pair identity, baseline, or independent mean predictor enters J.
pair_features <- function(a, b, u) {
  ha <- tanh(a %*% u)
  hb <- tanh(b %*% u)
  h <- matrix(0, nrow(a), nrow(latent_pairs))
  for (j in seq_len(nrow(latent_pairs))) {
    p <- latent_pairs[j, 1]
    q <- latent_pairs[j, 2]
    h[, j] <- if (p == q) ha[, p] * hb[, p] else
      (ha[, p] * hb[, q] + ha[, q] * hb[, p]) / sqrt(2)
  }
  cbind(tanh(a) * tanh(b), h)
}
pair_back <- function(a, b, gradient, par) {
  g <- as.numeric(gradient)
  ha <- tanh(a %*% par$u)
  hb <- tanh(b %*% par$u)
  ga <- gb <- matrix(0, nrow(a), decoder_rank)
  for (j in seq_len(nrow(latent_pairs))) {
    p <- latent_pairs[j, 1]
    q <- latent_pairs[j, 2]
    c <- g * par$interaction[40L + j]
    if (p == q) {
      ga[, p] <- ga[, p] + c * hb[, p]
      gb[, p] <- gb[, p] + c * ha[, p]
    } else {
      ga[, p] <- ga[, p] + c * hb[, q] / sqrt(2)
      ga[, q] <- ga[, q] + c * hb[, p] / sqrt(2)
      gb[, p] <- gb[, p] + c * ha[, q] / sqrt(2)
      gb[, q] <- gb[, q] + c * ha[, p] / sqrt(2)
    }
  }
  a_gradient <- g * sweep((1 - tanh(a)^2) * tanh(b), 2L, par$interaction[1:40], '*') +
    (ga * (1 - ha^2)) %*% t(par$u)
  b_gradient <- g * sweep((1 - tanh(b)^2) * tanh(a), 2L, par$interaction[1:40], '*') +
    (gb * (1 - hb^2)) %*% t(par$u)
  list(a = a_gradient, b = b_gradient,
    interaction = as.numeric(crossprod(pair_features(a, b, par$u), g)))
}
decoder_vector <- function(par) c(par$w, par$a, par$deviation, par$interaction)
decoder_unpack <- function(theta, u) {
  list(w = theta[1:40], a = theta[40 + seq_len(nonlinear_features)], u = u,
    deviation = matrix(theta[40L + nonlinear_features + seq_len(120L)], 40L, 3L),
    interaction = tail(theta, interaction_features))
}
ndecoder <- 40L + nonlinear_features + 120L + interaction_features

# Warm start uses one shared regularized decoder. Measured and inferred RNA
# supervision have equal total weights, with within-treatment contrasts weighted
# more than absolute effects. No context-specific calibration heads are fitted.
initialize_decoder <- function(p, family) {
  par <- list(w = rep(0, 40L), a = rep(0, nonlinear_features), u = p$u,
    deviation = matrix(0, 40L, 3L), interaction = rep(0, interaction_features))
  basis <- cbind(p$initial_response, phi(p$initial_response, par$u))
  native_design <- lapply(seq_along(assays), function(a)
    cbind(basis, do.call(cbind, lapply(2:length(assays), function(j)
      p$initial_response * as.numeric(j == a)))))
  native_design <- cbind(do.call(rbind, native_design), matrix(0, p$nquery * length(assays), interaction_features))
  pair_a <- p$initial_response[p$pair_points$a, , drop = FALSE]
  pair_b <- p$initial_response[p$pair_points$b, , drop = FALSE]
  pair_design <- cbind(matrix(0, nrow(pair_a), ndecoder - interaction_features),
    pair_features(pair_a, pair_b, par$u))
  xs <- as.matrix(p$ms %*% native_design)
  xa <- as.matrix(p$ma %*% native_design)
  xb <- as.matrix(p$mb %*% native_design)
  xj <- as.matrix(p$mp %*% pair_design)
  x <- rbind(xs, xa + xb + xj, xj)
  observed_x <- cbind(p$observed_response, phi(p$observed_response, par$u),
    matrix(0, nrow(p$observed_response), 120L + interaction_features))
  views <- function(x, geo) {
    if (!length(geo$y)) return(list(x = x, y = numeric(), weight = numeric()))
    xm <- rowsum(x * geo$weight, geo$group, reorder = FALSE) / geo$group_weight
    ym <- as.numeric(rowsum(geo$y * geo$weight, geo$group, reorder = FALSE)) / geo$group_weight
    list(x = rbind(x - xm[geo$group, , drop = FALSE], x), y = c(geo$y - ym[geo$group], geo$y),
      weight = c(geo$weight * within_function_weight / geo$within_scale,
        geo$weight * (1 - within_function_weight) / geo$absolute_scale))
  }
  observed <- views(observed_x, p$observed_geometry)
  inferred <- views(x, p$function_geometry)
  fitting <- list(x = rbind(inferred$x, observed$x), y = c(inferred$y, observed$y),
    weight = c(inferred$weight, observed$weight * observed_function_weight))
  fit <- glmnet(fitting$x, fitting$y, weights = fitting$weight, alpha = 0,
    lambda = .003, intercept = FALSE, standardize = TRUE, thresh = 1e-7, maxit = 100000L)
  stopifnot(fit$jerr == 0L)
  coefficient <- as.numeric(coef(fit))[-1L]
  par$w <- coefficient[seq_len(40L)]
  par$a <- coefficient[40L + seq_len(nonlinear_features)]
  par$deviation <- matrix(coefficient[40L + nonlinear_features + seq_len(120L)], 40L, 3L)
  par$interaction <- tail(coefficient, interaction_features)
  par
}

fit_model <- function(p, family, job_tag, check_gradient = FALSE) {
  linked <- TRUE
  par0 <- initialize_decoder(p, family)
  initial <- decoder_vector(par0)
  na <- baseline_rank * 40L * p$ngroup
  if (linked) initial <- c(initial, as.numeric(p$initial_coefficient))
  evaluations <- 0L
  cache_theta <- cache_result <- NULL
  objective <- function(theta) {
    if (!is.null(cache_theta) && identical(theta, cache_theta)) return(cache_result)
    par <- decoder_unpack(theta[seq_len(ndecoder)], p$u)
    if (linked) {
      coefficient <- array(theta[ndecoder + seq_len(na)], c(baseline_rank, 40L, p$ngroup))
      response <- matrix(0, p$nquery, 40L)
      for (g in seq_len(p$ngroup)) {
        ii <- p$group_rows[[g]]
        response[ii, ] <- sweep(p$s[ii, , drop = FALSE] %*% coefficient[, , g],
          2L, p$group_response_mean[g, ], '+')
      }
    } else response <- p$initial_response
    n <- decoder(response, par)
    pair_a <- response[p$pair_points$a, , drop = FALSE]
    pair_b <- response[p$pair_points$b, , drop = FALSE]
    pair_h <- pair_features(pair_a, pair_b, par$u)
    single_prediction <- as.numeric(p$ms %*% as.vector(n$output))
    a_prediction <- as.numeric(p$ma %*% as.vector(n$output))
    b_prediction <- as.numeric(p$mb %*% as.vector(n$output))
    interaction_prediction <- as.numeric(p$mp %*% (pair_h %*% par$interaction))
    predicted <- c(single_prediction, a_prediction + b_prediction + interaction_prediction, interaction_prediction)
    fl <- function_loss(predicted, p$function_geometry)
    obs_node <- decoder(p$observed_response, par, rep(1L, nrow(p$observed_response)))
    ol <- function_loss(obs_node$output, p$observed_geometry)
    ns <- p$ns
    np <- p$np
    gs <- fl$gradient[seq_len(ns)]
    gp <- fl$gradient[ns + seq_len(np)]
    gi <- fl$gradient[ns + np + seq_len(np)]
    native_gradient <- as.numeric(crossprod(p$ms, gs) + crossprod(p$ma, gp) + crossprod(p$mb, gp))
    native_gradient <- matrix(native_gradient, p$nquery, length(assays))
    pair_gradient <- as.numeric(crossprod(p$mp, gp + gi))
    ng <- decoder_back(n, native_gradient, par)
    pg <- pair_back(pair_a, pair_b, pair_gradient, par)
    og <- decoder_back(obs_node, ol$gradient * observed_function_weight, par, rep(1L, nrow(p$observed_response)))
    grad_decoder <- list(w = ng$w + og$w, a = ng$a + og$a,
      deviation = ng$deviation + og$deviation, interaction = pg$interaction)
    penalty_value <- decoder_penalty * (sum((par$w * p$decoder_rms)^2) + sum((par$a * p$nonlinear_rms)^2) +
      8 * sum((par$deviation * p$decoder_rms)^2)) +
      interaction_penalty * sum((par$interaction * p$interaction_rms)^2)
    grad_decoder$w <- grad_decoder$w + 2 * decoder_penalty * par$w * p$decoder_rms^2
    grad_decoder$a <- grad_decoder$a + 2 * decoder_penalty * par$a * p$nonlinear_rms^2
    grad_decoder$deviation <- grad_decoder$deviation + 16 * decoder_penalty * par$deviation * p$decoder_rms^2
    grad_decoder$interaction <- grad_decoder$interaction +
      2 * interaction_penalty * par$interaction * p$interaction_rms^2
    value <- fl$value + observed_function_weight * ol$value + penalty_value
    gradient <- decoder_vector(grad_decoder)
    if (linked) {
      grad_response <- ng$x
      for (which_input in c('a', 'b')) {
        index <- p$pair_points[[which_input]]
        contribution <- rowsum(pg[[which_input]], index, reorder = FALSE)
        ii <- as.integer(rownames(contribution))
        grad_response[ii, ] <- grad_response[ii, , drop = FALSE] + contribution
      }
      grad_coefficient <- coefficient * (2 * molecular_weight * molecular_ridge / (p$ngroup * 40L))
      for (g in seq_len(p$ngroup)) {
        ii <- p$group_rows[[g]]
        mi <- p$group_molecular_rows[[g]]
        qi <- p$molecular_index[mi]
        molecular_error <- response[qi, , drop = FALSE] - p$molecular_response[mi, , drop = FALSE]
        value <- value + molecular_weight / (p$ngroup * 40L) *
          (sum(p$molecular_weights[mi, , drop = FALSE] * molecular_error^2) +
            molecular_ridge * sum(coefficient[, , g]^2))
        molecular_gradient <- 2 * molecular_weight * p$molecular_weights[mi, , drop = FALSE] *
          molecular_error / (p$ngroup * 40L)
        contribution <- rowsum(molecular_gradient, qi, reorder = FALSE)
        unique_qi <- as.integer(rownames(contribution))
        grad_response[unique_qi, ] <- grad_response[unique_qi, , drop = FALSE] + contribution
        grad_coefficient[, , g] <- grad_coefficient[, , g] +
          crossprod(p$s[ii, , drop = FALSE], grad_response[ii, , drop = FALSE])
      }
      gradient <- c(gradient, grad_coefficient)
    }
    stopifnot(is.finite(value), all(is.finite(gradient)))
    evaluations <<- evaluations + 1L
    cache_theta <<- theta
    cache_result <<- list(value = value, gradient = gradient)
    cache_result
  }
  if (check_gradient && linked) {
    probe <- initial
    probe[ndecoder + 1L] <- .013
    positions <- unique(c(1L, 41L, 49L, 100L, ndecoder - c(0L, 2L, 10L),
      ndecoder + c(1L, 40L, 159L, 2301L), length(probe) - c(0L, 39L, 120L)))
    analytical <- objective(probe)$gradient[positions]
    numerical <- vapply(positions, function(j) {
      step <- 1e-5 * max(1, abs(probe[j]))
      plus <- minus <- probe
      plus[j] <- plus[j] + step
      minus[j] <- minus[j] - step
      (objective(plus)$value - objective(minus)$value) / (2 * step)
    }, numeric(1))
    relative <- abs(numerical - analytical) / pmax(1, abs(numerical), abs(analytical))
    write_table(data.table(position = positions, analytical, numerical, relative),
      paste0('analytical_gradient_check_', gsub('[^A-Za-z0-9]+', '_', family)))
    stopifnot(max(relative) < 1e-5)
    a <- p$initial_response[p$pair_points$a[1:20], , drop = FALSE]
    b <- p$initial_response[p$pair_points$b[1:20], , drop = FALSE]
    stopifnot(max(abs(pair_features(a, b, p$u) - pair_features(b, a, p$u))) < 1e-12,
      max(abs(pair_features(a, a * 0, p$u))) == 0,
      max(abs(pair_features(b * 0, b, p$u))) == 0)
  }
  status(paste(job_tag, '| optimizing', family))
  # Program-coordinate units differ substantially. Parameter scaling conditions
  # the same objective; it does not change the representation or supervision.
  decoder_scale <- c(1 / p$decoder_rms, 1 / p$nonlinear_rms, rep(1 / p$decoder_rms, 3L), 1 / p$interaction_rms)
  parameter_scale <- if (linked) c(decoder_scale, rep(sqrt(p$ngroup / molecular_weight), na)) else decoder_scale
  opt <- optim(initial, function(theta) objective(theta)$value, function(theta) objective(theta)$gradient,
    method = 'L-BFGS-B', control = list(maxit = maximum_iterations, factr = 1e8,
      pgtol = 1e-4, lmm = 7, parscale = parameter_scale))
  if (opt$convergence == 1L && max(abs(objective(opt$par)$gradient * parameter_scale)) > 1e-3) {
    status(paste(job_tag, '| continuing optimization of', family))
    opt <- optim(opt$par, function(theta) objective(theta)$value, function(theta) objective(theta)$gradient,
      method = 'L-BFGS-B', control = list(maxit = maximum_iterations, factr = 1e8,
        pgtol = 1e-4, lmm = 7, parscale = parameter_scale))
  }
  par <- decoder_unpack(opt$par[seq_len(ndecoder)], p$u)
  model <- list(family = family, decoder = par, context = p$context,
    molecular_source = p$molecular_source, partial_observation_operator = lincs$observer,
    molecular_coefficient = if (linked) array(opt$par[ndecoder + seq_len(na)],
      c(baseline_rank, 40L, p$ngroup)) else p$initial_coefficient,
    initial_coefficient = p$initial_coefficient, group_response_mean = p$group_response_mean,
    selected_drugs = selected_drugs,
    functional_features = P, functional_baseline_inputs = FALSE,
    training_cells = unique(p$metadata$depmap_id),
    training_molecular_cells = p$context$reference_cells,
    optimization = list(convergence = opt$convergence, counts = opt$counts,
      maximum_gradient = max(abs(objective(opt$par)$gradient)),
      maximum_scaled_gradient = max(abs(objective(opt$par)$gradient * parameter_scale)), value = opt$value))
  status(paste(job_tag, '| finished', family, '| solver status', opt$convergence))
  model
}

# Native-dose molecular inference is the sole baseline-dependent operation.
predict_displacement <- function(model, query) {
  baseline <- getmatrix(PREP$raw, query$depmap_id, P)
  s <- sweep(sweep(sweep(baseline, 2L, model$context$baseline_mean),
    2L, model$context$baseline_sd, '/') %*% model$context$baseline_projection,
    2L, model$context$baseline_projection_scale, '/')
  group <- (match(query$drug, model$selected_drugs) - 1L) * 3L + match(query$dose, native_doses)
  stopifnot(!anyNA(group))
  response <- matrix(0, nrow(query), 40L)
  for (g in unique(group)) {
    ii <- which(group == g)
    centered <- sweep(s[ii, , drop = FALSE], 2L, model$context$group_baseline_mean[g, ])
    response[ii, ] <- sweep(centered %*% model$molecular_coefficient[, , g],
      2L, model$group_response_mean[g, ], '+')
  }
  response
}
predict_single <- function(model, metadata) {
  prediction <- numeric(nrow(metadata))
  w <- dose_weights(metadata$dose)
  for (i in 1:3) {
    q <- metadata[, .(depmap_id, drug, dose = native_doses[i])]
    response <- predict_displacement(model, q)
    prediction <- prediction + decoder(response, model$decoder, match(metadata$assay, assays))$output * w[, i]
  }
  prediction
}
predict_pairs <- function(model, metadata) {
  w_a <- dose_weights(metadata$doseA)
  w_b <- dose_weights(metadata$doseB)
  qa <- qb <- list()
  for (i in 1:3) {
    qa[[i]] <- predict_displacement(model, metadata[, .(depmap_id, drug = drugA, dose = native_doses[i])])
    qb[[i]] <- predict_displacement(model, metadata[, .(depmap_id, drug = drugB, dose = native_doses[i])])
  }
  assay <- match(paste0('nci_', metadata$SCREENER), assays)
  total <- a <- b <- numeric(nrow(metadata))
  for (i in 1:3) {
    a <- a + decoder(qa[[i]], model$decoder, assay)$output * w_a[, i]
    b <- b + decoder(qb[[i]], model$decoder, assay)$output * w_b[, i]
    for (j in 1:3) total <- total + as.numeric(pair_features(qa[[i]], qb[[j]], model$decoder$u) %*%
      model$decoder$interaction) * w_a[, i] * w_b[, j]
  }
  list(total = a + b + total, a = a, b = b, departure = total)
}

ref_a <- paste(d$SCREENER, d$depmap_id, d$drugA, sprintf('%.10g', d$doseA), sep = '|')
ref_b <- paste(d$SCREENER, d$depmap_id, d$drugB, sprintf('%.10g', d$doseB), sep = '|')
jobs <- c(lapply(1:3, function(k) list(tag = paste0('recipient_', k), regime = 'recipient', fold = k,
  held_cells = names(recipient_fold)[recipient_fold == k], held_partner = NA_character_)),
  lapply(setdiff(selected_drugs, anchor), function(dr)
    list(tag = paste0('partner_', gsub('[^A-Za-z0-9]+', '_', dr)), regime = 'partner', fold = NA_integer_,
      held_cells = character(), held_partner = dr)),
  list(list(tag = 'final', regime = 'final', fold = NA_integer_, held_cells = character(), held_partner = NA_character_)))
single_predictions <- pair_predictions <- molecular_predictions <- oracle_predictions <- audit <- list()
run_contract <- digest(list(digest(script_path, file = TRUE, algo = 'sha256'), families, selected_drugs,
  baseline_rank, molecular_ridge, molecular_weight, within_function_weight, decoder_penalty,
  interaction_penalty, maximum_iterations), algo = 'sha256')
for (job in jobs) {
  done_path <- file.path(root, 'models', paste0(job$tag, '_completed.qs'))
  if (file.exists(done_path)) {
    result <- qread(done_path)
    stopifnot(identical(result$run_contract, run_contract))
    status(paste(job$tag, '| restored completed models'))
  } else {
    status(paste(job$tag, '| preparing joint RNA and functional training data'))
    si <- which(!SINGLE$depmap_id %in% job$held_cells)
    pi <- which(!d$depmap_id %in% job$held_cells)
    test_pair <- if (job$regime == 'recipient') which(d$depmap_id %in% job$held_cells) else
      if (job$regime == 'partner') which(d$partner == job$held_partner & !is.na(d$partner)) else integer()
    if (job$regime == 'partner') {
      pi <- pi[!d$drugA[pi] %in% job$held_partner & !d$drugB[pi] %in% job$held_partner]
      held_refs <- unique(c(ref_a[test_pair], ref_b[test_pair]))
      pi <- pi[!ref_a[pi] %in% held_refs & !ref_b[pi] %in% held_refs]
      si <- si[!SINGLE$control_key[si] %in% held_refs]
    }
    out_single <- out_pair <- out_molecular <- out_oracle <- out_audit <- models <- list()
    fitting_families <- families
    for (family in fitting_families) {
      p <- prepare(si, pi, union(AU, job$held_cells), family)
      model <- fit_model(p, family, job$tag, check_gradient = job$tag == 'recipient_1')
      model$held_cells <- job$held_cells
      model$held_partner <- job$held_partner
      model$training_single_rows <- si
      model$training_pair_rows <- pi
      stopifnot(!any(model$training_cells %in% job$held_cells), !any(model$training_molecular_cells %in% job$held_cells),
        !AU %in% c(model$training_cells, model$training_molecular_cells))
      if (!is.na(job$held_partner)) stopifnot(!job$held_partner %in% c(d$drugA[pi], d$drugB[pi]))
      qsave(model, file.path(root, 'models', paste0(job$tag, '__', gsub(' ', '_', family), '.qs')))
      models[[family]] <- model
      if (job$regime == 'recipient') {
        test <- evaluation[depmap_id %in% job$held_cells]
        out <- test[, .(evaluation_id, depmap_id, drug, dose, observed_log2 = y)]
        out[, `:=`(family = family, fold = job$fold, predicted_log2 = predict_single(model, test))]
        out_single[[family]] <- out
        measured <- as.data.table(PREP$m)[depmap_id %in% job$held_cells]
        predicted <- predict_displacement(model, measured)
        truth <- sweep(as.matrix(measured[, P, with = FALSE]), 2L, model$context$response_sd, '/')
        mol <- measured[, .(depmap_id, drug, dose)]
        mol[, `:=`(family = family, fold = job$fold)]
        mol <- cbind(mol, as.data.table(predicted))
        setnames(mol, names(mol)[5L + seq_len(40L)], P)
        out_molecular[[family]] <- mol
        observed <- merge(measured[, .(depmap_id, drug, dose, molecular_row = .I)],
          SINGLE[assay == 'secondary', .(depmap_id, drug, dose, observed_log2 = y)],
          by = c('depmap_id', 'drug', 'dose'))
        observed[, `:=`(family = family, fold = job$fold,
          measured_displacement_prediction = decoder(truth[molecular_row, , drop = FALSE], model$decoder,
            rep(1L, .N))$output,
          predicted_displacement_prediction = decoder(predicted[molecular_row, , drop = FALSE], model$decoder,
            rep(1L, .N))$output)]
        out_oracle[[family]] <- observed
      }
      if (length(test_pair)) {
        out <- d[test_pair, .(row_id, depmap_id, pair, drugA, drugB, doseA, doseB, SCREENER,
          observed_log2 = y, observed_departure_log2 = I)]
        prediction <- predict_pairs(model, d[test_pair])
        out[, `:=`(family = family, regime = job$regime, excluded_partner = job$held_partner,
          predicted_log2 = prediction$total, predicted_departure_log2 = prediction$departure)]
        out_pair[[family]] <- out
      }
      out_audit[[family]] <- data.table(job = job$tag, family, convergence = model$optimization$convergence,
        function_evaluations = model$optimization$counts[1], maximum_gradient = model$optimization$maximum_gradient,
        training_functional_cells = length(model$training_cells), training_molecular_cells = length(model$training_molecular_cells))
    }
    result <- list(run_contract = run_contract, single = rbindlist(out_single), pair = rbindlist(out_pair),
      molecular = rbindlist(out_molecular), oracle = rbindlist(out_oracle), audit = rbindlist(out_audit),
      models = if (job$regime == 'final') models else NULL)
    qsave(result, done_path)
    rm(p, model)
    gc(verbose = FALSE)
  }
  single_predictions[[job$tag]] <- result$single
  pair_predictions[[job$tag]] <- result$pair
  molecular_predictions[[job$tag]] <- result$molecular
  oracle_predictions[[job$tag]] <- result$oracle
  audit[[job$tag]] <- result$audit
  write_table(rbindlist(single_predictions), 'held_single_predictions_frozen')
  write_table(rbindlist(pair_predictions), 'held_combination_predictions_frozen')
  write_table(rbindlist(molecular_predictions), 'held_molecular_predictions_frozen')
  write_table(rbindlist(oracle_predictions), 'held_measured_displacement_function_frozen')
  write_table(rbindlist(audit), 'fit_audit')
  if (job$regime == 'final') final_models <- result$models
}

status('Applying linked and frozen decoders to observed HGF vectors')
hgf <- fread(file.path(reference, 'data/hgf_transported_inputs.tsv'))
h <- as.matrix(hgf[, P, with = FALSE])
hgf_predictions <- list()
for (family in families) {
  model <- final_models[[family]]
  hz <- sweep(h, 2L, model$context$response_sd, '/')
  fh <- decoder(hz, model$decoder, rep(1L, nrow(hz)))$output
  for (dose in c(.15625, .625, 1, 2.5)) {
    w <- dose_weights(dose)
    fl <- combined <- numeric(nrow(hz))
    for (i in 1:3) {
      l <- predict_displacement(model, data.table(depmap_id = AU, drug = anchor, dose = native_doses[i]))
      l <- l[rep(1L, nrow(hz)), , drop = FALSE]
      fl <- fl + decoder(l, model$decoder, rep(1L, nrow(hz)))$output * w[1, i]
      combined <- combined + (decoder(l, model$decoder, rep(1L, nrow(hz)))$output + fh +
        as.numeric(pair_features(l, hz, model$decoder$u) %*% model$decoder$interaction)) * w[1, i]
    }
    hgf_predictions[[paste(family, dose)]] <- data.table(family, dose, hgf_context = hgf$context,
      lapatinib_log2 = fl, hgf_log2 = fh, combination_log2 = combined,
      lapatinib_abundance_pct = 100 * 2^fl, hgf_abundance_pct = 100 * 2^fh,
      combination_abundance_pct = 100 * 2^combined, interaction_log2 = combined - fl - fh)
  }
}
write_table(rbindlist(hgf_predictions), 'au565_hgf_predictions_frozen')
prediction_paths <- file.path(root, 'data', c('held_single_predictions_frozen.tsv',
  'held_combination_predictions_frozen.tsv', 'held_molecular_predictions_frozen.tsv',
  'held_measured_displacement_function_frozen.tsv', 'au565_hgf_predictions_frozen.tsv'))
write_table(data.table(path = prediction_paths, sha256 = vapply(prediction_paths, digest,
  character(1), file = TRUE, algo = 'sha256')), 'prediction_freeze_hashes')
qsave(list(run_contract = run_contract, families = families, selected_drugs = selected_drugs,
  functional_drugs = functional_drugs, native_doses = native_doses, decoder_rank = decoder_rank,
  baseline_rank = baseline_rank, molecular_ridge = molecular_ridge,
  molecular_weight = molecular_weight, within_function_weight = within_function_weight,
  functional_baseline_inputs = FALSE, gene_program_projection_frozen = TRUE,
  fixed_configuration = TRUE, retrospective = TRUE), file.path(root, 'models/settings.qs'))

status('Evaluating fixed-treatment correlations and molecular anchors')
single <- rbindlist(single_predictions)
pair <- rbindlist(pair_predictions)
single_metrics <- single[, .(cells = uniqueN(depmap_id), pearson_log2 = pearson(observed_log2, predicted_log2),
  observed_below_50pct = sum(observed_log2 <= -1), predicted_below_50pct = sum(predicted_log2 <= -1)), by = .(family, drug, dose)]
write_table(single_metrics, 'fixed_drug_dose_single_metrics')
write_table(single_metrics[drug == anchor & dose == 2.5], 'b1_441_cell_comparison')
stopifnot(all(single_metrics[drug == anchor & dose == 2.5]$cells == 441L))
contexts <- pair[, .(cells = uniqueN(depmap_id), pearson_abundance = pearson(observed_log2, predicted_log2),
  pearson_departure = pearson(observed_departure_log2, predicted_departure_log2)),
  by = .(family, regime, SCREENER, pair, doseA, doseB)]
write_table(contexts, 'fixed_pair_dose_metrics')
combo_summary <- contexts[, .(contexts = .N, median_pearson_abundance = median(pearson_abundance, na.rm = TRUE),
  median_pearson_departure = median(pearson_departure, na.rm = TRUE)), by = .(family, regime)]
write_table(combo_summary, 'fixed_pair_dose_summary')
oracle <- rbindlist(oracle_predictions)
oracle_metrics <- oracle[, .(cells = uniqueN(depmap_id),
  measured_displacement_pearson = pearson(observed_log2, measured_displacement_prediction),
  predicted_displacement_pearson = pearson(observed_log2, predicted_displacement_prediction)), by = .(family, drug, dose)]
write_table(oracle_metrics, 'measured_vs_predicted_displacement_function')
mol <- rbindlist(molecular_predictions)
truth <- as.data.table(PREP$m)[, c('depmap_id', 'drug', 'dose', P), with = FALSE]
mol_metrics <- list()
for (family in families) for (k in 1:3) {
  model <- qread(file.path(root, 'models', paste0('recipient_', k, '__', gsub(' ', '_', family), '.qs')))
  family_value <- family
  z <- mol[family == family_value & fold == k]
  true_z <- merge(z[, .(depmap_id, drug, dose)], truth, by = c('depmap_id', 'drug', 'dose'), sort = FALSE)
  match_index <- match(key(z$depmap_id, z$drug, z$dose), key(true_z$depmap_id, true_z$drug, true_z$dose))
    observed <- as.matrix(true_z[match_index, P, with = FALSE])
    predicted_physical <- sweep(as.matrix(z[, P, with = FALSE]), 2L, model$context$response_sd, '*')
  for (j in seq_along(P)) mol_metrics[[paste(family, k, j)]] <- data.table(family, fold = k,
    depmap_id = z$depmap_id, drug = z$drug, dose = z$dose, program = P[j],
    observed = observed[, j], predicted = predicted_physical[, j])
}
molecular_long <- rbindlist(mol_metrics)
molecular_metrics <- molecular_long[, .(cells = uniqueN(depmap_id), pearson = pearson(observed, predicted)),
  by = .(family, drug, dose, program)]
write_table(molecular_metrics, 'held_program_displacement_correlations')
mol_summary <- molecular_metrics[, .(contexts = .N, median_pearson = median(pearson, na.rm = TRUE)), by = .(family, drug)]
write_table(mol_summary, 'held_molecular_summary')

# All three measured PRISM endpoints are shown; folds remain visible to expose
# calibration instability rather than interpreting artificial bands as biology.
plot_drugs <- selected_drugs
plots <- list()
for (family_value in families) for (dose_value in c(.15625, .625, 2.5)) {
  z <- single[family == family_value & drug %in% plot_drugs & abs(dose - dose_value) < 1e-9]
  z[, fold_label := factor(fold)]
  label <- single_metrics[family == family_value & drug %in% plot_drugs & abs(dose - dose_value) < 1e-9,
    .(drug, label = sprintf('r = %.3f; n = %d', pearson_log2, cells))]
  p <- ggplot(z, aes(observed_log2, predicted_log2)) + geom_point(aes(color = fold_label), size = .6, alpha = .5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .45, color = 'black') +
    geom_text(data = label, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.05, vjust = 1.3, size = 2.6) + facet_wrap(~drug, ncol = 3, scales = 'free') +
    scale_color_manual(values = c('#2378a0', '#bd5a28', '#72419d'), name = 'Validation fold') +
    theme_bw(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom', strip.text = element_text(size = 8)) +
    labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
      title = sprintf('%s: joint program rules at %g µM; common nine-drug panel', family_value, dose_value))
  plot_name <- paste(gsub('[^A-Za-z0-9]+', '_', family_value), gsub('\\.', '_', format(dose_value)), sep = '_')
  plots[[plot_name]] <- p
  ggsave(file.path(root, 'figures', paste0(plot_name, 'uM.png')), p, width = 9, height = 9, dpi = 140)
}
pdf(file.path(root, 'figures/all_measured_doses.pdf'), width = 9, height = 9)
for (p in plots) print(p)
dev.off()
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
input_paths <- c(script_path, file.path(original, 'modular_pair_pilot/bank_operator.R'),
  file.path(original, 'modular_pair_pilot/inputs.rds'),
  file.path(previous, 'data/individual_training_outcomes.tsv'), file.path(previous, 'data/combination_training_outcomes.tsv'),
  file.path(root, 'data/lincs_program_observations.qs'))
write_table(data.table(path = input_paths, sha256 = vapply(input_paths, digest, character(1), file = TRUE, algo = 'sha256')),
  'source_manifest')
capture.output(cat('V83: joint rules across the common nine drugs; topotecan absent in local LINCS.\n',
  'Three molecular-supervision sources: Tahoe, LINCS, and both; same functional rows and exclusions.\n',
  'LINCS anchors partial observations of full programs under a fixed gene-loading-pattern approximation.\n',
  'No LINCS partial value directly enters the functional decoder as measured full displacement.\n',
  'Baseline occurs only upstream; eight RNA-derived baseline directions predict the fixed 40 programs.\n',
  'Measured drug/dose RNA means are fixed. RNA ridge and functional supervision update the same slope coefficients.\n',
  'Pairs use F(Da) + F(Db) + J(Da,Db), with shared symmetric J and J(D,0)=0; no combined RNA sum.\n',
  'Measured displacement labels constrain imputation; actual measured displacements also train the decoder.\n',
  'Fixed configuration; no selection on held labels. Whole-cell and whole-partner exclusions.\n',
  'Within-treatment function contributes80% of functional loss; absolute calibration uses the same decoder.\n',
  'Prior evaluation labels are retrospective; HGF/Wilson functional outcomes never fitted.\n'),
  print(single_metrics[drug == anchor & dose == 2.5]), print(combo_summary),
  print(mol_summary[drug == anchor]), print(rbindlist(audit)),
  file = file.path(root, 'RESULTS.txt'))
status('Complete')
