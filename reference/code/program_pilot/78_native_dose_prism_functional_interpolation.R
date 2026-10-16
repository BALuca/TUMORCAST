#!/usr/bin/env Rscript
# Align PRISM function to native Tahoe doses, then jointly fit singles and pairs.
# Off-grid function is interpolated after F; molecular states are not interpolated.
# Molecular operators remain trained on molecular data only; Wilson is not fitted.
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

# Fixed, small iteration. The molecular shortlist and five nonlinear program
# dimensions are inherited from75; no response-variation or synergy filtering.
root <- '/work/results/program_pilot/78_native_dose_prism_functional_interpolation'
prior <- '/work/results/program_pilot/75_hgf_pattern_partner_transfer'
reference <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query'
original <- file.path(reference, 'inputs/original')
script_path <- '/work/src/program_pilot/78_native_dose_prism_functional_interpolation.R'
prism_dir <- '/work/cache/tahoe_prism_full_chain/prism'
native_doses <- c(.05, .5, 5)
baseline_rank <- 4L
lambda_grid <- 10^seq(0, -3, length.out = 12L)
families <- c('Joint singles and combinations', 'Singles only')
seed <- 77012L
anchor <- 'Lapatinib ditosylate'
set.seed(seed)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
}
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
options(error = function() {
  status('FAILED; see run log')
  traceback(5)
  quit(save = 'no', status = 1L)
})
status('Preparing all10 drugs, matched individual controls, and all combinations')
setwd(original)
source('modular_pair_pilot/model.R')
settings75 <- qread(file.path(prior, 'models/settings.qs'))
selected_drugs <- c(anchor, settings75$partners$drug)
support <- settings75$programs
stopifnot(length(selected_drugs) == 10L, length(support) == 5L)
DRUGS <- DRUGS[DRUGS %in% selected_drugs]
stopifnot(setequal(unname(DRUGS), selected_drugs))
PREP$m <- PREP$m[PREP$m$drug %in% selected_drugs, ]
independent <- as.data.table(SINGLE)[drug %in% selected_drugs & assay %in% c('secondary', 'primary')]
independent[, `:=`(source = paste0('prism_', assay), control_key = paste(assay, depmap_id, drug, dose, sep = '|'))]

# Preserve measured PRISM endpoints for evaluation. Training targets are local
# log-dose/log2-response interpolations of full eight-dose curves within screens.
# Screen curves are then combined in the same replicate-weighted log2 endpoint
# as the archived model; this iteration does not change the screen-selection rule.
evaluation <- copy(independent[assay == 'secondary'])
evaluation[, evaluation_id := .I]
meta_path <- file.path(prism_dir, 'secondary-screen-replicate-treatment-info.csv')
raw_path <- file.path(prism_dir, 'secondary-screen-logfold-change.csv')
treatment <- fread(meta_path)[tolower(name) %in% names(DRUGS) & perturbation_type == 'experimental_treatment']
treatment[, drug := unname(DRUGS[tolower(name)])]
stopifnot(!anyDuplicated(treatment$column_name), setequal(treatment$drug, selected_drugs))
header <- names(fread(raw_path, nrows = 0))
stopifnot(all(treatment$column_name %in% header))
raw <- fread(raw_path, select = c(header[1], treatment$column_name))
setnames(raw, 1L, 'depmap_id')
raw <- raw[depmap_id %in% evaluation$depmap_id]
long <- melt(raw, id.vars = 'depmap_id', variable.name = 'column_name', value.name = 'y', variable.factor = FALSE)
long <- merge(long, treatment[, .(column_name, drug, dose, screen_id, broad_id)], by = 'column_name')
eligible <- unique(evaluation[, .(depmap_id, drug)])
long <- merge(long[is.finite(y)], eligible, by = c('depmap_id', 'drug'))
stopifnot(!AU %in% long$depmap_id)
points <- long[, .(y = mean(y), replicates = .N), by = .(depmap_id, drug, screen_id, broad_id, dose)]

# Confirm that raw labels reproduce the previous endpoints before changing the
# alignment. All original441 lapatinib labels stay measured, not interpolated.
measured <- long[, .(reconstructed_y = mean(y)), by = .(depmap_id, drug, dose)]
measured[, dose_key := sprintf('%.7g', dose)]
evaluation[, dose_key := sprintf('%.7g', dose)]
label_check <- merge(evaluation, measured[, .(depmap_id, drug, dose_key, reconstructed_y)],
  by = c('depmap_id', 'drug', 'dose_key'), all.x = TRUE)
stopifnot(nrow(label_check) == nrow(evaluation), !anyNA(label_check$reconstructed_y),
  max(abs(label_check$y - label_check$reconstructed_y)) < 1e-8)
write_table(points, 'prism_full_curve_points')
write_table(label_check, 'measured_endpoint_reconstruction')
native <- points[, {
  order <- order(dose)
  x <- log(dose[order])
  estimate <- approx(x, y[order], xout = log(native_doses), rule = 1)$y
  effective_replicates <- approx(x, replicates[order], xout = log(native_doses), rule = 1)$y
  .(dose = native_doses, y = estimate, weight = effective_replicates,
    measured_doses = .N, minimum_measured_dose = min(dose), maximum_measured_dose = max(dose))
}, by = .(depmap_id, drug, screen_id, broad_id)]
write_table(native, 'prism_native_dose_screen_estimates')
write_table(native[, .(targets = .N, unbracketed_targets = sum(!is.finite(y))), by = .(drug, dose)],
  'prism_interpolation_coverage')
native <- native[is.finite(y) & is.finite(weight)]
native <- native[, .(y = weighted.mean(y, weight)), by = .(depmap_id, drug, dose)]
native[, `:=`(assay = 'secondary', source = 'prism_secondary',
  control_key = paste('secondary', depmap_id, drug, dose, sep = '|'))]
independent <- rbindlist(list(native, independent[assay == 'primary']), use.names = TRUE)
rm(raw, long, measured)
gc(verbose = FALSE)
d <- copy(qread(file.path(prior, 'data/prepared_inputs.qs'))$metadata)
stopifnot(nrow(d) == 5640L, uniqueN(d$pair) == 42L, uniqueN(d$depmap_id) == 39L,
  !AU %in% c(d$depmap_id, independent$depmap_id, PREP$m$depmap_id))

# The archive carries matched constituent outcomes but lacks physical control IDs.
# Deduplicate identical outcomes; retain distinct outcomes at a reused setting,
# weighting the setting once. Purging uses the conservative screen/cell/drug/dose key.
control_key <- function(screen, cell, drug, dose) {
  paste(screen, cell, drug, sprintf('%.10g', dose), sep = '|')
}
ref_a <- control_key(d$SCREENER, d$depmap_id, d$drugA, d$doseA)
ref_b <- control_key(d$SCREENER, d$depmap_id, d$drugB, d$doseB)
constituents <- unique(rbind(
  d[, .(depmap_id, drug = drugA, dose = doseA, assay = paste0('nci_', SCREENER), y = yA,
    control_key = ref_a)],
  d[, .(depmap_id, drug = drugB, dose = doseB, assay = paste0('nci_', SCREENER), y = yB,
    control_key = ref_b)]))
constituents[, source := 'nci_single']
SINGLE <- rbindlist(list(independent, constituents), use.names = TRUE)
SINGLE[, single_id := .I]
SINGLE[, group := paste(assay, drug, sprintf('%.10g', dose), sep = '|')]
PAIR <- as.data.frame(d)
SINGLE <- as.data.frame(SINGLE)
stopifnot(all(is.finite(SINGLE$y)), all(SINGLE$dose >= .05 - 1e-10), all(SINGLE$dose <= 5 + 1e-10),
  setequal(unique(SINGLE$drug), selected_drugs), !AU %in% SINGLE$depmap_id)
write_table(as.data.table(SINGLE)[, .(observations = .N, cells = uniqueN(depmap_id), doses = uniqueN(dose)),
  by = .(drug, assay)], 'all_drug_single_coverage')
write_table(as.data.table(SINGLE), 'individual_training_outcomes')
write_table(d, 'combination_training_outcomes')

# Three recipient exclusions retain the exact441 B1 lapatinib recipients/folds.
# Other selected-drug recipients are assigned deterministically to those folds.
original_oof <- fread(file.path(reference, 'inputs/results/predictions.tsv'))[model == 'full' & dose == 2.5]
stopifnot(nrow(original_oof) == 441L, !anyDuplicated(original_oof$depmap_id))
all_ids <- sort(unique(c(SINGLE$depmap_id, d$depmap_id)))
recipient_fold <- foldmap(all_ids, 3L, seed)
recipient_fold[original_oof$depmap_id] <- original_oof$fold
stopifnot(!anyNA(recipient_fold))
write_table(data.table(depmap_id = names(recipient_fold), fold = unname(recipient_fold)), 'recipient_folds')

# Every molecular query is at one of the three native Tahoe concentrations.
# Off-grid predictions average evaluations of F, not molecular vectors before F.
entities <- unique(rbind(as.data.table(SINGLE)[, .(depmap_id, drug)],
  d[, .(depmap_id, drug = drugA)], d[, .(depmap_id, drug = drugB)], evaluation[, .(depmap_id, drug)]))
native_query <- entities[, .(dose = native_doses), by = .(depmap_id, drug)]
query_key <- function(cell, drug, dose) paste(cell, drug, sprintf('%.10g', dose), sep = '|')
native_key <- query_key(native_query$depmap_id, native_query$drug, native_query$dose)
query_indices <- function(cell, drug) {
  index <- vapply(native_doses, function(dose)
    match(query_key(cell, drug, dose), native_key), integer(length(cell)))
  stopifnot(!anyNA(index))
  index
}
# Log-dose interpolation weights are nonnegative, sum to one, and never extrapolate.
dose_weights <- function(dose) {
  stopifnot(all(dose >= min(native_doses) - 1e-10), all(dose <= max(native_doses) + 1e-10))
  dose <- pmax(min(native_doses), pmin(max(native_doses), dose))
  weights <- matrix(0, length(dose), length(native_doses))
  for (i in seq_along(dose)) {
    lo <- max(which(native_doses <= dose[i]))
    hi <- min(which(native_doses >= dose[i]))
    if (lo == hi) weights[i, lo] <- 1 else {
      fraction <- log(dose[i] / native_doses[lo]) / log(native_doses[hi] / native_doses[lo])
      weights[i, lo] <- 1 - fraction
      weights[i, hi] <- fraction
    }
  }
  stopifnot(max(abs(rowSums(weights) - 1)) < 1e-12)
  weights
}
single_index <- query_indices(SINGLE$depmap_id, SINGLE$drug)
pair_a_index <- query_indices(d$depmap_id, d$drugA)
pair_b_index <- query_indices(d$depmap_id, d$drugB)
evaluation_index <- query_indices(evaluation$depmap_id, evaluation$drug)
single_weights <- dose_weights(SINGLE$dose)
pair_a_weights <- dose_weights(d$doseA)
pair_b_weights <- dose_weights(d$doseB)
evaluation_weights <- dose_weights(evaluation$dose)
write_table(native_query, 'native_molecular_queries')

# Shared F has all40 linear displacement terms plus15 signed quadratic products.
# Four baseline-program PCs condition every term through tanh gates. Assay
# deviations affect linear terms only; the nonlinear coefficients are shared.
ij <- which(upper.tri(matrix(0, 5L, 5L), diag = TRUE), arr.ind = TRUE)
nonlinear_names <- paste(support[ij[, 1]], support[ij[, 2]], sep = '*')
basis_names <- c(P, nonlinear_names)
deviation_assays <- c('primary', 'nci_FF', 'nci_FG')
gate_names <- c('constant', paste0('baseline_PC', 1:baseline_rank))
common_names <- unlist(lapply(gate_names, function(g) paste(g, basis_names, sep = ':')))
linear_columns <- unlist(lapply(0:baseline_rank, function(g) g * length(basis_names) + seq_len(40L)))
feature_names <- c(common_names, unlist(lapply(deviation_assays, function(a)
  paste(a, common_names[linear_columns], sep = ':'))))
base_penalty <- c(rep(1, 40L), rep(4, 15L))
common_penalty <- unlist(lapply(0:baseline_rank, function(g) base_penalty * if (g == 0) 1 else 4))
penalty <- c(common_penalty, rep(common_penalty[linear_columns] * 8, length(deviation_assays)))
stopifnot(length(feature_names) == length(penalty))

# Row metadata explicitly distinguish absolute combination abundance from the
# matched departure contrast. Both constrain the very same F coefficients.
ns <- nrow(SINGLE)
np <- nrow(d)
row_metadata <- rbindlist(list(
  as.data.table(SINGLE)[, .(depmap_id, assay, source, group, control_key, kind = 'single',
    item_id = single_id, y)],
  d[, .(depmap_id, assay = paste0('nci_', SCREENER), source = 'combination_total', group = pair,
    control_key = paste0('combination:', row_id), kind = 'combination', item_id = row_id, y)],
  d[, .(depmap_id, assay = paste0('nci_', SCREENER), source = 'combination_departure', group = pair,
    control_key = paste0('departure:', row_id), kind = 'departure', item_id = row_id, y = I)]))
row_metadata[, observation_id := .I]
stopifnot(nrow(row_metadata) == ns + 2L * np)
source_totals <- c(prism_secondary = .50, prism_primary = .10, nci_single = .15,
  combination_total = .125, combination_departure = .125)

# Source totals and drug/dose or pair balancing prevent large panels and repeated
# control reuse from dominating the fit. This function returns normalized weights.
row_weights <- function(rows, family) {
  md <- row_metadata[rows]
  if (family == 'Singles only') stopifnot(all(md$kind == 'single'))
  w <- numeric(nrow(md))
  for (src in unique(md$source)) {
    ii <- which(md$source == src)
    z <- md[ii]
    unit_weight <- 1 / as.numeric(table(z$control_key)[z$control_key])
    group_weight <- tapply(unit_weight, z$group, sum)
    unit_weight <- unit_weight / as.numeric(group_weight[z$group]) / uniqueN(z$group)
    w[ii] <- source_totals[src] * unit_weight
  }
  stopifnot(all(is.finite(w)), all(w > 0))
  w / sum(w)
}

# Response scales use training molecular recipients only. Baseline PCs use only
# training functional recipients. No functional outcome selects a program or PC.
make_context <- function(training_ids, excluded) {
  raw <- getmatrix(PREP$raw, training_ids, P)
  pc <- prcomp(raw, center = TRUE, scale. = TRUE, rank. = baseline_rank)
  molecular <- as.data.table(PREP$m)[!depmap_id %in% excluded]
  centered <- molecular[, lapply(.SD, function(x) x - mean(x)), by = .(drug, dose), .SDcols = P]
  scale <- vapply(centered[, P, with = FALSE], sd, numeric(1))
  stopifnot(all(is.finite(scale)), all(scale > 0), all(pc$sdev[1:baseline_rank] > 0))
  list(pc = pc, response_scale = scale, training_ids = training_ids,
    excluded_molecular_recipients = excluded, baseline_rank = baseline_rank)
}
response_basis <- function(response, context) {
  z <- sweep(response, 2L, context$response_scale, '/')
  selected <- z[, support, drop = FALSE]
  nonlinear <- vapply(seq_len(nrow(ij)), function(k)
    selected[, ij[k, 1]] * selected[, ij[k, 2]], numeric(nrow(z)))
  result <- cbind(z, nonlinear)
  colnames(result) <- basis_names
  result
}
# Design depends on molecular responses, baseline programs and assay only.
# Drug identity and measured constituent function are never predictor features.
make_design <- function(b, baseline, assay, context) {
  score <- predict(context$pc, newdata = baseline)[, 1:baseline_rank, drop = FALSE]
  gates <- cbind(1, tanh(sweep(score, 2L, context$pc$sdev[1:baseline_rank], '/')))
  common <- do.call(cbind, lapply(seq_len(ncol(gates)), function(g) b * gates[, g]))
  result <- cbind(common, do.call(cbind, lapply(deviation_assays,
    function(a) common[, linear_columns, drop = FALSE] * as.numeric(assay == a))))
  colnames(result) <- feature_names
  Matrix(result, sparse = TRUE)
}

interpolate_basis <- function(basis, index, weights) {
  result <- basis[index[, 1], , drop = FALSE] * weights[, 1]
  for (k in 2:3) result <- result + basis[index[, k], , drop = FALSE] * weights[, k]
  result
}
# Algebraically equivalent to bilinear interpolation of F on the nine native
# dose-pair evaluations. Quadratic within-agent moments remain intact; cross
# moments factorize because the dose interpolation weights form a tensor product.
sum_basis <- function(a, b) {
  result <- a + b
  for (k in seq_len(nrow(ij))) {
    i <- match(support[ij[k, 1]], P)
    j <- match(support[ij[k, 2]], P)
    result[, 40L + k] <- result[, 40L + k] + a[, i] * b[, j] + a[, j] * b[, i]
  }
  result
}

# Descriptors are cross-fitted by molecular recipient using the preserved
# original operator-selection code, with outer/inner held cells removed.
prepare <- function(single_train, pair_train, excluded) {
  training_ids <- sort(unique(c(SINGLE$depmap_id[single_train], d$depmap_id[pair_train])))
  stopifnot(!any(training_ids %in% excluded), !AU %in% training_ids)
  context <- make_context(training_ids, excluded)
  desc <- build_descriptors(PREP, as.data.frame(native_query), as.data.frame(d[0]),
    excluded = excluded, crossfit = TRUE)
  basis <- response_basis(desc$X, context)
  ss <- getmatrix(PREP$raw, SINGLE$depmap_id, P)
  sp <- getmatrix(PREP$raw, d$depmap_id, P)
  bs <- interpolate_basis(basis, single_index, single_weights)
  ba <- interpolate_basis(basis, pair_a_index, pair_a_weights)
  bb <- interpolate_basis(basis, pair_b_index, pair_b_weights)
  combined <- sum_basis(ba, bb)
  xs <- make_design(bs, ss, SINGLE$assay, context)
  xa <- make_design(ba, sp, paste0('nci_', d$SCREENER), context)
  xb <- make_design(bb, sp, paste0('nci_', d$SCREENER), context)
  xsum <- make_design(combined, sp, paste0('nci_', d$SCREENER), context)
  xe <- make_design(interpolate_basis(basis, evaluation_index, evaluation_weights),
    getmatrix(PREP$raw, evaluation$depmap_id, P), evaluation$assay, context)
  departure <- xsum - xa - xb
  # On real conditions, verify the shortcut against all nine native state sums.
  direct <- matrix(0, 3L, 55L)
  for (i in 1:3) for (j in 1:3) {
    response <- desc$X[pair_a_index[1:3, i], , drop = FALSE] +
      desc$X[pair_b_index[1:3, j], , drop = FALSE]
    direct <- direct + response_basis(response, context) * pair_a_weights[1:3, i] * pair_b_weights[1:3, j]
  }
  stopifnot(max(abs(departure[, linear_columns])) < 1e-8,
    max(abs(direct - combined[1:3, ])) < 1e-8,
    max(abs(make_design(bs[1:3, , drop = FALSE] * 0, ss[1:3, , drop = FALSE],
      SINGLE$assay[1:3], context))) == 0,
    all(vapply(desc$operators, function(op) !any(op$training_ids %in% excluded), logical(1))))
  list(x = rbind(xs, xsum, departure), xa = xa, xb = xb, xsum = xsum, xe = xe,
    context = context, operators = desc$operators)
}

# Ridge paths update shared linear, nonlinear and baseline-dependent coefficients
# jointly. The identical architecture without combination rows is the comparator.
fit_path <- function(prepared, rows, family) {
  if (family == 'Singles only') rows <- rows[row_metadata$kind[rows] == 'single']
  w <- row_weights(rows, family)
  x <- prepared$x[rows, , drop = FALSE]
  rms <- sqrt(as.numeric(crossprod(w, x^2)))
  keep <- is.finite(rms) & rms > 1e-10
  x <- x[, keep, drop = FALSE] %*% Diagonal(x = 1 / rms[keep])
  y <- row_metadata$y[rows]
  yscale <- sqrt(sum(w * y^2))
  stopifnot(is.finite(yscale), yscale > 0)
  fit <- glmnet(x, y / yscale, alpha = 0, lambda = lambda_grid, weights = w,
    type.gaussian = 'covariance',
    intercept = FALSE, standardize = FALSE, penalty.factor = penalty[keep],
    thresh = 1e-7, maxit = 100000L)
  stopifnot(fit$jerr == 0L)
  coefficients <- matrix(0, ncol(prepared$x), length(lambda_grid), dimnames = list(feature_names, NULL))
  coefficients[keep, ] <- sweep(as.matrix(coef(fit, s = lambda_grid))[-1, , drop = FALSE] * yscale,
    1L, rms[keep], '/')
  list(coefficients = coefficients, training_rows = rows, rms = rms, yscale = yscale)
}
training_rows <- function(si, pi) c(si, ns + pi, ns + np + pi)

# Same-screen matched controls of held combination conditions are purged from
# both other combinations and constituent training rows. Independent PRISM
# single-drug measurements remain available in held-partner evaluations.
outer_jobs <- c(lapply(1:3, function(k) list(regime = 'recipient', tag = paste0('recipient_', k), fold = k,
  held_cells = names(recipient_fold)[recipient_fold == k], held_partner = NA_character_)),
  lapply(settings75$partners$drug, function(dr) list(regime = 'partner',
    tag = paste0('partner_', gsub('[^A-Za-z0-9]+', '_', dr)), fold = NA_integer_,
    held_cells = character(), held_partner = dr)),
  list(list(regime = 'final', tag = 'final', fold = NA_integer_, held_cells = character(), held_partner = NA_character_)))
tuning <- predictions_single <- predictions_pair <- audits <- molecular_audits <- list()
final_models <- list()

# Resume fully saved jobs from this iteration after a solver/process restart.
# Coefficients, exclusions, feature geometry and selected lambdas must match.
completed_path <- file.path(root, 'data/fit_audit.tsv')
completed <- if (file.exists(completed_path)) fread(completed_path) else data.table()
saved_single <- if (file.exists(file.path(root, 'data/held_single_predictions_frozen.tsv')))
  fread(file.path(root, 'data/held_single_predictions_frozen.tsv')) else data.table()
saved_pair <- if (file.exists(file.path(root, 'data/held_combination_predictions_frozen.tsv')))
  fread(file.path(root, 'data/held_combination_predictions_frozen.tsv')) else data.table()
saved_tuning <- if (file.exists(file.path(root, 'data/inner_parameter_selection.tsv')))
  fread(file.path(root, 'data/inner_parameter_selection.tsv')) else data.table()
saved_molecular <- if (file.exists(file.path(root, 'data/molecular_operator_audit.tsv')))
  fread(file.path(root, 'data/molecular_operator_audit.tsv')) else data.table()

for (job in outer_jobs) {
  if (nrow(completed) && sum(completed$job == job$tag) == length(families)) {
    for (family in families) {
      tag <- paste(job$tag, gsub(' ', '_', family), sep = '__')
      model <- qread(file.path(root, 'models', paste0(tag, '.qs')))
      stopifnot(identical(model$feature_names, feature_names), identical(model$support, support),
        identical(model$baseline_rank, baseline_rank), identical(model$held_cells, job$held_cells),
        identical(model$held_partner, job$held_partner), model$lambda %in% lambda_grid)
      completed_index <- which(completed$job == job$tag & completed$family == family)
      audits[[tag]] <- completed[completed_index]
      if (job$regime == 'recipient') {
        saved_single_index <- which(saved_single$fold == job$fold & saved_single$family == family)
        pair_index <- which(saved_pair$regime == 'recipient' & saved_pair$family == family &
          saved_pair$depmap_id %in% job$held_cells)
        predictions_single[[tag]] <- saved_single[saved_single_index]
        predictions_pair[[tag]] <- saved_pair[pair_index]
      } else if (job$regime == 'partner') {
        pair_index <- which(saved_pair$regime == 'partner' & saved_pair$family == family &
          saved_pair$excluded_partner == job$held_partner)
        predictions_pair[[tag]] <- saved_pair[pair_index]
      } else final_models[[family]] <- model
    }
    tuning_index <- which(saved_tuning$job == job$tag)
    molecular_index <- which(saved_molecular$job == job$tag)
    tuning[[job$tag]] <- saved_tuning[tuning_index]
    molecular_audits[[job$tag]] <- saved_molecular[molecular_index]
    status(paste(job$tag, '| restored complete saved job'))
    next
  }
  status(paste(job$tag, '| nested recipient selection of the shared bridge'))
  si <- which(!SINGLE$depmap_id %in% job$held_cells)
  pi <- which(!d$depmap_id %in% job$held_cells)
  test_pair <- if (job$regime == 'recipient') which(d$depmap_id %in% job$held_cells) else
    if (job$regime == 'partner') which(d$partner == job$held_partner & !is.na(d$partner)) else integer()
  if (job$regime == 'partner') {
    pi <- pi[!d$drugA[pi] %in% job$held_partner & !d$drugB[pi] %in% job$held_partner]
    held_refs <- unique(c(ref_a[test_pair], ref_b[test_pair]))
    pi <- pi[!ref_a[pi] %in% held_refs & !ref_b[pi] %in% held_refs]
    si <- si[!SINGLE$control_key[si] %in% held_refs]
    stopifnot(!any(job$held_partner %in% c(d$drugA[pi], d$drugB[pi])),
      !length(intersect(c(ref_a[pi], ref_b[pi], SINGLE$control_key[si]), held_refs)))
  }
  excluded <- union(AU, job$held_cells)
  ids <- sort(unique(c(SINGLE$depmap_id[si], d$depmap_id[pi])))
  inner_fold <- foldmap(ids, 3L, seed + 1L)
  loss <- matrix(0, length(lambda_grid), length(families), dimnames = list(NULL, families))
  for (k in 1:3) {
    held <- names(inner_fold)[inner_fold == k]
    inner_si <- si[!SINGLE$depmap_id[si] %in% held]
    inner_pi <- pi[!d$depmap_id[pi] %in% held]
    val_si <- si[SINGLE$depmap_id[si] %in% held]
    val_pi <- pi[d$depmap_id[pi] %in% held]
    prepared <- prepare(inner_si, inner_pi, union(excluded, held))
    for (family in families) {
      fitted <- fit_path(prepared, training_rows(inner_si, inner_pi), family)
      val <- if (family == 'Singles only') val_si else training_rows(val_si, val_pi)
      w <- row_weights(val, family)
      predicted <- as.matrix(prepared$x[val, , drop = FALSE] %*% fitted$coefficients)
      loss[, family] <- loss[, family] + colSums((predicted - row_metadata$y[val])^2 * w) / 3
    }
    rm(prepared, fitted)
    gc(verbose = FALSE)
  }
  selected <- vapply(families, function(family) which(loss[, family] <= min(loss[, family]) * 1.01)[1], integer(1))
  tuning[[job$tag]] <- rbindlist(lapply(families, function(family)
    data.table(job = job$tag, family = family, lambda = lambda_grid,
      internal_selection_objective = loss[, family], selected = seq_along(lambda_grid) == selected[family])))
  prepared <- prepare(si, pi, excluded)
  for (family in families) {
    fitted <- fit_path(prepared, training_rows(si, pi), family)
    beta <- fitted$coefficients[, selected[family]]
    model <- list(family = family, coefficients = beta, context = prepared$context,
      operators = prepared$operators, lambda = lambda_grid[selected[family]], support = support,
      native_doses = native_doses, input_rule = 'Interpolate functional evaluations on native molecular dose grid',
      training_single_rows = si, training_pair_rows = if (family == 'Singles only') integer() else pi,
      regime = job$regime, held_cells = job$held_cells, held_partner = job$held_partner,
      feature_names = feature_names, baseline_rank = baseline_rank)
    tag <- paste(job$tag, gsub(' ', '_', family), sep = '__')
    qsave(model, file.path(root, 'models', paste0(tag, '.qs')))
    # Fresh serialized coefficients must reproduce a real prediction vector.
    replay <- qread(file.path(root, 'models', paste0(tag, '.qs')))
    stopifnot(max(abs(as.numeric(prepared$x[1:5, ] %*% beta) -
      as.numeric(prepared$x[1:5, ] %*% replay$coefficients))) < 1e-10)
    audits[[tag]] <- data.table(job = job$tag, family = family, lambda = model$lambda,
      training_single_rows = length(si), training_pair_rows = length(model$training_pair_rows),
      held_combination_conditions = length(test_pair), held_cells = length(job$held_cells),
      nonlinear_coefficient_norm = sqrt(sum(beta[setdiff(seq_along(common_names), linear_columns)]^2)))
    if (job$regime == 'recipient') {
      test <- which(evaluation$depmap_id %in% job$held_cells)
      out <- evaluation[test, .(single_id = evaluation_id, depmap_id, drug, dose, observed_log2 = y)]
      out[, `:=`(fold = job$fold, family = family,
        predicted_log2 = as.numeric(prepared$xe[test, , drop = FALSE] %*% beta))]
      predictions_single[[tag]] <- out
    }
    if (length(test_pair)) {
      out <- d[test_pair, .(row_id, depmap_id, CELLNAME, pair, drugA, drugB, doseA, doseB, SCREENER,
        partner, observed_log2 = y, observed_departure_log2 = I, measured_singleA_log2 = yA, measured_singleB_log2 = yB)]
      out[, `:=`(regime = job$regime, excluded_partner = job$held_partner, family = family,
        predicted_log2 = as.numeric(prepared$xsum[test_pair, , drop = FALSE] %*% beta),
        predicted_singleA_log2 = as.numeric(prepared$xa[test_pair, , drop = FALSE] %*% beta),
        predicted_singleB_log2 = as.numeric(prepared$xb[test_pair, , drop = FALSE] %*% beta))]
      out[, predicted_departure_log2 := predicted_log2 - predicted_singleA_log2 - predicted_singleB_log2]
      predictions_pair[[tag]] <- out
    }
    if (job$regime == 'final') final_models[[family]] <- model
  }
  molecular_audits[[job$tag]] <- rbindlist(lapply(names(prepared$operators), function(dr) {
    op <- prepared$operators[[dr]]
    data.table(job = job$tag, drug = dr, molecular_operator = op$kind, lambda = op$lambda,
      training_molecular_cells = length(op$training_ids), shared_held_cells = sum(op$training_ids %in% excluded))
  }))
  write_table(rbindlist(tuning), 'inner_parameter_selection')
  write_table(rbindlist(audits), 'fit_audit')
  write_table(rbindlist(molecular_audits), 'molecular_operator_audit')
  if (length(predictions_single)) write_table(rbindlist(predictions_single), 'held_single_predictions_frozen')
  if (length(predictions_pair)) write_table(rbindlist(predictions_pair), 'held_combination_predictions_frozen')
  rm(prepared, fitted)
  gc(verbose = FALSE)
}

status('Querying HGF through the same jointly trained F, with untreated control as100%')
hgf <- fread(file.path(reference, 'data/hgf_transported_inputs.tsv'))
stopifnot(nrow(hgf) == 13L, sum(hgf$context == 'Equal-context mean') == 1L)
h <- as.matrix(hgf[, P, with = FALSE])
b <- getmatrix(PREP$raw, AU, P)[rep(1L, nrow(h)), , drop = FALSE]
queries <- list()
for (family in families) {
  model <- final_models[[family]]
  for (dose in c(1, 2.5)) {
    weights <- dose_weights(dose)
    lap_basis <- matrix(0, nrow(h), 55L)
    for (i in 1:3) {
      l <- uppred(PREP, model$operators[[anchor]], AU, native_doses[i])
      l <- l[rep(1L, nrow(h)), , drop = FALSE]
      lap_basis <- lap_basis + response_basis(l, model$context) * weights[1, i]
    }
    hgf_basis <- response_basis(h, model$context)
    fl <- as.numeric(make_design(lap_basis, b, rep('secondary', nrow(h)), model$context) %*% model$coefficients)
    fh <- as.numeric(make_design(hgf_basis, b, rep('secondary', nrow(h)), model$context) %*% model$coefficients)
    flh <- as.numeric(make_design(sum_basis(lap_basis, hgf_basis), b, rep('secondary', nrow(h)), model$context) %*% model$coefficients)
    stopifnot(max(abs(make_design(lap_basis * 0, b, rep('secondary', nrow(h)), model$context))) == 0,
      max(abs(sum_basis(lap_basis, hgf_basis) - sum_basis(hgf_basis, lap_basis))) < 1e-10)
    queries[[paste(family, dose)]] <- data.table(family = family, dose = dose, hgf_context = hgf$context,
      lapatinib_log2 = fl, hgf_log2 = fh, combination_log2 = flh,
      lapatinib_abundance_pct = 100 * 2^fl, hgf_abundance_pct = 100 * 2^fh,
      combination_abundance_pct = 100 * 2^flh,
      interaction_log2 = flh - fl - fh, absolute_hgf_effect_pp = 100 * (2^flh - 2^fl))
  }
  write_table(data.table(feature = names(model$coefficients), coefficient = model$coefficients),
    paste0('final_coefficients_', gsub(' ', '_', family)))
}
query <- rbindlist(queries)
write_table(query, 'au565_hgf_predictions_frozen')
frozen_paths <- file.path(root, 'data', c('held_single_predictions_frozen.tsv',
  'held_combination_predictions_frozen.tsv', 'au565_hgf_predictions_frozen.tsv'))
write_table(data.table(path = frozen_paths,
  sha256 = vapply(frozen_paths, digest, character(1), file = TRUE, algo = 'sha256')), 'prediction_freeze_hashes')
qsave(list(selected_drugs = selected_drugs, nonlinear_programs = support, baseline_rank = baseline_rank,
  lambda_grid = lambda_grid, source_totals = source_totals, source_assay_deviations = deviation_assays,
  native_doses = native_doses,
  functional_rule = 'Evaluate F on native molecular dose grid; interpolate log2 function; pair departure is difference of the same F evaluations',
  phenotype = 'log2 abundance relative to untreated control', molecular_fitting = 'Molecular outcomes only; recipient cross-fitting',
  family = families, seed = seed), file.path(root, 'models/settings.qs'))

# Scientific evaluation begins only after predictions are frozen. Pearson axes
# distinguish fixed drug/dose across cells from mixtures of drugs/doses/pairs.
pearson <- function(x, y) if (length(x) >= 4L && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
single <- rbindlist(predictions_single)
stopifnot(!anyDuplicated(single[, .(family, single_id)]))
single_metrics <- single[, .(cells = uniqueN(depmap_id), pearson_log2 = pearson(observed_log2, predicted_log2),
  observed_below_50pct = sum(observed_log2 <= -1), predicted_below_50pct = sum(predicted_log2 <= -1)),
  by = .(family, drug, dose)]
write_table(single_metrics, 'fixed_drug_dose_single_metrics')
b1 <- single_metrics[drug == anchor & dose == 2.5]
stopifnot(all(b1$cells == 441L))
write_table(b1, 'b1_441_cell_comparison')
pair <- rbindlist(predictions_pair)
pair[, `:=`(observed_departure_pp = 100 * (2^observed_log2 - 2^(measured_singleA_log2 + measured_singleB_log2)),
  predicted_departure_pp = 100 * (2^predicted_log2 - 2^(predicted_singleA_log2 + predicted_singleB_log2)))]
write_table(pair, 'held_combination_predictions_evaluated')
pooled <- pair[, .(conditions = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
  pearson_abundance_log2 = pearson(observed_log2, predicted_log2),
  pearson_predicted_independent = pearson(observed_log2, predicted_singleA_log2 + predicted_singleB_log2),
  pearson_measured_independent = pearson(observed_log2, measured_singleA_log2 + measured_singleB_log2),
  pearson_departure_log2 = pearson(observed_departure_log2, predicted_departure_log2),
  pearson_departure_pp = pearson(observed_departure_pp, predicted_departure_pp)), by = .(family, regime)]
write_table(pooled, 'pooled_combination_metrics')
within <- pair[, .(cells = uniqueN(depmap_id),
  pearson_abundance_log2 = pearson(observed_log2, predicted_log2),
  pearson_departure_log2 = pearson(observed_departure_log2, predicted_departure_log2)),
  by = .(family, regime, pair, doseA, doseB, SCREENER)]
write_table(within, 'fixed_pair_dose_across_cell_metrics')
within_summary <- within[cells >= 25, .(contexts = .N,
  finite_abundance = sum(is.finite(pearson_abundance_log2)), median_pearson_abundance = median(pearson_abundance_log2, na.rm = TRUE),
  finite_departure = sum(is.finite(pearson_departure_log2)), median_pearson_departure = median(pearson_departure_log2, na.rm = TRUE)),
  by = .(family, regime)]
write_table(within_summary, 'fixed_pair_dose_summary')
write_table(pair[regime == 'partner', .(conditions = .N, cells = uniqueN(depmap_id),
  pearson_abundance_log2 = pearson(observed_log2, predicted_log2),
  pearson_departure_log2 = pearson(observed_departure_log2, predicted_departure_log2)),
  by = .(family, excluded_partner)], 'held_partner_metrics')
write_table(pair[, .(conditions = .N, cells = uniqueN(depmap_id),
  pearson_abundance_log2 = pearson(observed_log2, predicted_log2),
  pearson_departure_log2 = pearson(observed_departure_log2, predicted_departure_log2)),
  by = .(family, regime, pair, drugA, drugB)], 'all_pair_ground_truth_correlations')

# All-drug across-cell comparisons at the same measured dose, including the
# previous architecture with molecular rather than functional interpolation.
previous <- fread('/work/results/program_pilot/77_joint_single_combination_bridge/data/fixed_drug_dose_single_metrics.tsv')
previous <- previous[, .(family, drug, dose, previous_pearson = pearson_log2)]
all_drug_comparison <- merge(single_metrics, previous, by = c('family', 'drug', 'dose'), all.x = TRUE)
write_table(all_drug_comparison, 'all_drug_measured_dose_comparison')

plot_all <- single[family == families[1] & dose == 2.5]
all_labels <- single_metrics[family == families[1] & dose == 2.5]
all_labels[, panel := paste0(sub(' \\(.*', '', drug), '\nn=', cells, '; r=', sprintf('%.2f', pearson_log2))]
plot_all <- merge(plot_all, all_labels[, .(drug, panel)], by = 'drug')
gp <- ggplot(plot_all, aes(observed_log2, predicted_log2)) +
  geom_point(size = .6, alpha = .4, color = '#2166AC') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .45) +
  facet_wrap(~panel, ncol = 3, scales = 'free') + theme_classic(base_size = 9) +
  theme(aspect.ratio = 1, strip.background = element_blank(), strip.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0, size = 8)) +
  labs(x = 'Measured PRISM log2 abundance', y = 'Held-out predicted log2 abundance',
    title = 'All10 drugs at2.5µM: correlations across cell lines',
    caption = 'Every point is one excluded recipient. Observed outcomes are measured at2.5µM.\nTraining uses native Tahoe doses and interpolated PRISM function; evaluation uses functional interpolation.')
ggsave(file.path(root, 'figures/all_drugs_2_5uM_across_cells.pdf'), gp, width = 8.2, height = 11.2, device = cairo_pdf)
ggsave(file.path(root, 'figures/all_drugs_2_5uM_across_cells.png'), gp, width = 8.2, height = 11.2, dpi = 150)

# Ground truth is measured ALMANAC endpoint abundance and departure from matched
# measured singles. Mixture correlations are kept separate from fixed-dose metrics.
pg <- pair[family == families[1]]
cg <- rbind(pg[, .(regime, observed = observed_log2, predicted = predicted_log2, endpoint = 'Total abundance')],
  pg[, .(regime, observed = observed_departure_log2, predicted = predicted_departure_log2, endpoint = 'Departure from separate effects')])
cg[, panel := paste(fifelse(regime == 'recipient', 'Excluded cell lines', 'Excluded partner combinations'), endpoint, sep = '\n')]
cg_labels <- cg[, .(label = sprintf('n=%d; r=%.3f', .N, pearson(observed, predicted))), by = panel]
cp <- ggplot(cg, aes(observed, predicted)) +
  geom_point(size = .45, alpha = .25, color = '#2166AC') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
  geom_text(data = cg_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.2, size = 3) +
  facet_wrap(~panel, ncol = 2, scales = 'free') + theme_classic(base_size = 9) +
  theme(aspect.ratio = 1, strip.background = element_blank(), strip.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0, size = 8)) +
  labs(x = 'Measured ALMANAC outcome (log2)', y = 'Held-out prediction (log2)',
    title = 'Combination predictions compared with ground truth',
    caption = 'Each panel pools different pairs and doses. Fixed-pair/dose across-cell correlations are reported separately.\nBoth individual effects and combination effects constrain the same F; there is no independent combination decoder.')
ggsave(file.path(root, 'figures/combinations_vs_ground_truth.pdf'), cp, width = 6.7, height = 7.1, device = cairo_pdf)
ggsave(file.path(root, 'figures/combinations_vs_ground_truth.png'), cp, width = 6.7, height = 7.1, dpi = 160)

# Two readable views. Each AU565 treatment has one final prediction. Both bars
# use untreated control as100%; no fit-specific points or changed denominators.
plot_single <- single[drug == anchor & dose == 2.5]
annotations <- b1[, .(family, label = sprintf('441 cells; r = %.3f', pearson_log2))]
sp <- ggplot(plot_single, aes(observed_log2, predicted_log2)) +
  geom_point(size = .8, alpha = .45, color = '#2166AC') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
  geom_text(data = annotations, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.3, size = 3) + facet_wrap(~family, scales = 'free') +
  theme_classic(base_size = 10) + theme(aspect.ratio = 1, strip.background = element_blank()) +
  labs(x = 'Observed log2 abundance', y = 'Predicted log2 abundance', title = 'Lapatinib 2.5 µM: held-out recipients')
ggsave(file.path(root, 'figures/b1_441_cell_joint_bridge.pdf'), sp, width = 6.8, height = 3.8, device = cairo_pdf)
ggsave(file.path(root, 'figures/b1_441_cell_joint_bridge.png'), sp, width = 6.8, height = 3.8, dpi = 160)

# Published Wilson values are opened after model fitting, selection and freezing.
wilson_path <- '/work/almanac_combination_training_2026-10-03/results/wilson_observed_source.tsv'
observed <- fread(wilson_path)[series == 'Wilson et al.' & condition %in% c('Lapatinib', 'Lapatinib + HGF')]
mean_query <- query[hgf_context == 'Equal-context mean']
auplot <- rbind(mean_query[, .(family, dose, condition = 'Lapatinib', abundance = lapatinib_abundance_pct)],
  mean_query[, .(family, dose, condition = 'Lapatinib + HGF', abundance = combination_abundance_pct)])
auplot[, series := 'Model prediction']
wp <- observed[, .(family = families[1], dose = 1, condition, abundance, sd, series = 'Wilson measurement')]
write_table(merge(auplot[dose == 1], observed[, .(condition, observed_abundance = abundance, observed_sd = sd)],
  by = 'condition'), 'wilson_1uM_comparison')
ap <- ggplot(auplot[family == families[1] & dose == 1], aes(condition, abundance, group = series, color = series)) +
  geom_hline(yintercept = 100, color = 'grey65', linetype = 2, linewidth = .35) +
  geom_line(linewidth = .75) + geom_point(size = 2.8) +
  geom_errorbar(data = wp, aes(ymin = abundance - sd, ymax = abundance + sd), width = .1, linewidth = .6) +
  geom_point(data = wp, shape = 18, size = 3) +
  geom_text(aes(label = sprintf('%.1f%%', abundance)), vjust = -1, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c('Model prediction' = '#2166AC', 'Wilson measurement' = '#B2182B')) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .14))) +
  theme_classic(base_size = 11) + theme(legend.position = 'bottom', plot.caption = element_text(hjust = 0, size = 8)) +
  labs(x = NULL, y = 'Cell abundance (% untreated control)', color = NULL, title = 'AU565: lapatinib 1 µM',
    caption = 'One final joint-model prediction per condition. 100%=untreated control for both conditions.\nRed: Wilson’s published mean±SD. HGF uses the mean of12 transported molecular profiles.\nHGF function was not fitted; HGF dose and assay-duration transfer remain uncalibrated.')
ggsave(file.path(root, 'figures/au565_joint_bridge_wilson.pdf'), ap, width = 5.7, height = 4.6, device = cairo_pdf)
ggsave(file.path(root, 'figures/au565_joint_bridge_wilson.png'), ap, width = 5.7, height = 4.6, dpi = 160)

inputs <- c(file.path(original, 'modular_pair_pilot/inputs.rds'), file.path(prior, 'models/settings.qs'),
  file.path(prior, 'data/prepared_inputs.qs'), file.path(reference, 'data/hgf_transported_inputs.tsv'),
  meta_path, raw_path, script_path)
write_table(data.table(path = inputs,
  sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')), 'source_manifest')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
writeLines(c('Version78: native-dose molecular inputs and interpolated functional targets; same architecture as77.',
  'Full PRISM eight-dose curves are interpolated locally in log-dose/log2 abundance within each screen; no extrapolated targets.',
  'Screen estimates use replicate-count weighting to preserve the archived screen-mixture endpoint. This is not provider-preferred screen selection.',
  'Measured original PRISM labels were independently reconstructed from raw data before interpolation.',
  'No variation/synergy filtering. All441 original lapatinib recipients and original outer folds retained.',
  'Single predictions interpolate evaluations of F at native molecular doses. Pair predictions bilinearly interpolate F at nine native-dose state sums.',
  'All40 linear programs;15 quadratic terms in the five molecularly nominated HGF dimensions; four training-only baseline PCs.',
  'Shared nonlinear coefficients; heavily penalized linear assay adjustments for primary PRISM and each ALMANAC lab.',
  'Drug identity and measured function are not predictor features. Molecular operators use molecular labels only and recipient cross-fitting.',
  'Both absolute combination outcomes and matched departures constrain F. Their labels are correlated and not independent experiments.',
  'Whole receiving-cell exclusions for singles/combination prediction; whole-partner exclusions and conservative shared-control purging for combination transfer.',
  'Nested three-recipient-fold ridge selection; same architecture trained on singles alone is the comparator.',
  'Wilson/HGF functional outcomes not used in fitting or selection. Historical shortlist/program bank and revealed outcomes make this retrospective development.',
  'HGF-alone and combined output use the same final F. All displayed AU565 values have untreated control as denominator.',
  capture.output(print(b1)), capture.output(print(pooled)), capture.output(print(within_summary)),
  capture.output(print(mean_query))
), file.path(root, 'RESULTS.txt'))
print(b1)
print(pooled)
print(within_summary)
print(mean_query)
status('Complete')
