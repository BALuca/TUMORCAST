#!/usr/bin/env Rscript
# Learn f from singles only, freeze it, then learn g without modifying f.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(glmnet)
  library(kernlab)
  library(digest)
})
root <- '/work/results/program_pilot/69_sequential_single_combination_learning'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
single_path <- '/work/results/program_pilot/62_separate_response_interaction_bridge/data/single_predictions.tsv.gz'
baseline_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
script_path <- '/work/src/program_pilot/69_sequential_single_combination_learning.R'
f_lambda <- 10^seq(0, -5, length.out = 24L)
f_width <- c(.25, 1, 4)
centers_per_assay <- 256L
g_lambda <- 10^seq(2, -3, length.out = 28L)
g_alpha <- c(.1, .5, 1)
set.seed(69011L)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn = 1)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
}
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
options(error = function() {
  status('FAILED; see run log')
  traceback(6)
  quit(save = 'no', status = 1L)
})
f_model_path <- file.path(root, 'models/frozen_single_function.qs')
if (file.exists(f_model_path)) {
  # An in-place g revision must use the exact previously frozen f. Validate its
  # upstream files, rather than refitting singles or changing their predictions.
  status('Verifying and reusing frozen f; only g will be refitted')
  cached_manifest <- fread(file.path(root, 'data/input_manifest.tsv.gz'))
  current_hashes <- vapply(cached_manifest$path, digest, character(1), file = TRUE, algo = 'sha256')
  stopifnot(identical(unname(current_hashes), cached_manifest$sha256))
  frozen_inputs <- qread(file.path(root, 'data/frozen_inputs.qs'))
  f_hash <- digest(file = f_model_path, algo = 'sha256')
  stopifnot(identical(f_hash, frozen_inputs$f_sha256))
  f_profiles <- frozen_inputs$f_predictions
  bt_scores <- frozen_inputs$D[frozen_inputs$assay == 'BT20', , drop = FALSE]
  input_paths <- cached_manifest$path
  rm(frozen_inputs)
} else {
status('Preparing full Tahoe original40 scores and native BT20 measurements')
panel <- qread(panel_path)
bank <- qread(file.path(prior, 'data/representation_scores.qs'))
W <- qread(file.path(prior, 'models/original40_weights.qs'))
stopifnot(all(rownames(W) %in% panel$genes), nrow(W) == 8842L, ncol(W) == 40L)
labels <- copy(panel$sensitivity[is.finite(prism_abundance_pct) & prism_abundance_pct > 0])
tp_parts <- list()
input_paths <- character()
for (i in seq_along(panel$input_paths)) {
  id <- panel$cell_map$depmap_id[i]
  ids <- labels[depmap_id == id, observation_id]
  if (!length(ids)) next
  input_paths <- c(input_paths, panel$input_paths[i])
  obj <- qread(panel$input_paths[i])
  idx <- match(ids, rownames(obj$response))
  stopifnot(!anyNA(idx), all(rownames(W) %in% colnames(obj$response)))
  tp_parts[[id]] <- as.matrix(obj$response[idx, rownames(W), drop = FALSE] %*% W)
}
tp_scores <- do.call(rbind, tp_parts)
rownames(tp_scores) <- unlist(lapply(tp_parts, rownames), use.names = FALSE)
tp_scores <- tp_scores[labels$observation_id, , drop = FALSE]
bt_profiles <- bank$metadata[time_h == 3, profile_id]
bt_scores <- bank$original40[bt_profiles, , drop = FALSE]
stopifnot(identical(colnames(bt_scores), colnames(tp_scores)), all(is.finite(tp_scores)))
rm(obj, tp_parts)
singles <- fread(single_path)[representation == 'Original40' & time_h == 3,
  .(sample_id, source, drug, nominal_dose_uM, biological_replicate, profile_id, ratio)]
stopifnot(!anyDuplicated(paste(singles$source, singles$sample_id)), all(singles$ratio > 0))

# Baseline uses the same CCLE assay in every cell, including BT20 (ACH-000536).
# Collapse duplicated gene symbols, report coverage, and keep signed W unchanged.
header <- names(fread(baseline_path, nrows = 0L))
symbols <- sub(' \\(.*\\)$', '', header)
baseline_genes <- intersect(rownames(W), symbols)
selected <- which(symbols %chin% baseline_genes)
raw <- fread(baseline_path, select = c(header[1], header[selected]))
ids <- c(unique(labels$depmap_id), 'ACH-000536')
idx <- match(ids, raw[[1]])
stopifnot(!anyNA(idx), length(baseline_genes) > 8500L)
raw_matrix <- as.matrix(raw[idx, -1, with = FALSE])
collapsed <- rowsum(t(raw_matrix), symbols[selected], reorder = FALSE)
collapsed <- collapsed / as.numeric(table(symbols[selected])[rownames(collapsed)])
baseline <- as.matrix(t(collapsed[baseline_genes, , drop = FALSE]) %*% W[baseline_genes, ])
rownames(baseline) <- ids
stopifnot(all(is.finite(baseline)))
rm(raw, raw_matrix, collapsed)
write_table(data.table(program = colnames(W), total_genes = colSums(W != 0),
  Tahoe_genes = colSums(W != 0), baseline_genes = colSums(W[baseline_genes, ] != 0)), 'full_program_coverage')

# Available profiles are distinct from the repeated functional measurements.
profile_id <- c(labels$observation_id, bt_profiles)
D <- rbind(tp_scores, bt_scores)
H <- baseline[c(labels$depmap_id, rep('ACH-000536', length(bt_profiles))), , drop = FALSE]
assay <- c(rep('PRISM', nrow(tp_scores)), rep('BT20', nrow(bt_scores)))
rownames(D) <- rownames(H) <- profile_id
md <- rbind(labels[, .(sample_id = observation_id, source = 'Tahoe_PRISM', drug,
    cell = depmap_id, profile_id = observation_id, ratio = prism_abundance_pct / 100,
    dose = dose, biological_replicate = NA_integer_)],
  singles[, .(sample_id, source, drug, cell = 'ACH-000536', profile_id, ratio,
    dose = nominal_dose_uM, biological_replicate)])
md[, row := match(profile_id, rownames(D))]
stopifnot(!anyNA(md$row), !anyDuplicated(profile_id), all(md$ratio > 0))
tp_cells <- sort(unique(labels$depmap_id))
bt_drugs <- sort(unique(singles$drug))
tp_fold <- setNames(sample(rep(1:5, length.out = length(tp_cells))), tp_cells)
bt_fold <- setNames(sample(rep(1:5, length.out = length(bt_drugs))), bt_drugs)
md[, fold := fifelse(source == 'Tahoe_PRISM', tp_fold[cell], bt_fold[drug])]
stopifnot(!anyNA(md$fold), md[source != 'Tahoe_PRISM', uniqueN(fold), by = drug][, all(V1 == 1L)])
write_table(md, 'single_training_labels_and_folds')
source_mass <- c(Tahoe_PRISM = .5, 'Individual drug panel' = .25, 'Same-screen single controls' = .25)
single_weights <- function(m) {
  w <- numeric(nrow(m))
  for (ss in unique(m$source)) {
    ii <- which(m$source == ss)
    group <- paste(m$cell[ii], m$drug[ii])
    ww <- 1 / as.numeric(table(group)[group])
    w[ii] <- source_mass[ss] * ww / sum(ww)
  }
  w / sum(w)
}
target_f <- log2(md$ratio)

# A baseline kernel times a response kernel anchored at zero allows f(h,0)=0.
# Nyström centers come solely from fitting profiles, balanced across assays.
# The BT20 deviation uses the same nonlinear features with fourfold shrinkage.
conditional_kernel <- function(h, d, hc, dc, width) {
  kd <- as.matrix(kernelMatrix(rbfdot(sigma = width / (2 * ncol(d))), d, dc))
  zd <- exp(-width * rowSums(d^2) / (2 * ncol(d)))
  zc <- exp(-width * rowSums(dc^2) / (2 * ncol(dc)))
  anchored <- kd - outer(zd, rep(1, nrow(dc))) - outer(rep(1, nrow(d)), zc) + 1
  kh <- as.matrix(kernelMatrix(rbfdot(sigma = width / (2 * ncol(h))), h, hc))
  kh * anchored
}
make_basis <- function(rows, width) {
  fitting_profiles <- sort(unique(md$row[rows]))
  drms <- lapply(c('PRISM', 'BT20'), function(ss) sqrt(colMeans(D[fitting_profiles[assay[fitting_profiles] == ss], ]^2)))
  names(drms) <- c('PRISM', 'BT20')
  scaled_d <- D
  for (ss in names(drms)) {
    stopifnot(all(drms[[ss]] > 0))
    scaled_d[assay == ss, ] <- sweep(D[assay == ss, , drop = FALSE], 2, drms[[ss]], '/')
  }
  hb <- unique(H[fitting_profiles, , drop = FALSE])
  hmean <- colMeans(hb)
  hsd <- apply(hb, 2, sd)
  stopifnot(all(hsd > 0))
  scaled_h <- sweep(sweep(H, 2, hmean, '-'), 2, hsd, '/')
  set.seed(69012L)
  centers <- unlist(lapply(c('PRISM', 'BT20'), function(ss) {
    possible <- fitting_profiles[assay[fitting_profiles] == ss]
    sample(possible, min(centers_per_assay, length(possible)))
  }), use.names = FALSE)
  hc <- scaled_h[centers, , drop = FALSE]
  dc <- scaled_d[centers, , drop = FALSE]
  kc <- conditional_kernel(hc, dc, hc, dc, width)
  eig <- eigen((kc + t(kc)) / 2, symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-8
  stopifnot(min(eig$values) > -1e-7, sum(keep) > 20L)
  projection <- sweep(eig$vectors[, keep, drop = FALSE], 2, sqrt(eig$values[keep]), '/')
  shared <- conditional_kernel(scaled_h, scaled_d, hc, dc, width) %*% projection
  local <- shared * as.numeric(assay == 'BT20') / 2
  features <- cbind(shared, local)
  colnames(features) <- c(sprintf('shared_%03d', seq_len(ncol(shared))), sprintf('BT20_%03d', seq_len(ncol(shared))))
  zero <- conditional_kernel(scaled_h[1:3, , drop = FALSE], scaled_d[1:3, , drop = FALSE] * 0, hc, dc, width)
  stopifnot(max(abs(zero)) < 1e-10)
  list(features = features, response_rms = drms, baseline_mean = hmean, baseline_sd = hsd,
    center_h = hc, center_d = dc, projection = projection, centers = profile_id[centers], width = width)
}
fit_f_path <- function(x, y, w) {
  yscale <- sqrt(sum(w * y^2))
  fit <- glmnet(x, y / yscale, weights = w, alpha = 0, lambda = f_lambda,
    standardize = FALSE, intercept = FALSE, thresh = 1e-7, maxit = 30000L)
  stopifnot(fit$jerr == 0L)
  as.matrix(coef(fit, s = f_lambda))[-1, , drop = FALSE] * yscale
}
status('Selecting f using single-drug folds only: held Tahoe cells and held BT20 drugs')
f_oof <- vector('list', length(f_width))
f_losses <- matrix(0, length(f_lambda), length(f_width))
for (j in seq_along(f_width)) {
  pred <- matrix(NA_real_, nrow(md), length(f_lambda))
  for (fold in 1:5) {
    tr <- which(md$fold != fold)
    va <- which(md$fold == fold)
    basis <- make_basis(tr, f_width[j])
    beta <- fit_f_path(basis$features[md$row[tr], , drop = FALSE], target_f[tr], single_weights(md[tr]))
    pred[va, ] <- basis$features[md$row[va], , drop = FALSE] %*% beta
  }
  stopifnot(all(is.finite(pred)))
  f_losses[, j] <- colSums((pred - target_f)^2 * single_weights(md))
  f_oof[[j]] <- pred
  status(paste('f single-data CV completed for kernel width', f_width[j]))
}
choice <- which(f_losses == min(f_losses), arr.ind = TRUE)[1, ]
final_basis <- make_basis(seq_len(nrow(md)), f_width[choice[2]])
f_beta <- fit_f_path(final_basis$features[md$row, , drop = FALSE], target_f, single_weights(md))[, choice[1]]
f_profiles <- as.numeric(final_basis$features %*% f_beta)
names(f_profiles) <- profile_id
f_model_path <- file.path(root, 'models/frozen_single_function.qs')
qsave(list(beta = f_beta, basis = final_basis[setdiff(names(final_basis), 'features')],
  lambda = f_lambda[choice[1]], width = f_width[choice[2]], training_ids = md$sample_id,
  outcome_sources = unique(md$source), full_Tahoe_genes = rownames(W), BT20_measurement = 'Native partial L1000 projection',
  original_weight_sha256 = bank$original_weight_sha256), f_model_path)
f_hash <- digest(file = f_model_path, algo = 'sha256')
qsave(list(D = D, H = H, assay = assay, profile_ids = profile_id, baseline = baseline,
  f_predictions = f_profiles, f_sha256 = f_hash), file.path(root, 'data/frozen_inputs.qs'))
fs <- copy(md)
fs[, `:=`(predicted_fitted_log2 = f_profiles[profile_id],
  predicted_cv_log2 = f_oof[[choice[2]]][, choice[1]])]
write_table(fs, 'single_predictions')
f_grid <- CJ(width_index = seq_along(f_width), lambda_index = seq_along(f_lambda))
f_grid[, `:=`(width = f_width[width_index], lambda = f_lambda[lambda_index], internal_objective = as.vector(f_losses),
  selected = width_index == choice[2] & lambda_index == choice[1])]
write_table(f_grid, 'single_tuning')
write_table(data.table(profile_id = profile_id, predicted_f_log2 = f_profiles), 'frozen_single_profile_predictions')
status(paste('f frozen; now loading combination outcomes. Selected width', f_width[choice[2]], 'lambda', f_lambda[choice[1]]))
rm(f_oof, basis, final_basis, beta, pred)
}

# No operation below changes f. Combination labels first enter here.
obs <- fread(file.path(prior, 'data/matched_molecular_links.tsv.gz'))[time_h == 3 & molecular_supported == TRUE]
setorder(obs, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)
stopifnot(nrow(obs) == 496L, !anyDuplicated(obs$observation_id))
f_a <- unname(f_profiles[obs$profile_A])
f_b <- unname(f_profiles[obs$profile_B])
offset <- f_a + f_b
target_g <- log2(obs$ratio_AB) - offset
bt_rms <- sqrt(colMeans(bt_scores^2))
scaled_bt <- sweep(bt_scores, 2, bt_rms, '/')
smooth_basis <- function(x, f) {
  sat <- tanh(x)
  colnames(sat) <- paste0(colnames(x), '_sat')
  cbind(x, sat, learned_single_effect = f)
}
a <- smooth_basis(scaled_bt[obs$profile_A, , drop = FALSE], f_a)
b <- smooth_basis(scaled_bt[obs$profile_B, , drop = FALSE], f_b)
cross_features <- function(a, b) {
  ij <- which(upper.tri(matrix(0, ncol(a), ncol(a)), diag = TRUE), arr.ind = TRUE)
  x <- vapply(seq_len(nrow(ij)), function(k) {
    i <- ij[k, 1]
    j <- ij[k, 2]
    (a[, i] * b[, j] + a[, j] * b[, i]) / 2
  }, numeric(nrow(a)))
  colnames(x) <- paste(colnames(a)[ij[, 1]], colnames(a)[ij[, 2]], sep = '__')
  x
}
xg <- cross_features(a, b)
stopifnot(ncol(xg) == 3321L, max(abs(xg - cross_features(b, a))) == 0,
  max(abs(cross_features(a * 0, b))) == 0, all(is.finite(target_g)))
partners <- sort(unique(obs$partner))
balance <- function(m) {
  group <- paste(m$anchor, m$partner, m$biological_replicate)
  w <- 1 / as.numeric(table(group)[group])
  w / sum(w)
}
center_context <- function(x, m) {
  out <- x
  group <- interaction(m$anchor, m$nominal_A_uM, m$biological_replicate, drop = TRUE)
  for (idx in split(seq_len(nrow(m)), group)) {
    if (is.matrix(x)) out[idx, ] <- sweep(x[idx, , drop = FALSE], 2, colMeans(x[idx, , drop = FALSE]), '-')
    else out[idx] <- x[idx] - mean(x[idx])
  }
  out
}
objective <- function(y, p, m) {
  .5 * colSums((p - y)^2 * balance(m)) +
    .5 * colSums((center_context(p, m) - center_context(y, m))^2 * balance(m))
}
cache <- new.env(parent = emptyenv())
fit_g <- function(excluded, alpha) {
  key <- paste(paste(sort(excluded), collapse = '+'), alpha, sep = '|')
  if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
  tr <- which(!obs$partner %chin% excluded)
  w <- balance(obs[tr])
  rms <- sqrt(colSums(xg[tr, , drop = FALSE]^2 * w))
  keep <- is.finite(rms) & rms > 1e-10
  xx <- sweep(xg[tr, keep, drop = FALSE], 2, rms[keep], '/')
  yy <- target_g[tr]
  xx <- rbind(xx, center_context(xx, obs[tr]))
  yy <- c(yy, center_context(yy, obs[tr]))
  ww <- c(w * .5, w * .5)
  yscale <- sqrt(sum(ww * yy^2))
  fit <- glmnet(xx, yy / yscale, weights = ww, alpha = alpha, lambda = g_lambda,
    standardize = FALSE, intercept = FALSE, thresh = 1e-7, maxit = 30000L)
  stopifnot(fit$jerr == 0L)
  beta <- matrix(0, ncol(xg), length(g_lambda), dimnames = list(colnames(xg), NULL))
  beta[keep, ] <- sweep(as.matrix(coef(fit, s = g_lambda))[-1, , drop = FALSE] * yscale, 1, rms[keep], '/')
  ans <- list(coefficients = beta, predictions = xg %*% beta)
  assign(key, ans, cache)
  ans
}
predictions <- audits <- selections <- coefficients <- list()
for (held in partners) {
  status(paste('Learning g; all', held, 'combinations excluded'))
  training_partners <- setdiff(partners, held)
  # Select across independent intervention groups, not individual dose rows.
  # The zero correction is a real candidate. Within one SE of the best grouped
  # validation objective, choose the sparsest model, then strongest shrinkage.
  fold_losses <- array(NA_real_, c(length(training_partners), length(g_lambda), length(g_alpha)))
  fold_active <- fold_losses
  zero_losses <- numeric(length(training_partners))
  for (ii in seq_along(training_partners)) {
    va <- which(obs$partner == training_partners[ii])
    zero_losses[ii] <- objective(target_g[va], matrix(0, length(va), 1L), obs[va])
  }
  for (j in seq_along(g_alpha)) for (ii in seq_along(training_partners)) {
    inner <- training_partners[ii]
    candidate <- fit_g(c(held, inner), g_alpha[j])
    va <- which(obs$partner == inner)
    fold_losses[ii, , j] <- objective(target_g[va], candidate$predictions[va, , drop = FALSE], obs[va])
    fold_active[ii, , j] <- colSums(abs(candidate$coefficients) > 1e-8)
  }
  losses <- apply(fold_losses, c(2, 3), mean)
  standard_errors <- apply(fold_losses, c(2, 3), sd) / sqrt(length(training_partners))
  complexity <- apply(fold_active, c(2, 3), mean)
  grid <- CJ(alpha_index = seq_along(g_alpha), lambda_index = seq_along(g_lambda))
  grid[, `:=`(fold = held, alpha = g_alpha[alpha_index], lambda = g_lambda[lambda_index],
    internal_objective = as.vector(losses), objective_se = as.vector(standard_errors),
    mean_active_terms = as.vector(complexity), is_zero = FALSE)]
  grid <- rbind(grid, data.table(alpha_index = NA_integer_, lambda_index = NA_integer_, fold = held,
    alpha = NA_real_, lambda = Inf, internal_objective = mean(zero_losses),
    objective_se = sd(zero_losses) / sqrt(length(zero_losses)), mean_active_terms = 0, is_zero = TRUE))
  best_row <- which.min(grid$internal_objective)
  threshold <- grid$internal_objective[best_row] + grid$objective_se[best_row]
  grid[, eligible_one_se := internal_objective <= threshold]
  chosen <- grid[eligible_one_se == TRUE][order(mean_active_terms, -lambda, -alpha)][1]
  grid[, selected := (is_zero & chosen$is_zero) | (!is_zero & !chosen$is_zero &
    alpha_index == chosen$alpha_index & lambda_index == chosen$lambda_index)]
  grid[is.na(selected), selected := FALSE]
  stopifnot(sum(grid$selected) == 1L)
  if (chosen$is_zero) {
    g_beta <- setNames(rep(0, ncol(xg)), colnames(xg))
  } else {
    fit <- fit_g(held, chosen$alpha)
    g_beta <- fit$coefficients[, chosen$lambda_index]
  }
  gpred <- as.numeric(xg %*% g_beta)
  tr <- which(obs$partner != held)
  te <- which(obs$partner == held)
  ids <- function(x) unique(unlist(strsplit(x, ';', fixed = TRUE)))
  overlap <- intersect(ids(obs$source_sample_ids_AB[tr]), ids(obs$source_sample_ids_AB[te]))
  stopifnot(length(overlap) == 0L, identical(f_hash, digest(file = f_model_path, algo = 'sha256')))
  model_path <- file.path(root, 'models', paste0('g_exclude_', held, '.qs'))
  qsave(list(beta = g_beta, lambda = chosen$lambda, alpha = chosen$alpha,
    selection = 'sparsest within one SE across inner excluded partners, including g=0',
    one_se_threshold = threshold, frozen_f_sha256 = f_hash, response_rms = bt_rms,
    training_ids = obs$observation_id[tr]), model_path)
  replay <- as.numeric(xg %*% qread(model_path)$beta)
  stopifnot(max(abs(replay - gpred)) < 1e-10, all(is.finite(2^(offset + gpred))))
  out <- obs[, .(observation_id, condition_id, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)]
  out[, `:=`(fold = held, split = fifelse(partner == held, 'held', 'training'),
    predicted_A_log2 = f_a, predicted_B_log2 = f_b, predicted_g = gpred, predicted_log2 = offset + gpred,
    additive_log2 = offset, frozen_f_sha256 = f_hash)]
  write_table(out, paste0('predictions_exclude_', held))
  prediction_hash <- digest(out, algo = 'sha256')
  out[, `:=`(ratio_AB = obs$ratio_AB, ratio_A = obs$ratio_A, ratio_B = obs$ratio_B)]
  predictions[[held]] <- out
  audits[[held]] <- data.table(fold = held, source_overlap = length(overlap),
    replay_difference = max(abs(replay - gpred)), frozen_f_sha256 = f_hash,
    prediction_sha256 = prediction_hash, active_terms = sum(g_beta != 0),
    lambda = chosen$lambda, alpha = chosen$alpha, zero_selected = chosen$is_zero,
    training_observations = length(tr))
  fold_diagnostics <- rbindlist(lapply(seq_along(g_alpha), function(j) {
    rbindlist(lapply(seq_along(training_partners), function(ii) {
      data.table(outer_held = held, inner_held = training_partners[ii], alpha = g_alpha[j], lambda = g_lambda,
        internal_objective = fold_losses[ii, , j], active_terms = fold_active[ii, , j], zero_objective = zero_losses[ii])
    }))
  }))
  write_table(fold_diagnostics, paste0('inner_partner_selection_', held))
  selections[[held]] <- grid
  coefficients[[held]] <- data.table(fold = held, term = names(g_beta), coefficient = g_beta)
  write_table(rbindlist(predictions), 'predictions')
  write_table(rbindlist(audits), 'validation')
  write_table(rbindlist(selections), 'interaction_tuning')
  write_table(rbindlist(coefficients), 'interaction_coefficients')
}
inputs <- unique(c(panel_path, single_path, baseline_path, input_paths,
  file.path(prior, c('data/representation_scores.qs', 'models/original40_weights.qs', 'data/matched_molecular_links.tsv.gz'))))
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')), 'input_manifest')
write_table(data.table(path = script_path, sha256 = digest(file = script_path, algo = 'sha256')), 'script_manifest')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
writeLines(system('git -C /work rev-parse HEAD', intern = TRUE), file.path(root, 'git_revision.txt'))
status('f and all six excluded-partner g models complete; ready for figures')
