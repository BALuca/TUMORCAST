#!/usr/bin/env Rscript
# Retrospective MEK+PI3K family exclusion: can separate RNA responses predict
# additional combination inhibition beyond observed constituent activities?
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(Matrix)
  library(digest)
})
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/60_defensible_response_composition'
set.seed(6011)
lambda_grid <- 10^seq(3, -2, length.out = 24L)
n_folds <- 4L
gene_components <- 12L
write_table <- function(x, file) fwrite(x, file.path(root, 'data', file), sep = '\t', na = 'NA')

links_path <- file.path(root, 'data/transfer_molecular_links.tsv.gz')
molecular_path <- file.path(root, 'data/transfer_molecular.qs')
eligibility_path <- file.path(root, 'data/representation_eligibility.tsv')
stopifnot(file.exists(links_path), file.exists(molecular_path), file.exists(eligibility_path))
all_obs <- fread(links_path)
obs <- all_obs[molecular_supported == TRUE]
setorder(obs, observation_id)
stopifnot(!anyDuplicated(obs$observation_id), all(is.finite(obs$viability_A)),
  all(is.finite(obs$viability_B)), all(is.finite(obs$viability_AB)))
obs[, focal := pair == 'pictilisib + trametinib']
train <- which(!obs$heldout)
test <- which(obs$heldout)
stopifnot(length(test) > 0L, length(train) > 100L,
  !any(obs$pair[train] %in% obs$pair[test]), all(obs[focal == TRUE, heldout]))

# Outcomes are separated from predictor construction. Public retrospective
# exclusions are not claimed as prospective blindness.
y_train <- obs$excess_inhibition[train]
obs[, c('viability_AB', 'excess_inhibition') := NULL]
write_table(obs[, .(observation_id, cell, pair, heldout, focal, A_matching_status,
  B_matching_status, dose_A_uM, dose_B_uM)], 'transfer_analysis_manifest.tsv.gz')
pairs <- sort(unique(obs$pair[train]))
stopifnot(length(pairs) >= 12L)
folds <- data.table(pair = sample(pairs), fold = rep(seq_len(n_folds), length.out = length(pairs)))
fold_id <- folds$fold[match(obs$pair[train], folds$pair)]
write_table(folds, 'transfer_inner_pair_folds.tsv')
obs[, fitting_weight := 0]
wt <- obs[train, .(observation_id, cell, pair)]
wt[, pair_rows := .N, by = .(cell, pair)]
wt[, cell_pairs := uniqueN(pair), by = cell]
wt[, weight := 1 / (pair_rows * cell_pairs)]
weights_train <- wt$weight / mean(wt$weight)
obs$fitting_weight[train] <- weights_train

# The same measured single-agent and dose nuisance features enter every fit.
# Exchange symmetry is explicit, so reversing drug A/B cannot change a prediction.
va <- obs$viability_A / 100
vb <- obs$viability_B / 100
da <- log10(obs$dose_A_uM)
db <- log10(obs$dose_B_uM)
single_x <- cbind(single_min = pmin(va, vb), single_max = pmax(va, vb),
  single_product = va * vb, single_min_sq = pmin(va, vb)^2,
  single_max_sq = pmax(va, vb)^2, logdose_min = pmin(da, db), logdose_max = pmax(da, db))
cell_x <- model.matrix(~cell - 1, obs)
molecular <- qread(molecular_path)
stopifnot(all(is.finite(molecular$displacement)), all(is.finite(molecular$baseline)))
baseline_pca <- prcomp(t(molecular$baseline), center = TRUE, scale. = FALSE, rank. = 2L)
baseline_scores <- baseline_pca$x[, 1:2, drop = FALSE]
baseline_x <- baseline_scores[match(obs$cell, rownames(baseline_scores)), , drop = FALSE]
colnames(baseline_x) <- paste0('baseline_PC', seq_len(ncol(baseline_x)))
context_x <- cbind(single_x, cell_x, baseline_x,
  baseline_PC1_single_min = baseline_x[, 1] * pmin(va, vb),
  baseline_PC1_single_max = baseline_x[, 1] * pmax(va, vb),
  baseline_PC2_single_min = baseline_x[, 2] * pmin(va, vb),
  baseline_PC2_single_max = baseline_x[, 2] * pmax(va, vb))
qsave(baseline_pca, file.path(root, 'models/transfer_baseline_pca.qs'))

# Learn the outcome-free gene reference from constituent profiles used in
# training combinations. No joint RNA or excluded combination outcome enters.
profile_cols <- c('A_profile_low', 'A_profile_high', 'B_profile_low', 'B_profile_high')
training_profiles <- unique(unlist(obs[train, ..profile_cols], use.names = FALSE))
training_profiles <- intersect(training_profiles, colnames(molecular$displacement))
gene_pca <- prcomp(t(molecular$displacement[, training_profiles, drop = FALSE]),
  center = TRUE, scale. = FALSE, rank. = gene_components)
gene_scores <- t(predict(gene_pca, newdata = t(molecular$displacement))[, seq_len(gene_components), drop = FALSE])
qsave(gene_pca, file.path(root, 'models/transfer_gene_pca.qs'))
membership <- fread(file.path(root, 'data/representation_membership.tsv'))
eligibility <- fread(eligibility_path)
representation_scores <- list()
feature_manifest <- list()
for (rep in c('bank16', 'bank42', 'bank43', 'biology')) {
  features <- eligibility[representation == rep & eligible_LINCS == TRUE, feature]
  stopifnot(length(features) >= 1L)
  zz <- membership[representation == rep & feature %chin% features & gene %chin% rownames(molecular$displacement)]
  scores <- t(vapply(features, function(ff) {
    genes <- zz[feature == ff, gene]
    stopifnot(length(genes) >= 5L)
    colMeans(molecular$displacement[genes, , drop = FALSE])
  }, numeric(ncol(molecular$displacement))))
  dimnames(scores) <- list(features, colnames(molecular$displacement))
  representation_scores[[rep]] <- scores
  feature_manifest[[rep]] <- zz[, .(representation, feature, gene)]
}
representation_scores[['genes']] <- gene_scores
# Direct mechanistic RNA anchors remain individually measured features. They
# are not called complete pathway scores or phosphorylation measurements.
anchor_registry <- data.table(gene = c('DUSP6', 'DUSP4', 'SPRY2', 'ETV4', 'ETV5',
  'PHLDA1', 'FOS', 'EGR1', 'MYC', 'CCND1', 'FOXO3', 'AKT1', 'EGFR', 'ERBB2',
  'ERBB3', 'MET', 'RPS6KB1', 'EIF4EBP1'))
anchor_registry[, measured := gene %chin% rownames(molecular$displacement)]
write_table(anchor_registry, 'transfer_direct_RNA_anchor_registry.tsv')
representation_scores[['anchors']] <- molecular$displacement[anchor_registry[measured == TRUE, gene], , drop = FALSE]
write_table(rbindlist(feature_manifest), 'transfer_measured_feature_membership.tsv')
qsave(representation_scores, file.path(root, 'models/transfer_constituent_scores.qs'))

# Shared interpolation and interaction design, used for every representation.
interpolate_arm <- function(scores, arm) {
  low <- obs[[paste0(arm, '_profile_low')]]
  high <- obs[[paste0(arm, '_profile_high')]]
  w <- obs[[paste0(arm, '_high_weight')]]
  stopifnot(all(low %in% colnames(scores)), all(high %in% colnames(scores)), all(w >= 0 & w <= 1))
  t(sweep(scores[, low, drop = FALSE], 2, 1 - w, '*') + sweep(scores[, high, drop = FALSE], 2, w, '*'))
}
response_design <- function(a, b) {
  stopifnot(identical(dim(a), dim(b)))
  sums <- a + b
  colnames(sums) <- paste0('sum_', colnames(a))
  ij <- which(upper.tri(matrix(0, ncol(a), ncol(a)), diag = TRUE), arr.ind = TRUE)
  products <- vapply(seq_len(nrow(ij)), function(k) {
    i <- ij[k, 1]; j <- ij[k, 2]
    (a[, i] * b[, j] + a[, j] * b[, i]) / 2
  }, numeric(nrow(a)))
  colnames(products) <- paste0('cross_', colnames(a)[ij[, 1]], '__', colnames(a)[ij[, 2]])
  list(sum = sums, product = products)
}
designs <- list(Single_calibrated = single_x, Context = context_x)
for (rep in names(representation_scores)) {
  a <- interpolate_arm(representation_scores[[rep]], 'A')
  b <- interpolate_arm(representation_scores[[rep]], 'B')
  response <- response_design(a, b)
  reverse <- response_design(b, a)
  stopifnot(max(abs(response$product - reverse$product)) < 1e-10)
  bx <- do.call(cbind, lapply(seq_len(ncol(baseline_x)), function(j) {
    z <- response$sum * baseline_x[, j]
    colnames(z) <- paste0(colnames(baseline_x)[j], '__', colnames(z)); z
  }))
  designs[[paste0(rep, '_response')]] <- cbind(single_x, cell_x, response$sum, response$product)
  designs[[paste0(rep, '_conditioned')]] <- cbind(context_x, response$sum, response$product, bx)
}
stopifnot(all(vapply(designs, function(z) all(is.finite(z)), logical(1))))

# Standardization is fitted separately within every pair-exclusion fold.
fit_ridge <- function(x, y, w, lambda) {
  center <- colSums(x * w) / sum(w)
  scale <- sqrt(colSums(sweep(x, 2, center, '-')^2 * w) / sum(w))
  keep <- is.finite(scale) & scale > 1e-8
  z <- sweep(sweep(x[, keep, drop = FALSE], 2, center[keep], '-'), 2, scale[keep], '/')
  fit <- glmnet(z, y, weights = w, alpha = 0, lambda = lambda,
    standardize = FALSE, intercept = TRUE)
  list(fit = fit, center = center[keep], scale = scale[keep], keep = keep)
}
predict_ridge <- function(object, x, lambda) {
  z <- sweep(sweep(x[, object$keep, drop = FALSE], 2, object$center, '-'), 2, object$scale, '/')
  as.matrix(predict(object$fit, newx = z, s = lambda))
}
pred_parts <- list()
tuning_parts <- list()
coefficient_parts <- list()
capacity_parts <- list()
cor_diagnostic <- function(a, b) if (sd(a) > 1e-8 && sd(b) > 1e-8) cor(a, b) else NA_real_
for (model_name in names(designs)) {
  message('Fitting ', model_name, ': ', ncol(designs[[model_name]]), ' inputs; ', length(pairs), ' training pairs')
  x <- designs[[model_name]]
  x_train <- x[train, , drop = FALSE]
  oof <- matrix(NA_real_, length(train), length(lambda_grid))
  fold_loss <- matrix(NA_real_, n_folds, length(lambda_grid))
  for (fold in seq_len(n_folds)) {
    fit_idx <- which(fold_id != fold)
    val_idx <- which(fold_id == fold)
    ff <- fit_ridge(x_train[fit_idx, , drop = FALSE], y_train[fit_idx], weights_train[fit_idx], lambda_grid)
    pp <- predict_ridge(ff, x_train[val_idx, , drop = FALSE], lambda_grid)
    oof[val_idx, ] <- pp
    fold_loss[fold, ] <- colSums(abs(pp - y_train[val_idx]) * weights_train[val_idx]) / sum(weights_train[val_idx])
  }
  loss_mean <- colMeans(fold_loss)
  loss_se <- apply(fold_loss, 2, sd) / sqrt(n_folds)
  best <- which.min(loss_mean)
  selected <- min(which(loss_mean <= loss_mean[best] + loss_se[best]))
  lambda <- lambda_grid[selected]
  tuning_parts[[model_name]] <- data.table(model = model_name, lambda = lambda_grid,
    internal_pair_cv_absolute_loss = loss_mean, fold_se = loss_se,
    selected = seq_along(lambda_grid) == selected, rule = 'one_standard_error_largest_lambda')
  final_fit <- fit_ridge(x_train, y_train, weights_train, lambda_grid)
  fit_path <- predict_ridge(final_fit, x_train, lambda_grid)
  capacity_parts[[model_name]] <- rbindlist(lapply(seq_along(lambda_grid), function(j) {
    data.table(model = model_name, lambda = lambda_grid[j], selected = j == selected,
      training_extra_effect_r = cor_diagnostic(y_train, fit_path[, j]),
      pair_OOF_extra_effect_r = cor_diagnostic(y_train, oof[, j]),
      observed_extra_at_least10pp = sum(y_train >= 10),
      training_predicted_extra_at_least10pp = sum(fit_path[, j] >= 10),
      pair_OOF_predicted_extra_at_least10pp = sum(oof[, j] >= 10),
      training_prediction_p05 = quantile(fit_path[, j], .05),
      training_prediction_p95 = quantile(fit_path[, j], .95))
  }))
  pp <- as.numeric(predict_ridge(final_fit, x, lambda))
  qsave(list(model = final_fit, selected_lambda = lambda, inputs = colnames(x),
    train_observation_id = obs$observation_id[train], training_pairs = pairs, folds = folds,
    representation_is_outcome_frozen = TRUE),
    file.path(root, 'models', paste0('transfer_', model_name, '.qs')))
  cc <- as.matrix(coef(final_fit$fit, s = lambda))
  coefficient_parts[[model_name]] <- data.table(model = model_name, term = rownames(cc), coefficient = cc[, 1])
  pred_parts[[model_name]] <- rbind(
    data.table(observation_id = obs$observation_id, model = model_name,
      split = ifelse(obs$heldout, 'held_family', 'training_fit'), predicted_extra_pp = pp),
    data.table(observation_id = obs$observation_id[train], model = model_name,
      split = 'training_pair_OOF', predicted_extra_pp = oof[, selected]))
  # Flush every model and its predictions: interrupted work remains usable.
  write_table(pred_parts[[model_name]], paste0('transfer_predictions_', model_name, '.tsv.gz'))
  write_table(rbindlist(tuning_parts), 'transfer_internal_tuning.tsv')
  write_table(rbindlist(capacity_parts), 'transfer_training_capacity_diagnostic.tsv')
}
pred_parts[['Bliss']] <- rbind(
  data.table(observation_id = obs$observation_id, model = 'Bliss',
    split = ifelse(obs$heldout, 'held_family', 'training_fit'), predicted_extra_pp = 0),
  data.table(observation_id = obs$observation_id[train], model = 'Bliss',
    split = 'training_pair_OOF', predicted_extra_pp = 0))
predictions <- rbindlist(pred_parts)
write_table(predictions, 'transfer_predictions_before_outcome_join.tsv.gz')
writeLines(digest(predictions, algo = 'sha256'), file.path(root, 'models/transfer_prediction_sha256.txt'))
write_table(rbindlist(coefficient_parts), 'transfer_standardized_coefficients.tsv')

# Only after all prespecified predictions have been saved, reattach the held
# measured outcomes. Never change representation or tuning based on these.
observed <- fread(links_path)[observation_id %chin% obs$observation_id,
  .(observation_id, cell, pair, drug_A, dose_A_uM, drug_B, dose_B_uM,
    curve_id, plate_id, source_release, assay_day, viability_A, viability_B,
    viability_AB, observed_extra_pp = excess_inhibition,
    bliss_viability = bliss_expected_viability, heldout,
    focal = pair == 'pictilisib + trametinib')]
predictions <- merge(predictions, observed, by = 'observation_id', all.x = TRUE, sort = FALSE)
predictions[, predicted_viability := bliss_viability - predicted_extra_pp]
write_table(predictions, 'transfer_predictions.tsv.gz')

# Independent units remain cell/pair/curve contexts, never individual technical
# wells or the repeated dose points portrayed as independent sample sizes.
condition <- predictions[, .(observed_extra_pp = mean(observed_extra_pp),
  predicted_extra_pp = mean(predicted_extra_pp), viability_AB = mean(viability_AB),
  predicted_viability = mean(predicted_viability), viability_A = mean(viability_A),
  viability_B = mean(viability_B), n_source_wells = .N),
  by = .(model, split, cell, pair, drug_A, dose_A_uM, drug_B, dose_B_uM, focal)]
cor_safe <- function(a, b) if (length(a) >= 4 && sd(a) > 1e-8 && sd(b) > 1e-8) cor(a, b) else NA_real_
summary_table <- condition[, .(dose_contexts = .N, cells = uniqueN(cell), pairs = uniqueN(pair),
  mean_observed_extra_pp = mean(observed_extra_pp), mean_predicted_extra_pp = mean(predicted_extra_pp),
  observed_extra_at_least_10pp = sum(observed_extra_pp >= 10),
  predicted_extra_at_least_10pp = sum(predicted_extra_pp >= 10),
  direction_correct_at_observed_abs10pp = if (any(abs(observed_extra_pp) >= 10))
    mean(sign(predicted_extra_pp[abs(observed_extra_pp) >= 10]) == sign(observed_extra_pp[abs(observed_extra_pp) >= 10])) else NA_real_,
  viability_r = cor_safe(viability_AB, predicted_viability),
  extra_effect_r = cor_safe(observed_extra_pp, predicted_extra_pp)), by = .(model, split)]
write_table(condition, 'transfer_condition_predictions.tsv.gz')
write_table(summary_table, 'transfer_biological_summary.tsv')
source_files <- c(links_path, molecular_path, eligibility_path,
  '/work/src/program_pilot/60_combination_transfer.R')
write_table(data.table(path = source_files, sha256 = vapply(source_files, digest, character(1), file = TRUE, algo = 'sha256')),
  'transfer_model_source_manifest.tsv')
writeLines(c(
  'All MEK+PI3K family combinations are excluded in all four cells, doses, releases and orientations.',
  'Four seen cells: this is retrospective intervention-family transfer, not unseen-cell validation.',
  'Target is measured-single Bliss viability minus observed combination viability, in percentage points.',
  'Observed constituent function is an input. This does not test a fully molecular-only prediction bridge.',
  'Ridge penalty selected by whole-pair four-fold CV, using internal absolute loss and the one-standard-error rule.',
  'Pairs receive equal weight within cells; cells receive equal total weight. Wells are not independent biological replicates.',
  'Prediction features are exchange-symmetric constituent response sums and cross-products, with optional baseline-PC modulation.',
  'Two baseline PCs across four seen cells cannot establish transferable baseline conditioning.',
  'All frozen eligible representations are reported; none is selected on the excluded-family result.',
  'Only molecular-supported dose brackets are predicted; all other observed focal doses remain in transfer_molecular_links.tsv.gz.',
  'The 10 percentage-point effect threshold is descriptive, not a replicate-derived significance test.',
  paste('Training unordered pairs:', length(pairs), '; training wells:', length(train), '; reserved wells:', length(test)),
  paste('Finished', Sys.time())
), file.path(root, 'TRANSFER_MODEL_DESIGN.txt'))
writeLines(capture.output(sessionInfo()), file.path(root, 'transfer_model_sessionInfo.txt'))
print(summary_table[split == 'held_family'])
