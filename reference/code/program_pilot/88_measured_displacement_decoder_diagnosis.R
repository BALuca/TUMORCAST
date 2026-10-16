#!/usr/bin/env Rscript
# Locate v87 failures using a shared decoder trained only on measured RNA.
# This is a diagnostic of cross-cell sensitivity, not a replacement dose model.
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(qs)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/88_measured_displacement_decoder_diagnosis'
v84 <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
v87 <- '/work/results/program_pilot/87_selected_dose_joint_program_bridge'
programs <- sprintf('P%02d', 1:40)
lambda <- .1
dir.create(file.path(root, 'data'), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
pearson <- function(x, y) if (length(x) >= 4L && sd(x) > 1e-12 && sd(y) > 1e-12) cor(x, y) else NA_real_

selected <- fread(file.path(v87, 'data/selected_doses.tsv'))
points <- fread(file.path(v84, 'data/observed_program_viability_points.tsv.gz'))
points <- merge(points, selected[, .(drug, dose)], by = c('drug', 'dose'))
wide <- dcast(points, depmap_id + drug + dose + observed_log2_abundance ~ program,
  value.var = 'observed_displacement')
folds <- fread(file.path(v87, 'data/recipient_folds.tsv'))
wide <- merge(wide, folds[, .(depmap_id, fold)], by = 'depmap_id', all.x = TRUE)
stopifnot(!anyNA(wide$fold), !anyDuplicated(wide[, .(depmap_id, drug)]),
  nrow(wide) == sum(selected$cells), all(is.finite(as.matrix(wide[, ..programs]))))

# The same 130 functional features as v87: 40 linear coordinates, 40 tanh terms,
# 40 saturating square terms, and ten products of four fixed latent coordinates.
latent_pairs <- which(upper.tri(matrix(0, 4L, 4L), diag = TRUE), arr.ind = TRUE)
nonlinear_features <- function(x, u) {
  latent <- tanh(x %*% u)
  cbind(x, tanh(x), x^2 / (1 + x^2),
    latent[, latent_pairs[, 1], drop = FALSE] * latent[, latent_pairs[, 2], drop = FALSE])
}

# Fit one common coefficient vector across all drugs. Training features and
# outcomes are centered per drug to isolate cell-to-cell variation. Predictions
# below are functional scores, without a learned drug-specific calibration head.
predictions <- coefficients <- list()
for (k in 1:3) {
  model <- qread(file.path(v87, 'models', paste0('recipient_', k, '__Tahoe_+_LINCS.qs')))
  train <- which(wide$fold != k)
  test <- which(wide$fold == k)
  stopifnot(!length(intersect(wide$depmap_id[train], wide$depmap_id[test])))
  x <- sweep(as.matrix(wide[, ..programs]), 2L, model$context$response_sd, '/')
  design <- list(linear = x, nonlinear = nonlinear_features(x, model$decoder$u))
  for (family in names(design)) {
    xt <- design[[family]][train, , drop = FALSE]
    yt <- wide$observed_log2_abundance[train]
    train_drug <- wide$drug[train]
    groups <- unique(train_drug)
    for (drug in groups) {
      ii <- which(train_drug == drug)
      xt[ii, ] <- sweep(xt[ii, , drop = FALSE], 2L, colMeans(xt[ii, , drop = FALSE]))
      yt[ii] <- yt[ii] - mean(yt[ii])
    }
    weight <- 1 / as.numeric(table(train_drug)[train_drug])
    fit <- glmnet(xt, yt, weights = weight, alpha = 0, lambda = lambda,
      intercept = FALSE, standardize = TRUE, thresh = 1e-10)
    stopifnot(fit$jerr == 0L)
    beta <- as.numeric(coef(fit))[-1L]
    out <- wide[test, .(depmap_id, drug, dose, fold, observed_log2 = observed_log2_abundance)]
    out[, `:=`(decoder = family,
      predicted_function_score = as.numeric(design[[family]][test, , drop = FALSE] %*% beta))]
    predictions[[paste(k, family)]] <- out
    coefficients[[paste(k, family)]] <- data.table(fold = k, decoder = family,
      feature = seq_along(beta), coefficient = beta)
  }
}
pred <- rbindlist(predictions)
pred[, observed_within_fold := observed_log2 - mean(observed_log2), by = .(drug, decoder, fold)]
pred[, predicted_within_fold := predicted_function_score - mean(predicted_function_score),
  by = .(drug, decoder, fold)]
metrics <- pred[, .(cells = uniqueN(depmap_id),
  within_fold_r = pearson(observed_within_fold, predicted_within_fold)), by = .(drug, dose, decoder)]
write_table(pred, 'measured_only_shared_decoder_predictions')
write_table(rbindlist(coefficients), 'measured_only_shared_decoder_coefficients')
write_table(metrics, 'measured_only_shared_decoder_correlations')

# Compare the actual v87 decoder on the identical measured profiles, removing
# fold offsets from both axes for the same sensitivity comparison.
joint <- fread(file.path(v87, 'data/held_measured_displacement_function_frozen.tsv'))
joint[, observed_within_fold := observed_log2 - mean(observed_log2), by = .(drug, fold)]
joint[, predicted_within_fold := measured_displacement_prediction - mean(measured_displacement_prediction),
  by = .(drug, fold)]
joint_metrics <- joint[, .(joint_decoder_r = pearson(observed_within_fold, predicted_within_fold)),
  by = .(drug, dose)]
comparison <- merge(dcast(metrics, drug + dose + cells ~ decoder, value.var = 'within_fold_r'),
  joint_metrics, by = c('drug', 'dose'))

# Track whether the visibly informative program is actually predictable from
# baseline. These are separate comparison axes, explicitly named in the table.
associations <- fread(file.path(v84, 'data/observed_program_viability_each_dose_correlations.tsv'))
associations <- merge(associations, selected[, .(drug, dose)], by = c('drug', 'dose'))
strongest <- associations[, .SD[which.max(abs(r))], by = drug]
molecular <- fread(file.path(v87, 'data/held_program_displacement_correlations.tsv'))
chain <- merge(strongest[, .(drug, dose, program, observed_program_to_viability_r = r)],
  molecular[, .(drug, dose, program, baseline_to_program_r = pearson)], by = c('drug', 'dose', 'program'))
write_table(chain, 'strongest_observed_program_and_baseline_transfer')
write_table(comparison, 'decoder_diagnosis_comparison')
print(comparison)
print(chain)
writeLines(c('Retrospective diagnosis: doses were selected using all observed functional outcomes upfront.',
  'All fitting excludes the recipient being scored, including every drug observation from that recipient.',
  'Shared measured-only ridge decoders use fixed lambda=0.1; no tuning on these diagnostic correlations.',
  'Scores have no abundance calibration. Pearson compares cell variation within a drug/dose and coefficient-fitting fold.',
  'The linear and nonlinear diagnostic models use measured RNA only; no imputed RNA, combinations or direct baseline inputs.'),
  file.path(root, 'DIAGNOSTIC_SCOPE.txt'))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
