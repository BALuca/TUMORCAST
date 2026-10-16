#!/usr/bin/env Rscript
# Preserve51's fitted bridge. Learn its single-effect calibration around vehicle
# zero using development-cell OOF scores, without refitting any molecular model.
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
options(warn = 1)
root <- '/work/results/program_pilot/52_vehicle_anchored_program_bridge'
parent <- '/work/results/program_pilot/51_restored_program_baseline_function_combination_bridge'
decoder_root <- '/work/results/program_pilot/50_restore_calibrated_svr_current_program_bank'
families <- c('ERBB', 'MEK')
held_ids <- qread(file.path(parent, 'settings.qs'))$held_ids
rank <- 34L
zero <- matrix(0, 1L, rank)
for (p in file.path(root, c('data', 'models', 'figures'))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', name), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(6); quit(save = 'no', status = 1L)})
status('replaying frozen bridge; training vehicle-anchored calibration only')
fp <- fread(file.path(parent, 'data/function_predictions.tsv.gz'))
dp <- fread(file.path(parent, 'data/combination_predictions.tsv.gz'))
fp[, version := '51 affine calibration']
dp[, version := '51 affine calibration']
frozen50 <- fread(file.path(decoder_root, 'data/evaluated_predictions.tsv.gz'))
fp_new <- copy(fp)
fp_new[, version := '52 vehicle anchored']
boundaries <- models <- list()
predict_anchored_effect <- function(fitted, slope, full_zero, x) {
  slope * (predict_program_svr(fitted, x) - full_zero)
}

# Each calibration observation subtracts its OWN held-cell fold model's zero.
# Both slope and scoring zero are saved; no evaluation outcomes enter this fit.
for (family in families) {
  fitted <- qread(file.path(decoder_root, 'models', paste0(family, '_current_decoder.qs')))
  old_cal <- qread(file.path(decoder_root, 'models', paste0(family, '_current_calibration.qs')))
  family_key <- family
  train <- frozen50[family == family_key & bank == '43 current / 34 programs' &
    calibration == 'Raw' & evaluation == 'Own fitting data']
  train <- train[match(old_cal$training_ids, observation_id)]
  check(identical(train$observation_id, old_cal$training_ids) && !any(train$depmap_id %chin% held_ids),
    'PRISM calibration training identities differ')
  folds <- fread(file.path(decoder_root, 'data', paste0(family, '_original_folds.tsv')))
  fold <- folds$fold[match(train$depmap_id, folds$depmap_id)]
  check(!anyNA(fold), 'Missing PRISM whole-cell fold')
  fold_zero <- numeric(nrow(train))
  for (k in 1:3) {
    ff <- qread(file.path(decoder_root, 'models', paste0(family, '_current_fold_', k, '.qs')))
    fold_zero[fold == k] <- predict_program_svr(ff$fitted, zero)
  }
  train[, centered_oof := old_cal$oof - fold_zero]
  calibration <- lm(observed_log2 ~ 0 + centered_oof, data = train, weights = observation_weights(train))
  slope <- unname(coef(calibration))
  check(length(slope) == 1L && is.finite(slope) && slope > 0, 'Invalid PRISM anchored slope')
  full_zero <- predict_program_svr(fitted$fitted, zero)
  neutral <- 100 * 2^predict_anchored_effect(fitted$fitted, slope, full_zero, zero)
  check(abs(neutral - 100) < 1e-10, 'PRISM vehicle boundary failed')
  old_coefficients <- unname(coef(fitted$calibration))
  rows <- which(fp_new$decoder == family & fp_new$model %chin% c('Measured response', 'Predicted response', 'Mean response'))
  raw <- (fp_new$predicted_log2[rows] - old_coefficients[1]) / old_coefficients[2]
  set(fp_new, rows, 'predicted_log2', slope * (raw - full_zero))
  set(fp_new, rows, 'predicted_pct', 100 * 2^fp_new$predicted_log2[rows])
  models[[family]] <- list(calibration = calibration, slope = slope, full_zero = full_zero,
    training_ids = train$observation_id, source_model = file.path(decoder_root, 'models', paste0(family, '_current_decoder.qs')))
  qsave(models[[family]], file.path(root, 'models', paste0(family, '_anchored_calibration.qs')))
  boundaries[[family]] <- data.table(assay = 'PRISM', head = family,
    before_single_zero_pct = 100 * 2^(old_coefficients[1] + old_coefficients[2] * full_zero),
    after_single_zero_pct = neutral, before_two_zeros_pct = 100 * 2^(2 * (old_coefficients[1] + old_coefficients[2] * full_zero)),
    after_two_zeros_pct = 100 * (neutral / 100)^2, slope = slope)
}

# The interaction predictions are reused verbatim. Only matched ALMANAC single
# effects are calibrated differently; all51 controls retain their predictions.
dp_new <- copy(dp)
dp_new[, version := '52 vehicle anchored']
original_combo <- as.data.table(readRDS('/work/almanac_broad9_2026-10-07/results/prepared.rds')$d)
for (s in unique(dp$split)) {
  fitted <- qread(file.path(parent, 'models', paste0('almanac_single_decoder_', s, '.qs')))
  train <- copy(fitted$training_singles)
  check(!any(train$depmap_id %chin% held_ids), 'Held ALMANAC cell in anchored calibration')
  fold_zero <- numeric(nrow(train))
  for (cell in sort(unique(train$depmap_id))) {
    ff <- qread(file.path(parent, 'models', paste0('almanac_single_', s, '_', cell, '.qs')))
    check(!cell %in% ff$training_ids, 'Single-function calibration fold leakage')
    fold_zero[train$depmap_id == cell] <- predict_program_svr(ff$fitted, zero)
  }
  train[, centered_oof := fitted$oof - fold_zero]
  calibration <- lm(observed ~ 0 + centered_oof, data = train, weights = observation_weights(train))
  slope <- unname(coef(calibration))
  check(length(slope) == 1L && is.finite(slope) && slope > 0, 'Invalid ALMANAC anchored slope')
  full_zero <- predict_program_svr(fitted$fitted, zero)
  neutral <- 100 * 2^predict_anchored_effect(fitted$fitted, slope, full_zero, zero)
  check(abs(neutral - 100) < 1e-10, 'ALMANAC vehicle boundary failed')
  cf <- unname(coef(fitted$calibration))
  singles <- fread(file.path(parent, 'data', paste0('almanac_single_predictions_', s, '.tsv.gz')))
  key <- function(cell, drug, dose) paste(cell, drug, dose, sep = '__')
  singles[, coordinate := key(depmap_id, drug, dose)]
  check(all(singles[, max(almanac_matched) - min(almanac_matched), by = coordinate]$V1 < 1e-10),
    'Identical molecular queries have different predictions')
  singles <- unique(singles, by = 'coordinate')
  rows <- which(dp_new$split == s & dp_new$input %chin%
    c('Measured responses | ALMANAC head', 'Baseline only | ALMANAC head'))
  cd <- original_combo[match(dp_new$source_row[rows], source_row)]
  qa <- match(key(cd$depmap_id, cd$idA, cd$doseA), singles$coordinate)
  qb <- match(key(cd$depmap_id, cd$idB, cd$doseB), singles$coordinate)
  check(!anyNA(qa) && !anyNA(qb), 'Single queries missing in frozen52 composition')
  baseline_rows <- dp_new$input[rows] == 'Baseline only | ALMANAC head'
  single_a <- ifelse(baseline_rows, singles$baseline_almanac[qa], singles$almanac_matched[qa])
  single_b <- ifelse(baseline_rows, singles$baseline_almanac[qb], singles$almanac_matched[qb])
  raw_a <- (single_a - cf[1]) / cf[2]
  raw_b <- (single_b - cf[1]) / cf[2]
  additive <- slope * (raw_a + raw_b - 2 * full_zero)
  predicted <- additive + dp_new$predicted_interaction[rows]
  set(dp_new, rows, 'predicted_log2', predicted)
  set(dp_new, rows, 'predicted_pct', 100 * 2^predicted)
  set(dp_new, rows, 'predicted_bliss_pp', 100 * (2^additive - 2^predicted))
  qsave(list(calibration = calibration, slope = slope, full_zero = full_zero,
    source_model = file.path(parent, 'models', paste0('almanac_single_decoder_', s, '.qs')),
    training_ids = unique(train$depmap_id), training_pairs = fitted$training_pairs),
    file.path(root, 'models', paste0('ALMANAC_', s, '_anchored_calibration.qs')))
  # Validate a representative molecular partner and the exact neutral limits.
  f <- qread(file.path(parent, 'models', paste0('combination_', s, '_molecular_products.qs')))$model
  interaction_zero <- predict_ridge(f, matrix(0, 1L, nrow(coef(f$fit, s = f$lambda)) - 1L))
  other_single <- slope * (raw_a[1L] - full_zero)
  check(abs(interaction_zero) < 1e-12 && abs((0 + other_single + interaction_zero) - other_single) < 1e-12,
    'Zero-partner composition identity failed')
  boundaries[[s]] <- data.table(assay = 'ALMANAC', head = s,
    before_single_zero_pct = 100 * 2^(cf[1] + cf[2] * full_zero), after_single_zero_pct = neutral,
    before_two_zeros_pct = 100 * 2^(2 * (cf[1] + cf[2] * full_zero)),
    after_two_zeros_pct = 100 * (neutral / 100)^2, slope = slope)
  status(paste('saved anchored calibration:', s))
}
check(identical(dp_new$source_row, dp$source_row) &&
  max(abs(dp_new$predicted_interaction - dp$predicted_interaction)) == 0, 'Interaction predictions changed')
check(all(is.finite(fp_new$predicted_log2)) && all(is.finite(dp_new$predicted_log2)), 'Nonfinite anchored predictions')
responses <- qread(file.path(parent, 'data/bridge_responses.qs'))
replay <- list()
for (family in families) {
  row <- fp_new[decoder == family & model == 'Measured response' & split == 'Held cells'][1L]
  index <- match(row$observation_id, responses$metadata$observation_id)
  f <- qread(models[[family]]$source_model)
  predicted <- predict_anchored_effect(f$fitted, models[[family]]$slope,
    models[[family]]$full_zero, responses$actual[index, , drop = FALSE])
  check(abs(predicted - row$predicted_log2) < 1e-10, 'Frozen feature-level anchored prediction replay failed')
  replay[[family]] <- data.table(head = family, observation_id = row$observation_id,
    absolute_prediction_difference = abs(predicted - row$predicted_log2))
}
write_table(rbindlist(replay), 'direct_frozen_feature_replay.tsv')
write_table(fp_new, 'function_predictions_frozen.tsv.gz')
write_table(dp_new, 'combination_predictions_frozen.tsv.gz')
prediction_paths <- file.path(root, 'data', c('function_predictions_frozen.tsv.gz', 'combination_predictions_frozen.tsv.gz'))
write_table(data.table(path = prediction_paths,
  sha256 = vapply(prediction_paths, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'prediction_hashes.tsv')
write_table(rbindlist(boundaries), 'neutral_response_boundary.tsv')

# Report matched51/52 results after freezing predictions. Correlation shape is
# retained within a fitted head; absolute effect calibration and calls can change.
function_summary <- function(z) {
  context <- paste(z$drug, z$dose)
  oy <- z$observed_log2 - ave(z$observed_log2, context, FUN = mean)
  py <- z$predicted_log2 - ave(z$predicted_log2, context, FUN = mean)
  cross <- z[, .(r = cor_safe(observed_log2, predicted_log2)), by = .(depmap_id, dose)]
  data.table(rows = nrow(z), cells = uniqueN(z$depmap_id), pooled_r = cor_safe(z$observed_log2, z$predicted_log2),
    across_cells_fixed_treatment_r = if (sd(py) > 1e-10) cor_safe(oy, py) else NA_real_,
    median_across_drugs_same_cell_dose_r = median(cross$r, na.rm = TRUE),
    sensitive = sum(z$measured_pct <= 50), recovered = sum(z$measured_pct <= 50 & z$predicted_pct <= 50),
    resistant = sum(z$measured_pct > 50), false_sensitive = sum(z$measured_pct > 50 & z$predicted_pct <= 50))
}
f_all <- rbind(fp, fp_new)
d_all <- rbind(dp, dp_new)
fm <- f_all[, function_summary(.SD), by = .(version, decoder, model, split)]
cm <- d_all[, .(rows = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
  endpoint_r = cor_safe(observed_log2, predicted_log2), interaction_r = cor_safe(observed_interaction, predicted_interaction),
  sensitive = sum(measured_pct <= 50), recovered = sum(measured_pct <= 50 & predicted_pct <= 50),
  resistant = sum(measured_pct > 50), false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50)),
  by = .(version, scope, input, model)]
write_table(fm, 'single_drug_calibration_comparison.tsv')
write_table(cm, 'combination_calibration_comparison.tsv')
write_table(f_all, 'function_predictions.tsv.gz')
write_table(d_all, 'combination_predictions.tsv.gz')

# Compact comparison with actual measured outcomes; independently scaled square
# panels, percentage guides and Pearson/trend are retained for interpretation.
plot_data <- d_all[scope == 'Held cells + pairs' & model == 'molecular_products' &
  input %chin% c('Measured responses | ALMANAC head', 'Baseline only | ALMANAC head')]
stats <- plot_data[, .(r = cor_safe(observed_log2, predicted_log2)), by = .(version, input)]
p <- ggplot(plot_data, aes(observed_log2, predicted_log2)) + geom_point(alpha = .2, size = .45) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
  geom_vline(xintercept = -1, linetype = 2, color = 'grey50') + geom_hline(yintercept = -1, linetype = 2, color = 'grey50') +
  facet_wrap(~version + input, scales = 'free', ncol = 2) +
  geom_text(data = stats, aes(x = -Inf, y = Inf, label = sprintf('r=%.2f', r)),
    hjust = -.1, vjust = 1.3, inherit.aes = FALSE) +
  labs(title = 'Physical single-effect calibration | held recipients and pairs',
    x = 'Measured combination log2(% vehicle/100)', y = 'Predicted combination log2(% vehicle/100)') +
  theme_classic(base_size = 10) + theme(aspect.ratio = 1, strip.text = element_text(size = 9))
ggsave(file.path(root, 'figures/vehicle_anchored_bridge.pdf'), p, width = 9, height = 8.5)
ggsave(file.path(root, 'figures/vehicle_anchored_bridge.png'), p, width = 9, height = 8.5, dpi = 130)
checks <- data.table(check = c('calibration excludes held cells', 'whole-cell fold zeros used for OOF calibration',
  'original source rows retained', 'interaction predictions unchanged', 'neutral single effect100%',
  'two neutral responses100%', 'zero partner equals remaining single effect'), passed = TRUE)
write_table(checks, 'acceptance.tsv')
input_paths <- c(file.path(parent, 'settings.qs'), file.path(parent, 'data', c('function_predictions.tsv.gz',
  'combination_predictions.tsv.gz')), list.files(file.path(parent, 'models'), pattern = 'almanac_single.*qs$', full.names = TRUE),
  list.files(file.path(decoder_root, 'models'), pattern = 'current_.*qs$', full.names = TRUE))
write_table(data.table(path = input_paths,
  sha256 = vapply(input_paths, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'input_hashes.tsv')
write_table(data.table(path = '/work/src/program_pilot/52_vehicle_anchored_program_bridge.R',
  sha256 = digest::digest(file = '/work/src/program_pilot/52_vehicle_anchored_program_bridge.R', algo = 'sha256')), 'code_hash.tsv')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
report <- c('52: vehicle-anchored calibration of the frozen51 whole program bridge',
  'No SVR, baseline-response, program bank or interaction refitting. Development OOF calibration only.',
  'Neutral program input is defined as vehicle-like for composition; responses remain signed and unclipped.',
  'This boundary requirement does not prove the34-program representation retains all function-relevant information.',
  'Retrospective43 discovery and earlier exposed evaluation outcomes remain unchanged limitations.',
  'Boundary correction does not restore predictive sensitivity: the anchored ALMANAC gain shrinks strongly.',
  'Use the counts below to assess this failure; physically valid composition alone is insufficient grant evidence.',
  '', 'Neutral-response boundaries:', capture.output(print(rbindlist(boundaries))),
  '', 'Held focal function comparison:', capture.output(print(fm[split == 'Held cells' & decoder == 'ERBB'])),
  '', 'Held-cell/held-pair combination comparison:',
  capture.output(print(cm[scope == 'Held cells + pairs' & input %chin%
    c('Measured responses | ALMANAC head', 'Baseline only | ALMANAC head', 'Treatment-average singles')])))
writeLines(report, file.path(root, 'REPORT.txt'))
status('COMPLETE')
