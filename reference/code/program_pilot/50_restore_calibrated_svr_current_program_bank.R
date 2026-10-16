#!/usr/bin/env Rscript
# Restore the earlier calibrated RBF SVR; change only its molecular program bank.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(patchwork)
})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(20261010L)
options(warn = 1)

root <- '/work/results/program_pilot/50_restore_calibrated_svr_current_program_bank'
input_path <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr/data/frozen_comparison_inputs.qs'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
families <- c('ERBB', 'MEK')
workers <- 3L
for (path in file.path(root, c('data', 'models', 'figures'))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('restoring original settings, folds and calibrated SVR formulation')
input <- qread(input_path)
panel <- qread(panel_path)
current <- input$scores[['43 averaged / 34 programs']]
check(identical(rownames(current), input$metadata$observation_id) && ncol(current) == 34L,
  'Current43 score identities or dimensions changed')
check(all(is.finite(current)), 'Nonfinite43 scores')
results <- settings_tables <- checks <- list()

# Every summary uses the same assay scale as the restored decoder. Whole-cell CV
# centering includes fold to avoid comparing predictions from different fits.
summarize <- function(z, evaluation_name) {
  z <- copy(z)
  if (evaluation_name == 'Whole-cell validation') {
    z[, context := paste(drug, dose, fold, sep = '__')]
  } else z[, context := paste(drug, dose, sep = '__')]
  z[, `:=`(yc = observed_log2 - mean(observed_log2), pc = predicted_log2 - mean(predicted_log2)), by = context]
  per_cell <- z[, .(r = cor_safe(observed_log2, predicted_log2)), by = depmap_id]
  s <- z$measured_pct <= 50
  data.table(observations = nrow(z), cells = uniqueN(z$depmap_id), drugs = uniqueN(z$drug),
    pooled_r = cor_safe(z$observed_log2, z$predicted_log2),
    median_within_cell_r = median(per_cell$r, na.rm = TRUE),
    across_cell_same_treatment_r = cor_safe(z$yc, z$pc),
    sensitive = sum(s), recovered_sensitive = sum(s & z$predicted_pct <= 50),
    resistant = sum(!s), false_sensitive = sum(!s & z$predicted_pct <= 50),
    measured_sensitive_median_pct = median(z$measured_pct[s]),
    predicted_sensitive_median_pct = median(z$predicted_pct[s]),
    measured_resistant_median_pct = median(z$measured_pct[!s]),
    predicted_resistant_median_pct = median(z$predicted_pct[!s]))
}

for (family_name in families) {
  parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family_name, '_holdout')
  paths <- file.path(parent, c('settings.qs', 'data/prepared_scores.qs',
    'models/nonlinear_RBF_SVR_D.qs', 'models/nonlinear_RBF_SVR_D_calibration.qs',
    'data/nonlinear_training_folds.tsv', 'data/nonlinear_predictions.tsv.gz'))
  settings <- qread(paths[1L])
  prepared <- qread(paths[2L])
  old <- qread(paths[3L])
  old_cal <- qread(paths[4L])
  folds <- fread(paths[5L])
  cached <- fread(paths[6L])
  labels <- copy(panel$sensitivity[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
    depmap_id %chin% rownames(prepared$baseline)])
  labels[, `:=`(measured_pct = prism_abundance_pct, observed_log2 = log2(prism_abundance_pct / 100),
    fold = folds$fold[match(depmap_id, folds$depmap_id)],
    drug_set = fifelse(drug %chin% settings$excluded_drugs, 'Excluded family', 'Seen drugs'))]
  train <- labels[!depmap_id %chin% settings$held_ids & !drug %chin% settings$excluded_drugs]
  test <- labels[depmap_id %chin% settings$held_ids]
  check(uniqueN(train$depmap_id) == 26L && uniqueN(test$depmap_id) == 9L &&
    setequal(unique(train$depmap_id), old$training_ids), 'Original cell split changed')
  check(setequal(settings$held_ids, input$held_ids) && !anyNA(train$fold) && all(is.na(test$fold)),
    'Original folds or outer-held cells changed')
  check(!anyDuplicated(labels$observation_id) && !any(train$drug %chin% settings$excluded_drugs),
    'Duplicate identities or excluded-family fitting leakage')
  original_train <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  original_test <- prepared$measured[match(test$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  check(ncol(original_train) == 40L && all(is.finite(original_train)) && all(is.finite(original_test)),
    'Original40 input alignment failed')
  check(max(abs(colMeans(original_train) - old$model$x_mean)) < 1e-12 &&
    max(abs(apply(original_train, 2L, sd) - old$model$x_sd)) < 1e-12,
    'Reconstructed original fitting order/scales differ')
  check(length(old$oof) == nrow(train) && all(is.finite(old$oof)), 'Missing original OOF predictions')
  ids <- c(train$observation_id, test$observation_id)
  ii <- match(ids, rownames(current))
  check(!anyNA(ii), 'Original comparison rows missing in43')
  assay_rows <- match(ids, input$metadata$observation_id)
  check(max(abs(input$metadata$measured_pct[assay_rows] - c(train$measured_pct, test$measured_pct))) < 1e-10,
    'Assay outcomes changed between old and current bank')
  new_train <- current[ii[seq_len(nrow(train))], , drop = FALSE]
  new_test <- current[ii[nrow(train) + seq_len(nrow(test))], , drop = FALSE]
  cost <- old$selected$cost
  multiplier <- old$selected$gamma_multiplier
  check(length(cost) == 1L && length(multiplier) == 1L && multiplier == 1 &&
    cost == if (family_name == 'MEK') 1 else 10, 'Original successful setting changed')
  manifest_paths <- c(paths, input_path, panel_path,
    '/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R')
  manifest <- data.table(path = manifest_paths,
    sha256 = vapply(manifest_paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
  signature <- digest::digest(list(manifest, ids, folds, cost, multiplier, 'restore20-calibrated-D-only'), algo = 'sha256')
  fwrite(manifest, file.path(root, 'data', paste0(family_name, '_input_hashes.tsv')), sep = '\t')
  settings_tables[[family_name]] <- data.table(family = family_name, cost, gamma_multiplier = multiplier,
    old_programs = 40L, current_programs = 34L, folds = 3L, calibration = 'equal-drug weighted OOF affine')
  fwrite(rbindlist(settings_tables), file.path(root, 'data/restored_settings.tsv'), sep = '\t')
  fwrite(folds, file.path(root, 'data', paste0(family_name, '_original_folds.tsv')), sep = '\t')

  # Flush each fold from its worker. Return paths, never unsaved model objects.
  status(paste(family_name, 'fitting3 folds with original C', cost, 'gamma multiplier', multiplier))
  fold_paths <- parallel::mclapply(1:3, function(k) {
    tr <- which(train$fold != k)
    va <- which(train$fold == k)
    check(!any(train$depmap_id[tr] %chin% train$depmap_id[va]), 'Fold cells overlap')
    path <- file.path(root, 'models', sprintf('%s_current_fold_%d.qs', family_name, k))
    if (file.exists(path)) {
      saved <- qread(path)
      check(identical(saved$signature, signature) && identical(saved$validation_ids, train$observation_id[va]),
        'Stale restoration fold checkpoint')
    } else {
      fitted <- fit_program_svr(new_train[tr, , drop = FALSE], train$observed_log2[tr],
        'radial', cost, multiplier / ncol(new_train))
      predicted <- predict_program_svr(fitted, new_train[va, , drop = FALSE])
      qsave(list(signature = signature, fitted = fitted, validation_ids = train$observation_id[va],
        prediction = predicted, fitting_ids = train$observation_id[tr],
        fitting_cells = unique(train$depmap_id[tr]), validation_cells = unique(train$depmap_id[va])), path)
    }
    message(family_name, ' fold', k, ' MODEL SAVED')
    path
  }, mc.cores = workers, mc.preschedule = FALSE)
  check(all(vapply(fold_paths, function(p) is.character(p) && length(p) == 1L && file.exists(p), logical(1))),
    'Restoration fold worker failed')
  oof <- rep(NA_real_, nrow(train))
  for (path in unlist(fold_paths)) {
    saved <- qread(path)
    oof[match(saved$validation_ids, train$observation_id)] <- saved$prediction
  }
  check(all(is.finite(oof)), 'Incomplete restoration OOF predictions')
  weights <- observation_weights(train)
  calibration <- lm(train$observed_log2 ~ oof, weights = weights)
  coefficients <- unname(coef(calibration))
  check(length(coefficients) == 2L && all(is.finite(coefficients)), 'Invalid OOF calibration')
  qsave(list(signature = signature, calibration = calibration, oof = oof,
    training_ids = train$observation_id), file.path(root, 'models', paste0(family_name, '_current_calibration.qs')))
  final_path <- file.path(root, 'models', paste0(family_name, '_current_decoder.qs'))
  status(paste(family_name, 'fitting original selected setting to current bank'))
  if (file.exists(final_path)) {
    saved <- qread(final_path)
    check(identical(saved$signature, signature), 'Stale restoration final model')
    fitted <- saved$fitted
  } else {
    fitted <- fit_program_svr(new_train, train$observed_log2, 'radial', cost, multiplier / ncol(new_train))
    qsave(list(signature = signature, fitted = fitted, calibration = calibration,
      training_ids = train$observation_id, excluded_cells = settings$held_ids,
      excluded_drugs = settings$excluded_drugs, cost = cost, gamma_multiplier = multiplier), final_path)
  }

  # Predict and hash held outcomes separately before reporting their labels.
  raw_held <- predict_program_svr(fitted, new_test)
  pred_path <- file.path(root, 'data', paste0(family_name, '_held_predictions_frozen.tsv'))
  fwrite(data.table(observation_id = test$observation_id, raw_log2 = raw_held,
    calibrated_log2 = coefficients[1] + coefficients[2] * raw_held), pred_path, sep = '\t')
  fwrite(data.table(path = pred_path, sha256 = digest::digest(file = pred_path, algo = 'sha256')),
    file.path(root, 'data', paste0(family_name, '_prediction_hash.tsv')), sep = '\t')
  raw_original <- predict_program_svr(old$model, original_test)
  original_coefficients <- unname(coef(old_cal$calibration))
  for (variant in c('Raw', 'Calibrated')) {
    source_model <- if (variant == 'Raw') 'RBF_SVR_D' else 'calibrated_RBF_SVR_D'
    expected <- cached[model == source_model][match(test$observation_id, observation_id)]
    replay <- if (variant == 'Raw') raw_original else original_coefficients[1] + original_coefficients[2] * raw_original
    check(identical(expected$observation_id, test$observation_id) &&
      max(abs(expected$predicted - replay)) < 1e-10, 'Original held predictions do not replay')
    check(max(abs(expected$observed - test$observed_log2)) < 1e-10, 'Original held labels changed')
  }
  own_current <- predict_program_svr(fitted, new_train)
  own_original <- predict_program_svr(old$model, original_train)
  prediction_rows <- list()
  for (bank_name in c('20 original / 40 programs', '43 current / 34 programs')) {
    original_bank <- startsWith(bank_name, '20')
    cal <- if (original_bank) original_coefficients else coefficients
    own <- if (original_bank) own_original else own_current
    held <- if (original_bank) raw_original else raw_held
    cv <- if (original_bank) old$oof else oof
    for (evaluation in c('Own fitting data', 'Whole-cell validation', 'Held cells')) {
      metadata <- copy(if (evaluation == 'Held cells') test else train)
      prediction <- switch(evaluation, 'Own fitting data' = own, 'Whole-cell validation' = cv, 'Held cells' = held)
      variants <- if (evaluation == 'Whole-cell validation') 'Raw' else c('Raw', 'Calibrated')
      for (variant in variants) {
        z <- copy(metadata)
        z[, `:=`(family = family_name, bank = bank_name, calibration = variant, evaluation = evaluation,
          predicted_log2 = if (variant == 'Raw') prediction else cal[1] + cal[2] * prediction)]
        z[, predicted_pct := 100 * 2^predicted_log2]
        check(all(is.finite(z$predicted_pct)), 'Nonfinite restored predictions')
        prediction_rows[[paste(bank_name, evaluation, variant)]] <- z
      }
    }
  }
  results[[family_name]] <- rbindlist(prediction_rows)
  predictions <- rbindlist(results)
  fwrite(predictions, file.path(root, 'data/evaluated_predictions.tsv.gz'), sep = '\t')
  summary <- predictions[, summarize(.SD, .BY$evaluation),
    by = .(family, bank, calibration, evaluation, drug_set)]
  fwrite(summary, file.path(root, 'data/comparison_summary.tsv'), sep = '\t')
  print(summary[family == family_name & evaluation == 'Held cells'])
  checks[[family_name]] <- data.table(family = family_name, training_observations = nrow(train),
    held_observations = nrow(test), original_prediction_replay = TRUE,
    current_feature_dimensions = ncol(new_train), restored_cost = cost,
    calibration_intercept = coefficients[1], calibration_slope = coefficients[2])
  fwrite(rbindlist(checks), file.path(root, 'data/acceptance.tsv'), sep = '\t')
  status(paste(family_name, 'COMPLETE; model and predictions saved'))
}

# Compact consolidated views retain every dose and distinguish comparison axes.
per_cell <- predictions[evaluation == 'Held cells' & calibration == 'Calibrated',
  .(conditions = .N, r = cor_safe(observed_log2, predicted_log2)), by = .(family, bank, drug_set, recipient)]
fwrite(per_cell, file.path(root, 'data/held_results_by_cell.tsv'), sep = '\t')
per_drug <- predictions[evaluation == 'Held cells' & calibration == 'Calibrated', {
  z <- copy(.SD)
  z[, drug := .BY$drug]
  summarize(z, 'Held cells')
}, by = .(family, bank, drug_set, drug)]
fwrite(per_drug, file.path(root, 'data/all_held_drug_results.tsv'), sep = '\t')
fwrite(per_drug[drug_set == 'Excluded family'], file.path(root, 'data/all_excluded_family_drug_results.tsv'), sep = '\t')
overview <- melt(summary[calibration == 'Calibrated' & evaluation == 'Held cells'],
  id.vars = c('family', 'bank', 'drug_set'),
  measure.vars = c('median_within_cell_r', 'across_cell_same_treatment_r'),
  variable.name = 'comparison_axis', value.name = 'r')
overview[, comparison_axis := factor(comparison_axis,
  levels = c('median_within_cell_r', 'across_cell_same_treatment_r'),
  labels = c('Across drugs/doses within held cells', 'Across held cells at the same drug/dose'))]
colours <- c('20 original / 40 programs' = '#666666', '43 current / 34 programs' = '#247F8A')
p_overview <- ggplot(overview, aes(r, drug_set, colour = bank)) +
  geom_vline(xintercept = 0, colour = 'grey70', linetype = 2) +
  geom_point(position = position_dodge(width = .4), size = 3) +
  facet_grid(family ~ comparison_axis) + scale_colour_manual(values = colours) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1)) +
  labs(title = 'Original calibrated SVR restored: only the program bank changes',
    subtitle = 'Same26 fitting cells, nine held cells, family exclusions, original settings and folds; all doses retained.',
    x = 'Pearson r', y = NULL, colour = NULL,
    caption = 'Bank43 discovery used full-panel outcomes. This is a retrospective decoder comparison; baseline-to-response is not fitted.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom')
ggsave(file.path(root, 'figures/restored_model_overview.pdf'), p_overview, width = 9, height = 5)
ggsave(file.path(root, 'figures/restored_model_overview.png'), p_overview, width = 9, height = 5, dpi = 160)

# Show fitting, cell validation and held data on separate panels. Plot a fixed
# sample of fitting rows for readability; all quantitative summaries use all rows.
scatter_data <- predictions[bank == '43 current / 34 programs' &
  (calibration == 'Calibrated' | evaluation == 'Whole-cell validation')]
scatter_data <- scatter_data[, .SD[unique(round(seq(1, .N, length.out = min(.N, 1500L))))],
  by = .(family, evaluation, drug_set)]
scatter_stats <- summary[bank == '43 current / 34 programs' &
  (calibration == 'Calibrated' | evaluation == 'Whole-cell validation')]
scatter_stats[, label := sprintf('r=%.2f | recovered %d/%d sensitive\nfalse sensitive %d/%d resistant',
  pooled_r, recovered_sensitive, sensitive, false_sensitive, resistant)]
p_scatter <- ggplot(scatter_data, aes(measured_pct, predicted_pct)) +
  geom_point(alpha = .25, size = .7, colour = '#247F8A') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#B24C35', linewidth = .5) +
  geom_vline(xintercept = 50, linetype = 2, colour = 'grey65') +
  geom_hline(yintercept = 50, linetype = 2, colour = 'grey65') +
  geom_text(data = scatter_stats, aes(label = label), x = -Inf, y = Inf,
    hjust = -.05, vjust = 1.1, size = 2.6, inherit.aes = FALSE) +
  facet_wrap(family + evaluation + drug_set ~ ., scales = 'free', ncol = 4) +
  scale_x_continuous(transform = 'log2') + scale_y_continuous(transform = 'log2') +
  labs(title = 'Current bank with the original SVR and calibration',
    subtitle = 'CV panels use raw out-of-fold predictions; calibration uses training OOF outcomes for fitting/held panels.',
    x = 'Measured % vehicle', y = 'Predicted % vehicle',
    caption = 'Panels display up to1500 observations; statistics use every row. Pearson correlations use log2 abundance.') +
  theme_classic(base_size = 9) + theme(aspect.ratio = 1, strip.text = element_text(size = 7))
ggsave(file.path(root, 'figures/restored_model_training_validation.pdf'), p_scatter, width = 11, height = 6.5)
ggsave(file.path(root, 'figures/restored_model_training_validation.png'), p_scatter, width = 11, height = 6.5, dpi = 160)

family_data <- predictions[evaluation == 'Held cells' & calibration == 'Calibrated' & drug_set == 'Excluded family']
family_stats <- per_drug[drug_set == 'Excluded family']
family_stats[, label := sprintf('r=%.2f\n%d/%d sensitive recovered', pooled_r, recovered_sensitive, sensitive)]
p_family <- ggplot(family_data, aes(measured_pct, predicted_pct, colour = factor(dose))) +
  geom_point(size = 1.1, alpha = .7) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, colour = 'grey35', linewidth = .4) +
  geom_vline(xintercept = 50, linetype = 2, colour = 'grey65') +
  geom_hline(yintercept = 50, linetype = 2, colour = 'grey65') +
  facet_grid(drug ~ bank, scales = 'free') +
  scale_x_continuous(transform = 'log2') + scale_y_continuous(transform = 'log2') +
  labs(title = 'Every assayed drug in the excluded families: same held cells, all doses',
    x = 'Measured % vehicle', y = 'Predicted % vehicle', colour = 'Dose (uM)') +
  theme_classic(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom')
family_plots <- list()
for (family_name in families) {
  pp <- (p_family %+% family_data[family == family_name]) +
    geom_text(data = family_stats[family == family_name], aes(label = label), x = -Inf, y = Inf,
      hjust = -.05, vjust = 1.1, size = 2.6, inherit.aes = FALSE)
  h <- if (family_name == 'ERBB') 15 else 7.5
  ggsave(file.path(root, 'figures', paste0(family_name, '_all_drug_restoration.pdf')), pp, width = 7, height = h)
  ggsave(file.path(root, 'figures', paste0(family_name, '_all_drug_restoration.png')), pp, width = 7, height = h, dpi = 160)
  family_plots[[family_name]] <- pp
}
grDevices::cairo_pdf(file.path(root, 'figures/consolidated_restoration_report.pdf'), width = 11, height = 8.5)
print(p_overview)
print(p_scatter)
grDevices::dev.off()
report <- c('Restored original calibrated displacement-only RBF SVR. Only bank changes.',
  'Start with figures/restored_model_overview.pdf and restored_model_training_validation.pdf.',
  'ERBB_all_drug_restoration.pdf and MEK_all_drug_restoration.pdf show all assayed excluded-family drugs.',
  'comparison_summary.tsv separates own-fit, raw whole-cell CV and held cells; seen drugs versus excluded family.',
  'Calibrated CV is deliberately not reported as independent validation.',
  'Original20 held raw/calibrated predictions replay exactly; original models were not refitted.',
  'Current43 discovery used full-panel outcomes. Results are retrospective and do not test baseline-to-response.',
  '', 'MATCHED HELD-CELL RESULTS (calibrated; same9 cells)',
  capture.output(print(summary[calibration == 'Calibrated' & evaluation == 'Held cells',
    .(family, bank, drug_set, median_within_cell_r, across_cell_same_treatment_r,
      recovered_sensitive, sensitive, false_sensitive, resistant)])),
  '', 'Interpretation',
  'The earlier familiar-drug signal reproduces: median within-held-cell r is',
  'ERBB0.658 in both banks; MEK0.611 original versus0.581 current. With the',
  'same learner, current43 improves familiar-treatment across-held-cell r:',
  'ERBB0.122->0.198; MEK0.189->0.258. This distinguishes bank effects from',
  'later changes to the decoder setting/selection/calibration.',
  'Magnitude remains compressed: familiar-drug sensitive observations have',
  'median measured viability about21%, versus current calibrated predictions',
  'about58%. Sensitive recovery increases slightly, with more false calls.',
  'Whole-cell CV same-treatment r remains0.062 ERBB/0.033 MEK (raw), so the',
  'specific9 held cells do not establish general transfer to arbitrary cells.',
  'Excluded-family results remain heterogeneous. Across-cell centeredr is',
  '0.119 ERBB and-0.110 MEK. ERBB sensitive recovery7/17, false calls14/172;',
  'MEK recovery16/47, false calls4/34. This does not establish the full bridge.',
  'All10 assayed excluded-family drugs are shown separately, retaining all',
  'doses, and exported in all_excluded_family_drug_results.tsv.')
writeLines(report, file.path(root, 'REPORT.txt'))
code_path <- '/work/src/program_pilot/50_restore_calibrated_svr_current_program_bank.R'
fwrite(data.table(path = code_path, sha256 = digest::digest(file = code_path, algo = 'sha256')),
  file.path(root, 'data/run_code_hash.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('COMPLETE: original formulation restored on current bank')
