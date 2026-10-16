#!/usr/bin/env Rscript
# Choose a smoother molecular-to-function SVR using whole-cell validation.
# Reuse the frozen, identically scored bank42/bank43 inputs from script44.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(patchwork)})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(20261010L)
options(warn = 1)

root <- '/work/results/program_pilot/45_whole_cell_regularized_program_svr'
input_path <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr/data/frozen_comparison_inputs.qs'
reference_path <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr/data/comparison_summary.tsv'
# Three folds each hold8–9 entire training cells, matching the9-cell test size.
# Lower C penalizes complex fits; lower gamma broadens the RBF kernel.
grid <- CJ(cost = c(.1, 1, 10), gamma_multiplier = c(.1, 1))
grid[, candidate := .I]
fold_count <- 3L
workers <- 3L
for (p in file.path(root, c('data', 'models', 'figures'))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading the frozen44 comparison; held-cell outcomes excluded from selection')
input <- qread(input_path)
md <- copy(input$metadata)
designs <- input$scores
train <- which(md$split == 'Training cells')
held <- which(md$split == 'Held cells')
check(uniqueN(md$depmap_id[train]) == 26L && uniqueN(md$depmap_id[held]) == 9L, 'Unexpected cell split')
check(!any(md$depmap_id[train] %chin% md$depmap_id[held]), 'Decoder split has overlapping cells')
folds <- data.table(depmap_id = sample(sort(unique(md$depmap_id[train]))),
  fold = rep(seq_len(fold_count), length.out = 26L))
md[, fold := folds$fold[match(depmap_id, folds$depmap_id)]]
check(all(!is.na(md$fold[train])) && all(is.na(md$fold[held])), 'Held cells entered validation folds')
fwrite(folds, file.path(root, 'data/whole_cell_folds.tsv'), sep = '\t')
fwrite(grid, file.path(root, 'data/regularization_grid.tsv'), sep = '\t')
signature <- digest::digest(list(input$input_hashes, input$held_ids, md$observation_id, folds, grid), algo = 'sha256')
fwrite(data.table(path = c(input_path, reference_path),
  sha256 = vapply(c(input_path, reference_path), digest::digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/input_hashes.tsv'), sep = '\t')

# All association measurements use log2(%vehicle/100), matching the fitted target.
# Center only for evaluation; every raw dose/outcome remains in fitting.
measure <- function(metadata, prediction) {
  z <- copy(metadata)
  z[, predicted_log2 := prediction]
  check(all(is.finite(z$predicted_log2)), 'Nonfinite predictions')
  z[, `:=`(yc = observed_log2 - mean(observed_log2), pc = predicted_log2 - mean(predicted_log2)),
    by = evaluation_context]
  safe_cor <- function(a, b) if (length(a) >= 3L && sd(a) > 1e-12 && sd(b) > 1e-12) cor(a, b) else NA_real_
  per_context <- z[, .(r = safe_cor(observed_log2, predicted_log2)), by = evaluation_context]
  pair_scores <- z[, {
    r <- which(measured_pct >= 60)
    s <- which(measured_pct <= 40)
    if (length(r) && length(s)) {
      valid <- outer(measured_pct[r], measured_pct[s], '-') >= 25
      gap <- outer(predicted_log2[r], predicted_log2[s], '-')
      accuracy <- if (any(valid)) mean((gap[valid] > 0) + .5 * (gap[valid] == 0)) else NA_real_
    } else accuracy <- NA_real_
    .(accuracy)
  }, by = evaluation_context]
  data.table(observations = nrow(z), cells = uniqueN(z$depmap_id), drugs = uniqueN(z$drug),
    sensitive = sum(z$measured_pct <= 50), resistant = sum(z$measured_pct > 50),
    pooled_pearson_log2 = safe_cor(z$observed_log2, z$predicted_log2),
    same_drug_dose_centered_pearson_log2 = safe_cor(z$yc, z$pc),
    median_per_drug_dose_pearson_log2 = median(per_context$r, na.rm = TRUE),
    sensitive_resistant_correct_order_fraction = mean(pair_scores$accuracy, na.rm = TRUE))
}
cv_metrics <- selected_settings <- final_metrics <- final_predictions <- list()

# Each candidate holds whole cells out. Fold checkpoints are flushed immediately;
# fork workers avoid PSOCK serialization and reuse read-only score matrices.
for (name in names(designs)) {
  x <- designs[[name]]
  bank_id <- substr(name, 1, 2)
  bank_metrics <- list()
  for (j in seq_len(nrow(grid))) {
    setting <- grid[j]
    status(paste('validating bank', bank_id, 'candidate', j, '/', nrow(grid),
      '| C', setting$cost, '| gamma', setting$gamma_multiplier, '/p'))
    results <- parallel::mclapply(seq_len(fold_count), function(fold_id) {
      RhpcBLASctl::blas_set_num_threads(1L)
      fit_rows <- train[md$fold[train] != fold_id]
      validation_rows <- train[md$fold[train] == fold_id]
      path <- file.path(root, 'models', sprintf('%s_candidate_%02d_fold_%d.qs', bank_id, j, fold_id))
      if (file.exists(path)) {
        saved <- qread(path)
        check(identical(saved$signature, signature), 'Fold checkpoint inputs changed')
        return(saved)
      }
      check(!any(md$depmap_id[fit_rows] %chin% md$depmap_id[validation_rows]) &&
        !any(md$depmap_id[c(fit_rows, validation_rows)] %chin% input$held_ids), 'Invalid validation cell split')
      fitted <- fit_program_svr(x[fit_rows, , drop = FALSE], md$observed_log2[fit_rows],
        'radial', setting$cost, setting$gamma_multiplier / ncol(x))
      predicted <- predict_program_svr(fitted, x[validation_rows, , drop = FALSE])
      saved <- list(signature = signature, fitted = fitted, setting = setting, bank = name,
        fitting_cells = sort(unique(md$depmap_id[fit_rows])), validation_cells = sort(unique(md$depmap_id[validation_rows])),
        validation_ids = md$observation_id[validation_rows], validation_prediction = predicted)
      qsave(saved, path, preset = 'high')
      saved
    }, mc.cores = workers, mc.preschedule = FALSE)
    check(all(vapply(results, function(z) is.list(z) && !inherits(z, 'try-error'), logical(1))), 'A CV worker failed')
    oof <- rep(NA_real_, length(train))
    for (saved in results) oof[match(saved$validation_ids, md$observation_id[train])] <- saved$validation_prediction
    check(all(is.finite(oof)), 'Incomplete whole-cell validation predictions')
    validation <- copy(md[train])
    validation[, evaluation_context := paste(condition, fold, sep = '__fold_')]
    score <- cbind(data.table(bank = name, candidate = j, cost = setting$cost,
      gamma_multiplier = setting$gamma_multiplier), measure(validation, oof))
    bank_metrics[[j]] <- score
    cv_metrics[[paste(bank_id, j)]] <- score
    qsave(list(signature = signature, validation_ids = validation$observation_id, prediction = oof),
      file.path(root, 'data', sprintf('%s_candidate_%02d_oof.qs', bank_id, j)), preset = 'high')
    fwrite(rbindlist(cv_metrics), file.path(root, 'data/whole_cell_validation_metrics.tsv'), sep = '\t')
    message('CV READY: bank', bank_id, ', C=', setting$cost, ', gamma=', setting$gamma_multiplier,
      '/p, same-treatment r=', round(score$same_drug_dose_centered_pearson_log2, 3))
  }
  # Selection rewards cell differences at fixed treatment, not pooled drug averages.
  candidates <- rbindlist(bank_metrics)
  check(all(is.finite(candidates$same_drug_dose_centered_pearson_log2)), 'Unscorable CV candidate')
  chosen <- candidates[order(-same_drug_dose_centered_pearson_log2, cost, gamma_multiplier)][1]
  selected_settings[[name]] <- chosen
  fwrite(rbindlist(selected_settings), file.path(root, 'data/selected_settings.tsv'), sep = '\t')
  status(paste('SELECTED bank', bank_id, '| C', chosen$cost, '| gamma', chosen$gamma_multiplier, '/p'))
  model_path <- file.path(root, 'models', paste0(bank_id, '_selected_svr.qs'))
  if (file.exists(model_path)) {
    saved <- qread(model_path)
    check(identical(saved$signature, signature) && identical(saved$chosen$candidate, chosen$candidate),
      'Selected decoder checkpoint changed')
    fitted <- saved$fitted
  } else {
    fitted <- fit_program_svr(x[train, , drop = FALSE], md$observed_log2[train],
      'radial', chosen$cost, chosen$gamma_multiplier / ncol(x))
    qsave(list(signature = signature, chosen = chosen, fitted = fitted,
      training_cells = sort(unique(md$depmap_id[train])), excluded_cells = input$held_ids), model_path, preset = 'high')
  }
  prediction <- predict_program_svr(fitted, x)
  outcome_free <- md[, .(observation_id, depmap_id, recipient, drug, dose, split)]
  outcome_free[, `:=`(bank = name, predicted_log2 = prediction)]
  fwrite(outcome_free, file.path(root, 'data', paste0(bank_id, '_outcome_free_predictions.tsv')), sep = '\t')
  fwrite(data.table(sha256 = digest::digest(outcome_free, algo = 'sha256')),
    file.path(root, 'data', paste0(bank_id, '_prediction_hash.tsv')), sep = '\t')
  oof <- qread(file.path(root, 'data', sprintf('%s_candidate_%02d_oof.qs', bank_id, chosen$candidate)))$prediction
  for (evaluation in c('Training fit', 'Whole-cell CV', 'Held cells')) {
    ii <- if (evaluation == 'Held cells') held else train
    z <- copy(md[ii])
    z[, evaluation_context := if (evaluation == 'Whole-cell CV') paste(condition, fold, sep = '__fold_') else condition]
    pred <- if (evaluation == 'Whole-cell CV') oof else prediction[ii]
    score <- cbind(data.table(bank = name, evaluation = evaluation, cost = chosen$cost,
      gamma_multiplier = chosen$gamma_multiplier), measure(z, pred))
    final_metrics[[paste(name, evaluation)]] <- score
    z[, `:=`(bank = name, evaluation = evaluation, predicted_log2 = pred)]
    z[, `:=`(observed_context = observed_log2 - mean(observed_log2),
      predicted_context = predicted_log2 - mean(predicted_log2)), by = evaluation_context]
    final_predictions[[paste(name, evaluation)]] <- z
  }
  fwrite(rbindlist(final_metrics), file.path(root, 'data/comparison_summary.tsv'), sep = '\t')
  print(rbindlist(final_metrics)[bank == name])
  status(paste('RESULTS SAVED: bank', bank_id))
}
z <- rbindlist(final_predictions)
fwrite(z, file.path(root, 'data/evaluated_predictions.tsv'), sep = '\t')
drug_summary <- z[, .(observations = .N, cells = uniqueN(depmap_id),
  pooled_pearson_log2 = cor(observed_log2, predicted_log2),
  same_drug_dose_centered_pearson_log2 = cor(observed_context, predicted_context)), by = .(bank, evaluation, drug)]
fwrite(drug_summary, file.path(root, 'data/per_drug_comparison.tsv'), sep = '\t')
reference <- fread(reference_path)[bank %chin% names(designs)]
fwrite(reference, file.path(root, 'data/44_flexible_svr_reference.tsv'), sep = '\t')
status('selected predictions and metrics saved; rendering figures')

for (view in c('pooled', 'same_drug_dose')) {
  plots <- list()
  for (name in names(designs)) for (evaluation_name in c('Training fit', 'Whole-cell CV', 'Held cells')) {
    d <- z[bank == name & evaluation == evaluation_name]
    d[, `:=`(actual = if (view == 'pooled') observed_log2 else observed_context,
      prediction = if (view == 'pooled') predicted_log2 else predicted_context)]
    plots[[paste(name, evaluation_name)]] <- ggplot(d, aes(actual, prediction)) +
      geom_point(size = .6, alpha = .18, colour = '#27818A') +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .6, colour = '#B24C35') +
      labs(title = paste(name, evaluation_name, sep = '\n'),
        subtitle = sprintf('Pearson r = %.2f | n = %d', cor(d$actual, d$prediction), nrow(d)),
        x = if (view == 'pooled') 'Measured log2(% vehicle / 100)' else 'Measured deviation within drug/dose (log2)',
        y = if (view == 'pooled') 'Predicted log2(% vehicle / 100)' else 'Predicted deviation within drug/dose (log2)') +
      theme_classic(base_size = 9) + theme(aspect.ratio = 1, plot.title = element_text(size = 10))
  }
  p <- wrap_plots(plots, ncol = 3) + plot_annotation(title = 'Regularized molecular-to-function SVR',
    subtitle = 'Complexity selected by3 whole-cell folds within26 training cells;9 test cells excluded from selection\nBoth banks used all cells during discovery: retrospective decoder comparison')
  ggsave(file.path(root, 'figures', paste0(view, '_prediction_comparison.pdf')), p, width = 12, height = 8)
  ggsave(file.path(root, 'figures', paste0(view, '_prediction_comparison.png')), p, width = 12, height = 8, dpi = 160)
}
curve <- rbindlist(cv_metrics)
p <- ggplot(curve, aes(cost, same_drug_dose_centered_pearson_log2, colour = factor(gamma_multiplier))) +
  geom_line() + geom_point(size = 2) + scale_x_log10(breaks = c(.1, 1, 10)) + facet_wrap(~bank, nrow = 1) +
  labs(title = 'Select regularization using whole-cell validation', x = 'SVR cost C',
    y = 'Validation r within drug/dose', colour = 'Gamma multiplier / p') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom')
ggsave(file.path(root, 'figures/regularization_selection.pdf'), p, width = 8, height = 4)
ggsave(file.path(root, 'figures/regularization_selection.png'), p, width = 8, height = 4, dpi = 160)
writeLines(c('Same frozen hard-mean molecular inputs and26training/9held cell split as44.',
  'RBF-SVR regularization selected by3 whole-cell folds within the26training cells.',
  'Grid C={0.1,1,10}; gamma={0.1,1}/feature_count; epsilon0.1; train-only scaling.',
  'Selection: highest out-of-fold Pearson r after centering within each drug/dose/fold.',
  'All doses remain independent and unadjusted in fitting. No baseline/drug/dose predictors.',
  'Training fit is in-sample; CV holds whole training cells out; held9cells are outside decoder selection.',
  'Program discovery used all cells, and earlier test results were seen: all results are retrospective.',
  'Original flexible44SVR metrics saved for comparison. Baseline-to-displacement not tested.'), file.path(root, 'README.txt'))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('complete')
