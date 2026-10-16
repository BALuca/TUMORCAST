#!/usr/bin/env Rscript
# Test nonlinear readout capacity while rewarding absolute and same-treatment sensitivity.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

family <- 'ERBB'
root <- '/work/results/program_pilot/38_nonlinear_absolute_and_recipient_sensitivity_decoder'
score_path <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps/data/ERBB_all_training_inputs.qs'
decoder_root <- '/work/results/program_pilot/33_balanced_cell_and_dose_contrast_decoder'
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout'
grid <- data.table(candidate = 1:3, cost = c(10, 100, 100), gamma_multiplier = c(1, 1, 4))
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading unchanged training programs and whole-cell folds')
input <- qread(score_path)
prepared_path <- file.path(parent, 'data/prepared_scores.qs')
prepared <- qread(prepared_path)
ridge_path <- file.path(decoder_root, 'models/ERBB_balanced_decoder.qs')
ridge <- qread(ridge_path)
fold_paths <- file.path(decoder_root, 'models', paste0('ERBB_fold_', 1:5, '.qs'))
fold_models <- lapply(fold_paths, qread)
md <- copy(input$annotations)
d <- input$displacement
programs <- colnames(d)
y <- log2(md$measured_pct / 100)
stopifnot(nrow(md) == 6260L, ncol(d) == 40L, identical(rownames(d), md$observation_id),
  identical(prepared$signature, ridge$source_signature), uniqueN(md$depmap_id) == 26L,
  !any(md$depmap_id %in% ridge$excluded_ids), !any(md$drug %in% ridge$excluded_drugs))

# Scaling uses independent fitting recipients, rather than their repeated drug rows.
learn_scaling <- function(ids) {
  b <- prepared$baseline[sort(unique(ids)), , drop = FALSE]
  s <- list(mean = colMeans(b), sd = apply(b, 2L, sd), ids = sort(unique(ids)))
  stopifnot(all(s$sd > 0), all(is.finite(s$sd)))
  s
}
make_design <- function(ids, displacement, scaling) {
  b <- sweep(sweep(prepared$baseline[ids, , drop = FALSE], 2L, scaling$mean, '-'), 2L, scaling$sd, '/')
  x <- cbind(displacement, displacement * b)
  colnames(x) <- c(paste0('D_', programs), paste0('DB_', programs))
  stopifnot(ncol(x) == 80L, all(is.finite(x)))
  x
}

# Same fixed threshold for every model. Pair scores average settings equally,
# not the dependent cross-products of endpoints as independent samples.
measure <- function(metadata, prediction, context_keys) {
  z <- copy(metadata)
  z[, `:=`(observed = log2(measured_pct / 100), predicted = prediction, predicted_pct = 100 * 2^prediction)]
  stopifnot(all(is.finite(z$predicted_pct)))
  cell <- z[condition %in% context_keys]
  cell[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)), by = condition]
  pair_scores <- cell[, {
    r <- which(measured_pct >= 60)
    s <- which(measured_pct <= 40)
    if (!length(r) || !length(s)) list(pairs = 0L, correct = NA_real_) else {
      pp <- CJ(r = r, s = s)
      pp <- pp[measured_pct[r] - measured_pct[s] >= 25]
      list(pairs = nrow(pp), correct = if (nrow(pp)) mean(predicted[pp$r] > predicted[pp$s]) else NA_real_)
    }
  }, by = condition]
  n_sensitive <- sum(z$measured_pct <= 50)
  n_resistant <- sum(z$measured_pct > 50)
  sensitive_called <- sum(z$measured_pct <= 50 & z$predicted_pct <= 50)
  resistant_called <- sum(z$measured_pct > 50 & z$predicted_pct > 50)
  balanced_accuracy <- (sensitive_called / n_sensitive + resistant_called / n_resistant) / 2
  pair_accuracy <- mean(pair_scores$correct, na.rm = TRUE)
  data.table(observations = nrow(z), cells = uniqueN(z$depmap_id), global_r = cor_safe(z$observed, z$predicted),
    fixed_treatment_observations = nrow(cell), fixed_treatment_r = cor_safe(cell$yc, cell$pc),
    sensitive = n_sensitive, sensitive_called = sensitive_called, resistant = n_resistant,
    resistant_called = resistant_called, false_sensitive = n_resistant - resistant_called,
    balanced_accuracy_pct = 100 * balanced_accuracy, clear_cell_pair_order_pct = 100 * pair_accuracy,
    pair_contexts = sum(is.finite(pair_scores$correct)), cell_pairs = sum(pair_scores$pairs),
    selection_reward = (balanced_accuracy + pair_accuracy) / 2,
    measured_sensitive_median_pct = median(z$measured_pct[z$measured_pct <= 50]),
    measured_resistant_median_pct = median(z$measured_pct[z$measured_pct > 50]),
    predicted_sensitive_median_pct = median(z$predicted_pct[z$measured_pct <= 50]),
    predicted_resistant_median_pct = median(z$predicted_pct[z$measured_pct > 50]))
}
full_scaling <- learn_scaling(md$depmap_id)
x <- make_design(md$depmap_id, d, full_scaling)
stopifnot(max(abs(predict_ridge(ridge$model, x) - log2(md$predicted_pct / 100))) < 1e-10)
sources <- c(score_path, prepared_path, ridge_path, fold_paths,
  '/work/src/program_pilot/38_nonlinear_absolute_and_recipient_sensitivity_decoder.R',
  '/work/src/program_pilot/response_models.R', '/work/src/program_pilot/nonlinear_response.R')
manifest <- data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256'))
signature <- digest::digest(list(manifest, grid), algo = 'sha256')
settings_path <- file.path(root, 'settings.qs')
if (file.exists(settings_path)) {
  # Resume models after reporting fixes only when all fitting data, helpers,
  # folds and settings match the originally recorded run. The script revision
  # is recorded separately; no fitting contract changed in the report repair.
  previous <- qread(settings_path)
  previous_manifest <- fread(file.path(root, 'data', 'input_hashes.tsv'))
  fitting_paths <- sources[sources != '/work/src/program_pilot/38_nonlinear_absolute_and_recipient_sensitivity_decoder.R']
  stopifnot(isTRUE(all.equal(as.data.frame(previous$candidates), as.data.frame(grid))), identical(previous$training_ids, ridge$training_ids),
    identical(previous_manifest[match(fitting_paths, path), sha256], manifest[match(fitting_paths, path), sha256]))
  signature <- previous$signature
}
save_table(manifest, 'input_hashes')
qsave(list(signature = signature, family = family, candidates = grid,
  training_ids = ridge$training_ids, excluded_ids = ridge$excluded_ids, excluded_drugs = ridge$excluded_drugs,
  input = '40 D and 40 D times corresponding standardized baseline program',
  training_objective = 'absolute log2 abundance; epsilon-insensitive RBF SVR',
  selection = '50% balanced sensitive/resistant accuracy + 50% same-treatment clear-pair ordering',
  anchoring = 'native fitted absolute predictions; no post-fit zero subtraction'), file.path(root, 'settings.qs'))

training_metrics <- list(cbind(data.table(model = 'Current ridge', candidate = 0L),
  measure(md, log2(md$predicted_pct / 100), ridge$screen$cell_keys)))
training_predictions <- list(cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
  data.table(model = 'Current ridge', candidate = 0L, predicted = log2(md$predicted_pct / 100))))
oof <- matrix(NA_real_, nrow(md), nrow(grid))
oof_ridge <- rep(NA_real_, nrow(md))
cv_context <- rep(FALSE, nrow(md))
cv_fold <- rep(NA_integer_, nrow(md))
validation <- list()
final_models <- vector('list', nrow(grid))
for (j in seq_len(nrow(grid))) {
  spec <- grid[j]
  status(paste('candidate', j, '| fitting full training data; C=', spec$cost, 'gamma multiplier=', spec$gamma_multiplier))
  model_path <- file.path(root, 'models', paste0('candidate_', j, '_full.qs'))
  if (file.exists(model_path)) {
    cached <- qread(model_path)
    stopifnot(identical(cached$signature, signature), identical(cached$scaling, full_scaling),
      isTRUE(all.equal(as.data.frame(cached$selected_setting), as.data.frame(spec))), identical(cached$source_signature, prepared$signature))
    fit <- cached$model
  } else {
    fit <- fit_program_svr(x, y, 'radial', spec$cost, spec$gamma_multiplier / ncol(x))
    qsave(list(signature = signature, model = fit, scaling = full_scaling, selected_setting = spec,
      training_ids = ridge$training_ids, source_signature = prepared$signature), model_path)
  }
  final_models[[j]] <- fit
  prediction <- predict_program_svr(fit, x)
  zero_pct <- 100 * 2^predict_program_svr(fit, matrix(0, 1L, ncol(x)))
  training_metrics[[j + 1L]] <- cbind(data.table(model = paste0('SVR candidate ', j), candidate = j),
    measure(md, prediction, ridge$screen$cell_keys), data.table(zero_displacement_pct_vehicle = zero_pct))
  training_predictions[[j + 1L]] <- cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(model = paste0('SVR candidate ', j), candidate = j, predicted = prediction))
  save_table(rbindlist(training_metrics, fill = TRUE), 'training_metrics')
  save_table(rbindlist(training_predictions), 'training_predictions')
  status(paste('candidate', j, '| training fixed-treatment r=',
    round(training_metrics[[j + 1L]]$fixed_treatment_r, 3), '| sensitive called=',
    training_metrics[[j + 1L]]$sensitive_called, '/', training_metrics[[j + 1L]]$sensitive))

  for (k in 1:5) {
    saved <- fold_models[[k]]
    tr <- which(md$depmap_id %in% saved$training_ids)
    va <- which(md$depmap_id %in% saved$validation_ids)
    stopifnot(!length(intersect(md$depmap_id[tr], md$depmap_id[va])), length(tr) + length(va) == nrow(md),
      identical(saved$source_signature, prepared$signature))
    scaling <- learn_scaling(md$depmap_id[tr])
    xx <- make_design(md$depmap_id[tr], d[tr, , drop = FALSE], scaling)
    vx <- make_design(md$depmap_id[va], d[va, , drop = FALSE], scaling)
    status(paste('candidate', j, '| whole-cell selection fold', k, '/5'))
    checkpoint_path <- file.path(root, 'models', paste0('candidate_', j, '_fold_', k, '.qs'))
    if (file.exists(checkpoint_path)) {
      cached <- qread(checkpoint_path)
      stopifnot(identical(cached$signature, signature), identical(cached$scaling, scaling),
        isTRUE(all.equal(as.data.frame(cached$setting), as.data.frame(spec))), identical(cached$validation_rows, va),
        identical(cached$training_ids, saved$training_ids), identical(cached$validation_ids, saved$validation_ids))
      pp <- cached$prediction
    } else {
      ff <- fit_program_svr(xx, y[tr], 'radial', spec$cost, spec$gamma_multiplier / ncol(xx))
      pp <- predict_program_svr(ff, vx)
      # Commit the fitted object and predictions immediately after each fold.
      qsave(list(signature = signature, model = ff, scaling = scaling, setting = spec,
        training_ids = saved$training_ids, validation_ids = saved$validation_ids,
        validation_rows = va, prediction = pp, screen = saved$screen), checkpoint_path)
    }
    oof[va, j] <- pp
    cv_context[va] <- md$condition[va] %in% saved$screen$cell_keys
    cv_fold[va] <- k
    validation[[paste(j, k)]] <- cbind(data.table(candidate = j, fold = k),
      measure(md[va], pp, saved$screen$cell_keys))
    if (j == 1L) {
      rr <- saved$models$balanced
      rr$lambda <- ridge$model$lambda
      oof_ridge[va] <- predict_ridge(rr, x[va, , drop = FALSE])
    }
    save_table(rbindlist(validation), 'whole_cell_validation_metrics')
    qsave(list(oof = oof, ridge = oof_ridge, eligible = cv_context, observation_id = md$observation_id),
      file.path(root, 'data', 'development_prediction_checkpoint.qs'))
  }
}
stopifnot(all(is.finite(oof)), all(is.finite(oof_ridge)), !anyNA(cv_fold))
selection <- rbindlist(validation)[, .(selection_reward = mean(selection_reward),
  balanced_accuracy_pct = mean(balanced_accuracy_pct), clear_cell_pair_order_pct = mean(clear_cell_pair_order_pct)),
  by = candidate]
setorder(selection, -selection_reward, candidate)
best <- selection$candidate[1L]
selection[, selected := candidate == best]
save_table(selection, 'selected_setting')
selected <- qread(file.path(root, 'models', paste0('candidate_', best, '_full.qs')))
qsave(selected, file.path(root, 'models', 'selected_decoder.qs'))
status(paste('selected candidate', best, '| freezing original held-cell predictions'))

# Predict before loading held outcomes. No drug/cell identity enters the readout.
settings <- qread(file.path(parent, 'settings.qs'))
query <- copy(as.data.table(prepared$metadata)[depmap_id %in% ridge$excluded_ids & drug %in% md$drug])
query[, condition := paste(drug, dose, sep = '__')]
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
cell_lookup <- as.data.table(qread(panel_path)$cell_map)
query[, recipient := cell_lookup$depmap_name[match(depmap_id, cell_lookup$depmap_id)]]
stopifnot(!anyNA(query$recipient))
hd <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
hx <- make_design(query$depmap_id, hd, full_scaling)
query[, `:=`(svr_prediction = predict_program_svr(selected$model, hx), ridge_prediction = predict_ridge(ridge$model, hx))]
outcome_free <- query[, .(observation_id, depmap_id, recipient, drug, dose, condition, svr_prediction, ridge_prediction)]
save_table(outcome_free, 'outcome_free_held_predictions')
prediction_path <- file.path(root, 'data', 'outcome_free_held_predictions.tsv')
save_table(data.table(path = prediction_path, sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')

labels <- as.data.table(qread(panel_path)$sensitivity)
query[, measured_pct := labels$prism_abundance_pct[match(observation_id, labels$observation_id)]]
save_table(query[!is.finite(measured_pct) | measured_pct <= 0,
  .(observation_id, depmap_id, drug, dose, measured_pct)], 'unsupported_held_endpoints')
query <- query[is.finite(measured_pct) & measured_pct > 0]
stopifnot(uniqueN(query$depmap_id) == 9L, !any(query$depmap_id %in% md$depmap_id),
  !any(query$drug %in% ridge$excluded_drugs))
tab <- rbindlist(list(
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Training fit', model = 'Current ridge', predicted = log2(md$predicted_pct / 100),
      cell_context = md$condition %in% ridge$screen$cell_keys, fold = 0L)),
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Training fit', model = 'Selected SVR', predicted = predict_program_svr(selected$model, x),
      cell_context = md$condition %in% ridge$screen$cell_keys, fold = 0L)),
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Development cell folds', model = 'Current ridge', predicted = oof_ridge,
      cell_context = cv_context, fold = cv_fold)),
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Development cell folds', model = 'Selected SVR', predicted = oof[, best],
      cell_context = cv_context, fold = cv_fold)),
  cbind(query[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Original held cells', model = 'Current ridge', predicted = query$ridge_prediction,
      cell_context = query$condition %in% ridge$screen$cell_keys, fold = 0L)),
  cbind(query[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Original held cells', model = 'Selected SVR', predicted = query$svr_prediction,
      cell_context = query$condition %in% ridge$screen$cell_keys, fold = 0L))))
tab[, `:=`(observed = log2(measured_pct / 100), predicted_pct = 100 * 2^predicted)]
tab[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)), by = .(cohort, model, fold, condition)]
save_table(tab, 'evaluated_predictions')
summary <- tab[, .(observations = .N, cells = uniqueN(depmap_id), global_r = cor_safe(observed, predicted),
  fixed_treatment_observations = sum(cell_context), fixed_treatment_r = cor_safe(yc[cell_context], pc[cell_context]),
  sensitive = sum(measured_pct <= 50), sensitive_called = sum(measured_pct <= 50 & predicted_pct <= 50),
  resistant = sum(measured_pct > 50), resistant_called = sum(measured_pct > 50 & predicted_pct > 50),
  false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50),
  measured_sensitive_median_pct = median(measured_pct[measured_pct <= 50]),
  measured_resistant_median_pct = median(measured_pct[measured_pct > 50]),
  predicted_sensitive_median_pct = median(predicted_pct[measured_pct <= 50]),
  predicted_resistant_median_pct = median(predicted_pct[measured_pct > 50])), by = .(cohort, model)]
save_table(summary, 'comparison_summary')
save_table(tab[cell_context == TRUE, .(observations = .N, cells = uniqueN(depmap_id),
  r = cor_safe(observed, predicted)), by = .(cohort, model, fold, drug, dose)], 'per_treatment_correlations')
stopifnot(max(abs(predict_program_svr(qread(file.path(root, 'models', 'selected_decoder.qs'))$model, x) -
  tab[cohort == 'Training fit' & model == 'Selected SVR', predicted])) < 1e-10)
save_table(data.table(candidate = best,
  zero_displacement_pct_vehicle = 100 * 2^predict_program_svr(selected$model, matrix(0, 1L, ncol(x)))), 'zero_input_audit')

status('exporting compact readout diagnostics')
colors <- c('Current ridge' = '#777777', 'Selected SVR' = '#0072B2')
summary[, cohort := factor(cohort, levels = c('Training fit', 'Development cell folds', 'Original held cells'))]
metrics <- melt(summary, id.vars = c('cohort', 'model'), measure.vars = c('global_r', 'fixed_treatment_r'),
  variable.name = 'axis', value.name = 'r')
metrics[, axis := factor(axis, levels = c('global_r', 'fixed_treatment_r'),
  labels = c('All treatments combined', 'Cells at fixed drug/dose'))]
p <- ggplot(metrics, aes(model, r, fill = model)) + geom_col(width = .65) +
  geom_text(aes(label = sprintf('%.2f', r)), vjust = -.4, size = 3) +
  facet_grid(axis ~ cohort) + scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(min(0, min(metrics$r) - .05), 1), expand = expansion(mult = c(.02, .06))) +
  theme_classic(base_size = 10) + theme(legend.position = 'none', axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = NULL, y = 'Pearson r (log2 abundance)', title = 'Same programs: linear versus nonlinear response decoder',
    caption = 'Development folds validate the readout conditional on the outcome-informed bank; original held results are retrospective.')
grDevices::cairo_pdf(file.path(root, 'figures', 'readout_correlations.pdf'), width = 10, height = 6)
print(p)
dev.off()
ggsave(file.path(root, 'figures', 'readout_correlations.png'), p, width = 10, height = 6, dpi = 180)
strip <- tab[, .(observation_id, cohort, model, actual_class = ifelse(measured_pct <= 50, 'Measured sensitive', 'Measured resistant'),
  measured_pct, predicted_pct)]
strip[, cohort := factor(cohort, levels = levels(summary$cohort))]
p <- ggplot(strip, aes(actual_class, predicted_pct, fill = model)) +
  geom_boxplot(position = position_dodge(.8), width = .65, outlier.shape = NA) +
  geom_hline(yintercept = 50, linetype = 2) + facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = colors) + scale_y_log10() + theme_classic(base_size = 10) +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = 'Predicted abundance (% vehicle; log axis)', fill = NULL,
    title = 'Does the decoder separate measured sensitive and resistant observations?',
    caption = 'Dashed line: 50% vehicle cutoff. Boxes show median/IQR; outlying points omitted for readability, not removed from metrics.')
grDevices::cairo_pdf(file.path(root, 'figures', 'sensitive_resistant_separation.pdf'), width = 10, height = 4.5)
print(p)
dev.off()
ggsave(file.path(root, 'figures', 'sensitive_resistant_separation.png'), p, width = 10, height = 4.5, dpi = 180)
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(c('NONLINEAR DECODER WITH ABSOLUTE AND RECIPIENT SENSITIVITY SELECTION', '',
  'Same frozen ERBB 40-program bank and 6,260 training observations; no program discovery.',
  '80 features: D plus D times corresponding standardized baseline program. No identity predictors.',
  'Absolute log2 abundance SVR; three fixed capacities. Whole-cell selection rewards correct class calls and clear-pair ordering equally.',
  'Model objects and fold predictions are saved immediately as produced.',
  'Training fit, readout development folds and original nine held recipients are reported separately.',
  'Development folds condition on a bank learned from all development cells; they do not independently validate program discovery.',
  'Original held results are retrospective. Selection/predictions are frozen before reading their outcomes.',
  'No post-fit zero subtraction. zero_input_audit.tsv reports any failure to predict 100% at zero displacement.',
  'This tests the measured displacement-to-function readout, not baseline-to-displacement accuracy or combinations.'),
  file.path(root, 'READING_GUIDE.txt'))
print(summary)
status('COMPLETE')
