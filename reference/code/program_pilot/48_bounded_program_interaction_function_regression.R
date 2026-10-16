#!/usr/bin/env Rscript
# Recover a bounded surviving fraction with regularized program-pair effects.
# Reuse fixed43 scores and exact46 family/cell splits; preserve all dose rows.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(20261010L)
options(warn = 1)

root <- '/work/results/program_pilot/48_bounded_program_interaction_function_regression'
reference <- '/work/results/program_pilot/46_family_holdout_context_signal_pilot'
input_path <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr/data/frozen_comparison_inputs.qs'
fold_path <- '/work/results/program_pilot/45_whole_cell_regularized_program_svr/data/whole_cell_folds.tsv'
lambda_grid <- 10^seq(1, -6, length.out = 22L)
interaction_penalty <- 10
workers <- 3L
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading unchanged43 program scores and46 family projections')
input <- qread(input_path)
md <- copy(input$metadata)
x <- input$scores[['43 averaged / 34 programs']]
projection <- qread(file.path(reference, 'data/family_projection_inputs.qs'))
fm <- copy(projection$metadata)
fx <- projection$scores
family_map <- fread(file.path(reference, 'data/excluded_families.tsv'))
folds <- fread(fold_path)
md[, fold := folds$fold[match(depmap_id, folds$depmap_id)]]
check(identical(rownames(x), md$observation_id) && identical(colnames(x), colnames(fx)), 'Score identifiers differ')
check(ncol(x) == 34L && all(is.finite(x)) && all(is.finite(fx)), 'Invalid fixed program inputs')
pairs <- which(upper.tri(matrix(0, ncol(x), ncol(x)), diag = TRUE), arr.ind = TRUE)
feature_names <- c(colnames(x), paste(colnames(x)[pairs[, 1L]], colnames(x)[pairs[, 2L]], sep = ':'))
penalty_factors <- c(rep(1, ncol(x)), rep(interaction_penalty, nrow(pairs)))
check(length(feature_names) == 629L, 'Unexpected program-interaction design')
sources <- c(input_path, fold_path, file.path(reference, 'data/family_projection_inputs.qs'),
  file.path(reference, 'data/excluded_families.tsv'),
  file.path(reference, 'models', paste0(c('ERBB', 'MEK'), '_selected_decoder.qs')))
hashes <- data.table(path = sources,
  sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256'))
signature <- digest::digest(list(hashes, lambda_grid, interaction_penalty, pairs,
  'fractional-binomial-target; capped-vehicle-MAE; minimum-whole-cell-CV'), algo = 'sha256')
fwrite(hashes, file.path(root, 'data/input_hashes.tsv'), sep = '\t')
fwrite(data.table(feature = feature_names, penalty_factor = penalty_factors),
  file.path(root, 'data/feature_definition.tsv'), sep = '\t')
qsave(list(signature = signature, lambda_grid = lambda_grid, interaction_penalty = interaction_penalty,
  folds = folds, held_cells = input$held_ids, family_map = family_map,
  score_columns = colnames(x), target = 'min(%vehicle/100,1), continuous surviving fraction',
  selection = 'minimum drug/class-balanced absolute error in capped %vehicle'), file.path(root, 'settings.qs'))

# One design is used for CV, full fits and future molecular queries.
# Program products let coherent combinations affect function; ridge constrains them.
make_design <- function(scores, scaling) {
  z <- sweep(sweep(scores, 2L, scaling$mean, '-'), 2L, scaling$sd, '/')
  design <- cbind(z, z[, pairs[, 1L], drop = FALSE] * z[, pairs[, 2L], drop = FALSE])
  colnames(design) <- feature_names
  check(all(is.finite(design)), 'Nonfinite program-interaction design')
  design
}
fit_decoder <- function(scores, metadata) {
  scaling <- list(mean = colMeans(scores), sd = apply(scores, 2L, sd))
  check(all(scaling$sd > 1e-12), 'Constant fitting program')
  design <- make_design(scores, scaling)
  # Give each drug equal total weight; all its cells and doses remain observations.
  weights <- 1 / as.numeric(table(metadata$drug)[metadata$drug])
  # Fractional binomial regression uses the full continuous survival response,
  # not a sensitive/resistant training label. Above-vehicle growth is no killing.
  survived <- pmin(metadata$measured_pct / 100, 1)
  response <- cbind(1 - survived, survived)
  fit <- glmnet(design, response, alpha = 0, family = 'binomial',
    lambda = lambda_grid, weights = weights, standardize = TRUE, intercept = TRUE,
    penalty.factor = penalty_factors, thresh = 1e-7, maxit = 100000L)
  check(length(fit$lambda) == length(lambda_grid), 'Incomplete ridge path')
  list(fit = fit, scaling = scaling)
}
predict_path <- function(model, scores, penalty = lambda_grid) {
  fraction <- as.matrix(predict(model$fit, newx = make_design(scores, model$scaling),
    s = penalty, type = 'response'))
  # Numerical floor only, to represent saturated logistic predictions on log axes.
  p <- log2(pmax(fraction, .Machine$double.xmin))
  check(all(is.finite(p)), 'Nonfinite log2 predictions')
  p
}
safe_cor <- function(a, b) {
  if (length(a) >= 3L && sd(a) > 1e-12 && sd(b) > 1e-12) cor(a, b) else NA_real_
}
validation_error <- function(metadata, prediction) {
  z <- metadata[, .(drug, measured_pct)]
  z[, `:=`(sensitive = measured_pct <= 50, error_pct = abs(pmin(measured_pct, 100) - 100 * 2^prediction))]
  check(all(is.finite(z$error_pct)), 'Nonfinite assay-unit validation errors')
  per_drug_class <- z[, .(error_pct = mean(error_pct)), by = .(drug, sensitive)]
  mean(per_drug_class[, .(error_pct = mean(error_pct)), by = sensitive]$error_pct)
}
summarize <- function(metadata, prediction) {
  z <- copy(metadata)
  z[, `:=`(observed_log2 = log2(measured_pct / 100), predicted_log2 = prediction,
    predicted_pct = 100 * 2^prediction)]
  z[, `:=`(yc = observed_log2 - mean(observed_log2), pc = predicted_log2 - mean(predicted_log2)),
    by = .(drug, dose)]
  s <- z$measured_pct <= 50
  data.table(observations = nrow(z), cells = uniqueN(z$depmap_id), drugs = uniqueN(z$drug),
    sensitive = sum(s), resistant = sum(!s), sensitive_recovered = sum(s & z$predicted_pct <= 50),
    false_sensitive = sum(!s & z$predicted_pct <= 50),
    sensitive_recall = if (any(s)) mean(z$predicted_pct[s] <= 50) else NA_real_,
    resistant_specificity = if (any(!s)) mean(z$predicted_pct[!s] > 50) else NA_real_,
    pooled_r = safe_cor(z$observed_log2, prediction), within_drug_dose_r = safe_cor(z$yc, z$pc),
    measured_sensitive_median_pct = median(z$measured_pct[s]),
    predicted_sensitive_median_pct = median(z$predicted_pct[s]),
    measured_resistant_median_pct = median(z$measured_pct[!s]),
    predicted_resistant_median_pct = median(z$predicted_pct[!s]))
}

all_predictions <- cv_tables <- selected_settings <- acceptance <- list()
for (family_name in c('ERBB', 'MEK')) {
  ref <- qread(file.path(reference, 'models', paste0(family_name, '_selected_decoder.qs')))
  train <- match(ref$training_ids, md$observation_id)
  excluded <- family_map[family == family_name, drug]
  check(!anyNA(train) && all(!is.na(md$fold[train])) && uniqueN(md$depmap_id[train]) == 26L,
    'Invalid training split')
  check(!any(md$depmap_id[train] %chin% input$held_ids) && !any(md$drug[train] %chin% excluded),
    'Excluded family or held cell entered fitting')
  status(paste(family_name, 'fitting three whole-cell ridge paths'))
  fold_results <- parallel::mclapply(1:3, function(k) {
    tr <- train[md$fold[train] != k]
    va <- train[md$fold[train] == k]
    path <- file.path(root, 'models', sprintf('%s_fold_%d.qs', family_name, k))
    if (file.exists(path)) {
      saved <- qread(path)
      check(identical(saved$signature, signature), 'Stale fold checkpoint')
      return(saved)
    }
    model <- fit_decoder(x[tr, , drop = FALSE], md[tr])
    prediction <- predict_path(model, x[va, , drop = FALSE])
    saved <- list(signature = signature, model = model, prediction = prediction,
      fitting_ids = md$observation_id[tr], validation_ids = md$observation_id[va],
      fitting_cells = unique(md$depmap_id[tr]), validation_cells = unique(md$depmap_id[va]),
      excluded_drugs = excluded, excluded_cells = input$held_ids)
    qsave(saved, path, preset = 'high')
    message(family_name, ' FOLD ', k, ' SAVED')
    saved
  }, mc.cores = workers, mc.preschedule = FALSE)
  check(all(vapply(fold_results, function(z) is.list(z) && !inherits(z, 'try-error'), logical(1))), 'Fold worker failed')
  oof <- matrix(NA_real_, length(train), length(lambda_grid))
  reference_oof <- rep(NA_real_, length(train))
  fold_errors <- matrix(NA_real_, 3, length(lambda_grid))
  for (k in 1:3) {
    saved <- fold_results[[k]]
    ii <- match(saved$validation_ids, md$observation_id[train])
    check(!any(saved$fitting_cells %chin% saved$validation_cells), 'Overlapping CV recipients')
    oof[ii, ] <- saved$prediction
    for (j in seq_along(lambda_grid)) fold_errors[k, j] <- validation_error(md[train[ii]], saved$prediction[, j])
    old <- qread(file.path(reference, 'models', sprintf('%s_candidate_%02d_fold_%d.qs',
      family_name, ref$chosen$candidate, k)))
    check(setequal(old$validation_ids, saved$validation_ids), 'Reference validation recipients differ')
    reference_oof[match(old$validation_ids, md$observation_id[train])] <- old$prediction
  }
  check(all(is.finite(oof)) && all(is.finite(reference_oof)), 'Incomplete development predictions')
  cv <- data.table(family = family_name, lambda = lambda_grid, mean_error_pct = colMeans(fold_errors),
    se_error_pct = apply(fold_errors, 2, sd) / sqrt(3))
  best <- which.min(cv$mean_error_pct)
  chosen <- best
  cv[, selected := .I == chosen]
  cv_tables[[family_name]] <- cv
  selected_settings[[family_name]] <- cv[chosen]
  fwrite(rbindlist(cv_tables), file.path(root, 'data/whole_cell_validation.tsv'), sep = '\t')
  fwrite(rbindlist(selected_settings), file.path(root, 'data/selected_settings.tsv'), sep = '\t')
  status(paste(family_name, 'selected penalty', lambda_grid[chosen], '| fitting full decoder'))
  path <- file.path(root, 'models', paste0(family_name, '_selected_decoder.qs'))
  if (file.exists(path)) {
    full <- qread(path)
    check(identical(full$signature, signature) && identical(full$lambda, lambda_grid[chosen]), 'Stale full fit')
  } else {
    full <- list(signature = signature, model = fit_decoder(x[train, , drop = FALSE], md[train]),
      lambda = lambda_grid[chosen], fitting_ids = md$observation_id[train],
      fitting_cells = unique(md$depmap_id[train]), excluded_drugs = excluded, excluded_cells = input$held_ids)
    qsave(full, path, preset = 'high')
  }
  coef_table <- data.table(feature = c('(Intercept)', feature_names),
    coefficient = as.numeric(coef(full$model$fit, s = full$lambda)))
  fwrite(coef_table, file.path(root, 'data', paste0(family_name, '_coefficients.tsv')), sep = '\t')
  for (model_name in c('46 SVR', '48 Program interactions')) {
    for (stage in c('Training fit', 'Whole-cell validation')) {
      d <- copy(md[train, .(observation_id, depmap_id, recipient, drug, dose, measured_pct)])
      pred <- if (model_name == '46 SVR') {
        if (stage == 'Training fit') predict_program_svr(ref$fitted, x[train, , drop = FALSE]) else reference_oof
      } else {
        if (stage == 'Training fit') drop(predict_path(full$model, x[train, , drop = FALSE], full$lambda)) else oof[, chosen]
      }
      d[, `:=`(family = family_name, model = model_name, evaluation = stage, predicted_log2 = pred)]
      all_predictions[[paste(family_name, model_name, stage)]] <- d
    }
  }
  ii <- which(fm$family == family_name)
  pred <- drop(predict_path(full$model, fx[ii, , drop = FALSE], full$lambda))
  outcome_free <- fm[ii, .(observation_id, depmap_id, recipient, drug, dose, recipient_set)]
  outcome_free[, `:=`(family = family_name, predicted_log2 = pred)]
  fwrite(outcome_free, file.path(root, 'data', paste0(family_name, '_outcome_free_predictions.tsv')), sep = '\t')
  fwrite(data.table(sha256 = digest::digest(outcome_free, algo = 'sha256')),
    file.path(root, 'data', paste0(family_name, '_prediction_hash.tsv')), sep = '\t')
  old <- fread(file.path(reference, 'data', paste0(family_name, '_outcome_free_predictions.tsv')))
  check(identical(outcome_free$observation_id, old$observation_id), 'Family prediction rows differ')
  for (model_name in c('46 SVR', '48 Program interactions')) {
    d <- copy(outcome_free)
    d[, `:=`(measured_pct = fm$measured_pct[ii], model = model_name,
      predicted_log2 = if (model_name == '46 SVR') old$predicted_log2 else pred,
      evaluation = fifelse(recipient_set == 'Existing cells', 'Excluded family / familiar cells',
        fifelse(recipient_set == 'Held cells', 'Excluded family / held cells', 'Excluded family / no assay cells')))]
    all_predictions[[paste(family_name, model_name, 'family')]] <- d
  }
  interim <- rbindlist(all_predictions, fill = TRUE)
  fwrite(interim, file.path(root, 'data/evaluated_predictions.tsv'), sep = '\t')
  summary <- interim[is.finite(measured_pct) & measured_pct > 0,
    summarize(.SD, predicted_log2), by = .(family, model, evaluation)]
  fwrite(summary, file.path(root, 'data/summary.tsv'), sep = '\t')
  print(summary[family == family_name])
  acceptance[[family_name]] <- data.table(family = family_name, fitting_cells = uniqueN(md$depmap_id[train]),
    fold_models = 3L, held_cells_absent = !any(md$depmap_id[train] %chin% input$held_ids),
    family_absent = !any(md$drug[train] %chin% excluded), family_columns = length(ii))
  fwrite(rbindlist(acceptance), file.path(root, 'data/split_acceptance.tsv'), sep = '\t')
  status(paste(family_name, 'MODELS AND RESULTS SAVED'))
}

# Consolidated views: all drugs, native calibration, actual predictions and fitting/CV.
z <- rbindlist(all_predictions, fill = TRUE)
z[, predicted_pct := 100 * 2^predicted_log2]
check(all(is.finite(z$predicted_pct)), 'Nonfinite native-unit prediction')
fwrite(z, file.path(root, 'data/evaluated_predictions.tsv'), sep = '\t')
assayed <- z[is.finite(measured_pct) & measured_pct > 0]
per_drug <- assayed[, {
  d <- copy(.SD)
  d[, drug := .BY$drug]
  summarize(d, d$predicted_log2)
}, by = .(family, model, evaluation, drug)]
fwrite(per_drug, file.path(root, 'data/all_drug_results.tsv'), sep = '\t')
status('all family results saved; rendering consolidated report')
colours <- c('46 SVR' = '#B24C35', '48 Program interactions' = '#27818A')
family_results <- per_drug[grepl('Excluded family', evaluation)]
family_results[, cohort := ifelse(grepl('familiar', evaluation), 'Familiar cells', 'Held-out cells')]
family_results[, label := paste(family, drug, sep = ' | ')]
family_results[, label := factor(label, levels = rev(unique(family_results[order(family, drug), label])))]
p_overview <- ggplot(family_results, aes(pooled_r, label, colour = model, shape = cohort)) +
  geom_vline(xintercept = 0, colour = 'grey65', linetype = 2) +
  geom_point(position = position_dodge(width = .45), size = 2.5) +
  scale_colour_manual(values = colours) + scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
  labs(title = 'Full response across all doses: every assayed excluded-family drug',
    subtitle = 'Same measured34 programs;46 SVR versus regularized48 program interactions.\nBoth exclude the evaluated family and nine held cells from fitting and selection.',
    x = 'Pearson r: measured versus predicted log2 viability', y = NULL, colour = NULL, shape = NULL,
    caption = 'Seven other compounds lack matched PRISM. Fixed bank discovery used all outcomes: retrospective comparison.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', axis.text.y = element_text(size = 9),
    plot.caption = element_text(hjust = 0))
ggsave(file.path(root, 'figures/all_drugs_response_ordering.pdf'), p_overview, width = 9, height = 5.4)
ggsave(file.path(root, 'figures/all_drugs_response_ordering.png'), p_overview, width = 9, height = 5.4, dpi = 160)

calibration <- copy(family_results)
calibration[, row_label := sprintf('%s | %s | %s  (S/R %d/%d)', family, drug, cohort, sensitive, resistant)]
calibration[, row_label := factor(row_label,
  levels = rev(unique(calibration[order(family, drug, cohort), row_label])))]
observed <- unique(calibration[, .(row_label, measured_sensitive_median_pct)])
p_calibration <- ggplot(calibration, aes(predicted_sensitive_median_pct, row_label, colour = model)) +
  geom_vline(xintercept = 50, colour = 'grey65', linetype = 2) +
  geom_point(position = position_dodge(width = .4), size = 2.3) +
  geom_point(data = observed, aes(measured_sensitive_median_pct, row_label),
    inherit.aes = FALSE, shape = 4, size = 2.5, colour = 'black') +
  scale_colour_manual(values = colours) +
  labs(title = 'Does the model recover the measured killing?',
    subtitle = 'Median viability among observations actually sensitive (<=50% vehicle), across all doses.\nBlack cross: measured median. Coloured dots: predicted median for the same observations.',
    x = 'PRISM % vehicle', y = NULL, colour = NULL,
    caption = 'S/R denotes sensitive/resistant dose-cell observations. Empty rows have no sensitive observations.\nFalse sensitivity calls are reported alongside sensitivity recovery in data/all_drug_results.tsv.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', axis.text.y = element_text(size = 8),
    plot.caption = element_text(hjust = 0, size = 8))
ggsave(file.path(root, 'figures/all_drugs_killing_calibration.pdf'), p_calibration, width = 11, height = 7.5)
ggsave(file.path(root, 'figures/all_drugs_killing_calibration.png'), p_calibration, width = 11, height = 7.5, dpi = 160)

dose_colours <- c('0.05' = '#009E73', '0.5' = '#E69F00', '5' = '#0072B2')
family_pages <- list()
for (family_name in c('ERBB', 'MEK')) {
  d <- copy(assayed[family == family_name & model == '48 Program interactions' & grepl('Excluded family', evaluation)])
  drugs <- sort(unique(d$drug))
  stats <- family_results[family == family_name & model == '48 Program interactions']
  stats[, key := paste(drug, cohort)]
  stats[, panel := sprintf('%s | %s\nr=%.2f | S recovered %d/%d; false S %d/%d',
    drug, cohort, pooled_r, sensitive_recovered, sensitive, false_sensitive, resistant)]
  panel_order <- as.vector(rbind(paste(drugs, 'Familiar cells'), paste(drugs, 'Held-out cells')))
  d[, cohort := ifelse(grepl('familiar', evaluation), 'Familiar cells', 'Held-out cells')]
  d[, panel := factor(stats$panel[match(paste(drug, cohort), stats$key)],
    levels = stats$panel[match(panel_order, stats$key)])]
  p <- ggplot(d, aes(measured_pct, predicted_pct)) +
    geom_vline(xintercept = 50, colour = 'grey80', linetype = 2) +
    geom_hline(yintercept = 50, colour = 'grey80', linetype = 2) +
    geom_point(aes(colour = factor(dose)), size = 1.5, alpha = .8) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = 'grey35', linewidth = .5) +
    scale_x_continuous(trans = 'log2', breaks = c(1, 5, 10, 25, 50, 100, 200)) +
    scale_y_continuous(trans = 'log2', breaks = c(1, 5, 10, 25, 50, 100, 200)) +
    scale_colour_manual(values = dose_colours) +
    facet_wrap(~panel, ncol = if (family_name == 'ERBB') 4 else 2, scales = 'free') +
    labs(title = paste(family_name, '|48 program interaction regression'),
      subtitle = 'All doses together; family excluded from fitting and whole-cell penalty selection.\nS: <=50% vehicle; false S: measured resistant but predicted sensitive.',
      x = 'Measured PRISM % vehicle (log2 scale)', y = 'Predicted % vehicle (log2 scale)', colour = 'Dose (uM)',
      caption = 'Independent panel axes. Fixed43 discovery used all outcomes. Measured displacement-to-function test only.') +
    theme_classic(base_size = 8) + theme(aspect.ratio = 1, legend.position = 'bottom', strip.text = element_text(size = 7),
      plot.caption = element_text(hjust = 0))
  family_pages[[family_name]] <- p
  ggsave(file.path(root, 'figures', paste0(family_name, '_all_drug_predictions.pdf')), p,
    width = if (family_name == 'ERBB') 12 else 8, height = if (family_name == 'ERBB') 12 else 10)
  ggsave(file.path(root, 'figures', paste0(family_name, '_all_drug_predictions.png')), p,
    width = if (family_name == 'ERBB') 12 else 8, height = if (family_name == 'ERBB') 12 else 10, dpi = 160)
}
fit_data <- assayed[evaluation %chin% c('Training fit', 'Whole-cell validation')]
fit_stats <- summary[evaluation %chin% c('Training fit', 'Whole-cell validation')]
fit_stats[, key := paste(family, model, evaluation)]
fit_stats[, panel := sprintf('%s | %s | %s\nr=%.2f; S recovered %d/%d; false S %d/%d',
  family, model, evaluation, pooled_r, sensitive_recovered, sensitive, false_sensitive, resistant)]
fit_data[, panel := fit_stats$panel[match(paste(family, model, evaluation), fit_stats$key)]]
p_fit <- ggplot(fit_data, aes(log2(measured_pct / 100), predicted_log2)) +
  geom_vline(xintercept = -1, colour = 'grey80', linetype = 2) +
  geom_hline(yintercept = -1, colour = 'grey80', linetype = 2) +
  geom_point(size = .5, alpha = .1, colour = '#27818A') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#B24C35', linewidth = .5) +
  facet_wrap(~panel, ncol = 4, scales = 'free') +
  labs(title = 'Where does the model lose the response?',
    subtitle = 'Non-family drugs only. Training fit is in-sample; validation excludes whole familiar cells.',
    x = 'Measured log2(% vehicle /100)', y = 'Predicted log2(% vehicle /100)') +
  theme_classic(base_size = 8) + theme(aspect.ratio = 1, strip.text = element_text(size = 7))
ggsave(file.path(root, 'figures/training_and_whole_cell_validation.pdf'), p_fit, width = 12, height = 7)
ggsave(file.path(root, 'figures/training_and_whole_cell_validation.png'), p_fit, width = 12, height = 7, dpi = 160)
grDevices::cairo_pdf(file.path(root, 'figures/combined_model_report.pdf'), width = 12, height = 12)
print(p_overview)
print(p_calibration)
for (p in family_pages) print(p)
print(p_fit)
grDevices::dev.off()
writeLines(c('48: same fixed43 measured programs and exact46 family/cell exclusions.',
  'Fractional binomial ridge:34 main terms,595 pair/quadratic terms; interaction penalty multiplier10.',
  'All dose rows retained; equal total weight per drug; train-only standardization.',
  '22 ridge penalties; minimum whole-cell validation error selects regularization.',
  'Target: continuous surviving fraction; measurements above100% vehicle treated as no killing.',
  'Validation: absolute error in capped%vehicle, averaged by drug and sensitive/resistant class.',
  'Predictions constrained to0-100% vehicle. Original assay values retained in reports and plots.',
  'No evaluated family or held-cell labels used to select the penalty.',
  '46 SVR reference unchanged. Every assayed original family drug reported.',
  'Seven compounds have no matched PRISM; predictions saved but performance not inferred.',
  'No molecular map redrawn. Baseline-to-displacement not evaluated.',
  'Bank discovery used all outcomes; previous test outcomes seen. Retrospective comparison.'), file.path(root, 'README.txt'))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('complete')
