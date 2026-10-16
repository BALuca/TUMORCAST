#!/usr/bin/env Rscript
# Focus the frozen molecular bridge on training-defined informative focal doses.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

# Scope is selected from training phenotype variation, never held prediction success.
focal_drugs <- c('lapatinib', 'crizotinib', 'alpelisib', 'trametinib', 'gefitinib', 'everolimus', 'bortezomib')
minimum_cells <- 18L
minimum_each_class <- 4L
minimum_iqr_pct <- 25
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_MEK_holdout'
root <- '/work/results/program_pilot/30_targeted_grant_signal_pilots'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time()), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading frozen bridge and defining informative settings')
settings <- qread(file.path(parent, 'settings.qs'))
prepared <- qread(file.path(parent, 'data/prepared_scores.qs'))
decoder <- qread(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'))
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
folds <- fread(file.path(parent, 'data/nonlinear_training_folds.tsv'))
inner_map <- setNames(folds$fold, folds$depmap_id)
lambda_grid <- settings$lambda_grid
programs <- colnames(prepared$measured)
training <- labels[depmap_id %in% decoder$training_ids & !drug %in% settings$excluded_drugs &
  is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
stopifnot(!anyDuplicated(training$observation_id), uniqueN(training$depmap_id) == 26L,
  !any(decoder$training_ids %in% settings$held_ids), ncol(prepared$measured) == 40L)

sources <- c(panel_path, file.path(parent, c('settings.qs', 'data/prepared_scores.qs',
  'models/nonlinear_RBF_SVR_D.qs', 'data/nonlinear_training_folds.tsv',
  'data/nonlinear_predictions.tsv.gz', 'data/contextual_predictions.tsv.gz')),
  '/work/src/program_pilot/30_targeted_grant_signal_pilots.R',
  '/work/src/program_pilot/response_models.R', '/work/src/program_pilot/nonlinear_response.R',
  '/work/ERC_StG_2027_B2_v20.docx')
hashes <- data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256'))
save_table(hashes, 'input_hashes')
coverage <- training[, .(training_cells = uniqueN(depmap_id),
  sensitive = sum(prism_abundance_pct <= 50), resistant = sum(prism_abundance_pct > 50),
  q25_pct = as.numeric(quantile(prism_abundance_pct, .25)),
  q75_pct = as.numeric(quantile(prism_abundance_pct, .75)), iqr_pct = IQR(prism_abundance_pct)), by = .(drug, dose)]
coverage[, eligible := training_cells >= minimum_cells & sensitive >= minimum_each_class &
  resistant >= minimum_each_class & iqr_pct >= minimum_iqr_pct]
coverage[, reason := fcase(training_cells < minimum_cells, 'too few matched training cells',
  sensitive < minimum_each_class, 'too few sensitive training cells',
  resistant < minimum_each_class, 'too few resistant training cells',
  iqr_pct < minimum_iqr_pct, 'small abundance variation', default = 'eligible')]
save_table(coverage, 'all_training_settings')
focal_coverage <- merge(CJ(drug = focal_drugs, dose = c(.05, .5, 5)), coverage, by = c('drug', 'dose'), all.x = TRUE)
focal_coverage[drug %in% settings$excluded_drugs, `:=`(eligible = FALSE, reason = 'excluded MEK family')]
stopifnot(!anyNA(focal_coverage$reason))
save_table(focal_coverage, 'focal_training_settings')
selected <- coverage[eligible & drug %in% focal_drugs][order(drug, dose)]
stopifnot(nrow(selected) == 2L, identical(selected$drug, c('bortezomib', 'everolimus')), all(selected$dose == .05))
save_table(selected, 'selected_settings')
qsave(list(selected = selected, held_ids = settings$held_ids, training_ids = decoder$training_ids,
  minimum_cells = minimum_cells, minimum_each_class = minimum_each_class, minimum_iqr_pct = minimum_iqr_pct,
  functional_model = 'unchanged shared displacement-only RBF SVR', retrospective = TRUE), file.path(root, 'settings.qs'))

# The small upstream fit changes only baseline -> D; the shared D -> function map is frozen.
prior <- fread(file.path(parent, 'data/nonlinear_predictions.tsv.gz'))[model == 'RBF_SVR_D']
controls <- fread(file.path(parent, 'data/contextual_predictions.tsv.gz'))[
  model %in% c('baseline_reference', 'context_average')]
predictions <- molecular <- checks <- fits <- list()
for (j in seq_len(nrow(selected))) {
  drug_name <- selected$drug[j]
  dose_value <- selected$dose[j]
  status(paste('fitting baseline -> displacement:', drug_name, dose_value, 'uM'))
  fitting <- training[drug == drug_name & dose == dose_value]
  query <- prior[drug == drug_name & dose == dose_value,
    .(observation_id, depmap_id, recipient = cell, drug, dose)]
  stopifnot(nrow(fitting) == 26L, nrow(query) == 9L, !anyDuplicated(fitting$depmap_id),
    !anyDuplicated(query$depmap_id), !any(query$depmap_id %in% fitting$depmap_id))
  measured_fit <- prepared$measured[match(fitting$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  measured_query <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  baseline_fit <- prepared$baseline[fitting$depmap_id, , drop = FALSE]
  baseline_query <- prepared$baseline[query$depmap_id, , drop = FALSE]
  upstream <- fit_ridge(baseline_fit, measured_fit, fitting$depmap_id)
  inferred <- predict_ridge(upstream, baseline_query)
  mean_response <- matrix(rep(colMeans(measured_fit), each = nrow(query)), nrow(query),
    dimnames = list(NULL, programs))
  qsave(list(upstream = upstream, query = query, inferred = inferred,
    source_hashes = hashes, training_ids = fitting$depmap_id),
    file.path(root, 'models', paste0(drug_name, '_baseline_to_displacement.qs')))
  fits[[drug_name]] <- upstream
  response_inputs <- list(Measured_D = measured_query, Baseline_predicted_D = inferred, Treatment_mean_D = mean_response)
  for (name in names(response_inputs)) {
    tab <- copy(query)
    tab[, `:=`(model = name, predicted = predict_program_svr(decoder$model, response_inputs[[name]]))]
    predictions[[paste(drug_name, name)]] <- tab
  }
  ct <- controls[drug == drug_name & dose == dose_value]
  ct[, model := ifelse(model == 'baseline_reference', 'Baseline_drug_dose', 'Baseline_with_mean_D')]
  predictions[[paste0(drug_name, '_controls')]] <- ct
  for (k in seq_along(programs)) {
    molecular[[paste(drug_name, k)]] <- data.table(drug = drug_name, dose = dose_value,
      depmap_id = query$depmap_id, recipient = query$recipient, program = programs[k],
      measured = measured_query[, k], predicted = inferred[, k],
      training_mean = mean(measured_fit[, k]), training_sd = sd(measured_fit[, k]))
  }
  replay <- prior[drug == drug_name & dose == dose_value]
  difference <- max(abs(predictions[[paste(drug_name, 'Measured_D')]]$predicted -
    replay$predicted[match(query$observation_id, replay$observation_id)]))
  saved <- qread(file.path(root, 'models', paste0(drug_name, '_baseline_to_displacement.qs')))
  upstream_difference <- max(abs(predict_ridge(saved$upstream, baseline_query) - inferred))
  stopifnot(difference < 1e-10, upstream_difference < 1e-10, all(is.finite(inferred)))
  checks[[j]] <- data.table(drug = drug_name, dose = dose_value, fitting_cells = nrow(fitting), held_cells = nrow(query),
    held_training_overlap = length(intersect(fitting$depmap_id, query$depmap_id)),
    selected_lambda = upstream$lambda, shared_decoder_replay_difference = difference,
    upstream_replay_difference = upstream_difference)
  save_table(rbindlist(predictions, use.names = TRUE), 'outcome_free_predictions')
  save_table(rbindlist(checks), 'verification')
}
prediction_path <- file.path(root, 'data/outcome_free_predictions.tsv')
save_table(data.table(path = prediction_path, sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')

# Reveal labels only after freezing all targeted predictions. Every model uses the same nine cells.
prediction <- rbindlist(predictions, use.names = TRUE)
stopifnot(!anyDuplicated(prediction[, .(observation_id, model)]))
prediction[, observed_pct := labels$prism_abundance_pct[match(observation_id, labels$observation_id)]]
prediction[, predicted_pct := 100 * 2^predicted]
stopifnot(all(is.finite(prediction$observed_pct)), all(is.finite(prediction$predicted_pct)))
save_table(prediction, 'functional_predictions')
metrics <- prediction[, .(held_cells = .N, pearson_pct = cor_safe(observed_pct, predicted_pct),
  pearson_log2 = cor_safe(log2(observed_pct / 100), predicted),
  sensitive = sum(observed_pct <= 50), resistant = sum(observed_pct > 50),
  sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50),
  observed_resistant_minus_sensitive_pct = mean(observed_pct[observed_pct > 50]) - mean(observed_pct[observed_pct <= 50]),
  predicted_resistant_minus_sensitive_pct = mean(predicted_pct[observed_pct > 50]) - mean(predicted_pct[observed_pct <= 50])),
  by = .(drug, dose, model)]
save_table(metrics, 'functional_metrics')
# Resample recipients with frozen fits: uncertainty reflects this small held cohort.
uncertainty <- prediction[, {
  bootstrap_r <- replicate(1000L, {
    ix <- sample.int(.N, .N, replace = TRUE)
    cor_safe(observed_pct[ix], predicted_pct[ix])
  })
  loo_r <- vapply(seq_len(.N), function(i) cor_safe(observed_pct[-i], predicted_pct[-i]), numeric(1))
  valid <- bootstrap_r[is.finite(bootstrap_r)]
  valid_loo <- loo_r[is.finite(loo_r)]
  list(bootstrap_valid = length(valid),
    bootstrap_low = if (length(valid)) as.numeric(quantile(valid, .025)) else NA_real_,
    bootstrap_high = if (length(valid)) as.numeric(quantile(valid, .975)) else NA_real_,
    leave_one_cell_out_min = if (length(valid_loo)) min(valid_loo) else NA_real_,
    leave_one_cell_out_max = if (length(valid_loo)) max(valid_loo) else NA_real_)
}, by = .(drug, dose, model)]
save_table(uncertainty, 'functional_uncertainty')
mol <- rbindlist(molecular)
save_table(mol, 'molecular_predictions')
mol_metrics <- mol[, .(held_cells = .N, pearson = cor_safe(measured, predicted),
  measured_sd = sd(measured), predicted_sd = sd(predicted)), by = .(drug, dose, program)]
save_table(mol_metrics, 'all_program_molecular_metrics')

# Display training phenotype coverage, including doses without a useful class contrast.
plot_rows <- training[drug %in% focal_drugs]
plot_rows[, dose_label := factor(dose, levels = c(.05, .5, 5), labels = c('0.05', '0.5', '5'))]
coverage_plot <- ggplot(plot_rows, aes(dose_label, prism_abundance_pct)) +
  geom_hline(yintercept = 50, linetype = 2, color = '#888888') +
  geom_jitter(width = .15, height = 0, size = 1, alpha = .65, color = '#2166AC') +
  facet_wrap(~drug, ncol = 3) + theme_classic(base_size = 10) + theme(aspect.ratio = 1) +
  labs(x = 'Dose (µM)', y = 'PRISM abundance (% vehicle)',
    title = 'Training cells: which focal doses have an informative contrast?',
    subtitle = 'Dashed line: 50% vehicle. MEK drugs excluded by the frozen training partition.')
ggsave(file.path(root, 'figures/focal_training_dose_coverage.pdf'), coverage_plot, width = 8, height = 6.2)
ggsave(file.path(root, 'figures/focal_training_dose_coverage.png'), coverage_plot, width = 8, height = 6.2, dpi = 150)

# Square, independent-scale panels show both correlation and absolute call calibration.
display_models <- c('Measured_D', 'Baseline_predicted_D', 'Baseline_drug_dose')
plot_rows <- prediction[model %in% display_models]
annotation <- metrics[model %in% display_models]
annotation[, label := sprintf('r = %.2f; sensitive: %d/%d\nfalse sensitive: %d/%d',
  pearson_pct, sensitive_called, sensitive, false_sensitive, resistant)]
plot_rows[, model := factor(model, levels = display_models)]
annotation[, model := factor(model, levels = display_models)]
function_plot <- ggplot(plot_rows, aes(observed_pct, predicted_pct)) +
  geom_vline(xintercept = 50, linetype = 2, color = '#CCCCCC') +
  geom_hline(yintercept = 50, linetype = 2, color = '#CCCCCC') +
  geom_point(aes(color = recipient), size = 2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#2166AC', linewidth = .5) +
  geom_text(data = annotation, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.1, size = 2.5) +
  facet_wrap(vars(drug, model), ncol = 3, scales = 'free') + theme_classic(base_size = 10) +
  theme(aspect.ratio = 1, legend.position = 'bottom', strip.text = element_text(size = 9)) +
  labs(x = 'Measured PRISM abundance (% vehicle)', y = 'Predicted abundance (% vehicle)',
    title = 'Targeted molecular bridge | 0.05 µM | nine held cells per setting',
    subtitle = 'Unchanged shared decoder; scope selected from training variation. Retrospective evaluation.', color = NULL)
ggsave(file.path(root, 'figures/targeted_chain.pdf'), function_plot, width = 9, height = 7.1)
ggsave(file.path(root, 'figures/targeted_chain.png'), function_plot, width = 9, height = 7.1, dpi = 150)

report <- c('TARGETED GRANT SIGNAL PILOTS', '',
  'Training-only eligibility: >=18 cells, >=4 sensitive and >=4 resistant, abundance IQR >=25 percentage points.',
  'Both eligible focal settings assessed; no selection by held-out prediction success.',
  'Existing 40 hard-mean programs and shared displacement-only functional SVR retained.',
  'Only baseline -> displacement is newly fitted; no functional refit or global architecture search.',
  'Correlations compare different cells at the same drug and dose. All held evaluation is retrospective.', '',
  'FUNCTIONAL RESULTS', capture.output(print(metrics)), '',
  'MOLECULAR LINK (all 40 programs, no favorable-program selection)',
  capture.output(print(mol_metrics[, .(median_program_r = median(pearson, na.rm = TRUE),
    positive_programs = sum(pearson > 0, na.rm = TRUE)), by = drug])), '',
  'CAVEATS',
  'Strong correlation with all samples predicted sensitive establishes ordering, not sensitive/resistant calibration.',
  'The existing inner folds condition on a bank discovered on all outer-training cells; only outer cells are fully excluded.',
  'The grant historical lapatinib 441-cell cohort differs from these matched molecular cells.',
  'Training variance eligibility does not by itself solve or explain the broad model failure.')
writeLines(report, file.path(root, 'REPORT.txt'))
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
status('COMPLETE')
print(metrics)
