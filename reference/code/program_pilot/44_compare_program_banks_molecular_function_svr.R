#!/usr/bin/env Rscript
# Compare frozen molecular banks with the same SVR and whole-cell decoder split.
# Discovery used all cells: held cells evaluate decoder transfer retrospectively.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
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

# Fixed earlier SVR choice; no new tuning or baseline/drug/dose predictors.
root <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr'
bank42_path <- '/work/results/program_pilot/42_global_pooled_signature_programs_ward/models/program_bank.qs'
bank43_path <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
split_path <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout/settings.qs'
setting_path <- '/work/results/program_pilot/38_nonlinear_absolute_and_recipient_sensitivity_decoder/models/selected_decoder.qs'
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading the two frozen banks and existing decoder split')
bank42 <- qread(bank42_path)
bank43 <- qread(bank43_path)
panel <- qread(panel_path)
held_ids <- qread(split_path)$held_ids
setting <- qread(setting_path)$selected_setting
check(identical(as.numeric(setting$cost), 100) && identical(as.numeric(setting$gamma_multiplier), 4),
  'Earlier selected SVR setting changed')
check(length(unique(bank42$partition)) == 60L && length(unique(bank43$partition)) == 34L,
  'Expected the existing 60- and 34-program banks')
check(!anyDuplicated(bank43$metadata$observation_id) && !anyDuplicated(panel$sensitivity$observation_id),
  'Duplicate molecular/assay identities')
md <- copy(bank43$metadata)
ann <- panel$sensitivity[match(md$observation_id, observation_id)]
check(identical(ann$observation_id, md$observation_id), 'Molecular and assay rows do not align')
md[, measured_pct := ann$prism_abundance_pct]
md[, exclusion := fifelse(!observation_id %chin% bank42$metadata$observation_id, 'absent from bank42',
  fifelse(!is.finite(measured_pct), 'no matched PRISM value',
    fifelse(measured_pct <= 0, 'nonpositive PRISM value cannot be log transformed', NA_character_)))]
fwrite(md[!is.na(exclusion), .(observation_id, depmap_id, recipient, drug, dose, measured_pct, exclusion)],
  file.path(root, 'data/excluded_observations.tsv'), sep = '\t')
md <- md[is.na(exclusion)]
md[, `:=`(split = fifelse(depmap_id %chin% held_ids, 'Held cells', 'Training cells'),
  condition = paste(drug, sprintf('%.12g', dose), sep = '__'), observed_log2 = log2(measured_pct / 100))]
check(all(is.finite(md$observed_log2)) && all(held_ids %chin% md$depmap_id), 'Invalid outcome or split')
training <- which(md$split == 'Training cells')
check(uniqueN(md[training, depmap_id]) == 26L && uniqueN(md[split == 'Held cells', depmap_id]) == 9L,
  'Expected the established 26-training/9-held cell split')
fwrite(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, split, measured_pct)],
  file.path(root, 'data/observations_and_split.tsv'), sep = '\t')
inputs <- c(bank42_path, bank43_path, panel_path, split_path, setting_path, panel$input_paths)
hashes <- data.table(path = inputs,
  sha256 = vapply(inputs, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(hashes, file.path(root, 'data/input_hashes.tsv'), sep = '\t')

# Score bank42 with hard means too, so the bank comparison has identical scoring.
status('computing hard assigned-gene means for bank42 on the common observations')
genes42 <- bank42$genes
k42 <- max(bank42$partition)
weights42 <- sparseMatrix(i = seq_along(genes42), j = bank42$partition[genes42],
  x = 1 / tabulate(bank42$partition[genes42], k42)[bank42$partition[genes42]],
  dims = c(length(genes42), k42), dimnames = list(genes42, sprintf('P%02d', seq_len(k42))))
scores42 <- matrix(NA_real_, nrow(md), k42, dimnames = list(md$observation_id, colnames(weights42)))
for (id in sort(unique(md$depmap_id))) {
  molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
  ii <- which(md$depmap_id == id)
  xx <- molecular$response[md$observation_id[ii], genes42, drop = FALSE]
  check(all(is.finite(xx)), 'Nonfinite bank42 molecular input')
  scores42[ii, ] <- as.matrix(xx %*% weights42)
}
rm(molecular, xx)
scores43 <- bank43$scores[match(md$observation_id, bank43$metadata$observation_id), , drop = FALSE]
rownames(scores43) <- md$observation_id
check(all(is.finite(scores42)) && all(is.finite(scores43)) &&
  identical(rownames(scores42), rownames(scores43)), 'Program score inputs do not align')
designs <- list('42 pooled / 60 programs' = scores42, '43 averaged / 34 programs' = scores43)
qsave(list(metadata = md, scores = designs, held_ids = held_ids, scoring = 'hard assigned-gene means',
  selected_setting = setting, input_hashes = hashes), file.path(root, 'data/frozen_comparison_inputs.qs'), preset = 'high')
signature <- digest::digest(list(hashes, md$observation_id, held_ids, setting, scoring = 'hard means'), algo = 'sha256')
rm(bank42, bank43, panel)
invisible(gc())

# Each model is saved as soon as fitted. Held outcomes never enter SVR fitting.
predictions <- list()
for (name in names(designs)) {
  x <- designs[[name]]
  model_path <- file.path(root, 'models', paste0(substr(name, 1, 2), '_svr.qs'))
  if (file.exists(model_path)) {
    saved <- qread(model_path)
    check(identical(saved$signature, signature), 'Cached decoder inputs changed')
    fitted <- saved$fitted
    status(paste('reusing saved SVR:', name))
  } else {
    status(paste('fitting SVR:', name, '|', length(training), 'training signatures'))
    fitted <- fit_program_svr(x[training, , drop = FALSE], md$observed_log2[training],
      kernel = 'radial', cost = setting$cost, gamma = setting$gamma_multiplier / ncol(x))
    qsave(list(fitted = fitted, signature = signature, bank = name, training_ids = md$observation_id[training],
      held_cells = held_ids, feature_count = ncol(x), setting = setting), model_path, preset = 'high')
    status(paste('MODEL SAVED:', name))
  }
  pred <- predict_program_svr(fitted, x)
  z <- md[, .(observation_id, depmap_id, recipient, drug, dose, condition, split)]
  z[, `:=`(bank = name, predicted_log2 = pred, predicted_pct = 100 * 2^pred)]
  check(all(is.finite(z$predicted_pct)), 'Nonfinite SVR predictions')
  fwrite(z, file.path(root, 'data', paste0(substr(name, 1, 2), '_outcome_free_predictions.tsv')), sep = '\t')
  predictions[[name]] <- z
  status(paste('PREDICTIONS SAVED:', name))
}
prediction_table <- rbindlist(predictions)
fwrite(data.table(sha256 = digest::digest(prediction_table, algo = 'sha256')),
  file.path(root, 'data/prediction_hash.tsv'), sep = '\t')

# Drug/dose means learned in training are a transparent comparator, not features.
means <- md[training, .(predicted_log2 = mean(observed_log2)), by = condition]
mean_prediction <- means$predicted_log2[match(md$condition, means$condition)]
check(all(is.finite(mean_prediction)), 'A held drug/dose lacks training assay support')
mean_rows <- md[, .(observation_id, depmap_id, recipient, drug, dose, condition, split)]
mean_rows[, `:=`(bank = 'Training drug/dose mean', predicted_log2 = mean_prediction,
  predicted_pct = 100 * 2^mean_prediction)]
z <- rbind(prediction_table, mean_rows)
z[, `:=`(observed_log2 = md$observed_log2[match(observation_id, md$observation_id)],
  measured_pct = md$measured_pct[match(observation_id, md$observation_id)])]
# Centering is diagnostic only: training retains every original dose/outcome.
z[, `:=`(observed_context = observed_log2 - mean(observed_log2),
  predicted_context = predicted_log2 - mean(predicted_log2),
  measured_context_pct = measured_pct - mean(measured_pct),
  predicted_context_pct = predicted_pct - mean(predicted_pct)), by = .(bank, split, condition)]
z[, `:=`(observed_cell = observed_log2 - mean(observed_log2),
  predicted_cell = predicted_log2 - mean(predicted_log2)), by = .(bank, split, depmap_id)]
correlation <- function(a, b) {
  if (length(a) < 3L || sd(a) < 1e-12 || sd(b) < 1e-12) return(NA_real_)
  cor(a, b)
}
contexts <- z[, .(cells = uniqueN(depmap_id), observations = .N,
  sensitive = sum(measured_pct <= 50), resistant = sum(measured_pct > 50),
  pearson_log2 = correlation(observed_log2, predicted_log2),
  pearson_pct = correlation(measured_pct, predicted_pct)), by = .(bank, split, drug, dose, condition)]
# Clearly sensitive/resistant pairs within a drug/dose; equal weight per context.
pairs <- z[, {
  resistant <- which(measured_pct >= 60)
  sensitive <- which(measured_pct <= 40)
  if (length(resistant) && length(sensitive)) {
    gap <- outer(measured_pct[resistant], measured_pct[sensitive], '-')
    valid <- gap >= 25
    predicted_gap <- outer(predicted_log2[resistant], predicted_log2[sensitive], '-')
    count <- sum(valid)
    accuracy <- if (count) mean((predicted_gap[valid] > 0) + .5 * (predicted_gap[valid] == 0)) else NA_real_
  } else {count <- 0L; accuracy <- NA_real_}
  .(pairs = count, correct_order_fraction = accuracy)
}, by = .(bank, split, drug, dose, condition)]
fwrite(contexts, file.path(root, 'data/per_drug_dose_correlations.tsv'), sep = '\t')
fwrite(pairs, file.path(root, 'data/per_drug_dose_sensitive_resistant_ordering.tsv'), sep = '\t')
summary <- z[, .(observations = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug),
  sensitive = sum(measured_pct <= 50), resistant = sum(measured_pct > 50),
  pooled_pearson_log2 = correlation(observed_log2, predicted_log2),
  pooled_pearson_pct = correlation(measured_pct, predicted_pct),
  same_drug_dose_centered_pearson_log2 = correlation(observed_context, predicted_context),
  same_drug_dose_centered_pearson_pct = correlation(measured_context_pct, predicted_context_pct),
  across_drugs_doses_within_cell_pearson_log2 = correlation(observed_cell, predicted_cell)), by = .(bank, split)]
summary <- merge(summary, contexts[, .(median_per_drug_dose_pearson_log2 = median(pearson_log2, na.rm = TRUE),
  contexts_with_correlation = sum(is.finite(pearson_log2))), by = .(bank, split)], by = c('bank', 'split'))
summary <- merge(summary, pairs[, .(clearly_separated_contexts = sum(pairs > 0),
  sensitive_resistant_correct_order_fraction = mean(correct_order_fraction, na.rm = TRUE)),
  by = .(bank, split)], by = c('bank', 'split'))
fwrite(summary, file.path(root, 'data/comparison_summary.tsv'), sep = '\t')
fwrite(z, file.path(root, 'data/evaluated_predictions.tsv'), sep = '\t')
drug_summary <- z[, .(observations = .N, cells = uniqueN(depmap_id),
  pooled_pearson_log2 = correlation(observed_log2, predicted_log2),
  same_drug_dose_centered_pearson_log2 = correlation(observed_context, predicted_context)),
  by = .(bank, split, drug)]
fwrite(drug_summary, file.path(root, 'data/per_drug_comparison.tsv'), sep = '\t')
print(summary)
status('all prediction metrics saved; rendering comparison figures')

# A pooled panel and a separate same-treatment panel prevent misleading pooling.
for (view in c('pooled', 'same_drug_dose')) {
  plots <- list()
  for (name in names(designs)) for (split_name in c('Training cells', 'Held cells')) {
    d <- z[bank == name & split == split_name]
    d[, `:=`(actual = if (view == 'pooled') observed_log2 else observed_context,
      prediction = if (view == 'pooled') predicted_log2 else predicted_context)]
    r <- correlation(d$actual, d$prediction)
    plots[[paste(name, split_name)]] <- ggplot(d, aes(actual, prediction)) +
      geom_point(size = .65, alpha = .18, colour = '#27818A') +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .65, colour = '#B24C35') +
      labs(title = paste(name, split_name, sep = '\n'),
        subtitle = sprintf('Pearson r = %.2f | n = %s', r, nrow(d)),
        x = if (view == 'pooled') 'Measured log2(% vehicle / 100)' else 'Measured deviation within drug/dose (log2)',
        y = if (view == 'pooled') 'Predicted log2(% vehicle / 100)' else 'Predicted deviation within drug/dose (log2)') +
      theme_classic(base_size = 10) + theme(aspect.ratio = 1, plot.title = element_text(size = 10))
  }
  plot <- wrap_plots(plots, ncol = 2) + plot_annotation(
    title = if (view == 'pooled') 'Molecular displacement to PRISM viability' else 'Differences between cells receiving the same drug and dose',
    subtitle = 'Same hard-mean scoring and frozen SVR setting; 26 training / 9 held cells\nBank discovery used all cells: retrospective decoder-transfer comparison')
  ggsave(file.path(root, 'figures', paste0(view, '_prediction_comparison.pdf')), plot, width = 9, height = 9)
  ggsave(file.path(root, 'figures', paste0(view, '_prediction_comparison.png')), plot, width = 9, height = 9, dpi = 160)
}
# The previously focal seven drugs plus afatinib show variation behind the aggregate.
focus <- c('afatinib', local_drugs)
shown <- drug_summary[drug %chin% focus & bank %chin% names(designs)]
wide <- dcast(shown, drug + split ~ bank, value.var = 'same_drug_dose_centered_pearson_log2')
setnames(wide, names(designs), c('r42', 'r43'))
plot <- ggplot(shown, aes(x = same_drug_dose_centered_pearson_log2, y = drug, colour = bank)) +
  geom_segment(data = wide, aes(x = r42, xend = r43, y = drug, yend = drug),
    inherit.aes = FALSE, colour = 'grey65', linewidth = .6) +
  geom_vline(xintercept = 0, colour = 'grey80', linetype = 2) +
  geom_point(size = 2.3) + facet_wrap(~split, nrow = 1) +
  scale_colour_manual(values = setNames(c('#27818A', '#B24C35'), names(designs))) +
  labs(title = 'Previously focal drugs and afatinib',
    subtitle = 'Predicting differences between cells receiving the same drug and dose\nTraining is in-sample; held cells evaluate retrospective decoder transfer',
    x = 'Pearson r after centering within drug/dose', y = NULL, colour = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', plot.title = element_text(size = 12))
ggsave(file.path(root, 'figures/focal_drug_comparison.pdf'), plot, width = 8, height = 4.6)
ggsave(file.path(root, 'figures/focal_drug_comparison.png'), plot, width = 8, height = 4.6, dpi = 160)
writeLines(c('Frozen bank42 vs bank43, with hard assigned-gene mean displacement inputs for both.',
  'RBF epsilon-SVR: cost100, gamma4/feature_count, epsilon0.1; all scaling learned in training.',
  'Input: measured molecular displacement only; no baseline, drug identity or dose predictor.',
  'Target: log2(PRISM %vehicle/100). All positive matched observations; doses retained independently.',
  'Existing whole-cell split: 26 training and9 held; no drug-family exclusion.',
  'Training predictions are in-sample. Held cells are excluded from decoder fitting only.',
  'Both banks were discovered with all cell outcomes, so this is retrospective and not independent validation.',
  'Same-treatment centering is performed only after predictions, separately within each split.',
  'Sensitive/resistant counts use50%vehicle. Ordering uses<=40 vs>=60, minimum25percentage-point separation.',
  'Drug/dose-mean comparator uses training outcomes only. Baseline→displacement is not tested.'),
  file.path(root, 'README.txt'))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('complete')
