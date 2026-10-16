#!/usr/bin/env Rscript
# Diagnose fitting capacity and ixazomib-specific baseline/response information.
# Frozen program banks and held recipients are reused; no discovery is repeated.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
})
source('/work/src/program_pilot/nonlinear_response.R')
source('/work/src/program_pilot/response_models.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, grepl('^--family-transfer=(MEK|ERBB)$', args))
family <- sub('--family-transfer=', '', args, fixed = TRUE)
root <- paste0('/work/results/program_pilot/24_ixazomib_decoder_capacity_and_molecular_signal_', family, '_holdout')
parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
recent <- paste0('/work/results/program_pilot/22_compressed_baseline_mutation_conditioning_', family, '_holdout')
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', family, ' | ', stage)
  fwrite(data.table(stage = stage, time = as.character(Sys.time())), file.path(root, 'run_status.tsv'), sep = '\t')
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 140)
}
describe <- function(z, groups) {
  z[, .(n = .N, pearson_r = cor_safe(observed_pct, predicted_pct),
    observed_median_pct = median(observed_pct), predicted_median_pct = median(predicted_pct),
    sensitive = sum(observed_pct <= 50), sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
    resistant = sum(observed_pct > 50), false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)), by = groups]
}
status('loading original training data and frozen predictions')
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
prepared <- qread(file.path(parent, 'data/prepared_scores.qs'))
settings <- qread(file.path(parent, 'settings.qs'))
base <- qread(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'))
paths <- c(panel_path, file.path(parent, c('settings.qs', 'data/prepared_scores.qs', 'models/nonlinear_RBF_SVR_D.qs')),
  '/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R',
  '/work/src/program_pilot/24_ixazomib_decoder_capacity_and_molecular_signal.R')
manifest <- data.table(path = paths, sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
save_table(manifest, 'input_hashes')
signature <- digest::digest(manifest, algo = 'sha256')
labels <- as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
train <- labels[depmap_id %in% base$training_ids & !drug %in% settings$excluded_drugs]
held <- labels[depmap_id %in% settings$held_ids]
stopifnot(!length(intersect(train$depmap_id, held$depmap_id)), !anyDuplicated(train$observation_id),
  !anyDuplicated(held$observation_id))
md <- rbindlist(list(copy(train)[, cohort := 'Fitting cells'], copy(held)[, cohort := 'Held cells']))
md[, observed_pct := prism_abundance_pct]
d <- prepared$measured[match(md$observation_id, prepared$metadata$observation_id), , drop = FALSE]
b <- prepared$baseline[md$depmap_id, , drop = FALSE]
colnames(b) <- paste0('B_', colnames(b))
stopifnot(all(is.finite(d)), all(is.finite(b)))
tr <- which(md$cohort == 'Fitting cells')
fold_table <- fread(file.path(parent, 'data/nonlinear_training_folds.tsv'))
fold <- fold_table$fold[match(md$depmap_id, fold_table$depmap_id)]
stopifnot(!anyNA(fold[tr]))
ix <- md[drug == 'ixazomib']
save_table(ix, 'ixazomib_observations')
save_table(ix[, .(cells = .N, sensitive = sum(observed_pct <= 50), resistant = sum(observed_pct > 50),
  median_abundance_pct = median(observed_pct)), by = .(cohort, dose)], 'dose_counts')

# Exact previously frozen predictions, plus a treatment-only training reference.
old <- fread(file.path(parent, 'data/nonlinear_predictions.tsv.gz'))[drug == 'ixazomib',
  .(observation_id, depmap_id, dose, model, predicted_pct)]
context <- fread(file.path(parent, 'data/contextual_predictions.tsv.gz'))[drug == 'ixazomib',
  .(observation_id, depmap_id, dose, model, predicted_pct = 100 * 2^predicted)]
latest <- fread(file.path(recent, 'data/predictions.tsv.gz'))[drug == 'ixazomib' & model == 'pca_mutation_actual',
  .(observation_id, depmap_id, dose, model = 'compressed_plus_mutations', predicted_pct = 100 * 2^predicted)]
reference <- ix[cohort == 'Fitting cells', .(predicted_pct = median(observed_pct)), by = dose]
reference <- merge(ix[cohort == 'Held cells', .(observation_id, depmap_id, dose)], reference, by = 'dose')
reference[, model := 'training_treatment_median']
prior <- rbindlist(list(old, context, latest, reference), use.names = TRUE)
prior[, `:=`(observed_pct = md$observed_pct[match(observation_id, md$observation_id)],
  cell = md$recipient[match(observation_id, md$observation_id)])]
save_table(prior, 'earlier_model_predictions')
save_table(describe(prior, c('dose', 'model')), 'earlier_model_summary')

# Independent models at each dose ask whether cell differences are present in
# baseline or displacement. Regularization uses only the 26 fitting cells.
lambda_grid <- 10^seq(3, -8, length.out = 45L)
local_predictions <- local_selection <- associations <- list()
for (dose_value in c(.05, .5, 5)) {
  rows <- which(md$drug == 'ixazomib' & md$dose == dose_value)
  fitting <- which(md$cohort[rows] == 'Fitting cells')
  eval_rows <- which(md$cohort[rows] == 'Held cells')
  matrices <- list(Baseline = b[rows, , drop = FALSE], Displacement = d[rows, , drop = FALSE],
    Both = cbind(b[rows, , drop = FALSE], d[rows, , drop = FALSE]))
  y <- log2(md$observed_pct[rows] / 100)
  for (kind in names(matrices)) {
    status(paste('ixazomib-specific ridge', dose_value, kind))
    x <- matrices[[kind]]
    cv <- cv.glmnet(x[fitting, , drop = FALSE], y[fitting], alpha = 0, lambda = lambda_grid,
      foldid = fold[rows[fitting]], type.measure = 'mae', standardize = TRUE, keep = TRUE, thresh = 1e-8)
    qsave(list(signature = signature, fit = cv, training_ids = md$depmap_id[rows[fitting]],
      kind = kind, dose = dose_value), file.path(root, 'models', paste0('ixazomib_', dose_value, '_', kind, '.qs')))
    selected <- as.numeric(predict(cv, newx = x, s = 'lambda.min'))
    flexible <- as.numeric(predict(cv$glmnet.fit, newx = x, s = min(lambda_grid)))
    for (choice in c('CV selected', 'Weak-penalty capacity')) {
      pp <- if (choice == 'CV selected') selected else flexible
      z <- copy(md[rows, .(observation_id, depmap_id, recipient, dose, cohort, observed_pct)])
      z[, `:=`(model = kind, choice = choice, predicted_pct = 100 * 2^pp)]
      local_predictions[[paste(dose_value, kind, choice)]] <- z
    }
    local_selection[[paste(dose_value, kind)]] <- data.table(dose = dose_value, model = kind,
      selected_lambda = cv$lambda.min, fitting_cells = length(fitting), held_cells = length(eval_rows))
    if (kind != 'Both') {
      aa <- rbindlist(lapply(seq_len(ncol(x)), function(j) {
        test <- cor.test(x[fitting, j], y[fitting])
        data.table(dose = dose_value, input = kind, program = colnames(x)[j],
          fitting_r_log_abundance = unname(test$estimate), fitting_p = test$p.value,
          held_r_log_abundance = cor_safe(x[eval_rows, j], y[eval_rows]))
      }))
      aa[, fitting_q := p.adjust(fitting_p, 'BH')]
      associations[[paste(dose_value, kind)]] <- aa
    }
  }
}
local <- rbindlist(local_predictions)
save_table(local, 'ixazomib_specific_predictions')
save_table(rbindlist(local_selection), 'ixazomib_specific_selection')
save_table(describe(local, c('dose', 'model', 'choice', 'cohort')), 'ixazomib_specific_summary')
save_table(rbindlist(associations), 'program_associations')

# A broad dose-response plot makes the sensitive/resistant coverage explicit.
p <- ggplot(ix, aes(dose, observed_pct, group = recipient, color = cohort)) +
  geom_line(alpha = .5) + geom_point(size = 1) + scale_x_log10(breaks = c(.05, .5, 5)) +
  geom_hline(yintercept = 50, linetype = 2, color = 'grey50') + theme_classic(base_size = 10) +
  labs(x = 'Ixazomib dose (uM)', y = 'PRISM abundance (% vehicle)', color = NULL,
    title = 'Ixazomib response across the same 35 cell lines')
save_plot(p, 'observed_dose_response', 6, 4)
selected_local <- local[choice == 'CV selected']
pdf(file.path(root, 'figures/ixazomib_specific_models.pdf'), width = 9, height = 6)
for (dose_value in c(.05, .5, 5)) {
  values <- selected_local[dose == dose_value]
  a <- describe(values, c('model', 'cohort'))
  a[, label := sprintf('r=%.2f; sensitive %d/%d', pearson_r, sensitive_called, sensitive)]
  p <- ggplot(values, aes(observed_pct, predicted_pct)) + geom_point(aes(color = observed_pct <= 50), size = 1.7) +
    geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
    geom_hline(yintercept = 50, linetype = 2, color = 'grey70') + geom_vline(xintercept = 50, linetype = 2, color = 'grey70') +
    geom_text(data = a, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 2.5) +
    facet_wrap(vars(cohort, model), ncol = 3, scales = 'free') + theme_classic(base_size = 10) +
    theme(aspect.ratio = 1, legend.position = 'none') + labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)',
      title = paste(family, '| ixazomib', dose_value, 'uM | drug/dose-specific ridge'),
      subtitle = 'Penalty selected using fitting-cell folds; held outcomes are retrospective')
  print(p)
  if (dose_value == .05) save_plot(p, 'ixazomib_005_specific_models', 9, 6)
}
dev.off()
earlier_labels <- c(original_ridge = 'Earlier conditioned ridge',
  baseline_reference = 'Drug-specific baseline reference',
  compressed_plus_mutations = 'Current conditioned SVR')
earlier_values <- copy(prior[dose == .05 & model %in% names(earlier_labels)])
earlier_values[, model_label := factor(earlier_labels[model], levels = earlier_labels)]
a <- describe(earlier_values, 'model_label')
a[, label := sprintf('r=%.2f; sensitive %d/%d', pearson_r, sensitive_called, sensitive)]
p <- ggplot(earlier_values, aes(observed_pct, predicted_pct)) +
  geom_point(aes(color = observed_pct <= 50), size = 2) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_hline(yintercept = 50, linetype = 2, color = 'grey70') +
  geom_vline(xintercept = 50, linetype = 2, color = 'grey70') +
  geom_text(data = a, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.2, size = 2.5) +
  facet_wrap(~model_label, nrow = 1, scales = 'free') + theme_classic(base_size = 9) +
  theme(aspect.ratio = 1, legend.position = 'none') +
  labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)',
    title = paste(family, '| ixazomib 0.05 uM | earlier versus current models'),
    subtitle = 'Same nine held cells; three sensitive and six resistant; previously frozen predictions')
save_plot(p, 'earlier_models_005', 9, 3.7)
case_labels <- c(original_ridge = 'Earlier conditioned ridge',
  calibrated_RBF_SVR_D = 'Earlier calibrated SVR', baseline_reference = 'Drug-specific baseline',
  compressed_plus_mutations = 'Current conditioned SVR', training_treatment_median = 'Training drug/dose median')
case <- copy(prior[cell == 'CFPAC-1' & model %in% names(case_labels)])
case[, curve := unname(case_labels[model])]
observed <- unique(case[, .(dose, predicted_pct = observed_pct)])
observed[, curve := 'Observed CFPAC-1']
case <- rbindlist(list(case[, .(dose, predicted_pct, curve)], observed))
p <- ggplot(case, aes(dose, predicted_pct, color = curve, group = curve)) +
  geom_line(linewidth = .7) + geom_point(size = 2) +
  scale_x_log10(breaks = c(.05, .5, 5)) + geom_hline(yintercept = 50, linetype = 2, color = 'grey70') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom') +
  labs(x = 'Ixazomib dose (uM)', y = 'Abundance (% vehicle)', color = NULL,
    title = paste(family, '| CFPAC-1 dose response: measured and predicted'))
save_plot(p, 'cfpac1_earlier_models', 8, 4.6)
status('ixazomib-specific analyses saved; testing full shared SVR capacity')

# The old grid never exceeded gamma=1/p. These fixed probes deliberately test
# substantially more local response functions, without selecting on held data.
capacity <- data.table(model = c('Original SVR', 'Local SVR', 'Very local SVR'),
  cost = c(base$selected$cost, 100, 1000), gamma_multiplier = c(base$selected$gamma_multiplier, 16, 256))
global_predictions <- capacity_settings <- list()
for (j in seq_len(nrow(capacity))) {
  spec <- capacity[j]
  status(paste('full fitting capacity:', spec$model))
  path <- file.path(root, 'models', paste0('capacity_', j, '.qs'))
  if (j == 1L) {
    fit <- base$model
  } else if (file.exists(path)) {
    saved <- qread(path)
    stopifnot(identical(saved$signature, signature))
    fit <- saved$model
  } else {
    fit <- fit_program_svr(d[tr, , drop = FALSE], log2(md$observed_pct[tr] / 100),
      'radial', spec$cost, spec$gamma_multiplier / ncol(d))
    qsave(list(signature = signature, model = fit, spec = spec, training_ids = base$training_ids), path)
  }
  pp <- predict_program_svr(fit, d)
  stopifnot(all(is.finite(pp)))
  z <- copy(md[, .(observation_id, depmap_id, recipient, drug, dose, cohort, observed_pct)])
  z[, `:=`(model = spec$model, predicted_pct = 100 * 2^pp)]
  global_predictions[[spec$model]] <- z
  capacity_settings[[spec$model]] <- cbind(spec, data.table(zero_input_abundance_pct =
    100 * 2^predict_program_svr(fit, matrix(0, 1L, ncol(d)))))
  gp <- rbindlist(global_predictions)
  fwrite(gp, file.path(root, 'data/shared_capacity_predictions.tsv.gz'), sep = '\t')
  save_table(describe(gp[drug == 'ixazomib'], c('model', 'cohort', 'dose')), 'shared_capacity_ixazomib_summary')
  save_table(describe(gp, c('model', 'cohort')), 'shared_capacity_global_summary')
  save_table(rbindlist(capacity_settings), 'capacity_settings')
}
pdf(file.path(root, 'figures/shared_decoder_capacity.pdf'), width = 9, height = 6)
for (dose_value in c(.05, .5, 5)) {
  values <- gp[drug == 'ixazomib' & dose == dose_value]
  values[, model := factor(model, levels = capacity$model)]
  a <- describe(values, c('model', 'cohort'))
  a[, label := sprintf('r=%.2f; sensitive %d/%d', pearson_r, sensitive_called, sensitive)]
  p <- ggplot(values, aes(observed_pct, predicted_pct)) + geom_point(aes(color = observed_pct <= 50), size = 1.7) +
    geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
    geom_hline(yintercept = 50, linetype = 2, color = 'grey70') + geom_vline(xintercept = 50, linetype = 2, color = 'grey70') +
    geom_text(data = a, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 2.5) +
    facet_wrap(vars(cohort, model), ncol = 3, scales = 'free') + theme_classic(base_size = 10) +
    theme(aspect.ratio = 1, legend.position = 'none') + labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)',
      title = paste(family, '| ixazomib', dose_value, 'uM | shared displacement SVR capacity'),
      subtitle = 'Fits use all training drugs; fixed capacity probes; no post-hoc zero anchoring')
  print(p)
  if (dose_value == .5) save_plot(p, 'shared_decoder_capacity_05', 9, 6)
}
dev.off()
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('COMPLETE')
print(describe(local[choice == 'CV selected'], c('dose', 'model', 'cohort')))
print(describe(gp[drug == 'ixazomib'], c('model', 'cohort', 'dose')))
