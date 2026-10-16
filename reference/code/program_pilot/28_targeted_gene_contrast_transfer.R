#!/usr/bin/env Rscript
# Test whether gene contrasts missed by the local bank recur outside the index case.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)
root <- '/work/results/program_pilot/28_targeted_gene_contrast_transfer'
for (directory in file.path(root, c('data', 'models', 'figures'))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_MEK_holdout'
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
settings <- qread(file.path(parent, 'settings.qs'))
original <- qread(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'))
folds <- fread(file.path(parent, 'data/nonlinear_training_folds.tsv'))
eligible_ids <- c(original$training_ids, settings$held_ids)
stopifnot(length(original$training_ids) == 26L, length(settings$held_ids) == 9L,
  !anyDuplicated(eligible_ids))
target <- labels[depmap_id %in% eligible_ids & is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
  ((drug %in% c('anastrozole', 'ixazomib') & dose == .5) | (drug == 'raltitrexed' & dose == .05))]
target[, cohort := ifelse(depmap_id %in% original$training_ids, 'Fitting cells', 'Held cells')]
stopifnot(!anyDuplicated(target[, .(depmap_id, drug)]))
source_paths <- c(panel_path, file.path(parent, c('settings.qs', 'models/nonlinear_RBF_SVR_D.qs',
  'data/nonlinear_training_folds.tsv')), '/work/src/program_pilot/28_targeted_gene_contrast_transfer.R')
responses <- null_variances <- list()
message(format(Sys.time()), ' | reading only three target signatures per cell')
for (id in unique(target$depmap_id)) {
  path <- panel$input_paths[match(id, panel$cell_map$depmap_id)]
  source_paths <- c(source_paths, path)
  molecular <- qread(path)
  rows <- target[depmap_id == id, observation_id]
  stopifnot(all(rows %in% rownames(molecular$response)))
  responses[[id]] <- molecular$response[rows, , drop = FALSE]
  null_variances[[id]] <- matrixStats::colVars(molecular$null[[id]])
  if (id == 'ACH-000138') cfpac_null <- molecular$null[[id]]
}
response <- do.call(rbind, responses)
response <- response[target$observation_id, , drop = FALSE]
genes <- colnames(response)
stopifnot(all(is.finite(response)))
qsave(list(metadata = target, response = response, null_variances = null_variances),
  file.path(root, 'data/targeted_molecular_inputs.qs'))

# Learn a molecular drug-contrast direction using fitting cells only, without
# outcome-guided gene selection. Test its recurrence in the independent case.
pair <- merge(target[drug == 'ixazomib'], target[drug == 'anastrozole'],
  by = c('depmap_id', 'recipient', 'cohort'), suffixes = c('_ixazomib', '_anastrozole'))
coverage <- data.table(depmap_id = eligible_ids,
  cohort = ifelse(eligible_ids %in% original$training_ids, 'Fitting cells', 'Held cells'))
coverage[, paired_outcome_available := depmap_id %in% pair$depmap_id]
coverage[, reason := ifelse(paired_outcome_available, 'paired observations available',
  'no positive finite functional measurement for both target drugs')]
fwrite(coverage, file.path(root, 'data/paired_drug_coverage.tsv'), sep = '\t')
delta <- response[pair$observation_id_ixazomib, , drop = FALSE] -
  response[pair$observation_id_anastrozole, , drop = FALSE]
fit_rows <- which(pair$cohort == 'Fitting cells')
mean_direction <- colMeans(delta[fit_rows, , drop = FALSE])
training_sign_fraction <- colMeans(sweep(delta[fit_rows, , drop = FALSE], 2L, sign(mean_direction), '*') > 0)
index <- match('ACH-000138', pair$depmap_id)
stopifnot(!is.na(index), index %in% which(pair$cohort == 'Held cells'))
gene_table <- data.table(gene = genes, training_mean_difference = mean_direction,
  training_direction_fraction = training_sign_fraction, cfpac_difference = delta[index, ])
old_gene_path <- '/work/results/program_pilot/27_program_compression_failure_trace/data/gene_contrasts.tsv.gz'
source_paths <- c(source_paths, old_gene_path)
old <- fread(old_gene_path)[case == 'C01' & is.finite(difference_vehicle_sd)][order(-difference_vehicle_sd)]
gene_table[, cfpac_top50 := gene %in% head(old$gene, 50L)]
fwrite(gene_table, file.path(root, 'data/paired_drug_gene_recurrence.tsv'), sep = '\t')
top <- gene_table[cfpac_top50 == TRUE]
recurrence <- data.table(training_cells = length(fit_rows), held_cells = sum(pair$cohort == 'Held cells'),
  top50_same_training_mean_direction = sum(sign(top$training_mean_difference) == sign(top$cfpac_difference)),
  top50_median_training_direction_fraction = median(top$training_direction_fraction))
fwrite(recurrence, file.path(root, 'data/paired_drug_recurrence_summary.tsv'), sep = '\t')

# Decompose the shared gene direction into the part retained by hard means and
# the within-program part they discard. Compare the index projection to its
# own empirical vehicle contrasts; this tests molecular recurrence, not function.
projection_tables <- projection_observations <- list()
for (family in c('MEK', 'ERBB')) {
  bank_path <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_',
    family, '_holdout/models/program_bank.qs')
  bank <- qread(bank_path)
  source_paths <- c(source_paths, bank_path)
  direction <- mean_direction[bank$genes]
  retained <- ave(direction, bank$partition, FUN = mean)
  discarded <- direction - retained
  stopifnot(abs(sum(retained * discarded)) < 1e-10)
  for (part in c('Full genes in bank', 'Retained by hard means', 'Discarded within programs')) {
    axis <- switch(part, 'Full genes in bank' = direction,
      'Retained by hard means' = retained, 'Discarded within programs' = discarded)
    axis <- axis / sqrt(sum(axis^2))
    projection <- as.numeric(delta[, bank$genes, drop = FALSE] %*% axis)
    projection_observations[[paste(family, part)]] <- data.table(family = family, part = part,
      depmap_id = pair$depmap_id, recipient = pair$recipient, cohort = pair$cohort,
      molecular_projection = projection,
      observed_log2_ratio = log2(pair$prism_abundance_pct_ixazomib / pair$prism_abundance_pct_anastrozole))
    null_projection <- as.numeric(cfpac_null[, bank$genes, drop = FALSE] %*% axis)
    background <- as.vector(dist(matrix(null_projection)))
    projection_tables[[paste(family, part)]] <- data.table(family = family, part = part,
      cfpac_projection = projection[index], training_median_projection = median(projection[fit_rows]),
      held_median_projection = median(projection[pair$cohort == 'Held cells']),
      median_vehicle_projection_gap = median(background),
      fraction_vehicle_gaps_below_cfpac = mean(background < abs(projection[index])),
      training_direction_energy_retained_by_means = sum(retained^2) / sum(direction^2))
  }
}
fwrite(rbindlist(projection_tables), file.path(root, 'data/paired_drug_projection_decomposition.tsv'), sep = '\t')
projection_data <- rbindlist(projection_observations)
fwrite(projection_data, file.path(root, 'data/paired_drug_projection_observations.tsv'), sep = '\t')
fwrite(projection_data[, .(cells = .N, r_log2_ratio = cor(molecular_projection, observed_log2_ratio)),
  by = .(family, part, cohort)], file.path(root, 'data/paired_drug_projection_function_association.tsv'), sep = '\t')
plot <- ggplot(top, aes(training_mean_difference, cfpac_difference)) +
  geom_point(size = 1.5, color = '#2166AC') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey40') +
  theme_classic(base_size = 10) + theme(aspect.ratio = 1) +
  labs(x = sprintf('Mean gene contrast in %d fitting cells', length(fit_rows)), y = 'Gene contrast in held CFPAC-1',
    title = 'Do CFPAC-1 gene differences recur in fitting cells?',
    subtitle = sprintf('Ixazomib minus anastrozole, 0.5 uM | top 50 case contrasts\nPearson r = %.2f',
      cor(top$training_mean_difference, top$cfpac_difference)))
ggsave(file.path(root, 'figures/paired_drug_gene_recurrence.pdf'), plot, width = 5.7, height = 5)
ggsave(file.path(root, 'figures/paired_drug_gene_recurrence.png'), plot, width = 5.7, height = 5, dpi = 140)
message('Paired-drug recurrence and decomposition results saved')
print(recurrence)
print(rbindlist(projection_tables))

# Test predictive transfer of these representations, with the original recipient
# folds and training-only regularization. Paired contrast predicts the difference
# between outcomes; it is not a new baseline-only or combination architecture.
predictions <- summaries <- selections <- list()
for (task in c('Ixazomib minus anastrozole, 0.5 uM', 'Raltitrexed, 0.05 uM')) {
  if (task == 'Ixazomib minus anastrozole, 0.5 uM') {
    md <- pair[, .(depmap_id, recipient, cohort,
      observed_log2 = log2(prism_abundance_pct_ixazomib / prism_abundance_pct_anastrozole))]
    x <- delta
  } else {
    md <- target[drug == 'raltitrexed', .(depmap_id, recipient, cohort,
      observed_log2 = log2(prism_abundance_pct / 100))]
    x <- response[target[drug == 'raltitrexed', observation_id], , drop = FALSE]
  }
  fitting <- which(md$cohort == 'Fitting cells')
  fold <- folds$fold[match(md$depmap_id[fitting], folds$depmap_id)]
  stopifnot(length(fitting) >= 12L, !anyNA(fold), length(unique(fold)) == 3L,
    !any(md$depmap_id[fitting] %in% settings$held_ids))
  for (family in c('MEK', 'ERBB')) {
    parent_family <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
    bank <- qread(file.path(parent_family, 'models/program_bank.qs'))
    baseline_path <- file.path(parent_family, 'data/prepared_scores.qs')
    source_paths <- c(source_paths, baseline_path)
    prepared <- qread(baseline_path)
    inputs <- list('All measured genes' = x,
      'Bank genes before averaging' = x[, bank$genes, drop = FALSE],
      '40 hard program means' = as.matrix(x[, bank$genes, drop = FALSE] %*% bank$weights),
      '40 baseline programs' = prepared$baseline[md$depmap_id, , drop = FALSE])
    for (representation in names(inputs)) {
      xx <- inputs[[representation]]
      model <- cv.glmnet(xx[fitting, , drop = FALSE], md$observed_log2[fitting], alpha = 0,
        lambda = 10^seq(3, -8, length.out = 45), foldid = fold, type.measure = 'mae', standardize = TRUE)
      pp <- as.numeric(predict(model, newx = xx, s = 'lambda.min'))
      name <- gsub('[^A-Za-z0-9]+', '_', paste(task, family, representation))
      qsave(list(model = model, training_ids = md$depmap_id[fitting], features = colnames(xx),
        task = task, family = family, representation = representation), file.path(root, 'models', paste0(name, '.qs')))
      p <- copy(md)
      p[, `:=`(task = task, family = family, representation = representation,
        predicted_log2 = pp, observed_pct = 100 * 2^observed_log2, predicted_pct = 100 * 2^pp)]
      predictions[[name]] <- p
      frozen_prediction_path <- file.path(root, 'data/outcome_free_predictions.tsv')
      fwrite(rbindlist(predictions)[, .(depmap_id, recipient, cohort, task, family, representation,
        predicted_log2, predicted_pct)], frozen_prediction_path, sep = '\t')
      fwrite(data.table(path = frozen_prediction_path,
        sha256 = digest::digest(file = frozen_prediction_path, algo = 'sha256')),
        file.path(root, 'data/prediction_hash.tsv'), sep = '\t')
      selections[[name]] <- data.table(task = task, family = family, representation = representation,
        lambda = model$lambda.min)
      summaries[[name]] <- p[, .(cells = .N, r_log2 = cor(observed_log2, predicted_log2),
        r_pct = cor(observed_pct, predicted_pct),
        observed_below_half = sum(observed_pct <= 50), called_below_half = sum(observed_pct <= 50 & predicted_pct <= 50),
        observed_above_half = sum(observed_pct > 50), false_below_half = sum(observed_pct > 50 & predicted_pct <= 50)),
        by = .(task, family, representation, cohort)]
      fwrite(rbindlist(predictions), file.path(root, 'data/predictions.tsv'), sep = '\t')
      fwrite(rbindlist(summaries), file.path(root, 'data/representation_summary.tsv'), sep = '\t')
      fwrite(rbindlist(selections), file.path(root, 'data/regularization_selection.tsv'), sep = '\t')
    }
  }
  message(task, ': fitted models and predictions saved')
  print(rbindlist(summaries))
}
prediction_table <- rbindlist(predictions)
for (target_task in unique(prediction_table$task)) {
  show <- prediction_table[task == target_task & cohort == 'Held cells']
  plot <- ggplot(show, aes(observed_pct, predicted_pct)) + geom_point(size = 1.5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey40') +
    facet_grid(family ~ representation, scales = 'free') + theme_classic(base_size = 9) +
    theme(aspect.ratio = 1) + labs(x = if (grepl('minus', target_task))
      'Measured ixazomib / anastrozole abundance (%)' else 'Measured abundance (% vehicle)',
      y = 'Predicted (%)', title = target_task,
      subtitle = sprintf('%d held cells; regularization selected in %d fitting cells',
        uniqueN(show$depmap_id), uniqueN(prediction_table[task == target_task & cohort == 'Fitting cells', depmap_id])))
  association <- show[, .(r = cor(observed_pct, predicted_pct)), by = .(family, representation)]
  association[, label := sprintf('r = %.2f', r)]
  plot <- plot + geom_text(data = association, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.1, vjust = 1.2, size = 3)
  name <- gsub('[^A-Za-z0-9]+', '_', target_task)
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), plot, width = 11, height = 5.5)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), plot, width = 11, height = 5.5, dpi = 140)
}
source_paths <- unique(source_paths)
fwrite(data.table(path = source_paths, sha256 = vapply(source_paths, digest::digest,
  character(1), file = TRUE, algo = 'sha256')), file.path(root, 'data/input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
message('COMPLETE')
