#!/usr/bin/env Rscript
# Report the original frozen program-to-function models on all four split combinations.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})
source('/work/src/program_pilot/response_models.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

# Preserve the preceding figures; this version broadens reporting, without refitting.
root <- '/work/results/program_pilot/32_training_and_held_drug_transfer_figures'
families <- c('ERBB', 'MEK')
bootstrap_replicates <- 1000L
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
for (path in file.path(root, c('data', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time()), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(6); quit(save = 'no', status = 1L)})
status('replaying frozen models on fitted and excluded recipients')
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
cell_map <- as.data.table(panel$cell_map)
stopifnot(!anyDuplicated(labels$observation_id),
  all(c('observation_id', 'depmap_id', 'drug', 'dose', 'prism_abundance_pct') %in% names(labels)))
comparisons <- data.table(model = c('measured_D_conditioned', 'measured_D', 'response_magnitude',
  'measured_D_conditioned', 'measured_D_conditioned'),
  response_input = c('measured', 'measured', 'measured', 'training_mean', 'swapped_recipient'),
  comparison = c('Programs + baseline', 'Programs alone', 'Magnitude + baseline',
    'Mean response + baseline', 'Swapped response + baseline'))
method_colors <- setNames(c('#0072B2', '#009E73', '#D55E00', '#777777', '#CC79A7'), comparisons$comparison)
prediction_parts <- original_intervals <- audit_parts <- availability_parts <- list()
sources <- c(panel_path, '/work/src/program_pilot/32_training_and_held_drug_transfer_figures.R',
  '/work/src/program_pilot/response_models.R', '/work/src/R/lib_heatmap.R')

for (family_name in families) {
  parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family_name, '_holdout')
  files <- file.path(parent, c('settings.qs', 'models/functional_decoder.qs', 'models/program_bank.qs',
    'data/prepared_scores.qs', 'data/functional_predictions.tsv.gz', 'data/functional_by_drug_uncertainty.tsv'))
  stopifnot(all(file.exists(files)))
  sources <- c(sources, files)
  settings <- qread(files[1L])
  decoder <- qread(files[2L])
  bank <- qread(files[3L])
  prepared <- qread(files[4L])
  original <- fread(files[5L])
  intervals <- fread(files[6L])
  intervals[, family := family_name]
  original_intervals[[family_name]] <- intervals
  programs <- colnames(bank$weights)
  lambda_grid <- settings$lambda_grid
  stopifnot(length(programs) == 40L, identical(settings$signature, prepared$signature),
    identical(settings$signature, decoder$signature), identical(settings$signature, bank$signature),
    identical(bank$excluded_drugs, settings$excluded_drugs),
    !any(bank$training_ids %in% settings$held_ids), !any(decoder$training_ids %in% settings$held_ids),
    !any(decoder$training_drugs %in% settings$excluded_drugs),
    identical(prepared$metadata$observation_id, rownames(prepared$measured)))

  # The development cohort is the exact 26 recipients with functional fit data.
  # Four additional molecular-discovery cells have no paired functional endpoints.
  query <- copy(labels[depmap_id %in% c(decoder$training_ids, settings$held_ids) &
    drug %in% c(decoder$training_drugs, settings$excluded_drugs) &
    is.finite(prism_abundance_pct) & prism_abundance_pct > 0])
  query[, `:=`(family = family_name, cell = cell_map$depmap_name[match(depmap_id, cell_map$depmap_id)],
    cell_cohort = fifelse(depmap_id %in% settings$held_ids, 'Held-out cells', 'Training cells'),
    drug_cohort = fifelse(drug %in% settings$excluded_drugs, 'Held-out drugs', 'Training drugs'),
    observed_pct = prism_abundance_pct, observed = log2(prism_abundance_pct / 100))]
  setorder(query, drug, dose, depmap_id)
  ix <- match(query$observation_id, prepared$metadata$observation_id)
  stopifnot(!anyNA(ix), !anyNA(query$cell))
  d <- prepared$measured[ix, , drop = FALSE]
  b <- sweep(sweep(prepared$baseline[query$depmap_id, , drop = FALSE], 2L,
    decoder$baseline_mean, '-'), 2L, decoder$baseline_sd, '/')
  designs <- list(measured_D = d, measured_D_conditioned = cbind(d, d * b),
    response_magnitude = cbind(magnitude = sqrt(rowMeans(d^2)), b * sqrt(rowMeans(d^2))))
  fitted <- rbindlist(lapply(names(designs), function(name) {
    z <- query[, .(observation_id, family, depmap_id, cell, cell_cohort, drug_cohort, drug, dose,
      observed, observed_pct, prism_status, prism_lower_dose, prism_upper_dose)]
    z[, `:=`(model = name, response_input = 'measured', predicted = predict_ridge(decoder$models[[name]], designs[[name]]))]
    z
  }))
  matched <- merge(fitted, original[response_input == 'measured', .(observation_id, model, saved = predicted)],
    by = c('observation_id', 'model'))
  replay_error <- max(abs(matched$predicted - matched$saved))
  stopifnot(nrow(matched) == nrow(original[response_input == 'measured']), replay_error < 1e-10)

  # Retain original saved mean and swap predictions for excluded recipients.
  controls <- original[model == 'measured_D_conditioned' & response_input != 'measured']
  controls[, c('cell', 'observed_pct', 'predicted_pct') := NULL]
  controls <- merge(controls, query[, .(observation_id, family, cell, cell_cohort, drug_cohort,
    observed_pct, prism_status, prism_lower_dose, prism_upper_dose)], by = 'observation_id')

  # Apply analogous descriptive controls to development recipients on excluded drugs.
  # The mean reference is the original molecular development cohort (30 cells),
  # including the query development cell; this is not a leave-one-cell-out control.
  tr <- which(query$drug_cohort == 'Held-out drugs' & query$cell_cohort == 'Training cells')
  reference <- t(vapply(tr, function(i) {
    ii <- prepared$metadata$drug == query$drug[i] & prepared$metadata$dose == query$dose[i] &
      prepared$metadata$depmap_id %in% bank$training_ids
    stopifnot(any(ii))
    colMeans(prepared$measured[ii, , drop = FALSE])
  }, numeric(40L)))
  permutation <- seq_along(tr)
  for (key in unique(paste(query$drug[tr], query$dose[tr]))) {
    ii <- which(paste(query$drug[tr], query$dose[tr]) == key)
    permutation[ii] <- sample(ii)
  }
  for (input in c('training_mean', 'swapped_recipient')) {
    dd <- if (input == 'training_mean') reference else d[tr[permutation], , drop = FALSE]
    z <- copy(fitted[model == 'measured_D_conditioned' & drug_cohort == 'Held-out drugs' &
      cell_cohort == 'Training cells'])
    stopifnot(identical(z$observation_id, query$observation_id[tr]))
    z[, `:=`(response_input = input, predicted = predict_ridge(decoder$models$measured_D_conditioned,
      cbind(dd, dd * b[tr, , drop = FALSE])))]
    controls <- rbindlist(list(controls, z), use.names = TRUE)
  }
  prediction_parts[[family_name]] <- rbindlist(list(fitted, controls), use.names = TRUE)
  audit_parts[[family_name]] <- data.table(family = family_name, discovery_cells = length(bank$training_ids),
    functional_training_cells = length(decoder$training_ids), held_cells = length(settings$held_ids),
    fitted_drugs = length(decoder$training_drugs), genes = length(bank$genes), programs = length(programs),
    training_rows = nrow(query[cell_cohort == 'Training cells' & drug_cohort == 'Training drugs']),
    maximum_original_prediction_difference = replay_error)
  availability <- data.table(family = family_name, drug = settings$excluded_drugs)
  availability <- merge(availability, query[drug_cohort == 'Held-out drugs',
    .(paired_rows = .N, paired_cells = uniqueN(depmap_id)), by = drug], by = 'drug', all.x = TRUE)
  availability[is.na(paired_rows), `:=`(paired_rows = 0L, paired_cells = 0L)]
  availability[, status := fifelse(paired_rows > 0L, 'Shown', 'No paired positive finite PRISM endpoint')]
  availability_parts[[family_name]] <- availability
}

# Checkpoint the unchanged frozen-model predictions before evaluating them.
predictions <- merge(rbindlist(prediction_parts), comparisons, by = c('model', 'response_input'))
predictions[, predicted_pct := 100 * 2^predicted]
stopifnot(all(is.finite(predictions$predicted)),
  !anyDuplicated(predictions[, .(family, observation_id, comparison)]))
save_table(predictions, 'frozen_predictions')
save_table(rbindlist(audit_parts), 'split_and_replay_audit')
save_table(rbindlist(availability_parts), 'excluded_compound_availability')
hashes <- data.table(path = unique(sources), sha256 = vapply(unique(sources), digest::digest,
  character(1), file = TRUE, algo = 'sha256'))
save_table(hashes, 'input_hashes')
status('calculating condition-specific associations and sensitivity counts')

# Center both endpoints separately in each drug/dose AND cell cohort.
# Pooling those residuals compares recipients, removing drug and dose mean effects.
predictions[, `:=`(observed_centered = observed - mean(observed),
  predicted_centered = predicted - mean(predicted)), by = .(family, comparison, cell_cohort, drug, dose)]
save_table(predictions, 'scatter_plot_values')
condition_metrics <- predictions[, .(cells = uniqueN(depmap_id),
  pearson_log2 = cor_safe(observed, predicted), pearson_pct = cor_safe(observed_pct, predicted_pct),
  sensitive = sum(observed_pct <= 50), resistant = sum(observed_pct > 50),
  predicted_sensitive = sum(predicted_pct <= 50),
  sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  resistant_called = sum(observed_pct > 50 & predicted_pct > 50),
  measured_min_pct = min(observed_pct), measured_max_pct = max(observed_pct),
  predicted_min_pct = min(predicted_pct), predicted_max_pct = max(predicted_pct),
  measured_sd_pct = sd(observed_pct), predicted_sd_pct = sd(predicted_pct)),
  by = .(family, cell_cohort, drug_cohort, drug, dose, model, response_input, comparison)]
save_table(condition_metrics, 'fixed_dose_metrics')
drug_metrics <- predictions[, .(cells = uniqueN(depmap_id), observations = .N,
  pearson = cor_safe(observed_centered, predicted_centered),
  measured_centered_sd = sd(observed_centered), predicted_centered_sd = sd(predicted_centered)),
  by = .(family, cell_cohort, drug_cohort, drug, model, response_input, comparison)]
drug_metrics[, spread_ratio := predicted_centered_sd / measured_centered_sd]
drug_metrics[, correlation_status := fifelse(is.finite(pearson), 'Evaluable',
  'Insufficient recipient contrast: fewer than three rows or zero centered variance')]
save_table(drug_metrics[!is.finite(pearson)], 'unevaluable_drug_correlations')

# Recipient bootstrap keeps every dose for a sampled recipient together. Original
# excluded-recipient intervals are reused exactly; new intervals are descriptive.
bootstrap_recipient_r <- function(z) {
  ids <- unique(z$depmap_id)
  ri <- match(z$depmap_id, ids)
  groups <- match(z$dose, unique(z$dose))
  values <- replicate(bootstrap_replicates, {
    multiplicity <- tabulate(sample.int(length(ids), length(ids), replace = TRUE), length(ids))
    w <- multiplicity[ri]
    mass <- as.numeric(rowsum(w, groups, reorder = FALSE))
    stopifnot(all(mass > 0))
    ox <- z$observed - as.numeric(rowsum(w * z$observed, groups, reorder = FALSE))[groups] / mass[groups]
    px <- z$predicted - as.numeric(rowsum(w * z$predicted, groups, reorder = FALSE))[groups] / mass[groups]
    denom <- sqrt(sum(w * ox^2) * sum(w * px^2))
    if (denom > 0) sum(w * ox * px) / denom else NA_real_
  })
  quantile(values[is.finite(values)], c(.025, .975), names = FALSE)
}
interval_parts <- list()
for (family_name in families) {
  for (drug_name in sort(unique(predictions[family == family_name & drug_cohort == 'Held-out drugs', drug]))) {
    z <- predictions[family == family_name & drug == drug_name & cell_cohort == 'Training cells']
    for (method in comparisons$comparison) {
      ci <- bootstrap_recipient_r(z[comparison == method])
      interval_parts[[paste(family_name, drug_name, method)]] <- data.table(family = family_name,
        drug = drug_name, cell_cohort = 'Training cells', comparison = method, ci_low = ci[1L], ci_high = ci[2L])
    }
  }
}
saved_intervals <- merge(rbindlist(original_intervals), comparisons, by = c('model', 'response_input'))
saved_intervals <- saved_intervals[, .(family, drug, cell_cohort = 'Held-out cells', comparison,
  ci_low = cell_bootstrap_lower, ci_high = cell_bootstrap_upper, expected_r = within_dose_pearson)]
replayed <- merge(drug_metrics, saved_intervals, by = c('family', 'drug', 'cell_cohort', 'comparison'))
stopifnot(all(abs(replayed$pearson - replayed$expected_r) < 1e-12))
ci_table <- rbindlist(list(rbindlist(interval_parts), saved_intervals[, !'expected_r']), use.names = TRUE)
drug_metrics <- merge(drug_metrics, ci_table, by = c('family', 'drug', 'cell_cohort', 'comparison'), all.x = TRUE)
save_table(drug_metrics, 'per_drug_metrics')

# Report both the residual-pool association and the distribution across drugs and
# fixed-dose settings. A pooled correlation is not the number of drugs that work.
aggregate <- predictions[response_input == 'measured', .(cells = uniqueN(depmap_id), drugs = uniqueN(drug),
  observations = .N, pooled_condition_centered_r = cor_safe(observed_centered, predicted_centered)),
  by = .(family, cell_cohort, drug_cohort, comparison)]
aggregate <- merge(aggregate, drug_metrics[response_input == 'measured',
  .(evaluable_drugs = sum(is.finite(pearson)), unevaluable_drugs = sum(!is.finite(pearson)),
    median_per_drug_r = median(pearson[is.finite(pearson)]), positive_drugs = sum(pearson > 0, na.rm = TRUE),
    drugs_above_point_three = sum(pearson > .3, na.rm = TRUE),
    median_spread_ratio = median(spread_ratio[is.finite(pearson)])),
  by = .(family, cell_cohort, drug_cohort, comparison)],
  by = c('family', 'cell_cohort', 'drug_cohort', 'comparison'))
aggregate <- merge(aggregate, condition_metrics[response_input == 'measured',
  .(dose_settings = .N, median_fixed_dose_r = median(pearson_log2, na.rm = TRUE),
    settings_with_both_outcomes = sum(sensitive > 0 & resistant > 0),
    sensitive = sum(sensitive), sensitive_called = sum(sensitive_called),
    resistant = sum(resistant), resistant_called = sum(resistant_called)),
  by = .(family, cell_cohort, drug_cohort, comparison)],
  by = c('family', 'cell_cohort', 'drug_cohort', 'comparison'))
save_table(aggregate, 'aggregate_metrics')

# The fitting target includes variation between interventions and doses. Contrast
# that mixed readout with within-recipient and fixed-treatment recipient readouts.
within_recipient <- predictions[response_input == 'measured',
  .(observations = .N, drugs = uniqueN(drug), pearson = cor_safe(observed, predicted)),
  by = .(family, cell_cohort, drug_cohort, comparison, cell)]
save_table(within_recipient, 'within_recipient_across_drug_dose_metrics')
axis_metrics <- predictions[response_input == 'measured',
  .(mixed_r = cor_safe(observed, predicted)), by = .(family, cell_cohort, drug_cohort, comparison)]
axis_metrics <- merge(axis_metrics, within_recipient[,
  .(median_within_recipient_r = median(pearson[is.finite(pearson)]),
    evaluable_recipients = sum(is.finite(pearson))), by = .(family, cell_cohort, drug_cohort, comparison)],
  by = c('family', 'cell_cohort', 'drug_cohort', 'comparison'))
axis_metrics <- merge(axis_metrics, aggregate[, .(family, cell_cohort, drug_cohort, comparison,
  pooled_condition_centered_r)], by = c('family', 'cell_cohort', 'drug_cohort', 'comparison'))
save_table(axis_metrics, 'comparison_axis_metrics')
status('exporting aggregate figures')

# Export each figure immediately, with its source table already saved above.
figures <- list()
save_plot <- function(p, stem, width, height, caption) {
  ggsave(file.path(root, 'figures', paste0(stem, '.pdf')), p, width = width, height = height,
    device = grDevices::cairo_pdf, bg = 'white')
  ggsave(file.path(root, 'figures', paste0(stem, '.png')), p, width = width, height = height, dpi = 150, bg = 'white')
  figures[[stem]] <<- list(object = p, kind = 'plot', caption = caption)
}
scatter_theme <- theme_classic(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom',
  strip.background = element_rect(fill = 'grey95', color = NA), plot.title = element_text(face = 'bold', size = 12))
scope_caption <- paste('Frozen original ridge models; measured hard-mean displacement inputs.',
  'Each family has its own 40-program bank and decoder. Shared recipients make the experiments dependent.',
  'Development discovery uses permitted functional labels; inspected evaluations are retrospective.',
  'Training cells on training drugs are in-sample. Excluded-family outcomes never fitted either decoder.')

full <- predictions[comparison == 'Programs + baseline']
full[, split := factor(paste(drug_cohort, cell_cohort, sep = '\n'), levels = c(
  'Training drugs\nTraining cells', 'Training drugs\nHeld-out cells',
  'Held-out drugs\nTraining cells', 'Held-out drugs\nHeld-out cells'))]
overview_labels <- aggregate[comparison == 'Programs + baseline']
overview_labels[, `:=`(split = factor(paste(drug_cohort, cell_cohort, sep = '\n'), levels = levels(full$split)),
  label = sprintf('Pooled r = %.2f\nMedian drug r = %.2f\n%d cells; %d evaluable drugs',
    pooled_condition_centered_r, median_per_drug_r, cells, evaluable_drugs))]
overview <- ggplot(full, aes(observed_centered, predicted_centered)) +
  geom_hline(yintercept = 0, color = 'grey90') + geom_vline(xintercept = 0, color = 'grey90') +
  geom_point(size = .75, alpha = .35, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#222222', linewidth = .5) +
  geom_text(data = overview_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.04, vjust = 1.15, size = 3) +
  facet_wrap(vars(family, split), ncol = 4, scales = 'free') + scatter_theme +
  labs(x = 'Measured log2 abundance: drug/dose mean removed',
    y = 'Predicted log2 abundance: drug/dose mean removed',
    title = 'Where does the shared program-to-function map transfer?',
    subtitle = 'Four split combinations; each mean is computed within its own cell cohort. Axes retain prediction compression.')
save_plot(overview, 'aggregate_01_four_split_combinations', 14, 8,
  paste('Points pool residuals after separate drug/dose centering. Training-drug panels include all fitted drugs;',
    'held-drug panels include all evaluable excluded family drugs. These are across-recipient comparisons.', scope_caption))

held_metrics <- drug_metrics[drug_cohort == 'Held-out drugs']
held_metrics[, cell_cohort := factor(cell_cohort, levels = c('Training cells', 'Held-out cells'))]
drug_levels <- rev(c(sort(unique(held_metrics[family == 'ERBB', drug])), sort(unique(held_metrics[family == 'MEK', drug]))))
held_metrics[, drug_label := factor(drug, levels = drug_levels)]
core <- held_metrics[response_input == 'measured']
forest <- ggplot(core, aes(pearson, drug_label, color = comparison)) +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = 'y', width = .15,
    position = position_dodge(width = .55), linewidth = .35) +
  geom_point(position = position_dodge(width = .55), size = 1.9) +
  facet_grid(family ~ cell_cohort, scales = 'free_y', space = 'free_y') +
  scale_color_manual(values = method_colors) + theme_classic(base_size = 10) +
  labs(x = 'Across-cell Pearson r after centering each dose', y = NULL, color = NULL,
    title = 'Every held-out drug: training recipients versus held-out recipients',
    subtitle = 'Intervals resample recipients with all their doses together; 25–26 training cells or nine held-out cells') +
  theme(legend.position = 'bottom', strip.background = element_rect(fill = 'grey95', color = NA))
save_plot(forest, 'aggregate_02_every_held_drug', 11, 7,
  paste('All seven paired ERBB drugs and all three paired MEK drugs, without success-based selection.',
    'A fitted recipient on an excluded drug is an intervention-transfer test, not an in-sample drug fit.', scope_caption))

# A compact fixed-dose matrix places correlation beside the available biological
# contrast (observed sensitive count / recipient count), without hiding doses.
fixed <- condition_metrics[drug_cohort == 'Held-out drugs' & comparison == 'Programs + baseline']
fixed[, row_id := paste(family, drug, sep = ' | ')]
fixed[, col_id := paste(cell_cohort, paste0(dose, ' µM'), sep = ' | ')]
rows <- c(paste('ERBB', sort(unique(fixed[family == 'ERBB', drug])), sep = ' | '),
  paste('MEK', sort(unique(fixed[family == 'MEK', drug])), sep = ' | '))
cols <- unlist(lapply(c('Training cells', 'Held-out cells'), function(cohort)
  paste(cohort, paste0(sort(unique(fixed$dose)), ' µM'), sep = ' | ')))
stopifnot(nrow(fixed) == length(rows) * length(cols), !anyDuplicated(fixed[, .(row_id, col_id)]))
map <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
map_text <- matrix('', length(rows), length(cols), dimnames = dimnames(map))
ij <- cbind(match(fixed$row_id, rows), match(fixed$col_id, cols))
map[ij] <- fixed$pearson_log2
map_text[ij] <- sprintf('r %.2f\n%d/%d sensitive', fixed$pearson_log2, fixed$sensitive, fixed$cells)
stopifnot(all(is.finite(map)), identical(dimnames(map), dimnames(map_text)))
column_annotations <- data.frame(Cohort = rep(c('Training cells', 'Held-out cells'), each = 3), row.names = cols)
stopifnot(identical(colnames(map), rownames(column_annotations)))
hm <- draw_heatmap(map, col_anno_df = column_annotations, col_anno_cols = 'Cohort',
  ann_top_args = list(Cohort = ComplexHeatmap::anno_block(
    labels = c('Training cells', 'Held-out cells'), gp = grid::gpar(fill = c('grey90', '#DCEAF2'), col = NA),
    labels_gp = grid::gpar(fontsize = 11), height = grid::unit(6, 'mm')), show_annotation_name = FALSE),
  scale_rows = FALSE, value_range = c(-1, 0, 1),
  value_colors = c('#B2182B', 'white', '#2166AC'), hm_name = 'Pearson r',
  show_row_names = TRUE, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = factor(sub(' \\| .*', '', rows), levels = families),
  column_split = factor(rep(c('Training cells', 'Held-out cells'), each = 3), levels = c('Training cells', 'Held-out cells')),
  column_labels = rep(paste0(sort(unique(fixed$dose)), ' µM'), 2), column_names_rot = 0,
  row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 10),
  title = 'Fixed-dose prediction and observed sensitivity contrast',
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(map_text[i, j], x, y, gp = grid::gpar(fontsize = 8,
      col = if (abs(map[i, j]) > .65) 'white' else '#222222'))
  })
grDevices::cairo_pdf(file.path(root, 'figures/aggregate_03_fixed_dose_limits.pdf'), width = 11, height = 6.5)
ComplexHeatmap::draw(hm)
dev.off()
png(file.path(root, 'figures/aggregate_03_fixed_dose_limits.png'), width = 11, height = 6.5, units = 'in', res = 150)
ComplexHeatmap::draw(hm)
dev.off()
figures[['aggregate_03_fixed_dose_limits']] <- list(object = hm, kind = 'heatmap',
  caption = paste('Each tile: fixed-dose across-cell log2 correlation plus measured samples at or below 50% vehicle.',
    'Correlation is not sensitivity-classification accuracy. Native-unit predictions and missed sensitive cases are shown per drug.', scope_caption))

known <- drug_metrics[drug_cohort == 'Training drugs' & response_input == 'measured' & is.finite(pearson)]
known[, comparison := factor(comparison, levels = comparisons$comparison[1:3])]
distribution <- ggplot(known, aes(comparison, pearson, color = cell_cohort)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_boxplot(position = position_dodge(width = .65), outlier.shape = NA, width = .55, linewidth = .45) +
  geom_point(position = position_jitterdodge(jitter.width = .1, dodge.width = .65), size = 1, alpha = .45) +
  facet_wrap(vars(family), ncol = 2) + theme_classic(base_size = 10) +
  scale_color_manual(values = c('Training cells' = '#777777', 'Held-out cells' = '#0072B2')) +
  labs(x = NULL, y = 'Across-cell Pearson r: dose-centered, one point per drug', color = NULL,
    title = 'Known drugs: does the relationship survive recipient holdout?',
    subtitle = 'All 90 ERBB-experiment or 94 MEK-experiment training drugs; grey points are in-sample fits') +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 15, hjust = 1))
save_plot(distribution, 'aggregate_04_known_drug_transfer', 11, 5.5,
  paste('One point per training drug per decoder, pooling only centered within-dose recipient differences.',
    'The box plots give equal weight to drugs. Vincristine has only one paired held-cell observation and is explicitly excluded',
    'from correlation summaries (unevaluable_drug_correlations.tsv). Held-cell outcomes did not fit these models.', scope_caption))

axis_plot_values <- melt(axis_metrics[comparison == 'Programs + baseline'],
  id.vars = c('family', 'cell_cohort', 'drug_cohort', 'comparison', 'evaluable_recipients'),
  measure.vars = c('mixed_r', 'median_within_recipient_r', 'pooled_condition_centered_r'),
  variable.name = 'axis', value.name = 'pearson')
axis_plot_values[, axis := factor(axis, levels = c('mixed_r', 'median_within_recipient_r', 'pooled_condition_centered_r'),
  labels = c('Mixed cells,\ndrugs and doses', 'Within a cell,\nacross drugs/doses', 'Across cells,\nfixed drug/dose'))]
save_table(axis_plot_values, 'comparison_axis_plot_values')
axis_plot <- ggplot(axis_plot_values, aes(axis, pearson, color = cell_cohort, group = cell_cohort)) +
  geom_hline(yintercept = 0, color = 'grey70') + geom_line(linewidth = .55) + geom_point(size = 2.3) +
  geom_text(aes(label = sprintf('%.2f', pearson)), vjust = -1, size = 3, show.legend = FALSE) +
  facet_grid(family ~ drug_cohort) + theme_classic(base_size = 10) +
  scale_color_manual(values = c('Training cells' = '#777777', 'Held-out cells' = '#0072B2')) +
  scale_y_continuous(limits = c(-.2, .85), breaks = c(0, .2, .4, .6, .8)) +
  labs(x = NULL, y = 'Pearson correlation', color = NULL,
    title = 'Which biological comparison supplies the apparent prediction signal?',
    subtitle = 'Middle: median of recipient-specific correlations. Right: pooled within-drug/dose residual correlation.') +
  theme(legend.position = 'bottom', strip.background = element_rect(fill = 'grey95', color = NA))
save_plot(axis_plot, 'aggregate_05_comparison_axes', 11, 7,
  paste('The left readout mixes intervention, dose and recipient differences. The middle compares drug/dose responses',
    'within each cell then takes the median cell correlation. The right removes each drug/dose mean and compares cells.',
    'These distinct questions must not be treated as equivalent validation.', scope_caption))

# Three clear pages per excluded drug: centered associations, native calibration,
# and all excluded-recipient dose curves. Do not select plots by performance.
detail_names <- character()
cell_colors <- setNames(c('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00',
  '#56B4E9', '#6A3D9A', '#666666', '#B15928'), sort(unique(full[cell_cohort == 'Held-out cells', cell])))
for (family_name in families) {
  for (drug_name in sort(unique(held_metrics[family == family_name, drug]))) {
    status(paste('exporting', family_name, drug_name))
    stem <- paste0(tolower(family_name), '_', gsub('[^a-z0-9]+', '_', drug_name))
    drug_rows <- copy(predictions[family == family_name & drug == drug_name])
    drug_rows[, cell_cohort := factor(cell_cohort, levels = c('Training cells', 'Held-out cells'))]
    dm <- copy(held_metrics[family == family_name & drug == drug_name])
    dm[, cell_cohort := factor(cell_cohort, levels = levels(drug_rows$cell_cohort))]
    dm[, label := sprintf('r = %.2f [%.2f, %.2f]\n%d recipients; %d observations',
      pearson, ci_low, ci_high, cells, observations)]
    measured <- drug_rows[response_input == 'measured']
    measured[, comparison := factor(comparison, levels = comparisons$comparison[1:3])]
    centered <- ggplot(measured, aes(observed_centered, predicted_centered)) +
      geom_point(aes(shape = factor(dose)), color = '#0072B2', size = 1.65, alpha = .75) +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#333333', linewidth = .45) +
      geom_text(data = dm[response_input == 'measured'], aes(x = -Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = -.05, vjust = 1.1, size = 2.8) +
      facet_wrap(vars(cell_cohort, comparison), ncol = 3, scales = 'free') + scatter_theme +
      labs(x = 'Measured log2 abundance: dose mean removed', y = 'Predicted log2 abundance: dose mean removed',
        shape = 'Dose (µM)', title = paste(family_name, '|', drug_name, '| transfer by recipient cohort'),
        subtitle = 'Drug outcomes excluded from model fitting in both rows; only the bottom row also excludes recipients.')
    centered_name <- paste0(stem, '_01_centered_comparisons')
    save_plot(centered, centered_name, 11, 8.5,
      paste('Three independently fitted original models evaluated on the same paired observations.',
        'Axes have independent ranges and show prediction compression. Intervals are conditional on frozen models.', scope_caption))

    native <- drug_rows[comparison == 'Programs + baseline']
    native[, `:=`(dose_label = factor(paste0(dose, ' µM'), levels = paste0(sort(unique(dose)), ' µM')),
      sensitivity = fifelse(observed_pct <= 50, 'Sensitive (≤50%)', 'Resistant (>50%)'),
      point_label = fifelse(cell_cohort == 'Held-out cells', cell, ''))]
    nm <- condition_metrics[family == family_name & drug == drug_name & comparison == 'Programs + baseline']
    nm[, `:=`(cell_cohort = factor(cell_cohort, levels = levels(drug_rows$cell_cohort)),
      dose_label = factor(paste0(dose, ' µM'), levels = levels(native$dose_label)),
      label = sprintf('Native r = %.2f; log2 r = %.2f\nSensitive: %d/%d measured, %d/%d predicted',
        pearson_pct, pearson_log2, sensitive, cells, predicted_sensitive, cells))]
    calibration <- ggplot(native, aes(observed_pct, predicted_pct)) +
      geom_vline(xintercept = 50, color = 'grey65', linetype = 2) +
      geom_hline(yintercept = 50, color = 'grey65', linetype = 2) +
      geom_point(aes(color = sensitivity), size = 1.8) +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#333333', linewidth = .45) +
      geom_text_repel(aes(label = point_label), size = 2.5, max.overlaps = Inf, seed = 20261010L,
        box.padding = .3, min.segment.length = 0, max.time = .3) +
      geom_text(data = nm, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -.04, vjust = 1.1, size = 2.7) +
      facet_wrap(vars(cell_cohort, dose_label), ncol = 3, scales = 'free') + scatter_theme +
      scale_color_manual(values = c('Sensitive (≤50%)' = '#D55E00', 'Resistant (>50%)' = '#0072B2')) +
      labs(x = 'Measured abundance (% vehicle)', y = 'Predicted abundance (% vehicle)', color = NULL,
        title = paste(family_name, '|', drug_name, '| fixed-dose native calibration'),
        subtitle = 'Dashed lines: 50% vehicle. Independent axes preserve native ranges; each facet compares cells at one dose.')
    calibration_name <- paste0(stem, '_02_fixed_dose_calibration')
    save_plot(calibration, calibration_name, 11, 8.5,
      paste('Full baseline-conditioned decoder at each dose. Sensitive = abundance ≤50% vehicle.',
        'No cross-dose or cross-drug pooling in these panels. All held-cell points are labeled.', scope_caption))

    curves <- melt(drug_rows[cell_cohort == 'Held-out cells'],
      id.vars = c('cell', 'dose', 'comparison'), measure.vars = c('observed_pct', 'predicted_pct'),
      variable.name = 'source', value.name = 'abundance_pct')
    observed_curve <- unique(curves[source == 'observed_pct', .(cell, dose, abundance_pct)])
    observed_curve[, curve := 'Measured']
    predicted_curve <- curves[source == 'predicted_pct', .(cell, dose, abundance_pct, curve = comparison)]
    curves <- rbindlist(list(observed_curve, predicted_curve))
    curve_colors <- c('Measured' = '#111111', method_colors)
    curve_plot <- ggplot(curves, aes(dose, abundance_pct, color = curve, group = curve)) +
      geom_hline(yintercept = 50, color = 'grey60', linetype = 2) +
      geom_line(linewidth = .6) + geom_point(size = 1.55) + facet_wrap(vars(cell), ncol = 3) +
      scale_x_log10(breaks = sort(unique(curves$dose))) + scale_color_manual(values = curve_colors) +
      theme_classic(base_size = 10) + theme(legend.position = 'bottom') +
      guides(color = guide_legend(ncol = 2)) + labs(x = 'Dose (µM)', y = 'Abundance (% vehicle)', color = NULL,
        title = paste(family_name, '|', drug_name, '| all nine held-out recipients'),
        subtitle = 'Same abundance scale across cells; inspect dose trend, sensitive calls and comparator separation.')
    curve_name <- paste0(stem, '_03_held_recipient_dose_curves')
    save_plot(curve_plot, curve_name, 10, 9,
      paste('Every held recipient and all five original comparisons; measured endpoint plotted once.',
        'Swapped response is one saved descriptive shuffle, not a permutation test.',
        'Training-mean response contains no held recipient molecular data.', scope_caption))
    selected <- c(centered_name, calibration_name, curve_name)
    detail_names <- c(detail_names, selected)
    grDevices::cairo_pdf(file.path(root, paste0(stem, '_drug_report.pdf')), width = 11, height = 9.5, onefile = TRUE)
    for (name in selected) print(figures[[name]]$object)
    dev.off()
  }
}

status('assembling aggregate and detailed collections')
aggregate_names <- paste0('aggregate_0', 1:5, c('_four_split_combinations', '_every_held_drug',
  '_fixed_dose_limits', '_known_drug_transfer', '_comparison_axes'))
for (collection in c('aggregate', 'all_held_drugs')) {
  selected <- if (collection == 'aggregate') aggregate_names else detail_names
  grDevices::cairo_pdf(file.path(root, paste0(collection, '_figure_collection.pdf')),
    width = if (collection == 'aggregate') 14 else 11, height = if (collection == 'aggregate') 8.5 else 9.5, onefile = TRUE)
  for (name in selected) {
    if (figures[[name]]$kind == 'plot') print(figures[[name]]$object) else ComplexHeatmap::draw(figures[[name]]$object)
  }
  dev.off()
}
index <- rbindlist(lapply(names(figures), function(name) data.table(figure = name,
  pdf = file.path(root, 'figures', paste0(name, '.pdf')),
  png = file.path(root, 'figures', paste0(name, '.png')), caption = figures[[name]]$caption)))
fwrite(index, file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
writeLines(unlist(lapply(names(figures), function(name) c(name, figures[[name]]$caption, ''))),
  file.path(root, 'figure_captions.txt'))
qsave(list(families = families, bootstrap_replicates = bootstrap_replicates, hashes = hashes,
  models_refitted = FALSE, retrospective = TRUE, aggregate_order = aggregate_names, detail_order = detail_names),
  file.path(root, 'settings.qs'))
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(system('git -C /work rev-parse HEAD', intern = TRUE), file.path(root, 'git_revision.txt'))
writeLines(c('TRAINING / EXCLUDED-DRUG TRANSFER FIGURES', '',
  'Start with aggregate_figure_collection.pdf (five pages).',
  'all_held_drugs_figure_collection.pdf: three pages for every paired excluded drug (30 pages).',
  'Individual three-page *_drug_report.pdf files sit directly in this numbered result root.',
  'FIGURE_INDEX.tsv links compact individual figures, PNG previews and full captions.', '',
  'FOUR TESTS',
  'Training cells + training drugs: apparent in-sample fit, not validation.',
  'Held-out cells + training drugs: recipient transfer.',
  'Training cells + held-out drugs: intervention transfer in familiar recipients.',
  'Held-out cells + held-out drugs: joint recipient/intervention transfer.', '',
  'READOUT',
  'Across-cell log2 correlation removes drug/dose means separately within each recipient cohort.',
  'Aggregate panels also report median per-drug correlation; a pooled correlation can hide failures.',
  'Vincristine has one paired held-cell observation; it cannot supply an across-recipient correlation.',
  'It is retained in prediction tables and explicitly excluded from correlation summaries.',
  'Fixed-dose native plots retain percentage vehicle and show observed/predicted sensitive counts (≤50%).',
  'A good ordering can coexist with compressed prediction amplitude and missed sensitive recipients.',
  'All endpoints are PRISM, including bracket-interpolated exposures documented in frozen_predictions.tsv.', '',
  'SCOPE',
  'Original family-excluded ridge models; no refitting and no later SVR substitutions.',
  'Each family has a different permitted bank. Their shared recipients are not independent replicates.',
  'Discovery used development outcomes. Fitted recipients are not held out from discovery.',
  'Mean-response controls for training recipients include the query recipient in the molecular reference.',
  'Single saved response swaps are descriptive. Previously inspected tests are retrospective.',
  'Baseline conditions measured responses. Baseline→response, combinations and HGF rescue are not tested here.', '',
  'FULL-MODEL AGGREGATE', capture.output(print(aggregate[comparison == 'Programs + baseline'])), '',
  'FULL-MODEL EXCLUDED DRUGS', capture.output(print(held_metrics[comparison == 'Programs + baseline',
    .(family, drug, cell_cohort, cells, pearson, ci_low, ci_high, spread_ratio)]))),
  file.path(root, 'READING_GUIDE.txt'))
stopifnot(nrow(index) == 35L, all(file.exists(index$pdf)), all(file.exists(index$png)),
  all(file.info(index$pdf)$size > 0), all(file.info(index$png)$size > 0))
status('COMPLETE')
message('Collections: ', root)
