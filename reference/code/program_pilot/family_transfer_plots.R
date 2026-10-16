# Shared reporting for saved family-transfer predictions; no model fitting.
# Direct invocation: Rscript src/program_pilot/family_transfer_plots.R
plot_family_transfer <- function(pred, figure_dir, data_dir) {
  rows <- data.table::copy(pred[model == 'measured_D_conditioned' & response_input == 'measured'])
  stopifnot(nrow(rows) > 0L, !anyDuplicated(rows[, .(drug, dose, cell)]),
    all(is.finite(rows$observed_pct)), all(is.finite(rows$predicted_pct)))
  data.table::fwrite(rows[, .(drug, dose, cell, observed_pct, predicted_pct)],
    file.path(data_dir, 'fixed_dose_plot_values.tsv'), sep = '\t')
  colours <- c('Measured PRISM' = '#0072B2', 'Predicted from programs' = '#D55E00')
  for (drug_name in unique(rows$drug)) {
    d <- data.table::copy(rows[drug == drug_name])
    # Keep cell order identical across doses, sorted by high-dose measured abundance.
    cell_order <- d[dose == max(dose)][order(observed_pct), as.character(cell)]
    stopifnot(setequal(cell_order, d$cell))
    d[, cell := factor(cell, levels = rev(cell_order))]
    d[, dose_label := factor(paste(dose, 'uM'), levels = paste(sort(unique(dose)), 'uM'))]
    points <- data.table::melt(d, id.vars = c('drug', 'dose', 'dose_label', 'cell'),
      measure.vars = c('observed_pct', 'predicted_pct'), variable.name = 'series', value.name = 'abundance')
    points[, series := factor(series, levels = c('observed_pct', 'predicted_pct'), names(colours))]
    p <- ggplot2::ggplot(d, ggplot2::aes(y = cell)) +
      ggplot2::geom_vline(xintercept = 100, colour = 'grey65', linetype = 'dashed', linewidth = .35) +
      ggplot2::geom_segment(ggplot2::aes(x = observed_pct, xend = predicted_pct, yend = cell),
        colour = 'grey65', linewidth = .55) +
      ggplot2::geom_point(data = points, ggplot2::aes(x = abundance, colour = series, shape = series), size = 2.4) +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::scale_shape_manual(values = c('Measured PRISM' = 16, 'Predicted from programs' = 17)) +
      ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(.02, .06))) +
      ggplot2::facet_wrap(~dose_label, nrow = 1) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(legend.position = 'bottom', strip.background = ggplot2::element_rect(fill = 'grey95'),
        panel.grid.major.y = ggplot2::element_line(colour = 'grey93', linewidth = .3),
        plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
      ggplot2::labs(x = 'Abundance (% vehicle): lower means more sensitive', y = NULL,
        colour = NULL, shape = NULL, title = paste0(drug_name, ': measured versus predicted sensitivity'),
        subtitle = 'Same drug and dose within each panel; each row is one held-out cell line',
        caption = paste('Dashed line: 100% vehicle. Connecting lines show prediction errors; all panels share the same scale.',
          'Prediction uses measured program displacement + baseline conditioning; drug family and cells excluded from training.',
          'PRISM values interpolated to molecular doses. Cell order follows measured abundance at the highest dose.', sep = '\n'))
    stem <- gsub('[^a-z0-9]+', '_', drug_name)
    ggplot2::ggsave(file.path(figure_dir, paste0(stem, '_fixed_dose_comparison.pdf')), p, width = 10, height = 4.7)
    ggplot2::ggsave(file.path(figure_dir, paste0(stem, '_fixed_dose_comparison.png')), p, width = 10, height = 4.7, dpi = 150)
  }
}

plot_family_diagnostics <- function(root, panel) {
  figure_dir <- file.path(root, 'figures')
  data_dir <- file.path(root, 'data')
  prepared <- qs::qread(file.path(data_dir, 'prepared_scores.qs'))
  decoder <- qs::qread(file.path(root, 'models/functional_decoder.qs'))
  settings <- qs::qread(file.path(root, 'settings.qs'))
  labels <- data.table::as.data.table(panel$sensitivity)[
    depmap_id %in% decoder$training_ids & !drug %in% settings$excluded_drugs &
      is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  displacement <- prepared$measured[match(labels$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  baseline <- sweep(sweep(prepared$baseline[labels$depmap_id, , drop = FALSE],
    2, decoder$baseline_mean), 2, decoder$baseline_sd, '/')
  labels[, observed := log2(prism_abundance_pct / 100)]
  labels[, predicted := as.numeric(predict(decoder$models$measured_D_conditioned$fit,
    newx = cbind(displacement, displacement * baseline), s = decoder$models$measured_D_conditioned$lambda))]
  centered <- data.table::copy(labels)
  centered[, `:=`(observed = observed - mean(observed), predicted = predicted - mean(predicted)), by = .(drug, dose)]
  raw_r <- cor(labels$observed, labels$predicted)
  centered_r <- cor(centered$observed, centered$predicted)
  panels <- data.table::rbindlist(list(
    labels[, .(observed, predicted, panel = sprintf('All training observations\nPearson r = %.3f', raw_r))],
    centered[, .(observed, predicted, panel = sprintf('Differences between cells at fixed drug/dose\nPearson r = %.3f', centered_r))]))
  data.table::fwrite(panels, file.path(data_dir, 'diagnostic_training_plot_values.tsv'), sep = '\t')
  p <- ggplot2::ggplot(panels, ggplot2::aes(observed, predicted)) +
    ggplot2::geom_point(size = .45, alpha = .12, colour = '#0072B2') +
    ggplot2::geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#D55E00', linewidth = .65) +
    ggplot2::facet_wrap(~panel, scales = 'free', nrow = 1) +
    ggplot2::theme_classic(base_size = 10) + ggplot2::theme(aspect.ratio = 1, plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(x = 'Measured log2 abundance relative to vehicle', y = 'Predicted log2 abundance relative to vehicle',
      title = paste(settings$family, '| What the fitted model captures in its training data'),
      caption = 'Right: drug-dose means removed from both axes. Each dot is one cell-drug-dose observation.\nThese are fitted training values, not held-out validation. Axes are scaled independently.')
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_training_signal.pdf'), p, width = 8, height = 4.8)
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_training_signal.png'), p, width = 8, height = 4.8, dpi = 150)

  audit <- data.table::fread(file.path(data_dir, 'diagnostic_hard_mean_audit.tsv'))[
    is.finite(admitted_outcome_r) & is.finite(full_outcome_r)]
  audit[, `:=`(admitted_strength = abs(admitted_outcome_r), full_strength = abs(full_outcome_r))]
  p <- ggplot2::ggplot(audit, ggplot2::aes(admitted_strength, full_strength)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = 'grey65', linetype = 'dashed') +
    ggplot2::geom_point(size = 1, alpha = .35, colour = '#0072B2') +
    ggplot2::geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#D55E00', linewidth = .65) +
    ggplot2::theme_classic(base_size = 10) + ggplot2::theme(aspect.ratio = 1, plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::labs(x = 'Association using genes selected in that cell (|r|)', y = 'Association using every gene in the program (|r|)',
      title = paste(settings$family, '| Effect of averaging the full global program'),
      subtitle = sprintf('Median |r|: %.2f selected genes to %.2f full program; Pearson r between axes = %.2f',
        median(audit$admitted_strength), median(audit$full_strength), cor(audit$admitted_strength, audit$full_strength)),
      caption = paste(nrow(audit), 'cell-program pairs; same global memberships.\nBelow dashed line: weaker association across drugs/doses within a training cell.\nSelection used these outcomes; this is a descriptive comparison.'))
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_gene_averaging.pdf'), p, width = 6.4, height = 6)
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_gene_averaging.png'), p, width = 6.4, height = 6, dpi = 150)

  comparison <- data.table::fread(file.path(data_dir, 'diagnostic_decoder_comparison.tsv'))[
    diagnostic %in% c('existing', 'centered_training', 'drug_specific_displacement')]
  comparison[, diagnostic := factor(diagnostic,
    levels = c('existing', 'centered_training', 'drug_specific_displacement'),
    labels = c('Original shared model', 'Fit cell differences directly', 'Separate fit for each drug'))]
  comparison[, spread_fraction := predicted_sd / observed_sd]
  values <- data.table::melt(comparison, id.vars = c('drug', 'diagnostic'),
    measure.vars = c('within_dose_r', 'spread_fraction'), variable.name = 'metric', value.name = 'value')
  values[, metric := factor(metric, levels = c('within_dose_r', 'spread_fraction'),
    labels = c('Across-cell Pearson correlation\nHigher is better; zero = no linear association',
      'Predicted / measured across-cell spread\n1 = matched spread; near 0 = almost constant'))]
  p <- ggplot2::ggplot(values, ggplot2::aes(value, drug, colour = diagnostic, shape = diagnostic)) +
    ggplot2::geom_vline(xintercept = 0, colour = 'grey70', linewidth = .35) +
    ggplot2::geom_vline(data = data.frame(metric = levels(values$metric)[2], reference = 1),
      ggplot2::aes(xintercept = reference), colour = 'grey60', linetype = 'dashed', linewidth = .35) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = .55), size = 2.3) +
    ggplot2::facet_wrap(~metric, scales = 'free_x', nrow = 1) +
    ggplot2::scale_colour_manual(values = c('#0072B2', '#D55E00', '#009E73')) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(legend.position = 'bottom', legend.box = 'vertical', plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), shape = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, shape = NULL,
      title = paste(settings$family, '| Do the diagnostic refits recover sensitivity differences?'),
      caption = 'Nine held-out cells; three doses, centered within drug/dose. Spread is SD in log2 abundance units.\nSeparate drug fits use same-drug outcomes in training cells; they test cell transfer, not family transfer.')
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_refits.pdf'), p, width = 8.5, height = 4.8)
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_refits.png'), p, width = 8.5, height = 4.8, dpi = 150)
}

plot_within_cell_models <- function(root) {
  rows <- data.table::fread(file.path(root, 'data/diagnostic_within_cell_model_comparison.tsv'))
  ordering <- rows[model == 'cell_specific_D'][order(pearson), cell]
  rows[, cell := factor(cell, levels = ordering)]
  rows[, model := factor(model, levels = c('shared_D', 'shared_conditioned', 'cell_specific_D'),
    labels = c('Shared equation: displacement', 'Shared equation: displacement + baseline interactions', 'Separate equation for each cell'))]
  p <- ggplot2::ggplot(rows, ggplot2::aes(pearson, cell, colour = model, shape = model)) +
    ggplot2::geom_vline(xintercept = 0, colour = 'grey65', linewidth = .35) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = .55), size = 2) +
    ggplot2::scale_colour_manual(values = c('#0072B2', '#D55E00', '#009E73')) +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(legend.position = 'bottom',
      legend.box = 'vertical', plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1), shape = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::labs(x = 'Pearson correlation across held-out drugs/doses within each cell', y = NULL, colour = NULL, shape = NULL,
      title = 'Can the same 40 displacement scores reproduce the signal within a cell?',
      subtitle = 'Identical whole-drug holdouts for shared and cell-specific equations; nested training-only penalty selection',
      caption = 'All three models use the same fixed program scores. Green allows a separate mapping for each training cell.\nProgram discovery used these development outcomes: this tests mapping flexibility conditional on that bank.\nIt does not test prediction in unseen cell lines or independently validate program discovery.')
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_within_cell_models.pdf'), p, width = 8.5, height = 9)
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_within_cell_models.png'), p, width = 8.5, height = 9, dpi = 120)
}

plot_held_cell_generalization <- function(root, panel) {
  data_dir <- file.path(root, 'data')
  settings <- qs::qread(file.path(root, 'settings.qs'))
  prepared <- qs::qread(file.path(data_dir, 'prepared_scores.qs'))
  decoder <- qs::qread(file.path(root, 'models/functional_decoder.qs'))
  rows <- data.table::as.data.table(panel$sensitivity)[depmap_id %in% settings$held_ids &
    is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  d <- prepared$measured[match(rows$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  b <- sweep(sweep(prepared$baseline[rows$depmap_id, , drop = FALSE],
    2, decoder$baseline_mean), 2, decoder$baseline_sd, '/')
  rows[, predicted := as.numeric(predict(decoder$models$measured_D_conditioned$fit,
    newx = cbind(d, d * b), s = decoder$models$measured_D_conditioned$lambda))]
  rows[, `:=`(observed = log2(prism_abundance_pct / 100), predicted_pct = 100 * 2^predicted,
    cell = panel$cell_map$depmap_name[match(depmap_id, panel$cell_map$depmap_id)],
    drug_set = ifelse(drug %in% settings$excluded_drugs, 'Excluded drug family', 'Drugs seen in training'))]
  rows[, outcome_class := factor(ifelse(prism_abundance_pct <= 50, 'Sensitive: <=50% vehicle', 'Resistant: >50% vehicle'),
    levels = c('Sensitive: <=50% vehicle', 'Resistant: >50% vehicle'))]
  summary <- rows[, .(drugs = data.table::uniqueN(drug), conditions = .N,
    sensitive = sum(prism_abundance_pct <= 50), resistant = sum(prism_abundance_pct > 50),
    sensitive_called_sensitive = sum(prism_abundance_pct <= 50 & predicted_pct <= 50),
    resistant_called_resistant = sum(prism_abundance_pct > 50 & predicted_pct > 50),
    pearson_log2 = cor(observed, predicted), pearson_pct = cor(prism_abundance_pct, predicted_pct),
    predicted_sensitive_median_pct = if(any(prism_abundance_pct <= 50)) median(predicted_pct[prism_abundance_pct <= 50]) else NA_real_,
    predicted_resistant_median_pct = if(any(prism_abundance_pct > 50)) median(predicted_pct[prism_abundance_pct > 50]) else NA_real_),
    by = .(cell, drug_set)]
  data.table::fwrite(summary, file.path(data_dir, 'diagnostic_held_cell_across_drugs.tsv'), sep = '\t')
  data.table::fwrite(rows[, .(observation_id, cell, drug, dose, drug_set, prism_abundance_pct, predicted_pct, outcome_class)],
    file.path(data_dir, 'diagnostic_held_cell_across_drug_predictions.tsv.gz'), sep = '\t')
  for (drug_set_name in unique(rows$drug_set)) {
    selected <- rows[drug_set == drug_set_name]
    annotation <- summary[drug_set == drug_set_name]
    annotation[, panel := sprintf('%s | Pearson r = %.2f\n%d sensitive, %d resistant conditions',
      cell, pearson_pct, sensitive, resistant)]
    selected[, panel := annotation$panel[match(cell, annotation$cell)]]
    p <- ggplot2::ggplot(selected, ggplot2::aes(prism_abundance_pct, predicted_pct)) +
      ggplot2::geom_vline(xintercept = 50, colour = 'grey65', linetype = 'dashed', linewidth = .3) +
      ggplot2::geom_hline(yintercept = 50, colour = 'grey65', linetype = 'dashed', linewidth = .3) +
      ggplot2::geom_point(ggplot2::aes(colour = outcome_class), size = 1.1, alpha = .6) +
      ggplot2::geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#D55E00', linewidth = .5) +
      ggplot2::facet_wrap(~panel, scales = 'free', ncol = 3) +
      ggplot2::scale_colour_manual(values = c('#0072B2', 'grey60'), drop = FALSE) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom', plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
      ggplot2::labs(x = 'Measured abundance (% vehicle)', y = 'Predicted abundance (% vehicle)', colour = NULL,
        title = paste(settings$family, '| Sensitive versus resistant treatments within held-out cells'),
        subtitle = paste(drug_set_name, '| each dot is one drug-dose condition; all nine cells excluded from training'),
        caption = 'Lower abundance means greater sensitivity. Dashed lines mark 50% vehicle; blue dots are observed sensitive conditions.\nPrediction uses measured 40-program displacement and baseline interactions. PRISM interpolated to molecular doses.\nPearson correlations use percent-vehicle values; source table also reports log2 correlations. Panel axes are independent.')
    stem <- if(drug_set_name == 'Excluded drug family') 'diagnostic_held_cells_new_family' else 'diagnostic_held_cells_known_drugs'
    ggplot2::ggsave(file.path(root, 'figures', paste0(stem, '.pdf')), p, width = 10.5, height = 10)
    ggplot2::ggsave(file.path(root, 'figures', paste0(stem, '.png')), p, width = 10.5, height = 10, dpi = 120)
  }
  print(summary[, .(median_within_cell_r_log2 = median(pearson_log2), sensitive = sum(sensitive),
    sensitive_called_sensitive = sum(sensitive_called_sensitive), resistant = sum(resistant),
    resistant_called_resistant = sum(resistant_called_resistant)), by = drug_set])
}

plot_nonlinear_response <- function(root) {
  summary <- data.table::fread(file.path(root, 'data/nonlinear_summary.tsv'))
  names <- c('original_ridge', 'calibrated_ridge', 'ridge_DB', 'linear_SVR_DB',
    'RBF_SVR_D', 'calibrated_RBF_SVR_D', 'RBF_SVR_DB', 'calibrated_RBF_SVR_DB')
  descriptions <- c('Original ridge', 'Calibrated original ridge', 'Ridge: displacement + baseline',
    'Linear SVR: displacement + baseline', 'RBF SVR: displacement only',
    'Calibrated RBF SVR: displacement only', 'RBF SVR: displacement + baseline',
    'Calibrated RBF SVR: displacement + baseline')
  summary[, model := factor(model, levels = rev(names), labels = rev(descriptions))]
  summary[, `:=`(sensitive_pct = 100 * sensitive_correct / sensitive,
    resistant_pct = 100 * resistant_correct / resistant)]
  long <- data.table::melt(summary, id.vars = c('model', 'drug_set', 'sensitive',
    'sensitive_correct', 'resistant', 'resistant_correct'),
    measure.vars = c('median_within_cell_r_log2', 'sensitive_pct', 'resistant_pct'),
    variable.name = 'metric', value.name = 'value')
  long[, label := ifelse(metric == 'median_within_cell_r_log2', sprintf('%.2f', value),
    ifelse(metric == 'sensitive_pct', paste0(sensitive_correct, '/', sensitive),
      paste0(resistant_correct, '/', resistant)))]
  long[, metric := factor(metric, levels = c('median_within_cell_r_log2', 'sensitive_pct', 'resistant_pct'),
    labels = c('Median within-cell correlation\n(log2 abundance)', 'Sensitive conditions\ncorrectly called (%)',
      'Resistant conditions\ncorrectly called (%)'))]
  p <- ggplot2::ggplot(long, ggplot2::aes(value, model)) +
    ggplot2::geom_point(colour = '#0072B2', size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -.2, size = 2.4) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.08, .4))) +
    ggplot2::facet_grid(drug_set ~ metric, scales = 'free_x') +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::labs(x = NULL, y = NULL, title = 'Do nonlinear models recover sensitivity in the nine excluded cells?',
      subtitle = 'Identical frozen programs and recipient/family exclusions; settings and calibration use training cells only',
      caption = 'Sensitive means abundance <=50% vehicle; resistant means >50%. Within-cell correlations compare drug-dose conditions.\nHigher sensitive-call rates must be assessed alongside resistant-call rates. Previously inspected outcomes: retrospective evaluation.')
  ggplot2::ggsave(file.path(root, 'figures/nonlinear_model_comparison.pdf'), p, width = 11, height = 6.5)
  ggplot2::ggsave(file.path(root, 'figures/nonlinear_model_comparison.png'), p, width = 11, height = 6.5, dpi = 130)
  pred <- data.table::fread(file.path(root, 'data/nonlinear_predictions.tsv.gz'))
  # Compare drugs at a fixed dose within each recipient, alongside assay-scale
  # separation of measured sensitive and resistant conditions.
  by_dose <- pred[, .(conditions = .N, pearson_log2 = cor(observed, predicted)),
    by = .(model, drug_set, cell, dose)]
  data.table::fwrite(by_dose, file.path(root, 'data/nonlinear_by_cell_dose.tsv'), sep = '\t')
  magnitude <- pred[, .(conditions = .N, observed_median_pct = median(observed_pct),
    predicted_median_pct = median(predicted_pct), predicted_lower_decile_pct = quantile(predicted_pct, .1),
    predicted_upper_decile_pct = quantile(predicted_pct, .9)),
    by = .(model, drug_set, sensitive = observed_pct <= 50)]
  data.table::fwrite(magnitude, file.path(root, 'data/nonlinear_assay_magnitude.tsv'), sep = '\t')
  pred <- pred[model %in% c('original_ridge', 'calibrated_RBF_SVR_D')]
  pred[, model := factor(model, levels = c('original_ridge', 'calibrated_RBF_SVR_D'),
    labels = c('Original ridge', 'Calibrated RBF SVR: displacement'))]
  # Pool the nine held recipients in one figure, retaining cell identity and
  # separate panels for familiar drugs versus the excluded family.
  within_cell <- pred[, .(r = cor(observed_pct, predicted_pct), r_log2 = cor(observed, predicted)),
    by = .(drug_set, model, cell)]
  combined <- pred[, .(pooled_r = cor(observed_pct, predicted_pct), conditions = .N,
    sensitive = sum(observed_pct <= 50),
    sensitive_correct = sum(observed_pct <= 50 & predicted_pct <= 50)), by = .(drug_set, model)]
  combined <- merge(combined, within_cell[, .(median_within_cell_r = median(r),
    median_within_cell_r_log2 = median(r_log2)), by = .(drug_set, model)],
    by = c('drug_set', 'model'))
  combined[, label := sprintf(paste0('Median within-cell r: log2 = %.2f; %% vehicle = %.2f\n',
    'Pooled r (%% vehicle) = %.2f\nSensitive called: %d/%d; conditions: %d'),
    median_within_cell_r_log2, median_within_cell_r, pooled_r, sensitive_correct, sensitive, conditions)]
  data.table::fwrite(combined, file.path(root, 'data/nonlinear_combined_correlations.tsv'), sep = '\t')
  family <- if (grepl('_MEK_holdout$', root)) 'MEK' else 'ERBB'
  combined_plot <- ggplot2::ggplot(pred, ggplot2::aes(observed_pct, predicted_pct, colour = cell)) +
    ggplot2::geom_vline(xintercept = 50, colour = 'grey65', linetype = 'dashed', linewidth = .3) +
    ggplot2::geom_hline(yintercept = 50, colour = 'grey65', linetype = 'dashed', linewidth = .3) +
    ggplot2::geom_point(size = .8, alpha = .45) +
    ggplot2::geom_smooth(ggplot2::aes(group = 1), colour = 'black', method = 'lm',
      formula = y ~ x, se = FALSE, linewidth = .6) +
    ggplot2::geom_text(data = combined, ggplot2::aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -.03, vjust = 1.1, size = 2.6) +
    ggplot2::facet_wrap(~drug_set + model, ncol = 2, scales = 'free') +
    ggplot2::scale_colour_brewer(palette = 'Set1') +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.05, .23))) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom',
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, override.aes = list(alpha = 1, size = 2))) +
    ggplot2::labs(x = 'Measured abundance (% vehicle)', y = 'Predicted abundance (% vehicle)', colour = 'Held-out cell',
      title = paste0(family, ' exclusion | All nine held-out cells'),
      subtitle = 'Original ridge versus calibrated displacement-only RBF SVR',
      caption = 'Panels separate drugs seen in training from the excluded family. Dashed lines mark 50% vehicle.\nAxes use percent vehicle; within-cell correlations show both log2 and percent-vehicle scales. Calibration uses training cells only.\nMeasured molecular displacements; PRISM interpolated to molecular doses. Retrospective evaluation; independent axes.')
  ggplot2::ggsave(file.path(root, 'figures/nonlinear_combined.pdf'), combined_plot, width = 8.5, height = 9)
  ggplot2::ggsave(file.path(root, 'figures/nonlinear_combined.png'), combined_plot, width = 8.5, height = 9, dpi = 140)
  for (set_name in unique(pred$drug_set)) {
    rows <- pred[drug_set == set_name]
    correlations <- rows[, .(r = cor(observed_pct, predicted_pct)), by = .(cell, model)]
    correlations[, label := sprintf('%s: r = %.2f', model, r)]
    annotation <- correlations[, .(label = paste(label, collapse = '\n')), by = cell]
    p <- ggplot2::ggplot(rows, ggplot2::aes(observed_pct, predicted_pct, colour = model)) +
      ggplot2::geom_vline(xintercept = 50, colour = 'grey65', linetype = 'dashed', linewidth = .3) +
      ggplot2::geom_hline(yintercept = 50, colour = 'grey65', linetype = 'dashed', linewidth = .3) +
      ggplot2::geom_point(size = .8, alpha = .4) +
      ggplot2::geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .55) +
      ggplot2::geom_text(data = annotation, ggplot2::aes(x = -Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = -.03, vjust = 1.1, size = 2.3) +
      ggplot2::facet_wrap(~cell, ncol = 3, scales = 'free') +
      ggplot2::scale_colour_manual(values = c('#0072B2', '#D55E00')) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom', plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
      ggplot2::labs(x = 'Measured abundance (% vehicle)', y = 'Predicted abundance (% vehicle)', colour = NULL,
        title = paste('Original versus nonlinear prediction |', set_name),
        subtitle = 'Same drug-dose conditions in each excluded cell; lower abundance means greater sensitivity',
        caption = 'Dashed lines mark 50% vehicle. RBF SVR uses measured displacement; calibration uses training-cell predictions.\nPRISM values interpolated to molecular doses. Axes are independent. Source table provides all model predictions.')
    stem <- if(set_name == 'Seen drugs') 'nonlinear_known_drugs' else 'nonlinear_excluded_family'
    ggplot2::ggsave(file.path(root, 'figures', paste0(stem, '.pdf')), p, width = 10, height = 9.5)
    ggplot2::ggsave(file.path(root, 'figures', paste0(stem, '.png')), p, width = 10, height = 9.5, dpi = 120)
  }
}

# Report frozen contextual predictions; outcomes are attached only here.
report_contextual_response <- function(root, panel) {
  pred <- data.table::fread(file.path(root, 'data/contextual_predictions.tsv.gz'))
  settings <- qs::qread(file.path(root, 'settings.qs'))
  labels <- data.table::as.data.table(panel$sensitivity)
  ii <- match(pred$observation_id, labels$observation_id)
  stopifnot(!anyNA(ii), !anyDuplicated(pred[, .(observation_id, model)]), all(is.finite(pred$predicted)))
  pred[, `:=`(observed_pct = labels$prism_abundance_pct[ii], predicted_pct = 100 * 2^predicted,
    cell = recipient, drug_set = ifelse(drug %in% settings$excluded_drugs, 'Excluded family', 'Seen drugs'))]
  pred[, observed := log2(observed_pct / 100)]
  stopifnot(all(is.finite(pred$predicted_pct)))
  save <- function(x, name) data.table::fwrite(x,
    file.path(root, 'data', paste0('contextual_', name, '.tsv')), sep = '\t')
  save(pred, 'evaluation')
  by_cell <- pred[, .(conditions = .N, pearson_log2 = cor_safe(observed, predicted),
    pearson_pct = cor_safe(observed_pct, predicted_pct), sensitive = sum(observed_pct <= 50),
    sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50), resistant = sum(observed_pct > 50),
    false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)), by = .(model, drug_set, cell)]
  save(by_cell, 'by_cell')
  by_cell_dose <- pred[, .(conditions = .N, pearson_log2 = cor_safe(observed, predicted)),
    by = .(model, drug_set, cell, dose)]
  save(by_cell_dose, 'by_cell_dose')
  pred[, `:=`(observed_centered = observed - mean(observed), predicted_centered = predicted - mean(predicted),
    observed_pp_centered = observed_pct - mean(observed_pct),
    predicted_pp_centered = predicted_pct - mean(predicted_pct)), by = .(model, drug, dose)]
  within <- pred[, .(within_treatment_r_log2 = cor_safe(observed_centered, predicted_centered),
    within_treatment_r_pct = cor_safe(observed_pp_centered, predicted_pp_centered),
    observed_spread_pp = sd(observed_pp_centered), predicted_spread_pp = sd(predicted_pp_centered)),
    by = .(model, drug_set)]
  save(within, 'within_treatment')
  summary <- by_cell[, .(median_within_cell_r_log2 = median(pearson_log2, na.rm = TRUE),
    sensitive = sum(sensitive), sensitive_called = sum(sensitive_called),
    resistant = sum(resistant), false_sensitive = sum(false_sensitive)), by = .(model, drug_set)]
  summary <- merge(summary, by_cell_dose[, .(median_fixed_dose_within_cell_r = median(pearson_log2, na.rm = TRUE)),
    by = .(model, drug_set)], by = c('model', 'drug_set'))
  summary <- merge(summary, within, by = c('model', 'drug_set'))
  save(summary, 'summary')
  # Show why shrinkage was selected using development predictions only. These
  # curves never evaluate alternative penalties on the outer held-out cells.
  context_settings <- qs::qread(file.path(root, 'contextual_settings.qs'))
  training <- labels[depmap_id %in% context_settings$training_ids &
    !drug %in% settings$excluded_drugs & is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  cv_parts <- list()
  for (variant in c('actual', 'average')) {
    for (k in 1:3) {
      saved <- qs::qread(file.path(root, 'models', sprintf('contextual_predictions_%s_fold%d.qs', variant, k)))
      md <- training[saved$validation_rows]
      cv_parts[[paste(variant, k)]] <- data.table::data.table(variant = variant,
        lambda = rep(context_settings$penalties, each = nrow(md)),
        observed = rep(log2(md$prism_abundance_pct / 100), times = length(context_settings$penalties)),
        predicted = as.vector(saved$prediction), cell = rep(md$recipient, times = length(context_settings$penalties)),
        drug = rep(md$drug, times = length(context_settings$penalties)),
        dose = rep(md$dose, times = length(context_settings$penalties)))
    }
  }
  cv <- data.table::rbindlist(cv_parts)
  cv[, `:=`(yy = observed - mean(observed), pp = predicted - mean(predicted)), by = .(variant, lambda, drug, dose)]
  cv_treatment <- cv[, .(across_cells_r = cor_safe(yy, pp)), by = .(variant, lambda)]
  cv_cell <- cv[, .(r = cor_safe(observed, predicted)), by = .(variant, lambda, cell, dose)]
  cv_metrics <- merge(cv_treatment, cv_cell[, .(across_drugs_r = median(r, na.rm = TRUE)),
    by = .(variant, lambda)], by = c('variant', 'lambda'))
  save(cv_metrics, 'validation_axes')
  per_treatment <- pred[, .(cells = .N, pearson_pct = cor_safe(observed_pct, predicted_pct),
    pearson_log2 = cor_safe(observed, predicted), observed_min_pct = min(observed_pct),
    observed_max_pct = max(observed_pct), predicted_min_pct = min(predicted_pct),
    predicted_max_pct = max(predicted_pct), sensitive = sum(observed_pct <= 50),
    sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
    false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)), by = .(model, drug, dose)]
  save(per_treatment, 'by_treatment')
  print(summary)
  descriptions <- c(anchored_actual = 'Anchored SVR', context_actual = 'Context + individual response',
    context_average = 'Context + average response', baseline_reference = 'Baseline + drug/dose',
    anchored_average = 'Anchored SVR: average response', previous_SVR = 'Previous calibrated SVR')
  colours <- c(anchored_actual = '#0072B2', context_actual = '#D55E00',
    context_average = '#009E73', baseline_reference = '#CC79A7')
  figure_dir <- file.path(root, 'figures')
  plot_save <- function(p, stem, width, height) {
    ggplot2::ggsave(file.path(figure_dir, paste0('contextual_', stem, '.pdf')), p, width = width, height = height)
    ggplot2::ggsave(file.path(figure_dir, paste0('contextual_', stem, '.png')), p, width = width, height = height, dpi = 140)
  }
  cv_long <- data.table::melt(cv_metrics, id.vars = c('variant', 'lambda'),
    measure.vars = c('across_cells_r', 'across_drugs_r'), variable.name = 'axis', value.name = 'correlation')
  p <- ggplot2::ggplot(cv_long, ggplot2::aes(lambda, correlation, colour = variant)) +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~axis, ncol = 2, scales = 'free_y') +
    ggplot2::scale_x_log10() + ggplot2::scale_colour_manual(values = c(actual = '#D55E00', average = '#009E73')) +
    ggplot2::theme_classic(base_size = 10) + ggplot2::theme(legend.position = 'bottom') +
    ggplot2::labs(x = 'Context penalty (larger = stronger shrinkage)', y = 'Training-cell validation Pearson r',
      colour = 'Response input', title = 'Do stronger baseline interactions generalize between training cells?',
      subtitle = 'Whole-cell validation predictions; frozen bank and previously selected SVR parameters',
      caption = 'Penalty selection uses balanced absolute and contrast errors; these correlations show its consequences.\nNo alternative penalties are screened on outer evaluation outcomes.')
  plot_save(p, 'validation_penalty', 8, 4.5)
  main <- pred[model %in% names(colours)]
  main[, model_label := factor(descriptions[model], levels = descriptions[names(colours)])]
  # Native abundance is essential for seeing calibration and compressed ranges.
  make_scatter <- function(rows, facets, title, subtitle) {
    ggplot2::ggplot(rows, ggplot2::aes(observed_pct, predicted_pct, colour = cell)) +
      ggplot2::geom_point(size = 1.2, alpha = .8) +
      ggplot2::geom_smooth(ggplot2::aes(group = 1), method = 'lm', formula = y ~ x,
        se = FALSE, colour = 'black', linewidth = .55) +
      ggplot2::geom_vline(xintercept = 50, linetype = 'dashed', colour = 'grey65') +
      ggplot2::geom_hline(yintercept = 50, linetype = 'dashed', colour = 'grey65') +
      ggplot2::facet_grid(facets, scales = 'free') +
      ggplot2::scale_colour_brewer(palette = 'Set1') +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom',
        plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
      ggplot2::labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)',
        colour = 'Held-out cell', title = title, subtitle = subtitle,
        caption = 'Each point is one held-out cell. Dashed lines mark 50% vehicle. No clipping or profile normalization.\nPrograms, SVR and context fitting exclude these cells; previously inspected outcomes are retrospective. PRISM doses are interpolated.')
  }
  focal_drugs <- c('lapatinib', 'crizotinib', 'alpelisib', 'trametinib', 'gefitinib', 'everolimus', 'bortezomib')
  grDevices::pdf(file.path(figure_dir, 'contextual_focal_drugs.pdf'), width = 12, height = 8.8)
  for (drug_name in focal_drugs) {
    rows <- main[drug == drug_name]
    stopifnot(nrow(rows) > 0L)
    annotation <- rows[, .(label = sprintf('n=%d; Pearson r=%.2f\nSensitive called %d/%d', .N,
      cor_safe(observed_pct, predicted_pct), sum(observed_pct <= 50 & predicted_pct <= 50),
      sum(observed_pct <= 50))), by = .(dose, model_label)]
    p <- make_scatter(rows, dose ~ model_label, paste('Held cells |', drug_name),
      'Columns compare decoders; rows separate doses (uM)') +
      ggplot2::geom_text(data = annotation, ggplot2::aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = 1.03, vjust = 1.15, size = 2.2)
    print(p)
    plot_save(p, paste0('focal_', drug_name), 12, 8.8)
  }
  grDevices::dev.off()
  # Across drugs, within each cell: retain the dose and cell axes separately.
  rows <- main[drug_set == 'Seen drugs']
  p <- ggplot2::ggplot(rows, ggplot2::aes(observed_pct, predicted_pct, colour = model)) +
    ggplot2::geom_point(size = .55, alpha = .3) +
    ggplot2::geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .6) +
    ggplot2::facet_wrap(~cell, ncol = 3, scales = 'free') +
    ggplot2::scale_colour_manual(values = colours, labels = descriptions[names(colours)]) +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom') +
    ggplot2::labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)', colour = NULL,
      title = 'Across familiar drugs within each held-out cell')
  plot_save(p, 'by_cell', 10, 10)
  # Center within held-out drug-dose groups only for descriptive visualization;
  # these evaluation means never enter model fitting or prediction.
  p <- ggplot2::ggplot(main, ggplot2::aes(observed_pp_centered, predicted_pp_centered, colour = cell)) +
    ggplot2::geom_point(size = .65, alpha = .35) +
    ggplot2::geom_smooth(ggplot2::aes(group = 1), method = 'lm', formula = y ~ x,
      se = FALSE, colour = 'black', linewidth = .6) +
    ggplot2::facet_grid(drug_set ~ model_label, scales = 'free') +
    ggplot2::scale_colour_brewer(palette = 'Set1') +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom') +
    ggplot2::labs(x = 'Observed difference from drug-dose mean (percentage points)',
      y = 'Predicted difference from drug-dose mean (percentage points)', colour = 'Held-out cell',
      title = 'Across cells receiving the same drug and dose',
      subtitle = 'Evaluation-only centering removes differences in average drug potency')
  plot_save(p, 'within_treatment', 12, 6.8)
  # Show whether the learned adjustment moves predictions in the needed direction.
  comp <- data.table::fread(file.path(root, 'data/contextual_components.tsv'))
  jj <- match(comp$observation_id, labels$observation_id)
  comp[, `:=`(needed_pp = labels$prism_abundance_pct[jj] - 100 * 2^shared,
    applied_pp = 100 * 2^predicted - 100 * 2^shared,
    cell = recipient, drug_set = ifelse(drug %in% settings$excluded_drugs, 'Excluded family', 'Seen drugs'))]
  save(comp, 'component_evaluation')
  p <- ggplot2::ggplot(comp, ggplot2::aes(needed_pp, applied_pp, colour = cell)) +
    ggplot2::geom_point(size = .65, alpha = .4) +
    ggplot2::geom_smooth(ggplot2::aes(group = 1), method = 'lm', formula = y ~ x,
      se = FALSE, colour = 'black', linewidth = .6) +
    ggplot2::facet_grid(drug_set ~ variant, scales = 'free') +
    ggplot2::scale_colour_brewer(palette = 'Set1') +
    ggplot2::theme_classic(base_size = 10) + ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom') +
    ggplot2::labs(x = 'Adjustment needed from shared SVR contribution (percentage points)',
      y = 'Adjustment supplied by baseline interactions (percentage points)', colour = 'Held-out cell',
      title = 'Does baseline conditioning correct the remaining errors?',
      subtitle = 'Positive adjustments predict more remaining abundance; negative adjustments predict greater sensitivity')
  plot_save(p, 'adjustment', 8.5, 8)
  # Full coefficient matrices retain their magnitudes and share a common order.
  source('/work/src/R/lib_heatmap.R')
  coefficients <- data.table::fread(file.path(root, 'data/contextual_coefficients_actual.tsv'))
  programs <- sort(unique(coefficients$baseline_program))
  matrices <- lapply(c('actual', 'average'), function(variant) {
    tab <- data.table::fread(file.path(root, 'data', paste0('contextual_coefficients_', variant, '.tsv')))
    mat <- matrix(NA_real_, length(programs), length(programs), dimnames = list(programs, programs))
    mat[cbind(match(tab$baseline_program, programs), match(tab$displaced_program, programs))] <- tab$coefficient
    stopifnot(all(is.finite(mat)))
    mat
  })
  ro <- stats::hclust(stats::dist(matrices[[1L]]), method = 'ward.D2')$order
  co <- stats::hclust(stats::dist(t(matrices[[1L]])), method = 'ward.D2')$order
  limit <- max(abs(unlist(matrices)))
  if (limit == 0) limit <- 1
  grDevices::pdf(file.path(figure_dir, 'contextual_interaction_weights.pdf'), width = 7.8, height = 7)
  for (k in seq_along(matrices)) {
    h <- draw_heatmap(matrices[[k]][ro, co], scale_rows = FALSE, show_row_names = TRUE,
      show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
      value_range = c(-limit, 0, limit), value_colors = c('#2166AC', '#FFFFFF', '#B2182B'),
      hm_name = 'Log2 effect', title = paste(c('Individual response', 'Average response')[k],
        '| rows: baseline; columns: displacement'))
    ComplexHeatmap::draw(h)
  }
  grDevices::dev.off()
  # A compact table-like figure keeps the two predictive axes distinct.
  long <- data.table::melt(summary, id.vars = c('model', 'drug_set'),
    measure.vars = c('median_fixed_dose_within_cell_r', 'within_treatment_r_pct'),
    variable.name = 'axis', value.name = 'correlation')
  long[, `:=`(model_label = factor(descriptions[model], levels = rev(descriptions)),
    axis = factor(axis, levels = c('median_fixed_dose_within_cell_r', 'within_treatment_r_pct'),
      labels = c('Across drugs within a cell, fixed dose\nMedian Pearson r (log2 abundance)',
        'Across cells, fixed drug and dose\nPearson r (percentage-point contrasts)')))]
  long[, `:=`(label_x = ifelse(is.finite(correlation), correlation, 0),
    label = ifelse(is.finite(correlation), sprintf('%.2f', correlation), 'No cell differences'))]
  p <- ggplot2::ggplot(long, ggplot2::aes(correlation, model_label)) +
    ggplot2::geom_point(data = long[is.finite(correlation)], size = 2, colour = '#0072B2') +
    ggplot2::geom_text(ggplot2::aes(x = label_x, label = label), hjust = -.2, size = 2.6) +
    ggplot2::facet_grid(drug_set ~ axis, scales = 'free_x') +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.1, .3))) +
    ggplot2::theme_classic(base_size = 9) + ggplot2::labs(x = NULL, y = NULL,
      title = 'Contextual decoder | two distinct predictive comparisons')
  plot_save(p, 'summary', 10.5, 6)
}

# Direct pairwise scatterplots for inspecting patterns by cell and drug.
plot_response_distance_scatter <- function(root) {
  pairs <- data.table::fread(file.path(root, 'data/diagnostic_geometry_pairs.tsv.gz'))
  focal_drugs <- c('lapatinib', 'crizotinib', 'alpelisib', 'trametinib', 'gefitinib', 'everolimus', 'bortezomib')
  pairs[, dose_label := factor(dose, levels = sort(unique(dose)))]
  colours <- c('#0072B2', '#D55E00', '#009E73')
  stopifnot(length(levels(pairs$dose_label)) <= length(colours))
  figure_dir <- file.path(root, 'figures')
  correlations <- pairs[, .(pairs = .N,
    pearson_r = stats::cor(displacement_distance, outcome_gap_pp)),
    by = .(comparison, cohort, unit, dose)]
  data.table::fwrite(correlations, file.path(root, 'data/diagnostic_geometry_scatter_correlations.tsv'), sep = '\t')
  # Reused for cell and drug panels; every point is an actual eligible pair.
  make_plot <- function(rows, title, subtitle, columns = 3L) {
    annotation <- rows[, .(label = sprintf('n = %s; Pearson r = %.2f',
      format(.N, big.mark = ','), stats::cor(displacement_distance, outcome_gap_pp))), by = unit]
    ggplot2::ggplot(rows, ggplot2::aes(displacement_distance, outcome_gap_pp, colour = dose_label)) +
      ggplot2::geom_point(size = .5, alpha = .35, stroke = 0) +
      ggplot2::geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .65) +
      ggplot2::geom_text(data = annotation, ggplot2::aes(x = -Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = -.03, vjust = 1.3, size = 2.5) +
      ggplot2::facet_wrap(~unit, ncol = columns, scales = 'free') +
      ggplot2::scale_colour_manual(values = colours) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.03, .16))) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom',
        strip.background = ggplot2::element_rect(fill = 'grey95', colour = NA),
        plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
      ggplot2::labs(x = 'Distance between measured displacement profiles',
        y = 'Absolute sensitivity difference (percentage points of vehicle)',
        colour = 'Dose (uM)', title = title, subtitle = subtitle,
        caption = paste('Each dot is a pair at the same dose. All points retained; axes vary independently by panel.',
          'Distance: Euclidean / sqrt(40), after dividing each program by its training SD. No profile normalization.',
          'Lines: linear trends by dose. Pearson r in each panel pools doses. Pairs share observations; descriptive, retrospective.', sep = '\n'))
  }
  held <- pairs[comparison == 'same_cell_different_drugs' & cohort == 'Held cells']
  p <- make_plot(held, 'Within each cell: molecular distance versus sensitivity difference',
    'Nine held-out cell lines | each dot compares two different drugs at the same dose')
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_geometry_scatter_cells.pdf'), p, width = 10, height = 10)
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_geometry_scatter_cells.png'), p, width = 10, height = 10, dpi = 150)
  # All cell lines, with held and development cells explicitly identified.
  cells <- pairs[comparison == 'same_cell_different_drugs']
  cells[, unit := paste0(unit, ifelse(cohort == 'Held cells', ' [held]', ' [development]'))]
  cell_names <- sort(unique(cells$unit))
  grDevices::pdf(file.path(figure_dir, 'diagnostic_geometry_scatter_all_cells.pdf'), width = 10, height = 10)
  for (start in seq.int(1L, length(cell_names), by = 9L)) {
    names_page <- cell_names[seq.int(start, min(start + 8L, length(cell_names)))]
    print(make_plot(cells[unit %in% names_page], 'Within each cell: different drugs, matched dose',
      'All eligible pairs | colour and fitted line indicate dose'))
  }
  grDevices::dev.off()
  drugs <- pairs[comparison == 'same_treatment_different_cells']
  focal <- drugs[drug_i %in% focal_drugs]
  stopifnot(nrow(focal) > 0L)
  p <- make_plot(focal, 'Within each drug: molecular distance versus sensitivity difference',
    'Seven focal drugs | each dot compares two cells at the same drug and dose')
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_geometry_scatter_focal_drugs.pdf'), p, width = 10, height = 10)
  ggplot2::ggsave(file.path(figure_dir, 'diagnostic_geometry_scatter_focal_drugs.png'), p, width = 10, height = 10, dpi = 150)
  drug_names <- sort(unique(drugs$unit))
  grDevices::pdf(file.path(figure_dir, 'diagnostic_geometry_scatter_all_drugs.pdf'), width = 10, height = 10)
  for (start in seq.int(1L, length(drug_names), by = 9L)) {
    names_page <- drug_names[seq.int(start, min(start + 8L, length(drug_names)))]
    print(make_plot(drugs[unit %in% names_page], 'Within each drug: different cells, matched dose',
      'Pairs stay within held/development cohorts | colour and fitted line indicate dose'))
  }
  grDevices::dev.off()
}

# Systematic molecular-distance audit. All summaries are descriptive; pairs
# share observations. Fixed drug-dose/cell-dose strata receive equal weight.
plot_systematic_response_geometry <- function(root) {
  trend <- data.table::fread(file.path(root, 'data/diagnostic_geometry_distance_trend.tsv'))
  nearest <- data.table::fread(file.path(root, 'data/diagnostic_geometry_nearest_summary.tsv'))
  units <- data.table::fread(file.path(root, 'data/diagnostic_geometry_by_cell_drug.tsv'))
  dose <- data.table::fread(file.path(root, 'data/diagnostic_geometry_distance_by_dose.tsv'))
  comparison_names <- c(same_treatment_different_cells = 'Same drug + dose, different cells',
    same_cell_different_drugs = 'Same cell + dose, different drugs')
  representation_names <- c(displacement_distance = 'Displacement', post_state_distance = 'Vehicle baseline + displacement')
  for (tab in list(trend, nearest, units, dose)) {
    tab[, comparison_label := factor(comparison_names[comparison], levels = comparison_names)]
  }
  colours <- c(displacement_distance = '#0072B2', post_state_distance = '#D55E00')
  family <- sub('.*_(MEK|ERBB)_holdout$', '\\1', root)
  p <- ggplot2::ggplot(trend, ggplot2::aes(decile, median_gap_pp, colour = representation)) +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_grid(comparison_label ~ cohort) +
    ggplot2::scale_colour_manual(values = colours, labels = representation_names) +
    ggplot2::scale_x_continuous(breaks = c(1, 5, 10), labels = c('Closest 10%', 'Middle', 'Farthest 10%')) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom',
      strip.text.y = ggplot2::element_text(angle = 0),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::labs(x = 'Molecular distance within each fixed comparison',
      y = 'Median sensitivity difference\n(percentage points of vehicle)', colour = NULL,
      title = paste(family, 'bank | Do closer molecular profiles have more similar sensitivity?'),
      subtitle = 'All eligible pairs; distance uses all 40 programs, preserving direction and amplitude',
      caption = paste('Each drug-dose or cell-dose stratum contributes equally. Lines are medians of stratum medians.',
        'Orange uses same-assay mean vehicle expression + displacement, an approximate post-treatment state.',
        'Within a cell, adding the same baseline leaves distances unchanged: the two lines overlap exactly.',
        'PRISM is interpolated to molecular doses. Pairs reuse observations; this is a retrospective descriptive audit.', sep = '\n'))
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_systematic.pdf'), p, width = 10.5, height = 8)
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_systematic.png'), p, width = 10.5, height = 8, dpi = 140)
  q <- ggplot2::ggplot(nearest, ggplot2::aes(representation, median_gap_pp, fill = representation)) +
    ggplot2::geom_col(width = .6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf('%.1f pp\n%d/%d discordant',
      median_gap_pp, discordant, observations)), vjust = -.2, size = 2.7) +
    ggplot2::facet_grid(comparison_label ~ cohort) +
    ggplot2::scale_fill_manual(values = colours, guide = 'none') +
    ggplot2::scale_x_discrete(labels = c(displacement_distance = 'Displacement',
      post_state_distance = 'Baseline +\ndisplacement')) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .3))) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
    ggplot2::labs(x = NULL, y = 'Median sensitivity difference (percentage points)',
      title = 'Sensitivity difference from the closest molecular neighbour',
      subtitle = 'One neighbour per observation, selected using molecular measurements only',
      caption = 'Discordant: one sample at or below 50% vehicle and the other above 50%. Counts are directed observations.\nThe same pair can serve as neighbour for both samples. Baseline + displacement is an approximate same-assay state.')
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_state_comparison.pdf'), q, width = 10.5, height = 7)
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_state_comparison.png'), q, width = 10.5, height = 7, dpi = 140)
  # Named breakdowns include every unit with at least one matched-profile pair.
  # Pages are ordered alphabetically rather than selecting on outcome differences.
  matched <- units[selection == 'Matched profiles']
  for (comparison_name in unique(matched$comparison)) {
    tab <- matched[comparison == comparison_name]
    data.table::setorder(tab, cohort, unit)
    tab[, row_label := sprintf('%s [%s]  n=%d; discordant=%d', unit,
      ifelse(cohort == 'Held cells', 'held', 'development'), pairs, discordant)]
    stem <- if (comparison_name == 'same_cell_different_drugs') 'by_cell' else 'by_drug'
    grDevices::pdf(file.path(root, 'figures', paste0('diagnostic_geometry_', stem, '.pdf')), width = 9.5, height = 8)
    for (start in seq.int(1L, nrow(tab), by = 30L)) {
      page <- data.table::copy(tab[seq.int(start, min(start + 29L, nrow(tab)))])
      page[, row_label := factor(row_label, levels = rev(row_label))]
      chart <- ggplot2::ggplot(page, ggplot2::aes(100 * discordant_fraction, row_label)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 100 * discordant_fraction, yend = row_label),
          colour = 'grey75', linewidth = .4) +
        ggplot2::geom_point(colour = '#0072B2', size = 2) +
        ggplot2::scale_x_continuous(limits = c(0, 100)) +
        ggplot2::theme_classic(base_size = 9) +
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8)) +
        ggplot2::labs(x = 'Similar-profile pairs with opposite sensitive/resistant classes (%)', y = NULL,
          title = comparison_names[comparison_name],
          subtitle = 'Direction cosine >= 0.9 and amplitude ratio >= 0.8; doses compared separately',
          caption = 'n is the number of matched-profile pairs across doses; pairs share observations.\nSensitive <=50% vehicle. Units with no matched-profile pairs are absent; full denominators are in the source table.')
      print(chart)
    }
    grDevices::dev.off()
  }
  # Dose-specific curves show whether the pooled pattern is driven by one dose.
  p_dose <- ggplot2::ggplot(dose[representation == 'displacement_distance'],
    ggplot2::aes(decile, median_gap_pp, colour = factor(dose))) +
    ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 1) +
    ggplot2::facet_grid(comparison_label ~ cohort) +
    ggplot2::theme_classic(base_size = 10) + ggplot2::theme(aspect.ratio = 1, legend.position = 'bottom') +
    ggplot2::labs(x = 'Displacement distance decile (1 = closest)', y = 'Median sensitivity difference (percentage points)',
      colour = 'Dose (uM)', title = 'Systematic comparison separated by dose')
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_by_dose.pdf'), p_dose, width = 10.5, height = 7)
}

# Display the actual frozen displacement profiles alongside the assay outcomes.
# The appendix contains every audited example, including held-cell examples.
plot_response_geometry <- function(root) {
  stopifnot(requireNamespace('patchwork', quietly = TRUE))
  examples <- data.table::fread(file.path(root, 'data/diagnostic_geometry_example_resolution.tsv'))
  prepared <- qs::qread(file.path(root, 'data/prepared_scores.qs'))
  response <- prepared$measured
  identifiers <- prepared$metadata$observation_id
  stopifnot(!anyDuplicated(identifiers), nrow(response) == length(identifiers))
  values <- data.table::rbindlist(lapply(seq_len(nrow(examples)), function(k) {
    pair <- examples[k]
    ii <- match(c(pair$observation_i, pair$observation_j), identifiers)
    stopifnot(!anyNA(ii))
    data.table::data.table(example = k, program = rep(colnames(response), 2L),
      member = rep(c('A', 'B'), each = ncol(response)),
      observation_id = rep(identifiers[ii], each = ncol(response)),
      displacement = as.vector(t(response[ii, , drop = FALSE])))
  }))
  data.table::fwrite(values, file.path(root, 'data/diagnostic_geometry_plot_values.tsv'), sep = '\t')
  colours <- c(A = '#0072B2', B = '#D55E00')
  panels <- lapply(seq_len(nrow(examples)), function(k) {
    pair <- examples[k]
    profile <- values[example == k]
    # Order programs by their average displacement; both traces use that order.
    order <- profile[, .(average = mean(displacement)), by = program][order(average)]$program
    profile[, position := match(program, order)]
    title <- sprintf('%s | %s', pair$example_type, pair$cohort)
    labels <- c(sprintf('A: %s / %s', pair$cell_i, pair$drug_i),
      sprintf('B: %s / %s', pair$cell_j, pair$drug_j))
    p <- ggplot2::ggplot(profile, ggplot2::aes(position, displacement, colour = member)) +
      ggplot2::geom_hline(yintercept = 0, colour = 'grey75', linewidth = .3) +
      ggplot2::geom_line(linewidth = .6) + ggplot2::geom_point(size = 1) +
      ggplot2::scale_colour_manual(values = colours, labels = labels) +
      ggplot2::scale_x_continuous(breaks = seq_along(order), labels = order) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(legend.position = 'bottom', legend.text = ggplot2::element_text(size = 8),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5, size = 6)) +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
      ggplot2::labs(x = 'All 40 programs (same order for A and B)',
        y = 'Measured gene displacement\n(hard cluster mean)', title = title,
        subtitle = sprintf('Dose: %g uM | direction cosine %.3f | amplitude ratio %.3f | raw-gene Pearson r %.3f',
          pair$dose, pair$direction_cosine, pair$magnitude_ratio, pair$raw_gene_pearson), colour = NULL)
    outcomes <- data.table::data.table(member = factor(c('A', 'B'), levels = c('A', 'B')),
      abundance = c(pair$abundance_i, pair$abundance_j))
    q <- ggplot2::ggplot(outcomes, ggplot2::aes(member, abundance, fill = member)) +
      ggplot2::geom_col(width = .6) +
      ggplot2::geom_hline(yintercept = 50, linetype = 'dashed', colour = 'grey40') +
      ggplot2::geom_text(ggplot2::aes(label = sprintf('%.1f%%', abundance)), vjust = -.5, size = 3) +
      ggplot2::scale_fill_manual(values = colours, guide = 'none') +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .18))) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::labs(x = NULL, y = 'PRISM abundance (% vehicle)', title = 'Functional outcome')
    patchwork::wrap_plots(p, q, nrow = 1, widths = c(4, 1))
  })
  # Three explicit illustrations, not a random sample or validation estimate.
  selected <- c(
    which(examples$cell_i == 'SK-MEL-2' & examples$cell_j == 'SW480' & examples$drug_i == 'clofarabine'),
    which(examples$cell_i == 'AsPC-1' & examples$drug_i == 'everolimus' & examples$drug_j == 'gemcitabine'),
    which(examples$cell_i == 'NCI-H1573' & examples$drug_i == 'homoharringtonine' & examples$drug_j == 'methotrexate'))
  stopifnot(length(selected) == 3L)
  combined <- patchwork::wrap_plots(panels[selected], ncol = 1) +
    patchwork::plot_annotation(caption = paste(
      'Retrospective examples; not evidence of causality or prediction accuracy. Lower PRISM abundance means greater sensitivity.',
      'Profiles are expression changes, not functional loadings. No per-profile normalization. Program order varies by pair.',
      'Cosine and amplitude ratio use training-SD-scaled coordinates. PRISM values are interpolated to molecular doses.',
      sep = '\n'))
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_examples.pdf'), combined, width = 11, height = 9)
  ggplot2::ggsave(file.path(root, 'figures/diagnostic_geometry_examples.png'), combined, width = 11, height = 9, dpi = 140)
  grDevices::pdf(file.path(root, 'figures/diagnostic_geometry_all_examples.pdf'), width = 11, height = 3.5)
  for (p in panels) print(p)
  grDevices::dev.off()
}

# Frozen contextual decoder audit: named discordant cases, source baseline,
# molecular resolution and training-versus-excluded-cell behavior.
plot_contextual_audit <- function(root) {
  read_table <- function(name) data.table::fread(file.path(root, 'data', paste0('contextual_audit_', name, '.tsv')))
  cases <- read_table('cases')
  programs <- read_table('case_programs')
  genes <- read_table('case_response_genes')
  baseline_genes <- read_table('case_baseline_genes')
  neighbors <- read_table('case_training_neighbors')
  metrics <- read_table('fold_metrics')
  figure_dir <- file.path(root, 'figures')
  palette <- c('i' = '#0072B2', 'j' = '#D55E00')
  theme_audit <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = 'bottom')
  # Plot matched model paths on fitting recipients and recipients excluded from
  # that fold. This reveals overfitting without rescreening any held-cell penalty.
  path <- data.table::melt(metrics, id.vars = c('fold', 'lambda', 'cohort'),
    measure.vars = c('across_cells_r', 'context_residual_r'), variable.name = 'axis', value.name = 'r')
  path[, axis := factor(axis, levels = c('across_cells_r', 'context_residual_r'),
    labels = c('Same drug/dose: cell differences', 'Context contribution versus required correction'))]
  p <- ggplot2::ggplot(path, ggplot2::aes(lambda, r, color = cohort, group = interaction(cohort, fold))) +
    ggplot2::geom_hline(yintercept = 0, color = 'grey70') + ggplot2::geom_line(alpha = .65) +
    ggplot2::geom_point(size = 1.3) + ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~axis, nrow = 1) +
    ggplot2::labs(x = 'Interaction penalty (larger = more shrinkage)', y = 'Pearson r', color = NULL,
      title = 'Why the context terms shrink: three whole-cell development folds') +
    theme_audit + ggplot2::theme(aspect.ratio = 1)
  ggplot2::ggsave(file.path(figure_dir, 'contextual_audit_training_transfer.pdf'), p, width = 8.4, height = 4.6)
  ggplot2::ggsave(file.path(figure_dir, 'contextual_audit_training_transfer.png'), p, width = 8.4, height = 4.6, dpi = 150)
  agreement <- read_table('baseline_source_agreement')
  agreement[, program := factor(program, levels = program[order(across_cell_ccle_vehicle_r)])]
  p <- ggplot2::ggplot(agreement, ggplot2::aes(across_cell_ccle_vehicle_r, program)) +
    ggplot2::geom_vline(xintercept = 0, color = 'grey70') +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = across_cell_ccle_vehicle_r, yend = program), color = '#0072B2') +
    ggplot2::geom_point(color = '#0072B2', size = 1.5) +
    ggplot2::labs(x = 'Pearson r across training cells', y = NULL,
      title = 'Do CCLE program baselines agree with Tahoe vehicle baselines?',
      subtitle = 'Same frozen genes and hard means; independent expression sources') + theme_audit
  ggplot2::ggsave(file.path(figure_dir, 'contextual_audit_baseline_source_agreement.pdf'), p, width = 6.7, height = 7.4)
  # Retain every pre-existing case in one PDF, including the convergent examples.
  # These are selected illustrations, not new validation estimates.
  grDevices::pdf(file.path(figure_dir, 'contextual_audit_cases.pdf'), width = 10, height = 7.8)
  for (case_id in cases$case) {
    one <- cases[case == case_id]
    pp <- programs[case == case_id]
    endpoint_names <- c(i = paste(one$cell_i, one$drug_i), j = paste(one$cell_j, one$drug_j))
    values <- data.table::data.table(endpoint = rep(c('i', 'j'), each = 3L),
      measurement = rep(c('Measured\nPRISM', 'Response\nprediction', '+ baseline\ncontext'), 2L),
      abundance = c(one$observed_i, one$shared_i, one$predicted_i,
        one$observed_j, one$shared_j, one$predicted_j))
    values[, measurement := factor(measurement,
      levels = c('Measured\nPRISM', 'Response\nprediction', '+ baseline\ncontext'))]
    a <- ggplot2::ggplot(values, ggplot2::aes(measurement, abundance, fill = endpoint)) +
      ggplot2::geom_col(position = 'dodge', width = .75, show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = 50, linetype = 2, color = 'grey40') +
      ggplot2::scale_fill_manual(values = palette, labels = endpoint_names) +
      ggplot2::labs(x = NULL, y = 'Abundance (% vehicle)', fill = NULL,
        title = sprintf('%s | %g uM | %s', case_id, one$dose, one$cohort),
        subtitle = if (one$cohort == 'Held cells') 'Frozen predictions in held cells' else 'Whole-cell development OOF') +
      theme_audit + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7))
    dd <- data.table::melt(pp, id.vars = 'program',
      measure.vars = c('displacement_i', 'displacement_j'), variable.name = 'endpoint', value.name = 'value')
    dd[, endpoint := sub('displacement_', '', endpoint)]
    dd[, index := match(program, pp$program)]
    b <- ggplot2::ggplot(dd, ggplot2::aes(index, value, color = endpoint)) +
      ggplot2::geom_hline(yintercept = 0, color = 'grey75') + ggplot2::geom_line(show.legend = FALSE) +
      ggplot2::scale_color_manual(values = palette, labels = endpoint_names) +
      ggplot2::labs(x = 'Frozen program (P01-P40)', y = 'Raw displacement', color = NULL,
        title = sprintf('Difference = %.3g; vehicle median = %.3g',
          one$program_difference, one$vehicle_difference_median)) + theme_audit
    bb <- data.table::melt(pp, id.vars = 'program',
      measure.vars = c('baseline_ccle_z_i', 'baseline_ccle_z_j'), variable.name = 'endpoint', value.name = 'value')
    bb[, endpoint := sub('baseline_ccle_z_', '', endpoint)]
    bb[, index := match(program, pp$program)]
    c <- ggplot2::ggplot(bb, ggplot2::aes(index, value, color = endpoint)) +
      ggplot2::geom_hline(yintercept = 0, color = 'grey75') + ggplot2::geom_line() +
      ggplot2::scale_color_manual(values = palette, labels = endpoint_names) +
      ggplot2::labs(x = 'Frozen program (P01-P40)', y = 'CCLE baseline (training-cell z)', color = NULL,
        title = sprintf('Baseline difference = %.2f training SD', one$ccle_baseline_rms_z)) + theme_audit
    # Per-gene vehicle SD is a background scale, not a significance test or
    # treated-replicate error. Display the largest differences explicitly.
    gg <- genes[case == case_id][1:15]
    gg[, label := paste0(gene, ifelse(in_bank, paste0(' [', program, ']'), ' [outside bank]'))]
    gg[, label := factor(label, levels = rev(label))]
    e <- ggplot2::ggplot(gg, ggplot2::aes(difference_vehicle_sd, label)) +
      ggplot2::geom_vline(xintercept = c(-2, 0, 2), color = 'grey70', linetype = c(2, 1, 2)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = difference_vehicle_sd, yend = label), color = '#0072B2') +
      ggplot2::geom_point(color = '#0072B2', size = 1.2) +
      ggplot2::labs(x = 'Response difference / pooled vehicle SD', y = NULL,
        title = sprintf('%d/%d genes exceed 2 vehicle SD', one$genes_above_two_vehicle_sd, one$assayed_common_genes)) +
      theme_audit
    page <- (patchwork::wrap_plots(a, b, c, e, ncol = 2, guides = 'collect') &
      ggplot2::theme(legend.position = 'bottom')) +
      patchwork::plot_annotation(caption = paste('Raw program values; no row scaling. Gene ranking is descriptive.',
        'Vehicle contrasts are background controls, not treated replicates.'))
    print(page)
    if (case_id %in% c('C01', 'C17', 'C19', 'C25', 'C34')) {
      ggplot2::ggsave(file.path(figure_dir, paste0('contextual_audit_', case_id, '.png')),
        page, width = 10, height = 7.8, dpi = 150)
    }
  }
  grDevices::dev.off()
  # Named context genes retain their source expression scale and exact pair order.
  source('/work/src/R/lib_heatmap.R')
  selected_cases <- c('C01', 'C17', 'C19', 'C25', 'C34')
  grDevices::pdf(file.path(figure_dir, 'contextual_audit_baseline_genes.pdf'), width = 9.5, height = 7.2)
  for (case_id in selected_cases) {
    one <- cases[case == case_id]
    tab <- baseline_genes[case == case_id]
    matrix <- as.matrix(tab[, .(ccle_i, ccle_j, vehicle_i, vehicle_j)])
    rownames(matrix) <- tab$gene
    order <- stats::hclust(stats::dist(matrix), method = 'ward.D2')$order
    matrix <- matrix[order, , drop = FALSE]
    colnames(matrix) <- c(one$cell_i, one$cell_j, one$cell_i, one$cell_j)
    maps <- lapply(1:2, function(source_index) {
      mm <- matrix[, (source_index - 1L) * 2L + 1:2, drop = FALSE]
      # Independent source scales avoid conflating TPM and CP10K expression.
      draw_heatmap(mm, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = source_index == 1L, show_column_names = TRUE,
        hm_name = if (source_index == 1L) 'log2(TPM+1)' else 'log2(CP10K+1)',
        value_range = c(min(mm), mean(range(mm)), max(mm)),
        value_colors = c('#f7fbff', '#6baed6', '#08306b'),
        title = paste(case_id, if (source_index == 1L) '| CCLE baseline' else '| Tahoe vehicle baseline'),
        row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 8))
    })
    ComplexHeatmap::draw(maps[[1L]] + maps[[2L]])
  }
  grDevices::dev.off()
  # A concrete omitted-context example: plot cached BRAF hotspot flags without
  # fitting a genotype predictor. Hotspot presence is not an exact variant call.
  driver <- read_table('vemurafenib_braf_context')
  driver[, label := paste0(recipient, ' [', ifelse(cohort == 'Training cells', 'train', 'held'), ']')]
  driver[, label := factor(label, levels = label[order(prism_abundance_pct)])]
  p <- ggplot2::ggplot(driver, ggplot2::aes(prism_abundance_pct, label,
    color = factor(hotspot_flag), shape = cohort)) +
    ggplot2::geom_vline(xintercept = 50, linetype = 2, color = 'grey65') +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c('0' = '#0072B2', '1' = '#D55E00'),
      labels = c('0' = 'No hotspot flag', '1' = 'Hotspot flag')) +
    ggplot2::labs(x = 'PRISM abundance (% vehicle)', y = NULL, color = 'BRAF', shape = NULL,
      title = 'Vemurafenib 0.5 uM: genotype context absent from the decoder',
      subtitle = 'Descriptive audit; no genotype model fitted') + theme_audit
  ggplot2::ggsave(file.path(figure_dir, 'contextual_audit_vemurafenib_context.pdf'), p, width = 7.5, height = 7)
  ggplot2::ggsave(file.path(figure_dir, 'contextual_audit_vemurafenib_context.png'), p, width = 7.5, height = 7, dpi = 150)
  data.table::fwrite(neighbors[case %in% selected_cases],
    file.path(root, 'data/contextual_audit_focus_neighbors.tsv'), sep = '\t')
}

if (sys.nframe() == 0L) {
  suppressPackageStartupMessages({library(data.table); library(glmnet)})
  held_cells_only <- identical(commandArgs(TRUE), '--held-cells')
  stopifnot(length(commandArgs(TRUE)) == 0L || held_cells_only)
  panel <- qs::qread('/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
  for (family in c('MEK', 'ERBB')) {
    root <- file.path('/work/results/program_pilot',
      paste0('20_multidrug_baseline_response_function_', family, '_holdout'))
    if (held_cells_only) {
      plot_held_cell_generalization(root, panel)
    } else {
      pred <- fread(file.path(root, 'data/functional_predictions.tsv.gz'))
      plot_family_transfer(pred, file.path(root, 'figures'), file.path(root, 'data'))
      plot_family_diagnostics(root, panel)
    }
    message('Updated comparisons: ', family)
  }
}
