#!/usr/bin/env Rscript
# Plot variants of v79: all three measured PRISM endpoints, with validation folds
# visible. This is a diagnostic view of the existing fit, not another model.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
root <- '/work/results/program_pilot/79_displacement_only_within_treatment_bridge'
single <- fread(file.path(root, 'data/held_single_predictions_frozen.tsv'))
metrics <- fread(file.path(root, 'data/fixed_drug_dose_single_metrics.tsv'))
points <- copy(single[family == 'Selected from training'])
points[, fold_label := factor(fold, levels = 1:3, labels = paste('Fold', 1:3))]
doses <- sort(unique(points$dose))
stopifnot(length(doses) == 3L, !anyDuplicated(points$single_id),
  all(points[drug == 'Lapatinib ditosylate', .N, by = dose]$N == 441L))
fold_summary <- points[, .(cells = .N, mean_prediction_log2 = mean(predicted_log2),
  average_component_log2 = unique(average_log2), sensitivity_sd_log2 = sd(sensitivity_log2)),
  by = .(drug, dose, fold)]
fwrite(fold_summary, file.path(root, 'data/all_dose_fold_calibration.tsv'), sep = '\t')
fwrite(dcast(metrics, drug + dose + cells ~ family, value.var = 'pearson_log2'),
  file.path(root, 'data/all_dose_model_correlations.tsv'), sep = '\t')
fwrite(points, file.path(root, 'data/all_dose_plot_points.tsv'), sep = '\t')

plots <- list()
for (dose_value in doses) {
  plot_data <- points[dose == dose_value]
  labels <- metrics[family == 'Selected from training' & dose == dose_value,
    .(drug, label = sprintf('Pooled r = %.3f; n = %d', pearson_log2, cells))]
  p <- ggplot(plot_data, aes(observed_log2, predicted_log2)) +
    geom_point(aes(color = fold_label), size = .65, alpha = .55) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .45, color = 'black') +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -.05, vjust = 1.3, size = 2.65) +
    facet_wrap(~drug, ncol = 5, scales = 'free') +
    scale_color_manual(values = c('#2378a0', '#bd5a28', '#72419d'), name = 'Validation fit') +
    theme_bw(base_size = 9) + theme(aspect.ratio = 1, strip.text = element_text(size = 8),
      legend.position = 'bottom', plot.subtitle = element_text(size = 9)) +
    labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
      title = sprintf('Measured PRISM endpoint: %g µM', dose_value),
      subtitle = 'Architecture selected inside each training fold; coloured bands expose calibration artifacts')
  name <- paste0('all_drugs_', gsub('\\.', '_', format(dose_value, scientific = FALSE)), 'uM_by_fold')
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = 12, height = 5.8, dpi = 150)
  plots[[as.character(dose_value)]] <- p
}
pdf(file.path(root, 'figures/all_drugs_all_measured_doses_by_fold.pdf'), width = 12, height = 5.8)
for (p in plots) print(p)
dev.off()
cat('Three dose pages saved; all441 lapatinib points retained on every page.\n')
