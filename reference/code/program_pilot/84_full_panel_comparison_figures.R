#!/usr/bin/env Rscript
# Present the frozen three-source predictions; this script never fits a model.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setDTthreads(1L)
root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
reference <- '/work/results/program_pilot/82_joint_program_rules_and_pair_interactions'
families <- c('Tahoe only', 'Tahoe + LINCS')
wilson <- fread(file.path(root, 'data/wilson_1uM_comparison.tsv'))
old_queries <- fread(file.path(reference, 'data/au565_hgf_predictions_frozen.tsv'))[
  family == 'Linked molecular bridge' & hgf_context == 'Equal-context mean' & dose == 1]
old_wilson <- rbindlist(list(
  old_queries[, .(family = 'Tahoe only', condition = 'Lapatinib', dose,
    predicted_abundance = lapatinib_abundance_pct, observed_abundance = 36.5, observed_sd = 4.9)],
  old_queries[, .(family = 'Tahoe only', condition = 'Lapatinib + HGF', dose,
    predicted_abundance = combination_abundance_pct, observed_abundance = 94, observed_sd = 5.7)]))
wilson <- rbindlist(list(old_wilson, wilson), use.names = TRUE)
fwrite(wilson, file.path(root, 'data/full_panel_wilson_source_comparison.tsv'), sep = '\t')
predicted <- wilson[, .(series = family, condition, abundance = predicted_abundance, sd = NA_real_)]
observed <- unique(wilson[, .(condition, abundance = observed_abundance, sd = observed_sd)])
observed[, series := 'Wilson observed']
plot_data <- rbindlist(list(predicted, observed), use.names = TRUE)
plot_data[, condition := factor(condition, levels = c('Lapatinib', 'Lapatinib + HGF'))]
plot_data[, series := factor(series, levels = c('Wilson observed', families))]
p <- ggplot(plot_data, aes(condition, abundance, color = series, group = series)) +
  geom_line(linewidth = .6) + geom_point(size = 2.3) +
  geom_errorbar(data = plot_data[series == 'Wilson observed'],
    aes(ymin = abundance - sd, ymax = abundance + sd), width = .08) +
  scale_color_manual(values = c('Wilson observed' = '#b32631', 'Tahoe only' = '#5574a6',
    'LINCS only' = '#a77917', 'Tahoe + LINCS' = '#198273'), name = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = 'bottom') +
  labs(x = NULL, y = 'Abundance (% untreated control)',
    title = 'AU565: lapatinib at 1 µM, with and without HGF',
    caption = 'Wilson: observed mean ± SD. Models: one frozen prediction per condition.\nHGF uses the mean transported molecular vector; Wilson outcomes were excluded from training.\nLINCS uses partial program observations; assay and HGF-dose transfer remain uncalibrated.')
ggsave(file.path(root, 'figures/wilson_1uM_comparison.pdf'), p, width = 6.2, height = 4.9)
ggsave(file.path(root, 'figures/wilson_1uM_comparison.png'), p, width = 6.2, height = 4.9, dpi = 160)
fwrite(plot_data, file.path(root, 'data/wilson_plot_source.tsv'), sep = '\t')

# The lapatinib figure contains every one of the 441 cells at every measured
# endpoint. Each source/dose panel has independent x and y scales.
single <- fread(file.path(root, 'data/held_single_predictions_frozen.tsv'))[
  drug == 'Lapatinib ditosylate']
old_single <- fread(file.path(reference, 'data/held_single_predictions_frozen.tsv'))[
  drug == 'Lapatinib ditosylate' & family == 'Linked molecular bridge']
old_single[, family := 'Tahoe only']
single <- rbindlist(list(old_single, single), use.names = TRUE)
single[, family := factor(family, levels = families)]
single[, dose_label := factor(dose, levels = c(.15625, .625, 2.5),
  labels = c('0.15625 µM', '0.625 µM', '2.5 µM'))]
metrics <- single[, .(cells = .N, r = cor(observed_log2, predicted_log2)), by = .(family, dose_label)]
stopifnot(all(metrics$cells == 441L))
metrics[, label := sprintf('r = %.3f; n = %d', r, cells)]
single[, panel := paste(family, dose_label, sep = ': ')]
metrics[, panel := paste(family, dose_label, sep = ': ')]
levels <- unlist(lapply(families, function(f) paste(f, levels(single$dose_label), sep = ': ')))
single[, panel := factor(panel, levels = levels)]
metrics[, panel := factor(panel, levels = levels)]
p <- ggplot(single, aes(observed_log2, predicted_log2)) +
  geom_point(aes(color = factor(fold)), size = .65, alpha = .5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = 'black') +
  geom_text(data = metrics, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.3, size = 2.8) + facet_wrap(~panel, ncol = 3, scales = 'free') +
  scale_color_manual(values = c('#2378a0', '#bd5a28', '#72419d'), name = 'Held-cell fold') +
  theme_bw(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom',
    strip.text = element_text(size = 8)) +
  labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
    title = 'Lapatinib: matched 441-cell comparison of molecular training sources')
ggsave(file.path(root, 'figures/lapatinib_all_sources_and_doses.pdf'), p, width = 9, height = 6.1)
ggsave(file.path(root, 'figures/lapatinib_all_sources_and_doses.png'), p, width = 9, height = 6.1, dpi = 150)
fwrite(metrics, file.path(root, 'data/lapatinib_source_comparison.tsv'), sep = '\t')

all_metrics <- fread(file.path(root, 'data/fixed_drug_dose_single_metrics.tsv'))
combination <- fread(file.path(root, 'data/fixed_pair_dose_summary.tsv'))
RNA <- fread(file.path(root, 'data/held_molecular_summary.tsv'))
fits <- fread(file.path(root, 'data/fit_audit.tsv'))
capture.output(cat('V84 assessment: all correlations are across cells receiving the same treatment and dose.\n',
  'Original ten drugs retained. Topotecan uses Tahoe; nine other drugs also receive available LINCS.\n',
  'All 441 lapatinib cells preserved at each endpoint; v82 is the original full-panel Tahoe comparator.\n',
  'RNA source comparison uses the same PRISM/ALMANAC functional data and shared F/J architecture.\n',
  'Eight source-training-RNA baseline directions predict 40 fixed original program coordinates.\n',
  'LINCS partial observations use a frozen loading-pattern approximation; inferred full scores are not measurements.\n',
  'One shared F and one shared symmetric J. Baseline occurs only in the molecular mapping.\n',
  'Wilson predictions at 1 µM were frozen before attaching observed values; HGF/Wilson never fitted.\n',
  'Previously revealed held outcomes remain retrospective. No selection or tuning on these outcomes.\n'),
  print(all_metrics[drug == 'Lapatinib ditosylate']), print(all_metrics[dose == 2.5]), print(combination),
  print(RNA[drug == 'Lapatinib ditosylate']), print(wilson),
  print(fits[, .(fits = .N, converged = sum(convergence == 0), iteration_limited = sum(convergence != 0)), by = family]),
  file = file.path(root, 'ASSESSMENT.txt'))
cat('Wilson and source-comparison figures saved.\n')

# Correct the full-panel overview captions and retain all ten drugs at all three
# measured endpoints. These are plot views of the same frozen v84 predictions.
full_prediction <- fread(file.path(root, 'data/held_single_predictions_frozen.tsv'))
overview <- list()
for (dose_value in c(.15625, .625, 2.5)) {
  rows <- full_prediction[abs(dose - dose_value) < 1e-9]
  labels <- all_metrics[abs(dose - dose_value) < 1e-9,
    .(drug, label = sprintf('r = %.3f; n = %d', pearson_log2, cells))]
  p <- ggplot(rows, aes(observed_log2, predicted_log2)) +
    geom_point(aes(color = factor(fold)), size = .6, alpha = .5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .45, color = 'black') +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.05, vjust = 1.3, size = 2.6) + facet_wrap(~drug, ncol = 5, scales = 'free') +
    scale_color_manual(values = c('#2378a0', '#bd5a28', '#72419d'), name = 'Held-cell fold') +
    theme_bw(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom',
      strip.text = element_text(size = 8)) +
    labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
      title = sprintf('Tahoe + LINCS at %g µM: original ten-drug panel', dose_value))
  overview[[as.character(dose_value)]] <- p
  ggsave(file.path(root, 'figures', paste0('Tahoe_LINCS_', gsub('\\.', '_', format(dose_value)), 'uM.png')),
    p, width = 12, height = 5.8, dpi = 140)
}
pdf(file.path(root, 'figures/all_measured_doses.pdf'), width = 12, height = 5.8)
for (p in overview) print(p)
dev.off()
