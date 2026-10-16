#!/usr/bin/env Rscript
# A view of frozen v104 predictions; no model refitting or selection by outcomes.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
root <- '/work/results/program_pilot/104_erbb2_met_lapatinib_shared_decoder/linear'
measured <- fread(file.path(root, 'data/measured_predictions_frozen.tsv'))
predictions <- fread(file.path(root, 'data/baseline_predictions_frozen.tsv'))
matched_ids <- measured[evaluation == 'held_cell' & drug == 'Lapatinib ditosylate', depmap_id]
plot_data <- predictions[evaluation == 'held_cell' & drug == 'Lapatinib ditosylate' & depmap_id %in% matched_ids]
stopifnot(nrow(plot_data) == 30L, !anyDuplicated(plot_data$depmap_id), all(plot_data$dose == 5))
plot_data[, measured_sensitivity := -observed_log2]
plot_data[, predicted_sensitivity := -predicted_log2]
r <- cor(plot_data$measured_sensitivity, plot_data$predicted_sensitivity)
stopifnot(abs(r - 0.5894183) < 1e-6)
fwrite(plot_data, file.path(root, 'data/lapatinib_matched_receptor_bridge_scatter.tsv'), sep = '\t')
plot <- ggplot(plot_data, aes(measured_sensitivity, predicted_sensitivity)) +
  geom_point(color = '#0072B2', size = 2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey30', linewidth = .6) +
  annotate('text', x = -Inf, y = Inf, hjust = -.1, vjust = 1.2,
           label = sprintf('Pearson r = %.3f\nn = 30 held-out cells', r), size = 3.4) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(), plot.title = element_text(size = 11)) +
  labs(title = 'v104 | Lapatinib 5 µM: matched 30 cells',
       x = 'Measured sensitivity (−log2 abundance)',
       y = 'Predicted sensitivity (−log2 abundance)',
       caption = 'ERBB2/HER2 + MET → predicted displacements → shared decoder.\nWhole-cell held-out predictions; previously revealed outcomes are retrospective.\nSubset defined by molecular data availability, not prediction performance.')
ggsave(file.path(root, 'figures/lapatinib_matched_receptor_bridge.pdf'), plot,
       width = 4.6, height = 5.1, device = cairo_pdf)
ggsave(file.path(root, 'figures/lapatinib_matched_receptor_bridge.png'), plot,
       width = 4.6, height = 5.1, dpi = 180)
cat(sprintf('Saved matched 30-cell scatter: r = %.7f\n', r))
