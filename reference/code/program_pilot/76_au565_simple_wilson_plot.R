#!/usr/bin/env Rscript
# A readable view of saved predictions: one average per treatment condition.
# The three training fits are retained in the source table, not plotted as cells.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
root <- '/work/results/program_pilot/76_variable_context_and_wilson_comparison'
prediction_path <- file.path(root, 'data/wilson_dose_absolute_predictions.tsv')
observed_path <- file.path(root, 'data/wilson_observed_reference.tsv')
stopifnot(file.exists(prediction_path), file.exists(observed_path))
prediction <- fread(prediction_path)
observed <- fread(observed_path)
stopifnot(nrow(prediction) == 6L, nrow(observed) == 2L,
  !anyDuplicated(prediction[, .(original_fold, condition)]),
  all(is.finite(prediction$abundance)), all(is.finite(observed$abundance)))

# Arithmetic mean of the three saved abundance predictions. These are different
# training splits for the same AU565 cell line, not three biological replicates.
average <- prediction[, .(abundance = mean(abundance), training_fits = .N), by = condition]
stopifnot(all(average$training_fits == 3L))
average[, `:=`(series = 'Model prediction', sd = NA_real_)]
observed[, series := 'Wilson measurement']
plot_data <- rbind(average[, .(condition, series, abundance, sd)],
  observed[, .(condition, series, abundance, sd)])
plot_data[, condition := factor(condition, levels = c('Lapatinib', 'Lapatinib + HGF'))]
plot_data[, series := factor(series, levels = c('Model prediction', 'Wilson measurement'))]
fwrite(plot_data, file.path(root, 'data/au565_simple_wilson_plot.tsv'), sep = '\t')

# Both series use untreated control as 100%. Only experimental SD is displayed;
# differences among training fits are not an experimental uncertainty interval.
p <- ggplot(plot_data, aes(condition, abundance, fill = series)) +
  geom_hline(yintercept = 100, color = 'grey65', linetype = 2, linewidth = .35) +
  geom_col(position = position_dodge(width = .72), width = .62) +
  geom_errorbar(aes(ymin = abundance - sd, ymax = abundance + sd),
    position = position_dodge(width = .72), width = .13, na.rm = TRUE, linewidth = .55) +
  geom_text(aes(label = sprintf('%.1f%%', abundance),
    y = abundance + fifelse(is.na(sd), 0, sd) + 4),
    position = position_dodge(width = .72), size = 3.5) +
  scale_fill_manual(values = c('Model prediction' = '#2166AC', 'Wilson measurement' = '#B2182B')) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, .01))) +
  scale_x_discrete(labels = c('Lapatinib', 'Lapatinib\n+ HGF')) +
  theme_classic(base_size = 11) +
  theme(legend.position = 'bottom', plot.caption = element_text(hjust = 0, size = 8),
    plot.title = element_text(size = 12)) +
  labs(x = NULL, y = 'Cell abundance (% untreated control)', fill = NULL,
    title = 'AU565 treated with 1 µM lapatinib',
    caption = '100% = untreated control. Lower bars mean less cell abundance.\nBlue: average of three saved model fits. Red: Wilson’s mean ± SD.\nThis compares absolute abundance; the version75 HGF plot used HGF-only control.\nHGF-alone prediction and transfer across assays remain unvalidated.')
ggsave(file.path(root, 'figures/au565_simple_wilson_comparison.pdf'), p,
  width = 5.25, height = 4.8, device = cairo_pdf)
ggsave(file.path(root, 'figures/au565_simple_wilson_comparison.png'), p,
  width = 5.25, height = 4.8, dpi = 170)

# Confirm that simplifying the display preserves all original predictions.
verification <- merge(average, prediction[, .(replayed_average = mean(abundance)), by = condition], by = 'condition')
stopifnot(max(abs(verification$abundance - verification$replayed_average)) == 0)
manifest <- data.table(path = c(prediction_path, observed_path,
  '/work/src/program_pilot/76_au565_simple_wilson_plot.R'))
manifest[, sha256 := vapply(path, digest, character(1), file = TRUE, algo = 'sha256')]
fwrite(manifest, file.path(root, 'data/au565_simple_wilson_plot_sources.tsv'), sep = '\t')
print(plot_data)
