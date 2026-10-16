#!/usr/bin/env Rscript
# Display measured departures from the declared independent-single reference.
# This is an observed-data view of workflow61; no fitting or model selection.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
data_dir <- file.path(root, 'data')
figure_dir <- file.path(root, 'figures')
obs <- fread(file.path(data_dir, 'matched_molecular_links.tsv.gz'))[
  time_h == 3 & molecular_supported == TRUE]
stopifnot(nrow(obs) == 496L, !anyDuplicated(obs$observation_id))
conditions <- obs[, .(expected_pct = 100 * mean(ratio_A * ratio_B),
  observed_pct = 100 * mean(ratio_AB), extra_pp = mean(excess_inhibition),
  biological_repeats = uniqueN(biological_replicate), positive_every_repeat = all(excess_inhibition > 0)),
  by = .(anchor, partner, nominal_A_uM, nominal_B_uM)]
stopifnot(nrow(conditions) == 198L,
  max(abs(conditions$extra_pp - (conditions$expected_pct - conditions$observed_pct))) < 1e-8)
fwrite(conditions, file.path(data_dir, 'figure_observed_single_expectation_departures.tsv'), sep = '\t')
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom'))
p <- ggplot(conditions, aes(expected_pct, observed_pct)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey40', linewidth = .7) +
  geom_point(aes(color = anchor), size = 2, alpha = .75) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .65) +
  scale_color_manual(values = c(alpelisib = '#0072B2', dactolisib = '#D55E00')) +
  annotate('text', x = -Inf, y = Inf, hjust = -.03, vjust = 1.15, size = 3.8,
    label = sprintf('Pearson r = %.3f; n = 198\n128 conditions: >=10 pp below expectation\nand same direction in every repeat',
      cor(conditions$expected_pct, conditions$observed_pct))) +
  labs(title = 'Measured combinations depart from the individual-drug expectation',
    subtitle = 'Below the dashed equality line = fewer cells than expected',
    x = 'Independent-single expectation (% of vehicle)', y = 'Measured combination cell count (% of vehicle)',
    color = 'First drug', caption = 'No prediction model. Dashed line: equality; solid line: fitted trend. Each point averages 2–3 biological experiments.') +
  theme(aspect.ratio = 1)
ggsave(file.path(figure_dir, 'observed_combination_departures.pdf'), p, width = 7, height = 6.6, device = cairo_pdf)
ggsave(file.path(figure_dir, 'observed_combination_departures.png'), p, width = 7, height = 6.6, dpi = 160)

# Illustrative named case; the complete dataset remains in the preceding plot.
case <- obs[anchor == 'alpelisib' & partner == 'trametinib' & nominal_A_uM == 1.11 & nominal_B_uM == .37]
stopifnot(nrow(case) == 3L, uniqueN(case$biological_replicate) == 3L)
z <- rbindlist(list(
  case[, .(biological_replicate, group = 'Alpelisib\nalone', cell_count_pct = 100 * ratio_A)],
  case[, .(biological_replicate, group = 'Trametinib\nalone', cell_count_pct = 100 * ratio_B)],
  case[, .(biological_replicate, group = 'Expected\ncombination', cell_count_pct = 100 * ratio_A * ratio_B)],
  case[, .(biological_replicate, group = 'Measured\ncombination', cell_count_pct = 100 * ratio_AB)]))
z[, group := factor(group, levels = c('Alpelisib\nalone', 'Trametinib\nalone',
  'Expected\ncombination', 'Measured\ncombination'))]
means <- z[, .(mean_pct = mean(cell_count_pct), label_y = max(cell_count_pct) + 5), by = group]
fwrite(z, file.path(data_dir, 'figure_observed_alpelisib_trametinib_example.tsv'), sep = '\t')
q <- ggplot(means, aes(group, mean_pct)) +
  geom_col(aes(fill = group), width = .65) +
  geom_point(data = z, aes(group, cell_count_pct), inherit.aes = FALSE, size = 2,
    position = position_jitter(width = .045, height = 0, seed = 61)) +
  geom_text(aes(y = label_y, label = sprintf('%.1f%%', mean_pct)), vjust = 0, size = 4) +
  scale_fill_manual(values = c('grey80', 'grey80', 'grey55', '#0072B2'), guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0, .14))) +
  labs(title = 'Actual example: the combination leaves fewer cells than expected',
    subtitle = 'Alpelisib 1.11 µM + trametinib 0.37 µM; BT20, 72 h',
    x = NULL, y = 'Cell count (% of vehicle)',
    caption = 'Bars: biological means. Dots: three biological experiments. Expectation is calculated within each experiment, then averaged.')
ggsave(file.path(figure_dir, 'observed_combination_example.pdf'), q, width = 7.7, height = 4.5, device = cairo_pdf)
ggsave(file.path(figure_dir, 'observed_combination_example.png'), q, width = 7.7, height = 4.5, dpi = 160)
message('Measured departures displayed without fitting any model.')
