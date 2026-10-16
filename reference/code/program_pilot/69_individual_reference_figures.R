#!/usr/bin/env Rscript
# Descriptive views of measured individual effects, not another model fit.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
root <- '/work/results/program_pilot/69_sequential_single_combination_learning'
# Stream gzip directly, avoiding temporary decompression on the full disk.
z <- as.data.table(read.delim(gzfile(file.path(root, 'data/predictions.tsv.gz')), check.names = FALSE))[split == 'held']
stopifnot(nrow(z) == 496L, uniqueN(z$condition_id) == 198L)
m <- z[, .(observed = 100 * mean(ratio_AB), expected = 100 * mean(ratio_A * ratio_B)),
  by = .(condition_id, anchor, partner, nominal_A_uM, nominal_B_uM)]
m[, `:=`(observed_centered = observed - mean(observed), expected_centered = expected - mean(expected)),
  by = .(anchor, nominal_A_uM, nominal_B_uM)]
r_raw <- cor(m$expected, m$observed)
r_centered <- cor(m$expected_centered, m$observed_centered)
stopifnot(abs(r_centered - .7687249) < 1e-6)
m[, group := paste(anchor, nominal_A_uM, '+ partner', nominal_B_uM, 'uM')]
fwrite(m, file.path(root, 'data/individual_reference_scatter_values.tsv'), sep = '\t')
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom'))
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 140)
}
colors <- c(alpelisib = '#0072B2', dactolisib = '#D55E00')
left <- ggplot(m, aes(expected, observed, color = anchor)) + geom_point(size = 1.5, alpha = .7) +
  geom_abline(slope = 1, intercept = 0, color = 'grey60', linetype = 2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .6) +
  scale_color_manual(values = colors) + theme(aspect.ratio = 1) +
  labs(title = sprintf('Actual cell counts: Pearson r = %.3f', r_raw),
    subtitle = 'Dashed line: exact agreement with independence',
    x = 'Independent expectation (% vehicle)', y = 'Observed combination (% vehicle)', color = 'First drug')
right <- ggplot(m, aes(expected_centered, observed_centered, color = anchor)) +
  geom_hline(yintercept = 0, color = 'grey85') + geom_vline(xintercept = 0, color = 'grey85') +
  geom_point(size = 1.5, alpha = .7) +
  geom_abline(slope = 1, intercept = 0, color = 'grey60', linetype = 2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .6) +
  scale_color_manual(values = colors) + theme(aspect.ratio = 1) +
  labs(title = sprintf('Within-dose differences: Pearson r = %.3f', r_centered),
    subtitle = 'Both doses and first-drug identity fixed; compare partners',
    x = 'Expected minus dose-group mean (percentage points)',
    y = 'Observed minus dose-group mean (percentage points)', color = 'First drug')
p <- (left | right) + plot_layout(guides = 'collect') + plot_annotation(
  title = 'What the individual-effects correlation actually measures',
  caption = paste('198 combination conditions in BT20; each point averages biological repeats. No fitted combination model is used.',
    'Right: expected and observed means are subtracted separately within each of45 fixed-dose groups.',
    'Observed group means use combination measurements for descriptive evaluation. Correlation is not percent accuracy.', sep = '\n'))
save_plot(p, 'individual_reference_raw_and_within_dose', 11.5, 6.3)

# Fixed example chosen before reading its outcomes: both drugs at1.11uM.
example <- m[anchor == 'alpelisib' & nominal_A_uM == 1.11 & nominal_B_uM == 1.11]
stopifnot(nrow(example) == 6L)
example[, partner := factor(partner, levels = rev(partner))]
long <- melt(example, id.vars = 'partner', measure.vars = c('expected', 'observed'),
  variable.name = 'quantity', value.name = 'count')
p <- ggplot(example, aes(y = partner)) +
  geom_segment(aes(x = observed, xend = expected, yend = partner), color = 'grey70', linewidth = 1) +
  geom_point(data = long, aes(x = count, color = quantity), size = 3) +
  geom_text(data = long, aes(x = count, label = sprintf('%.1f', count)), nudge_y = .19, size = 3.3) +
  scale_color_manual(values = c(expected = '#777777', observed = '#0072B2'),
    labels = c(expected = 'Individual-effect expectation', observed = 'Measured combination')) +
  scale_x_continuous(expand = expansion(mult = c(.05, .08))) +
  labs(title = 'Large departures remain at exactly the same doses',
    subtitle = 'BT20: alpelisib1.11uM plus each partner1.11uM',
    x = 'Cell count (% vehicle)', y = NULL, color = NULL,
    caption = 'Expectation = product of measured single-drug relative counts. No fitted combination model.')
save_plot(p, 'fixed_dose_individual_expectation_vs_observed', 8.5, 4.8)

groups <- unique(m$group)
cairo_pdf(file.path(root, 'figures/all_fixed_dose_individual_reference.pdf'), width = 10, height = 10)
for (start in seq(1, length(groups), by = 9L)) {
  d <- m[group %chin% groups[start:min(start + 8L, length(groups))]]
  ann <- d[, .(label = sprintf('n=%d; r=%.2f', .N, cor(expected, observed))), by = group]
  p <- ggplot(d, aes(expected, observed, color = partner)) + geom_point(size = 1.8) +
    geom_abline(slope = 1, intercept = 0, color = 'grey65', linetype = 2) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .4) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.04, vjust = 1.2, size = 2.8) +
    facet_wrap(vars(group), ncol = 3, scales = 'free') + theme(aspect.ratio = 1, strip.text = element_text(size = 8)) +
    labs(title = 'Every fixed-dose comparison: measured singles versus combination',
      x = 'Independent expectation (% vehicle)', y = 'Observed combination (% vehicle)', color = 'Partner',
      caption = 'Both doses fixed within each panel; BT20 only. Small-panel correlations are descriptive.')
  print(p)
}
dev.off()
cat('Raw r:', r_raw, '; within-dose r:', r_centered, '\n')
