#!/usr/bin/env Rscript
# Describe the completed functional-only predictions; no fitting or selection.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
setDTthreads(2L)
d <- fread(file.path(root, 'data/predictions.tsv.gz'))[
  time_h == 3 & split == 'held_partner' & model == 'Singles_SVR']
stopifnot(nrow(d) == 496L, !anyDuplicated(d$observation_id))
z <- d[, .(observed = 100 * mean(relative_AB), single_A = 100 * mean(relative_A),
  single_B = 100 * mean(relative_B), product = 100 * mean(relative_A * relative_B),
  prediction = 100 * mean(predicted_relative_AB),
  measured_extra = mean(observed_extra_pp), predicted_extra = mean(predicted_extra_pp)),
  by = .(anchor, partner, nominal_A_uM, nominal_B_uM)]
stopifnot(nrow(z) == 198L)
summary <- z[, .(n = .N, r_single_A = cor(single_A, observed),
  r_single_B = cor(single_B, observed), r_product = cor(product, observed),
  r_prediction = cor(prediction, observed), r_extra_effect = cor(measured_extra, predicted_extra),
  single_A_min_pct = min(single_A), single_A_max_pct = max(single_A),
  observed_min_pct = min(observed), observed_max_pct = max(observed),
  observed_sd_pp = sd(observed)), by = anchor]
overall <- z[, .(anchor = 'Pooled', n = .N, r_single_A = cor(single_A, observed),
  r_single_B = cor(single_B, observed), r_product = cor(product, observed),
  r_prediction = cor(prediction, observed), r_extra_effect = cor(measured_extra, predicted_extra),
  single_A_min_pct = min(single_A), single_A_max_pct = max(single_A),
  observed_min_pct = min(observed), observed_max_pct = max(observed), observed_sd_pp = sd(observed))]
fwrite(rbind(summary, overall), file.path(root, 'data/no_RNA_prediction_drivers_summary.tsv'), sep = '\t')
fwrite(z, file.path(root, 'data/no_RNA_prediction_drivers_conditions.tsv'), sep = '\t')

# Remove pair means to distinguish between-pair variation from dose variation.
centered <- copy(z)
centered[, c('observed', 'product', 'prediction') := lapply(.SD, function(x) x - mean(x)),
  by = .(anchor, partner), .SDcols = c('observed', 'product', 'prediction')]
variance_summary <- data.table(
  fraction_observed_variation_between_pair_means = 1 - sum(centered$observed^2) / sum((z$observed - mean(z$observed))^2),
  within_pair_product_r = cor(centered$product, centered$observed),
  within_pair_prediction_r = cor(centered$prediction, centered$observed))
fwrite(variance_summary, file.path(root, 'data/no_RNA_prediction_drivers_variance.tsv'), sep = '\t')
pair_summary <- z[, .(n = .N, r_product = cor(product, observed),
  r_prediction = cor(prediction, observed)), by = .(anchor, partner)]
fwrite(pair_summary, file.path(root, 'data/no_RNA_prediction_drivers_pairs.tsv'), sep = '\t')

long <- melt(z, id.vars = c('anchor', 'partner', 'nominal_A_uM', 'nominal_B_uM', 'observed'),
  measure.vars = c('single_A', 'product', 'prediction'), variable.name = 'input', value.name = 'value')
long[, input := factor(input, levels = c('single_A', 'product', 'prediction'),
  labels = c('Measured first drug alone', 'Measured singles multiplied\nNo model trained', 'Functional-only SVR\nNo RNA supplied'))]
long[, anchor := factor(anchor, levels = c('alpelisib', 'dactolisib'))]
stats <- long[, .(label = sprintf('r = %.3f; n = %d', cor(value, observed), .N)), by = .(anchor, input)]
p <- ggplot(long, aes(value, observed)) +
  geom_point(aes(color = nominal_A_uM), size = 1.6, alpha = .75) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .55) +
  geom_text(data = stats, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.3, size = 3.2) +
  facet_wrap(vars(anchor, input), ncol = 3, scales = 'free') +
  scale_color_viridis_c(trans = 'log10', breaks = c(.04, .37, 10)) +
  theme_bw(base_size = 10) + theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
    legend.position = 'bottom', strip.text = element_text(size = 9), plot.title.position = 'plot') +
  labs(title = 'The high pooled correlation is driven by known single-drug functional dose responses',
    subtitle = 'BT20 only; combinations containing each partner excluded from training; 198 dose conditions',
    x = 'Supplied single-drug count or predicted combination count (% of vehicle)',
    y = 'Measured combination cell count (% of vehicle)', color = 'First-drug dose (µM)',
    caption = paste('Each point averages 2–3 biological experiments. Individual function and combined function are measured at 72 h.',
      '\nTop: alpelisib combinations. Bottom: dactolisib combinations. Axes scale independently; lines are fitted trends.'))
ggsave(file.path(root, 'figures/no_RNA_prediction_drivers.pdf'), p, width = 10.4, height = 7.5, device = cairo_pdf)
ggsave(file.path(root, 'figures/no_RNA_prediction_drivers.png'), p, width = 10.4, height = 7.5, dpi = 160)
print(rbind(summary, overall))
print(variance_summary)
