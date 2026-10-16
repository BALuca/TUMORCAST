#!/usr/bin/env Rscript
# Reuse68 comparison and complete-dose-curve views for frozen f followed by g.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(qs)
  library(digest)
})
root <- '/work/results/program_pilot/69_sequential_single_combination_learning'
setDTthreads(2L)
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom'))
save_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 140)
}
safe_cor <- function(x, y) {
  if (length(x) < 4L || sd(x) < 1e-10 || sd(y) < 1e-10) return(NA_real_)
  cor(x, y)
}
new <- fread(file.path(root, 'data/predictions.tsv.gz'))
stopifnot(uniqueN(new$fold) == 6L, new[split == 'held', .N, by = observation_id][, all(N == 1L)])
new[, method := 'Frozen f + learned g']
additive <- copy(new)
additive[, `:=`(method = 'Same frozen f; g = 0', predicted_log2 = additive_log2, predicted_g = 0)]
independent <- copy(new)
independent[, `:=`(method = 'Measured individual effects', predicted_log2 = log2(ratio_A * ratio_B), predicted_g = 0)]
z <- rbind(new, additive, independent)
means <- z[, .(actual_count = 100 * mean(ratio_AB), predicted_count = 100 * mean(2^predicted_log2),
  observed_g = mean(log2(ratio_AB / (ratio_A * ratio_B))), predicted_g = mean(predicted_g),
  residual_to_f = mean(log2(ratio_AB) - additive_log2),
  actual_log2 = mean(log2(ratio_AB)), predicted_log2 = mean(predicted_log2),
  low_count = 100 * min(ratio_AB), high_count = 100 * max(ratio_AB)),
  by = .(method, split, anchor, partner, nominal_A_uM, nominal_B_uM, condition_id)]
means[, `:=`(centered_actual = actual_count - mean(actual_count),
  centered_predicted = predicted_count - mean(predicted_count)), by = .(method, split, anchor, nominal_A_uM)]
summary <- means[, .(conditions = .N, combination_count_r = safe_cor(actual_count, predicted_count),
  combination_log2_r = safe_cor(actual_log2, predicted_log2), empirical_interaction_r = safe_cor(observed_g, predicted_g),
  residual_to_f_r = safe_cor(residual_to_f, predicted_g),
  within_anchor_dose_count_r = safe_cor(centered_actual, centered_predicted)), by = .(method, split)]
by_partner <- means[split == 'held', .(conditions = .N, combination_count_r = safe_cor(actual_count, predicted_count),
  empirical_interaction_r = safe_cor(observed_g, predicted_g), residual_to_f_r = safe_cor(residual_to_f, predicted_g)),
  by = .(method, partner)]
by_pair <- means[split == 'held', .(conditions = .N, combination_count_r = safe_cor(actual_count, predicted_count),
  empirical_interaction_r = safe_cor(observed_g, predicted_g)), by = .(method, anchor, partner)]
save_table(means, 'figure_condition_predictions')
save_table(summary, 'correlation_summary')
save_table(by_partner, 'held_correlations_by_partner')
save_table(by_pair, 'held_correlations_by_pair')
print(summary)
print(by_partner)

methods <- c('Measured individual effects', 'Same frozen f; g = 0', 'Frozen f + learned g')
colors <- setNames(c('#777777', '#009E73', '#D55E00'), methods)
d <- copy(means)
d[, `:=`(method = factor(method, levels = methods), split = factor(split, levels = c('training', 'held')))]
ann <- d[, .(label = sprintf('Pearson r = %.2f', safe_cor(actual_count, predicted_count))), by = .(method, split)]
p <- ggplot(d, aes(actual_count, predicted_count)) + geom_point(aes(color = anchor), size = 1.25, alpha = .6) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 3) +
  facet_wrap(vars(split, method), ncol = 3, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'Learn individual effects first, then combination modifications',
    subtitle = 'The same frozen single-drug function is used in both learned models',
    x = 'Measured cell count (% vehicle)', y = 'Predicted cell count (% vehicle)', color = 'First drug',
    caption = 'Held: every combination involving that partner excluded; individual experiments remain available.\nTraining points average fitted predictions across outer models. Previously inspected outcomes: retrospective development.')
save_plot(p, 'sequential_counts_training_vs_held', 10.5, 7.5)

# Separate g's actual training target from empirical departure from measured
# singles; disagreement between these two diagnoses imperfect single calibration.
d0 <- means[method == 'Frozen f + learned g']
d <- rbind(d0[, .(anchor, split, measured = residual_to_f, predicted = predicted_g,
    endpoint = 'Residual to learned f (training target)')],
  d0[, .(anchor, split, measured = observed_g, predicted = predicted_g,
    endpoint = 'Departure from measured individual effects')])
d[, split := factor(split, levels = c('training', 'held'))]
ann <- d[, .(label = sprintf('Pearson r = %.2f', safe_cor(measured, predicted))), by = .(endpoint, split)]
p <- ggplot(d, aes(measured, predicted)) + geom_point(aes(color = anchor), size = 1.4, alpha = .65) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2) +
  facet_wrap(vars(split, endpoint), ncol = 2, scales = 'free') + theme(aspect.ratio = 1, strip.text = element_text(size = 9)) +
  labs(title = 'What is the interaction component learning?',
    x = 'Measured log2 departure', y = 'Predicted g(a,b)', color = 'First drug')
save_plot(p, 'interaction_target_vs_empirical_departure', 9, 8)

d <- copy(means[split == 'held' & method != 'Measured individual effects'])
d[, method := factor(method, levels = methods[-1])]
ann <- d[, .(label = sprintf('r = %.2f', safe_cor(actual_count, predicted_count))), by = .(method, partner)]
p <- ggplot(d, aes(actual_count, predicted_count)) + geom_point(aes(color = anchor), size = 1.15, alpha = .7) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .45) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 2.8) +
  facet_wrap(vars(method, partner), ncol = 6, scales = 'free') +
  theme(aspect.ratio = 1, strip.text = element_text(size = 7)) +
  labs(title = 'Combination predictions for every excluded partner',
    x = 'Measured cell count (% vehicle)', y = 'Predicted cell count (% vehicle)', color = 'First drug')
save_plot(p, 'held_counts_each_partner', 14, 5.8)

d <- copy(means[split == 'held' & method == 'Frozen f + learned g'])
ann <- d[, .(label = sprintf('r = %.2f', safe_cor(observed_g, predicted_g))), by = partner]
zero_partners <- d[, .(is_zero = max(abs(predicted_g)) < 1e-10), by = partner][is_zero == TRUE, partner]
p <- ggplot(d, aes(observed_g, predicted_g)) + geom_point(aes(color = anchor), size = 1.4, alpha = .7) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2) +
  facet_wrap(vars(partner), ncol = 3, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'Actual interaction departures for every excluded partner',
    subtitle = if (length(zero_partners)) paste('Zero correction; r undefined:', paste(zero_partners, collapse = ', ')) else
      'Every partner excluded under both anchors and all doses',
    x = 'Measured log2[combination / individual-product]', y = 'Predicted g(a,b)', color = 'First drug')
save_plot(p, 'held_empirical_interactions_each_partner', 10.5, 7.5)
save_table(d[, .(measured_mean_departure = mean(observed_g), predicted_mean_g = mean(predicted_g),
  mean_residual_to_f = mean(residual_to_f)), by = partner], 'partner_interaction_offsets')

cairo_pdf(file.path(root, 'figures/all12_combination_dose_curves.pdf'), width = 10.4, height = 6.3)
for (aa in sort(unique(means$anchor))) for (bb in sort(unique(means$partner))) {
  d <- means[split == 'held' & anchor == aa & partner == bb]
  dots <- d[method == 'Measured individual effects']
  lines <- rbind(d[, .(nominal_A_uM, nominal_B_uM, method, value = predicted_count, endpoint = 'Cell count (% vehicle)')],
    d[, .(nominal_A_uM, nominal_B_uM, method, value = predicted_g, endpoint = 'Pair modification (log2)')])
  observed <- rbind(dots[, .(nominal_A_uM, nominal_B_uM, value = actual_count,
    low = low_count, high = high_count, endpoint = 'Cell count (% vehicle)')],
    z[split == 'held' & method == 'Measured individual effects' & anchor == aa & partner == bb,
      .(value = mean(log2(ratio_AB / (ratio_A * ratio_B))), low = min(log2(ratio_AB / (ratio_A * ratio_B))),
        high = max(log2(ratio_AB / (ratio_A * ratio_B))), endpoint = 'Pair modification (log2)'),
      by = .(nominal_A_uM, nominal_B_uM)])
  p <- ggplot(lines, aes(nominal_A_uM, value, color = method)) + geom_line(linewidth = .7) +
    geom_linerange(data = observed, aes(x = nominal_A_uM, ymin = low, ymax = high), inherit.aes = FALSE, color = 'grey55') +
    geom_point(data = observed, aes(nominal_A_uM, value), inherit.aes = FALSE, size = 1.6) +
    facet_grid(endpoint ~ nominal_B_uM, scales = 'free_y',
      labeller = labeller(nominal_B_uM = function(x) paste(bb, x, 'µM'))) +
    scale_x_log10(breaks = sort(unique(d$nominal_A_uM))) + scale_color_manual(values = colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 8)) +
    labs(title = paste(aa, '+', bb), subtitle = 'All outcomes for this partner excluded; black points/ranges are biological measurements.',
      x = paste(aa, 'dose (µM)'), y = NULL, color = NULL)
  print(p)
}
dev.off()

s <- fread(file.path(root, 'data/single_predictions.tsv.gz'))
sl <- rbind(s[, .(sample_id, source, drug, cell, dose, observed = log2(ratio), predicted = predicted_fitted_log2, evaluation = 'Fitted')],
  s[, .(sample_id, source, drug, cell, dose, observed = log2(ratio), predicted = predicted_cv_log2, evaluation = 'Single-data tuning CV')])
sl[, evaluation := factor(evaluation, levels = c('Fitted', 'Single-data tuning CV'))]
ann <- sl[, .(n = .N, r = safe_cor(observed, predicted)), by = .(source, evaluation)]
save_table(ann, 'single_correlations')
save_table(sl[, .(n = .N, r = safe_cor(observed, predicted)), by = .(source, evaluation, drug)], 'single_correlations_by_drug')
ann[, label := sprintf('r = %.2f; n = %s', r, format(n, big.mark = ','))]
p <- ggplot(sl, aes(observed, predicted)) + geom_point(size = .65, alpha = .2, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 2.9) +
  facet_wrap(vars(evaluation, source), ncol = 3, scales = 'free') + theme(aspect.ratio = 1, strip.text = element_text(size = 8)) +
  labs(title = 'Single-drug function learned before seeing combinations',
    subtitle = 'CV excludes whole Tahoe cells and whole BT20 drugs; these folds select f and are not independent validation',
    x = 'Measured log2 relative abundance', y = 'Predicted individual effect f(h,x)')
save_plot(p, 'single_drug_fitted_and_cv', 11, 7.8)

# Calibration of the exact eight drugs entering the combination study.
s <- s[source == 'Same-screen single controls']
sm <- s[, .(observed = 100 * mean(ratio), fitted = 100 * mean(2^predicted_fitted_log2),
  cv = 100 * mean(2^predicted_cv_log2), low = 100 * min(ratio), high = 100 * max(ratio)), by = .(drug, dose)]
save_table(sm, 'combination_constituent_single_curves')
lines <- melt(sm, id.vars = c('drug', 'dose', 'observed', 'low', 'high'), measure.vars = c('fitted', 'cv'),
  variable.name = 'prediction', value.name = 'value')
p <- ggplot(lines, aes(dose, value, color = prediction)) + geom_line(linewidth = .7) +
  geom_linerange(data = sm, aes(x = dose, ymin = low, ymax = high), inherit.aes = FALSE, color = 'grey55') +
  geom_point(data = sm, aes(dose, observed), inherit.aes = FALSE, size = 1.5) +
  facet_wrap(vars(drug), ncol = 4, scales = 'free_x') + scale_x_log10() +
  labs(title = 'Individual effects of the drugs used in combinations',
    subtitle = 'Black: measured biological mean/range. Green: frozen f used for combination prediction.',
    x = 'Dose (µM)', y = 'Cell count (% vehicle)', color = NULL) +
  scale_color_manual(values = c(fitted = '#009E73', cv = '#0072B2'),
    labels = c(fitted = 'Frozen f (singles available)', cv = 'Single-data tuning CV'))
save_plot(p, 'constituent_single_drug_dose_curves', 10.5, 6)

writeLines(c('Version69: learn f from singles, freeze it, then learn g from combinations.',
  'Control is the exact same frozen f with g=0. Measured individual independence is an additional reference.',
  'Tahoe uses all8842 original40 genes, not the L1000 subset. BT20 remains a partial L1000 projection.',
  'f: baseline-conditioned Gaussian response kernel, Nyström approximation, shared function plus BT20 adjustment.',
  'Baseline original40 scores use matched CCLE expression. Source-specific response RMS preserves biological zero.',
  'f settings use only single-data CV: whole Tahoe cells and BT20 drugs. Single CV reported is tuning CV.',
  'g:3321 symmetric products of program/saturating-program/frozen-f features; elastic net with inner partner exclusions.',
  'In-place correction: sparsest candidate within one SE of best grouped CV objective, including explicit g=0.',
  'f reused unchanged after upstream hashes verified. Combination outcomes cannot update the frozen function.',
  'Six whole-partner combination exclusions; individual measurements for those partners are legitimately available.',
  'g learns residual to predicted f. This differs from empirical interaction when f is inaccurate; both are shown.',
  'All molecular responses are measured, not baseline-predicted. No held-cell combination transfer or combination RNA input.',
  'Assay-specific adjustment does not establish cross-assay biological equivalence; missing BT20 genes are not imputed.',
  'Previously inspected outcomes remain retrospective. PRISM dose interpolation retained from the established panel.',
  '', capture.output(print(summary)), '', capture.output(print(ann))), file.path(root, 'REPORT.txt'))
save_table(data.table(path = '/work/src/program_pilot/69_sequential_learning_figures.R',
  sha256 = digest(file = '/work/src/program_pilot/69_sequential_learning_figures.R', algo = 'sha256')), 'figure_script_manifest')
fwrite(data.table(path = list.files(file.path(root, 'figures'), full.names = TRUE)), file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
writeLines('f, all six g models, and figures complete.', file.path(root, 'STATUS.txt'))
