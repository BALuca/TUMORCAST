#!/usr/bin/env Rscript
# Matched comparisons for joint learning: all conditions, all excluded partners.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(qs)
  library(digest)
})
root <- '/work/results/program_pilot/68_joint_tahoe_prism_combination_learning'
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
z <- fread(file.path(root, 'data/predictions.tsv.gz'))
stopifnot(uniqueN(z$fold) == 6L, uniqueN(z$model) == 3L,
  z[split == 'held', .N, by = .(model, observation_id)][, all(N == 1L)])
names_map <- c(Joint_Tahoe_PRISM = 'Joint + Tahoe/PRISM', Additive_Tahoe_PRISM = 'Additive + Tahoe/PRISM',
  Joint_BT20_only = 'Joint: BT20 only')
z[, method := unname(names_map[model])]
independent <- copy(z[model == 'Joint_Tahoe_PRISM'])
independent[, `:=`(method = 'Measured individual effects', predicted_log2 = log2(ratio_A * ratio_B), predicted_g = 0)]
z <- rbind(z, independent)
means <- z[, .(actual_count = 100 * mean(ratio_AB), predicted_count = 100 * mean(2^predicted_log2),
  observed_g = mean(log2(ratio_AB / (ratio_A * ratio_B))), predicted_g = mean(predicted_g),
  actual_log2 = mean(log2(ratio_AB)), predicted_log2 = mean(predicted_log2),
  low_count = 100 * min(ratio_AB), high_count = 100 * max(ratio_AB)),
  by = .(method, split, anchor, partner, nominal_A_uM, nominal_B_uM, condition_id)]
means[, `:=`(centered_actual = actual_count - mean(actual_count),
  centered_predicted = predicted_count - mean(predicted_count)), by = .(method, split, anchor, nominal_A_uM)]
summary <- means[, .(conditions = .N, combination_count_r = safe_cor(actual_count, predicted_count),
  combination_log2_r = safe_cor(actual_log2, predicted_log2), interaction_r = safe_cor(observed_g, predicted_g),
  within_anchor_dose_count_r = safe_cor(centered_actual, centered_predicted)), by = .(method, split)]
by_partner <- means[split == 'held', .(conditions = .N, combination_count_r = safe_cor(actual_count, predicted_count),
  interaction_r = safe_cor(observed_g, predicted_g)), by = .(method, partner)]
by_pair <- means[split == 'held', .(conditions = .N, combination_count_r = safe_cor(actual_count, predicted_count),
  interaction_r = safe_cor(observed_g, predicted_g)), by = .(method, anchor, partner)]
save_table(means, 'figure_condition_predictions')
save_table(summary, 'correlation_summary')
save_table(by_partner, 'held_correlations_by_partner')
save_table(by_pair, 'held_correlations_by_pair')
print(summary)
print(by_partner)

methods <- c('Measured individual effects', 'Additive + Tahoe/PRISM', 'Joint + Tahoe/PRISM', 'Joint: BT20 only')
colors <- setNames(c('#777777', '#009E73', '#D55E00', '#0072B2'), methods)
d <- copy(means)
d[, `:=`(method = factor(method, levels = methods), split = factor(split, levels = c('training', 'held')))]
ann <- d[, .(label = sprintf('Pearson r = %.2f', safe_cor(actual_count, predicted_count))), by = .(method, split)]
p <- ggplot(d, aes(actual_count, predicted_count)) + geom_point(aes(color = anchor), size = 1.25, alpha = .6) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 3) +
  facet_wrap(vars(split, method), ncol = 4, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'Jointly learned individual and combination effects',
    subtitle = 'Same combinations and exclusions; additive control has no pair-interaction terms',
    x = 'Measured cell count (% vehicle)', y = 'Predicted cell count (% vehicle)', color = 'First drug',
    caption = 'Held: every combination involving that partner excluded. Its single-drug experiments remain available.\nTraining points average fitted predictions across outer models. Retrospective development, not independent confirmation.')
save_plot(p, 'joint_counts_training_vs_held', 13, 7.6)

d <- copy(means[method %chin% c('Joint + Tahoe/PRISM', 'Joint: BT20 only')])
d[, split := factor(split, levels = c('training', 'held'))]
ann <- d[, .(label = sprintf('Pearson r = %.2f', safe_cor(observed_g, predicted_g))), by = .(method, split)]
p <- ggplot(d, aes(observed_g, predicted_g)) + geom_point(aes(color = anchor), size = 1.4, alpha = .65) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2) +
  facet_wrap(vars(split, method), ncol = 2, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'Does the learned interaction recover departure from individual effects?',
    subtitle = 'f and g learned jointly; measured singles used only to calculate the observed departure here',
    x = 'Measured log2[combination / individual-product]', y = 'Learned pair component g(a,b)', color = 'First drug')
save_plot(p, 'joint_interactions_training_vs_held', 8.3, 8)

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
save_plot(p, 'held_counts_each_partner', 14, 8)

# All 12 combinations are displayed, including biological repeat ranges.
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

# Single-effect fitting is shown explicitly; these are training observations.
s <- fread(file.path(root, 'data/single_predictions.tsv.gz'))
s[, method := unname(names_map[model])]
sm <- s[, .(observed = mean(log2(ratio)), predicted = mean(predicted_log2)),
  by = .(method, source, sample_id, drug, nominal_dose_uM)]
ann <- sm[, .(n = .N, r = safe_cor(observed, predicted)), by = .(method, source)]
save_table(ann, 'single_fitted_correlations')
save_table(sm, 'single_fitted_predictions')
ann[, label := sprintf('Fitted r = %.2f; n = %s', r, format(n, big.mark = ','))]
p <- ggplot(sm, aes(observed, predicted)) + geom_point(size = .65, alpha = .2, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 2.7) +
  facet_wrap(vars(source, method), ncol = 3, scales = 'free') + theme(aspect.ratio = 1, strip.text = element_text(size = 8)) +
  labs(title = 'Single-drug effects learned inside the joint model', subtitle = 'Training fits, not held-out single-drug predictions',
    x = 'Measured log2 relative abundance', y = 'Predicted individual effect f(x)')
save_plot(p, 'single_drug_fitting_by_source', 11, 10)

writeLines(c('Version68: simultaneous single and combination functional learning.',
  'Primary comparison: Joint + Tahoe/PRISM versus Additive + Tahoe/PRISM.',
  'Both learn singles and combinations together. Only the full model has g(a,b).',
  'Third model omits external Tahoe/PRISM singles, retaining BT20 singles and combinations.',
  'Original40 frozen weights; common measured genes; source-specific outcome-free RMS scaling;',
  'shared single coefficients plus penalized BT20 deviations. f(0)=0; symmetric g vanishes if either input is zero.',
  'PRISM labels are interpolated within measured dose brackets, without extrapolation, from the existing panel.',
  'Equal source totals prevent PRISM row count from dominating. All single experiments are permitted, including excluded partners.',
  'All combinations containing each predicted partner excluded under both anchors; tuning repeats partner exclusions internally.',
  'Molecular input is measured response. Baseline-to-response learning and cross-cell combination transfer are not tested.',
  'BT20 RNA:3h; function:72h imaged-cell count. PRISM and Tahoe use their own assay/exposure measurements.',
  'Joint fitting does not ensure learned g equals empirical interaction if f is inaccurate; both are plotted.',
  'Previously inspected outcomes: retrospective development. Partial program coverage and assay/time differences remain limitations.',
  '', capture.output(print(summary)), '', capture.output(print(by_partner))), file.path(root, 'REPORT.txt'))
save_table(data.table(path = '/work/src/program_pilot/68_joint_learning_figures.R',
  sha256 = digest(file = '/work/src/program_pilot/68_joint_learning_figures.R', algo = 'sha256')), 'figure_script_manifest')
fwrite(data.table(path = list.files(file.path(root, 'figures'), full.names = TRUE)),
  file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
writeLines('Joint fits and figures complete.', file.path(root, 'STATUS.txt'))
