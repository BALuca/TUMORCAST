#!/usr/bin/env Rscript
# Interpret the prespecified combination exclusion using assay-scale effects.
# Every frozen representation is shown; held outcomes do not select a winner.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
setDTthreads(2L)
root <- '/work/results/program_pilot/60_defensible_response_composition'
path <- file.path(root, 'data/transfer_predictions.tsv.gz')
stopifnot(file.exists(path))
dd <- fread(path)
required <- c('observation_id', 'cell', 'pair', 'drug_A', 'dose_A_uM', 'drug_B', 'dose_B_uM',
  'curve_id', 'plate_id', 'source_release', 'viability_A', 'viability_B', 'viability_AB',
  'observed_extra_pp', 'bliss_viability', 'model', 'predicted_extra_pp', 'predicted_viability', 'split', 'focal')
stopifnot(all(required %in% names(dd)), nrow(dd) > 0L,
  all(is.finite(dd$observed_extra_pp)), all(is.finite(dd$predicted_extra_pp)))
for (folder in c('data', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0('transfer_plot_', name, '.tsv.gz')), sep = '\t')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0('transfer_', name, '.pdf')), p,
    device = cairo_pdf, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0('transfer_', name, '.png')), p,
    width = width, height = height, dpi = 130)
}
cor_safe <- function(x, y) {
  good <- is.finite(x) & is.finite(y)
  if (sum(good) < 4L || sd(x[good]) < 1e-10 || sd(y[good]) < 1e-10) return(NA_real_)
  cor(x[good], y[good])
}
model_order <- c('Bliss', 'Single_calibrated', 'Context', 'bank16_response', 'bank16_conditioned',
  'bank42_response', 'bank42_conditioned', 'bank43_response', 'bank43_conditioned',
  'biology_response', 'biology_conditioned', 'genes_response', 'genes_conditioned',
  'anchors_response', 'anchors_conditioned')
stopifnot(all(dd$model %in% model_order))
labels <- c(Bliss = 'Independent single effects', Single_calibrated = 'Calibrated singles', Context = 'Singles + baseline',
  bank16_response = 'Bank16 response', bank16_conditioned = 'Bank16 + baseline',
  bank42_response = 'Bank42 response', bank42_conditioned = 'Bank42 + baseline',
  bank43_response = 'Bank43 response', bank43_conditioned = 'Bank43 + baseline',
  biology_response = 'Biological readouts', biology_conditioned = 'Biological readouts + baseline',
  genes_response = 'Gene reference (12 PCs)', genes_conditioned = 'Gene PCs (12) + baseline',
  anchors_response = 'Named RNA markers', anchors_conditioned = 'Named RNA markers + baseline')

# Summarize exact drug/cell/dose duplicates without declaring separate database
# releases to be independent biology. Retain min/max spread beside every mean.
group_keys <- c('cell', 'pair', 'drug_A', 'dose_A_uM', 'drug_B', 'dose_B_uM', 'model', 'split', 'focal')
conditions <- dd[, .(viability_A = mean(viability_A), viability_B = mean(viability_B),
  viability_AB = mean(viability_AB), observed_extra_pp = mean(observed_extra_pp),
  bliss_viability = mean(bliss_viability), predicted_extra_pp = mean(predicted_extra_pp),
  predicted_viability = mean(predicted_viability),
  observed_min = min(viability_AB), observed_max = max(viability_AB),
  predicted_min = min(predicted_viability), predicted_max = max(predicted_viability),
  extra_min = min(observed_extra_pp), extra_max = max(observed_extra_pp),
  original_rows = .N, releases = paste(sort(unique(source_release)), collapse = ';'),
  plates = paste(sort(unique(plate_id)), collapse = ';')), by = group_keys]
conditions[, model_label := factor(labels[model], levels = labels[model_order])]
conditions[, split_label := factor(split,
  levels = c('training_fit', 'training_pair_OOF', 'held_family'),
  labels = c('Training fit', 'Training: excluded pairs', 'Excluded MEK-PI3K family'))]
save_table(conditions, 'condition_means_and_ranges')

# A 10 percentage-point additional effect is a prespecified descriptive size,
# not a significance threshold or a claim about assay precision.
summary <- conditions[, {
  clear <- abs(observed_extra_pp) >= 10
  .(conditions = .N, cells = uniqueN(cell), pairs = uniqueN(pair),
    informative_conditions = sum(clear), beneficial = sum(observed_extra_pp >= 10),
    antagonistic = sum(observed_extra_pp <= -10),
    correct_direction = sum(sign(predicted_extra_pp[clear]) == sign(observed_extra_pp[clear])),
    recovered_large_direction = sum(abs(predicted_extra_pp[clear]) >= 10 &
      sign(predicted_extra_pp[clear]) == sign(observed_extra_pp[clear])),
    correct_direction_fraction = if (sum(clear) > 0L)
      mean(sign(predicted_extra_pp[clear]) == sign(observed_extra_pp[clear])) else NA_real_,
    observed_extra_median_pp = median(observed_extra_pp), predicted_extra_median_pp = median(predicted_extra_pp),
    median_signed_difference_pp = median(predicted_extra_pp - observed_extra_pp),
    response_range_r = cor_safe(observed_extra_pp, predicted_extra_pp))
}, by = .(model, model_label, split, split_label, focal)]
save_table(summary, 'scientific_summary')
groups <- conditions[, .(conditions = .N, measured_extra_median_pp = median(observed_extra_pp),
  predicted_extra_median_pp = median(predicted_extra_pp),
  clear_conditions = sum(abs(observed_extra_pp) >= 10),
  clear_directions_correct = sum(abs(observed_extra_pp) >= 10 &
    sign(observed_extra_pp) == sign(predicted_extra_pp))),
  by = .(cell, pair, model, split, focal)]
save_table(groups, 'cell_pair_results')

# Compare fits and genuinely excluded pairs, with the same scientific yardstick.
ss <- summary[, .(conditions = sum(conditions), informative_conditions = sum(informative_conditions),
  correct_direction_fraction = sum(correct_direction) / sum(informative_conditions),
  large_direction_fraction = sum(recovered_large_direction) / sum(informative_conditions),
  beneficial = sum(beneficial), antagonistic = sum(antagonistic)),
  by = .(model, model_label, split, split_label)]
ss_long <- melt(ss, id.vars = c('model', 'model_label', 'split', 'split_label', 'conditions',
  'informative_conditions', 'beneficial', 'antagonistic'),
  measure.vars = c('correct_direction_fraction', 'large_direction_fraction'),
  variable.name = 'criterion', value.name = 'fraction')
ss_long[, criterion := factor(criterion, levels = c('correct_direction_fraction', 'large_direction_fraction'),
  labels = c('Correct sign (even a tiny prediction)', 'Correct sign AND predicted effect >=10 pp'))]
p <- ggplot(ss_long, aes(fraction, model_label, color = split_label)) +
  geom_point(position = position_dodge(width = .45), size = 2) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  facet_wrap(~criterion, nrow = 1) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom') +
  labs(x = 'Correct direction among measured additional effects of at least 10 percentage points',
    y = NULL, color = NULL, title = 'Getting the sign right is insufficient if predictions miss the added effect',
    caption = 'Each exact cell-pair-dose condition contributes once; repeated releases are averaged.\nPositive extra effect = fewer surviving cells than expected from independent single effects. The 10 pp cutoff is descriptive.')
save_plot(p, 'direction_summary', 13, 6)

# All conditioned representations are retained, with the baseline-only context
# model. Association panels use square geometry and independently scaled axes.
conditioned <- c('Context', 'bank16_conditioned', 'bank42_conditioned',
  'bank43_conditioned', 'biology_conditioned', 'genes_conditioned', 'anchors_conditioned')
scatter <- conditions[model %in% conditioned]
ann <- scatter[, .(label = sprintf('r = %.2f; n = %d', cor_safe(observed_extra_pp, predicted_extra_pp), .N)),
  by = .(model_label, split_label)]
p <- ggplot(scatter, aes(observed_extra_pp, predicted_extra_pp)) +
  geom_hline(yintercept = 0, color = 'grey80') + geom_vline(xintercept = 0, color = 'grey80') +
  geom_point(aes(color = cell), size = .45, alpha = .3) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey25', linewidth = .45) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.25, size = 2.6) +
  facet_grid(model_label ~ split_label, scales = 'free') +
  theme_classic(base_size = 8.5) + theme(aspect.ratio = 1, legend.position = 'bottom') +
  labs(x = 'Measured additional effect (percentage points)', y = 'Predicted additional effect (percentage points)',
    color = 'Cell', title = 'Training fit, excluded training pairs and excluded target family',
    caption = 'Each dot is one cell-pair-dose condition. Correlations span doses and drugs, not fixed-dose differences between cells.\nAll conditioned representations shown; held outcomes do not choose a model.')
save_plot(p, 'training_exclusion_sanity', 10.5, 19)

# The focal combination is displayed in complete dose curves, keeping the
# anchor concentration fixed in each page and cell. No favored dose is selected.
focal <- conditions[focal == TRUE & split == 'held_family' & model %in% conditioned]
stopifnot(nrow(focal) > 0L)
focal[, context := paste(cell, drug_A, paste0(format(dose_A_uM, trim = TRUE), ' uM'), sep = ' | ')]
focal[, variable_axis := paste0(drug_B, ' dose (uM)')]
source_lines <- list()
for (kind in c('viability_A', 'viability_B', 'viability_AB', 'bliss_viability', 'predicted_viability')) {
  z <- copy(focal)
  z[, `:=`(value = get(kind), quantity = kind)]
  source_lines[[kind]] <- z
}
lines <- rbindlist(source_lines)
line_labels <- c(viability_A = 'Anchor drug alone', viability_B = 'Varied drug alone',
  viability_AB = 'Measured combination', bliss_viability = 'Independent single effects',
  predicted_viability = 'Predicted combination')
line_colors <- c(viability_A = '#999999', viability_B = '#B99C58',
  viability_AB = '#222222', bliss_viability = '#7B68A6', predicted_viability = '#087F8C')
pdf_path <- file.path(root, 'figures/transfer_focal_dose_curves.pdf')
cairo_pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
for (ctx in unique(focal$context)) {
  zz <- lines[context == ctx]
  observed_range <- focal[context == ctx]
  pp <- ggplot(zz, aes(dose_B_uM, value, color = quantity, group = quantity)) +
    geom_hline(yintercept = 100, color = 'grey85') + geom_line(linewidth = .55) +
    geom_point(size = 1.1) +
    geom_errorbar(data = observed_range, aes(x = dose_B_uM, ymin = observed_min, ymax = observed_max),
      inherit.aes = FALSE, width = 0, color = 'grey30') +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = .001)) +
    scale_color_manual(values = line_colors, labels = line_labels, name = NULL) +
    facet_wrap(~model_label, ncol = 3) + theme_classic(base_size = 9) +
    theme(legend.position = 'bottom') +
    labs(x = unique(zz$variable_axis)[1], y = 'Viability (% vehicle)',
      title = paste('Excluded trametinib + pictilisib:', ctx),
      subtitle = 'Same observations in every panel; each panel shows a different frozen molecular representation.',
      caption = 'Dots average exact cell/drug/dose matches across records. Vertical bars show the observed range across repeats/releases, not confidence intervals.\nAll RNA-supported doses are shown; unsupported low doses are excluded. The x axis uses a log-like scale including zero.')
  print(pp)
}
dev.off()

# A single-page focal summary retains the paired observations, without ranking
# representations by their held-family success.
fs <- focal[, .(measured_median_pp = median(observed_extra_pp),
  predicted_median_pp = median(predicted_extra_pp),
  measured_lo = quantile(observed_extra_pp, .25), measured_hi = quantile(observed_extra_pp, .75),
  predicted_lo = quantile(predicted_extra_pp, .25), predicted_hi = quantile(predicted_extra_pp, .75),
  conditions = .N), by = .(cell, model_label)]
save_table(fs, 'focal_effect_summary')
pl <- rbind(fs[, .(cell, model_label, value = measured_median_pp, lo = measured_lo, hi = measured_hi,
  quantity = 'Measured')], fs[, .(cell, model_label, value = predicted_median_pp, lo = predicted_lo,
  hi = predicted_hi, quantity = 'Predicted')])
p <- ggplot(pl, aes(value, model_label, color = quantity)) +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, position = position_dodge(width = .4), orientation = 'y') +
  geom_point(position = position_dodge(width = .4), size = 1.8) +
  facet_wrap(~cell, ncol = 2) + theme_classic(base_size = 9) + theme(legend.position = 'bottom') +
  scale_color_manual(values = c(Measured = '#222222', Predicted = '#087F8C')) +
  labs(x = 'Additional combination effect beyond independent singles (percentage points)', y = NULL, color = NULL,
    title = 'Excluded trametinib + pictilisib: does each model recover the added benefit?',
    caption = 'Points: median across all RNA-supported dose contexts. Bars: middle half of dose-context values, not uncertainty intervals.\nThe detailed PDF retains the supported dose curves and constituent single-agent activities.')
save_plot(p, 'focal_effect_summary', 11, 7)

# Existing path fits only: compare the selected penalty with the least penalty
# that was already fitted. This diagnostic never reads held-family outcomes.
capacity_path <- file.path(root, 'data/transfer_training_capacity_diagnostic.tsv')
stopifnot(file.exists(capacity_path))
capacity <- fread(capacity_path)
stopifnot(all(c('model', 'lambda', 'selected', 'training_extra_effect_r', 'pair_OOF_extra_effect_r',
  'training_prediction_p05', 'training_prediction_p95') %in% names(capacity)))
selected_path <- copy(capacity[selected == TRUE])
selected_path[, penalty := 'Selected by training CV']
minimum_path <- capacity[order(lambda), .SD[1L], by = model]
minimum_path[, penalty := 'Lowest already-fitted penalty']
capacity_compare <- rbind(selected_path, minimum_path)
capacity_compare[, model_label := factor(labels[model], levels = labels[model_order])]
save_table(capacity_compare, 'training_capacity_comparison')
p_range <- ggplot(capacity_compare, aes(training_prediction_p95, model_label, color = penalty)) +
  geom_vline(xintercept = 10, color = 'grey70', linetype = 2) +
  geom_point(position = position_dodge(width = .45), size = 1.8) +
  theme_classic(base_size = 9) + theme(legend.position = 'bottom') +
  labs(x = '95th percentile of predicted training extra effect (pp)', y = NULL, color = NULL,
    title = 'How much did regularization suppress predicted effects?')
p_oof <- ggplot(capacity_compare, aes(pair_OOF_extra_effect_r, model_label, color = penalty)) +
  geom_point(position = position_dodge(width = .45), size = 1.8) +
  theme_classic(base_size = 9) + theme(legend.position = 'bottom') +
  labs(x = 'Extra-effect correlation across excluded training-pair conditions', y = NULL, color = NULL,
    title = 'What happens within the training pool when the penalty is reduced?')
combined_capacity <- (p_range | p_oof) + plot_layout(guides = 'collect') +
  plot_annotation(caption = paste('Training pool only; all values come from the already-fitted regularization path.',
    'No target-family outcomes choose these comparisons or alter the frozen predictions.',
    'A wider prediction range is not proof of accurate effects; the right panel evaluates excluded pairs within training.', sep = '\n'))
combined_capacity <- combined_capacity & theme(legend.position = 'bottom')
save_plot(combined_capacity, 'training_capacity_diagnostic', 12, 7)
save_table(data.table(path = path, sha256 = digest(file = path, algo = 'sha256')), 'prediction_manifest')
capture.output(sessionInfo(), file = file.path(root, 'transfer_plot_sessionInfo.txt'))
message('Completed transfer figures: ', uniqueN(focal$context), ' focal cell/anchor pages.')
