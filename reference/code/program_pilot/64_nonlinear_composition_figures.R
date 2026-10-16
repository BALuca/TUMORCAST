#!/usr/bin/env Rscript
# Separate gains from nonlinear individual calibration and nonlinear composition.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/64_nonlinear_response_composition'
previous <- '/work/results/program_pilot/62_separate_response_interaction_bridge'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
setDTthreads(2L)
stopifnot(grepl('complete', readLines(file.path(root, 'STATUS.txt'))))
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
save_plot <- function(x, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), x, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), x, width = width, height = height, dpi = 150)
}
safe_cor <- function(x, y) {
  if (length(x) < 4 || sd(x) < 1e-10 || sd(y) < 1e-10) return(NA_real_)
  cor(x, y)
}
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
  legend.position = 'bottom', plot.title.position = 'plot'))
audit <- fread(file.path(root, 'data/validation.tsv.gz'))
stopifnot(nrow(audit) == 168L, all(audit$source_overlap == 0), all(audit$replay_error < 1e-10),
  all(audit$zero_input_error < 1e-10))
old <- fread(file.path(previous, 'data/predictions.tsv.gz'))[
  model %chin% c('Original40_interaction', 'Original40_bridge')]
new <- fread(file.path(root, 'data/predictions.tsv.gz'))
old[, method := fifelse(model == 'Original40_interaction', 'Linear RNA interaction', 'Linear molecular bridge')]
new[, method := fcase(model == 'Nonlinear_control', 'Nonlinear functional control',
  model == 'Nonlinear_conditioned', 'Nonlinear RNA interaction', model == 'Nonlinear_g',
  'Nonlinear interaction only', default = 'Nonlinear molecular bridge')]
comparison <- merge(new[model == 'Nonlinear_g'], old[model == 'Original40_bridge'],
  by = c('time_h', 'fold_id', 'observation_id'), suffixes = c('_new', '_old'))
stopifnot(nrow(comparison) == nrow(new[model == 'Nonlinear_g']), all(comparison$split_new == comparison$split_old),
  max(abs(comparison$ratio_AB_new - comparison$ratio_AB_old)) == 0,
  max(abs(comparison$predicted_single_sum_log2_new - comparison$predicted_single_sum_log2_old)) < 1e-12)
z <- rbind(old, new)
additive <- copy(new[model == 'Nonlinear_bridge'])
additive[, `:=`(method = 'Nonlinear singles only', model = 'Nonlinear_additive',
  predicted_ratio_AB = 2^predicted_single_sum_log2, predicted_interaction_log2 = 0)]
z <- rbind(z, additive)
svr <- fread(file.path(prior, 'data/predictions.tsv.gz'))[model == 'Singles_SVR' & split == 'held_partner']
svr <- svr[, .(time_h, fold_id = held_partner, fold_type = 'single_partner',
  method = 'Previous functional SVR', model = 'Singles_SVR', split = 'held',
  observation_id, condition_id, anchor, partner, nominal_A_uM, nominal_B_uM,
  biological_replicate, ratio_A = relative_A, ratio_B = relative_B, ratio_AB = relative_AB,
  predicted_ratio_AB = predicted_relative_AB,
  predicted_interaction_log2 = log2(predicted_relative_AB / (relative_A * relative_B)))]
z <- rbind(z, svr, fill = TRUE)
stopifnot(all(z$predicted_ratio_AB > 0), all(is.finite(z$predicted_ratio_AB)))
z[, `:=`(observed_log_interaction = log2(ratio_AB / (ratio_A * ratio_B)),
  reconstructed_log_interaction = log2(predicted_ratio_AB / (ratio_A * ratio_B)),
  observed_extra = 100 * (ratio_A * ratio_B - ratio_AB),
  predicted_extra = 100 * (ratio_A * ratio_B - predicted_ratio_AB))]

# Primary whole-partner comparison. Biological means match62's reporting.
means <- z[fold_type == 'single_partner', .(observed_pct = 100 * mean(ratio_AB),
  predicted_pct = 100 * mean(predicted_ratio_AB), low_pct = 100 * min(ratio_AB), high_pct = 100 * max(ratio_AB),
  observed_log_interaction = mean(observed_log_interaction),
  predicted_log_interaction = mean(reconstructed_log_interaction),
  fitted_g = mean(predicted_interaction_log2),
  observed_extra = mean(observed_extra), predicted_extra = mean(predicted_extra)),
  by = .(time_h, method, split, anchor, partner, nominal_A_uM, nominal_B_uM, condition_id)]
means[, `:=`(observed_centered = observed_pct - mean(observed_pct),
  predicted_centered = predicted_pct - mean(predicted_pct)),
  by = .(time_h, method, split, anchor, nominal_A_uM)]
metrics <- means[, .(n = .N, interaction_r = safe_cor(observed_log_interaction, predicted_log_interaction),
  count_r = safe_cor(observed_pct, predicted_pct),
  within_first_drug_dose_r = safe_cor(observed_centered, predicted_centered),
  median_absolute_count_difference_pp = median(abs(observed_pct - predicted_pct))),
  by = .(time_h, method, split)]
by_anchor <- means[, .(n = .N, interaction_r = safe_cor(observed_log_interaction, predicted_log_interaction),
  within_first_drug_dose_r = safe_cor(observed_centered, predicted_centered),
  median_absolute_count_difference_pp = median(abs(observed_pct - predicted_pct))),
  by = .(time_h, method, split, anchor)]
by_pair <- means[, .(n = .N, interaction_r = safe_cor(observed_log_interaction, predicted_log_interaction),
  count_r = safe_cor(observed_pct, predicted_pct), median_absolute_count_difference_pp = median(abs(observed_pct - predicted_pct))),
  by = .(time_h, method, split, anchor, partner)]
save_table(means, 'primary_condition_predictions')
save_table(metrics, 'model_summary')
save_table(by_anchor, 'model_by_first_drug')
save_table(by_pair, 'model_by_pair')

main_methods <- c('Linear RNA interaction', 'Nonlinear functional control', 'Nonlinear RNA interaction')
colors <- c('Linear RNA interaction' = '#0072B2', 'Nonlinear functional control' = '#777777',
  'Nonlinear RNA interaction' = '#D55E00', 'Linear molecular bridge' = '#0072B2',
  'Nonlinear interaction only' = '#CC79A7', 'Nonlinear molecular bridge' = '#D55E00',
  'Nonlinear singles only' = '#A6761D', 'Previous functional SVR' = '#009E73')
for (tt in c(3, 24)) {
  points <- means[time_h == tt & method %chin% main_methods]
  points[, method := factor(method, levels = main_methods)]
  points[, split := factor(split, levels = c('training', 'held'))]
  labels <- points[, .(label = sprintf('r = %.2f; n = %d',
    safe_cor(observed_log_interaction, predicted_log_interaction), .N)), by = .(method, split)]
  plot <- ggplot(points, aes(observed_log_interaction, predicted_log_interaction)) +
    geom_point(aes(color = anchor), size = 1.25, alpha = .55) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .55) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.04, vjust = 1.2, size = 3.3) +
    facet_wrap(vars(split, method), ncol = 3, scales = 'free') + theme(aspect.ratio = 1) +
    scale_color_manual(values = c(alpelisib = '#0072B2', dactolisib = '#D55E00')) +
    labs(title = paste('Does nonlinear response composition transfer better?', tt, 'h RNA'),
      subtitle = 'Individual function supplied; identical original40 inputs and complete-partner exclusions',
      x = 'Measured log2(combination / single-product)', y = 'Predicted log2 modification', color = 'First drug',
      caption = 'Zero = independent singles. -1 = half the expected cell count. Kernel width and regularization use inner training partners only.')
  save_plot(plot, paste0('interaction_training_vs_held_', tt, 'h'), 10.8, 8)

  bridge_methods <- c('Linear molecular bridge', 'Nonlinear interaction only', 'Nonlinear molecular bridge')
  bridge_points <- means[time_h == tt & split == 'held' & method %chin% bridge_methods]
  bridge_points[, method := factor(method, levels = bridge_methods)]
  bridge_labels <- bridge_points[, .(label = sprintf('r = %.2f', safe_cor(observed_centered, predicted_centered))), by = method]
  plot <- ggplot(bridge_points, aes(observed_centered, predicted_centered)) +
    geom_point(aes(color = anchor), alpha = .6, size = 1.5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .6) +
    geom_text(data = bridge_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.04, vjust = 1.2, size = 4) +
    facet_wrap(~method, nrow = 1, scales = 'free') + theme(aspect.ratio = 1) +
    scale_color_manual(values = c(alpelisib = '#0072B2', dactolisib = '#D55E00')) +
    labs(title = paste('Where nonlinearity improves the complete molecular bridge:', tt, 'h'),
      subtitle = 'Excluded partners; each first-drug dose mean removed only for evaluation',
      x = 'Measured combination deviation (percentage points)', y = 'Predicted combination deviation (percentage points)',
      color = 'First drug', caption = 'Individual outcomes can train the single-response decoder; combination outcomes involving the target partner are withheld.')
  save_plot(plot, paste0('complete_bridge_', tt, 'h'), 11, 4.8)

  for (scope in c('conditional', 'full_bridge')) {
    use <- if (scope == 'conditional') c('Linear RNA interaction', 'Nonlinear RNA interaction',
      'Nonlinear functional control', 'Previous functional SVR') else
      c('Linear molecular bridge', 'Nonlinear singles only', 'Nonlinear molecular bridge', 'Previous functional SVR')
    cairo_pdf(file.path(root, 'figures', paste0('all12_', scope, '_curves_', tt, 'h.pdf')), width = 10.5, height = 6.6)
    for (aa in c('alpelisib', 'dactolisib')) {
      for (bb in sort(unique(means$partner))) {
        d <- means[time_h == tt & split == 'held' & anchor == aa & partner == bb & method %chin% use]
        observed <- d[method == use[1]]
        line_data <- rbind(d[, .(nominal_A_uM, nominal_B_uM, method, value = predicted_pct, endpoint = 'Combination count (% vehicle)')],
          d[, .(nominal_A_uM, nominal_B_uM, method, value = predicted_extra, endpoint = 'Additional inhibition (pp)')])
        dot_data <- rbind(observed[, .(nominal_A_uM, nominal_B_uM, value = observed_pct,
          low = low_pct, high = high_pct, endpoint = 'Combination count (% vehicle)')],
          observed[, .(nominal_A_uM, nominal_B_uM, value = observed_extra,
            low = observed_extra - (high_pct - observed_pct), high = observed_extra + (observed_pct - low_pct),
            endpoint = 'Additional inhibition (pp)')])
        # For the extra-effect row, biological intervals must be computed from
        # paired per-experiment extra effects, not inferred from count ranges.
        ranges <- z[time_h == tt & fold_type == 'single_partner' & split == 'held' &
          anchor == aa & partner == bb & method == use[1],
          .(lo = min(observed_extra), hi = max(observed_extra)), by = .(nominal_A_uM, nominal_B_uM)]
        dot_data <- merge(dot_data, ranges, by = c('nominal_A_uM', 'nominal_B_uM'))
        dot_data[endpoint == 'Additional inhibition (pp)', `:=`(low = lo, high = hi)]
        order <- c('Combination count (% vehicle)', 'Additional inhibition (pp)')
        line_data[, endpoint := factor(endpoint, levels = order)]
        dot_data[, endpoint := factor(endpoint, levels = order)]
        plot <- ggplot(line_data, aes(nominal_A_uM, value, color = method)) +
          geom_line(linewidth = .7) +
          geom_linerange(data = dot_data, aes(x = nominal_A_uM, ymin = low, ymax = high),
            inherit.aes = FALSE, color = 'grey50', linewidth = .45) +
          geom_point(data = dot_data, aes(nominal_A_uM, value), inherit.aes = FALSE, color = 'black', size = 1.8) +
          facet_grid(endpoint ~ nominal_B_uM, scales = 'free_y',
            labeller = labeller(nominal_B_uM = function(x) paste(bb, x, 'µM'))) +
          scale_x_log10(breaks = sort(unique(d$nominal_A_uM))) + scale_color_manual(values = colors) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = paste(aa, '+', bb, '|', tt, 'h RNA'),
            subtitle = if (scope == 'conditional') 'Interaction comparison with measured individual effects supplied' else
              'Complete molecular bridge; additive nonlinear singles separate from learned joint modification',
            x = paste(aa, 'dose (µM)'), y = NULL, color = NULL,
            caption = 'Black points/ranges: measured biological mean/range. All combinations containing this partner were excluded.')
        print(plot)
        if (tt == 3 && aa == 'alpelisib' && bb == 'lapatinib') {
          ggsave(file.path(root, 'figures', paste0(scope, '_curve_example.png')), plot, width = 10.5, height = 6.6, dpi = 150)
        }
      }
    }
    dev.off()
  }
}

# Same-model contrasts retain the stricter two-partner exclusion from62.
h <- z[fold_type == 'two_partners' & split == 'held']
keys <- c('time_h', 'method', 'fold_id', 'anchor', 'nominal_A_uM', 'nominal_B_uM', 'biological_replicate')
c <- merge(h, h, by = keys, suffixes = c('_1', '_2'), allow.cartesian = TRUE)[partner_1 < partner_2]
c[, `:=`(count_difference = 100 * (ratio_AB_1 - ratio_AB_2),
  predicted_count_difference = 100 * (predicted_ratio_AB_1 - predicted_ratio_AB_2),
  extra_difference = observed_extra_1 - observed_extra_2,
  predicted_extra_difference = predicted_extra_1 - predicted_extra_2)]
cm <- c[, .(observed = mean(count_difference), predicted = mean(predicted_count_difference),
  observed_extra = mean(extra_difference), predicted_extra = mean(predicted_extra_difference),
  clear_count = abs(mean(count_difference)) >= 10 && all(sign(count_difference) == sign(mean(count_difference))),
  clear_extra = abs(mean(extra_difference)) >= 10 && all(sign(extra_difference) == sign(mean(extra_difference)))),
  by = .(time_h, method, fold_id, anchor, nominal_A_uM, nominal_B_uM)]
contrast_summary <- cm[, .(n = .N, count_r = safe_cor(observed, predicted), extra_r = safe_cor(observed_extra, predicted_extra),
  clear_count_n = sum(clear_count), count_directions_correct = sum(clear_count & sign(observed) == sign(predicted)),
  clear_extra_n = sum(clear_extra), extra_directions_correct = sum(clear_extra & sign(observed_extra) == sign(predicted_extra))),
  by = .(time_h, method)]
save_table(cm, 'common_model_partner_contrasts')
save_table(contrast_summary, 'common_model_contrast_summary')

# Decompose the apparent gain: can f reproduce the known individual responses,
# and does g improve held combinations beyond these calibrated individual effects?
single <- fread(file.path(root, 'data/single_predictions.tsv.gz'))
single_long <- melt(single, id.vars = c('time_h', 'source', 'drug', 'nominal_dose_uM', 'ratio'),
  measure.vars = c('predicted_fitted_log2', 'predicted_nonlinear_fitted_log2'),
  variable.name = 'method', value.name = 'prediction_log2')
single_long[, method := fifelse(method == 'predicted_fitted_log2', 'Linear single decoder', 'Nonlinear single decoder')]
single_means <- single_long[, .(observed_pct = 100 * mean(ratio), low = 100 * min(ratio),
  high = 100 * max(ratio), predicted_pct = 100 * mean(2^prediction_log2)),
  by = .(time_h, source, drug, nominal_dose_uM, method)]
single_metrics <- single_means[, .(conditions = .N, count_r = safe_cor(observed_pct, predicted_pct),
  median_absolute_count_difference_pp = median(abs(observed_pct - predicted_pct))),
  by = .(time_h, source, method)]
save_table(single_means, 'single_effect_calibration')
save_table(single_metrics, 'single_effect_calibration_summary')
contrast_by_anchor <- cm[, .(n = .N, count_r = safe_cor(observed, predicted),
  clear_count_n = sum(clear_count), count_directions_correct = sum(clear_count & sign(observed) == sign(predicted))),
  by = .(time_h, method, anchor)]
save_table(contrast_by_anchor, 'common_model_contrasts_by_first_drug')

for (tt in c(3, 24)) {
  d <- single_means[time_h == tt & source == 'Same-screen single controls']
  dots <- d[method == 'Linear single decoder']
  p <- ggplot(d, aes(nominal_dose_uM, predicted_pct, color = method)) +
    geom_line(linewidth = .75) +
    geom_linerange(data = dots, aes(x = nominal_dose_uM, ymin = low, ymax = high),
      inherit.aes = FALSE, color = 'grey50') +
    geom_point(data = dots, aes(nominal_dose_uM, observed_pct), inherit.aes = FALSE, size = 1.8) +
    facet_wrap(~drug, ncol = 4, scales = 'free_x') + scale_x_log10() +
    scale_color_manual(values = c('Linear single decoder' = '#0072B2', 'Nonlinear single decoder' = '#D55E00')) +
    labs(title = paste('Does the molecular decoder recover individual effects?', tt, 'h'),
      subtitle = 'Known single-drug experiments used to train f; these are fitted values, not held-out single-drug predictions',
      x = 'Individual drug dose (µM)', y = 'Viable-cell count (% vehicle)', color = NULL,
      caption = 'Black points/ranges: biological mean/range. The same f is used for every combination exclusion.')
  save_plot(p, paste0('single_effect_calibration_', tt, 'h'), 10, 6.3)

  use <- c('Linear molecular bridge', 'Nonlinear singles only', 'Nonlinear molecular bridge')
  d <- means[time_h == tt & split == 'held' & method %chin% use]
  d[, method := factor(method, levels = use)]
  labs <- d[, .(label = sprintf('r = %.2f\nmedian error = %.1f pp',
    safe_cor(observed_centered, predicted_centered), median(abs(observed_pct - predicted_pct)))),
    by = .(method, anchor)]
  p <- ggplot(d, aes(observed_centered, predicted_centered)) +
    geom_point(aes(color = partner), alpha = .75, size = 1.5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .6) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.04, vjust = 1.1, size = 3.2) +
    facet_wrap(vars(anchor, method), ncol = 3, scales = 'free') + theme(aspect.ratio = 1) +
    labs(title = paste('Separate individual-effect calibration from learned interaction:', tt, 'h'),
      subtitle = 'All combinations containing the predicted partner excluded; all six partners shown',
      x = 'Measured deviation at fixed first-drug dose (pp)', y = 'Predicted deviation at fixed first-drug dose (pp)',
      color = 'Excluded partner',
      caption = 'Correlations compare partners and partner doses within one BT20 cell line. Error uses uncentered cell counts.')
  save_plot(p, paste0('bridge_ablation_by_first_drug_', tt, 'h'), 10.8, 8.2)

  use <- c('Nonlinear functional control', 'Nonlinear RNA interaction', 'Nonlinear molecular bridge')
  d <- cm[time_h == tt & method %chin% use]
  d[, method := factor(method, levels = use)]
  labs <- d[, .(label = sprintf('r = %.2f\n%d/%d clear directions', safe_cor(observed, predicted),
    sum(clear_count & sign(observed) == sign(predicted)), sum(clear_count))), by = .(method, anchor)]
  p <- ggplot(d, aes(observed, predicted)) +
    geom_vline(xintercept = 0, color = 'grey80') + geom_hline(yintercept = 0, color = 'grey80') +
    geom_point(aes(color = clear_count), size = 1.4, alpha = .55) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .55) +
    geom_text(data = labs, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.04, vjust = 1.1, size = 3.2) +
    facet_wrap(vars(anchor, method), ncol = 3, scales = 'free') + theme(aspect.ratio = 1) +
    scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = '#0072B2'),
      labels = c('Below 10 pp or replicate-discordant', 'At least 10 pp, consistent replicate direction')) +
    labs(title = paste('Can one model distinguish two unseen partners?', tt, 'h'),
      subtitle = 'Both partners excluded together; matched first-drug dose, partner dose and biological experiment',
      x = 'Measured difference in combination cell count (pp)', y = 'Predicted difference in combination cell count (pp)',
      color = NULL, caption = 'All contrasts shown. Correlation uses all points; direction count uses only the predeclared clear subset.')
  save_plot(p, paste0('two_excluded_partner_contrasts_', tt, 'h'), 10.8, 8)
}

selection <- fread(file.path(root, 'data/selection.tsv.gz'))
save_table(selection[, .(folds = .N, selected_zero = sum(zero_interaction)),
  by = .(time_h, fold_type, model)], 'selected_zero_interactions')
coverage <- new[split == 'held', .N, by = .(time_h, model, fold_type, observation_id)]
stopifnot(all(coverage[fold_type == 'single_partner', N] == 1L),
  all(coverage[fold_type == 'two_partners', N] == 5L),
  nrow(coverage) == 2L * 4L * 2L * 496L)
save_table(coverage[, .(observations = .N, min_appearances = min(N), max_appearances = max(N)),
  by = .(time_h, model, fold_type)], 'prediction_coverage')
stopifnot(!anyNA(fread(file.path(root, 'data/single_selection.tsv.gz'))))

report <- c('WORKFLOW64: NONLINEAR SEPARATE-RESPONSE COMPOSITION', '',
  'Question: does nonlinear calibration improve the original40 molecular bridge?',
  'Same original Tahoe+Cytosig40 bank; same data and partner exclusions as62.',
  'BT20 only: two first drugs (alpelisib/dactolisib), six partners,198 dose conditions.',
  'The six partners are lapatinib, neratinib, NVP-TAE684, saracatinib, selumetinib, trametinib.',
  'Each primary fold excludes ALL combinations involving one partner: both first drugs,',
  'every dose and every biological replicate. Ten distinct pairs remain in training.',
  'Single-drug molecular AND functional observations remain available, including target partners.',
  'No combination RNA, excluded combination outcomes, drug identity or replicate identity enters fitting.',
  'This is transfer from known individual experiments to unseen joint outcomes, not unseen-drug singles.', '',
  'MODEL: log2(relative count AB) = f(xA) + f(xB) + g(xA,xB).',
  'f is origin-anchored Gaussian kernel ridge; g uses its symmetric pair kernel.',
  'This permits smooth nonlinear effects; f(0)=0 and g(0,x)=g(x,0)=0.',
  'Three kernel widths and ridge regularization chosen on inner training exclusions.',
  'Individual f uses drug-grouped tuning on individual measurements only.',
  'Interaction g uses equal raw/within-first-drug-dose residual objective, as in62.',
  'Kernel inputs and scales use available individual inputs, never held joint outcomes.',
  'Full molecular bridge uses RNA inputs at prediction; f was calibrated on known individual function.',
  'Conditional RNA test supplies measured singles, then adds RNA+function/dose interaction kernels.',
  'Equally flexible functional control supplies the same singles with function/dose kernel only.', '',
  'PRIMARY3H RESULTS:', capture.output(print(metrics[time_h == 3])), '',
  'HELD RESULTS BY FIRST DRUG:', capture.output(print(by_anchor[time_h == 3 & split == 'held'])), '',
  'SINGLE-EFFECT FITTED CALIBRATION (not a held-out performance claim):',
  capture.output(print(single_metrics[time_h == 3])), '',
  'STRICTER TWO-PARTNER EXCLUSION:', capture.output(print(contrast_summary[time_h == 3])), '',
  'INTERPRETATION:',
  'Nonlinear f fixes much of the attenuation in the linear molecular bridge.',
  'The alpelisib combination differences are captured substantially better.',
  'Dactolisib partner differences remain weak; pooling the first drugs hides this failure.',
  'Calibrated nonlinear singles already explain the improved ranking. Adding g improves absolute',
  'joint-count calibration, but does not improve within-dose ranking over nonlinear singles alone.',
  'RNA interactions do not clearly outperform the equally flexible functional-only control.',
  'The conditional RNA training/held gap persists. This is not yet evidence that the molecular',
  'representation contributes distinctive information about unseen functional interactions.', '',
  'HOW TO READ:',
  'bridge_ablation_by_first_drug_3h.pdf: start here; isolates f versus g and both first drugs.',
  'single_effect_calibration_3h.pdf: which known individual effects were lost by linear calibration.',
  'interaction_training_vs_held_3h.pdf: test molecular value against the nonlinear functional control.',
  'all12_full_bridge_curves_3h.pdf: one page per pair, actual counts and extra inhibition, biological ranges.',
  'two_excluded_partner_contrasts_3h.pdf: same model predicts both excluded partners at matched doses.',
  '24h versions and tables include every pair, both successful and unsuccessful.', '',
  'MEASUREMENT NOTES:',
  'Count = viable cells as percentage of vehicle; pp = percentage points.',
  'Within-first-drug-dose correlations remove dose means only for evaluation, not to hide dose effects;',
  'the accompanying curves retain absolute dose responses. Variation is across partners and partner doses,',
  'not across cell lines. Biological replicates are averaged for correlations;198 conditions are not198 drugs.',
  'Training summaries average a condition across fits where it was in training, as in62.',
  'interaction_r uses predicted count relative to the MEASURED single-product for all models.',
  'For full molecular bridges this combines f error and g; shared references can inflate it.',
  'Interpret actual-count predictions, paired contrasts and the additive ablation together.',
  'Clear contrasts: absolute observed difference >=10 pp and same direction in all biological replicates.',
  'All outcomes have been inspected before. Retrospective development, not independent confirmation.',
  'One cell line cannot establish baseline conditioning, new-cell transfer or HGF rescue.',
  'The original40 bank has partial L1000 coverage; unchanged from61/62.', '',
  'VALIDATION:168 models; source-well exclusions, swap symmetry, zero inputs, PSD tolerance,',
  'saved prediction replay, frozen linear offsets, complete held coverage and finite tuning tables checked.',
  'Every glmnet fit returned jerr=0. Predictions saved and hashed before outcome evaluation.')
writeLines(report, file.path(root, 'REPORT.txt'))
writeLines('Models, comparisons and figures complete.', file.path(root, 'STATUS.txt'))
fwrite(data.table(path = list.files(file.path(root, 'figures'), pattern = '\\.pdf$', full.names = TRUE)),
  file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
scripts <- file.path('/work/src/program_pilot', c('64_nonlinear_response_composition.R', '64_nonlinear_composition_figures.R'))
save_table(data.table(path = scripts, sha256 = vapply(scripts, digest, character(1), file = TRUE, algo = 'sha256')), 'script_manifest')
print(metrics[time_h == 3])
print(by_anchor[time_h == 3 & split == 'held'])
print(single_metrics[time_h == 3])
print(contrast_summary[time_h == 3])

