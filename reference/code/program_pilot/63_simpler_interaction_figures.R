#!/usr/bin/env Rscript
# Direct complexity comparison on identical saved data, offsets and exclusions.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/63_simpler_program_interactions'
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
stopifnot(nrow(audit) == 84L, all(audit$source_overlap == 0), all(audit$molecular_terms == 40),
  all(audit$frozen_single_offset_error < 1e-12), all(audit$replay_error < 1e-10),
  all(audit$symmetry_error < 1e-12), all(audit$zero_input_error < 1e-12))
old <- fread(file.path(previous, 'data/predictions.tsv.gz'))[
  model %chin% c('Functional_control', 'Original40_interaction', 'Original40_bridge')]
new <- fread(file.path(root, 'data/predictions.tsv.gz'))
old[, method := fcase(model == 'Functional_control', 'Functional control',
  model == 'Original40_interaction', '820 cross-program terms', default = '820 terms: full bridge')]
new[, method := fifelse(model == 'Original40_interaction', '40 same-program terms', '40 terms: full bridge')]
comparison <- merge(new, old[model != 'Functional_control'],
  by = c('time_h', 'fold_id', 'model', 'observation_id'), suffixes = c('_new', '_old'))
stopifnot(nrow(comparison) == nrow(new), all(comparison$split_new == comparison$split_old),
  max(abs(comparison$ratio_AB_new - comparison$ratio_AB_old)) == 0,
  max(abs(comparison$predicted_single_sum_log2_new - comparison$predicted_single_sum_log2_old)) < 1e-12)
z <- rbind(old, new)
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
  observed_extra = 100 * (ratio_A * ratio_B - ratio_AB),
  predicted_extra = 100 * (ratio_A * ratio_B - predicted_ratio_AB))]

# Primary whole-partner comparison. Biological means match62's reporting.
means <- z[fold_type == 'single_partner', .(observed_pct = 100 * mean(ratio_AB),
  predicted_pct = 100 * mean(predicted_ratio_AB), low_pct = 100 * min(ratio_AB), high_pct = 100 * max(ratio_AB),
  observed_log_interaction = mean(observed_log_interaction),
  predicted_log_interaction = mean(predicted_interaction_log2),
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
save_table(metrics, 'complexity_summary')
save_table(by_anchor, 'complexity_by_first_drug')
save_table(by_pair, 'complexity_by_pair')

main_methods <- c('Functional control', '820 cross-program terms', '40 same-program terms')
colors <- c('Functional control' = '#777777', '820 cross-program terms' = '#0072B2',
  '40 same-program terms' = '#D55E00', '820 terms: full bridge' = '#0072B2',
  '40 terms: full bridge' = '#D55E00', 'Previous functional SVR' = '#009E73')
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
    labs(title = paste('Does a smaller interaction model transfer better?', tt, 'h RNA'),
      subtitle = 'Same data, original40 bank, single-effect offset, ridge selection and complete-partner exclusions',
      x = 'Measured log2(combination / single-product)', y = 'Predicted log2 modification', color = 'First drug',
      caption = 'Zero = independent singles. -1 = half the expected cell count. 40-term model allows only same-program interactions.')
  save_plot(plot, paste0('complexity_training_vs_held_', tt, 'h'), 10.8, 8)

  for (scope in c('conditional', 'full_bridge')) {
    use <- if (scope == 'conditional') c('820 cross-program terms', '40 same-program terms', 'Previous functional SVR') else
      c('820 terms: full bridge', '40 terms: full bridge', 'Previous functional SVR')
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
              'Complete molecular bridge; identical frozen single-response decoder',
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

# Every program coefficient remains visible. Sign consistency across outer fits
# is algorithmic stability, not independent biological evidence.
coef <- fread(file.path(root, 'data/interaction_coefficients.tsv.gz'))[
  fold_type == 'single_partner' & model == 'Original40_interaction' & !grepl('^functional_', term)]
coef[, program := sub('__.*', '', term)]
bank <- qread(file.path(prior, 'data/representation_scores.qs'))
annotations <- bank$program_annotations
for (tt in c(3, 24)) {
  wide <- dcast(coef[time_h == tt], program ~ fold_id, value.var = 'raw_coefficient')
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$program
  stopifnot(nrow(mat) == 40L, ncol(mat) == 6L, all(is.finite(mat)))
  mat <- mat[hclust(dist(mat), method = 'ward.D2')$order, hclust(dist(t(mat)), method = 'ward.D2')$order, drop = FALSE]
  labels <- substr(annotations$display_name[match(rownames(mat), annotations$program)], 1, 58)
  limit <- max(abs(mat))
  if (limit < 1e-10) limit <- 1e-10
  hm <- draw_heatmap(mat, scale_rows = FALSE, value_range = c(-limit, 0, limit),
    value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'Interaction\ncoefficient',
    show_row_names = TRUE, show_column_names = TRUE, row_labels = labels,
    row_names_gp = grid::gpar(fontsize = 7.5), column_names_gp = grid::gpar(fontsize = 8), column_names_rot = 45,
    title = paste('Same-program interaction coefficients:', tt, 'h | columns identify excluded partners'))
  cairo_pdf(file.path(root, 'figures', paste0('program_coefficient_consistency_', tt, 'h.pdf')), width = 9, height = 8)
  ComplexHeatmap::draw(hm)
  dev.off()
}
selection <- fread(file.path(root, 'data/interaction_selection.tsv.gz'))[selected == TRUE]
save_table(selection[, .(folds = .N, selected_zero = sum(is.infinite(lambda))),
  by = .(time_h, fold_type, model)], 'selected_zero_interactions')
report <- c('WORKFLOW63: SIMPLER SAME-PROGRAM INTERACTIONS', '',
  '62 used ridge regression with820 molecular interaction coefficients.',
  '63 uses40 same-program products, keeping the frozen40-coefficient single decoder.',
  'Equation: log2 count = beta dot (xA+xB) + sum_k gamma_k*xA_k*xB_k.',
  'The conditional test replaces beta dot (xA+xB) with measured individual log2 effects',
  'and retains exactly the same six function/dose control terms as62.',
  'Each primary fold trains on10 distinct pairs and excludes the target partner entirely.',
  'All84 models completed; input-order symmetry, zero inputs, source exclusion,',
  'frozen-single-offset identity and saved model replay checks passed.', '',
  'PRIMARY3H: TRAINING AND EXCLUDED-PARTNER INTERACTION:',
  capture.output(print(metrics[time_h == 3 & method %chin% main_methods])), '',
  'HELD PREDICTIONS BY FIRST DRUG (same first-drug dose comparisons):',
  capture.output(print(by_anchor[time_h == 3 & split == 'held'])), '',
  'STRICTER TWO-PARTNER CONTRASTS:', capture.output(print(contrast_summary[time_h == 3])), '',
  '24H TEMPORAL SENSITIVITY:', capture.output(print(metrics[time_h == 24 & method %chin% main_methods])), '',
  'The original single-response calibration is unchanged. This run isolates g complexity.',
  'A smaller training/test gap alone is not success; excluded-partner prediction must improve.',
  'Coefficient consistency across folds is descriptive; six partners are not198 independent interventions.',
  'All outcomes have been inspected before this analysis. This is retrospective model development.',
  'High derived extra-effect correlations for full bridges can share the individual-product reference;',
  'actual joint-count contrasts and strong functional comparators remain necessary.',
  'Figures: complexity_training_vs_held_3h.pdf first; then all12_conditional_curves_3h.pdf',
  'and all12_full_bridge_curves_3h.pdf. All12 pairs and both times remain visible.')
writeLines(report, file.path(root, 'REPORT.txt'))
writeLines('Models, comparisons and figures complete.', file.path(root, 'STATUS.txt'))
fwrite(data.table(path = list.files(file.path(root, 'figures'), pattern = '\\.pdf$', full.names = TRUE)),
  file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
print(metrics[time_h == 3])
print(contrast_summary[time_h == 3])
