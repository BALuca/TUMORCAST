#!/usr/bin/env Rscript
# Evaluate frozen predictions using whole curves and same-model partner contrasts.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/62_separate_response_interaction_bridge'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data'
setDTthreads(2L)
options(warn = 1)
stopifnot(grepl('complete', readLines(file.path(root, 'STATUS.txt'))))
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 150)
}
safe_cor <- function(a, b) {
  if (length(a) < 4L || sd(a) < 1e-10 || sd(b) < 1e-10) return(NA_real_)
  cor(a, b)
}
model_labels <- c(Functional_control = 'Individual function + dose',
  Original40_interaction = 'Same inputs + original 40 RNA',
  Original40_bridge = 'Original 40: molecular f + g',
  GenePC12_interaction = 'Same inputs + gene-response PCs',
  GenePC12_bridge = 'Gene PCs: molecular f + g',
  Previous_functional_SVR = 'Previous functional SVR')
main_models <- c('Functional_control', 'Original40_interaction', 'GenePC12_interaction')
anchor_colors <- c(alpelisib = '#0072B2', dactolisib = '#D55E00')
model_colors <- c('Individual function + dose' = '#777777',
  'Same inputs + original 40 RNA' = '#0072B2',
  'Same inputs + gene-response PCs' = '#009E73',
  'Independent singles' = '#A6761D', 'Original 40: molecular f + g' = '#CC79A7',
  'Previous functional SVR' = '#E69F00')
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
  legend.position = 'bottom', plot.title.position = 'plot'))

p <- fread(file.path(root, 'data/predictions.tsv.gz'))
audit <- fread(file.path(root, 'data/validation.tsv.gz'))
stopifnot(nrow(audit) == 210L, all(audit$source_overlap == 0L),
  all(audit$zero_input_error < 1e-12), all(audit$replay_error < 1e-10),
  all(is.finite(p$predicted_ratio_AB)), !anyDuplicated(p[, .(time_h, fold_id, model, observation_id)]))
p[, `:=`(expected_pct = 100 * ratio_A * ratio_B, measured_pct = 100 * ratio_AB,
  predicted_pct = 100 * predicted_ratio_AB,
  measured_extra_pp = 100 * (ratio_A * ratio_B - ratio_AB),
  predicted_extra_pp = 100 * (ratio_A * ratio_B - predicted_ratio_AB),
  measured_interaction_log2 = log2(ratio_AB / (ratio_A * ratio_B)))]

# One model supplies both sides: neither partner occurred in its training pairs.
held <- p[split == 'held' & fold_type == 'two_partners']
keys <- c('time_h', 'fold_id', 'model', 'anchor', 'nominal_A_uM', 'nominal_B_uM', 'biological_replicate')
contrasts <- merge(held, held, by = keys, suffixes = c('_1', '_2'), allow.cartesian = TRUE)[partner_1 < partner_2]
stopifnot(nrow(contrasts) > 0L, all(contrasts$fold_id == paste(contrasts$partner_1, contrasts$partner_2, sep = '__')))
contrasts[, `:=`(measured_count_difference = measured_pct_1 - measured_pct_2,
  predicted_count_difference = predicted_pct_1 - predicted_pct_2,
  measured_extra_difference = measured_extra_pp_1 - measured_extra_pp_2,
  predicted_extra_difference = predicted_extra_pp_1 - predicted_extra_pp_2,
  single_expectation_difference = expected_pct_1 - expected_pct_2)]
save_table(contrasts, 'matched_partner_contrast_replicates')
contrast_means <- contrasts[, .(n_bio = .N,
  measured_count_difference = mean(measured_count_difference),
  predicted_count_difference = mean(predicted_count_difference),
  measured_count_low = min(measured_count_difference), measured_count_high = max(measured_count_difference),
  measured_extra_difference = mean(measured_extra_difference),
  predicted_extra_difference = mean(predicted_extra_difference),
  measured_extra_low = min(measured_extra_difference), measured_extra_high = max(measured_extra_difference),
  count_direction_repeated = all(sign(measured_count_difference) == sign(mean(measured_count_difference))),
  extra_direction_repeated = all(sign(measured_extra_difference) == sign(mean(measured_extra_difference))),
  single_expectation_difference = mean(single_expectation_difference)),
  by = .(time_h, fold_id, model, anchor, nominal_A_uM, nominal_B_uM, partner_1, partner_2)]
contrast_means[, `:=`(clear_count = abs(measured_count_difference) >= 10 & count_direction_repeated & n_bio >= 2,
  clear_extra = abs(measured_extra_difference) >= 10 & extra_direction_repeated & n_bio >= 2,
  similar_individual_expectation = abs(single_expectation_difference) <= 5)]
save_table(contrast_means, 'matched_partner_contrast_means')

summarize_contrasts <- function(d) {
  cc <- d$clear_count
  ee <- d$clear_extra
  list(n = nrow(d), count_r = safe_cor(d$measured_count_difference, d$predicted_count_difference),
    extra_r = safe_cor(d$measured_extra_difference, d$predicted_extra_difference),
    clear_count_n = sum(cc), clear_count_correct = sum(cc & sign(d$measured_count_difference) == sign(d$predicted_count_difference)),
    clear_extra_n = sum(ee), clear_extra_correct = sum(ee & sign(d$measured_extra_difference) == sign(d$predicted_extra_difference)),
    clear_extra_correct_and_10pp = sum(ee & sign(d$measured_extra_difference) == sign(d$predicted_extra_difference) &
      abs(d$predicted_extra_difference) >= 10),
    median_abs_extra_error_pp = median(abs(d$measured_extra_difference - d$predicted_extra_difference)),
    median_abs_measured_extra_pp = median(abs(d$measured_extra_difference)),
    median_abs_predicted_extra_pp = median(abs(d$predicted_extra_difference)))
}
metrics <- contrast_means[, summarize_contrasts(.SD), by = .(time_h, model)]
anchor_metrics <- contrast_means[, summarize_contrasts(.SD), by = .(time_h, model, anchor)]
pair_metrics <- contrast_means[, summarize_contrasts(.SD), by = .(time_h, model, anchor, fold_id)]
matched_metrics <- contrast_means[similar_individual_expectation == TRUE,
  summarize_contrasts(.SD), by = .(time_h, model, anchor)]
save_table(metrics, 'contrast_summary')
save_table(anchor_metrics, 'contrast_summary_by_anchor')
save_table(pair_metrics, 'contrast_summary_by_partner_pair')
save_table(matched_metrics, 'similar_single_effect_contrast_summary')

# Keep the previously stronger no-RNA comparator visible on its exact original
# single-partner exclusions. It is not a two-partner exclusion comparator.
previous <- fread(file.path(prior, 'predictions.tsv.gz'))[model == 'Singles_SVR' & split == 'held_partner']
stopifnot(all(previous$predicted_relative_AB > 0), nrow(previous) == 992L)
previous <- previous[, .(time_h, model = 'Previous_functional_SVR', fold_id = held_partner,
  fold_type = 'single_partner', split = 'held', observation_id, condition_id,
  anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate,
  expected_pct = 100 * relative_A * relative_B, measured_pct = 100 * relative_AB,
  predicted_pct = 100 * predicted_relative_AB,
  measured_extra_pp = 100 * (relative_A * relative_B - relative_AB),
  predicted_extra_pp = 100 * (relative_A * relative_B - predicted_relative_AB),
  measured_interaction_log2 = log2(relative_AB / (relative_A * relative_B)),
  predicted_interaction_log2 = log2(predicted_relative_AB / (relative_A * relative_B)))]
p <- rbind(p, previous, fill = TRUE)

# The primary plot retains small/inconsistent effects. Repeated >=10pp effects
# are emphasized descriptively; they never select features or fitted parameters.
for (tt in c(3, 24)) {
  z <- contrast_means[time_h == tt & model %chin% main_models]
  z[, model_label := factor(model, levels = main_models, labels = model_labels[main_models])]
  stat <- z[, .(label = sprintf('r = %.2f; %d contrasts\nRepeated >=10 pp: %d/%d directions correct',
    safe_cor(measured_extra_difference, predicted_extra_difference), .N,
    sum(clear_extra & sign(measured_extra_difference) == sign(predicted_extra_difference)), sum(clear_extra))),
    by = .(model_label, anchor)]
  plot <- ggplot(z, aes(measured_extra_difference, predicted_extra_difference)) +
    geom_hline(yintercept = 0, color = 'grey75', linewidth = .3) +
    geom_vline(xintercept = 0, color = 'grey75', linewidth = .3) +
    geom_point(aes(color = anchor, alpha = clear_extra), size = 1.5) +
    scale_alpha_manual(values = c('FALSE' = .18, 'TRUE' = .85), guide = 'none') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
    geom_text(data = stat, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.03, vjust = 1.15, size = 2.9) +
    facet_wrap(vars(anchor, model_label), ncol = 3, scales = 'free') +
    scale_color_manual(values = anchor_colors) + theme(aspect.ratio = 1) +
    labs(title = paste0('Does separate RNA predict which partner changes the interaction? ', tt, ' h RNA'),
      subtitle = 'Both partners excluded together; matched first drug and both doses; 2–3 biological experiments',
      x = 'Measured difference in additional inhibition (percentage points)',
      y = 'Predicted difference in additional inhibition (percentage points)', color = 'First drug',
      caption = paste('Positive means the alphabetically first partner produces more additional inhibition.',
        '\nFaint points: small or direction-inconsistent effects. All points enter correlations; contrasts share observations.'))
  save_plot(plot, paste0('partner_modification_', tt, 'h'), 11, 8.2)
}

# Average repeated outer predictions first, then biological experiments. The
# between-model spread is model sensitivity, not biological replication.
per_bio <- p[fold_type == 'single_partner', .(predicted_pct = mean(predicted_pct), predicted_min_pct = min(predicted_pct),
  predicted_max_pct = max(predicted_pct), measured_pct = unique(measured_pct),
  expected_pct = unique(expected_pct), measured_extra_pp = unique(measured_extra_pp),
  predicted_extra_pp = mean(predicted_extra_pp), measured_interaction_log2 = unique(measured_interaction_log2),
  predicted_interaction_log2 = mean(predicted_interaction_log2)),
  by = .(time_h, model, split, observation_id, condition_id, anchor, partner,
    nominal_A_uM, nominal_B_uM, biological_replicate)]
condition <- per_bio[, .(measured_pct = mean(measured_pct), predicted_pct = mean(predicted_pct),
  measured_low = min(measured_pct), measured_high = max(measured_pct), expected_pct = mean(expected_pct),
  measured_extra_pp = mean(measured_extra_pp), predicted_extra_pp = mean(predicted_extra_pp),
  extra_low = min(measured_extra_pp), extra_high = max(measured_extra_pp),
  predicted_interaction_log2 = mean(predicted_interaction_log2),
  measured_interaction_log2 = mean(measured_interaction_log2),
  extra_consistent = all(sign(measured_extra_pp) == sign(mean(measured_extra_pp)))),
  by = .(time_h, model, split, condition_id, anchor, partner, nominal_A_uM, nominal_B_uM)]
save_table(per_bio, 'predictions_by_biological_experiment')
save_table(condition, 'predictions_by_condition')
curve_metrics <- condition[, .(conditions = .N, count_r = safe_cor(measured_pct, predicted_pct),
  interaction_r = safe_cor(measured_interaction_log2, predicted_interaction_log2),
  additional_inhibition_r = safe_cor(measured_extra_pp, predicted_extra_pp),
  median_count_difference_pp = median(abs(measured_pct - predicted_pct)),
  mean_observed_extra_pp = mean(measured_extra_pp), mean_predicted_extra_pp = mean(predicted_extra_pp)),
  by = .(time_h, model, split, anchor, partner)]
save_table(curve_metrics, 'whole_pair_metrics')

primary <- copy(condition[split == 'held'])
primary[, `:=`(centered_count = measured_pct - mean(measured_pct),
  centered_prediction = predicted_pct - mean(predicted_pct),
  centered_extra = measured_extra_pp - mean(measured_extra_pp),
  centered_predicted_extra = predicted_extra_pp - mean(predicted_extra_pp)),
  by = .(time_h, model, anchor, nominal_A_uM)]
primary_summary <- primary[, .(n = .N,
  count_r = safe_cor(measured_pct, predicted_pct),
  within_anchor_dose_count_r = safe_cor(centered_count, centered_prediction),
  within_anchor_dose_extra_r = safe_cor(centered_extra, centered_predicted_extra),
  median_absolute_difference_pp = median(abs(measured_pct - predicted_pct))),
  by = .(time_h, model, anchor)]
save_table(primary_summary, 'primary_single_partner_summary')

for (tt in c(3, 24)) {
  # All12 pairs appear, with identical layout and complete supported dose ranges.
  cairo_pdf(file.path(root, 'figures', paste0('all12_excluded_combination_curves_', tt, 'h.pdf')), width = 11, height = 7)
  for (aa in c('alpelisib', 'dactolisib')) {
    for (bb in sort(unique(condition$partner))) {
      z <- condition[time_h == tt & split == 'held' & anchor == aa & partner == bb &
        model %chin% c('Original40_interaction', 'GenePC12_interaction', 'Original40_bridge', 'Previous_functional_SVR')]
      observed <- z[model == 'Original40_interaction']
      curves <- rbind(
        z[, .(nominal_A_uM, nominal_B_uM, label = unname(model_labels[model]),
          value = predicted_pct, endpoint = 'Combination cell count (% vehicle)')],
        z[, .(nominal_A_uM, nominal_B_uM, label = unname(model_labels[model]),
          value = predicted_extra_pp, endpoint = 'Additional inhibition (percentage points)')],
        observed[, .(nominal_A_uM, nominal_B_uM, label = 'Independent singles',
          value = expected_pct, endpoint = 'Combination cell count (% vehicle)')],
        observed[, .(nominal_A_uM, nominal_B_uM, label = 'Independent singles',
          value = 0, endpoint = 'Additional inhibition (percentage points)')])
      dots <- rbind(
        observed[, .(nominal_A_uM, nominal_B_uM, value = measured_pct, low = measured_low,
          high = measured_high, endpoint = 'Combination cell count (% vehicle)')],
        observed[, .(nominal_A_uM, nominal_B_uM, value = measured_extra_pp, low = extra_low,
          high = extra_high, endpoint = 'Additional inhibition (percentage points)')])
      endpoint_order <- c('Combination cell count (% vehicle)', 'Additional inhibition (percentage points)')
      curves[, endpoint := factor(endpoint, levels = endpoint_order)]
      dots[, endpoint := factor(endpoint, levels = endpoint_order)]
      fig <- ggplot(curves, aes(nominal_A_uM, value, color = label, group = label)) +
        geom_line(linewidth = .7) +
        geom_linerange(data = dots, aes(x = nominal_A_uM, ymin = low, ymax = high),
          inherit.aes = FALSE, color = 'grey40', linewidth = .45) +
        geom_point(data = dots, aes(nominal_A_uM, value), inherit.aes = FALSE, color = 'black', size = 1.9) +
        facet_grid(endpoint ~ nominal_B_uM, scales = 'free_y',
          labeller = labeller(nominal_B_uM = function(x) paste(bb, x, 'µM'))) +
        scale_x_log10(breaks = sort(unique(z$nominal_A_uM))) + scale_color_manual(values = model_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 8)) +
        guides(color = guide_legend(nrow = 2)) +
        labs(title = paste(aa, '+', bb, '|', tt, 'h individual RNA'),
          subtitle = 'Black: measured combination, range across biological experiments. Lines: predictions for excluded combinations.',
          x = paste(aa, 'dose (µM)'), y = NULL, color = NULL,
          caption = paste('Columns show', bb, 'dose (µM). Every combination involving this partner was withheld.',
            '\nThe molecular f + g line uses RNA-derived single effects; blue/green models receive measured individual function.'))
      print(fig)
      if (tt == 3 && aa == 'alpelisib' && bb == 'lapatinib') {
        ggsave(file.path(root, 'figures', 'alpelisib_lapatinib_curve_example.png'), fig, width = 11, height = 7, dpi = 150)
      }
    }
  }
  dev.off()
}

# Single-effect diagnostics distinguish calibration of known individual
# experiments from prediction of entirely excluded individual drugs.
single <- fread(file.path(root, 'data/single_predictions.tsv.gz'))
single_means <- single[, .(measured_pct = 100 * mean(ratio),
  fitted_pct = 100 * mean(2^predicted_fitted_log2),
  excluded_drug_pct = 100 * mean(2^predicted_excluded_drug_log2)),
  by = .(time_h, representation, source, drug, nominal_dose_uM)]
single_long <- melt(single_means, id.vars = c('time_h', 'representation', 'source', 'drug', 'nominal_dose_uM', 'measured_pct'),
  measure.vars = c('fitted_pct', 'excluded_drug_pct'), variable.name = 'stage', value.name = 'predicted_pct')
single_metrics <- single_long[, .(n = .N, r = safe_cor(measured_pct, predicted_pct),
  median_absolute_difference_pp = median(abs(measured_pct - predicted_pct))),
  by = .(time_h, representation, source, stage)]
save_table(single_metrics, 'single_decoder_summary')
source_agreement <- dcast(single_means[representation == 'Original40'],
  time_h + drug + nominal_dose_uM ~ source, value.var = 'measured_pct')
source_agreement <- source_agreement[complete.cases(source_agreement)]
source_agreement[, absolute_difference_pp := abs(`Individual drug panel` - `Same-screen single controls`)]
save_table(source_agreement, 'individual_source_agreement')
single_drug_summary <- single_means[source == 'Same-screen single controls',
  .(doses = .N, measured_min_pct = min(measured_pct), measured_max_pct = max(measured_pct),
    predicted_min_pct = min(fitted_pct), predicted_max_pct = max(fitted_pct),
    median_absolute_difference_pp = median(abs(measured_pct - fitted_pct))),
  by = .(time_h, representation, drug)]
save_table(single_drug_summary, 'single_decoder_by_drug')
for (tt in c(3, 24)) {
  z <- single_long[time_h == tt & representation == 'Original40']
  stat <- single_metrics[time_h == tt & representation == 'Original40']
  stat[, label := sprintf('r = %.2f; n = %d', r, n)]
  fig <- ggplot(z, aes(measured_pct, predicted_pct)) +
    geom_point(alpha = .35, size = 1.1, color = '#0072B2') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .55) +
    geom_text(data = stat, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.03, vjust = 1.2, size = 3.4) +
    facet_wrap(vars(source, stage), ncol = 2, scales = 'free') + theme(aspect.ratio = 1) +
    labs(title = paste('Original 40: can the single-response component recover individual effects?', tt, 'h'),
      subtitle = 'Fitted = available single experiments; excluded_drug = all individual measurements of that drug withheld',
      x = 'Measured individual cell count (% vehicle)', y = 'RNA-predicted individual cell count (% vehicle)',
      caption = 'Combination outcomes never enter the single-response decoder. Sources have equal total fitting weight.')
  save_plot(fig, paste0('single_response_decoder_', tt, 'h'), 8.4, 8.6)

  selected_drugs <- unique(single_means[source == 'Same-screen single controls', drug])
  single_curves <- rbind(
    single_means[time_h == tt & representation == 'Original40' & drug %chin% selected_drugs,
      .(drug, nominal_dose_uM, value = measured_pct, series = source)],
    single_means[time_h == tt & representation == 'Original40' & source == 'Same-screen single controls',
      .(drug, nominal_dose_uM, value = fitted_pct, series = 'Original 40: fitted single effect')])
  fig <- ggplot(single_curves, aes(nominal_dose_uM, value, color = series)) +
    geom_line(linewidth = .65) + geom_point(size = 1.3) +
    facet_wrap(~drug, ncol = 4, scales = 'free_x') + scale_x_log10() +
    scale_color_manual(values = c('Individual drug panel' = '#AAAAAA',
      'Same-screen single controls' = '#111111', 'Original 40: fitted single effect' = '#0072B2')) +
    labs(title = paste('Where single-response calibration loses information:', tt, 'h RNA'),
      subtitle = 'All individual experiments are available; these are fitted single effects, before predicting any combination',
      x = 'Individual drug dose (µM)', y = 'Cell count (% of vehicle)', color = NULL,
      caption = 'Black: singles from the combination experiment. Grey: the separate individual-drug screen. Blue: molecular decoder.')
  save_plot(fig, paste0('single_effect_dose_curves_', tt, 'h'), 11, 5.9)
}

# Capacity/transfer is shown on the actual multiplicative departure, avoiding
# total-count dose correlations. Complete bridge departures are relative to f.
capacity <- condition[model %chin% c('Functional_control', 'Original40_interaction', 'GenePC12_interaction')]
capacity[, model_label := factor(model, levels = main_models, labels = model_labels[main_models])]
capacity_stats <- capacity[, .(r = safe_cor(measured_interaction_log2, predicted_interaction_log2),
  n = .N), by = .(time_h, model_label, split)]
save_table(capacity_stats, 'interaction_capacity_summary')
for (tt in c(3, 24)) {
  stat <- capacity_stats[time_h == tt]
  stat[, label := sprintf('r = %.2f; n = %d', r, n)]
  fig <- ggplot(capacity[time_h == tt], aes(measured_interaction_log2, predicted_interaction_log2)) +
    geom_point(aes(color = anchor), alpha = .5, size = 1.1) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
    geom_text(data = stat, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.03, vjust = 1.2, size = 3.2) +
    facet_wrap(vars(split, model_label), ncol = 3, scales = 'free') +
    scale_color_manual(values = anchor_colors) + theme(aspect.ratio = 1) +
    labs(title = paste('Learning the actual departure from individual effects:', tt, 'h'),
      subtitle = 'Training and primary whole-partner exclusions; original 40 and gene PCs are conditional on known singles',
      x = 'Measured log2(combination / single-product)', y = 'Predicted log2 modification', color = 'First drug',
      caption = 'Zero = independent product; -1 = half its cell count; +1 = twice its cell count. Biological and outer-fit means.')
  save_plot(fig, paste0('interaction_training_and_exclusion_', tt, 'h'), 11, 8)
}

# Display all40 inputs, with explicit partial-measurement labels and unchanged
# signed scores. Drug blocks keep adjacent doses; programs are clustered once.
bank <- qread(file.path(prior, 'representation_scores.qs'))
md <- copy(bank$metadata)[drug %chin% unique(c(p$anchor, p$partner))]
setorder(md, time_h, drug, dose_uM)
mm <- t(bank$original40[md$profile_id, , drop = FALSE])
row_order <- hclust(dist(mm), method = 'ward.D2')$order
mm <- mm[row_order, , drop = FALSE]
annotation <- bank$program_annotations
row_labels <- annotation$display_name[match(rownames(mm), annotation$program)]
row_labels <- substr(row_labels, 1, 61)
lim <- as.numeric(quantile(abs(mm), .98))
cairo_pdf(file.path(root, 'figures/original40_separate_response_inputs.pdf'), width = 12.5, height = 8.7)
for (tt in c(3, 24)) {
  meta <- md[time_h == tt]
  mat <- mm[, meta$profile_id, drop = FALSE]
  col_anno <- data.frame(Drug = meta$drug, row.names = colnames(mat))
  stopifnot(identical(rownames(col_anno), colnames(mat)))
  hm <- draw_heatmap(mat, col_anno_df = col_anno, col_anno_cols = 'Drug', scale_rows = FALSE,
    value_range = c(-lim, 0, lim), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'Weighted\nresponse', show_row_names = TRUE, show_column_names = TRUE,
    row_labels = row_labels, column_labels = paste(meta$dose_uM, 'µM'), column_names_rot = 70,
    column_split = factor(meta$drug, levels = unique(meta$drug)), cluster_column_slices = FALSE,
    row_names_gp = grid::gpar(fontsize = 7.5), column_names_gp = grid::gpar(fontsize = 6),
    title = paste('Original 40 separate responses at', tt, 'h | measured978-gene contribution; missing members unknown'))
  ComplexHeatmap::draw(hm)
}
dev.off()

selection <- fread(file.path(root, 'data/interaction_selection.tsv.gz'))[selected == TRUE]
save_table(selection[, .(outer_fits = .N, selected_zero_interaction = sum(is.infinite(lambda)),
  median_selected_lambda = median(lambda)), by = .(time_h, fold_type, model)], 'selected_model_summary')
coefficients <- fread(file.path(root, 'data/interaction_coefficients.tsv.gz'))
coefficient_summary <- coefficients[, .(mean_raw_coefficient = mean(raw_coefficient),
  median_raw_coefficient = median(raw_coefficient), mean_absolute_scaled = mean(abs(scaled_coefficient)),
  positive_fits = sum(raw_coefficient > 0), negative_fits = sum(raw_coefficient < 0), zero_fits = sum(raw_coefficient == 0)),
  by = .(time_h, fold_type, model, term)]
setorder(coefficient_summary, time_h, fold_type, model, -mean_absolute_scaled)
save_table(coefficient_summary, 'interaction_coefficient_consistency')

report <- c('WORKFLOW62: separate responses and joint functional modification',
  'Completed retrospective test. See memo72 for the design frozen before fitting.',
  'Primary3h; secondary24h. Six whole-partner exclusions are primary.',
  'Fifteen two-partner exclusions additionally support same-model matched comparisons.',
  'All12 drug pairs retained in BT20. These are not held-out cell lines or studies.',
  'Known individual functional experiments train f; no combination trains f.',
  'Functional_control: measured singles offset + function/dose interaction.',
  'Original40_interaction: same inputs + symmetric original40 program interactions.',
  'Original40_bridge: RNA f(xA)+f(xB)+g(xA,xB), with no single function at inference.',
  'GenePC12 analogues diagnose compression; they are not an exhaustive gene-information test.',
  'Baseline conditioning is fixed because this dataset contains one cell line.',
  '', 'COMMON-MODEL PARTNER CONTRASTS (dose contexts are dependent):',
  capture.output(print(metrics)), '', 'SEPARATE FIRST-DRUG RESULTS:',
  capture.output(print(anchor_metrics)), '', 'SIMILAR INDIVIDUAL EXPECTATION (<=5pp difference):',
  capture.output(print(matched_metrics)), '', 'SINGLE-RESPONSE DECODER:',
  capture.output(print(single_metrics)), '',
  'Clear effects mean >=10pp and the same direction in every available biological experiment.',
  'This labels observed reproducibility, not a prospective signal selection or significance test.',
  'All small/inconsistent effects remain in tables, plots and all-context correlations.',
  'Raw count or derived extra-effect correlations alone cannot prove added RNA information.',
  'For molecular f+g, derived extra effects share the measured single-product reference with observations.',
  'Its apparently high extra-effect correlation must be checked against actual combination-count contrasts.',
  'The matched control comparison and common-model partner contrasts are the central tests.',
  'Previously inspected outcomes remain retrospective; no environmental-rescue claim is established.',
  'Original bank projection measures only the observed978-gene contribution; programs may be incomplete.',
  '', 'PRIMARY WHOLE-PARTNER EXCLUSIONS, INCLUDING THE PREVIOUS FUNCTIONAL SVR:',
  capture.output(print(primary_summary)))
writeLines(report, file.path(root, 'REPORT.txt'))
findings <- c('WORKFLOW62: WHAT THE COMPLETED TEST SHOWS', '',
  'The agreed architecture was implemented and evaluated. It does not yet provide convincing evidence for the grant bridge.',
  'Original40 is the primary frozen bank; 3h is primary and24h is a temporal sensitivity analysis.',
  'Training includes other combinations. Each test partner has no combination outcomes in training.',
  'The separate individual RNA and functional experiments remain available, as intended.', '',
  '1. Main test: six complete partner exclusions. The molecular f+g bridge predicts',
  'single effects from RNA, then adds a learned symmetric program interaction.',
  'Conditional versions receive measured single effects to isolate interaction learning.',
  'The previous functional-only SVR is retained as a strong comparator, not replaced by a weak baseline.', '',
  '2. Original40 fits the log interaction target in training, but transfers poorly:',
  capture.output(print(capacity_stats[time_h == 3 & model_label == model_labels[['Original40_interaction']]])), '',
  '3. Primary prediction of partner/dose differences at the SAME first-drug dose:',
  capture.output(print(primary_summary[time_h == 3 & model %chin%
    c('Previous_functional_SVR', 'Original40_interaction', 'Original40_bridge'),
    .(model, anchor, within_anchor_dose_count_r, median_absolute_difference_pp)])), '',
  '4. The original40 single decoder already distorts some known individual effects:',
  capture.output(print(single_drug_summary[time_h == 3 & representation == 'Original40' &
    drug %chin% c('alpelisib', 'lapatinib', 'nvp-tae684')])),
  'The source screens themselves differ despite mapping to the same RNA profiles:',
  capture.output(print(source_agreement[, .(matched_conditions = .N,
    median_absolute_source_difference_pp = median(absolute_difference_pp),
    maximum_absolute_source_difference_pp = max(absolute_difference_pp)), by = time_h])), '',
  '5. A stricter diagnostic excludes TWO partners together, so the same fitted model',
  'predicts both sides of each matched contrast. It finds modest added molecular',
  'information in some alpelisib comparisons, but not a consistent win across anchors.',
  'Only six partners and two first drugs constrain the interaction map. Training fit',
  'does not establish that these statistical program interactions transfer to other drugs.',
  'The gene-PC alternative also fails to establish the complete bridge; it is not an exhaustive gene test.', '',
  'READ THE FIGURES IN THIS ORDER:',
  'figures/all12_excluded_combination_curves_3h.pdf: actual joint outcomes and all comparators, every pair.',
  'figures/interaction_training_and_exclusion_3h.pdf: where interaction fitting loses transfer.',
  'figures/single_effect_dose_curves_3h.pdf: where the individual-effect decoder loses amplitude and ordering.',
  'figures/partner_modification_3h.pdf: stricter, same-model differences between two excluded partners.',
  'figures/original40_separate_response_inputs.pdf: complete original40 molecular inputs and biological annotations.', '',
  'These are retrospective BT20 results. They do not establish HGF rescue or cross-cell transfer.',
  'The functional-only model remains stronger overall. The central molecular-composition claim is not established by this run.')
writeLines(findings, file.path(root, 'READ_ME.txt'))
figures <- list.files(file.path(root, 'figures'), pattern = '\\.pdf$', full.names = TRUE)
fwrite(data.table(path = figures), file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
writeLines('Models, validation, figures and report complete.', file.path(root, 'STATUS.txt'))
artifacts <- c(list.files(root, recursive = TRUE, full.names = TRUE),
  '/work/src/program_pilot/62_separate_response_interaction_bridge.R',
  '/work/src/program_pilot/62_interaction_bridge_figures.R',
  '/work/memos/72_separate_response_interaction_bridge.txt')
artifacts <- artifacts[!file.info(artifacts)$isdir & basename(artifacts) != 'artifact_manifest.tsv']
fwrite(data.table(path = artifacts, bytes = file.info(artifacts)$size,
  sha256 = vapply(artifacts, digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'artifact_manifest.tsv'), sep = '\t')
print(metrics)
print(anchor_metrics)
