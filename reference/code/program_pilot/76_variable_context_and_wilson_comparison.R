#!/usr/bin/env Rscript
# Unfiltered diagnostic of saved predictions; no coefficients are refitted.
# All 1300 combination conditions are retained, as requested by the user.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
setDTthreads(1L)

# Measure variation of whole pair/dose/laboratory contexts without filtering.
prior <- '/work/results/program_pilot/75_hgf_pattern_partner_transfer'
root <- '/work/results/program_pilot/76_variable_context_and_wilson_comparison'
minimum_cells <- 25L
wilson_dose_uM <- 1
script <- '/work/src/program_pilot/76_variable_context_and_wilson_comparison.R'
inputs <- c(file.path(prior, 'data/held_predictions_evaluated.tsv'),
  file.path(prior, 'data/hgf_predictions_frozen.tsv'),
  '/work/almanac_combination_training_2026-10-03/results/wilson_observed_source.tsv')
stopifnot(all(file.exists(inputs)))
for (folder in c('data', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
}
# Undefined correlations are retained when a saved model predicts a constant.
pearson <- function(x, y) {
  if (length(x) >= 4L && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
}

held <- fread(inputs[1])
stopifnot(all(c('row_id', 'family', 'excluded_partner', 'SCREENER', 'doseA', 'doseB',
  'I', 'predicted_g', 'observed_departure_pp', 'predicted_departure_pp') %in% names(held)),
  !anyDuplicated(held[, .(row_id, family)]),
  all(is.finite(held$observed_departure_pp)), all(is.finite(held$predicted_g)))
held[, lapatinib_uM := fifelse(drugA == 'Lapatinib ditosylate', doseA, doseB)]
held[, partner_uM := fifelse(drugA == 'Lapatinib ditosylate', doseB, doseA)]
context_keys <- c('excluded_partner', 'SCREENER', 'lapatinib_uM', 'partner_uM')
observed <- held[family == 'Program interactions']
stopifnot(nrow(observed) == 1300L,
  !anyDuplicated(observed[, c(context_keys, 'depmap_id'), with = FALSE]))

# Measure variation of departures, rather than variation of total drug killing.
# Every cell in every context is retained, including its near-zero effects.
contexts <- observed[, .(cells = uniqueN(depmap_id),
  observed_departure_sd_pp = sd(observed_departure_pp),
  observed_interaction_sd_log2 = sd(I),
  observed_departure_min_pp = min(observed_departure_pp),
  observed_departure_max_pp = max(observed_departure_pp),
  single_agent_below_10pct = sum(pmin(vA, vB) < .1),
  repeated_combination_conditions = sum(replicates > 1)), by = context_keys]
contexts[, retained := TRUE]
contexts[, reason := 'All contexts restored; no variation cutoff']
setorderv(contexts, context_keys)
write_table(contexts, 'context_variation_and_selection')
held <- merge(held, contexts[, c(context_keys, 'retained', 'observed_departure_sd_pp'), with = FALSE],
  by = context_keys, all.x = TRUE, sort = FALSE)
stopifnot(nrow(held) == 2600L, !anyNA(held$retained))
write_table(held[retained == TRUE], 'retained_frozen_predictions')

# The primary comparison is across cells receiving the same pair and doses.
within <- held[, .(cells = uniqueN(depmap_id),
  pearson_interaction_log2 = pearson(I, predicted_g),
  pearson_departure_pp = pearson(observed_departure_pp, predicted_departure_pp)),
  by = c('family', context_keys, 'retained', 'observed_departure_sd_pp')]
write_table(within, 'fixed_treatment_correlations')
within_summary <- within[cells >= minimum_cells, .(contexts = .N,
  finite_log2_correlations = sum(is.finite(pearson_interaction_log2)),
  median_pearson_log2 = if (any(is.finite(pearson_interaction_log2))) median(pearson_interaction_log2, na.rm = TRUE) else NA_real_,
  finite_pp_correlations = sum(is.finite(pearson_departure_pp)),
  median_pearson_pp = if (any(is.finite(pearson_departure_pp))) median(pearson_departure_pp, na.rm = TRUE) else NA_real_),
  by = .(family, retained)]
write_table(within_summary, 'fixed_treatment_summary')
pooled <- held[, .(conditions = .N, cells = uniqueN(depmap_id), partners = uniqueN(excluded_partner),
  pearson_interaction_log2 = pearson(I, predicted_g),
  pearson_departure_pp = pearson(observed_departure_pp, predicted_departure_pp)), by = .(family, retained)]
write_table(pooled, 'pooled_context_correlations')

# Show the complete panel grouped by partner, including negative associations.
plot_rows <- held[retained == TRUE & family == 'Program interactions']
plot_metrics <- plot_rows[, .(cells = uniqueN(depmap_id), conditions = .N,
  pearson_interaction_log2 = pearson(I, predicted_g)), by = excluded_partner]
if (nrow(plot_rows)) {
  plot_metrics[, panel := paste0(sub(' \\(.*', '', excluded_partner), '\nn=', conditions,
    '; r=', fifelse(is.finite(pearson_interaction_log2), sprintf('%.2f', pearson_interaction_log2), 'undefined'))]
  plot_rows <- merge(plot_rows, plot_metrics[, .(excluded_partner, panel)], by = 'excluded_partner')
  p <- ggplot(plot_rows, aes(I, predicted_g)) +
    geom_point(size = 1.3, alpha = .7, color = '#2166AC') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
    facet_wrap(~panel, ncol = min(3L, nrow(plot_metrics)), scales = 'free') +
    theme_classic(base_size = 9) +
    theme(aspect.ratio = 1, strip.background = element_blank(), strip.text = element_text(size = 8),
      plot.caption = element_text(hjust = 0, size = 8)) +
    labs(x = 'Observed departure from separate effects (log2)', y = 'Predicted departure (log2)',
      title = 'Lapatinib combinations: all 1,300 conditions restored',
      caption = 'No variation filtering; all partners, doses and cells shown. Each partner excluded from interaction training.\nFrozen version75 predictions. Retrospective development evaluation conditional on the molecular shortlist.')
  plot_width <- 2.7 * min(3L, nrow(plot_metrics)) + .6
  plot_height <- 2.9 * ceiling(nrow(plot_metrics) / 3) + 1.25
} else {
  p <- ggplot() + annotate('text', x = 0, y = 0,
    label = 'No combination conditions available.') + theme_void()
  plot_width <- 6.5
  plot_height <- 2
}
ggsave(file.path(root, 'figures/variable_context_predictions.pdf'), p,
  width = plot_width, height = plot_height, device = cairo_pdf)
ggsave(file.path(root, 'figures/variable_context_predictions.png'), p,
  width = plot_width, height = plot_height, dpi = 160)

# Wilson's published abundance is relative to untreated control. Use the saved
# absolute combination extrapolation here; conditional HGF-only normalization
# cannot be compared directly to Wilson's untreated-normalized measurements.
hgf <- fread(inputs[2])[family == 'Program interactions' & dose == wilson_dose_uM &
  hgf_context == 'Equal-context mean']
stopifnot(nrow(hgf) == 3L, !anyDuplicated(hgf$original_fold))
model <- rbind(hgf[, .(original_fold, condition = 'Lapatinib', abundance = lapatinib_alone_pct)],
  hgf[, .(original_fold, condition = 'Lapatinib + HGF', abundance = combined_extrapolated_pct)])
published <- fread(inputs[3])[series == 'Wilson et al.' & condition %in% model$condition]
stopifnot(nrow(published) == 2L, !anyDuplicated(published$condition))
write_table(model, 'wilson_dose_absolute_predictions')
write_table(published, 'wilson_observed_reference')
comparison <- model[, .(predicted_min = min(abundance), predicted_max = max(abundance)), by = condition]
comparison <- merge(comparison, published[, .(condition, observed_mean = abundance, observed_sd = sd)], by = 'condition')
write_table(comparison, 'wilson_dose_comparison')
condition_order <- c('Lapatinib', 'Lapatinib + HGF')
model[, condition := factor(condition, levels = condition_order)]
published[, condition := factor(condition, levels = condition_order)]
wp <- ggplot() +
  geom_hline(yintercept = 100, color = 'grey70', linetype = 2, linewidth = .35) +
  geom_line(data = model, aes(condition, abundance, group = factor(original_fold)),
    color = '#2166AC', alpha = .5, linewidth = .55) +
  geom_point(data = model, aes(condition, abundance, color = 'Current model'),
    size = 2, position = position_nudge(x = -.07)) +
  geom_errorbar(data = published, aes(condition, ymin = abundance - sd, ymax = abundance + sd),
    color = '#B2182B', width = .12, position = position_nudge(x = .07), linewidth = .65) +
  geom_point(data = published, aes(condition, abundance, color = 'Wilson observed'),
    shape = 18, size = 3, position = position_nudge(x = .07)) +
  scale_color_manual(values = c('Current model' = '#2166AC', 'Wilson observed' = '#B2182B')) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  theme_classic(base_size = 10) +
  theme(legend.position = 'bottom', plot.caption = element_text(hjust = 0, size = 8)) +
  labs(x = NULL, y = 'Abundance (% untreated control)', color = NULL,
    title = 'AU565 at Wilson’s lapatinib dose: 1 µM',
    caption = 'Blue points: three frozen B1 fits plus the same learned interaction modifier. Red: published mean ± SD.\nHGF uses the mean of 12 transported profiles; no matched HGF dose calibration.\nAbsolute HGF-alone effect and interaction transfer across assays remain extrapolations. Assay durations differ.')
ggsave(file.path(root, 'figures/wilson_1uM_comparison.pdf'), wp, width = 6.4, height = 4.6, device = cairo_pdf)
ggsave(file.path(root, 'figures/wilson_1uM_comparison.png'), wp, width = 6.4, height = 4.6, dpi = 160)

# Ensure that this diagnostic did not modify any prior input or predictions.
manifest <- data.table(path = c(inputs, script),
  sha256 = vapply(c(inputs, script), digest, character(1), file = TRUE, algo = 'sha256'))
write_table(manifest, 'source_manifest')
freeze <- fread(file.path(prior, 'data/prediction_freeze_hashes.tsv'))
stopifnot(identical(digest(inputs[2], file = TRUE, algo = 'sha256'),
  freeze[grepl('hgf_predictions_frozen.tsv$', path), sha256]))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
writeLines(c(
  'Unfiltered diagnostic of version75 frozen held-new-partner predictions. No refitting.',
  'User retracted the variation filter. All1300 conditions,34 fixed-dose contexts and9 partners restored.',
  'Within-treatment summaries use at least25 cells, as in version75; every current context satisfies this count.',
  'Withdrawn filtered outputs preserved in the sibling76_variable_context_and_wilson_comparison_withdrawn_filter folder.',
  'None of the1300 lapatinib combination conditions has repeated combination records in the prepared archive.',
  'Wilson reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC3724525/ Figure3c; lapatinib1uM, crizotinib0.5uM. Only the two modeled lapatinib/HGF conditions are compared.',
  'Current queries used1 and2.5uM lapatinib. The older selected Wilson reconstruction used5uM molecular descriptors and2.5uM functional calibration as exposure proxies.',
  'Wilson comparison uses absolute abundance relative to untreated control; conditional abundance relative to HGF-only control is not substituted.',
  'Wilson outcomes are displayed only after frozen predictions; no fitting or selection on their agreement.',
  capture.output(print(within_summary)), capture.output(print(pooled)), capture.output(print(comparison))
), file.path(root, 'RESULTS.txt'))
print(contexts[retained == TRUE])
print(within_summary)
print(pooled)
print(comparison)
