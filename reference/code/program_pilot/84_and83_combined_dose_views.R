#!/usr/bin/env Rscript
# Descriptive dose pooling of already-frozen held-cell predictions. No fitting,
# selection, response filtering, or modification of earlier predictions occurs.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(digest)})
setDTthreads(1L)
versions <- c('83_lincs_and_tahoe_joint_program_bridge', '84_full_panel_lincs_tahoe_joint_bridge')
doses <- c(.15625, .625, 2.5)
anchor <- 'Lapatinib ditosylate'
correlation <- function(x, y) {
  if (length(x) < 4L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y)
}
lapatinib_points <- lapatinib_metrics <- combination_points <- combination_metrics <- list()
for (version in versions) {
  root <- file.path('/work/results/program_pilot', version)
  source_paths <- file.path(root, 'data', c('held_single_predictions_frozen.tsv', 'held_combination_predictions_frozen.tsv'))
  frozen <- fread(file.path(root, 'data/prediction_freeze_hashes.tsv'))
  hash <- vapply(source_paths, digest, character(1), file = TRUE, algo = 'sha256')
  stopifnot(identical(unname(hash), frozen$sha256[match(source_paths, frozen$path)]))
  single <- fread(source_paths[1L])
  stopifnot(!anyDuplicated(single[, .(family, evaluation_id)]),
    all(is.finite(single$observed_log2)), all(is.finite(single$predicted_log2)),
    max(abs(sort(unique(single$dose)) - doses)) < 1e-9)
  single[, `:=`(observed_centered = observed_log2 - mean(observed_log2),
    predicted_centered = predicted_log2 - mean(predicted_log2)), by = .(family, drug, dose)]
  metrics <- single[, .(cells = uniqueN(depmap_id), cell_dose_points = .N,
    pooled_r = correlation(observed_log2, predicted_log2),
    within_dose_r = correlation(observed_centered, predicted_centered)), by = .(family, drug)]
  fwrite(metrics, file.path(root, 'data/combined_dose_single_metrics.tsv'), sep = '\t')
  lap <- single[drug == anchor]
  stopifnot(all(lap[, .N, by = .(family, depmap_id)]$N == 3L),
    all(lap[, .(cells = uniqueN(depmap_id)), by = family]$cells == 441L))
  averaged <- lap[, .(observed = mean(observed_log2), predicted = mean(predicted_log2)),
    by = .(family, depmap_id)]
  average_metric <- averaged[, .(dose_averaged_cell_r = correlation(observed, predicted)), by = family]
  lap_metric <- merge(metrics[drug == anchor], average_metric, by = 'family')
  lap[, version_label := substr(version, 1L, 2L)]
  lap_metric[, version_label := substr(version, 1L, 2L)]
  lapatinib_points[[version]] <- lap
  lapatinib_metrics[[version]] <- lap_metric
  fwrite(lap_metric, file.path(root, 'data/combined_dose_lapatinib_metrics.tsv'), sep = '\t')

  # Single-drug panels retain all available cell/dose rows. Some drugs have
  # missing released endpoints; they are not imputed or newly filtered here.
  figures <- list()
  for (family_value in unique(single$family)) {
    rows <- single[family == family_value]
    labels <- metrics[family == family_value,
      .(drug, label = sprintf('pooled r = %.3f\nwithin-dose r = %.3f', pooled_r, within_dose_r))]
    p <- ggplot(rows, aes(observed_log2, predicted_log2)) +
      geom_point(aes(color = factor(dose)), size = .5, alpha = .4) +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .45, color = 'black') +
      geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -.05, vjust = 1.15, size = 2.4) +
      facet_wrap(~drug, ncol = if (uniqueN(rows$drug) == 10L) 5 else 3, scales = 'free') +
      scale_color_manual(values = c('#1879a4', '#df9b32', '#ae3769'), name = 'Dose (µM)') +
      theme_bw(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom',
        strip.text = element_text(size = 8)) +
      labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
        title = paste('v', substr(version, 1L, 2L), ': ', family_value, '; measured doses pooled', sep = ''))
    figures[[family_value]] <- p
  }
  width <- if (uniqueN(single$drug) == 10L) 12 else 9
  height <- if (uniqueN(single$drug) == 10L) 6 else 9
  pdf(file.path(root, 'figures/combined_dose_all_drugs.pdf'), width = width, height = height)
  for (p in figures) print(p)
  dev.off()

  combination <- fread(source_paths[2L])[regime == 'recipient']
  stopifnot(!anyDuplicated(combination[, .(family, row_id)]))
  combination[, `:=`(observed_centered = observed_log2 - mean(observed_log2),
    predicted_centered = predicted_log2 - mean(predicted_log2),
    interaction_observed_centered = observed_departure_log2 - mean(observed_departure_log2),
    interaction_predicted_centered = predicted_departure_log2 - mean(predicted_departure_log2)),
    by = .(family, SCREENER, pair, doseA, doseB)]
  combination[, `:=`(observed_pair_centered = observed_log2 - mean(observed_log2),
    predicted_pair_centered = predicted_log2 - mean(predicted_log2)), by = .(family, SCREENER, pair)]
  cm <- combination[, .(cells = uniqueN(depmap_id), pairs = uniqueN(pair), cell_pair_dose_points = .N,
    pooled_total_r = correlation(observed_log2, predicted_log2),
    pair_centered_total_r = correlation(observed_pair_centered, predicted_pair_centered),
    within_pair_dose_total_r = correlation(observed_centered, predicted_centered),
    pooled_departure_r = correlation(observed_departure_log2, predicted_departure_log2),
    within_pair_dose_departure_r = correlation(interaction_observed_centered, interaction_predicted_centered)), by = family]
  per_pair <- combination[, .(cells = uniqueN(depmap_id), cell_dose_points = .N,
    pooled_total_r = correlation(observed_log2, predicted_log2),
    within_dose_total_r = correlation(observed_centered, predicted_centered)), by = .(family, SCREENER, pair)]
  fwrite(cm, file.path(root, 'data/combined_dose_combination_summary.tsv'), sep = '\t')
  fwrite(per_pair, file.path(root, 'data/combined_dose_per_pair_metrics.tsv'), sep = '\t')
  combination[, version_label := substr(version, 1L, 2L)]
  cm[, version_label := substr(version, 1L, 2L)]
  combination_points[[version]] <- combination
  combination_metrics[[version]] <- cm
  fwrite(data.table(path = source_paths, sha256 = unname(hash)),
    file.path(root, 'data/combined_dose_view_inputs.tsv'), sep = '\t')
}

root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
lap <- rbindlist(lapatinib_points)
lm <- rbindlist(lapatinib_metrics)
lap[, panel := paste0('v', version_label, ': ', family)]
lm[, panel := paste0('v', version_label, ': ', family)]
lm[, label := sprintf('pooled r = %.3f\nwithin-dose r = %.3f', pooled_r, within_dose_r)]
p <- ggplot(lap, aes(observed_log2, predicted_log2)) +
  geom_point(aes(color = factor(dose)), size = .65, alpha = .45) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = 'black') +
  geom_text(data = lm, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.15, size = 2.6) + facet_wrap(~panel, ncol = 4, scales = 'free') +
  scale_color_manual(values = c('#1879a4', '#df9b32', '#ae3769'), name = 'Dose (µM)') +
  theme_bw(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom', strip.text = element_text(size = 8)) +
  labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
    title = 'Lapatinib: three measured doses pooled',
    caption = 'Each panel: 1,323 cell–dose points from 441 held cells. Within-dose r removes each dose mean from both axes.\nv83: nine-drug training panel; v84: original ten-drug panel. No predictions were refitted or recalibrated.')
ggsave(file.path(root, 'figures/combined_dose_lapatinib_v83_v84.pdf'), p, width = 11.5, height = 3.8)
ggsave(file.path(root, 'figures/combined_dose_lapatinib_v83_v84.png'), p, width = 11.5, height = 3.8, dpi = 150)
fwrite(lm, file.path(root, 'data/combined_dose_lapatinib_v83_v84.tsv'), sep = '\t')
fwrite(rbindlist(combination_metrics), file.path(root, 'data/combined_dose_combinations_v83_v84.tsv'), sep = '\t')
print(lm[, .(version_label, family, cells, cell_dose_points, pooled_r, within_dose_r, dose_averaged_cell_r)])
print(rbindlist(combination_metrics))
