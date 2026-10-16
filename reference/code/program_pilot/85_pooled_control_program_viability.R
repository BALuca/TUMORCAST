#!/usr/bin/env Rscript
# Reference replacement diagnostic. Preserve treated-profile normalization,
# frozen program weights, PRISM responses, and all 867 previously plotted rows.
# A pooled control is the arithmetic mean of raw DMSO profiles from all fourteen
# plates for that same line, normalized using its full-transcriptome total.
# The old per-profile trimmed scaling remains fixed to isolate reference effects.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(digest)
})
setDTthreads(1L)
root <- '/work/results/program_pilot/85_pooled_control_program_viability'
prior <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
archive_path <- paste0('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/',
  'inputs/original/modular_pair_pilot/inputs.rds')
reference_path <- '/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge/artifacts/reference_genes.rds'
points_path <- file.path(prior, 'data/observed_program_viability_points.tsv.gz')
for (folder in c('data', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
programs <- sprintf('P%02d', 1:40)
doses <- c(.05, .5, 5)
colors <- c('0.05' = '#1879a4', '0.5' = '#df9b32', '5' = '#ae3769')
pearson <- function(x, y) {
  if (length(x) < 3L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y)
}

archive <- readRDS(archive_path)
old <- fread(points_path)
profiles <- unique(old[, .(depmap_id, drug, dose, contributing_plates)])
stopifnot(nrow(profiles) == 867L, nrow(old) == 867L * 40L,
  !anyDuplicated(old[, .(depmap_id, drug, dose, program)]))
W <- Matrix(as.matrix(archive$prep$W), sparse = TRUE)
stopifnot(identical(colnames(W), programs))
references <- readRDS(reference_path)
cm <- as.data.table(references$cm)
pm <- as.data.table(references$pcm)
genes <- rownames(W)
stopifnot(all(genes %in% colnames(references$C)), all(genes %in% colnames(references$PC)),
  !anyDuplicated(cm[, .(depmap_id, plate)]), !anyDuplicated(pm[, .(depmap_id, plate)]))

# Public and cached plate controls overlap exactly, so all fourteen cached plates
# can be used consistently for every plotted recipient without mixing estimates.
public_index <- match(paste(pm$depmap_id, pm$plate), paste(cm$depmap_id, cm$plate))
common <- which(!is.na(public_index))
reference_error <- max(abs(references$PC[common, genes, drop = FALSE] -
  references$C[public_index[common], genes, drop = FALSE]))
stopifnot(reference_error < 1e-8, all(profiles$depmap_id %in% cm$depmap_id))
selected <- which(cm$depmap_id %in% profiles$depmap_id)
C <- references$C[selected, genes, drop = FALSE]
cm <- cm[selected]
stopifnot(all(cm$total > 0), all(cm[, .N, by = depmap_id]$N == 14L),
  all(cm[, .(complete = setequal(plate, 1:14)), by = depmap_id]$complete))
ids <- sort(unique(cm$depmap_id))
pooled_log <- t(vapply(ids, function(id) {
  ii <- which(cm$depmap_id == id)
  log2(1 + colSums(C[ii, , drop = FALSE]) / sum(cm$total[ii]) * 10000)
}, numeric(length(genes))))
rownames(pooled_log) <- ids
plate_log <- log2(1 + C / cm$total * 10000)
control_scores <- as.matrix((plate_log - pooled_log[match(cm$depmap_id, ids), ]) %*% W)
control_wide <- cbind(cm[, .(depmap_id, plate)], as.data.table(control_scores))
controls <- melt(control_wide, id.vars = c('depmap_id', 'plate'), variable.name = 'program',
  value.name = 'dmso_minus_pooled', variable.factor = FALSE)

# D_new = D_old + log(control_on_treatment_plates) - log(pooled_control).
# This identity changes only the reference, including when a condition averages
# more than one plate. Averaging is done on raw counts before log transformation.
adjustments <- matrix(NA_real_, nrow(profiles), 40L, dimnames = list(NULL, programs))
for (i in seq_len(nrow(profiles))) {
  plate_ids <- as.integer(trimws(strsplit(profiles$contributing_plates[i], ',', fixed = TRUE)[[1]]))
  ii <- which(cm$depmap_id == profiles$depmap_id[i] & cm$plate %in% plate_ids)
  stopifnot(length(ii) == length(plate_ids))
  matched_log <- log2(1 + colSums(C[ii, , drop = FALSE]) / sum(cm$total[ii]) * 10000)
  adjustments[i, ] <- as.numeric((matched_log - pooled_log[profiles$depmap_id[i], ]) %*% W)
}
adjustment_long <- melt(cbind(profiles, as.data.table(adjustments)),
  id.vars = names(profiles), variable.name = 'program', value.name = 'reference_adjustment',
  variable.factor = FALSE)
pooled <- merge(old, adjustment_long[, .(depmap_id, drug, dose, program, reference_adjustment)],
  by = c('depmap_id', 'drug', 'dose', 'program'))
stopifnot(nrow(pooled) == nrow(old), all(is.finite(pooled$reference_adjustment)))
setnames(pooled, 'observed_displacement', 'matched_control_displacement')
pooled[, observed_displacement := matched_control_displacement + reference_adjustment]
pooled[, dose_label := factor(as.character(dose), levels = names(colors))]
old[, reference := 'Matched plate control']
pooled[, reference := 'Pooled 14-plate control']
comparison <- rbindlist(list(old[, .(depmap_id, drug, dose, program, observed_displacement,
  observed_log2_abundance, reference)], pooled[, .(depmap_id, drug, dose, program,
  observed_displacement, observed_log2_abundance, reference)]))
comparison[, `:=`(centered_displacement = observed_displacement - mean(observed_displacement),
  centered_response = observed_log2_abundance - mean(observed_log2_abundance)),
  by = .(reference, drug, dose, program)]
metrics <- comparison[, .(cells = uniqueN(depmap_id), points = .N,
  pooled_r = pearson(observed_displacement, observed_log2_abundance),
  within_dose_r = pearson(centered_displacement, centered_response)), by = .(reference, drug, program)]
dose_metrics <- comparison[, .(cells = .N, r = pearson(observed_displacement, observed_log2_abundance),
  mean_displacement = mean(observed_displacement), sd_displacement = sd(observed_displacement)),
  by = .(reference, drug, dose, program)]
control_metrics <- controls[, .(cells = .N, mean_control_shift = mean(dmso_minus_pooled),
  sd_control_shift = sd(dmso_minus_pooled)), by = .(plate, program)]
fwrite(pooled, file.path(root, 'data/pooled_control_scatter_points.tsv.gz'), sep = '\t')
fwrite(controls, file.path(root, 'data/dmso_plate_vs_pooled_programs.tsv.gz'), sep = '\t')
fwrite(metrics, file.path(root, 'data/matched_vs_pooled_correlations.tsv'), sep = '\t')
fwrite(dose_metrics, file.path(root, 'data/matched_vs_pooled_each_dose.tsv'), sep = '\t')
fwrite(control_metrics, file.path(root, 'data/dmso_plate_program_summary.tsv'), sep = '\t')
fwrite(cm[, .(depmap_id, plate, total)], file.path(root, 'data/pooled_control_membership.tsv'), sep = '\t')

# Reproduce the same forty program panels for each drug, in their frozen order.
drugs <- c('Lapatinib ditosylate', setdiff(unique(pooled$drug), 'Lapatinib ditosylate'))
pages <- lap_pages <- list()
for (drug_value in drugs) {
  for (page in 1:4) {
    pp <- programs[(page - 1L) * 10L + 1:10]
    rows <- pooled[drug == drug_value & program %in% pp]
    rows[, program := factor(program, levels = pp)]
    labels <- metrics[reference == 'Pooled 14-plate control' & drug == drug_value & program %in% pp]
    labels[, program := factor(program, levels = pp)]
    labels[, label := sprintf('Pooled r = %.2f\nWithin-dose r = %.2f', pooled_r, within_dose_r)]
    p <- ggplot(rows, aes(observed_displacement, observed_log2_abundance)) +
      geom_point(aes(color = dose_label), size = .8, alpha = .65) +
      geom_smooth(aes(color = dose_label), method = 'lm', formula = y ~ x, se = FALSE, linewidth = .45) +
      geom_label(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -.02, vjust = 1.03, size = 2.4, linewidth = 0, fill = '#FFFFFFCC') +
      facet_wrap(~program, ncol = 5, scales = 'free_x') +
      scale_color_manual(values = colors, name = 'Measured RNA dose (µM)', drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(.04, .16))) +
      theme_bw(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom',
        panel.grid.minor = element_blank()) +
      labs(x = 'Measured program displacement relative to pooled DMSO',
        y = 'Observed PRISM log2 abundance / vehicle',
        title = paste0(drug_value, ': pooled-control program displacements'),
        subtitle = 'Same measured cell–dose profiles and PRISM responses; only the subtracted RNA reference changes.',
        caption = paste0('Each line uses its own DMSO control averaged across all 14 plates. ',
          'Treated-profile scaling remains fixed. No fitting or filtering.\n',
          'PRISM responses are interpolated to Tahoe doses. 0 = 100%, −1 = 50%, −2 = 25% of vehicle abundance.'))
    pages[[paste(drug_value, page)]] <- p
    if (drug_value == 'Lapatinib ditosylate') {
      lap_pages[[page]] <- p
      ggsave(file.path(root, 'figures', paste0('lapatinib_pooled_control_page_', page, '.png')),
        p, width = 12.5, height = 6.5, dpi = 150)
    }
  }
}
pdf(file.path(root, 'figures/all_ten_drugs_all_40_programs_pooled_control.pdf'),
  width = 12.5, height = 6.5, useDingbats = FALSE)
for (p in pages) print(p)
dev.off()
pdf(file.path(root, 'figures/lapatinib_all_40_programs_pooled_control.pdf'),
  width = 12.5, height = 6.5, useDingbats = FALSE)
for (p in lap_pages) print(p)
dev.off()

# Actual DMSO control shifts, independent of drug-response measurements.
pdf(file.path(root, 'figures/dmso_controls_all_40_programs.pdf'),
  width = 12.5, height = 6.5, useDingbats = FALSE)
for (page in 1:4) {
  pp <- programs[(page - 1L) * 10L + 1:10]
  rows <- controls[program %in% pp]
  rows[, program := factor(program, levels = pp)]
  rows[, plate_group := factor(ifelse(plate <= 3L, as.character(plate), 'Other'),
    levels = c('1', '2', '3', 'Other'))]
  p <- ggplot(rows, aes(factor(plate), dmso_minus_pooled)) +
    geom_hline(yintercept = 0, color = '#888888', linewidth = .3) +
    geom_boxplot(aes(fill = plate_group), outlier.shape = NA, linewidth = .3) +
    geom_point(position = position_jitter(width = .15, height = 0, seed = 85), size = .4, alpha = .35) +
    facet_wrap(~program, ncol = 5, scales = 'free_y') +
    scale_fill_manual(values = c('1' = unname(colors[1]), '2' = unname(colors[2]),
      '3' = unname(colors[3]), 'Other' = '#CCCCCC')) +
    theme_bw(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'none',
      panel.grid.minor = element_blank(), axis.text.x = element_text(size = 7)) +
    labs(x = 'DMSO plate', y = 'Control program shift relative to same-line pooled control',
      title = 'DMSO controls: plate differences in the same 40 molecular programs',
      subtitle = '30 matched cell lines, 14 plates each; no treated cells or PRISM responses enter this calculation.',
      caption = 'Plates 1, 2 and 3 are highlighted because they provide the lapatinib dose references. Each point is one cell line.')
  print(p)
  if (page %in% c(2L, 4L)) ggsave(file.path(root, 'figures', paste0('dmso_controls_page_', page, '.png')),
    p, width = 12.5, height = 6.5, dpi = 150)
}
dev.off()

# Input identities and the exact reference-only interpretation remain reviewable.
input_paths <- c(archive_path, reference_path, points_path,
  '/work/src/program_pilot/85_pooled_control_program_viability.R')
fwrite(data.table(path = input_paths,
  sha256 = vapply(input_paths, digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'data/session.txt'))
writeLines(c('Reference-only replacement: treated RNA and original trimmed scaling are fixed.',
  'Pooled reference: raw DMSO expression averaged across all 14 plates of the same cell line, then normalized to CP10K.',
  'All 867 original cell-dose profiles and all 40 frozen programs retained; no model refitting.',
  paste('Cached/public overlapping control profiles:', length(common), '; maximum difference:', reference_error)),
  file.path(root, 'data/normalization_definition.txt'))
print(dose_metrics[drug == 'Lapatinib ditosylate' & program %in% c('P11', 'P16', 'P35', 'P40')])
print(metrics[drug == 'Lapatinib ditosylate' & program %in% c('P11', 'P16', 'P35', 'P40')])
cat('Completed pooled-control scatters and independent DMSO plate audit.\n')
