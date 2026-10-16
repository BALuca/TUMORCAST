#!/usr/bin/env Rscript
# Locate molecular information before program compression in the MEK examples.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/58_molecular_function_information_audit'
prior <- '/work/results/program_pilot/57_erbb_met_hgf_program_biology/data'
drugs <- c('trametinib', 'cobimetinib')
paths <- c(
  '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs',
  '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs',
  file.path(prior, 'drug_observations.tsv.gz'),
  file.path(prior, 'biological_readout_membership.tsv.gz'),
  '/work/results/program_pilot/56_drug_program_biology_reliability/data/mek_genotype_biology.tsv.gz')
stopifnot(all(file.exists(paths)))
for (folder in c('data', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0('mek_', name, '.tsv.gz')), sep = '\t')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0('mek_', name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0('mek_', name, '.png')), p, width = width, height = height, dpi = 150)
}
cor_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L || sd(x[ok]) < 1e-10 || sd(y[ok]) < 1e-10) return(NA_real_)
  cor(x[ok], y[ok])
}
panel <- qread(paths[1])
bank <- qread(paths[2])
md <- fread(paths[3])[drug %in% drugs & is.finite(observed_log2)]
genotype <- unique(fread(paths[5])[, .(depmap_id, PIK3CA, RAS_RAF_hotspot)])
stopifnot(!anyDuplicated(genotype$depmap_id), !anyDuplicated(md$observation_id))
md <- merge(md, genotype, by = 'depmap_id')
setorder(md, drug, dose, depmap_id)
members <- fread(paths[4])
selected_sets <- c('ERK_output', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'MTORC1_SIGNALING',
  'UNFOLDED_PROTEIN_RESPONSE', 'APOPTOSIS')
members <- members[readout %in% selected_sets & in_Tahoe_input == TRUE]

# The same measured gene universe is used for both sister compounds.
delta_parts <- baseline_parts <- null_parts <- list()
for (id in unique(md$depmap_id)) {
  obj <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
  ids <- md[depmap_id == id, observation_id]
  stopifnot(identical(obj$normalizer, 'paired_trimmed_median_log_ratio_cp10k_v1'),
    identical(colnames(obj$response), obj$genes), all(ids %in% rownames(obj$response)))
  delta_parts[[id]] <- obj$response[ids, , drop = FALSE]
  baseline_parts[[id]] <- matrix(obj$baseline[[id]], nrow = 1L, dimnames = list(id, obj$genes))
  null_parts[[id]] <- obj$null[[id]]
}
delta <- do.call(rbind, delta_parts)[md$observation_id, , drop = FALSE]
baseline <- do.call(rbind, baseline_parts)
stopifnot(identical(rownames(delta), md$observation_id), identical(colnames(delta), colnames(baseline)))
baseline_observations <- baseline[md$depmap_id, , drop = FALSE]
rownames(baseline_observations) <- md$observation_id
# A cell's mean vehicle baseline differs from each treatment's matched control.
# This proxy asks whether retaining baseline level restores information, without
# claiming to reconstruct the exact treated profile or kinase phosphorylation.
treated_proxy <- baseline_observations + delta
modalities <- list(Baseline = baseline_observations, Displacement = delta,
  `Baseline + displacement proxy` = treated_proxy)
save_table(md, 'observations')
save_table(members, 'fixed_readout_membership')
save_table(data.table(path = paths, sha256 = vapply(paths, function(p) digest(file = p, algo = 'sha256'), character(1))),
  'input_manifest')
message('Loaded ', nrow(md), ' molecular/function observations in ', uniqueN(md$depmap_id), ' cells.')

# Fixed biological sets, independent of the current outcome association.
readout_list <- list()
control_list <- list()
for (term_name in selected_sets) {
  genes <- members[readout == term_name, gene]
  stopifnot(length(genes) > 0L, all(genes %in% colnames(delta)))
  for (mode in names(modalities)) {
    readout_list[[paste(term_name, mode)]] <- data.table(observation_id = md$observation_id,
      readout = term_name, modality = mode, value = rowMeans(modalities[[mode]][, genes, drop = FALSE]))
  }
  for (id in names(null_parts)) {
    control_list[[paste(term_name, id)]] <- data.table(depmap_id = id, readout = term_name,
      vehicle_abs95 = quantile(abs(rowMeans(null_parts[[id]][, genes, drop = FALSE])), .95, names = FALSE))
  }
}
readouts <- merge(rbindlist(readout_list), md, by = 'observation_id')
readouts <- merge(readouts, rbindlist(control_list), by = c('depmap_id', 'readout'))
save_table(readouts, 'readout_states')
old <- fread(file.path(prior, 'drug_biological_readouts.tsv.gz'))[drug %in% drugs & readout %in% selected_sets]
check_scores <- merge(readouts[modality == 'Displacement', .(observation_id, readout, value)],
  old[, .(observation_id, readout, measured_genes)], by = c('observation_id', 'readout'))
stopifnot(nrow(check_scores) == nrow(readouts[modality == 'Displacement']),
  max(abs(check_scores$value - check_scores$measured_genes)) < 1e-10)

# Per-cell full-dose means preserve threshold transitions without treating doses
# as independent cells. These descriptive associations are not predictions.
cell_states <- readouts[, .(value = mean(value), observed_log2 = mean(observed_log2),
  abundance_geomean = 100 * 2^mean(observed_log2), doses = uniqueN(dose)),
  by = .(drug, depmap_id, recipient, PIK3CA, RAS_RAF_hotspot, readout, modality)]
save_table(cell_states, 'cell_curve_readouts')
associations <- cell_states[, {
  r <- cor_safe(value, observed_log2)
  lo <- vapply(seq_len(.N), function(j) cor_safe(value[-j], observed_log2[-j]), numeric(1))
  .(cells = .N, r = r, omit_one_min_r = min(lo, na.rm = TRUE), omit_one_max_r = max(lo, na.rm = TRUE))
}, by = .(drug, readout, modality)]
save_table(associations, 'readout_function_associations')
group_summary <- readouts[RAS_RAF_hotspot == TRUE, .(cells = uniqueN(depmap_id),
  median_value = median(value), median_abundance = median(measured_pct),
  displacement_above_vehicle95 = if (modality[1] == 'Displacement') mean(abs(value) > vehicle_abs95) else NA_real_),
  by = .(drug, dose, readout, modality, PIK3CA)]
save_table(group_summary, 'hotspot_readout_summary')
noise_summary <- readouts[RAS_RAF_hotspot == TRUE & modality == 'Displacement',
  .(cells = .N, above_vehicle95 = sum(abs(value) > vehicle_abs95),
    median_vehicle95 = median(vehicle_abs95), median_displacement = median(value)),
  by = .(drug, dose, readout)]
save_table(noise_summary, 'response_versus_vehicle_variation')
noise_rows <- readouts[RAS_RAF_hotspot == TRUE & modality == 'Displacement' & dose == .5]
stopifnot(all(noise_rows$vehicle_abs95 > 0))
noise_rows[, relative_to_vehicle := value / vehicle_abs95]
noise_rows[, readout := factor(readout, levels = selected_sets)]
np <- ggplot(noise_rows, aes(readout, relative_to_vehicle)) +
  geom_hline(yintercept = c(-1, 1), linetype = 2, color = 'grey50') +
  geom_hline(yintercept = 0, color = 'grey80') +
  geom_boxplot(outlier.shape = NA, width = .5, color = 'grey60') +
  geom_point(aes(color = factor(PIK3CA)), position = position_jitter(width = .13, height = 0, seed = 1), size = 1.7) +
  facet_wrap(~drug, nrow = 1) +
  scale_color_manual(values = c('#247F8A', '#B95720'), name = 'PIK3CA hotspot') +
  scale_x_discrete(labels = c('ERK', 'E2F', 'MYC', 'mTORC1', 'UPR', 'Apoptosis')) +
  labs(x = 'Fixed gene readout', y = 'Displacement / vehicle 95% magnitude',
    title = 'Clear MEK engagement; downstream changes mostly within vehicle variation',
    subtitle = '0.5 uM; 25 RAS/RAF-hotspot cells per compound; points retain cell-specific control variation.',
    caption = 'Dashed lines: each cell/readout 95th percentile of absolute vehicle pseudo-responses.\nA control reference, not a treated-replicate significance test or a measure of kinase activity.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(np, 'response_versus_vehicle_variation', 10, 5)

# Check all measured genes without choosing markers using outcome significance.
# Correlations have one observation per cell, averaging its entire measured curve.
gene_tests <- list()
for (drug_name in drugs) {
  z <- md[drug == drug_name]
  cell_ids <- sort(unique(z$depmap_id))
  counts <- table(z$depmap_id)
  stopifnot(length(unique(as.integer(counts))) == 1L)
  yy <- z[, .(y = mean(observed_log2)), by = depmap_id][match(cell_ids, depmap_id), y]
  for (mode in names(modalities)) {
    x <- modalities[[mode]][z$observation_id, , drop = FALSE]
    means <- rowsum(x, z$depmap_id)[cell_ids, , drop = FALSE] / as.numeric(counts[cell_ids])
    xc <- sweep(means, 2L, colMeans(means), '-')
    yc <- yy - mean(yy)
    denominator <- sqrt(colSums(xc^2) * sum(yc^2))
    rr <- as.numeric(crossprod(yc, xc)) / denominator
    rr[denominator < 1e-12] <- NA_real_
    tt <- rr * sqrt((length(yy) - 2) / pmax(1e-15, 1 - rr^2))
    pp <- 2 * pt(-abs(tt), df = length(yy) - 2)
    gene_tests[[paste(drug_name, mode)]] <- data.table(drug = drug_name, modality = mode,
      gene = colnames(x), cells = length(yy), r = rr, p_value = pp, q_value = p.adjust(pp, 'BH'))
  }
}
gene_tests <- rbindlist(gene_tests)
save_table(gene_tests, 'all_gene_function_associations')
sister <- merge(gene_tests[drug == drugs[1]], gene_tests[drug == drugs[2]], by = c('gene', 'modality'),
  suffixes = c('_trametinib', '_cobimetinib'))
sister[, same_direction := sign(r_trametinib) == sign(r_cobimetinib)]
save_table(sister, 'sister_drug_gene_associations')
gene_summary <- gene_tests[, .(tested = sum(is.finite(r)), q05 = sum(q_value <= .05, na.rm = TRUE),
  q10 = sum(q_value <= .10, na.rm = TRUE)), by = .(drug, modality)]
save_table(gene_summary, 'gene_signal_counts')
save_table(sister[, .(genes = .N, association_agreement_r = cor_safe(r_trametinib, r_cobimetinib),
  same_direction_q05_both = sum(q_value_trametinib <= .05 & q_value_cobimetinib <= .05 & same_direction, na.rm = TRUE)),
  by = modality], 'sister_drug_summary')
save_table(merge(sister, members[, .(gene, readout)], by = 'gene', allow.cartesian = TRUE), 'named_gene_associations')

# Plot every recipient in the hotspot comparison, preserving the sensitive exception.
plot_rows <- readouts[RAS_RAF_hotspot == TRUE & dose == .5 &
  readout %in% c('ERK_output', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'MTORC1_SIGNALING')]
plot_rows[, genotype := factor(PIK3CA, levels = 0:1, labels = c('No PIK3CA hotspot', 'PIK3CA hotspot'))]
plot_rows[, modality := factor(modality, levels = names(modalities))]
p <- ggplot(plot_rows, aes(genotype, value, color = genotype)) +
  geom_boxplot(outlier.shape = NA, color = 'grey60', width = .55) +
  geom_point(position = position_jitter(width = .12, height = 0, seed = 1), size = 1.4) +
  facet_wrap(vars(readout, modality), scales = 'free_y', ncol = 3) +
  scale_color_manual(values = c('#247F8A', '#B95720')) +
  labs(x = NULL, y = 'Mean gene expression score',
    title = 'MEK inhibition: starting level, displacement and retained-state proxy',
    subtitle = '0.5 uM; 25 RAS/RAF-hotspot cells; each cell contributes one point per compound.',
    caption = 'Trametinib and cobimetinib share cells. Baseline + delta uses a mean vehicle reference, not exact treated RNA.') +
  theme_classic(base_size = 9) + theme(legend.position = 'none', axis.text.x = element_text(size = 7))
# Separate compounds avoid displaying repeat cells as independent evidence.
for (drug_name in drugs) save_plot(p %+% plot_rows[drug == drug_name] + labs(subtitle = paste(drug_name, '| 0.5 uM | n=25 cells')),
  paste0(drug_name, '_baseline_response'), 10, 10)

sp <- cell_states[readout %in% c('ERK_output', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'MTORC1_SIGNALING')]
sp[, modality := factor(modality, levels = names(modalities))]
for (drug_name in drugs) {
  a <- associations[drug == drug_name & readout %in% unique(sp$readout)]
  p <- ggplot(sp[drug == drug_name], aes(value, observed_log2)) +
    geom_point(aes(color = factor(PIK3CA)), size = 1.5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey35', linewidth = .5) +
    geom_text(data = a, aes(x = -Inf, y = Inf, label = sprintf('r = %.2f', r)),
      inherit.aes = FALSE, hjust = -.1, vjust = 1.2, size = 3) +
    facet_wrap(vars(readout, modality), scales = 'free', ncol = 3) +
    scale_color_manual(values = c('#247F8A', '#B95720'), name = 'PIK3CA hotspot') +
    labs(x = 'Mean molecular score over measured doses', y = 'Mean log2 abundance / vehicle',
      title = paste(drug_name, '| molecular information across cells'),
      subtitle = 'One point per cell, all measured doses included. Descriptive association, not a fitted prediction.') +
    theme_classic(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom')
  save_plot(p, paste0(drug_name, '_molecular_function'), 10, 11)
}

# Fixed named genes expose heterogeneity hidden by broad pathway means.
marker_sets <- list(ERK = c('PHLDA1','SPRY2','DUSP6','ETV4','ETV5'),
  Growth = c('MYC','CCND1','CCND2','CDK2','CCNE1','CCNA2','MKI67','PCNA','MCM2','MCM3','MCM5','MCM6','CDC20','CCNB1','BIRC5'),
  Fate = c('CDKN1A','BCL2L1','MCL1','BCL2','BCL2L11','BBC3','DDIT3','ATF4','ASNS','XBP1','TRIB3'))
markers <- rbindlist(lapply(names(marker_sets), function(nm) data.table(gene = marker_sets[[nm]], role = nm)))
markers <- markers[gene %in% colnames(delta)]
save_table(markers, 'displayed_genes')
for (drug_name in drugs) {
  z <- md[drug == drug_name & dose == .5][order(measured_pct)]
  ids <- z$observation_id
  anno <- data.frame(PIK3CA = factor(z$PIK3CA), RAS_RAF = factor(z$RAS_RAF_hotspot), row.names = ids)
  ht <- NULL
  for (mode in names(modalities)) {
    mm <- t(modalities[[mode]][ids, markers$gene, drop = FALSE])
    stopifnot(identical(colnames(mm), ids))
    # Expression levels vary greatly between genes: baseline/proxy are row-centered;
    # displacement retains its absolute signed scale. No row variance scaling.
    if (mode != 'Displacement') mm <- sweep(mm, 1L, rowMeans(mm), '-')
    h <- draw_heatmap(mm, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = factor(markers$role, levels = names(marker_sets)), cluster_row_slices = FALSE,
      col_anno_df = if (mode == 'Baseline') anno else NULL,
      col_anno_cols = if (mode == 'Baseline') names(anno) else NULL,
      palette_overrides = list(PIK3CA = c('0'='#247F8A','1'='#B95720'),
        RAS_RAF = c('FALSE'='grey80','TRUE'='#705B95')),
      ann_top_args = list(`Abundance %` = ComplexHeatmap::anno_barplot(z$measured_pct,
        gp = grid::gpar(fill = '#666666'), border = FALSE)),
      show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7),
      show_column_names = TRUE, column_labels = z$recipient, column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = 6),
      value_range = c(-.6, 0, .6), value_colors = c('#2166AC','white','#B2182B'),
      heatmap_legend_param = list(at = c(-.6, 0, .6)),
      hm_name = mode, title = paste(mode, if (mode != 'Displacement') '(gene-centered)' else '(original delta)'))
    ht <- if (is.null(ht)) h else ht %v% h
  }
  name <- paste0('mek_', drug_name, '_genes_baseline_response')
  pdf(file.path(root, 'figures', paste0(name, '.pdf')), width = 13, height = 18)
  figure_title <- paste(drug_name, '0.5 uM | same cells ordered by PRISM abundance\n',
    'Top: baseline | middle: displacement | bottom: baseline + displacement proxy')
  ComplexHeatmap::draw(ht, column_title = figure_title, column_title_gp = grid::gpar(fontsize = 11),
    merge_legends = TRUE, padding = grid::unit(rep(6,4),'mm'))
  dev.off()
  png(file.path(root, 'figures', paste0(name, '.png')), width = 1560, height = 2160, res = 120)
  ComplexHeatmap::draw(ht, column_title = figure_title, column_title_gp = grid::gpar(fontsize = 11),
    merge_legends = TRUE, padding = grid::unit(rep(6,4),'mm'))
  dev.off()
}
writeLines(capture.output(sessionInfo()), file.path(root, 'mek_session_info.txt'))
save_table(data.table(script = '/work/src/program_pilot/58_molecular_function_information_audit.R',
  sha256 = digest(file = '/work/src/program_pilot/58_molecular_function_information_audit.R', algo = 'sha256')),
  'script_provenance')
print(associations)
print(gene_summary)
message('MEK information audit complete. No new program discovery or outcome predictor fitted.')
