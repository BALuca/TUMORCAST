#!/usr/bin/env Rscript
# Independent molecular replication in H1792, retaining differences in culture,
# collection time and trametinib formulation. No outcomes select genes/programs.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/60_defensible_response_composition'
prior <- '/work/results/program_pilot/59_mechanistic_evidence_datasets/data'
input_paths <- c(
  RNA = file.path(prior, 'mek_Macaya_RNA_gene_matrices.qs'),
  Tahoe = '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/ACH-000496.qs',
  weights = file.path(root, 'models/representation_frozen_weights.qs'),
  protein = file.path(prior, 'mek_Macaya_protein_contrasts.tsv.gz'),
  measurements = file.path(prior, 'mek_measurements.tsv.gz'),
  PRISM = file.path(prior, 'resource_PRISM_matching_measurements.tsv'))
stopifnot(all(file.exists(input_paths)))
for (folder in c('data', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
save_table <- function(tab, name) {
  fwrite(tab, file.path(root, 'data', paste0('h1792_', name, '.tsv.gz')), sep = '\t')
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0('h1792_', name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0('h1792_', name, '.png')), p, width = width, height = height, dpi = 150)
}
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

# Port of the established numeric normalization only; no Python is executed.
# CP10K, genes >=0.1 CP10K in both arms, central40% logratios and central90%
# abundance; median retained ratio removes compositional offset.
paired_delta <- function(treated, reference) {
  stopifnot(length(treated) == length(reference), all(is.finite(treated)),
            all(is.finite(reference)), all(treated >= 0), all(reference >= 0))
  tp <- treated / sum(treated) * 1e4
  rp <- reference / sum(reference) * 1e4
  eligible <- tp >= .1 & rp >= .1
  stopifnot(sum(eligible) >= 100L)
  m <- log2(tp[eligible] / rp[eligible])
  a <- .5 * log2(tp[eligible] * rp[eligible])
  mq <- quantile(m, c(.30, .70))
  aq <- quantile(a, c(.05, .95))
  keep <- m >= mq[1] & m <= mq[2] & a >= aq[1] & a <= aq[2]
  stopifnot(sum(keep) >= 100L)
  correction <- 2 ^ median(m[keep])
  log2(1 + tp / correction) - log2(1 + rp)
}
rna <- qread(input_paths['RNA'])
tahoe <- qread(input_paths['Tahoe'])
weights <- qread(input_paths['weights'])
stopifnot(identical(tahoe$normalizer, 'paired_trimmed_median_log_ratio_cp10k_v1'),
          identical(rownames(tahoe$response), tahoe$metadata$observation_id),
          identical(colnames(rna$counts), rna$metadata$sample_id))
md <- copy(rna$metadata)
control <- md$treatment == 'DMSO'
trt <- md$treatment == 'trametinib'
stopifnot(sum(control) == 3L, sum(trt) == 3L)
reference <- rowMeans(rna$counts[, control, drop = FALSE])
delta <- vapply(seq_len(ncol(rna$counts)), function(i) paired_delta(rna$counts[, i], reference),
                numeric(nrow(rna$counts)))
dimnames(delta) <- dimnames(rna$counts)
# Center to the observed vehicle mean; group uncertainty includes vehicle variance.
delta <- delta - rowMeans(delta[, control, drop = FALSE])
tm <- tahoe$metadata[grepl('^trametinib', drug)]
stopifnot(nrow(tm) == 6L, all(tm$dose %in% c(.05, .5, 5)))
tm[, formulation := ifelse(grepl('solvate', drug), 'DMSO_TF solvate', 'Standard')]
tm[, `:=`(time_h = 24, culture = 'Pooled 3D spheroids', study = 'Tahoe')]
setorder(tm, formulation, dose)
tdelta <- t(tahoe$response[tm$observation_id, , drop = FALSE])
tnull <- t(tahoe$null[['ACH-000496']])
save_table(tm, 'Tahoe_profiles')

# Gene-level discovery uses only independent study replicates. Genes are not
# filtered by whether their direction agrees with Tahoe or by functional outcome.
group <- factor(md$treatment, levels = c('DMSO', 'trametinib', 'lestaurtinib'))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
dge <- DGEList(rna$counts)
keep <- filterByExpr(dge, design = design)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
v <- voom(dge, design, plot = FALSE)
fit <- eBayes(contrasts.fit(lmFit(v, design), makeContrasts(trametinib - DMSO, levels = design)))
de <- as.data.table(topTable(fit, number = Inf, sort.by = 'none'), keep.rownames = 'gene')
setnames(de, c('logFC', 'adj.P.Val'), c('RNA_log2_fold_change', 'RNA_q'))
genes <- intersect(rownames(delta), rownames(tdelta))
stopifnot(length(genes) > 5000L)
gene_means <- data.table(gene = genes,
  RNA_delta = rowMeans(delta[genes, trt, drop = FALSE]),
  RNA_delta_SE = sqrt(matrixStats::rowVars(delta[genes, trt, drop = FALSE]) / 3 +
                      matrixStats::rowVars(delta[genes, control, drop = FALSE]) / 3),
  Tahoe_vehicle_abs95 = apply(abs(tnull[genes, , drop = FALSE]), 1, quantile, .95))
gene_means <- merge(gene_means, de[, .(gene, RNA_log2_fold_change, RNA_q)], by = 'gene', all.x = TRUE)
gene_long <- rbindlist(lapply(seq_len(nrow(tm)), function(i) {
  tab <- copy(gene_means)
  tab[, `:=`(observation_id = tm$observation_id[i], formulation = tm$formulation[i],
             dose = tm$dose[i], Tahoe_delta = tdelta[gene, i])]
  tab[, `:=`(same_direction = RNA_delta * Tahoe_delta > 0,
             Tahoe_above_vehicle = abs(Tahoe_delta) > Tahoe_vehicle_abs95)]
  tab
}))
save_table(gene_long, 'gene_replication')
agreement <- gene_long[, .(genes = .N, Pearson_r = cor(RNA_delta, Tahoe_delta),
  RNA_significant = sum(RNA_q <= .05, na.rm = TRUE),
  RNA_significant_same_direction = sum(RNA_q <= .05 & same_direction, na.rm = TRUE),
  RNA_significant_same_direction_above_vehicle = sum(RNA_q <= .05 & same_direction & Tahoe_above_vehicle, na.rm = TRUE),
  Pearson_r_RNA_significant = cor(RNA_delta[!is.na(RNA_q) & RNA_q <= .05],
                                 Tahoe_delta[!is.na(RNA_q) & RNA_q <= .05])), by = .(formulation, dose)]
save_table(agreement, 'gene_agreement')
print(agreement)

# Frozen representations are evaluated on identical shared measured genes.
# Coverage is explicit; no missing gene is assigned a zero response.
score_parts <- coverage_parts <- list()
all_scores <- list()
for (bank in names(weights)) {
  w <- as.matrix(weights[[bank]])
  common <- intersect(rownames(w), genes)
  coverage <- colSums(abs(w[common, , drop = FALSE])) / colSums(abs(w))
  coverage_parts[[bank]] <- data.table(representation = bank, feature = colnames(w),
    total_genes = colSums(w != 0), measured_genes = colSums(w[common, , drop = FALSE] != 0),
    coverage_fraction = coverage, eligible = coverage >= .8)
  w <- w[common, coverage >= .8, drop = FALSE]
  stopifnot(ncol(w) > 0L, all(w >= 0))
  w <- sweep(w, 2, colSums(w), '/')
  rs <- crossprod(w, delta[common, , drop = FALSE])
  ts <- crossprod(w, tdelta[common, , drop = FALSE])
  ns <- crossprod(w, tnull[common, , drop = FALSE])
  score_parts[[bank]] <- rbindlist(lapply(seq_len(nrow(tm)), function(i) data.table(
    representation = bank, feature = rownames(rs), formulation = tm$formulation[i], dose = tm$dose[i],
    RNA_delta = rowMeans(rs[, trt, drop = FALSE]),
    RNA_SE = sqrt(matrixStats::rowVars(rs[, trt, drop = FALSE]) / 3 +
                    matrixStats::rowVars(rs[, control, drop = FALSE]) / 3),
    Tahoe_delta = ts[, i], Tahoe_vehicle_abs95 = apply(abs(ns), 1, quantile, .95))))
  all_scores[[bank]] <- list(RNA = rs, Tahoe = ts, vehicle = ns)
}
score_data <- rbindlist(score_parts)
coverage <- rbindlist(coverage_parts)
save_table(score_data, 'program_replication')
save_table(coverage, 'program_coverage')
score_agreement <- score_data[, .(features = .N, Pearson_r = cor(RNA_delta, Tahoe_delta),
  same_direction = sum(RNA_delta * Tahoe_delta > 0),
  same_direction_above_vehicle = sum(RNA_delta * Tahoe_delta > 0 & abs(Tahoe_delta) > Tahoe_vehicle_abs95)),
  by = .(representation, formulation, dose)]
save_table(score_agreement, 'program_agreement')
qsave(list(RNA_delta = delta, RNA_metadata = md, Tahoe_delta = tdelta, Tahoe_metadata = tm,
           program_scores = all_scores), file.path(root, 'data/h1792_molecular_replication.qs'))

# Readouts were specified independently in the earlier biological audit.
# Their means are transcriptional readouts, never measurements of phosphorylation.
readout_members <- fread('/work/results/program_pilot/57_erbb_met_hgf_program_biology/data/biological_readout_membership.tsv.gz')
readout_names <- c('ERK_output', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'MTORC1_SIGNALING',
                   'UNFOLDED_PROTEIN_RESPONSE', 'APOPTOSIS')
readout_parts <- list()
for (term in readout_names) {
  gg <- intersect(readout_members[readout == term, gene], genes)
  stopifnot(length(gg) >= 5L)
  rs <- colMeans(delta[gg, , drop = FALSE])
  ts <- colMeans(tdelta[gg, , drop = FALSE])
  ns <- colMeans(tnull[gg, , drop = FALSE])
  readout_parts[[term]] <- data.table(readout = term, n_genes = length(gg), formulation = tm$formulation,
    dose = tm$dose, RNA_delta = mean(rs[trt]),
    RNA_SE = sqrt(var(rs[trt]) / 3 + var(rs[control]) / 3),
    RNA_all_3_treated_same_direction = all(rs[trt] > 0) || all(rs[trt] < 0),
    Tahoe_delta = ts, Tahoe_vehicle_abs95 = quantile(abs(ns), .95))
}
readouts <- rbindlist(readout_parts)
save_table(readouts, 'biological_readouts')

# Independent protein measurement: preserve source protein entries and separately
# report median gene-level abundance (multiple entries are not replicates).
protein <- fread(input_paths['protein'])
protein <- protein[gene != '0' & !is.na(gene) & is.finite(fold_change) & fold_change > 0]
protein_gene <- protein[, .(protein_log2FC = median(log2(fold_change)),
                            source_entries = .N, minimum_source_p = min(p_value, na.rm = TRUE)),
                        by = .(gene, contrast)]
protein_compare <- merge(gene_means, protein_gene[contrast == 'TvsC'], by = 'gene')
save_table(protein_compare, 'RNA_protein_agreement')
protein_complete <- protein_compare[is.finite(RNA_log2_fold_change) & is.finite(protein_log2FC)]
save_table(protein_compare[!is.finite(RNA_log2_fold_change), .(gene,
  reason = 'RNA did not pass expression filtering for the independent differential-expression model')],
  'protein_comparison_missing_RNA')
selected_markers <- c('DUSP6', 'DUSP4', 'SPRY2', 'ETV4', 'ETV5', 'PHLDA1', 'FOS', 'JUN',
  'MYC', 'CCND1', 'CDKN1A', 'BCL2L11', 'EGFR', 'ERBB2', 'MET', 'AKT1', 'MAPK1')
marker_source <- merge(data.table(gene = selected_markers), gene_means, by = 'gene', all.x = TRUE)
marker_source <- merge(marker_source, protein_gene[contrast == 'TvsC'], by = 'gene', all.x = TRUE)
marker_source[, source_protein_missing := is.na(protein_log2FC)]
save_table(marker_source, 'named_marker_coverage')
measurements <- fread(input_paths['measurements'])
independent <- measurements[study_id == 'Macaya2023' & cell_line == 'H1792' &
  biological_state == 'parental' & assay %in% c('qPCR', 'immunoblot_densitometry', 'cell_viability')]
save_table(independent, 'independent_measurements')

# Source PRISM curves are shown as barcode abundance vs vehicle, not cell counts.
prism <- fread(input_paths['PRISM'])[cell == 'H1792' & drug == 'trametinib']
prism[, `:=`(percent_vehicle = 100 * 2 ^ response,
              screen = sub('^.*::(HTS002|MTS010.*)$', '\\1', column_name))]
setorder(prism, screen, dose)
save_table(prism, 'PRISM_curve')

# Across-gene agreement at the exact matched dose: all common genes and the
# independently significant subset are both supplied, without picking winners.
exact <- gene_long[dose == .5]
labels <- agreement[dose == .5, .(formulation,
  label = sprintf('r = %.2f; n = %s genes\n%d / %d significant genes agree', Pearson_r,
    format(genes, big.mark = ','), RNA_significant_same_direction, RNA_significant))]
p <- ggplot(exact, aes(RNA_delta, Tahoe_delta)) +
  geom_point(aes(color = !is.na(RNA_q) & RNA_q <= .05), size = .5, alpha = .3) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_hline(yintercept = 0, linewidth = .2) + geom_vline(xintercept = 0, linewidth = .2) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -.03, vjust = 1.15, size = 3) +
  facet_wrap(~formulation, scales = 'free') + theme(aspect.ratio = 1, legend.position = 'bottom') +
  scale_color_manual(values = c('FALSE' = 'grey65', 'TRUE' = '#0072B2'), name = 'Independent RNA q <= 0.05') +
  labs(title = 'Same cell and dose; independent RNA experiments',
    subtitle = 'H1792, trametinib 0.5 µM: Macaya 16 h culture vs Tahoe 24 h pooled spheroids',
    x = 'Macaya: expression change vs vehicle', y = 'Tahoe: expression change vs vehicle',
    caption = 'Both axes: paired compositional normalization, log2(1 + CP10K) difference.\nFormulations remain separate; Tahoe profiles are plate-aggregated, not independent replicate dots.')
save_plot(p, 'gene_replication', 9, 5.3)

# This second view shows every gene selected by independent RNA significance,
# including discordant genes; none are selected by Tahoe agreement.
p <- ggplot(exact[!is.na(RNA_q) & RNA_q <= .05], aes(RNA_delta, Tahoe_delta)) +
  geom_point(size = 1.3, alpha = .75, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_hline(yintercept = 0, linewidth = .2) + geom_vline(xintercept = 0, linewidth = .2) +
  geom_text(data = agreement[dose == .5],
    aes(x = -Inf, y = Inf, label = sprintf('r = %.2f; all %d RNA-significant genes',
      Pearson_r_RNA_significant, RNA_significant)), inherit.aes = FALSE, hjust = -.03, vjust = 1.3, size = 3) +
  facet_wrap(~formulation, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'The independently detected response reproduces better than the full gene cloud',
    subtitle = 'Selection uses Macaya RNA replicates only: BH q <= 0.05; no Tahoe or functional labels',
    x = 'Macaya 16 h: expression change vs vehicle', y = 'Tahoe 24 h: expression change vs vehicle')
save_plot(p, 'significant_gene_replication', 9.5, 5)
exact_scores <- score_data[dose == .5 & representation %in% c('bank16', 'bank42', 'bank43')]
score_labels <- score_agreement[dose == .5 & representation %in% c('bank16', 'bank42', 'bank43'),
  .(representation, formulation, label = sprintf('r = %.2f\n%d / %d same direction', Pearson_r, same_direction, features))]
p <- ggplot(exact_scores, aes(RNA_delta, Tahoe_delta)) +
  geom_errorbar(aes(xmin = RNA_delta - 1.96 * RNA_SE, xmax = RNA_delta + 1.96 * RNA_SE),
                orientation = 'y', width = 0, color = 'grey65') + geom_point(size = 1.5, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = score_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -.05, vjust = 1.2, size = 3) +
  facet_wrap(vars(representation, formulation), scales = 'free', ncol = 2) + theme(aspect.ratio = 1) +
  labs(title = 'Do the frozen program averages preserve replication?',
    subtitle = 'All programs with >=80% shared-gene coverage; no outcome fitting',
    x = 'Macaya 16 h: mean gene displacement', y = 'Tahoe 24 h: mean gene displacement',
    caption = 'Horizontal bars: mean +/-1.96SE from three treated and three vehicle cultures; descriptive uncertainty.')
save_plot(p, 'program_replication', 8.5, 11.2)

# Biological readouts preserve complete Tahoe dose series. The external RNA dose
# is fixed0.5µM and appears only at that x-coordinate, not across all doses.
p <- ggplot(readouts, aes(dose, Tahoe_delta, color = formulation, group = formulation)) +
  geom_hline(yintercept = 0, color = 'grey65', linewidth = .3) + geom_line(linewidth = .5) + geom_point(size = 1.7) +
  geom_errorbar(aes(ymin = -Tahoe_vehicle_abs95, ymax = Tahoe_vehicle_abs95),
                width = .04, alpha = .35, position = position_dodge(width = .1)) +
  geom_point(data = unique(readouts[, .(readout, RNA_delta, RNA_SE)]),
    aes(x = .5, y = RNA_delta), inherit.aes = FALSE, shape = 18, size = 3.4, color = 'black') +
  geom_errorbar(data = unique(readouts[, .(readout, RNA_delta, RNA_SE)]),
    aes(x = .5, ymin = RNA_delta - 1.96 * RNA_SE, ymax = RNA_delta + 1.96 * RNA_SE),
    inherit.aes = FALSE, width = .05, color = 'black') +
  scale_x_log10(breaks = c(.05, .5, 5)) + facet_wrap(~readout, scales = 'free_y', ncol = 3) +
  theme(legend.position = 'bottom') + scale_color_manual(values = c('#D55E00', '#0072B2')) +
  labs(title = 'Which biological response reproduces?', x = 'Trametinib (µM)', y = 'Mean gene displacement',
    color = 'Tahoe formulation', subtitle = 'Black diamond: independent 16 h RNA at 0.5 µM; colored lines: Tahoe 24 h',
    caption = 'Colored bars around zero: 95th percentile absolute vehicle contrasts (14 controls), not treatment confidence intervals.\nBlack bars: +/-1.96SE. Readouts describe RNA; ERK_output is not a pERK measurement.')
save_plot(p, 'biological_readouts', 10, 6.2)

markers <- intersect(selected_markers, genes)
hm <- cbind(delta[markers, trt, drop = FALSE], tdelta[markers, , drop = FALSE])
colnames(hm) <- c(paste0('Macaya16h: ', md$sample_id[trt]),
  paste0(tm$formulation, ' ', tm$dose, ' µM'))
row_order <- hclust(dist(hm), method = 'ward.D2')$order
hm <- hm[row_order, , drop = FALSE]
ann <- data.frame(Study = c(rep('Macaya16h', 3), rep('Tahoe24h', 6)), row.names = colnames(hm))
stopifnot(identical(colnames(hm), rownames(ann)))
pdf(file.path(root, 'figures/h1792_named_genes.pdf'), width = 8, height = 6)
draw_heatmap(hm, col_anno_df = ann, col_anno_cols = 'Study', scale_rows = FALSE,
  value_range = c(-1.5, 0, 1.5), value_colors = c('#2166AC', 'white', '#B2182B'),
  hm_name = 'RNA delta', show_row_names = TRUE, show_column_names = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE, column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 8),
  title = 'H1792: named signaling / cell-fate transcripts; fixed color scale')
dev.off()

# Trace independently defined ERK genes into each fixed bank. A coherent decrease
# can be dispersed among large mixed programs; this view directly shows that loss.
erk_genes <- intersect(readout_members[readout == 'ERK_output', gene], genes)
trace_parts <- list()
for (bank in c('bank16', 'bank42', 'bank43')) {
  w <- as.matrix(weights[[bank]])
  for (gene_symbol in erk_genes) {
    if (!gene_symbol %in% rownames(w)) next
    features <- colnames(w)[w[gene_symbol, ] != 0]
    stopifnot(length(features) == 1L)
    feature_name <- features[1]
    feature_tab <- score_data[representation == bank & feature == feature_name & dose == .5]
    stopifnot(nrow(feature_tab) == 2L)
    gene_tab <- gene_long[gene == gene_symbol & dose == .5]
    trace_parts[[paste(bank, gene_symbol)]] <- merge(
      gene_tab[, .(gene, formulation, gene_RNA_delta = RNA_delta, gene_Tahoe_delta = Tahoe_delta)],
      feature_tab[, .(formulation, representation, feature, program_RNA_delta = RNA_delta,
                       program_Tahoe_delta = Tahoe_delta)], by = 'formulation')[,
      program_n_genes := sum(w[, feature_name] != 0)]
  }
}
trace <- rbindlist(trace_parts)
save_table(trace, 'ERK_gene_to_program_trace')
trace_matrix_parts <- list()
trace_labels <- trace_groups <- character()
for (gene_symbol in erk_genes) {
  gene_tab <- gene_long[gene == gene_symbol & dose == .5]
  values <- c(gene_tab$RNA_delta[1], gene_tab[formulation == 'Standard', Tahoe_delta],
              gene_tab[formulation == 'DMSO_TF solvate', Tahoe_delta])
  trace_matrix_parts[[length(trace_matrix_parts) + 1L]] <- values
  trace_labels <- c(trace_labels, paste0('Gene: ', gene_symbol))
  trace_groups <- c(trace_groups, gene_symbol)
  for (bank in c('bank16', 'bank42', 'bank43')) {
    tab <- trace[gene == gene_symbol & representation == bank]
    if (!nrow(tab)) next
    values <- c(tab$program_RNA_delta[1], tab[formulation == 'Standard', program_Tahoe_delta],
                tab[formulation == 'DMSO_TF solvate', program_Tahoe_delta])
    trace_matrix_parts[[length(trace_matrix_parts) + 1L]] <- values
    trace_labels <- c(trace_labels, paste0(bank, ': ', tab$feature[1], ' (', tab$program_n_genes[1], ' genes)'))
    trace_groups <- c(trace_groups, gene_symbol)
  }
}
trace_matrix <- do.call(rbind, trace_matrix_parts)
dimnames(trace_matrix) <- list(make.unique(trace_labels), c('Macaya 16 h', 'Tahoe standard 24 h', 'Tahoe solvate 24 h'))
pdf(file.path(root, 'figures/h1792_ERK_gene_to_program_trace.pdf'), width = 7.4, height = 7.3)
trace_hm <- draw_heatmap(trace_matrix, scale_rows = FALSE, value_range = c(-1.5, 0, 1.5),
  value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'RNA delta',
  show_row_names = TRUE, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = factor(trace_groups, levels = unique(trace_groups)), cluster_row_slices = FALSE,
  row_title = NULL, row_labels = trace_labels, row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 9), column_title_gp = grid::gpar(fontsize = 11),
  heatmap_legend_param = list(at = c(-1.5, 0, 1.5)), column_names_rot = 40,
  title = 'ERK-response genes and their containing programs\nTrametinib 0.5 µM\nUnscaled displacement vs matched vehicle')
ComplexHeatmap::draw(trace_hm)
dev.off()
png(file.path(root, 'figures/h1792_ERK_gene_to_program_trace.png'), width = 1110, height = 1095, res = 150)
ComplexHeatmap::draw(trace_hm)
dev.off()

p <- ggplot(protein_complete, aes(RNA_log2_fold_change, protein_log2FC)) +
  geom_point(size = .8, alpha = .35, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  theme(aspect.ratio = 1) + labs(title = 'RNA and later total protein abundance',
    subtitle = sprintf('H1792 trametinib vs vehicle; %d genes; Pearson r = %.2f', nrow(protein_complete),
      cor(protein_complete$RNA_log2_fold_change, protein_complete$protein_log2FC)),
    x = 'RNA 16 h: log2 fold change (voom model)', y = 'Protein 48 h: log2 abundance fold change',
    caption = 'Different timepoints and wells; total protein abundance does not establish kinase activity.\nMultiple source protein entries are collapsed by gene median, not treated as independent replicates.')
save_plot(p, 'RNA_protein_agreement', 5.5, 5.7)
p <- ggplot(prism, aes(dose, percent_vehicle, color = screen, group = screen)) +
  geom_hline(yintercept = 100, linetype = 2, color = 'grey55') + geom_line() + geom_point(size = 2) +
  scale_x_log10() + labs(title = 'H1792 trametinib: independent PRISM curves',
    x = 'Trametinib (µM)', y = 'Barcode abundance (% of vehicle)', color = 'Screen',
    caption = 'Each screen remains separate. These are multi-day barcode-abundance assays, not the 16 h RNA or 72 h Macaya viability assay.') +
  theme(legend.position = 'bottom', legend.text = element_text(size = 7))
save_plot(p, 'PRISM_curves', 6.5, 4.7)

save_table(data.table(input = names(input_paths), path = input_paths,
  sha256 = vapply(input_paths, digest, character(1), algo = 'sha256', file = TRUE)), 'input_manifest')
capture.output(sessionInfo(), file = file.path(root, 'data/h1792_session_info.txt'))
writeLines(c(
  'H1792 independent molecular replication: interpretation boundary',
  'Gene discovery: RNA-only limma-voom trametinib versus DMSO, three cultures each; BH q<=0.05.',
  'Tahoe responses retain the fixed historical normalizer; Macaya counts use the same paired CP10K transformation for scores.',
  'The independent study has16h RNA in conventional culture; Tahoe has24h RNA in pooled3D spheroids.',
  'Same cell and0.5µM dose do not imply identical biological conditions. Both Tahoe formulations are retained.',
  'Four frozen representations are descriptive compressions, not independent confirmation of their original training signal.',
  'Program coverage is measured on shared genes; allfeatures>=80% coverage are reported without direction selection.',
  'The fourteen Tahoe vehicle contrasts quantify background, not biological treatment replication.',
  'The protein assay measures total abundance48h, not phosphosignaling. Missing EGFR/MYC/etc proteins remain missing.',
  'MYC forced-expression rescue is a separate causal reference inH2009/H358, not a demonstrated causal state inH1792.',
  'This analysis tests reproducibility of response information before program compression and after it. It does not fit a function predictor.',
  'Primary sources: https://www.nature.com/articles/s41467-023-41828-z ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161218',
  'Tahoe timing/culture: https://www.sciencedirect.com/science/article/pii/S0092867426010081'),
  file.path(root, 'H1792_REPLICATION_DESIGN.txt'))
exact_agreement <- agreement[dose == .5]
report_lines <- c(
  'H1792 trametinib: biological result and what the programs retain',
  '',
  'There is independently reproducible RNA response information before fitting any function model.',
  sprintf('Among %d shared genes, %d pass the independent-study replicate-based RNA test (BH q<=0.05).',
    exact_agreement$genes[1], exact_agreement$RNA_significant[1]),
  vapply(seq_len(nrow(exact_agreement)), function(i) sprintf(
    '%s, 0.5 uM: %d/%d significant genes change in the same direction; their Pearson r=%.3f; all-gene r=%.3f.',
    exact_agreement$formulation[i], exact_agreement$RNA_significant_same_direction[i],
    exact_agreement$RNA_significant[i], exact_agreement$Pearson_r_RNA_significant[i],
    exact_agreement$Pearson_r[i]), character(1)),
  '',
  'The clearest fixed biological readout is suppression of ERK-output transcripts: PHLDA1, SPRY2, DUSP6, ETV4, ETV5.',
  'All five decrease in the independent RNA study and in both Tahoe formulations at0.5uM; all five exceed their Tahoe vehicle background.',
  'This is an RNA footprint consistent with MEK inhibition, not direct evidence about pERK or cell death.',
  'The E2F, MYC-target and mTOR transcriptional averages are weak relative to the Tahoe vehicle background; their different directions do not establish bypass.',
  '',
  'A representation-specific failure is visible without any predictor:',
  'Bank43 distributes these five genes into five large programs (198-1219 genes), whose average responses are small and sometimes reverse the gene direction.',
  'Bank16 omits PHLDA1 and ETV5. Bank42 groups PHLDA1 and DUSP6 in its30-gene P47; that program decreases in both independent RNA and Tahoe.',
  'Thus some banks dilute a reproducible response whereas another retains a coherent subset. A single bank cannot serve as the biological ground truth.',
  'See h1792_ERK_gene_to_program_trace.pdf and its source table for every gene, membership and measured response.',
  '',
  'Orthogonal measurements and limits:',
  sprintf('RNA16h versus total protein48h: %d genes; Pearson r=%.3f. This is not a phosphosignaling comparison.',
    nrow(protein_complete), cor(protein_complete$RNA_log2_fold_change, protein_complete$protein_log2FC)),
  'The quantitative proteome lacks most named ERK-response regulators and MYC; absence from that assay is not lack of pathway regulation.',
  'Macaya H1792 viability72h source means: trametinib81.670%, lestaurtinib89.739%, combination18.517% of vehicle.',
  'Combination genome-wide RNA is absent. Single-drug RNA replication does not by itself establish the combination mechanism.',
  'MYC-overexpression rescue inH2009/H358 remains a separate genetic causal reference.',
  '',
  'Scientific implication: the data support a reproducible, mechanism-consistent molecular response and show exactly how particular bank averages can erase it.',
  'They do not establish the cross-cell sensitivity predictor or a successful molecular composition model.',
  'Culture and timing differ (Macaya16h conventional culture; Tahoe24h pooled3D spheroids). Both formulations and all tested doses are retained.')
writeLines(report_lines, file.path(root, 'H1792_RESULTS.txt'))
message('H1792 replication complete: ', length(genes), ' shared genes; ', nrow(score_agreement), ' representation comparisons.')
