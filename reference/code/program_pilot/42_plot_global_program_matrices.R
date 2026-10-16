#!/usr/bin/env Rscript
# Plot the frozen script-42 bank; no discovery, regression or clustering is fitted.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2)})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(8L)
output_dir <- '/work/results/program_pilot/42_global_pooled_signature_programs_ward'
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
bank_path <- file.path(output_dir, 'models/program_bank.qs')
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
# Adjacent columns are averaged only for the global display, at approximately
# one column per horizontal pixel. Bins never cross program or assay-match blocks.
display_columns <- 2200L
bank <- qread(bank_path)
stopifnot(ncol(bank$weights) == 60L, length(unique(bank$partition)) == 60L)
bank_hash <- digest::digest(bank, algo = 'sha256')
panel <- qread(panel_path)
md <- copy(bank$metadata)
sensitivity <- copy(panel$sensitivity)
responses <- list()
message('Loading saved molecular inputs for ', length(bank$genes), ' genes')
for (id in panel$cell_map$depmap_id) {
  molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
  responses[[id]] <- molecular$response[md[depmap_id == id, observation_id], bank$genes, drop = FALSE]
}
response <- do.call(rbind, responses)[md$observation_id, , drop = FALSE]
rm(responses, molecular)
invisible(gc())
stopifnot(identical(rownames(response), md$observation_id), all(is.finite(response)))
bank <- restore_functional_bank(bank, response)
programs <- colnames(bank$weights)
palette <- list(Program = setNames(hcl.colors(length(programs), 'Dark 3'), programs))
ro <- order(bank$partition, match(seq_along(bank$genes), bank$tree$order))
row_meta <- data.frame(Program = factor(programs[bank$partition[ro]], levels = programs),
  row.names = bank$genes[ro])
legend_args <- list(annotation_legend_param = list(Program = list(ncol = 3, labels_gp = grid::gpar(fontsize = 7))))
fwrite(data.table(order = seq_along(ro), gene = bank$genes[ro], program = as.character(row_meta$Program)),
  file.path(data_dir, 'global_gene_order_rank60.tsv'), sep = '\t')

# These are the exact two inputs and their mean, checked against saved fit hashes
# by restore_functional_bank. Identical gene order and color scale on all maps.
functional_correlation <- gene_correlation(t(bank$functional_support$loadings),
  rep(1, ncol(bank$functional_support$loadings)))
stopifnot(max(abs((bank$correlation + functional_correlation) / 2 - bank$clustering_correlation)) < 1e-12)
matrices <- list(gene_correlation = bank$correlation,
  gene_functional_loading_correlation = functional_correlation,
  gene_mean_clustering_correlation = bank$clustering_correlation)
labels <- c('Gene displacement correlation', 'Gene functional-loading correlation',
  'Mean of displacement and functional-loading correlations: actual Ward input')
qsave(list(genes = bank$genes, gene_order = ro, partition = bank$partition, matrices = matrices),
  file.path(data_dir, 'global_gene_relationships_rank60.qs'))
for (j in seq_along(matrices)) {
  mat <- matrices[[j]][ro, ro, drop = FALSE]
  dimnames(mat) <- list(bank$genes[ro], bank$genes[ro])
  hm <- draw_heatmap(mat, col_anno_df = row_meta, col_anno_cols = 'Program',
    row_anno_df = row_meta, row_anno_cols = 'Program', palette_overrides = palette,
    ann_row_args = legend_args, ann_top_args = list(show_legend = FALSE),
    hm_name = 'Pearson r', scale_rows = FALSE, value_range = c(-1, 0, 1),
    value_colors = c('#2166AC', 'white', '#B2182B'), raster_resize_mat = FALSE,
    title = paste0(labels[j], '\n', nrow(mat), ' genes; 60 fixed programs; identical row/column order'),
    column_title_gp = grid::gpar(fontsize = 10))
  save_heatmap(hm, names(matrices)[j], 10, 8, figure_dir)
  message('READY: ', file.path(figure_dir, paste0(names(matrices)[j], '.pdf')))
}

# Match the existing overview order: dominant program, PRISM sensitivity,
# ALMANAC-only sensitivity, then no outcome data. Original scores and scale.
annotation <- sensitivity[match(md$observation_id, observation_id)]
stopifnot(identical(annotation$observation_id, md$observation_id))
dominant <- max.col(abs(sweep(bank$scores, 2, bank$score_sd, '/')), ties.method = 'first')
source_rank <- fifelse(!is.na(annotation$prism_abundance_pct), 1L,
  fifelse(!is.na(annotation$almanac_abundance_pct), 2L, 3L))
outcome <- fifelse(source_rank == 1L, annotation$prism_abundance_pct, annotation$almanac_abundance_pct)
co <- order(dominant, source_rank, outcome, annotation$almanac_abundance_pct, na.last = TRUE)
md <- md[co]
annotation <- annotation[co]
dominant <- dominant[co]
# Contiguous bins include every signature exactly once, with no cross-block
# averaging and no missing outcome replaced with zero in the exported values.
block <- rleid(dominant, !is.na(annotation$prism_abundance_pct), !is.na(annotation$almanac_abundance_pct))
bin_width <- ceiling(nrow(md) / display_columns)
bins <- data.table(position = seq_len(nrow(md)), block = block)
bins[, local_bin := ceiling(seq_len(.N) / bin_width), by = block]
bins[, display_bin := rleid(block, local_bin)]
nb <- max(bins$display_bin)
counts <- tabulate(bins$display_bin, nb)
projection <- Matrix::sparseMatrix(i = seq_len(nrow(md)), j = bins$display_bin,
  x = 1 / counts[bins$display_bin], dims = c(nrow(md), nb))
first <- match(seq_len(nb), bins$display_bin)
cm <- data.frame(Program = factor(programs[dominant[first]], levels = programs),
  row.names = paste0('bin_', seq_len(nb)))
order_table <- cbind(data.table(order = seq_len(nrow(md)), display_bin = bins$display_bin,
  dominant_program = programs[dominant]), md)
order_table[, prism_abundance_pct := annotation$prism_abundance_pct]
order_table[, almanac_abundance_pct := annotation$almanac_abundance_pct]
fwrite(order_table, file.path(data_dir, 'global_condition_order_rank60.tsv'), sep = '\t')
bar_data <- order_table[, .(signatures = .N, PRISM = mean(prism_abundance_pct),
  ALMANAC = mean(almanac_abundance_pct)), by = display_bin]
fwrite(bar_data, file.path(data_dir, 'global_display_bins_rank60.tsv'), sep = '\t')

# Reused across displacement and assay loading views, preserving all ordering.
top_annotation <- function() {
  out <- list(show_legend = FALSE, annotation_name_gp = grid::gpar(fontsize = 8))
  for (assay in c('PRISM', 'ALMANAC')) {
    values <- bar_data[[assay]]
    cap <- max(125, ceiling(max(values, na.rm = TRUE) / 25) * 25)
    out[[paste(assay, '% vehicle')]] <- ComplexHeatmap::anno_barplot(replace(values, is.na(values), 0),
      gp = grid::gpar(fill = if (assay == 'PRISM') '#27818A' else '#B24C35', col = NA),
      ylim = c(0, cap), height = grid::unit(10, 'mm'), axis_param = list(at = c(0, 100, cap)))
    out[[paste(assay, 'value')]] <- sensitivity_track(values)
  }
  out
}
# Display reduction is explicit; the full exact source values remain in the
# frozen model and original input files. No reduction enters any analysis.
displacement <- t(as.matrix(Matrix::crossprod(projection, response[co, ro, drop = FALSE])))
dimnames(displacement) <- list(bank$genes[ro], rownames(cm))
loading_limit <- max(.1, as.numeric(quantile(abs(bank$functional_support$loadings), .98)))
plot_matrices <- list(displacement_overview = displacement)
for (assay in c('PRISM', 'ALMANAC')) {
  # data.table lookup uses explicit parent value to avoid ambiguous column names.
  selected <- which(bank$functional_support$observations$assay == assay)
  obs <- bank$functional_support$observations[selected]
  positions <- match(obs$observation_id, md$observation_id)
  stopifnot(!anyNA(positions))
  assay_counts <- tabulate(bins$display_bin[positions], nb)
  assay_projection <- Matrix::sparseMatrix(i = seq_along(positions), j = bins$display_bin[positions],
    x = 1 / assay_counts[bins$display_bin[positions]], dims = c(length(positions), nb))
  loading <- as.matrix(bank$functional_support$loadings[ro, selected, drop = FALSE] %*% assay_projection)
  loading[, assay_counts == 0] <- NA_real_
  dimnames(loading) <- dimnames(displacement)
  plot_matrices[[paste0('functional_loadings_', tolower(assay))]] <- loading
}
plot_matrices$program_scores <- t(as.matrix(Matrix::crossprod(projection, bank$scores[co, , drop = FALSE])))
colnames(plot_matrices$program_scores) <- rownames(cm)
qsave(list(matrices = plot_matrices, bars = bar_data, column_annotation = cm, row_annotation = row_meta,
  note = 'Display-only adjacent-column means; exact signature-to-bin map in global_condition_order_rank60.tsv'),
  file.path(data_dir, 'global_display_matrices_rank60.qs'))
for (name in names(plot_matrices)) {
  mat <- plot_matrices[[name]]
  is_score <- name == 'program_scores'
  is_loading <- grepl('functional_loadings', name)
  rm <- if (is_score) data.frame(Program = factor(programs, levels = programs), row.names = programs) else row_meta
  limit <- if (is_loading) loading_limit else max(.25, as.numeric(quantile(abs(mat), .98, na.rm = TRUE)))
  stopifnot(identical(rownames(mat), rownames(rm)), identical(colnames(mat), rownames(cm)))
  hm <- draw_heatmap(mat, col_anno_df = cm, col_anno_cols = 'Program', row_anno_df = rm,
    row_anno_cols = 'Program', palette_overrides = palette, ann_row_args = legend_args,
    ann_top_args = top_annotation(), scale_rows = FALSE, raster_resize_mat = FALSE,
    hm_name = if (is_loading) 'Signed gain' else 'log2 delta', na_col = '#D9D9D9',
    value_range = c(-limit, 0, limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    heatmap_legend_param = list(at = c(-limit, 0, limit), labels = sprintf('%.2f', c(-limit, 0, limit))),
    row_split = if (is_score) NULL else rm$Program, cluster_row_slices = FALSE,
    row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 5), row_gap = grid::unit(1.5, 'mm'),
    show_row_names = is_score, row_names_gp = grid::gpar(fontsize = 7),
    title = paste0(gsub('_', ' ', name), ' | 60 fixed global programs | 42,679 signatures\n',
      'Within dominant program: PRISM sensitivity, ALMANAC-only sensitivity, then no data\n',
      'Display only: adjacent-column means (maximum ', bin_width, ' signatures); no row scaling; grey = no outcome'),
    column_title_gp = grid::gpar(fontsize = 10))
  save_heatmap(hm, name, 16, if (is_score) 10 else 13.6, figure_dir, list(sensitivity_legend))
  message('READY: ', file.path(figure_dir, paste0(name, '.pdf')))
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'global_plot_sessionInfo.txt'))
fwrite(data.table(rank = 60L, genes = length(ro), signatures = nrow(md), display_columns = nb,
  max_signatures_per_bin = max(counts), mean_correlation_verified = TRUE,
  saved_bank_unchanged = identical(bank_hash, digest::digest(qread(bank_path), algo = 'sha256'))),
  file.path(data_dir, 'global_plot_checks_rank60.tsv'), sep = '\t')
message('COMPLETE: all global rank-60 maps')
writeLines(c('Script 42: fixed global bank with 60 programs, 2378 genes and 42679 signatures.',
  'Outcome-informed exploratory dictionary; displayed outcomes contributed to gene selection.',
  'Global figures: displacement_overview.pdf; functional_loadings_prism.pdf; functional_loadings_almanac.pdf;',
  'program_scores.pdf; gene_correlation.pdf; gene_functional_loading_correlation.pdf; gene_mean_clustering_correlation.pdf.',
  'All gene maps use the same frozen gene/program order; displacement and loading maps share column order.',
  'Global signature maps display adjacent-column means, at most 20 signatures per bin, without row scaling.',
  'Bins never cross dominant-program or assay-availability blocks. All signatures contribute.',
  'Gene-by-gene matrices retain every exact pair value; their verified 50/50 mean is the actual Ward input.',
  'Source values and orders: data/global_gene_relationships_rank60.qs; data/global_display_matrices_rank60.qs;',
  'data/global_gene_order_rank60.tsv; data/global_condition_order_rank60.tsv; data/global_display_bins_rank60.tsv.',
  'Per-drug maps in drug_*/figures retain individual signatures and use the same global bank.',
  'Plotting entry point: /work/src/program_pilot/42_plot_global_program_matrices.R'),
  file.path(output_dir, 'README.txt'))
