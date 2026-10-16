#!/usr/bin/env Rscript
# Exact script-15 gene-pair means, original-program comparison, fixed global Ward programs.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(cluster); library(patchwork)})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
arguments <- commandArgs(trailingOnly = TRUE)
check(length(arguments) == 0L || identical(arguments, '--aggregate-only') || identical(arguments, '--matrices-only'),
  'Usage: Rscript 16_global_mean_relationships_fixed_ward_programs.R [--aggregate-only|--matrices-only]')
aggregate_only <- identical(arguments, '--aggregate-only')
matrices_only <- identical(arguments, '--matrices-only')
local_rank_grid <- 2:80
global_rank_override <- 40L
local_minimum_rank <- 2L
local_resamples <- 0L
outcome_p_threshold <- 0.05
outcome_minimum_abs_correlation <- 0
functional_input_transform <- 'Exact script-15 pair means across cells where used; fixed global Ward.D2 partition'
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
data.table::setDTthreads(2L)
options(warn = 1)
theme_set(theme_classic(base_size = 10))
parent_dir <- file.path(work_dir, 'results/program_pilot/15_mean_functional_displacement_correlation_programs_ward')
prep_path <- file.path(work_dir, 'results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
output_dir <- file.path(work_dir, 'results/program_pilot/16_global_mean_relationships_fixed_ward_programs')
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
model_dir <- file.path(output_dir, 'models')
for (path in c(output_dir, data_dir, figure_dir, model_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
panel <- qread(prep_path)
cell_map <- panel$cell_map
sensitivity <- copy(panel$sensitivity)
completed <- fread(file.path(parent_dir, 'data/panel_completed.tsv'))
check(setequal(completed$depmap_id, cell_map$depmap_id) && all(completed$status == 'complete'), 'Script 15 must be complete')
scopes <- paste0('cell_', gsub('[^A-Za-z0-9]+', '_', cell_map$depmap_name))
bank_paths <- file.path(parent_dir, scopes, 'models/program_bank.qs')
banks <- lapply(bank_paths, qread)
genes <- sort(unique(unlist(lapply(banks, `[[`, 'genes'))))
source_paths <- c(prep_path, bank_paths, panel$input_paths)
hashes <- data.table(path = source_paths, sha256 = vapply(source_paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(hashes, file.path(data_dir, 'input_hashes.tsv'), sep = '\t')
signature <- digest::digest(list(hashes, aggregation = 'mean each pair only across script-15 cells where used'), algo = 'sha256')
aggregate_path <- file.path(model_dir, 'averaged_script15_relationships.qs')
if (file.exists(aggregate_path)) {
  aggregate <- qread(aggregate_path)
  check(identical(aggregate$input_signature, signature), 'Source script-15 matrices changed')
} else {
  totals <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  counts <- matrix(0L, length(genes), length(genes), dimnames = list(genes, genes))
  audit <- list()
  for (i in seq_along(banks)) {
    molecular <- qread(panel$input_paths[i])
    bank <- restore_functional_bank(banks[[i]], molecular$response)
    indices <- match(bank$genes, genes)
    totals[indices, indices] <- totals[indices, indices] + bank$clustering_correlation
    counts[indices, indices] <- counts[indices, indices] + 1L
    audit[[i]] <- data.table(recipient = cell_map$depmap_name[i], genes = length(bank$genes),
      input_pairs = choose(length(bank$genes), 2), relationship_sha256 = digest::digest(bank$clustering_correlation, algo = 'sha256'))
    fwrite(rbindlist(audit), file.path(data_dir, 'cell_relationship_contributions.tsv'), sep = '\t')
    fwrite(data.table(status = 'running', phase = 'exact pair aggregation', cells_ready = i, cells = length(banks)), file.path(output_dir, 'run_status.tsv'), sep = '\t')
    message('Exact script-15 matrix aggregated: ', cell_map$depmap_name[i], ' (', i, '/', length(banks), ')')
    rm(molecular, bank)
    invisible(gc())
  }
  average <- totals / pmax(counts, 1L)
  average[counts == 0L] <- NA_real_
  aggregate <- list(genes = genes, clustering_correlation = average, pair_cell_counts = counts,
    input_signature = signature, contribution_audit = rbindlist(audit))
  temporary <- tempfile(pattern = 'pair_means.', tmpdir = model_dir)
  qsave(aggregate, temporary, preset = 'high')
  check(file.rename(temporary, aggregate_path), 'Could not commit exact pair means')
  indices <- which(upper.tri(counts) & counts > 0L, arr.ind = TRUE)
  fwrite(data.table(gene_a = genes[indices[, 1L]], gene_b = genes[indices[, 2L]], mean_relationship = average[indices], contributing_cells = counts[indices]),
    file.path(data_dir, 'averaged_gene_pairs.tsv.gz'), sep = '\t')
  rm(totals, counts, average, indices)
  invisible(gc())
}
# Original script-15 seed membership. Zero is non-membership, not an imputed relationship.
comparison_path <- file.path(model_dir, 'original_program_comparison.qs')
comparison_figures <- file.path(figure_dir, paste0('original_programs_gene_membership.', c('pdf', 'png')))
reuse_comparison <- file.exists(comparison_path) && all(file.exists(comparison_figures)) &&
  identical(qread(comparison_path)$input_signature, signature)
if (!reuse_comparison) {
membership_parts <- column_parts <- list()
for (i in seq_along(banks)) {
  bank <- banks[[i]]
  if (is.null(bank$weights)) next
  programs <- colnames(bank$weights)
  part <- matrix(0L, length(genes), length(programs), dimnames = list(genes, paste(scopes[i], programs, sep = ':')))
  part[cbind(match(bank$genes, genes), bank$partition)] <- 1L
  membership_parts[[i]] <- part
  column_parts[[i]] <- data.table(column = colnames(part), cell_line = cell_map$depmap_name[i], program = programs)
}
membership <- do.call(cbind, membership_parts)
columns <- rbindlist(column_parts)
check(identical(colnames(membership), columns$column), 'Program comparison alignment failed')
fwrite(columns, file.path(data_dir, 'original_program_comparison_columns.tsv'), sep = '\t')
fwrite(cbind(data.table(gene = genes), as.data.table(membership)), file.path(data_dir, 'original_gene_program_membership.tsv.gz'), sep = '\t')
row_tree <- hclust(dist(membership), method = 'ward.D2')
column_tree <- hclust(dist(t(membership)), method = 'ward.D2')
column_annotation <- data.frame(Cell = factor(columns$cell_line, levels = cell_map$depmap_name), row.names = columns$column)
colors <- setNames(hcl.colors(nrow(cell_map), 'Dark 3'), cell_map$depmap_name)
map <- draw_heatmap(membership, scale_rows = FALSE, cluster_rows = row_tree, cluster_columns = column_tree,
  row_dend_reorder = FALSE, column_dend_reorder = FALSE, show_row_names = FALSE, show_column_names = TRUE,
  column_labels = paste(columns$cell_line, columns$program), column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 3.5),
  col_anno_df = column_annotation, col_anno_cols = 'Cell', palette_overrides = list(Cell = colors),
  ann_top_args = list(annotation_legend_param = list(Cell = list(ncol = 2))), value_range = c(0, 1), value_colors = c('white', '#2166AC'),
  hm_name = 'Seed member', title = 'Original script-15 programs across cell lines | Ward.D2 on genes and programs\n1: seed member; 0: non-member; unresolved cells have no program columns')
save_heatmap(map, 'original_programs_gene_membership', 17, 10)
qsave(list(row_tree = row_tree, column_tree = column_tree, genes = genes, columns = columns, input_signature = signature), file.path(model_dir, 'original_program_comparison.qs'))
message('Original gene-by-program comparison ready: ', nrow(membership), ' genes x ', ncol(membership), ' programs')
rm(membership_parts, membership, map)
} else message('Reusing saved original gene-by-program comparison')
if (!matrices_only) fwrite(data.table(status = 'aggregated', phase = 'exact pair union and original-program comparison ready', cells = length(banks)), file.path(output_dir, 'run_status.tsv'), sep = '\t')
writeLines(c('Exact script-15 relationships only. No expanded loadings or newly fitted regressions.',
  'Each pair is averaged only across cells where it was used; counts and observed-pair union are exported.',
  'Unused pairs stay NA in raw averages; the finite Ward input encodes their correlation as zero.',
  'Per-cell plots retain each cell\'s original script-15 gene set; global rank is explicitly fixed at 40.',
  'original_programs_gene_membership.pdf/png compares original seed programs, clustered on both axes.'), file.path(output_dir, 'README.txt'))
if (aggregate_only) quit(save = 'no', status = 0L)
# The union is a sparse relationship graph. Only observed pairs carry a mean.
# No-edge entries are zero in the finite Ward input; raw averages stay NA there,
# and these no-edge zeros never enter the observed-pair means or their counts.
clustering_input <- aggregate$clustering_correlation
clustering_input[aggregate$pair_cell_counts == 0L] <- 0
distance <- as.dist(sqrt(2 * (1 - clustering_input)))
check(all(is.finite(distance)), 'Invalid aggregated relationship distances')
global <- list(scope = '.', orientation = 'pooled', genes = genes, metadata = panel$metadata,
  display_name = 'Global observed-pair mean of script-15 relationships', clustering_correlation = clustering_input,
  distance = distance, tree = hclust(distance, method = 'ward.D2'), path = file.path(model_dir, 'global_gene_relationships.qs'),
  input_signature = signature, nmf_input_signature = digest::digest(clustering_input, algo = 'sha256'),
  functional_support = list(settings = list(clustering_algorithm = 'ward.D2', rank_selection_rule = 'silhouette_elbow')))
qsave(global, global$path, preset = 'high')
fits <- fit_ward_cuts(distance, global$tree, local_rank_grid, '.', 'pooled', global$nmf_input_signature, model_dir)
global$rank_override <- global_rank_override
global$ranks <- select_silhouette(rbindlist(lapply(fits, `[[`, 'diagnostics')), rank_override = global_rank_override, selection_rule = 'silhouette_elbow')
global$all_partitions <- lapply(fits, `[[`, 'partition')
global$resampling <- banks[[1L]]$resampling
selected <- global$ranks[selected == TRUE, rank]
check(length(selected) == 1L && selected == global_rank_override, 'Requested global rank unavailable')
global$partition <- global$all_partitions[[as.character(selected)]]
global$diagnostics <- data.table(scope = '.', orientation = 'pooled', conditions = nrow(panel$metadata), discovery_units = nrow(cell_map),
  genes = length(genes), selected_rank = selected, selection_status = 'selected', raw_tree_cophenetic = cor(as.numeric(distance), as.numeric(cophenetic(global$tree))))
qsave(global, global$path, preset = 'high')
render_diagnostics(global, rank_plots_only = TRUE)
if (!matrices_only) fwrite(data.table(status = 'running', phase = 'fixed-global per-cell plotting', selected_rank = selected, cells_ready = 0L, cells = length(banks)), file.path(output_dir, 'run_status.tsv'), sep = '\t')
fwrite(data.table(gene = genes, global_program = sprintf('P%02d', global$partition)), file.path(data_dir, 'global_gene_assignments.tsv'), sep = '\t')
# Separate the exact pairwise averages into their two original script-15 views.
component_path <- file.path(model_dir, 'averaged_relationship_components.qs')
if (file.exists(component_path)) {
  components <- qread(component_path)
  check(identical(components$input_signature, signature), 'Component inputs changed')
} else {
  functional_totals <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  for (i in seq_along(banks)) {
    loading <- banks[[i]]$functional_support$loadings
    component <- gene_correlation(t(loading), rep(1, ncol(loading)))
    indices <- match(banks[[i]]$genes, genes)
    functional_totals[indices, indices] <- functional_totals[indices, indices] + component
  }
  functional_average <- functional_totals / pmax(aggregate$pair_cell_counts, 1L)
  functional_average[aggregate$pair_cell_counts == 0L] <- NA_real_
  # For each observed pair, mean(displacement) = 2*mean(mixed) - mean(functional).
  components <- list(functional_correlation = functional_average,
    displacement_correlation = 2 * aggregate$clustering_correlation - functional_average,
    input_signature = signature)
  qsave(components, component_path, preset = 'high')
  rm(functional_totals, functional_average, component)
}
gene_order <- genes[global$tree$order]
programs <- sprintf('P%02d', seq_len(selected))
annotation <- data.frame(Program = factor(sprintf('P%02d', global$partition[gene_order]), levels = programs), row.names = gene_order)
# Precompute ordinary raster averages in R instead of repeating millions of
# callbacks inside the graphics device. Bins never cross a program boundary.
bin_size <- ceiling(length(genes) / 1200L)
within_program <- ave(seq_along(gene_order), annotation$Program, FUN = seq_along)
gene_bins <- paste0(annotation$Program, ':', ceiling(within_program / bin_size))
bin_ids <- unique(gene_bins)
bin_sizes <- as.numeric(table(factor(gene_bins, levels = bin_ids)))
bin_annotation <- annotation[match(bin_ids, gene_bins), , drop = FALSE]
rownames(bin_annotation) <- bin_ids
observed <- !is.na(aggregate$clustering_correlation[gene_order, gene_order, drop = FALSE])
observed_bin_counts <- t(rowsum(t(rowsum(observed * 1L, gene_bins, reorder = FALSE)), gene_bins, reorder = FALSE))
rm(observed)
palette <- list(Program = setNames(hcl.colors(selected, 'Dark 3'), programs))
matrix_views <- list(global_ward_input = global$clustering_correlation,
  averaged_gene_relationships = aggregate$clustering_correlation,
  averaged_functional_loading_correlation = components$functional_correlation,
  averaged_displacement_correlation = components$displacement_correlation)
view_titles <- c('Actual Ward input | unobserved edges = 0',
  'Observed gene-pair means | gray = unobserved',
  'Average functional-loading correlation | gray = unobserved',
  'Average displacement correlation | gray = unobserved')
fwrite(data.table(position = seq_along(gene_order), gene = gene_order, raster_bin = gene_bins,
  global_program = as.character(annotation$Program)), file.path(data_dir, 'global_heatmap_gene_order.tsv'), sep = '\t')
for (i in seq_along(matrix_views)) {
  values <- matrix_views[[i]][gene_order, gene_order, drop = FALSE]
  bin_sum <- t(rowsum(t(rowsum(values, gene_bins, reorder = FALSE, na.rm = TRUE)), gene_bins, reorder = FALSE, na.rm = TRUE))
  bin_count <- if (i == 1L) outer(bin_sizes, bin_sizes) else observed_bin_counts
  raster_values <- bin_sum / bin_count
  raster_values[bin_count == 0L] <- NA_real_
  map <- draw_heatmap(raster_values,
    col_anno_df = bin_annotation, col_anno_cols = 'Program', row_anno_df = bin_annotation, row_anno_cols = 'Program',
    palette_overrides = palette, ann_top_args = list(show_legend = FALSE),
    ann_row_args = list(show_annotation_name = FALSE, annotation_legend_param = list(Program = list(ncol = 2))),
    scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'), na_col = '#D9D9D9',
    raster_resize_mat = FALSE,
    hm_name = 'Pearson r', title = paste0(view_titles[i], '\nRank 40 | ', length(genes), ' genes | same order | raster means'))
  save_heatmap(map, names(matrix_views)[i], 8.8, 7.6)
  message('Rank-40 matrix heatmap ready: ', names(matrix_views)[i])
  rm(values, bin_sum, raster_values)
}
rm(components, matrix_views, map)
invisible(gc())
if (matrices_only) quit(save = 'no', status = 0L)
selection <- fread(file.path(parent_dir, 'data/gene_selection.tsv.gz'))
gene_selection <- outcome_test_tables <- outcome_test_summary <- list()
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm = TRUE) / 25) * 25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm = TRUE) / 25) * 25)
artifacts <- list()
for (i in seq_along(banks)) {
  id <- cell_map$depmap_id[i]
  scope <- scopes[i]
  parent <- banks[[i]]
  molecular <- qread(panel$input_paths[i])
  prepared <- list(response = molecular$response)
  gene_selection[[id]] <- selection[depmap_id == id]
  outcome_test_tables[[id]] <- fread(file.path(parent_dir, scope, 'data/gene_outcome_correlations.tsv.gz'))
  outcome_test_summary[[id]] <- fread(file.path(parent_dir, scope, 'data/gene_outcome_test_summary.tsv'))
  for (folder in c('models', 'data', 'figures')) dir.create(file.path(output_dir, scope, folder), recursive = TRUE, showWarnings = FALSE)
  bank <- discover_local(molecular$response, molecular$null[[id]], parent$metadata, parent$genes, 'drug', scope, 'cell',
    support = parent$functional_support, sensitivity_sign = TRUE, rank_override = global_rank_override, shared_structure = global)
  save_functional_bank(bank, file.path(output_dir, scope, 'models/program_bank.qs'))
  render_diagnostics(bank)
  render_program_score_clustered_maps(bank, molecular$response)
  render_local(bank)
  render_functional_loadings(bank)
  names <- c('silhouette.pdf', 'diagnostics.pdf', 'gene_correlation.pdf', 'displacement_overview.pdf', 'displacement_pages.pdf',
    'program_scores.pdf', 'focal_panel.pdf', 'functional_loadings.pdf', 'functional_loadings_focal.pdf',
    'displacement_program_score_clustered.pdf', 'functional_loadings_program_score_clustered.pdf')
  artifacts[[i]] <- data.table(scope = scope, artifact = names, path = file.path(output_dir, scope, 'figures', names))
  fwrite(rbindlist(artifacts), file.path(output_dir, 'artifact_index.tsv'), sep = '\t')
  fwrite(data.table(status = 'running', phase = 'fixed-global per-cell plotting', selected_rank = selected, cells_ready = i, cells = length(banks)), file.path(output_dir, 'run_status.tsv'), sep = '\t')
  message('Fixed-global figures ready: ', scope)
  rm(molecular, prepared, bank)
  invisible(gc())
}
write_manifest('16_global_mean_relationships_fixed_ward_programs', source_paths, list(pair_aggregation = 'only cells using each pair',
  fixed_global_rank = selected, original_cell_genes_preserved = TRUE))
fwrite(data.table(status = 'complete', cells = length(banks), selected_rank = selected), file.path(output_dir, 'run_status.tsv'), sep = '\t')
