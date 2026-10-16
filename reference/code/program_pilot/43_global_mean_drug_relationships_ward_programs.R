#!/usr/bin/env Rscript
# Average original script-41 relationships; project every drug through fixed global programs.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(Matrix)})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(4L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn = 1)
parent_dir <- '/work/results/program_pilot/41_drug_first_functional_displacement_programs_ward'
output_dir <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs'
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
model_dir <- file.path(output_dir, 'models')
rank <- 34L # User requested the low-30s region, before the sharp silhouette drop at35.
plot_workers <- 4L
for (path in c(output_dir, data_dir, figure_dir, model_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
cell_map <- panel$cell_map
sensitivity <- copy(panel$sensitivity)
coverage <- fread(file.path(parent_dir, 'data/drug_assay_coverage.tsv'))
plot_drugs <- coverage[eligible == TRUE, drug]
plot_drugs <- c(intersect(local_drugs, plot_drugs), sort(setdiff(plot_drugs, local_drugs)))
source_paths <- sort(list.files(parent_dir, pattern = 'program_bank.qs$', recursive = TRUE, full.names = TRUE))
check(length(source_paths) == 102L, 'Expected all 102 completed script-41 matrices')
banks <- lapply(source_paths, qread)
genes <- sort(unique(unlist(lapply(banks, `[[`, 'genes'))))
check(!anyDuplicated(genes) && all(genes %chin% panel$genes), 'Invalid global gene union')
source_drugs <- vapply(banks, function(b) unique(b$metadata$drug), character(1))
check(!anyDuplicated(source_drugs) && all(source_drugs %chin% plot_drugs), 'Duplicate/unknown source drugs')
metadata <- copy(panel$metadata[drug %chin% plot_drugs])
# All original columns and doses remain separate; load only the union genes and eligible drugs.
responses <- list()
for (id in cell_map$depmap_id) {
  molecular <- qread(panel$input_paths[match(id, cell_map$depmap_id)])
  ids <- metadata[depmap_id == id, observation_id]
  responses[[id]] <- molecular$response[ids, genes, drop = FALSE]
}
response <- do.call(rbind, responses)[metadata$observation_id, , drop = FALSE]
rm(responses, molecular)
invisible(gc())
check(identical(rownames(response), metadata$observation_id) && all(is.finite(response)), 'Invalid projection inputs')
input_hashes <- data.table(path = source_paths,
  sha256 = vapply(source_paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(input_hashes, file.path(data_dir, 'source_bank_hashes.tsv'), sep = '\t')
input_signature <- digest::digest(list(input_hashes, genes, rule = 'pair means only over drugs using each pair'), algo = 'sha256')
aggregate_path <- file.path(model_dir, 'averaged_script41_relationships.qs')
if (file.exists(aggregate_path)) {
  aggregate <- qread(aggregate_path)
  check(identical(aggregate$input_signature, input_signature), 'Script-41 source banks changed')
} else {
  mixed_sum <- functional_sum <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  counts <- matrix(0L, length(genes), length(genes), dimnames = list(genes, genes))
  audit <- list()
  for (i in seq_along(banks)) {
    b <- restore_functional_bank(banks[[i]], response)
    ii <- match(b$genes, genes)
    functional <- gene_correlation(t(b$functional_support$loadings), rep(1, ncol(b$functional_support$loadings)))
    check(max(abs((b$correlation + functional) / 2 - b$clustering_correlation)) < 1e-12,
      'Original script-41 matrix components differ')
    mixed_sum[ii, ii] <- mixed_sum[ii, ii] + b$clustering_correlation
    functional_sum[ii, ii] <- functional_sum[ii, ii] + functional
    counts[ii, ii] <- counts[ii, ii] + 1L
    audit[[i]] <- data.table(drug = source_drugs[i], genes = length(ii), pairs = choose(length(ii), 2),
      relationship_sha256 = digest::digest(b$clustering_correlation, algo = 'sha256'))
    fwrite(rbindlist(audit), file.path(data_dir, 'drug_relationship_contributions.tsv'), sep = '\t')
    message('AGGREGATED ', i, '/', length(banks), ': ', source_drugs[i])
    rm(b, functional)
  }
  mixed_mean <- mixed_sum / pmax(counts, 1L)
  functional_mean <- functional_sum / pmax(counts, 1L)
  mixed_mean[counts == 0L] <- NA_real_
  functional_mean[counts == 0L] <- NA_real_
  aggregate <- list(genes = genes, clustering_correlation = mixed_mean,
    functional_correlation = functional_mean,
    displacement_correlation = 2 * mixed_mean - functional_mean,
    pair_drug_counts = counts, contribution_audit = rbindlist(audit), input_signature = input_signature)
  qsave(aggregate, aggregate_path, preset = 'high')
  rm(mixed_sum, functional_sum, mixed_mean, functional_mean, counts)
  invisible(gc())
}
message('SAVED exact drug-pair means: ', length(genes), ' genes')
# Raw unobserved pairs remain NA. Neutral correlation0 supplies finite distances,
# exactly as in script16; these zeros never enter the observed-pair averages.
clustering_input <- aggregate$clustering_correlation
clustering_input[aggregate$pair_drug_counts == 0L] <- 0
clustering_input[clustering_input > 1] <- 1
clustering_input[clustering_input < -1] <- -1
distance <- as.dist(sqrt(2 * (1 - clustering_input)))
check(all(is.finite(distance)), 'Nonfinite global Ward distance')
bank_path <- file.path(model_dir, 'program_bank.qs')
previous_bank <- if (file.exists(bank_path)) qread(bank_path) else NULL
if (!is.null(previous_bank)) {
  check(identical(previous_bank$input_signature, input_signature) && identical(previous_bank$genes, genes),
    'Saved hierarchy inputs changed')
  tree <- previous_bank$tree
} else tree <- hclust(distance, method = 'ward.D2')
partition <- cutree(tree, k = rank)
programs <- sprintf('P%02d', seq_len(rank))
# Hard means use every assigned gene and no gene from another program.
weights <- sparseMatrix(i = seq_along(genes), j = partition,
  x = 1 / tabulate(partition, rank)[partition], dims = c(length(genes), rank),
  dimnames = list(genes, programs))
scores <- as.matrix(response %*% weights)
score_sd <- matrixStats::colSds(scores)
check(all(score_sd > 0), 'Constant global program score')
global <- list(genes = genes, tree = tree, partition = partition, weights = weights,
  metadata = metadata, scores = scores, score_sd = score_sd,
  diagnostics = data.table(selected_rank = rank), input_signature = input_signature,
  aggregate_path = aggregate_path, scoring = 'mean of all and only assigned genes')
if (!is.null(previous_bank$ranks)) {
  global$ranks <- select_silhouette(previous_bank$ranks, rank_override = rank)
  global$all_partitions <- previous_bank$all_partitions
  fwrite(global$ranks, file.path(data_dir, 'rank_diagnostics.tsv'), sep = '\t')
}
qsave(global, bank_path, preset = 'high')
missing_pairs <- sum(aggregate$pair_drug_counts == 0L) / (length(genes)^2 - length(genes))
fwrite(data.table(rank = rank, genes = length(genes), contributing_drugs = length(banks),
  projected_drugs = length(plot_drugs), signatures = nrow(metadata), unobserved_pair_fraction = missing_pairs),
  file.path(data_dir, 'aggregate_summary.tsv'), sep = '\t')
fwrite(data.table(gene = genes, program = programs[partition], contributing_drugs = diag(aggregate$pair_drug_counts)),
  file.path(data_dir, 'global_gene_assignments.tsv'), sep = '\t')
message('GLOBAL BANK SAVED: ', rank, ' programs; unobserved pairs ', round(100 * missing_pairs, 1), '%')

# Global matrix overviews reuse script16's display-only adjacent-gene reduction.
ro <- order(partition, match(seq_along(genes), tree$order))
row_meta <- data.frame(Program = factor(programs[partition[ro]], levels = programs), row.names = genes[ro])
palette <- list(Program = setNames(hcl.colors(rank, 'Dark 3'), programs))
within_program <- ave(seq_along(ro), row_meta$Program, FUN = seq_along)
gene_bins <- paste0(row_meta$Program, ':', ceiling(within_program / ceiling(length(genes) / 1200L)))
bin_ids <- unique(gene_bins)
bin_meta <- row_meta[match(bin_ids, gene_bins), , drop = FALSE]
rownames(bin_meta) <- bin_ids
bin_sizes <- as.numeric(table(factor(gene_bins, levels = bin_ids)))
observed <- aggregate$pair_drug_counts[ro, ro, drop = FALSE] > 0L
observed_counts <- t(rowsum(t(rowsum(observed * 1L, gene_bins, reorder = FALSE)), gene_bins, reorder = FALSE))
rm(observed)
fwrite(data.table(order = seq_along(ro), gene = genes[ro], program = as.character(row_meta$Program), display_bin = gene_bins),
  file.path(data_dir, 'global_heatmap_gene_order.tsv'), sep = '\t')
views <- list(averaged_gene_relationships = aggregate$clustering_correlation,
  averaged_displacement_correlation = aggregate$displacement_correlation,
  averaged_functional_loading_correlation = aggregate$functional_correlation,
  global_ward_input = clustering_input)
for (name in names(views)) {
  value <- views[[name]][ro, ro, drop = FALSE]
  sums <- t(rowsum(t(rowsum(value, gene_bins, reorder = FALSE, na.rm = TRUE)), gene_bins, reorder = FALSE, na.rm = TRUE))
  denominator <- if (name == 'global_ward_input') outer(bin_sizes, bin_sizes) else observed_counts
  plotted <- sums / denominator
  plotted[denominator == 0] <- NA_real_
  hm <- draw_heatmap(plotted, scale_rows = FALSE, col_anno_df = bin_meta, col_anno_cols = 'Program',
    row_anno_df = bin_meta, row_anno_cols = 'Program', palette_overrides = palette,
    ann_top_args = list(show_legend = FALSE), ann_row_args = list(show_annotation_name = FALSE,
      annotation_legend_param = list(Program = list(ncol = 3, labels_gp = grid::gpar(fontsize = 7)))),
    value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'), na_col = '#D9D9D9',
    raster_resize_mat = FALSE, hm_name = 'Pearson r', column_title_gp = grid::gpar(fontsize = 10),
    title = paste0(gsub('_', ' ', name), ' | mean across contributing script-41 drugs\n',
      length(genes), ' genes; ', rank, ' programs; same order; adjacent-gene means for display only\n',
      if (name == 'global_ward_input') 'Unobserved pairs = neutral0 in Ward input' else 'Grey = unobserved pairs'))
  save_heatmap(hm, name, 10, 8, figure_dir)
  message('GLOBAL MATRIX READY: ', name)
  rm(value, sums, plotted, hm)
}
rm(views, aggregate, clustering_input, distance, observed_counts)
invisible(gc())
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm = TRUE) / 25) * 25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm = TRUE) / 25) * 25)

# Each independent rendering job shares read-only molecular inputs. No clustering
# of genes, gene selection or regression is repeated per drug.
render_drug <- function(drug_name) {
  RhpcBLASctl::blas_set_num_threads(1L)
  scope <- paste0('drug_', gsub('[^A-Za-z0-9]+', '_', drug_name))
  dirs <- file.path(output_dir, scope, c('data', 'models', 'figures'))
  for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  idx <- which(metadata$drug == drug_name)
  md <- metadata[idx]
  regression_path <- file.path(parent_dir, scope, 'models/functional_regression.qs')
  check(file.exists(regression_path), paste('Missing original fitted regression:', drug_name))
  support <- qread(regression_path)$support
  # Expand the existing fitted loading rows to the full gene union for display.
  # Genes never fitted for this drug are explicit NA, never zero/imputed.
  original_loading <- support$loadings
  loading <- matrix(NA_real_, length(genes), ncol(original_loading), dimnames = list(genes, colnames(original_loading)))
  fitted_genes <- intersect(genes, rownames(original_loading))
  loading[fitted_genes, ] <- original_loading[fitted_genes, , drop = FALSE]
  support$loadings <- loading
  support$raw_prediction_gain <- NULL
  view <- list(scope = scope, orientation = 'drug', genes = genes, tree = tree, partition = partition,
    weights = weights, metadata = md, scores = scores[idx, , drop = FALSE], score_sd = score_sd,
    functional_support = support, global_structure = list(path = bank_path, program_colors = palette$Program))
  render_program_score_clustered_maps(view, response)
  # Sensitivity-ordered full-gene map: all genes, without top-gene sampling.
  ann <- sensitivity[match(md$observation_id, observation_id)]
  source <- fifelse(!is.na(ann$prism_abundance_pct), 1L, fifelse(!is.na(ann$almanac_abundance_pct), 2L, 3L))
  abundance <- fifelse(source == 1L, ann$prism_abundance_pct, ann$almanac_abundance_pct)
  co <- order(source, abundance, ann$almanac_abundance_pct, md$recipient, md$dose, na.last = TRUE)
  cm <- data.frame(Dose = factor(sprintf('%g', md$dose[co])), Cell = factor(md$recipient[co]), row.names = md$observation_id[co])
  av <- ann$almanac_abundance_pct[co]
  pv <- ann$prism_abundance_pct[co]
  top <- list(show_legend = c(Dose = TRUE, Cell = FALSE),
    `PRISM % vehicle` = ComplexHeatmap::anno_barplot(replace(pv, is.na(pv), 0),
      ylim = c(0, prism_bar_max), height = grid::unit(11, 'mm'), gp = grid::gpar(fill = '#27818A', col = NA)),
    `ALMANAC % vehicle` = ComplexHeatmap::anno_barplot(replace(av, is.na(av), 0),
      ylim = c(0, almanac_bar_max), height = grid::unit(11, 'mm'), gp = grid::gpar(fill = '#B24C35', col = NA)),
    `PRISM value` = sensitivity_track(pv), `ALMANAC value` = sensitivity_track(av))
  displacement <- t(response[md$observation_id[co], genes[ro], drop = FALSE])
  limit <- max(.25, as.numeric(quantile(abs(displacement), .98)))
  hm <- draw_heatmap(displacement, scale_rows = FALSE, row_anno_df = row_meta, row_anno_cols = 'Program',
    col_anno_df = cm, col_anno_cols = c('Dose', 'Cell'), ann_top_args = top, palette_overrides = palette,
    ann_row_args = list(show_annotation_name = FALSE, annotation_legend_param = list(Program = list(ncol = 3))),
    row_split = row_meta$Program, cluster_row_slices = FALSE, row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 5), row_gap = grid::unit(2, 'mm'),
    raster_resize_mat = FALSE, value_range = c(-limit, 0, limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'log2 delta', column_title_gp = grid::gpar(fontsize = 10),
    title = paste0(drug_name, ' | all ', length(genes), ' genes in ', rank, ' fixed global programs\n',
      'PRISM sensitivity, ALMANAC-only sensitivity, then no data; original displacements; no row scaling'))
  save_heatmap(hm, 'displacement_sensitivity_ordered', 16, 13.6, dirs[3], list(sensitivity_legend))
  score_matrix <- t(scores[idx[co], , drop = FALSE])
  sm <- data.frame(Program = factor(programs, levels = programs), row.names = programs)
  score_limit <- max(.25, as.numeric(quantile(abs(score_matrix), .98)))
  hm <- draw_heatmap(score_matrix, scale_rows = FALSE, row_anno_df = sm, row_anno_cols = 'Program',
    col_anno_df = cm, col_anno_cols = c('Dose', 'Cell'), ann_top_args = top, palette_overrides = palette,
    ann_row_args = list(show_annotation_name = FALSE, annotation_legend_param = list(Program = list(ncol = 3))),
    show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7), raster_resize_mat = FALSE,
    value_range = c(-score_limit, 0, score_limit), value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'Mean delta',
    title = paste0(drug_name, ' | means of all assigned genes | same sensitivity order'))
  save_heatmap(hm, 'program_scores_sensitivity_ordered', 16, 10, dirs[3], list(sensitivity_legend))
  fwrite(cbind(md[co, .(observation_id, recipient, drug, dose)], ann[co, .(prism_abundance_pct, almanac_abundance_pct)]),
    file.path(dirs[1], 'condition_order_sensitivity.tsv'), sep = '\t')
  fwrite(data.table(gene = genes[ro], program = as.character(row_meta$Program), order = seq_along(ro)),
    file.path(dirs[1], 'gene_order.tsv'), sep = '\t')
  fwrite(cbind(md[, .(observation_id, recipient, drug, dose)], as.data.table(scores[idx, , drop = FALSE])),
    file.path(dirs[1], 'program_scores.tsv'), sep = '\t')
  qsave(list(global_bank = bank_path, genes = genes, partition = partition, metadata = md,
    scores = scores[idx, , drop = FALSE], functional_support = support), file.path(dirs[2], 'projection.qs'))
  data.table(drug = drug_name, status = 'plots complete', genes = length(genes), programs = rank,
    signatures = length(idx), fitted_loading_genes = length(fitted_genes), figure_dir = dirs[3])
}
status <- list()
message('PROJECTING ', length(plot_drugs), ' drugs through the same fixed global programs')
for (batch in split(plot_drugs, ceiling(seq_along(plot_drugs) / plot_workers))) {
  results <- parallel::mclapply(batch, render_drug, mc.cores = plot_workers, mc.preschedule = FALSE)
  check(all(vapply(results, is.data.table, logical(1))), 'A drug plotting worker failed; inspect the log')
  status <- c(status, results)
  fwrite(rbindlist(status), file.path(output_dir, 'drug_status.tsv'), sep = '\t')
  message('DRUG PLOTS READY: ', length(status), '/', length(plot_drugs))
}
writeLines(c('Global Ward.D2 programs from averaged original script41 drug-specific gene relationships.',
  'Mean each observed pair only over drugs using that pair; no new global regression or gene screen.',
  'All8429 global genes appear in every drug displacement heatmap; no top-gene subsampling.',
  'Program scores are hard means of all and only assigned genes. Doses remain separate.',
  'Functional maps reuse original script41 fitted regressions; never-fitted genes are grey NA.',
  'Raw unobserved gene pairs remain NA; actual finite Ward input uses neutral correlation0.',
  'Outcome-informed retrospective discovery; no independent validation claimed.'), file.path(output_dir, 'README.txt'))
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
fwrite(data.table(status = 'complete', contributing_drugs = length(banks), projected_drugs = length(status), genes = length(genes), rank = rank),
  file.path(output_dir, 'run_status.tsv'), sep = '\t')
message('COMPLETE: ', output_dir)
