#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
# Reuse the earlier sensitivity-annotated heatmap style for the exact v89 data.
# Molecular values alone determine row/column clusters; outcomes annotate them.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
root <- '/work/results/program_pilot/89_shared_program_interaction_decoder_mock'
figures <- file.path(root, 'figures/clustered_program_sensitivity/lapatinib_crizotinib_5uM')
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
points_path <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/data/observed_program_viability_points.tsv.gz'
selection_path <- file.path(root, 'data/selected_doses.tsv')
programs <- sprintf('P%02d', 1:40)
drug_names <- c('Lapatinib ditosylate', 'crizotinib')
drug_labels <- c('Lapatinib', 'Crizotinib')
drug_colors <- c(Lapatinib = '#B07AA1', Crizotinib = '#59A14F')
selected <- data.table(drug = drug_names, dose = c(5, 5))
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
lap <- as.data.table(archive$prep$m)[drug == 'Lapatinib ditosylate' & dose == 5]
outcomes <- fread('/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/data/individual_training_outcomes.tsv')[drug == 'Lapatinib ditosylate' & dose == 5 & source == 'prism_secondary', .(depmap_id, drug, dose, observed_log2_abundance = y)]
lap <- merge(lap, outcomes, by = c('depmap_id', 'drug', 'dose'))
lap <- merge(lap, as.data.table(archive$meta)[, .(depmap_id = ModelID, cell_line = CellLineName)], by = 'depmap_id')
lap[, observed_percent_vehicle := 100 * 2^observed_log2_abundance]
lap_points <- melt(lap, id.vars = c('depmap_id', 'drug', 'dose', 'cell_line', 'contributing_plates', 'observed_log2_abundance', 'observed_percent_vehicle'), measure.vars = programs, variable.name = 'program', value.name = 'observed_displacement', variable.factor = FALSE)
criz_points <- fread('/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/figures/observed_program_viability/crizotinib/observed_program_viability_points.tsv.gz')[dose == 5]
points <- rbindlist(list(lap_points, criz_points), use.names = TRUE, fill = TRUE)
profiles <- unique(points[, .(depmap_id, cell_line, drug, dose, contributing_plates,
  observed_log2_abundance, observed_percent_vehicle)])
setorder(profiles, drug, depmap_id)
profiles[, profile_id := paste(depmap_id, drug, dose, sep = '|')]
stopifnot(nrow(profiles) == 60L, uniqueN(profiles$depmap_id) == 30L,
  !anyDuplicated(profiles$profile_id))
wide <- dcast(points, depmap_id + drug + dose ~ program, value.var = 'observed_displacement')
wide[, profile_id := paste(depmap_id, drug, dose, sep = '|')]
wide <- wide[match(profiles$profile_id, profile_id)]
d <- as.matrix(wide[, ..programs])
rownames(d) <- profiles$profile_id
stopifnot(all(is.finite(d)), identical(rownames(d), profiles$profile_id),
  all(is.finite(profiles$observed_percent_vehicle)))

# The second display isolates molecular cell variation within each treatment.
# Both scalings are descriptive display transforms, not corrected RNA inputs.
within <- d
for (drug in unique(profiles$drug)) {
  ii <- which(profiles$drug == drug)
  within[ii, ] <- sweep(d[ii, , drop = FALSE], 2L, colMeans(d[ii, , drop = FALSE]))
}
z_global <- sweep(sweep(d, 2L, colMeans(d)), 2L, apply(d, 2L, sd), '/')
z_within <- sweep(within, 2L, apply(within, 2L, sd), '/')
stopifnot(all(is.finite(z_global)), all(is.finite(z_within)))
views <- list(global_displacement_patterns = t(z_global), within_drug_displacement_variation = t(z_within))
bar_max <- max(100, ceiling(max(profiles$observed_percent_vehicle) / 25) * 25)
book <- file.path(figures, 'clustered_40_programs_with_measured_sensitivity.pdf')
cairo_pdf(book, width = 11.5, height = 9.5, onefile = TRUE)
book_device <- dev.cur()
orders <- list()
for (view in names(views)) {
  m <- views[[view]]
  row_tree <- hclust(as.dist(sqrt(pmax(2 * (1 - cor(t(m))), 0))), method = 'ward.D2')
  column_tree <- hclust(dist(t(m)), method = 'ward.D2')
  ro <- row_tree$order
  co <- column_tree$order
  row_cluster <- cutree(row_tree, k = 5L)
  column_cluster <- cutree(column_tree, k = 6L)
  ordered <- profiles[co]
  mat <- m[ro, co, drop = FALSE]
  annotations <- data.frame(
    Drug = drug_labels[match(ordered$drug, drug_names)],
    Dose = paste0(ordered$dose, ' µM'),
    Sensitivity = ifelse(ordered$observed_percent_vehicle <= 50, '<=50% vehicle', '>50% vehicle'),
    row.names = ordered$profile_id, check.names = FALSE)
  stopifnot(identical(colnames(mat), rownames(annotations)),
    identical(rownames(mat), programs[ro]))
  row_split <- factor(paste0('R', row_cluster[ro]), levels = unique(paste0('R', row_cluster[ro])))
  column_split <- factor(paste0('C', column_cluster[co]), levels = unique(paste0('C', column_cluster[co])))
  heading <- if (view == 'global_displacement_patterns') {
    'Observed program patterns | row z-scores across all profiles'
  } else {
    'Within-drug program variation | drug/dose means removed; pooled within-drug SD units'
  }
  render <- function() {
    h <- draw_heatmap(mat, col_anno_df = annotations, col_anno_cols = names(annotations),
      palette_overrides = list(Drug = drug_colors,
        Dose = c('0.05 µM' = '#DCEAF2', '0.5 µM' = '#6BAED6', '5 µM' = '#08519C'),
        Sensitivity = c('<=50% vehicle' = '#D55E00', '>50% vehicle' = '#0072B2')),
      ann_top_args = list(
        `Measured %vehicle` = ComplexHeatmap::anno_barplot(ordered$observed_percent_vehicle,
          ylim = c(0, bar_max), baseline = 0, height = grid::unit(18, 'mm'), border = FALSE,
          gp = grid::gpar(fill = ifelse(ordered$observed_percent_vehicle <= 50, '#D55E00', '#0072B2'), col = NA),
          axis_param = list(at = unique(c(0, 50, 100, bar_max)), gp = grid::gpar(fontsize = 7))),
        simple_anno_size = grid::unit(2.5, 'mm'), annotation_name_gp = grid::gpar(fontsize = 8),
        annotation_legend_param = list(Drug = list(ncol = 2, labels_gp = grid::gpar(fontsize = 8)))),
      scale_rows = FALSE, value_range = c(-3, 0, 3), value_colors = c('#2166AC', 'white', '#B2182B'),
      hm_name = 'Program score\n(display SD units)',
      heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = row_split, column_split = column_split,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_gap = grid::unit(1.5, 'mm'), column_gap = grid::unit(2, 'mm'),
      show_row_names = TRUE, show_column_names = FALSE, row_names_gp = grid::gpar(fontsize = 8),
      row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 9),
      title = paste0(heading, '\n40 programs x 60 measured profiles | lapatinib and crizotinib | both 5 µM'))
    drawn <- ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
    actual_columns <- unlist(ComplexHeatmap::column_order(drawn), use.names = FALSE)
    actual_rows <- unlist(ComplexHeatmap::row_order(drawn), use.names = FALSE)
    stopifnot(identical(as.integer(actual_columns), seq_len(ncol(mat))),
      identical(as.integer(actual_rows), seq_len(nrow(mat))))
  }
  dev.set(book_device)
  render()
  png(file.path(figures, paste0(view, '.png')), width = 11.5, height = 9.5,
    units = 'in', res = 160, type = 'cairo')
  render()
  dev.off()
  ordered[, `:=`(display_column = .I, molecular_cluster = column_cluster[co], view = view)]
  orders[[view]] <- ordered
  qsave(list(raw_displacements = d, display_matrix = mat, row_tree = row_tree, column_tree = column_tree,
    row_order = ro, column_order = co, row_cluster = row_cluster, column_cluster = column_cluster),
    file.path(figures, paste0('heatmap_', view, '_inputs.qs')))
  fwrite(data.table(program = programs[ro], display_row = seq_along(ro), cluster = row_cluster[ro], view),
    file.path(figures, paste0('heatmap_', view, '_program_order.tsv')), sep = '\t')
  cat('Saved ', view, '\n', sep = '')
}
dev.set(book_device)
dev.off()
fwrite(rbindlist(orders), 