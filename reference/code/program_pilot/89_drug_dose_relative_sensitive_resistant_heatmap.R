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
figures <- file.path(root, 'figures/clustered_program_sensitivity/drug_dose_relative')
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
points_path <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/data/observed_program_viability_points.tsv.gz'
selection_path <- file.path(root, 'data/selected_doses.tsv')
programs <- sprintf('P%02d', 1:40)
drug_names <- c('5-Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan',
  'Cytarabine (hydrochloride)', 'Lapatinib ditosylate', 'Topotecan (hydrochloride)',
  'Vandetanib', 'Vismodegib')
drug_labels <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
  'Lapatinib', 'Topotecan', 'Vandetanib', 'Vismodegib')
drug_colors <- setNames(c('#4E79A7', '#F28E2B', '#59A14F', '#E15759', '#76B7B2',
  '#B07AA1', '#9C755F', '#B89800', '#79706E'), drug_labels)

selected <- fread(selection_path)
points <- merge(fread(points_path), selected[, .(drug, dose)], by = c('drug', 'dose'))
profiles <- unique(points[, .(depmap_id, cell_line, drug, dose, contributing_plates,
  observed_log2_abundance, observed_percent_vehicle)])
setorder(profiles, drug, depmap_id)
profiles[, profile_id := paste(depmap_id, drug, dose, sep = '|')]
stopifnot(nrow(profiles) == 256L, uniqueN(profiles$depmap_id) == 30L,
  !anyDuplicated(profiles$profile_id))
wide <- dcast(points, depmap_id + drug + dose ~ program, value.var = 'observed_displacement')
wide[, profile_id := paste(depmap_id, drug, dose, sep = '|')]
wide <- wide[match(profiles$profile_id, profile_id)]
d <- as.matrix(wide[, ..programs])
rownames(d) <- profiles$profile_id
stopifnot(all(is.finite(d)), identical(rownames(d), profiles$profile_id),
  all(is.finite(profiles$observed_percent_vehicle)))

# Exactly the scatter transformation: mean 0 / SD 1 separately per drug/dose/program.
z_relative <- d
for (dr in unique(profiles$drug)) {
  ii <- which(profiles$drug == dr)
  mu <- colMeans(d[ii, , drop = FALSE])
  sc <- apply(d[ii, , drop = FALSE], 2L, sd)
  stopifnot(all(sc > 1e-12))
  z_relative[ii, ] <- sweep(sweep(d[ii, , drop = FALSE], 2L, mu), 2L, sc, '/')
}
stopifnot(all(is.finite(z_relative)))
views <- list(drug_dose_relative_sensitivity_split = t(z_relative))
bar_max <- max(100, ceiling(max(profiles$observed_percent_vehicle) / 25) * 25)
book <- file.path(figures, 'clustered_40_programs_drug_dose_relative_sensitive_vs_resistant.pdf')
cairo_pdf(book, width = 15, height = 9.5, onefile = TRUE)
book_device <- dev.cur()
orders <- list()
for (view in names(views)) {
  m <- views[[view]]
  row_tree <- hclust(as.dist(sqrt(pmax(2 * (1 - cor(t(m))), 0))), method = 'ward.D2')
  ro <- row_tree$order
  row_cluster <- cutree(row_tree, k = 5L)
  sensitivity_group <- ifelse(profiles$observed_percent_vehicle <= 50, 'Sensitive (<=50% vehicle)', 'Resistant (>50% vehicle)')
  group_levels <- c('Sensitive (<=50% vehicle)', 'Resistant (>50% vehicle)')
  group_indices <- lapply(group_levels, function(g) which(sensitivity_group == g))
  stopifnot(all(lengths(group_indices) >= 2L))
  column_trees <- lapply(group_indices, function(ii) hclust(dist(t(m[, ii, drop = FALSE])), method = 'ward.D2'))
  co <- unlist(lapply(seq_along(group_indices), function(j) group_indices[[j]][column_trees[[j]]$order]))
  column_cluster <- match(sensitivity_group, group_levels)
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
  column_split <- factor(sensitivity_group[co], levels = group_levels)
  heading <- 'Drug/dose-relative program variation | sensitivity split, molecular clustering within each group'
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
      title = paste0(heading, '\n40 programs x 256 measured profiles | 30 cell lines | nine selected drug/dose contexts'))
    drawn <- ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
    actual_columns <- unlist(ComplexHeatmap::column_order(drawn), use.names = FALSE)
    actual_rows <- unlist(ComplexHeatmap::row_order(drawn), use.names = FALSE)
    stopifnot(identical(as.integer(actual_columns), seq_len(ncol(mat))),
      identical(as.integer(actual_rows), seq_len(nrow(mat))))
  }
  dev.set(book_device)
  render()
  png(file.path(figures, paste0(view, '.png')), width = 15, height = 9.5,
    units = 'in', res = 160, type = 'cairo')
  render()
  dev.off()
  ordered[, `:=`(display_column = .I, molecular_cluster = column_cluster[co], view = view)]
  orders[[view]] <- ordered
  qsave(list(raw_displacements = d, display_matrix = mat, row_tree = row_tree, column_trees = column_trees,
    row_order = ro, column_order = co, row_cluster = row_cluster, column_cluster = column_cluster),
    file.path(figures, paste0('heatmap_', view, '_inputs.qs')))
  fwrite(data.table(program = programs[ro], display_row = seq_along(ro), cluster = row_cluster[ro], view),
    file.path(figures, paste0('heatmap_', view, '_program_order.tsv')), sep = '\t')
  cat('Saved ', view, '\n', sep = '')
}
dev.set(book_device)
dev.off()
fwrite(rbindlist(orders), file.path(figures, 'heatmap_sensitive_vs_resistant_column_annotations.tsv'), sep = '\t')
writeLines(c('Each program: (displacement - own drug/dose mean) / own drug/dose sample SD, exactly as the relative scatterplots.',
  '40 programs, 256 measured profiles, 30 cell lines, nine fixed selected drug/dose contexts. No crizotinib in this requested v89 view.',
  'Observed PRISM first splits columns at <=50% vehicle; molecular values then cluster columns independently within each group. Rows share a global molecular clustering order.',
  'Top bars retain observed abundance percent of vehicle; lower means greater sensitivity. Orange <=50%; blue >50%.',
  'Heatmap colors saturated at +/-3 SD for display only. No model fit or batch correction.'),
  file.path(figures, 'SENSITIVITY_SPLIT_DISPLAY_DEFINITION.txt'))
cat('Complete: per-drug/dose z-score heatmap matching the scatter transformation.
')
