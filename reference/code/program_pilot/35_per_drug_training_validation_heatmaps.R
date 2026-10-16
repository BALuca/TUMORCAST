#!/usr/bin/env Rscript
# Per-drug views of the existing ERBB-bank displacement data; no model refitting.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps'
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout'
input_paths <- c(
  training = file.path(root, 'data/ERBB_all_training_inputs.qs'),
  scores = file.path(parent, 'data/prepared_scores.qs'),
  bank = file.path(parent, 'models/program_bank.qs'),
  labels = '/work/results/program_pilot/40_constrained_program_sensitivity_classifier/data/evaluated_scores.tsv')
stopifnot(all(file.exists(input_paths)))
training <- qread(input_paths['training'])
prepared <- qread(input_paths['scores'])
bank <- qread(input_paths['bank'])
labels <- fread(input_paths['labels'])
labels <- labels[cohort %in% c('Training fit', 'Original held cells')]
labels[, split := factor(ifelse(cohort == 'Training fit', 'Training', 'Validation'),
  levels = c('Training', 'Validation'))]
setorder(labels, drug, split, measured_pct, recipient, dose, observation_id)
stopifnot(!anyDuplicated(labels$observation_id), nrow(labels) == 8489L,
  all(is.finite(labels$measured_pct)), uniqueN(labels$drug) == 90L,
  uniqueN(labels[split == 'Training', depmap_id]) == 26L,
  uniqueN(labels[split == 'Validation', depmap_id]) == 9L,
  !any(labels[split == 'Training', depmap_id] %in% labels[split == 'Validation', depmap_id]),
  identical(prepared$signature, bank$signature),
  identical(prepared$metadata$observation_id, rownames(prepared$measured)))
d <- prepared$measured[match(labels$observation_id, rownames(prepared$measured)), , drop = FALSE]
stopifnot(identical(rownames(d), labels$observation_id), ncol(d) == 40L, all(is.finite(d)),
  identical(colnames(d), colnames(training$displacement)),
  isTRUE(all.equal(d[rownames(training$displacement), ], training$displacement, tolerance = 1e-12)))
# Preserve training-derived scaling and program order across both cohorts and all drugs.
z <- sweep(sweep(d, 2L, training$program_mean, '-'), 2L, training$program_sd, '/')
ro <- training$row_order
row_labels <- paste0(colnames(d)[ro], ' (', as.integer(colSums(bank$weights > 0))[ro], ' genes)')
dose_colors <- c('0.05 µM' = '#DCEAF2', '0.5 µM' = '#6BAED6', '5 µM' = '#08519C')
stopifnot(all(paste0(labels$dose, ' µM') %in% names(dose_colors)))
drugs <- sort(unique(labels$drug))
slugs <- gsub('[^A-Za-z0-9_-]+', '_', drugs)
stopifnot(!anyDuplicated(slugs))

# Reused for individual PDFs, complete collections and previews. Column trees
# are fitted separately inside Training and Validation and never use outcomes.
render_drug <- function(drug_labels, drug_name, clustered) {
  mat <- t(z[drug_labels$observation_id, ro, drop = FALSE])
  sensitivity <- ifelse(drug_labels$measured_pct <= 50, 'Sensitive <=50%', 'Resistant >50%')
  annotations <- data.frame(Dose = paste0(drug_labels$dose, ' µM'),
    Sensitivity = sensitivity, row.names = drug_labels$observation_id, check.names = FALSE)
  stopifnot(identical(colnames(mat), rownames(annotations)))
  bar_max <- max(100, ceiling(max(drug_labels$measured_pct) / 25) * 25)
  counts <- drug_labels[, .(n = .N, cells = uniqueN(depmap_id), sensitive = sum(measured_pct <= 50)), by = split]
  split_titles <- paste0(counts$split, ': ', counts$cells, ' cells | ', counts$n,
    ' observations | ', counts$sensitive, ' sensitive')
  h <- draw_heatmap(mat, col_anno_df = annotations, col_anno_cols = colnames(annotations),
    palette_overrides = list(Dose = dose_colors,
      Sensitivity = c('Sensitive <=50%' = '#D55E00', 'Resistant >50%' = '#0072B2')),
    ann_top_args = list(
      `Measured %vehicle` = ComplexHeatmap::anno_barplot(drug_labels$measured_pct,
        ylim = c(0, bar_max), baseline = 0, height = grid::unit(22, 'mm'), border = FALSE,
        gp = grid::gpar(fill = ifelse(drug_labels$measured_pct <= 50, '#D55E00', '#0072B2'), col = NA),
        axis_param = list(at = unique(c(0, 50, 100, bar_max)), gp = grid::gpar(fontsize = 7))),
      simple_anno_size = grid::unit(2.5, 'mm'), annotation_name_gp = grid::gpar(fontsize = 8)),
    scale_rows = FALSE, value_range = c(-3, 0, 3), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'Displacement\n(training SD; capped +/-3)',
    heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)),
    cluster_rows = FALSE, row_labels = row_labels, show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_split = drug_labels$split, cluster_column_slices = FALSE,
    cluster_columns = clustered, clustering_distance_columns = 'euclidean',
    clustering_method_columns = 'ward.D2', column_dend_reorder = FALSE,
    column_gap = grid::unit(6, 'mm'), show_column_names = TRUE,
    column_labels = paste0(drug_labels$recipient, ' | ', drug_labels$dose, ' µM'),
    column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 6.5),
    column_names_max_height = grid::unit(40, 'mm'),
    title = split_titles, column_title_gp = grid::gpar(fontsize = 10))
  method <- if (clustered) 'Ward.D2 clustering separately within each cohort' else
    'Ordered from most sensitive to most resistant within each cohort'
  drawn <- ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right',
    column_title = paste0(drug_name, ' | 40 program displacements\n', method,
      ' | all measured doses'), column_title_gp = grid::gpar(fontsize = 12))
  orders <- ComplexHeatmap::column_order(drawn)
  stopifnot(length(orders) == 2L)
  for (j in seq_along(orders)) {
    block <- drug_labels[orders[[j]]]
    stopifnot(all(as.character(block$split) == levels(drug_labels$split)[j]))
    if (!clustered) stopifnot(all(diff(block$measured_pct) >= 0))
    if (clustered) {
      base <- which(drug_labels$split == levels(drug_labels$split)[j])
      expected <- if (length(base) == 1L) base else
        base[hclust(dist(z[drug_labels$observation_id[base], , drop = FALSE]), method = 'ward.D2')$order]
      stopifnot(identical(as.integer(orders[[j]]), as.integer(expected)))
    }
  }
  drug_labels[unlist(orders, use.names = FALSE), .(observation_id, recipient, depmap_id,
    drug, dose, split, measured_pct)][, `:=`(view = if (clustered) 'clustered' else 'sensitivity_ordered',
      display_order = seq_len(.N))]
}

# Open two collections; each drug's own two-page PDF is closed immediately.
collection_paths <- file.path(root, 'figures', c('ERBB_per_drug_sensitivity_ordered.pdf',
  'ERBB_per_drug_clustered_within_cohort.pdf'))
grDevices::cairo_pdf(collection_paths[1L], width = 15, height = 10, onefile = TRUE)
ordered_device <- dev.cur()
grDevices::cairo_pdf(collection_paths[2L], width = 15, height = 10, onefile = TRUE)
clustered_device <- dev.cur()
index <- orders <- list()
for (i in seq_along(drugs)) {
  drug_name <- drugs[i]
  subset <- labels[drug == drug_name]
  pdf_path <- file.path(root, 'figures', paste0('ERBB_per_drug_', slugs[i], '_training_validation.pdf'))
  grDevices::cairo_pdf(pdf_path, width = 15, height = 10, onefile = TRUE)
  individual_device <- dev.cur()
  for (clustered in c(FALSE, TRUE)) {
    dev.set(individual_device)
    order_table <- render_drug(subset, drug_name, clustered)
    orders[[length(orders) + 1L]] <- order_table
    dev.set(if (clustered) clustered_device else ordered_device)
    collection_order <- render_drug(subset, drug_name, clustered)
    stopifnot(identical(collection_order, order_table))
  }
  dev.off(individual_device)
  index[[i]] <- subset[, .(drug = drug_name, collection_page = i, individual_pdf = pdf_path,
    training_n = sum(split == 'Training'), validation_n = sum(split == 'Validation'),
    training_sensitive = sum(split == 'Training' & measured_pct <= 50),
    validation_sensitive = sum(split == 'Validation' & measured_pct <= 50))]
  fwrite(rbindlist(index), file.path(root, 'PER_DRUG_FIGURE_INDEX.tsv'), sep = '\t')
  message(format(Sys.time()), ' | ', i, '/', length(drugs), ' | ', drug_name, ' | saved ', pdf_path)
}
dev.off(ordered_device)
dev.off(clustered_device)
fwrite(rbindlist(orders), file.path(root, 'data/ERBB_per_drug_displayed_columns.tsv'), sep = '\t')

# Inspectable raster previews, produced from exactly the same draw function.
for (drug_name in c('bortezomib', 'everolimus')) {
  stopifnot(drug_name %in% drugs)
  for (clustered in c(FALSE, TRUE)) {
    png_path <- file.path(root, 'figures', paste0('ERBB_per_drug_', drug_name,
      if (clustered) '_clustered.png' else '_sensitivity_ordered.png'))
    png(png_path, width = 15, height = 10, units = 'in', res = 120, type = 'cairo')
    render_drug(labels[drug == drug_name], drug_name, clustered)
    dev.off()
  }
}
fwrite(data.table(path = c(input_paths, '/work/src/program_pilot/35_per_drug_training_validation_heatmaps.R'),
  sha256 = vapply(c(input_paths, '/work/src/program_pilot/35_per_drug_training_validation_heatmaps.R'),
    digest::digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/ERBB_per_drug_input_hashes.tsv'), sep = '\t')
writeLines(capture.output(sessionInfo()), file.path(root, 'per_drug_sessionInfo.txt'))
writeLines(c('90 familiar training drugs; frozen ERBB 40-program bank.',
  'Training: 26 cells / 6260 observations. Validation: original 9 held cells / 2229 observations.',
  'Validation labels are retrospective; used only for display, never program scaling or clustering.',
  'All doses retained. Columns are cell line x dose; bars show measured PRISM abundance (% vehicle).',
  'Orange <=50%: sensitive; blue >50%: resistant. Lower bars mean greater sensitivity.',
  'Individual PDFs: page 1 sensitivity-ordered, page 2 Ward.D2 clustered independently within each cohort.',
  'Two collection PDFs: one page per drug, alphabetical; page lookup in PER_DRUG_FIGURE_INDEX.tsv.',
  'Program rows and scaling are identical to script 35; color saturates at +/-3 training SD.',
  'All drug/cohort column orders and identifiers are exported in data/ERBB_per_drug_displayed_columns.tsv.',
  paste('Git revision:', system2('git', c('-C', '/work', 'rev-parse', 'HEAD'), stdout = TRUE))),
  file.path(root, 'PER_DRUG_README.txt'))
stopifnot(length(index) == 90L, all(file.info(c(collection_paths, rbindlist(index)$individual_pdf))$size > 0))
message('COMPLETE: 90 individual two-page PDFs and both 90-page collections.')
