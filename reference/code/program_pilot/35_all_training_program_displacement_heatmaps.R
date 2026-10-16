#!/usr/bin/env Rscript
# Display every original fitting observation together, without the script-33 signal screen.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
})
source('/work/src/R/lib_heatmap.R')
source('/work/src/program_pilot/response_models.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps'
families <- c('ERBB', 'MEK')
z_color_limit <- 3
for (path in file.path(root, c('data', 'figures'))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(warn = 1, error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading complete training endpoints')
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
labels <- as.data.table(qread(panel_path)$sensitivity)
sources <- c(panel_path, '/work/src/program_pilot/35_all_training_program_displacement_heatmaps.R',
  '/work/src/R/lib_heatmap.R', '/work/src/program_pilot/response_models.R')
audits <- index <- list()

for (family_name in families) {
  parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family_name, '_holdout')
  paths <- file.path(parent, c('settings.qs', 'data/prepared_scores.qs', 'models/program_bank.qs',
    'models/functional_decoder.qs'))
  decoder_path <- file.path('/work/results/program_pilot/33_balanced_cell_and_dose_contrast_decoder/models',
    paste0(family_name, '_balanced_decoder.qs'))
  sources <- c(sources, paths, decoder_path)
  settings <- qread(paths[1L])
  prepared <- qread(paths[2L])
  bank <- qread(paths[3L])
  old <- qread(paths[4L])
  decoder <- qread(decoder_path)
  stopifnot(identical(prepared$signature, bank$signature), identical(bank$signature, old$signature),
    identical(bank$signature, decoder$source_signature), length(colnames(bank$weights)) == 40L,
    identical(prepared$metadata$observation_id, rownames(prepared$measured)),
    !any(old$training_ids %in% settings$held_ids), !any(old$training_drugs %in% settings$excluded_drugs))
  training <- copy(labels[depmap_id %in% old$training_ids & drug %in% old$training_drugs &
    is.finite(prism_abundance_pct) & prism_abundance_pct > 0])
  training[, `:=`(measured_pct = prism_abundance_pct, condition = paste(drug, dose, sep = '__'))]
  setorder(training, drug, dose, recipient, observation_id)
  stopifnot(!anyDuplicated(training$observation_id), uniqueN(training$depmap_id) == 26L,
    nrow(training) == if (family_name == 'ERBB') 6260L else 6572L)
  d <- prepared$measured[match(training$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  programs <- colnames(d)
  stopifnot(identical(rownames(d), training$observation_id), identical(programs, colnames(bank$weights)),
    all(is.finite(d)))
  b <- sweep(sweep(prepared$baseline[training$depmap_id, , drop = FALSE],
    2L, decoder$baseline_mean, '-'), 2L, decoder$baseline_sd, '/')
  training[, predicted_pct := 100 * 2^predict_ridge(decoder$model, cbind(d, d * b))]
  mean_d <- colMeans(d)
  sd_d <- apply(d, 2L, sd)
  stopifnot(all(is.finite(sd_d)), all(sd_d > 0), all(is.finite(training$predicted_pct)))
  z <- sweep(sweep(d, 2L, mean_d, '-'), 2L, sd_d, '/')

  # All columns are clustered together using molecular data alone. Responses do
  # not set the order. Centering is within each program across the full training set.
  status(paste(family_name, '| clustering all', nrow(d), 'training profiles together'))
  row_tree <- hclust(as.dist(sqrt(pmax(2 * (1 - cor(d)), 0))), method = 'ward.D2')
  column_tree <- hclust(dist(z), method = 'ward.D2')
  ro <- row_tree$order
  co <- column_tree$order
  ordered <- training[co]
  row_labels <- paste0(programs[ro], ' (', as.integer(colSums(bank$weights > 0))[ro], ' genes)')
  annotations <- data.frame(
    Drug = ordered$drug,
    Dose = paste0(ordered$dose, ' µM'),
    Cell = ordered$recipient,
    Sensitivity = ifelse(ordered$measured_pct <= 50, 'Sensitive <=50%', 'Resistant >50%'),
    row.names = ordered$observation_id, check.names = FALSE)
  bar_max <- ceiling(max(c(ordered$measured_pct, ordered$predicted_pct)) / 100) * 100
  bar_ticks <- function(maximum) unique(c(0, 50, 100, maximum))
  palettes <- list(
    Drug = setNames(grDevices::hcl.colors(uniqueN(ordered$drug), 'Dark 3'), sort(unique(ordered$drug))),
    Cell = setNames(grDevices::hcl.colors(uniqueN(ordered$recipient), 'Dynamic'), sort(unique(ordered$recipient))),
    Dose = c('0.05 µM' = '#DCEAF2', '0.5 µM' = '#6BAED6', '5 µM' = '#08519C'),
    Sensitivity = c('Sensitive <=50%' = '#D55E00', 'Resistant >50%' = '#0072B2'))
  matrices <- list(row_zscore = t(z[co, ro, drop = FALSE]), raw = t(d[co, ro, drop = FALSE]))
  status(paste(family_name, '| exporting combined training heatmap'))
  for (representation in names(matrices)) {
    mat <- matrices[[representation]]
    stopifnot(identical(colnames(mat), rownames(annotations)), identical(rownames(mat), programs[ro]))
    limit <- if (representation == 'row_zscore') z_color_limit else max(abs(mat))
    hm_name <- if (representation == 'row_zscore') 'Program z-score\n(colors saturate at +/-3)' else 'Signed displacement'
    stem <- paste0(family_name, '_all_training_', representation)
    pdf_path <- file.path(root, 'figures', paste0(stem, '.pdf'))
    png_path <- file.path(root, 'figures', paste0(stem, '.png'))
    # Measure annotations/layout afresh on each export device, avoiding stale
    # viewport units when moving between PDF and raster graphics.
    for (format in c('pdf', 'png')) {
      if (format == 'pdf') {
        grDevices::cairo_pdf(pdf_path, width = 20, height = 12)
      } else {
        png(png_path, width = 20, height = 12, units = 'in', res = 180, type = 'cairo')
      }
      grid::grid.newpage()
      annotation_args <- list(
        `Measured %vehicle` = ComplexHeatmap::anno_barplot(ordered$measured_pct,
          ylim = c(0, bar_max), baseline = 0, height = grid::unit(16, 'mm'), border = FALSE,
          gp = grid::gpar(fill = ifelse(ordered$measured_pct <= 50, '#D55E00', '#0072B2'), col = NA),
          axis_param = list(at = bar_ticks(bar_max), gp = grid::gpar(fontsize = 7))),
        `Predicted %vehicle` = ComplexHeatmap::anno_barplot(ordered$predicted_pct,
          ylim = c(0, bar_max), baseline = 0, height = grid::unit(16, 'mm'), border = FALSE,
          gp = grid::gpar(fill = '#777777', col = NA),
          axis_param = list(at = bar_ticks(bar_max), gp = grid::gpar(fontsize = 7))),
        simple_anno_size = grid::unit(2.5, 'mm'), annotation_name_gp = grid::gpar(fontsize = 8),
        annotation_legend_param = list(
          Drug = list(ncol = 3, labels_gp = grid::gpar(fontsize = 6)),
          Cell = list(ncol = 2, labels_gp = grid::gpar(fontsize = 7))))
      h <- draw_heatmap(mat, col_anno_df = annotations, col_anno_cols = colnames(annotations),
      ann_top_args = annotation_args, palette_overrides = palettes, scale_rows = FALSE,
      value_range = c(-limit, 0, limit), value_colors = c('#2166AC', 'white', '#B2182B'),
      hm_name = hm_name,
      heatmap_legend_param = if (representation == 'row_zscore') list(at = c(-3, -1.5, 0, 1.5, 3)) else list(),
      cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
      show_column_names = FALSE, row_labels = row_labels, row_names_gp = grid::gpar(fontsize = 9),
      title = paste0('ALL TRAINING DATA | ', family_name, ' bank | 40 programs × ', nrow(training),
        ' observations | 26 cells | ', uniqueN(training$drug), ' drugs | all doses'),
      column_title_gp = grid::gpar(fontsize = 12))
      ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right', newpage = FALSE)
      dev.off()
    }
    index[[stem]] <- data.table(family = family_name, representation, pdf = pdf_path, png = png_path,
      observations = nrow(training), cells = uniqueN(training$depmap_id), drugs = uniqueN(training$drug))
  }
  qsave(list(displacement = d, program_zscore = z, annotations = training,
    program_mean = mean_d, program_sd = sd_d, row_tree = row_tree, column_tree = column_tree,
    row_order = ro, column_order = co), file.path(root, 'data', paste0(family_name, '_all_training_inputs.qs')))
  fwrite(ordered[, display_order := seq_len(.N)],
    file.path(root, 'data', paste0(family_name, '_ordered_annotations.tsv')), sep = '\t')
  fwrite(data.table(program = programs, genes = as.integer(colSums(bank$weights > 0)),
    training_mean = mean_d, training_sd = sd_d, display_order = match(seq_along(programs), ro)),
    file.path(root, 'data', paste0(family_name, '_program_scaling_and_order.tsv')), sep = '\t')
  audits[[family_name]] <- data.table(family = family_name, observations = nrow(training),
    cells = uniqueN(training$depmap_id), drugs = uniqueN(training$drug), contexts = uniqueN(training$condition),
    sensitive = sum(training$measured_pct <= 50), resistant = sum(training$measured_pct > 50),
    molecular_programs = ncol(d), held_cells = sum(training$depmap_id %in% settings$held_ids),
    excluded_drugs = sum(training$drug %in% settings$excluded_drugs),
    displayed_zscore_saturation_pct = 100 * mean(abs(z) > z_color_limit))
}
fwrite(rbindlist(audits), file.path(root, 'data', 'all_training_counts.tsv'), sep = '\t')
figure_index <- rbindlist(index)
fwrite(figure_index, file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
fwrite(data.table(path = unique(sources), sha256 = vapply(unique(sources), digest::digest,
  character(1), file = TRUE, algo = 'sha256')), file.path(root, 'data', 'input_hashes.tsv'), sep = '\t')
revision <- system2('git', c('-C', '/work', 'rev-parse', 'HEAD'), stdout = TRUE)
qsave(list(families = families, all_original_fitting_observations = TRUE, no_signal_screen = TRUE,
  rows = '40 signed hard-mean molecular displacement programs', columns = 'all cells/drugs/doses clustered together',
  clustering = 'Ward.D2 on per-program standardized displacements; rows by Pearson program correlation',
  row_zscore_color_saturation = z_color_limit, raw_color_saturation = FALSE,
  outcomes_used_for_order = FALSE, models_refitted = FALSE, git_revision = revision), file.path(root, 'settings.qs'))
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(c('ALL ORIGINAL TRAINING DATA, COMBINED', '',
  'Start with figures/ERBB_all_training_row_zscore.pdf: the bank whose screened training r was 0.165.',
  'All 6,260 original fitting observations are together: 26 fitting cells, 90 drugs, all available doses.',
  'No informative-context filter is applied. No held cells or excluded-family drugs are included.',
  'MEK_all_training_row_zscore.pdf provides its separate 40-program bank (6,572 observations, 94 drugs).',
  'The two banks have different gene memberships; their program numbers are not interchangeable.',
  'Each matrix shows actual expression displacement hard-cluster means, not functional loadings or predictions.',
  'Rows and columns are clustered using molecular data alone; neither order uses sensitivity.',
  'Row z-score centers/scales each program across its whole training set. Color saturates at ±3 SD, labeled in legend.',
  'Raw PDFs retain signed displacement values and full color range with exactly the same row/column order.',
  'Top bars: measured PRISM abundance (% vehicle), orange <=50%, blue >50%. Grey bars: script-33 prediction.',
  'Measured and predicted bars share identical axes/heights; annotations show drug, dose, cell and sensitivity.',
  'Use the complete combined view to inspect association across drugs/doses/cells.',
  'It does not alone establish discrimination across cells receiving the SAME drug and dose.',
  'Full source matrices and dendrograms are in data/*_all_training_inputs.qs; ordered identities are exported as TSV.'),
  file.path(root, 'READING_GUIDE.txt'))
stopifnot(nrow(figure_index) == 4L, all(file.exists(figure_index$pdf)), all(file.info(figure_index$pdf)$size > 0))
status('COMPLETE')
