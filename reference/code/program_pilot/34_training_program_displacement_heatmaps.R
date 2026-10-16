#!/usr/bin/env Rscript
# Inspect exact frozen program displacements behind the script-33 training fit.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet)})
source('/work/src/R/lib_heatmap.R')
source('/work/src/program_pilot/response_models.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/34_training_program_displacement_heatmaps'
parent <- '/work/results/program_pilot/33_balanced_cell_and_dose_contrast_decoder'
families <- c('ERBB', 'MEK')
preview_drugs <- c('bortezomib', 'everolimus', 'irinotecan')
for (path in file.path(root, c('data', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading completed decoder predictions')
stopifnot(fread(file.path(parent, 'run_status.tsv'))$stage == 'COMPLETE')
predictions <- fread(file.path(parent, 'data/evaluated_predictions.tsv'))
sources <- c(file.path(parent, c('data/evaluated_predictions.tsv', 'settings.qs')),
  '/work/src/program_pilot/34_training_program_displacement_heatmaps.R',
  '/work/src/R/lib_heatmap.R', '/work/src/program_pilot/response_models.R')
index <- overview_figures <- verification <- list()

# Bars retain actual assay units and a shared range for measured/predicted values.
# Color marks <=50% vehicle. The 40-60% band describes fitting support only.
make_map <- function(mat, annotations, row_labels, title, hm_name, limit, overview = FALSE) {
  stopifnot(identical(colnames(mat), annotations$observation_id),
    identical(rownames(mat), names(row_labels)))
  df <- data.frame(Drug = annotations$drug, Dose = paste0(annotations$dose, ' µM'), Cell = annotations$cell,
    Sensitivity = ifelse(annotations$observed_pct <= 50, 'Sensitive ≤50%', 'Resistant >50%'),
    Support = ifelse(annotations$observed_pct <= 40, 'Clear sensitive ≤40%',
      ifelse(annotations$observed_pct >= 60, 'Clear resistant ≥60%', '40–60% buffer')),
    row.names = annotations$observation_id, check.names = FALSE)
  stopifnot(identical(colnames(mat), rownames(df)))
  bar_max <- ceiling(max(c(annotations$observed_pct, annotations$predicted_pct)) / 50) * 50
  bar_fill <- ifelse(annotations$observed_pct <= 50, '#D55E00', '#0072B2')
  args <- list(
    `Measured %vehicle` = ComplexHeatmap::anno_barplot(annotations$observed_pct,
      ylim = c(0, bar_max), baseline = 0, height = grid::unit(16, 'mm'), border = FALSE,
      gp = grid::gpar(fill = bar_fill, col = NA), axis_param = list(at = unique(c(0, 50, 100, bar_max)),
        gp = grid::gpar(fontsize = 7))),
    `Predicted %vehicle` = ComplexHeatmap::anno_barplot(annotations$predicted_pct,
      ylim = c(0, bar_max), baseline = 0, height = grid::unit(16, 'mm'), border = FALSE,
      gp = grid::gpar(fill = '#777777', col = NA), axis_param = list(at = unique(c(0, 50, 100, bar_max)),
        gp = grid::gpar(fontsize = 7))),
    simple_anno_size = grid::unit(2.5, 'mm'), annotation_name_gp = grid::gpar(fontsize = 8))
  show <- if (overview) colnames(df) else c('Sensitivity', 'Support')
  if (overview) args$annotation_legend_param <- list(
    Drug = list(ncol = 2, labels_gp = grid::gpar(fontsize = 7)),
    Cell = list(ncol = 2, labels_gp = grid::gpar(fontsize = 7)))
  draw_heatmap(mat, col_anno_df = df, col_anno_cols = show, ann_top_args = args,
    palette_overrides = list(
      Drug = setNames(grDevices::hcl.colors(length(unique(df$Drug)), 'Dark 3'), sort(unique(df$Drug))),
      Cell = setNames(grDevices::hcl.colors(length(unique(df$Cell)), 'Dynamic'), sort(unique(df$Cell))),
      Dose = c('0.05 µM' = '#DCEAF2', '0.5 µM' = '#6BAED6', '5 µM' = '#08519C'),
      Sensitivity = c('Sensitive ≤50%' = '#D55E00', 'Resistant >50%' = '#0072B2'),
      Support = c('Clear sensitive ≤40%' = '#D55E00', 'Clear resistant ≥60%' = '#0072B2', '40–60% buffer' = '#CCCCCC')),
    scale_rows = FALSE, value_range = c(-limit, 0, limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = hm_name, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = !overview,
    row_labels = unname(row_labels), column_labels = annotations$cell,
    row_names_gp = grid::gpar(fontsize = if (overview) 9 else 8),
    column_names_gp = grid::gpar(fontsize = 8), column_names_rot = 90,
    title = title, column_title_gp = grid::gpar(fontsize = 11))
}

for (family_name in families) {
  status(paste(family_name, 'clustering actual training program scores'))
  source_root <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family_name, '_holdout')
  paths <- c(file.path(source_root, c('data/prepared_scores.qs', 'models/program_bank.qs', 'models/functional_decoder.qs')),
    file.path(parent, 'models', paste0(family_name, '_balanced_decoder.qs')))
  sources <- c(sources, paths)
  prepared <- qread(paths[1L])
  bank <- qread(paths[2L])
  original <- qread(paths[3L])
  decoder <- qread(paths[4L])
  stopifnot(identical(decoder$source_signature, bank$signature),
    identical(prepared$signature, bank$signature), setequal(decoder$training_ids, original$training_ids))
  query <- predictions[family == family_name & cohort == 'Fitting cells' &
    model == 'Balanced cell + dose fit' & cell_context == TRUE]
  setorder(query, drug, dose, cell)
  stopifnot(!anyDuplicated(query$observation_id), uniqueN(query$depmap_id) == 26L)
  programs <- colnames(bank$weights)
  d <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  all_training <- prepared$metadata$depmap_id %in% original$training_ids &
    prepared$metadata$drug %in% original$training_drugs
  # Use exactly the available functional fitting observations, not unsupported RNA conditions.
  train_ids <- predictions[family == family_name & cohort == 'Fitting cells' & model == 'Balanced cell + dose fit', observation_id]
  # Full original fitting IDs are reconstructed from the same paired endpoint table.
  panel <- qread('/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
  fitting <- as.data.table(panel$sensitivity)[depmap_id %in% original$training_ids & drug %in% original$training_drugs &
    is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  td <- prepared$measured[match(fitting$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  sd_program <- apply(td, 2L, sd)
  stopifnot(length(programs) == 40L, all(sd_program > 0), all(is.finite(d)),
    identical(colnames(d), programs), identical(rownames(d), query$observation_id))
  b <- sweep(sweep(prepared$baseline[query$depmap_id, , drop = FALSE], 2L, decoder$baseline_mean, '-'),
    2L, decoder$baseline_sd, '/')
  replay_error <- max(abs(predict_ridge(decoder$model, cbind(d, d * b)) - query$predicted))
  stopifnot(replay_error < 1e-10)
  verification[[family_name]] <- data.table(family = family_name, programs = 40L, cells = 26L,
    observations = nrow(query), contexts = uniqueN(query$condition), prediction_replay_difference = replay_error,
    score_source = paths[1L])
  scaled <- sweep(d, 2L, sd_program, '/')
  row_order <- hclust(as.dist(sqrt(pmax(2 * (1 - cor(d)), 0))), method = 'ward.D2')$order
  column_order <- hclust(dist(scaled), method = 'ward.D2')$order
  row_labels <- setNames(paste0(programs[row_order], ' (', as.integer(colSums(bank$weights > 0))[row_order], ' genes)'),
    programs[row_order])
  save_table(data.table(program = programs, genes = as.integer(colSums(bank$weights > 0)),
    training_sd = sd_program, row_order = match(seq_along(programs), row_order)), paste0(family_name, '_programs_and_scaling'))
  save_table(data.table(observation_id = query$observation_id[column_order], order = seq_along(column_order)),
    paste0(family_name, '_overview_column_order'))
  save_table(query, paste0(family_name, '_annotations'))
  qsave(list(displacement = d, sd_scaled_displacement = scaled, annotations = query,
    program_order = row_order, overview_column_order = column_order, training_program_sd = sd_program),
    file.path(root, 'data', paste0(family_name, '_heatmap_inputs.qs')))
  save_table(melt(data.table(observation_id = query$observation_id, d), id.vars = 'observation_id',
    variable.name = 'program', value.name = 'signed_displacement'), paste0(family_name, '_program_displacements'))

  for (representation in c('raw', 'training_sd')) {
    mat <- t((if (representation == 'raw') d else scaled)[column_order, row_order, drop = FALSE])
    limit <- max(abs(mat))
    h <- make_map(mat, query[column_order], row_labels,
      paste(family_name, '| all informative training contexts |', ncol(mat), 'observations'),
      if (representation == 'raw') 'Signed D' else 'D / training SD', limit, TRUE)
    stem <- paste0(family_name, '_overview_', representation)
    grDevices::cairo_pdf(file.path(root, 'figures', paste0(stem, '.pdf')), width = 16, height = 10)
    ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
    dev.off()
    png(file.path(root, 'figures', paste0(stem, '.png')), width = 16, height = 10, units = 'in', res = 150)
    ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
    dev.off()
    overview_figures[[stem]] <- h
    index[[stem]] <- data.table(family = family_name, kind = 'Overview', figure = stem,
      pdf = file.path(root, 'figures', paste0(stem, '.pdf')), png = file.path(root, 'figures', paste0(stem, '.png')),
      caption = paste('All 40 programs; exact signed hard-mean measured displacements.',
        'Columns clustered by SD-scaled molecular profiles; rows by program correlation. Neither order uses outcomes.',
        if (representation == 'raw') 'Colors use raw displacement units.' else 'Each program divided by its original fitting SD, without centering.',
        'Bars show measured and balanced-decoder-predicted PRISM abundance; sensitivity cutoff 50% vehicle.',
        'Drug, dose, cell and clear/buffer support labels are synchronized. These are fitting observations.'))
  }
  status(paste(family_name, 'exporting every fixed drug/dose'))
  # Two synchronized views expose the shared treatment response and the recipient
  # deviations around it. They share the same molecular ordering and color scale.
  contexts <- unique(query[, .(drug, dose, condition)])[order(drug, dose)]
  context_maps <- list()
  context_order <- list()
  for (j in seq_len(nrow(contexts))) {
    ii <- which(query$condition == contexts$condition[j])
    co <- hclust(dist(scaled[ii, , drop = FALSE]), method = 'ward.D2')$order
    jj <- ii[co]
    absolute <- t(scaled[jj, row_order, drop = FALSE])
    deviation <- sweep(absolute, 1L, rowMeans(absolute), '-')
    stopifnot(identical(colnames(absolute), query$observation_id[jj]), identical(dimnames(absolute), dimnames(deviation)))
    limit <- max(abs(c(absolute, deviation)))
    h1 <- make_map(absolute, query[jj], row_labels, 'Signed program displacement', 'D / training SD', limit)
    h2 <- make_map(deviation, query[jj], row_labels, 'Recipient deviation from drug/dose mean', 'D deviation / SD', limit)
    h <- h1 + h2
    stem <- paste0(family_name, '_', gsub('[^a-z0-9]+', '_', contexts$drug[j]), '_', contexts$dose[j], 'uM')
    title <- paste(family_name, '|', contexts$drug[j], '|', contexts$dose[j], 'µM |', length(ii), 'training cells')
    grDevices::cairo_pdf(file.path(root, 'figures', paste0(stem, '.pdf')), width = 15, height = 9)
    ComplexHeatmap::draw(h, column_title = title, merge_legends = TRUE,
      heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
    dev.off()
    preview <- contexts$drug[j] %in% preview_drugs || j == 1L
    if (preview) {
      png(file.path(root, 'figures', paste0(stem, '.png')), width = 15, height = 9, units = 'in', res = 150)
      ComplexHeatmap::draw(h, column_title = title, merge_legends = TRUE,
        heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
      dev.off()
    }
    context_maps[[stem]] <- list(heatmap = h, title = title)
    context_order[[stem]] <- data.table(condition = contexts$condition[j], observation_id = query$observation_id[jj],
      cell = query$cell[jj], order = seq_along(jj))
    index[[stem]] <- data.table(family = family_name, kind = 'Fixed drug/dose', figure = stem,
      pdf = file.path(root, 'figures', paste0(stem, '.pdf')),
      png = if (preview) file.path(root, 'figures', paste0(stem, '.png')) else NA_character_,
      caption = paste(title, '| Left: signed D divided by training program SD; right: same scores minus context mean.',
        'Identical row/column order, measured/predicted abundance bars and 50% sensitivity annotation.',
        'Column clustering uses molecular data only. Every retained context is exported, without performance selection.'))
  }
  save_table(rbindlist(context_order), paste0(family_name, '_fixed_context_column_orders'))
  grDevices::cairo_pdf(file.path(root, paste0(family_name, '_all_fixed_drug_dose_heatmaps.pdf')),
    width = 15, height = 9, onefile = TRUE)
  for (item in context_maps) ComplexHeatmap::draw(item$heatmap, column_title = item$title, merge_legends = TRUE,
    heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
  dev.off()
}
status('assembling overviews and verifying exports')
grDevices::cairo_pdf(file.path(root, 'training_program_overviews.pdf'), width = 16, height = 10, onefile = TRUE)
for (h in overview_figures) ComplexHeatmap::draw(h, heatmap_legend_side = 'bottom', annotation_legend_side = 'right')
dev.off()
figure_index <- rbindlist(index)
fwrite(figure_index, file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
save_table(rbindlist(verification), 'score_and_prediction_replay')
sources <- c(sources, '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
save_table(data.table(path = unique(sources), sha256 = vapply(unique(sources), digest::digest,
  character(1), file = TRUE, algo = 'sha256')), 'input_hashes')
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
qsave(list(families = families, source_predictions = file.path(parent, 'data/evaluated_predictions.tsv'),
  row_clustering = 'program Pearson correlation; Ward.D2', column_clustering = 'SD-scaled measured D; Ward.D2',
  matrices_recentered = FALSE, models_refitted = FALSE), file.path(root, 'settings.qs'))
writeLines(c('TRAINING PROGRAM DISPLACEMENT HEATMAPS', '',
  'training_program_overviews.pdf: two pages per bank, raw and training-SD-scaled displacements.',
  'ERBB_all_fixed_drug_dose_heatmaps.pdf: all 34 contexts from the model with training r=0.165.',
  'MEK_all_fixed_drug_dose_heatmaps.pdf: all 29 contexts from the model with training r=0.283.',
  'Each fixed-context page shows the full program response and its across-cell deviation side by side.',
  'The measured/predicted bars above each matrix follow exactly the same column order.',
  'Column identities are printed on fixed-context pages. Drug, dose and cell colors identify overview columns.',
  'Blue displacement is negative, red is positive. No outcome-based orientation or clustering is used.',
  'Bar colors identify measured sensitivity: orange <=50% vehicle, blue >50%.',
  'Support annotation identifies the 40-60% buffer excluded from direct fitting contrasts.',
  'Each bank has distinct memberships; identical program IDs across banks do not identify the same genes.',
  'Input qs objects, long TSV scores, annotations, scaling and exact orders are exported in data/.',
  'FIGURE_INDEX.tsv lists every individual PDF, available PNG preview and caption.', '',
  'WHAT THIS CAN SHOW',
  'The overview shows structure across treatments and cells. That alone does not establish recipient discrimination.',
  'Fixed drug/dose pages isolate recipient structure, without changing programs or choosing successful examples.',
  'All figures use measured expression displacements; no functional loadings or outcomes are substituted for scores.'),
  file.path(root, 'READING_GUIDE.txt'))
stopifnot(all(file.exists(figure_index$pdf)), all(file.info(figure_index$pdf)$size > 0),
  nrow(figure_index) == 67L)
status('COMPLETE')
