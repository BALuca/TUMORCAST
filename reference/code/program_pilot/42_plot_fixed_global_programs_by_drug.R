#!/usr/bin/env Rscript
# Render the saved script-42 bank per drug without repeating any fitting stages.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2)})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(8L)
output_dir <- '/work/results/program_pilot/42_global_pooled_signature_programs_ward'
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
model_dir <- file.path(output_dir, 'models')
bank_path <- file.path(model_dir, 'program_bank.qs')
bank <- qread(bank_path)
check(bank$diagnostics$selected_rank == 60L && ncol(bank$weights) == 60L,
  'This plotting run requires the requested saved 60-program bank')
signature <- digest::digest(list(bank$genes, bank$partition, bank$weights, bank$score_sd), algo = 'sha256')
prep_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(prep_path)
cell_map <- panel$cell_map
metadata <- copy(bank$metadata)
sensitivity <- copy(panel$sensitivity)
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm = TRUE) / 25) * 25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm = TRUE) / 25) * 25)
# Load only the admitted genes. No screen, regression or hierarchy is rerun.
responses <- list()
for (id in cell_map$depmap_id) {
  molecular <- qread(panel$input_paths[match(id, cell_map$depmap_id)])
  ids <- metadata[depmap_id == id, observation_id]
  responses[[id]] <- molecular$response[ids, bank$genes, drop = FALSE]
}
prepared <- list(response = do.call(rbind, responses))
prepared$response <- prepared$response[metadata$observation_id, , drop = FALSE]
check(identical(rownames(prepared$response), metadata$observation_id) && all(is.finite(prepared$response)),
  'Plot source rows do not match the saved global bank')
rm(responses, molecular)
invisible(gc())
bank <- restore_functional_bank(bank, prepared$response)
coverage <- fread(file.path(data_dir, 'drug_assay_coverage.tsv'))
plot_drugs <- coverage[eligible == TRUE, drug]
plot_drugs <- c(intersect(local_drugs, plot_drugs), sort(setdiff(plot_drugs, local_drugs)))
message('Rendering ', length(plot_drugs), ' drugs with the same saved 60-program bank')
view_status <- list()
for (view_drug in plot_drugs) {
  idx <- which(bank$metadata$drug == view_drug)
  view <- bank
  view$scope <- paste0('drug_', gsub('[^A-Za-z0-9]+', '_', view_drug))
  view$orientation <- 'drug'
  view$metadata <- bank$metadata[idx]
  view$scores <- if (is.null(bank$scores)) NULL else bank$scores[idx, , drop = FALSE]
  cols <- which(bank$functional_support$observations$observation_id %chin% view$metadata$observation_id)
  view$functional_support$observations <- bank$functional_support$observations[cols]
  view$functional_support$loadings <- bank$functional_support$loadings[, cols, drop = FALSE]
  view$functional_support$raw_prediction_gain <- bank$functional_support$raw_prediction_gain[, cols, drop = FALSE]
  view$functional_support$folds <- bank$functional_support$folds[observation_id %chin% view$metadata$observation_id]
  # This small reference suppresses repeated export of the global quadratic
  # correlation input while preserving globally fixed program colors.
  if (!is.null(bank$weights)) view$global_structure <- list(path = bank_path,
    program_colors = setNames(hcl.colors(ncol(bank$weights), 'Dark 3'), colnames(bank$weights)))
  view_data <- file.path(output_dir, view$scope, 'data')
  dir.create(view_data, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, view$scope, 'models'), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, view$scope, 'figures'), recursive = TRUE, showWarnings = FALSE)
  if (is.null(view$weights)) render_unresolved(view) else {
    render_program_score_clustered_maps(view, prepared$response)
    render_local(view)
  }
  render_functional_loadings(view)
  check(identical(view$genes, bank$genes) && identical(view$partition, bank$partition) &&
    identical(view$weights, bank$weights), 'Per-drug view changed the global programs')
  writeLines(c(paste('Fixed global bank:', bank_path),
    'No gene screen, regression, rank selection or clustering of genes was repeated for this drug.'),
    file.path(output_dir, view$scope, 'README.txt'))
  view_status[[view_drug]] <- data.table(drug = view_drug, status = 'plots complete',
    genes = length(bank$genes), rank = bank$diagnostics$selected_rank, signatures = length(idx))
  fwrite(rbindlist(view_status), file.path(output_dir, 'per_drug_plot_status.tsv'), sep = '\t')
  message('FIXED GLOBAL PROGRAM FIGURES READY: ', file.path(output_dir, view$scope, 'figures'))
  rm(view)
} 
check(identical(signature, digest::digest(list(bank$genes, bank$partition, bank$weights, bank$score_sd), algo = 'sha256')), 'Plotting changed the global bank')
fwrite(data.table(stage = 'per-drug plots complete', rank = 60L, drugs = length(plot_drugs)), file.path(output_dir, 'run_status.tsv'), sep = '\t')
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'per_drug_plot_sessionInfo.txt'))
message('COMPLETE: all per-drug plots use the same global 60-program bank')
