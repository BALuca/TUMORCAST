#!/usr/bin/env Rscript
# Redraw existing diagnostic figures from frozen tables. No model fitting.
entry <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
source(file.path(dirname(normalizePath(entry)), 'config.R'))
check_packages(c('data.table', 'ggplot2', 'ComplexHeatmap', 'circlize', 'dplyr',
                 'matrixStats', 'RColorBrewer', 'viridisLite'))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setDTthreads(1L)
source(file.path(bundle_root, 'code/lib_heatmap.R'))
source(file.path(bundle_root, 'code/single_drug_prediction_functions.R'))
root <- file.path(regenerated_root, 'single_drug_bridge/linear')
dir.create(file.path(root, 'figures'), recursive = TRUE, showWarnings = FALSE)
measured <- fread(file.path(result_root, 'single_drug_bridge/linear/data/measured_predictions_frozen.tsv'))
replay <- fread(file.path(result_root, 'single_drug_bridge/linear/data/baseline_predictions_frozen.tsv'))
combo <- fread(file.path(result_root, 'single_drug_bridge/linear/data/combination_predictions_frozen.tsv'))
bank_selection <- fread(file.path(input_root, 'training_outcomes/selected_doses.tsv'))
drugs <- bank_selection$drug
selection <- fread(file.path(result_root, 'single_drug_bridge/linear/data/selected_doses.tsv'))
included_indices <- match(selection$drug, drugs)
short <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
           'Lapatinib', 'Topotecan', 'Vandetanib', 'Vismodegib', 'Crizotinib')

# Exact plotting expressions retained from single-drug bridge, including its original dimensions.
source(file.path(bundle_root, 'code/single_drug_figures_from_tables.R'))
model_paths <- list.files(file.path(result_root, 'single_drug_bridge/linear/models'), pattern = '\\.rds$', full.names = TRUE)
model_paths <- c(model_paths[basename(model_paths) != 'final.rds'], model_paths[basename(model_paths) == 'final.rds'])
book <- file.path(root, 'figures/training_correlations_and_collapsed_signatures_LOO.pdf')
cairo_pdf(book, width = 9, height = 9, onefile = TRUE)
for (path in model_paths) {
  model <- readRDS(path)
  label <- sub('\\.rds$', '', basename(path))
  render_heatmap(model$template, paste(label, paste(model$held_cells, collapse = ', ')))
}
dev.off()
for (label in c('loo_01', 'loo_30', 'final')) {
  model <- readRDS(file.path(result_root, 'single_drug_bridge/linear/models', paste0(label, '.rds')))
  png(file.path(root, 'figures', paste0(label, '_correlations_and_signatures.png')),
      width = 9, height = 9, units = 'in', res = 150, type = 'cairo')
  render_heatmap(model$template, label)
  dev.off()
}

root <- file.path(regenerated_root, 'ligand_rescue')
dir.create(file.path(root, 'figures'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, 'data'), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
metrics <- fread(file.path(result_root, 'ligand_rescue/data/combination_validation_correlations.tsv'))
pred <- fread(file.path(result_root, 'ligand_rescue/data/combination_predictions_frozen.tsv'))
answer <- fread(file.path(result_root, 'ligand_rescue/data/AU565_HGF_predictions_frozen.tsv'))
# Preserves the historical Wilson 12-day comparison; its mismatch is documented.
source(file.path(bundle_root, 'code/ligand_rescue_figures_from_tables.R'))
writeLines(capture.output(sessionInfo()), file.path(regenerated_root, 'figure_session_info.txt'))
cat('Existing single-drug bridge/ligand rescue figures regenerated in regenerated/. No new panel arrangements or refits.\n')
