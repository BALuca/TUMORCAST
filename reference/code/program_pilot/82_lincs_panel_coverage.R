#!/usr/bin/env Rscript
# Read-only assessment of local LINCS coverage for the existing ten-drug panel.
# LINCS expression is not pooled with Tahoe or used to refit functional labels.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(rhdf5)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
setDTthreads(1L)
root <- '/work/results/program_pilot/82_joint_program_rules_and_pair_interactions'
source_root <- '/work/v20_pilot_design_2026-10-09/source'
metadata_path <- file.path(source_root, 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz')
expression_path <- file.path(source_root, 'GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx')
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
archive <- readRDS(archive_path)
panel <- qread(file.path(root, 'models/settings.qs'))$selected_drugs
aliases <- data.table(drug = c('Lapatinib ditosylate', '5-Azacytidine', 'Bendamustine', 'Bortezomib',
  'Busulfan', 'Cytarabine (hydrochloride)', 'Topotecan (hydrochloride)', 'Vandetanib', 'Vemurafenib', 'Vismodegib'),
  pert_iname = c('lapatinib', 'azacitidine', 'bendamustine', 'bortezomib', 'busulfan', 'cytarabine',
    'topotecan', 'vandetanib', 'vemurafenib', 'vismodegib'))
stopifnot(setequal(aliases$drug, panel))
metadata <- fread(metadata_path)
metadata[, pert_iname := tolower(pert_iname)]
metadata <- merge(metadata, aliases, by = 'pert_iname')
normalize_cell <- function(x) toupper(gsub('[^A-Za-z0-9]', '', x))
cell_map <- as.data.table(archive$meta)
cell_map[, cell_key := normalize_cell(CellLineName)]
duplicates <- cell_map[, .N, by = cell_key][N > 1, cell_key]
stopifnot(!any(normalize_cell(metadata$cell_id) %in% duplicates))
metadata[, cell_key := normalize_cell(cell_id)]
metadata[cell_id == 'LNCAP', cell_key := 'LNCAPCLONEFGC']
metadata <- merge(metadata, cell_map[, .(cell_key, depmap_id = ModelID)], by = 'cell_key', all.x = TRUE)
metadata[, `:=`(in_CCLE = depmap_id %in% archive$prep$raw$depmap_id,
  in_Tahoe = depmap_id %in% archive$prep$m$depmap_id)]
columns <- h5read(expression_path, '/0/META/COL/id')
metadata[, in_matrix := inst_id %in% columns]
stopifnot(all(metadata$in_matrix))
metadata[, context := paste(cell_id, pert_dose, pert_time, sep = '|')]
coverage <- metadata[, .(replicate_wells = .N, cell_dose_time_contexts = uniqueN(context), cells = uniqueN(cell_id),
  CCLE_cells = uniqueN(depmap_id[in_CCLE]),
  additional_CCLE_cells_beyond_Tahoe = uniqueN(depmap_id[in_CCLE & !in_Tahoe]),
  time_hours = paste(sort(unique(pert_time)), collapse = ','),
  dose_min_uM = min(pert_dose), dose_max_uM = max(pert_dose),
  cells_names = paste(sort(unique(cell_id)), collapse = ',')), by = .(drug, pert_iname)]
coverage <- merge(aliases, coverage, by = c('drug', 'pert_iname'), all.x = TRUE)
for (column in c('replicate_wells', 'cell_dose_time_contexts', 'cells', 'CCLE_cells',
  'additional_CCLE_cells_beyond_Tahoe')) set(coverage, which(is.na(coverage[[column]])), column, 0L)
entrez <- h5read(expression_path, '/0/META/ROW/id')
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(entrez), keytype = 'ENTREZID',
  column = 'SYMBOL', multiVals = 'first')
W <- archive$prep$W
measured <- rownames(W) %in% na.omit(symbols)
program_coverage <- data.table(program = colnames(W), bank_genes = nrow(W),
  measured_bank_genes = sum(measured),
  fraction_squared_loading_on_measured_genes = colSums(W[measured, , drop = FALSE]^2) / colSums(W^2))
fwrite(coverage, file.path(root, 'data/lincs_selected_panel_coverage.tsv'), sep = '\t')
fwrite(metadata, file.path(root, 'data/lincs_selected_panel_instances.tsv'), sep = '\t')
fwrite(program_coverage, file.path(root, 'data/lincs_landmark_program_coverage.tsv'), sep = '\t')
capture.output(cat('LINCS assessment: existing GSE70138 Level-2 measured-landmark expression.\n',
  'No LINCS values have entered the v82 fit. Dose/time/cell and platform alignment are required before augmentation.\n',
  'Replicate wells are not independent cells. Restricted landmark loading coverage does not measure program accuracy.\n'),
  print(coverage[, .(drug, replicate_wells, cell_dose_time_contexts, cells, CCLE_cells,
    additional_CCLE_cells_beyond_Tahoe, time_hours)]),
  print(program_coverage[, .(bank_genes = unique(bank_genes), measured_bank_genes = unique(measured_bank_genes),
    smallest_loading_fraction = min(fraction_squared_loading_on_measured_genes),
    median_loading_fraction = median(fraction_squared_loading_on_measured_genes),
    largest_loading_fraction = max(fraction_squared_loading_on_measured_genes))]),
  file = file.path(root, 'LINCS_ASSESSMENT.txt'))
print(coverage[, .(drug, replicate_wells, cells, CCLE_cells, additional_CCLE_cells_beyond_Tahoe)])
print(program_coverage[, .(bank_genes = unique(bank_genes), measured_bank_genes = unique(measured_bank_genes),
  median_loading_fraction = median(fraction_squared_loading_on_measured_genes))])
