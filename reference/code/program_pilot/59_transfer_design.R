#!/usr/bin/env Rscript
# Feasibility count for the dataset59 transfer design. No model or label-based
# selection; matches are drug identity and cell identity, not dose/time matches.
suppressPackageStartupMessages(library(data.table))
root <- '/work/results/program_pilot/59_mechanistic_evidence_datasets'
key <- function(x) tolower(gsub('[^a-zA-Z0-9]', '', x))
cell_key <- function(x) toupper(key(x))
lincs_path <- '/work/v20_pilot_design_2026-10-09/source/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz'
lincs <- fread(lincs_path)[pert_type == 'trt_cp']
lincs[, `:=`(cell = cell_key(cell_id), drug = key(pert_iname))]
lincs[drug == 'gdc0941', drug := 'pictilisib']
jaaks_paths <- file.path('/work/jaaks_2022/raw', paste0(c('Original','Validation'),
  ' screen_All tissues_fitted data.csv'))
jaaks <- rbindlist(lapply(jaaks_paths, function(path) {
  x <- fread(path, select = c('CELL_LINE_NAME','ANCHOR_NAME','LIBRARY_NAME','ANCHOR_TARGET','LIBRARY_TARGET'))
  x[, source_file := path]
  x
}))
jaaks[, `:=`(cell = cell_key(CELL_LINE_NAME), A = key(ANCHOR_NAME), B = key(LIBRARY_NAME))]
jaaks <- jaaks[cell %chin% c('MCF7','BT20','HS578T','MDAMB231')]
available <- unique(paste(lincs$cell, lincs$drug))
jaaks[, both_molecular := paste(cell,A) %chin% available & paste(cell,B) %chin% available]
jaaks[, pair := fifelse(A < B, paste(A,B,sep=' + '), paste(B,A,sep=' + '))]
selected <- unique(jaaks[both_molecular == TRUE,
  .(cell,pair,ANCHOR_NAME,LIBRARY_NAME,ANCHOR_TARGET,LIBRARY_TARGET)])
selected[, primary_held_pair := pair == 'pictilisib + trametinib']
summary <- selected[, .(unique_pairs = uniqueN(pair),
  other_pairs_for_training = uniqueN(pair[!primary_held_pair]),
  focal_pair_present = any(primary_held_pair)), by = cell]
summary[, limitation := 'Identity-level availability; dose/time/replicate QC and program measurement validation precede model fitting']
stopifnot(all(summary$focal_pair_present), all(summary$other_pairs_for_training > 0))
fwrite(selected, file.path(root,'data/transfer_design_available_Jaaks_LINCS_pairs.tsv'),sep='\t')
fwrite(summary, file.path(root,'data/transfer_design_training_coverage.tsv'),sep='\t')
writeLines(c('One explicit alias: LINCS GDC-0941 = Jaaks Pictilisib. Remaining names match after punctuation/case normalization.',
  'Count unique unordered drug pairs within each cell; reciprocal anchor/library layouts are not new pairs.',
  'Focal held pair: trametinib+pictilisib, every cell and dose, including both screen releases and orientations.',
  'Counts are not independent sample sizes and are not an assertion of dose-matched RNA.',
  'LINCS landmark restriction must be validated against full RNA before calling a measured subset a bank43 program.'),
  file.path(root,'TRANSFER_DESIGN_SCOPE.txt'))
print(summary)
