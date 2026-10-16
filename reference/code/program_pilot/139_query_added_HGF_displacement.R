#!/usr/bin/env Rscript
# Counterfactual only: f(D_lapatinib + D_HGF), using fixed lapatinib calibration.
# No new model, g term, ligand-specific normalization or outcome fitting.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/139_raw_ligand_rescue_transfer/linear'
source_root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
model <- qread(file.path(source_root,'models/final.qs'))
programs <- sprintf('P%02d',1:40)
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
drugs <- bank_selection$drug
contexts <- fread(file.path(source_root,'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
required <- c('context_key','context_index','predict_rna','signature_node','predict_single')
found <- character()
for(e in parse('/work/src/program_pilot/138_pooled_direction_constrained_rescue.R')) {
 if(is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% required) {
   eval(e); found <- c(found,as.character(e[[2]]))
 }
}
stopifnot(setequal(found,required))
H <- fread(file.path(source_root,'data/time_matched_HGF_program_inputs.tsv'))[CellLine=='AU565']
rows <- list(); contributions <- list()
for(dose in c(.5,5)) {
 meta <- data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=dose)
 lap <- predict_single(model,meta)
 repeated <- meta[rep(1,nrow(H))]
 raw_H <- as.matrix(H[,..programs])
 raw_L <- lap$D[rep(1,nrow(H)),,drop=FALSE]
 combined <- predict_single(model,repeated,measured_D=raw_L+raw_H)
 null <- predict_single(model,meta,measured_D=lap$D+0)
 stopifnot(abs(null$sensitivity-lap$sensitivity)<1e-12)
 g <- context_index(meta$drug,meta$dose)
 effective_beta <- lap$amplitude*model$context$W[,lap$route]/model$context$sigma[g,]
 by_program <- -sweep(raw_H,2,effective_beta,'*')
 increment <- lap$sensitivity-combined$sensitivity
 stopifnot(max(abs(rowSums(by_program)-increment))<1e-10)
 rows[[as.character(dose)]] <- data.table(cell='AU565',dose=dose,source_profile=H$source_profile,
   lap_pct_control=100*2^(-lap$sensitivity),summed_displacements_pct_control=100*2^(-combined$sensitivity),
   HGF_increment_log2_abundance=increment,route=lap$route,
   interpretation='f(D_lapatinib + D_HGF) with unchanged lapatinib context calibration; no learned interaction g')
 contributions[[as.character(dose)]] <- data.table(dose=dose,source_profile=H$source_profile,by_program)
}
answer <- rbindlist(rows)
fwrite(answer,file.path(root,'data/added_HGF_displacement_query.tsv'),sep='\t')
fwrite(rbindlist(contributions),file.path(root,'data/added_HGF_displacement_program_contributions.tsv'),sep='\t')
print(answer[source_profile=='Primary_study_mean'])
print(answer[dose==5 & !grepl('Primary|Leave_out',source_profile),.(source_profile,summed_displacements_pct_control)])
