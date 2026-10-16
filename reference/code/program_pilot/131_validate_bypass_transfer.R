#!/usr/bin/env Rscript
# Validate zero/symmetry, unseen-compound transfer, and HGF source-study dependence.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(glmnet);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/131_precision_weighted_HGF_bridge/linear'
parent <- '/work/results/program_pilot/107_within_pair_departure_interaction/linear'
programs <- sprintf('P%02d',1:40)
source_weight <- .2
interaction_penalty <- .05
assay_penalty <- .01
external_within_weight <- .8
primary_endpoint <- 'normalized_cell_count'
challenge_genes <- c('HGF','NRG1')
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
pearson <- function(x,y)if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10)cor(x,y) else NA_real_
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
drugs <- bank_selection$drug
contexts <- fread(file.path(root,'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
context_key <- function(drug,dose)paste(drug,sprintf('%.10g',dose),sep='|')
context_index <- function(drug,dose)match(context_key(drug,dose),context_key(contexts$drug,contexts$dose))
ligands <- fread(file.path(root,'data/time_matched_mema_molecular_inputs.tsv'))
external_ligands <- fread(file.path(root,'data/external_ligand_programs.tsv'))
external_molecular <- fread(file.path(root,'data/external_exact_drug_molecular_observations.tsv'))
external <- fread(file.path(root,'data/external_combination_outcomes.tsv'))
mutation <- fread('/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv',
  select=c('V1','BRAF (673)'))
setnames(mutation,c('depmap_id','BRAF_hotspot'))
training_labels <- fread(file.path(root,'data/mema_nonchallenge_training_outcomes.tsv'))
pairs <- fread(file.path(root,'data/combination_training_outcomes.tsv'))
original_final <- qread(file.path(root,'models/final_normalized_cell_count.qs'))
auxiliary_root<-file.path(root,'data')
auxiliary_molecular<-fread(file.path(auxiliary_root,'tahoe_HGF_coordinate_molecular_observations.tsv'))
auxiliary_external<-fread(file.path(auxiliary_root,'external_exact_drug_HGF_coordinate.tsv'))
auxiliary_ligands<-fread(file.path(auxiliary_root,'ligand_HGF_coordinate_molecular_inputs.tsv'))
auxiliary_ligand_scale<-sd(auxiliary_ligands$HGF_RNA_coordinate[!grepl('^Leave_out_',auxiliary_ligands$input)])
source_views <- fread(file.path(root,'data/time_matched_HGF_program_inputs.tsv'))[
  !grepl('^Leave_out_|^Primary_',source_profile)]

# Reuse exact fitted function definitions, without executing any pipeline IO/fit.
required <- c('predict_rna','signature_node','pair_features','predict_single','predict_pair',
  'mema_features','center_geometry','external_features','fit_interaction','select_penalty','fit_selected',
  'predict_external_departure','auxiliary_key','ligand_coordinate','fit_auxiliary_rna','predict_auxiliary_coordinate')
expressions <- parse('/work/src/program_pilot/131_precision_weighted_HGF_bridge.R')
found <- character()
for(expression in expressions) {
  if(is.call(expression) && identical(expression[[1]],as.name('<-')) &&
     is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% required) {
    eval(expression,envir=.GlobalEnv)
    found <- c(found,as.character(expression[[2]]))
  }
}
stopifnot(setequal(required,found))
penalty_audit <- list()

query_hgf <- function(model,Dh,profile,auxiliary_input='HGF') {
  ctx <- model$context
  lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=5))
  Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
  scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
    crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
  minimum <- matrix(apply(ctx$sigma[1:10,],2,min),nrow(Dh),40,byrow=TRUE)
  maximum <- matrix(apply(ctx$sigma[1:10,],2,max),nrow(Dh),40,byrow=TRUE)
  hscale <- pmin(pmax(exp(Qh %*% scale_map),minimum),maximum)
  H40 <- Dh/hscale
  H <- cbind(H40,H40 %*% ctx$W,ligand_coordinate(rep('AU565',nrow(Dh)),auxiliary_input))
  L <- lap$program_activity[rep(1L,nrow(Dh)),,drop=FALSE]
  product <- pair_features(L,H,ctx$W,model$interaction_basis)
  g <- -as.numeric(sweep(product,2,ctx$interaction_scale,'/') %*% model$par$Gamma)
  standalone <- -((H40 %*% ctx$W)*as.numeric(exp(Qh %*% model$par$log_amplitude)))[,1]*model$ligand_single_conversion
  data.table(source_profile=profile,lapatinib_pct=100*2^(-lap$sensitivity),
    HGF_increment_pct=100*2^standalone,interaction_log2=g,
    conditional_drug_pct_HGF_control=100*2^(-lap$sensitivity+g),
    combination_pct_vehicle=100*2^(-lap$sensitivity+standalone+g))
}

# Structural checks on actual molecular inputs and the fitted g.
mf <- mema_features(original_final)
A <- mf$A
B <- mf$B
P <- pair_features(A,B,original_final$context$W,original_final$interaction_basis)
swapped <- pair_features(B,A,original_final$context$W,original_final$interaction_basis)
zero <- pair_features(A,B*0,original_final$context$W,original_final$interaction_basis)
stopifnot(max(abs(P-swapped))<1e-12,max(abs(zero))<1e-12)
old <- qread(file.path(parent,'models/final.qs'))
stopifnot(identical(old$rna,original_final$rna),identical(old$par$level,original_final$par$level),
  identical(old$par$log_amplitude,original_final$par$log_amplitude),identical(old$template,original_final$template))
checks <- data.table(check=c('symmetric_g','zero_displacement_g','fixed_original_single_drug_chain'),
  maximum_difference=c(max(abs(P-swapped)),max(abs(zero)),0))
write_table(checks,'structural_validation')

# HGF source sensitivity: first query a frozen fit, then refit g after excluding
# the source study from HGF RNA in EVERY training receiver. No target Y is fitted.
source_profiles <- fread(file.path(root,'data/time_matched_HGF_program_inputs.tsv'))
au <- source_profiles[CellLine=='AU565']
frozen_source <- query_hgf(original_final,as.matrix(au[,..programs]),au$source_profile,
  ifelse(au$source_profile=='Primary_study_mean','HGF',au$source_profile))
write_table(frozen_source,'HGF_source_study_queries_frozen_model')
original_ligands <- copy(ligands)
original_external_ligands <- copy(external_ligands)
original_source_views <- copy(source_views)
original_auxiliary_ligands <- copy(auxiliary_ligands)
refits <- list()
for(profile in unique(source_profiles$source_profile)[grepl('^Leave_out_',unique(source_profiles$source_profile))]) {
  ligands <- copy(original_ligands)
  external_ligands <- copy(original_external_ligands)
  source_views <- original_source_views[source_profile!=sub('^Leave_out_','',profile)]
  auxiliary_ligands<-copy(original_auxiliary_ligands)
  for(cell in unique(auxiliary_ligands$CellLine)) {
    value<-auxiliary_ligands[CellLine==cell & input==profile,HGF_RNA_coordinate]
    stopifnot(length(value)==1L)
    auxiliary_ligands[CellLine==cell & input=='HGF',HGF_RNA_coordinate:=value]
  }
  replacement <- source_profiles[source_profile==profile]
  for(cell in unique(replacement$CellLine)) {
    values <- replacement[CellLine==cell,..programs]
    if(cell %in% ligands$CellLine)ligands[CellLine==cell & gene=='HGF',(programs):=values]
    else external_ligands[CellLine==cell & gene=='HGF',(programs):=values]
  }
  fit <- fit_selected(original_final,mema_features(original_final),primary_endpoint)
  D <- as.matrix(replacement[CellLine=='AU565',..programs])
  query <- query_hgf(fit$model,D,profile,profile)
  query[,selected_penalty:=fit$model$interaction_training$interaction_penalty]
  refits[[profile]] <- query
  query[,selection_run:=profile]
  write_table(rbindlist(penalty_audit)[,selection_run:=profile],paste0('penalty_audit_',profile))
  qsave(fit$model,file.path(root,'models',paste0('RNA_source_',profile,'.qs')))
  cat('RNA source exclusion completed:',profile,'\n')
}
write_table(rbindlist(refits),'HGF_source_study_exclusion_refits')
ligands <- copy(original_ligands)
external_ligands <- copy(original_external_ligands)
source_views <- copy(original_source_views)
auxiliary_ligands <- copy(original_auxiliary_ligands)

# Whole compound exclusion: RNA of that characterized compound is permitted,
# but ALL its ligand-combination functional outcomes are excluded from g fitting.
all_external <- copy(external)
query_external <- external_features(original_final)
predictions <- list()
for(held_drug in unique(all_external$drug)) {
  external <- all_external[all_external$drug!=held_drug]
  stopifnot(!any(external$drug==held_drug))
  fit <- fit_selected(original_final,mema_features(original_final),primary_endpoint)
  rows <- which(all_external$drug==held_drug)
  products <- pair_features(query_external$A[rows,,drop=FALSE],query_external$B[rows,,drop=FALSE],
    fit$model$context$W,fit$model$interaction_basis)
  predicted <- -as.numeric(sweep(products,2,fit$model$context$interaction_scale,'/') %*% fit$model$par$Gamma)
  predictions[[held_drug]] <- all_external[rows,.(CellLine,depmap_id,gene,drug,dose)]
  predictions[[held_drug]][,`:=`(predicted_departure_log2=predicted,
    selected_penalty=fit$model$interaction_training$interaction_penalty)]
  stopifnot(!held_drug %in% external$drug)
  write_table(rbindlist(penalty_audit)[,selection_run:=held_drug],paste0('penalty_audit_hold_',held_drug))
}
external <- all_external
frozen <- rbindlist(predictions)
write_table(frozen,'external_whole_compound_predictions_frozen')
write_table(data.table(path=file.path(root,'data/external_whole_compound_predictions_frozen.tsv'),
  sha256=digest(file.path(root,'data/external_whole_compound_predictions_frozen.tsv'),file=TRUE,algo='sha256')),
  'external_whole_compound_prediction_hash')
scored <- merge(frozen,all_external[,.(CellLine,depmap_id,gene,drug,dose,observed_departure_log2)],
  by=c('CellLine','depmap_id','gene','drug','dose'),all.x=TRUE)
write_table(scored,'external_whole_compound_predictions_scored')
metrics <- scored[,.(ligands=.N,r_departure=pearson(observed_departure_log2,predicted_departure_log2)),
  by=.(CellLine,drug,dose)]
write_table(metrics,'external_whole_compound_correlations')
write_table(rbindlist(penalty_audit),'source_and_compound_validation_penalty_audit')
print(frozen_source)
print(rbindlist(refits))
print(metrics)
