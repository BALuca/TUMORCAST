#!/usr/bin/env Rscript
# Audit functional-gain extrapolation on new molecular inputs, without choosing by HGF outcomes.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/112_external_signature_anchored_interaction/linear/data'
model <- qread('/work/results/program_pilot/107_within_pair_departure_interaction/linear/models/final.qs')
ctx <- model$context
programs <- sprintf('P%02d',1:40)
drug <- fread(file.path(root,'external_exact_drug_molecular_observations.tsv'))
ligand <- fread(file.path(root,'external_ligand_programs.tsv'))
drug <- drug[,lapply(.SD,mean),by=drug,.SDcols=programs]
drug[,source:='Added exact drug (LINCS)']
ligand[,`:=`(drug=gene,source='Added ligand (CytoSig)')]
raw <- rbind(drug[,c('drug','source',programs),with=FALSE],
  ligand[,c('drug','source',programs),with=FALSE])
D <- as.matrix(raw[,..programs])
Q <- cbind(1,sweep(sweep(sweep(D,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
raw[,`:=`(log_single_effect_gain=as.numeric(Q %*% model$par$log_amplitude),
  footprint_coordinate_distance=sqrt(rowSums(Q[,-1,drop=FALSE]^2)))]
raw[,single_effect_gain:=exp(log_single_effect_gain)]
fwrite(raw[,.(drug,source,log_single_effect_gain,single_effect_gain,footprint_coordinate_distance)],
  file.path(root,'interaction_input_gain_extrapolation.tsv'),sep='\t')
print(raw[source=='Added exact drug (LINCS)',.(drug,single_effect_gain,footprint_coordinate_distance)])
print(raw[,.(minimum_gain=min(single_effect_gain),median_gain=median(single_effect_gain),
  maximum_gain=max(single_effect_gain)),by=source])
print(raw[drug=='HGF',.(single_effect_gain,footprint_coordinate_distance)])
cat('Original drug footprint gain range:',range(exp(ctx$Q[1:10,] %*% model$par$log_amplitude)),'\n')
