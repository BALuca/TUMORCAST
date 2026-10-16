#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-p-linux-gnu','/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/112_external_signature_anchored_interaction/linear/data'
model <- qread('/work/results/program_pilot/112_external_signature_anchored_interaction/linear/models/final.qs')
P <- sprintf('P%02d',1:40)
ctx <- model$context
pred <- fread(file.path(root,'external_predicted_drug_programs.tsv'))
rows <- list()
for(dr in unique(pred$drug)) {
  map <- model$external_rna[[dr]]
  D <- as.matrix(pred[drug==dr,..P])
  scores <- (sweep(D,2,map$sigma,'/') %*% ctx$W)
  rows[[dr]] <- data.table(drug=dr,sig1=mean(scores[,1]),sig2=mean(scores[,2]))
}
# Parent lapatinib receptor forecast at native5, no functions rerun or Y inputs.
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
receptors <- as.data.table(a$prep$receptors)
lap <- model$rna$lap_receptors
X <- as.matrix(receptors[depmap_id=='ACH-000248',.(ERBB2,MET)])
D <- sweep(sweep(sweep(sweep(X,2,lap$xmean),2,lap$xsd,'/') %*% lap$B,2,lap$sigma,'*'),2,lap$mu,'+')
j <- 6L
stopifnot(fread(file.path(root,'molecular_contexts.tsv'))$drug[j]=='Lapatinib ditosylate')
scores <- (D/ctx$sigma[j,,drop=FALSE]) %*% ctx$W
rows[['lap']] <- data.table(drug='Lapatinib5_AU565',sig1=scores[1,1],sig2=scores[1,2])
out <- rbindlist(rows)
fwrite(out,file.path(root,'drug_signature_sign_transfer.tsv'),sep='\t')
print(out)
