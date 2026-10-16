#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/113_nonlinear_molecular_bypass_interaction/linear/data'
mutation <- fread('/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv')
print(names(mutation)[grepl('BRAF|DepMap|Model|depmap|^V1$',names(mutation),ignore.case=TRUE)])
setnames(mutation,1,'depmap_id')
setnames(mutation,'BRAF (673)','BRAF')
stopifnot('BRAF' %in% names(mutation))
rna <- fread(file.path(root,'external_exact_drug_molecular_observations.tsv'))
P <- sprintf('P%02d',1:40)
rna <- merge(rna,mutation[,.(depmap_id,BRAF)],by='depmap_id')
print(rna[,.(depmap_id,drug,BRAF,P21,P35,P39)])
model <- qread('/work/results/program_pilot/113_nonlinear_molecular_bypass_interaction/linear/models/final.qs')
drugs <- rna[,lapply(.SD,mean),by=.(drug,BRAF),.SDcols=P]
scores <- list()
for(i in seq_len(nrow(drugs))) {
  map <- model$external_rna[[drugs$drug[i]]]
  D <- as.matrix(drugs[i,..P])
  score <- (D/map$sigma) %*% model$context$W
  scores[[i]] <- data.table(drug=drugs$drug[i],BRAF=drugs$BRAF[i],sig1=score[1,1],sig2=score[1,2])
}
out <- rbindlist(scores)
fwrite(out,file.path(root,'BRAF_mutant_molecular_signature_diagnostic.tsv'),sep='\t')
print(out)
