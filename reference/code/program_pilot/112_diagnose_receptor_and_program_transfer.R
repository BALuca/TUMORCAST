#!/usr/bin/env Rscript
# Diagnose measured HGF departures against receptor baseline and ligand input.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/112_external_signature_anchored_interaction/linear/data'
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
receptors <- as.data.table(a$prep$receptors)
outcomes <- fread(file.path(root,'external_combination_outcomes.tsv'))
z <- merge(outcomes[gene=='HGF'],receptors[,.(depmap_id,ERBB2,MET)],by='depmap_id')
print(z[,.(CellLine,drug,MET,observed_departure_log2)])
print(z[,.(cells=.N,r_MET_departure=cor(MET,observed_departure_log2)),by=drug])
fwrite(z,file.path(root,'external_HGF_receptor_association.tsv'),sep='\t')
# Correlation to observed molecular ligand footprints in each source study is not
# receptor signaling competence. Compare the transported recipient directions.
lig <- fread(file.path(root,'external_ligand_programs.tsv'))[gene=='HGF']
P <- sprintf('P%02d',1:40)
D <- as.matrix(lig[,..P])
cosine <- tcrossprod(D)/sqrt(outer(rowSums(D^2),rowSums(D^2)))
rownames(cosine) <- colnames(cosine) <- lig$CellLine
fwrite(as.data.table(cosine,keep.rownames='CellLine'),file.path(root,'external_HGF_receiver_cosines.tsv'),sep='\t')
cat('HGF receiver cosine range:',range(cosine[lower.tri(cosine)]),'\n')
