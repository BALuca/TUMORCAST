#!/usr/bin/env Rscript
# Localize information loss in transported ligand and predicted drug representations.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/110_signature_anchored_program_interaction/linear/data'
programs <- sprintf('P%02d',1:40)
ligands <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/recipient_ligand_program_scores.tsv')
labels <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/all_gene_endpoints.tsv')
genes <- sort(unique(ligands$gene))
A <- as.matrix(ligands[CellLine=='AU565'][match(genes,gene),..programs])
B <- as.matrix(ligands[CellLine=='HCC1954'][match(genes,gene),..programs])
geometry <- data.table(gene=genes,
  receiver_vector_cosine=rowSums(A*B)/sqrt(rowSums(A^2)*rowSums(B^2)),
  receiver_difference_fraction=sqrt(rowSums((A-B)^2)/rowSums(A^2)),
  AU565_program_norm=sqrt(rowSums(A^2)),HCC1954_program_norm=sqrt(rowSums(B^2)))
fwrite(geometry,file.path(root,'receiver_ligand_geometry.tsv'),sep='\t')
print(geometry[gene %in% c('HGF','NRG1','EGF','AREG','FGF2')])
# Nonchallenge associations characterize the training signal; they do not select HGF predictions.
train <- merge(ligands,labels[!gene %in% c('HGF','NRG1')],by=c('CellLine','gene'))
associations <- rbindlist(lapply(programs,function(p)train[,.(program=p,
  r_program_response=cor(get(p),y)),by=.(CellLine,endpoint)]))
fwrite(associations,file.path(root,'nonchallenge_program_response_associations.tsv'),sep='\t')
print(associations[,.(programs_r_above_05=sum(abs(r_program_response)>.5),
  max_absolute_r=max(abs(r_program_response))),by=.(CellLine,endpoint)])
models <- list(v108='/work/results/program_pilot/108_mema_calibrated_shared_interaction/linear',
  v109='/work/results/program_pilot/109_joint_conditional_interaction/linear',
  v110='/work/results/program_pilot/110_signature_anchored_program_interaction/linear')
effects <- rbindlist(lapply(names(models),function(version) {
  input <- fread(file.path(models[[version]],'data/mema_predictions_scored.tsv'))
  z <- input[evaluation=='fitted' & !gene %in% c('HGF','NRG1')]
  z[,.(version=version,r_fitted=cor(observed_log2,predicted_log2),
    observed_sd=sd(observed_log2),predicted_sd=sd(predicted_log2)),by=.(CellLine,endpoint)]
}))
fwrite(effects,file.path(root,'interaction_capacity_diagnostic.tsv'),sep='\t')
print(effects)
