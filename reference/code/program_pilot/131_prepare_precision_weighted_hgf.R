#!/usr/bin/env Rscript
# Reconstruct one noisy HGF source from actual matched GEO control/treated data.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(digest)})
root<-'/work/results/program_pilot/131_precision_weighted_HGF_bridge/linear'
for(folder in c('data','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table<-function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
original<-'/work/results/program_pilot/118_time_matched_HGF_bypass_transfer/linear/data'
auxiliary<-'/work/results/program_pilot/124_HER2_RNA_coordinate_pair_decoder/linear/data'
audit<-'/work/results/program_pilot/130_HGF_source_platform_audit/linear/data'
definition<-qread(file.path(auxiliary,'HGF_consensus_coordinate_definition.qs'))
measurements<-fread(file.path(audit,'reconstructed_platform_gene_responses.tsv'))[
  (study=='GSE26732'&preprocessing=='native_RMA_log2')|
  (study=='GSE26735'&preprocessing=='normexp_quantile_log2')]
precision<-fread(file.path(audit,'platform_replicate_reliability.tsv'))[
  (study=='GSE26732'&preprocessing=='native_RMA_log2')|
  (study=='GSE26735'&preprocessing=='normexp_quantile_log2')]
precision[,inverse_variance:=1/expected_mean_contrast_noise_SD^2]
precision[,platform_weight:=inverse_variance/sum(inverse_variance)]
write_table(precision,'HGF_source_platform_precision_weights')
measurements[,weight:=precision$inverse_variance[match(study,precision$study)]]
new_mda<-measurements[,.(logFC=sum(weight*logFC)/sum(weight),platforms=.N),by=gene]
SX<-definition$SX
# Missing native annotation is explicit, and does not become a measured zero.
SX[,'GSE26736']<-new_mda$logFC[match(rownames(SX),new_mda$gene)]
write_table(data.table(gene=rownames(SX),previous_MDA_logFC=definition$SX[,'GSE26736'],
  corrected_MDA_logFC=SX[,'GSE26736'],native_platforms=new_mda$platforms[match(rownames(SX),new_mda$gene)]),
  'corrected_HGF_source_gene_response')
studies<-colnames(SX)
profiles<-cbind(Primary_study_mean=rowMeans(SX,na.rm=TRUE),SX,
  vapply(studies,function(st)rowMeans(SX[,studies!=st,drop=FALSE],na.rm=TRUE),numeric(nrow(SX))))
colnames(profiles)<-c('Primary_study_mean',studies,paste0('Leave_out_',studies))
profiles[!is.finite(profiles)]<-0
archive<-readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
W<-archive$prep$W;P<-colnames(W)
mema<-fread(file.path(original,'time_matched_mema_molecular_inputs.tsv'))
external<-fread('/work/results/program_pilot/123_HER2_ligand_combination_transfer/linear/data/expanded_external_ligand_programs.tsv')
coordinate<-fread(file.path(auxiliary,'ligand_HGF_coordinate_molecular_inputs.tsv'))
cells<-unique(coordinate[,.(CellLine,depmap_id)])
ccle<-'/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
connection<-file(ccle,'r');header<-readLines(connection,1);rows<-list()
repeat {
  lines<-readLines(connection,100,warn=FALSE);if(!length(lines))break
  ids<-gsub('"','',sub(',.*','',lines))
  for(j in which(ids %in% cells$depmap_id))rows[[ids[j]]]<-lines[j]
  if(length(rows)==nrow(cells))break
}
close(connection)
outputs<-list()
for(i in seq_len(nrow(cells))) {
  cell<-cells$CellLine[i];id<-cells$depmap_id[i]
  b<-fread(text=paste(header,rows[[id]],sep='\n'))
  genes<-sub(' \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$','',names(b)[-1],perl=TRUE)
  raw<-data.table(gene=genes,u=pmax(2^as.numeric(b[,-1])-1,0))[,.(u=mean(u)),by=gene]
  u<-raw$u;reference<-u/sum(u)*1e4
  index<-match(raw$gene,rownames(profiles));gi<-match(rownames(W),raw$gene)
  hi<-match(definition$genes,raw$gene)
  stopifnot(!anyNA(gi),!anyNA(hi))
  Z<-matrix(NA_real_,ncol(profiles),41)
  for(k in seq_len(ncol(profiles))) {
    fc<-profiles[index,k];fc[!is.finite(fc)]<-0
    treated<-u*2^fc;treated<-treated/sum(treated)*1e4
    ok<-treated>=.1&reference>=.1;mv<-log2(treated[ok]/reference[ok]);av<-(log2(treated[ok])+log2(reference[ok]))/2
    ml<-quantile(mv,c(.3,.7));al<-quantile(av,c(.05,.95));keep<-mv>=ml[1]&mv<=ml[2]&av>=al[1]&av<=al[2]
    factor<-2^median(mv[keep]);displacement<-log2(1+treated/factor)-log2(1+reference)
    Z[k,1:40]<-as.numeric(displacement[gi] %*% W)
    Z[k,41]<-sum(displacement[hi]*definition$weights)
  }
  colnames(Z)<-c(P,'HGF_RNA_coordinate')
  output<-cbind(cells[rep(i,nrow(Z))],data.table(gene='HGF',source_profile=colnames(profiles)),as.data.table(Z))
  outputs[[cell]]<-output
  if(cell %in% mema$CellLine)mema[CellLine==cell&gene=='HGF',(P):=output[source_profile=='Primary_study_mean',..P]]
  if(cell %in% external$CellLine)external[CellLine==cell&gene=='HGF',(P):=output[source_profile=='Primary_study_mean',..P]]
  for(label in colnames(profiles)) {
    input<-if(label=='Primary_study_mean')'HGF' else label
    value<-output[source_profile==label,HGF_RNA_coordinate]
    if(any(coordinate$CellLine==cell&coordinate$input==input))
      coordinate[CellLine==cell&coordinate$input==input,HGF_RNA_coordinate:=value]
    else coordinate<-rbind(coordinate,data.table(CellLine=cell,depmap_id=id,input,HGF_RNA_coordinate=value))
  }
}
write_table(rbindlist(outputs),'time_matched_HGF_program_inputs')
write_table(mema,'time_matched_mema_molecular_inputs')
write_table(external,'expanded_external_ligand_programs')
write_table(coordinate,'ligand_HGF_coordinate_molecular_inputs')
for(name in c('HGF_consensus_coordinate_definition.qs','HGF_consensus_coordinate_weights.tsv',
  'tahoe_HGF_coordinate_molecular_observations.tsv','external_exact_drug_HGF_coordinate.tsv'))
  stopifnot(file.copy(file.path(auxiliary,name),file.path(root,'data',name),overwrite=TRUE))
for(name in c('HGF_time_context_eligibility.tsv','time_matched_HGF_source_hashes.tsv'))
  stopifnot(file.copy(file.path(original,name),file.path(root,'data',name),overwrite=TRUE))
qsave(list(original_gene_profiles=definition$SX,corrected_gene_profiles=SX,precision=precision,
  missing_annotation_policy='primary mean over available study observations; source-only absent annotation treated as unmeasured'),
  file.path(root,'data/HGF_input_reconstruction.qs'))
paths<-c('/work/src/program_pilot/131_prepare_precision_weighted_hgf.R',file.path(audit,'reconstructed_platform_gene_responses.tsv'),
  file.path(audit,'platform_replicate_reliability.tsv'),file.path(auxiliary,'HGF_consensus_coordinate_definition.qs'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'HGF_source_reconstruction_hashes')
print(precision)
cat('HGF source reconstruction complete; model and programme loadings unchanged.\n')
