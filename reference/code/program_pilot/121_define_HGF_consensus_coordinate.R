.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
root <- '/work/results/program_pilot/121_HGF_consensus_coordinate_bridge/linear/data'
dir.create(root,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(qs);library(rhdf5);library(AnnotationDbi);library(org.Hs.eg.db)})
meta <- fread('/work/results/program_pilot/118_time_matched_HGF_bypass_transfer/linear/data/HGF_time_context_eligibility.tsv')[selected==TRUE]
cy <- fread('/work/data/01_CytoSig_data/diff.merge.gz',select=c('V1',meta$profile_id),showProgress=FALSE)
setnames(cy,1,'gene')
studies <- sort(unique(meta$study_family))
SX <- vapply(studies,function(st) {
  contexts <- unique(meta[study_family==st,context_key])
  CX <- vapply(contexts,function(key) {
    columns <- meta[context_key==key,profile_id]
    values <- rowMeans(as.matrix(cy[,..columns]),na.rm=TRUE)
    values[!is.finite(values)] <- 0
    values
  },numeric(nrow(cy)))
  rowMeans(CX)
},numeric(nrow(cy)))
rownames(SX) <- cy$gene
model <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
W <- model$prep$W
entrez <- as.character(h5read('/work/v20_pilot_design_2026-10-09/source/GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx','/0/META/ROW/id'))
landmark <- unname(mapIds(org.Hs.eg.db,keys=entrez,keytype='ENTREZID',column='SYMBOL',multiVals='first'))
genes <- intersect(intersect(rownames(W),rownames(SX)),landmark)
X <- SX[genes,,drop=FALSE]
median_response <- apply(X,1,median)
agreement <- rowMeans(sign(X)==sign(median_response))
v <- median_response*(agreement>=.75 & abs(median_response)>=.25)
WL <- W[genes,,drop=FALSE]
projected <- WL %*% solve(crossprod(WL),crossprod(WL,v))
residual <- v-as.numeric(projected)
cat('Genes:',length(genes),'; consistent responsive landmarks:',sum(v!=0),
  '; fraction captured by existing40:',sum(projected^2)/sum(v^2),'\n')
stopifnot(sum(v!=0)>=10L,sum(residual^2)>0)
weights <- residual/sqrt(sum(residual^2))
table <- data.table(gene=genes,median_HGF_logFC=median_response,sign_agreement=agreement,
  consistent_loading=v,residual_loading=weights)
fwrite(table,file.path(root,'HGF_consensus_coordinate_weights.tsv'),sep='\t')
print(table[order(-abs(residual_loading))][1:20])
print(data.table(study=studies,untransported_residual_score=as.numeric(crossprod(weights,X))))
qsave(list(genes=genes,weights=weights,SX=SX,metadata=meta,all_gene_names=cy$gene,
  fraction_captured=sum(projected^2)/sum(v^2),old_weights=W),file.path(root,'HGF_consensus_coordinate_definition.qs'))
