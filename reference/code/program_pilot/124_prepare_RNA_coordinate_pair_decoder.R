#!/usr/bin/env Rscript
# One RNA-only coordinate outside the existing40 means; no phenotype-based gene choice.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(rhdf5);library(AnnotationDbi);
  library(org.Hs.eg.db);library(limma);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/124_HER2_RNA_coordinate_pair_decoder/linear'
for(folder in c('data','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
definition <- qread(file.path(root,'data/HGF_consensus_coordinate_definition.qs'))
# The definition was made solely from independent HGF RNA studies. Copy its
# explicit fixed rule and source weights, rather than any outcome-fitted score.
weights <- definition$weights
genes <- definition$genes
programs <- sprintf('P%02d',1:40)
archive <- readRDS("/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds")
W <- archive$prep$W
# All ligand coordinates use full gene movements, the same transport and fixed
# loading as their40 scores. Existing HGF source views are kept individually.
meta <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/cytosig_full_metadata.tsv')
meta[,treatment:=toupper(trimws(treatment))]
ext <- fread('/work/results/program_pilot/123_HER2_ligand_combination_transfer/linear/data/expanded_external_ligand_programs.tsv')
mema <- fread('/work/results/program_pilot/118_time_matched_HGF_bypass_transfer/linear/data/time_matched_mema_molecular_inputs.tsv')
targets <- unique(rbind(ext[,.(CellLine,depmap_id,gene)],mema[,.(CellLine,depmap_id,gene)]))
needed <- sort(unique(targets$gene))
meta[treatment=='INSULIN',treatment:='INS']
meta[treatment=='SCF',treatment:='KITLG']
meta[treatment=='FL',treatment:='FLT3LG']
meta <- meta[treatment %in% needed & !grepl('AU[ _-]?565|HCC[ _-]?1954',paste(profile_id,condition),ignore.case=TRUE)]
meta[is.na(condition)|condition=='',condition:='Unspecified context']
meta[,context_key:=paste(treatment,study,condition,sep='::')]
cy <- fread('/work/data/01_CytoSig_data/diff.merge.gz',select=c('V1',meta$profile_id),showProgress=FALSE)
setnames(cy,1,'gene')
FC <- vapply(needed,function(lg) {
  contexts <- unique(meta[treatment==lg,context_key])
  CX <- vapply(contexts,function(key) {
    columns <- meta[context_key==key,profile_id]
    v <- rowMeans(as.matrix(cy[,..columns]),na.rm=TRUE)
    v[!is.finite(v)] <- 0
    v
  },numeric(nrow(cy)))
  studies <- meta[match(contexts,context_key),study]
  SX <- vapply(unique(studies),function(st)rowMeans(CX[,studies==st,drop=FALSE]),numeric(nrow(cy)))
  rowMeans(SX)
},numeric(nrow(cy)))
rownames(FC) <- cy$gene
FC[,'HGF'] <- rowMeans(definition$SX[rownames(FC),,drop=FALSE])
study_names <- colnames(definition$SX)
leave_study_profiles <- vapply(study_names,function(st)
  rowMeans(definition$SX[rownames(FC),study_names!=st,drop=FALSE]),numeric(nrow(FC)))
colnames(leave_study_profiles) <- paste0('Leave_out_',study_names)
all_profiles <- cbind(FC,definition$SX[rownames(FC),,drop=FALSE],leave_study_profiles)
rownames(all_profiles) <- cy$gene
rm(cy)
ccle <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
connection <- file(ccle,'r')
header <- readLines(connection,1)
rows <- list()
repeat {
  lines <- readLines(connection,100,warn=FALSE)
  if(!length(lines))break
  ids <- gsub('"','',sub(',.*','',lines))
  for(i in which(ids %in% targets$depmap_id))rows[[ids[i]]] <- lines[i]
  if(length(rows)==uniqueN(targets$depmap_id))break
}
close(connection)
outputs <- list()
for(cell in unique(targets$CellLine)) {
  id <- targets[CellLine==cell,depmap_id][1]
  b <- fread(text=paste(header,rows[[id]],sep='\n'))
  names_gene <- sub(' \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$','',names(b)[-1],perl=TRUE)
  raw <- data.table(gene=names_gene,u=pmax(2^as.numeric(b[,-1])-1,0))[,.(u=mean(u)),by=gene]
  u <- raw$u
  reference <- u/sum(u)*1e4
  index <- match(raw$gene,rownames(all_profiles))
  gi <- match(genes,raw$gene)
  stopifnot(!anyNA(gi))
  convert <- function(fc) {
    fc <- fc[index]
    fc[!is.finite(fc)] <- 0
    treated <- u*2^fc
    treated <- treated/sum(treated)*1e4
    ok <- treated>=.1 & reference>=.1
    mv <- log2(treated[ok]/reference[ok])
    av <- (log2(treated[ok])+log2(reference[ok]))/2
    ml <- quantile(mv,c(.3,.7));al <- quantile(av,c(.05,.95))
    keep <- mv>=ml[1]&mv<=ml[2]&av>=al[1]&av<=al[2]
    factor <- 2^median(mv[keep])
    sum((log2(1+treated[gi]/factor)-log2(1+reference[gi]))*weights)
  }
  scores <- vapply(seq_len(ncol(all_profiles)),function(k)convert(all_profiles[,k]),numeric(1))
  outputs[[cell]] <- data.table(CellLine=cell,depmap_id=id,input=colnames(all_profiles),HGF_RNA_coordinate=scores)
}
write_table(rbindlist(outputs),'ligand_HGF_coordinate_molecular_inputs')
write_table(data.table(genes=length(genes),original_responsive_genes=sum(fread(file.path(root,'data/HGF_consensus_coordinate_weights.tsv'))$consistent_loading!=0),
  fraction_captured_by_original40=definition$fraction_captured,orthogonality=max(abs(crossprod(W[genes,,drop=FALSE],weights)))),
 'HGF_coordinate_definition_audit')
capture.output(sessionInfo(),file=file.path(root,'coordinate_preparation_session.txt'))
cat("Expanded RNA coordinate preparation complete.\n")
