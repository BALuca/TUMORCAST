#!/usr/bin/env Rscript
# Outcome-free HGF source alignment: cancer-cell RNA nearest24h, independent studies equally.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(digest)})
root <- '/work/results/program_pilot/118_time_matched_HGF_bypass_transfer/linear'
for(folder in c('data','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
meta <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/cytosig_full_metadata.tsv')[toupper(treatment)=='HGF']
meta[,hours:=suppressWarnings(as.numeric(sub('[^0-9.].*$','',duration)))]
# GSE135552 RNA methods specify16h; CytoSig labels48h. The source correction is
# independently recorded in v57 and is applied before selecting a time window.
meta[study=='GSE135552',hours:=16]
meta[,study_family:=fifelse(study %in% c('GSE26732','GSE26735'),'GSE26736',study)]
meta[,cancer_cell:=grepl('cancer|H1975|Detroit|Fadu|Scc9|Pc3m|Mda.Mb435',condition,ignore.case=TRUE)]
meta[,eligible:=cancer_cell & is.finite(hours) & hours>=12 & hours<=48 &
  !grepl('AU[ _-]?565|HCC[ _-]?1954',paste(profile_id,condition),ignore.case=TRUE)]
meta[,context_key:=paste(study_family,condition,sep='|')]
selected <- meta[eligible==TRUE]
selected[,distance_to_24h:=abs(log(hours/24))]
selected <- selected[, .SD[distance_to_24h<=min(distance_to_24h)+1e-10],by=context_key]
meta[,selected:=profile_id %in% selected$profile_id]
write_table(meta,'HGF_time_context_eligibility')
stopifnot(uniqueN(selected$study_family)==4L,uniqueN(selected$context_key)==6L,nrow(selected)>=6L)
cy <- fread('/work/data/01_CytoSig_data/diff.merge.gz',select=c('V1',selected$profile_id),showProgress=FALSE)
setnames(cy,1,'gene')
contexts <- unique(selected$context_key)
CX <- vapply(contexts,function(key) {
  cols <- selected[context_key==key,profile_id]
  v <- rowMeans(as.matrix(cy[,..cols]),na.rm=TRUE)
  v[!is.finite(v)] <- 0
  v
},numeric(nrow(cy)))
studies <- selected[match(contexts,context_key),study_family]
study_names <- sort(unique(studies))
SX <- vapply(study_names,function(st)rowMeans(CX[,studies==st,drop=FALSE]),numeric(nrow(cy)))
profiles <- cbind(Primary_study_mean=rowMeans(SX),SX,
  vapply(study_names,function(st)rowMeans(SX[,study_names!=st,drop=FALSE]),numeric(nrow(cy))))
colnames(profiles) <- c('Primary_study_mean',study_names,paste0('Leave_out_',study_names))
rownames(profiles) <- cy$gene
rm(cy)
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
W <- a$prep$W
P <- colnames(W)
mema <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/recipient_ligand_program_scores.tsv')
ext <- fread('/work/results/program_pilot/111_external_bypass_combination_transfer/linear/data/external_ligand_programs.tsv')
cells <- unique(rbind(mema[,.(CellLine,depmap_id)],ext[,.(CellLine,depmap_id)]))
path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
connection <- file(path,'r')
header <- readLines(connection,1)
rows <- list()
repeat {
  lines <- readLines(connection,100,warn=FALSE)
  if(!length(lines))break
  ids <- gsub('"','',sub(',.*','',lines))
  for(i in which(ids %in% cells$depmap_id))rows[[ids[i]]] <- lines[i]
  if(length(rows)==nrow(cells))break
}
close(connection)
stopifnot(length(rows)==8L)
outputs <- list()
for(i in seq_len(nrow(cells))) {
  b <- fread(text=paste(header,rows[[cells$depmap_id[i]]],sep='\n'))
  genes <- sub(' \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$','',names(b)[-1],perl=TRUE)
  raw <- data.table(gene=genes,u=pmax(2^as.numeric(b[,-1])-1,0))[,.(u=mean(u)),by=gene]
  u <- raw$u
  reference <- u/sum(u)*1e4
  index <- match(raw$gene,rownames(profiles))
  gi <- match(rownames(W),raw$gene)
  stopifnot(!anyNA(gi))
  convert <- function(fc) {
    fc <- fc[index]
    fc[!is.finite(fc)] <- 0
    treated <- u*2^fc
    treated <- treated/sum(treated)*1e4
    ok <- treated>=.1 & reference>=.1
    mv <- log2(treated[ok]/reference[ok])
    av <- (log2(treated[ok])+log2(reference[ok]))/2
    ml <- quantile(mv,c(.3,.7))
    al <- quantile(av,c(.05,.95))
    keep <- mv>=ml[1]&mv<=ml[2]&av>=al[1]&av<=al[2]
    factor <- 2^median(mv[keep])
    as.numeric((log2(1+treated[gi]/factor)-log2(1+reference[gi])) %*% W)
  }
  Z <- t(vapply(seq_len(ncol(profiles)),function(k)convert(profiles[,k]),numeric(40)))
  colnames(Z) <- P
  output <- cbind(cells[rep(i,nrow(Z))],data.table(gene='HGF',source_profile=colnames(profiles)),as.data.table(Z))
  outputs[[cells$CellLine[i]]] <- output
  primary <- output[source_profile=='Primary_study_mean',..P]
  if(cells$CellLine[i] %in% mema$CellLine) {
    ii <- which(mema$CellLine==cells$CellLine[i] & mema$gene=='HGF')
    mema[ii,(P):=primary]
  } else {
    ii <- which(ext$CellLine==cells$CellLine[i] & ext$gene=='HGF')
    ext[ii,(P):=primary]
  }
}
write_table(rbindlist(outputs),'time_matched_HGF_program_inputs')
write_table(mema,'time_matched_mema_molecular_inputs')
write_table(ext,'time_matched_external_ligand_programs')
write_table(data.table(path=c('/work/src/program_pilot/118_prepare_time_matched_hgf.R',
 '/work/data/01_CytoSig_data/diff.merge.gz'),sha256=vapply(c('/work/src/program_pilot/118_prepare_time_matched_hgf.R',
 '/work/data/01_CytoSig_data/diff.merge.gz'),digest,character(1),file=TRUE,algo='sha256')),'time_matched_HGF_source_hashes')
print(selected[,.(profile_id,condition,hours,study_family)])
