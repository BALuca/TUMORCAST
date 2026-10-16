#!/usr/bin/env Rscript
# Check released Level3/Level4 ligand/control annotations independently of outcomes.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages(library(data.table))
input <- '/work/grant_rebuild_2026-10-04/transfer/mema_raw'
output <- '/work/results/program_pilot/108_mema_calibrated_shared_interaction/linear/data'
scratch <- '/work/tmp/mema_control_check'
dir.create(scratch,recursive=TRUE,showWarnings=FALSE)
tables <- list()
for(id in c('LDS-1467','LDS-1468','LDS-1471','LDS-1472','LDS-1475','LDS-1476','LDS-1479','LDS-1480')) {
  path <- file.path(input,paste0(id,'.tar.gz'))
  members <- untar(path,list=TRUE)
  member <- members[grepl('/Data/.+\\.tsv$',members)]
  stopifnot(length(member)==1L)
  untar(path,files=member,exdir=scratch)
  columns <- names(fread(file.path(scratch,member),nrows=0))
  wanted <- intersect(c('CellLine','Drug','Ligand','Ligand1','Barcode'),columns)
  data <- fread(file.path(scratch,member),select=wanted)
  data[,dataset:=id]
  tables[[id]] <- data[,.(spots=.N),by=c('dataset',setdiff(wanted,'Barcode'))]
  unlink(file.path(scratch,member))
}
summary <- rbindlist(tables,fill=TRUE)
fwrite(summary,file.path(output,'released_ligand_control_annotations.tsv'),sep='\t')
print(summary[is.na(Ligand)|grepl('PBS|DMSO|none|control|blank|FBS',Ligand,ignore.case=TRUE)])
print(summary[,.(ligand_annotations=uniqueN(Ligand),missing_ligand_rows=sum(spots[is.na(Ligand)])),by=dataset])
