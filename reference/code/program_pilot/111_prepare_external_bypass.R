#!/usr/bin/env Rscript
# Prepare exact-compound RNA and independent ligand/RAF/MEK departures.
# No lapatinib+HGF or Wilson outcome is used. No substitute kinase inhibitor.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(readxl);library(rhdf5);
  library(AnnotationDbi);library(org.Hs.eg.db);library(limma);library(glmnet);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/111_external_bypass_combination_transfer/linear'
for(folder in c('data','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
source <- '/work/v20_pilot_design_2026-10-09/source'
metadata_path <- file.path(source,'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz')
expression_path <- file.path(source,'GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx')
url <- 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11183/MediaObjects/41586_2012_BFnature11183_MOESM12_ESM.zip'
zip_path <- file.path(root,'data/straussman_2012_original_tables.zip')
download.file(url,zip_path,quiet=TRUE,mode='wb')
unzip(zip_path,exdir=file.path(root,'data/source'))
workbook <- file.path(root,'data/source/nature11183-s2/2011-11-14225B-Supplementary table 8 - Rescue effect of 22 cytokines.xlsx')
x <- as.data.table(read_excel(workbook,col_names=FALSE,col_types='text',.name_repair='minimal'))
setnames(x,paste0('V',seq_len(ncol(x))))
cells <- c('A2058','C32','COLO829','G361','SKMEL28','SKMEL5')
annotation <- fread('/work/outputs/dataset_inventory/cache/depmap/Model.csv')
ids <- annotation[match(cells,StrippedCellLineName),ModelID]
stopifnot(length(ids)==6L,!anyNA(ids),!anyDuplicated(ids))
ligand_names <- x$V1
for(i in seq_len(nrow(x)))if(i>1L && is.na(ligand_names[i]))ligand_names[i] <- ligand_names[i-1L]
# Exact single-protein mappings. NRG1 alpha is excluded because the cached NRG1
# transcriptome averages isoforms. Collagen II has no soluble-ligand match.
mapping <- data.table(name=c('Angiopoietin-1','BDNF','EGF','FGF1','flt-3 ligand (FL)','Gas6','GDNF','HGF',
  'IGF-I','Insulin','MSP','NGF','NT3','PDGF-BB','Pleiotrophin','SCF','VEGF-A','VEGF-C','Wnt1'),
  gene=c('ANGPT1','BDNF','EGF','FGF1','FLT3LG','GAS6','GDNF','HGF','IGF1','INS','MST1','NGF','NTF3',
    'PDGFB','PTN','KITLG','VEGFA','VEGFC','WNT1'))
meta <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/cytosig_full_metadata.tsv')
meta[,treatment:=toupper(trimws(treatment))]
meta[treatment=='INSULIN',treatment:='INS']
meta[treatment=='SCF',treatment:='KITLG']
meta[treatment=='FL',treatment:='FLT3LG']
mapping[,molecular_support:=gene %in% meta$treatment]
write_table(mapping,'external_ligand_mapping')
supported <- mapping[molecular_support==TRUE]
outcomes <- list()
for(i in which(is.finite(suppressWarnings(as.numeric(x$V2))))) {
  gene <- supported$gene[match(ligand_names[i],supported$name)]
  if(is.na(gene))next
  dose <- as.numeric(x$V2[i])
  # Fixed 25 ng/mL query: it is present in each supported protein's tested grid.
  # No functional response is used to select a ligand dose.
  if(abs(dose-25)>1e-10)next
  for(j in seq_along(cells)) {
    columns <- 3L+3L*(j-1L)+0:2
    fold <- as.numeric(unlist(x[i,columns,with=FALSE]))
    stopifnot(all(is.finite(fold)),all(fold>0))
    for(k in 1:2)outcomes[[paste(i,j,k)]] <- data.table(CellLine=cells[j],depmap_id=ids[j],gene=gene,
      drug=c('PD-184352','PLX-4720')[k],dose=c(1,2)[k],ligand_dose_ng_mL=dose,
      observed_ligand_log2=log2(fold[1]),observed_drug_arm_increment_log2=log2(fold[k+1]),
      observed_departure_log2=log2(fold[k+1])-log2(fold[1]))
  }
}
outcomes <- rbindlist(outcomes)
stopifnot(nrow(outcomes)>60L,'HGF' %in% outcomes$gene,!anyNA(outcomes))
write_table(outcomes,'external_combination_outcomes')
cat('External outcomes:',nrow(outcomes),'profiles;',uniqueN(outcomes$gene),'ligands; six melanoma cells.\n')

# RNA inputs use all independent CytoSig studies equally; never weighted by rescue.
needed <- sort(unique(outcomes$gene))
meta <- meta[treatment %in% needed & !grepl('AU[ _-]?565|HCC[ _-]?1954',paste(profile_id,condition),ignore.case=TRUE)]
meta[is.na(condition)|condition=='',condition:='Unspecified context']
meta[,context_key:=paste(treatment,study,condition,sep='::')]
cy <- fread('/work/data/01_CytoSig_data/diff.merge.gz',select=c('V1',meta$profile_id),showProgress=FALSE)
setnames(cy,1,'gene')
FC <- vapply(needed,function(lg) {
  contexts <- unique(meta[treatment==lg,context_key])
  cm <- vapply(contexts,function(key) {
    columns <- meta[context_key==key,profile_id]
    v <- rowMeans(as.matrix(cy[,..columns]),na.rm=TRUE)
    v[!is.finite(v)] <- 0
    v
  },numeric(nrow(cy)))
  studies <- meta[match(contexts,context_key),study]
  sm <- vapply(unique(studies),function(st)rowMeans(cm[,studies==st,drop=FALSE]),numeric(nrow(cy)))
  rowMeans(sm)
},numeric(nrow(cy)))
rownames(FC) <- cy$gene
rm(cy)
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
W <- archive$prep$W
programs <- colnames(W)
baseline <- as.data.table(archive$prep$raw)
original <- qread('/work/results/program_pilot/107_within_pair_departure_interaction/linear/models/final.qs')
stopifnot(all(ids %in% baseline$depmap_id))
expression_ccle <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
connection <- file(expression_ccle,'r')
header <- readLines(connection,1)
rows <- list()
repeat {
  lines <- readLines(connection,100,warn=FALSE)
  if(!length(lines))break
  rowids <- gsub('"','',sub(',.*','',lines))
  for(j in which(rowids %in% ids))rows[[rowids[j]]] <- lines[j]
  if(length(rows)==length(ids))break
}
close(connection)
stopifnot(length(rows)==6L)
ligand_inputs <- list()
for(i in seq_along(cells)) {
  b <- fread(text=paste(header,rows[[ids[i]]],sep='\n'))
  genes <- sub(' \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$','',names(b)[-1],perl=TRUE)
  raw <- data.table(gene=genes,u=pmax(2^as.numeric(b[,-1])-1,0))[,.(u=mean(u)),by=gene]
  u <- raw$u
  reference <- u/sum(u)*1e4
  index <- match(raw$gene,rownames(FC))
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
  D <- t(vapply(needed,function(lg)convert(FC[,lg]),numeric(40)))
  colnames(D) <- programs
  ligand_inputs[[cells[i]]] <- cbind(data.table(CellLine=cells[i],depmap_id=ids[i],gene=needed),as.data.table(D))
}
write_table(rbindlist(ligand_inputs),'external_ligand_programs')
write_table(meta,'external_ligand_source_profiles')
rm(FC)

# Exact PLX-4720 and PD-184352 measurements at 24 h. Read own-plate vehicles and
# unselected wells only for quantile normalization; molecular maps have no Y input.
raw_metadata <- fread(metadata_path)
selected <- raw_metadata[pert_iname %in% c('PLX-4720','PD-184352') & pert_time==24 & pert_type=='trt_cp']
normalized_name <- gsub('[^A-Z0-9]','',toupper(selected$cell_id))
selected[,depmap_id:=annotation$ModelID[match(normalized_name,annotation$StrippedCellLineName)]]
selected <- selected[!is.na(depmap_id) & depmap_id %in% baseline$depmap_id]
metadata <- raw_metadata[det_plate %in% selected$det_plate]
qc <- metadata[,.(vehicles=sum(pert_type=='ctl_vehicle'),wells=.N),by=det_plate]
metadata <- metadata[det_plate %in% qc[vehicles>=5,det_plate]]
selected <- selected[det_plate %in% metadata$det_plate]
metadata[,column_index:=match(inst_id,h5read(expression_path,'/0/META/COL/id'))]
stopifnot(!anyNA(metadata$column_index))
setorder(metadata,column_index)
entrez <- as.character(h5read(expression_path,'/0/META/ROW/id'))
symbols <- unname(mapIds(org.Hs.eg.db,keys=entrez,keytype='ENTREZID',column='SYMBOL',multiVals='first'))
stopifnot(!anyNA(symbols),!anyDuplicated(symbols))
expression <- matrix(NA_real_,978,nrow(metadata),dimnames=list(symbols,metadata$inst_id))
cat('Reading LINCS:',nrow(metadata),'wells across',uniqueN(metadata$det_plate),'plates.\n')
for(start in seq(1L,nrow(metadata),by=512L)) {
  keep <- start:min(start+511L,nrow(metadata))
  ix <- metadata$column_index[keep]
  block <- h5read(expression_path,'/0/DATA/0/matrix',start=c(1L,min(ix)),count=c(978L,max(ix)-min(ix)+1L))
  expression[,keep] <- block[,ix-min(ix)+1L,drop=FALSE]
}
stopifnot(all(is.finite(expression)),all(expression>=0))
expression <- log2(1+expression)
for(plate in unique(metadata$det_plate)) {
  ii <- which(metadata$det_plate==plate)
  expression[,ii] <- limma::normalizeQuantiles(expression[,ii,drop=FALSE])
  vehicle <- ii[metadata$pert_type[ii]=='ctl_vehicle']
  expression[,ii] <- sweep(expression[,ii,drop=FALSE],1,rowMeans(expression[,vehicle,drop=FALSE]),'-')
}
genes <- intersect(rownames(W),symbols)
fraction <- colSums(W[genes,,drop=FALSE]^2)/diag(crossprod(W))
partial <- t(expression[genes,,drop=FALSE]) %*% W[genes,,drop=FALSE]
full <- sweep(partial,2,fraction,'/')
observed <- cbind(selected[,.(inst_id,depmap_id,drug=pert_iname,dose=pert_dose,det_plate)],
  as.data.table(full[match(selected$inst_id,metadata$inst_id),,drop=FALSE]))
grid <- c(.04,.12,.37,1.11,3.33,10)
for(dose in grid)observed[abs(log(observed$dose/dose))<=.035,dose:=dose]
curves <- observed[,lapply(.SD,mean),by=.(depmap_id,drug,dose,det_plate),.SDcols=programs]
curves <- curves[,lapply(.SD,mean),by=.(depmap_id,drug,dose),.SDcols=programs]
native <- curves[,{
  target <- if(drug[1]=='PLX-4720')2 else 1
  stopifnot(min(dose)<=target,max(dose)>=target)
  values <- vapply(programs,function(p)approx(log(dose),get(p),xout=log(target))$y,numeric(1))
  cbind(data.table(dose=target),as.data.table(as.list(values)))
},by=.(depmap_id,drug)]
write_table(native,'external_exact_drug_molecular_observations')
write_table(qc,'external_lincs_plate_coverage')
write_table(data.table(program=programs,landmark_loading_fraction=fraction),'external_lincs_observation_operator')
qsave(list(baseline=baseline,original_rna=original$rna,programs=programs),file.path(root,'data/external_molecular_map_inputs.qs'))
paths <- c(zip_path,workbook,metadata_path,expression_path,
 '/work/src/program_pilot/111_prepare_external_bypass.R')
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'external_source_hashes')
capture.output(sessionInfo(),file=file.path(root,'external_preparation_session.txt'))
print(native[,.(molecular_cells=uniqueN(depmap_id)),by=drug])
print(outcomes[,.(rows=.N,cells=uniqueN(depmap_id)),by=gene])
