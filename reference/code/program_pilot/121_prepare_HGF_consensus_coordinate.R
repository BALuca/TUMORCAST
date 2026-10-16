#!/usr/bin/env Rscript
# One RNA-only coordinate outside the existing40 means; no phenotype-based gene choice.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(rhdf5);library(AnnotationDbi);
  library(org.Hs.eg.db);library(limma);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/121_HGF_consensus_coordinate_bridge/linear'
for(folder in c('data','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
definition <- qread(file.path(root,'data/HGF_consensus_coordinate_definition.qs'))
# The definition was made solely from independent HGF RNA studies. Copy its
# explicit fixed rule and source weights, rather than any outcome-fitted score.
qsave(definition,file.path(root,'data/HGF_consensus_coordinate_definition.qs'))
weights <- definition$weights
genes <- definition$genes
programs <- sprintf('P%02d',1:40)
bank <- '/work/almanac_conditionals/state_conditioned_bridge/program_first'
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
W <- archive$prep$W
molecular <- as.data.table(archive$prep$m)
contexts <- fread('/work/results/program_pilot/107_within_pair_departure_interaction/linear/data/molecular_contexts.tsv')
keys <- function(cell,drug,dose)paste(cell,drug,sprintf('%.10g',dose),sep='|')
selected <- merge(molecular,contexts,by=c('drug','dose'))
cache_dir <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs'
source_key <- function(drug) {
  out <- tolower(trimws(drug))
  out[out %in% c('5-azacitidine','5-azacytidine')] <- 'azacitidine'
  out
}
input_parts <- list()
cache_coverage <- data.table(depmap_id=sort(unique(selected$depmap_id)))
cache_coverage[,gene_cache:=file.path(cache_dir,paste0(depmap_id,'.qs'))]
cache_coverage[,available:=file.exists(gene_cache)]
cache_coverage[,disposition:=ifelse(available,'full gene RNA cache available','no full gene cache; not fabricated')]
write_table(cache_coverage,'full_gene_cache_coverage')
for(id in sort(unique(selected$depmap_id))) {
  path <- file.path(cache_dir,paste0(id,'.qs'))
  # Exact full-gene coverage is explicit; absent caches are not fabricated.
  if(!file.exists(path))next
  x <- qread(path)
  stopifnot(identical(x$normalizer,'paired_trimmed_median_log_ratio_cp10k_v1'),
    all(genes %in% colnames(x$response)),all(rownames(W) %in% colnames(x$response)))
  md <- as.data.table(x$metadata)
  md[,`:=`(drug_key=source_key(source_drug),rna_row=.I)]
  wanted <- contexts[,.(drug,drug_key=source_key(drug),dose)]
  take <- merge(md,wanted,by=c('drug_key','dose'),suffixes=c('_cache',''))
  stopifnot(!anyDuplicated(take[,.(depmap_id,drug,dose)]))
  D <- x$response[take$rna_row,,drop=FALSE]
  replayed <- D[,rownames(W),drop=FALSE] %*% W
  stored_index <- match(keys(take$depmap_id,take$drug,take$dose),keys(selected$depmap_id,selected$drug,selected$dose))
  stopifnot(!anyNA(stored_index))
  stored <- as.matrix(selected[stored_index,..programs])
  take[,`:=`(HGF_RNA_coordinate=as.numeric(D[,genes,drop=FALSE] %*% weights),
    maximum_old_program_replay_difference=apply(abs(replayed-stored),1,max),gene_cache=path)]
  input_parts[[id]] <- take[,.(depmap_id,drug,dose,HGF_RNA_coordinate,maximum_old_program_replay_difference,gene_cache)]
}
source_rows <- rbindlist(input_parts)
stopifnot(nrow(source_rows)>350L,!anyDuplicated(source_rows[,.(depmap_id,drug,dose)]))
write_table(source_rows[,.(depmap_id,drug,dose,HGF_RNA_coordinate,maximum_old_program_replay_difference)],
  'tahoe_HGF_coordinate_molecular_observations')
write_table(source_rows[,.(profiles=.N,max_difference=max(maximum_old_program_replay_difference)),by=drug],
  'tahoe_old_program_replay_audit')
write_table(source_rows[,.(RNA_profiles=.N),by=.(drug,dose)],'auxiliary_coordinate_RNA_context_coverage')

# Exact added compounds: new coordinate uses only measured landmark genes.
# The existing40 LINCS estimates are unchanged; no imputed genes enter this axis.
source_root <- '/work/v20_pilot_design_2026-10-09/source'
path <- file.path(source_root,'GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx')
metadata <- fread(file.path(source_root,'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz'))
annotation <- fread('/work/outputs/dataset_inventory/cache/depmap/Model.csv')
baseline <- as.data.table(archive$prep$raw)
selected <- metadata[pert_iname %in% c('PLX-4720','PD-184352') & pert_time==24 & pert_type=='trt_cp']
selected[,depmap_id:=annotation$ModelID[match(gsub('[^A-Z0-9]','',toupper(cell_id)),annotation$StrippedCellLineName)]]
selected <- selected[!is.na(depmap_id) & depmap_id %in% baseline$depmap_id]
metadata <- metadata[det_plate %in% selected$det_plate]
plates <- metadata[,.(vehicles=sum(pert_type=='ctl_vehicle')),by=det_plate]
metadata <- metadata[det_plate %in% plates[vehicles>=5,det_plate]]
selected <- selected[det_plate %in% metadata$det_plate]
metadata[,column_index:=match(inst_id,h5read(path,'/0/META/COL/id'))]
setorder(metadata,column_index)
entrez <- as.character(h5read(path,'/0/META/ROW/id'))
symbols <- unname(mapIds(org.Hs.eg.db,keys=entrez,keytype='ENTREZID',column='SYMBOL',multiVals='first'))
stopifnot(all(genes %in% symbols),!anyNA(metadata$column_index))
expression <- matrix(NA_real_,978,nrow(metadata),dimnames=list(symbols,metadata$inst_id))
for(start in seq(1L,nrow(metadata),by=512L)) {
  keep <- start:min(start+511L,nrow(metadata))
  ix <- metadata$column_index[keep]
  block <- h5read(path,'/0/DATA/0/matrix',start=c(1L,min(ix)),count=c(978L,max(ix)-min(ix)+1L))
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
scores <- as.numeric(crossprod(weights,expression[genes,,drop=FALSE]))
observed <- selected[,.(depmap_id,drug=pert_iname,dose=pert_dose,det_plate,inst_id)]
observed[,HGF_RNA_coordinate:=scores[match(inst_id,metadata$inst_id)]]
for(dose in c(.04,.12,.37,1.11,3.33,10))observed[abs(log(observed$dose/dose))<=.035,dose:=dose]
curves <- observed[,.(HGF_RNA_coordinate=mean(HGF_RNA_coordinate)),by=.(depmap_id,drug,dose,det_plate)]
curves <- curves[,.(HGF_RNA_coordinate=mean(HGF_RNA_coordinate)),by=.(depmap_id,drug,dose)]
native <- curves[,{
  target <- if(drug[1]=='PLX-4720')2 else 1
  stopifnot(min(dose)<=target,max(dose)>=target)
  data.table(dose=target,HGF_RNA_coordinate=approx(log(dose),HGF_RNA_coordinate,xout=log(target))$y)
},by=.(depmap_id,drug)]
write_table(native,'external_exact_drug_HGF_coordinate')

# All ligand coordinates use full gene movements, the same transport and fixed
# loading as their40 scores. Existing HGF source views are kept individually.
meta <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/cytosig_full_metadata.tsv')
meta[,treatment:=toupper(trimws(treatment))]
ext <- fread('/work/results/program_pilot/111_external_bypass_combination_transfer/linear/data/external_ligand_programs.tsv')
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
paths <- c('/work/src/program_pilot/121_prepare_HGF_consensus_coordinate.R',
 file.path(root,'data/HGF_consensus_coordinate_definition.qs'),path,
 unique(source_rows$gene_cache))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'coordinate_source_hashes')
cat('Additional RNA coordinate prepared; existing40 and all functional outputs are unchanged.\n')
