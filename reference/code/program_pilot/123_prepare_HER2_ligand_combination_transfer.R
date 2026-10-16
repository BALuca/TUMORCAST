#!/usr/bin/env Rscript
# Independent paired drug/ligand phenotypes. No AU565 outcomes enter preparation.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(xml2);library(qs);library(digest)})
setDTthreads(1L)
root <- '/work/results/program_pilot/123_HER2_ligand_combination_transfer/linear'
for(folder in c('data/source','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
names_source <- c('elife-32271-fig1-data1-v1.pzfx','elife-32271-fig1-figsupp1-data1-v1.pzfx')
paths <- file.path(root,'data/source',names_source)
for(i in seq_along(paths))download.file(paste0('https://cdn.elifesciences.org/articles/32271/',names_source[i]),paths[i],quiet=TRUE,mode='wb')
documents <- lapply(paths,read_xml)
# Only unmodified parental2D assays: exclude72h AU565, transfected cells,3D,
# and additional inhibitors. Same study's MCF7 duplicate table is read only once.
selected <- data.table(file=c(1,1,2,2,2,2,2,2),table=c('Table11','Table11','Table13','Table17','Table27','Table29','Table19','Table21'),
  CellLine=c('SKBR3','MCF7','BT474','ZR751','HCC1569','HCC1419','BT474','SKBR3'),
  gene=c(rep('NRG1',6),'EGF','EGF'),
  drug_column=c('SKBR3 LAP','MCF7 LAP','LAP','LAP','HCC1569 Lapatinib','Lapatinib','LAP','LAP'),
  combination_column=c('SKBR3 LAP NRG','MCF7 LAP NRG','LAP NRG','LAP NRG','HCC1569 LAP NRG','Lapatinib + Heregulin','LAP EGF','LAP EGF'))
annotation <- fread('/work/outputs/dataset_inventory/cache/depmap/Model.csv')
selected[,depmap_id:=annotation$ModelID[match(CellLine,annotation$StrippedCellLineName)]]
stopifnot(!anyNA(selected$depmap_id),!any(selected$depmap_id=='ACH-000248'))
parts <- list()
for(i in seq_len(nrow(selected))) {
  t <- xml_find_first(documents[[selected$file[i]]],paste0('//Table[@ID="',selected$table[i],'"]'))
  doses <- as.numeric(xml_text(xml_find_all(t,'./XColumn/Subcolumn/d')))
  # First abscissa0.0001 is the study's graph-only vehicle placeholder; its
  # matcheddrug/ligand-free measurements define0-dose control, not0.0001uM drug.
  stopifnot(doses[1]==.0001)
  get_matrix <- function(title) {
    y <- xml_find_first(t,paste0('./YColumn[Title="',title,'"]'))
    stopifnot(!inherits(y,'xml_missing'))
    z <- lapply(xml_find_all(y,'./Subcolumn'),function(s)as.numeric(xml_text(xml_find_all(s,'./d'))))
    stopifnot(all(lengths(z)==8L),length(z)>=3L)
    do.call(cbind,z)
  }
  a <- get_matrix(selected$drug_column[i]);b <- get_matrix(selected$combination_column[i])
  stopifnot(all(a>0),all(b>0))
  for(ds in c(.5,5)) {
    # Match both measured outcome curves to native Tahoe doses in log dose and
    # log abundance. Endpoints and molecular vectors are never extrapolated.
    stopifnot(ncol(a)==ncol(b))
    aa <- vapply(seq_len(ncol(a)),function(k)approx(log(doses[2:8]),log2(a[2:8,k]),xout=log(ds))$y,numeric(1))
    bb <- vapply(seq_len(ncol(b)),function(k)approx(log(doses[2:8]),log2(b[2:8,k]),xout=log(ds))$y,numeric(1))
    growth <- log2(b[1,]/a[1,])
    departure <- bb-aa-growth
    parts[[paste(i,ds)]] <- data.table(CellLine=selected$CellLine[i],depmap_id=selected$depmap_id[i],
      gene=selected$gene[i],drug='Lapatinib ditosylate',dose=ds,
      observed_departure_log2=mean(departure),replicate_departure_SD=sd(departure),
      observed_ligand_log2=mean(growth),observed_drug_arm_increment_log2=mean(bb-aa),
      observed_drug_pct_vehicle=mean(100*2^aa/a[1,]),
      observed_combination_pct_vehicle=mean(100*2^bb/a[1,]),
      bio_replicates=ncol(a),study='elife32271',ligand_nM=10,time_h=72,
      source_table=selected$table[i],source_file=paths[selected$file[i]],
      outcome_match='within-replicate log-abundance/log-dose interpolation; nativeTahoe0.5/5uM')
  }
}
new_outcomes <- rbindlist(parts)
stopifnot(nrow(new_outcomes)==16L,!anyDuplicated(new_outcomes[,.(depmap_id,gene,drug,dose)]))
write_table(new_outcomes,'HER2_ligand_independent_outcomes')
old <- fread('/work/results/program_pilot/111_external_bypass_combination_transfer/linear/data/external_combination_outcomes.tsv')
old[,study:='Straussman2012']
write_table(rbindlist(list(old,new_outcomes),fill=TRUE),'expanded_external_combination_outcomes')

# Transfer the same pre-existingCytoSig gene-logFC recipes into NEW receivers.
# This is a characterized-ligand prediction, not observed combination RNA.
meta <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/cytosig_full_metadata.tsv')
meta[,treatment:=toupper(trimws(treatment))]
meta <- meta[treatment %in% c('NRG1','EGF') & !grepl('AU[ _-]?565|HCC[ _-]?1954',paste(profile_id,condition),ignore.case=TRUE)]
meta[is.na(condition)|condition=='',condition:='Unspecified context']
meta[,context_key:=paste(treatment,study,condition,sep='::')]
cy <- fread('/work/data/01_CytoSig_data/diff.merge.gz',select=c('V1',meta$profile_id),showProgress=FALSE)
setnames(cy,1,'gene')
FC <- vapply(c('NRG1','EGF'),function(lg) {
  keys <- unique(meta[treatment==lg,context_key])
  CX <- vapply(keys,function(key) {
    columns <- meta[context_key==key,profile_id]
    v <- rowMeans(as.matrix(cy[,..columns]),na.rm=TRUE);v[!is.finite(v)]<-0;v
  },numeric(nrow(cy)))
  studies <- meta[match(keys,context_key),study]
  SX <- vapply(unique(studies),function(st)rowMeans(CX[,studies==st,drop=FALSE]),numeric(nrow(cy)))
  rowMeans(SX)
},numeric(nrow(cy)))
rownames(FC)<-cy$gene
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
W <- archive$prep$W
P <- colnames(W)
path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
connection <- file(path,'r');header<-readLines(connection,1);rows<-list()
repeat {
  lines<-readLines(connection,100,warn=FALSE);if(!length(lines))break
  ids<-gsub('"','',sub(',.*','',lines))
  for(j in which(ids %in% selected$depmap_id))rows[[ids[j]]]<-lines[j]
  if(length(rows)==uniqueN(selected$depmap_id))break
}
close(connection)
inputs <- list()
for(cell in unique(selected$CellLine)) {
  id<-selected[CellLine==cell,depmap_id][1]
  b<-fread(text=paste(header,rows[[id]],sep='\n'))
  genes<-sub(' \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$','',names(b)[-1],perl=TRUE)
  raw<-data.table(gene=genes,u=pmax(2^as.numeric(b[,-1])-1,0))[,.(u=mean(u)),by=gene]
  u<-raw$u;reference<-u/sum(u)*1e4;index<-match(raw$gene,rownames(FC));gi<-match(rownames(W),raw$gene)
  stopifnot(!anyNA(gi))
  for(lg in unique(selected[CellLine==cell,gene])) {
    fc<-FC[index,lg];fc[!is.finite(fc)]<-0
    treated<-u*2^fc;treated<-treated/sum(treated)*1e4
    ok<-treated>=.1&reference>=.1;mv<-log2(treated[ok]/reference[ok]);av<-(log2(treated[ok])+log2(reference[ok]))/2
    ml<-quantile(mv,c(.3,.7));al<-quantile(av,c(.05,.95));keep<-mv>=ml[1]&mv<=ml[2]&av>=al[1]&av<=al[2]
    factor<-2^median(mv[keep]);z<-as.numeric((log2(1+treated[gi]/factor)-log2(1+reference[gi])) %*% W)
    inputs[[paste(cell,lg)]]<-cbind(data.table(CellLine=cell,depmap_id=id,gene=lg),as.data.table(as.list(setNames(z,P))))
  }
}
new_inputs<-rbindlist(inputs)
old_inputs<-fread('/work/results/program_pilot/118_time_matched_HGF_bypass_transfer/linear/data/time_matched_external_ligand_programs.tsv')
write_table(rbindlist(list(old_inputs,new_inputs),fill=TRUE),'expanded_external_ligand_programs')
write_table(new_inputs,'HER2_ligand_transported_programs')
write_table(meta,'HER2_ligand_RNA_provenance')
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'HER2_source_hashes')
print(new_outcomes[,.(CellLine,gene,dose,observed_departure_log2,observed_combination_pct_vehicle)])
cat('16 independent combinations prepared; all AU565 outcomes excluded.\n')
