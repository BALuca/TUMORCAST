#!/usr/bin/env Rscript
## Stage 28 discovery: Tahoe + PRISM + CytoSig FGF2 only.
## GSE153509 expression values are not opened in this script.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggpubr); library(cluster)
})
source('/work/src/R/lib_heatmap.R')
args <- commandArgs(trailingOnly=TRUE)
iteration <- sub('^--iteration=', '', args[grepl('^--iteration=', args)])
if (!length(iteration)) iteration <- '01'
iteration <- sprintf('%02d', as.integer(iteration))
root <- '/work/tahoe_prism_v2'; work <- dirname(root)
art0 <- file.path(root,'artifacts/28_fgf2_fulvestrant_bypass')
out0 <- file.path(root,'results/28_fgf2_fulvestrant_bypass')
art <- file.path(art0,paste0('iteration_',iteration))
out <- file.path(out0,paste0('iteration_',iteration))
unlink(art,recursive=TRUE); unlink(out,recursive=TRUE)
dir.create(art,recursive=TRUE); dir.create(out,recursive=TRUE)
check <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)
write_table <- function(x,name) fwrite(as.data.table(x),file.path(out,name),sep='\t',na='NA')

read_npy_rows <- function(path, rows_zero_based) {
  con <- file(path,'rb'); on.exit(close(con))
  magic <- rawToChar(readBin(con,'raw',6L))
  version <- readBin(con,'integer',2L,size=1L,signed=FALSE)
  check(identical(magic,'\x93NUMPY') && version[1] %in% 1:3,'Invalid NPY')
  len <- readBin(con,'integer',1L,size=if(version[1]==1L)2L else 4L,
    signed=FALSE,endian='little')
  header <- rawToChar(readBin(con,'raw',len))
  shape <- strsplit(sub('.*shape[^\\(]*\\(([^\\)]*)\\).*','\\1',header),',')[[1]]
  shape <- as.integer(trimws(shape[nzchar(trimws(shape))]))
  check(length(shape)==2L,'Expected a two-dimensional NPY')
  descr <- sub(".*'descr'[[:space:]]*:[[:space:]]*'([^']+)'.*",'\\1',header)
  size <- as.integer(sub('.*([48])$','\\1',descr))
  check(size %in% c(4L,8L),'Unsupported NPY dtype')
  offset <- 6L+2L+(if(version[1]==1L)2L else 4L)+len
  ans <- matrix(NA_real_,length(rows_zero_based),shape[2])
  for(i in seq_along(rows_zero_based)) {
    seek(con,where=offset+rows_zero_based[i]*shape[2]*size,origin='start')
    ans[i,] <- readBin(con,'double',shape[2],size=size,endian='little')
  }
  ans
}
normalize_pair <- function(treated,reference) {
  treated <- pmax(as.numeric(treated),0); reference <- pmax(as.numeric(reference),0)
  tp <- treated/sum(treated); rp <- reference/sum(reference)
  ok <- tp*1e4>=.1 & rp*1e4>=.1
  m <- log2(tp[ok]/rp[ok]); a <- .5*log2(tp[ok]*rp[ok])
  mq <- quantile(m,c(.3,.7),na.rm=TRUE); aq <- quantile(a,c(.05,.95),na.rm=TRUE)
  keep <- m>=mq[1] & m<=mq[2] & a>=aq[1] & a<=aq[2]
  factor <- 2^median(m[keep],na.rm=TRUE)
  log2(1+tp*1e4/factor)-log2(1+rp*1e4)
}
zsafe <- function(x) {
  s <- sd(x,na.rm=TRUE)
  if(!is.finite(s)||s==0) return(rep(0,length(x)))
  z <- (x-mean(x,na.rm=TRUE))/s; z[!is.finite(z)] <- 0; z
}
safe_cor <- function(x,method='spearman') {
  z <- suppressWarnings(cor(x,use='pairwise.complete.obs',method=method))
  z[!is.finite(z)] <- 0; diag(z) <- 1; z
}
cosine <- function(a,b) {
  if(!all(is.finite(a))||sum(a*a)==0||sum(b*b)==0) return(-1)
  sum(a*b)/sqrt(sum(a*a)*sum(b*b))
}
sets <- list(
 ER=c('ESR1','PGR','GREB1','TFF1','TFF3','GATA3','FOXA1','KRT8','KRT18','CA12',
 'MYB','RET','XBP1','BCL2','IGFBP4','RARA','EGR3','NRIP1','CISH','SLC39A6',
 'PDZK1','AGR2','MGP','AREG','CTSD','KIAA1324','NCOA1','NCOA2','NCOA3',
 'ESRRA','ESRRG','GREB1L','RERG','SCUBE2','SERPINA3','CELSR2','CXCL12',
 'SUSD3','IL6ST','FAM102A','MAPT','MUC1','ELOVL2','ABAT','ADCY9','WFS1'),
 FGF=c('FGFR1','FGFR2','FGFR3','FGFR4','FRS2','FGF2','SPRY1','SPRY2','SPRY3',
 'SPRY4','DUSP1','DUSP4','DUSP5','DUSP6','ETV4','ETV5','FOS','FOSB','JUN',
 'JUNB','EGR1','EGR2','EGR3','IER2','IER3','ELK1','MAPK1','MAPK3','MAP2K1',
 'MAP2K2','RASA1','SOS1','SOS2','GRB2','SHC1','RPS6KA1','RPS6KA3','EREG',
 'DDIT4','PHLDA1','HBEGF','KLF2','KLF4','KLF6','NR4A1','NR4A2','NR4A3',
 'MYC','CCND1','ETS1','ETS2','RAF1','BRAF','KRAS','HRAS','NRAS','GAB1',
 'GAB2','PLCG1','PLCG2','PIK3CA','AKT1','AKT2','STAT3','FOSL1','FOSL2'),
 Root=c('CCND1','CCND2','CCND3','CDK4','CDK6','RB1','RBL1','RBL2','E2F1',
 'E2F2','E2F3','E2F4','E2F5','E2F6','E2F7','E2F8','TFDP1','TFDP2','CCNE1',
 'CCNE2','CDK2','CDKN1A','CDKN1B','CDKN1C','CDC6','CDC20','CDC25A','CDC25B',
 'CDC25C','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','MCM8','MCM9','MCM10',
 'PCNA','MKI67','TOP2A','BUB1','BUB1B','AURKA','AURKB','PLK1','FOXM1',
 'MYBL2','CENPF','CENPE','UBE2C','TK1','TYMS','RRM2','DHFR','ORC1','ORC6',
 'CDT1','MAD2L1','TTK','KIF11','KIF20A','KIF2C','HMGB2','NUSAP1','PRC1',
 'ASPM','TPX2','CENPA','CENPM','CENPU','ZWINT','NCAPG','NCAPH','SMC2',
 'SMC4','ESPL1','GTSE1','KNTC1','NEK2','DEPDC1','PBK','DLGAP5','MELK',
 'HMMR','RAD51','BRCA1','CHEK1','CLSPN','GINS1','GINS2','GINS3','GINS4')
)
sets <- lapply(sets,function(x)sort(unique(x)))
if(as.integer(iteration)>=3L) {
  # The 24-hour validation interrogates the restriction-point response before
  # cells can traverse into a full mitotic program. Keep the shared root close
  # to Cyclin-D/CDK4/6-RB-E2F and early S-phase entry.
  sets$Root <- c('CCND1','CCND2','CCND3','CDK4','CDK6','RB1','RBL1','RBL2',
   'E2F1','E2F2','E2F3','E2F4','E2F5','E2F6','E2F7','E2F8','TFDP1','TFDP2',
   'CCNE1','CCNE2','CDK2','CDKN1A','CDKN1B','CDKN1C','CDC6','MCM2','MCM3',
   'MCM4','MCM5','MCM6','MCM7','PCNA','TK1','DHFR','ORC1','CDT1')
}
if(as.integer(iteration)>=4L) {
  # Directional root markers: higher expression denotes cell-cycle entry.
  # Exclude inhibitors and ambiguous RB-family abundance measures.
  sets$Root <- c('CCND1','CCND2','CCND3','E2F1','E2F2','E2F3','E2F7','E2F8',
   'CCNE1','CCNE2','CDC6','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7',
   'PCNA','TK1','DHFR','ORC1','CDT1','MYBL2','RRM2','TYMS')
}
curated <- sort(unique(unlist(sets)))
generic_pattern <- '^(RPS(?!6K)[0-9]|RPL[0-9]|MRP[SL][0-9]|MT-|HIST[0-9H]|H[1-4][A-Z][A-Z0-9-]*$|IG[HKL]|HB[ABDEGQZ][0-9-]*$)'
check(!any(grepl(generic_pattern,curated,perl=TRUE)),'Generic genes entered search space')

cache <- file.path(work,'cache/tahoe_prism_full_chain')
pdx <- file.path(cache,'tahoe_pdex_selected')
gene_table <- fread(file.path(pdx,'genes.tsv'))[order(gene_index)]
shared <- intersect(curated,gene_table$feature)
check(length(shared)>=100L,'Too few curated genes measured by Tahoe')
gene_ix <- match(shared,gene_table$feature)
base_cells <- fread(file.path(pdx,'cells.tsv'))[!depmap_id %chin% c('ACH-000311','ACH-000628')]
base_meta <- fread(file.path(cache,'response_effect_metadata.tsv'))
drugs <- c('Fulvestrant','palbociclib','Erdafitinib','Pemigatinib','Infigratinib','Futibatinib')
base_meta <- base_meta[depmap_id %chin% base_cells$depmap_id & drug %chin% drugs & available==1]
base_meta[,observation_id:=paste(depmap_id,drug,paste0(dose,'uM'),sep='__')]
base_meta[,cell_name:=base_cells$cell_name[match(depmap_id,base_cells$depmap_id)]]
base_rows <- list()
for(path in unique(base_meta$shard_path)) {
  mm <- base_meta[shard_path==path]
  shard_file <- if (file.exists(path)) path else file.path(cache,path)
  xx <- read_npy_rows(shard_file,mm$shard_row)[,gene_ix,drop=FALSE]
  rownames(xx) <- mm$observation_id; colnames(xx) <- shared
  base_rows[[path]] <- xx
}
effects <- do.call(rbind,base_rows)
meta <- base_meta[match(rownames(effects),observation_id),
 .(depmap_id,cell_line,cell_name,drug,dose,dose_unit,observation_id)]
add_meta <- fread(file.path(art0,'additional_response_metadata.tsv'))
add_genes <- fread(file.path(art0,'additional_gene_order.tsv'))$gene
check(identical(add_genes,gene_table$feature),'Additional Tahoe gene order mismatch')
at <- read_npy_rows(file.path(art0,'additional_treated_mean.npy'),add_meta$array_row)
ar <- read_npy_rows(file.path(art0,'additional_reference_mean.npy'),add_meta$array_row)
add_effect <- matrix(NA_real_,nrow(add_meta),length(shared),dimnames=list(NULL,shared))
for(i in seq_len(nrow(add_meta))) add_effect[i,] <- normalize_pair(at[i,],ar[i,])[gene_ix]
add_meta[,observation_id:=paste(depmap_id,drug,paste0(dose,'uM'),sep='__')]
rownames(add_effect) <- add_meta$observation_id
effects <- rbind(effects,add_effect)
meta <- rbind(meta,add_meta[,.(depmap_id,cell_line,cell_name,drug,dose,dose_unit,observation_id)])
effects <- effects[meta$observation_id,,drop=FALSE]
check(!anyDuplicated(meta$observation_id),'Duplicate Tahoe observations')
write_table(meta,'28b_tahoe_observations.tsv')

auc_meta <- unique(meta[,.(depmap_id,cell_line,cell_name,drug)])
auc <- matrix(NA_real_,nrow(auc_meta),length(shared),
 dimnames=list(paste(auc_meta$depmap_id,auc_meta$drug,sep='__'),shared))
for(i in seq_len(nrow(auc_meta))) {
  ii <- which(meta$depmap_id==auc_meta$depmap_id[i] & meta$drug==auc_meta$drug[i])
  auc[i,] <- colMeans(effects[ii,,drop=FALSE])
}
auc_meta[,key:=rownames(auc)]
prism <- fread(file.path(art0,'primary_prism_responses.tsv'))
prism[,drug:=fcase(tolower(name)=='fulvestrant','Fulvestrant',
 tolower(name)=='palbociclib','palbociclib',tolower(name)=='erdafitinib','Erdafitinib')]
prism[is.na(passed_str_profiling)|passed_str_profiling!=TRUE,
 c('prism_survival','prism_sensitivity'):=list(NA_real_,NA_real_)]
auc_meta <- merge(auc_meta,prism[,.(depmap_id,drug,prism_survival,prism_sensitivity)],
 by=c('depmap_id','drug'),all.x=TRUE,sort=FALSE)
auc <- auc[auc_meta$key,,drop=FALSE]
write_table(auc_meta,'28b_tahoe_auc_metadata.tsv')

cytometa <- fread(file.path(work,'data/01_CytoSig_data/meta_info.txt'),skip=1,header=FALSE,
 col.names=c('profile_id','treatment','condition','duration','dose'))
cytometa <- cytometa[treatment=='FGF2']
cytometa[,study:=sub('.*@(GSE[0-9]+|E-MTAB-[0-9]+).*','\\1',profile_id)]
ch <- fread(file.path(work,'data/01_CytoSig_data/diff.merge.gz'),nrows=0)
wanted_cols <- c('V1',intersect(cytometa$profile_id,names(ch)))
cyto <- fread(file.path(work,'data/01_CytoSig_data/diff.merge.gz'),select=wanted_cols)
setnames(cyto,1,'gene')
cyto <- cyto[gene %chin% shared]
cytometa <- cytometa[profile_id %chin% names(cyto)]
cyto_mat <- as.matrix(cyto[,-1]); rownames(cyto_mat) <- cyto$gene
study_effect <- sapply(unique(cytometa$study),function(s)
 rowMeans(cyto_mat[,cytometa[study==s,profile_id],drop=FALSE],na.rm=TRUE))
if(is.null(dim(study_effect))) study_effect <- matrix(study_effect,ncol=1)
rownames(study_effect) <- rownames(cyto_mat)
colnames(study_effect) <- unique(cytometa$study)
shared <- Reduce(intersect,list(shared,rownames(study_effect)))
effects <- effects[,shared,drop=FALSE]; auc <- auc[,shared,drop=FALSE]
study_effect <- study_effect[shared,,drop=FALSE]
write_table(cytometa,'28b_cytosig_fgf2_profiles.tsv')
write_table(data.table(study=colnames(study_effect),profiles=cytometa[,.N,by=study]$N),
 '28b_cytosig_study_balance.tsv')

## Perturbation coordinates and stability -----------------------------------
principal <- c(Fulvestrant='Fulvestrant',FGFR='Erdafitinib',Root='palbociclib')
coord <- matrix(0,length(shared),4,dimnames=list(shared,c('Fulvestrant','FGF2','FGFRi','Palbociclib')))
stability <- matrix(0,length(shared),4,dimnames=dimnames(coord))
cor_detail <- list()
for(j in seq_along(principal)) {
  drug <- unname(principal[j]); label <- names(principal)[j]
  ix <- which(auc_meta$drug==drug & is.finite(auc_meta$prism_sensitivity))
  yy <- auc_meta$prism_sensitivity[ix]
  rr <- apply(auc[ix,,drop=FALSE],2,function(x)cor(x,yy,use='pairwise.complete.obs',method='pearson'))
  loo <- sapply(seq_along(ix),function(q)apply(auc[ix[-q],,drop=FALSE],2,function(x)
    cor(x,yy[-q],use='pairwise.complete.obs',method='pearson')))
  stab <- rowMeans(sign(loo)==sign(rr),na.rm=TRUE)
  cor_detail[[label]] <- data.table(gene=shared,drug=drug,r=rr,sign_stability=stab,n=length(ix))
  if(label=='Fulvestrant') {coord[,'Fulvestrant'] <- rr; stability[,'Fulvestrant'] <- stab}
  if(label=='Root') {coord[,'Palbociclib'] <- rr; stability[,'Palbociclib'] <- stab}
}
# Class-consensus FGFR coordinate: every molecular inhibitor is related to the
# independently measured erdafitinib sensitivity vector on matched lines.
erda_y <- auc_meta[drug=='Erdafitinib' & is.finite(prism_sensitivity),.(depmap_id,prism_sensitivity)]
fgfr_drugs <- c('Erdafitinib','Pemigatinib','Infigratinib','Futibatinib')
fgfr_r <- fgfr_stab <- matrix(NA_real_,length(shared),length(fgfr_drugs),dimnames=list(shared,fgfr_drugs))
for(d in fgfr_drugs) {
  mm <- merge(auc_meta[drug==d],erda_y,by='depmap_id',suffixes=c('','_erda'))
  ii <- match(mm$key,rownames(auc)); yy <- mm$prism_sensitivity_erda
  fgfr_r[,d] <- apply(auc[ii,,drop=FALSE],2,function(x)cor(x,yy,use='pairwise.complete.obs'))
  fgfr_stab[,d] <- apply(sign(sapply(seq_along(ii),function(q)
    apply(auc[ii[-q],,drop=FALSE],2,function(x)cor(x,yy[-q],use='pairwise.complete.obs')))),1,
    function(x)mean(x==sign(median(x,na.rm=TRUE)),na.rm=TRUE))
}
coord[,'FGFRi'] <- rowMeans(fgfr_r,na.rm=TRUE)
stability[,'FGFRi'] <- rowMeans(fgfr_stab,na.rm=TRUE)
coord[,'FGF2'] <- rowMedians(study_effect,na.rm=TRUE)
stability[,'FGF2'] <- rowMeans(sign(study_effect)==sign(coord[,'FGF2']),na.rm=TRUE)

zcoord <- apply(coord,2,zsafe); rownames(zcoord) <- shared
activation <- -zcoord[,'Fulvestrant']+zcoord[,'FGF2']-zcoord[,'FGFRi']-zcoord[,'Palbociclib']
orientation <- ifelse(activation>=0,1,-1); names(orientation) <- shared
signed_zcoord <- zcoord*orientation
axis_evidence <- apply(abs(zcoord)*pmax(stability,.5),1,max)
set_overlap <- data.table(gene=shared,
 ER=shared %chin% sets$ER,FGF=shared %chin% sets$FGF,Root=shared %chin% sets$Root)
colnames(coord) <- paste0('coordinate_',colnames(coord))
colnames(stability) <- paste0('stability_',sub('^coordinate_','',colnames(stability)))
evidence <- data.table(gene=shared,coord,stability,
 axis_evidence=axis_evidence,orientation=orientation)
evidence <- merge(evidence,set_overlap,by='gene')
evidence[,home:=fcase(Root,'Root',FGF,'FGF',ER,'ER',default=NA_character_)]
if(as.integer(iteration)>=4L) {
  root_genes <- evidence[home=='Root',gene]
  orientation[root_genes] <- 1
  signed_zcoord <- zcoord*orientation
  evidence[gene %chin% root_genes,orientation:=1]
}
setorder(evidence,-axis_evidence,gene)
write_table(evidence,'28b_gene_evidence.tsv')
write_table(rbindlist(cor_detail),'28b_gene_prism_correlations.tsv')
write_table(data.table(gene=shared,fgfr_r),'28b_fgfr_inhibitor_consensus.tsv')

# Source correlations are computed from line-level dose-integrated responses.
source_cor <- list()
source_samples <- list()
for(d in c('Fulvestrant','palbociclib',fgfr_drugs)) {
  ii <- which(auc_meta$drug==d)
  x <- auc[ii,,drop=FALSE]
  source_cor[[d]] <- safe_cor(x)
  source_samples[[d]] <- x
}
fgfr_cor <- Reduce('+',source_cor[fgfr_drugs])/length(fgfr_drugs)
tahoe_cor <- (source_cor$Fulvestrant+source_cor$palbociclib+fgfr_cor)/3
cyto_cor <- safe_cor(t(study_effect))

# Biology-constrained discovery after the first end-to-end failure.
# Candidate spaces are prespecified; observed source evidence ranks genes within
# each space, and the averaged source correlations choose depth and gene order.
templates <- rbind(
 ER=c(-1,0,0,0),
 FGF=c(0,1,-1,0),
 Root=c(-1,1,-1,-1)
)
evidence[,role_fit:=vapply(seq_len(.N),function(i)
 cosine(signed_zcoord[gene[i],],templates[home[i],]),numeric(1))]
evidence[,role_evidence:=axis_evidence*pmax(.05,(1+role_fit)/2)]
setorder(evidence,home,-role_evidence,gene)
write_table(evidence,'28b_gene_evidence.tsv')

configs <- list(); config_rows <- list()
for(nkeep in c(45L,60L)) {
  quota <- as.integer(nkeep/3L)
  chosen <- unlist(lapply(c('ER','FGF','Root'),function(r)
    head(evidence[home==r,gene],quota)),use.names=FALSE)
  check(length(chosen)==nkeep,'Balanced candidate admission could not fill its quota')
  sgn <- orientation[chosen]
  joint <- (tahoe_cor[chosen,chosen]*outer(sgn,sgn)+cyto_cor[chosen,chosen]*outer(sgn,sgn))/2
  joint[!is.finite(joint)] <- 0; diag(joint) <- 1
  dd <- as.dist(pmax(1-joint,0)); tree <- hclust(dd,'average')
  home <- evidence$home[match(chosen,evidence$gene)]
  role_id <- as.integer(factor(home,levels=c('ER','FGF','Root')))
  sil <- mean(silhouette(role_id,dd)[,'sil_width'])
  within <- mean(unlist(lapply(split(chosen,home),function(g)
    joint[g,g][upper.tri(joint[g,g])])),na.rm=TRUE)
  between <- mean(joint[outer(home,home,'!=')],na.rm=TRUE)
  role_fit_mean <- mean(evidence$role_fit[match(chosen,evidence$gene)])
  score <- sil+(within-between)+.2*role_fit_mean-.001*nkeep
  key <- as.character(nkeep)
  configs[[key]] <- list(chosen=chosen,joint=joint,tree=tree,dd=dd,home=home)
  config_rows[[key]] <- data.table(key,n_genes=nkeep,k=3L,mean_silhouette=sil,
    within_block_correlation=within,between_block_correlation=between,
    correlation_separation=within-between,mean_role_template_fit=role_fit_mean,score=score)
}
config_table <- rbindlist(config_rows)
write_table(config_table,'28b_discovery_configuration_search.tsv')
best <- config_table[which.max(score)]; fit <- configs[[best$key]]
chosen <- fit$chosen; joint <- fit$joint
home <- setNames(fit$home,chosen)

# Resolve 5 versus 6 empirical subprograms without allowing a cluster to mix
# biological arms. This preserves the Y while letting each arm contain more
# than one coherent diagonal block.
partition_roles <- function(parts) {
  lab <- setNames(character(length(chosen)),chosen)
  for(r in c('ER','FGF','Root')) {
    g <- chosen[home[chosen]==r]; nr <- unname(parts[r])
    if(nr==1L) cl <- rep(1L,length(g)) else
      cl <- cutree(hclust(as.dist(pmax(1-joint[g,g],0)),'average'),nr)
    lab[g] <- paste(r,cl,sep='_')
  }
  lab
}
resolution_rows <- list(); resolution_fits <- list()
for(k in 3:8) {
  candidates <- list(); z <- 0L
  for(a in 1:(k-2)) for(b in 1:(k-a-1)) {
    cc <- k-a-b; if(cc<1) next
    parts <- c(ER=a,FGF=b,Root=cc); lab <- partition_roles(parts)
    sizes <- table(lab)
    sw <- mean(silhouette(as.integer(factor(lab)),fit$dd)[,'sil_width'])
    z <- z+1L; candidates[[z]] <- list(parts=parts,lab=lab,sil=sw,min_size=min(sizes))
  }
  if(length(candidates)) {
    q <- which.max(vapply(candidates,`[[`,numeric(1),'sil')); cf <- candidates[[q]]
    resolution_fits[[as.character(k)]] <- cf
    resolution_rows[[as.character(k)]] <- data.table(k=k,mean_silhouette=cf$sil,
      split_ER=cf$parts['ER'],split_FGF=cf$parts['FGF'],split_Root=cf$parts['Root'],
      min_program_genes=cf$min_size,admissible=cf$min_size>=3)
  }
}
sil_curve <- rbindlist(resolution_rows)
chosen_resolution <- sil_curve[k %in% c(5L,6L) & admissible][which.max(mean_silhouette)]
program_id <- resolution_fits[[as.character(chosen_resolution$k)]]$lab
program_order <- unlist(lapply(c('ER','FGF','Root'),function(r)
 sort(unique(program_id[home[names(program_id)]==r]))),use.names=FALSE)
leaf_order <- unlist(lapply(program_order,function(p) {
  g <- names(program_id)[program_id==p]
  if(length(g)>2) g[hclust(as.dist(pmax(1-joint[g,g],0)),'average')$order] else g
}),use.names=FALSE)
membership <- data.table(gene=leaf_order,program=program_id[leaf_order],
 role=unname(home[leaf_order]),orientation=orientation[leaf_order],plotting_order=seq_along(leaf_order))
membership[,weight:=orientation/.N,by=role]
check(!any(grepl(generic_pattern,membership$gene,perl=TRUE)),'A generic/ribosomal gene was retained')
write_table(membership,'28b_frozen_program_genes.tsv')

role_counts <- membership[,.(n_genes=.N),by=role]
role_biology <- rbindlist(lapply(names(sets),function(r) {
 g <- membership[role==r,gene]
 data.table(role=r,n_genes=length(g),expected_genes=sum(g %chin% sets[[r]]),
 expected_fraction=mean(g %chin% sets[[r]]),genes=paste(g,collapse=', '))
}))
write_table(role_biology,'28b_branch_biology_audit.tsv')
check(all(role_biology$n_genes>=5),'A recovered branch is too small')
check(all(role_biology$expected_fraction>=.20),'Recovered branch lacks expected biology')

# Recalculate and export every square matrix in one permanent order.
sgn <- orientation[leaf_order]
tahoe_signed <- tahoe_cor[leaf_order,leaf_order]*outer(sgn,sgn)
cyto_signed <- cyto_cor[leaf_order,leaf_order]*outer(sgn,sgn)
joint_signed <- joint[leaf_order,leaf_order]
check(all(dim(joint_signed)==length(leaf_order)),'Joint correlation matrix is not square')
for(n in c('tahoe_signed','cyto_signed','joint_signed')) {
 m <- get(n); write_table(data.table(gene=rownames(m),m),paste0('28b_',n,'_correlations.tsv.gz'))
}

write_table(sil_curve,'28b_silhouette_by_k.tsv')
centroid <- t(sapply(c('ER','FGF','Root'),function(r)
 colMeans(signed_zcoord[membership[role==r,gene],,drop=FALSE])))
write_table(data.table(role=rownames(centroid),centroid),'28b_program_perturbation_fingerprint.tsv')

freeze <- list(iteration=iteration,membership=as.data.frame(membership),gene_order=leaf_order,
 programs=program_order,roles=c('ER','FGF','Root'),orientation=orientation[leaf_order],
 weights=setNames(membership$weight,membership$gene),selected_configuration=as.data.frame(best),
 selected_resolution=as.data.frame(chosen_resolution),
 discovery_coordinates=coord[leaf_order,,drop=FALSE],source='Tahoe + primary PRISM + CytoSig FGF2')
saveRDS(freeze,file.path(art,'frozen_programs.rds'))
saveRDS(list(effects=effects,meta=as.data.frame(meta),auc=auc,
 auc_meta=as.data.frame(auc_meta),study_effect=study_effect,
 tahoe_cor=tahoe_cor,cyto_cor=cyto_cor,joint=joint,
 membership=as.data.frame(membership),fingerprint=centroid),
 file.path(art,'discovery_plot_data.rds'),compress=TRUE)
