#!/usr/bin/env Rscript
## Stage 27: repair BRAF--ERBB discovery around the actual bypass hypothesis.
## BRAF target engagement may occur in BOTH sensitive and bypassing cells.
## Its markers must not be required to correlate individually with viability.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cluster)})
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(arg))dirname(dirname(normalizePath(sub('^--file=','',arg[1])))) else '/work/tahoe_prism_v2'
work<-dirname(root);input_art<-file.path(root,'artifacts/26_braf_erbb_bypass')
art<-file.path(root,'artifacts/27_braf_erbb_hierarchy');out<-file.path(root,'results/27_braf_erbb_hierarchy')
dir.create(art,recursive=TRUE,showWarnings=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(work,'src/R/lib_heatmap.R'))
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
write_table<-function(x,name)fwrite(as.data.table(x),file.path(out,name),sep='\t',na='NA')
EXCLUDED_IDS<-c('ACH-000311','ACH-000628');N_PER_ROLE<-40L;K<-6L


read_npy_rows <- function(path, rows_zero_based) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  stopifnot(identical(magic, "\x93NUMPY"), version[1] %in% 1:3)
  len <- readBin(con, "integer", 1L, size = if (version[1] == 1L) 2L else 4L,
    signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", len))
  stopifnot(grepl("'descr': '<f4'", header, fixed = TRUE),
    grepl("'fortran_order': False", header, fixed = TRUE))
  shape <- as.integer(trimws(strsplit(sub(".*'shape': \\(([^)]*)\\).*", "\\1", header), ",")[[1L]]))
  shape <- shape[!is.na(shape)]
  stopifnot(length(shape) == 2L, all(rows_zero_based >= 0 & rows_zero_based < shape[1]))
  offset <- seek(con)
  answer <- matrix(NA_real_, length(rows_zero_based), shape[2])
  for (i in seq_along(rows_zero_based)) {
    seek(con, offset + 4 * shape[2] * rows_zero_based[i], origin = "start")
    v <- readBin(con, "numeric", n = shape[2], size = 4L, endian = "little")
    stopifnot(length(v) == shape[2]); answer[i, ] <- v
  }
  answer
}
## Apply the existing paired, sample-local composition normalization to the
## saved processed pseudobulks. All 62,710 genes remain in the expression output.
## This is a measurement transformation, not a gene or cell-line eligibility gate.
normalize_pair <- function(treated,reference) {
  tol <- max(1e-6,32*2^-23*max(abs(c(treated,reference)),1))
  check(min(c(treated,reference))>=-tol,'Material negative reconstructed abundance')
  treated <- pmax(treated,0);reference <- pmax(reference,0)
  tt <- sum(treated); rt <- sum(reference)
  check(tt>0 && rt>0,'Empty pseudobulk')
  tp <- treated/tt; rp <- reference/rt
  eligible_scale <- tp*10000>=.1 & rp*10000>=.1
  m <- log2(tp[eligible_scale]/rp[eligible_scale]);a <- .5*log2(tp[eligible_scale]*rp[eligible_scale])
  mq <- quantile(m,c(.3,.7)); aq <- quantile(a,c(.05,.95))
  kept <- m>=mq[1]&m<=mq[2]&a>=aq[1]&a<=aq[2]
  check(sum(kept)>0,'Cannot estimate sample-local composition scaling')
  factor <- 2^median(m[kept])
  list(effect=log2(1+tp*10000/factor)-log2(1+rp*10000),
    treated_total=tt,reference_total=rt,factor=factor)
}

plot_save <- function(p,name,w=7,h=5) {
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='white')
  ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=140,bg='white')
}

heat <- function(x,name,title,rows=NULL,cols=NULL,limit=2.5,width=10,height=NULL,show_cols=TRUE) {
  if(is.null(height))height<-max(5,nrow(x)*.14+2)
  if(!is.null(rows))check(length(rows)==nrow(x),'Heatmap row annotations misaligned')
  if(!is.null(cols))check(length(cols)==ncol(x),'Heatmap column annotations misaligned')
  h<-draw_heatmap(x,scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    cluster_rows=FALSE,cluster_columns=FALSE,row_split=rows,column_split=cols,
    cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(2,'mm'),column_gap=grid::unit(2,'mm'),
    show_row_names=TRUE,show_column_names=show_cols,row_names_gp=grid::gpar(fontsize=9),column_names_gp=grid::gpar(fontsize=8),
    column_names_rot=90,hm_name=if(grepl('correlation',name))'Spearman rho' else 'Value',heatmap_legend_param=list(at=c(-limit,0,limit),labels=round(c(-limit,0,limit),2)),title=title)
  cairo_pdf(file.path(out,paste0(name,'.pdf')),width=width,height=height)
  drawn<-ComplexHeatmap::draw(h,heatmap_legend_side='bottom')
  check(identical(as.integer(unlist(ComplexHeatmap::row_order(drawn))),seq_len(nrow(x))),'Heatmap changed frozen gene order')
  check(identical(as.integer(unlist(ComplexHeatmap::column_order(drawn))),seq_len(ncol(x))),'Heatmap changed frozen sample order')
  invisible(dev.off())
  png(file.path(out,paste0(name,'.png')),width=width,height=height,units='in',res=120,type='cairo')
  ComplexHeatmap::draw(h,heatmap_legend_side='bottom');invisible(dev.off())
}

source_heat <- function(x,name,title) {
  x<-x[leaf_order,,drop=FALSE]
  c<-cor(x,method='spearman',use='pairwise.complete.obs');check(all(is.finite(c)),'Non-estimable sample correlation')
  tr<-hclust(as.dist(pmax(1-c,0)),method='average')
  z<-t(scale(t(x)));z[!is.finite(z)]<-NA_real_
  signed_z<-z*gene_sign[rownames(z)]
  activity<-sapply(program_order,function(p)colMeans(signed_z[membership[program==p,gene],,drop=FALSE],na.rm=TRUE))
  winner<-max.col(activity,ties.method='first')
  ord<-unlist(lapply(seq_along(program_order),function(j)tr$order[winner[tr$order]==j]))
  cb<-factor(program_order[winner[ord]],levels=program_order)
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),dominant_program=as.character(cb)),paste0(name,'_sample_order.tsv'))
  heat(z[,ord,drop=FALSE],name,title,row_blocks,cb,limit=2.5,width=12)
  heat(signed_z[,ord,drop=FALSE],paste0(name,'_signed'),paste0(title,' | signed activity'),row_blocks,cb,limit=2.5,width=12)
  invisible(ord)
}

get_soft_samples <- function(s) {
  starts <- grep('^\\^SAMPLE = ',s)
  rbindlist(lapply(seq_along(starts),function(i) {
    block <- s[starts[i]:(if(i==length(starts)) length(s) else starts[i+1]-1L)]
    field <- function(prefix) {
      x <- block[startsWith(block,prefix)]
      check(length(x)==1L,paste('Ambiguous SOFT field',prefix))
      substring(x,nchar(prefix)+1L)
    }
    data.table(gsm=field('^SAMPLE = '),title=field('!Sample_title = '),
      characteristics=paste(block[startsWith(block,'!Sample_characteristics_ch1 = ')],collapse=' | '),
      protocol=paste(block[startsWith(block,'!Sample_treatment_protocol_ch1 = ')],collapse=' | '))
  }))
}

## 1. Extend the molecular context using processed data, not raw cells --------
required<-file.path(art,c('additional_treated_mean.npy','additional_reference_mean.npy','additional_response_metadata.tsv','primary_prism_responses.tsv'))
check(all(file.exists(required)),'Run scripts/27_prepare_braf_context.py first to retrieve the two processed melanoma profiles')
a<-readRDS(file.path(input_art,'discovery_measurements.rds'))
effects<-a$effects;meta<-a$meta;cells<-a$cells;genes<-a$genes
reference_cp10k<-a$reference_cp10k;treated_cp10k<-a$treated_cp10k
cyto<-a$cyto;cytometa<-a$cytometa
new_meta<-fread(required[3]);new_treated<-read_npy_rows(required[1],new_meta$array_row);new_reference<-read_npy_rows(required[2],new_meta$array_row)
check(identical(fread(file.path(art,'additional_gene_order.tsv'))$gene,colnames(effects)),'Additional gene order differs')
new_effect<-new_treated;new_ref_cp<-new_treated;new_treat_cp<-new_treated;new_qc<-list()
for(i in seq_len(nrow(new_meta))) {
  z<-normalize_pair(new_treated[i,],new_reference[i,]);new_effect[i,]<-z$effect
  new_ref_cp[i,]<-new_reference[i,]/z$reference_total*10000
  new_treat_cp[i,]<-new_treated[i,]/z$treated_total*10000/z$factor
  new_qc[[i]]<-data.table(depmap_id=new_meta$depmap_id[i],cell_name=new_meta$cell_name[i],drug=new_meta$drug[i],dose=new_meta$dose[i],
    composition_factor=z$factor,treated_total=z$treated_total,reference_total=z$reference_total,
    all_gene_nonzero_fraction=mean(pmax(new_ref_cp[i,],new_treat_cp[i,])>0),
    shared_gene_nonzero_fraction=mean(pmax(new_ref_cp[i,match(a$universe,colnames(effects))],new_treat_cp[i,match(a$universe,colnames(effects))])>0))
}
write_table(rbindlist(new_qc),'additional_line_quality.tsv')
model<-fread(file.path(work,'outputs/dataset_inventory/cache/depmap/Model.csv'))
new_cells<-unique(new_meta[,.(depmap_id,cell_line,cell_name)])
new_cells<-merge(new_cells,model[,.(depmap_id=ModelID,rrid=RRID,tissue=OncotreeLineage,disease=OncotreePrimaryDisease)],by='depmap_id')
check(all(new_cells$rrid==new_cells$cell_line) && all(new_cells$tissue=='Skin'),'Additional melanoma identity mismatch')
# Genotype is independently documented in the cached original Tahoe metadata.
full_cells<-as.data.table(nanoparquet::read_parquet(file.path(work,'outputs/dataset_inventory/cache/tahoe/cell_line_metadata.parquet')))
v600_ids<-unique(full_cells[Driver_Gene_Symbol=='BRAF'&grepl('V600',Driver_ProtEffect_or_CdnaEffect),Cell_ID_DepMap])
check(all(new_cells$depmap_id %chin% v600_ids),'Additional BRAF-V600 annotation absent')
new_cells[,`:=`(braf_v600_confirmed=TRUE,split='Development')]
cells<-rbindlist(list(cells,new_cells),fill=TRUE)
cells[,`:=`(analysis_included=!depmap_id %chin% EXCLUDED_IDS,
  melanoma_v600=braf_v600_confirmed & tissue=='Skin')]
check(sum(cells$melanoma_v600 & cells$analysis_included)==3L,'Expected three usable V600 melanoma lines')
new_meta[,`:=`(compound=c('Vemurafenib'='vemurafenib','Lapatinib ditosylate'='lapatinib','Cobimetinib'='cobimetinib')[drug],
  perturbation_id=paste0(drug,'__',dose,'uM'))]
new_meta[,observation_id:=paste(depmap_id,perturbation_id,sep='__')]
for(n in c('new_effect','new_ref_cp','new_treat_cp')) {
  x<-get(n);rownames(x)<-new_meta$observation_id;colnames(x)<-colnames(effects);assign(n,x)
}
effects<-rbind(effects,new_effect);reference_cp10k<-rbind(reference_cp10k,new_ref_cp);treated_cp10k<-rbind(treated_cp10k,new_treat_cp)
meta<-rbindlist(list(meta[,.(depmap_id,compound,drug,dose,dose_unit,perturbation_id,observation_id,n_plates)],
  new_meta[,.(depmap_id,compound,drug,dose,dose_unit,perturbation_id,observation_id,n_plates)]))
meta<-merge(meta,cells[,.(depmap_id,cell_name,tissue,split,braf_v600_confirmed,melanoma_v600,analysis_included)],by='depmap_id')
setorder(meta,compound,depmap_id,dose);effects<-effects[meta$observation_id,,drop=FALSE]
reference_cp10k<-reference_cp10k[meta$observation_id,,drop=FALSE];treated_cp10k<-treated_cp10k[meta$observation_id,,drop=FALSE]
# Use primary PRISM consistently for everyone; do not mix it with secondary AUC.
prism<-fread(required[4]);check(!anyDuplicated(prism[,.(depmap_id,compound)]),'Duplicate primary PRISM response')
prism[passed_str_profiling==FALSE|is.na(passed_str_profiling),`:=`(prism_survival=NA_real_,prism_sensitivity=NA_real_)]
meta<-merge(meta,prism[,.(depmap_id,compound,prism_survival,prism_sensitivity,functional_dose_uM,passed_str_profiling)],by=c('depmap_id','compound'),all.x=TRUE)
setorder(meta,compound,depmap_id,dose)
# A matched reference can carry a common perturbation response. Estimate it
# independently in each cell/plate from the other compounds, never these three
# compounds or any PRISM outcome. Retain raw displacement for the full audit.
background_meta<-fread(file.path(art,'plate_common_response_metadata.tsv'))
background<-read_npy_rows(file.path(art,'plate_common_response.npy'),background_meta$array_row)
colnames(background)<-colnames(effects)
plate_map<-rbindlist(list(
  fread(file.path(work,'cache/tahoe_prism_full_chain/tahoe_pdex_selected/pdex_plate_audit.tsv.gz'))[,.(depmap_id,drug,dose,plate)],
  fread(file.path(art,'additional_plate_audit.tsv'))[,.(depmap_id,drug,dose,plate)]))
raw_effects<-effects
background_audit<-list()
for(i in which(meta$analysis_included)) {
  pp<-plate_map[depmap_id==meta$depmap_id[i] & drug==meta$drug[i] & dose==meta$dose[i],unique(plate)]
  jj<-which(background_meta$depmap_id==meta$depmap_id[i] & background_meta$plate %in% pp)
  check(length(jj)==length(pp) & length(jj)>0,'Missing matched plate calibration')
  effects[i,]<-raw_effects[i,]-colMeans(background[jj,,drop=FALSE])
  background_audit[[length(background_audit)+1L]]<-data.table(observation_id=meta$observation_id[i],depmap_id=meta$depmap_id[i],
    compound=meta$compound[i],dose=meta$dose[i],plates=paste(pp,collapse=','),background_profiles=sum(background_meta$n_profiles[jj]))
}
write_table(rbindlist(background_audit),'plate_calibration_applied.tsv')
write_table(background_meta,'plate_calibration_coverage.tsv')
meta[,dose_weight:=c(.25,.5,.25)[match(dose,c(.05,.5,5))]]
line_drug<-unique(meta[,.(depmap_id,compound,cell_name,tissue,split,braf_v600_confirmed,melanoma_v600,analysis_included,prism_survival,prism_sensitivity,functional_dose_uM)])
line_drug[,`:=`(key=paste(depmap_id,compound,sep='__'),eligible=is.finite(prism_survival))];setorder(line_drug,compound,depmap_id)
auc<-matrix(NA_real_,nrow(line_drug),ncol(effects),dimnames=list(line_drug$key,colnames(effects)))
for(i in seq_len(nrow(line_drug))) {
  ii<-which(meta$depmap_id==line_drug$depmap_id[i]&meta$compound==line_drug$compound[i]);check(length(ii)==3L,'Incomplete molecular dose series')
  auc[i,]<-colSums(effects[ii,,drop=FALSE]*meta$dose_weight[ii])
}
write_table(cells,'experiment_cell_lines.tsv');write_table(meta,'tahoe_primary_prism_observations.tsv')
write_table(line_drug,'dose_integrated_molecular_metadata.tsv')
write_table(line_drug[,.(lines=.N,functional_responses=sum(eligible)),by=.(compound,split,analysis_included)],'input_coverage.tsv')
write_table(data.table(path=c(file.path(input_art,'discovery_measurements.rds'),required),
  role=c('Original processed measurements; no previous gene/program/model inputs','Processed treated means for two melanomas','Matched processed reference means','Exact molecular profile metadata','One primary-screen outcome per line and drug at 2.5 uM')),'input_manifest.tsv')
write_table(rbind(fread(file.path(out,'input_manifest.tsv')),data.table(path=file.path(art,c('plate_common_response.npy','plate_common_response_metadata.tsv','plate_background_observations.tsv')),role=c('Per-cell/plate common response across all genes','Independent calibration coverage for retained lines','Non-target compounds that estimated each common response'))),'input_manifest.tsv')
rm(a,new_treated,new_reference,new_effect,new_ref_cp,new_treat_cp);invisible(gc())

## 2. Study-balanced ligand evidence and measurable gene universe ------------
cytometa[,study:=sub('.*@((?:GSE[0-9]+|E-MTAB-[0-9]+)).*','\\1',profile_id,perl=TRUE)]
cytometa[,study_ligand:=paste(study,treatment,sep=' | ')]
study_meta<-unique(cytometa[,.(study_ligand,study,treatment)])
cyto_full<-as.matrix(cyto[,-1,with=FALSE]);rownames(cyto_full)<-cyto$gene
study_effect<-sapply(study_meta$study_ligand,function(k)matrixStats::rowMedians(cyto_full[,cytometa[study_ligand==k,profile_id],drop=FALSE],na.rm=TRUE))
rownames(study_effect)<-cyto$gene
universe<-intersect(colnames(effects),rownames(study_effect)[rowSums(is.finite(study_effect))>=ceiling(ncol(study_effect)/2)])
write_table(cytometa,'cytosig_profile_audit.tsv')
write_table(cytometa[,.(profiles=.N,conditions=paste(unique(condition),collapse='; '),durations=paste(unique(duration),collapse='; ')),by=.(study,treatment)],'cytosig_study_context.tsv')

## 3. Perturbation contrasts select target engagement, receptor and root -----
## BRAF transcript suppression must be allowed in viable bypassing cells.
## Its coordinate is the median molecular AUC in the THREE V600 melanomas.
## Lapatinib/cobimetinib use response-weighted displacement over the retained
## primary-PRISM development lines, rather than ranking noisy gene correlations.
## This is hypothesis-guided marker selection; the discovery labels are NOT
## independent confirmation of the Y. Prediction comes from external tests.
drug_map<-data.table(compound=c('vemurafenib','lapatinib','cobimetinib'))
coordinate<-matrix(NA_real_,length(universe),3,dimnames=list(universe,drug_map$compound))
minimum_effect<-coordinate;measured<-coordinate;selection_detail<-list()
for(drug in drug_map$compound) {
  ix<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included &
    (if(drug=='vemurafenib')line_drug$melanoma_v600 else line_drug$eligible))
  xx<-auc[ix,universe,drop=FALSE]
  if(drug=='vemurafenib') {
    estimate<-matrixStats::colMedians(xx)
    loo<-sapply(seq_along(ix),function(j)matrixStats::colMedians(xx[-j,,drop=FALSE]))
  } else {
    yy<-line_drug$prism_sensitivity[ix]
    estimate<-colSums(xx*yy)/sum(abs(yy))
    loo<-sapply(seq_along(ix),function(j)colSums(xx[-j,,drop=FALSE]*yy[-j])/sum(abs(yy[-j])))
  }
  ii<-which(meta$compound==drug & meta$depmap_id %chin% line_drug$depmap_id[ix])
  ab<-pmax(reference_cp10k[ii,universe,drop=FALSE],treated_cp10k[ii,universe,drop=FALSE])
  coordinate[,drug]<-estimate;minimum_effect[,drug]<-apply(loo*sign(estimate),1,min)
  measured[,drug]<-colMeans(ab>=.1)
  selection_detail[[drug]]<-data.table(compound=drug,gene=universe,n_lines=length(ix),coordinate=estimate,
    minimum_leave_line_signed_effect=minimum_effect[,drug],expressed_fraction=measured[,drug],
    coordinate_type=if(drug=='vemurafenib')'Median V600-melanoma target-engagement displacement' else 'Primary-PRISM-response-weighted displacement')
}
egf<-study_effect[universe,study_meta$treatment=='EGF',drop=FALSE]
egf_est<-matrixStats::rowMedians(egf,na.rm=TRUE)
egf_loo<-sapply(seq_len(ncol(egf)),function(j)matrixStats::rowMedians(egf[,-j,drop=FALSE],na.rm=TRUE))
egf_min<-apply(egf_loo*sign(egf_est),1,min)
fp<-cbind(coordinate,EGF=egf_est);coordinate_sd<-apply(fp,2,sd,na.rm=TRUE)
zfp<-sweep(fp,2,coordinate_sd,'/') # retain biological zero, no centering
v<-zfp[,'vemurafenib'];l<-zfp[,'lapatinib'];c<-zfp[,'cobimetinib'];e<-zfp[,'EGF']
matched_coordinate<-sapply(drug_map$compound,function(drug){
  ix<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included & line_drug$melanoma_v600)
  matrixStats::colMedians(auc[ix,universe,drop=FALSE])
})
rownames(matched_coordinate)<-universe
# B specificity compares the SAME melanoma contexts, with one common drug-effect
# scale. Different cohort composition must not masquerade as a branch difference.
matched_z<-matched_coordinate/coordinate_sd['vemurafenib']
bs<-abs(matched_z[,'vemurafenib'])-pmax(abs(matched_z[,'lapatinib']),abs(matched_z[,'cobimetinib']),abs(e))
write_table(data.table(gene=universe,matched_coordinate),'matched_melanoma_coordinates.tsv.gz')

es<-pmin(abs(e),abs(l))-pmax(abs(v),abs(c));es[e*l>=0]<-NA_real_
rs<-pmin(abs(v),abs(c),abs(e));rs[e*v>=0|e*c>=0]<-NA_real_
# Stability is evaluated at the independent cell-line / ligand-study unit.
# It is not a nominal significance claim from three melanoma lines.
b_ok<-measured[,'vemurafenib']>=.5 & minimum_effect[,'vemurafenib']>0
e_ok<-measured[,'lapatinib']>=.5 & minimum_effect[,'lapatinib']>0 & egf_min>0
r_ok<-measured[,'vemurafenib']>=.5 & measured[,'cobimetinib']>=.5 & minimum_effect[,'vemurafenib']>0 & minimum_effect[,'cobimetinib']>0 & egf_min>0
admission_scores<-data.table(gene=universe,B=bs,E=es,R=rs,B_eligible=b_ok,E_eligible=e_ok,R_eligible=r_ok)
# If a gene supports several contrasts, assign its strongest eligible evidence
# before ranking. This deterministic discovery-only rule prevents double use.
eligible_scores<-as.matrix(admission_scores[,.(B,E,R)])
for(j in 1:3)eligible_scores[!admission_scores[[paste0(c('B','E','R')[j],'_eligible')]]|!is.finite(eligible_scores[,j]),j]<--Inf
admission_scores[,preferred_role:=c('B','E','R')[max.col(eligible_scores,ties.method='first')]]
admissions<-rbindlist(lapply(c('B','E','R'),function(role) {
  z<-admission_scores[preferred_role==role & get(paste0(role,'_eligible'))==TRUE & is.finite(get(role)) & get(role)>0]
  setorderv(z,c(role,'gene'),c(-1,1));check(nrow(z)>=5L,paste('Too few measured/stable candidates for',role,'; inspect coordinate evidence'))
  z<-head(z,N_PER_ROLE);data.table(gene=z$gene,admission_role=role,selection_score=z[[role]],rank=seq_len(nrow(z)))
}))
chosen<-sort(unique(admissions$gene));tests<-rbindlist(selection_detail);selected<-merge(tests,admissions,by='gene',allow.cartesian=TRUE)
write_table(tests,'gene_perturbation_evidence.tsv.gz');write_table(admission_scores,'gene_role_selection_scores.tsv.gz')
write_table(admissions,'selected_gene_admissions.tsv');write_table(data.table(gene=universe,fp),'raw_perturbation_coordinates.tsv.gz')
write_table(data.table(coordinate=names(coordinate_sd),scale=coordinate_sd),'coordinate_scaling.tsv')
write_table(data.table(gene=chosen,zfp[chosen,,drop=FALSE]),'selected_perturbation_coordinates.tsv')

## 4. Six subprograms from real sample correlations, not from templates ------
# Drug-balanced correlations prevent the BRAF response in three melanomas
# from being swamped by unrelated tumor contexts. All three doses are used.
source_rows<-lapply(drug_map$compound,function(drug)which(meta$compound==drug & meta$split=='Development' & meta$analysis_included &
  (drug!='vemurafenib'|meta$melanoma_v600)))
names(source_rows)<-drug_map$compound
gene_sign<-setNames(ifelse(admissions$admission_role=='B',-sign(fp[admissions$gene,'vemurafenib']),sign(fp[admissions$gene,'EGF'])),admissions$gene)
check(all(abs(gene_sign)==1),'Undefined discovery gene orientation')
by_drug_cor<-lapply(source_rows,function(ii)cor(effects[ii,chosen,drop=FALSE],method='spearman'))
ct<-Reduce('+',by_drug_cor)/length(by_drug_cor)
cyto_x<-study_effect[chosen,,drop=FALSE];cc<-cor(t(cyto_x),method='spearman',use='pairwise.complete.obs')
raw_correlation<-(ct+cc)/2
# Opposite-direction genes describe the same activity when their discovery signs
# are reversed. Keep raw correlations as well as activity-oriented structure.
orientation_matrix<-outer(gene_sign[chosen],gene_sign[chosen])
correlation<-raw_correlation*orientation_matrix;check(all(is.finite(correlation)),'Non-estimable source correlation')
distance<-as.dist(pmax(1-correlation,0));tree<-hclust(distance,'average')
# A global hclust height cut repeatedly produced a huge block plus singleton
# branches. PAM partitions around representative genes and avoids that failure;
# hclust still supplies the deterministic within-block display order and tree.
check(!anyDuplicated(admissions$gene),'Overlapping admission pools require an explicit multi-membership design')
# Preserve the hypothesis-defined family contrasts. Unrestricted clustering
# merged receptor-specific markers into the larger shared response. Split each
# family into two empirical subprograms; do not infer the three-family topology
# from the very contrasts used to select genes.
global_partition<-cluster::pam(distance,k=K,diss=TRUE)
family<-setNames(admissions$admission_role,admissions$gene)[chosen]
write_table(data.table(gene=chosen,admission_family=family,unrestricted_program=global_partition$clustering),
  'unrestricted_clustering_family_mixing.tsv')
groups<-setNames(integer(length(chosen)),chosen);partition<-list()
for(i in seq_along(c('B','E','R'))) {
  r<-c('B','E','R')[i];gg<-chosen[family==r]
  part<-cluster::pam(as.dist(1-correlation[gg,gg]),k=2,diss=TRUE)
  partition[[r]]<-part;groups[gg]<-part$clustering+2L*(i-1L)
}
leaf0<-chosen[tree$order];programs<-paste0('P',seq_len(K));k<-K
leaf_order<-unlist(lapply(seq_len(K),function(j)leaf0[groups[leaf0]==j]),use.names=FALSE)
membership<-data.table(gene=leaf_order,program=paste0('P',groups[leaf_order]),plotting_order=seq_along(leaf_order))
colors<-setNames(grDevices::hcl.colors(K,'Dark 3'),programs);row_blocks<-factor(membership$program,levels=programs)
curve<-rbindlist(lapply(1:3,function(q) {
  gg<-setNames(integer(length(chosen)),chosen)
  for(i in seq_along(c('B','E','R'))) {
    ix<-chosen[family==c('B','E','R')[i]]
    gg[ix]<-(if(q==1L)rep(1L,length(ix)) else pam(as.dist(1-correlation[ix,ix]),q,diss=TRUE)$clustering)+q*(i-1L)
  }
  si<-silhouette(gg,distance)
  data.table(k=3*q,mean_silhouette=mean(si[,'sil_width']),minimum_size=min(table(gg)),negative_fraction=mean(si[,'sil_width']<0))
}))
write_table(curve,'silhouette_by_k.tsv')
write_table(data.table(gene=chosen,program=paste0('P',groups),silhouette=silhouette(groups,distance)[,'sil_width']),'gene_silhouettes.tsv')
write_table(data.table(method=c('hclust height cut','Unrestricted PAM','Two subprograms per candidate family'),
  minimum_size=c(min(table(cutree(tree,K))),min(table(global_partition$clustering)),min(table(groups))),
  largest_size=c(max(table(cutree(tree,K))),max(table(global_partition$clustering)),max(table(groups)))), 'partition_failure_repair.tsv')
coherence<-rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];rbindlist(lapply(c('Tahoe','CytoSig'),function(src){m<-if(src=='Tahoe')ct[g,g] else cc[g,g];vv<-m[upper.tri(m)];data.table(program=p,source=src,n_genes=length(g),mean_rho=mean(vv),median_rho=median(vv),positive_fraction=mean(vv>0))}))}))
write_table(coherence,'program_source_coherence.tsv')
context_cor<-c(by_drug_cor,list(CytoSig=cc,Tahoe_balanced=ct))
context_coherence<-rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];rbindlist(lapply(names(context_cor),function(src){
  m<-context_cor[[src]][g,g];sg<-m*outer(gene_sign[g],gene_sign[g]);
  data.table(program=p,source=src,n_genes=length(g),raw_mean_rho=mean(m[upper.tri(m)]),signed_mean_rho=mean(sg[upper.tri(sg)]))
}))}))
write_table(context_coherence,'program_context_coherence.tsv')
write_table(data.table(gene=chosen,correlation),'joint_gene_correlations.tsv.gz')
write_table(data.table(gene=chosen,tcrossprod(is.finite(cyto_x)*1L)),'cytosig_pairwise_study_counts.tsv.gz')

## 5. Candidate hierarchy: classify observed subprogram centroids ------------
signed_fp<-zfp[chosen,,drop=FALSE]*gene_sign[chosen]
interpretation_fp<-copy(zfp)
interpretation_fp[admissions[admission_role=='B',gene],1:3]<-matched_z[admissions[admission_role=='B',gene],]
signed_fp<-interpretation_fp[chosen,,drop=FALSE]*gene_sign[chosen]
centroid<-t(sapply(programs,function(p)colMeans(signed_fp[membership[program==p,gene],,drop=FALSE])))
v<-centroid[,'vemurafenib'];l<-centroid[,'lapatinib'];c<-centroid[,'cobimetinib'];e<-centroid[,'EGF']
b_specificity<-abs(v)-pmax(abs(l),abs(c),abs(e));e_specificity<-pmin(abs(e),abs(l))-pmax(abs(v),abs(c))
observed_roles<-ifelse(b_specificity>0,'B',ifelse(e_specificity>0 & e*l<0,'E',ifelse(e*v<0 & e*c<0,'R','Unresolved')))
roles<-rep(c('B','E','R'),each=2L)
orient<-gene_sign
membership[,`:=`(role=roles[match(program,programs)],orientation=gene_sign[gene])]
membership[,orientation_basis:=ifelse(role=='B','Opposite discovery vemurafenib effect','Discovery EGF effect')]
membership[,weight:=orientation/.N,by=program]
program_order<-programs[order(match(roles,c('B','E','R','Unresolved')),seq_along(programs))]
leaf_order<-unlist(lapply(program_order,function(p)membership[program==p,gene]),use.names=FALSE)
membership<-membership[match(leaf_order,gene)];membership[,plotting_order:=seq_len(.N)]
row_blocks<-factor(membership$program,levels=program_order)
role_audit<-data.table(program=programs,role=roles,observed_centroid_role=observed_roles,centroid_matches_family=roles==observed_roles,n_genes=as.integer(table(factor(membership$program,levels=programs))),centroid,
  braf_specificity=b_specificity,erbb_specificity=e_specificity)
write_table(role_audit,'program_hierarchy_evidence.tsv');write_table(membership,'frozen_program_genes.tsv')
hierarchy_resolved<-all(roles==observed_roles)
# Render the Y as a hypothesis-guided candidate hierarchy. The edges are NOT
# inferred causal gene edges, and centroid matching is not validation evidence.
node<-data.table(role=c('B','E','R'),x=c(0,2,1),y=c(1,1,0),
  title=c('BRAFi-specific response','ERBB-associated response','Shared downstream response'))
node[,label:=paste0(title,'\n',vapply(role,function(r)paste(role_audit[role==r,program],collapse=', '),character(1)))]
p<-ggplot()+geom_segment(data=data.table(x=c(.24,1.76),y=.84,xend=c(.86,1.14),yend=.18),aes(x=x,y=y,xend=xend,yend=yend),linewidth=.7,linetype=2,
  arrow=grid::arrow(length=grid::unit(2.5,'mm')))+geom_label(data=node,aes(x,y,label=label,fill=role),size=4,label.padding=grid::unit(.35,'lines'))+
  scale_fill_manual(values=c(B='#E3B88B',E='#A0D6D7',R='#C4AFD8'))+theme_void()+theme(legend.position='none',plot.margin=margin(10,12,10,12))+
  scale_x_continuous(limits=c(-.8,2.8))+scale_y_continuous(limits=c(-.4,1.4))+
  labs(title='Candidate program hierarchy | perturbation-informed selection',subtitle='Family topology is hypothesis-defined; subprograms are data-derived; validation is separate')
plot_save(p,'27h_candidate_hierarchy',10,5.5)

## 6. Freeze before evaluating primary-PRISM holdout or external outcomes -----
score<-function(x)sapply(programs,function(p){g<-membership[program==p,gene];colMeans(x[g,,drop=FALSE]*gene_sign[g])})
hierarchy_score<-function(x,partial=FALSE)sapply(c('B','E','R'),function(r){g<-membership[role==r,gene];if(!length(g))return(rep(NA_real_,ncol(x)))
  present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)]
  if(!length(present)||(!partial&&length(present)!=length(g)))return(rep(NA_real_,ncol(x)))
  colMeans(x[present,,drop=FALSE]*membership$orientation[match(present,membership$gene)])})
freeze<-list(membership=as.data.frame(membership),admissions=admissions,orientation=orient,roles=roles,partition=partition,
  coordinate_scales=coordinate_sd,selection='Target-engagement/specificity contrasts; NOT discovery proof of the architecture',
  validation='Primary PRISM at 2.5 uM; frozen B/E/R and six subprogram scores; no selection using external outcomes')
saveRDS(freeze,file.path(art,'frozen_programs.rds'))
write_table(data.table(prediction=c('BRAF inhibition suppresses B','Bypass retains E and R while B remains suppressed',
  'Adding lapatinib lowers the receptor-supported root in susceptible contexts','EGF activates receptor-associated and shared-response programs',
  'Receptor-driven ERK may resist RAF/MEK inhibition: estimate cobimetinib effects without assuming a 50% root-off threshold'),
  basis=c(rep('Candidate hierarchy; test in held-out expression',4),'Gerosa et al. report differential resistance of receptor-driven MAPK signaling')),'external_prediction_freeze.tsv')
## 7. Render the actual structure and the repaired selection mechanism -------
cor_limit<-as.numeric(quantile(abs(c(ct[upper.tri(ct)],cc[upper.tri(cc)])),.95))
write_table(data.table(limit=cor_limit,tahoe_clipped=mean(abs(ct[upper.tri(ct)])>cor_limit),cytosig_clipped=mean(abs(cc[upper.tri(cc)])>cor_limit)),'correlation_display_scale.tsv')
for(src in c('Tahoe','CytoSig','Joint')) {
  m<-switch(src,Tahoe=ct,CytoSig=cc,Joint=raw_correlation)[leaf_order,leaf_order];diag(m)<-NA_real_
  heat(m,paste0('27f_',tolower(src),'_gene_correlations'),paste(src,'| six subprograms | fixed gene order'),row_blocks,row_blocks,cor_limit,13)
}
for(src in c('Tahoe','CytoSig','Joint')) {
  m<-(switch(src,Tahoe=ct,CytoSig=cc,Joint=raw_correlation)*orientation_matrix)[leaf_order,leaf_order];diag(m)<-NA_real_
  heat(m,paste0('27f_',tolower(src),'_signed_correlations'),paste(src,'| activity-oriented correlations | frozen gene signs'),row_blocks,row_blocks,cor_limit,13)
}
heat(zfp[leaf_order,,drop=FALSE],'27g_perturbation_fingerprint','Selection/interpretation coordinates | not the clustering input',row_blocks,width=8)
tahoe_x<-t(effects[unlist(source_rows),chosen,drop=FALSE])
source_heat(tahoe_x,'27g_tahoe_gene_samples','Tahoe | three V600 melanomas for BRAFi | all doses | fixed genes')
source_heat(cyto_x,'27g_cytosig_gene_samples','CytoSig | EGF/NRG1 studies | same genes and order')
p<-ggplot(curve,aes(k,mean_silhouette))+geom_line()+geom_point(size=2)+geom_vline(xintercept=K,lty=2,color='#B64039')+
  theme_classic(base_size=12)+labs(x='Number of subprograms',y='Mean silhouette',title='PAM resolution diagnostic | six subprograms retained')
plot_save(p,'27e_rank_silhouette',6,4.5)
for(ext in c('pdf','png')) {
  if(ext=='pdf')cairo_pdf(file.path(out,'27h_gene_hierarchy.pdf'),width=13,height=7) else png(file.path(out,'27h_gene_hierarchy.png'),width=13,height=7,units='in',res=140,type='cairo')
  par(mar=c(8,5,3,1));plot(tree,labels=paste0(chosen,' [',paste0('P',groups),']'),cex=.6,hang=-1,
    main='Gene hierarchy | subprogram labels from PAM',xlab='',sub='Average linkage on joint correlation distance; within-block display order derives from this tree')
  invisible(dev.off())
}
# A separate empirical program-similarity tree uses real-sample gene distances.
# It is not an upstream/downstream inference; the Y states the directional hypothesis.
program_distance<-outer(programs,programs,Vectorize(function(p,q) {
  if(p==q)return(0)
  mean(1-correlation[membership[program==p,gene],membership[program==q,gene],drop=FALSE])
}))
dimnames(program_distance)<-list(programs,programs)
program_tree<-hclust(as.dist(program_distance),'average')
write_table(data.table(program=programs,program_distance),'empirical_program_distances.tsv')
for(ext in c('pdf','png')) {
  if(ext=='pdf')cairo_pdf(file.path(out,'27h_empirical_program_tree.pdf'),width=8,height=4.8) else png(file.path(out,'27h_empirical_program_tree.png'),width=8,height=4.8,units='in',res=140,type='cairo')
  par(mar=c(4,3,3,8));display_tree<-program_tree
  display_tree$labels<-paste0(display_tree$labels,' | ',roles[match(display_tree$labels,programs)])
  dd<-as.dendrogram(display_tree)
  plot(dd,horiz=TRUE,main='Empirical subprogram similarity',xlab='Mean gene-pair distance; no causal direction',cex=1.2)
  invisible(dev.off())
}
# A small literature-defined diagnostic panel is traced, never manually added.
markers<-intersect(c('EGFR','ERBB2','ERBB3','DUSP4','DUSP6','SPRY2','SPRY4','ERRFI1','EGR1','FOS','FOSL1','ETV4','ETV5'),universe)
marker_trace<-rbindlist(lapply(drug_map$compound,function(drug)rbindlist(lapply(c(.05,.5,5),function(d){
  ii<-which(meta$compound==drug & meta$dose==d & meta$melanoma_v600 & meta$analysis_included)
  data.table(compound=drug,dose=d,gene=markers,n_melanomas=length(ii),median_effect=matrixStats::colMedians(effects[ii,markers,drop=FALSE]),
    sign_agreement=abs(colMeans(sign(effects[ii,markers,drop=FALSE]))),selected=markers %chin% chosen)
}))))
write_table(marker_trace,'mechanism_marker_trace.tsv')
old_tests<-fread(file.path(root,'results/26_braf_erbb_bypass/gene_selection_tests.tsv.gz'))
marker_loss<-merge(old_tests[gene %chin% markers,.(gene,compound,old_cross_tumor_rho=rho_auc,old_selected=selected)],
  admissions[gene %chin% markers],by='gene',all.x=TRUE,allow.cartesian=TRUE)
write_table(marker_loss,'why_original_selection_lost_pathway_genes.tsv')
p<-ggplot(marker_trace,aes(factor(dose),median_effect,group=compound,color=compound))+geom_hline(yintercept=0,color='grey75')+geom_line()+geom_point(size=1.7)+
  facet_wrap(~gene,ncol=4,scales='free_y')+theme_classic(base_size=10)+labs(x='Dose (uM)',y='Median signed molecular effect',title='Pathway-marker diagnostic | three V600 melanomas',subtitle='Markers are literature-based checks; admission follows the written contrast rule')
plot_save(p,'27c_marker_dose_trace',11,10)

# Decompose the old B signature into raw response and common plate response.
# This archived membership is a diagnostic comparator only; it never admits genes.
old_membership_path<-file.path(art,'archive/program_level_signs/frozen_program_genes.tsv')
if(file.exists(old_membership_path)) {
  old_m<-fread(old_membership_path)[role=='B'];old_g<-old_m$gene
  diag_rows<-which(meta$melanoma_v600 & meta$analysis_included)
  plate_failure<-rbindlist(lapply(diag_rows,function(ii) {
    raw<-mean(raw_effects[ii,old_g]*old_m$orientation)
    adjusted<-mean(effects[ii,old_g]*old_m$orientation)
    data.table(cell_line=meta$cell_name[ii],compound=meta$compound[ii],dose=meta$dose[ii],
      raw_B=raw,common_response_B=raw-adjusted,compound_specific_B=adjusted)
  }))
  write_table(plate_failure,'old_B_plate_confound_decomposition.tsv')
  pp<-melt(plate_failure,id.vars=c('cell_line','compound','dose'),variable.name='component',value.name='score')
  p<-ggplot(pp,aes(factor(dose),score,color=component,group=component))+geom_hline(yintercept=0,color='grey75')+
    geom_line()+geom_point(size=1.5)+facet_grid(cell_line~compound)+theme_classic(base_size=11)+
    labs(x='Dose (uM)',y='Old B signature score',color='Component',title='Why the previous B branch failed',
      subtitle='The matched-plate common response explains much of the apparent BRAF-specific displacement')
  plot_save(p,'27b_plate_confound_diagnosis',12,8)
}
p<-ggplot(context_coherence,aes(source,signed_mean_rho,fill=program))+geom_col(position='dodge')+
  facet_wrap(~program,ncol=3)+theme_classic(base_size=10)+theme(legend.position='none',axis.text.x=element_text(angle=45,hjust=1))+
  labs(x=NULL,y='Mean signed within-program Spearman correlation',title='Source-specific coherence | same frozen genes and signs')
plot_save(p,'27f_context_coherence',10,6)

## 8. Functional transfer: primary PRISM, unchanged seven-line holdout --------
scores<-score(t(auc[,chosen,drop=FALSE]));hscores<-hierarchy_score(t(auc[,chosen,drop=FALSE]))
fit_rows<-which(line_drug$split=='Development'&line_drug$analysis_included&line_drug$eligible)
held<-which(line_drug$split=='Holdout'&line_drug$analysis_included&line_drug$eligible)
train_meta<-copy(line_drug[fit_rows]);base<-lm(prism_survival~compound,data=train_meta)
fit_results<-list();held_points<-list()
for(j in seq_along(programs)) {
  train<-copy(train_meta);train[,program_score:=scores[fit_rows,j]];fit<-lm(prism_survival~compound+program_score,data=train)
  tt<-copy(line_drug[held]);tt[,program_score:=scores[held,j]]
  tt[,`:=`(program=programs[j],predicted_baseline=predict(base,tt),predicted_program=predict(fit,tt))]
  tt[,residual_survival:=prism_survival-predicted_baseline]
  fit_results[[j]]<-data.table(program=programs[j],role=roles[j],n_lines=uniqueN(tt$depmap_id),n_observations=nrow(tt),
    baseline_rmse=sqrt(mean((tt$prism_survival-tt$predicted_baseline)^2)),program_rmse=sqrt(mean((tt$prism_survival-tt$predicted_program)^2)),
    pooled_residual_spearman=cor(tt$program_score,tt$residual_survival,method='spearman'))
  held_points[[j]]<-tt
}
fit_results<-rbindlist(fit_results);held_points<-rbindlist(held_points)
write_table(fit_results,'subprogram_heldout_model_comparison.tsv');write_table(held_points,'subprogram_heldout_predictions.tsv')
train<-copy(train_meta);train[,shared_root_score:=hscores[fit_rows,'R']]
root_fit<-lm(prism_survival~compound+shared_root_score,data=train)
root_test<-copy(line_drug[held]);root_test[,shared_root_score:=hscores[held,'R']]
root_test[,`:=`(predicted_baseline=predict(base,root_test),predicted_root=predict(root_fit,root_test))]
root_test[,residual_survival:=prism_survival-predicted_baseline]
root_comparison<-data.table(n_lines=uniqueN(root_test$depmap_id),n_observations=nrow(root_test),
  baseline_rmse=sqrt(mean((root_test$prism_survival-root_test$predicted_baseline)^2)),
  root_rmse=sqrt(mean((root_test$prism_survival-root_test$predicted_root)^2)),
  residual_spearman=cor(root_test$shared_root_score,root_test$residual_survival,method='spearman'))
rank_x<-resid(lm(rank(shared_root_score)~compound,data=root_test))
rank_y<-resid(lm(rank(prism_survival)~compound,data=root_test))
root_comparison[,partial_spearman_controlling_drug:=cor(rank_x,rank_y)]
root_test[,`:=`(within_drug_score=shared_root_score-mean(shared_root_score),
  within_drug_survival=prism_survival-mean(prism_survival)),by=compound]
write_table(root_comparison,'root_heldout_model_comparison.tsv');write_table(root_test,'root_heldout_predictions.tsv')
write_table(root_test[,.(n_lines=.N,rho=cor(shared_root_score,prism_survival,method='spearman')),by=compound],'root_within_drug_holdout_correlations.tsv')
p<-ggplot(held_points,aes(program_score,residual_survival,color=compound))+geom_point(size=1.8)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+facet_wrap(~program,scales='free')+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+labs(x='Frozen subprogram score',y='Primary-PRISM survival minus training drug mean',title='All six subprograms | held-out functional association')
plot_save(p,'27i_subprogram_prism_associations',11,8)
p<-ggplot(root_test,aes(within_drug_score,within_drug_survival,color=compound))+geom_point(size=2)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+theme_classic(base_size=12)+theme(aspect.ratio=1)+
  labs(x='Shared-root score, centered within drug',y='Primary-PRISM survival, centered within drug',title='Shared root | six held-out lines',subtitle='Both variables adjusted for drug; repeated lines')
plot_save(p,'27i_shared_root_prism',7,5.5)
# Save compact development state for review; it does not select on holdout outcomes.
saveRDS(list(freeze=freeze,source_correlations=list(Tahoe=ct,CytoSig=cc),cohort=cells,metadata=line_drug),file.path(art,'discovery_objects.rds'))

## 9. External metadata: include all 96 GSE127988 wells in identity audit -----
s684<-fread(file.path(root,'results/26_braf_erbb_bypass/gse68453_sample_audit.tsv'))
factorial_lines<-c('A375','DM331','HT144','SKMEL24','SKMEL28')
external_lines<-c('A375','SKMEL24','HT144','DM331','SKMEL28','DM13','VMM5A','VMM18','DM6','VMM15','VMM12','VMM17')
s127_all<-get_soft_samples(readLines(gzfile(file.path(input_art,'gse127988/GSE127988_family.soft.gz')),warn=FALSE))
number<-function(x,label)as.numeric(sub(paste0('.*',label,': ([0-9.]+).*'),'\\1',x))
s127_all[,`:=`(vem_uM=number(characteristics,'vemurafenib for 24 hrs \\(um\\)'),cobi_uM=number(characteristics,'cobimetinib for 24 hrs \\(um\\)'),
  egf_h=number(characteristics,'egf at 100 ng/ml \\(h\\)'),biological_set=number(characteristics,'set id'),sample_id=as.integer(sub('.*\\(S([0-9]+)\\).*','\\1',title)))]
s127_all[,`:=`(cell_line='A375',arm=ifelse(vem_uM==0,'DMSO',ifelse(cobi_uM==0,'Vem','Vem + Cobi')))]
s127<-s127_all[vem_uM==1 & biological_set %in% 3:4]
check(nrow(s127_all)==96L & nrow(s127)==56L,'GSE127988 design mismatch')
write_table(s127_all,'gse127988_complete_sample_audit.tsv')
write_table(s127_all[,.(wells=.N),by=.(biological_set,arm,vem_uM,cobi_uM,egf_h)],'gse127988_complete_design.tsv')
## 10. External gene identity and platform coverage ---------------------------
## Gene identifiers are repaired using approved/unique historical HGNC symbols.
## External values cannot change membership, orientation or gene weights.
hgnc_path<-file.path(work,'outputs/dataset_inventory/cache/cellfie/hgnc_complete_set.txt')
hgnc<-fread(hgnc_path);hgnc<-hgnc[status=='Approved']
lookup<-rbindlist(lapply(seq_len(nrow(hgnc)),function(i) {
  a<-unique(c(hgnc$symbol[i],unlist(strsplit(paste(hgnc$prev_symbol[i],hgnc$alias_symbol[i],sep='|'),'\\|'))))
  data.table(input=a[nzchar(a)&!is.na(a)&a!='NA'],symbol=hgnc$symbol[i])
}))
lookup<-lookup[,if(uniqueN(symbol)==1L).(symbol=symbol[1]) else .(symbol=NA_character_),by=input]
canonical<-function(v){ans<-lookup$symbol[match(v,lookup$input)];exact<-v %chin% hgnc$symbol;ans[exact]<-v[exact];ans}
frozen_symbols<-canonical(leaf_order)
check(!anyNA(frozen_symbols)&&!anyDuplicated(frozen_symbols),'Frozen gene aliases are missing or ambiguous')
write_table(data.table(gene=leaf_order,official_symbol=frozen_symbols,hgnc_name=hgnc$name[match(frozen_symbols,hgnc$symbol)],
  locus_type=hgnc$locus_type[match(frozen_symbols,hgnc$symbol)]),'selected_gene_identities.tsv')
# Read the platform annotation embedded in the downloaded SOFT file.
platform_path<-file.path(input_art,'GSE68453_platform_annotation.tsv')
if(!file.exists(platform_path)) {
  con<-gzfile(file.path(input_art,'GSE68453_family.soft.gz'),'rt');buffer<-character();reading<-FALSE;done<-FALSE
  repeat {
    lines<-readLines(con,n=10000,warn=FALSE);if(!length(lines))break
    for(l in lines) {
      if(startsWith(l,'!platform_table_end')) {done<-TRUE;break}
      if(reading)buffer<-c(buffer,l)
      if(startsWith(l,'!platform_table_begin'))reading<-TRUE
    }
    if(done)break
  };close(con);check(done,'Platform annotation not found in SOFT');writeLines(buffer,platform_path)
}
platform<-fread(platform_path,quote='');platform[,official_symbol:=canonical(Symbol)]
matfile<-file.path(input_art,'GSE68453_series_matrix.txt.gz')
m684<-suppressWarnings(fread(matfile,skip='"ID_REF"',header=TRUE,check.names=FALSE,select=c('ID_REF',s684$gsm)))
m684<-m684[ID_REF %chin% platform$ID]
x684<-as.matrix(m684[,s684$gsm,with=FALSE]);storage.mode(x684)<-'double'
check(all(is.finite(x684)),'Non-finite GSE68453 measurements')
# GEO supplies normalized log2 array expression. Assert the scale, do not log twice.
check(quantile(x684,.99)<30 && quantile(x684,.5)>0,'Unexpected GSE68453 expression scale')
probe<-data.table(probe=m684$ID_REF,official_symbol=platform$official_symbol[match(m684$ID_REF,platform$ID)],
  iqr=matrixStats::rowIQRs(x684),array_row=seq_len(nrow(x684)))
probe<-probe[!is.na(official_symbol)];setorder(probe,official_symbol,-iqr,probe)
probe[,chosen_probe:=seq_len(.N)==1L,by=official_symbol]
write_table(probe,'gse68453_probe_mapping.tsv.gz')
keep<-probe[chosen_probe==TRUE];x684<-x684[keep$array_row,,drop=FALSE];rownames(x684)<-keep$official_symbol
# Individual RPM columns are selected by sample ID and independently checked
# against GEO drug doses, exposure times and biological day. Summary columns
# are explicitly excluded; they are not additional replicates.
rpm_path<-file.path(input_art,'gse127988/GSE127988_RPM_Gerosa_et_al.csv.gz')
rh<-fread(rpm_path,header=FALSE,nrows=6,check.names=FALSE)
header<-as.matrix(rh[,-c(1,2),with=FALSE])
sample_ids<-as.integer(sub('.*\\(S([0-9]+)\\).*','\\1',s127_all$title))
ci<-vapply(sample_ids,function(id){z<-which(header[4,]=='Individual'&header[5,]==as.character(id));check(length(z)==1L,'Ambiguous RPM sample ID');z},integer(1))
day_mismatch<-which(as.numeric(header[6,ci])!=s127_all$biological_set)
if(length(day_mismatch)) {
  check(setequal(s127_all$sample_id[day_mismatch],c(15L,16L)) &
    all(s127_all$biological_set[day_mismatch]==2L) & all(as.numeric(header[6,ci[day_mismatch]])==1L),
    'Unexpected GEO/author-matrix day discrepancy')
  write_table(data.table(gsm=s127_all$gsm[day_mismatch],sample_id=s127_all$sample_id[day_mismatch],
    GEO_day=s127_all$biological_set[day_mismatch],author_matrix_day=as.numeric(header[6,ci[day_mismatch]]),
    decision='Use author processed-matrix set ID; S15/S16 complete day-1 untreated controls'), 'gse127988_metadata_repairs.tsv')
  s127_all[day_mismatch,biological_set:=as.numeric(header[6,ci[day_mismatch]])]
  write_table(s127_all,'gse127988_complete_sample_audit.tsv')
  write_table(s127_all[,.(wells=.N),by=.(biological_set,arm,vem_uM,cobi_uM,egf_h)],'gse127988_complete_design.tsv')
}
check(all(as.numeric(header[1,ci])==s127_all$vem_uM)&all(as.numeric(header[2,ci])==s127_all$cobi_uM)&
  all(as.numeric(header[3,ci])==s127_all$egf_h)&all(as.numeric(header[6,ci])==s127_all$biological_set),'RPM/GEO design mismatch')
rpm<-fread(rpm_path,skip=7,header=FALSE,select=c(1L,2L,ci+2L),check.names=FALSE)
numeric_rpm<-as.matrix(rpm[,-c(1,2),with=FALSE]);storage.mode(numeric_rpm)<-'double'
check(all(is.finite(numeric_rpm)&numeric_rpm>=0),'Invalid RPM measurements')
rna_symbol<-hgnc$symbol[match(sub('\\..*','',rpm[[1]]),hgnc$ensembl_gene_id)]
fallback<-canonical(rpm[[2]]);rna_symbol[is.na(rna_symbol)]<-fallback[is.na(rna_symbol)]
rna_map<-data.table(ensembl=rpm[[1]],reported_symbol=rpm[[2]],official_symbol=rna_symbol)
write_table(rna_map,'gse127988_gene_mapping.tsv.gz')
# Several Ensembl rows mapping to one gene are aggregated on the RPM scale.
rna<-data.table(official_symbol=rna_symbol,numeric_rpm)
rna<-rna[!is.na(official_symbol),lapply(.SD,sum),by=official_symbol]
x127<-log2(1+as.matrix(rna[,-1,with=FALSE]));rownames(x127)<-rna$official_symbol;colnames(x127)<-s127_all$gsm
coverage<-rbindlist(lapply(c('GSE68453','GSE127988'),function(ds) {
  xx<-if(ds=='GSE68453')x684 else x127
  data.table(dataset=ds,gene=leaf_order,official_symbol=frozen_symbols,program=membership$program,
    measured=frozen_symbols %chin% rownames(xx),reason=ifelse(frozen_symbols %chin% rownames(xx),'Mapped','No uniquely mapped platform/RNA measurement'))
}))
write_table(coverage,'external_gene_coverage.tsv')
hierarchy_coverage<-merge(coverage,membership[,.(gene,role)],by='gene')[,.(n_frozen=.N,n_measured=sum(measured),partial=!all(measured)),by=.(dataset,role)]
write_table(hierarchy_coverage,'external_hierarchy_coverage.tsv')
write_table(data.table(gene=leaf_order,official_symbol=frozen_symbols,ensembl=hgnc$ensembl_gene_id[match(frozen_symbols,hgnc$symbol)],
  raw_symbol_present=frozen_symbols %chin% rpm[[2]],raw_ensembl_present=hgnc$ensembl_gene_id[match(frozen_symbols,hgnc$symbol)] %chin% rpm[[1]]),
  'gse127988_raw_identifier_audit.tsv')
extract_genes<-function(x) {
  m<-matrix(NA_real_,length(leaf_order),ncol(x),dimnames=list(leaf_order,colnames(x)))
  ok<-frozen_symbols %chin% rownames(x);m[ok,]<-x[frozen_symbols[ok],,drop=FALSE];m
}
x684<-extract_genes(x684);x127_full<-extract_genes(x127);x127<-x127_full[,s127$gsm,drop=FALSE]
# Report complete and observed-gene scores separately; no zeros or replacement genes.
external_score<-function(x) {
  # Missingness-only adaptation: keep the frozen membership and orientation,
  # but expose an observed-gene mean separately from the complete score. A
  # missing gene is never zero. Counts/coverage identify every partial score.
  sapply(programs,function(p) {
    g<-membership[program==p,gene]
    present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)]
    if(!length(present))return(rep(NA_real_,ncol(x)))
    colMeans(x[present,,drop=FALSE]*gene_sign[present])
  })
}
write_table(data.table(dataset=coverage$dataset,program=coverage$program,measured=coverage$measured)[,
  .(n_genes=.N,n_measured=sum(measured),complete=all(measured)),by=.(dataset,program)],'external_program_coverage.tsv')

## 11. GSE68453 factorial displacement and independent efficacy association --
# Controls are a within-line reference. Replicate labels are not assumed to be
# matched biological pairs across treatment arms.
d684<-x684
for(l in unique(s684$cell_line)) {
  jj<-which(s684$cell_line==l);ref<-rowMeans(x684[,s684$cell_line==l&s684$arm=='control',drop=FALSE])
  d684[,jj]<-x684[,jj,drop=FALSE]-ref
}
s684scores<-external_score(d684)
hierarchy_values684<-hierarchy_score(d684,partial=TRUE)
hierarchy_long684<-melt(data.table(s684[,.(gsm,cell_line,arm,replicate)],hierarchy_values684),id.vars=c('gsm','cell_line','arm','replicate'),variable.name='role',value.name='score')
hierarchy_means684<-hierarchy_long684[,.(score=mean(score)),by=.(cell_line,arm,role)]
write_table(hierarchy_long684,'gse68453_hierarchy_scores.tsv')
hierarchy_predictors<-dcast(hierarchy_means684[arm=='plx'],cell_line~role,value.var='score')
hierarchy_predictors[,BypassIndex:=pmin(as.numeric(scale(E)),as.numeric(scale(R)))]
write_table(hierarchy_predictors,'gse68453_frozen_bypass_predictors.tsv')

long684<-melt(data.table(s684[,.(gsm,cell_line,arm,replicate)],s684scores),id.vars=c('gsm','cell_line','arm','replicate'),variable.name='program',value.name='score')
long684<-merge(long684,coverage[dataset=='GSE68453',.(n_frozen=.N,n_measured=sum(measured)),by=program],by='program',sort=FALSE)
long684[,complete_score:=ifelse(n_frozen==n_measured,score,NA_real_)]
write_table(long684,'gse68453_program_scores.tsv')
means684<-long684[,.(score=mean(score)),by=.(cell_line,arm,program)]
factorial<-dcast(means684[cell_line %chin% factorial_lines],cell_line+program~arm,value.var='score')
factorial[,`:=`(plx_minus_control=plx-control,lap_minus_control=lapatinib-control,
  combo_minus_control=plx_lap-control,adding_lapatinib=plx_lap-plx,interaction=plx_lap-plx-lapatinib+control)]
write_table(factorial,'gse68453_factorial_contrasts.tsv')
# Frozen genes keep their permanent order in both external representations.
external_heat<-function(x,metadata,prefix,semantic_order) {
  finite_genes<-rowSums(is.finite(x))==ncol(x)
  check(sum(finite_genes)>1,'Too few measured frozen genes for an external heatmap')
  c<-cor(x[finite_genes,,drop=FALSE],method='spearman');check(all(is.finite(c)),'Non-estimable external sample correlations')
  tr<-hclust(as.dist(pmax(1-c,0)),method='average');ord<-tr$order
  blocks<-cutree(tr,min(4,ncol(x)-1));bl<-factor(blocks[ord],levels=unique(blocks[ord]))
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),cluster=blocks[ord]),paste0(prefix,'_sample_order.tsv'))
  sample_limit<-as.numeric(quantile(abs(c[upper.tri(c)]),.95));diag(c)<-NA_real_
  heat(c[ord,ord],paste0(prefix,'_sample_correlations'),paste(prefix,'| sample correlations; 95% off-diagonal scale'),bl,bl,limit=sample_limit,width=12,height=11,show_cols=FALSE)
  heat(x[leaf_order,ord,drop=FALSE],paste0(prefix,'_gene_samples_clustered'),paste(prefix,'| log-expression displacement | frozen genes'),row_blocks,bl,width=13)
  heat(x[leaf_order,semantic_order,drop=FALSE],paste0(prefix,'_gene_samples_semantic'),paste(prefix,'| experimental order | frozen genes'),row_blocks,width=13)
}
external_heat(d684,s684,'27k_gse68453',order(match(s684$cell_line,external_lines),match(s684$arm,c('control','plx','lapatinib','plx_lap')),s684$replicate))
p<-ggplot(means684[cell_line %chin% factorial_lines],aes(factor(arm,levels=c('control','plx','lapatinib','plx_lap')),score,color=program,group=program))+
  geom_hline(yintercept=0,color='grey70')+geom_line()+geom_point(size=2)+facet_wrap(~cell_line,ncol=3)+
  scale_color_manual(values=colors)+theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=30,hjust=1))+
  labs(x=NULL,y='Frozen program displacement',title='GSE68453 | five complete factorials',subtitle='Same genes and weights; empirical program labels')
plot_save(p,'27l_factorial_programs',10,7)
# Write expression-only predictors BEFORE attaching digitized efficacy.
predictors684<-means684[arm=='plx',.(cell_line,program,plx_displacement=score)]
write_table(predictors684,'gse68453_frozen_expression_predictors.tsv')
bliss<-data.table(cell_line=external_lines,bliss=c(5.7,29.3,16.3,27.3,49.7,39.3,43,31.7,-2.7,-8.7,3.7,-1))
benefit<-merge(predictors684,bliss,by='cell_line')
benefit_stats<-benefit[,{
  ok<-is.finite(plx_displacement)&is.finite(bliss)
  if(sum(ok)>=3) {tt<-cor.test(plx_displacement[ok],bliss[ok],method='spearman',exact=FALSE);.(n=sum(ok),rho=unname(tt$estimate),p=tt$p.value)}
  else .(n=sum(ok),rho=NA_real_,p=NA_real_)
},by=program]
benefit_stats[,bh_q:=p.adjust(p,'BH')]
write_table(benefit,'gse68453_expression_bliss_join.tsv');write_table(benefit_stats,'gse68453_exploratory_program_bliss.tsv')
p<-ggplot(benefit,aes(plx_displacement,bliss))+geom_point(color='#287D99',size=2)+geom_text(aes(label=cell_line),size=3,vjust=-.6,check_overlap=TRUE)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='#B64039')+facet_wrap(~program,scales='free_x')+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+labs(x='Frozen PLX program displacement',y='Digitized Bliss benefit',title='Exploratory external efficacy association',subtitle='Subprogram-level associations; the frozen hierarchy BypassIndex is reported separately')
plot_save(p,'27m_external_efficacy',9,5)

## GSE127988 ligand displacement, respecting biological day --------------
d127<-x127
for(arm in unique(s127$arm))for(day in unique(s127$biological_set)) {
  jj<-which(s127$arm==arm&s127$biological_set==day)
  ref<-rowMeans(x127[,s127$arm==arm&s127$biological_set==day&s127$egf_h==0,drop=FALSE])
  d127[,jj]<-x127[,jj,drop=FALSE]-ref
}
s127scores<-external_score(d127)
long127<-melt(data.table(s127[,.(gsm,arm,egf_h,biological_set)],s127scores),id.vars=c('gsm','arm','egf_h','biological_set'),variable.name='program',value.name='score')
long127<-merge(long127,coverage[dataset=='GSE127988',.(n_frozen=.N,n_measured=sum(measured)),by=program],by='program',sort=FALSE)
long127[,complete_score:=ifelse(n_frozen==n_measured,score,NA_real_)]
write_table(long127,'gse127988_program_scores.tsv')
day127<-long127[,.(score=mean(score)),by=.(arm,egf_h,biological_set,program)]
write_table(day127,'gse127988_biological_day_means.tsv')
auc127<-day127[order(egf_h),.(auc=sum(diff(egf_h)*(head(score,-1)+tail(score,-1))/2)/8),by=.(arm,biological_set,program)]
write_table(auc127,'gse127988_auc_by_day.tsv')
attenuation<-dcast(auc127,program+biological_set~arm,value.var='auc')
setnames(attenuation,c('Vem','Vem + Cobi'),c('vem_auc','vem_cobi_auc'))
attenuation[,`:=`(difference=vem_cobi_auc-vem_auc,attenuation_fraction=ifelse(vem_auc>0,1-vem_cobi_auc/vem_auc,NA_real_))]
write_table(attenuation,'gse127988_auc_attenuation.tsv')
external_heat(d127,s127,'27n_gse127988',order(s127$cobi_uM,s127$egf_h,s127$biological_set,s127$gsm))
p<-ggplot(day127,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+
  geom_hline(yintercept=0,color='grey70')+geom_line(aes(linetype=factor(biological_set)))+geom_point(size=1.8)+facet_wrap(~program)+
  theme_classic(base_size=12)+labs(x='EGF exposure (h)',y='Frozen program displacement from own time zero',color='Drug background',linetype='Biological day',
  title='GSE127988 | frozen programs under EGF',subtitle='Technical wells averaged within day; two biological days shown separately')
plot_save(p,'27o_egf_time_courses',9,5)
p<-ggplot(auc127,aes(arm,auc,color=program,group=interaction(program,biological_set)))+geom_hline(yintercept=0,color='grey70')+
  geom_line(aes(linetype=factor(biological_set)))+geom_point(size=2)+facet_wrap(~program)+scale_color_manual(values=colors)+
  theme_classic(base_size=12)+labs(x=NULL,y='Normalized 0–8 h AUC',linetype='Biological day',title='EGF program activation with and without MEK inhibition')
plot_save(p,'27p_egf_auc',8,4.5)


## 12. Evaluate the frozen hierarchy, including absolute B suppression -------
# These predictors were written before the digitized Bliss outcomes were joined.
hier_benefit<-merge(hierarchy_predictors,bliss,by='cell_line')
by_ok<-is.finite(hier_benefit$BypassIndex)
by_test<-if(sum(by_ok)>=3)cor.test(hier_benefit$BypassIndex[by_ok],hier_benefit$bliss[by_ok],method='spearman',exact=FALSE) else NULL
pred_top<-hier_benefit[order(-BypassIndex),head(cell_line,3)]
observed_top<-hier_benefit[order(-bliss),head(cell_line,4)]
bypass_result<-data.table(n=sum(by_ok),spearman=if(is.null(by_test))NA_real_ else unname(by_test$estimate),
  nominal_p=if(is.null(by_test))NA_real_ else by_test$p.value,top3_in_top4=sum(pred_top %chin% observed_top),
  prediction='min[z(E PLX-C), z(R PLX-C)]',calibration='No Bliss-fitted weights or line removal')
write_table(hier_benefit,'gse68453_bypass_prediction_vs_bliss.tsv');write_table(bypass_result,'gse68453_bypass_test.tsv')
hi<-dcast(hierarchy_means684,cell_line+role~arm,value.var='score')
hi<-hi[cell_line %chin% factorial_lines]
hi[,`:=`(plx_minus_control=plx-control,adding_lapatinib=plx_lap-plx,interaction=plx_lap-plx-lapatinib+control)]
write_table(hi,'gse68453_hierarchy_factorial.tsv')
p<-ggplot(hierarchy_means684[cell_line %chin% factorial_lines],aes(factor(arm,levels=c('control','plx','lapatinib','plx_lap')),score,group=role,color=role))+
  geom_hline(yintercept=0,color='grey75')+geom_line()+geom_point(size=2)+facet_wrap(~cell_line,ncol=3)+
  scale_color_manual(values=c(B='#BA6B29',E='#188991',R='#8351A2'))+theme_classic(base_size=11)+
  theme(axis.text.x=element_text(angle=25,hjust=1))+labs(x=NULL,y='Frozen hierarchy score',title='B / E / R across the five complete factorials')
plot_save(p,'27q_hierarchy_factorial',10,7)
p<-ggplot(hier_benefit,aes(BypassIndex,bliss))+geom_point(size=2,color='#188991')+
  geom_text(aes(label=cell_line),size=3,vjust=-.6,check_overlap=TRUE)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#8351A2')+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+labs(x='Frozen BypassIndex',y='Digitized Bliss benefit',title='Hierarchy-based combination prediction',
  subtitle=paste0('Spearman rho = ',round(bypass_result$spearman,2),'; top-three/top-four overlap = ',bypass_result$top3_in_top4))
plot_save(p,'27r_bypass_prediction',7,5.7)
# GSE127988 actually contains untreated controls on both time-course days.
# Use all those controls to assess whether B remains suppressed relative to
# untreated cells. A positive EGF increment alone does not mean B was restored.
full_time<-s127_all[biological_set %in% 3:4]
absolute127<-x127_full[,full_time$gsm,drop=FALSE]
for(day in unique(full_time$biological_set)) {
  jj<-which(full_time$biological_set==day)
  refs<-full_time[biological_set==day & vem_uM==0 & cobi_uM==0 & egf_h==0,gsm]
  check(length(refs)==4L,'Unexpected untreated controls per biological day')
  absolute127[,jj]<-absolute127[,jj,drop=FALSE]-rowMeans(x127_full[,refs,drop=FALSE])
}
h_abs<-hierarchy_score(absolute127,partial=TRUE)
h_inc<-hierarchy_score(d127,partial=TRUE)
h127abs<-melt(data.table(full_time[,.(gsm,arm,egf_h,biological_set)],h_abs),id.vars=c('gsm','arm','egf_h','biological_set'),variable.name='role',value.name='score')
h127inc<-melt(data.table(s127[,.(gsm,arm,egf_h,biological_set)],h_inc),id.vars=c('gsm','arm','egf_h','biological_set'),variable.name='role',value.name='score')
write_table(h127abs,'gse127988_hierarchy_vs_untreated.tsv');write_table(h127inc,'gse127988_hierarchy_egf_increments.tsv')
habsday<-h127abs[,.(score=mean(score)),by=.(role,arm,egf_h,biological_set)]
hincday<-h127inc[,.(score=mean(score)),by=.(role,arm,egf_h,biological_set)]
hauc<-hincday[order(egf_h),.(auc=sum(diff(egf_h)*(head(score,-1)+tail(score,-1))/2)/8),by=.(role,arm,biological_set)]
write_table(hauc,'gse127988_hierarchy_auc_by_day.tsv')
hcontrast<-dcast(hauc,role+biological_set~arm,value.var='auc');setnames(hcontrast,c('Vem','Vem + Cobi'),c('vem_auc','vem_cobi_auc'))
hcontrast[,`:=`(cobimetinib_change=vem_cobi_auc-vem_auc,fraction_retained=ifelse(vem_auc>0,vem_cobi_auc/vem_auc,NA_real_))]
write_table(hcontrast,'gse127988_hierarchy_cobimetinib_effect.tsv')
bsupp<-habsday[role=='B'&arm!='DMSO',.(baseline_B=score[egf_h==0],maximum_B_after_EGF=max(score[egf_h>0]),
  B_remains_below_untreated=all(score<0)),by=.(arm,biological_set)]
write_table(bsupp,'gse127988_B_suppression_vs_untreated.tsv')
p<-ggplot(habsday[arm!='DMSO'],aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+
  geom_hline(yintercept=0,color='grey60',linetype=2)+geom_line(aes(linetype=factor(biological_set)))+geom_point(size=1.8)+facet_wrap(~role,scales='free_y')+
  theme_classic(base_size=12)+labs(x='EGF exposure (h)',y='Hierarchy score relative to untreated controls',linetype='Biological day',color='Background',
  title='Does B remain suppressed while receptor/shared responses return?',subtitle=paste('Same-day untreated zero | measured genes:',paste(hierarchy_coverage[dataset=='GSE127988',paste0(role,' ',n_measured,'/',n_frozen)],collapse=', '),'| days separate'))
plot_save(p,'27s_hierarchy_absolute_time_courses',11,5)
p<-ggplot(hincday,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+geom_hline(yintercept=0,color='grey70')+
  geom_line(aes(linetype=factor(biological_set)))+geom_point(size=1.8)+facet_wrap(~role,scales='free_y')+theme_classic(base_size=12)+
  labs(x='EGF exposure (h)',y='Increment from own drug-background time zero',color='Background',linetype='Biological day',
  title='EGF-induced hierarchy displacement | no assumed MEK-inhibitor root-off state')
plot_save(p,'27s_hierarchy_egf_increments',11,5)

# Preserve gene-level contributions so scoring/transfer failures can be localized.
contributions<-rbindlist(lapply(3:4,function(day) {
  ids<-full_time[biological_set==day & arm=='Vem' & egf_h==0,gsm]
  raw<-rowMeans(absolute127[,ids,drop=FALSE])
  data.table(gene=leaf_order,role=membership$role,program=membership$program,biological_day=day,
    discovery_V=fp[leaf_order,'vemurafenib'],discovery_EGF=fp[leaf_order,'EGF'],
    frozen_sign=gene_sign[leaf_order],external_V_minus_DMSO=raw,signed_contribution=raw*gene_sign[leaf_order])
}))
write_table(contributions,'gse127988_gene_transfer_diagnostic.tsv')
saveRDS(list(membership=membership,fp=fp[chosen,,drop=FALSE],gene_sign=gene_sign,
  effects=effects[,chosen,drop=FALSE],raw_effects=raw_effects[,chosen,drop=FALSE],
  background=background[,chosen,drop=FALSE],background_meta=background_meta,auc=auc[,chosen,drop=FALSE],meta=meta,line_drug=line_drug,
  study_effect=cyto_x,source_rows=source_rows,by_drug_cor=by_drug_cor,
  absolute127=absolute127,d127=d127,d684=d684,s127=s127,full_time=full_time,s684=s684),file.path(art,'verification_inputs.rds'))

for(n in c('effects','raw_effects','auc','background')) {
  xx<-get(n)[,chosen,drop=FALSE]
  fwrite(data.table(id=if(is.null(rownames(xx)))seq_len(nrow(xx)) else rownames(xx),xx),file.path(art,paste0('verify_',n,'.tsv.gz')),sep='\t',na='NA')
}
for(n in c('cyto_x','d127','absolute127','d684')) {
  xx<-get(n);fwrite(data.table(gene=rownames(xx),xx),file.path(art,paste0('verify_',n,'.tsv.gz')),sep='\t',na='NA')
}
write_table(root_test[,.(compound,cell_name,within_drug_score,within_drug_survival)],'root_adjustment_diagnostic.tsv')
# Localize the combination-prediction mismatch without fitting new Bliss weights.
combination_diagnostic<-rbindlist(lapply(c('B','E','R'),function(r){
  xx<-hier_benefit[[r]]
  data.table(role=r,n_lines=length(xx),positive_displacements=sum(xx>0),
    spearman_with_Bliss=cor(xx,hier_benefit$bliss,method='spearman'),minimum=min(xx),maximum=max(xx))
}))
write_table(combination_diagnostic,'combination_prediction_failure_diagnostic.tsv')
p<-ggplot(melt(hier_benefit,id.vars=c('cell_line','bliss','BypassIndex'),measure.vars=c('B','E','R'),variable.name='role',value.name='displacement'),
  aes(displacement,bliss))+geom_vline(xintercept=0,color='grey75')+geom_point(size=2,color='#287D99')+
  geom_text(aes(label=cell_line),size=2.7,vjust=-.7,check_overlap=TRUE)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#B64039')+
  facet_wrap(~role,scales='free_x')+theme_classic(base_size=11)+theme(aspect.ratio=1)+
  labs(x='Frozen PLX-minus-control displacement',y='Digitized Bliss benefit',title='Which part limits combination prediction?',
    subtitle='Diagnostic only: no refitting, gene changes or line exclusions')
plot_save(p,'27r_combination_failure_diagnosis',11,4.5)

## 13. Save a clean, inspectable experiment ----------------------------------
check(identical(readRDS(file.path(art,'frozen_programs.rds'))$membership,as.data.frame(membership)),'Frozen genes/weights changed during validation')
check(!any(line_drug$depmap_id[c(fit_rows,held)] %chin% EXCLUDED_IDS),'Excluded outlier entered analysis')
check(!length(intersect(line_drug$depmap_id[fit_rows],line_drug$depmap_id[held])),'Development and holdout share a line')
write_table(data.table(check=c('Two confirmed outliers excluded','Three V600 melanomas contribute molecular target-engagement evidence',
  'One primary-PRISM assay used consistently','Gene and program weights frozen before validation','Whole-line holdout separation',
  'Full external identifier and missingness audit','Untreated GSE127988 controls retained for B suppression','PAM membership and frozen gene order explicit'),passed=TRUE),'verification.tsv')
write_table(data.table(step=c('Source extension','Discovery','Hierarchy','Primary PRISM holdout','GSE68453 factorial and efficacy','GSE127988 controls and time courses'),
  status=c('Two processed melanoma lines added',paste(length(chosen),'genes;',K,'subprograms'),
    if(hierarchy_resolved)'All centroid patterns agree with their hypothesis-defined families; external evaluation remains separate' else 'Some candidate roles unresolved; inspect the centroid evidence',rep('Executed',3))),'execution_status.tsv')
decisions<-data.table(failure=c('Unnumbered pilots and mixed outputs','Two previously excluded technical outliers restored',
  'Only one BRAF-mutant melanoma in the secondary-PRISM intersection','Gene-wise viability correlation deletes target-engagement markers',
  'Strong common response obscures weaker branches','A global hclust cut leaves singleton branches','NRG1 profiles dominated by one study',
  'Cobimetinib was assumed to force the shared root off','Positive B change was confused with restored B activity'),
  repair=c('Numbered stage-27 R/Python scripts with a fresh matching results/artifacts folder','Both remain excluded; new-line QC is exported',
  'Retrieve processed C32/RPMI-7951 profiles; use one locally cached primary-PRISM assay at 2.5 uM',
  'B markers use consistent melanoma target engagement; functional weighting enters receptor/root coordinates and is tested at program level',
  'Select explicit B-specific, E-specific and shared-response contrasts, 40 candidates per role; never claim this prior-guided selection independently proves the Y',
  'Preserve the three admission families and split each into two PAM subprograms; show unrestricted mixing and centroid agreement separately','Use study-ligand medians and retain source context/timing in the audit',
  'Gerosa describes receptor-driven RAF/MEK resistance: measure the response instead of imposing 50% attenuation',
  'Use the available same-day untreated controls to determine whether B remains below baseline'))
decisions<-rbind(decisions,data.table(
  failure=c('B markers respond to many unrelated compounds on the same plates','B specificity compared different cellular contexts',
    'One sign per program cancels induced and repressed member genes','Pooled functional association adjusted only one variable for drug',
    'Transparent hierarchy PNG hid edges against a dark background','Early root displacement does not rank combination benefit'),
  repair=c('Subtract each cell/plate median of non-target compound responses; retain raw measurements and calibration audits. This removes common response, not a proven purely technical batch effect.',
    'Compare B against L/C in the same three melanomas using a common drug-effect scale; assign overlapping candidate genes by strongest eligible contrast.',
    'Freeze gene-level signs from discovery V or EGF effects; show raw and signed correlations separately.',
    'Report partial Spearman after residualizing both ranked variables for drug; preserve the original marginal statistic and held-out prediction errors.',
    'Save a white background and inspect the rendered Y and heatmaps.',
    'Report E/R contributions separately: distinguish early EGF transcription, basal receptor competence and three-day Bliss synergy. Do not relabel weak efficacy prediction as failure of the molecular Y.')))
write_table(decisions,'failure_and_repair_log.tsv')
write_table(data.table(figure=sort(list.files(out,pattern='^27.*\\.pdf$'))),'FIGURE_INDEX.tsv')
writeLines(c('# Stage 27: BRAF–ERBB hierarchy repair','',
  sprintf('Current result: B remains below untreated levels in all %d drug-background/day combinations. EGF-induced E and R are positive on both biological days. Root AUC is attenuated %.1f–%.1f%%; E retains %.1f–%.1f%%.',
    sum(bsupp$B_remains_below_untreated),100*min(1-hcontrast[role=='R',fraction_retained]),100*max(1-hcontrast[role=='R',fraction_retained]),
    100*min(hcontrast[role=='E',fraction_retained]),100*max(hcontrast[role=='E',fraction_retained])),
  sprintf('Primary PRISM: drug-adjusted partial Spearman %.3f; RMSE %.3f versus %.3f for drug alone. Combination-benefit prediction remains weaker: Spearman %.3f across twelve lines.',
    root_comparison$partial_spearman_controlling_drug,root_comparison$root_rmse,root_comparison$baseline_rmse,bypass_result$spearman),
  'These are retrospective molecular checks. Gene coverage is reported for each dataset and family; missing genes produce explicit observed-gene partial scores.','',
  '[Results viewer](index.html) | [R analysis](../../scripts/27_repair_braf_erbb_hierarchy.R) | [Processed-data preparation](../../scripts/27_prepare_braf_context.py)','',
  'Run from /work:','',
  '```sh','python3 tahoe_prism_v2/scripts/27_prepare_braf_context.py','python3 tahoe_prism_v2/scripts/27_calibrate_plate_background.py','Rscript tahoe_prism_v2/scripts/27_repair_braf_erbb_hierarchy.R','python3 tahoe_prism_v2/scripts/27_verify_braf_erbb_hierarchy.py','```','',
  '## What was repaired','',
  'The old ranking criterion was incompatible with the bypass hypothesis: BRAF target-engagement markers can be suppressed in both viable bypassing cells and sensitive cells, so they need not correlate with viability individually. This stage uses target engagement for the BRAF coordinate and tests functional relevance at the root/program level.',
  'The original secondary-PRISM intersection had one V600 melanoma, LOX-IMVI, plus two V600 colorectal lines. C32 and RPMI-7951 have primary-PRISM responses and processed Tahoe measurements. Only the three requested drug/dose profiles were retrieved by Parquet byte ranges (no raw single-cell data). NCI-H2122 and NCI-H596 remain excluded.',
  'All functional outcomes now come from primary PRISM at 2.5 uM, for both original and additional lines. Molecular AUC still uses all Tahoe doses with 0.25/0.50/0.25 weights. These are distinct quantities; the functional outcome is not a dose AUC.',
  'All retained cell lines receive sample-local calibration: subtract the median response to other compounds on the same cell/plate before dose integration. The three target compounds are excluded from calibration. This removes a common perturbation response that may contain both technical and generic biological components; it is not a claim that all of it is technical noise.',
  'The BRAF coordinate is median molecular AUC across the three V600 melanomas. Lapatinib and cobimetinib coordinates are mean molecular AUC weighted by signed primary-PRISM sensitivity, divided by the sum of absolute sensitivity. EGF is the median of study-level effects. Coordinates are divided by across-gene SD without subtracting their means, retaining biological zero.',
  'Genes must be measurable in the relevant drug cohort (>=0.1 counts per 10,000 in at least half of observations), present in at least half the CytoSig study-ligand summaries, and retain the relevant response direction under leave-one-line/study deletion. The abundance rule is a measurement heuristic, not a calibrated assay detection limit.',
  'B specificity = |V| − max(|L|,|C|,|EGF|), comparing V/L/C medians in the same three melanomas on a common drug-effect SD scale. E specificity = min(|EGF|,|L|) − max(|V|,|C|), requiring opposite EGF/lapatinib directions. Shared strength = min(|V|,|C|,|EGF|), requiring EGF opposite V and C. Each gene first takes its strongest eligible contrast; up to 40 genes per family enter the union, without duplicate membership. This explicitly uses the biological hypothesis; matching discovery centroids cannot independently validate the hypothesis.',
  'Programs use actual sample correlations, oriented using each gene’s discovery-derived direction. Both raw and signed matrices are displayed with the same permanent order. Within-drug Tahoe matrices (three V600 melanomas for vemurafenib; retained development lines for other drugs) receive equal weight; their average and the CytoSig study-level gene correlation matrix each get one half. Each hypothesis-defined admission family is partitioned into two PAM subprograms. A global hclust cut in the initial uncalibrated analysis produced a large block plus singleton branches, while unrestricted PAM merged receptor markers into the shared response. Hierarchical clustering supplies the within-block order and a separately rendered gene tree. The three-family topology is imposed by the biological hypothesis, not discovered anew by clustering.',
  'Scores are equal-weight signed gene means. A B gene is oriented opposite its own discovery V effect; E/R genes follow their own discovery EGF effect. B is an operational BRAFi-specific transcriptional response, including translation/metabolic genes; it is not a direct BRAF kinase-activity assay.',
  'Subprogram family membership is inherited from the admission contrast. Observed centroid roles are computed separately, and their agreement is exported without relabeling discordant subprograms. The Y is a candidate, perturbation-informed hierarchy, not a causal gene tree. Coherence and external responses are the checks beyond the selection criteria.','',
  '## Validation and interpretation','',
  'Genes, directions, weights and hierarchy assignments are written before primary-PRISM holdout evaluation or either external expression matrix is loaded. The same seven nominal holdout cell lines remain together across drugs; six have usable primary-PRISM measurements, because SW48 lacks this assay outcome. Root and subprogram models are all reported; external outcomes never select weights or clusters. These public data have been inspected in earlier analyses, so results are retrospective, not a prospectively untouched validation.',
  'External probes map through approved/unique historical HGNC symbols. The largest-IQR array probe is chosen before efficacy is joined. RNA Ensembl IDs are preferred over names; duplicated gene rows sum on the RPM scale before log2(1+RPM). Complete scores and observed-gene partial scores are distinguished, and every missing gene is listed.',
  'The expression-only BypassIndex is written before joining the locked digitized Bliss values. It uses min[z(E),z(R)] across the twelve PLX-versus-control profiles. Neither weights, cell lines nor doses are optimized against Bliss.',
  'GSE127988 has 96 individual wells, not just the 56 selected treatment time-course wells. The extra same-day untreated controls allow direct assessment of B suppression; a positive increment from vemurafenib time zero does not by itself establish recovery to untreated B activity. Technical wells are averaged within two biological days, which remain separate in plots.',
  'Gerosa reports receptor-driven MAPK resistance under RAF/MEK inhibition. In these RNA data, however, the shared transcript response is strongly attenuated by cobimetinib. Kinase pulses and time-integrated transcription are distinct readouts; report the observed attenuation, not an assumed complete pathway-off state.',
  'Combination prediction: all twelve early R displacements are negative. Across-line z scoring describes relative retention, not absolute root activation. E shows more association with Bliss than R; the minimum combines a receptor signal with an early root readout that does not rank three-day synergy. Capaldo also reports basal NRG1-related resistance and distinct single-agent response phenotypes. Testing receptor competence/basal state and matched later-time function requires a separately evaluated predictor; those data are not retroactively fitted here.',
  'No expression program alone proves kinase activity, survival causality, apoptosis or a universal pathway hierarchy. The study provides independent mechanistic support; these analyses test transfer of the molecular programs and their association with functional outcomes.','',
  '## Sources','',
  '- [Gerosa et al., Cell Systems](https://pubmed.ncbi.nlm.nih.gov/33113355/)',
  '- [Capaldo et al., PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0138210)',
  '- [GSE127988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988)',
  '- [GSE68453](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68453)'),file.path(out,'README.md'))
html_escape<-function(x){x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
html_table<-function(x)paste0('<table><tr>',paste0('<th>',names(x),'</th>',collapse=''),'</tr>',paste(apply(as.data.frame(x),1,function(r)paste0('<tr>',paste0('<td>',html_escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table>')
figures<-sort(list.files(out,pattern='^27.*\\.pdf$'));tables<-sort(list.files(out,pattern='\\.tsv(\\.gz)?$'))
# Show hierarchy and the most informative validation panels first.
priority<-c('27h_candidate_hierarchy.pdf','27f_joint_signed_correlations.pdf','27q_hierarchy_factorial.pdf','27s_hierarchy_absolute_time_courses.pdf','27r_bypass_prediction.pdf')
figures<-c(intersect(priority,figures),setdiff(figures,priority))
cards<-vapply(figures,function(f){png<-sub('\\.pdf$','.png',f);paste0('<article><h3>',html_escape(sub('\\.pdf$','',f)),'</h3>',
  if(file.exists(file.path(out,png)))paste0('<a href="',png,'"><img loading="lazy" src="',png,'"></a>') else '',
  '<p><a href="',f,'">Open PDF</a></p></article>')},character(1))
code_path<-file.path(root,'scripts/27_repair_braf_erbb_hierarchy.R')
writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Stage 27 — hierarchy repair</title>',
  '<style>body{font:16px/1.5 system-ui;background:#f5f7f9;color:#203040;margin:0}main{max-width:1350px;margin:auto;padding:28px}a{color:#075b8e}.status{background:#e9f0f6;border-left:5px solid #287D99;padding:16px}.grid{display:grid;grid-template-columns:1fr;gap:18px}article{background:white;padding:16px;border:1px solid #d5dfe6;overflow:hidden}h3{font-size:15px;overflow-wrap:anywhere}img{display:block;max-width:100%;height:auto;margin:auto}table{border-collapse:collapse;width:100%;font-size:14px;background:white}th,td{padding:8px;border:1px solid #ddd;text-align:left}pre{background:#182530;color:#e4edf5;padding:16px;overflow:auto;max-height:700px}.tables{columns:2}summary{cursor:pointer;font-weight:bold}</style><main>',
  '<h1>Stage 27 · BRAF–ERBB hierarchy repair</h1><p><a href="../../scripts/27_repair_braf_erbb_hierarchy.R">R analysis</a> · <a href="../../scripts/27_prepare_braf_context.py">Processed-data preparation</a> · <a href="README.md">Methods and limits</a></p>',
  paste0('<div class="status"><strong>',length(chosen),' genes · six subprograms · three V600 melanoma contexts</strong><p>The hierarchy is now tested using target engagement, branch specificity and a shared response. Its discovery labels are hypothesis-guided; external prediction and perturbation behavior are reported separately below.</p></div>'),
  '<h2>Why the previous implementation lost the hierarchy</h2>',html_table(decisions),
  '<h2>Candidate subprogram assignments</h2>',html_table(role_audit[,.(program,role,observed_centroid_role,centroid_matches_family,n_genes,braf_specificity,erbb_specificity)]),
  '<h2>Actual source coherence (signed activity)</h2>',html_table(context_coherence),
  '<h2>External score coverage</h2>',html_table(hierarchy_coverage),
  '<h2>EGF and cobimetinib displacement</h2>',html_table(hcontrast),
  '<h2>Primary-PRISM root transfer</h2>',html_table(root_comparison),
  '<h2>Frozen combination prediction</h2>',html_table(bypass_result),html_table(combination_diagnostic),
  '<h2>Does B remain below untreated levels?</h2>',html_table(bsupp),
  '<h2>Figures</h2><div class="grid">',cards,'</div>',
  '<h2>Tables and audits</h2><ul class="tables">',paste0('<li><a href="',tables,'">',tables,'</a></li>'),'</ul>',
  '<p><a href="../../scripts/27_calibrate_plate_background.py">Plate-calibration source</a> · <a href="../../scripts/27_verify_braf_erbb_hierarchy.py">Independent verification source</a> · <a href="independent_verification.tsv">Independent verification results</a></p>',
  '<h2>Code</h2><details><summary>Expand the complete R analysis</summary><pre><code>',paste(html_escape(readLines(code_path,warn=FALSE)),collapse='\n'),'</code></pre></details>',
  '<details><summary>Expand the processed-data preparation script</summary><pre><code>',paste(html_escape(readLines(file.path(root,'scripts/27_prepare_braf_context.py'),warn=FALSE)),collapse='\n'),'</code></pre></details></main></html>'),file.path(out,'index.html'))
message('Stage 27 completed: candidate hierarchy, source structure and external tests saved to ',out)
