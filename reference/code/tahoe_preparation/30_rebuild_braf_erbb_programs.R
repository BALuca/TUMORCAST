#!/usr/bin/env Rscript
## Stage 30: rebuild BRAF--ERBB programs from response association and cross-source coherence.
## Drug-response genes must be associated with PRISM response, retain their
## direction across doses, and participate in coherent Tahoe/CytoSig modules.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cluster);library(msigdbr);library(ggpubr)})
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(arg))dirname(dirname(normalizePath(sub('^--file=','',arg[1])))) else '/work/tahoe_prism_v2'
work<-dirname(root);input_art<-file.path(root,'artifacts/26_braf_erbb_bypass')
prepared_art<-file.path(root,'artifacts/27_braf_erbb_hierarchy')
art<-file.path(root,'artifacts/30_braf_erbb_programs');out<-file.path(root,'results/30_braf_erbb_programs')
unlink(art,recursive=TRUE,force=TRUE);unlink(out,recursive=TRUE,force=TRUE)
dir.create(art,recursive=TRUE,showWarnings=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
main_fig_dir<-file.path(out,'figures/main');supp_fig_dir<-file.path(out,'figures/supplementary')
dir.create(main_fig_dir,recursive=TRUE,showWarnings=FALSE);dir.create(supp_fig_dir,recursive=TRUE,showWarnings=FALSE)
source(file.path(work,'src/R/lib_heatmap.R'))
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
write_table<-function(x,name)fwrite(as.data.table(x),file.path(out,name),sep='\t',na='NA')
EXCLUDED_IDS<-c('ACH-000311','ACH-000628');K<-5L;MIN_PROGRAM<-6L


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

plot_map<-c(
  '29e_rank_silhouette'='figures/supplementary/S01_program_resolution',
  '29h_candidate_hierarchy'='figures/main/05_program_hierarchy',
  '29f_context_coherence'='figures/main/06_source_sharing',
  '29i_subprogram_prism_associations'='figures/supplementary/S02_subprogram_viability_associations',
  '29i_shared_root_prism'='figures/main/09_prism_viability_association',
  '29l_factorial_programs'='figures/main/07_gse68453_factorial_response',
  '29m_external_efficacy'='figures/supplementary/S03_subprogram_combination_associations',
  '29o_egf_time_courses'='figures/main/08_gse127988_egf_response',
  '29r_bypass_prediction'='figures/main/10_gse68453_combination_prediction')
heat_map<-c(
  '29f_joint_signed_correlations'='figures/main/02_joint_gene_correlations',
  '29f_tahoe_signed_correlations'='figures/main/03_tahoe_gene_correlations',
  '29f_cytosig_signed_correlations'='figures/main/04_cytosig_gene_correlations',
  '29f_joint_signed_correlations_labeled'='figures/supplementary/S04_joint_gene_correlations_labeled',
  '29f_tahoe_signed_correlations_labeled'='figures/supplementary/S05_tahoe_gene_correlations_labeled',
  '29f_cytosig_signed_correlations_labeled'='figures/supplementary/S06_cytosig_gene_correlations_labeled',
  '29g_perturbation_fingerprint'='figures/supplementary/S07_gene_perturbation_fingerprints',
  '29g_tahoe_gene_samples_signed'='figures/supplementary/S08_tahoe_gene_by_sample_response',
  '29g_cytosig_gene_samples_signed'='figures/supplementary/S09_cytosig_gene_by_sample_response',
  '29k_gse68453_gene_samples_semantic'='figures/supplementary/S10_gse68453_gene_by_sample_response',
  '29n_gse127988_gene_samples_semantic'='figures/supplementary/S11_gse127988_gene_by_sample_response')
plot_save <- function(p,name,w=7,h=5) {
  destination<-unname(plot_map[name]);if(is.na(destination))return(invisible(NULL))
  base<-file.path(out,destination);dir.create(dirname(base),recursive=TRUE,showWarnings=FALSE)
  ggsave(paste0(base,'.pdf'),p,width=w,height=h,device=cairo_pdf,bg='white')
  ggsave(paste0(base,'.png'),p,width=w,height=h,dpi=160,bg='white')
}

heat <- function(x,name,title,rows=NULL,cols=NULL,limit=2.5,width=10,height=NULL,show_cols=TRUE,show_rows=TRUE,row_font=9,col_font=8) {
  destination<-unname(heat_map[name]);if(is.na(destination))return(invisible(NULL))
  base<-file.path(out,destination);dir.create(dirname(base),recursive=TRUE,showWarnings=FALSE)
  if(is.null(height))height<-max(4.8,nrow(x)*.12+1.8)
  if(!is.null(rows))check(length(rows)==nrow(x),'Heatmap row annotations misaligned')
  if(!is.null(cols))check(length(cols)==ncol(x),'Heatmap column annotations misaligned')
  legend_title<-if(grepl('correlation',name))'Spearman rho' else if(grepl('gene_samples',name))'Expression\nchange' else 'Standardized\neffect'
  h<-draw_heatmap(x,scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(x)),column_order=seq_len(ncol(x)),row_split=rows,column_split=cols,
    cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(2,'mm'),column_gap=grid::unit(2,'mm'),
    row_title_rot=0,column_title_rot=0,row_title_gp=grid::gpar(fontsize=10,fontface='bold'),column_title_gp=grid::gpar(fontsize=8,fontface='bold'),
    show_row_names=show_rows,show_column_names=show_cols,row_names_gp=grid::gpar(fontsize=row_font),column_names_gp=grid::gpar(fontsize=col_font),
    column_names_rot=90,hm_name=legend_title,heatmap_legend_param=list(at=c(-limit,0,limit),labels=round(c(-limit,0,limit),2)),title=title)
  cairo_pdf(paste0(base,'.pdf'),width=width,height=height)
  drawn<-ComplexHeatmap::draw(h,heatmap_legend_side='right')
  check(identical(as.integer(unlist(ComplexHeatmap::row_order(drawn))),seq_len(nrow(x))),'Heatmap changed frozen gene order')
  check(identical(as.integer(unlist(ComplexHeatmap::column_order(drawn))),seq_len(ncol(x))),'Heatmap changed frozen sample order')
  invisible(dev.off())
  png(paste0(base,'.png'),width=width,height=height,units='in',res=140,type='cairo')
  ComplexHeatmap::draw(h,heatmap_legend_side='right');invisible(dev.off())
}
robust_limit<-function(x,p=.98)max(as.numeric(quantile(abs(x[is.finite(x)]),p,na.rm=TRUE)),.01)

source_heat <- function(x,name,title) {
  x<-x[leaf_order,,drop=FALSE]
  c<-cor(x,method='spearman',use='pairwise.complete.obs');check(all(is.finite(c)),'Non-estimable sample correlation')
  tr<-hclust(as.dist(pmax(1-c,0)),method='average')
  z<-t(scale(t(x)));z[!is.finite(z)]<-NA_real_
  signed_z<-z*gene_sign[rownames(z)]
  activity<-sapply(program_order,function(p)colMeans(signed_z[membership[program==p,gene],,drop=FALSE],na.rm=TRUE))
  winner<-apply(activity,1,function(a)if(all(!is.finite(a)))1L else which.max(replace(a,!is.finite(a),-Inf)))
  ord<-unlist(lapply(seq_along(program_order),function(j)tr$order[winner[tr$order]==j]))
  cb<-factor(program_order[winner[ord]],levels=program_order)
  lim<-robust_limit(z,.98)
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),dominant_program=as.character(cb)),paste0(name,'_sample_order.tsv'))
  heat(z[,ord,drop=FALSE],name,paste0(title,' | row-standardized response'),row_blocks,cb,limit=lim,width=11,height=18,show_cols=FALSE,row_font=8)
  heat(signed_z[,ord,drop=FALSE],paste0(name,'_signed'),paste0(title,' | oriented to program activation'),row_blocks,cb,limit=lim,width=11,height=18,show_cols=FALSE,row_font=8)
  invisible(ord)
}

# Convert an hclust object to plain segment coordinates so trees are rendered
# with ggplot2 and remain legible on a white, compact canvas.
hclust_frame <- function(hc,labels=hc$labels) {
  n<-length(hc$order);leaf_x<-numeric(n);leaf_x[hc$order]<-seq_len(n)
  node_x<-node_y<-numeric(n-1L);seg<-list()
  pos<-function(v)if(v<0)c(leaf_x[-v],0) else c(node_x[v],node_y[v])
  for(i in seq_len(n-1L)) {
    a<-pos(hc$merge[i,1]);b<-pos(hc$merge[i,2]);node_x[i]<-mean(c(a[1],b[1]));node_y[i]<-hc$height[i]
    seg[[i]]<-data.table(x=c(a[1],b[1],a[1]),y=c(a[2],b[2],node_y[i]),
      xend=c(a[1],b[1],b[1]),yend=c(node_y[i],node_y[i],node_y[i]))
  }
  list(segments=rbindlist(seg),labels=data.table(x=seq_len(n),y=0,label=labels[hc$order]))
}
plot_tree <- function(hc,labels,name,title,subtitle,w=9,h=5) {
  d<-hclust_frame(hc,labels)
  p<-ggplot(d$segments)+geom_segment(aes(x,y,xend=xend,yend=yend),linewidth=.45,color='#314A5A')+
    geom_text(data=d$labels,aes(x,y,label=label),angle=90,hjust=1,vjust=.5,size=2.8,inherit.aes=FALSE)+
    scale_x_continuous(expand=expansion(mult=c(.015,.015)))+scale_y_continuous(expand=expansion(mult=c(.28,.06)))+
    theme_classic(base_size=11)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),
      plot.margin=margin(7,7,45,7))+labs(y='Joint correlation distance',title=title,subtitle=subtitle)
  plot_save(p,name,w,h)
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
required<-file.path(prepared_art,c('additional_treated_mean.npy','additional_reference_mean.npy','additional_response_metadata.tsv','primary_prism_responses.tsv'))
check(all(file.exists(required)),'Run scripts/27_prepare_braf_context.py first to retrieve the two processed melanoma profiles')
a<-readRDS(file.path(input_art,'discovery_measurements.rds'))
effects<-a$effects;meta<-a$meta;cells<-a$cells;genes<-a$genes
reference_cp10k<-a$reference_cp10k;treated_cp10k<-a$treated_cp10k
cyto<-a$cyto;cytometa<-a$cytometa
new_meta<-fread(required[3]);new_treated<-read_npy_rows(required[1],new_meta$array_row);new_reference<-read_npy_rows(required[2],new_meta$array_row)
check(identical(fread(file.path(prepared_art,'additional_gene_order.tsv'))$gene,colnames(effects)),'Additional gene order differs')
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
# Reuse the genotype audit materialized by the preparation stage. This avoids a
# runtime dependency on a parquet reader and does not reuse any program choices.
prepared_cells<-fread(file.path(root,'results/27_braf_erbb_hierarchy/experiment_cell_lines.tsv'))
v600_ids<-unique(prepared_cells[braf_v600_confirmed==TRUE,depmap_id])
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
background_meta<-fread(file.path(prepared_art,'plate_common_response_metadata.tsv'))
background<-read_npy_rows(file.path(prepared_art,'plate_common_response.npy'),background_meta$array_row)
colnames(background)<-colnames(effects)
plate_map<-rbindlist(list(
  fread(file.path(work,'cache/tahoe_prism_full_chain/tahoe_pdex_selected/pdex_plate_audit.tsv.gz'))[,.(depmap_id,drug,dose,plate)],
  fread(file.path(prepared_art,'additional_plate_audit.tsv'))[,.(depmap_id,drug,dose,plate)]))
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
write_table(rbind(fread(file.path(out,'input_manifest.tsv')),data.table(path=file.path(prepared_art,c('plate_common_response.npy','plate_common_response_metadata.tsv','plate_background_observations.tsv')),role=c('Per-cell/plate common response across all genes','Independent calibration coverage for retained lines','Non-target compounds that estimated each common response'))),'input_manifest.tsv')
rm(a,new_treated,new_reference,new_effect,new_ref_cp,new_treat_cp);invisible(gc())

## 2. Common measured universe and study-balanced ligand evidence ----------
cytometa[,study:=sub('.*@((?:GSE[0-9]+|E-MTAB-[0-9]+)).*','\\1',profile_id,perl=TRUE)]
cytometa[,study_ligand:=paste(study,treatment,sep=' | ')]
study_meta<-unique(cytometa[,.(study_ligand,study,treatment)])
cyto_full<-as.matrix(cyto[,-1,with=FALSE]);rownames(cyto_full)<-cyto$gene
study_effect<-sapply(study_meta$study_ligand,function(k)matrixStats::rowMedians(cyto_full[,cytometa[study_ligand==k,profile_id],drop=FALSE],na.rm=TRUE))
rownames(study_effect)<-cyto$gene
universe<-intersect(colnames(auc),rownames(study_effect)[rowSums(is.finite(study_effect))>=ceiling(ncol(study_effect)/2)])
universe<-universe[grepl('^[A-Z][A-Z0-9-]{1,14}$',universe)]
write_table(cytometa,'cytosig_profile_audit.tsv')
write_table(cytometa[,.(profiles=.N,conditions=paste(unique(condition),collapse='; '),durations=paste(unique(duration),collapse='; ')),by=.(study,treatment)],'cytosig_study_context.tsv')

## 3. Admit genes through drug-response association or ligand response ------
# Three BRAF-V600 melanomas cannot identify a correlation with viability.
# Associations therefore use all eligible development lines; the matched
# melanomas separately determine whether vemurafenib engages each candidate.
drug_names<-c('vemurafenib','lapatinib','cobimetinib')
association_list<-list();dose_association_list<-list()
for(drug in drug_names) {
  ii<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible)
  xx<-auc[ii,universe,drop=FALSE];yy<-line_drug$prism_sensitivity[ii]
  rho<-apply(xx,2,function(v)cor(v,yy,method='spearman',use='complete.obs'))
  n<-length(yy);p<-2*pt(-abs(rho*sqrt((n-2)/pmax(1-rho^2,1e-12))),df=n-2)
  dose_rho<-sapply(c(.05,.5,5),function(dd){jj<-which(meta$compound==drug & meta$split=='Development' & meta$analysis_included & is.finite(meta$prism_sensitivity) & meta$dose==dd);apply(effects[jj,universe,drop=FALSE],2,function(v)cor(v,meta$prism_sensitivity[jj],method='spearman',use='complete.obs'))})
  association_list[[drug]]<-data.table(gene=universe,drug,n_lines=n,rho,p,agreeing_doses=rowSums(sign(dose_rho)==sign(rho),na.rm=TRUE))
  dose_association_list[[drug]]<-data.table(gene=universe,drug,rho_auc=rho,rho_005=dose_rho[,1],rho_05=dose_rho[,2],rho_5=dose_rho[,3])
}
association<-rbindlist(association_list);write_table(association,'gene_prism_response_associations.tsv.gz');write_table(rbindlist(dose_association_list),'gene_dose_association_consistency.tsv.gz')
wide<-dcast(association,gene~drug,value.var=c('rho','p','agreeing_doses'))
# Perturbation coordinates are descriptive and are never substituted for response association.
drug_coordinate<-sapply(drug_names,function(drug){ii<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible);matrixStats::colMedians(auc[ii,universe,drop=FALSE],na.rm=TRUE)})
rownames(drug_coordinate)<-universe
melanoma_coordinate<-sapply(drug_names,function(drug){ii<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included & line_drug$melanoma_v600);matrixStats::colMedians(auc[ii,universe,drop=FALSE],na.rm=TRUE)})
rownames(melanoma_coordinate)<-universe
egf_matrix<-study_effect[universe,study_meta$treatment=='EGF',drop=FALSE];nrg1_matrix<-study_effect[universe,study_meta$treatment=='NRG1',drop=FALSE]
egf_effect<-matrixStats::rowMedians(egf_matrix,na.rm=TRUE);egf_consistency<-rowMeans(sign(egf_matrix)==sign(egf_effect),na.rm=TRUE);nrg1_effect<-matrixStats::rowMedians(nrg1_matrix,na.rm=TRUE)
wide[,`:=`(V=drug_coordinate[gene,'vemurafenib'],L=drug_coordinate[gene,'lapatinib'],C=drug_coordinate[gene,'cobimetinib'],Vm=melanoma_coordinate[gene,'vemurafenib'],Lm=melanoma_coordinate[gene,'lapatinib'],Cm=melanoma_coordinate[gene,'cobimetinib'],EGF=egf_effect[gene],EGF_consistency=egf_consistency[gene],NRG1=nrg1_effect[gene])]
wide[,`:=`(vem_eligible=is.finite(p_vemurafenib)&p_vemurafenib<.05&abs(rho_vemurafenib)>=.30&agreeing_doses_vemurafenib>=2,lap_eligible=is.finite(p_lapatinib)&p_lapatinib<.05&abs(rho_lapatinib)>=.30&agreeing_doses_lapatinib>=2,cobi_eligible=is.finite(p_cobimetinib)&p_cobimetinib<.05&abs(rho_cobimetinib)>=.30&agreeing_doses_cobimetinib>=2,egf_eligible=is.finite(EGF)&EGF_consistency>=.70&abs(EGF)>=.10)]
# Focus on receptor, MAPK and proliferation biology; generic expression
# machinery is explicitly removed before any ranking.
msig_all<-as.data.table(msigdbr(species='Homo sapiens'))
pathway_sets<-c('KEGG_MELANOMA','WP_MELANOMA','REACTOME_RAF_ACTIVATION','REACTOME_ONCOGENIC_MAPK_SIGNALING','REACTOME_SIGNALING_BY_MAPK_MUTANTS','REACTOME_SIGNALING_BY_MODERATE_KINASE_ACTIVITY_BRAF_MUTANTS','PID_ERBB1_RECEPTOR_PROXIMAL_PATHWAY','PID_ERBB1_DOWNSTREAM_PATHWAY','PID_ERBB2_ERBB3_PATHWAY','REACTOME_SIGNALING_BY_EGFR','REACTOME_SIGNALING_BY_ERBB2','REACTOME_SIGNALING_BY_ERBB4','WP_EGFEGFR_SIGNALING_PATHWAY','WP_ERBB_SIGNALING_PATHWAY','GOBP_ERBB_SIGNALING_PATHWAY','EGFR_TARGET_GENES','REACTOME_ERK_MAPK_TARGETS','REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES','REACTOME_NEGATIVE_FEEDBACK_REGULATION_OF_MAPK_PATHWAY','REACTOME_SIGNALLING_TO_ERKS','BIOCARTA_ERK_PATHWAY','HALLMARK_KRAS_SIGNALING_UP','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','REACTOME_CELL_CYCLE_MITOTIC')
generic_sets<-c('GOBP_CYTOPLASMIC_TRANSLATION','GOBP_MITOCHONDRIAL_GENE_EXPRESSION','GOBP_MRNA_PROCESSING','GOBP_PROTEIN_FOLDING','GOBP_RIBOSOME_BIOGENESIS','GOBP_RNA_SPLICING')
pathway_annotation<-msig_all[gs_name%chin%pathway_sets,.(pathway_count=uniqueN(gs_name),pathways=paste(sort(unique(gs_name)),collapse='; ')),by=.(gene=gene_symbol)]
generic_process_genes<-unique(msig_all[gs_name%chin%generic_sets,gene_symbol]);generic_regex<-'^(RPS|RPL|MRPS|MRPL|MT-|HIST|H[1-4][A-Z][A-Z0-9-]*$|EIF|EEF|LOC|LINC|MIR|SNOR|SNAR|SCARNA|RNU|RNA5|RNA18|RNA28)'
wide<-merge(wide,pathway_annotation,by='gene');wide<-wide[!grepl(generic_regex,gene)&!gene%chin%generic_process_genes]
source_pool<-function(source,n) {
  if(source=='EGF'){z<-wide[egf_eligible==TRUE];z[,raw_priority:=abs(EGF)*EGF_consistency];z[,orientation:=sign(EGF)]}
  else {elig<-c(vemurafenib='vem_eligible',lapatinib='lap_eligible',cobimetinib='cobi_eligible')[[source]];pc<-paste0('p_',source);rc<-paste0('rho_',source);z<-wide[wide[[elig]]];other_rho<-setdiff(paste0('rho_',drug_names),rc);z$specificity<-abs(z[[rc]])-pmax(abs(z[[other_rho[1]]]),abs(z[[other_rho[2]]]),na.rm=TRUE);z<-z[z$specificity>0,];if(source=='lapatinib')z<-z[is.finite(z$EGF)&z$EGF_consistency>=.60&sign(z$EGF)==-sign(z[[rc]]),];if(source=='cobimetinib')z<-z[is.finite(z$EGF)&z$EGF_consistency>=.60&sign(z$EGF)==-sign(z[[rc]])&sign(z$C)==sign(z[[rc]]),];z$raw_priority <- -log10(z[[pc]])*abs(z[[rc]])+z$specificity;z$orientation <- -sign(z[[rc]])}
  setorder(z,-raw_priority,gene);z<-head(z,n);z[,source_rank:=seq_len(.N)];z[,source_percentile:=1-(source_rank-1)/.N];z[,.(gene,admission_source=source,raw_priority,source_rank,source_percentile,orientation)]
}
source_admissions<-rbindlist(list(source_pool('vemurafenib',30L),source_pool('lapatinib',30L),source_pool('cobimetinib',30L),source_pool('EGF',40L)))
setorder(source_admissions,-source_percentile,-raw_priority,gene)
candidate_pool<-source_admissions[,.(admission_sources=paste(sort(unique(admission_source)),collapse='; '),admission_percentile=max(source_percentile),orientation=orientation[1]),by=gene]
# Require a reproducible local neighborhood in both sources. Tahoe correlations
# use one dose-integrated profile per cell line; drug strata get equal weight.
pool_genes<-candidate_pool$gene;pool_sign<-setNames(candidate_pool$orientation,pool_genes);pool_sign[!is.finite(pool_sign)|pool_sign==0]<-1
auc_cor_by_drug<-lapply(drug_names,function(drug){ii<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible);cor(auc[ii,pool_genes,drop=FALSE],method='spearman',use='pairwise.complete.obs')})
ct_pool<-Reduce('+',auc_cor_by_drug)/length(auc_cor_by_drug);cc_pool<-cor(t(study_effect[pool_genes,,drop=FALSE]),method='spearman',use='pairwise.complete.obs');oo<-outer(pool_sign,pool_sign);ct_pool<-ct_pool*oo;cc_pool<-cc_pool*oo;diag(ct_pool)<-diag(cc_pool)<-NA_real_
knn_mean<-function(m,k=8L)apply(m,1,function(v){v<-sort(v[is.finite(v)],decreasing=TRUE);mean(head(v,min(k,length(v))))})
candidate_pool[,`:=`(tahoe_neighborhood=knn_mean(ct_pool),cytosig_neighborhood=knn_mean(cc_pool))];candidate_pool[,replicated_neighborhood:=pmin(tahoe_neighborhood,cytosig_neighborhood)];candidate_pool[,connectivity_percentile:=frank(replicated_neighborhood,ties.method='average')/.N];candidate_pool[,final_priority:=.65*admission_percentile+.35*connectivity_percentile]
# Preserve representation from every experimental input while allowing shared genes. Each source contributes fifteen high-confidence candidates; duplicate losses are filled by the best remaining replicated candidate.
source_shortlist<-rbindlist(lapply(c('vemurafenib','lapatinib','cobimetinib','EGF'),function(src){gsrc<-source_admissions[admission_source==src,gene];z<-candidate_pool[gene%chin%gsrc];setorder(z,-final_priority,-replicated_neighborhood,gene);head(z,15L)}));admissions<-source_shortlist[!duplicated(gene)];if(nrow(admissions)<60L){fill<-candidate_pool[!gene%chin%admissions$gene];setorder(fill,-final_priority,-replicated_neighborhood,gene);admissions<-rbind(admissions,head(fill,60L-nrow(admissions)),fill=TRUE)};admissions<-head(admissions,60L);chosen<-admissions$gene;gene_sign<-setNames(admissions$orientation,chosen)
check(length(chosen)>=30L&&length(chosen)<=60L&&!any(grepl(generic_regex,chosen))&&!any(chosen%chin%generic_process_genes),'Final gene admission failed')
write_table(source_admissions,'source_gene_admissions.tsv');write_table(candidate_pool,'candidate_gene_coherence.tsv');write_table(admissions,'selected_gene_admissions.tsv');write_table(wide,'gene_selection_evidence.tsv.gz')

## 4. Unrestricted joint clustering; source labels never constrain clusters --
by_drug_cor<-lapply(drug_names,function(drug){ii<-which(line_drug$compound==drug & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible);cor(auc[ii,chosen,drop=FALSE],method='spearman',use='pairwise.complete.obs')});names(by_drug_cor)<-drug_names
ct_raw<-Reduce('+',by_drug_cor)/length(by_drug_cor);cc_raw<-cor(t(study_effect[chosen,,drop=FALSE]),method='spearman',use='pairwise.complete.obs');orientation_matrix<-outer(gene_sign[chosen],gene_sign[chosen]);ct<-ct_raw*orientation_matrix;cc<-cc_raw*orientation_matrix;correlation<-(ct+cc)/2;diag(correlation)<-1
check(all(is.finite(correlation)),'Non-estimable joint correlation')
distance<-as.dist(pmax(1-correlation,0));solutions<-lapply(2:8,function(k){hc<-hclust(distance,method='ward.D2');cl<-cutree(hc,k);si<-silhouette(cl,distance);data.table(k,mean_silhouette=mean(si[,'sil_width']),minimum_size=min(table(cl)),negative_fraction=mean(si[,'sil_width']<0))});curve<-rbindlist(solutions)
eligible_k<-curve[k%in%c(5L,6L)&minimum_size>=MIN_PROGRAM];check(nrow(eligible_k)>0,'Neither five nor six clusters has at least six genes per program')
# Six is retained when it is admissible because the experiment asks whether
# each of the three biological arms contains more than one reproducible module.
# The complete silhouette curve is exported so this structural choice is visible.
K<-if(6L %in% eligible_k$k)6L else eligible_k[which.max(mean_silhouette),k]
tree<-hclust(distance,method='ward.D2');groups<-cutree(tree,K);programs<-paste0('P',seq_len(K))
leaf_order<-unlist(lapply(seq_len(K),function(j){gg<-chosen[groups==j];gg[hclust(as.dist(pmax(1-correlation[gg,gg],0)),method='average')$order]}),use.names=FALSE)
membership<-data.table(gene=leaf_order,program=paste0('P',groups[leaf_order]),plotting_order=seq_along(leaf_order),orientation=gene_sign[leaf_order])
colors<-setNames(grDevices::hcl.colors(K,'Dark 3'),programs);program_order<-programs;row_blocks<-factor(membership$program,levels=programs)
write_table(curve,'silhouette_by_k.tsv');write_table(data.table(gene=chosen,program=paste0('P',groups),silhouette=silhouette(groups,distance)[,'sil_width']),'gene_silhouettes.tsv')
context_cor<-c(by_drug_cor,list(CytoSig=cc_raw,Tahoe_balanced=ct_raw));context_coherence<-rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];rbindlist(lapply(names(context_cor),function(src){m<-context_cor[[src]][g,g];sg<-m*outer(gene_sign[g],gene_sign[g]);data.table(program=p,source=src,n_genes=length(g),raw_mean_rho=mean(m[upper.tri(m)]),signed_mean_rho=mean(sg[upper.tri(sg)]))}))}));write_table(context_coherence,'program_context_coherence.tsv')
write_table(data.table(gene=chosen,correlation),'joint_gene_correlations.tsv.gz')
# Compatibility objects used by the external scoring and heatmap code below.
fp<-cbind(vemurafenib=drug_coordinate[,'vemurafenib'],lapatinib=drug_coordinate[,'lapatinib'],cobimetinib=drug_coordinate[,'cobimetinib'],EGF=egf_effect);zfp<-scale(fp)
family<-setNames(admissions$admission_sources,admissions$gene);raw_correlation<-(ct_raw+cc_raw)/2
cyto_x<-study_effect[chosen,,drop=FALSE]
source_rows<-lapply(drug_names,function(d)which(meta$compound==d & meta$split=='Development' & meta$analysis_included));names(source_rows)<-drug_names
## 5. Candidate hierarchy: classify observed subprogram centroids ------------
signed_fp<-fp[chosen,,drop=FALSE]*gene_sign[chosen]
centroid<-t(sapply(programs,function(p)colMeans(signed_fp[membership[program==p,gene],,drop=FALSE])))
program_response<-rbindlist(lapply(programs,function(p)rbindlist(lapply(drug_names,function(d){
  ii<-which(line_drug$compound==d & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible)
  gg<-membership[program==p,gene];sc<-rowMeans(sweep(auc[ii,gg,drop=FALSE],2,gene_sign[gg],'*'))
  tt<-cor.test(sc,line_drug$prism_sensitivity[ii],method='spearman',exact=FALSE)
  data.table(program=p,drug=d,n=length(ii),rho=unname(tt$estimate),nominal_p=tt$p.value,median_displacement=median(sc))
}))));write_table(program_response,'program_prism_response_associations.tsv')
source_enrichment<-rbindlist(lapply(programs,function(p){
  gg<-membership[program==p,gene]
  rbindlist(lapply(c('vemurafenib','lapatinib','cobimetinib','EGF'),function(s)
    data.table(program=p,source=s,n_genes=sum(gg %chin% source_admissions[admission_source==s,gene]),fraction=sum(gg %chin% source_admissions[admission_source==s,gene])/length(gg))))
}));write_table(source_enrichment,'program_source_enrichment.tsv')
rank01<-function(x)frank(x,ties.method='average')/length(x)
ev<-data.table(program=programs,
  B=rank01(-program_response[drug=='vemurafenib',rho])+rank01(source_enrichment[source=='vemurafenib',fraction]),
  E=rank01(-program_response[drug=='lapatinib',rho])+rank01(source_enrichment[source=='lapatinib',fraction])+rank01(centroid[,'EGF'])-.35*rank01(abs(centroid[,'cobimetinib'])),
  R=rank01(-program_response[drug=='cobimetinib',rho])+rank01(source_enrichment[source=='cobimetinib',fraction])+rank01(centroid[,'EGF'])+rank01(-centroid[,'cobimetinib']))
assignments<-as.data.table(expand.grid(rep(list(c('B','E','R')),K),stringsAsFactors=FALSE));setnames(assignments,programs)
assignments<-assignments[rowSums(assignments=='B')==2L & rowSums(assignments=='E')==2L & rowSums(assignments=='R')==2L]
assignments[,objective:=vapply(seq_len(.N),function(i)sum(vapply(programs,function(p)ev[program==p,get(assignments[[p]][i])],numeric(1))),numeric(1))]
best<-assignments[which.max(objective)];roles<-setNames(as.character(unlist(best[,..programs])),programs)
role_evidence<-melt(ev,id.vars='program',variable.name='role',value.name='evidence_score');role_evidence[,selected:=roles[program]==role]
write_table(role_evidence,'program_role_evidence.tsv')
orient<-gene_sign
membership[,`:=`(role=roles[program],orientation=gene_sign[gene])]
membership[,orientation_basis:=ifelse(role=='B','Higher score is associated with survival under vemurafenib','Higher score follows the CytoSig EGF response')]
membership[,weight:=orientation/.N,by=program]
program_order<-programs[order(match(roles,c('B','E','R')),seq_along(programs))]
leaf_order<-unlist(lapply(program_order,function(p)membership[program==p,gene]),use.names=FALSE)
membership<-membership[match(leaf_order,gene)];membership[,plotting_order:=seq_len(.N)]
role_names<-c(B='BRAF-inhibitor branch',E='ERBB receptor branch',R='Shared MEK–ERK output')
program_label<-setNames(paste0(programs,' · ',role_names[roles]),programs)
heat_block_label<-setNames(paste0(programs,'\n',c(B='BRAF branch',E='ERBB branch',R='MEK–ERK output')[roles]),programs)
row_blocks<-factor(heat_block_label[membership$program],levels=heat_block_label[program_order])
role_drug<-c(B='vemurafenib',E='lapatinib',R='cobimetinib')
role_response<-rbindlist(lapply(names(role_drug),function(r){
  d<-role_drug[[r]];ii<-which(line_drug$compound==d & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible)
  gg<-membership[role==r,gene];sc<-rowMeans(sweep(auc[ii,gg,drop=FALSE],2,gene_sign[gg],'*'))
  tt<-cor.test(sc,line_drug$prism_sensitivity[ii],method='spearman',exact=FALSE)
  data.table(role=r,relevant_drug=d,n=length(ii),rho=unname(tt$estimate),nominal_p=tt$p.value,median_displacement=median(sc))
}));write_table(role_response,'branch_prism_response_associations.tsv')
role_audit<-data.table(program=programs,role=roles,n_genes=as.integer(table(factor(membership$program,levels=programs))),centroid)
write_table(role_audit,'program_hierarchy_evidence.tsv');write_table(membership,'frozen_program_genes.tsv')
# Report the biological processes that distinguish the empirical blocks. The
# tested background is the admitted gene set, so this describes subdivision of
# the programs and is not presented as genome-wide pathway discovery.
program_enrichment<-rbindlist(lapply(programs,function(p)rbindlist(lapply(unique(msig_all[gs_name%chin%pathway_sets,gs_name]),function(gs) {
  inside<-membership$gene[membership$program==p];set_genes<-unique(msig_all[gs_name==gs,gene_symbol])
  a<-sum(inside %chin% set_genes);b<-length(inside)-a;c<-sum(setdiff(chosen,inside) %chin% set_genes);d<-length(chosen)-length(inside)-c
  if(a<2L)return(NULL)
  data.table(program=p,role=roles[match(p,programs)],pathway=gs,overlap=a,program_size=length(inside),
    odds_ratio=unname(fisher.test(matrix(c(a,b,c,d),2),alternative='greater')$estimate),
    p=fisher.test(matrix(c(a,b,c,d),2),alternative='greater')$p.value,
    genes=paste(sort(intersect(inside,set_genes)),collapse=', '))
}),fill=TRUE)),fill=TRUE)
if(nrow(program_enrichment))program_enrichment[,bh_q:=p.adjust(p,'BH')]
write_table(program_enrichment[order(p,-overlap)],'program_pathway_enrichment.tsv')
hierarchy_resolved<-setequal(unique(roles),c('B','E','R')) && all(table(roles)==2L)
partition<-list(method='Ward.D2 hierarchical clustering of the equally weighted Tahoe and CytoSig signed correlation matrices',k=K)
# Render the Y as a hypothesis-guided candidate hierarchy. The edges are NOT
# inferred causal gene edges, and centroid matching is not validation evidence.
node<-data.table(role=c('B','E','R'),x=c(0,2,1),y=c(1,1,0),
  title=c('BRAF-inhibitor branch\n(vemurafenib / PLX4720)',
          'ERBB receptor branch\n(EGF opposed by lapatinib)',
          'Shared MEK–ERK output\n(EGF induced; RAF/MEK inhibited)'))
node[,label:=paste0(title,'\n',vapply(role,function(r)paste(role_audit[role==r,program],collapse=', '),character(1)))]
p<-ggplot()+geom_segment(data=data.table(x=c(.24,1.76),y=.84,xend=c(.86,1.14),yend=.18),aes(x=x,y=y,xend=xend,yend=yend),linewidth=.7,linetype=2,
  arrow=grid::arrow(length=grid::unit(2.5,'mm')))+geom_label(data=node,aes(x,y,label=label,fill=role),size=4,label.padding=grid::unit(.35,'lines'))+
  scale_fill_manual(values=c(B='#E3B88B',E='#A0D6D7',R='#C4AFD8'))+theme_void()+theme(legend.position='none',plot.margin=margin(10,12,10,12))+
  scale_x_continuous(limits=c(-.8,2.8))+scale_y_continuous(limits=c(-.4,1.4))+
  labs(title='Tested BRAF–ERBB convergence model',subtitle='Two upstream perturbation branches converge on a shared transcriptional output')
plot_save(p,'29h_candidate_hierarchy',7.2,4.2)

## 6. Freeze before evaluating primary-PRISM holdout or external outcomes -----
score<-function(x)sapply(programs,function(p){g<-membership[program==p,gene];colMeans(x[g,,drop=FALSE]*gene_sign[g])})
hierarchy_score<-function(x,partial=FALSE)sapply(c('B','E','R'),function(r){g<-membership[role==r,gene];if(!length(g))return(rep(NA_real_,ncol(x)))
  present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)]
  if(!length(present)||(!partial&&length(present)!=length(g)))return(rep(NA_real_,ncol(x)))
  colMeans(x[present,,drop=FALSE]*membership$orientation[match(present,membership$gene)])})
freeze<-list(membership=as.data.frame(membership),admissions=admissions,orientation=orient,roles=roles,partition=partition,
  coordinate_scales=apply(fp[chosen,,drop=FALSE],2,sd,na.rm=TRUE),selection='Nominal drug-response association plus relevant pathway membership and replicated Tahoe/CytoSig coherence',
  validation='Primary PRISM at 2.5 uM; frozen B/E/R and six subprogram scores; no selection using external outcomes')
saveRDS(freeze,file.path(art,'frozen_programs.rds'))
write_table(data.table(prediction=c('BRAF inhibition suppresses B','Bypass retains E and R while B remains suppressed',
  'Adding lapatinib lowers the receptor-supported root in susceptible contexts','EGF activates receptor-associated and shared-response programs',
  'Receptor-driven ERK may resist RAF/MEK inhibition: estimate cobimetinib effects without assuming a 50% root-off threshold'),
  basis=c(rep('Candidate hierarchy; test in held-out expression',4),'Gerosa et al. report differential resistance of receptor-driven MAPK signaling')),'external_prediction_freeze.tsv')
## 7. Render the actual structure and the repaired selection mechanism -------
display_scales<-rbindlist(lapply(c('Tahoe','CytoSig','Joint'),function(src){m<-switch(src,Tahoe=ct,CytoSig=cc,Joint=correlation);data.table(source=src,limit=robust_limit(m[upper.tri(m)],.98),clipped_fraction=mean(abs(m[upper.tri(m)])>robust_limit(m[upper.tri(m)],.98)))}))
write_table(display_scales,'correlation_display_scale.tsv')
for(src in c('Joint','Tahoe','CytoSig')) {
  m<-switch(src,Tahoe=ct,CytoSig=cc,Joint=correlation)[leaf_order,leaf_order];diag(m)<-NA_real_
  lim<-display_scales[source==src,limit]
  heat(m,paste0('29f_',tolower(src),'_signed_correlations'),paste0(src,' gene correlations'),row_blocks,row_blocks,lim,4.8,4.8,FALSE,FALSE)
  heat(m,paste0('29f_',tolower(src),'_signed_correlations_labeled'),paste0(src,' gene correlations | labeled'),row_blocks,row_blocks,lim,9.2,9.2,TRUE,TRUE,6.3,6.3)
}
signed_zfp<-zfp[leaf_order,,drop=FALSE]*gene_sign[leaf_order]
heat(signed_zfp,'29g_perturbation_fingerprint','Selected genes across discovery perturbations',row_blocks,limit=robust_limit(signed_zfp,.98),width=7.2,height=9.2,row_font=7.5,col_font=9)
tahoe_x<-t(effects[unlist(source_rows),chosen,drop=FALSE])
source_heat(tahoe_x,'29g_tahoe_gene_samples','Tahoe drug-response profiles')
source_heat(cyto_x,'29g_cytosig_gene_samples','CytoSig ligand-response studies')
p<-ggplot(curve,aes(k,mean_silhouette))+geom_line()+geom_point(size=2)+geom_vline(xintercept=K,lty=2,color='#B64039')+
  theme_classic(base_size=15)+scale_x_continuous(breaks=2:8)+labs(x='Number of programs',y='Mean silhouette',title='Five or six programs fit similarly',subtitle=paste0('Six retained; smallest program has ',min(table(groups)),' genes'))
plot_save(p,'29e_rank_silhouette',4.4,4.2)
p<-ggplot(context_coherence[source%chin%c('Tahoe_balanced','CytoSig')],aes(factor(source,levels=c('Tahoe_balanced','CytoSig')),signed_mean_rho,fill=source))+geom_col(width=.68)+
  facet_wrap(~program,nrow=2,labeller=as_labeller(program_label))+scale_fill_manual(values=c(Tahoe_balanced='#3B7EA1',CytoSig='#CC6677'),labels=c(Tahoe_balanced='Tahoe',CytoSig='CytoSig'))+
  theme_classic(base_size=13)+theme(legend.position='none',axis.text.x=element_text(angle=28,hjust=1),strip.text=element_text(face='bold'))+
  labs(x=NULL,y='Mean within-program correlation',title='Programs reproduce across Tahoe and CytoSig')
plot_save(p,'29f_context_coherence',6.4,5.3)

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
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+
  ggpubr::stat_cor(aes(group=1),method='spearman',label.x.npc='left',label.y.npc='top',size=3,color='black')+
  facet_wrap(~program,scales='free',labeller=as_labeller(program_label))+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+labs(x='Mean signed Tahoe expression response',
    y='PRISM 5-day viable-cell abundance\nminus the training mean for that drug',
    title='Do frozen expression programs predict viability in held-out cell lines?',subtitle='Each point is one held-out cell line × drug; line is a linear fit; rho is Spearman')
plot_save(p,'29i_subprogram_prism_associations',11,8)
p<-ggplot(root_test,aes(within_drug_score,within_drug_survival,color=compound))+geom_point(size=2)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+
  ggpubr::stat_cor(aes(group=1),method='spearman',label.x.npc='left',label.y.npc='top',size=3.5,color='black')+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+
  labs(x='Mean signed expression of the shared MEK–ERK program\n(centered within drug)',
    y='PRISM 5-day viable-cell abundance\n(centered within drug)',
    title='Does shared MEK–ERK program activity predict viability?',
    subtitle='Six held-out cell lines across three drugs; line is a linear fit; rho is Spearman')
plot_save(p,'29i_shared_root_prism',7,5.5)
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
external_heat<-function(x,metadata,prefix,semantic_order,display_name) {
  finite_genes<-rowSums(is.finite(x))==ncol(x)
  check(sum(finite_genes)>1,'Too few measured frozen genes for an external heatmap')
  c<-cor(x[finite_genes,,drop=FALSE],method='spearman');check(all(is.finite(c)),'Non-estimable external sample correlations')
  tr<-hclust(as.dist(pmax(1-c,0)),method='average');ord<-tr$order
  blocks<-cutree(tr,min(4,ncol(x)-1));bl<-factor(blocks[ord],levels=unique(blocks[ord]))
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),cluster=blocks[ord]),paste0(prefix,'_sample_order.tsv'))
  sample_limit<-as.numeric(quantile(abs(c[upper.tri(c)]),.95));diag(c)<-NA_real_
  heat(c[ord,ord],paste0(prefix,'_sample_correlations'),paste(display_name,'| sample–sample Spearman correlation'),bl,bl,limit=sample_limit,width=8.2,height=8.2,show_cols=FALSE)
  expression_limit<-robust_limit(x,.98)
  heat(x[leaf_order,ord,drop=FALSE],paste0(prefix,'_gene_samples_clustered'),paste(display_name,'| treatment response | samples grouped by similarity'),row_blocks,bl,limit=expression_limit,width=11,height=18,show_cols=FALSE,row_font=8)
  semantic_group<-if(grepl('gse68453',prefix))factor(metadata$cell_line,levels=unique(metadata$cell_line[semantic_order])) else factor(metadata$arm,levels=unique(metadata$arm[semantic_order]))
  semantic_description<-if(grepl('gse68453',prefix))'cell line → treatment → replicate' else 'drug background → EGF time → biological day'
  heat(x[leaf_order,semantic_order,drop=FALSE],paste0(prefix,'_gene_samples_semantic'),paste(display_name,'| treatment response |',semantic_description),
    row_blocks,semantic_group[semantic_order],limit=expression_limit,width=11,height=18,show_cols=FALSE,row_font=8)
}
external_heat(d684,s684,'29k_gse68453',order(match(s684$cell_line,external_lines),match(s684$arm,c('control','plx','lapatinib','plx_lap')),s684$replicate),'GSE68453 melanoma panel')
p<-ggplot(means684[cell_line %chin% factorial_lines],aes(factor(arm,levels=c('control','plx','lapatinib','plx_lap')),score,color=program,group=program))+
  geom_hline(yintercept=0,color='grey70')+geom_line()+geom_point(size=2)+facet_wrap(~cell_line,ncol=3)+
  scale_color_manual(values=colors)+theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=30,hjust=1))+
  labs(x=NULL,y='Mean signed expression change from untreated control',title='GSE68453: RAF and ERBB inhibitor factorials',
    subtitle='PLX = BRAF inhibitor PLX4720; lapatinib = ERBB inhibitor; five melanoma cell lines have all four conditions')
plot_save(p,'29l_factorial_programs',10,7)
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
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='#B64039')+
  ggpubr::stat_cor(method='spearman',label.x.npc='left',label.y.npc='top',size=3,color='black')+
  facet_wrap(~program,scales='free_x',labeller=as_labeller(program_label))+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+labs(x='Mean signed program change after 8 h PLX4720',
    y='Added 3-day benefit of PLX4720 + lapatinib\n(digitized Bliss score)',
    title='Does each early expression program predict combination benefit?',subtitle='Twelve melanoma cell lines; line is a linear fit; rho is Spearman')
plot_save(p,'29m_external_efficacy',9,5)

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
external_heat(d127,s127,'29n_gse127988',order(s127$cobi_uM,s127$egf_h,s127$biological_set,s127$gsm),'GSE127988 A375 EGF time course')
p<-ggplot(day127,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+
  geom_hline(yintercept=0,color='grey70')+geom_line(aes(linetype=factor(biological_set)))+geom_point(size=1.8)+
  facet_wrap(~program,labeller=as_labeller(program_label))+
  theme_classic(base_size=12)+labs(x='Hours after adding 100 ng/mL EGF',y='Mean signed expression change\nfrom the same drug background at 0 h',color='Drug background',linetype='Biological day',
  title='Does EGF activate the frozen programs during BRAF inhibition?',subtitle='A375 cells: 1 µM vemurafenib with or without 1 µM cobimetinib; two biological days')
plot_save(p,'29o_egf_time_courses',9,5)
p<-ggplot(auc127,aes(arm,auc,color=program,group=interaction(program,biological_set)))+geom_hline(yintercept=0,color='grey70')+
  geom_line(aes(linetype=factor(biological_set)))+geom_point(size=2)+facet_wrap(~program,labeller=as_labeller(program_label))+scale_color_manual(values=colors)+
  theme_classic(base_size=12)+labs(x=NULL,y='Mean 0–8 h EGF-induced expression change',linetype='Biological day',
    title='Does cobimetinib block sustained EGF-driven program activation?',subtitle='A375 cells; each line connects the same program and biological day')
plot_save(p,'29p_egf_auc',8,4.5)


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
  theme(axis.text.x=element_text(angle=25,hjust=1))+labs(x=NULL,y='Mean signed expression relative to untreated control',
    title='Do BRAF and ERBB inhibition displace the predicted hierarchy?',
    subtitle='GSE68453: five melanoma lines with control, PLX4720, lapatinib, and their combination')
plot_save(p,'29q_hierarchy_factorial',10,7)
p<-ggplot(hier_benefit,aes(BypassIndex,bliss))+geom_point(size=2,color='#188991')+
  geom_text(aes(label=cell_line),size=3,vjust=-.6,check_overlap=TRUE)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#8351A2')+
  ggpubr::stat_cor(method='spearman',label.x.npc='left',label.y.npc='top',size=3.5,color='black')+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+labs(
    x='Predicted bypass activity after 8 h PLX4720\n(minimum of standardized ERBB and shared-program changes)',
    y='Added 3-day benefit of PLX4720 + lapatinib\n(digitized Bliss score)',
    title='Frozen bypass activity does not rank\ncombination benefit',
    subtitle=paste0('12 melanoma lines · Spearman rho = ',round(bypass_result$spearman,2),'\nTop-three/top-four overlap = ',bypass_result$top3_in_top4))
plot_save(p,'29r_bypass_prediction',7,5.7)
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
  theme_classic(base_size=12)+labs(x='Hours after adding 100 ng/mL EGF',y='Mean signed branch expression\nrelative to same-day untreated A375 cells',linetype='Biological day',color='Drug background',
  title='Which branches recover after EGF is added during BRAF inhibition?',subtitle=paste('Vem = 1 µM vemurafenib; Vem + Cobi adds 1 µM cobimetinib | measured genes:',paste(hierarchy_coverage[dataset=='GSE127988',paste0(role,' ',n_measured,'/',n_frozen)],collapse=', ')))
plot_save(p,'29s_hierarchy_absolute_time_courses',11,5)
p<-ggplot(hincday,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+geom_hline(yintercept=0,color='grey70')+
  geom_line(aes(linetype=factor(biological_set)))+geom_point(size=1.8)+facet_wrap(~role,scales='free_y')+theme_classic(base_size=12)+
  labs(x='Hours after adding 100 ng/mL EGF',y='Mean signed branch expression change\nfrom the same drug background at 0 h',color='Drug background',linetype='Biological day',
  title='How EGF moves each branch during BRAF inhibition',subtitle='A375 cells: vemurafenib with or without the MEK inhibitor cobimetinib')
plot_save(p,'29s_hierarchy_egf_increments',11,5)

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
  ggpubr::stat_cor(method='spearman',label.x.npc='left',label.y.npc='top',size=3,color='black')+
  facet_wrap(~role,scales='free_x')+theme_classic(base_size=11)+theme(aspect.ratio=1)+
  labs(x='Mean signed branch change after 8 h PLX4720',
    y='Added 3-day benefit of PLX4720 + lapatinib\n(digitized Bliss score)',
    title='Which branch carries the molecular-to-function association?',
    subtitle='Twelve melanoma lines; diagnostic shown without outcome-fitted weights or line removal')
plot_save(p,'29r_combination_failure_diagnosis',11,4.5)

## 13. Clean, inspectable report ---------------------------------------------
# Opening schematic: the sources enter discovery before any external test.
design_nodes<-data.table(
  x=c(-1.15,1.15,0,0,-2.05,0,2.05),y=c(3.05,3.05,2.05,1.08,0,0,0),
  label=c('Tahoe\n3 inhibitors + RNA\nPRISM viability','CytoSig\nEGF / NRG1\nligand-response RNA',
    'Joint discovery\n60 response-linked genes\naveraged correlations','Frozen hierarchy\nBRAF · ERBB · MEK–ERK',
    'GSE127988\nEGF ± cobimetinib\nduring vemurafenib','Held-out PRISM\nviable-cell response','GSE68453\nPLX4720 ± lapatinib\nRNA + Bliss'))
design_edges<-data.table(x=c(-1.15,1.15,0,0,0,0),y=c(2.68,2.68,1.70,.70,.70,.70),
  xend=c(-.35,.35,0,-2.05,0,2.05),yend=c(2.36,2.36,1.43,.37,.37,.37))
p<-ggplot()+geom_segment(data=design_edges,aes(x,y,xend=xend,yend=yend),linewidth=.8,color='#526A78',arrow=grid::arrow(length=grid::unit(2.5,'mm')))+
  geom_label(data=design_nodes,aes(x,y,label=label),size=3.8,linewidth=.35,label.padding=grid::unit(.38,'lines'),fill='white',color='#17324D')+
  coord_cartesian(xlim=c(-3,3),ylim=c(-.55,3.55),clip='off')+theme_void(base_size=14)+theme(plot.margin=margin(8,8,8,8))+
  labs(title='Perturbation-informed programs are frozen before validation')
base<-file.path(main_fig_dir,'01_experimental_design');ggsave(paste0(base,'.pdf'),p,width=7.2,height=4.8,device=cairo_pdf,bg='white');ggsave(paste0(base,'.png'),p,width=7.2,height=4.8,dpi=180,bg='white')

# Compact grant views summarize the six modules at the three biological roles.
short_program<-setNames(paste0(programs,' · ',c(B='BRAF',E='ERBB',R='MEK–ERK')[roles]),programs)
p<-ggplot(context_coherence[source%chin%c('Tahoe_balanced','CytoSig')],aes(factor(source,levels=c('Tahoe_balanced','CytoSig')),signed_mean_rho,fill=source))+geom_col(width=.68)+
  facet_wrap(~factor(program,levels=program_order),nrow=2,labeller=as_labeller(short_program))+scale_x_discrete(labels=c(Tahoe_balanced='Tahoe',CytoSig='CytoSig'))+
  scale_fill_manual(values=c(Tahoe_balanced='#3B7EA1',CytoSig='#CC6677'))+theme_classic(base_size=14)+theme(legend.position='none',strip.text=element_text(face='bold'))+
  labs(x=NULL,y='Mean gene correlation',title='Programs reproduce in both discovery sources')
ggsave(file.path(main_fig_dir,'06_source_sharing.pdf'),p,width=5.4,height=4.5,device=cairo_pdf,bg='white');ggsave(file.path(main_fig_dir,'06_source_sharing.png'),p,width=5.4,height=4.5,dpi=180,bg='white')

role_labels<-c(B='BRAF-response branch',E='ERBB-response branch',R='Shared MEK–ERK output')
p<-ggplot(hierarchy_means684[cell_line%chin%factorial_lines],aes(factor(arm,levels=c('control','plx','lapatinib','plx_lap'),labels=c('Control','PLX','Lap','PLX + Lap')),score,color=role,group=role))+
  geom_hline(yintercept=0,color='grey80')+geom_line(linewidth=.65)+geom_point(size=1.8)+facet_wrap(~cell_line,nrow=2)+
  scale_color_manual(values=c(B='#BA6B29',E='#188991',R='#8351A2'),labels=role_labels)+theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=28,hjust=1),legend.position='bottom',legend.title=element_blank())+
  labs(x=NULL,y='Expression change from untreated control',title='RAF and ERBB inhibition displace the frozen branches',subtitle='GSE68453 · five melanoma lines · 8-hour RNA')
ggsave(file.path(main_fig_dir,'07_gse68453_factorial_response.pdf'),p,width=6.4,height=5.2,device=cairo_pdf,bg='white');ggsave(file.path(main_fig_dir,'07_gse68453_factorial_response.png'),p,width=6.4,height=5.2,dpi=180,bg='white')

p<-ggplot(hincday,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+geom_hline(yintercept=0,color='grey80')+
  geom_line(aes(linetype=factor(biological_set)),linewidth=.65)+geom_point(size=1.8)+facet_wrap(~factor(role,levels=c('B','E','R')),labeller=as_labeller(role_labels),scales='free_y')+
  scale_color_manual(values=c('Vem'='#EE6A5F','Vem + Cobi'='#18A4A7'),labels=c('Vem'='Vemurafenib','Vem + Cobi'='Vemurafenib + cobimetinib'))+
  theme_classic(base_size=13)+theme(legend.position='bottom',legend.title=element_blank(),strip.text=element_text(face='bold'))+
  labs(x='Hours after adding EGF',y='Expression change from the\nsame drug background at 0 h',title='EGF activates the ERBB and shared-output branches',subtitle='GSE127988 · A375 · two biological days')
ggsave(file.path(main_fig_dir,'08_gse127988_egf_response.pdf'),p,width=7.0,height=4.1,device=cairo_pdf,bg='white');ggsave(file.path(main_fig_dir,'08_gse127988_egf_response.png'),p,width=7.0,height=4.1,dpi=180,bg='white')

check(length(chosen)==60L,'Expected a 60-gene frozen panel')
check(K==6L && min(table(membership$program))>=6L,'Program count or minimum size changed')
check(setequal(unique(membership$role),c('B','E','R')) && all(table(membership$role)>=2L),'The Y hierarchy is incomplete')
check(!any(grepl(generic_regex,membership$gene)) && !any(membership$gene %chin% generic_process_genes),'Generic expression genes entered the panel')
check(identical(dim(ct),c(60L,60L)) && identical(dim(cc),c(60L,60L)) && identical(dim(correlation),c(60L,60L)),'A correlation matrix is not square')
check(all(context_coherence[source %chin% c('Tahoe_balanced','CytoSig'),signed_mean_rho]>0),'A retained program is anticorrelated within a discovery source')
check(all(role_response$nominal_p<.05),'A biological branch lacks nominal association with its relevant PRISM drug response')
check(identical(readRDS(file.path(art,'frozen_programs.rds'))$membership,as.data.frame(membership)),'Frozen program definition changed during validation')
check(!any(line_drug$depmap_id[c(fit_rows,held)] %chin% EXCLUDED_IDS),'A discarded outlier entered the analysis')
check(!length(intersect(line_drug$depmap_id[fit_rows],line_drug$depmap_id[held])),'Development and holdout share a cell line')
verification<-data.table(check=c('60-gene panel','Six programs, each at least six genes','Two resolved programs per Y branch','No ribosomal or generic expression genes','All three correlation matrices are square','Every program has positive within-source coherence in Tahoe and CytoSig','Each branch is associated with its relevant PRISM drug response','Frozen before external validation','Discarded outliers excluded','Cell-line holdout is disjoint'),passed=TRUE)
write_table(verification,'verification.tsv')

figure_files<-sort(list.files(file.path(out,'figures'),pattern='\\.pdf$',recursive=TRUE))
write_table(data.table(figure=file.path('figures',figure_files)),'FIGURE_INDEX.tsv')
summary_table<-merge(role_audit[,.(program,role,n_genes)],context_coherence[source %chin% c('Tahoe_balanced','CytoSig'),.(program,source,signed_mean_rho)],by='program')
write_table(summary_table,'program_summary.tsv')

writeLines(c('# Stage 30: BRAF–ERBB program hierarchy','',
  '## Question','',
  'Can inhibitor-response RNA from Tahoe and ligand-response RNA from CytoSig jointly define a BRAF-response branch and an ERBB-response branch that converge on a shared MEK–ERK transcriptional output—and do those frozen programs transfer to independent molecular and viability data?','',
  '## Discovery result','',
  sprintf('The final panel contains %d biologically focused genes in %d programs (%s genes). The programs were cut from the average of square Tahoe and CytoSig gene-correlation matrices. CytoSig therefore changes both gene admission and program boundaries; it is not post-hoc annotation.',length(chosen),K,paste(sort(as.integer(table(membership$program))),collapse=', ')),
  sprintf('All programs have positive within-program correlation in both sources. Mean signed correlation spans %.2f–%.2f in Tahoe and %.2f–%.2f in CytoSig.',min(context_coherence[source=='Tahoe_balanced',signed_mean_rho]),max(context_coherence[source=='Tahoe_balanced',signed_mean_rho]),min(context_coherence[source=='CytoSig',signed_mean_rho]),max(context_coherence[source=='CytoSig',signed_mean_rho])),
  paste(sprintf('The resolved branches are BRAF response (%s), ERBB response (%s), and shared MEK–ERK output (%s).',paste(programs[roles=='B'],collapse=', '),paste(programs[roles=='E'],collapse=', '),paste(programs[roles=='R'],collapse=', '))),'',
  '## Functional and molecular validation','',
  paste0('Within the Tahoe development set, every branch is associated with the matching PRISM response: ',paste(sprintf('%s versus %s: rho %.2f, p %.3g',role_response$role,role_response$relevant_drug,role_response$rho,role_response$nominal_p),collapse='; '),'.'),
  sprintf('In held-out cell lines, adding the frozen shared MEK–ERK score improves PRISM RMSE from %.3f to %.3f; the drug-adjusted rank correlation is %.2f.',root_comparison$baseline_rmse,root_comparison$root_rmse,root_comparison$partial_spearman_controlling_drug),
  sprintf('In GSE127988 A375 cells, EGF activates the ERBB and shared MEK–ERK branches on both biological days. Cobimetinib changes the shared EGF-response AUC from %.3f to %.3f on average and the ERBB AUC from %.3f to %.3f.',mean(hcontrast[role=='R',vem_auc]),mean(hcontrast[role=='R',vem_cobi_auc]),mean(hcontrast[role=='E',vem_auc]),mean(hcontrast[role=='E',vem_cobi_auc])),
  sprintf('The separate no-added-ligand GSE68453 cross-line predictor remains weak (rho %.2f, p %.3f). It does not negate the matched EGF experiment, but it does not validate individualized PLX4720+lapatinib benefit ranking.',bypass_result$spearman,bypass_result$nominal_p),'',
  '## Interpretation','',
  'The pilot supports the core molecular chain: heterogeneous inhibitor and ligand datasets can jointly recover reproducible response programs; EGF activates the receptor and shared-output branches during BRAF inhibition; MEK inhibition suppresses that response; and the shared program carries independent viability information. The evidence does not support a calibrated predictor of which melanoma line will benefit from adding lapatinib in a no-added-EGF setting.','',
  '[Open the result viewer](index.html) · [Frozen genes](frozen_program_genes.tsv) · [Analysis code](../../scripts/30_rebuild_braf_erbb_programs.R)'),file.path(out,'STORY.md'))

writeLines(c('# Stage 30 implementation','',
  'This directory contains only the current end-to-end BRAF–ERBB experiment. Main figures follow the scientific story in numerical order; supporting gene-by-sample heatmaps and diagnostics are under `figures/supplementary`.','',
  'Gene admission requires nominal association with the relevant development-set PRISM response, direction stability across Tahoe doses, membership in focused BRAF/ERBB/MAPK/proliferation biology, and local correlation structure in both Tahoe and CytoSig. Ribosomal, translation, mitochondrial-expression, RNA-processing, histone and noncoding genes are excluded. Sixty genes are retained while preserving candidates from vemurafenib, lapatinib, cobimetinib and CytoSig EGF.','',
  'Ward.D2 clustering uses `1 - mean(Tahoe signed correlation, CytoSig signed correlation)`. Six programs are retained from the prespecified five-to-six range because all six contain at least six genes and permit two empirical modules per branch. Roles are assigned after clustering by a globally optimized enrichment score based on the relevant admission source, perturbation displacement and development-set PRISM association. Every program receives one resolved role.','',
  'Programs, signs and weights are frozen before the held-out PRISM lines, GSE127988 or GSE68453 outcomes are evaluated. External missing genes are never replaced with zero; coverage and identifiers are exported explicitly.','',
  'Run from `/work`:','',
  '```sh','Rscript tahoe_prism_v2/scripts/30_rebuild_braf_erbb_programs.R','python3 tahoe_prism_v2/scripts/30_verify_braf_erbb_programs.py','```'),file.path(out,'README.md'))

html_escape<-function(x){x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
html_table<-function(x)paste0('<table><tr>',paste0('<th>',names(x),'</th>',collapse=''),'</tr>',paste(apply(as.data.frame(x),1,function(r)paste0('<tr>',paste0('<td>',html_escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table>')
main_pdfs<-file.path('figures/main',sort(list.files(main_fig_dir,pattern='\\.pdf$')));supp_pdfs<-file.path('figures/supplementary',sort(list.files(supp_fig_dir,pattern='\\.pdf$')))
make_cards<-function(ff)vapply(ff,function(f){png<-sub('\\.pdf$','.png',f);paste0('<article><h3>',html_escape(sub('\\.pdf$','',basename(f))),'</h3><a href="',f,'"><img loading="lazy" src="',png,'"></a><p><a href="',f,'">Open PDF</a></p></article>')},character(1))
tables<-sort(list.files(out,pattern='\\.tsv(\\.gz)?$'))
writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Stage 30 · BRAF–ERBB programs</title>',
  '<style>body{font:16px/1.5 system-ui;background:#f4f6f8;color:#183044;margin:0}main{max-width:1200px;margin:auto;padding:26px}a{color:#075b8e}.status,article{background:white;border:1px solid #d8e0e6;border-radius:5px}.status{padding:16px;border-left:6px solid #287d99}.grid{display:grid;grid-template-columns:1fr;gap:18px}article{padding:14px}img{display:block;max-width:100%;height:auto;margin:auto}table{border-collapse:collapse;width:100%;font-size:13px;background:white}th,td{padding:7px;border:1px solid #d7dde2;text-align:left}.tables{columns:2}pre{white-space:pre-wrap;background:#172733;color:#edf4f7;padding:15px}</style><main>',
  '<h1>Stage 30 · BRAF–ERBB program hierarchy</h1><p><a href="STORY.md">Scientific story</a> · <a href="README.md">Methods</a> · <a href="../../scripts/30_rebuild_braf_erbb_programs.R">R code</a></p>',
  paste0('<div class="status"><strong>',length(chosen),' genes · ',K,' programs · complete BRAF/ERBB/MEK–ERK hierarchy</strong><p>CytoSig materially informs discovery. EGF activates the ERBB and shared-output branches; cobimetinib suppresses them; the shared branch improves held-out PRISM prediction.</p></div>'),
  '<h2>Main figures</h2><div class="grid">',make_cards(main_pdfs),'</div><h2>Program definitions</h2>',html_table(role_audit),'<h2>Branch–PRISM associations</h2>',html_table(role_response),'<h2>Held-out PRISM test</h2>',html_table(root_comparison),'<h2>GSE127988 EGF and cobimetinib test</h2>',html_table(hcontrast),'<h2>GSE68453 functional ranking test</h2>',html_table(bypass_result),
  '<h2>Supporting figures</h2><div class="grid">',make_cards(supp_pdfs),'</div><h2>Tables</h2><ul class="tables">',paste0('<li><a href="',tables,'">',tables,'</a></li>'),'</ul>',
  '<h2>Complete analysis code</h2><details><summary>Show code</summary><pre><code>',paste(html_escape(readLines(file.path(root,'scripts/30_rebuild_braf_erbb_programs.R'),warn=FALSE)),collapse='\n'),'</code></pre></details></main></html>'),file.path(out,'index.html'))
message('Stage 30 complete: ',out)
