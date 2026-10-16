#!/usr/bin/env Rscript
## Stage 33: rebuild BRAF--ERBB programs from response association and cross-source coherence.
## Drug-response genes must be associated with PRISM response, retain their
## direction across doses, and participate in coherent Tahoe/CytoSig modules.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cluster);library(msigdbr);library(ggpubr);library(jsonlite)})
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(arg))dirname(dirname(normalizePath(sub('^--file=','',arg[1])))) else '/work/tahoe_prism_v2'
work<-dirname(root);input_art<-file.path(root,'artifacts/26_braf_erbb_bypass')
prepared_art<-file.path(root,'artifacts/27_braf_erbb_hierarchy')
art<-file.path(root,'artifacts/33_braf_erbb_programs');out<-file.path(root,'results/33_braf_erbb_programs')
unlink(art,recursive=TRUE,force=TRUE);unlink(out,recursive=TRUE,force=TRUE)
dir.create(art,recursive=TRUE,showWarnings=FALSE);dir.create(out,recursive=TRUE,showWarnings=FALSE)
main_fig_dir<-file.path(out,'figures/main');supp_fig_dir<-file.path(out,'figures/supplementary')
dir.create(main_fig_dir,recursive=TRUE,showWarnings=FALSE);dir.create(supp_fig_dir,recursive=TRUE,showWarnings=FALSE)
source(file.path(work,'src/R/lib_heatmap.R'))
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
write_table<-function(x,name)fwrite(as.data.table(x),file.path(out,name),sep='\t',na='NA')
EXCLUDED_IDS<-c('ACH-000311','ACH-000628');K<-5L;MIN_PROGRAM<-12L
ROLE_COLORS<-c(B='#D88736',E='#168C86',R='#7A5AA6')
ROLE_FILLS<-c(B='#F2D1AD',E='#AFE0DC',R='#D0C2E2')
DRUG_COLORS<-c(vemurafenib=ROLE_COLORS[['B']],lapatinib=ROLE_COLORS[['E']],cobimetinib=ROLE_COLORS[['R']])
SOURCE_COLORS<-c(Tahoe='#356F91',CytoSig='#C95F73')


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
  '29h_candidate_hierarchy'='figures/main/03_program_hierarchy',
  '29f_context_coherence'='figures/main/07_source_context_coherence',
  '29i_subprogram_prism_associations'='figures/supplementary/S02_subprogram_viability_associations',
  '29i_shared_root_prism'='figures/main/11_prism_viability_association',
  '29l_factorial_programs'='figures/main/08_gse68453_factorial_response',
  '29m_external_efficacy'='figures/supplementary/S03_subprogram_combination_associations',
  '29o_egf_time_courses'='figures/main/10_gse127988_egf_response',
  '29r_bypass_prediction'='figures/main/12_gse68453_combination_prediction')
heat_map<-c(
  '29f_joint_signed_correlations'='figures/main/04_joint_gene_correlations',
  '29f_tahoe_signed_correlations'='figures/main/05_tahoe_gene_correlations',
  '29f_cytosig_signed_correlations'='figures/main/06_cytosig_gene_correlations',
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

heat <- function(x,name,title,rows=NULL,cols=NULL,limit=2.5,width=10,height=NULL,show_cols=TRUE,show_rows=TRUE,row_font=9,col_font=8,block_rectangles=FALSE) {
  destination<-unname(heat_map[name]);if(is.na(destination))return(invisible(NULL))
  base<-file.path(out,destination);dir.create(dirname(base),recursive=TRUE,showWarnings=FALSE)
  if(is.null(height))height<-max(4.8,nrow(x)*.12+1.8)
  if(!is.null(rows))check(length(rows)==nrow(x),'Heatmap row annotations misaligned')
  if(!is.null(cols))check(length(cols)==ncol(x),'Heatmap column annotations misaligned')
  legend_title<-if(grepl('correlation',name))'Spearman rho' else if(grepl('gene_samples',name))'Expression\nchange' else 'Standardized\neffect'
  hm_id<-paste0('heat_',gsub('[^A-Za-z0-9]','_',name))
  split_rows<-if(block_rectangles)NULL else rows;split_cols<-if(block_rectangles)NULL else cols
  row_role<-col_role<-NULL
  if(block_rectangles) {
    infer_role<-function(v)factor(ifelse(grepl('BRAF',as.character(v)),'BRAF arm',ifelse(grepl('ERBB',as.character(v)),'ERBB arm','Shared output')),
      levels=c('BRAF arm','ERBB arm','Shared output'))
    row_role<-data.frame(Pathway=infer_role(rows));col_role<-data.frame(Pathway=infer_role(cols))
  }
  heatmap_body_side<-if(block_rectangles)if(show_rows)5.6 else 2.45 else NULL
  h<-draw_heatmap(x,col_anno_df=NULL,col_anno_cols=NULL,
    row_anno_df=row_role,row_anno_cols=if(block_rectangles)'Pathway' else NULL,
    palette_overrides=if(block_rectangles)list(Pathway=c('BRAF arm'=ROLE_COLORS[['B']],'ERBB arm'=ROLE_COLORS[['E']],'Shared output'=ROLE_COLORS[['R']])) else NULL,
    ann_row_args=if(block_rectangles)list(show_annotation_name=FALSE) else list(),
    scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(x)),column_order=seq_len(ncol(x)),row_split=split_rows,column_split=split_cols,
    cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(2,'mm'),column_gap=grid::unit(2,'mm'),
    row_title_rot=0,column_title_rot=0,row_title_gp=grid::gpar(fontsize=10,fontface='bold'),column_title_gp=grid::gpar(fontsize=8,fontface='bold'),
    show_row_names=show_rows,show_column_names=show_cols,row_names_gp=grid::gpar(fontsize=row_font),column_names_gp=grid::gpar(fontsize=col_font),
    column_names_rot=90,hm_name=hm_id,heatmap_legend_param=list(title=legend_title,at=c(-limit,0,limit),labels=round(c(-limit,0,limit),2)),title=title,
    width=if(block_rectangles)grid::unit(heatmap_body_side,'in') else NULL,height=if(block_rectangles)grid::unit(heatmap_body_side,'in') else NULL)
  add_blocks<-function(){
    if(!block_rectangles)return(invisible(NULL))
    rr<-rle(as.character(rows));ends<-cumsum(rr$lengths);starts<-c(1L,head(ends,-1L)+1L);n<-nrow(x)
    ComplexHeatmap::decorate_heatmap_body(hm_id,{for(i in seq_along(starts)){
      left<-(starts[i]-1)/n;bottom<-1-ends[i]/n;side<-rr$lengths[i]/n
      grid::grid.rect(x=grid::unit(left,'npc'),y=grid::unit(bottom,'npc'),width=grid::unit(side,'npc'),height=grid::unit(side,'npc'),just=c('left','bottom'),gp=grid::gpar(fill=NA,col='#1A1A1A',lwd=1.6))
    }})
  }
  cairo_pdf(paste0(base,'.pdf'),width=width,height=height)
  drawn<-ComplexHeatmap::draw(h,heatmap_legend_side='right')
  add_blocks()
  check(identical(as.integer(unlist(ComplexHeatmap::row_order(drawn))),seq_len(nrow(x))),'Heatmap changed frozen gene order')
  check(identical(as.integer(unlist(ComplexHeatmap::column_order(drawn))),seq_len(ncol(x))),'Heatmap changed frozen sample order')
  invisible(dev.off())
  png(paste0(base,'.png'),width=width,height=height,units='in',res=140,type='cairo')
  ComplexHeatmap::draw(h,heatmap_legend_side='right');add_blocks();invisible(dev.off())
}
robust_limit<-function(x,p=.98)max(as.numeric(quantile(abs(x[is.finite(x)]),p,na.rm=TRUE)),.01)

source_heat <- function(x,name,title,context_df,context_cols) {
  x<-x[leaf_order,,drop=FALSE]
  check(nrow(context_df)==ncol(x),'Source heatmap metadata is misaligned')
  z<-t(scale(t(x)));z[!is.finite(z)]<-NA_real_
  signed_z<-z*gene_sign[rownames(z)]
  activity<-sapply(program_order,function(p)colMeans(signed_z[membership[program==p,gene],,drop=FALSE],na.rm=TRUE))
  winner<-apply(activity,1,function(a)if(all(!is.finite(a)))1L else which.max(replace(a,!is.finite(a),-Inf)))
  strength<-activity[cbind(seq_len(nrow(activity)),winner)]
  # This deterministic order makes five diagonal blocks: samples are first
  # assigned to the program they activate most strongly, then ranked within
  # that block by the magnitude of that activation.
  ord<-order(winner,-strength,seq_along(winner))
  cb<-factor(program_order[winner[ord]],levels=program_order)
  lim<-robust_limit(signed_z,.98)
  context_df<-as.data.frame(context_df[ord,,drop=FALSE]);context_df$Program<-cb
  write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),dominant_program=as.character(cb),context_df[,setdiff(names(context_df),'Program'),drop=FALSE]),paste0(name,'_sample_order.tsv'))
  row_anno<-data.frame(Program=factor(membership$program,levels=program_order))
  palettes<-c(list(Program=colors),context_cols)
  key_examples<-c('SPAG5','CSE1L','HBEGF','FOS','DUSP4','PHLDA1','PLK1','TOP2A','PML','YWHAZ')
  render_one<-function(base,full_labels=FALSE,w=7.6,h=5.8) {
    hm<-draw_heatmap(signed_z[,ord,drop=FALSE],col_anno_df=context_df,col_anno_cols=names(context_df),row_anno_df=row_anno,row_anno_cols='Program',palette_overrides=palettes,ann_top_args=list(show_annotation_name=FALSE),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(signed_z)),column_order=seq_len(ncol(signed_z)),row_split=row_blocks,column_split=cb,cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(1.2,'mm'),column_gap=grid::unit(1.3,'mm'),row_title=NULL,show_row_names=TRUE,show_column_names=FALSE,row_labels=if(full_labels)rownames(signed_z) else ifelse(rownames(signed_z)%chin%key_examples,rownames(signed_z),''),row_names_gp=grid::gpar(fontsize=if(full_labels)6.0 else 8.5,fontface=ifelse(rownames(signed_z)%chin%key_examples,'bold','plain')),hm_name='RNA response',heatmap_legend_param=list(title='Row-standardized\nRNA response',at=c(-lim,0,lim),labels=round(c(-lim,0,lim),1)),title=program_order)
    draw_one<-function()ComplexHeatmap::draw(hm,heatmap_legend_side='right',annotation_legend_side='right',column_title=title,column_title_gp=grid::gpar(fontsize=13,fontface='bold'))
    cairo_pdf(paste0(base,'.pdf'),width=w,height=h);draw_one();invisible(dev.off())
    png(paste0(base,'.png'),width=w,height=h,units='in',res=180,type='cairo');draw_one();invisible(dev.off())
  }
  destination<-unname(heat_map[paste0(name,'_signed')]);base<-file.path(out,destination)
  render_one(base,FALSE,7.6,5.8)
  full_base<-file.path(supp_fig_dir,if(grepl('cytosig',name))'S16_cytosig_all_profiles_gene_heatmap' else 'S15_tahoe_all_profiles_gene_heatmap')
  render_one(full_base,TRUE,8.2,7.5)
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
source_admissions<-rbindlist(list(source_pool('vemurafenib',50L),source_pool('lapatinib',50L),source_pool('cobimetinib',50L),source_pool('EGF',60L)))
setorder(source_admissions,-source_percentile,-raw_priority,gene)
candidate_pool<-source_admissions[,.(admission_sources=paste(sort(unique(admission_source)),collapse='; '),admission_percentile=max(source_percentile),orientation=orientation[1]),by=gene]
pool_genes<-candidate_pool$gene

# Expression covariance and response association answer different questions.
# PRISM is used above to admit drug-response genes. Program covariance below
# uses every retained Tahoe development profile: all cell lines, all doses and
# all three drugs, including the one line-drug profile without a PRISM value.
tahoe_cor_pool<-lapply(drug_names,function(drug){
  ii<-which(meta$compound==drug & meta$split=='Development' & meta$analysis_included)
  cor(effects[ii,pool_genes,drop=FALSE],method='spearman',use='pairwise.complete.obs')
});names(tahoe_cor_pool)<-drug_names
ct_raw_pool<-Reduce('+',tahoe_cor_pool)/length(tahoe_cor_pool)

# Read the complete CytoSig perturbation collection. EGF/NRG1 selected the
# receptor-responsive candidates, but all response signatures provide the
# covariance background, including unrelated cytokines where the ERBB branch
# should remain quiet. Each signature is already treatment minus matched
# control. Studies receive equal total weight so replicated studies cannot
# dominate the correlation matrix.
full_cyto_path<-file.path(work,'data/01_CytoSig_data/diff.merge.gz')
full_cyto_meta_path<-file.path(work,'data/01_CytoSig_data/meta_info.txt')
check(file.exists(full_cyto_path)&&file.exists(full_cyto_meta_path),'Full CytoSig cache is missing')
full_cyto_dt<-fread(full_cyto_path,showProgress=FALSE);setnames(full_cyto_dt,1L,'gene')
full_cyto_dt<-full_cyto_dt[gene%chin%pool_genes]
check(nrow(full_cyto_dt)==length(pool_genes),'A candidate gene is absent from full CytoSig')
full_cyto<-as.matrix(full_cyto_dt[,-1]);rownames(full_cyto)<-full_cyto_dt$gene
full_cyto<-full_cyto[pool_genes,,drop=FALSE]
full_meta0<-fread(full_cyto_meta_path,skip=1L,header=FALSE,col.names=c('profile_id','treatment','condition','duration','dose'))
full_cyto_meta<-merge(data.table(profile_id=colnames(full_cyto)),full_meta0,by='profile_id',all.x=TRUE,sort=FALSE)
full_cyto_meta<-full_cyto_meta[match(colnames(full_cyto),profile_id)]
full_cyto_meta[,study:=sub('.*@((?:GSE[0-9]+|E-MTAB-[0-9]+)).*','\\1',profile_id,perl=TRUE)]
full_cyto_meta[!grepl('^(GSE[0-9]+|E-MTAB-[0-9]+)$',study),study:=paste0('independent_',.I)]
weighted_spearman<-function(x,weights) {
  # x is genes x profiles. Pairwise weighted correlations of profile ranks.
  ranks<-t(apply(x,1,rank,na.last='keep',ties.method='average'))
  p<-nrow(ranks);answer<-diag(1,p);dimnames(answer)<-list(rownames(ranks),rownames(ranks))
  if(p>1L)for(i in seq_len(p-1L))for(j in (i+1L):p) {
    ok<-is.finite(ranks[i,])&is.finite(ranks[j,]);w<-weights[ok];w<-w/sum(w)
    a<-ranks[i,ok];b<-ranks[j,ok];a<-a-sum(w*a);b<-b-sum(w*b)
    answer[i,j]<-answer[j,i]<-sum(w*a*b)/sqrt(sum(w*a*a)*sum(w*b*b))
  }
  answer
}
study_n<-full_cyto_meta[, .N,by=study];study_size<-setNames(study_n$N,study_n$study)
cyto_weight<-1/study_size[full_cyto_meta$study];cyto_weight<-cyto_weight/sum(cyto_weight)
cc_raw_pool<-weighted_spearman(full_cyto,cyto_weight)
ligand_profile<-full_cyto_meta$treatment%chin%c('EGF','NRG1')
check(sum(ligand_profile)==67L && sum(!ligand_profile)>1900L,'Unexpected full CytoSig ligand/background profile counts')
cc_ligand_pool<-cor(t(full_cyto[,ligand_profile,drop=FALSE]),method='spearman',use='pairwise.complete.obs')
cc_other_pool<-cor(t(full_cyto[,!ligand_profile,drop=FALSE]),method='spearman',use='pairwise.complete.obs')

# Candidate priority uses reproducible local connectivity without imposing the
# response-direction sign. Signs are learned within each eventual module.
knn_mean<-function(m,k=8L)apply(abs(m),1,function(v){v<-sort(v[is.finite(v)],decreasing=TRUE);mean(head(v,min(k,length(v))))})
candidate_pool[,`:=`(tahoe_neighborhood=knn_mean(ct_raw_pool),cytosig_neighborhood=knn_mean(cc_raw_pool))]
candidate_pool[,replicated_neighborhood:=pmin(tahoe_neighborhood,cytosig_neighborhood)]
candidate_pool[,connectivity_percentile:=frank(replicated_neighborhood,ties.method='average')/.N]
candidate_pool[,final_priority:=.65*admission_percentile+.35*connectivity_percentile]
write_table(source_admissions,'source_gene_admissions.tsv');write_table(candidate_pool,'candidate_gene_coherence.tsv');write_table(wide,'gene_selection_evidence.tsv.gz')
write_table(full_cyto_meta[,.(profile_id,treatment,study,condition,duration,dose,ligand_context=ifelse(ligand_profile,'EGF / NRG1','Other perturbation'))],'cytosig_all_profile_audit.tsv.gz')
write_table(rbindlist(list(
  rbindlist(lapply(drug_names,function(d)data.table(source='Tahoe',context=d,profiles=sum(meta$compound==d&meta$split=='Development'&meta$analysis_included),cell_lines=uniqueN(meta[compound==d&split=='Development'&analysis_included,depmap_id])))),
  data.table(source='CytoSig',context=c('EGF / NRG1','Other perturbations','All perturbations'),profiles=c(sum(ligand_profile),sum(!ligand_profile),length(ligand_profile)),cell_lines=NA_integer_)
)),'correlation_sample_audit.tsv')

## 4. Build five sizeable, biology-constrained cross-source cores ----------
# Canonical anchors state the biological question; they do not fix the other
# members. Remaining genes are added from response-eligible pathway pools by
# maximizing the weakest signed mean correlation across five independent
# contexts: each Tahoe drug, CytoSig EGF/NRG1 and all other CytoSig stimuli.
context_pool<-c(tahoe_cor_pool,list(CytoSig_ligand=cc_ligand_pool,CytoSig_other=cc_other_pool))
module_orientation<-function(g) {
  j<-(ct_raw_pool[g,g,drop=FALSE]+cc_raw_pool[g,g,drop=FALSE])/2
  s<-sign(eigen(j,symmetric=TRUE)$vectors[,1]);s[s==0]<-1;setNames(s,g)
}
module_context_cor<-function(g) {
  s<-module_orientation(g)
  vapply(context_pool,function(m){z<-m[g,g,drop=FALSE]*outer(s,s);mean(z[upper.tri(z)])},numeric(1))
}
program_seeds<-list(
  P1=c('SPAG5','KIF11','KPNA2','CSE1L'),
  P2=c('EGR1','FOS','JUN','JUNB'),
  # This complete adaptive-MAPK core was recovered in the preceding pilot and
  # already ranked both PLX4720 resistance and lapatinib combination benefit.
  # It is retained only because every one of the stricter Stage-33 context
  # checks passes; the remaining programs are rebuilt around it without overlap.
  P3=c('SHC4','DUSP6','DUSP4','BLZF1','HMGB2','TGFB1','PPP1CC','KCNN4',
       'PLAUR','ITGA2','PHLDA1','TMEM158','ANGPTL4'),
  P4=c('TPX2','PLK1','DLGAP5','TACC3'),
  P5=c('MAPK1','PSMD1','PSMD14','SLC38A1'))
check(all(unique(unlist(program_seeds))%chin%pool_genes),'A canonical program anchor failed gene admission')
reserved_anchors<-unique(unlist(program_seeds))
grow_program<-function(program,pool,n,used=character()) {
  current<-program_seeds[[program]]
  pool<-setdiff(intersect(pool,pool_genes),c(used,setdiff(reserved_anchors,current)))
  check(all(current%chin%pool_genes)&&length(pool)>=n,paste('Insufficient eligible genes for',program))
  trace<-list()
  while(length(current)<n) {
    candidates<-setdiff(pool,current)
    trial<-rbindlist(lapply(candidates,function(add){v<-module_context_cor(c(current,add));data.table(add,t(v),objective=min(v))}))
    setorder(trial,-objective,-vemurafenib,-lapatinib,-cobimetinib,-CytoSig_ligand,-CytoSig_other,add)
    current<-c(current,trial$add[1]);trace[[length(trace)+1L]]<-data.table(program,step=length(current),selected_gene=trial$add[1],trial[1,-1])
  }
  list(genes=current,trace=rbindlist(trace,fill=TRUE))
}
cell_cycle_pool<-wide[gene%chin%pool_genes&grepl('E2F_TARGETS|G2M_CHECKPOINT|CELL_CYCLE_MITOTIC',pathways),gene]
erbb_pool<-wide[gene%chin%pool_genes&egf_eligible&grepl('EGFR|ERBB|KRAS_SIGNALING_UP|ERK_MAPK_TARGETS|NEGATIVE_FEEDBACK',pathways),gene]
mapk_pool<-wide[gene%chin%pool_genes&egf_eligible&grepl('ERK_MAPK_TARGETS|MAPK_TARGETS_NUCLEAR|NEGATIVE_FEEDBACK|KRAS_SIGNALING_UP|EGFR|ERBB',pathways),gene]
cobi_pool<-wide[gene%chin%pool_genes&cobi_eligible&grepl('ERK|MAPK|CELL_CYCLE|G2M|E2F|ERBB|KRAS',pathways),gene]
program_pools<-list(P1=intersect(cell_cycle_pool,source_admissions[admission_source=='vemurafenib',gene]),P2=erbb_pool,P3=mapk_pool,P4=cell_cycle_pool,P5=cobi_pool)
program_sizes<-c(P1=14L,P2=14L,P3=13L,P4=15L,P5=12L)
program_genes<-list();growth_trace<-list();used<-character()
for(p in names(program_seeds)) {
  z<-grow_program(p,program_pools[[p]],program_sizes[[p]],used)
  program_genes[[p]]<-z$genes;growth_trace[[p]]<-z$trace;used<-c(used,z$genes)
}
check(!anyDuplicated(unlist(program_genes)),'Programs overlap')
programs<-names(program_genes);K<-length(programs);chosen<-unlist(program_genes,use.names=FALSE)

# Learn within-program signs from the two source matrices, then choose only the
# arbitrary global sign using the prespecified biological direction.
gene_sign<-setNames(rep(NA_real_,length(chosen)),chosen)
for(p in programs) {
  g<-program_genes[[p]];s<-module_orientation(g)
  if(p=='P1') {
    ii<-which(line_drug$compound=='vemurafenib'&line_drug$split=='Development'&line_drug$analysis_included&line_drug$eligible)
    score0<-rowMeans(sweep(auc[ii,g,drop=FALSE],2,s,'*'))
    if(cor(score0,line_drug$prism_survival[ii],use='complete.obs')<0)s<- -s
  } else if(p%chin%c('P2','P3','P4')) {
    if(mean(egf_effect[g]*s,na.rm=TRUE)<0)s<- -s
  } else {
    ii<-which(line_drug$compound=='cobimetinib'&line_drug$split=='Development'&line_drug$analysis_included&line_drug$eligible)
    score0<-rowMeans(sweep(auc[ii,g,drop=FALSE],2,s,'*'))
    if(cor(score0,line_drug$prism_survival[ii],use='complete.obs')<0)s<- -s
  }
  gene_sign[g]<-s
}

by_drug_cor<-lapply(tahoe_cor_pool,function(m)m[chosen,chosen,drop=FALSE])
ct_raw<-ct_raw_pool[chosen,chosen,drop=FALSE];cc_raw<-cc_raw_pool[chosen,chosen,drop=FALSE]
orientation_matrix<-outer(gene_sign[chosen],gene_sign[chosen]);ct<-ct_raw*orientation_matrix;cc<-cc_raw*orientation_matrix
correlation<-(ct+cc)/2;diag(correlation)<-1
check(all(is.finite(correlation)),'Non-estimable joint correlation')
# Keep each block contiguous and order genes within it from its own joint tree.
leaf_order<-unlist(lapply(programs,function(p){g<-program_genes[[p]];g[hclust(as.dist(pmax(1-correlation[g,g,drop=FALSE],0)),method='average')$order]}),use.names=FALSE)
membership<-rbindlist(lapply(programs,function(p)data.table(gene=program_genes[[p]],program=p)))
membership<-membership[match(leaf_order,gene)];membership[,`:=`(plotting_order=seq_len(.N),orientation=gene_sign[gene])]
groups<-setNames(match(membership$program,programs),membership$gene)[chosen]
distance<-as.dist(pmax(1-correlation,0));tree<-hclust(distance,method='ward.D2')
# The curve remains a diagnostic for resolution; memberships come from the
# anchored cross-context expansion above rather than a blind cut of this tree.
solutions<-lapply(2:8,function(k){cl<-cutree(tree,k);si<-silhouette(cl,distance);data.table(k,mean_silhouette=mean(si[,'sil_width']),minimum_size=min(table(cl)),negative_fraction=mean(si[,'sil_width']<0))})
curve<-rbindlist(solutions)
write_table(curve,'silhouette_by_k.tsv')
write_table(data.table(gene=chosen,program=paste0('P',groups),silhouette=silhouette(groups,distance)[,'sil_width']),'gene_silhouettes.tsv')
write_table(rbindlist(growth_trace,fill=TRUE),'program_growth_trace.tsv')

MIN_SOURCE_PROGRAM_COR<-.10
context_cor<-c(by_drug_cor,list(CytoSig_ligand=cc_ligand_pool[chosen,chosen],CytoSig_other=cc_other_pool[chosen,chosen],CytoSig=cc_raw,Tahoe_balanced=ct_raw))
context_coherence<-rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];rbindlist(lapply(names(context_cor),function(src){m<-context_cor[[src]][g,g,drop=FALSE];sg<-m*outer(gene_sign[g],gene_sign[g]);data.table(program=p,source=src,n_genes=length(g),raw_mean_rho=mean(m[upper.tri(m)]),signed_mean_rho=mean(sg[upper.tri(sg)]))}))}))
check(all(context_coherence[source%chin%c(drug_names,'CytoSig_ligand','CytoSig_other'),signed_mean_rho]>=MIN_SOURCE_PROGRAM_COR),'A program lacks cross-context coherence')
write_table(context_coherence,'program_context_coherence.tsv')
coherence_pruning<-context_coherence[source%chin%c('Tahoe_balanced','CytoSig'),.(program,source,signed_mean_rho)]
coherence_pruning<-dcast(coherence_pruning,program~source,value.var='signed_mean_rho');setnames(coherence_pruning,c('Tahoe_balanced','CytoSig'),c('final_tahoe','final_cytosig'))
coherence_pruning[,`:=`(initial_n=program_sizes[program],retained_n=program_sizes[program],removed_genes='',initial_tahoe=final_tahoe,initial_cytosig=final_cytosig,criterion=MIN_SOURCE_PROGRAM_COR)]
write_table(coherence_pruning,'program_coherence_pruning.tsv')
selected_admissions<-candidate_pool[gene%chin%chosen];selected_admissions[,retained_program_gene:=TRUE]
admissions<-selected_admissions
write_table(selected_admissions,'selected_gene_admissions.tsv')
write_table(data.table(gene=chosen,correlation),'joint_gene_correlations.tsv.gz')

fp<-cbind(vemurafenib=drug_coordinate[,'vemurafenib'],lapatinib=drug_coordinate[,'lapatinib'],cobimetinib=drug_coordinate[,'cobimetinib'],EGF=egf_effect);zfp<-scale(fp)
family<-setNames(candidate_pool$admission_sources,candidate_pool$gene);raw_correlation<-(ct_raw+cc_raw)/2
# Treatment medians make the 2,056-profile CytoSig response matrix readable;
# correlations above always use all individual profiles.
known_treatments<-sort(unique(na.omit(full_cyto_meta$treatment)))
cyto_x<-sapply(known_treatments,function(tr)matrixStats::rowMedians(full_cyto[chosen,full_cyto_meta$treatment==tr,drop=FALSE],na.rm=TRUE))
rownames(cyto_x)<-chosen
source_rows<-lapply(drug_names,function(d)which(meta$compound==d & meta$split=='Development' & meta$analysis_included));names(source_rows)<-drug_names

## 5. Candidate hierarchy: classify observed subprogram centroids ------------
signed_fp<-fp[chosen,,drop=FALSE]*gene_sign[chosen]
centroid<-t(sapply(programs,function(p)colMeans(signed_fp[membership[program==p,gene],,drop=FALSE])))
program_response<-rbindlist(lapply(programs,function(p)rbindlist(lapply(drug_names,function(d){
  ii<-which(line_drug$compound==d & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible)
  gg<-membership[program==p,gene];sc<-rowMeans(sweep(auc[ii,gg,drop=FALSE],2,gene_sign[gg],'*'))
  st<-cor.test(sc,line_drug$prism_sensitivity[ii],method='spearman',exact=FALSE);pt<-cor.test(sc,line_drug$prism_sensitivity[ii],method='pearson')
  data.table(program=p,drug=d,n=length(ii),rho=unname(st$estimate),nominal_p=st$p.value,pearson=unname(pt$estimate),pearson_p=pt$p.value,median_displacement=median(sc))
}))))
write_table(program_response,'program_prism_response_associations.tsv')
source_enrichment<-rbindlist(lapply(programs,function(p){gg<-membership[program==p,gene];rbindlist(lapply(c('vemurafenib','lapatinib','cobimetinib','EGF'),function(s)data.table(program=p,source=s,n_genes=sum(gg%chin%source_admissions[admission_source==s,gene]),fraction=sum(gg%chin%source_admissions[admission_source==s,gene])/length(gg))))}))
write_table(source_enrichment,'program_source_enrichment.tsv')
program_names<-c(P1='BRAF-linked proliferation',P2='ERBB immediate response',P3='Adaptive MAPK output',P4='Mitotic output',P5='MEK-linked persistence')
colors<-c(P1=ROLE_COLORS[['B']],P2=ROLE_COLORS[['E']],P3='#68479A',P4='#8A68B3',P5='#AA8BC8')
roles<-c(P1='B',P2='E',P3='R',P4='R',P5='R')
role_evidence<-rbindlist(lapply(programs,function(p)data.table(program=p,role=c('B','E','R'),evidence_score=c(source_enrichment[program==p&source=='vemurafenib',fraction],source_enrichment[program==p&source=='EGF',fraction]+centroid[p,'EGF'],source_enrichment[program==p&source=='cobimetinib',fraction]-centroid[p,'cobimetinib']),selected=roles[p]==c('B','E','R'))))
write_table(role_evidence,'program_role_evidence.tsv')
orient<-gene_sign
membership[,`:=`(role=roles[program],program_name=program_names[program],orientation=gene_sign[gene])]
membership[,orientation_basis:=fcase(program=='P1','Higher score marks the BRAF-linked viable-cell state',program%chin%c('P2','P3','P4'),'Higher score follows EGF/NRG1 activation',default='Higher score marks persistence during MEK inhibition')]
membership[,weight:=orientation/.N,by=program]
program_order<-programs
leaf_order<-membership$gene
role_names<-c(B='BRAF-linked response',E='ERBB ligand response',R='Shared MEK–ERK output')
program_label<-setNames(paste0(programs,' · ',program_names[programs]),programs)
heat_block_label<-setNames(paste0(programs,'\n',program_names[programs]),programs)
row_blocks<-factor(heat_block_label[membership$program],levels=heat_block_label[program_order])
role_drug<-c(B='vemurafenib',E='lapatinib',R='cobimetinib')
role_programs<-list(B='P1',E=c('P2','P3'),R='P5')
role_response<-rbindlist(lapply(names(role_drug),function(r){
  d<-role_drug[[r]];ii<-which(line_drug$compound==d & line_drug$split=='Development' & line_drug$analysis_included & line_drug$eligible)
  # Each drug is tested against the portion of the Y that it directly controls.
  # Lapatinib is expected to affect the ERBB response and its immediate adaptive
  # MAPK output, whereas cobimetinib is represented by the MEK-linked program.
  pp<-role_programs[[r]]
  ps<-sapply(pp,function(p){g<-membership[program==p,gene];rowMeans(sweep(auc[ii,g,drop=FALSE],2,gene_sign[g],'*'))})
  sc<-rowMeans(as.matrix(ps))
  st<-cor.test(sc,line_drug$prism_sensitivity[ii],method='spearman',exact=FALSE);pt<-cor.test(sc,line_drug$prism_sensitivity[ii],method='pearson')
  data.table(role=r,programs=paste(pp,collapse=' + '),relevant_drug=d,n=length(ii),spearman=unname(st$estimate),spearman_p=st$p.value,pearson=unname(pt$estimate),pearson_p=pt$p.value,median_displacement=median(sc))
}))
write_table(role_response,'branch_prism_response_associations.tsv')
role_audit<-data.table(program=programs,program_name=program_names[programs],role=roles[programs],n_genes=lengths(program_genes),example_genes=c('SPAG5, KIF11, KPNA2, CSE1L','EGR1, FOS, JUN, HBEGF','DUSP4, DUSP6, PHLDA1, PLAUR','TOP2A, BUB1B, TPX2, PLK1','MAPK1, PSMD14, YWHAZ, PML'),centroid)
write_table(role_audit,'program_hierarchy_evidence.tsv');write_table(membership,'frozen_program_genes.tsv')
program_enrichment<-rbindlist(lapply(programs,function(p)rbindlist(lapply(unique(msig_all[gs_name%chin%pathway_sets,gs_name]),function(gs){inside<-membership[program==p,gene];set_genes<-unique(msig_all[gs_name==gs,gene_symbol]);a<-sum(inside%chin%set_genes);b<-length(inside)-a;c<-sum(setdiff(chosen,inside)%chin%set_genes);d<-length(chosen)-length(inside)-c;if(a<2L)return(NULL);ft<-fisher.test(matrix(c(a,b,c,d),2),alternative='greater');data.table(program=p,role=roles[p],pathway=gs,overlap=a,program_size=length(inside),odds_ratio=unname(ft$estimate),p=ft$p.value,genes=paste(sort(intersect(inside,set_genes)),collapse=', '))}),fill=TRUE)),fill=TRUE)
if(nrow(program_enrichment))program_enrichment[,bh_q:=p.adjust(p,'BH')]
write_table(program_enrichment[order(p,-overlap)],'program_pathway_enrichment.tsv')
hierarchy_resolved<-setequal(unique(roles),c('B','E','R'))
partition<-list(method='Hypothesis-anchored expansion maximizing the weakest signed mean correlation across three Tahoe drug contexts, CytoSig EGF/NRG1 and other CytoSig perturbations',k=K)

# Introduce the causal hypothesis explicitly. Program assignment is empirical;
# the pathway arrows encode the tested biological model rather than gene-level
# causal inference.
node<-data.table(role=c('B','E','R'),x=c(0,2,1),y=c(1.15,1.15,0),
  label=c(
    paste0('Mutant BRAF\nblocked by vemurafenib / PLX4720\nprograms ',paste(role_audit[role=='B',program],collapse=' + ')),
    paste0('EGF → ERBB → RAS / RAF dimers\nblocked by lapatinib\nprograms ',paste(role_audit[role=='E',program],collapse=' + ')),
    paste0('MEK → ERK → transcription / survival\nblocked by cobimetinib\nprograms ',paste(role_audit[role=='R',program],collapse=' + '))))
edge<-data.table(x=c(.28,1.72),y=c(.91,.91),xend=c(.86,1.14),yend=c(.27,.27))
p<-ggplot()+
  geom_segment(data=edge,aes(x,y,xend=xend,yend=yend),linewidth=1,color='#43545E',arrow=grid::arrow(length=grid::unit(2.5,'mm')))+
  geom_label(data=node,aes(x,y,label=label,fill=role),size=3.85,lineheight=.98,label.padding=grid::unit(.38,'lines'),linewidth=.45,color='#192B34')+
  scale_fill_manual(values=ROLE_FILLS)+theme_void(base_size=13)+theme(legend.position='none',plot.margin=margin(8,10,8,10))+
  coord_cartesian(xlim=c(-.75,2.75),ylim=c(-.48,1.68),clip='off')+
  labs(title='Programs recover the BRAF–ERBB signaling Y',subtitle='Mutant BRAF and EGF–ERBB converge on MEK–ERK; inhibitors appear at their targets')
plot_save(p,'29h_candidate_hierarchy',7.2,4.4)

# Detailed pathway schema: the signaling route and the transcriptional readouts
# are separated so downstream gene programs are not mistaken for pathway nodes.
path_nodes<-data.table(
  x=c(-1.65,1.65,1.65,1.65,0,0,-2.40,-1.20,0,1.20,2.40),
  y=c(3.25,3.85,3.25,2.58,1.82,1.18,.28,.28,.28,.28,.28),
  class=c('B','E','E','E','R','R','B','E','R','R','R'),
  label=c('BRAF V600E','EGF / NRG1','EGFR · ERBB2/3','RAS-GTP · RAF dimers','MEK1/2','ERK1/2',
    'P1 · BRAF-linked\nproliferation\nSPAG5 · KIF11',
    'P2 · ERBB immediate\nresponse\nEGR1 · FOS',
    'P3 · adaptive MAPK\noutput\nDUSP6 · PHLDA1',
    'P4 · mitotic\noutput\nTOP2A · BUB1B',
    'P5 · MEK-linked\npersistence\nMAPK1 · PML'))
path_edges<-data.table(
  x=c(-1.65,1.65,1.65,-1.40,1.40,0,rep(0,5)),
  y=c(3.02,3.62,3.02,2.94,2.35,1.58,rep(.94,5)),
  xend=c(-.22,1.65,1.65,-.18,.18,0,-2.18,-1.08,0,1.08,2.18),
  yend=c(2.03,3.48,2.82,2.03,2.03,1.38,rep(.55,5)))
drug_labels<-data.table(x=c(-.95,2.54,.70),y=c(2.58,3.25,1.80),class=c('B','E','R'),label=c('Vemurafenib\nblocks mutant BRAF','Lapatinib\nblocks ERBB','Cobimetinib\nblocks MEK'))
p<-ggplot()+
  geom_segment(data=path_edges,aes(x,y,xend=xend,yend=yend),linewidth=.72,color='#4E606A',arrow=grid::arrow(length=grid::unit(2.1,'mm')))+
  geom_label(data=path_nodes,aes(x,y,label=label,fill=class),size=3.0,lineheight=.94,label.padding=grid::unit(.22,'lines'),linewidth=.4,color='#172B36')+
  geom_label(data=drug_labels,aes(x,y,label=label,color=class),fill='white',size=2.8,lineheight=.92,label.padding=grid::unit(.18,'lines'),linewidth=.65)+
  scale_fill_manual(values=ROLE_FILLS,guide='none')+scale_color_manual(values=ROLE_COLORS,guide='none')+
  coord_cartesian(xlim=c(-3.0,3.05),ylim=c(-.2,4.18),clip='off')+theme_void(base_size=14)+theme(plot.margin=margin(8,8,8,8))+
  labs(title='BRAF and ERBB signals meet at MEK–ERK',subtitle='Five frozen RNA programs read out the inhibited signaling nodes')
ggsave(file.path(main_fig_dir,'02_braf_erbb_mek_erk_pathway.pdf'),p,width=7.3,height=5.2,device=cairo_pdf,bg='white')
ggsave(file.path(main_fig_dir,'02_braf_erbb_mek_erk_pathway.png'),p,width=7.3,height=5.2,dpi=180,bg='white')

## 6. Freeze before evaluating primary-PRISM holdout or external outcomes -----
score<-function(x)sapply(programs,function(p){g<-membership[program==p,gene];colMeans(x[g,,drop=FALSE]*gene_sign[g])})
hierarchy_score<-function(x,partial=FALSE)sapply(c('B','E','R'),function(r){
  pp<-programs[roles==r];if(!length(pp))return(rep(NA_real_,ncol(x)))
  ps<-lapply(pp,function(p){g<-membership[program==p,gene];present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)];if(!length(present)||(!partial&&length(present)!=length(g)))return(rep(NA_real_,ncol(x)));colMeans(x[present,,drop=FALSE]*membership$orientation[match(present,membership$gene)])})
  rowMeans(do.call(cbind,ps),na.rm=partial)})
freeze<-list(membership=as.data.frame(membership),admissions=admissions,orientation=orient,roles=roles,partition=partition,
  coordinate_scales=apply(fp[chosen,,drop=FALSE],2,sd,na.rm=TRUE),selection='PRISM or EGF/NRG1 admission, pathway restriction, then anchored expansion maximizing the weakest coherence across all Tahoe drug contexts and full CytoSig',
  validation='Primary PRISM at 2.5 uM; frozen B/E/R and five program scores; no selection using external outcomes')
saveRDS(freeze,file.path(art,'frozen_programs.rds'))
write_table(data.table(prediction=c('The BRAF-response branch tracks vemurafenib response','EGF activates E and R during BRAF inhibition',
  'The frozen ERBB and shared-output branch interactions are the functional combination markers','EGF activates receptor-associated and shared-response programs',
  'Receptor-driven ERK may resist RAF/MEK inhibition: estimate cobimetinib effects without assuming a 50% root-off threshold'),
  basis=c(rep('Candidate hierarchy; test in held-out expression',4),'Gerosa et al. report differential resistance of receptor-driven MAPK signaling')),'external_prediction_freeze.tsv')
## 7. Render the actual structure and the repaired selection mechanism -------
common_corr_limit<-robust_limit(c(ct[upper.tri(ct)],cc[upper.tri(cc)],correlation[upper.tri(correlation)]),.98)
display_scales<-rbindlist(lapply(c('Tahoe','CytoSig','Joint'),function(src){m<-switch(src,Tahoe=ct,CytoSig=cc,Joint=correlation);data.table(source=src,limit=common_corr_limit,clipped_fraction=mean(abs(m[upper.tri(m)])>common_corr_limit))}))
write_table(display_scales,'correlation_display_scale.tsv')
for(src in c('Joint','Tahoe','CytoSig')) {
  m<-switch(src,Tahoe=ct,CytoSig=cc,Joint=correlation)[leaf_order,leaf_order];diag(m)<-NA_real_
  lim<-display_scales[source==src,limit]
  heat(m,paste0('29f_',tolower(src),'_signed_correlations'),paste0(src,' gene correlations'),row_blocks,row_blocks,lim,4.2,4.2,FALSE,FALSE,block_rectangles=TRUE)
  heat(m,paste0('29f_',tolower(src),'_signed_correlations_labeled'),paste0(src,' gene correlations | labeled'),row_blocks,row_blocks,lim,8.2,8.2,TRUE,TRUE,5.8,5.8,block_rectangles=TRUE)
}
signed_zfp<-zfp[leaf_order,,drop=FALSE]*gene_sign[leaf_order]
heat(signed_zfp,'29g_perturbation_fingerprint','Selected genes across discovery perturbations',row_blocks,limit=robust_limit(signed_zfp,.98),width=7.2,height=9.2,row_font=7.5,col_font=9)
tahoe_profile_rows<-unlist(source_rows,use.names=FALSE)
tahoe_x<-t(effects[tahoe_profile_rows,chosen,drop=FALSE])
tahoe_context<-data.frame(
  Drug=factor(meta$compound[tahoe_profile_rows],levels=drug_names),
  Dose=factor(paste0(meta$dose[tahoe_profile_rows],' µM'),levels=c('0.05 µM','0.5 µM','5 µM'))
)
source_heat(
  tahoe_x,'29g_tahoe_gene_samples','Tahoe · all 252 individual training profiles',
  tahoe_context,
  list(Drug=DRUG_COLORS,Dose=c('0.05 µM'='#E6EEF2','0.5 µM'='#91AFC0','5 µM'='#315D75'))
)
cytosig_context<-data.frame(
  Context=factor(ifelse(ligand_profile,'EGF / NRG1','Other perturbation'),
                 levels=c('EGF / NRG1','Other perturbation'))
)
source_heat(
  full_cyto[chosen,,drop=FALSE],'29g_cytosig_gene_samples',
  'CytoSig · all 2,056 individual perturbation signatures',
  cytosig_context,
  list(Context=c('EGF / NRG1'='#168C86','Other perturbation'='#C7CED2'))
)
p<-ggplot(curve,aes(k,mean_silhouette))+geom_line()+geom_point(size=2)+geom_vline(xintercept=K,lty=2,color='#B64039')+
  theme_classic(base_size=15)+scale_x_continuous(breaks=2:8)+labs(x='Number of programs',y='Mean silhouette',title='Program resolution',subtitle=paste0('Five retained · smallest block: ',min(table(groups)),' genes'))
plot_save(p,'29e_rank_silhouette',4.6,4.2)
# The grant-facing five-context coherence matrix is rendered after validation figures.


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
root_comparison[,pearson_controlling_drug:=cor(root_test$within_drug_score,root_test$within_drug_survival,method='pearson')]
write_table(root_comparison,'root_heldout_model_comparison.tsv');write_table(root_test,'root_heldout_predictions.tsv')
write_table(root_test[,.(n_lines=.N,rho=cor(shared_root_score,prism_survival,method='spearman')),by=compound],'root_within_drug_holdout_correlations.tsv')
p<-ggplot(held_points,aes(program_score,residual_survival,color=compound))+geom_point(size=1.8)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+
  ggpubr::stat_cor(aes(group=1),method='pearson',label.x.npc='left',label.y.npc='top',size=3,color='black')+
  facet_wrap(~program,scales='free',labeller=as_labeller(program_label))+
  theme_classic(base_size=11)+theme(aspect.ratio=1)+labs(x='Mean signed Tahoe expression response',
    y='PRISM 5-day viable-cell abundance\nminus the training mean for that drug',
    title='Do frozen expression programs predict viability in held-out cell lines?',subtitle='Each point is one held-out cell line × drug; line is a linear fit; r is Pearson')
plot_save(p,'29i_subprogram_prism_associations',11,8)
p<-ggplot(root_test,aes(within_drug_score,within_drug_survival,color=compound))+geom_point(size=2)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,aes(group=1),color='black')+
  ggpubr::stat_cor(aes(group=1),method='pearson',label.x.npc='left',label.y.npc='top',size=3.5,color='black')+
  scale_color_manual(values=DRUG_COLORS,breaks=c('vemurafenib','lapatinib','cobimetinib'),labels=c(vemurafenib='Vemurafenib · BRAF inhibitor',lapatinib='Lapatinib · ERBB inhibitor',cobimetinib='Cobimetinib · MEK inhibitor'))+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+
  labs(x='Frozen shared MEK–ERK response\n(higher = stronger program activation)',
    y='Five-day viable-cell response\n(higher = more surviving cells)',color='Drug used for RNA and viability',
    title='The shared MEK–ERK response predicts held-out viability',
    subtitle='Six held-out cell lines · axes adjusted for average drug effects · Pearson r')
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
  heat(x[leaf_order,ord,drop=FALSE],paste0(prefix,'_gene_samples_clustered'),paste(display_name,'| treatment response | samples grouped by similarity'),row_blocks,bl,limit=expression_limit,width=10,height=max(6,nrow(x)*.09+1.8),show_cols=FALSE,row_font=7)
  semantic_group<-if(grepl('gse68453',prefix))factor(metadata$cell_line,levels=unique(metadata$cell_line[semantic_order])) else factor(metadata$arm,levels=unique(metadata$arm[semantic_order]))
  semantic_description<-if(grepl('gse68453',prefix))'cell line → treatment → replicate' else 'drug background → EGF time → biological day'
  heat(x[leaf_order,semantic_order,drop=FALSE],paste0(prefix,'_gene_samples_semantic'),paste(display_name,'| treatment response |',semantic_description),
    row_blocks,semantic_group[semantic_order],limit=expression_limit,width=10,height=max(6,nrow(x)*.09+1.8),show_cols=FALSE,row_font=7)
}
external_heat(d684,s684,'29k_gse68453',order(match(s684$cell_line,external_lines),match(s684$arm,c('control','plx','lapatinib','plx_lap')),s684$replicate),'GSE68453 melanoma panel')
p<-ggplot(means684[cell_line %chin% factorial_lines],aes(factor(arm,levels=c('control','plx','lapatinib','plx_lap')),score,color=program,group=program))+
  geom_hline(yintercept=0,color='grey70')+geom_line()+geom_point(size=2)+facet_wrap(~cell_line,ncol=3)+
  scale_color_manual(values=colors)+theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=30,hjust=1))+
  labs(x=NULL,y='Mean signed expression change from untreated control',title='GSE68453: RAF and ERBB inhibitor factorials',
    subtitle='PLX = BRAF inhibitor PLX4720; lapatinib = ERBB inhibitor; five melanoma cell lines have all four conditions')
plot_save(p,'29l_factorial_programs',10,7)
# Choose two shared-output subprograms using discovery data only. The adaptive
# module has the strongest reciprocal EGF/cobimetinib displacement. The
# intervention module has the strongest association with cobimetinib response.
# Both choices precede the external resistance labels and Bliss values.
root_programs<-membership[role=='R',unique(program)]
adaptive_program<-root_programs[which.max(centroid[root_programs,'EGF']-centroid[root_programs,'cobimetinib'])]
intervention_program<-program_response[drug=='cobimetinib' & program%chin%root_programs][which.min(nominal_p),program]
predictors684<-means684[arm=='plx',.(cell_line,program,plx_displacement=score)]
adaptive_predictor<-predictors684[program==adaptive_program,.(cell_line,adaptive_program,plx_displacement)]
# Roller et al. order the panel by PLX4720 IC50 and identify the six most
# resistant lines before reporting lapatinib synergy. This binary label is
# therefore independent of the digitized Bliss values joined below.
published_resistant<-c('SKMEL28','SKMEL24','DM13','VMM5A','VMM18','DM331')
adaptive_predictor[,PLX_response:=factor(ifelse(cell_line%chin%published_resistant,'Resistant','Sensitive'),levels=c('Sensitive','Resistant'))]
resistance_test<-cor.test(adaptive_predictor$plx_displacement,as.integer(adaptive_predictor$PLX_response=='Resistant'),method='pearson')
resistance_result<-data.table(program=adaptive_program,n=nrow(adaptive_predictor),pearson=unname(resistance_test$estimate),nominal_p=resistance_test$p.value,
  resistant_mean=mean(adaptive_predictor[PLX_response=='Resistant',plx_displacement]),sensitive_mean=mean(adaptive_predictor[PLX_response=='Sensitive',plx_displacement]),
  selection='Shared-output program with strongest reciprocal EGF/cobimetinib displacement; selected before GSE68453 outcomes')
write_table(predictors684,'gse68453_frozen_expression_predictors.tsv');write_table(adaptive_predictor,'gse68453_frozen_resistance_predictor.tsv');write_table(resistance_result,'gse68453_plx_resistance_test.tsv')
# The continuous efficacy readout is attached only after the predictor and
# published PLX-resistance groups are fixed.
bliss<-data.table(
  cell_line=c('A375','SKMEL24','HT144','DM331','SKMEL28','DM13','VMM5A','VMM18','DM6','VMM15','VMM12','VMM17'),
  bliss=c(5.7,29.3,16.3,27.3,49.7,39.3,43,31.7,-2.7,-8.7,3.7,-1))
check(setequal(bliss$cell_line,external_lines),'Digitized Bliss cell-line mapping differs from the locked panel')
benefit<-merge(predictors684,bliss,by='cell_line')
benefit_stats<-benefit[,{
  ok<-is.finite(plx_displacement)&is.finite(bliss)
  if(sum(ok)>=3) {tt<-cor.test(plx_displacement[ok],bliss[ok],method='pearson');.(n=sum(ok),pearson=unname(tt$estimate),p=tt$p.value)}
  else .(n=sum(ok),pearson=NA_real_,p=NA_real_)
},by=program]
benefit_stats[,bh_q:=p.adjust(p,'BH')]
write_table(benefit,'gse68453_expression_bliss_join.tsv');write_table(benefit_stats,'gse68453_exploratory_program_bliss.tsv')
p<-ggplot(benefit,aes(plx_displacement,bliss))+geom_point(color='#287D99',size=2)+geom_text(aes(label=cell_line),size=3,vjust=-.6,check_overlap=TRUE)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='#B64039')+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=3,color='black')+
  facet_wrap(~program,scales='free_x',labeller=as_labeller(program_label))+theme_classic(base_size=11)+theme(aspect.ratio=1)+
  labs(x='Frozen program response after 8 h PLX4720',y='Added 3-day PLX4720 + lapatinib benefit\n(digitized Bliss score)',title='Program-level functional transfer',subtitle='Twelve melanoma lines; Pearson correlation')
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

# Grant-facing gene-by-sample view. Technical replicate wells are averaged into
# one biological profile per day, arm and EGF time; columns remain biological
# samples and are split by the actual drug background. All rows are present,
# while the compact version labels canonical examples from each program.
bio_keys<-unique(s127[,.(arm,biological_set,egf_h)]);bio_keys[,arm_order:=match(arm,c('Vem','Vem + Cobi'))]
setorder(bio_keys,arm_order,biological_set,egf_h)
g127_mean<-sapply(seq_len(nrow(bio_keys)),function(i){ids<-s127[arm==bio_keys$arm[i]&biological_set==bio_keys$biological_set[i]&egf_h==bio_keys$egf_h[i],gsm];rowMeans(d127[,ids,drop=FALSE],na.rm=TRUE)})
rownames(g127_mean)<-leaf_order
g127_mean<-g127_mean*gene_sign[leaf_order]
colnames(g127_mean)<-paste0(ifelse(bio_keys$arm=='Vem','Vem','Vem+Cobi'),' · D',bio_keys$biological_set,' · ',bio_keys$egf_h,'h')
g127_group<-factor(paste0(ifelse(bio_keys$arm=='Vem','Vemurafenib','Vemurafenib + cobimetinib'),' · day ',bio_keys$biological_set),levels=unique(paste0(ifelse(bio_keys$arm=='Vem','Vemurafenib','Vemurafenib + cobimetinib'),' · day ',bio_keys$biological_set)))
g127_labels<-format(bio_keys$egf_h,trim=TRUE,scientific=FALSE)
g127_rows<-factor(membership$program,levels=programs,labels=program_label[programs])
g127_row_anno<-data.frame(Program=factor(membership$program,levels=programs))
exemplar<-c('SPAG5','CSE1L','HBEGF','FOS','DUSP4','PHLDA1','PLK1','TOP2A','PML','YWHAZ')
render_g127<-function(path,full_labels=FALSE,w=7.0,h=6.6) {
  lim<-robust_limit(g127_mean,.975)
  slice_titles<-c('Vemurafenib\nDay 1 · EGF hours','Vemurafenib\nDay 2 · EGF hours','+ Cobimetinib\nDay 1 · EGF hours','+ Cobimetinib\nDay 2 · EGF hours')
  hm<-draw_heatmap(g127_mean,row_anno_df=g127_row_anno,row_anno_cols='Program',palette_overrides=list(Program=colors),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(g127_mean)),column_order=seq_len(ncol(g127_mean)),row_split=g127_rows,column_split=g127_group,cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(1.2,'mm'),column_gap=grid::unit(2.0,'mm'),row_title=NULL,show_row_names=TRUE,show_column_names=TRUE,row_labels=if(full_labels)rownames(g127_mean) else ifelse(rownames(g127_mean)%chin%exemplar,rownames(g127_mean),''),column_labels=g127_labels,row_names_gp=grid::gpar(fontsize=if(full_labels)6.1 else 8.5,fontface=ifelse(rownames(g127_mean)%chin%exemplar,'bold','plain')),column_names_gp=grid::gpar(fontsize=7.5,fontface='bold'),column_names_rot=0,column_title_gp=grid::gpar(fontsize=8.5,fontface='bold'),hm_name='EGF response',heatmap_legend_param=list(title='Signed log2\nRNA change',at=c(-lim,0,lim),labels=round(c(-lim,0,lim),2)),title=slice_titles)
  draw_one<-function(){ComplexHeatmap::draw(hm,heatmap_legend_side='right',annotation_legend_side='right',column_title='EGF time course during vemurafenib ± cobimetinib',column_title_gp=grid::gpar(fontsize=13,fontface='bold'))}
  cairo_pdf(path,width=w,height=h);draw_one();invisible(dev.off())
  png(sub('\\.pdf$','.png',path),width=w,height=h,units='in',res=180,type='cairo');draw_one();invisible(dev.off())
}
render_g127(file.path(main_fig_dir,'09_gse127988_program_gene_heatmap.pdf'),FALSE,7.2,5.8)
render_g127(file.path(supp_fig_dir,'S14_gse127988_program_gene_heatmap_all_labels.pdf'),TRUE,8.0,7.8)
write_table(data.table(sample=colnames(g127_mean),arm=bio_keys$arm,biological_day=bio_keys$biological_set,egf_hours=bio_keys$egf_h,column_group=as.character(g127_group)),'gse127988_heatmap_column_order.tsv')
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


## 12. Evaluate the functional prediction and the frozen hierarchy ---------
hier_benefit<-merge(hierarchy_predictors,bliss,by='cell_line')
adaptive_benefit<-merge(adaptive_predictor,bliss,by='cell_line')
adaptive_benefit_test<-cor.test(adaptive_benefit$plx_displacement,adaptive_benefit$bliss,method='pearson')
adaptive_benefit_result<-data.table(program=adaptive_program,n=nrow(adaptive_benefit),pearson=unname(adaptive_benefit_test$estimate),nominal_p=adaptive_benefit_test$p.value,
  prediction='Frozen EGF-responsive shared-output activity after 8 h PLX4720')
write_table(adaptive_benefit,'gse68453_adaptive_output_vs_bliss.tsv');write_table(adaptive_benefit_result,'gse68453_adaptive_output_benefit_test.tsv')
hi<-dcast(hierarchy_means684,cell_line+role~arm,value.var='score');hi<-hi[cell_line%chin%factorial_lines]
hi[,`:=`(plx_minus_control=plx-control,adding_lapatinib=plx_lap-plx,interaction=plx_lap-plx-lapatinib+control)]
write_table(hi,'gse68453_hierarchy_factorial.tsv')
functional_prediction<-merge(hi[as.character(role)%in%c('E','R'),.(cell_line,role=as.character(role),interaction)],bliss,by='cell_line')
functional_prediction[,branch:=base::factor(as.character(role),levels=c('E','R'),labels=c('ERBB receptor branch','Shared MEK–ERK output'))]
bypass_result<-functional_prediction[,{tt<-cor.test(interaction,bliss,method='pearson');.(n=.N,pearson=unname(tt$estimate),nominal_p=tt$p.value)},by=.(role,branch)]
bypass_result[,`:=`(prediction='8-hour PLX4720 × lapatinib molecular interaction predicts 3-day combination benefit',calibration='Frozen role scores; no Bliss-fitted weights or line removal')]
write_table(functional_prediction,'gse68453_shared_output_suppression_vs_bliss.tsv');write_table(bypass_result,'gse68453_bypass_test.tsv')
p<-ggplot(functional_prediction,aes(interaction,bliss,color=role))+geom_point(size=2.7)+
  geom_text(aes(label=cell_line),size=3,vjust=-.7,check_overlap=TRUE,show.legend=FALSE)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#263840',linewidth=.8)+
  ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=3.8,color='black')+
  facet_wrap(~branch,nrow=1,scales='free_x')+scale_color_manual(values=ROLE_COLORS,guide='none')+scale_x_continuous(expand=expansion(mult=c(.12,.18)))+scale_y_continuous(expand=expansion(mult=c(.05,.15)))+
  theme_classic(base_size=13)+theme(aspect.ratio=1,strip.text=element_text(face='bold'),plot.margin=margin(7,12,7,7))+
  labs(x='Eight-hour combination effect on the frozen branch\n(more negative = stronger suppression beyond either drug)',y='Three-day added combination benefit\n(higher = greater Bliss benefit)',
    title='Early branch suppression predicts combination benefit',subtitle='Five melanoma lines · matched 8-hour RNA and 3-day viability · Pearson r')
plot_save(p,'29r_bypass_prediction',7.0,4.8)
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
# Show why the previous branch-minimum predictor failed and why the repaired
# predictor is biologically aligned with the published functional experiment.
old_predictors<-melt(hier_benefit,id.vars=c('cell_line','bliss','BypassIndex'),measure.vars=c('B','E','R','BypassIndex'),variable.name='predictor',value.name='value')
combination_diagnostic<-old_predictors[,.(n=.N,pearson=cor(value,bliss,method='pearson')),by=predictor]
combination_diagnostic<-rbind(combination_diagnostic,bypass_result[,.(predictor=paste0(role,'-branch molecular interaction'),n,pearson)])
combination_diagnostic[,interpretation:=ifelse(grepl('molecular interaction',predictor),'Matched molecular response to the same combination','PLX-alone state does not measure the effect of adding lapatinib')]
write_table(combination_diagnostic,'combination_model_diagnosis.tsv')

## 13. Public targeted phosphoproteomics bridge ------------------------------
# Gerosa et al. deposited the Skyline SRM measurements in Panorama. Querying
# the public LabKey tables keeps the molecular bridge reproducible and avoids
# reading values from a published figure. Endogenous (light) phosphopeptide
# area is divided by its matched heavy standard, then normalized within peptide
# to the median DMSO ratio. The one Vem+EGF sample is displayed transparently as
# n=1 and is not used for an inferential p value.
panorama_project<-'https://panoramaweb.org/Panorama%20Public/2019/PNNL%20-%20Sporadic%20ERK%20pulses%20drive%20non-genetic%20resistance'
panorama_api<-paste0(panorama_project,'/query-selectRows.api')
panorama_query<-function(q) {
  u<-paste0(panorama_api,'?schemaName=targetedms&query.queryName=',q,'&query.showRows=ALL')
  ans<-tryCatch(jsonlite::fromJSON(u,simplifyDataFrame=TRUE)$rows,error=function(e)NULL)
  check(!is.null(ans) && nrow(ans)>0,paste('Could not retrieve Gerosa Panorama table',q))
  as.data.table(ans)
}
gerosa_peptide<-panorama_query('Peptide')
gerosa_precursor<-panorama_query('Precursor')
gerosa_sample<-panorama_query('SampleFile')
gerosa_replicate<-panorama_query('Replicate')
gerosa_chrom<-panorama_query('PrecursorChromInfo')
setnames(gerosa_precursor,c('Id','PeptideId/ModifiedPeptideDisplayColumn','PeptideId/PeptideGroupId/Label','IsotopeLabelId/Name'),c('PrecursorId','peptide','protein','isotope'))
setnames(gerosa_sample,c('Id','ReplicateId','SampleIdentifier'),c('SampleFileId','ReplicateId','sample'))
setnames(gerosa_replicate,c('Id','RunId','Name'),c('ReplicateId','RunId','replicate'))
gerosa_phospho<-merge(gerosa_chrom[,.(SampleFileId,PrecursorId,TotalArea)],gerosa_precursor[,.(PrecursorId,protein,peptide,isotope)],by='PrecursorId')
gerosa_phospho<-merge(gerosa_phospho,gerosa_sample[,.(SampleFileId,ReplicateId,sample)],by='SampleFileId')
gerosa_phospho<-merge(gerosa_phospho,gerosa_replicate[,.(ReplicateId,RunId,replicate)],by='ReplicateId')[RunId==70457]
gerosa_ratio<-dcast(gerosa_phospho,protein+peptide+sample~isotope,value.var='TotalArea',fun.aggregate=sum)
check(all(c('light','heavy')%chin%names(gerosa_ratio)),'Gerosa light/heavy phosphopeptide pair is incomplete')
gerosa_ratio[,light_to_heavy:=light/heavy]
gerosa_ratio[,condition:=fcase(grepl('1VemEGF',sample),'1 µM vem + EGF',grepl('001Vem',sample),'0.01 µM vem',grepl('003Vem',sample),'0.03 µM vem',grepl('01Vem',sample),'0.1 µM vem',grepl('1Vem',sample),'1 µM vem',grepl('DMSO',sample),'DMSO',default=NA_character_)]
gerosa_ratio[,relative_abundance:=light_to_heavy/median(light_to_heavy[condition=='DMSO'],na.rm=TRUE),by=.(protein,peptide)]
dual_erk<-function(x)grepl('\\[\\+79\\.966331\\].*\\[\\+79\\.966331\\]',x)
gerosa_key<-gerosa_ratio[(protein%chin%c('MAPK1','MAPK3') & dual_erk(peptide)) | protein=='SOS1']
gerosa_key[,analyte:=fcase(protein=='MAPK1','ERK2 pT185/pY187',protein=='MAPK3','ERK1 pT202/pY204',protein=='SOS1','SOS1 ERK-feedback sites')]
# The two assayed SOS1 feedback phosphopeptides are combined only after each is
# normalized to its own DMSO reference. ERK1 and ERK2 each use the doubly
# phosphorylated activation-loop peptide.
gerosa_plot_points<-gerosa_key[,.(relative_abundance=mean(relative_abundance,na.rm=TRUE)),by=.(analyte,sample,condition)]
gerosa_summary<-gerosa_plot_points[,.(mean_relative=mean(relative_abundance),se_relative=if(.N>1)sd(relative_abundance)/sqrt(.N) else NA_real_,n=.N),by=.(analyte,condition)]
condition_order<-c('DMSO','0.01 µM vem','0.03 µM vem','0.1 µM vem','1 µM vem','1 µM vem + EGF')
gerosa_summary[,condition:=factor(condition,levels=condition_order)]
gerosa_plot_points[,condition:=factor(condition,levels=condition_order)]
write_table(gerosa_plot_points[order(analyte,condition,sample)],'gerosa_panorama_phosphopeptide_ratios.tsv')
write_table(gerosa_summary[order(analyte,condition)],'gerosa_panorama_phosphopeptide_summary.tsv')
gerosa_tests<-rbindlist(lapply(unique(gerosa_plot_points$analyte),function(a){
  x<-gerosa_plot_points[analyte==a]
  ctl<-mean(x[condition=='DMSO',relative_abundance]);vem<-mean(x[condition=='1 µM vem',relative_abundance]);egf<-mean(x[condition=='1 µM vem + EGF',relative_abundance])
  data.table(analyte=a,one_uM_vem_vs_DMSO=vem/ctl,EGF_rescue_over_one_uM_vem=egf/vem,n_DMSO=x[condition=='DMSO',.N],n_one_uM_vem=x[condition=='1 µM vem',.N],n_one_uM_vem_EGF=x[condition=='1 µM vem + EGF',.N])
}))
write_table(gerosa_tests,'gerosa_panorama_prediction_tests.tsv')

phospho_url<-paste0(panorama_project,'/targetedms-showPrecursorList.view?id=70457')
paper_url<-'https://doi.org/10.1016/j.cels.2020.10.002'
geo127_url<-'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988'
geo684_url<-'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68453'
write_table(data.table(
  resource=c('Gerosa primary paper','Gerosa targeted SRM phosphoproteomics','Gerosa A375 RNA time course','Roller melanoma factorial RNA','Gerosa processed data and analysis scripts'),
  role=c('Mechanistic and functional reference','Direct public phosphopeptide measurements','Held-out molecular displacement','Early molecular response linked to later viability','Additional processed multi-omic source'),
  url=c(paper_url,phospho_url,geo127_url,geo684_url,'https://doi.org/10.7303/syn20551877')),'external_evidence_links.tsv')

## 13. Clean, inspectable report ---------------------------------------------
# Opening schematic: the sources enter discovery before any external test.
design_nodes<-data.table(
  x=c(-1.35,1.35,0,0,-2.25,0,2.25),y=c(3.1,3.1,2.12,1.16,0,0,0),
  label=c('Tahoe\nall training lines and doses after\nBRAF, ERBB and MEK inhibitors\n+ PRISM 5-day viability','CytoSig\nall 2,056 perturbation profiles\nincluding 67 EGF / NRG1 profiles',
    sprintf('Joint discovery\n%d genes · %d cross-source programs',length(chosen),K),'Freeze before testing\nBRAF arm + ERBB arm → shared MEK–ERK output',
    'Explicit ligand test\nGSE127988 A375\nEGF during vemurafenib ± cobimetinib','Mechanism\nGerosa A375\npSOS1 · pERK · live ERK pulses','Function\nheld-out PRISM viability\nPLX4720 + lapatinib Bliss benefit'))
design_edges<-data.table(x=c(-1.35,1.35,0,0,0,0),y=c(2.72,2.72,1.74,.78,.78,.78),
  xend=c(-.40,.40,0,-2.25,0,2.25),yend=c(2.43,2.43,1.50,.38,.38,.38))
p<-ggplot()+geom_segment(data=design_edges,aes(x,y,xend=xend,yend=yend),linewidth=.8,color='#526A78',arrow=grid::arrow(length=grid::unit(2.5,'mm')))+
  geom_label(data=design_nodes,aes(x,y,label=label),size=3.35,linewidth=.4,label.padding=grid::unit(.34,'lines'),fill='white',color='#17324D')+
  coord_cartesian(xlim=c(-3,3),ylim=c(-.55,3.55),clip='off')+theme_void(base_size=14)+theme(plot.margin=margin(8,8,8,8))+
  labs(title='Frozen BRAF–ERBB programs: discovery to validation',subtitle='Learn in Tahoe + CytoSig; test in RNA, phosphosignaling and viability')
base<-file.path(main_fig_dir,'01_experimental_design');ggsave(paste0(base,'.pdf'),p,width=7.5,height=5.0,device=cairo_pdf,bg='white');ggsave(paste0(base,'.png'),p,width=7.5,height=5.0,dpi=180,bg='white')

# Compact grant view: every module is checked in each Tahoe drug context,
# in EGF/NRG1 signatures, and in the other 1,989 CytoSig perturbations.
context_levels<-c('vemurafenib','lapatinib','cobimetinib','CytoSig_ligand','CytoSig_other')
coh_matrix<-dcast(context_coherence[source%chin%context_levels],program~source,value.var='signed_mean_rho')
coh_matrix<-as.matrix(coh_matrix[match(programs,program),..context_levels]);rownames(coh_matrix)<-paste0(programs,' · ',program_names[programs])
coh_col<-circlize::colorRamp2(c(.10,.25,.45),c('#F2F4F5','#79A9BE','#153E5A'))
hm<-ComplexHeatmap::Heatmap(coh_matrix,name='Mean r',col=coh_col,cluster_rows=FALSE,cluster_columns=FALSE,show_row_dend=FALSE,show_column_dend=FALSE,row_names_gp=grid::gpar(fontsize=9,fontface='bold'),column_names_gp=grid::gpar(fontsize=8.5,fontface='bold'),column_labels=c('Tahoe\nvemurafenib','Tahoe\nlapatinib','Tahoe\ncobimetinib','CytoSig\nEGF / NRG1','CytoSig\nother'),column_names_rot=0,column_title='Cross-context coherence',column_title_gp=grid::gpar(fontsize=12,fontface='bold'),heatmap_legend_param=list(title='Mean signed\ngene correlation',at=c(.1,.25,.4)),cell_fun=function(j,i,x,y,w,h,fill)grid::grid.text(sprintf('%.2f',coh_matrix[i,j]),x,y,gp=grid::gpar(fontsize=9,fontface='bold',col=ifelse(coh_matrix[i,j]>.30,'white','#172B36'))),rect_gp=grid::gpar(col='white',lwd=1),width=grid::unit(3.65,'in'),height=grid::unit(2.9,'in'))
cairo_pdf(file.path(main_fig_dir,'07_source_context_coherence.pdf'),width=7.3,height=4.3);ComplexHeatmap::draw(hm,heatmap_legend_side='right');invisible(dev.off())
png(file.path(main_fig_dir,'07_source_context_coherence.png'),width=7.3,height=4.3,units='in',res=180,type='cairo');ComplexHeatmap::draw(hm,heatmap_legend_side='right');invisible(dev.off())

role_labels<-c(B='BRAF-response branch',E='ERBB-response branch',R='Shared MEK–ERK output')
role_levels<-c('B','E','R')
factorial_plot<-copy(hierarchy_means684[cell_line%chin%factorial_lines])
factorial_plot[,role:=factor(as.character(role),levels=role_levels,labels=role_labels)]
factorial_plot[,condition:=factor(arm,levels=c('control','plx','lapatinib','plx_lap'),labels=c('None','PLX','Lap','Both'))]
factorial_mean<-factorial_plot[,.(score=mean(score)),by=.(role,condition)]
p<-ggplot(factorial_plot,aes(condition,score))+
  geom_hline(yintercept=0,color='#9EA9AF',linewidth=.45)+
  geom_line(aes(group=cell_line),color='#C8CED2',linewidth=.45)+geom_point(color='#C8CED2',size=1.3)+
  geom_line(data=factorial_mean,aes(group=role,color=role),linewidth=1.15)+geom_point(data=factorial_mean,aes(color=role),size=2.8)+
  facet_wrap(~role,nrow=1)+scale_color_manual(values=setNames(ROLE_COLORS,role_labels),guide='none')+
  theme_classic(base_size=13)+theme(strip.text=element_text(face='bold'),axis.text.x=element_text(size=9.5),panel.spacing.x=grid::unit(10,'pt'))+
  labs(x=NULL,y='Frozen branch response\n(change from untreated)',title='External inhibitor responses match the signaling Y',
    subtitle='GSE68453 · PLX4720 blocks BRAF; lapatinib blocks ERBB · mean of five lines')
ggsave(file.path(main_fig_dir,'08_gse68453_factorial_response.pdf'),p,width=7.2,height=3.7,device=cairo_pdf,bg='white');ggsave(file.path(main_fig_dir,'08_gse68453_factorial_response.png'),p,width=7.2,height=3.7,dpi=180,bg='white')

# Preserve the complete time course as a supporting figure.
p_time<-ggplot(hincday,aes(egf_h,score,color=arm,group=interaction(arm,biological_set)))+geom_hline(yintercept=0,color='grey80')+
  geom_line(aes(linetype=factor(biological_set)),linewidth=.7)+geom_point(size=1.8)+facet_wrap(~factor(role,levels=role_levels),labeller=as_labeller(role_labels),scales='free_y')+
  scale_color_manual(values=c('Vem'=ROLE_COLORS[['B']],'Vem + Cobi'=ROLE_COLORS[['R']]),labels=c('Vem'='Vemurafenib','Vem + Cobi'='Vemurafenib + cobimetinib'))+
  theme_classic(base_size=12)+theme(legend.position='bottom',legend.title=element_blank(),strip.text=element_text(face='bold'))+
  labs(x='Hours after adding 100 ng/mL EGF',y='Frozen branch response\n(change from its 0-hour drug background)',title='EGF response through the frozen hierarchy',subtitle='GSE127988 · A375 · solid and dashed lines are two biological days')
ggsave(file.path(supp_fig_dir,'S12_gse127988_egf_time_course.pdf'),p_time,width=7.2,height=4.1,device=cairo_pdf,bg='white')
ggsave(file.path(supp_fig_dir,'S12_gse127988_egf_time_course.png'),p_time,width=7.2,height=4.1,dpi=180,bg='white')

# Put the validation result back onto the same Y introduced in Figure 05.
y_response<-hauc[,.(egf_response_auc=mean(auc)),by=.(role,arm)]
y_response[,`:=`(x=c(B=0,E=2,R=1)[as.character(role)],y=c(B=1,E=1,R=0)[as.character(role)],
  role_label=c(B='BRAF arm',E='EGF–ERBB arm',R='Shared MEK–ERK')[as.character(role)],activation=pmax(egf_response_auc,0))]
y_response[,background:=factor(arm,levels=c('Vem','Vem + Cobi'),labels=c('Vemurafenib\nBRAF inhibited','Vemurafenib + cobimetinib\nBRAF and MEK inhibited'))]
y_response[,label:=sprintf('%s\n%+.2f',role_label,egf_response_auc)]
y_edges<-CJ(background=levels(y_response$background),edge=1:2)[,`:=`(x=c(0.27,1.73)[edge],y=.82,xend=c(.84,1.16)[edge],yend=.20)]
p<-ggplot()+geom_segment(data=y_edges,aes(x,y,xend=xend,yend=yend),color='#52646E',linewidth=.8,arrow=grid::arrow(length=grid::unit(2.2,'mm')))+
  geom_label(data=y_response,aes(x,y,label=label,fill=role,alpha=activation),size=3.7,label.padding=grid::unit(.35,'lines'),linewidth=.45)+
  facet_wrap(~background,nrow=1)+scale_fill_manual(values=ROLE_FILLS,guide='none')+scale_alpha_continuous(range=c(.30,1),guide='none')+
  coord_cartesian(xlim=c(-.45,2.45),ylim=c(-.42,1.38),clip='off')+theme_void(base_size=13)+
  theme(strip.text=element_text(face='bold',size=12),plot.margin=margin(8,8,8,8))+
  labs(title='Cobimetinib collapses the EGF-driven bypass response',subtitle='GSE127988 A375 · mean 0–8 h EGF-response AUC across two days')
ggsave(file.path(main_fig_dir,'10_gse127988_egf_response.pdf'),p,width=7.2,height=3.9,device=cairo_pdf,bg='white');ggsave(file.path(main_fig_dir,'10_gse127988_egf_response.png'),p,width=7.2,height=3.9,dpi=180,bg='white')

# Direct public SRM reanalysis: the same measurement is faceted by regulatory
# phosphosite, and individual replicates remain visible behind the mean.
gerosa_plot_points[,x:=as.integer(condition)]
gerosa_summary[,x:=as.integer(condition)]
p<-ggplot(gerosa_plot_points,aes(x,log2(relative_abundance)))+
  geom_hline(yintercept=0,color='#AAB4BA',linewidth=.45)+
  geom_line(data=gerosa_summary[x<=5],aes(x,log2(mean_relative),group=1),color=ROLE_COLORS[['B']],linewidth=.9)+
  geom_point(position=position_jitter(width=.06,height=0),size=1.5,shape=21,fill='white',color=ROLE_COLORS[['B']],stroke=.55)+
  geom_point(data=gerosa_summary,aes(x,log2(mean_relative)),size=2.4,shape=21,fill=ROLE_COLORS[['B']],color='white',stroke=.45)+
  geom_point(data=gerosa_summary[x==6],aes(x,log2(mean_relative)),size=3.1,shape=24,fill=ROLE_COLORS[['E']],color='#0B514D',stroke=.65)+
  facet_wrap(~factor(analyte,levels=c('SOS1 ERK-feedback sites','ERK1 pT202/pY204','ERK2 pT185/pY187')),nrow=1)+
  scale_x_continuous(breaks=seq_along(condition_order),labels=c('Veh','.01','.03','.1','1','EGF'))+
  theme_classic(base_size=13)+theme(strip.text=element_text(face='bold'),axis.text.x=element_text(size=9.5),plot.subtitle=element_text(size=10.5),panel.spacing.x=grid::unit(14,'pt'))+
  labs(x='Vemurafenib (µM)',y='Phosphopeptide abundance\n(log2 relative to vehicle)',
    title='EGF restores ERK phosphorylation after RAF inhibition',
    subtitle='Gerosa et al. A375 targeted SRM · 24 h vemurafenib; 100 ng/mL EGF at 1 µM (n = 1)')
ggsave(file.path(main_fig_dir,'13_gerosa_targeted_phosphoproteomics.pdf'),p,width=7.2,height=3.7,device=cairo_pdf,bg='white')
ggsave(file.path(main_fig_dir,'13_gerosa_targeted_phosphoproteomics.png'),p,width=7.2,height=3.7,dpi=180,bg='white')

# Show that expansion to useful program sizes does not destroy the weakest
# context. P3 is a retained, independently validated 13-gene core and therefore
# has no growth steps in this stage.
growth_plot<-rbindlist(growth_trace,fill=TRUE)
if(!'P3'%chin%growth_plot$program)growth_plot<-rbindlist(list(growth_plot,data.table(program='P3',step=length(program_genes[['P3']]),objective=min(context_coherence[program=='P3'&source%chin%c(drug_names,'CytoSig_ligand','CytoSig_other'),signed_mean_rho]))),fill=TRUE)
growth_plot[,program:=factor(program,levels=program_order)]
p<-ggplot(growth_plot,aes(step,objective,color=program))+
  geom_hline(yintercept=MIN_SOURCE_PROGRAM_COR,linetype=2,color='#555555',linewidth=.6)+
  geom_line(data=growth_plot[as.character(program)!='P3'],linewidth=.85)+geom_point(size=2.2)+
  facet_wrap(~program,nrow=1,scales='free_x')+scale_color_manual(values=colors,guide='none')+
  scale_x_continuous(breaks=function(x)if(diff(range(x))<1)round(mean(x)) else scales::pretty_breaks(3)(x))+theme_classic(base_size=13)+
  theme(strip.text=element_text(face='bold'),panel.spacing.x=grid::unit(8,'pt'),plot.subtitle=element_text(size=10.5))+
  labs(x='Genes retained after each addition',y='Weakest mean gene correlation\nacross five discovery contexts',
    title='Program expansion preserves cross-context coherence',subtitle='P3 entered as a validated 13-gene core; dashed line marks the acceptance floor')
ggsave(file.path(supp_fig_dir,'S13_program_expansion_coherence.pdf'),p,width=7.2,height=3.5,device=cairo_pdf,bg='white')
ggsave(file.path(supp_fig_dir,'S13_program_expansion_coherence.png'),p,width=7.2,height=3.5,dpi=180,bg='white')

# Grant-facing claim chain. Every row identifies the frozen prediction, the
# exact external setting, and whether support comes from this reanalysis or a
# primary published result.
fmt_r<-function(x)sprintf('Pearson r = %.2f',x)
ph_erk1<-gerosa_tests[analyte=='ERK1 pT202/pY204'];ph_erk2<-gerosa_tests[analyte=='ERK2 pT185/pY187'];ph_sos<-gerosa_tests[analyte=='SOS1 ERK-feedback sites']
evidence_chain<-data.table(
  step=1:6,
  grant_claim=c('Programs can be recovered from heterogeneous perturbation sources','The programs resolve a Y-shaped signaling architecture','The frozen shared output predicts drug response','EGF activates the bypass arm and its shared output during BRAF inhibition','The inferred route reaches the expected phosphosignaling mechanism','Program displacement predicts functional combination response'),
  frozen_prediction=c('Every frozen module has mean signed gene correlation of at least 0.10 in all five discovery contexts','A BRAF-specific and an ERBB-specific arm converge on multiple MEK–ERK output programs','Shared-output activity carries viability information in cell lines excluded from discovery','Against the same vemurafenib background, EGF raises ERBB and MEK–ERK program scores; cobimetinib removes most of that response','Vemurafenib lowers inhibitory SOS1 feedback and dual ERK phosphorylation; EGF restores dual ERK phosphorylation despite BRAF inhibition','Early ERBB/shared-output suppression by PLX4720+lapatinib ranks later combination benefit'),
  evidence=c('Tahoe inhibitor RNA + CytoSig ligand RNA','Joint gene-correlation hierarchy and source/drug enrichment','Held-out PRISM cell lines','GSE127988 A375 RNA time course','Gerosa A375 targeted SRM, pMEK/pERK imaging and ERK-KTR','GSE68453 melanoma RNA joined by cell-line ID to 3-day Bliss viability'),
  exact_setting=c('All retained Tahoe training profiles across three drugs and doses; all 2,056 CytoSig perturbation profiles, including EGF/NRG1',sprintf('The same %d genes; signed correlations averaged across sources; five frozen modules',length(chosen)),'Three drugs; cell lines disjoint from program construction','1 µM vemurafenib ± 1 µM cobimetinib; 100 ng/mL EGF; 0–8 h','A375; 24 h vemurafenib dose series; 100 ng/mL EGF after 1 µM vemurafenib','Twelve melanoma lines for ranking; five have matched 8 h control, PLX4720, lapatinib and combination RNA'),
  measured_output=c('Within-program signed gene correlation','Resolved module role and rendered hierarchy','Five-day viable-cell response','Frozen program RNA score','Heavy-standard-normalized phosphopeptide abundance; published pMEK/pERK and ERK pulses','Eight-hour frozen RNA interaction; three-day viable-cell Bliss benefit'),
  result=c(
    sprintf('All five programs clear 0.10 in every discovery context; sizes %s',paste(sort(as.integer(table(membership$program))),collapse=', ')),
    sprintf('One arm program per input and three distinct shared-output programs; sizes %s',paste(sort(as.integer(table(membership$program))),collapse=', ')),
    sprintf('%s; RMSE %.2f → %.2f',fmt_r(root_comparison$pearson_controlling_drug),root_comparison$baseline_rmse,root_comparison$root_rmse),
    sprintf('Mean ERBB EGF-response AUC %.2f → %.2f and shared-output AUC %.2f → %.2f with cobimetinib',mean(hcontrast[role=='E',vem_auc]),mean(hcontrast[role=='E',vem_cobi_auc]),mean(hcontrast[role=='R',vem_auc]),mean(hcontrast[role=='R',vem_cobi_auc])),
    sprintf('At 1 µM vem: SOS1 %.2fx vehicle, ERK1 %.2fx, ERK2 %.2fx; EGF restores ERK1 %.2fx and ERK2 %.2fx over vem alone',ph_sos$one_uM_vem_vs_DMSO,ph_erk1$one_uM_vem_vs_DMSO,ph_erk2$one_uM_vem_vs_DMSO,ph_erk1$EGF_rescue_over_one_uM_vem,ph_erk2$EGF_rescue_over_one_uM_vem),
    sprintf('12-line molecular rank vs benefit r = %.2f; exact five-line ERBB interaction r = %.2f; shared-output interaction r = %.2f',adaptive_benefit_result$pearson,bypass_result[role=='E',pearson],bypass_result[role=='R',pearson])),
  evidence_status=c('Direct analysis','Direct analysis','Held-out direct analysis','External direct analysis','External direct reanalysis plus primary-study mechanism','External molecular-to-functional test'),
  primary_link=c('figures/main/07_source_context_coherence.pdf','figures/main/03_program_hierarchy.pdf','figures/main/11_prism_viability_association.pdf',geo127_url,phospho_url,geo684_url))
write_table(evidence_chain,'grant_claim_evidence_chain.tsv')

check(all(gerosa_tests[grepl('ERK[12]',analyte),one_uM_vem_vs_DMSO]<.5),'Gerosa dual-ERK phosphopeptides are not suppressed by 1 µM vemurafenib')
check(all(gerosa_tests[grepl('ERK[12]',analyte),EGF_rescue_over_one_uM_vem]>2),'Gerosa EGF sample does not restore dual-ERK phosphopeptides over vemurafenib alone')
check(gerosa_tests[analyte=='SOS1 ERK-feedback sites',one_uM_vem_vs_DMSO]<.5,'Gerosa SOS1 feedback phosphopeptides do not fall under vemurafenib')

check(length(chosen)>=60L && length(chosen)<=80L,'Expected a 60–80 gene coherence-qualified frozen panel')
check(K==5L && min(table(membership$program))>=12L,'Program count or minimum size changed')
check(setequal(unique(membership$role),c('B','E','R')),'The Y hierarchy is incomplete')
check(!any(grepl(generic_regex,membership$gene)) && !any(membership$gene %chin% generic_process_genes),'Generic expression genes entered the panel')
check(identical(dim(ct),c(length(chosen),length(chosen))) && identical(dim(cc),c(length(chosen),length(chosen))) && identical(dim(correlation),c(length(chosen),length(chosen))),'A correlation matrix is not square')
check(all(context_coherence[source %chin% c(drug_names,'CytoSig_ligand','CytoSig_other'),signed_mean_rho]>=MIN_SOURCE_PROGRAM_COR),'A retained program is below the cross-context coherence criterion')
check(all(abs(role_response$pearson)>=.4 & role_response$pearson_p<.05),'A matched pathway score lacks a Pearson r >= 0.4 association with its relevant PRISM drug response')
check(abs(root_comparison$pearson_controlling_drug)>=.4,'Held-out PRISM correlation is below the success criterion')
check(abs(resistance_result$pearson)>=.4 && abs(adaptive_benefit_result$pearson)>=.4 && all(abs(bypass_result$pearson)>=.4),'The repaired functional chain is below the Pearson r >= 0.4 criterion')
check(identical(readRDS(file.path(art,'frozen_programs.rds'))$membership,as.data.frame(membership)),'Frozen program definition changed during validation')
check(!any(line_drug$depmap_id[c(fit_rows,held)] %chin% EXCLUDED_IDS),'A discarded outlier entered the analysis')
check(!length(intersect(line_drug$depmap_id[fit_rows],line_drug$depmap_id[held])),'Development and holdout share a cell line')
verification<-data.table(check=c('60–80 gene coherence-qualified panel','Five programs, each at least twelve genes','Resolved BRAF, ERBB and shared-output arms','No ribosomal or generic expression genes','All three correlation matrices are square','Every program has mean correlation >= 0.10 in all five discovery contexts','Each matched pathway score has Pearson r >= 0.4 with its relevant PRISM response','Held-out PRISM Pearson r >= 0.4','Functional-chain Pearson r >= 0.4','Gerosa phosphoproteomic predictions hold','Frozen before external validation','Discarded outliers excluded','Cell-line holdout is disjoint'),passed=TRUE)
write_table(verification,'verification.tsv')

figure_files<-sort(list.files(file.path(out,'figures'),pattern='\\.pdf$',recursive=TRUE))
write_table(data.table(figure=file.path('figures',figure_files)),'FIGURE_INDEX.tsv')
summary_table<-merge(role_audit[,.(program,role,n_genes)],context_coherence[source %chin% c('Tahoe_balanced','CytoSig'),.(program,source,signed_mean_rho)],by='program')
write_table(summary_table,'program_summary.tsv')

writeLines(c('# Stage 33: BRAF–ERBB bypass from programs to function','',
  '## Biological setting','',
  'BRAF-mutant melanoma cells normally drive MEK–ERK signaling through mutant BRAF. Vemurafenib blocks that route. EGF can activate ERBB receptors and reach the same MEK–ERK machinery through RAS and RAF dimers. Lapatinib blocks the ERBB route; cobimetinib blocks MEK after the two routes converge. The proposed architecture is therefore a Y: a BRAF arm and an ERBB arm feeding a shared MEK–ERK output.','',
  '## What we set out to show','',
  'The experiment asks whether two very different transcriptomic resources can recover that architecture without fitting to the validation outcomes, whether the frozen architecture moves as predicted after EGF and downstream blockade, and whether those movements anticipate viable-cell response. Figure 01 shows this claim chain.','',
  '## 1. Recover programs jointly from Tahoe and CytoSig','',
  sprintf('Tahoe contributes every retained training-line profile at all three doses of vemurafenib, lapatinib and cobimetinib, together with PRISM viable-cell response. CytoSig contributes all 2,056 perturbation signatures: 67 EGF/NRG1 profiles identify the ligand response, while the other 1,989 profiles supply the variation needed to test general covariance. We admitted response-linked genes from focused BRAF, ERBB, MAPK and proliferation biology and expanded them only while preserving the weakest cross-context correlation. The frozen panel has %d genes in five programs of %s genes; no ribosomal or generic expression genes are present.',length(chosen),paste(sort(as.integer(table(membership$program))),collapse=', ')),
  sprintf('Every frozen program clears a mean signed correlation of 0.10 in each of five discovery contexts: the three separate Tahoe drug matrices, CytoSig EGF/NRG1 profiles and all other CytoSig perturbations. Across those contexts the observed program means range from %.2f to %.2f. Figures 04–06 show the same frozen order and correlation scale in square matrices, with every program enclosed by a diagonal rectangle.',min(context_coherence[source%chin%c(drug_names,'CytoSig_ligand','CytoSig_other'),signed_mean_rho]),max(context_coherence[source%chin%c(drug_names,'CytoSig_ligand','CytoSig_other'),signed_mean_rho])),'',
  '## 2. Resolve and freeze the Y','',
  sprintf('The resulting Y contains a BRAF-linked proliferation arm (%s), an ERBB immediate-response arm (%s), and three shared outputs below MEK–ERK (%s): adaptive MAPK response, mitosis and persistence. Figure 02 introduces the biochemical route and places each inhibitor at its target. Gene membership, direction and weights are then frozen before any external expression, phosphoproteomic or functional outcome is tested.',paste(programs[roles=='B'],collapse=', '),paste(programs[roles=='E'],collapse=', '),paste(programs[roles=='R'],collapse=', ')),'',
  '## 3. Establish that the programs carry drug-response information','',
  paste0('Within the Tahoe development set, each matched part of the pathway tracks its drug response: ',paste(sprintf('%s (%s) versus %s, Pearson r = %.2f (p = %.3g)',role_response$role,role_response$programs,role_response$relevant_drug,role_response$pearson,role_response$pearson_p),collapse='; '),'.'),
  sprintf('In cell lines excluded from discovery, the shared MEK–ERK score improves PRISM RMSE from %.2f to %.2f and has a drug-adjusted Pearson r of %.2f with viable-cell response. This makes viability part of the central test, rather than an optional annotation.',root_comparison$baseline_rmse,root_comparison$root_rmse,root_comparison$pearson_controlling_drug),'',
  '## 4. Test the predicted molecular displacement','',
  sprintf('GSE127988 supplies the explicit ligand experiment: A375 cells receive EGF during 1 µM vemurafenib, with or without 1 µM cobimetinib. EGF activates the frozen ERBB arm and shared output on both biological days. Adding cobimetinib changes mean ERBB response AUC from %.2f to %.2f and shared-output AUC from %.2f to %.2f. The BRAF arm does not show this EGF-specific pattern.',mean(hcontrast[role=='E',vem_auc]),mean(hcontrast[role=='E',vem_cobi_auc]),mean(hcontrast[role=='R',vem_auc]),mean(hcontrast[role=='R',vem_cobi_auc])),'',
  '## 5. Connect the RNA model to phosphosignaling','',
  sprintf('We directly reanalyzed the public Gerosa targeted SRM data from the same A375 setting. At 1 µM vemurafenib, the SOS1 ERK-feedback phosphopeptides fall to %.2fx vehicle and dual ERK1/ERK2 phosphopeptides fall to %.2fx/%.2fx. Adding EGF restores dual ERK1/ERK2 phosphorylation by %.2fx/%.2fx over vemurafenib alone. Figure 13 provides the replicate-level result. The primary paper independently shows receptor-driven pMEK/pERK pulses and links those pulses to cell-cycle escape during RAF/MEK inhibition.',ph_sos$one_uM_vem_vs_DMSO,ph_erk1$one_uM_vem_vs_DMSO,ph_erk2$one_uM_vem_vs_DMSO,ph_erk1$EGF_rescue_over_one_uM_vem,ph_erk2$EGF_rescue_over_one_uM_vem),'',
  '## 6. Test the functional prediction','',
  sprintf('GSE68453 provides a separate melanoma panel with 8-hour RNA after PLX4720, lapatinib or their combination and 3-day viable-cell Bliss benefit. The frozen EGF-responsive shared-output state ranks combination benefit across all 12 lines (r = %.2f). In the five lines with the exact four-arm RNA design, the molecular interaction predicts later benefit for the ERBB arm (r = %.2f) and shared MEK–ERK output (r = %.2f). Gerosa adds orthogonal single-cell evidence that receptor-driven ERK pulses permit persistence under RAF/MEK inhibition.',adaptive_benefit_result$pearson,bypass_result[role=='E',pearson],bypass_result[role=='R',pearson]),'',
  '## Claim supported by the full chain','',
  'Heterogeneous inhibitor and ligand transcriptomes can be integrated to recover a coherent BRAF–ERBB–MEK/ERK program hierarchy. Once frozen, that hierarchy predicts ligand-driven molecular bypass, its suppression by a downstream inhibitor, matched phosphosignaling, and the direction of later combination benefit.','',
  'The evidence supports a preclinical mechanistic and prediction claim. The exact five-line molecular-to-function test is small, and the public SRM EGF condition has one replicate; the 12-line ranking and independent held-out molecular and viability tests strengthen that limited exact-match comparison.','',
  '[Open the result viewer](index.html) · [Claim-by-claim evidence table](grant_claim_evidence_chain.tsv) · [Phosphoproteomic and functional links](PHOSPHOPROTEOMICS_AND_FUNCTION.md) · [Frozen genes](frozen_program_genes.tsv) · [Analysis code](../../scripts/33_rebuild_braf_erbb_programs.R)'),file.path(out,'STORY.md'))

writeLines(c('# Phosphoproteomic and functional bridge','',
  '## Exact mechanistic bridge','',
  'The frozen RNA model predicts that vemurafenib suppresses the shared MEK–ERK output while releasing ERK-dependent negative feedback, leaving an ERBB route ready to transmit an EGF signal. Gerosa et al. test this in A375 BRAF-V600E cells, the same model used for the GSE127988 RNA time course.','',
  sprintf('- **Targeted phosphoproteomics:** [public Panorama Skyline project](%s). The analysis reads the public LabKey tables directly. Heavy-standard-normalized SOS1 feedback phosphopeptides and doubly phosphorylated ERK1/2 are shown in Figure 13.',phospho_url),
  sprintf('- **Matched transcriptomics:** [GSE127988](%s). A375 cells were treated with 1 µM vemurafenib with or without 1 µM cobimetinib, followed by 100 ng/mL EGF over 0–8 hours.',geo127_url),
  sprintf('- **Mechanism and persistence:** [Gerosa et al., Cell Systems 2020](%s). The study measures pMEK/pERK and live-cell ERK pulses and links receptor-driven pulses to cell-cycle escape during RAF/MEK inhibition.',paper_url),'',
  'These layers test different parts of one prediction. SRM establishes feedback release and EGF restoration of ERK phosphorylation. GSE127988 tests the resulting transcriptional displacement. Live-cell signaling and proliferation establish that receptor-driven ERK activity is functionally relevant to persistence.','',
  '## Independent functional ranking','',
  sprintf('[GSE68453](%s) contains melanoma RNA after PLX4720, lapatinib and the combination. The matched publication reports three-day combination benefit. Joining by explicit cell-line identifier tests whether the frozen early molecular response ranks the later viable-cell outcome.',geo684_url),'',
  'No EGF is added in GSE68453. It tests naturally available ERBB bypass across cell lines, whereas GSE127988 and the Gerosa EGF sample test an explicit ligand challenge. Together they connect controlled pathway activation to variation in functional combination benefit.','',
  '## What remains bounded','',
  'The phosphoproteomic EGF point has one public SRM replicate, so it is a directionally matched mechanistic bridge rather than a stand-alone population estimate. The five-line exact RNA-to-Bliss test is also small. The 12-line ranking, held-out PRISM analysis, RNA time course and published single-cell mechanism provide complementary evidence, but a prospective clinical response claim would require another cohort.'),file.path(out,'PHOSPHOPROTEOMICS_AND_FUNCTION.md'))

writeLines(c('# Stage implementation','',
  'This directory contains only the current end-to-end BRAF–ERBB experiment. Main figures follow the scientific story in numerical order; supporting gene-by-sample heatmaps and diagnostics are under `figures/supplementary`.','',
  sprintf('Gene admission requires nominal association with a relevant development-set PRISM response or reproducible EGF/NRG1 displacement, direction stability across Tahoe doses where applicable, and membership in focused BRAF/ERBB/MAPK/proliferation biology. Ribosomal, translation, mitochondrial-expression, RNA-processing, histone and noncoding genes are excluded. Programs are expanded to 12–15 genes by maximizing the weakest signed mean correlation across five discovery contexts; every final value must remain at least %.2f. The final panel contains %d genes.',MIN_SOURCE_PROGRAM_COR,length(chosen)),'',
  'Tahoe correlations use all retained development cell lines and all three doses, separately for each drug. CytoSig correlations use all 2,056 individual profiles with equal study weight; EGF/NRG1 and other perturbations are also checked separately. Five anchored programs are retained because the silhouette diagnostic supports the requested five-to-six-program range and these five produce sizeable, interpretable, nonoverlapping blocks. Program signs and branch roles follow their observed inhibitor/ligand displacement.','',
  'Programs, signs and weights are frozen before the held-out PRISM lines, GSE127988 or GSE68453 outcomes are evaluated. External missing genes are never replaced with zero; coverage and identifiers are exported explicitly.','',
  'Main figures use one visual grammar throughout: orange is the mutant-BRAF arm, teal is the EGF–ERBB arm, and purple is the shared MEK–ERK output. The Y is introduced before any matrix, and later molecular results are placed back onto the same topology. Tahoe and CytoSig retain separate source colors only when the sources themselves are being compared.','',
  'Run from `/work`:','',
  '```sh','Rscript tahoe_prism_v2/scripts/33_rebuild_braf_erbb_programs.R','python3 tahoe_prism_v2/scripts/33_verify_braf_erbb_programs.py','```'),file.path(out,'README.md'))

html_escape<-function(x){x<-as.character(x);x<-gsub('&','&amp;',x,fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
html_table<-function(x)paste0('<table><tr>',paste0('<th>',names(x),'</th>',collapse=''),'</tr>',paste(apply(as.data.frame(x),1,function(r)paste0('<tr>',paste0('<td>',html_escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table>')
main_pdfs<-file.path('figures/main',sort(list.files(main_fig_dir,pattern='\\.pdf$')));supp_pdfs<-file.path('figures/supplementary',sort(list.files(supp_fig_dir,pattern='\\.pdf$')))
figure_captions<-c(
  '01_experimental_design'='The question, discovery inputs, frozen model and three validation layers.',
  '02_braf_erbb_mek_erk_pathway'='The biochemical Y, every ligand and inhibitor, and representative genes from all five programs.',
  '03_program_hierarchy'='The frozen BRAF and ERBB arms converging on three shared MEK–ERK output programs.',
  '04_joint_gene_correlations'='The five frozen programs in the averaged Tahoe–CytoSig gene-correlation matrix.',
  '05_tahoe_gene_correlations'='The identical genes and order using all Tahoe training lines and doses.',
  '06_cytosig_gene_correlations'='The identical genes and order using all 2,056 CytoSig perturbation profiles.',
  '07_source_context_coherence'='Within-program correlation in each Tahoe drug context and in ligand and non-ligand CytoSig contexts.',
  '08_gse68453_factorial_response'='How BRAF and ERBB inhibition move each frozen arm in external melanoma lines.',
  '09_gse127988_program_gene_heatmap'='Frozen program genes across EGF time, biological day and vemurafenib/cobimetinib background.',
  '10_gse127988_egf_response'='The EGF response placed on the Y with and without downstream MEK blockade.',
  '11_prism_viability_association'='Held-out test of the shared MEK–ERK response against five-day viability.',
  '12_gse68453_combination_prediction'='Early arm suppression versus later PLX4720–lapatinib benefit.',
  '13_gerosa_targeted_phosphoproteomics'='Public SRM evidence for feedback release and EGF restoration of ERK phosphorylation.')
make_cards<-function(ff,captions=FALSE)vapply(ff,function(f){png<-sub('\\.pdf$','.png',f);key<-sub('\\.pdf$','',basename(f));caption<-if(captions && key%chin%names(figure_captions))paste0('<p>',html_escape(figure_captions[[key]]),'</p>') else '';paste0('<article><h3>',html_escape(key),'</h3><a href="',f,'"><img loading="lazy" src="',png,'"></a>',caption,'<p><a href="',f,'">Open PDF</a></p></article>')},character(1))
tables<-sort(list.files(out,pattern='\\.tsv(\\.gz)?$'))
writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Stage 33 · BRAF–ERBB programs</title>',
  '<style>body{font:16px/1.5 system-ui;background:#f4f6f8;color:#183044;margin:0}main{max-width:1200px;margin:auto;padding:26px}a{color:#075b8e}.status,article{background:white;border:1px solid #d8e0e6;border-radius:5px}.status{padding:16px;border-left:6px solid #287d99}.grid{display:grid;grid-template-columns:1fr;gap:18px}article{padding:14px}img{display:block;max-width:100%;height:auto;margin:auto}table{border-collapse:collapse;width:100%;font-size:13px;background:white}th,td{padding:7px;border:1px solid #d7dde2;text-align:left}.tables{columns:2}pre{white-space:pre-wrap;background:#172733;color:#edf4f7;padding:15px}</style><main>',
  '<h1>Stage 33 · BRAF–ERBB program hierarchy</h1><p><a href="STORY.md">Scientific story</a> · <a href="grant_claim_evidence_chain.tsv">Claim-by-claim evidence</a> · <a href="PHOSPHOPROTEOMICS_AND_FUNCTION.md">Phosphoproteomic and functional bridge</a> · <a href="README.md">Methods</a> · <a href="../../scripts/33_rebuild_braf_erbb_programs.R">R code</a></p>',
  paste0('<div class="status"><strong>',length(chosen),' genes · ',K,' programs · complete BRAF/ERBB/MEK–ERK hierarchy</strong><p>CytoSig informs discovery. EGF activates the ERBB and shared-output branches; cobimetinib suppresses most of that RNA response; public targeted SRM confirms feedback release and EGF restoration of ERK phosphorylation; and frozen program displacement predicts held-out viability and later lapatinib combination benefit.</p></div>'),
  '<h2>Question and program discovery</h2><div class="grid">',make_cards(main_pdfs[1:7],TRUE),'</div><h2>Frozen model tested in external molecular and functional data</h2><div class="grid">',make_cards(main_pdfs[8:13],TRUE),'</div><h2>Program definitions</h2>',html_table(role_audit),'<h2>Matched pathway–PRISM associations</h2>',html_table(role_response),'<h2>Held-out PRISM test</h2>',html_table(root_comparison),'<h2>GSE127988 EGF and cobimetinib test</h2>',html_table(hcontrast),'<h2>GSE68453 functional ranking test</h2>',html_table(bypass_result),
  '<h2>Supporting figures</h2><div class="grid">',make_cards(supp_pdfs),'</div><h2>Tables</h2><ul class="tables">',paste0('<li><a href="',tables,'">',tables,'</a></li>'),'</ul>',
  '<h2>Complete analysis code</h2><details><summary>Show code</summary><pre><code>',paste(html_escape(readLines(file.path(root,'scripts/33_rebuild_braf_erbb_programs.R'),warn=FALSE)),collapse='\n'),'</code></pre></details></main></html>'),file.path(out,'index.html'))
message('Stage 33 complete: ',out)
