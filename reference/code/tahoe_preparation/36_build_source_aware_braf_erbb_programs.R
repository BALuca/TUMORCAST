#!/usr/bin/env Rscript
## Stage 36: source-aware RAF/ERBB/MEK programs for Aim 1.
## A source-specific RAF module is allowed to rely on Tahoe because CytoSig has
## no matched RAF-inhibitor perturbation. Shared modules must reproduce across
## Tahoe and CytoSig. Stages 34 and 35 remain immutable reference analyses.
suppressPackageStartupMessages({library(data.table);library(matrixStats);library(ggplot2);library(ggpubr);library(cluster)})
root<-'/work/tahoe_prism_v2';work<-dirname(root);input_art<-file.path(root,'artifacts/26_braf_erbb_bypass');prepared_art<-file.path(root,'artifacts/27_braf_erbb_hierarchy')
out<-file.path(root,'results/36_source_aware_braf_erbb_programs');art<-file.path(root,'artifacts/36_source_aware_braf_erbb_programs')
unlink(out,recursive=TRUE,force=TRUE);unlink(art,recursive=TRUE,force=TRUE)
dir.create(file.path(out,'figures/main'),recursive=TRUE);dir.create(file.path(out,'figures/supplementary'),recursive=TRUE);dir.create(art,recursive=TRUE)
source(file.path(work,'src/R/lib_heatmap.R'))
source(file.path(root,'src/R/plot_braf_erbb_pathway.R'))
main<-file.path(out,'figures/main');supp<-file.path(out,'figures/supplementary')
write_table<-function(x,n)fwrite(as.data.table(x),file.path(out,n),sep='\t',na='NA')
robust_limit<-function(x,p=.98)max(as.numeric(quantile(abs(x[is.finite(x)]),p,na.rm=TRUE)),.01)
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
read_npy_rows<-function(path,rows_zero_based){con<-file(path,'rb');on.exit(close(con));magic<-rawToChar(readBin(con,'raw',6L));version<-readBin(con,'integer',2L,size=1L,signed=FALSE);stopifnot(identical(magic,'\x93NUMPY'),version[1]%in%1:3);len<-readBin(con,'integer',1L,size=if(version[1]==1L)2L else 4L,signed=FALSE,endian='little');header<-rawToChar(readBin(con,'raw',len));shape<-as.integer(trimws(strsplit(sub(".*'shape': \\(([^)]*)\\).*",'\\1',header),',')[[1L]]));shape<-shape[!is.na(shape)];offset<-seek(con);answer<-matrix(NA_real_,length(rows_zero_based),shape[2]);for(i in seq_along(rows_zero_based)){seek(con,offset+4*shape[2]*rows_zero_based[i],origin='start');answer[i,]<-readBin(con,'numeric',n=shape[2],size=4L,endian='little')};answer}
normalize_pair<-function(treated,reference){treated<-pmax(treated,0);reference<-pmax(reference,0);tt<-sum(treated);rt<-sum(reference);tp<-treated/tt;rp<-reference/rt;eligible<-tp*10000>=.1&rp*10000>=.1;m<-log2(tp[eligible]/rp[eligible]);aa<-.5*log2(tp[eligible]*rp[eligible]);mq<-quantile(m,c(.3,.7));aq<-quantile(aa,c(.05,.95));kept<-m>=mq[1]&m<=mq[2]&aa>=aq[1]&aa<=aq[2];fac<-2^median(m[kept]);log2(1+tp*10000/fac)-log2(1+rp*10000)}
e<-fread(file.path(root,'results/33_braf_erbb_programs/gene_selection_evidence.tsv.gz'))
old<-fread(file.path(root,'results/34_aim1_cross_source_programs/frozen_program_genes.tsv'))
raf_pool<-e[vem_eligible==TRUE & abs(rho_vemurafenib)-pmax(abs(rho_lapatinib),abs(rho_cobimetinib))>=.15 & grepl('RAF|ERK|KRAS|MELANOMA|ERBB|EGFR',pathways),gene]
raf_pool<-unique(c('ARAF','RAF1','HRAS','IGF1R','KIT','RALBP1','GSK3B','PPP2R1B',raf_pool))
candidates<-unique(c(raf_pool,old$gene))
a<-readRDS(file.path(input_art,'discovery_measurements.rds'));idx<-match(candidates,colnames(a$effects));check(!anyNA(idx),'candidate absent')
effects<-a$effects[,idx,drop=FALSE];colnames(effects)<-candidates;meta<-a$meta;cells<-a$cells
new_meta<-fread(file.path(prepared_art,'additional_response_metadata.tsv'));nt<-read_npy_rows(file.path(prepared_art,'additional_treated_mean.npy'),new_meta$array_row);nr<-read_npy_rows(file.path(prepared_art,'additional_reference_mean.npy'),new_meta$array_row)
ne<-t(vapply(seq_len(nrow(nt)),function(i)normalize_pair(nt[i,],nr[i,]),numeric(ncol(nt))))[,idx,drop=FALSE];colnames(ne)<-candidates
new_cells<-unique(new_meta[,.(depmap_id,cell_line,cell_name)]);model<-fread(file.path(work,'outputs/dataset_inventory/cache/depmap/Model.csv'));new_cells<-merge(new_cells,model[,.(depmap_id=ModelID,rrid=RRID,tissue=OncotreeLineage,disease=OncotreePrimaryDisease)],by='depmap_id');new_cells[,`:=`(braf_v600_confirmed=TRUE,split='Development')]
cells<-rbindlist(list(cells,new_cells),fill=TRUE);cells[,`:=`(analysis_included=!depmap_id%chin%c('ACH-000311','ACH-000628'),melanoma_v600=braf_v600_confirmed&tissue=='Skin')]
new_meta[,`:=`(compound=c('Vemurafenib'='vemurafenib','Lapatinib ditosylate'='lapatinib','Cobimetinib'='cobimetinib')[drug],perturbation_id=paste0(drug,'__',dose,'uM'))];new_meta[,observation_id:=paste(depmap_id,perturbation_id,sep='__')];rownames(ne)<-new_meta$observation_id
effects<-rbind(effects,ne);meta<-rbindlist(list(meta[,.(depmap_id,compound,drug,dose,dose_unit,perturbation_id,observation_id,n_plates)],new_meta[,.(depmap_id,compound,drug,dose,dose_unit,perturbation_id,observation_id,n_plates)]));meta<-merge(meta,cells[,.(depmap_id,cell_name,tissue,split,braf_v600_confirmed,melanoma_v600,analysis_included)],by='depmap_id');setorder(meta,compound,depmap_id,dose);effects<-effects[meta$observation_id,,drop=FALSE]
prism<-fread(file.path(prepared_art,'primary_prism_responses.tsv'));prism[passed_str_profiling==FALSE|is.na(passed_str_profiling),`:=`(prism_survival=NA_real_,prism_sensitivity=NA_real_)];meta<-merge(meta,prism[,.(depmap_id,compound,prism_survival,prism_sensitivity,functional_dose_uM,passed_str_profiling)],by=c('depmap_id','compound'),all.x=TRUE);setorder(meta,compound,depmap_id,dose)
bm<-fread(file.path(prepared_art,'plate_common_response_metadata.tsv'));bg<-read_npy_rows(file.path(prepared_art,'plate_common_response.npy'),bm$array_row)[,idx,drop=FALSE];pm<-rbindlist(list(fread(file.path(work,'cache/tahoe_prism_full_chain/tahoe_pdex_selected/pdex_plate_audit.tsv.gz'))[,.(depmap_id,drug,dose,plate)],fread(file.path(prepared_art,'additional_plate_audit.tsv'))[,.(depmap_id,drug,dose,plate)]))
for(i in which(meta$analysis_included)){pp<-pm[depmap_id==meta$depmap_id[i]&drug==meta$drug[i]&dose==meta$dose[i],unique(plate)];jj<-which(bm$depmap_id==meta$depmap_id[i]&bm$plate%in%pp);effects[i,]<-effects[i,]-colMeans(bg[jj,,drop=FALSE])}
meta[,dose_weight:=c(.25,.5,.25)[match(dose,c(.05,.5,5))]];ld<-unique(meta[,.(depmap_id,compound,cell_name,tissue,split,braf_v600_confirmed,melanoma_v600,analysis_included,prism_survival,prism_sensitivity)]);ld[,`:=`(key=paste(depmap_id,compound,sep='__'),eligible=is.finite(prism_survival))];setorder(ld,compound,depmap_id)
auc<-matrix(NA_real_,nrow(ld),length(candidates),dimnames=list(ld$key,candidates));for(i in seq_len(nrow(ld))){ii<-which(meta$depmap_id==ld$depmap_id[i]&meta$compound==ld$compound[i]);auc[i,]<-colSums(effects[ii,,drop=FALSE]*meta$dose_weight[ii])}
## 2. Freeze four source-aware biological programs -----------------------
# The RAF network is chosen from nominally response-associated genes with a
# recognizable RAF/receptor/adaptor role. It is required to be coherent in
# Tahoe vemurafenib profiles and to retain ARAF and RAF1. CytoSig is not used
# as a gate because it contains ligand perturbations, not RAF inhibition.
program_genes<-list(
  P1=c('ARAF','RAF1','HRAS','IGF1R','KIT','RALBP1','PTK2B','ITCH','STAM2','TNFAIP3'),
  P2=setdiff(old[program=='P2',gene],c('SHC4','KCNN4')),
  P3=old[program=='P3',gene],
  P4=c('NUMA1','TFRC','HAUS5','KPNA2','PLK1','CDC25C','SPAG5','TACC3','BUB1B','KIF11','TPX2','DLGAP5','DTYMK','TUBG1','TUBA1C','TUBB4B')
)
programs<-names(program_genes)
program_names<-c(P1='RAF signaling network',P2='EGF-ERBB adaptive response',P3='MEK-ERK shared output',P4='Proliferation output')
program_drug<-c(P1='vemurafenib',P2='lapatinib',P3='cobimetinib',P4='vemurafenib')
program_colors<-c(P1='#E41A1C',P2='#377EB8',P3='#4DAF4A',P4='#984EA3')
program_fills<-c(P1='#F6B7B4',P2='#B8D0EA',P3='#B9DFB8',P4='#D8C1E5')
drug_colors<-c(vemurafenib='#D88736',lapatinib='#168C86',cobimetinib='#7656A5')
drug_labels<-c(vemurafenib='Vemurafenib (BRAF inhibitor)',lapatinib='Lapatinib (ERBB inhibitor)',cobimetinib='Cobimetinib (MEK inhibitor)')
keep<-unlist(program_genes,use.names=FALSE);check(length(keep)==58L&&!anyDuplicated(keep),'Programs must contain 58 unique genes')
check(all(keep%chin%candidates),'A frozen gene is absent from the prepared Tahoe matrix')

gene_sign<-setNames(old$orientation,old$gene)
gene_sign[program_genes$P1]<--sign(e[match(program_genes$P1,gene),rho_vemurafenib])
check(all(gene_sign[keep]%in%c(-1,1)),'A gene lacks a signed orientation')

## 3. Correlations use individual profiles from both sources -------------
drug_names<-c('vemurafenib','lapatinib','cobimetinib')
tahoe_raw<-lapply(drug_names,function(d){ii<-which(meta$compound==d&meta$split=='Development'&meta$analysis_included);cor(effects[ii,keep,drop=FALSE],method='spearman',use='pairwise.complete.obs')});names(tahoe_raw)<-drug_names
tahoe_all_raw<-Reduce('+',tahoe_raw)/length(tahoe_raw)

cyto_path<-file.path(work,'data/01_CytoSig_data/diff.merge.gz');cyto_meta_path<-file.path(work,'data/01_CytoSig_data/meta_info.txt')
cyto_dt<-fread(cyto_path,showProgress=FALSE);setnames(cyto_dt,1L,'gene');cyto_dt<-cyto_dt[gene%chin%keep]
cyto_x<-as.matrix(cyto_dt[match(keep,gene),-1]);rownames(cyto_x)<-keep
cm0<-fread(cyto_meta_path,skip=1L,header=FALSE,col.names=c('profile_id','treatment','condition','duration','dose'))
cm<-merge(data.table(profile_id=colnames(cyto_x)),cm0,by='profile_id',all.x=TRUE,sort=FALSE)[match(colnames(cyto_x),profile_id)]
cm[,study:=sub('.*@((?:GSE[0-9]+|E-MTAB-[0-9]+)).*','\\1',profile_id,perl=TRUE)]
cm[!grepl('^(GSE[0-9]+|E-MTAB-[0-9]+)$',study),study:=paste0('independent_',.I)]
study_n<-cm[,.N,by=study];study_size<-setNames(study_n$N,study_n$study);cyto_weight<-1/study_size[cm$study];cyto_weight<-cyto_weight/sum(cyto_weight)
weighted_spearman<-function(x,weights){ranks<-t(apply(x,1,rank,na.last='keep',ties.method='average'));p<-nrow(ranks);answer<-diag(1,p);dimnames(answer)<-list(rownames(ranks),rownames(ranks));if(p>1L)for(i in seq_len(p-1L))for(j in (i+1L):p){ok<-is.finite(ranks[i,])&is.finite(ranks[j,]);w<-weights[ok];w<-w/sum(w);a<-ranks[i,ok];b<-ranks[j,ok];a<-a-sum(w*a);b<-b-sum(w*b);answer[i,j]<-answer[j,i]<-sum(w*a*b)/sqrt(sum(w*a*a)*sum(w*b*b))};answer}
cyto_raw<-weighted_spearman(cyto_x,cyto_weight)
orient<-outer(gene_sign[keep],gene_sign[keep]);dimnames(orient)<-list(keep,keep)
tahoe_signed<-tahoe_all_raw*orient;cyto_signed<-cyto_raw*orient;tahoe_vem_signed<-tahoe_raw$vemurafenib*orient

# Source-aware square matrix: every pair involving the RAF module uses the
# matched Tahoe-vemurafenib covariance; all other pairs use the equal-source
# average. This avoids manufacturing absent RAF evidence in CytoSig.
source_aware<-0.5*(tahoe_signed+cyto_signed);raf_pair<-outer(keep%chin%program_genes$P1,keep%chin%program_genes$P1,'|');source_aware[raf_pair]<-tahoe_vem_signed[raf_pair];diag(source_aware)<-1

membership<-rbindlist(lapply(programs,function(p)data.table(gene=program_genes[[p]],program=p,program_name=program_names[p],orientation=gene_sign[program_genes[[p]]],support_rule=if(p=='P1')'Tahoe vemurafenib' else 'Tahoe + CytoSig')))
leaf_order<-unlist(lapply(programs,function(p){g<-membership[program==p,gene];m<-if(p=='P1')tahoe_vem_signed[g,g,drop=FALSE] else source_aware[g,g,drop=FALSE];g[hclust(as.dist(pmax(1-m,0)),method='average')$order]}),use.names=FALSE)
membership<-membership[match(leaf_order,gene)];membership[,plotting_order:=seq_len(.N)]
write_table(membership,'frozen_program_genes.tsv')

mean_upper<-function(m,g)mean(m[g,g,drop=FALSE][upper.tri(m[g,g,drop=FALSE])])
coherence<-rbindlist(lapply(programs,function(p){g<-program_genes[[p]];data.table(program=p,program_name=program_names[p],source=c('Tahoe','CytoSig'),profiles=c(if(p=='P1')84L else 252L,2056L),mean_signed_rho=c(if(p=='P1')mean_upper(tahoe_vem_signed,g) else mean_upper(tahoe_signed,g),mean_upper(cyto_signed,g)),required_for_admission=c(TRUE,p!='P1'),interpretation=c(if(p=='P1')'Matched vemurafenib context' else 'All three drug contexts',if(p=='P1')'Descriptive only: no RAF-inhibitor perturbation in CytoSig' else 'Required cross-source replication'))}))
check(all(coherence[required_for_admission==TRUE,mean_signed_rho]>=.10),'A required source-program coherence is below 0.10')
write_table(coherence,'program_source_coherence.tsv')
write_table(data.table(source=c('Tahoe, vemurafenib','Tahoe, lapatinib','Tahoe, cobimetinib','CytoSig','Source-aware joint'),individual_profiles=c(84L,84L,84L,2056L,NA_integer_),correlation_unit=c(rep('Individual cell-line by dose RNA profile',3),'Individual control-normalized perturbation signature','Gene pair'),method=c(rep('Spearman',3),'Study-balanced weighted Spearman','Tahoe-vemurafenib for RAF pairs; equal-source average otherwise'),treatment_medians_used=FALSE),'correlation_method_audit.tsv')

## 4. Each program must track its matched drug response ------------------
program_score<-function(rows,p)rowMeans(sweep(auc[rows,program_genes[[p]],drop=FALSE],2,gene_sign[program_genes[[p]]],'*'))
assoc<-rbindlist(lapply(programs,function(p){d<-program_drug[p];ii<-which(ld$compound==d&ld$split=='Development'&ld$analysis_included&ld$eligible);x<-program_score(ii,p);tt<-cor.test(x,ld$prism_survival[ii]);data.table(program=p,program_name=program_names[p],drug=d,n=length(ii),pearson=unname(tt$estimate),p_value=tt$p.value)}))
check(all(assoc$pearson>=.4&assoc$p_value<.05),'A matched development association fails the prespecified success criterion')
write_table(assoc,'matched_program_training_viability.tsv')

biology<-data.table(program=programs,program_name=program_names[programs],biological_role=c('RAF-family signaling and receptor/adaptor control','ERBB-induced immediate-early and feedback transcription','MEK-ERK signaling, feedback and signal routing','Mitotic and cell-cycle execution'),representative_genes=c('ARAF, RAF1, HRAS, IGF1R, KIT, RALBP1','EREG, HBEGF, FOS, JUN, DUSP4, DUSP6','MAPK1, SPRY2, ERRFI1, PRKCA, PPP2R1B','PLK1, CDC25C, BUB1B, KIF11, TPX2, DLGAP5'),n_genes=lengths(program_genes),source_logic=c('Tahoe-specific evidence retained','Joint Tahoe and CytoSig evidence','Joint Tahoe and CytoSig evidence','Joint Tahoe and CytoSig evidence'))
write_table(biology,'program_biological_content.tsv')

## 5. Freeze, then predict viability in disjoint cell lines --------------
# Only one program per drug enters the general holdout: distal proliferation
# for vemurafenib, ERBB adaptive response for lapatinib and MEK-ERK output for
# cobimetinib. The RAF module is tested in a BRAF-specific external setting.
held_program_drug<-c(P2='lapatinib',P3='cobimetinib',P4='vemurafenib')
scored<-rbindlist(lapply(seq_len(nrow(ld)),function(i){p<-names(held_program_drug)[held_program_drug==ld$compound[i]];data.table(row=i,program=p,program_score=program_score(i,p))}))
scored<-cbind(copy(ld),scored[,.(program,program_score)]);fit<-scored[split=='Development'&analysis_included&eligible];held<-scored[split=='Holdout'&analysis_included&eligible]
baseline<-lm(prism_survival~compound,data=fit);full_model<-lm(prism_survival~compound+program_score,data=fit)
held[,`:=`(baseline_prediction=predict(baseline,newdata=held),program_prediction=predict(full_model,newdata=held))]
held[,`:=`(within_drug_score=program_score-mean(program_score),within_drug_viability=prism_survival-mean(prism_survival)),by=compound]
ht<-cor.test(held$within_drug_score,held$within_drug_viability)
comparison<-data.table(n_lines=uniqueN(held$depmap_id),n_observations=nrow(held),pearson=unname(ht$estimate),p_value=ht$p.value,baseline_rmse=sqrt(mean((held$prism_survival-held$baseline_prediction)^2)),program_rmse=sqrt(mean((held$prism_survival-held$program_prediction)^2)))
check(comparison$pearson>=.4&&comparison$program_rmse<comparison$baseline_rmse,'Source-aware programs fail held-out viability')
write_table(held,'heldout_prism_predictions.tsv');write_table(comparison,'heldout_prism_test.tsv')

saveRDS(list(membership=membership,program_names=program_names,program_drug=program_drug,held_program_drug=held_program_drug,gene_sign=gene_sign[keep],leaf_order=leaf_order,training_model=full_model,source_aware_correlation=source_aware[leaf_order,leaf_order],tahoe_correlation=tahoe_signed[leaf_order,leaf_order],cytosig_correlation=cyto_signed[leaf_order,leaf_order]),file.path(art,'frozen_source_aware_programs.rds'))

## 6. Grant heatmaps ------------------------------------------------------
draw_correlation<-function(m,file,title,labels=FALSE){m<-m[leaf_order,leaf_order,drop=FALSE];lim<-max(abs(m[upper.tri(m)]));hm_id<-paste0('hm_',gsub('[^A-Za-z0-9]','_',file));ann<-data.frame(Program=factor(membership$program,levels=programs));h<-draw_heatmap(m,row_anno_df=ann,row_anno_cols='Program',palette_overrides=list(Program=program_colors),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(m)),column_order=seq_len(ncol(m)),show_row_names=labels,show_column_names=labels,row_names_gp=grid::gpar(fontsize=5.4),column_names_gp=grid::gpar(fontsize=5.4),column_names_rot=90,hm_name=hm_id,heatmap_legend_param=list(title='Signed\nSpearman ρ',at=round(c(-lim,0,lim),2)),title=title,width=grid::unit(if(labels)4.8 else 2.7,'in'),height=grid::unit(if(labels)4.8 else 2.7,'in'));add_boxes<-function(){rr<-rle(as.character(membership$program));en<-cumsum(rr$lengths);st<-c(1,head(en,-1)+1);n<-nrow(m);ComplexHeatmap::decorate_heatmap_body(hm_id,{for(i in seq_along(st))grid::grid.rect(x=grid::unit((st[i]-1)/n,'npc'),y=grid::unit(1-en[i]/n,'npc'),width=grid::unit(rr$lengths[i]/n,'npc'),height=grid::unit(rr$lengths[i]/n,'npc'),just=c('left','bottom'),gp=grid::gpar(fill=NA,col='#111111',lwd=1.7))})};base<-file.path(if(labels)supp else main,file);wh<-if(labels)7.0 else 4.15;render<-function(){ComplexHeatmap::draw(h,heatmap_legend_side='right');add_boxes()};cairo_pdf(paste0(base,'.pdf'),wh,wh);render();dev.off();png(paste0(base,'.png'),wh,wh,units='in',res=190,type='cairo');render();dev.off()}
draw_correlation(source_aware,'03_source_aware_gene_correlations','Source-aware gene correlations')
draw_correlation(tahoe_signed,'04_tahoe_gene_correlations','Tahoe, all 252 training profiles')
draw_correlation(cyto_signed,'05_cytosig_gene_correlations','CytoSig, all 2,056 signatures')
draw_correlation(source_aware,'S01_source_aware_gene_correlations_labeled','Source-aware correlations, genes labeled',TRUE)

tahoe_idx<-which(meta$split=='Development'&meta$analysis_included);tahoe_x<-t(effects[tahoe_idx,keep,drop=FALSE]);colnames(tahoe_x)<-meta$observation_id[tahoe_idx]
tahoe_meta<-meta[tahoe_idx,.(Drug=factor(compound,levels=drug_names),Dose=factor(paste0(dose,' µM'),levels=c('0.05 µM','0.5 µM','5 µM')))]
cm[,`ERBB ligand status`:=factor(ifelse(treatment%chin%c('EGF','NRG1'),'EGF / NRG1','Other perturbation'),levels=c('EGF / NRG1','Other perturbation'))]
check(ncol(tahoe_x)==252L&&ncol(cyto_x)==2056L,'Full profile heatmaps lost samples')

source_heat<-function(x,md,file,title,mode){x<-x[leaf_order,,drop=FALSE];z<-t(scale(t(x)));z[!is.finite(z)]<-0;signed<-z*gene_sign[leaf_order];activity<-sapply(programs,function(p)colMeans(signed[membership[program==p,gene],,drop=FALSE]));win<-max.col(activity,ties.method='first');set.seed(3601);rnd<-runif(ncol(x));ord<-if(mode=='Tahoe')order(win,md$Dose,rnd) else order(win,-activity[cbind(seq_len(ncol(x)),win)],rnd);cb<-factor(programs[win[ord]],levels=programs,labels=paste0(programs,', ',program_names[programs]));annr<-data.frame(Program=factor(membership$program,levels=programs));lim<-robust_limit(signed,.98);key<-c('RAF1','HRAS','FOS','DUSP6','HBEGF','MAPK1','ERRFI1','PLK1','TPX2');row_labels<-ifelse(leaf_order%chin%key,leaf_order,'');row_group<-factor(membership$program,levels=programs,labels=paste0(programs,', ',program_names[programs]));pal<-if(mode=='Tahoe')list(Program=program_colors,Drug=drug_colors,Dose=c('0.05 µM'='#E6EEF2','0.5 µM'='#91AFC0','5 µM'='#315D75')) else list(Program=program_colors,'ERBB ligand status'=c('EGF / NRG1'='#168C86','Other perturbation'='#C7CED2'));h<-draw_heatmap(signed[,ord,drop=FALSE],col_anno_df=as.data.frame(md[ord],check.names=FALSE),col_anno_cols=names(md),row_anno_df=annr,row_anno_cols='Program',palette_overrides=pal,ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(signed)),column_order=seq_len(ncol(signed)),row_split=row_group,column_split=cb,cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(.75,'mm'),column_gap=grid::unit(.8,'mm'),row_title=NULL,show_row_names=TRUE,show_column_names=FALSE,row_labels=row_labels,row_names_gp=grid::gpar(fontsize=8.2,fontface='bold'),hm_name=paste0('hm_',file),heatmap_legend_param=list(title='Signed row z-score',at=round(c(-lim,0,lim),1)),title=programs);draw1<-function()ComplexHeatmap::draw(h,heatmap_legend_side='right',annotation_legend_side='right',column_title=title,column_title_gp=grid::gpar(fontsize=13,fontface='bold'));base<-file.path(main,file);cairo_pdf(paste0(base,'.pdf'),7.2,5.1);draw1();dev.off();png(paste0(base,'.png'),7.2,5.1,units='in',res=190,type='cairo');draw1();dev.off();write_table(data.table(sample=colnames(x)[ord],plotting_order=seq_along(ord),dominant_program=cb,dominant_activity=activity[cbind(ord,win[ord])],as.data.table(md[ord])) ,paste0(file,'_sample_order.tsv'))}
source_heat(tahoe_x,tahoe_meta,'06_tahoe_all_profiles_gene_heatmap','Tahoe, all 252 individual training profiles','Tahoe')
source_heat(cyto_x,cm[,.(`ERBB ligand status`=`ERBB ligand status`)],'07_cytosig_all_profiles_gene_heatmap','CytoSig, strongest program activation first','CytoSig')

## 7. Compact ggplot2 figures --------------------------------------------
# Biological mechanism only: exact nodes and drug targets; inferred programs are separate.
p<-build_braf_erbb_pathway(show_plx4720=FALSE)
ggsave(file.path(main,'01_biological_pathway_and_drugs.pdf'),p,width=7.4,height=3.75,device=cairo_pdf,bg='white')
ggsave(file.path(main,'01_biological_pathway_and_drugs.png'),p,width=7.4,height=3.75,dpi=240,bg='white')

bio_plot<-copy(biology);bio_plot[,label:=factor(program,levels=rev(programs))];bio_plot[,detail:=paste0(program,', ',program_name,'\n',biological_role,'\n',representative_genes)]
p<-ggplot(bio_plot,aes(x=0,y=label,color=program))+geom_point(size=5)+geom_text(aes(x=.08,label=detail),hjust=0,size=3.35,lineheight=.93,color='#263840')+scale_color_manual(values=program_colors,guide='none')+coord_cartesian(xlim=c(-.05,1.22),clip='off')+theme_void(base_size=13)+theme(plot.title=element_text(size=15),plot.subtitle=element_text(size=10),plot.margin=margin(7,16,7,7))+labs(title='Each program captures biology at its pathway level',subtitle='RAF-specific genes use Tahoe; joint programs replicate in CytoSig')
ggsave(file.path(main,'02_program_biological_content.pdf'),p,width=7.2,height=4.1,device=cairo_pdf,bg='white');ggsave(file.path(main,'02_program_biological_content.png'),p,width=7.2,height=4.1,dpi=190,bg='white')

short<-c(P1='P1, RAF network',P2='P2, ERBB response',P3='P3, MEK-ERK output',P4='P4, Proliferation');train_plot<-rbindlist(lapply(programs,function(p){d<-program_drug[p];ii<-which(ld$compound==d&ld$split=='Development'&ld$analysis_included&ld$eligible);data.table(program=factor(p,levels=programs,labels=short[programs]),program_id=p,score=program_score(ii,p),viability=ld$prism_survival[ii])}))
p<-ggplot(train_plot,aes(score,viability,color=program_id))+geom_point(size=1.8)+geom_smooth(method='lm',formula=y~x,se=FALSE,linewidth=.8)+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=3.7,color='black')+facet_wrap(~program,nrow=1,scales='free')+scale_color_manual(values=program_colors,guide='none')+theme_classic(base_size=12)+theme(aspect.ratio=1,strip.text=element_text(face='bold',size=10.5),plot.title=element_text(size=15),plot.subtitle=element_text(size=10))+labs(x='Signed program RNA response',y='Five-day PRISM viability',title='Each pathway-level program tracks its matched inhibitor',subtitle='Development lines, axes are free within each program')
ggsave(file.path(main,'08_training_program_viability.pdf'),p,width=8.4,height=3.05,device=cairo_pdf,bg='white');ggsave(file.path(main,'08_training_program_viability.png'),p,width=8.4,height=3.05,dpi=190,bg='white')

p<-ggplot(held,aes(within_drug_score,within_drug_viability,color=compound))+geom_point(size=2.7)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#263840',linewidth=.9)+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=4.3,color='black')+scale_color_manual(values=drug_colors,breaks=names(drug_labels),labels=drug_labels)+theme_classic(base_size=13)+theme(aspect.ratio=1,legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=9),plot.title=element_text(size=13),plot.subtitle=element_text(size=9.5))+guides(color=guide_legend(nrow=3,byrow=TRUE))+labs(x='Matched program RNA response\n(centered within inhibitor)',y='Five-day PRISM viability\n(centered within inhibitor)',title='Frozen programs predict held-out viability',subtitle=sprintf('6 disjoint lines, Pearson r = %.2f, RMSE %.2f to %.2f',comparison$pearson,comparison$baseline_rmse,comparison$program_rmse))
ggsave(file.path(main,'09_heldout_viability_prediction.pdf'),p,width=5.0,height=5.5,device=cairo_pdf,bg='white');ggsave(file.path(main,'09_heldout_viability_prediction.png'),p,width=5.0,height=5.5,dpi=190,bg='white')

## 8. Verification and viewer -------------------------------------------
verification<-data.table(check=c('58 unique genes in four sizeable programs','ARAF and RAF1 retained','No ribosomal genes','All 252 Tahoe profiles used','All 2,056 CytoSig signatures used','Every required source coherence >= 0.10','Every training Pearson >= 0.40','Heldout Pearson >= 0.40','Heldout RMSE improves','No discovery/holdout line overlap','Source-aware correlation matrix square and symmetric'),passed=c(length(keep)==58L&&min(lengths(program_genes))>=10L,all(c('ARAF','RAF1')%chin%keep),!any(grepl('^RP[SL][0-9]',keep)),ncol(tahoe_x)==252L,ncol(cyto_x)==2056L,all(coherence[required_for_admission==TRUE,mean_signed_rho]>=.10),all(assoc$pearson>=.4),comparison$pearson>=.4,comparison$program_rmse<comparison$baseline_rmse,!length(intersect(fit$depmap_id,held$depmap_id)),nrow(source_aware)==ncol(source_aware)&&max(abs(source_aware-t(source_aware)))<1e-10))
check(all(verification$passed),'Stage 36 verification failed');write_table(verification,'verification.tsv')
figs<-sort(list.files(main,pattern='\\.pdf$'));write_table(data.table(figure=file.path('figures/main',figs)),'FIGURE_INDEX.tsv')
caption<-c('01_biological_pathway_and_drugs'='The biological pathway alone: BRAF and EGF-ERBB converge through MEK-ERK to proliferation, with each inhibitor at its target.','02_program_biological_content'='The separate data-derived representation names each program, its biological role and representative genes.','03_source_aware_gene_correlations'='Square signed-correlation matrix: Tahoe vemurafenib supports RAF pairs; the other blocks use the equal-source average.','04_tahoe_gene_correlations'='The identical genes and order across all 252 individual Tahoe development profiles.','05_cytosig_gene_correlations'='The identical genes and order across all 2,056 individual CytoSig signatures; RAF coherence is descriptive because CytoSig lacks RAF inhibition.','06_tahoe_all_profiles_gene_heatmap'='All Tahoe profiles; columns are grouped by dominant program, then dose, with deterministic within-dose order.','07_cytosig_all_profiles_gene_heatmap'='All CytoSig signatures; the strongest enrichment is first within each dominant-program block.','08_training_program_viability'='Each program is associated with viability under its matched inhibitor in development lines; every facet has free axes.','09_heldout_viability_prediction'='One matched program per inhibitor predicts viability in six disjoint cell lines; full inhibitor names and targets are shown.')
writeLines(c('# Source-aware RAF-ERBB program pilot','',sprintf('The new definition freezes %d genes in four programs of %s genes. ARAF and RAF1 are retained using matched Tahoe vemurafenib evidence; CytoSig replication is required only where it has relevant perturbational information.',length(keep),paste(lengths(program_genes),collapse=', ')),'',sprintf('Matched development correlations are %s.',paste(sprintf('%s/%s r = %.2f',assoc$program,assoc$drug,assoc$pearson),collapse='; ')),sprintf('The pre-frozen matched score predicts viability in six disjoint held-out lines: Pearson r = %.2f (p = %.3g), with RMSE %.2f to %.2f.',comparison$pearson,comparison$p_value,comparison$baseline_rmse,comparison$program_rmse),'','## Figure sequence','',paste(sprintf('%d. **%s** - %s',seq_along(figs),sub('^[0-9]+_','',sub('\\.pdf$','',figs)),caption[sub('\\.pdf$','',figs)]),collapse='\n'),'', '[Open the figure viewer](index.html), [Frozen genes](frozen_program_genes.tsv), [Code](../../scripts/36_build_source_aware_braf_erbb_programs.R)'),file.path(out,'README.md'))
cards<-paste0('<article><h3>',sub('\\.pdf$','',figs),'</h3><a href="figures/main/',figs,'"><img src="figures/main/',sub('\\.pdf$','.png',figs),'"></a><p>',caption[sub('\\.pdf$','',figs)],'</p><p><a href="figures/main/',figs,'">Open PDF</a></p></article>',collapse='')
writeLines(c('<!doctype html><meta charset="utf-8"><title>Source-aware RAF-ERBB programs</title><style>body{font:16px/1.45 system-ui;background:#f4f6f8;color:#183044;margin:auto;max-width:1100px;padding:24px}article{background:white;border:1px solid #d8e0e6;margin:18px 0;padding:14px}img{display:block;max-width:100%;margin:auto}a{color:#075b8e}</style><h1>Aim 1, source-aware RAF-ERBB programs</h1><p><a href="README.md">Result summary</a>, <a href="frozen_program_genes.tsv">Frozen genes</a>, <a href="../../scripts/36_build_source_aware_braf_erbb_programs.R">Analysis code</a></p>',cards),file.path(out,'index.html'))
message('Stage 36 complete: ',out)
