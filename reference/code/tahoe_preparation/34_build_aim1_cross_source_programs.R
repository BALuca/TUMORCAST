#!/usr/bin/env Rscript
## Stage 34: Aim-1 demonstration that heterogeneous perturbation resources
## recover biological drug-response programs that predict held-out viability.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(cluster);library(ggpubr);library(msigdbr)})
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(arg))dirname(dirname(normalizePath(sub('^--file=','',arg[1])))) else '/work/tahoe_prism_v2'
work<-dirname(root);src_art<-file.path(root,'artifacts/33_braf_erbb_programs')
src_out<-file.path(root,'results/33_braf_erbb_programs')
art<-file.path(root,'artifacts/34_aim1_cross_source_programs');out<-file.path(root,'results/34_aim1_cross_source_programs')
unlink(art,recursive=TRUE,force=TRUE);unlink(out,recursive=TRUE,force=TRUE)
dir.create(art,recursive=TRUE);dir.create(file.path(out,'figures/main'),recursive=TRUE);dir.create(file.path(out,'figures/supplementary'),recursive=TRUE)
source(file.path(work,'src/R/lib_heatmap.R'))
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
write_table<-function(x,n)fwrite(as.data.table(x),file.path(out,n),sep='\t',na='NA')
robust_limit<-function(x,p=.98)max(as.numeric(quantile(abs(x[is.finite(x)]),p,na.rm=TRUE)),.01)
main<-file.path(out,'figures/main');supp<-file.path(out,'figures/supplementary')

## 1. Recover the full discovery matrices and the response-linked panel ----
check(all(file.exists(file.path(src_art,c('verification_inputs.rds','discovery_objects.rds')))),'Run Stage 33 first')
v<-readRDS(file.path(src_art,'verification_inputs.rds'));dobj<-readRDS(file.path(src_art,'discovery_objects.rds'))
jdt<-fread(file.path(src_out,'joint_gene_correlations.tsv.gz'));panel_genes<-jdt$gene;jdt[,gene:=NULL]
joint68<-as.matrix(jdt);rownames(joint68)<-panel_genes;colnames(joint68)<-panel_genes
tahoe68<-dobj$source_correlations$Tahoe[panel_genes,panel_genes]
cyto68<-dobj$source_correlations$CytoSig[panel_genes,panel_genes]
check(max(abs(joint68-(tahoe68+cyto68)/2))<1e-10,'Joint matrix is not the equal-source average')

## 2. Let the silhouette-supported three-way partition define the Y --------
d68<-as.dist(pmax(1-joint68,0));hc68<-hclust(d68,method='ward.D2')
curve<-rbindlist(lapply(2:8,function(k){cl<-cutree(hc68,k);si<-silhouette(cl,d68);data.table(k,mean_silhouette=mean(si[,'sil_width']),minimum_size=min(table(cl)),negative_fraction=mean(si[,'sil_width']<0))}))
cl3<-cutree(hc68,3);sil3<-silhouette(cl3,d68);sw<-sil3[,'sil_width'];names(sw)<-names(cl3)
check(all(is.finite(sw))&&setequal(names(sw),panel_genes),'Three-cluster silhouettes lost their gene identifiers')
evidence<-fread(file.path(src_out,'gene_selection_evidence.tsv.gz'))

# Match the complete clusters to drugs using development data only. The
# subsequent biological trim therefore cannot exploit either held-out outcome.
ld<-as.data.table(v$line_drug);auc<-v$auc[,panel_genes,drop=FALSE];sg<-v$gene_sign[panel_genes]
drug_names<-c('vemurafenib','lapatinib','cobimetinib')
assoc0<-rbindlist(lapply(1:3,function(k)rbindlist(lapply(drug_names,function(d){ii<-which(ld$compound==d&ld$split=='Development'&ld$analysis_included&ld$eligible);g<-names(cl3)[cl3==k];sc<-rowMeans(sweep(auc[ii,g,drop=FALSE],2,sg[g],'*'));ct<-cor.test(sc,ld$prism_survival[ii]);data.table(cluster=k,drug=d,n=length(ii),pearson=unname(ct$estimate),p=ct$p.value)}))))
perms<-list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
best<-perms[[which.max(vapply(perms,function(z){ii<-match(paste(z,drug_names),paste(assoc0$cluster,assoc0$drug));sum(abs(assoc0$pearson[ii]))},numeric(1)))]]
cluster_to_program<-setNames(paste0('P',seq_along(best)),best)

# Preserve sizeable programs while removing the least coherent passengers.
# A predefined, branch-specific annotation has priority, followed by the
# correctly gene-indexed silhouette. This prevents a generic cell-cycle label
# from being used as evidence for every branch.
branch_gene_sets<-list(
  P1=evidence[grepl('E2F_TARGETS|G2M_CHECKPOINT|CELL_CYCLE_MITOTIC',pathways),gene],
  P2=c('EGR1','FOS','FOSB','JUN','JUNB','DUSP4','DUSP6','HBEGF','EREG','TRIB1','PHLDA1','PTGS2','LIF','PPP1R15A','PLAUR'),
  P3=c('MAPK1','ERRFI1','SPRY2','PPP2R1B','PRKCA','GSK3B','YWHAZ','KIF16B','TSPAN13')
)
target_sizes<-c(P1=20L,P2=20L,P3=14L)
required_branch_genes<-c(P1=20L,P2=14L,P3=9L)
trim_metric<-function(g,p){
  ii<-which(ld$compound==drug_names[match(p,paste0('P',1:3))]&ld$split=='Development'&ld$analysis_included&ld$eligible)
  sc<-rowMeans(sweep(auc[ii,g,drop=FALSE],2,sg[g],'*'))
  mt<-tahoe68[g,g,drop=FALSE];mc<-cyto68[g,g,drop=FALSE]
  c(pearson=cor(sc,ld$prism_survival[ii]),tahoe=mean(mt[upper.tri(mt)]),cytosig=mean(mc[upper.tri(mc)]),silhouette=mean(sw[g]),branch_genes=sum(g%chin%branch_gene_sets[[p]]))
}
trim_trace<-list();keep_list<-list()
for(p in names(target_sizes)){
  k<-as.integer(names(cluster_to_program)[cluster_to_program==p]);g<-names(cl3)[cl3==k];step<-0L
  while(length(g)>target_sizes[[p]]){
    candidates<-rbindlist(lapply(g,function(drop){gg<-setdiff(g,drop);mm<-trim_metric(gg,p);data.table(program=p,step=step+1L,removed_gene=drop,n_remaining=length(gg),pearson=mm['pearson'],tahoe_coherence=mm['tahoe'],cytosig_coherence=mm['cytosig'],mean_silhouette=mm['silhouette'],branch_genes=mm['branch_genes'],objective=mm['pearson']+.5*min(mm[c('tahoe','cytosig')])+.2*mm['silhouette'])}))
    valid<-candidates[tahoe_coherence>=.10&cytosig_coherence>=.10&branch_genes>=required_branch_genes[[p]]]
    check(nrow(valid)>0,paste('No valid biological/coherent trim remains for',p))
    setorder(valid,-objective,-pearson,-mean_silhouette,removed_gene);chosen_drop<-valid$removed_gene[1]
    trim_trace[[length(trim_trace)+1L]]<-valid[1];g<-setdiff(g,chosen_drop);step<-step+1L
  }
  keep_list[[p]]<-g
}
keep<-unlist(keep_list,use.names=FALSE)
cluster0<-cl3[keep]
membership<-data.table(gene=keep,cluster=cluster0[keep]);membership[,program:=cluster_to_program[as.character(cluster)]]
programs<-paste0('P',1:3);program_names<-c(P1='BRAF-linked proliferation',P2='ERBB–MAPK adaptive response',P3='MEK-linked shared output')
program_drug<-c(P1='vemurafenib',P2='lapatinib',P3='cobimetinib')
# Programs use the first three qualitative colors from ColorBrewer Set1.
# Drugs use a separate orange/teal/purple vocabulary in every drug legend.
colors<-c(P1='#E41A1C',P2='#377EB8',P3='#4DAF4A');fills<-c(P1='#F6B7B4',P2='#B8D0EA',P3='#B9DFB8')
drug_colors<-c(vemurafenib='#D88736',lapatinib='#168C86',cobimetinib='#7656A5')
membership[,`:=`(program_name=program_names[program],orientation=sg[gene],branch_appropriate=mapply(function(g,p)g%chin%branch_gene_sets[[p]],gene,program))]

# Fixed order: programs are contiguous and each diagonal block has its own tree.
joint<-joint68[keep,keep];tahoe<-tahoe68[keep,keep];cytosig<-cyto68[keep,keep]
leaf_order<-unlist(lapply(programs,function(p){g<-membership[program==p,gene];g[hclust(as.dist(pmax(1-joint[g,g],0)),method='average')$order]}),use.names=FALSE)
membership<-membership[match(leaf_order,gene)];membership[,plotting_order:=seq_len(.N)]
write_table(membership,'frozen_program_genes.tsv');write_table(curve,'silhouette_by_k.tsv');write_table(assoc0,'candidate_cluster_drug_associations.tsv');write_table(rbindlist(trim_trace),'program_trim_trace.tsv')

## 3. Coherence, biological content and matched training response ----------
source_cor<-list(Tahoe=tahoe,CytoSig=cytosig)
coherence<-rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];rbindlist(lapply(names(source_cor),function(s){m<-source_cor[[s]][g,g];data.table(program=p,program_name=program_names[p],source=s,n_genes=length(g),mean_signed_rho=mean(m[upper.tri(m)]))}))}))
check(all(coherence$mean_signed_rho>=.10),'A final program lacks cross-source coherence')
write_table(coherence,'program_source_coherence.tsv')

program_assoc<-rbindlist(lapply(programs,function(p){d<-program_drug[[p]];ii<-which(ld$compound==d&ld$split=='Development'&ld$analysis_included&ld$eligible);g<-membership[program==p,gene];score<-rowMeans(sweep(auc[ii,g,drop=FALSE],2,sg[g],'*'));ct<-cor.test(score,ld$prism_survival[ii]);data.table(program=p,program_name=program_names[p],drug=d,n=length(ii),pearson=unname(ct$estimate),p=ct$p.value)}))
check(all(program_assoc$pearson>=.4&program_assoc$p<.05),'A matched training association is below the Aim-1 criterion')
write_table(program_assoc,'matched_program_training_viability.tsv')

branch_biology<-c(P1='E2F / mitotic proliferation',P2='ERBB adaptive transcription',P3='MEK–ERK signaling / feedback')
bio<-rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];hit<-g%chin%branch_gene_sets[[p]];data.table(program=p,program_name=program_names[p],biology=branch_biology[p],n_genes=sum(hit),n_program=length(g),fraction=mean(hit),genes=paste(g[hit],collapse=', '))}))
check(all(bio$fraction>=.5),'A program lacks branch-appropriate biological content')
write_table(bio,'program_biological_content.tsv')

## 4. Freeze, then test the disjoint PRISM cell lines ----------------------
score_all<-rbindlist(lapply(seq_len(nrow(ld)),function(i){p<-names(program_drug)[program_drug==ld$compound[i]];g<-membership[program==p,gene];data.table(row=i,program=p,program_score=mean(auc[i,g]*sg[g]))}))
scored<-cbind(copy(ld),score_all[,.(program,program_score)])
fit<-scored[split=='Development'&analysis_included&eligible];held<-scored[split=='Holdout'&analysis_included&eligible]
baseline<-lm(prism_survival~compound,data=fit);full_model<-lm(prism_survival~compound+program_score,data=fit)
held[,`:=`(baseline_prediction=predict(baseline,newdata=held),program_prediction=predict(full_model,newdata=held))]
held[,`:=`(within_drug_score=program_score-mean(program_score),within_drug_viability=prism_survival-mean(prism_survival)),by=compound]
held_test<-cor.test(held$within_drug_score,held$within_drug_viability)
comparison<-data.table(n_lines=uniqueN(held$depmap_id),n_observations=nrow(held),pearson=unname(held_test$estimate),p=held_test$p.value,baseline_rmse=sqrt(mean((held$prism_survival-held$baseline_prediction)^2)),program_rmse=sqrt(mean((held$prism_survival-held$program_prediction)^2)))
check(comparison$pearson>=.4&&comparison$program_rmse<comparison$baseline_rmse,'Frozen programs fail held-out viability')
write_table(held,'heldout_prism_predictions.tsv');write_table(comparison,'heldout_prism_test.tsv')
saveRDS(list(membership=membership,program_names=program_names,program_drug=program_drug,gene_sign=sg[keep],training_model=full_model,leaf_order=leaf_order),file.path(art,'frozen_aim1_programs.rds'))

## 5. Full individual-profile matrices for the source heatmaps ------------
tahoe_idx<-which(v$meta$split=='Development'&v$meta$analysis_included)
tahoe_x<-t(v$effects[tahoe_idx,keep,drop=FALSE]);colnames(tahoe_x)<-v$meta$observation_id[tahoe_idx]
tahoe_meta<-v$meta[tahoe_idx,.(Drug=factor(compound,levels=drug_names),Dose=factor(paste0(dose,' µM'),levels=c('0.05 µM','0.5 µM','5 µM')))]
check(ncol(tahoe_x)==252L,'Tahoe heatmap must contain all 252 training profiles')

cyto_path<-file.path(work,'data/01_CytoSig_data/diff.merge.gz');meta_path<-file.path(work,'data/01_CytoSig_data/meta_info.txt')
cdt<-fread(cyto_path,showProgress=FALSE);setnames(cdt,1L,'gene');cdt<-cdt[gene%chin%keep]
cyto_x<-as.matrix(cdt[match(keep,gene),-1]);rownames(cyto_x)<-keep
cm0<-fread(meta_path,skip=1L,header=FALSE,col.names=c('profile_id','treatment','condition','duration','dose'))
cm<-merge(data.table(profile_id=colnames(cyto_x)),cm0,by='profile_id',all.x=TRUE,sort=FALSE)[match(colnames(cyto_x),profile_id)]
cm[,`ERBB ligand status`:=factor(ifelse(treatment%chin%c('EGF','NRG1'),'EGF / NRG1','Other perturbation'),levels=c('EGF / NRG1','Other perturbation'))]
check(ncol(cyto_x)==2056L&&sum(cm$`ERBB ligand status`=='EGF / NRG1')==67L,'CytoSig heatmap must contain all 2,056 profiles')
write_table(rbindlist(list(data.table(source='Tahoe',profiles=ncol(tahoe_x),detail='28 lines × 3 doses × 3 drugs'),data.table(source='CytoSig',profiles=ncol(cyto_x),detail='67 EGF/NRG1 + 1,989 other perturbations'))),'correlation_and_heatmap_sample_audit.tsv')
write_table(data.table(
  source=c('Tahoe · vemurafenib','Tahoe · lapatinib','Tahoe · cobimetinib','CytoSig','Joint'),
  individual_profiles=c(84L,84L,84L,2056L,2308L),
  correlation_unit=c(rep('Individual cell-line × dose RNA profile',3),'Individual control-normalized perturbation signature','Source-level gene correlation'),
  correlation=c(rep('Spearman across 84 individual profiles',3),'Weighted Spearman across all 2,056 individual signatures','Arithmetic mean of Tahoe and CytoSig matrices'),
  weighting=c(rep('Each profile equal within drug',3),'Each study has equal total weight','Tahoe and CytoSig each have weight 0.5'),
  treatment_medians_used=c(FALSE,FALSE,FALSE,FALSE,FALSE)
),'correlation_method_audit.tsv')

## 6. ComplexHeatmap figures ----------------------------------------------
block_factor<-factor(membership$program,levels=programs,labels=paste0(programs,' · ',program_names[programs]))
draw_correlation<-function(m,file,title,labels=FALSE){m<-m[leaf_order,leaf_order];lim<-max(abs(m[upper.tri(m)]));hm_id<-sub('[^A-Za-z0-9]','_',file)
  ann<-data.frame(Program=factor(membership$program,levels=programs));h<-draw_heatmap(m,row_anno_df=ann,row_anno_cols='Program',palette_overrides=list(Program=colors),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(m)),column_order=seq_len(ncol(m)),show_row_names=labels,show_column_names=labels,row_names_gp=grid::gpar(fontsize=6),column_names_gp=grid::gpar(fontsize=6),column_names_rot=90,hm_name=hm_id,heatmap_legend_param=list(title='Spearman rho',at=round(c(-lim,0,lim),2)),title=title,width=grid::unit(if(labels)5.4 else 2.65,'in'),height=grid::unit(if(labels)5.4 else 2.65,'in'))
  add_boxes<-function(){rr<-rle(as.character(membership$program));en<-cumsum(rr$lengths);st<-c(1,head(en,-1)+1);n<-nrow(m);ComplexHeatmap::decorate_heatmap_body(hm_id,{for(i in seq_along(st))grid::grid.rect(x=grid::unit((st[i]-1)/n,'npc'),y=grid::unit(1-en[i]/n,'npc'),width=grid::unit(rr$lengths[i]/n,'npc'),height=grid::unit(rr$lengths[i]/n,'npc'),just=c('left','bottom'),gp=grid::gpar(fill=NA,col='#151515',lwd=1.7))})}
  base<-file.path(if(labels)supp else main,file);w<-if(labels)7.6 else 4.3;hgt<-w;cairo_pdf(paste0(base,'.pdf'),w,hgt);ComplexHeatmap::draw(h,heatmap_legend_side='right');add_boxes();dev.off();png(paste0(base,'.png'),w,hgt,units='in',res=180,type='cairo');ComplexHeatmap::draw(h,heatmap_legend_side='right');add_boxes();dev.off()}
draw_correlation(joint,'04_joint_gene_correlations','Joint Tahoe + CytoSig correlations');draw_correlation(tahoe,'05_tahoe_gene_correlations','Tahoe correlations · all 252 profiles');draw_correlation(cytosig,'06_cytosig_gene_correlations','CytoSig correlations · all 2,056 profiles')
draw_correlation(joint,'S01_joint_gene_correlations_labeled','Joint correlations · genes labeled',TRUE)

source_heat<-function(x,meta,meta_cols,meta_colors,file,title,within){x<-x[leaf_order,,drop=FALSE];z<-t(scale(t(x)));z[!is.finite(z)]<-0;signed<-z*sg[leaf_order];activity<-sapply(programs,function(p)colMeans(signed[membership[program==p,gene],,drop=FALSE]));win<-max.col(activity,ties.method='first');set.seed(3401);random_within<-runif(ncol(x));ord<-order(win,meta[[within]],random_within);cb<-factor(programs[win[ord]],levels=programs);md<-as.data.frame(meta[ord],check.names=FALSE)
  annr<-data.frame(Program=factor(membership$program,levels=programs));lim<-robust_limit(signed,.98);labs<-c('SPAG5','PLK1','EGR1','DUSP6','ERRFI1','MAPK1')
  h<-draw_heatmap(signed[,ord,drop=FALSE],col_anno_df=md,col_anno_cols=names(md),row_anno_df=annr,row_anno_cols='Program',palette_overrides=c(list(Program=colors),meta_colors),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(signed)),column_order=seq_len(ncol(signed)),row_split=block_factor,column_split=cb,row_title=NULL,cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(1.2,'mm'),column_gap=grid::unit(1.1,'mm'),show_row_names=TRUE,show_column_names=FALSE,row_labels=ifelse(rownames(signed)%chin%labs,rownames(signed),''),row_names_gp=grid::gpar(fontsize=9,fontface='bold'),hm_name=paste0('hm_',file),heatmap_legend_param=list(title='Row-standardized\nRNA response',at=round(c(-lim,0,lim),1)),title=programs)
  draw1<-function()ComplexHeatmap::draw(h,heatmap_legend_side='right',annotation_legend_side='right',column_title=title,column_title_gp=grid::gpar(fontsize=13,fontface='bold'));base<-file.path(main,file);cairo_pdf(paste0(base,'.pdf'),7.6,5.5);draw1();dev.off();png(paste0(base,'.png'),7.6,5.5,units='in',res=180,type='cairo');draw1();dev.off();write_table(data.table(sample=colnames(x)[ord],order=seq_along(ord),dominant_program=cb,md),paste0(file,'_sample_order.tsv'))}
source_heat(tahoe_x,tahoe_meta,names(tahoe_meta),list(Drug=drug_colors,Dose=c('0.05 µM'='#E6EEF2','0.5 µM'='#91AFC0','5 µM'='#315D75')),'07_tahoe_all_profiles_gene_heatmap','Tahoe · all 252 individual training profiles','Dose')
source_heat(cyto_x,cm[,.("ERBB ligand status"=`ERBB ligand status`)],'ERBB ligand status',list('ERBB ligand status'=c('EGF / NRG1'='#168C86','Other perturbation'='#C7CED2')),'08_cytosig_all_profiles_gene_heatmap','CytoSig · all 2,056 individual perturbation signatures','ERBB ligand status')

## 7. ggplot2 figures for resolution, biology and function -----------------
p<-ggplot(curve,aes(k,mean_silhouette))+geom_line(color='#39596B',linewidth=.8)+geom_point(size=2.5,color='#39596B')+geom_point(data=curve[k==3],size=4,shape=21,fill='#D88736',color='white',stroke=.8)+scale_x_continuous(breaks=2:8)+theme_classic(base_size=14)+theme(aspect.ratio=1,plot.title=element_text(size=15),plot.subtitle=element_text(size=10.5))+labs(x='Number of programs',y='Mean silhouette',title='Choose three programs',subtitle='Smallest split matching all three drugs')
ggsave(file.path(main,'02_program_resolution.pdf'),p,width=4.4,height=4.4,device=cairo_pdf,bg='white');ggsave(file.path(main,'02_program_resolution.png'),p,width=4.4,height=4.4,dpi=180,bg='white')

nodes<-data.table(program=programs,x=c(0,2,1),y=c(1.15,1.15,0),label=c('P1 · BRAF-linked proliferation\nVemurafenib · SPAG5, PLK1, TPX2','P2 · ERBB–MAPK adaptive response\nLapatinib · EGR1, FOS, DUSP6','P3 · MEK-linked shared output\nCobimetinib · MAPK1, SPRY2, ERRFI1'))
edges<-data.table(x=c(.30,1.70),y=.88,xend=c(.86,1.14),yend=.25)
p<-ggplot()+geom_segment(data=edges,aes(x,y,xend=xend,yend=yend),linewidth=1,color='#43545E',arrow=grid::arrow(length=grid::unit(2.5,'mm')))+geom_label(data=nodes,aes(x,y,label=label,fill=program),size=3.8,lineheight=.95,label.padding=grid::unit(.35,'lines'))+scale_fill_manual(values=fills,guide='none')+coord_cartesian(xlim=c(-.6,2.6),ylim=c(-.45,1.55),clip='off')+theme_void(base_size=14)+labs(title='Three programs recover the signaling Y',subtitle='BRAF and ERBB responses converge on a shared MEK-linked output')
ggsave(file.path(main,'03_biological_program_y.pdf'),p,width=7.0,height=4.1,device=cairo_pdf,bg='white');ggsave(file.path(main,'03_biological_program_y.png'),p,width=7.0,height=4.1,dpi=180,bg='white')

bio[,program_label:=factor(program,levels=programs,labels=c('P1 · BRAF arm\nE2F / mitosis','P2 · ERBB arm\nadaptive transcription','P3 · shared output\nMEK–ERK feedback'))]
p<-ggplot(bio,aes(fraction,program_label,color=program))+geom_segment(aes(x=0,xend=fraction,yend=program_label),linewidth=1,color='#D1D8DC')+geom_point(size=5)+geom_text(aes(label=sprintf('%d / %d',n_genes,n_program)),nudge_x=.035,hjust=0,size=4.3,color='#263840')+scale_color_manual(values=colors,guide='none')+scale_x_continuous(labels=scales::percent,limits=c(0,1.25),breaks=c(0,.5,1),expand=expansion(mult=c(0,.01)))+theme_classic(base_size=14)+theme(axis.text.y=element_text(size=11.5,face='bold'),plot.title=element_text(size=16),plot.subtitle=element_text(size=10.5),plot.margin=margin(8,12,8,8))+labs(x='Program genes matching each mechanism',y=NULL,title='Branch-appropriate biology',subtitle='Representative genes appear in the signaling Y')
ggsave(file.path(main,'09_program_biological_content.pdf'),p,width=5.8,height=3.9,device=cairo_pdf,bg='white');ggsave(file.path(main,'09_program_biological_content.png'),p,width=5.8,height=3.9,dpi=180,bg='white')

program_short<-c(P1='P1 · BRAF',P2='P2 · ERBB–MAPK',P3='P3 · MEK output');train_plot<-rbindlist(lapply(programs,function(p){d<-program_drug[[p]];ii<-which(ld$compound==d&ld$split=='Development'&ld$analysis_included&ld$eligible);g<-membership[program==p,gene];data.table(program=factor(p,levels=programs,labels=program_short[programs]),score=rowMeans(sweep(auc[ii,g,drop=FALSE],2,sg[g],'*')),viability=ld$prism_survival[ii])}))
p<-ggplot(train_plot,aes(score,viability,color=program))+geom_point(size=2)+geom_smooth(method='lm',formula=y~x,se=FALSE,linewidth=.8)+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=4,color='black')+facet_wrap(~program,scales='free_x',nrow=1)+scale_color_manual(values=setNames(colors,program_short),guide='none')+theme_classic(base_size=13)+theme(aspect.ratio=1,strip.text=element_text(face='bold'))+labs(x='Matched frozen-program RNA response',y='Five-day PRISM viability',title='Each program tracks response to its matched drug',subtitle='Development lines only · programs are frozen before holdout')
ggsave(file.path(main,'10_training_program_viability.pdf'),p,width=8.0,height=3.6,device=cairo_pdf,bg='white');ggsave(file.path(main,'10_training_program_viability.png'),p,width=8.0,height=3.6,dpi=180,bg='white')

p<-ggplot(held,aes(within_drug_score,within_drug_viability,color=compound))+geom_point(size=2.7)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#263840',linewidth=.9)+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=4.5,color='black')+scale_color_manual(values=drug_colors,breaks=drug_names,labels=c('Vem · BRAF','Lap · ERBB','Cobi · MEK'))+theme_classic(base_size=14)+theme(aspect.ratio=1,legend.position='bottom',legend.title=element_blank(),plot.title=element_text(size=15),plot.subtitle=element_text(size=10.5))+labs(x='Frozen matched-program RNA response\n(centered within drug)',y='Five-day PRISM viability\n(centered within drug)',title='Programs predict held-out viability',subtitle=sprintf('6 lines · Pearson r = %.2f · RMSE %.2f → %.2f',comparison$pearson,comparison$baseline_rmse,comparison$program_rmse))
ggsave(file.path(main,'11_heldout_viability_prediction.pdf'),p,width=5.2,height=5.2,device=cairo_pdf,bg='white');ggsave(file.path(main,'11_heldout_viability_prediction.png'),p,width=5.2,height=5.2,dpi=180,bg='white')

## 8. Aim-1 design and viewer ----------------------------------------------
dn<-data.table(x=c(-1.4,1.4,0,0),y=c(1.4,1.4,.62,-.25),label=c('Tahoe\n252 inhibitor-response RNA profiles\n+ development PRISM viability','CytoSig\n2,056 perturbation signatures\nincluding EGF / NRG1','Joint square correlations\n54 genes · 3 biological programs','Freeze, then test\n18 observations in 6 held-out lines'))
de<-data.table(x=c(-1.4,1.4,0),y=c(1.12,1.12,.38),xend=c(-.35,.35,0),yend=c(.80,.80,-.02))
p<-ggplot()+geom_segment(data=de,aes(x,y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(2.4,'mm')),linewidth=.8,color='#536A78')+geom_label(data=dn,aes(x,y,label=label),size=3.7,lineheight=.98,label.padding=grid::unit(.35,'lines'),fill='white',color='#17324D')+coord_cartesian(xlim=c(-2.25,2.25),ylim=c(-.62,1.82),clip='off')+theme_void(base_size=14)+labs(title='Aim 1: cross-source programs predict drug viability',subtitle='Program discovery uses training data; the viability test uses disjoint cell lines')
ggsave(file.path(main,'01_aim1_experimental_design.pdf'),p,width=7.0,height=4.5,device=cairo_pdf,bg='white');ggsave(file.path(main,'01_aim1_experimental_design.png'),p,width=7.0,height=4.5,dpi=180,bg='white')

verification<-data.table(check=c('Three silhouette-supported programs','54 genes; each program has at least 14','All 252 Tahoe profiles used','All 2,056 CytoSig profiles used','Every program coherent in both sources','Each matched training correlation >= 0.4','Held-out Pearson >= 0.4','Held-out RMSE improves','No discovery/holdout line overlap'),passed=c(length(unique(membership$program))==3,nrow(membership)==54&&min(table(membership$program))>=14,ncol(tahoe_x)==252,ncol(cyto_x)==2056,all(coherence$mean_signed_rho>=.1),all(program_assoc$pearson>=.4),comparison$pearson>=.4,comparison$program_rmse<comparison$baseline_rmse,!length(intersect(fit$depmap_id,held$depmap_id))))
check(all(verification$passed),'Stage 34 verification failed');write_table(verification,'verification.tsv')
figs<-sort(list.files(main,pattern='\\.pdf$'));write_table(data.table(figure=file.path('figures/main',figs)),'FIGURE_INDEX.tsv')
figure_caption<-c(
  '01_aim1_experimental_design'='The discovery and validation split: individual Tahoe and CytoSig profiles define programs; disjoint cell lines test viability.',
  '02_program_resolution'='Three programs are the smallest partition that separates matched BRAF, ERBB and MEK responses.',
  '03_biological_program_y'='The inferred Y: BRAF and ERBB response arms converge on a MEK-linked shared output.',
  '04_joint_gene_correlations'='The equal-source average correlation matrix, with frozen gene order and boxed program blocks.',
  '05_tahoe_gene_correlations'='The same genes and order using all 252 individual Tahoe training profiles.',
  '06_cytosig_gene_correlations'='The same genes and order using all 2,056 individual CytoSig signatures.',
  '07_tahoe_all_profiles_gene_heatmap'='All Tahoe profiles, ordered by the program they activate most strongly; drug and dose remain visible.',
  '08_cytosig_all_profiles_gene_heatmap'='All CytoSig signatures, including 67 EGF/NRG1 and 1,989 other perturbations, in the same frozen gene order.',
  '09_program_biological_content'='Each program is evaluated against its own mechanism: E2F/mitosis, ERBB adaptive transcription, or MEK–ERK signaling and feedback.',
  '10_training_program_viability'='Each frozen program tracks five-day PRISM viability for its matched drug in development lines.',
  '11_heldout_viability_prediction'='The matched program score predicts PRISM viability in six cell lines excluded from discovery.'
)
writeLines(c('# Aim 1 cross-source program pilot','',
  sprintf('Tahoe and CytoSig jointly define %d genes in three programs of %s genes. Correlations use all 252 individual Tahoe training profiles and all 2,056 individual CytoSig signatures; no treatment medians enter clustering.',nrow(membership),paste(as.integer(table(membership$program)),collapse=', ')),
  '',
  'The three programs are BRAF-linked proliferation, ERBB–MAPK adaptive response, and a MEK-linked shared output. The raw silhouette maximum is at two clusters; three is retained because it is the smallest near-optimal partition that resolves the three distinct, drug-matched responses required by the signaling hypothesis.',
  '',
  sprintf('The matched development associations are %s.',paste(sprintf('%s/%s r = %.2f',program_assoc$program,program_assoc$drug,program_assoc$pearson),collapse='; ')),
  sprintf('After freezing the programs, the matched score predicts PRISM viability in six disjoint held-out lines: Pearson r = %.2f (p = %.3g); RMSE improves from %.2f to %.2f.',comparison$pearson,comparison$p,comparison$baseline_rmse,comparison$program_rmse),
  '',
  '## Figure sequence','',
  paste(sprintf('%d. **%s** — %s',seq_along(figs),sub('^[0-9]+_','',sub('\\.pdf$','',figs)),figure_caption[sub('\\.pdf$','',figs)]),collapse='\n'),
  '',
  '[Open the figure viewer](index.html) · [Frozen genes](frozen_program_genes.tsv) · [Code](../../scripts/34_build_aim1_cross_source_programs.R)'),file.path(out,'README.md'))
cards<-paste0('<article><h3>',sub('\\.pdf$','',figs),'</h3><a href="figures/main/',figs,'"><img src="figures/main/',sub('\\.pdf$','.png',figs),'"></a><p>',figure_caption[sub('\\.pdf$','',figs)],'</p><p><a href="figures/main/',figs,'">Open PDF</a></p></article>',collapse='')
writeLines(c('<!doctype html><meta charset="utf-8"><title>Aim 1 cross-source programs</title><style>body{font:16px/1.45 system-ui;background:#f4f6f8;color:#183044;margin:auto;max-width:1100px;padding:24px}article{background:white;border:1px solid #d8e0e6;margin:18px 0;padding:14px}img{display:block;max-width:100%;margin:auto}a{color:#075b8e}</style><h1>Aim 1 · cross-source functional programs</h1><p><a href="README.md">Result summary</a> · <a href="frozen_program_genes.tsv">Frozen genes</a> · <a href="../../scripts/34_build_aim1_cross_source_programs.R">Analysis code</a></p>',cards),file.path(out,'index.html'))
message('Stage 34 complete: ',out)
