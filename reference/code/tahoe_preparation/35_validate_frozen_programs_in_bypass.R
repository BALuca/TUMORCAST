#!/usr/bin/env Rscript
## Stage 35: transfer the exact Stage-34 programs into the bypass mechanism.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggpubr)})
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
root<-if(length(arg))dirname(dirname(normalizePath(sub('^--file=','',arg[1])))) else '/work/tahoe_prism_v2'
work<-dirname(root);src34<-file.path(root,'artifacts/34_aim1_cross_source_programs');src33<-file.path(root,'artifacts/33_braf_erbb_programs')
out<-file.path(root,'results/35_frozen_program_bypass_validation');art<-file.path(root,'artifacts/35_frozen_program_bypass_validation')
unlink(out,recursive=TRUE,force=TRUE);unlink(art,recursive=TRUE,force=TRUE);dir.create(file.path(out,'figures/main'),recursive=TRUE);dir.create(file.path(out,'figures/supplementary'),recursive=TRUE);dir.create(art,recursive=TRUE)
source(file.path(work,'src/R/lib_heatmap.R'))
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
write_table<-function(x,n)fwrite(as.data.table(x),file.path(out,n),sep='\t',na='NA')
main<-file.path(out,'figures/main');supp<-file.path(out,'figures/supplementary')
program_colors<-c(P1='#E41A1C',P2='#377EB8',P3='#4DAF4A')
drug_colors<-c('Vemurafenib'='#D88736','Vemurafenib + cobimetinib'='#7656A5')

## 1. Load the frozen Aim-1 definition without changing a gene or sign ------
frozen<-readRDS(file.path(src34,'frozen_aim1_programs.rds'));v<-readRDS(file.path(src33,'verification_inputs.rds'))
membership<-as.data.table(frozen$membership);setorder(membership,plotting_order)
programs<-c('P1','P2','P3');genes<-membership$gene;gene_sign<-frozen$gene_sign[genes]
check(identical(genes,as.character(frozen$leaf_order)),'Aim-1 gene order changed during transfer')
check(identical(as.character(membership$program),rep(programs,c(20L,20L,14L))),'Aim-1 membership changed during transfer')
write_table(membership,'frozen_program_genes_exact_copy.tsv')

score_programs<-function(x){
  sapply(programs,function(p){g<-membership[program==p,gene];present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)];check(length(present)>0,paste('No measured genes for',p));colMeans(sweep(x[present,,drop=FALSE],1,gene_sign[present],'*'))})
}
coverage_table<-function(x,dataset)rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];ok<-rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x);data.table(dataset,program=p,n_frozen=length(g),n_measured=sum(ok),missing_genes=paste(g[!ok],collapse=', '))}))

## 2. Explicit EGF bypass in GSE127988 --------------------------------------
x127<-v$d127[genes,,drop=FALSE];m127<-as.data.table(v$s127);m127<-m127[match(colnames(x127),gsm)]
coverage<-coverage_table(x127,'GSE127988');check(coverage[program=='P1',n_measured]==20L&&coverage[program=='P2',n_measured]==19L&&coverage[program=='P3',n_measured]==14L,'Unexpected GSE127988 frozen-gene coverage')
s127<-score_programs(x127)
long127<-melt(data.table(m127[,.(gsm,arm,egf_h,biological_set)],s127),id.vars=c('gsm','arm','egf_h','biological_set'),variable.name='program',value.name='score')
long127[,`:=`(program=factor(program,levels=programs),drug_background=factor(ifelse(arm=='Vem','Vemurafenib','Vemurafenib + cobimetinib'),levels=names(drug_colors)))]
day127<-long127[,.(score=mean(score),se=sd(score)/sqrt(.N)),by=.(drug_background,egf_h,biological_set,program)]
trapz<-function(t,y)sum(diff(t)*(head(y,-1)+tail(y,-1))/2)
auc127<-day127[order(egf_h),.(egf_auc=trapz(egf_h,score)),by=.(drug_background,biological_set,program)]
auc_wide<-dcast(auc127,program+biological_set~drug_background,value.var='egf_auc')
setnames(auc_wide,c('Vemurafenib','Vemurafenib + cobimetinib'),c('vem_auc','vem_cobi_auc'))
auc_wide[,`:=`(cobimetinib_change=vem_cobi_auc-vem_auc,fraction_retained=vem_cobi_auc/vem_auc)]
check(all(auc_wide[as.character(program)%in%c('P2','P3'),vem_auc]>0),'EGF fails to activate a frozen bypass program')
check(all(auc_wide[as.character(program)%in%c('P2','P3'),cobimetinib_change]<0),'Cobimetinib fails to attenuate a frozen bypass program')
check(all(auc_wide[program=='P1',vem_auc]<0),'EGF unexpectedly activates the BRAF-linked proliferation program')
write_table(long127,'gse127988_frozen_program_scores.tsv');write_table(day127,'gse127988_time_course_means.tsv');write_table(auc_wide,'gse127988_egf_auc_and_cobimetinib_effect.tsv')

## 3. Same programs in the matched external functional factorial ------------
x684<-v$d684[genes,,drop=FALSE];m684<-as.data.table(v$s684);m684<-m684[match(colnames(x684),gsm)]
coverage<-rbindlist(list(coverage,coverage_table(x684,'GSE68453')));write_table(coverage,'external_gene_coverage.tsv')
s684<-score_programs(x684)
long684<-melt(data.table(m684[,.(gsm,cell_line,arm,replicate)],s684),id.vars=c('gsm','cell_line','arm','replicate'),variable.name='program',value.name='score')
means684<-long684[,.(score=mean(score)),by=.(cell_line,arm,program)]
factorial684<-dcast(means684,cell_line+program~arm,value.var='score');factorial684[,interaction:=plx_lap-plx-lapatinib+control]
bliss<-data.table(cell_line=c('A375','SKMEL24','HT144','DM331','SKMEL28','DM13','VMM5A','VMM18','DM6','VMM15','VMM12','VMM17'),bliss=c(5.7,29.3,16.3,27.3,49.7,39.3,43,31.7,-2.7,-8.7,3.7,-1))
functional<-merge(factorial684[complete.cases(interaction)],bliss,by='cell_line');functional[,suppression:=-interaction]
functional_stats<-functional[as.character(program)%in%c('P2','P3'),{tt<-cor.test(suppression,bliss);.(n=.N,pearson=unname(tt$estimate),p=tt$p.value)},by=program]
check(all(functional_stats$pearson>=.4),'A frozen bypass program fails the external functional direction')
write_table(factorial684,'gse68453_frozen_program_factorial.tsv');write_table(functional,'gse68453_frozen_program_function.tsv');write_table(functional_stats,'gse68453_functional_prediction_stats.tsv')

## 4. Gerosa phosphoproteomics bridges RNA programs to pathway activity ------
gerosa<-fread(file.path(root,'results/33_braf_erbb_programs/gerosa_panorama_phosphopeptide_summary.tsv'))
gtest<-fread(file.path(root,'results/33_braf_erbb_programs/gerosa_panorama_prediction_tests.tsv'))
check(all(gtest[grepl('^ERK',analyte),EGF_rescue_over_one_uM_vem]>=3),'EGF does not restore ERK phosphorylation as expected')
write_table(gerosa,'gerosa_phosphoproteomic_summary.tsv');write_table(gtest,'gerosa_phosphoproteomic_tests.tsv')

## 5. Grant figures ----------------------------------------------------------
# Design: one frozen definition feeds both Aim 1 and Aim 2.
dn<-data.table(x=c(0,-1.45,1.45,0),y=c(1.55,.62,.62,-.35),program=c(NA,'P2','P3',NA),label=c('Same 54 frozen genes and signs\nfrom Aim 1','EGF → ERBB adaptive arm\nP2','MEK-linked shared output\nP3','External tests\nGSE127988 RNA · Gerosa pERK\nGSE68453 combination viability'))
de<-data.table(x=c(0,0,-1.05,1.05),y=c(1.25,1.25,.30,.30),xend=c(-1.15,1.15,0,0),yend=c(.87,.87,-.05,-.05))
p<-ggplot()+geom_segment(data=de,aes(x,y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(2.4,'mm')),linewidth=.9,color='#50636F')+geom_label(data=dn,aes(x,y,label=label,fill=program),size=4,lineheight=.98,label.padding=grid::unit(.35,'lines'),color='#183044')+scale_fill_manual(values=c(program_colors,`NA`='white'),na.value='white',guide='none')+coord_cartesian(xlim=c(-2.25,2.25),ylim=c(-.75,1.95),clip='off')+theme_void(base_size=14)+theme(plot.margin=margin(10,10,10,10),plot.title=element_text(size=17),plot.subtitle=element_text(size=11))+labs(title='Aim 2 uses the exact Aim 1 programs',subtitle='External outcomes do not change membership, sign or weight')
ggsave(file.path(main,'01_same_programs_aim2_design.pdf'),p,width=7,height=4.7,device=cairo_pdf,bg='white');ggsave(file.path(main,'01_same_programs_aim2_design.png'),p,width=7,height=4.7,dpi=180,bg='white')

# Individual GSE127988 profiles in biological order; EREG remains visible as NA.
set.seed(3501);rnd<-runif(nrow(m127));ord<-order(match(m127$arm,c('Vem','Vem + Cobi')),m127$biological_set,m127$egf_h,rnd)
signed127<-sweep(x127,1,gene_sign,'*');lim<-max(as.numeric(quantile(abs(signed127[is.finite(signed127)]),.975)),.01)
row_group<-factor(membership$program,levels=programs,labels=c('P1 · BRAF proliferation','P2 · ERBB adaptive','P3 · MEK output'))
col_group<-factor(paste(ifelse(m127$arm=='Vem','Vem','Vem + Cobi'),paste0('day ',m127$biological_set),sep=' · '),levels=c('Vem · day 3','Vem · day 4','Vem + Cobi · day 3','Vem + Cobi · day 4'))
ann<-data.frame('Drug background'=factor(ifelse(m127$arm=='Vem','Vemurafenib','Vemurafenib + cobimetinib'),levels=names(drug_colors)),'EGF time'=factor(paste0(m127$egf_h,' h'),levels=paste0(c(0,.5,1,2,3,4,8),' h')),check.names=FALSE)
annr<-data.frame(Program=factor(membership$program,levels=programs));example<-c('PLK1','SPAG5','FOS','DUSP6','HBEGF','EREG','ERRFI1','SPRY2','MAPK1')
render127<-function(base,all_labels=FALSE,w=7.5,hgt=6){
  hm<-draw_heatmap(signed127[,ord,drop=FALSE],col_anno_df=ann[ord,,drop=FALSE],col_anno_cols=names(ann),row_anno_df=annr,row_anno_cols='Program',palette_overrides=list(Program=program_colors,'Drug background'=drug_colors,'EGF time'=setNames(hcl.colors(7,'Blues 3'),levels(ann$`EGF time`))),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(signed127)),column_order=seq_len(ncol(signed127)),row_split=row_group,column_split=col_group[ord],cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(1.2,'mm'),column_gap=grid::unit(1.4,'mm'),row_title=NULL,show_row_names=TRUE,show_column_names=FALSE,row_labels=if(all_labels)genes else ifelse(genes%chin%example,genes,''),row_names_gp=grid::gpar(fontsize=if(all_labels)6.2 else 8.5,fontface=ifelse(genes%chin%example,'bold','plain')),hm_name='EGF response',na_col='#B9BEC1',heatmap_legend_param=list(title='Signed log2 RNA change\nfrom 0 h drug background',at=round(c(-lim,0,lim),1)))
  draw1<-function()ComplexHeatmap::draw(hm,heatmap_legend_side='right',annotation_legend_side='right',column_title='GSE127988 · EGF time course',column_title_gp=grid::gpar(fontsize=13,fontface='bold'))
  cairo_pdf(paste0(base,'.pdf'),w,hgt);draw1();dev.off();png(paste0(base,'.png'),w,hgt,units='in',res=180,type='cairo');draw1();dev.off()
}
render127(file.path(main,'02_gse127988_frozen_gene_heatmap'));render127(file.path(supp,'S01_gse127988_all_gene_labels'),TRUE,8.2,7.3)
write_table(data.table(gsm=m127$gsm[ord],plotting_order=seq_along(ord),drug_background=ann$`Drug background`[ord],egf_time=ann$`EGF time`[ord],biological_day=m127$biological_set[ord]),'gse127988_heatmap_column_order.tsv')

# Time courses: same visualization repeated for all three programs as facets.
program_labels<-c(P1='P1 · BRAF proliferation',P2='P2 · ERBB adaptive',P3='P3 · MEK output')
day127[,program_label:=factor(program,levels=programs,labels=program_labels)]
p<-ggplot(day127,aes(egf_h,score,color=drug_background,group=interaction(drug_background,biological_set)))+geom_hline(yintercept=0,color='#B8C1C6')+geom_line(aes(linetype=factor(biological_set)),linewidth=.9)+geom_point(size=2.2)+facet_wrap(~program_label,nrow=1,scales='free_y')+scale_color_manual(values=drug_colors)+scale_linetype_manual(values=c('3'='solid','4'='dashed'),labels=c('Day 3','Day 4'))+guides(color=guide_legend(order=1,nrow=1),linetype=guide_legend(order=2,nrow=1))+theme_classic(base_size=13)+theme(strip.text=element_text(face='bold'),legend.position='bottom',legend.box='vertical',legend.title=element_text(face='bold'),legend.margin=margin(0,0,0,0))+labs(x='Hours after 100 ng/mL EGF',y='Frozen program RNA response\n(change from 0 h drug background)',color='Drug background',linetype='Biological replicate',title='EGF activates the bypass programs',subtitle='Cobimetinib prevents sustained ERBB-adaptive and shared-output activation')
ggsave(file.path(main,'03_gse127988_program_time_courses.pdf'),p,width=8.2,height=4.25,device=cairo_pdf,bg='white');ggsave(file.path(main,'03_gse127988_program_time_courses.png'),p,width=8.2,height=4.25,dpi=180,bg='white')

# AUC places the two independent biological days side by side.
auc127[,program_label:=factor(program,levels=programs,labels=program_labels)]
p<-ggplot(auc127,aes(drug_background,egf_auc,color=program_label,group=interaction(program_label,biological_set)))+geom_hline(yintercept=0,color='#B8C1C6')+geom_line(linewidth=.8,color='#596A74')+geom_point(size=3)+facet_wrap(~program_label,nrow=1,scales='free_y')+scale_color_manual(values=setNames(program_colors,program_labels),guide='none')+scale_x_discrete(labels=c('Vemurafenib'='Vem','Vemurafenib + cobimetinib'='Vem + Cobi'))+theme_classic(base_size=13)+theme(strip.text=element_text(face='bold'),axis.text.x=element_text(size=10.5))+labs(x='Drug background',y='0–8 h EGF-response AUC',title='MEK inhibition blocks the EGF response',subtitle='Each line is one independent biological day')
ggsave(file.path(main,'04_gse127988_cobimetinib_attenuation.pdf'),p,width=7.8,height=3.7,device=cairo_pdf,bg='white');ggsave(file.path(main,'04_gse127988_cobimetinib_attenuation.png'),p,width=7.8,height=3.7,dpi=180,bg='white')

# Protein bridge in the same A375/vemurafenib/EGF setting.
keep_conditions<-c('DMSO','1 µM vem','1 µM vem + EGF');gp<-gerosa[condition%chin%keep_conditions];gp[,condition:=factor(condition,levels=keep_conditions,labels=c('DMSO','Vem','Vem + EGF'))];gp[,analyte:=factor(analyte,levels=c('ERK1 pT202/pY204','ERK2 pT185/pY187','SOS1 ERK-feedback sites'),labels=c('pERK1','pERK2','pSOS1 feedback'))]
p<-ggplot(gp,aes(condition,mean_relative,group=analyte,color=analyte))+geom_hline(yintercept=1,color='#B8C1C6',linetype='dashed')+geom_errorbar(aes(ymin=pmax(mean_relative-se_relative,0),ymax=mean_relative+se_relative),width=.08,na.rm=TRUE)+geom_line(linewidth=.8)+geom_point(size=3)+facet_wrap(~analyte,nrow=1)+scale_color_manual(values=c('#377EB8','#4DAF4A','#984EA3'),guide='none')+theme_classic(base_size=13)+theme(strip.text=element_text(face='bold'),axis.text.x=element_text(angle=25,hjust=1))+labs(x=NULL,y='Phosphopeptide abundance\n(relative to control)',title='EGF restores ERK phosphorylation after BRAF inhibition',subtitle='Gerosa A375 phosphoproteomics · EGF sample n = 1')
ggsave(file.path(main,'05_gerosa_phosphoproteomic_bridge.pdf'),p,width=7.6,height=3.7,device=cairo_pdf,bg='white');ggsave(file.path(main,'05_gerosa_phosphoproteomic_bridge.png'),p,width=7.6,height=3.7,dpi=180,bg='white')

# Functional prediction: an early molecular interaction predicts later benefit.
fp<-functional[as.character(program)%in%c('P2','P3')];fp[,program_label:=factor(program,levels=c('P2','P3'),labels=program_labels[c('P2','P3')])]
p<-ggplot(fp,aes(suppression,bliss,color=program_label))+geom_point(size=2.8)+geom_text(aes(label=cell_line),size=3,vjust=-.7,check_overlap=TRUE,show.legend=FALSE)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#283C46',linewidth=.9)+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=4,color='black')+facet_wrap(~program_label,nrow=1,scales='free_x')+scale_color_manual(values=setNames(program_colors[c('P2','P3')],program_labels[c('P2','P3')]),guide='none')+scale_x_continuous(expand=expansion(mult=c(.14,.25)))+scale_y_continuous(expand=expansion(mult=c(.06,.22)))+coord_cartesian(clip='off')+theme_classic(base_size=13)+theme(aspect.ratio=1,strip.text=element_text(face='bold'),plot.margin=margin(8,16,8,8))+labs(x='Eight-hour suppression by PLX4720 + lapatinib\n(beyond either drug alone)',y='Three-day added combination benefit\n(Bliss score)',title='Frozen bypass programs predict combination benefit',subtitle='Five melanoma lines with matched molecular and functional factorials')
ggsave(file.path(main,'06_gse68453_functional_prediction.pdf'),p,width=7.2,height=4.2,device=cairo_pdf,bg='white');ggsave(file.path(main,'06_gse68453_functional_prediction.png'),p,width=7.2,height=4.2,dpi=180,bg='white')

## 6. Verification and self-contained viewer --------------------------------
verification<-data.table(check=c('Exact Aim-1 genes and order','Exact Aim-1 signs','GSE127988 coverage explicit','P2 activated by EGF on both days','P3 activated by EGF on both days','Cobimetinib attenuates P2 on both days','Cobimetinib attenuates P3 on both days','P1 not EGF activated','Gerosa pERK rescue at least threefold','P2 functional Pearson at least 0.4','P3 functional Pearson at least 0.4'),passed=c(identical(genes,as.character(frozen$leaf_order)),identical(gene_sign,frozen$gene_sign[genes]),all(coverage[dataset=='GSE127988',n_measured]==c(20,19,14)),all(auc_wide[program=='P2',vem_auc]>0),all(auc_wide[program=='P3',vem_auc]>0),all(auc_wide[program=='P2',cobimetinib_change]<0),all(auc_wide[program=='P3',cobimetinib_change]<0),all(auc_wide[program=='P1',vem_auc]<0),all(gtest[grepl('^ERK',analyte),EGF_rescue_over_one_uM_vem]>=3),functional_stats[program=='P2',pearson]>=.4,functional_stats[program=='P3',pearson]>=.4))
check(all(verification$passed),'Aim-2 verification failed');write_table(verification,'verification.tsv')
saveRDS(list(source_frozen_programs=normalizePath(file.path(src34,'frozen_aim1_programs.rds')),membership=as.data.frame(membership),gene_sign=gene_sign,gse127988_auc=auc_wide,functional_stats=functional_stats),file.path(art,'frozen_program_bypass_transfer.rds'))
figs<-sort(list.files(main,pattern='\\.pdf$'));write_table(data.table(figure=file.path('figures/main',figs)),'FIGURE_INDEX.tsv')
caption<-c('01_same_programs_aim2_design'='The exact Aim 1 genes, signs and weights are carried into all Aim 2 tests.','02_gse127988_frozen_gene_heatmap'='All individual GSE127988 EGF time-course profiles, split by drug background and biological day. EREG is the single unavailable gene and is gray.','03_gse127988_program_time_courses'='EGF activates the ERBB-adaptive and MEK-linked programs during vemurafenib; cobimetinib prevents sustained activation.','04_gse127988_cobimetinib_attenuation'='The same effect summarized as 0–8 h AUC for each of two biological days.','05_gerosa_phosphoproteomic_bridge'='Matched A375 protein data show BRAF inhibition lowers pERK and EGF restores it.','06_gse68453_functional_prediction'='Suppression of the frozen bypass programs predicts later PLX4720 + lapatinib benefit in the five complete factorial cell lines.')
writeLines(c('# Aim 2 bypass validation with the frozen Aim 1 programs','', 'The exact 54 Aim 1 genes, memberships, signs and equal-gene weights are reused. No program is rediscovered from GSE127988, Gerosa or GSE68453.','',sprintf('In GSE127988, EGF activates P2 and P3 on both biological days during vemurafenib. Cobimetinib changes mean P2 EGF-response AUC from %.2f to %.2f and P3 from %.2f to %.2f.',mean(auc_wide[program=='P2',vem_auc]),mean(auc_wide[program=='P2',vem_cobi_auc]),mean(auc_wide[program=='P3',vem_auc]),mean(auc_wide[program=='P3',vem_cobi_auc])),sprintf('In the external five-line functional factorial, early suppression predicts three-day combination benefit: P2 Pearson r = %.2f and P3 r = %.2f.',functional_stats[program=='P2',pearson],functional_stats[program=='P3',pearson]),'', '[Open the figure viewer](index.html) · [Exact frozen genes](frozen_program_genes_exact_copy.tsv) · [Code](../../scripts/35_validate_frozen_programs_in_bypass.R)','', 'GSE127988 measures 53/54 frozen genes. EREG is absent from the processed matrix; it is displayed as missing and is never replaced with zero.'),file.path(out,'README.md'))
cards<-paste0('<article><h3>',sub('\\.pdf$','',figs),'</h3><a href="figures/main/',figs,'"><img src="figures/main/',sub('\\.pdf$','.png',figs),'"></a><p>',caption[sub('\\.pdf$','',figs)],'</p><p><a href="figures/main/',figs,'">Open PDF</a></p></article>',collapse='')
writeLines(c('<!doctype html><meta charset="utf-8"><title>Aim 2 frozen-program bypass validation</title><style>body{font:16px/1.45 system-ui;background:#f4f6f8;color:#183044;margin:auto;max-width:1100px;padding:24px}article{background:white;border:1px solid #d8e0e6;margin:18px 0;padding:14px}img{display:block;max-width:100%;margin:auto}a{color:#075b8e}</style><h1>Aim 2 · same frozen programs in bypass</h1><p><a href="README.md">Result summary</a> · <a href="frozen_program_genes_exact_copy.tsv">Exact frozen genes</a> · <a href="../../scripts/35_validate_frozen_programs_in_bypass.R">Analysis code</a></p>',cards),file.path(out,'index.html'))
message('Stage 35 complete: ',out)
