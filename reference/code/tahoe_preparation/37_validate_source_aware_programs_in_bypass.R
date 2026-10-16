#!/usr/bin/env Rscript
## Stage 37: use the exact Stage-36 programs in explicit EGF bypass,
## phosphoproteomic, and matched combination-viability tests.
suppressPackageStartupMessages({library(data.table);library(matrixStats);library(ggplot2);library(ggpubr)})
root<-'/work/tahoe_prism_v2';work<-dirname(root);input_art<-file.path(root,'artifacts/26_braf_erbb_bypass')
src36<-file.path(root,'artifacts/36_source_aware_braf_erbb_programs');src33<-file.path(root,'artifacts/33_braf_erbb_programs')
out<-file.path(root,'results/37_source_aware_bypass_validation');art<-file.path(root,'artifacts/37_source_aware_bypass_validation')
unlink(out,recursive=TRUE,force=TRUE);unlink(art,recursive=TRUE,force=TRUE);dir.create(file.path(out,'figures/main'),recursive=TRUE);dir.create(file.path(out,'figures/supplementary'),recursive=TRUE);dir.create(art,recursive=TRUE)
source(file.path(work,'src/R/lib_heatmap.R'));source(file.path(root,'src/R/plot_braf_erbb_pathway.R'));main<-file.path(out,'figures/main');supp<-file.path(out,'figures/supplementary')
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
write_table<-function(x,n)fwrite(as.data.table(x),file.path(out,n),sep='\t',na='NA')
robust_limit<-function(x,p=.98)max(as.numeric(quantile(abs(x[is.finite(x)]),p,na.rm=TRUE)),.01)

## 1. Transfer the exact frozen definition --------------------------------
frozen<-readRDS(file.path(src36,'frozen_source_aware_programs.rds'));membership<-as.data.table(frozen$membership);setorder(membership,plotting_order)
programs<-c('P1','P2','P3','P4');genes<-membership$gene;gene_sign<-frozen$gene_sign[genes]
program_names<-frozen$program_names;program_colors<-c(P1='#E41A1C',P2='#377EB8',P3='#4DAF4A',P4='#984EA3');program_fills<-c(P1='#F6B7B4',P2='#B8D0EA',P3='#B9DFB8',P4='#D8C1E5')
check(identical(genes,as.character(frozen$leaf_order))&&length(genes)==58L,'Stage-36 genes or order changed during transfer')
write_table(membership,'frozen_program_genes_exact_copy.tsv')
v<-readRDS(file.path(src33,'verification_inputs.rds'));m127<-as.data.table(v$s127);m684<-as.data.table(v$s684)

## 2. Re-extract all 58 genes from the raw processed external matrices -----
hgnc<-fread(file.path(work,'outputs/dataset_inventory/cache/cellfie/hgnc_complete_set.txt'));hgnc<-hgnc[status=='Approved']
lookup<-rbindlist(lapply(seq_len(nrow(hgnc)),function(i){a<-unique(c(hgnc$symbol[i],unlist(strsplit(paste(hgnc$prev_symbol[i],hgnc$alias_symbol[i],sep='|'),'\\|'))));data.table(input=a[nzchar(a)&!is.na(a)&a!='NA'],symbol=hgnc$symbol[i])}))
lookup<-lookup[,if(uniqueN(symbol)==1L).(symbol=symbol[1])else .(symbol=NA_character_),by=input]
canonical<-function(x){a<-lookup$symbol[match(x,lookup$input)];exact<-x%chin%hgnc$symbol;a[exact]<-x[exact];a}
official<-canonical(genes);check(!anyNA(official),'A frozen gene lacks an unambiguous HGNC symbol')

# GSE127988: author RPM matrix, individual columns only.
rpm_path<-file.path(input_art,'gse127988/GSE127988_RPM_Gerosa_et_al.csv.gz');rh<-fread(rpm_path,header=FALSE,nrows=6,check.names=FALSE);header<-as.matrix(rh[,-c(1,2),with=FALSE])
ci<-vapply(m127$sample_id,function(id){z<-which(header[4,]=='Individual'&header[5,]==as.character(id));check(length(z)==1L,'Ambiguous GSE127988 individual sample');z},integer(1))
rpm<-fread(rpm_path,skip=7,header=FALSE,select=c(1L,2L,ci+2L),check.names=FALSE);rpm_x<-as.matrix(rpm[,-c(1,2),with=FALSE]);storage.mode(rpm_x)<-'double'
rna_symbol<-hgnc$symbol[match(sub('\\..*','',rpm[[1]]),hgnc$ensembl_gene_id)];fallback<-canonical(rpm[[2]]);rna_symbol[is.na(rna_symbol)]<-fallback[is.na(rna_symbol)]
rna<-data.table(official_symbol=rna_symbol,rpm_x);rna<-rna[!is.na(official_symbol),lapply(.SD,sum),by=official_symbol]
x127_raw<-log2(1+as.matrix(rna[,-1]));rownames(x127_raw)<-rna$official_symbol;colnames(x127_raw)<-m127$gsm

# GSE68453: normalized log2 microarray expression, best-IQR probe per gene.
platform<-fread(file.path(input_art,'GSE68453_platform_annotation.tsv'),quote='');platform[,official_symbol:=canonical(Symbol)]
mat<-suppressWarnings(fread(file.path(input_art,'GSE68453_series_matrix.txt.gz'),skip='"ID_REF"',header=TRUE,check.names=FALSE,select=c('ID_REF',m684$gsm)));mat<-mat[ID_REF%chin%platform$ID]
array_x<-as.matrix(mat[,m684$gsm,with=FALSE]);storage.mode(array_x)<-'double';probe<-data.table(probe=mat$ID_REF,official_symbol=platform$official_symbol[match(mat$ID_REF,platform$ID)],iqr=rowIQRs(array_x),array_row=seq_len(nrow(array_x)));probe<-probe[!is.na(official_symbol)];setorder(probe,official_symbol,-iqr,probe);probe<-probe[,head(.SD,1),by=official_symbol]
x684_raw<-array_x[probe$array_row,,drop=FALSE];rownames(x684_raw)<-probe$official_symbol

extract_frozen<-function(x){m<-matrix(NA_real_,length(genes),ncol(x),dimnames=list(genes,colnames(x)));ok<-official%chin%rownames(x);m[ok,]<-x[official[ok],,drop=FALSE];m}
x127<-extract_frozen(x127_raw);x684<-extract_frozen(x684_raw)
coverage<-rbindlist(lapply(c('GSE127988','GSE68453'),function(ds){x<-if(ds=='GSE127988')x127 else x684;rbindlist(lapply(programs,function(p){g<-membership[program==p,gene];ok<-rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x);data.table(dataset=ds,program=p,program_name=program_names[p],n_frozen=length(g),n_measured=sum(ok),missing_genes=paste(g[!ok],collapse=', '))}))}))
write_table(coverage,'external_gene_coverage.tsv')

score_programs<-function(x)sapply(programs,function(p){g<-membership[program==p,gene];present<-g[rowSums(is.finite(x[g,,drop=FALSE]))==ncol(x)];check(length(present)>=.8*length(g),paste('Insufficient measured genes for',p));colMeans(sweep(x[present,,drop=FALSE],1,gene_sign[present],'*'))})

## 3. Explicit EGF bypass in A375 (GSE127988) -----------------------------
# Each time point is expressed relative to the matching 0-hour sample under
# the same drug background and biological day.
d127<-x127
for(a in unique(m127$arm))for(day in unique(m127$biological_set)){jj<-which(m127$arm==a&m127$biological_set==day);ref<-rowMeans(x127[,m127$arm==a&m127$biological_set==day&m127$egf_h==0,drop=FALSE],na.rm=TRUE);d127[,jj]<-x127[,jj,drop=FALSE]-ref}
s127<-score_programs(d127);long127<-melt(data.table(m127[,.(gsm,arm,egf_h,biological_set)],s127),id.vars=c('gsm','arm','egf_h','biological_set'),variable.name='program',value.name='score');long127[,program:=as.character(program)]
long127[,drug_background:=factor(ifelse(arm=='Vem','Vemurafenib','Vemurafenib + cobimetinib'),levels=c('Vemurafenib','Vemurafenib + cobimetinib'))]
day127<-long127[,.(score=mean(score),se=sd(score)/sqrt(.N)),by=.(drug_background,egf_h,biological_set,program)]
trapz<-function(t,y)sum(diff(t)*(head(y,-1)+tail(y,-1))/2)
auc127<-day127[order(egf_h),.(egf_response_auc=trapz(egf_h,score)),by=.(drug_background,biological_set,program)]
auc_wide<-dcast(auc127,program+biological_set~drug_background,value.var='egf_response_auc');setnames(auc_wide,c('Vemurafenib','Vemurafenib + cobimetinib'),c('vemurafenib_auc','vemurafenib_cobimetinib_auc'));auc_wide[,cobimetinib_change:=vemurafenib_cobimetinib_auc-vemurafenib_auc]
check(all(auc_wide[program%chin%c('P2','P3'),vemurafenib_auc]>0)&&all(auc_wide[program%chin%c('P2','P3'),cobimetinib_change]<0),'Explicit EGF bypass pattern failed')
write_table(long127,'gse127988_program_scores.tsv');write_table(day127,'gse127988_time_course_means.tsv');write_table(auc_wide,'gse127988_egf_response_auc.tsv')

## 4. Molecular-to-functional prediction in GSE68453 ---------------------
d684<-x684
for(line in unique(m684$cell_line)){jj<-which(m684$cell_line==line);ref<-rowMeans(x684[,m684$cell_line==line&m684$arm=='control',drop=FALSE],na.rm=TRUE);d684[,jj]<-x684[,jj,drop=FALSE]-ref}
s684<-score_programs(d684);long684<-melt(data.table(m684[,.(gsm,cell_line,arm,replicate)],s684),id.vars=c('gsm','cell_line','arm','replicate'),variable.name='program',value.name='score');long684[,program:=as.character(program)]
means684<-long684[,.(score=mean(score)),by=.(cell_line,arm,program)];factorial684<-dcast(means684,cell_line+program~arm,value.var='score');factorial684[,interaction:=plx_lap-plx-lapatinib+control]
bliss<-data.table(cell_line=c('A375','SKMEL24','HT144','DM331','SKMEL28','DM13','VMM5A','VMM18','DM6','VMM15','VMM12','VMM17'),bliss=c(5.7,29.3,16.3,27.3,49.7,39.3,43,31.7,-2.7,-8.7,3.7,-1))
functional<-merge(factorial684[is.finite(interaction)],bliss,by='cell_line');functional[,extra_program_reduction:=-interaction]
functional_stats<-functional[,{tt<-cor.test(extra_program_reduction,bliss);.(n=.N,pearson=unname(tt$estimate),p_value=tt$p.value)},by=program]
check(functional_stats[program=='P3',pearson]>=.4,'Shared MEK-ERK output fails the functional direction')
write_table(factorial684,'gse68453_program_factorial.tsv');write_table(functional,'gse68453_program_function.tsv');write_table(functional_stats,'gse68453_functional_prediction_stats.tsv')

## 5. Gerosa phosphoproteomics tests the same pathway edge ----------------
gerosa<-fread(file.path(root,'results/33_braf_erbb_programs/gerosa_panorama_phosphopeptide_summary.tsv'));gtest<-fread(file.path(root,'results/33_braf_erbb_programs/gerosa_panorama_prediction_tests.tsv'))
check(all(gtest[grepl('^ERK',analyte),EGF_rescue_over_one_uM_vem]>=3),'Gerosa EGF-to-pERK rescue failed')
write_table(gerosa,'gerosa_phosphoproteomic_summary.tsv');write_table(gtest,'gerosa_phosphoproteomic_tests.tsv')

## 6. Grant figures: pathway first, then each test -------------------------
# Figure 1: precise biological mechanism and exact drug targets only.
p<-build_braf_erbb_pathway(show_plx4720=TRUE)
ggsave(file.path(main,'01_biological_pathway_and_drugs.pdf'),p,width=7.4,height=3.75,device=cairo_pdf,bg='white')
ggsave(file.path(main,'01_biological_pathway_and_drugs.png'),p,width=7.4,height=3.75,dpi=240,bg='white')

# Figure 2: one complete, compact schema for Aim 2. The tests are parallel.
input_box<-data.table(xmin=3.35,xmax=6.57,ymin=3.84,ymax=4.48)
test_boxes<-data.table(xmin=c(.18,5.10),xmax=c(4.82,9.74),ymin=.66,ymax=3.34,fill=c('#EFF5F9','#F3EFF7'),border=c('#276C9B','#5B4775'))
arrows<-data.table(x=c(4.96,4.96),y=c(3.84,3.84),xend=c(2.50,7.42),yend=c(3.39,3.39))
left_steps<-data.table(x=c(1.38,3.62),label=c('Vemurafenib
BRAF inhibitor, 24 h','EGF
0-8 h'))
right_steps<-data.table(x=c(6.20,8.63),label=c('P3 RNA score
8 h','Viability
3 d'))
step_arrows<-data.table(x=c(1.94,6.86),y=c(2.20,1.83),xend=c(3.02,7.96),yend=c(2.20,1.83),col=c('#276C9B','#5B4775'))
p<-ggplot()+geom_rect(data=input_box,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='#E8EFEA',color='#24745E',linewidth=.8)+annotate('text',x=4.96,y=4.25,label='FROZEN FROM AIM 1',size=3,fontface='bold',color='#24745E')+annotate('text',x=4.96,y=4.00,label='Program genes, signs and weights',size=2.9,color='#213746')+geom_segment(data=arrows,aes(x,y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(1.8,'mm'),type='closed'),linewidth=.8,color='#526875')+geom_rect(data=test_boxes,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=test_boxes$fill,color=test_boxes$border,linewidth=.8)+annotate('text',x=2.50,y=3.08,label='LIGAND BYPASS TEST',size=3.15,fontface='bold',color='#276C9B')+annotate('text',x=2.50,y=2.77,label='A375, BRAF V600E, GSE127988',size=2.9,fontface='bold',color='#213746')+geom_label(data=left_steps,aes(x=x,y=2.20,label=label),size=2.75,fontface='bold',lineheight=.9,fill='white',color='#276C9B',label.padding=grid::unit(.22,'lines'),linewidth=.6)+geom_segment(data=step_arrows[1],aes(x,y,xend=xend,yend=yend),arrow=grid::arrow(length=grid::unit(1.6,'mm'),type='closed'),linewidth=.8,color='#276C9B')+annotate('text',x=3.62,y=1.75,label='repeat with cobimetinib (MEK inhibitor)',size=2.5,fontface='bold',color='#A33870')+annotate('text',x=2.50,y=1.38,label='RNA programs; ERK phosphorylation, Gerosa et al.',size=2.65,color='#213746')+annotate('text',x=2.50,y=1.06,label='Expected: P2 and P3 increase;\ncobimetinib blocks the increase',size=2.55,fontface='bold',lineheight=.9,color='#174C70')+annotate('text',x=7.42,y=3.08,label='FUNCTIONAL TEST',size=3.15,fontface='bold',color='#5B4775')+annotate('text',x=7.42,y=2.77,label='Five BRAF-mutant melanoma lines, GSE68453',size=2.9,fontface='bold',color='#213746')+annotate('text',x=7.42,y=2.35,label='PLX4720 (BRAF inhibitor)\nwith or without lapatinib (ERBB inhibitor)\nno EGF added',size=2.65,lineheight=.9,color='#213746')+geom_label(data=right_steps,aes(x=x,y=1.76,label=label),size=2.8,fontface='bold',lineheight=.9,fill='white',color='#5B4775',label.padding=grid::unit(.22,'lines'),linewidth=.6)+geom_segment(data=step_arrows[2],aes(x,y=1.76,xend=xend,yend=1.76),arrow=grid::arrow(length=grid::unit(1.6,'mm'),type='closed'),linewidth=.8,color='#5B4775')+annotate('text',x=7.42,y=1.24,label='Expected: more P3 suppression,\nmore combination benefit',size=2.55,fontface='bold',lineheight=.9,color='#4B3862')+annotate('rect',xmin=1.08,xmax=8.84,ymin=.08,ymax=.44,fill='white',color='#24745E',linewidth=.8)+annotate('text',x=4.96,y=.26,label='Shared readout: frozen P3 connects pathway reactivation to viability benefit',size=2.8,fontface='bold',color='#214E41')+coord_cartesian(xlim=c(0,9.92),ylim=c(-.02,4.68),clip='off')+theme_void(base_size=13)+theme(plot.title=element_text(size=15,face='bold'),plot.subtitle=element_text(size=10.2,color='#40535E'),plot.margin=margin(7,8,7,8))+labs(title='Aim 2: test receptor-driven bypass of BRAF inhibition',subtitle='One frozen program definition is evaluated in two complementary external experiments')
ggsave(file.path(main,'02_aim2_experimental_design.pdf'),p,width=7.4,height=4.15,device=cairo_pdf,bg='white');ggsave(file.path(main,'02_aim2_experimental_design.png'),p,width=7.4,height=4.15,dpi=240,bg='white')

# Figure 3: P2 and P3 genes in GSE127988; the complete four-program audit is supplementary.
aim2_programs<-c('P2','P3')
set.seed(3701);rnd<-runif(nrow(m127));ord<-order(match(m127$arm,c('Vem','Vem + Cobi')),m127$biological_set,m127$egf_h,rnd)
signed127<-sweep(d127,1,gene_sign,'*')
condition_levels<-c('Vemurafenib, day 3','Vemurafenib, day 4','Vemurafenib + cobimetinib, day 3','Vemurafenib + cobimetinib, day 4');condition_colors<-c('#F2C690','#D88736','#C5AFDA','#7656A5');names(condition_colors)<-condition_levels
ann<-data.frame('Drug condition'=factor(paste(ifelse(m127$arm=='Vem','Vemurafenib','Vemurafenib + cobimetinib'),paste0('day ',m127$biological_set),sep=', '),levels=condition_levels),'Hours after EGF'=factor(paste0(m127$egf_h,' h'),levels=paste0(c(0,.5,1,2,3,4,8),' h')),check.names=FALSE)
col_group<-factor(paste(ifelse(m127$arm=='Vem','Vem','Vem + Cobi'),paste0('day ',m127$biological_set),sep=', '),levels=c('Vem, day 3','Vem, day 4','Vem + Cobi, day 3','Vem + Cobi, day 4'))
render127<-function(base,show_all=FALSE){plot_genes<-if(show_all)genes else membership[program%chin%aim2_programs,gene];mat<-signed127[plot_genes,ord,drop=FALSE];mi<-membership[match(plot_genes,gene)];row_levels<-if(show_all)programs else aim2_programs;row_group<-factor(mi$program,levels=row_levels,labels=paste0(row_levels,': ',program_names[row_levels]));annr<-data.frame(Program=factor(mi$program,levels=programs));lim<-robust_limit(mat,.975);key<-c('FOS','EREG','DUSP6','HBEGF','MAPK1','ERRFI1','SPRY2');row_labels<-if(show_all)plot_genes else ifelse(plot_genes%chin%key,plot_genes,'');h<-draw_heatmap(mat,col_anno_df=ann[ord,,drop=FALSE],col_anno_cols=names(ann),row_anno_df=annr,row_anno_cols='Program',palette_overrides=list(Program=program_colors,'Drug condition'=condition_colors,'Hours after EGF'=setNames(hcl.colors(7,'Blues 3'),levels(ann$`Hours after EGF`))),ann_row_args=list(show_annotation_name=FALSE),scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=seq_len(nrow(mat)),column_order=seq_len(ncol(mat)),row_split=row_group,column_split=col_group[ord],cluster_row_slices=FALSE,cluster_column_slices=FALSE,row_gap=grid::unit(.75,'mm'),column_gap=grid::unit(.9,'mm'),row_title=NULL,show_row_names=TRUE,show_column_names=FALSE,row_labels=row_labels,row_names_gp=grid::gpar(fontsize=if(show_all)5.5 else 8.5,fontface='bold'),hm_name='Signed RNA response',na_col='#B9BEC1',heatmap_legend_param=list(title='Signed log2 RNA change
relative to 0 h',at=round(c(-lim,0,lim),1)));draw1<-function()ComplexHeatmap::draw(h,heatmap_legend_side='right',annotation_legend_side='right',column_title=if(show_all)'Complete frozen-program audit' else 'P2 and P3 genes after EGF in A375',column_title_gp=grid::gpar(fontsize=13,fontface='bold'));wh<-if(show_all)c(8,6.7) else c(7.2,4.6);cairo_pdf(paste0(base,'.pdf'),wh[1],wh[2]);draw1();dev.off();png(paste0(base,'.png'),wh[1],wh[2],units='in',res=220,type='cairo');draw1();dev.off()}
render127(file.path(main,'03_gse127988_bypass_gene_heatmap'));render127(file.path(supp,'S01_gse127988_all_program_genes'),TRUE)
write_table(data.table(gsm=m127$gsm[ord],plotting_order=seq_along(ord),drug_condition=ann$`Drug condition`[ord],egf_time=ann$`Hours after EGF`[ord],biological_day=m127$biological_set[ord]),'gse127988_heatmap_column_order.tsv')

# Figure 4: collapse the EGF time course to one mean response per day and drug background.
display_programs<-c('P1','P2','P3');collapsed127<-copy(auc127)[program%chin%display_programs];collapsed127[,mean_egf_response:=egf_response_auc/8];collapsed127[,condition:=factor(paste(drug_background,paste0('Day ',biological_set),sep='|'),levels=c('Vemurafenib|Day 3','Vemurafenib|Day 4','Vemurafenib + cobimetinib|Day 3','Vemurafenib + cobimetinib|Day 4'))]
write_table(collapsed127,'gse127988_collapsed_program_response.tsv')
wide_prog<-dcast(collapsed127,program~condition,value.var='mean_egf_response');program_matrix<-as.matrix(wide_prog[,-1]);rownames(program_matrix)<-c('P1: RAF-associated response','P2: EGF-EGFR response','P3: shared MEK-ERK output');colnames(program_matrix)<-rep(c('Day 3','Day 4'),2);plim<-max(abs(program_matrix));program_ann<-data.frame(Program=factor(display_programs,levels=programs));drug_split<-factor(rep(c('Vemurafenib','Vemurafenib + cobimetinib'),each=2),levels=c('Vemurafenib','Vemurafenib + cobimetinib'))
hprog<-draw_heatmap(program_matrix,row_anno_df=program_ann,row_anno_cols='Program',palette_overrides=list(Program=program_colors),ann_row_args=list(show_annotation_name=FALSE,show_legend=FALSE),scale_rows=FALSE,value_range=c(-plim,0,plim),value_colors=c('#2166AC','white','#B2182B'),cluster_rows=FALSE,cluster_columns=FALSE,row_order=1:3,column_order=1:4,column_split=drug_split,cluster_column_slices=FALSE,column_gap=grid::unit(1.2,'mm'),show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=10.5,fontface='bold'),column_names_gp=grid::gpar(fontsize=9.5,fontface='bold'),column_names_rot=0,title=c('Vemurafenib','Vemurafenib +\ncobimetinib'),column_title_gp=grid::gpar(fontsize=9.5,fontface='bold'),hm_name='Mean EGF response',heatmap_legend_param=list(title='Mean signed RNA response
after EGF, 0-8 h',at=round(c(-plim,0,plim),2)))
drawprog<-function()ComplexHeatmap::draw(hprog,heatmap_legend_side='right',annotation_legend_side='right',column_title='GSE127988, EGF-induced program response',column_title_gp=grid::gpar(fontsize=13,fontface='bold'))
cairo_pdf(file.path(main,'04_gse127988_collapsed_program_activity.pdf'),6.8,2.8);drawprog();dev.off();png(file.path(main,'04_gse127988_collapsed_program_activity.png'),6.8,2.8,units='in',res=240,type='cairo');drawprog();dev.off()

# Supplementary: retain kinetics for audit; main program figure collapses hours.
program_labels<-c(P1='P1: RAF network',P2='P2: EGF-EGFR response',P3='P3: MEK-ERK output',P4='P4: proliferation');time_plot<-day127[program%chin%aim2_programs];time_plot[,program_label:=factor(program,levels=aim2_programs,labels=program_labels[aim2_programs])]
p<-ggplot(time_plot,aes(egf_h,score,color=drug_background,group=interaction(drug_background,biological_set)))+geom_hline(yintercept=0,color='#B8C1C6')+geom_line(aes(linetype=factor(biological_set)),linewidth=.85)+geom_point(size=2)+facet_wrap(~program_label,nrow=1,scales='free_y')+scale_color_manual(values=c('Vemurafenib'='#D88736','Vemurafenib + cobimetinib'='#7656A5'))+scale_linetype_manual(values=c('3'='solid','4'='dashed'),labels=c('Day 3','Day 4'))+theme_classic(base_size=12)+theme(aspect.ratio=1,strip.text=element_text(face='bold',size=10),legend.position='bottom',legend.title=element_text(face='bold'),plot.title=element_text(size=14),plot.subtitle=element_text(size=10))+guides(color=guide_legend(nrow=1),linetype=guide_legend(nrow=1))+labs(x='Hours after 100 ng/mL EGF',y='Signed program RNA change
relative to 0 h',color='Drug background',linetype='Biological day',title='Time-course audit of the two bypass programs',subtitle='The main figure summarizes each curve over the full 0-8 h interval')
ggsave(file.path(supp,'S02_gse127988_bypass_program_time_courses.pdf'),p,width=6.3,height=3.7,device=cairo_pdf,bg='white');ggsave(file.path(supp,'S02_gse127988_bypass_program_time_courses.png'),p,width=6.3,height=3.7,dpi=220,bg='white')

# Figure 5: direct summary for the two bypass programs.
auc_plot<-auc127[program%chin%c('P2','P3')];auc_plot[,program_label:=factor(program,levels=c('P2','P3'),labels=program_labels[c('P2','P3')])]
p<-ggplot(auc_plot,aes(drug_background,egf_response_auc,group=biological_set,color=program))+geom_hline(yintercept=0,color='#B8C1C6')+geom_line(linewidth=.9,color='#53656F')+geom_point(size=3)+facet_wrap(~program_label,nrow=1,scales='free_y')+scale_color_manual(values=program_colors,guide='none')+scale_x_discrete(labels=c('Vemurafenib'='BRAF\ninhibitor','Vemurafenib + cobimetinib'='BRAF + MEK\ninhibitors'))+theme_classic(base_size=13)+theme(aspect.ratio=1,strip.text=element_text(face='bold'),axis.text.x=element_text(size=10),plot.title=element_text(size=14),plot.subtitle=element_text(size=10))+labs(x='Drug background before EGF addition',y='Total EGF-induced\nRNA response',title='Cobimetinib blocks both EGF-bypass programs',subtitle='The receptor response and shared MEK-ERK output fall on both biological days')
ggsave(file.path(main,'05_cobimetinib_blocks_egf_bypass.pdf'),p,width=6.8,height=3.8,device=cairo_pdf,bg='white');ggsave(file.path(main,'05_cobimetinib_blocks_egf_bypass.png'),p,width=6.8,height=3.8,dpi=190,bg='white')

# Figure 6: matched protein measurement of the EGF-to-ERK edge.
gp<-gerosa[condition%chin%c('DMSO','1 µM vem','1 µM vem + EGF')];gp[,condition:=factor(condition,levels=c('DMSO','1 µM vem','1 µM vem + EGF'),labels=c('Control','BRAF\ninhibitor','BRAF inhibitor\n+ EGF'))];gp[,analyte:=factor(analyte,levels=c('ERK1 pT202/pY204','ERK2 pT185/pY187'),labels=c('Activated ERK1','Activated ERK2'))]
p<-ggplot(gp[!is.na(analyte)],aes(condition,mean_relative,group=analyte,color=analyte))+geom_hline(yintercept=1,color='#B8C1C6',linetype='dashed')+geom_errorbar(aes(ymin=pmax(mean_relative-se_relative,0),ymax=mean_relative+se_relative),width=.08,na.rm=TRUE)+geom_line(linewidth=.85)+geom_point(size=3)+facet_wrap(~analyte,nrow=1)+scale_color_manual(values=c('#377EB8','#4DAF4A'),guide='none')+theme_classic(base_size=13)+theme(aspect.ratio=1,strip.text=element_text(face='bold'),axis.text.x=element_text(size=10),plot.title=element_text(size=15),plot.subtitle=element_text(size=10))+labs(x=NULL,y='Phosphopeptide abundance\n(relative to untreated control)',title='EGF restores ERK activity after BRAF inhibition',subtitle='Gerosa A375 phosphoproteomics, independent protein-level bridge')
ggsave(file.path(main,'06_gerosa_egf_restores_erk_activity.pdf'),p,width=6.6,height=3.8,device=cairo_pdf,bg='white');ggsave(file.path(main,'06_gerosa_egf_restores_erk_activity.png'),p,width=6.6,height=3.8,dpi=190,bg='white')

# Figure 7: early P3 molecular shutdown versus later functional benefit.
fp<-functional[program=='P3'];fs<-functional_stats[program=='P3']
p<-ggplot(fp,aes(extra_program_reduction,bliss))+geom_point(size=3,color=program_colors['P3'])+geom_text(aes(label=cell_line),size=3.2,vjust=-.75,check_overlap=TRUE)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#263840',linewidth=.9)+ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=4.2,color='black')+scale_x_continuous(expand=expansion(mult=c(.12,.20)))+scale_y_continuous(expand=expansion(mult=c(.07,.20)))+coord_cartesian(clip='off')+theme_classic(base_size=13)+theme(aspect.ratio=1,plot.title=element_text(size=13.5,lineheight=.95),plot.subtitle=element_text(size=9.2,lineheight=.95),plot.margin=margin(7,15,7,7))+labs(x='P3 suppression at 8 h beyond single-drug effects
positive = greater suppression than expected',y='Viability benefit at 3 d
positive Bliss score = lower viability than expected',title='P3 suppression ranks
combination benefit',subtitle='PLX4720 (BRAF inhibitor) + lapatinib (ERBB inhibitor)
five BRAF-mutant melanoma lines')
ggsave(file.path(main,'07_shared_output_predicts_combination_benefit.pdf'),p,width=5.3,height=5.3,device=cairo_pdf,bg='white');ggsave(file.path(main,'07_shared_output_predicts_combination_benefit.png'),p,width=5.3,height=5.3,dpi=190,bg='white')

## 7. Verification and viewer --------------------------------------------
verification<-data.table(check=c('Exact Stage-36 genes and order','Exact Stage-36 signs','GSE127988 missingness explicit','GSE68453 measures all 58 genes','P2 activated by EGF on both days','P3 activated by EGF on both days','Cobimetinib attenuates P2 on both days','Cobimetinib attenuates P3 on both days','Gerosa ERK rescue at least threefold','P3 functional Pearson at least 0.4','Every main figure has PDF and PNG'),passed=c(identical(genes,as.character(frozen$leaf_order)),identical(gene_sign,frozen$gene_sign[genes]),coverage[dataset=='GSE127988',sum(n_frozen-n_measured)]==2L,coverage[dataset=='GSE68453',sum(n_measured)]==58L,all(auc_wide[program=='P2',vemurafenib_auc]>0),all(auc_wide[program=='P3',vemurafenib_auc]>0),all(auc_wide[program=='P2',cobimetinib_change]<0),all(auc_wide[program=='P3',cobimetinib_change]<0),all(gtest[grepl('^ERK',analyte),EGF_rescue_over_one_uM_vem]>=3),fs$pearson>=.4,setequal(sub('\\.pdf$','',list.files(main,'\\.pdf$')),sub('\\.png$','',list.files(main,'\\.png$')))))
check(all(verification$passed),'Stage 37 verification failed');write_table(verification,'verification.tsv')
saveRDS(list(source_frozen_programs=normalizePath(file.path(src36,'frozen_source_aware_programs.rds')),membership=membership,gene_sign=gene_sign,gse127988_auc=auc_wide,functional_stats=functional_stats),file.path(art,'source_aware_bypass_transfer.rds'))
figs<-sort(list.files(main,pattern='\\.pdf$'));write_table(data.table(figure=file.path('figures/main',figs)),'FIGURE_INDEX.tsv')
caption<-c(
  '01_biological_pathway_and_drugs'='Biological model: EGF activates EGFR and receptor-driven RAS-RAF signaling, which rejoins the pathway at MEK-ERK after mutant BRAF is inhibited.',
  '02_aim2_experimental_design'='Frozen Aim 1 programs are evaluated in two parallel external tests: an added-EGF molecular experiment and a separate combination-viability experiment without added EGF.',
  '03_gse127988_bypass_gene_heatmap'='P2 and P3 genes across all 56 A375 EGF time-course samples. The full four-program audit is supplementary.',
  '04_gse127988_collapsed_program_activity'='Mean EGF-induced response over 0-8 hours. P1 is retained as context; P2 and P3 are the bypass programs supported by this experiment.',
  '05_cobimetinib_blocks_egf_bypass'='Adding MEK inhibition reduces both EGF-induced programs on both biological days.',
  '06_gerosa_egf_restores_erk_activity'='Independent protein measurements confirm that EGF restores activated ERK after BRAF inhibition.',
  '07_shared_output_predicts_combination_benefit'='Across five BRAF-mutant melanoma lines, stronger early suppression of P3 forecasts greater three-day combination benefit.'
)
writeLines(c('# Source-aware bypass validation','', 'The exact Stage 36 program definitions remain frozen. Aim 2 focuses on P2 and P3 because these are the programs supported by the EGF experiment. P4 is retained in the supplementary audit but excluded from the bypass claim because it does not increase after EGF.','',sprintf('In GSE127988, EGF activates the EGF-EGFR response (P2) and shared MEK-ERK output (P3) during vemurafenib on both biological days. Cobimetinib reduces both responses on both days.'),sprintf('Gerosa protein measurements independently confirm the same pathway edge: EGF restores activated ERK after BRAF inhibition. In GSE68453, early suppression of P3 predicts later benefit from PLX4720 (BRAF inhibitor) plus lapatinib (ERBB inhibitor) across five matched melanoma lines (Pearson r = %.2f).',fs$pearson),'','GSE127988 measures 56 of the 58 frozen genes. KIT and EREG are absent from the processed source matrix, shown as missing, and never replaced with zero. GSE68453 measures all 58 genes.','','## Figure sequence','',paste(sprintf('%d. **%s** - %s',seq_along(figs),sub('^[0-9]+_','',sub('\\.pdf$','',figs)),caption[sub('\\.pdf$','',figs)]),collapse='\n'),'', '[Open the figure viewer](index.html), [Exact frozen genes](frozen_program_genes_exact_copy.tsv), [Code](../../scripts/37_validate_source_aware_programs_in_bypass.R)'),file.path(out,'README.md'))
cards<-paste0('<article><h3>',sub('\\.pdf$','',figs),'</h3><a href="figures/main/',figs,'"><img src="figures/main/',sub('\\.pdf$','.png',figs),'"></a><p>',caption[sub('\\.pdf$','',figs)],'</p><p><a href="figures/main/',figs,'">Open PDF</a></p></article>',collapse='')
writeLines(c('<!doctype html><meta charset="utf-8"><title>Source-aware bypass validation</title><style>body{font:16px/1.45 system-ui;background:#f4f6f8;color:#183044;margin:auto;max-width:1100px;padding:24px}article{background:white;border:1px solid #d8e0e6;margin:18px 0;padding:14px}img{display:block;max-width:100%;margin:auto}a{color:#075b8e}</style><h1>Aim 2, source-aware bypass validation</h1><p><a href="README.md">Result summary</a>, <a href="frozen_program_genes_exact_copy.tsv">Exact frozen genes</a>, <a href="../../scripts/37_validate_source_aware_programs_in_bypass.R">Analysis code</a></p>',cards),file.path(out,'index.html'))
message('Stage 37 complete: ',out)
