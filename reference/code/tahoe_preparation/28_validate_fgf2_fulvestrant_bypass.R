#!/usr/bin/env Rscript
## Held-out validation. Discovery programs are loaded and frozen before counts.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggpubr); library(edgeR); library(ComplexHeatmap)
})
source('/work/src/R/lib_heatmap.R')
args <- commandArgs(trailingOnly=TRUE)
iteration <- sub('^--iteration=', '', args[grepl('^--iteration=', args)])
if(!length(iteration)) iteration <- '01'
iteration <- sprintf('%02d',as.integer(iteration))
root <- '/work/tahoe_prism_v2'
out <- file.path(root,'results/28_fgf2_fulvestrant_bypass',paste0('iteration_',iteration))
art <- file.path(root,'artifacts/28_fgf2_fulvestrant_bypass',paste0('iteration_',iteration))
fig <- file.path(out,'figures'); dir.create(fig,recursive=TRUE,showWarnings=FALSE)
frozen_file <- file.path(art,'frozen_programs.rds')
frozen_hash <- unname(tools::md5sum(frozen_file))
frozen <- readRDS(frozen_file)
memb <- as.data.table(frozen$membership)[order(plotting_order)]
genes <- memb$gene; orient <- setNames(memb$orientation,memb$gene)
role_levels <- c('ER','FGF','Root')
role_labels <- c(ER='ER branch',FGF='FGF receptor branch',Root='Shared proliferation root')
role_cols <- c(ER='#7A3E9D',FGF='#D97818',Root='#167C80')

meta <- fread(file.path(root,'artifacts/28_fgf2_fulvestrant_bypass/gse153509_parental_samples.tsv'))
raw <- fread(file.path(root,'artifacts/28_fgf2_fulvestrant_bypass/GSE153509_FGFR.Rescue.seq.counts.QC.csv.gz'))
setnames(raw,1,'gene'); raw <- raw[gene %chin% genes]
coverage <- data.table(gene=genes,measured=genes %chin% raw$gene)
fwrite(coverage,file.path(out,'28c_gse153509_gene_coverage.tsv'),sep='\t')
stopifnot(all(coverage$measured),!anyDuplicated(raw$gene))
sample_cols <- intersect(meta$matrix_column,names(raw)); meta <- meta[match(sample_cols,matrix_column)]
counts <- as.matrix(raw[,..sample_cols]); rownames(counts) <- raw$gene; storage.mode(counts) <- 'integer'
dge <- calcNormFactors(DGEList(counts=counts),method='TMM')
logcpm <- cpm(dge,log=TRUE,prior.count=2)
logcpm <- logcpm[genes,,drop=FALSE]

## Predefined replicate QC; only extreme, low-correlation profiles can be removed.
qc_cor <- cor(logcpm,use='pairwise.complete.obs')
qc <- rbindlist(lapply(seq_len(ncol(qc_cor)),function(i) {
 same <- which(meta$arm==meta$arm[i] & seq_len(nrow(meta))!=i)
 data.table(matrix_column=colnames(qc_cor)[i],arm=meta$arm[i],
   median_within_arm_cor=median(qc_cor[i,same],na.rm=TRUE))
}))
qc_med <- median(qc$median_within_arm_cor); qc_mad <- mad(qc$median_within_arm_cor)
qc[,extreme_outlier:=median_within_arm_cor < min(.85,qc_med-4*qc_mad)]
fwrite(qc,file.path(out,'28c_replicate_qc.tsv'),sep='\t')
keep <- !qc$extreme_outlier
meta <- meta[keep]; logcpm <- logcpm[,keep,drop=FALSE]

z <- t(scale(t(logcpm))); z[!is.finite(z)] <- 0
z <- sweep(z,1,orient[rownames(z)],'*')
scores <- sapply(role_levels,function(r) colMeans(z[memb[role==r,gene],,drop=FALSE]))
score_dt <- cbind(meta[,.(gsm,matrix_column,arm)],as.data.table(scores))
setnames(score_dt,role_levels,unname(role_labels[role_levels]))
fwrite(score_dt,file.path(out,'28c_sample_program_scores.tsv'),sep='\t')

contrasts <- data.table(
 contrast=c('Fulvestrant vs control','FGF2 vs control','FGF2 rescue under fulvestrant',
 'Palbociclib vs control','FGF2 under palbociclib','Palbociclib blocks FGF2 rescue',
 'FGFR inhibition blocks rescue','MEK inhibition blocks rescue',
 'FGFR inhibition after palbociclib','MEK inhibition after palbociclib',
 'FIIN-3 vs control','FGF2 under FIIN-3','FIIN-3 blocks FGF2',
 'FIIN-3 after palbociclib + FGF2'),
 treated=c('Fulvestrant','FGF2','Fulvestrant + FGF2','Palbociclib','Palbociclib + FGF2',
 'Fulvestrant + Palbociclib + FGF2','Fulvestrant + FGF2 + FIIN-3',
 'Fulvestrant + FGF2 + Trametinib','Fulvestrant + Palbociclib + FGF2 + FIIN-3',
 'Fulvestrant + Palbociclib + FGF2 + Trametinib','FIIN-3','FIIN-3 + FGF2',
 'FIIN-3 + FGF2','Palbociclib + FGF2 + FIIN-3'),
 reference=c('Control','Control','Fulvestrant','Control','Palbociclib',
 'Fulvestrant + FGF2','Fulvestrant + FGF2','Fulvestrant + FGF2',
 'Fulvestrant + Palbociclib + FGF2','Fulvestrant + Palbociclib + FGF2',
 'Control','FIIN-3','FGF2','Palbociclib + FGF2'))
res <- rbindlist(lapply(seq_len(nrow(contrasts)),function(i) rbindlist(lapply(seq_along(role_levels),function(j) {
 a <- scores[meta$arm==contrasts$treated[i],j]; b <- scores[meta$arm==contrasts$reference[i],j]
 tt <- t.test(a,b)
 data.table(contrast=contrasts$contrast[i],treated=contrasts$treated[i],reference=contrasts$reference[i],
  role=unname(role_labels[role_levels[j]]),effect=mean(a)-mean(b),ci_low=tt$conf.int[1],ci_high=tt$conf.int[2],
  p_value=tt$p.value,n_treated=length(a),n_reference=length(b))
})) ))
res[,q_value:=p.adjust(p_value,'BH')]
res[,contrast:=factor(contrast,levels=rev(contrasts$contrast))]
fwrite(res,file.path(out,'28c_program_contrasts.tsv'),sep='\t')

theme_grant <- theme_classic(base_size=12)+theme(plot.title=element_text(face='bold',size=13),
 axis.title=element_text(face='bold'),strip.background=element_rect(fill='#F2F2F2',colour=NA),
 strip.text=element_text(face='bold',size=10),legend.position='none')
save_plot <- function(p,name,w=7,h=4.5) {
 ggsave(file.path(fig,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)
 ggsave(file.path(fig,paste0(name,'.png')),p,width=w,height=h,dpi=220,bg='white')
}
save_ht <- function(ht,name,w=7,h=5.8) {
 cairo_pdf(file.path(fig,paste0(name,'.pdf')),w,h); draw(ht,heatmap_legend_side='right',annotation_legend_side='right'); dev.off()
 png(file.path(fig,paste0(name,'.png')),w,h,units='in',res=220,bg='white'); draw(ht,heatmap_legend_side='right',annotation_legend_side='right'); dev.off()
}
p <- ggplot(res,aes(effect,contrast,colour=role))+
 geom_vline(xintercept=0,colour='#777777',linewidth=.45)+
 geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',width=.18,linewidth=.7)+
 geom_point(aes(shape=q_value<.05),size=2.8)+facet_wrap(~role,nrow=1,scales='free_x')+
 scale_colour_manual(values=setNames(role_cols,unname(role_labels)))+scale_shape_manual(values=c(`TRUE`=19,`FALSE`=1))+
 labs(title='Held-out T47D program displacement',x='Difference in frozen program score',y=NULL)+theme_grant+
 theme(axis.text.y=element_text(size=8),aspect.ratio=1)
save_plot(p,'28c_01_program_contrasts',9.2,5.4)

arm_order <- c('Control','FGF2','FIIN-3','FIIN-3 + FGF2','Fulvestrant',
 'Fulvestrant + FGF2','Fulvestrant + FGF2 + FIIN-3','Fulvestrant + FGF2 + Trametinib',
 'Palbociclib','Palbociclib + FGF2','Palbociclib + FGF2 + FIIN-3',
 'Fulvestrant + Palbociclib + FGF2','Fulvestrant + Palbociclib + FGF2 + FIIN-3',
 'Fulvestrant + Palbociclib + FGF2 + Trametinib')
ord <- order(factor(meta$arm,levels=arm_order),scores[,'Root'])
colann <- data.frame(Condition=factor(meta$arm[ord],levels=arm_order))
rowann <- data.frame(Role=factor(role_labels[memb$role],levels=unname(role_labels)),
 Program=factor(memb$program,levels=unique(memb$program)))
ht <- draw_heatmap(z[,ord,drop=FALSE],col_anno_df=colann,col_anno_cols='Condition',
 row_anno_df=rowann,row_anno_cols=c('Role','Program'),palette_overrides=list(Role=setNames(role_cols,unname(role_labels))),
 scale_rows=FALSE,value_range=c(-2.5,0,2.5),value_colors=c('#2166AC','#FFFFFF','#B2182B'),hm_name='Gene z',
 show_row_names=TRUE,show_column_names=FALSE,cluster_rows=FALSE,cluster_columns=FALSE,
 row_split=rowann$Program,column_split=colann$Condition,row_gap=unit(1.2,'mm'),column_gap=unit(.6,'mm'),
 row_names_gp=grid::gpar(fontsize=6,fontface='bold'),title='GSE153509: frozen genes × parental T47D samples',border=TRUE)
save_ht(ht,'28c_02_gse_gene_sample',8.0,6.0)

arm_means <- sapply(arm_order,function(a) colMeans(scores[meta$arm==a,,drop=FALSE]))
rownames(arm_means) <- unname(role_labels[role_levels]); colnames(arm_means) <- arm_order
lim <- max(abs(quantile(arm_means,c(.02,.98),na.rm=TRUE)))
ht <- draw_heatmap(arm_means,scale_rows=FALSE,value_range=c(-lim,0,lim),value_colors=c('#2166AC','#FFFFFF','#B2182B'),
 hm_name='Mean score',show_row_names=TRUE,show_column_names=TRUE,cluster_rows=FALSE,cluster_columns=FALSE,
 row_names_gp=grid::gpar(fontsize=9,fontface='bold'),column_names_gp=grid::gpar(fontsize=7,fontface='bold'),
 column_names_rot=45,title='Mean frozen-program activity by treatment',border=TRUE)
save_ht(ht,'28c_03_condition_program_heatmap',7.4,3.8)

sd <- data.table(FGF=scores[,'FGF'],Root=scores[,'Root'],arm=meta$arm)
focal <- c('Control','Fulvestrant','FGF2','Fulvestrant + FGF2','Fulvestrant + Palbociclib + FGF2')
sd <- sd[arm %chin% focal]; sd[,arm:=factor(arm,levels=focal)]
p <- ggplot(sd,aes(FGF,Root,colour=arm))+geom_point(size=2.5,alpha=.85)+
 geom_smooth(aes(group=1),method='lm',se=TRUE,colour='#333333',fill='#C9C9C9',linewidth=.8)+
 ggpubr::stat_cor(aes(group=1),method='pearson',label.x.npc='left',label.y.npc='top',size=3.7,colour='#222222')+
 scale_colour_brewer(palette='Dark2')+labs(title='FGF branch and shared root move together',
 x='FGF receptor branch score',y='Shared proliferation root score',colour='Treatment')+
 theme_grant+theme(legend.position='right',aspect.ratio=1)
save_plot(p,'28c_04_fgf_root_coupling',5.3,4.5)

geteff <- function(con,rr) res[as.character(contrast)==con & role==role_labels[rr],effect]
nodes <- data.table(label=c('ER branch','FGF receptor branch','Shared proliferation root','Viable-cell output'),x=c(0,2,1,1),y=c(2,2,1,0))
edges <- data.table(x=c(0,2,1),y=c(2,2,1),xend=c(1,1,1),yend=c(1.12,1.12,.18))
p <- ggplot()+geom_segment(data=edges,aes(x,y,xend=xend,yend=yend),linewidth=1.2,colour='#4D4D4D',arrow=arrow(length=unit(3,'mm'),type='closed'))+
 geom_point(data=nodes,aes(x,y),shape=21,size=12,fill=c(role_cols,'#555555'),stroke=1.2,colour='white')+
 geom_label(data=nodes,aes(x,y,label=label),fontface='bold',size=3.4,linewidth=0,fill='white',alpha=.9)+
 annotate('label',x=2.15,y=1.62,label=sprintf('FGF2 rescue\nFGF %+0.2f | root %+0.2f',
  geteff('FGF2 rescue under fulvestrant','FGF'),geteff('FGF2 rescue under fulvestrant','Root')),size=3.4,fill='#FFF3E8',linewidth=0)+
 annotate('label',x=1.72,y=.58,label=sprintf('Add palbociclib\nroot %+0.2f',
  geteff('Palbociclib blocks FGF2 rescue','Root')),size=3.4,fill='#E7F3F2',linewidth=0)+
 coord_cartesian(xlim=c(-.55,2.85),ylim=c(-.2,2.5),clip='off')+labs(title='Held-out molecular test of the bypass hierarchy')+
 theme_void(base_size=12)+theme(plot.title=element_text(face='bold',hjust=.5,size=13),plot.margin=margin(6,10,6,10))
save_plot(p,'28c_05_validated_y_hierarchy',5.4,4.4)

criteria <- data.table(
 criterion=c('FGF branch activates during FGF2 rescue','Shared root restores during FGF2 rescue',
 'Palbociclib suppresses the restored root'),
 contrast=c('FGF2 rescue under fulvestrant','FGF2 rescue under fulvestrant','Palbociclib blocks FGF2 rescue'),
 role=c('FGF receptor branch','Shared proliferation root','Shared proliferation root'),
 expected_sign=c(1,1,-1))
criteria <- merge(criteria,res[,.(contrast=as.character(contrast),role,effect,q_value)],by=c('contrast','role'),all.x=TRUE)
criteria[,passed:=sign(effect)==expected_sign & q_value<.05]
fwrite(criteria,file.path(out,'28c_primary_hypothesis_tests.tsv'),sep='\t')
writeLines(c(sprintf('frozen_program_md5\t%s',frozen_hash),sprintf('n_parental_samples\t%d',nrow(meta)),
 sprintf('n_extreme_outliers_removed\t%d',sum(qc$extreme_outlier)),sprintf('all_primary_criteria_passed\t%s',all(criteria$passed))),
 file.path(out,'28c_validation_manifest.tsv'))
