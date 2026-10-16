#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggpubr); library(ComplexHeatmap)
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
x <- readRDS(file.path(art,'discovery_plot_data.rds'))
memb <- as.data.table(x$membership)[order(plotting_order)]
genes <- memb$gene; orient <- setNames(memb$orientation,memb$gene)
role_levels <- c('ER','FGF','Root')
role_labels <- c(ER='ER branch',FGF='FGF receptor branch',Root='Shared proliferation root')
role_cols <- c(ER='#7A3E9D',FGF='#D97818',Root='#167C80')
theme_grant <- theme_classic(base_size=12)+theme(
  plot.title=element_text(face='bold',size=13),axis.title=element_text(face='bold'),
  strip.background=element_rect(fill='#F2F2F2',colour=NA),strip.text=element_text(face='bold',size=10),
  legend.title=element_text(face='bold'),aspect.ratio=1)
save_plot <- function(p,name,w=4.8,h=4.2) {
  ggsave(file.path(fig,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf)
  ggsave(file.path(fig,paste0(name,'.png')),p,width=w,height=h,dpi=220,bg='white')
}
save_ht <- function(ht,name,w=6.8,h=6.8) {
  cairo_pdf(file.path(fig,paste0(name,'.pdf')),width=w,height=h)
  draw(ht,heatmap_legend_side='right',annotation_legend_side='right',padding=unit(c(3,3,3,3),'mm'))
  dev.off()
  png(file.path(fig,paste0(name,'.png')),width=w,height=h,units='in',res=220,bg='white')
  draw(ht,heatmap_legend_side='right',annotation_legend_side='right',padding=unit(c(3,3,3,3),'mm'))
  dev.off()
}
read_square <- function(name) {
  z <- fread(file.path(out,name)); rn <- z[[1]]; z[[1]] <- NULL
  a <- as.matrix(z); rownames(a) <- rn; a[genes,genes,drop=FALSE]
}
corr_ht <- function(a,title,name) {
  off <- abs(a[row(a)!=col(a)]); lim <- as.numeric(quantile(off,.98,na.rm=TRUE)); lim <- max(lim,.08)
  ann <- data.frame(Role=factor(role_labels[memb$role],levels=unname(role_labels)),
    Program=factor(memb$program,levels=unique(memb$program)))
  ht <- draw_heatmap(a,col_anno_df=ann,col_anno_cols=c('Role','Program'),row_anno_df=ann,row_anno_cols=c('Role','Program'),
    palette_overrides=list(Role=setNames(role_cols,unname(role_labels))),scale_rows=FALSE,
    ann_row_args=list(show_legend=FALSE,show_annotation_name=FALSE),
    value_range=c(-lim,0,lim),value_colors=c('#2166AC','#FFFFFF','#B2182B'),hm_name='Pearson r',
    show_row_names=TRUE,show_column_names=TRUE,cluster_rows=FALSE,cluster_columns=FALSE,title=title,
    row_split=ann$Program,column_split=ann$Program,row_gap=unit(1.2,'mm'),column_gap=unit(1.2,'mm'),
    row_names_gp=grid::gpar(fontsize=5.7,fontface='bold'),column_names_gp=grid::gpar(fontsize=5.7,fontface='bold'),
    column_names_rot=45,row_title_gp=grid::gpar(fontsize=6,fontface='bold'),
    column_title_gp=grid::gpar(fontsize=6,fontface='bold'),border=TRUE,width=unit(90,'mm'),height=unit(90,'mm'),
    heatmap_legend_param=list(at=c(-lim,0,lim),labels=sprintf('%.2f',c(-lim,0,lim))))
  save_ht(ht,name,8.8,5.9)
}
tahoe <- read_square('28b_tahoe_signed_correlations.tsv.gz')
cyto <- read_square('28b_cyto_signed_correlations.tsv.gz')
joint <- read_square('28b_joint_signed_correlations.tsv.gz')
corr_ht(joint,'Integrated Tahoe + CytoSig gene correlation','28b_01_joint_correlation')
corr_ht(tahoe,'Tahoe drug-response gene correlation','28b_02_tahoe_correlation')
corr_ht(cyto,'CytoSig FGF2-response gene correlation','28b_03_cytosig_correlation')

gene_sample_ht <- function(mat,colann,title,name,show_cols=FALSE,w=7.2,h=5.8) {
  mat <- mat[genes,,drop=FALSE]
  rowann <- data.frame(Role=factor(role_labels[memb$role],levels=unname(role_labels)),
    Program=factor(memb$program,levels=unique(memb$program)))
  ht <- draw_heatmap(mat,col_anno_df=colann,col_anno_cols=names(colann),row_anno_df=rowann,row_anno_cols='Program',
    scale_rows=TRUE,
    value_range=c(-2.5,0,2.5),value_colors=c('#2166AC','#FFFFFF','#B2182B'),hm_name='Row z',
    show_row_names=TRUE,show_column_names=show_cols,cluster_rows=FALSE,cluster_columns=FALSE,title=title,
    row_split=rowann$Program,column_split=colann$Dominant_program,
    row_gap=unit(1.2,'mm'),column_gap=unit(1.2,'mm'),
    row_names_gp=grid::gpar(fontsize=6,fontface='bold'),column_names_gp=grid::gpar(fontsize=6),
    column_names_rot=45,border=TRUE)
  save_ht(ht,name,w,h)
}
auc <- x$auc[,genes,drop=FALSE]
auc <- sweep(auc,2,orient[genes],'*')
am <- as.data.table(x$auc_meta)
rs <- sapply(role_levels,function(r) rowMeans(auc[,memb[role==r,gene],drop=FALSE],na.rm=TRUE))
dom <- role_levels[max.col(rs,ties.method='first')]
domval <- rs[cbind(seq_len(nrow(rs)),match(dom,role_levels))]
ord <- order(factor(dom,levels=role_levels),am$drug,-domval,am$cell_name)
ca <- data.frame(Dominant_program=factor(role_labels[dom[ord]],levels=unname(role_labels)),Drug=am$drug[ord])
gene_sample_ht(t(auc[ord,,drop=FALSE]),ca,'Tahoe molecular responses: genes × cell-line/drug profiles',
 '28b_04_tahoe_gene_sample',FALSE,7.4,5.8)

cm <- x$study_effect[genes,,drop=FALSE]
cm <- sweep(cm,1,orient[genes],'*')
crs <- sapply(role_levels,function(r) colMeans(cm[memb[role==r,gene],,drop=FALSE],na.rm=TRUE))
cdom <- role_levels[max.col(crs,ties.method='first')]
cval <- crs[cbind(seq_len(nrow(crs)),match(cdom,role_levels))]
cord <- order(factor(cdom,levels=role_levels),-cval)
cca <- data.frame(Dominant_program=factor(role_labels[cdom[cord]],levels=unname(role_labels)))
gene_sample_ht(cm[,cord,drop=FALSE],cca,'CytoSig FGF2 responses: genes × studies',
 '28b_05_cytosig_gene_sample',TRUE,6.3,5.8)

sil <- fread(file.path(out,'28b_silhouette_by_k.tsv'))
selected_k <- uniqueN(memb$program)
p <- ggplot(sil,aes(k,mean_silhouette))+geom_line(linewidth=.8,colour='#2D5F8B')+
  geom_point(aes(shape=admissible),size=2.8,colour='#2D5F8B')+
  geom_point(data=sil[k==selected_k],size=4,shape=21,fill='#F2B134')+
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1),name='Minimum 3 genes')+
  scale_x_continuous(breaks=sil$k)+labs(title='Resolution diagnostic',x='Number of clusters',y='Mean silhouette')+
  theme_grant+theme(legend.position='none')
save_plot(p,'28b_06_silhouette',4.2,4.0)

fp <- as.matrix(x$fingerprint); fp <- fp[role_levels,,drop=FALSE]
rownames(fp) <- unname(role_labels[rownames(fp)])
fr <- as.numeric(quantile(abs(fp),.98,na.rm=TRUE)); fr <- max(fr,.1)
ht <- draw_heatmap(fp,scale_rows=FALSE,value_range=c(-fr,0,fr),value_colors=c('#2166AC','#FFFFFF','#B2182B'),
 hm_name='Signed effect',show_row_names=TRUE,show_column_names=TRUE,cluster_rows=FALSE,cluster_columns=FALSE,
 title='Program perturbation fingerprints',row_names_gp=grid::gpar(fontsize=9,fontface='bold'),
 column_names_gp=grid::gpar(fontsize=8,fontface='bold'),column_names_rot=30,border=TRUE)
save_ht(ht,'28b_07_program_fingerprint',5.2,3.4)

parts <- list(); map <- c(ER='Fulvestrant',FGF='Erdafitinib',Root='palbociclib')
for(r in role_levels) {
  ii <- which(am$drug==map[r] & is.finite(am$prism_sensitivity))
  parts[[r]] <- data.table(role=unname(role_labels[r]),drug=map[r],
    molecular_score=rowMeans(auc[ii,memb[role==r,gene],drop=FALSE],na.rm=TRUE),
    prism_sensitivity=am$prism_sensitivity[ii],cell_line=am$cell_name[ii])
}
scatter <- rbindlist(parts); fwrite(scatter,file.path(out,'28b_prism_program_associations.tsv'),sep='\t')
p <- ggplot(scatter,aes(molecular_score,prism_sensitivity))+
  geom_point(size=2,alpha=.8,colour='#2D5F8B')+
  geom_smooth(method='lm',se=TRUE,linewidth=.8,colour='#B2182B',fill='#E7B7AD')+
  ggpubr::stat_cor(method='pearson',label.x.npc='left',label.y.npc='top',size=3.3)+
  facet_wrap(~role,scales='free')+labs(title='Molecular program displacement predicts PRISM response',
    x='Signed Tahoe program displacement',y='PRISM sensitivity')+theme_grant
save_plot(p,'28b_08_prism_association',8.6,3.3)

nodes <- data.table(label=c('ER branch','FGF receptor branch','Shared proliferation root','Viable-cell output'),
 x=c(0,2,1,1),y=c(2,2,1,0))
edges <- data.table(x=c(0,2,1),y=c(2,2,1),xend=c(1,1,1),yend=c(1.12,1.12,.18))
p <- ggplot()+geom_segment(data=edges,aes(x,y,xend=xend,yend=yend),linewidth=1.2,colour='#4D4D4D',
 arrow=arrow(length=unit(3,'mm'),type='closed'))+
  geom_point(data=nodes,aes(x,y,fill=label),shape=21,size=12,stroke=1.2,colour='white')+
  geom_label(data=nodes,aes(x,y,label=label),fontface='bold',size=3.4,linewidth=0,fill='white',alpha=.9)+
  annotate('text',x=0,y=2.45,label='Fulvestrant ↓',fontface='bold',size=4,colour=role_cols['ER'])+
  annotate('text',x=2,y=2.45,label='FGF2 ↑   |   FGFRi ↓',fontface='bold',size=4,colour=role_cols['FGF'])+
  annotate('text',x=1.58,y=.62,label='Palbociclib ↓',fontface='bold',size=4,colour=role_cols['Root'])+
  scale_fill_manual(values=c('ER branch'=role_cols['ER'],'FGF receptor branch'=role_cols['FGF'],
    'Shared proliferation root'=role_cols['Root'],'Viable-cell output'='#555555'))+
  coord_cartesian(xlim=c(-.6,2.6),ylim=c(-.2,2.65),clip='off')+
  labs(title='Tested Y-shaped bypass hierarchy')+theme_void(base_size=12)+
  theme(plot.title=element_text(face='bold',hjust=.5,size=13),legend.position='none',plot.margin=margin(6,10,6,10))
save_plot(p,'28b_09_y_hierarchy',4.8,4.2)
