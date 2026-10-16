#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
iteration <- sub('^--iteration=', '', args[grepl('^--iteration=',args)])
if(!length(iteration)) iteration <- '04'
iteration <- sprintf('%02d',as.integer(iteration))
root <- '/work/tahoe_prism_v2'
out <- file.path(root,'results/28_fgf2_fulvestrant_bypass',paste0('iteration_',iteration))
fig <- file.path(out,'figures'); dir.create(fig,recursive=TRUE,showWarnings=FALSE)

# Values are the maximum fold-change annotations printed over the two
# fulvestrant groups in Shee et al. 2018 Fig. S3H. The downloadable raw SRB
# replicates were not available, so no inferential statistic is reconstructed.
d <- data.table(
 treatment=c('Fulvestrant + FGF2','Fulvestrant + FGF2 + palbociclib'),
 palbociclib=c('Absent','Present'),max_viable_cell_fold_change=c(1.64,1.36),
 source_panel='Shee et al. 2018, Fig. S3H',assay='SRB viable-cell abundance, 5 days',
 provenance='Transcribed from published maximum fold-change labels; not raw measurements')
d[,excess_rescue_over_baseline:=max_viable_cell_fold_change-1]
attenuation <- 1-d[palbociclib=='Present',excess_rescue_over_baseline]/d[palbociclib=='Absent',excess_rescue_over_baseline]
fwrite(d,file.path(out,'28d_shee_functional_rescue.tsv'),sep='\t')

p <- ggplot(d,aes(palbociclib,max_viable_cell_fold_change,group=1))+
 geom_hline(yintercept=1,linetype=2,linewidth=.6,colour='#777777')+
 geom_line(linewidth=1.1,colour='#555555')+
 geom_point(aes(fill=palbociclib),shape=21,size=5,stroke=1,colour='white')+
 geom_text(aes(label=sprintf('%.2f×',max_viable_cell_fold_change)),vjust=-1.25,fontface='bold',size=4.2)+
 annotate('label',x=1.5,y=1.12,label=sprintf('Palbociclib removes %.0f%%\nof excess FGF2 rescue',100*attenuation),
   size=3.6,fontface='bold',fill='#F3F3F3',linewidth=.3)+
 scale_fill_manual(values=c(Absent='#D97818',Present='#167C80'))+
 scale_y_continuous(limits=c(.95,1.78),breaks=seq(1,1.8,.2),expand=expansion(mult=c(.01,.03)))+
 labs(title='Functional rescue in parental T47D',subtitle='FGF2 rescue during fulvestrant treatment',
  x='Palbociclib (1 µM)',y='Maximum viable-cell fold change\nwith FGF2')+
 theme_classic(base_size=12)+theme(plot.title=element_text(face='bold',size=13),
  plot.subtitle=element_text(size=10),axis.title=element_text(face='bold'),legend.position='none',aspect.ratio=1)
ggsave(file.path(fig,'28d_01_shee_functional_rescue.pdf'),p,width=4.3,height=4.1,device=cairo_pdf)
ggsave(file.path(fig,'28d_01_shee_functional_rescue.png'),p,width=4.3,height=4.1,dpi=220,bg='white')
