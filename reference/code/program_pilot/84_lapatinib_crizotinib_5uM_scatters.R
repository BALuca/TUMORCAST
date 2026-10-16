.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2)})
setDTthreads(1L)
root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/figures/observed_program_viability/lapatinib_crizotinib_5uM'
a <- qread(file.path(root,'heatmap_global_displacement_patterns_inputs.qs'))
metadata <- fread('/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/figures/observed_program_viability/lapatinib_crizotinib_5uM/heatmap_clustered_sensitivity_column_annotations.tsv')[view=='global_displacement_patterns']
metadata <- metadata[match(rownames(a$raw_displacements),profile_id)]
stopifnot(identical(metadata$profile_id,rownames(a$raw_displacements)),nrow(metadata)==60L,all(metadata$dose==5))
programs <- sprintf('P%02d',1:40)
wide <- cbind(metadata[,.(depmap_id,cell_line,drug,dose,observed_log2_abundance)],as.data.table(a$raw_displacements))
z <- melt(wide,id.vars=c('depmap_id','cell_line','drug','dose','observed_log2_abundance'),measure.vars=programs,
  variable.name='program',value.name='displacement',variable.factor=FALSE)
z[, Drug:=ifelse(drug=='crizotinib','Crizotinib','Lapatinib')]
colors <- c(Crizotinib='#59A14F',Lapatinib='#B07AA1')
metrics <- z[,.(n=.N,r=cor(displacement,observed_log2_abundance)),by=.(program,Drug)]
fwrite(metrics,file.path(root,'scatter_correlations_5uM.tsv'),sep='\t')
labels <- dcast(metrics,program~Drug,value.var='r')
labels[,label:=sprintf('Lapatinib r = %.2f\nCrizotinib r = %.2f',Lapatinib,Crizotinib)]
cairo_pdf(file.path(root,'lapatinib_crizotinib_5uM_all_40_program_scatters.pdf'),width=12.5,height=6.5)
book <- dev.cur()
for(page in 1:4) {
  pp <- programs[(page-1)*10+1:10]
  zz <- copy(z[program %in% pp])
  zz[,program:=factor(program,levels=pp)]
  ll <- copy(labels[program %in% pp])
  ll[,program:=factor(program,levels=pp)]
  p <- ggplot(zz,aes(displacement,observed_log2_abundance))+
    geom_point(aes(color=Drug),size=1.1,alpha=.8)+
    geom_smooth(aes(color=Drug),method='lm',formula=y~x,se=FALSE,linewidth=.55)+
    geom_label(data=ll,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.03,
      size=2.5,label.size=0,fill='#FFFFFFDD',label.padding=grid::unit(.06,'lines'))+
    facet_wrap(~program,ncol=5,scales='free_x')+
    scale_color_manual(values=colors,name='Drug; both 5 µM')+
    scale_y_continuous(expand=expansion(mult=c(.04,.16)))+
    theme_bw(base_size=10)+theme(aspect.ratio=1,legend.position='bottom',panel.grid.minor=element_blank())+
    labs(title='Lapatinib and crizotinib at 5 µM: measured program displacement versus PRISM abundance',
      subtitle='30 measured cell lines per drug; separate trend lines and Pearson correlations across cells.',
      x='Measured program displacement (signed RNA score)',y='Observed PRISM log2 abundance / vehicle',
      caption='PRISM curves interpolated to Tahoe 5 µM; separate experiments, not paired wells.\n0 = 100%, -1 = 50%, -2 = 25% vehicle abundance. No model fitted; no filtering.')
  dev.set(book)
  print(p)
  png(file.path(root,paste0('scatters_P',sprintf('%02d',min((page-1)*10+1:10)),'_P',sprintf('%02d',page*10),'.png')),
    width=12.5,height=6.5,units='in',res=100,type='cairo')
  print(p)
  dev.off()
}
dev.set(book)
dev.off()
cat('Saved all 40 program scatters, both drugs at 5 µM.\n')
