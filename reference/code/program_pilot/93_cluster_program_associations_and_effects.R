#!/usr/bin/env Rscript
# Cluster observed correlations and apply the same order to every model comparison.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2)})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
figures <- file.path(root,'figures/program_and_route_diagnostics/clustered')
dir.create(figures,recursive=TRUE,showWarnings=FALSE)
z <- fread(file.path(root,'data/program_association_vs_model_effect.tsv'))
selection <- fread(file.path(root,'data/selected_doses.tsv'))
programs <- sprintf('P%02d',1:40)
drugs <- selection$drug
short <- c('Azacytidine','Bendamustine','Bortezomib','Busulfan','Cytarabine','Lapatinib',
  'Topotecan','Vandetanib','Vismodegib','Crizotinib')
labels <- paste0(short,' | ',selection$dose,' µM')
matrices <- list()
for(field in c('observed_r','model_marginal_r','mean_sensitivity_effect_per_program_SD')) {
  wide <- dcast(z,program~drug,value.var=field)
  M <- as.matrix(wide[match(programs,program),drugs,with=FALSE])
  dimnames(M) <- list(programs,labels)
  stopifnot(all(is.finite(M)),nrow(M)==40L,ncol(M)==10L)
  matrices[[field]] <- M
}
# Euclidean distances retain correlation sign and strength; Ward groups patterns.
observed <- matrices$observed_r
row_tree <- hclust(dist(observed),method='ward.D2')
column_tree <- hclust(dist(t(observed)),method='ward.D2')
ro <- row_tree$order
co <- column_tree$order
row_groups <- cutree(row_tree,5)
col_groups <- cutree(column_tree,3)
row_split <- factor(paste0('P-group ',row_groups[ro]),levels=unique(paste0('P-group ',row_groups[ro])))
col_split <- factor(paste0('Drug-group ',col_groups[co]),levels=unique(paste0('Drug-group ',col_groups[co])))
cluster_names <- paste0('Pattern ',row_groups)
ncell <- z[,.(cells=unique(cells)),by=drug]
annotation <- data.frame(Dose=paste0(selection$dose[co],' µM'),
  Cells=as.character(ncell$cells[match(drugs[co],ncell$drug)]),row.names=labels[co])
qsave(list(matrices=matrices,row_tree=row_tree,column_tree=column_tree,row_order=ro,column_order=co,
  row_groups=row_groups,column_groups=col_groups),file.path(root,'data/correlation_cluster_inputs.qs'))
fwrite(data.table(program=programs[ro],display_row=1:40,cluster=row_groups[ro]),
  file.path(root,'data/correlation_program_order.tsv'),sep='\t')
fwrite(data.table(drug=drugs[co],dose=selection$dose[co],display_column=1:10,cluster=col_groups[co]),
  file.path(root,'data/correlation_drug_order.tsv'),sep='\t')
headings <- c(observed_r='Observed program association with sensitivity',
  model_marginal_r='Model-implied program association with sensitivity',
  mean_sensitivity_effect_per_program_SD='Effective model sensitivity beta per drug-specific program SD')
book <- file.path(figures,'clustered_program_associations_and_model_effects.pdf')
cairo_pdf(book,width=9,height=10,onefile=TRUE)
book_device <- dev.cur()
for(field in names(matrices)) {
  mat <- matrices[[field]][ro,co,drop=FALSE]
  stopifnot(identical(colnames(mat),rownames(annotation)))
  range <- if(field=='mean_sensitivity_effect_per_program_SD')max(abs(mat)) else 1
  render <- function() {
    h <- draw_heatmap(mat,col_anno_df=annotation,col_anno_cols=names(annotation),
      palette_overrides=list(Dose=c('0.05 µM'='#DCEAF2','0.5 µM'='#6BAED6','5 µM'='#08519C'),
        Cells=c('17'='#CECECE','29'='#969696','30'='#555555')),
      scale_rows=FALSE,value_range=c(-range,0,range),value_colors=c('#2166AC','white','#B2182B'),
      hm_name=if(field=='mean_sensitivity_effect_per_program_SD')'Conditional effect\n(log2 / program SD)' else 'Pearson r',
      cluster_rows=FALSE,cluster_columns=FALSE,row_split=row_split,column_split=col_split,
      cluster_row_slices=FALSE,cluster_column_slices=FALSE,
      row_gap=grid::unit(2,'mm'),column_gap=grid::unit(2,'mm'),
      show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=8),
      column_names_gp=grid::gpar(fontsize=8),column_names_rot=45,column_labels=short[co],
      row_title_gp=grid::gpar(fontsize=8),column_title_gp=grid::gpar(fontsize=8),
      title=paste0(headings[field],'\nSame ordering learned from observed correlations; positive = greater sensitivity'))
    drawn <- ComplexHeatmap::draw(h,heatmap_legend_side='bottom',annotation_legend_side='right')
    stopifnot(identical(as.integer(unlist(ComplexHeatmap::row_order(drawn))),1:40),
      identical(as.integer(unlist(ComplexHeatmap::column_order(drawn))),1:10))
  }
  dev.set(book_device)
  render()
  png(file.path(figures,paste0('clustered_',field,'.png')),width=9,height=10,units='in',res=150,type='cairo')
  render()
  dev.off()
}
dev.set(book_device)
dev.off()
# Reorder all barplots by exactly the observed matrix clusters.
z[,panel:=factor(labels[match(drug,drugs)],levels=labels[co])]
z[,program:=factor(program,levels=rev(programs[ro]))]
z[,pattern:=factor(cluster_names[match(as.character(program),programs)])]
bar_theme <- theme_bw(base_size=9)+theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
  strip.text=element_text(size=8),axis.text.y=element_text(size=6),legend.position='bottom')
paired <- melt(z,id.vars=c('drug','dose','program','panel'),measure.vars=c('observed_r','model_marginal_r'),
  variable.name='source',value.name='r')
p <- ggplot(paired,aes(r,program,fill=source))+geom_vline(xintercept=0,color='grey55',linewidth=.3)+
  geom_col(position=position_dodge(width=.85),width=.78)+facet_wrap(~panel,ncol=5)+
  scale_x_continuous(limits=c(-1,1),breaks=c(-1,0,1))+
  scale_fill_manual(values=c(observed_r='#39799D',model_marginal_r='#C27535'),
    labels=c(observed_r='Observed sensitivity',model_marginal_r='Fitted model sensitivity'),name=NULL)+bar_theme+
  labs(title='Clustered programs and drugs: observed versus fitted associations',
    x='Pearson r with sensitivity (-log2 abundance)',y=NULL,
    caption='Rows and drug panels have the same clustering order as the observed correlation heatmap.\nAll-data fitted model associations are descriptive; no model refitting or validation claim.')
ggsave(file.path(figures,'clustered_observed_vs_model_correlations_barplots.pdf'),p,width=15,height=13,device=cairo_pdf)
ggsave(file.path(figures,'clustered_observed_vs_model_correlations_barplots.png'),p,width=15,height=13,dpi=110)
beta <- ggplot(z,aes(mean_sensitivity_effect_per_program_SD,program,fill=pattern))+
  geom_vline(xintercept=0,color='grey55',linewidth=.3)+geom_col(width=.8)+facet_wrap(~panel,ncol=5,scales='free_x')+
  bar_theme+labs(title='Effective model betas in observed-correlation cluster order',
    x='Mean conditional sensitivity effect per program SD',y=NULL,fill='Observed association group',
    caption='Same program/drug order as the observed heatmap; these derivatives are conditional effects, not marginal correlations.')
ggsave(file.path(figures,'clustered_effective_beta_barplots.pdf'),beta,width=15,height=13,device=cairo_pdf)
ggsave(file.path(figures,'clustered_effective_beta_barplots.png'),beta,width=15,height=13,dpi=110)
# Program loadings use the same program order, exposing common versus opposing directions.
model <- qread(file.path(root,'models/final.qs'))
loadings <- data.table(program=rep(programs,2),route=factor(rep(1:2,each=40)),weight=as.numeric(model$par$W))
loadings[,program:=factor(program,levels=rev(programs[ro]))]
lp <- ggplot(loadings,aes(weight,program,fill=route))+geom_vline(xintercept=0,color='grey55',linewidth=.3)+
  geom_col(width=.8)+facet_wrap(~route,ncol=2)+scale_fill_manual(values=c('#39799D','#A96E9E'),guide='none')+
  theme_bw(base_size=10)+theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=8))+
  labs(title='Learned route loadings in observed-correlation cluster order',x='Preactivation weight W',y=NULL)
ggsave(file.path(figures,'clustered_route_loading_barplots.pdf'),lp,width=8,height=9,device=cairo_pdf)
ggsave(file.path(figures,'clustered_route_loading_barplots.png'),lp,width=8,height=9,dpi=130)
writeLines(c('Clusters are descriptive, derived from the signed observed 40-program x 10-drug Pearson correlation matrix.',
  'Programs and drugs use Euclidean distances and Ward.D2, preserving signed correlation strength.',
  'Five program groups and three drug groups are display cuts, not inferred biological route counts.',
  'All heatmap pages and barplots share the same ordering. Correlations are not row-standardized.',
  'Sensitivity = -log2 abundance; positive r/effect means greater sensitivity.',
  'Observed and fitted marginal correlations are comparable; conditional beta units differ.',
  'Model directions are all-data fitted descriptions; no model refitting or held-out selection.'),
  file.path(figures,'CLUSTER_DEFINITIONS.txt'))
cat('Saved three clustered heatmap pages and consistently reordered barplots.\n')
