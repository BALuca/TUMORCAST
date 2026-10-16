#!/usr/bin/env Rscript
# Compare the two completed v98 arms; no additional fitting or model selection.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/98_drug_mixtures_shared_quadratic_decoder'
for (folder in c('data','figures')) dir.create(file.path(root,folder),showWarnings=FALSE)
arms <- c('linear','quadratic')
metrics <- combos <- query <- mixture <- checks <- list()
for (arm in arms) {
  branch <- file.path(root,arm)
  measured <- fread(file.path(branch,'data/measured_decoder_correlations.tsv'))[evaluation=='held_cell']
  baseline <- fread(file.path(branch,'data/baseline_bridge_correlations.tsv'))[evaluation=='held_cell']
  metrics[[arm]] <- merge(measured[,.(drug,dose,measured_cells=cells,measured_r=r)],
    baseline[,.(drug,dose,baseline_cells=cells,baseline_r=r,observed_median_pct,predicted_median_pct)],
    by=c('drug','dose'))[,readout:=arm]
  combos[[arm]] <- fread(file.path(branch,'data/combination_correlations.tsv'))[evaluation=='held_cell'][,readout:=arm]
  query[[arm]] <- fread(file.path(branch,'data/au565_hgf_5uM_predictions_frozen.tsv'))[,readout:=arm]
  mixture[[arm]] <- fread(file.path(branch,'data/molecular_signature_mixtures.tsv'))[fit=='final'][,readout:=arm]
  audit <- fread(file.path(branch,'data/fit_audit.tsv'))
  gradient <- fread(file.path(branch,'data/analytic_gradient_check.tsv'))
  stopifnot(nrow(audit)==63L,all(audit$convergence==0L),all(audit$functional_leakage==0L),
    all(audit$RNA_leakage==0L),max(audit$zero_displacement_effect)<1e-12,
    max(gradient$relative_difference)<1e-5)
  checks[[arm]] <- data.table(readout=arm,fits=nrow(audit),solver_failures=sum(audit$convergence!=0L),
    zero_effect=max(audit$zero_displacement_effect),gradient_difference=max(gradient$relative_difference))
}
metrics <- rbindlist(metrics)
combos <- rbindlist(combos)
query <- rbindlist(query)
mixture <- rbindlist(mixture)
fwrite(metrics,file.path(root,'data/decoder_arm_comparison.tsv'),sep='\t')
fwrite(combos,file.path(root,'data/combination_arm_comparison.tsv'),sep='\t')
fwrite(query,file.path(root,'data/hgf_arm_comparison.tsv'),sep='\t')
fwrite(rbindlist(checks),file.path(root,'data/arm_validation.tsv'),sep='\t')

# Lapatinib scatter panels separate measured and baseline-predicted displacements.
points <- list()
for (arm in arms) {
  for (input in c('measured','baseline')) {
    p <- fread(file.path(root,arm,'data',paste0(input,'_predictions_frozen.tsv')))[
      evaluation=='held_cell' & drug=='Lapatinib ditosylate']
    p[,`:=`(readout=arm,input=if (input=='measured') 'Measured displacements (30 cells)' else 'Baseline-predicted displacements (441 cells)')]
    points[[length(points)+1L]] <- p
  }
}
points <- rbindlist(points)
stats <- points[,.(label=paste0('Held-out r = ',sprintf('%.2f',cor(observed_log2,predicted_log2)),'\nn = ',.N)),by=.(input,readout)]
p <- ggplot(points,aes(observed_log2,predicted_log2))+
  geom_point(color='#0072B2',alpha=.7,size=1.4)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey25',linewidth=.5)+
  geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.15,size=3)+
  facet_wrap(~input+readout,ncol=2,scales='free')+theme_bw(base_size=10)+theme(aspect.ratio=1)+
  labs(title='v98 | lapatinib 5 µM: linear versus quadratic shared readout',
    x='Measured absolute log2 PRISM abundance',y='Predicted absolute log2 PRISM abundance',
    caption='Pearson across held-out cells at the same drug/dose. 0 = 100% viability; -1 = 50%; -2 = 25%.')
ggsave(file.path(root,'figures/lapatinib_decoder_arm_comparison.pdf'),p,width=8,height=8.5,device=cairo_pdf)
ggsave(file.path(root,'figures/lapatinib_decoder_arm_comparison.png'),p,width=8,height=8.5,dpi=150)

# Different reference controls are deliberately displayed in separate panels.
h <- query[hgf_context=='Equal-context mean']
h[,conditional_rescue_percentage_points:=lapatinib_in_hgf_pct_hgf_control-lapatinib_pct_vehicle]
fwrite(h,file.path(root,'data/au565_hgf_equal_context_summary.tsv'),sep='\t')
vehicle <- melt(h,id.vars='readout',measure.vars=c('lapatinib_pct_vehicle','hgf_alone_pct_vehicle','combination_pct_vehicle'),
                variable.name='condition',value.name='viability')
vehicle[,reference:='Relative to untreated vehicle']
vehicle[,condition:=factor(condition,levels=c('lapatinib_pct_vehicle','hgf_alone_pct_vehicle','combination_pct_vehicle'),
                          labels=c('Lapatinib','HGF','Lapatinib + HGF'))]
conditional <- melt(h,id.vars='readout',measure.vars=c('lapatinib_pct_vehicle','lapatinib_in_hgf_pct_hgf_control'),
                    variable.name='condition',value.name='viability')
conditional[,reference:='Lapatinib effect relative to its matched control']
conditional[,condition:=factor(condition,levels=c('lapatinib_pct_vehicle','lapatinib_in_hgf_pct_hgf_control'),
                              labels=c('No HGF: vehicle control','With HGF: HGF-only control'))]
p <- ggplot(rbind(vehicle,conditional),aes(condition,viability,fill=readout))+
  geom_col(position=position_dodge(width=.8),width=.7)+
  geom_text(aes(label=sprintf('%.2f',viability)),position=position_dodge(width=.8),vjust=-.4,size=3)+
  facet_wrap(~reference,ncol=1,scales='free_x')+scale_y_continuous(limits=c(0,110),breaks=c(0,25,50,75,100))+
  scale_fill_manual(values=c(linear='#0072B2',quadratic='#D55E00'),name='Shared readout')+
  theme_bw(base_size=10)+theme(legend.position='bottom',axis.title.x=element_blank(),panel.grid.major.x=element_blank())+
  labs(title='v98 | AU565 lapatinib 5 µM + transported HGF query',y='Predicted viability (%)',
    caption='HGF mixture uses the shared footprint map; drug-specific correction = 0. No HGF/Wilson labels fitted.\nHGF receiver SD is inferred. Wilson reference was at 1 µM, so this is not a dose-matched validation.')
ggsave(file.path(root,'figures/au565_hgf_decoder_arm_comparison.pdf'),p,width=8,height=7,device=cairo_pdf)
ggsave(file.path(root,'figures/au565_hgf_decoder_arm_comparison.png'),p,width=8,height=7,dpi=150)

# Signed mixture maps are displayed in the observed-correlation clustering order.
model <- qread(file.path(root,'linear/models/final.qs'))
co <- model$template$drug_tree$order
drugs <- model$selected$drug
short <- c('Azacytidine','Bendamustine','Bortezomib','Busulfan','Cytarabine',
           'Lapatinib','Topotecan','Vandetanib','Vismodegib','Crizotinib')
cairo_pdf(file.path(root,'figures/learned_signed_mixture_maps.pdf'),width=8,height=3.5,onefile=TRUE)
for (arm in arms) {
  z <- mixture[readout==arm][match(drugs,drug)]
  M <- t(as.matrix(z[,.(signature_1,signature_2)]))[ ,co,drop=FALSE]
  dimnames(M) <- list(c('Signature 1','Signature 2'),short[co])
  hm <- draw_heatmap(M,scale_rows=FALSE,value_range=c(-1,0,1),
    value_colors=c('#2166AC','white','#B2182B'),hm_name='Signed mixture',
    show_row_names=TRUE,show_column_names=TRUE,column_names_rot=45,
    row_names_gp=grid::gpar(fontsize=9),column_names_gp=grid::gpar(fontsize=8),
    cluster_rows=FALSE,cluster_columns=FALSE,title=paste('v98',arm,'| learned drug mixtures'),
    height=grid::unit(15,'mm'))
  ComplexHeatmap::draw(hm,heatmap_legend_side='bottom')
}
dev.off()
cat('Saved comparison plots, HGF query and signed mixture heatmaps.\n')
print(metrics[drug=='Lapatinib ditosylate'])
print(h[,.(readout,lapatinib_pct_vehicle,combination_pct_vehicle,
           lapatinib_in_hgf_pct_hgf_control,conditional_rescue_percentage_points)])
