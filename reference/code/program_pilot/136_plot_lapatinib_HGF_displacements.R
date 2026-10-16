#!/usr/bin/env Rscript
# Paired raw displacements, preserving the program order in the v136 heatmaps.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
out <- file.path(root, 'figures/single_signatures_beside_combination_benefit')
data_root <- file.path(root, 'data/single_signatures_beside_combination_benefit')
input <- '/work/results/program_pilot/134_HGF_g_diagnostics/linear/data/AU565_program_level_opposition.tsv'
displacements <- fread(input)[source_profile == 'Primary_study_mean',
  .(program, lapatinib_displacement, HGF_displacement)]
ordering <- fread(file.path(data_root, 'program_order_from_single_drug_correlations.tsv'))
correlations <- fread(file.path(data_root, 'single_drug_sensitivity_correlations.tsv'))
stopifnot(nrow(displacements)==40L, !anyDuplicated(displacements$program),
          identical(sort(displacements$program), sort(ordering$program)))
displacements <- displacements[match(ordering$program, program)]
displacements[, lapatinib_sensitivity_r := correlations[['Lapatinib 5 µM']][match(program, correlations$program)]]
displacements[, strong_lapatinib_association := abs(lapatinib_sensitivity_r) >= .5]
displacements[, program_label := paste0(program, ifelse(strong_lapatinib_association, ' *', ''))]
displacements[, panel := factor(rep(c('First 20 programs', 'Next 20 programs'), each=20),
  levels=c('First 20 programs', 'Next 20 programs'))]
stopifnot(all(is.finite(as.matrix(displacements[, .(lapatinib_displacement, HGF_displacement)]))),
          !anyNA(displacements$lapatinib_sensitivity_r))
# Verify that the displayed HGF vector is the current, unchanged model input.
current_hgf <- fread(file.path(root, 'data/time_matched_HGF_program_inputs.tsv'))[
  CellLine=='AU565' & source_profile=='Primary_study_mean']
stopifnot(nrow(current_hgf)==1L,
  max(abs(as.numeric(current_hgf[,displacements$program,with=FALSE])-displacements$HGF_displacement))<1e-12)
long <- melt(displacements, id.vars=c('program','program_label','panel'),
  measure.vars=c('lapatinib_displacement','HGF_displacement'),
  variable.name='intervention', value.name='displacement')
long[, intervention := factor(intervention,
  levels=c('lapatinib_displacement','HGF_displacement'), labels=c('Lapatinib 5 µM','HGF'))]
long[, program_label := factor(program_label, levels=rev(displacements$program_label))]
p <- ggplot(long, aes(displacement, program_label, fill=intervention)) +
  geom_vline(xintercept=0, colour='#777777', linewidth=.3) +
  geom_col(position=position_dodge(width=.8), width=.72, orientation='y') +
  facet_wrap(~panel, ncol=2, scales='free_y') +
  scale_fill_manual(values=c('Lapatinib 5 µM'='#0072B2','HGF'='#D55E00')) +
  theme_classic(base_size=10) +
  theme(legend.position='top', legend.title=element_blank(),
    strip.background=element_blank(), strip.text=element_blank(),
    panel.spacing=grid::unit(12,'mm'), plot.caption=element_text(size=8,hjust=0)) +
  labs(x='Program displacement (original score units)', y=NULL,
    title='AU565 | lapatinib and HGF molecular displacements',
    caption=paste('Identical program order to the v136 heatmaps; common horizontal scale.',
      '* Strong lapatinib sensitivity association (|r| ≥ 0.5), irrespective of HGF direction.',
      'Lapatinib: predicted displacement at 5 µM. HGF: transported study-balanced profile.',
      sep='\n'))
name <- 'AU565_lapatinib_5uM_vs_HGF_raw_displacements'
ggsave(file.path(out, paste0(name,'.pdf')), p, width=8.5, height=6.1, device=cairo_pdf)
ggsave(file.path(out, paste0(name,'.png')), p, width=8.5, height=6.1, dpi=170)
fwrite(displacements, file.path(data_root,paste0(name,'.tsv')), sep='\t')
cat('Saved:', file.path(out,paste0(name,'.pdf')), '\n')
