#!/usr/bin/env Rscript
# Plot existing frozen v107 queries beside published Wilson measurements; no model fitting.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(qs)})
root <- '/work/results/program_pilot/107_within_pair_departure_interaction/linear'
query_path <- file.path(root,'data/au565_hgf_5uM_zero_reference_query.tsv')
wilson_path <- '/work/almanac_combination_training_repair_2026-10-03/results/wilson_observed_source.tsv'
query <- fread(query_path)[hgf_context=='Equal-context mean' & assumed_HGF_route==1]
wilson <- fread(wilson_path)[series=='Wilson et al.' & condition!='Lapatinib + HGF + crizotinib']
stopifnot(nrow(query)==1L,query$assumed_HGF_route==1,nrow(wilson)==2L,
          uniqueN(query$lapatinib_pct_zero_control)==1L,
          all(is.finite(wilson$abundance)),all(is.finite(wilson$sd)))
conditions <- c('Lapatinib','Lapatinib + HGF')
stopifnot(setequal(wilson$condition,conditions))
# PRISM-trained drug head already predicts abundance relative to the assay control.
# Its zero-displacement prediction is not a measured vehicle control and must not
# be used to renormalize that absolute drug prediction. HGF remains an incremental
# zero-referenced external response; add its effect and g to the fitted drug value.
model <- qread(file.path(root,'models/final.qs'))
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
receptors <- as.data.table(archive$prep$receptors)
lapmap <- model$rna$lap_receptors
receptor_values <- as.matrix(receptors[depmap_id=='ACH-000248',.(ERBB2,MET)])
lap_D <- as.numeric(sweep(sweep(sweep(sweep(receptor_values,2,lapmap$xmean),2,lapmap$xsd,'/') %*%
                              lapmap$B,2,lapmap$sigma,'*'),2,lapmap$mu,'+'))
contexts <- fread(file.path(root,'data/molecular_contexts.tsv'))
g <- which(contexts$drug=='Lapatinib ditosylate' & contexts$dose==5)
ctx <- model$context
lap_score <- sum(((lap_D-ctx$mu[g,])/ctx$sigma[g,])*ctx$W[,ctx$route[g]])
lap_level <- as.numeric(ctx$Q[g,,drop=FALSE] %*% model$par$level)
lap_amplitude <- exp(as.numeric(ctx$Q[g,,drop=FALSE] %*% model$par$log_amplitude))
lap_abundance <- 100*2^(-lap_level-lap_amplitude*lap_score)
hgf_increment <- query$hgf_alone_pct_zero_control/100 * 2^query$interaction_log2
combined_abundance <- lap_abundance*hgf_increment
prediction_panel <- 'Predicted\nLapatinib 5 µM'
observed_panel <- 'Wilson\nLapatinib 1 µM'
predictions <- rbind(
 data.table(panel=prediction_panel,condition='Lapatinib',x=1,
            series='Predicted',abundance=lap_abundance,sd=NA_real_),
 query[,.(panel=prediction_panel,condition='Lapatinib + HGF',
           x=2,series='Predicted',
           abundance=combined_abundance,sd=NA_real_)])
observed <- wilson[,.(panel=observed_panel,condition,x=match(condition,conditions),
                      series='Wilson',abundance,sd)]
plot_data <- rbind(predictions,observed)
plot_data[,panel:=factor(panel,levels=c(prediction_panel,observed_panel))]
plot_data[,label:=sprintf('%.1f',abundance)]
plot_data[,label_y:=abundance+ifelse(is.na(sd),2,sd+3)]
write_data <- copy(plot_data)
write_data[,`:=`(cell='AU565',reference='Fitted absolute drug abundance with zero-referenced HGF increment and g; no drug-head renormalization',
 model_lapatinib_uM=ifelse(panel==prediction_panel,5,NA_real_),
 observed_lapatinib_uM=ifelse(panel==observed_panel,1,NA_real_))]
fwrite(write_data,file.path(root,'data/predictions_beside_wilson_barplot.tsv'),sep='\t')
plot <- ggplot(plot_data,aes(x,abundance,fill=series))+
 geom_col(width=.55)+
 geom_errorbar(data=plot_data[!is.na(sd)],aes(ymin=abundance-sd,ymax=abundance+sd),
               width=.18,linewidth=.35)+
 geom_text(aes(y=label_y,label=label),size=2.4)+
 facet_wrap(~panel,nrow=1,scales='free_y')+
 scale_fill_manual(values=c('Predicted'='#D55E00','Wilson'='#0072B2'))+
 scale_x_continuous(breaks=1:2,labels=c('Lapatinib','Lapatinib\n+ HGF'),
                    limits=c(.5,2.5),expand=expansion(mult=0))+
 scale_y_continuous(breaks=scales::breaks_pretty(n=3),expand=expansion(mult=c(0,.07)))+
 theme_bw(base_size=7.5)+
 theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
       strip.background=element_rect(fill='grey95'),strip.text=element_text(size=7.5,margin=margin(2,2,2,2)),
       legend.position='none',panel.spacing=grid::unit(5,'pt'),
       axis.text.x=element_text(size=7),plot.caption=element_text(size=6.5,hjust=0),
       axis.title.y=element_text(margin=margin(r=3)),plot.margin=margin(3,4,3,3))+
 labs(x=NULL,y='Abundance (% of control)',
      caption='AU565 · HGF pattern 1; Wilson mean ± SD.')
ggsave(file.path(root,'figures/predictions_beside_wilson.pdf'),plot,width=3.4,height=2.4,device=cairo_pdf)
ggsave(file.path(root,'figures/predictions_beside_wilson.png'),plot,width=3.4,height=2.4,dpi=400)
print(plot_data[,.(panel,condition,series,abundance,sd)])
