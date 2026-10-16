#!/usr/bin/env Rscript
# Consolidated diagnostics from completed ALMANAC exclusions; no model changes.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/72_almanac_full_program_interactions'
setDTthreads(2L)
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
save_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
save_plot <- function(p,name,width,height) {
  ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=width,height=height,dpi=140)
}
safe_cor <- function(x,y) {
  if (length(x)<3L || sd(x)<1e-10 || sd(y)<1e-10) return(NA_real_)
  cor(x,y)
}
label_r <- function(x,y) {
  r <- safe_cor(x,y)
  if (is.na(r)) return('Pearson r undefined: constant prediction')
  sprintf('Pearson r = %.2f',r)
}
input <- qread(file.path(root,'data/prepared_inputs.qs'))
md <- input$metadata
# Read only folds recorded complete at the start, so plots can be produced while
# later exclusions are fitting. This changes reporting scope, never the model.
audit <- fread(file.path(root,'data/validation.tsv.gz'))
drug_ids <- sort(unique(c(md$drug_a,md$drug_b)))
tags <- sprintf('%02d_%s',match(audit$held_drug,drug_ids),gsub('[^A-Za-z0-9]+','_',audit$held_drug))
files <- file.path(root,'data',paste0('prediction_',tags,'.tsv.gz'))
stopifnot(all(file.exists(files)))
stopifnot(length(files)>0)
predictions <- rbindlist(lapply(files,fread))
z <- merge(predictions,md[,.(combination_id,observed_g,observed_ratio,expected_ratio,individual_log2)],
  by='combination_id',all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(z$observed_g),!anyDuplicated(z[,.(combination_id,held_drug,split)]))
z[, `:=`(observed_count=100*observed_ratio,independent_count=100*expected_ratio,
  predicted_count=100*2^(individual_log2+predicted_g))]
stopifnot(all(is.finite(z$predicted_count)))
summary <- z[,.(prediction_events=.N,conditions=uniqueN(combination_id),pairs=uniqueN(pair),
  interaction_r=safe_cor(observed_g,predicted_g),
  individual_reference_count_r=safe_cor(observed_count,independent_count),
  combination_count_r=safe_cor(observed_count,predicted_count)),by=split]
by_drug <- z[,.(conditions=.N,pairs=uniqueN(pair),cells=uniqueN(recipient),
  interaction_r=safe_cor(observed_g,predicted_g)),by=.(held_drug,split)]
by_cell <- z[split=='held',.(prediction_events=.N,conditions=uniqueN(combination_id),
  interaction_r=safe_cor(observed_g,predicted_g)),by=recipient]
by_drug_cell <- z[split=='held',.(conditions=.N,pairs=uniqueN(pair),
  interaction_r=safe_cor(observed_g,predicted_g)),by=.(held_drug,recipient)]
save_table(summary,'correlation_summary')
save_table(by_drug,'correlations_by_drug_complete')
save_table(by_cell,'correlations_by_cell')
save_table(by_drug_cell,'correlations_by_drug_cell')
save_table(z[split=='held'],'held_figure_values')
print(summary)
print(by_drug[split=='held'])

# Display actual predictions from each fitted model in both panels. Neither
# training nor held predictions are averaged into an ensemble.
training <- z[split=='training',.(combination_id,recipient,observed_g,predicted_g)]
training[,split:='training predictions']
held <- z[split=='held',.(combination_id,recipient,observed_g,predicted_g)]
held[,split:='held drug']
d <- rbind(training,held)
d[,split:=factor(split,levels=c('training predictions','held drug'))]
ann <- d[,.(label=label_r(observed_g,predicted_g)),by=split]
p <- ggplot(d,aes(observed_g,predicted_g))+geom_point(aes(color=recipient),size=.6,alpha=.25)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.6)+
  geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.3,size=3.3)+
  facet_wrap(vars(split),scales='free')+theme(aspect.ratio=1)+
  labs(title='ALMANAC: full explicit program interactions',
    subtitle=paste(length(files),'of 40 drug exclusions complete; original40 programs, bootstrap screen + elastic net'),
    x='Observed log2[combination / (single A × single B)]',y='Predicted log2 interaction',color='Cell line',
    caption='All combinations containing the test drug are excluded across cells and doses. Shared single references purged.\nPredictions from each fold shown separately; repeated conditions are dependent. Retrospective development; measured single responses supplied.')
save_plot(p,'interaction_training_vs_held',10,6)

drug_order <- sort(unique(z$held_drug))
cairo_pdf(file.path(root,'figures/all40_held_drugs.pdf'),width=10,height=7.6)
for (page in seq_len(ceiling(length(drug_order)/6))) {
  ids <- drug_order[seq((page-1)*6+1,min(page*6,length(drug_order)))]
  d <- z[split=='held' & held_drug %chin% ids]
  ann <- d[,.(label=label_r(observed_g,predicted_g)),by=held_drug]
  p <- ggplot(d,aes(observed_g,predicted_g))+geom_point(aes(color=recipient),size=.7,alpha=.35)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
    geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.2,size=2.8)+
    facet_wrap(vars(held_drug),ncol=3,scales='free')+theme(aspect.ratio=1)+
    labs(title=paste('Excluded drug interaction predictions — page',page),
      x='Observed log2 interaction',y='Predicted log2 interaction',color='Cell line')
  print(p)
}
dev.off()
ann <- z[split=='held',.(label=label_r(observed_g,predicted_g)),by=recipient]
p <- ggplot(z[split=='held'],aes(observed_g,predicted_g))+geom_point(size=.6,alpha=.2,color='#0072B2')+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
  geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.2,size=2.8)+
  facet_wrap(vars(recipient),ncol=4,scales='free')+theme(aspect.ratio=1)+
  labs(title='Held-drug predictions in every cell line',
    subtitle='Cell lines are represented in training; drug combinations containing the queried drug are excluded',
    x='Observed log2 interaction',y='Predicted log2 interaction')
save_plot(p,'held_interactions_by_cell',12,9)

stopifnot(all(audit$shared_reference_groups==0L),nrow(audit)==length(files))
funnel <- melt(audit[,.(held_drug,retained_terms,active_terms)],id.vars='held_drug',
  variable.name='stage',value.name='terms')
p <- ggplot(funnel,aes(terms,reorder(held_drug,terms),color=stage))+geom_point(size=1.5)+
  labs(title='How many of the820 interactions remain?',x='Number of terms',y=NULL,color=NULL)+
  theme(axis.text.y=element_text(size=7))
save_plot(p,'retained_interaction_counts',6,9)

co <- fread(file.path(root,'data/coefficients.tsv.gz'))[held_drug %chin% audit$held_drug]
terms <- co[,.(selected_folds=sum(coefficient!=0),median_coefficient=median(coefficient),
  positive_folds=sum(coefficient>0),negative_folds=sum(coefficient<0)),by=term]
setorder(terms,-selected_folds)
save_table(terms,'term_recurrence')
programs <- colnames(input$query_scores)
mat <- matrix(0,length(programs),length(programs),dimnames=list(programs,programs))
ab <- tstrsplit(terms$term,'__',fixed=TRUE)
for (i in seq_len(nrow(terms))) {
  v <- terms$median_coefficient[i]/if(ab[[1]][i]==ab[[2]][i]) 1 else 2
  mat[ab[[1]][i],ab[[2]][i]] <- mat[ab[[2]][i],ab[[1]][i]] <- v
}
order <- hclust(dist(mat),method='ward.D2')$order
mat <- mat[order,order,drop=FALSE]
limit <- max(abs(mat),1e-8)
hm <- draw_heatmap(mat,scale_rows=FALSE,value_range=c(-limit,0,limit),
  value_colors=c('#2166AC','white','#B2182B'),hm_name='Median\ncoefficient',
  show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=7),
  column_names_gp=grid::gpar(fontsize=7),column_names_rot=90,
  title='Explicit program interactions: median coefficient across drug exclusions')
cairo_pdf(file.path(root,'figures/program_interaction_coefficients.pdf'),width=8,height=7.5)
ComplexHeatmap::draw(hm)
dev.off()
report <- c('72: ALMANAC full820 original40 program interactions, bootstrap screened.',
  paste('Completed drug exclusions:',length(files),'of40.'),
  'Target: mean raw matched log2(combo/(singleA*singleB)); known separate molecular and functional measurements.',
  'Whole-drug exclusions across cells/doses; exact shared-single references purged in outer and inner folds.',
  '500 pair-block bootstrap draws;90% slope interval excludes0; elastic-net shrinkage; no nonlinear kernels.',
  capture.output(print(summary)),capture.output(print(by_drug[split=='held'])),
  'Each condition is evaluated twice, once per excluded constituent. These are dependent evaluations, not independent observations.',
  'Both table and main plot use actual per-fold predictions; no averaging into a training or held ensemble.',
  'Cell lines are not held out. Molecular exposure/culture differs from ALMANAC; retrospective development.',
  'No claim that high total-abundance correlation establishes interaction prediction.')
writeLines(report,file.path(root,'REPORT.txt'))
