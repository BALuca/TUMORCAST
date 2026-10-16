#!/usr/bin/env Rscript
# Primary comparison: independent singles versus the same singles plus shared g.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/73_bt20_almanac_program_interactions'
setDTthreads(2L)
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
save_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
save_plot <- function(p,name,width,height) {
  ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=width,height=height,dpi=150)
}
safe_cor <- function(x,y) if(length(x)<3L || sd(x)<1e-10 || sd(y)<1e-10) NA_real_ else cor(x,y)
z <- fread(file.path(root,'data/predictions.tsv.gz'))
stopifnot(uniqueN(z$fold_id)==6L,z[split=='held',.N,by=observation_id][,all(N==1L)])
z[,method:='Singles + program interactions']
baseline <- copy(z)
baseline[,`:=`(method='Independent singles',predicted_interaction_log2=0,predicted_ratio_AB=ratio_A*ratio_B)]
old <- fread('/work/results/program_pilot/71_bootstrap_screened_program_interactions/data/predictions.tsv.gz')[model=='Full']
old[,method:='BT20-only interactions']
stopifnot(identical(sort(old[split=='held',observation_id]),sort(z[split=='held',observation_id])))
shared <- intersect(names(z),names(old))
z <- rbind(z[,..shared],baseline[,..shared],old[,..shared])
# Average biological repeats within each condition and fit. Training predictions
# from different outer models are never averaged into an artificial ensemble.
means <- z[,.(observed=mean(log2(ratio_AB/(ratio_A*ratio_B))),predicted=mean(predicted_interaction_log2),
  actual_count=100*mean(ratio_AB),predicted_count=100*mean(predicted_ratio_AB),
  single_product=100*mean(ratio_A*ratio_B),low_count=100*min(ratio_AB),high_count=100*max(ratio_AB),
  low_g=min(log2(ratio_AB/(ratio_A*ratio_B))),high_g=max(log2(ratio_AB/(ratio_A*ratio_B)))),
  by=.(method,split,fold_id,anchor,partner,nominal_A_uM,nominal_B_uM,condition_id)]
stopifnot(means[split=='held',.N,by=method][,all(N==198L)])
means[,`:=`(centered_actual=actual_count-mean(actual_count),centered_predicted=predicted_count-mean(predicted_count)),
  by=.(method,split,fold_id,anchor,nominal_A_uM)]
summary <- means[,.(condition_predictions=.N,unique_conditions=uniqueN(condition_id),
  interaction_r=safe_cor(observed,predicted),combination_count_r=safe_cor(actual_count,predicted_count),
  median_absolute_count_error_pp=median(abs(actual_count-predicted_count)),
  within_10pp=sum(abs(actual_count-predicted_count)<=10)),by=.(method,split)]
by_partner <- means[split=='held',.(n=.N,interaction_r=safe_cor(observed,predicted),
  combination_count_r=safe_cor(actual_count,predicted_count),
  median_absolute_count_error_pp=median(abs(actual_count-predicted_count))),by=.(method,partner)]
by_pair <- means[split=='held',.(n=.N,interaction_r=safe_cor(observed,predicted),
  combination_count_r=safe_cor(actual_count,predicted_count)),by=.(method,anchor,partner)]
large <- means[split=='held',.(conditions=.N,median_observed_departure=median(observed),
  median_predicted_departure=median(predicted),median_observed_count=median(actual_count),
  median_predicted_count=median(predicted_count)),
  by=.(method,departure_group=fifelse(actual_count<single_product-10,'More than 10pp below independence',
    fifelse(actual_count>single_product+10,'More than 10pp above independence','Within 10pp of independence')))]
save_table(means,'figure_condition_predictions')
save_table(summary,'correlation_summary')
save_table(by_partner,'held_correlations_by_partner_all_models')
save_table(by_pair,'held_correlations_by_pair')
save_table(large,'large_departure_diagnostics')
colors <- c('Independent singles'='#777777','Singles + program interactions'='#D55E00','BT20-only interactions'='#0072B2')
primary <- c('Independent singles','Singles + program interactions')
# The central figure shows count predictions and the actual interaction target.
h <- means[split=='held' & method %chin% primary]
d <- rbind(h[,.(anchor,partner,x=actual_count,y=predicted_count,panel=method)],
  h[method==primary[2],.(anchor,partner,x=observed,y=predicted,panel='Departure from independence')])
d[,panel:=factor(panel,levels=c(primary,'Departure from independence'))]
labels <- d[,.(label=sprintf('Pearson r = %.3f\nn = %d conditions',safe_cor(x,y),.N)),by=panel]
p <- ggplot(d,aes(x,y))+geom_point(aes(color=anchor),size=1.4,alpha=.7)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.15,size=3)+
  facet_wrap(vars(panel),nrow=1,scales='free')+theme(aspect.ratio=1)+
  labs(title='BT20: excluded-partner combination predictions',
    subtitle='Left / middle: cell abundance (% vehicle). Right: log2[combination / single-drug product].',
    x='Measured',y='Predicted',color='First drug',
    caption='Shared g trained on BT20 + ALMANAC. Measured single effects supplied to both predictions.\nAll combinations containing each test partner excluded from both sources; retrospective evaluation.')
save_plot(p,'held_independent_vs_combined',10.5,4.3)

# Show all held partners rather than selecting successful exclusions.
d <- means[split=='held' & method==primary[2]]
labels <- d[,.(label=sprintf('r = %.3f; n = %d',safe_cor(observed,predicted),.N)),by=partner]
p <- ggplot(d,aes(observed,predicted))+geom_point(aes(color=anchor),size=1.5)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.15,size=3)+
  facet_wrap(vars(partner),ncol=3,scales='free')+theme(aspect.ratio=1)+
  labs(title='Does g recover departures from independent effects?',
    x='Measured log2 departure',y='Predicted log2 departure',color='First drug')
save_plot(p,'held_interactions_each_partner',9,7)

# Primary training/test diagnosis plus explicitly secondary supplementation check.
d <- means[method!='Independent singles']
d[,split:=factor(split,levels=c('training','held'))]
labels <- d[,.(label=sprintf('r = %.3f',safe_cor(observed,predicted))),by=.(method,split)]
p <- ggplot(d,aes(observed,predicted))+geom_point(aes(color=anchor),size=.9,alpha=.35)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.03,vjust=1.15)+
  facet_wrap(vars(split,method),ncol=2,scales='free')+theme(aspect.ratio=1)+
  labs(title='Interaction learning: training versus excluded partners',
    subtitle='Secondary comparison: BT20 alone versus BT20 supplemented with ALMANAC',
    x='Measured log2 departure',y='Predicted log2 departure',color='First drug',
    caption='Training: separate predictions from each outer fit; no averaging predictions across models.')
save_plot(p,'interaction_training_vs_held',8,8)

# Preserve every dose curve, including strong departures that the model misses.
cairo_pdf(file.path(root,'figures/all12_combination_curves.pdf'),width=10.4,height=6.4)
for(aa in c('alpelisib','dactolisib')) for(bb in sort(unique(means$partner))) {
  dd <- means[split=='held' & anchor==aa & partner==bb & method %chin% primary]
  dots <- dd[method=='Independent singles']
  lines <- rbind(dd[,.(nominal_A_uM,nominal_B_uM,method,value=predicted_count,endpoint='Cell count (% vehicle)')],
    dd[,.(nominal_A_uM,nominal_B_uM,method,value=predicted,endpoint='Log2 departure from singles')])
  dots <- rbind(dots[,.(nominal_A_uM,nominal_B_uM,value=actual_count,lo=low_count,hi=high_count,endpoint='Cell count (% vehicle)')],
    dots[,.(nominal_A_uM,nominal_B_uM,value=observed,lo=low_g,hi=high_g,endpoint='Log2 departure from singles')])
  p <- ggplot(lines,aes(nominal_A_uM,value,color=method))+geom_line(linewidth=.7)+
    geom_linerange(data=dots,aes(x=nominal_A_uM,ymin=lo,ymax=hi),inherit.aes=FALSE,color='grey55')+
    geom_point(data=dots,aes(nominal_A_uM,value),inherit.aes=FALSE,size=1.7)+
    facet_grid(endpoint~nominal_B_uM,scales='free_y',labeller=labeller(nominal_B_uM=function(x) paste(bb,x,'µM')))+
    scale_x_log10(breaks=sort(unique(dd$nominal_A_uM)))+scale_color_manual(values=colors)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    labs(title=paste(aa,'+',bb),subtitle='Black: measured biological mean/range. This partner excluded from both training sources.',
      x=paste(aa,'dose (µM)'),y=NULL,color=NULL)
  print(p)
}
dev.off()

coefficients <- fread(file.path(root,'data/coefficients.tsv.gz'))
cairo_pdf(file.path(root,'figures/program_interaction_coefficients.pdf'),width=8,height=8)
for(held in sort(unique(z$fold_id))) {
  dd <- coefficients[fold_id==held]
  terms <- tstrsplit(dd$term,'__',fixed=TRUE)
  programs <- sort(unique(c(terms[[1]],terms[[2]])))
  mat <- matrix(0,40,40,dimnames=list(programs,programs))
  for(i in seq_len(nrow(dd))) {
    value <- dd$coefficient[i]/if(terms[[1]][i]==terms[[2]][i]) 1 else 2
    mat[terms[[1]][i],terms[[2]][i]] <- value
    mat[terms[[2]][i],terms[[1]][i]] <- value
  }
  order <- hclust(dist(mat),method='ward.D2')$order
  mat <- mat[order,order,drop=FALSE]
  limit <- max(abs(mat),1e-8)
  hm <- draw_heatmap(mat,scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    hm_name='Coefficient',show_row_names=TRUE,show_column_names=TRUE,
    row_names_gp=grid::gpar(fontsize=6),column_names_gp=grid::gpar(fontsize=6),column_names_rot=90,
    title=paste('Shared BT20 + ALMANAC g; excluded',held))
  ComplexHeatmap::draw(hm)
}
dev.off()
# Validate persisted model predictions and exclusion manifests independently.
audit <- fread(file.path(root,'data/validation.tsv.gz'))
stopifnot(nrow(audit)==6L,all(audit$source_overlap==0),all(audit$replay_error<1e-8))
for(i in seq_len(nrow(audit))) {
  held <- audit$fold_id[i]
  fit <- qread(file.path(root,'models',paste0(held,'_Full.qs')))
  stopifnot(!held %chin% fit$training_partners,
    digest(file=file.path(root,'data',paste0('prediction_',held,'_Full.tsv.gz')),algo='sha256')==audit$prediction_sha256[i])
}
report <- c('73: BT20 PLUS ALMANAC, EXPLICIT PROGRAM INTERACTIONS',
  'Primary: independent single-drug product vs same product times 2^g, on held BT20 combinations.',
  'g target: log2[observed combination/(measured single A * measured single B)].',
  'BT20 and ALMANAC each half total training weight. Inner selection uses BT20 exclusions only.',
  'Each outer and inner partner excluded from both sources; ALMANAC sharing exact single references purged.',
  'Only lapatinib overlaps drug panels. Partial L1000 vs full Tahoe scores use separate outcome-free RMS scales.',
  'One shared coefficient vector; no drug identities or combination RNA as model features.',
  '',capture.output(print(summary)), '',capture.output(print(by_partner)), '',capture.output(print(large)),
  'Independent singles predict g=0; interaction correlation is undefined, not zero.',
  '198 held conditions / 496 biological observations. Training rows occur in five different fitted models.',
  'No outcome filtering. These repeatedly examined partner exclusions are retrospective, not independent confirmation.',
  'A high total count correlation driven by individual effects is not evidence of learned interactions.',
  'Cross-assay program comparability and transfer of cell-context relationships remain assumptions.')
writeLines(report,file.path(root,'REPORT.txt'))
writeLines('COMPLETE: joint training, all six exclusions, primary comparison and diagnostics saved.',file.path(root,'STATUS.txt'))
fwrite(data.table(path=list.files(file.path(root,'figures'),pattern='\\.pdf$',full.names=TRUE)),file.path(root,'FIGURE_INDEX.tsv'),sep='\t')
print(summary)
print(by_partner)
