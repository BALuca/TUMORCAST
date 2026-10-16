#!/usr/bin/env Rscript
# Inspect added predictive value beyond the SAME observed individual effects,
# along with each excluded-partner model's learned program-set interactions.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/66_bootstrapped_program_interaction_patterns'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
nonlinear <- '/work/results/program_pilot/64_nonlinear_response_composition'
setDTthreads(2L)
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
save_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
save_plot <- function(p,name,width,height) {
  ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=width,height=height,dpi=150)
}
safe_cor <- function(x,y) if(length(x)<4 || sd(x)<1e-10 || sd(y)<1e-10) NA_real_ else cor(x,y)
new <- fread(file.path(root,'data/predictions.tsv.gz'))
audit <- fread(file.path(root,'data/validation.tsv.gz'))
stopifnot(uniqueN(new$fold_id)==6L,all(audit$source_overlap==0),all(audit$swap_error<1e-10),
  all(audit$replay_error<1e-10),all(audit$zero_error==0))
stopifnot(new[split=='held',.N,by=observation_id][,all(N==1L)],nrow(new[split=='held'])==496L)
new[,method:='Bootstrap program interactions']
old <- fread(file.path(nonlinear,'data/predictions.tsv.gz'))[
  time_h==3 & fold_type=='single_partner' & model %chin% c('Nonlinear_control','Nonlinear_conditioned')]
old[,method:=fifelse(model=='Nonlinear_control','Function/dose control','Previous RNA kernel')]
shared <- intersect(names(new),names(old))
z <- rbind(new[,..shared],old[,..shared])
unbagged <- copy(new)
unbagged[,`:=`(method='Flexible single fit',predicted_interaction_log2=unbagged_interaction_log2,
  predicted_ratio_AB=2^(predicted_single_sum_log2+unbagged_interaction_log2))]
previous_sparse <- fread('/work/results/program_pilot/65_sparse_program_interaction_patterns/data/predictions.tsv.gz')
previous_sparse[,method:='Previous sparse fit']
z <- rbind(z,unbagged[,..shared],previous_sparse[,..shared])
baseline <- copy(new)
baseline[,`:=`(method='Independent singles',predicted_interaction_log2=0,predicted_ratio_AB=ratio_A*ratio_B)]
z <- rbind(z,baseline[,..shared])
z[,`:=`(measured_departure=log2(ratio_AB/(ratio_A*ratio_B)),
  independent_pct=100*ratio_A*ratio_B,observed_pct=100*ratio_AB,predicted_pct=100*predicted_ratio_AB)]
means <- z[,.(observed_pct=mean(observed_pct),predicted_pct=mean(predicted_pct),
  independent_pct=mean(independent_pct),observed_departure=mean(measured_departure),
  predicted_departure=mean(predicted_interaction_log2),low=min(observed_pct),high=max(observed_pct)),
  by=.(method,split,anchor,partner,nominal_A_uM,nominal_B_uM,condition_id)]
means[,`:=`(observed_centered=observed_pct-mean(observed_pct),
  predicted_centered=predicted_pct-mean(predicted_pct)),by=.(method,split,anchor,nominal_A_uM)]
means[,`:=`(absolute_error_pp=abs(observed_pct-predicted_pct),
  independent_error_pp=abs(observed_pct-independent_pct))]
metrics <- means[,.(n=.N,departure_r=safe_cor(observed_departure,predicted_departure),
  count_r=safe_cor(observed_pct,predicted_pct),within_first_drug_dose_r=safe_cor(observed_centered,predicted_centered),
  median_absolute_count_difference_pp=median(absolute_error_pp),
  median_gain_over_independent_pp=median(independent_error_pp-absolute_error_pp),
  conditions_closer_than_independent=sum(absolute_error_pp<independent_error_pp-1e-8)),by=.(method,split)]
by_anchor <- means[,.(n=.N,departure_r=safe_cor(observed_departure,predicted_departure),
  within_first_drug_dose_r=safe_cor(observed_centered,predicted_centered),
  median_absolute_count_difference_pp=median(absolute_error_pp),
  median_gain_over_independent_pp=median(independent_error_pp-absolute_error_pp)),by=.(method,split,anchor)]
by_pair <- means[split=='held',.(n=.N,departure_r=safe_cor(observed_departure,predicted_departure),
  median_gain_over_independent_pp=median(independent_error_pp-absolute_error_pp),
  median_absolute_count_difference_pp=median(absolute_error_pp)),by=.(method,anchor,partner)]
save_table(means,'figure_condition_predictions')
save_table(metrics,'model_summary')
save_table(by_anchor,'model_by_first_drug')
save_table(by_pair,'model_by_pair')
colors <- c('Bootstrap program interactions'='#D55E00','Independent singles'='#777777',
  'Function/dose control'='#009E73','Previous RNA kernel'='#0072B2',
  'Flexible single fit'='#CC79A7','Previous sparse fit'='#A6761D')
use <- c('Function/dose control','Flexible single fit','Bootstrap program interactions')
d <- means[method %chin% use]
d[,`:=`(method=factor(method,levels=use),split=factor(split,levels=c('training','held')))]
labels <- d[,.(label=sprintf('r = %.2f',safe_cor(observed_departure,predicted_departure))),by=.(method,split)]
p <- ggplot(d,aes(observed_departure,predicted_departure))+
  geom_point(aes(color=anchor),size=1.4,alpha=.6)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.55)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.2)+
  facet_wrap(vars(split,method),ncol=3,scales='free')+theme(aspect.ratio=1)+
  scale_color_manual(values=c(alpelisib='#0072B2',dactolisib='#D55E00'))+
  labs(title='Can bootstrapped program relationships predict departure from individual effects?',
    subtitle='Original40, 3 h RNA; combinations containing the excluded partner withheld under both first drugs',
    x='Measured log2[combination / individual-product]',y='Predicted log2 departure',color='First drug',
    caption='Negative: fewer cells than the independent expectation. Zero: independent individual effects. All198 conditions retained.')
save_plot(p,'interaction_training_and_excluded_partners',10.8,8)

d <- copy(by_pair[method %chin% c('Function/dose control','Flexible single fit','Bootstrap program interactions','Previous sparse fit')])
d[,partner:=factor(partner,levels=rev(sort(unique(partner))))]
p <- ggplot(d,aes(median_gain_over_independent_pp,partner,color=method))+
  geom_vline(xintercept=0,color='grey50')+geom_point(position=position_dodge(width=.45),size=2.5)+
  facet_wrap(~anchor,nrow=1)+scale_color_manual(values=colors)+
  labs(title='Does adding interaction improve on the SAME individual effects?',
    subtitle='Median across the complete held-out dose curve for each pair; positive values mean better predictions',
    x='Reduction in absolute cell-count error (percentage points)',y='Excluded partner',color=NULL)
save_plot(p,'gain_beyond_individual_effects',9,4.8)

use <- c('Independent singles','Function/dose control','Bootstrap program interactions')
d <- means[split=='held' & method %chin% use]
d[,method:=factor(method,levels=use)]
labels <- d[,.(label=sprintf('r = %.2f\nmedian error = %.1f pp',
  safe_cor(observed_centered,predicted_centered),median(absolute_error_pp))),by=.(method,anchor)]
p <- ggplot(d,aes(observed_centered,predicted_centered))+
  geom_point(aes(color=partner),size=1.4,alpha=.7)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.55)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.1,size=3.2)+
  facet_wrap(vars(anchor,method),ncol=3,scales='free')+theme(aspect.ratio=1)+
  labs(title='Combination differences after accounting for first-drug dose',
    subtitle='Same held-out combinations in every panel; measured individual effects provide the starting point',
    x='Measured deviation at fixed first-drug dose (pp)',y='Predicted deviation at fixed first-drug dose (pp)',
    color='Partner',caption='Comparisons are across partners and partner doses in BT20; absolute dose curves are provided separately.')
save_plot(p,'held_combinations_vs_individual_baseline',10.8,8.2)

# Preserve all twelve complete dose curves, including replicate ranges.
cairo_pdf(file.path(root,'figures/all12_combination_curves.pdf'),width=10.4,height=6.4)
for(aa in c('alpelisib','dactolisib')) for(bb in sort(unique(means$partner))) {
  d <- means[split=='held' & anchor==aa & partner==bb & method %chin% use]
  dots <- d[method=='Independent singles']
  lines <- rbind(d[,.(nominal_A_uM,nominal_B_uM,method,value=predicted_pct,endpoint='Cell count (% vehicle)')],
    d[,.(nominal_A_uM,nominal_B_uM,method,value=independent_pct-predicted_pct,endpoint='Extra inhibition (pp)')])
  dots <- rbind(dots[,.(nominal_A_uM,nominal_B_uM,value=observed_pct,lo=low,hi=high,endpoint='Cell count (% vehicle)')],
    z[split=='held' & anchor==aa & partner==bb & method=='Independent singles',
      .(value=mean(independent_pct-observed_pct),lo=min(independent_pct-observed_pct),
        hi=max(independent_pct-observed_pct),endpoint='Extra inhibition (pp)'),by=.(nominal_A_uM,nominal_B_uM)])
  p <- ggplot(lines,aes(nominal_A_uM,value,color=method))+geom_line(linewidth=.75)+
    geom_linerange(data=dots,aes(x=nominal_A_uM,ymin=lo,ymax=hi),inherit.aes=FALSE,color='grey50')+
    geom_point(data=dots,aes(nominal_A_uM,value),inherit.aes=FALSE,size=1.7)+
    facet_grid(endpoint~nominal_B_uM,scales='free_y',labeller=labeller(nominal_B_uM=function(x) paste(bb,x,'µM')))+
    scale_x_log10(breaks=sort(unique(d$nominal_A_uM)))+scale_color_manual(values=colors)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    labs(title=paste(aa,'+',bb),subtitle='Black: measured biological mean/range; the entire partner is excluded from combination learning',
      x=paste(aa,'dose (µM)'),y=NULL,color=NULL)
  print(p)
  if(aa=='alpelisib' && bb=='lapatinib') ggsave(file.path(root,'figures/lapatinib_curve_example.png'),p,width=10.4,height=6.4,dpi=150)
}
dev.off()

# W is identifiable only on the observed response span. U/V component axes have
# sign/order/rotation ambiguities; show each fitted fold separately, not an average.
annotations <- fread(file.path(prior,'data/original40_program_annotations.tsv'))
cairo_pdf(file.path(root,'figures/learned_program_interaction_patterns.pdf'),width=10,height=8.5)
for(held in sort(unique(new$fold_id))) {
  fit <- qread(file.path(root,'models',paste0('3h_exclude_',held,'.qs')))
  w <- fit$W
  order <- hclust(dist(w),method='ward.D2')$order
  mat <- do.call(cbind,lapply(seq_len(ncol(fit$u)),function(k) cbind(fit$u[,k],fit$v[,k])))
  rownames(mat) <- rownames(fit$u)
  colnames(mat) <- paste0('Pattern ',rep(seq_len(ncol(fit$u)),each=2),': ',rep(c('U','V'),ncol(fit$u)))
  mat <- mat[order,,drop=FALSE]
  labels <- substr(gsub('GOBP ','',annotations$display_name[match(rownames(mat),annotations$program)],fixed=TRUE),1,48)
  limit <- max(abs(mat),1e-8)
  hm <- draw_heatmap(mat,scale_rows=FALSE,value_range=c(-limit,0,limit),
    value_colors=c('#2166AC','white','#B2182B'),hm_name='Factor\nloading',
    show_row_names=TRUE,show_column_names=TRUE,row_labels=labels,row_names_gp=grid::gpar(fontsize=7),
    row_names_side='left',row_names_max_width=grid::unit(6,'cm'),
    column_names_rot=45,title=paste('Excluded',held,'| representative bootstrap member',fit$representative_member,'| U interacts with V'))
  ComplexHeatmap::draw(hm,heatmap_legend_side='bottom',annotation_legend_side='bottom')
  if(held=='lapatinib') {
    png(file.path(root,'figures/representative_program_patterns.png'),width=1500,height=1275,res=150)
    ComplexHeatmap::draw(hm,heatmap_legend_side='bottom',annotation_legend_side='bottom')
    dev.off()
  }
  ww <- w[order,order,drop=FALSE]
  limit <- max(abs(ww),1e-8)
  hm <- draw_heatmap(ww,scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    hm_name='Program-pair\ncoefficient',show_row_names=TRUE,show_column_names=TRUE,
    row_names_gp=grid::gpar(fontsize=6.5),column_names_gp=grid::gpar(fontsize=6.5),column_names_rot=90,
    title=paste('Excluded',held,'| ensemble mean program interaction matrix'))
  ComplexHeatmap::draw(hm)
  if(held=='lapatinib') {
    png(file.path(root,'figures/ensemble_program_interactions.png'),width=1500,height=1275,res=150)
    ComplexHeatmap::draw(hm)
    dev.off()
  }
  wc <- vapply(fit$members,function(f) as.vector(f$W),numeric(1600))
  stability <- matrix(pmax(rowMeans(wc>1e-8),rowMeans(wc< -1e-8)),40,40,dimnames=dimnames(w))
  hm <- draw_heatmap(stability[order,order],scale_rows=FALSE,value_range=c(0,.5,1),
    value_colors=c('white','#9ECAE1','#08519C'),hm_name='Bootstrap\nsign fraction',
    show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=6.5),
    column_names_gp=grid::gpar(fontsize=6.5),column_names_rot=90,
    title=paste('Excluded',held,'| fraction of all20 members supporting the predominant nonzero sign'))
  ComplexHeatmap::draw(hm)
}
dev.off()

# Apply EXACTLY the same learned interaction to64's frozen molecular individual
# decoder. No new g fit is allowed to absorb its single-effect prediction errors.
f <- fread(file.path(nonlinear,'data/predictions.tsv.gz'))[
  time_h==3 & fold_type=='single_partner' & model=='Nonlinear_bridge' & split=='held',
  .(observation_id,fold_id,molecular_offset=predicted_single_sum_log2)]
bridge <- merge(new[split=='held'],f,by=c('observation_id','fold_id'))
stopifnot(nrow(bridge)==496L)
bridge[,`:=`(individual_prediction_pct=100*2^molecular_offset,
  combination_prediction_pct=100*2^(molecular_offset+predicted_interaction_log2))]
save_table(bridge,'frozen_individual_decoder_plus_same_interaction')
report <- c('66: BOOTSTRAPPED PROGRAM-SET INTERACTIONS',
  'Same original40 bank, BT20 observations and whole-partner exclusions. 3h primary only for speed.',
  'g = sum_k [(a U_k)(b V_k)+(a V_k)(b U_k)]/2 in each bootstrap member.',
  'Ranks1,3,6 and penalties chosen using the ensemble prediction in inner partner exclusions.',
  'Grouped Bayesian bootstrap reweights entire partners: both anchors, all doses/replicates kept together.',
  'Six draws during tuning;20 at final fitting. Sparse penalty alpha0.5; penalties0.1,0.03,0.01,0.003.',
  'Flexible single fit uses the identical chosen settings without bootstrap reweighting.',
  'The averaged interaction matrix need not have rank<=6; that constraint applies to each member.',
  'Primary target is log2(measured AB / measured single-product). No individual-decoder errors enter g training.',
  'Zero interaction is an eligible candidate. All six excluded partners and all198 dose conditions retained.',
  'All methods use the identical measured individual effects; functional-only comparator reused from64.',
  '',capture.output(print(metrics)), '',capture.output(print(by_anchor[split=='held'])),
  '', 'Selected patterns:',capture.output(print(audit[,.(fold_id,rank,effective_patterns,lambda)])),
  '', 'Interpretation boundaries:',
  'Any gain over independent singles is distinct from added value over the function/dose control.',
  'Low rank constrains response interactions but does not establish pathway mechanism.',
  'U/V axes need not be unique. The member closest to the ensemble W is shown, with mean W and sign stability.',
  'Bootstrap stability measures training sensitivity, not biological validation or calibrated uncertainty.',
  'These outcomes have been inspected before: retrospective model development, one BT20 background.',
  'BT20 combination RNA is absent. Diaz joint-RNA diagnostics are separate and cannot validate BT20 mechanisms.',
  'Start with gain_beyond_individual_effects.pdf, interaction_training_and_excluded_partners.pdf,',
  'all12_combination_curves.pdf and learned_program_interaction_patterns.pdf.')
writeLines(report,file.path(root,'REPORT.txt'))
print(metrics)
print(audit[,.(fold_id,rank,effective_patterns,lambda)])
