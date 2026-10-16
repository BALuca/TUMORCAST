#!/usr/bin/env Rscript
# Correlation-first diagnostics on the identical held combinations and their
# measured individual baseline. Reuse prior plots' scales and dose-curve views.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
root <- '/work/results/program_pilot/67_full_shrunken_program_interactions'
previous <- '/work/results/program_pilot/64_nonlinear_response_composition'
setDTthreads(2L)
theme_set(theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))
save_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
save_plot <- function(p,name,width,height) {
  ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=width,height=height,dpi=150)
}
safe_cor <- function(x,y) if(length(x)<4 || sd(x)<1e-10 || sd(y)<1e-10) NA_real_ else cor(x,y)
new <- fread(file.path(root,'data/predictions.tsv.gz'))
stopifnot(uniqueN(new$fold_id)==6L,new[split=='held',.N,by=.(model,observation_id)][,all(N==1L)])
new[,method:=fifelse(model=='Full','Full interactions','Full + nonlinear')]
selected <- copy(new[selected_overall==TRUE])
selected[,method:='Training-selected variant']
baseline <- copy(new[model=='Full'])
baseline[,`:=`(method='Independent singles',predicted_interaction_log2=0,predicted_ratio_AB=ratio_A*ratio_B)]
old <- fread(file.path(previous,'data/predictions.tsv.gz'))[
  time_h==3 & fold_type=='single_partner' & model %chin% c('Nonlinear_control','Nonlinear_conditioned')]
old[,method:=fifelse(model=='Nonlinear_control','Function/dose control','Previous RNA + control')]
shared <- intersect(names(new),names(old))
z <- rbind(new[,..shared],selected[,..shared],baseline[,..shared],old[,..shared])
means <- z[,.(observed=mean(log2(ratio_AB/(ratio_A*ratio_B))),predicted=mean(predicted_interaction_log2),
  actual_count=100*mean(ratio_AB),predicted_count=100*mean(predicted_ratio_AB),
  single_product=100*mean(ratio_A*ratio_B),low_count=100*min(ratio_AB),high_count=100*max(ratio_AB)),
  by=.(method,split,anchor,partner,nominal_A_uM,nominal_B_uM,condition_id)]
means[,`:=`(centered_actual=actual_count-mean(actual_count),centered_predicted=predicted_count-mean(predicted_count)),
  by=.(method,split,anchor,nominal_A_uM)]
summary <- means[,.(n=.N,interaction_r=safe_cor(observed,predicted),
  combination_count_r=safe_cor(actual_count,predicted_count),
  within_first_drug_dose_count_r=safe_cor(centered_actual,centered_predicted)),by=.(method,split)]
by_partner <- means[split=='held',.(n=.N,interaction_r=safe_cor(observed,predicted),
  combination_count_r=safe_cor(actual_count,predicted_count)),by=.(method,partner)]
by_pair <- means[split=='held',.(n=.N,interaction_r=safe_cor(observed,predicted),
  combination_count_r=safe_cor(actual_count,predicted_count)),by=.(method,anchor,partner)]
save_table(means,'figure_condition_predictions')
save_table(summary,'correlation_summary')
save_table(by_partner,'held_correlations_by_partner_all_models')
save_table(by_pair,'held_correlations_by_pair')
colors <- c('Full interactions'='#D55E00','Full + nonlinear'='#0072B2','Independent singles'='#777777',
  'Function/dose control'='#009E73','Previous RNA + control'='#CC79A7','Training-selected variant'='#A6761D')
use <- c('Function/dose control','Full interactions','Full + nonlinear')
d <- copy(means[method %chin% use])
d[,`:=`(method=factor(method,levels=use),split=factor(split,levels=c('training','held')))]
labels <- d[,.(label=sprintf('Pearson r = %.2f',safe_cor(observed,predicted))),by=.(method,split)]
p <- ggplot(d,aes(observed,predicted))+geom_point(aes(color=anchor),alpha=.6,size=1.4)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.55)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.2)+
  facet_wrap(vars(split,method),ncol=3,scales='free')+theme(aspect.ratio=1)+
  scale_color_manual(values=c(alpelisib='#0072B2',dactolisib='#D55E00'))+
  labs(title='Full program interactions with shrinkage: training versus excluded partners',
    subtitle='Target is departure from measured individual effects; original40, BT20, 3h RNA',
    x='Measured log2[combination / individual-product]',y='Predicted log2 departure',color='First drug',
    caption='Every combination containing the predicted partner is excluded. Hyperparameters use inner partner exclusions only.')
save_plot(p,'interaction_training_vs_held',10.8,8)

d <- means[split=='held' & method %chin% c('Full interactions','Full + nonlinear')]
d[,method:=factor(method,levels=c('Full interactions','Full + nonlinear'))]
labels <- d[,.(label=sprintf('r = %.2f',safe_cor(observed,predicted))),by=.(method,partner)]
p <- ggplot(d,aes(observed,predicted))+geom_point(aes(color=anchor),size=1.3,alpha=.7)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.5)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.2,size=3)+
  facet_grid(method~partner,scales='free')+theme(aspect.ratio=1,strip.text=element_text(size=8))+
  scale_color_manual(values=c(alpelisib='#0072B2',dactolisib='#D55E00'))+
  labs(title='Every excluded partner: interaction predictions',
    x='Measured log2 departure',y='Predicted log2 departure',color='First drug')
save_plot(p,'held_interactions_each_partner',13,5.5)

use <- c('Independent singles','Full interactions','Full + nonlinear')
d <- means[split=='held' & method %chin% use]
d[,method:=factor(method,levels=use)]
labels <- d[,.(label=sprintf('r = %.2f',safe_cor(centered_actual,centered_predicted))),by=.(method,anchor)]
p <- ggplot(d,aes(centered_actual,centered_predicted))+geom_point(aes(color=partner),alpha=.7,size=1.4)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.55)+
  geom_text(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.04,vjust=1.2)+
  facet_wrap(vars(anchor,method),ncol=3,scales='free')+theme(aspect.ratio=1)+
  labs(title='Does the interaction improve combination differences beyond individual effects?',
    subtitle='Identical held-out observations; first-drug dose means removed only for this comparison',
    x='Measured deviation in cell count (pp)',y='Predicted deviation in cell count (pp)',color='Partner')
save_plot(p,'held_counts_vs_individual_baseline',10.8,8)

cairo_pdf(file.path(root,'figures/all12_combination_curves.pdf'),width=10.4,height=6.4)
for(aa in c('alpelisib','dactolisib')) for(bb in sort(unique(means$partner))) {
  d <- means[split=='held' & anchor==aa & partner==bb & method %chin% c(use,'Function/dose control')]
  dots <- d[method=='Independent singles']
  lines <- rbind(d[,.(nominal_A_uM,nominal_B_uM,method,value=predicted_count,endpoint='Cell count (% vehicle)')],
    d[,.(nominal_A_uM,nominal_B_uM,method,value=predicted,endpoint='Log2 departure from singles')])
  dots <- rbind(dots[,.(nominal_A_uM,nominal_B_uM,value=actual_count,lo=low_count,hi=high_count,endpoint='Cell count (% vehicle)')],
    z[split=='held' & method=='Independent singles' & anchor==aa & partner==bb,
      .(value=mean(log2(ratio_AB/(ratio_A*ratio_B))),lo=min(log2(ratio_AB/(ratio_A*ratio_B))),
        hi=max(log2(ratio_AB/(ratio_A*ratio_B))),endpoint='Log2 departure from singles'),by=.(nominal_A_uM,nominal_B_uM)])
  p <- ggplot(lines,aes(nominal_A_uM,value,color=method))+geom_line(linewidth=.7)+
    geom_linerange(data=dots,aes(x=nominal_A_uM,ymin=lo,ymax=hi),inherit.aes=FALSE,color='grey55')+
    geom_point(data=dots,aes(nominal_A_uM,value),inherit.aes=FALSE,size=1.7)+
    facet_grid(endpoint~nominal_B_uM,scales='free_y',labeller=labeller(nominal_B_uM=function(x) paste(bb,x,'µM')))+
    scale_x_log10(breaks=sort(unique(d$nominal_A_uM)))+scale_color_manual(values=colors)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    labs(title=paste(aa,'+',bb),subtitle='Black: measured biological mean/range. All combinations containing this partner excluded.',
      x=paste(aa,'dose (µM)'),y=NULL,color=NULL)
  print(p)
}
dev.off()

# Map learned coefficients directly; off-diagonal entries are halved because
# each design column is (a_i b_j + a_j b_i)/2. Nonlinear maps use80 basis terms.
coef <- fread(file.path(root,'data/coefficients.tsv.gz'))
cairo_pdf(file.path(root,'figures/full_interaction_coefficients.pdf'),width=9,height=8)
for(held in sort(unique(new$fold_id))) for(variant in c('Full','Smooth')) {
  d <- coef[fold_id==held & model==variant]
  terms <- tstrsplit(d$term,'__',fixed=TRUE)
  features <- unique(c(terms[[1]],terms[[2]]))
  mat <- matrix(0,length(features),length(features),dimnames=list(features,features))
  for(i in seq_len(nrow(d))) {
    v <- d$coefficient[i]/if(terms[[1]][i]==terms[[2]][i]) 1 else 2
    mat[terms[[1]][i],terms[[2]][i]] <- v
    mat[terms[[2]][i],terms[[1]][i]] <- v
  }
  order <- hclust(dist(mat),method='ward.D2')$order
  mat <- mat[order,order,drop=FALSE]
  limit <- max(abs(mat),1e-8)
  hm <- draw_heatmap(mat,scale_rows=FALSE,value_range=c(-limit,0,limit),value_colors=c('#2166AC','white','#B2182B'),
    hm_name='Interaction\ncoefficient',show_row_names=TRUE,show_column_names=TRUE,
    row_names_gp=grid::gpar(fontsize=if(variant=='Full') 6.5 else 5),
    column_names_gp=grid::gpar(fontsize=if(variant=='Full') 6.5 else 5),column_names_rot=90,
    title=paste('Excluded',held,'|',variant,'| _sat denotes tanh(program displacement)'))
  ComplexHeatmap::draw(hm)
  if(held=='lapatinib' && variant=='Full') {
    png(file.path(root,'figures/full_interaction_coefficient_example.png'),width=1350,height=1200,res=150)
    ComplexHeatmap::draw(hm)
    dev.off()
  }
}
dev.off()
audit <- fread(file.path(root,'data/validation.tsv.gz'))
stopifnot(nrow(audit)==12L,all(audit$source_overlap==0),all(audit$replay_error<1e-8))
for(i in seq_len(nrow(audit))) stopifnot(digest(file=file.path(root,'data',paste0('prediction_',audit$fold_id[i],'_',audit$model[i],'.tsv.gz')),
  algo='sha256')==audit$prediction_sha256[i])
report <- c('67: FULL INTERACTIONS WITH ELASTIC-NET SHRINKAGE',
  'All820 symmetric cross-program terms; no rank constraint or single-effect gene screen.',
  'Smooth extension uses all3240 products of x and tanh(x), permitting saturation.',
  'Alpha and lambda selected within training partners. Linear/smooth selection also uses inner exclusions only.',
  'Measured individual effects are fixed. Target=log2(AB/(A*B)); one BT20 context,3h measured separate RNA.',
  'No bootstrap, no joint RNA in predictions. Same198 conditions and6 partner exclusions as earlier versions.',
  '',capture.output(print(summary)), '',capture.output(print(by_partner)),
  '', 'All models are shown; choosing the best by these held results would reuse evaluation outcomes.',
  'Independent singles predict zero interaction, so their interaction correlation is undefined.',
  'Their actual-count correlations and complete dose curves are retained as the primary starting point.',
  'The function/dose control uses individual function and dose without RNA; Previous RNA + control adds RNA to that input.',
  'New Full/Smooth g models use only molecular cross-products, with the identical measured-single offset.',
  'Training summaries average predictions for cases included in multiple outer training sets, as before.',
  'This is retrospective development. Cross-partner transfer remains weak despite training fit.',
  '', 'Start with interaction_training_vs_held.pdf and held_interactions_each_partner.pdf.',
  'all12_combination_curves.pdf includes every pair; full_interaction_coefficients.pdf shows all fitted terms.')
writeLines(report,file.path(root,'REPORT.txt'))
scripts <- file.path('/work/src/program_pilot',c('67_full_shrunken_program_interactions.R','67_full_interaction_figures.R'))
save_table(data.table(path=scripts,sha256=vapply(scripts,digest,character(1),file=TRUE,algo='sha256')),'script_manifest')
fwrite(data.table(path=list.files(file.path(root,'figures'),pattern='\\.pdf$',full.names=TRUE)),file.path(root,'FIGURE_INDEX.tsv'),sep='\t')
writeLines('Full models, held correlations, diagnostic figures and validation complete.',file.path(root,'STATUS.txt'))
print(summary)
print(by_partner)
