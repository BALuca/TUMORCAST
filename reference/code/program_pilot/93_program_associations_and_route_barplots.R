#!/usr/bin/env Rscript
# Diagnostic barplots for frozen v93. No model refitting or program selection.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2)})
setDTthreads(1L)
root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
figures <- file.path(root,'figures/program_and_route_diagnostics')
dir.create(figures,recursive=TRUE,showWarnings=FALSE)
programs <- sprintf('P%02d',1:40)
selection <- fread(file.path(root,'data/selected_doses.tsv'))
model <- qread(file.path(root,'models/final.qs'))
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
m <- merge(as.data.table(archive$prep$m),selection,by=c('drug','dose'))
y <- fread(file.path(root,'data/individual_training_outcomes.tsv'))[assay=='secondary',.(depmap_id,drug,dose,observed_log2=y)]
profiles <- merge(m,y,by=c('depmap_id','drug','dose'))
setorder(profiles,drug,depmap_id)
drugs <- selection$drug
g <- match(profiles$drug,drugs)
context <- model$context
D <- as.matrix(profiles[,..programs])
X <- (D-context$mu[g,,drop=FALSE])/context$sigma[g,,drop=FALSE]
gain <- exp(context$Q%*%model$par$gain)
shift <- context$Q%*%model$par$shift
logits <- context$Q%*%model$par$mix
alpha <- exp(logits-apply(logits,1,max))
alpha <- alpha/rowSums(alpha)
calibrated <- X*gain[g,,drop=FALSE]+shift[g,,drop=FALSE]
latent <- calibrated%*%model$par$W
h <- tanh(latent)
hzero <- tanh(((-context$mu/context$sigma)*gain+shift)%*%model$par$W)
activity <- (h-hzero[g,,drop=FALSE])*alpha[g,,drop=FALSE]
profiles[,predicted_log2:=as.numeric(activity%*%model$par$C[1,])]
profiles[,`:=`(observed_sensitivity=-observed_log2,predicted_sensitivity=-predicted_log2)]
stopifnot(nrow(profiles)==286L,uniqueN(profiles$drug)==10L)
# Sensitivity is -log2 abundance: positive correlation means more displacement
# accompanies greater sensitivity. These signs reverse the previous abundance plots.
long <- melt(profiles,id.vars=c('depmap_id','drug','dose','observed_sensitivity','predicted_sensitivity'),
  measure.vars=programs,variable.name='program',value.name='displacement',variable.factor=FALSE)
correlations <- long[,.(cells=.N,observed_r=cor(displacement,observed_sensitivity),
  model_marginal_r=cor(displacement,predicted_sensitivity)),by=.(drug,dose,program)]
# Actual local model beta: derivative of sensitivity with respect to a program,
# scaled to +1 training RNA SD of that drug/program, averaged over observed cells.
# This includes output coefficients, drug mixture, gain and tanh saturation.
route_derivative <- -(1-h^2)*alpha[g,,drop=FALSE]*matrix(model$par$C[1,],nrow(profiles),2,byrow=TRUE)
effects_per_SD <- (route_derivative%*%t(model$par$W))*gain[g,,drop=FALSE]
raw_effects <- effects_per_SD/context$sigma[g,,drop=FALSE]
by_drug <- rbindlist(lapply(seq_along(drugs),function(j) {
  ii <- which(g==j)
  data.table(drug=drugs[j],dose=selection$dose[j],program=programs,
    mean_sensitivity_effect_per_program_SD=colMeans(effects_per_SD[ii,,drop=FALSE]),
    mean_sensitivity_effect_per_raw_unit=colMeans(raw_effects[ii,,drop=FALSE]))
}))
comparison <- merge(correlations,by_drug,by=c('drug','dose','program'))
fwrite(comparison,file.path(root,'data/program_association_vs_model_effect.tsv'),sep='\t')
# Numerical derivative checks on real fitted profiles validate signs and units.
recalculate <- function(Dnew) {
  xx <- (Dnew-context$mu[g,,drop=FALSE])/context$sigma[g,,drop=FALSE]
  hh <- tanh((xx*gain[g,,drop=FALSE]+shift[g,,drop=FALSE])%*%model$par$W)
  -as.numeric(((hh-hzero[g,,drop=FALSE])*alpha[g,,drop=FALSE])%*%model$par$C[1,])
}
for(j in c(1,21,35,40)) {
  step <- context$sigma[g,j]*1e-5
  plus <- minus <- D
  plus[,j] <- plus[,j]+step
  minus[,j] <- minus[,j]-step
  numerical <- (recalculate(plus)-recalculate(minus))/(2e-5)
  stopifnot(max(abs(numerical-effects_per_SD[,j]))<1e-7)
}
short <- c('Azacytidine','Bendamustine','Bortezomib','Busulfan','Cytarabine','Lapatinib',
  'Topotecan','Vandetanib','Vismodegib','Crizotinib')
panel_levels <- paste0(short,' | ',selection$dose,' µM')
add_panel <- function(z) {
  z[,panel:=factor(paste0(short[match(drug,drugs)],' | ',dose,' µM'),levels=panel_levels)]
  z[,program:=factor(program,levels=rev(programs))]
  z
}
comparison <- add_panel(comparison)
cal <- add_panel(fread(file.path(root,'data/learned_drug_program_calibration.tsv'))[
  ,dose:=selection$dose[match(drug,drugs)]])
# Ten facets contain 400 bars; horizontal labels keep all programs readable.
bar_theme <- theme_bw(base_size=9)+theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
  strip.text=element_text(size=8),axis.text.y=element_text(size=6),legend.position='bottom')
save_bar <- function(p,stem) {
  ggsave(file.path(figures,paste0(stem,'.pdf')),p,width=15,height=13,device=cairo_pdf)
  ggsave(file.path(figures,paste0(stem,'.png')),p,width=15,height=13,dpi=110)
}
p1 <- ggplot(comparison,aes(observed_r,program,fill=observed_r>0))+
  geom_vline(xintercept=0,color='grey55',linewidth=.3)+geom_col(width=.8)+
  facet_wrap(~panel,ncol=5)+scale_x_continuous(limits=c(-1,1),breaks=c(-1,0,1))+
  scale_fill_manual(values=c('#B45A68','#39799D'),guide='none')+bar_theme+
  labs(title='Observed program correlations with sensitivity, within each drug/dose',
    x='Pearson r with sensitivity (-log2 PRISM abundance)',y=NULL,
    caption='All 40 measured programs; 286 matched profiles / 30 cell lines. Positive = higher displacement associates with greater sensitivity.\nDescriptive observed associations, not held-out predictions. Values reverse signs from the earlier abundance scatters.')
save_bar(p1,'observed_program_sensitivity_correlations')
p2 <- ggplot(comparison,aes(mean_sensitivity_effect_per_program_SD,program,
  fill=mean_sensitivity_effect_per_program_SD>0))+
  geom_vline(xintercept=0,color='grey55',linewidth=.3)+geom_col(width=.8)+
  facet_wrap(~panel,ncol=5,scales='free_x')+
  scale_fill_manual(values=c('#B45A68','#39799D'),guide='none')+bar_theme+
  labs(title='v93 effective program effects: the actual nonlinear decoder',
    x='Mean local sensitivity change for +1 drug-specific program SD (log2 abundance units)',y=NULL,
    caption='All-data fitted model. Conditional local derivatives, holding the other programs and drug footprint fixed; not marginal correlations.\nIncludes route weights, drug mixing, learned gain and nonlinear saturation. Positive = greater predicted sensitivity.')
save_bar(p2,'effective_program_sensitivity_betas')
pairbars <- melt(comparison,id.vars=c('drug','dose','program','panel'),
  measure.vars=c('observed_r','model_marginal_r'),variable.name='source',value.name='r')
p3 <- ggplot(pairbars,aes(r,program,fill=source))+
  geom_vline(xintercept=0,color='grey55',linewidth=.3)+geom_col(position=position_dodge(width=.85),width=.78)+
  facet_wrap(~panel,ncol=5)+scale_x_continuous(limits=c(-1,1),breaks=c(-1,0,1))+
  scale_fill_manual(values=c(observed_r='#39799D',model_marginal_r='#C27535'),
    labels=c(observed_r='Observed sensitivity',model_marginal_r='Fitted model sensitivity'),name=NULL)+bar_theme+
  labs(title='Observed versus model-implied program associations: identical cells and coordinates',
    x='Pearson r with sensitivity (-log2 abundance)',y=NULL,
    caption='Both bars are marginal correlations across the same measured cells, making their directions directly comparable.\nThe model bars use the all-data fitted decoder, not held-out predictions; this explains what it fitted, not its generalization.')
save_bar(p3,'observed_vs_model_program_correlations')
loadings <- data.table(program=rep(programs,2),route=factor(rep(1:2,each=40)),
  weight=as.numeric(model$par$W))
loadings[,program:=factor(program,levels=rev(programs))]
p4 <- ggplot(loadings,aes(weight,program,fill=route))+geom_vline(xintercept=0,color='grey55',linewidth=.3)+
  geom_col(width=.8)+facet_wrap(~route,ncol=2)+scale_fill_manual(values=c('#39799D','#A96E9E'),guide='none')+
  theme_bw(base_size=10)+theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=8))+
  labs(title='v93: two learned program directions entering tanh',x='Learned preactivation weight W',y=NULL,
    caption='Each direction reads all 40 drug-calibrated programs. These are mathematical route loadings, not viability beta coefficients.\nOutput coefficients, drug mixtures and nonlinear saturation determine the functional effect.')
ggsave(file.path(figures,'two_route_program_loadings.pdf'),p4,width=8,height=9,device=cairo_pdf)
ggsave(file.path(figures,'two_route_program_loadings.png'),p4,width=8,height=9,dpi=130)
# Separate molecular anchors from learned adjustments: gain=1 and shift=0 mean
# no adjustment beyond the per-drug RNA mean/SD standardization.
pg <- ggplot(cal,aes(log2(learned_gain),program,fill=log2(learned_gain)>0))+
  geom_vline(xintercept=0,color='grey55',linewidth=.3)+geom_col(width=.8)+facet_wrap(~panel,ncol=5)+
  scale_fill_manual(values=c('#B45A68','#39799D'),guide='none')+bar_theme+
  labs(title='Learned drug/program gain adjustments beyond measured RNA scales',
    x='log2 learned gain (0 = gain 1; +1 = doubling)',y=NULL,
    caption='The starting scale is each drug/program measured RNA SD. These gains are jointly learned shared footprint-map adjustments.')
save_bar(pg,'learned_program_gain_adjustments')
ps <- ggplot(cal,aes(learned_shift,program,fill=learned_shift>0))+
  geom_vline(xintercept=0,color='grey55',linewidth=.3)+geom_col(width=.8)+facet_wrap(~panel,ncol=5)+
  scale_fill_manual(values=c('#B45A68','#39799D'),guide='none')+bar_theme+
  labs(title='Learned drug/program shift adjustments beyond measured RNA centering',
    x='Learned shift in drug-specific program SD units (0 = no added shift)',y=NULL,
    caption='Shifts act inside the functional calibration layer. They do not alter measured RNA or the molecular prediction targets.')
save_bar(ps,'learned_program_shift_adjustments')
pa <- ggplot(cal,aes(measured_SD,program))+geom_col(fill='#6C88A3',width=.8)+
  facet_wrap(~panel,ncol=5)+scale_x_log10()+bar_theme+
  labs(title='Starting molecular scale: measured program SD for each drug/dose',
    x='Measured RNA displacement SD (log axis)',y=NULL,
    caption='These molecular scales are fixed anchors in the all-data fit; the gain/shift panels show the additional learned adjustments.')
save_bar(pa,'measured_drug_program_scale_anchors')
summary <- cal[,.(gain_min=min(learned_gain),gain_max=max(learned_gain),
  shift_min=min(learned_shift),shift_max=max(learned_shift)),by=.(drug,dose)]
fwrite(summary,file.path(root,'data/drug_calibration_adjustment_summary.tsv'),sep='\t')
# Export raw direction loadings together with assay output coefficients.
fwrite(data.table(route=1:2,secondary_output_coefficient=model$par$C[1,],
  nci_FF_output_coefficient=model$par$C[2,],nci_FG_output_coefficient=model$par$C[3,]),
  file.path(root,'data/route_assay_output_coefficients.tsv'),sep='\t')
print(summary)
print(comparison[drug=='Lapatinib ditosylate'][order(-abs(observed_r))][1:8,
  .(program,observed_r,model_marginal_r,mean_sensitivity_effect_per_program_SD)])
print(model$par$C)
writeLines(c('Observed sensitivity = -log2 PRISM abundance; positive bars mean greater sensitivity.',
  'Observed correlations and model marginal correlations use identical measured profiles.',
  'Effective beta is the derivative of sensitivity with respect to one program, multiplied by its training RNA SD and averaged over profiles.',
  'Effective betas are conditional effects; marginal correlations can differ because other programs covary.',
  'Route weight bars show the two shared tanh preactivation directions W; readout C, mixture alpha and saturation are separate.',
  'Gain/shift bars are additions to the fixed per-drug/program measured mean/SD anchors.',
  'All-data fitted representation diagnostics are descriptive and do not replace held-out performance.',
  'Numerical central differences checked effective derivatives for P01/P21/P35/P40 on all 286 profiles.'),
  file.path(figures,'BARPLOT_DEFINITIONS.txt'))
cat('Complete: observed correlations, effective betas, fitted marginal associations, route loadings and drug calibration barplots.\n')
