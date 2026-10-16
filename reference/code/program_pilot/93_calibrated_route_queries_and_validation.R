#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest)})
root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
v87 <- '/work/results/program_pilot/87_selected_dose_joint_program_bridge'
programs <- sprintf('P%02d',1:40)
K <- 2L
assays <- c('secondary','nci_FF','nci_FG')
write_table <- function(x,n) fwrite(x,file.path(root,'data',paste0(n,'.tsv')),sep='\t')
pearson <- function(x,y) if(length(x)>=4L && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
archive <- readRDS(archive_path)
base <- as.data.table(archive$prep$raw)
selection <- fread(file.path(root,'data/selected_doses.tsv'))
drugs <- selection$drug
molecular <- merge(as.data.table(archive$prep$m),selection,by=c('drug','dose'))
getbaseline <- function(ids) as.matrix(base[match(ids,depmap_id),..programs])
route_node <- function(X,g,par,Q,mu,sigma) {
  gain <- exp(Q%*%par$gain)
  shift <- Q%*%par$shift
  logits <- Q%*%par$mix
  alpha <- exp(logits-apply(logits,1,max))
  alpha <- alpha/rowSums(alpha)
  zero <- -mu/sigma
  calibrated <- X*gain[g,,drop=FALSE]+shift[g,,drop=FALSE]
  calibrated0 <- (zero*gain+shift)[g,,drop=FALSE]
  h <- tanh(calibrated%*%par$W)
  h0 <- tanh(calibrated0%*%par$W)
  activity <- (h-h0)*alpha[g,,drop=FALSE]
  list(activity=activity,X=X,g=g,gain=gain,zero=zero,calibrated=calibrated,
    calibrated0=calibrated0,h=h,h0=h0,alpha=alpha)
}
pair_features <- function(A,B) cbind(A[,1]*B[,1],A[,1]*B[,2]+A[,2]*B[,1],A[,2]*B[,2])
predict_RNA <- function(model,cell,drug) {
  context <- model$context
  s <- sweep(sweep(sweep(getbaseline(cell),2,context$bmean),2,context$bsd,'/')%*%context$bv,2,context$bscale,'/')
  g <- match(drug,drugs)
  X <- matrix(0,length(cell),40)
  for(j in unique(g)) {
    ii <- which(g==j)
    X[ii,] <- sweep(s[ii,,drop=FALSE],2,context$bm[j,])%*%model$par$B[,,j]
  }
  list(X=X,D=X*context$sigma[g,,drop=FALSE]+context$mu[g,,drop=FALSE],g=g)
}
predict_single <- function(model,meta) {
  q <- predict_RNA(model,meta$depmap_id,meta$drug)
  context <- model$context
  node <- route_node(q$X,q$g,model$par,context$Q,context$mu,context$sigma)
  list(prediction=rowSums(node$activity*model$par$C[match(meta$assay,assays),,drop=FALSE]),
    activity=node$activity,D=q$D)
}
predict_pair <- function(model,meta) {
  aa <- predict_single(model,meta[,.(depmap_id,drug=drugA,assay=paste0('nci_',SCREENER))])
  bb <- predict_single(model,meta[,.(depmap_id,drug=drugB,assay=paste0('nci_',SCREENER))])
  J <- as.numeric(pair_features(aa$activity,bb$activity)%*%model$par$Gamma)
  list(total=aa$prediction+bb$prediction+J,departure=J)
}
# Verify all three held-cell models and molecular mean preservation by replay.
single <- fread(file.path(root,'data/single_predictions_frozen.tsv'))
measured <- fread(file.path(root,'data/measured_displacement_predictions_frozen.tsv'))
combo <- fread(file.path(root,'data/combination_predictions_frozen.tsv'))
checks <- list()
for(k in 1:3) {
  model <- qread(file.path(root,'models',paste0('fold_',k,'.qs')))
  ii <- which(single$evaluation=='held_cell' & single$fold==k)
  s <- single[ii]
  s[,assay:='secondary']
  replay <- predict_single(model,s)$prediction
  stopifnot(!any(s$depmap_id %in% model$training_single_cells),
    !any(s$depmap_id %in% model$context$training_RNA_cells))
  pi <- which(combo$evaluation=='held_cell' & combo$fold==k)
  p <- combo[pi]
  pr <- predict_pair(model,p)
  stopifnot(!any(p$depmap_id %in% model$training_pair_cells),
    !any(p$depmap_id %in% model$context$training_RNA_cells))
  context <- model$context
  zero <- -context$mu/context$sigma
  zn <- route_node(zero,1:10,model$par,context$Q,context$mu,context$sigma)
  zero_error <- max(abs(zn$activity))
  A <- matrix(c(.2,.4,.1,.3),2,2)
  B <- matrix(c(.1,.5,.3,.2),2,2)
  symmetry <- max(abs(pair_features(A,B)-pair_features(B,A)))
  training_m <- molecular[depmap_id %in% context$training_RNA_cells]
  rna <- predict_RNA(model,training_m$depmap_id,training_m$drug)$D
  means <- do.call(rbind,lapply(drugs,function(dr)colMeans(rna[training_m$drug==dr,,drop=FALSE])))
  mean_error <- max(abs(means-context$mu))
  si_error <- max(abs(replay-s$predicted_log2))
  co_error <- max(abs(pr$total-p$predicted_log2),abs(pr$departure-p$predicted_departure_log2))
  stopifnot(si_error<1e-10,co_error<1e-10,zero_error<1e-12,symmetry<1e-12,mean_error<1e-10)
  checks[[k]] <- data.table(fold=k,single_replay_difference=si_error,combination_replay_difference=co_error,
    zero_displacement_effect=zero_error,pair_symmetry_difference=symmetry,molecular_mean_difference=mean_error)
}
write_table(rbindlist(checks),'prediction_validation')

# Apply the fitted shared molecular maps to HGF, with no HGF functional fitting.
model <- qread(file.path(root,'models/final.qs'))
context <- model$context
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
hmean <- as.numeric(hgf[context=='Equal-context mean',..programs])
f <- (hmean-context$fm)/context$fs
qh <- cbind(1,sweep(matrix(f,nrow=1)%*%context$fv,2,context$fscale,'/'))
inferred_sigma <- as.numeric(exp(qh%*%context$scale_map))
used_sigma <- pmin(pmax(inferred_sigma,apply(context$sigma,2,min)),apply(context$sigma,2,max))
stopifnot(all(is.finite(used_sigma)),all(used_sigma>0))
Q <- rbind(context$Q,qh)
mu <- rbind(context$mu,hmean)
sigma <- rbind(context$sigma,used_sigma)
Xh <- sweep(sweep(as.matrix(hgf[,..programs]),2,hmean),2,used_sigma,'/')
hnode <- route_node(Xh,rep(11L,nrow(hgf)),model$par,Q,mu,sigma)
lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',assay='secondary'))
lact <- lap$activity[rep(1,nrow(hgf)),,drop=FALSE]
fh <- as.numeric(hnode$activity%*%model$par$C[1,])
J <- as.numeric(pair_features(lact,hnode$activity)%*%model$par$Gamma)
query <- data.table(dose=5,cell='AU565',hgf_context=hgf$context,
  lapatinib_log2=lap$prediction,hgf_log2=fh,combination_log2=lap$prediction+fh+J,
  lapatinib_pct_vehicle=100*2^lap$prediction,hgf_alone_pct_vehicle=100*2^fh,
  combination_pct_vehicle=100*2^(lap$prediction+fh+J),
  lapatinib_in_hgf_pct_hgf_control=100*2^(lap$prediction+J),interaction_log2=J,
  qualification='HGF receiver scale RNA-inferred, not measured; no HGF or Wilson functional labels fitted')
stopifnot(all(is.finite(query$interaction_log2)),!'ACH-000248' %in% model$training_single_cells,
  !'ACH-000248' %in% model$context$training_RNA_cells)
write_table(query,'au565_hgf_5uM_predictions_frozen')
write_table(data.table(program=programs,HGF_scale_inferred=inferred_sigma,HGF_scale_used=used_sigma,
  clipped=inferred_sigma!=used_sigma),'hgf_scale_transfer')
write_table(data.table(axis=1:3,HGF_descriptor=as.numeric(qh[,-1]),
  measured_drug_min=apply(context$Q[,-1,drop=FALSE],2,min),
  measured_drug_max=apply(context$Q[,-1,drop=FALSE],2,max)), 'hgf_footprint_domain')
paths <- file.path(root,'data',c('single_predictions_frozen.tsv','measured_displacement_predictions_frozen.tsv',
  'combination_predictions_frozen.tsv','au565_hgf_5uM_predictions_frozen.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')

# Evaluation uses Pearson at fixed drug/dose across excluded recipient cells.
metrics <- single[,.(cells=.N,r=pearson(observed_log2,predicted_log2),
  observed_sensitive=sum(observed_log2< -1),predicted_sensitive=sum(predicted_log2< -1)),by=.(drug,dose,evaluation)]
write_table(metrics,'single_correlations')
measured[,`:=`(yc=observed_log2-mean(observed_log2),pc=predicted_log2-mean(predicted_log2)),by=.(drug,evaluation,fold)]
write_table(measured[,.(cells=.N,r=pearson(observed_log2,predicted_log2),
  within_fold_r=pearson(yc,pc)),by=.(drug,dose,evaluation)],'measured_displacement_correlations')
pm <- combo[,.(cells=.N,r_total=pearson(observed_log2,predicted_log2),
  r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_departure_SD=sd(observed_departure_log2),predicted_departure_SD=sd(predicted_departure_log2)),
  by=.(pair,SCREENER,evaluation)]
write_table(pm,'combination_correlations')
old <- fread('/work/results/program_pilot/75_hgf_pattern_partner_transfer/data/protected_441_monotherapy_predictions.tsv')
common <- single[evaluation=='held_cell' & drug=='Lapatinib ditosylate' & depmap_id %in% old$depmap_id]
write_table(common[,.(cells=.N,r=pearson(observed_log2,predicted_log2))],'original_441_cell_lapatinib_5uM')
print(metrics[evaluation=='held_cell'])
print(pm[evaluation=='held_cell'])
print(query[hgf_context=='Equal-context mean'])

held <- single[evaluation=='held_cell']
held[,panel:=paste0(drug,' | ',dose,' µM')]
lab <- metrics[evaluation=='held_cell',.(panel=paste0(drug,' | ',dose,' µM'),label=sprintf('r = %.2f; n = %d',r,cells))]
p <- ggplot(held,aes(observed_log2,predicted_log2))+
  geom_point(aes(color=factor(fold)),size=.8,alpha=.6)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.45)+
  geom_label(data=lab,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=2.5,label.size=0)+
  facet_wrap(~panel,ncol=5,scales='free')+theme_bw(base_size=10)+
  theme(aspect.ratio=1,legend.position='bottom',strip.text=element_text(size=7))+
  labs(title='v93 joint calibrated sensitivity routes: held-cell PRISM predictions',
    x='Observed log2 abundance / vehicle',y='Predicted log2 abundance / vehicle',color='Held-cell fold',
    caption='Two learned nonlinear routes; per-drug RNA scale anchors; calibration, mixing, RNA bridge and interaction jointly trained.\nFixed selected doses; baseline reaches function only through predicted molecular programs.')
ggsave(file.path(root,'figures/held_cell_single_predictions.pdf'),p,width=13,height=7,device=cairo_pdf)
ggsave(file.path(root,'figures/held_cell_single_predictions.png'),p,width=13,height=7,dpi=120)
# Measured-RNA view localizes any failure in the displacement->function step.
mo <- measured[evaluation=='held_cell']
mo[,panel:=paste0(drug,' | ',dose,' µM')]
ml <- mo[,.(label=sprintf('r = %.2f; n = %d',pearson(observed_log2,predicted_log2),.N)),by=panel]
po <- ggplot(mo,aes(observed_log2,predicted_log2))+
  geom_point(aes(color=factor(fold)),size=1.4)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.45)+
  geom_label(data=ml,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=2.5,label.size=0)+
  facet_wrap(~panel,ncol=5,scales='free')+theme_bw(base_size=10)+
  theme(aspect.ratio=1,legend.position='bottom',strip.text=element_text(size=7))+
  labs(title='v93: held-cell functional prediction using measured displacements',
    x='Observed PRISM log2 abundance / vehicle',y='Predicted log2 abundance from measured RNA',color='Held-cell fold')
ggsave(file.path(root,'figures/held_cell_measured_displacement_predictions.pdf'),po,width=13,height=6.5,device=cairo_pdf)
ggsave(file.path(root,'figures/held_cell_measured_displacement_predictions.png'),po,width=13,height=6.5,dpi=120)
for(endpoint in c('total','departure')) {
  pc <- combo[evaluation=='held_cell']
  if(endpoint=='total') pc[,`:=`(observed=observed_log2,predicted=predicted_log2)] else
    pc[,`:=`(observed=observed_departure_log2,predicted=predicted_departure_log2)]
  labels <- pc[,.(label=sprintf('r = %.2f; n = %d',pearson(observed,predicted),.N)),by=.(pair,SCREENER)]
  pp <- ggplot(pc,aes(observed,predicted))+geom_point(aes(color=factor(fold)),size=1.5)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.45)+
    geom_label(data=labels,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=2.5,label.size=0)+
    facet_wrap(~SCREENER+pair,ncol=3,scales='free')+theme_bw(base_size=10)+
    theme(aspect.ratio=1,legend.position='bottom',strip.text=element_text(size=7))+
    labs(title=paste('v93 held-cell combination',endpoint),x=paste('Observed',endpoint,'log2 abundance'),
      y=paste('Predicted',endpoint,'log2 abundance'),color='Held-cell fold')
  ggsave(file.path(root,'figures',paste0('held_cell_combination_',endpoint,'.pdf')),pp,width=10,height=7,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0('held_cell_combination_',endpoint,'.png')),pp,width=10,height=7,dpi=120)
}
write_table(data.table(program=rep(programs,K),route=rep(1:K,each=40),weight=as.numeric(model$par$W)),'route_program_weights')
logits <- context$Q%*%model$par$mix
alpha <- exp(logits-apply(logits,1,max))
alpha <- alpha/rowSums(alpha)
write_table(data.table(drug=rep(drugs,K),route=rep(1:K,each=10),mixing_weight=as.numeric(alpha)),'drug_route_mixtures')
calibration <- data.table(drug=rep(drugs,40),program=rep(programs,each=10),
  measured_mean=as.numeric(context$mu),measured_SD=as.numeric(context$sigma),
  learned_gain=as.numeric(exp(context$Q%*%model$par$gain)),learned_shift=as.numeric(context$Q%*%model$par$shift))
write_table(calibration,'learned_drug_program_calibration')
# Display absolute and HGF-control-normalized queries as different endpoints.
mean_query <- query[hgf_context=='Equal-context mean']
qplot <- data.table(condition=factor(c('Lapatinib','HGF alone','Lapatinib + HGF'),
  levels=c('Lapatinib','HGF alone','Lapatinib + HGF')),
  percent=c(mean_query$lapatinib_pct_vehicle,mean_query$hgf_alone_pct_vehicle,mean_query$combination_pct_vehicle))
ph <- ggplot(qplot,aes(condition,percent))+geom_col(fill='#7373A5',width=.65)+
  geom_text(aes(label=sprintf('%.1f%%',percent)),vjust=-.3,size=3.5)+theme_bw(base_size=10)+
  scale_y_continuous(expand=expansion(mult=c(0,.12)))+
  labs(title='AU565: lapatinib 5 µM and transferred mean HGF',x=NULL,y='Predicted % untreated abundance',
    caption=paste0('Lapatinib in HGF, relative to HGF alone: ',round(mean_query$lapatinib_in_hgf_pct_hgf_control,1),'%.\n',
      'HGF receiver scale inferred, not measured; Wilson outcomes excluded. Wilson reference is at 1 µM.'))
ggsave(file.path(root,'figures/au565_hgf_query.pdf'),ph,width=5.5,height=4.5,device=cairo_pdf)
ggsave(file.path(root,'figures/au565_hgf_query.png'),ph,width=5.5,height=4.5,dpi=140)
source_paths <- c(archive_path,file.path(v87,'data/individual_training_outcomes.tsv'),
  file.path(v87,'data/combination_training_outcomes.tsv'),file.path(root,'data/selected_doses.tsv'),
  '/work/src/program_pilot/93_joint_calibrated_sensitivity_routes.R')
write_table(data.table(path=source_paths,sha256=vapply(source_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
writeLines(c('Fixed two-route joint mock; no held-out selection of architecture or hyperparameters.',
  'Nine v87 contexts plus crizotinib 5 uM. No extra drugs/doses or LINCS.',
  'Per-drug/program measured RNA means/SDs are anchors. Shared RNA-footprint maps learn small log-gain/shift adjustments and route mixtures.',
  'Baseline-conditioned molecular deviations, nonlinear routes, calibration, mixture weights and three symmetric interaction terms fit jointly.',
  'PRISM singles, matched NCI constituents, combination totals/departures and observed-RNA functional outcomes supervise the same objective.',
  'No baseline functional shortcut, free drug intercept or invented summed combination RNA. Zero displacement gives zero functional effect.',
  'Functional loss: 80% within-treatment variation, 20% absolute endpoint. Measured molecular means remain fixed.',
  'HGF mean transported from external contexts; receiver SD inferred by an RNA-only footprint->variance map and bounded to measured program SD ranges.',
  'No HGF/Wilson functional outcome fitted; 5 uM query is not dose-matched to Wilson 1 uM.',
  'Dose choices and reused outcomes are retrospective.',capture.output(sessionInfo())),file.path(root,'RUN_DETAILS.txt'))
writeLines('Complete: four converged fits, replayed predictions, held-cell evaluations and frozen HGF query.',file.path(root,'STATUS.txt'))
