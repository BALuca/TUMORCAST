#!/usr/bin/env Rscript
# Diagnose the saved v106 fit. No model changes, tuning or new training.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/106_lapatinib_supported_dose_interactions/linear'
programs <- sprintf('P%02d',1:40)
calibration_penalty <- .005
interaction_penalty <- .05
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
drugs <- bank_selection$drug
contexts <- fread(file.path(root,'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
# Replay the exact pure prediction functions, avoiding the pipeline runner.
expressions <- parse('/work/src/program_pilot/106_lapatinib_supported_dose_interactions.R')
needed <- c('context_key','context_index','predict_rna','signature_node','pair_features','predict_single','predict_pair','geometry')
for(e in expressions) if(is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% needed)eval(e)
model <- qread(file.path(root,'models/final.qs'))
pairs <- fread(file.path(root,'data/combination_training_outcomes.tsv'))
A <- predict_single(model,pairs[,.(depmap_id,drug=drugA,dose=doseA)])
B <- predict_single(model,pairs[,.(depmap_id,drug=drugB,dose=doseB)])
P <- pair_features(A$activity,B$activity)
F <- sweep(P,2,model$context$interaction_scale,'/')
assay <- exp(model$par$assay)[match(pairs$SCREENER,c('FF','FG'))]
X <- F*assay
s0 <- (A$sensitivity+B$sensitivity)*assay
J <- as.numeric(X %*% model$par$Gamma)
reported <- fread(file.path(root,'data/combination_predictions_frozen.tsv'))[evaluation=='fitted']
stopifnot(max(abs(J+reported$predicted_departure_log2))<1e-10)
pairs[,`:=`(activity_A=rowSums(A$activity),activity_B=rowSums(B$activity),
             route_A=A$route,route_B=B$route,product=rowSums(P),
             independent_prediction_log2=-s0,predicted_departure_log2=-J)]
fixed <- pairs[,.(cells=.N,route_A=unique(route_A),route_B=unique(route_B),
  observed_departure_mean=mean(I),predicted_departure_mean=mean(predicted_departure_log2),
  observed_departure_SD=sd(I),predicted_departure_SD=sd(predicted_departure_log2),
  product_mean=mean(product),product_SD=sd(product),product_relative_SD=sd(product)/abs(mean(product)),
  raw_product_r=cor(product,I),fixed_prediction_r=cor(predicted_departure_log2,I),
  independent_prediction_r=cor(independent_prediction_log2,y),
  activity_A_mean=mean(activity_A),activity_A_SD=sd(activity_A),
  activity_B_mean=mean(activity_B),activity_B_SD=sd(activity_B)),by=.(pair,SCREENER)]
held <- fread(file.path(root,'data/combination_predictions_frozen.tsv'))[evaluation=='held_cell',
 .(held_prediction_r=cor(observed_departure_log2,predicted_departure_log2),held_SD=sd(predicted_departure_log2)),by=.(pair,SCREENER)]
fixed <- merge(fixed,held,by=c('pair','SCREENER'))
fixed[,compression_fold:=observed_departure_SD/held_SD]
fwrite(fixed,file.path(root,'data/interaction_sign_and_compression_diagnostic.tsv'),sep='\t')
# Split each gradient into between-pair mean and within-pair variation terms.
geoP <- geometry(-pairs$y,pairs$pair)
geoI <- geometry(-pairs$I,pairs$pair)
rP <- s0+J+pairs$y
rI <- J+pairs$I
meanP <- ave(rP,pairs$pair,FUN=mean)
meanI <- ave(rI,pairs$pair,FUN=mean)
gPmean <- .4*as.numeric(crossprod(X,geoP$w*meanP))/geoP$within
gPwithin <- .4*as.numeric(crossprod(X,geoP$w*(rP-meanP)))/geoP$within
gImean <- .4*as.numeric(crossprod(X,geoI$w*meanI))/geoI$within
gIwithin <- .4*as.numeric(crossprod(X,geoI$w*(rI-meanI)))/geoI$within
gPenalty <- 2*interaction_penalty*model$par$Gamma
contributions <- data.table(term=c('route1_route1','cross_route','route2_route2'),
 coefficient_sensitivity=model$par$Gamma,total_mean_gradient=gPmean,total_within_gradient=gPwithin,
 departure_mean_gradient=gImean,departure_within_gradient=gIwithin,penalty_gradient=gPenalty,
 net_gradient=gPmean+gPwithin+gImean+gIwithin+gPenalty)
fwrite(contributions,file.path(root,'data/interaction_gradient_contributions.tsv'),sep='\t')
# Compare observed departure directions with model activities under fixed final fit.
print(fixed[,.(pair,product_relative_SD,raw_product_r,fixed_prediction_r,held_prediction_r,compression_fold)])
print(contributions)
# Hold the molecular activities fixed and inspect variation attributable solely
# to changing interaction coefficients/scales across LOO fits.
coeff <- rbindlist(lapply(seq_len(62),function(i) {
 m <- qread(file.path(root,'models',sprintf('loo_%02d.qs',i)))
 data.table(fold=i,assay=c('nci_FF','nci_FG'),
   route11=exp(m$par$assay)*(-m$par$Gamma[1])/m$context$interaction_scale[1],
   cross=exp(m$par$assay)*(-m$par$Gamma[2])/m$context$interaction_scale[2])
}))
fwrite(coeff,file.path(root,'data/held_interaction_effective_coefficients.tsv'),sep='\t')
# Isolate fold-coefficient effects while holding all molecular activities fixed.
oof <- fread(file.path(root,'data/combination_predictions_frozen.tsv'))[evaluation=='held_cell']
target <- pairs[,.(depmap_id,pair,SCREENER,product,I,route_A,route_B)]
target <- merge(target,oof[,.(depmap_id,pair,SCREENER,fold,predicted_departure_log2)],by=c('depmap_id','pair','SCREENER'))
target[,assay:=paste0('nci_',SCREENER)]
target <- merge(target,coeff,by=c('fold','assay'))
target[,coefficient_only_prediction:=product*ifelse(route_A==route_B,route11,cross)]
fold_effect <- target[,.(fixed_activity_changing_coefficient_r=cor(I,coefficient_only_prediction),
  actual_held_r=cor(I,predicted_departure_log2)),by=pair]
fwrite(fold_effect,file.path(root,'data/interaction_fold_coefficient_diagnostic.tsv'),sep='\t')
print(fold_effect)
# Quantify how much product variation lies within pairs, versus pair means.
variance <- rbindlist(lapply(which(model$context$interaction_supported),function(k) {
  xx <- F[,k]
  active <- abs(xx)>1e-12
  xx <- xx[active]; groups <- pairs$pair[active]
  centered <- xx-ave(xx,groups,FUN=mean)
  data.table(term=c('route1_route1','cross_route','route2_route2')[k],
    within_pair_fraction_of_centered_feature_variance=sum(centered^2)/sum((xx-mean(xx))^2))
}))
fwrite(variance,file.path(root,'data/interaction_feature_variation_diagnostic.tsv'),sep='\t')
print(variance)
