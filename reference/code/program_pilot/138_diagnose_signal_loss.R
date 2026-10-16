#!/usr/bin/env Rscript
# Diagnostic decomposition of the existing v138 model, not a new predictor.
x <- readLines('/work/src/program_pilot/138_pooled_direction_constrained_rescue.R')
cut <- which(grepl('^# Predictions are frozen',x))
eval(parse(text=x[seq_len(cut-1L)]))
model <- qread(file.path(root,'models/final.qs'))
f <- qread(file.path(root,'models/final_molecular_features.qs'))
out <- file.path(root,'data/signal_loss_diagnostics')
fig <- file.path(root,'figures/signal_loss_diagnostics')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
save_table <- function(x,n)fwrite(x,file.path(out,paste0(n,'.tsv')),sep='\t')

# Reconstruct exactly the fitted normal equations, including source-view weights.
inspect_fit <- fit_g
body(inspect_fit)[[length(body(inspect_fit))]] <- quote(list(model=model,z=z,scale=scale,
  gram=gram,rhs=rhs,Xc=Xc,Xm=Xm,variance_c=variance_c,variance_m=variance_m))
d <- inspect_fit(f,seq_len(nrow(external)),model$interaction_training$penalty)
stopifnot(max(abs(d$model$par$Gamma-model$par$Gamma))<1e-8)
mean_views <- function(X) {
  if(is.null(dim(X)))X <- matrix(X,ncol=1)
  ans <- rowsum(X,f$view_rows,reorder=FALSE)/as.numeric(table(f$view_rows))
  stopifnot(nrow(ans)==nrow(external))
  ans
}
A <- mean_views(f$Aviews[,1:40])
B <- mean_views(f$Bviews[,1:40])
P <- mean_views(f$products)
primary_B <- f$B[,1:40]
primary_P <- pair_features(f$A,f$B,model$context$W,basis)
ligand_key <- paste(external$depmap_id,external$gene)
raw_key <- paste(external_ligands$depmap_id,external_ligands$gene)
ix <- match(ligand_key,raw_key)
stopifnot(!anyNA(ix))
raw_B <- as.matrix(external_ligands[ix,..programs])
raw_H <- as.numeric(profiles[CellLine=='AU565' & source_profile=='Primary_study_mean',..programs])
ctx <- model$context
Dh <- matrix(raw_H,1)
Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
sm <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
hscale <- pmin(pmax(exp(Qh %*% sm),matrix(apply(ctx$sigma[1:10,],2,min),1)),
  matrix(apply(ctx$sigma[1:10,],2,max),1))
H <- as.numeric(Dh/hscale)
D <- copy(external)
D[,`:=`(raw_HGF_direction=as.numeric(raw_B %*% raw_H)/sqrt(sum(raw_H^2)),
  calibrated_HGF_direction=as.numeric(B %*% H)/sqrt(sum(H^2)),
  primary_HGF_direction=as.numeric(primary_B %*% H)/sqrt(sum(H^2)),
  reversal=P[,1],reinforcement=P[,2],primary_reversal=primary_P[,1],
  fitted_departure=predict_g(model,f))]
stages <- c('raw_HGF_direction','primary_HGF_direction','calibrated_HGF_direction',
  'primary_reversal','reversal','reinforcement','fitted_departure')
correlations <- rbindlist(lapply(c(FALSE,TRUE),function(exclude_HGF) {
  D[!exclude_HGF | gene!='HGF',lapply(.SD,function(z)pearson(z,observed_departure_log2)),
    by=.(CellLine,drug,dose),.SDcols=stages][,excluded_HGF:=exclude_HGF]
}))
save_table(correlations,'stage_correlations_with_observed_rescue')
summary <- melt(correlations,id.vars=c('CellLine','drug','dose','excluded_HGF'),variable.name='stage')[,
  .(median_r=median(value),min_r=min(value),max_r=max(value),positive_contexts=sum(value>0)),by=.(excluded_HGF,stage)]
save_table(summary,'stage_correlation_summary')
# Observed HGF rank and its rank at each representation stage within 13 ligands.
ranks <- D[,c(list(gene=gene),lapply(.SD,function(z)rank(-z,ties.method='average'))),
  by=.(CellLine,drug,dose),.SDcols=c('observed_departure_log2',stages)]
save_table(ranks[gene=='HGF'],'HGF_rank_at_each_stage')
save_table(D[,c('CellLine','drug','dose','gene','observed_departure_log2',stages),with=FALSE],
 'all_ligand_stage_scores')

# Locate the program contributions lost or amplified in the reversal statistic.
# R = sum_p mean_k(W_pk^2) * max(-A_p*B_p,0): a nonnegative magnitude gate.
weights <- rowMeans(ctx$W^2)
RV <- sweep(pmax(-f$Aviews[,1:40]*f$Bviews[,1:40],0),2,weights,'*')
stopifnot(max(abs(rowSums(RV)-f$products[,1]))<1e-12)
RP <- mean_views(RV)
program_tables <- list()
for(j in seq_len(40)) {
  program_tables[[j]] <- D[,.(CellLine,drug,dose,gene,observed_departure_log2)][,
    `:=`(program=programs[j],ligand_displacement=B[,j],drug_displacement=A[,j],
      reversal_contribution=RP[,j],signature_weight=weights[j])]
}
program_long <- rbindlist(program_tables)
program_assoc <- program_long[,.(r_ligand=pearson(ligand_displacement,observed_departure_log2),
  r_gated_reversal=pearson(reversal_contribution,observed_departure_log2)),by=.(CellLine,drug,dose,program)]
save_table(program_assoc,'per_program_before_after_reversal_gate')
contribution_summary <- program_long[,.(mean_reversal=mean(reversal_contribution)),by=.(gene,program)]
contribution_summary[,share:=mean_reversal/sum(mean_reversal),by=gene]
save_table(contribution_summary,'reversal_program_contributions_by_ligand')
# Which ligands the reversal score incorrectly ranks above HGF.
ligand_summary <- D[,.(mean_observed=mean(observed_departure_log2),mean_reversal=mean(reversal),
  mean_predicted=mean(fitted_departure)),by=gene]
setorder(ligand_summary,-mean_reversal)
save_table(ligand_summary,'ligand_ranking_by_reversal')

# Attribute shrinkage without changing the deployed fit: no-ridge sensitivity
# calculation, and removal of alternative-view uncertainty (one molecular view per observation).
no_penalty <- fit_g(f,seq_len(nrow(external)),0)
collapsed <- f
collapsed$products <- P
collapsed$view_rows <- seq_len(nrow(external))
collapsed_fit <- fit_g(collapsed,seq_len(nrow(external)),model$interaction_training$penalty)
# AU565 query features are exactly the previously saved f+HGF molecular inputs.
lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=5))
HP <- pair_features(lap$program_activity,cbind(matrix(H,1),matrix(H,1)%*%ctx$W),ctx$W,basis)
variants <- list(actual=model,no_ridge_diagnostic=no_penalty,mean_views_diagnostic=collapsed_fit)
checks <- rbindlist(lapply(names(variants),function(name) {
  m <- variants[[name]]
  g <- -sum(HP/m$context$interaction_scale*m$par$Gamma)
  data.table(diagnostic=name,reversal_coefficient=-m$par$Gamma[1],reinforcement_coefficient=m$par$Gamma[2],
    HGF_departure_log2=g,AU565_combination_pct=100*2^(-lap$sensitivity+g))
}))
save_table(checks,'regularization_and_molecular_uncertainty_diagnostics')
# Sufficient-statistic contributions to the fitted coefficient. These concern
# within-ligand ordering and group means, not an unconstrained baseline term.
save_table(data.table(term=terms,
  within_rhs=within_weight*as.numeric(crossprod(d$Xc,d$z$w*d$z$yc))/d$variance_c,
  mean_rhs=(1-within_weight)*as.numeric(crossprod(d$Xm,d$z$w*d$z$ym))/d$variance_m,
  fitted_sensitivity_coefficient=model$par$Gamma),'coefficient_evidence_by_objective_component')

# A direct diagnostic contrast: pooled mean across six cells for each inhibitor,
# 13 ligand points. These are fitted/descriptive, not cross-validation.
plot_data <- D[,lapply(.SD,mean),by=.(drug,gene),
 .SDcols=c('observed_departure_log2','raw_HGF_direction','calibrated_HGF_direction','reversal','fitted_departure')]
plot_data <- melt(plot_data,id.vars=c('drug','gene','observed_departure_log2'),variable.name='stage',value.name='score')
plot_data[,stage:=factor(stage,levels=c('raw_HGF_direction','calibrated_HGF_direction','reversal','fitted_departure'),
 labels=c('Raw ligand: HGF direction','Calibrated ligand: HGF direction','Pooled reversal score','Fitted g'))]
ann <- plot_data[,.(label=sprintf('r = %.2f',pearson(score,observed_departure_log2))),by=.(drug,stage)]
p <- ggplot(plot_data,aes(observed_departure_log2,score))+
 geom_point(aes(colour=gene=='HGF'),size=1.8)+scale_colour_manual(values=c('FALSE'='#555555','TRUE'='#D55E00'),guide='none')+
 geom_smooth(method='lm',formula=y~x,se=FALSE,colour='#0072B2',linewidth=.4)+
 geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=2.6)+
 facet_wrap(vars(drug,stage),ncol=4,scales='free')+theme_bw(base_size=8)+theme(aspect.ratio=1)+
 labs(x='Measured additional growth (log2)',y='Molecular score or fitted effect',
 title='Where the signal is lost | 13 ligands, averaged across six cells per inhibitor',
 caption='Orange = HGF. Descriptive/fitted comparisons, not held-out validation. Molecular scores have different units.')
ggsave(file.path(fig,'signal_loss_stages.pdf'),p,width=10,height=5.6,device=cairo_pdf)
ggsave(file.path(fig,'signal_loss_stages.png'),p,width=10,height=5.6,dpi=150)
print(summary)
print(ligand_summary)
print(checks)
print(ranks[gene=='HGF'])
