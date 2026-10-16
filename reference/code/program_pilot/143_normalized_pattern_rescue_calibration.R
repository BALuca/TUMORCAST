#!/usr/bin/env Rscript
# Fixed v106 single-drug function and seven program weights.
# Normalize the ligand's complete selected-program vector, preserving its shape.
# Primary calibration: all 12 NRG1 count departures, including negative values.
# g = drug_activity * (gain * relative_opposition - dose_specific_threshold).
# Three nonnegative coefficients: one gain and two thresholds. No penalty or tuning.
# EGF is an external control; an all-16-observation fit is a declared sensitivity check.
# HGF/Wilson outcomes never fit, select, or calibrate this model.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(glmnet);library(ggplot2);library(digest)})
RhpcBLASctl::blas_set_num_threads(1L)
setDTthreads(1L)
root <- '/work/results/program_pilot/143_normalized_pattern_rescue_calibration'
source_root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
programs <- sprintf('P%02d',1:40)
association_cutoff <- .5
for(folder in c('data','models','figures')) dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
pearson <- function(x,y) if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
model <- qread(file.path(source_root,'models/final.qs'))
reference <- qread('/work/results/program_pilot/106_lapatinib_supported_dose_interactions/linear/models/final.qs')
stopifnot(identical(model$par$level,reference$par$level),identical(model$par$log_amplitude,reference$par$log_amplitude),
          identical(model$context$W,reference$context$W),identical(model$rna$lap_receptors,reference$rna$lap_receptors))
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
drugs <- bank_selection$drug
contexts <- fread(file.path(source_root,'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
required <- c('context_key','context_index','predict_rna','signature_node','predict_single')
found <- character()
for(e in parse('/work/src/program_pilot/138_pooled_direction_constrained_rescue.R')) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% required) {
    eval(e)
    found <- c(found,as.character(e[[2]]))
  }
}
stopifnot(setequal(found,required))
outcomes <- fread(file.path(source_root,'data/external_combination_outcomes.tsv'))[study=='elife32271']
ligands <- fread(file.path(source_root,'data/external_ligand_programs.tsv'))
profiles <- fread(file.path(source_root,'data/time_matched_HGF_program_inputs.tsv'))
stopifnot(nrow(outcomes)==16,all(outcomes$gene %in% c('NRG1','EGF')),
          !any(outcomes$CellLine=='AU565'),all(outcomes$time_h==72),uniqueN(outcomes$depmap_id)==6)
outcomes[, observation_id := .I]
key <- function(x) paste(x$depmap_id,x$gene,sep='|')
stopifnot(!anyDuplicated(key(ligands)))
ix <- match(key(outcomes),key(ligands))
stopifnot(!anyNA(ix))
H <- as.matrix(ligands[ix,..programs])
single <- predict_single(model,outcomes[,.(depmap_id,drug,dose)])
association <- model$template$correlations[,'Lapatinib ditosylate']
keep <- abs(association)>=association_cutoff
stopifnot(sum(keep)==7)
write_table(data.table(program=programs,observed_single_drug_r=association,used_for_g=keep),
            'fixed_program_selection')

# Whole-vector normalization is performed after the frozen decoder transformation.
# The threshold acts on drug-sensitive molecular activity, not on baseline directly.
# Zero ligand has no direction and explicitly gives g=0. The normalized ligand
# representation describes nonzero response patterns, not a ligand-dose response.
make_features <- function(meta,L,H) {
  g <- context_index(meta$drug,meta$dose)
  route <- model$context$route[g]
  beta <- t(model$context$W[,route,drop=FALSE])/model$context$sigma[g,,drop=FALSE]
  beta[,!keep] <- 0
  d <- pmax(L*beta,0)
  h <- H*beta
  ligand_norm <- sqrt(rowSums(h^2))
  present <- ligand_norm>0
  normalized_h <- matrix(0,nrow(h),ncol(h),dimnames=dimnames(h))
  normalized_h[present,] <- h[present,,drop=FALSE]/ligand_norm[present]
  drug_activity <- rowSums(d)
  drug_activity[!present] <- 0
  opposed <- d*pmax(-normalized_h,0)
  opposition <- rowSums(opposed)
  relative_opposition <- ifelse(drug_activity>0,opposition/pmax(drug_activity,.Machine$double.xmin),0)
  X <- cbind(opposition=opposition,
    threshold_0.5=-drug_activity*(meta$dose==.5),threshold_5=-drug_activity*(meta$dose==5))
  stopifnot(all(meta$dose %in% c(.5,5)),all(relative_opposition>=0),all(relative_opposition<=1+1e-10))
  stopifnot(all(is.finite(X)))
  list(X=X,opposed=opposed,drug_activity=drug_activity,ligand_norm=ligand_norm,
       relative_opposition=relative_opposition,directed_drug=d,normalized_ligand=normalized_h)
}
features <- make_features(outcomes,single$D,H)
X <- features$X
y <- outcomes$observed_departure_log2
stopifnot(all(is.finite(y)),all(X[,1]>=0),all(X[,2:3]<=0))
primary_rows <- which(outcomes$gene=='NRG1')
control_rows <- which(outcomes$gene=='EGF')
stopifnot(length(primary_rows)==12,length(control_rows)==4)
write_table(cbind(outcomes,as.data.table(X),data.table(ligand_norm=features$ligand_norm,
  drug_activity=features$drug_activity,relative_opposition=features$relative_opposition,
  role=ifelse(outcomes$gene=='NRG1','primary calibration','external control'))),'training_features_and_outcomes')

# Fit absolute 72h departures. Equal weight per ligand/dose group. RMS scaling
# is numerical conditioning only: no penalty and no response standardization.
fit_g <- function(rows) {
  group <- paste(outcomes$gene[rows],outcomes$dose[rows],sep='|')
  weight <- 1/as.numeric(table(group)[group])/length(unique(group))
  scale <- sqrt(colSums(X[rows,,drop=FALSE]^2*weight))
  stopifnot(all(scale>1e-10))
  Z <- sweep(X[rows,,drop=FALSE],2,scale,'/')
  A <- Z*sqrt(weight)
  target <- y[rows]*sqrt(weight)
  stopifnot(qr(A)$rank==ncol(A))
  fit <- glmnet(A,target,alpha=1,lambda=0,intercept=FALSE,standardize=FALSE,
                lower.limits=0,thresh=1e-14,maxit=100000)
  coef <- as.numeric(coef(fit))[-1]
  gradient <- as.numeric(crossprod(Z,weight*(as.numeric(Z %*% coef)-y[rows])))
  active <- coef>1e-7
  kkt <- max(c(abs(gradient[active]),pmax(-gradient[!active],0),0))
  stopifnot(kkt<1e-5,all(coef>=0))
  list(coef=setNames(coef/scale,colnames(X)),kkt=kkt,training_rows=rows,
       condition_number=kappa(A),training_cells=unique(outcomes$depmap_id[rows]))
}
rows <- seq_len(nrow(outcomes))
final <- fit_g(primary_rows)
predictions <- list()
fits <- list(final=final)
for(cell in unique(outcomes$depmap_id)) {
  valid <- primary_rows[outcomes$depmap_id[primary_rows]==cell]
  fit <- fit_g(setdiff(primary_rows,valid))
  fits[[paste0('held_cell_',cell)]] <- fit
  predictions[[cell]] <- outcomes[valid][,`:=`(
    predicted_departure_log2=as.numeric(X[valid,,drop=FALSE] %*% fit$coef),evaluation='held_cell')]
}
predictions[['EGF']] <- outcomes[control_rows][,`:=`(
  predicted_departure_log2=as.numeric(X[control_rows,,drop=FALSE] %*% final$coef),evaluation='external_EGF')]
fits[['all_16_sensitivity']] <- fit_g(rows)
predictions[['fitted_NRG1']] <- outcomes[primary_rows][,`:=`(
  predicted_departure_log2=as.numeric(X[primary_rows,,drop=FALSE] %*% final$coef),evaluation='fitted')]
pred <- rbindlist(predictions)
write_table(pred,'combination_predictions_frozen')
qsave(list(single_drug_parent=source_root,g=final,folds=fits,programs=programs[keep],
           calibration='12 NRG1 observations; six recipients; negative departures retained',
           formula='g = gain * normalized_opposition - threshold[dose] * drug_activity',
           endpoint='72h log2 count departure',time_projection='5/3 from 72h to five days'),
      file.path(root,'models/final.qs'))
write_table(rbindlist(lapply(names(fits),function(label) data.table(fit=label,
  term=names(fits[[label]]$coef),coef=unname(fits[[label]]$coef),kkt=fits[[label]]$kkt,
  condition_number=fits[[label]]$condition_number))),
  'fitted_interaction_gains')
metrics <- pred[,.(cells=.N,r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2),
  observed_SD=sd(observed_departure_log2),predicted_SD=sd(predicted_departure_log2),
  observed_min=min(observed_departure_log2),observed_max=max(observed_departure_log2),
  predicted_min=min(predicted_departure_log2),predicted_max=max(predicted_departure_log2),
  observed_positive=sum(observed_departure_log2>0),predicted_positive=sum(predicted_departure_log2>0)),
  by=.(evaluation,gene,dose)]
write_table(metrics,'combination_validation_correlations')

# Query HGF only after the coefficients and selection are frozen. Five-day
# projection assumes constant relative growth-rate differences across assay times.
query <- list()
contributions <- list()
calibration_queries <- list()
molecular_support <- list()
Ha <- profiles[CellLine=='AU565']
for(dose in c(.5,5)) {
  meta <- data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=dose)
  lap <- predict_single(model,meta)
  repeated <- meta[rep(1,nrow(Ha))]
  hh <- as.matrix(Ha[,..programs])
  ll <- lap$D[rep(1,nrow(Ha)),,drop=FALSE]
  ff <- make_features(repeated,ll,hh)
  train_at_dose <- primary_rows[outcomes$dose[primary_rows]==dose]
  molecular_support[[as.character(dose)]] <- data.table(dose=dose,source_profile=Ha$source_profile,
    normalized_opposition=ff$X[,1],drug_activity=ff$drug_activity,relative_opposition=ff$relative_opposition,
    training_opposition_min=min(X[train_at_dose,1]),training_opposition_max=max(X[train_at_dose,1]),
    training_relative_min=min(features$relative_opposition[train_at_dose]),
    training_relative_max=max(features$relative_opposition[train_at_dose]),
    training_drug_activity_min=min(features$drug_activity[train_at_dose]),
    training_drug_activity_max=max(features$drug_activity[train_at_dose]))
  departure_72h <- as.numeric(ff$X %*% final$coef)
  departure_5day <- departure_72h*5/3
  query[[as.character(dose)]] <- data.table(CellLine='AU565',dose=dose,source_profile=Ha$source_profile,
    lapatinib_pct_control=100*2^(-lap$sensitivity),g_72h_log2=departure_72h,g_5day_log2=departure_5day,
    predicted_combination_pct_control=100*2^(-lap$sensitivity+departure_5day),
    time_qualification='72h count interaction projected to five-day PRISM assuming constant relative growth-rate difference',
    HGF_alone_qualification='Standalone HGF effect not learned; interaction shown with neutral HGF standalone reference')
  threshold <- final$coef[paste0('threshold_',dose)]
  by_program <- ff$opposed*final$coef[1]-ff$directed_drug*threshold
  by_program[ff$ligand_norm==0,] <- 0
  stopifnot(max(abs(rowSums(by_program)-departure_72h))<1e-8)
  contributions[[as.character(dose)]] <- cbind(data.table(dose=dose,source_profile=Ha$source_profile),
    as.data.table(by_program))
  for(label in names(fits)) {
    gg <- as.numeric(ff$X %*% fits[[label]]$coef)
    calibration_queries[[paste(dose,label)]] <- data.table(fit=label,dose=dose,source_profile=Ha$source_profile,
      g_72h_log2=gg,predicted_combination_pct_control=100*2^(-lap$sensitivity+gg*5/3))
  }
}
answer <- rbindlist(query)
write_table(answer,'AU565_HGF_predictions_frozen')
write_table(rbindlist(contributions),'AU565_HGF_program_contributions')
write_table(rbindlist(calibration_queries),'AU565_HGF_calibration_sensitivity')
write_table(rbindlist(molecular_support),'AU565_HGF_molecular_support')
stopifnot(all(make_features(outcomes,single$D*0,H)$X==0),
          all(make_features(outcomes,single$D,H*0)$X==0))
# Whole-vector scaling must leave shape-based features unchanged at nonzero input.
stopifnot(max(abs(make_features(outcomes,single$D,H*10)$X-X))<1e-10)
write_table(data.table(check=c('zero_drug_zero_g','zero_ligand_zero_g','ligand_global_scale_invariance',
  'single_drug_coefficients_unchanged','HGF_Wilson_excluded_from_calibration'),passed=TRUE),'structural_validation')
prediction_files <- file.path(root,'data',c('combination_predictions_frozen.tsv','AU565_HGF_predictions_frozen.tsv'))
write_table(data.table(path=prediction_files,sha256=vapply(prediction_files,digest,character(1),file=TRUE,algo='sha256')),
  'prediction_hashes')

# One compact figure reports the prediction and the measured-count validation.
plot_data <- answer[source_profile=='Primary_study_mean',.(dose,lapatinib=lapatinib_pct_control,
                                                         `Lapatinib + HGF`=predicted_combination_pct_control)]
plot_data <- melt(plot_data,id.vars='dose',variable.name='condition',value.name='pct')
p <- ggplot(plot_data,aes(condition,pct,fill=condition))+geom_col(width=.65)+
  facet_wrap(~dose,nrow=1,scales='free_y',labeller=labeller(dose=function(x)paste(x,'uM')))+
  scale_fill_manual(values=c('#2769A5','#D86A25'))+
  labs(x=NULL,y='Predicted abundance (% control)',subtitle='Unbounded five-day projection; HGF outcomes excluded')+
  theme_classic(base_size=9)+theme(legend.position='none',axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(root,'figures/AU565_HGF_prediction.pdf'),p,width=4.3,height=2.5)
ggsave(file.path(root,'figures/AU565_HGF_prediction.png'),p,width=4.3,height=2.5,dpi=180)
annotations <- metrics[evaluation=='held_cell' & gene=='NRG1'][,
  .(gene,dose,label=sprintf('r = %.2f; n = %d',r_departure,cells))]
s <- ggplot(pred[evaluation=='held_cell' & gene=='NRG1'],aes(observed_departure_log2,predicted_departure_log2))+
  geom_point(size=1.7,color='#2769A5')+geom_smooth(method='lm',formula=y~x,se=FALSE,linewidth=.4,color='black')+
  facet_wrap(~dose,scales='free',labeller=labeller(dose=function(x)paste(x,'uM')))+
  geom_text(data=annotations,aes(x=-Inf,y=Inf,label=label),
    inherit.aes=FALSE,hjust=-.05,vjust=1.4,size=2.8)+
  labs(x='Measured 72h log2 count departure',y='Held-cell predicted departure')+
  theme_classic(base_size=9)+theme(aspect.ratio=1)
ggsave(file.path(root,'figures/NRG1_held_cell_validation.pdf'),s,width=4.4,height=2.5)
ggsave(file.path(root,'figures/NRG1_held_cell_validation.png'),s,width=4.4,height=2.5,dpi=220)
# Compare recipient-held predictions on the same 12 NRG1 endpoints. Training sets
# differ: v142 also calibrated on EGF; v143 explicitly calibrates on NRG1 alone.
old_predictions <- fread('/work/results/program_pilot/142_HER2_anchored_program_rescue/data/combination_predictions_frozen.tsv')
old_predictions <- old_predictions[evaluation=='held_cell' & gene=='NRG1']
new_predictions <- pred[evaluation=='held_cell' & gene=='NRG1']
check <- merge(old_predictions[,.(CellLine,dose,observed_departure_log2)],
  new_predictions[,.(CellLine,dose,observed_departure_log2)],by=c('CellLine','dose'))
stopifnot(nrow(check)==12,max(abs(check$observed_departure_log2.x-check$observed_departure_log2.y))<1e-12)
old_predictions[, version := 'v142: raw pattern; NRG1 + EGF calibration']
new_predictions[, version := 'v143: normalized pattern + threshold; NRG1 calibration']
comparison_predictions <- rbindlist(list(old_predictions,new_predictions),use.names=TRUE,fill=TRUE)
comparison_metrics <- comparison_predictions[,.(cells=.N,
  r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_SD=sd(observed_departure_log2),predicted_SD=sd(predicted_departure_log2),
  observed_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2)),by=.(version,dose)]
write_table(comparison_metrics,'comparison_v142_v143_held_recipients')
wilson_path <- '/work/almanac_combination_training_repair_2026-10-03/results/wilson_observed_source.tsv'
wilson <- fread(wilson_path)[series=='Wilson et al.' & condition %in% c('Lapatinib','Lapatinib + HGF')]
stopifnot(nrow(wilson)==2,all(is.finite(wilson$abundance)),all(is.finite(wilson$sd)))
qa <- answer[dose==5 & source_profile=='Primary_study_mean']
comparison <- rbind(
  data.table(panel='Prediction\n5 uM, five-day projection',condition=c('Lapatinib','Lapatinib + HGF'),
    abundance=c(qa$lapatinib_pct_control,qa$predicted_combination_pct_control),sd=NA_real_),
  wilson[,.(panel='Wilson measurement\n1 uM, 12 days',condition,abundance,sd)])
comparison[, condition := factor(condition,levels=c('Lapatinib','Lapatinib + HGF'))]
write_table(comparison,'prediction_beside_Wilson')
cp <- ggplot(comparison,aes(condition,abundance,fill=condition))+geom_col(width=.65)+
  geom_errorbar(data=comparison[is.finite(sd)],aes(ymin=abundance-sd,ymax=abundance+sd),width=.15)+
  facet_wrap(~panel,scales='free_y')+scale_fill_manual(values=c('#2769A5','#D86A25'))+
  labs(x=NULL,y='Abundance (% control)')+theme_classic(base_size=8)+
  theme(legend.position='none',axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(root,'figures/prediction_beside_Wilson.pdf'),cp,width=3.6,height=2.3)
ggsave(file.path(root,'figures/prediction_beside_Wilson.png'),cp,width=3.6,height=2.3,dpi=220)
source_paths <- c('/work/src/program_pilot/143_normalized_pattern_rescue_calibration.R',
  '/work/src/program_pilot/138_pooled_direction_constrained_rescue.R',
  file.path(source_root,'models/final.qs'),file.path(source_root,'data/external_combination_outcomes.tsv'),
  file.path(source_root,'data/external_ligand_programs.tsv'),file.path(source_root,'data/time_matched_HGF_program_inputs.tsv'),
  wilson_path)
write_table(data.table(path=source_paths,sha256=vapply(source_paths,digest,character(1),file=TRUE,algo='sha256')),
  'source_hashes')
write_table(data.table(item=c('primary_calibration','external_control','sensitivity_check','learned_parameters',
  'normalization','evaluation_scope','time_projection','scope_limit'),value=c(
  'All 12 NRG1 observations; 6 cells; both doses; negative departures included',
  'Four EGF observations; no EGF labels used in primary model',
  'All 16 NRG1/EGF observations fitted separately; never used to select primary HGF result',
  'One nonnegative gain and two nonnegative dose thresholds; unpenalized constrained least squares',
  'One Euclidean norm over seven frozen weighted program coordinates, not a per-program rescaling',
  'Recipient-held combination calibration with frozen single-drug molecular maps; retrospective development',
  '72h g multiplied by 5/3 for five-day query; constant relative growth-rate assumption',
  'Nonzero ligand response-pattern transfer; does not model ligand dose or independently predict HGF-alone effects')),
  'model_contract')
writeLines(capture.output(sessionInfo()),file.path(root,'data/R_session_info.txt'))
print(answer[source_profile=='Primary_study_mean'])
print(metrics)
print(final$coef)
