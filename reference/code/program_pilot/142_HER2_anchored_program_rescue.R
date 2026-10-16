#!/usr/bin/env Rscript
# Keep v106 single-drug function fixed. Fit only two directional interaction gains
# to 16 measured 72h lapatinib-NRG1/EGF count departures; HGF/Wilson never fitted.
# Programs are fixed from strong single-drug associations, not combination outcomes.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(glmnet);library(ggplot2);library(digest)})
RhpcBLASctl::blas_set_num_threads(1L)
setDTthreads(1L)
root <- '/work/results/program_pilot/142_HER2_anchored_program_rescue'
source_root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
programs <- sprintf('P%02d',1:40)
penalty_grid <- c(0,.01,.1,1,10)
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

# A drug's positively directed sensitivity contributions identify the programs
# on which a ligand can act. Opposed/aligned ligand moves remain separate features.
# The sensitivity weights and program SD are frozen from the single-drug decoder.
make_features <- function(meta,L,H) {
  g <- context_index(meta$drug,meta$dose)
  route <- model$context$route[g]
  beta <- t(model$context$W[,route,drop=FALSE])/model$context$sigma[g,,drop=FALSE]
  beta[,!keep] <- 0
  d <- pmax(L*beta,0)
  h <- H*beta
  opposed <- d*pmax(-h,0)
  aligned <- d*pmax(h,0)
  X <- cbind(opposition=rowSums(opposed),reinforcement=-rowSums(aligned))
  stopifnot(all(is.finite(X)))
  list(X=X,opposed=opposed,aligned=aligned)
}
features <- make_features(outcomes,single$D,H)
X <- features$X
y <- outcomes$observed_departure_log2
stopifnot(all(is.finite(y)),all(X[,1]>=0),all(X[,2]<=0))
write_table(cbind(outcomes,as.data.table(X)),'training_features_and_outcomes')

# Equal total weight per ligand/dose group. Only two nonnegative coefficients;
# the signed columns preserve opposition => rescue and alignment => reinforcement.
fit_g <- function(rows,penalty) {
  group <- paste(outcomes$gene[rows],outcomes$dose[rows],sep='|')
  weight <- 1/as.numeric(table(group)[group])/length(unique(group))
  scale <- sqrt(colSums(X[rows,,drop=FALSE]^2*weight))
  stopifnot(all(scale>1e-10))
  Z <- sweep(X[rows,,drop=FALSE],2,scale,'/')
  A <- rbind(Z*sqrt(weight),diag(sqrt(penalty),2))
  target <- c(y[rows]*sqrt(weight),0,0)
  # lambda=0 on augmented rows gives the explicit constrained ridge objective.
  fit <- glmnet(A,target,alpha=1,lambda=0,intercept=FALSE,standardize=FALSE,
                lower.limits=0,thresh=1e-14,maxit=100000)
  coef <- as.numeric(coef(fit))[-1]
  gradient <- as.numeric(crossprod(Z,weight*(as.numeric(Z %*% coef)-y[rows])))+penalty*coef
  active <- coef>1e-7
  kkt <- max(c(abs(gradient[active]),pmax(-gradient[!active],0),0))
  stopifnot(kkt<1e-5,all(coef>=0))
  list(coef=setNames(coef/scale,colnames(X)),penalty=penalty,kkt=kkt,training_rows=rows)
}
select_g <- function(rows) {
  cells <- unique(outcomes$depmap_id[rows])
  pred <- matrix(NA_real_,length(rows),length(penalty_grid))
  for(cell in cells) {
    valid <- which(outcomes$depmap_id[rows]==cell)
    train <- rows[-valid]
    for(j in seq_along(penalty_grid)) {
      fit <- fit_g(train,penalty_grid[j])
      pred[valid,j] <- as.numeric(X[rows[valid],,drop=FALSE] %*% fit$coef)
    }
  }
  group <- paste(outcomes$gene[rows],outcomes$dose[rows],sep='|')
  weight <- 1/as.numeric(table(group)[group])/length(unique(group))
  score <- colSums((pred-y[rows])^2*weight)
  fit <- fit_g(rows,penalty_grid[which.min(score)])
  fit$selection <- data.table(penalty=penalty_grid,internal_objective=score,
                              selected=seq_along(score)==which.min(score))
  fit
}
rows <- seq_len(nrow(outcomes))
final <- select_g(rows)
predictions <- list()
fits <- list(final=final)
for(cell in unique(outcomes$depmap_id)) {
  valid <- which(outcomes$depmap_id==cell)
  fit <- select_g(setdiff(rows,valid))
  fits[[paste0('held_cell_',cell)]] <- fit
  predictions[[cell]] <- outcomes[valid][,`:=`(
    predicted_departure_log2=as.numeric(X[valid,,drop=FALSE] %*% fit$coef),evaluation='held_cell')]
}
for(gene in unique(outcomes$gene)) {
  valid <- which(outcomes$gene==gene)
  fit <- select_g(setdiff(rows,valid))
  fits[[paste0('held_ligand_',gene)]] <- fit
  predictions[[gene]] <- outcomes[valid][,`:=`(
    predicted_departure_log2=as.numeric(X[valid,,drop=FALSE] %*% fit$coef),evaluation='held_ligand')]
}
pred <- rbindlist(predictions)
write_table(pred,'combination_predictions_frozen')
qsave(list(single_drug_parent=source_root,g=final,programs=programs[keep],
           endpoint='72h log2 count departure',time_projection='5/3 from 72h to five days'),
      file.path(root,'models/final.qs'))
write_table(rbindlist(lapply(names(fits),function(label) data.table(fit=label,
  term=names(fits[[label]]$coef),coef=unname(fits[[label]]$coef),penalty=fits[[label]]$penalty,kkt=fits[[label]]$kkt))),
  'fitted_interaction_gains')
write_table(final$selection,'final_training_only_penalty_selection')
metrics <- pred[,.(cells=.N,r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2)),
  by=.(evaluation,gene,dose)]
write_table(metrics,'combination_validation_correlations')

# Query HGF only after the coefficients and selection are frozen. Five-day
# projection assumes constant relative growth-rate differences across assay times.
query <- list()
contributions <- list()
Ha <- profiles[CellLine=='AU565']
for(dose in c(.5,5)) {
  meta <- data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=dose)
  lap <- predict_single(model,meta)
  repeated <- meta[rep(1,nrow(Ha))]
  hh <- as.matrix(Ha[,..programs])
  ll <- lap$D[rep(1,nrow(Ha)),,drop=FALSE]
  ff <- make_features(repeated,ll,hh)
  departure_72h <- as.numeric(ff$X %*% final$coef)
  departure_5day <- departure_72h*5/3
  query[[as.character(dose)]] <- data.table(CellLine='AU565',dose=dose,source_profile=Ha$source_profile,
    lapatinib_pct_control=100*2^(-lap$sensitivity),g_72h_log2=departure_72h,g_5day_log2=departure_5day,
    predicted_combination_pct_control=100*2^(-lap$sensitivity+departure_5day),
    time_qualification='72h count interaction projected to five-day PRISM assuming constant relative growth-rate difference',
    HGF_alone_qualification='Standalone HGF effect not learned; interaction shown with neutral HGF standalone reference')
  by_program <- ff$opposed*final$coef[1]-ff$aligned*final$coef[2]
  stopifnot(max(abs(rowSums(by_program)-departure_72h))<1e-8)
  contributions[[as.character(dose)]] <- cbind(data.table(dose=dose,source_profile=Ha$source_profile),
    as.data.table(by_program))
}
answer <- rbindlist(query)
write_table(answer,'AU565_HGF_predictions_frozen')
write_table(rbindlist(contributions),'AU565_HGF_program_contributions')
stopifnot(all(make_features(outcomes,single$D*0,H)$X==0),
          all(make_features(outcomes,single$D,H*0)$X==0))

# One compact figure reports the prediction and the measured-count validation.
plot_data <- answer[source_profile=='Primary_study_mean',.(dose,lapatinib=lapatinib_pct_control,
                                                         `Lapatinib + HGF`=predicted_combination_pct_control)]
plot_data <- melt(plot_data,id.vars='dose',variable.name='condition',value.name='pct')
p <- ggplot(plot_data,aes(condition,pct,fill=condition))+geom_col(width=.65)+
  facet_wrap(~dose,nrow=1,labeller=labeller(dose=function(x)paste(x,'uM')))+
  scale_fill_manual(values=c('#2769A5','#D86A25'))+
  labs(x=NULL,y='Predicted abundance (% control)',subtitle='Five-day projection; HGF outcomes excluded')+
  theme_classic(base_size=9)+theme(legend.position='none',axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(root,'figures/AU565_HGF_prediction.pdf'),p,width=4.3,height=2.5)
ggsave(file.path(root,'figures/AU565_HGF_prediction.png'),p,width=4.3,height=2.5,dpi=180)
annotations <- metrics[evaluation=='held_cell' & gene=='NRG1'][,
  .(gene,dose,label=sprintf('r = %.2f; n = %d',r_departure,cells))]
s <- ggplot(pred[evaluation=='held_cell' & gene=='NRG1'],aes(observed_departure_log2,predicted_departure_log2))+
  geom_point(size=1.7,color='#2769A5')+geom_smooth(method='lm',formula=y~x,se=FALSE,linewidth=.4,color='black')+
  facet_wrap(~dose,scales='free')+geom_text(data=annotations,aes(x=-Inf,y=Inf,label=label),
    inherit.aes=FALSE,hjust=-.05,vjust=1.4,size=2.8)+
  labs(x='Measured 72h log2 count departure',y='Held-cell predicted departure')+
  theme_classic(base_size=9)+theme(aspect.ratio=1)
ggsave(file.path(root,'figures/NRG1_held_cell_validation.pdf'),s,width=4.4,height=2.5)
annotations_ligand <- metrics[evaluation=='held_ligand' & gene=='NRG1'][,
  .(gene,dose,label=sprintf('r = %.2f; n = %d',r_departure,cells))]
sl <- ggplot(pred[evaluation=='held_ligand' & gene=='NRG1'],
  aes(observed_departure_log2,predicted_departure_log2))+
  geom_point(size=1.7,color='#2769A5')+
  geom_smooth(method='lm',formula=y~x,se=FALSE,linewidth=.4,color='black')+
  facet_wrap(~dose,scales='free')+
  geom_text(data=annotations_ligand,aes(x=-Inf,y=Inf,label=label),
    inherit.aes=FALSE,hjust=-.05,vjust=1.4,size=2.8)+
  labs(x='Measured 72h log2 count departure',y='Predicted departure',
    subtitle='All NRG1 outcomes excluded; six breast cell lines')+
  theme_classic(base_size=9)+theme(aspect.ratio=1)
ggsave(file.path(root,'figures/NRG1_held_ligand_validation.pdf'),sl,width=4.4,height=2.7)
ggsave(file.path(root,'figures/NRG1_held_ligand_validation.png'),sl,width=4.4,height=2.7,dpi=220)
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
write_table(data.table(path=c('/work/src/program_pilot/142_HER2_anchored_program_rescue.R',
  file.path(source_root,'models/final.qs'),file.path(source_root,'data/external_combination_outcomes.tsv'),
  file.path(source_root,'data/external_ligand_programs.tsv')),sha256=vapply(c(
  '/work/src/program_pilot/142_HER2_anchored_program_rescue.R',file.path(source_root,'models/final.qs'),
  file.path(source_root,'data/external_combination_outcomes.tsv'),file.path(source_root,'data/external_ligand_programs.tsv')),
    digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
writeLines(capture.output(sessionInfo()),file.path(root,'data/R_session_info.txt'))
print(answer[source_profile=='Primary_study_mean'])
print(metrics)
print(final$coef)
