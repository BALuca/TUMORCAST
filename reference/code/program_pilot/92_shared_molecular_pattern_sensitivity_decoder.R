#!/usr/bin/env Rscript
# Measured displacement decoder: molecular patterns gate 1–3 regularized rules.
# Combination departures use the same learned sensitivity signatures. No PLS.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(glmnet);library(qs);library(ggplot2);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
root <- '/work/results/program_pilot/92_shared_molecular_pattern_sensitivity_decoder'
for(d in c('data','models','figures')) dir.create(file.path(root,d),recursive=TRUE,showWarnings=FALSE)
v87 <- '/work/results/program_pilot/87_selected_dose_joint_program_bridge'
programs <- sprintf('P%02d',1:40)
rank <- 6L
lambda_grid <- 10^seq(0,-2.5,length.out=14)
seed <- 92012L
write_table <- function(x,n) fwrite(x,file.path(root,'data',paste0(n,'.tsv')),sep='\t')
pearson <- function(x,y) if(length(x)>=4 && sd(x)>1e-12 && sd(y)>1e-12) cor(x,y) else NA_real_
status <- function(x) {cat(format(Sys.time()),x,'\n');flush.console();writeLines(x,file.path(root,'STATUS.txt'))}
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
selected <- rbind(fread(file.path(v87,'data/selected_doses.tsv'))[,.(drug,dose)],data.table(drug='crizotinib',dose=5))
molecular <- merge(as.data.table(archive$prep$m),selected,by=c('drug','dose'))[depmap_id!='ACH-000248']
baseline <- as.data.table(archive$prep$raw)
# PRISM target interpolation is inherited exactly from v87 and the crizotinib audit.
functional <- fread(file.path(v87,'data/individual_training_outcomes.tsv'))[source=='prism_secondary',.(depmap_id,drug,dose,y)]
criz <- fread('/work/results/program_pilot/91_lapatinib_crizotinib_5uM_joint_bridge/data/individual_training_outcomes.tsv')[drug=='crizotinib',.(depmap_id,drug,dose,y)]
functional <- rbind(functional,criz)
profiles <- merge(molecular,functional,by=c('depmap_id','drug','dose'))
setorder(profiles,drug,depmap_id)
folds <- fread(file.path(v87,'data/recipient_folds.tsv'))
profiles <- merge(profiles,folds[,.(depmap_id,fold)],by='depmap_id',all.x=TRUE)
setorder(profiles,drug,depmap_id)
pairs <- fread(file.path(v87,'data/combination_training_outcomes.tsv'))
pairs <- merge(pairs,folds[,.(depmap_id,fold)],by='depmap_id',all.x=TRUE)
setorder(pairs,row_id)
stopifnot(nrow(profiles)==286L,uniqueN(profiles$depmap_id)==30L,uniqueN(profiles$drug)==10L,
  nrow(pairs)==193L,uniqueN(pairs$depmap_id)==39L,!anyNA(profiles$fold),!anyNA(pairs$fold),
  !anyDuplicated(profiles[,.(depmap_id,drug,dose)]),!anyDuplicated(pairs$row_id),
  all(is.finite(as.matrix(profiles[,..programs]))),all(abs(pairs$I-(pairs$y-pairs$yA-pairs$yB))<1e-10))
write_table(selected,'selected_doses')
write_table(profiles[,.(depmap_id,drug,dose,fold,y)],'measured_profile_outcomes')
write_table(pairs,'combination_outcomes')
status(paste('Inputs:',nrow(profiles),'measured PRISM profiles;',nrow(pairs),'combination outcomes.'))

# RNA alone determines six axes of coordinated WITHIN-drug movement. Full signed
# scores (including treatment context) are projected on those axes. No viability
# label selects coordinates or molecular cluster centers.
fit_representation <- function(excluded,tag) {
  m <- molecular[!depmap_id %in% excluded]
  D <- as.matrix(m[,..programs])
  residual <- D
  for(dr in unique(m$drug)) {
    ii <- which(m$drug==dr)
    residual[ii,] <- sweep(D[ii,,drop=FALSE],2,colMeans(D[ii,,drop=FALSE]))
  }
  scale <- sqrt(colMeans(residual^2))
  stopifnot(all(scale>1e-10))
  norm <- sweep(residual,2,scale,'/')
  V <- svd(norm,nu=0,nv=rank)$v
  axis_scale <- sqrt(colMeans((norm%*%V)^2))
  transform <- sweep(sweep(V,1,scale,'/'),2,axis_scale,'/')
  Z <- D%*%transform
  gates <- list()
  for(K in 1:3) {
    if(K==1) gates[[K]] <- list(centers=matrix(colMeans(Z),1),radius=1) else {
      set.seed(seed+tag*10+K)
      km <- kmeans(Z,centers=K,nstart=20,iter.max=100)
      radius <- max(mean(km$withinss)/mean(km$size),.1)
      gates[[K]] <- list(centers=km$centers,radius=radius)
    }
  }
  # Baseline is used only to supply anchored program estimates where combination
  # participants lack measured RNA. These operators never see functional labels.
  ids <- sort(unique(m$depmap_id))
  B <- as.matrix(baseline[match(ids,depmap_id),..programs])
  bmean <- colMeans(B)
  bsd <- apply(B,2,sd)
  bnorm <- sweep(sweep(B,2,bmean),2,bsd,'/')
  BV <- svd(bnorm,nu=0,nv=8)$v
  bscale <- sqrt(colMeans((bnorm%*%BV)^2))
  operator <- list()
  for(dr in selected$drug) {
    rows <- which(m$drug==dr)
    xb <- sweep(sweep(as.matrix(baseline[match(m$depmap_id[rows],depmap_id),..programs]),2,bmean),2,bsd,'/')%*%BV
    xb <- sweep(xb,2,bscale,'/')
    xm <- colMeans(xb)
    ym <- colMeans(D[rows,,drop=FALSE])
    xb <- sweep(xb,2,xm)
    fit <- glmnet(xb,sweep(D[rows,,drop=FALSE],2,ym),family='mgaussian',alpha=0,lambda=1,
      intercept=FALSE,standardize=FALSE,thresh=1e-9)
    operator[[dr]] <- list(xmean=xm,ymean=ym,beta=vapply(coef(fit),function(x)as.numeric(x)[-1],numeric(8)))
  }
  list(transform=transform,gates=gates,training_cells=unique(m$depmap_id),
    operators=operator,bmean=bmean,bsd=bsd,BV=BV,bscale=bscale)
}
# Query programs remain measured whenever available. Missing RNA is supplied by
# the training-only baseline->program operator; never optimized to fit viability.
query_programs <- function(rep,cell,drug) {
  measured_key <- paste(molecular$depmap_id,molecular$drug)
  ix <- match(paste(cell,drug),measured_key)
  D <- matrix(NA_real_,length(cell),40)
  available <- which(!is.na(ix))
  D[available,] <- as.matrix(molecular[ix[available],..programs])
  missing <- which(is.na(ix))
  if(length(missing)) {
    B <- as.matrix(baseline[match(cell[missing],depmap_id),..programs])
    s <- sweep(sweep(sweep(B,2,rep$bmean),2,rep$bsd,'/')%*%rep$BV,2,rep$bscale,'/')
    for(dr in unique(drug[missing])) {
      ii <- which(drug[missing]==dr)
      op <- rep$operators[[dr]]
      D[missing[ii],] <- sweep(sweep(s[ii,,drop=FALSE],2,op$xmean)%*%op$beta,2,op$ymean,'+')
    }
  }
  stopifnot(all(is.finite(D)))
  list(D=D,measured=!is.na(ix))
}
# Soft Gaussian membership in observed molecular patterns; no drug identities.
pattern_design <- function(D,rep,K) {
  Z <- D%*%rep$transform
  if(K==1) return(list(X=Z,Z=Z,gates=matrix(1,nrow(D),1)))
  centers <- rep$gates[[K]]$centers
  distance <- vapply(1:K,function(k)rowSums(sweep(Z,2,centers[k,])^2),numeric(nrow(Z)))
  logg <- -distance/(2*rep$gates[[K]]$radius)
  logg <- logg-apply(logg,1,max)
  g <- exp(logg)
  g <- g/rowSums(g)
  # One common rule plus heavily penalized local deviations. Zero displacement
  # has zero score, regardless of molecular-pattern memberships.
  X <- cbind(Z,do.call(cbind,lapply(1:K,function(k)Z*g[,k])))
  list(X=X,Z=Z,gates=g)
}
# Both X and y are centered within treatment AFTER the gated features are made.
# Mean-only fits are ineligible; the function has no intercept or drug input.
fit_rule <- function(rep,K,rows) {
  design <- pattern_design(as.matrix(profiles[rows,..programs]),rep,K)
  X <- design$X
  y <- profiles$y[rows]
  drug <- profiles$drug[rows]
  for(dr in unique(drug)) {
    ii <- which(drug==dr)
    X[ii,] <- sweep(X[ii,,drop=FALSE],2,colMeans(X[ii,,drop=FALSE]))
    y[ii] <- y[ii]-mean(y[ii])
  }
  count <- table(drug)
  penalty <- c(rep(1,rank),if(K>1)rep(4,rank*K))
  fit <- glmnet(X,y,alpha=0,lambda=lambda_grid,weights=1/as.numeric(count[drug]),
    penalty.factor=penalty,intercept=FALSE,standardize=FALSE,thresh=1e-9,maxit=100000)
  stopifnot(fit$jerr==0)
  list(rep=rep,K=K,fit=fit,training_cells=unique(profiles$depmap_id[rows]))
}
score_predictions <- function(y,pred,group,fold) {
  if(is.null(dim(pred))) pred <- matrix(pred,ncol=1)
  residual <- y
  block <- paste(group,fold)
  for(g in unique(block)) {
    ii <- which(block==g)
    residual[ii] <- residual[ii]-mean(residual[ii])
    pred[ii,] <- sweep(pred[ii,,drop=FALSE],2,colMeans(pred[ii,,drop=FALSE]))
  }
  r <- vapply(unique(group),function(g) {
    ii <- which(group==g)
    vapply(seq_len(ncol(pred)),function(j)pearson(residual[ii],pred[ii,j]),numeric(1))
  },numeric(ncol(pred)))
  eligible <- rowSums(is.finite(r))==length(unique(group))
  score <- rep(-Inf,ncol(pred))
  score[eligible] <- rowMeans(r[eligible,,drop=FALSE])
  score
}
# Signature contributions are a low-dimensional, functionally learned summary
# of the molecular patterns. Pair inputs are separate and J is symmetric, zero
# when either displacement is zero. K=1/2/3 gives only 1/3/6 interaction terms.
contributions <- function(rule,D,beta) {
  design <- pattern_design(D,rule$rep,rule$K)
  K <- rule$K
  common <- as.numeric(design$Z%*%beta[1:rank])
  if(K==1) return(matrix(common,ncol=1))
  local <- vapply(1:K,function(k) {
    delta <- beta[rank+((k-1)*rank+1:rank)]
    design$gates[,k]*(common+as.numeric(design$Z%*%delta))
  },numeric(nrow(D)))
  stopifnot(max(abs(rowSums(local)-as.numeric(design$X%*%beta)))<1e-8)
  local
}
pair_design <- function(rule,beta,rows) {
  a <- query_programs(rule$rep,pairs$depmap_id[rows],pairs$drugA[rows])
  b <- query_programs(rule$rep,pairs$depmap_id[rows],pairs$drugB[rows])
  ca <- tanh(contributions(rule,a$D,beta))
  cb <- tanh(contributions(rule,b$D,beta))
  ij <- which(upper.tri(matrix(0,rule$K,rule$K),diag=TRUE),arr.ind=TRUE)
  X <- matrix(0,length(rows),nrow(ij))
  for(j in seq_len(nrow(ij))) {
    p <- ij[j,1];q <- ij[j,2]
    X[,j] <- if(p==q) ca[,p]*cb[,p] else (ca[,p]*cb[,q]+ca[,q]*cb[,p])/sqrt(2)
  }
  list(X=X,both_measured=a$measured & b$measured)
}
fit_interaction <- function(rule,beta,rows) {
  design <- pair_design(rule,beta,rows)
  X <- design$X
  y <- pairs$I[rows]
  group <- paste(pairs$SCREENER[rows],pairs$pair[rows])
  for(g in unique(group)) {
    ii <- which(group==g)
    X[ii,] <- sweep(X[ii,,drop=FALSE],2,colMeans(X[ii,,drop=FALSE]))
    y[ii] <- y[ii]-mean(y[ii])
  }
  count <- table(group)
  # glmnet requires >=2 columns; the zero pad adds no functional degree of freedom.
  if(ncol(X)==1) X <- cbind(X,0)
  fit <- glmnet(X,y,alpha=0,lambda=lambda_grid,weights=1/as.numeric(count[group]),
    intercept=FALSE,standardize=TRUE,thresh=1e-9)
  list(fit=fit,terms=ncol(design$X))
}

predictions <- interaction_predictions <- selections <- audits <- list()
for(job in 1:4) {
  label <- if(job==4)'all_data' else paste0('fold_',job)
  outer_excluded <- if(job==4) character() else folds[fold==job,depmap_id]
  outer_excluded <- union('ACH-000248',outer_excluded)
  train <- which(!profiles$depmap_id %in% outer_excluded)
  test <- if(job==4) train else which(profiles$depmap_id %in% outer_excluded)
  pair_train <- which(!pairs$depmap_id %in% outer_excluded)
  pair_test <- if(job==4) pair_train else which(pairs$depmap_id %in% outer_excluded)
  cells <- sort(unique(c(profiles$depmap_id[train],pairs$depmap_id[pair_train])))
  set.seed(seed+job)
  inner <- setNames(sample(rep(1:3,length.out=length(cells))),cells)
  inner_reps <- lapply(1:3,function(j)fit_representation(union(outer_excluded,cells[inner==j]),job*10+j))
  final_rep <- fit_representation(outer_excluded,job*10)
  chosen <- list()
  status(paste(label,'fitting 1–3 molecular rules'))
  for(K in 1:3) {
    held_prediction <- matrix(NA_real_,length(train),length(lambda_grid))
    for(j in 1:3) {
      held <- which(inner[profiles$depmap_id[train]]==j)
      rule <- fit_rule(inner_reps[[j]],K,train[-held])
      design <- pattern_design(as.matrix(profiles[train[held],..programs]),rule$rep,K)$X
      held_prediction[held,] <- predict(rule$fit,newx=design,s=lambda_grid)
      stopifnot(!any(profiles$depmap_id[train[held]] %in% rule$rep$training_cells))
    }
    score <- score_predictions(profiles$y[train],held_prediction,profiles$drug[train],inner[profiles$depmap_id[train]])
    stopifnot(any(is.finite(score)))
    lambda_index <- which.max(score)
    rule <- fit_rule(final_rep,K,train)
    beta <- as.numeric(coef(rule$fit,s=lambda_grid[lambda_index]))[-1]
    rule$beta <- beta
    rule$lambda <- lambda_grid[lambda_index]
    out <- profiles[test,.(depmap_id,drug,dose,fold,observed_log2=y)]
    design <- pattern_design(as.matrix(profiles[test,..programs]),final_rep,K)
    out[,`:=`(K=K,evaluation=if(job==4)'fitted' else 'held_cell',predicted_score=as.numeric(design$X%*%beta))]
    predictions[[paste(label,K)]] <- out
    # Select the tiny pair rule using the same excluded-cell inner partitions.
    pair_held <- matrix(NA_real_,length(pair_train),length(lambda_grid))
    for(j in 1:3) {
      ps <- which(inner[pairs$depmap_id[pair_train]]==j)
      fs <- train[inner[profiles$depmap_id[train]]!=j]
      irule <- fit_rule(inner_reps[[j]],K,fs)
      ibeta <- as.numeric(coef(irule$fit,s=rule$lambda))[-1]
      interaction <- fit_interaction(irule,ibeta,pair_train[-ps])
      X <- pair_design(irule,ibeta,pair_train[ps])$X
      if(ncol(X)==1) X <- cbind(X,0)
      pair_held[ps,] <- predict(interaction$fit,newx=X,s=lambda_grid)
    }
    pair_group <- paste(pairs$SCREENER[pair_train],pairs$pair[pair_train])
    pair_score <- score_predictions(pairs$I[pair_train],pair_held,pair_group,inner[pairs$depmap_id[pair_train]])
    stopifnot(any(is.finite(pair_score)))
    pj <- which.max(pair_score)
    interaction <- fit_interaction(rule,beta,pair_train)
    pair_beta <- as.numeric(coef(interaction$fit,s=lambda_grid[pj]))[-1]
    query <- pair_design(rule,beta,pair_test)
    X <- query$X
    if(ncol(X)==1) X <- cbind(X,0)
    pout <- pairs[pair_test,.(row_id,depmap_id,drugA,drugB,doseA,doseB,pair,SCREENER,fold,observed_departure_log2=I)]
    pout[,`:=`(K=K,evaluation=if(job==4)'fitted' else 'held_cell',predicted_departure_score=as.numeric(X%*%pair_beta),
      both_inputs_measured=query$both_measured)]
    interaction_predictions[[paste(label,K)]] <- pout
    chosen[[K]] <- data.table(job=label,K,lambda=rule$lambda,inner_mean_r=score[lambda_index],
      interaction_lambda=lambda_grid[pj],inner_interaction_r=pair_score[pj])
    rule$pair_beta <- pair_beta
    rule$interaction_lambda <- lambda_grid[pj]
    rule$training_pair_cells <- unique(pairs$depmap_id[pair_train])
    rule$held_cells <- outer_excluded
    rule$fit <- NULL
    qsave(rule,file.path(root,'models',paste0(label,'_K',K,'.qs')))
    status(paste(label,'K=',K,'inner r=',round(score[lambda_index],3),'pair r=',round(pair_score[pj],3)))
  }
  choice <- rbindlist(chosen)
  choice[,selected:=K==K[which.max(inner_mean_r)]]
  selections[[label]] <- choice
  write_table(rbindlist(predictions),'single_predictions_frozen')
  write_table(rbindlist(interaction_predictions),'interaction_predictions_frozen')
  write_table(rbindlist(selections),'model_selection')
}
pred <- rbindlist(predictions)
pred[,observed_variation:=observed_log2-mean(observed_log2),by=.(drug,K,evaluation,fold)]
pred[,predicted_variation:=predicted_score-mean(predicted_score),by=.(drug,K,evaluation,fold)]
metrics <- pred[,.(cells=.N,raw_r=pearson(observed_log2,predicted_score),within_fold_r=pearson(observed_variation,predicted_variation)),by=.(drug,dose,K,evaluation)]
write_table(metrics,'single_correlations')
choice <- rbindlist(selections)
outer_choice <- choice[selected & job!='all_data',.(fold=as.integer(sub('fold_','',job)),K)]
selected_pred <- merge(pred[evaluation=='held_cell'],outer_choice,by=c('fold','K'))
selected_metrics <- selected_pred[,.(cells=.N,raw_r=pearson(observed_log2,predicted_score),within_fold_r=pearson(observed_variation,predicted_variation)),by=.(drug,dose)]
write_table(selected_metrics,'selected_model_single_correlations')
print(selected_metrics)
pout <- rbindlist(interaction_predictions)
pout[,observed_variation:=observed_departure_log2-mean(observed_departure_log2),by=.(pair,SCREENER,K,evaluation,fold)]
pout[,predicted_variation:=predicted_departure_score-mean(predicted_departure_score),by=.(pair,SCREENER,K,evaluation,fold)]
write_table(pout[,.(cells=uniqueN(depmap_id),rows=.N,both_measured=sum(both_inputs_measured),
  raw_r=pearson(observed_departure_log2,predicted_departure_score),within_fold_r=pearson(observed_variation,predicted_variation)),
  by=.(pair,SCREENER,K,evaluation)],'interaction_correlations')
sp <- merge(pout[evaluation=='held_cell'],outer_choice,by=c('fold','K'))
pm <- sp[,.(cells=uniqueN(depmap_id),raw_r=pearson(observed_departure_log2,predicted_departure_score),
  within_fold_r=pearson(observed_variation,predicted_variation)),by=.(pair,SCREENER)]
write_table(pm,'selected_model_interaction_correlations')
print(pm)
# The selected plots show cell variation only; no absolute calibration is claimed.
selected_pred[,panel:=paste0(drug,' | ',dose,' µM')]
lab <- selected_metrics[,.(panel=paste0(drug,' | ',dose,' µM'),label=sprintf('r = %.2f; n = %d',within_fold_r,cells))]
p <- ggplot(selected_pred,aes(observed_variation,predicted_variation))+
  geom_point(aes(color=factor(fold)),size=1.3,alpha=.85)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.45)+
  geom_label(data=lab,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=2.5,label.size=0)+
  facet_wrap(~panel,ncol=5,scales='free')+theme_bw(base_size=10)+
  theme(aspect.ratio=1,legend.position='bottom',strip.text=element_text(size=7))+
  labs(title='Shared molecular-pattern sensitivity rules: training-selected 1–3 patterns',
    x='Observed PRISM log2 abundance variation',y='Predicted sensitivity-score variation',color='Held-cell fold',
    caption='286 measured profiles; whole cells excluded from RNA representation, operators, coefficient fitting and selection.\nBoth axes have drug/dose/fold means removed for sensitivity evaluation; scores are not calibrated viability.')
ggsave(file.path(root,'figures/selected_held_cell_singles.pdf'),p,width=13,height=6.8,device=cairo_pdf)
ggsave(file.path(root,'figures/selected_held_cell_singles.png'),p,width=13,height=6.8,dpi=120)
labp <- pm[,.(pair,SCREENER,label=sprintf('r = %.2f; n = %d',within_fold_r,cells))]
pp <- ggplot(sp,aes(observed_variation,predicted_variation))+
  geom_point(aes(color=factor(fold),shape=both_inputs_measured),size=1.6)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='black',linewidth=.45)+
  geom_label(data=labp,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.02,vjust=1.1,size=2.5,label.size=0)+
  facet_wrap(~SCREENER+pair,ncol=3,scales='free')+theme_bw(base_size=10)+
  theme(aspect.ratio=1,legend.position='bottom',strip.text=element_text(size=7))+
  labs(title='Shared signature interaction: departure from measured single-drug effects',
    x='Observed variation in combination departure (log2 abundance)',y='Predicted departure-score variation',
    color='Held-cell fold',shape='Both RNA inputs measured',
    caption='193 ALMANAC outcomes at supported v87 doses; measured RNA used where available, otherwise anchored baseline-to-RNA estimates.\nPair/screen/fold means removed for evaluation. No joint combination RNA vector is invented.')
ggsave(file.path(root,'figures/selected_held_cell_interactions.pdf'),pp,width=10,height=7,device=cairo_pdf)
ggsave(file.path(root,'figures/selected_held_cell_interactions.png'),pp,width=10,height=7,dpi=120)
final_choice <- choice[job=='all_data' & selected,K]
final <- qread(file.path(root,'models',paste0('all_data_K',final_choice,'.qs')))
physical_rule <- final$rep$transform%*%matrix(final$beta,nrow=rank)
write_table(data.table(program=rep(programs,ncol(physical_rule)),rule=rep(0:(ncol(physical_rule)-1),each=40),
  physical_program_coefficient=as.numeric(physical_rule)),'all_data_program_rule_weights')
paths <- file.path(root,'data',c('single_predictions_frozen.tsv','interaction_predictions_frozen.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')
writeLines(c('Measured-only functional decoder diagnostic, returning to v89 scope plus crizotinib 5 uM.',
  'No PLS; six RNA-only within-treatment covariance axes, molecular K-means centers and soft membership gates.',
  '1–3 shared rules, deviations penalized 4x relative to common rule; ridge strength selected in inner whole-cell folds.',
  'Sensitivity variation only: functional features and labels centered within treatment after feature creation.',
  'Separate symmetric pair interaction on 1/3/6 gated sensitivity contribution products; pair means removed for fitting.',
  'No drug identity or baseline input in F/J. Missing combination RNA supplied by training-only anchored baseline->program ridge operators.',
  'This is not a jointly optimized baseline bridge and not calibrated absolute viability or HGF rescue.',
  'Dose selection remains retrospective. Actual combinations at selected doses are only five supported pairs.',
  capture.output(sessionInfo())),file.path(root,'RUN_DETAILS.txt'))
status('Complete')
