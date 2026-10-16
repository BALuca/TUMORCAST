#!/usr/bin/env Rscript
# Sparse low-rank relationships between separate program responses. Individual
# measured effects anchor the target, so g cannot repair single-decoder errors.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(digest)})
root <- '/work/results/program_pilot/65_sparse_program_interaction_patterns'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
nonlinear <- '/work/results/program_pilot/64_nonlinear_response_composition'
ranks <- 1:3
penalties <- c(.3, .1, .03, .01)
alpha <- .8
max_iterations <- 250L
convergence_tolerance <- 1e-6
set.seed(65011L)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn = 1)
for (folder in c('data','models','figures')) dir.create(file.path(root,folder), recursive=TRUE, showWarnings=FALSE)
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
status <- function(x) {
  message(format(Sys.time(),tz='UTC',usetz=TRUE),' | ',x)
  writeLines(x,file.path(root,'STATUS.txt'))
}
inputs <- c(file.path(prior,'data/matched_molecular_links.tsv.gz'),
  file.path(prior,'data/representation_scores.qs'),file.path(nonlinear,'data/predictions.tsv.gz'))
stopifnot(all(file.exists(inputs)))
obs <- fread(inputs[1])[time_h==3 & molecular_supported==TRUE]
setorder(obs,anchor,partner,nominal_A_uM,nominal_B_uM,biological_replicate)
bank <- qread(inputs[2])
stopifnot(nrow(obs)==496L,uniqueN(obs$condition_id)==198L,!anyDuplicated(obs$observation_id))
profiles <- bank$metadata[time_h==3,profile_id]
response <- bank$original40[profiles,,drop=FALSE]
rms <- sqrt(colMeans(response^2))
stopifnot(ncol(response)==40L,all(is.finite(response)),all(rms>0))
response <- sweep(response,2,rms,'/')
a <- response[obs$profile_A,,drop=FALSE]
b <- response[obs$profile_B,,drop=FALSE]
stopifnot(all(obs$ratio_A>0),all(obs$ratio_B>0),all(obs$ratio_AB>0),!anyNA(a),!anyNA(b))
offset <- log2(obs$ratio_A*obs$ratio_B)
y <- log2(obs$ratio_AB)-offset
programs <- colnames(a)
partners <- sort(unique(obs$partner))
write_table(data.table(path=inputs,sha256=vapply(inputs,digest,character(1),file=TRUE,algo='sha256')),'input_manifest')
qsave(list(ranks=ranks,penalties=penalties,alpha=alpha,max_iterations=max_iterations,
  convergence_tolerance=convergence_tolerance,rms=rms,programs=programs,
  original_weight_sha256=bank$original_weight_sha256,primary_time_h=3),file.path(root,'models/settings.qs'))

# Identical weighting and raw/within-first-drug-dose objective as62 and64.
balanced_weights <- function(md) {
  group <- interaction(md$anchor,md$partner,md$biological_replicate,drop=TRUE)
  w <- 1/as.numeric(table(group)[as.character(group)])
  w/sum(w)
}
center_context <- function(x,md) {
  group <- interaction(md$anchor,md$nominal_A_uM,md$biological_replicate,drop=TRUE)
  out <- x
  for(idx in split(seq_len(nrow(md)),group)) {
    if(is.matrix(x)) out[idx,] <- sweep(x[idx,,drop=FALSE],2,colMeans(x[idx,,drop=FALSE]),'-')
    else out[idx] <- x[idx]-mean(x[idx])
  }
  out
}
objective <- function(y,pred,md) {
  w <- balanced_weights(md)
  .5*sum(w*(y-pred)^2)+.5*sum(w*(center_context(y,md)-center_context(pred,md))^2)
}
pattern_contributions <- function(a,b,u,v) .5*((a%*%u)*(b%*%v)+(a%*%v)*(b%*%u))
design_for_factor <- function(a,b,v) {
  av <- a%*%v
  bv <- b%*%v
  do.call(cbind,lapply(seq_len(ncol(v)),function(k) .5*(a*bv[,k]+b*av[,k])))
}
penalty_value <- function(u,v,lambda) {
  lambda*(alpha*(sum(abs(u))+sum(abs(v)))+(1-alpha)/2*(sum(u^2)+sum(v^2)))
}

# Each block is ordinary elastic-net regression solved by glmnet. Balancing a
# factor pair leaves predictions unchanged and minimizes its scale penalty.
fit_patterns <- function(idx,rank,lambda,start) {
  aa <- a[idx,,drop=FALSE]
  bb <- b[idx,,drop=FALSE]
  md <- obs[idx]
  yy <- y[idx]
  w <- balanced_weights(md)
  # glmnet's Gaussian path standardizes its response internally. Explicit unit
  # RMS on the augmented target makes its elastic-net penalty agree exactly
  # with the monitored block objective; scale the fitted factors back afterward.
  target_rms <- sqrt(.5*sum(w*yy^2)+.5*sum(w*center_context(yy,md)^2))
  stopifnot(target_rms>0)
  yy <- yy/target_rms
  yc <- center_context(yy,md)
  # A symmetric supervised initialization uses ONLY the current training labels.
  m <- crossprod(aa,bb*(w*yy))
  m <- (m+t(m))/2
  eig <- eigen(m,symmetric=TRUE)
  ii <- order(abs(eig$values),decreasing=TRUE)[seq_len(rank)]
  u <- eig$vectors[,ii,drop=FALSE]
  v <- sweep(u,2,sign(eig$values[ii]),'*')
  if(start==2L) {
    set.seed(65011L+rank)
    u <- u+matrix(rnorm(length(u),sd=.12),nrow(u))
    v <- v+matrix(rnorm(length(v),sd=.12),nrow(v))
  }
  initial <- rowSums(pattern_contributions(aa,bb,u,v))
  multiplier <- max(.01,sum(w*initial*yy)/max(sum(w*initial^2),1e-12))
  u <- u*sqrt(multiplier)
  v <- v*sqrt(multiplier)
  ww <- c(w,w)/2
  target <- c(yy,yc)
  evaluate <- function(u,v) {
    pred <- rowSums(pattern_contributions(aa,bb,u,v))
    .5*objective(yy,pred,md)+penalty_value(u,v,lambda)
  }
  previous <- evaluate(u,v)
  trace <- previous
  converged <- FALSE
  for(iter in seq_len(max_iterations)) {
    if(all(u==0) || all(v==0)) {
      u[,] <- 0
      v[,] <- 0
      current <- evaluate(u,v)
      trace <- c(trace,current)
      converged <- TRUE
      break
    }
    du <- design_for_factor(aa,bb,v)
    fit <- glmnet(rbind(du,center_context(du,md)),target,weights=ww,
      lambda=lambda,alpha=alpha,standardize=FALSE,intercept=FALSE,thresh=1e-10,maxit=100000L)
    stopifnot(fit$jerr==0L)
    u <- matrix(as.numeric(coef(fit))[-1],nrow=length(programs),ncol=rank)
    if(all(u==0)) {
      v[,] <- 0
      current <- evaluate(u,v)
      trace <- c(trace,current)
      converged <- TRUE
      break
    }
    dv <- design_for_factor(aa,bb,u)
    fit <- glmnet(rbind(dv,center_context(dv,md)),target,weights=ww,
      lambda=lambda,alpha=alpha,standardize=FALSE,intercept=FALSE,thresh=1e-10,maxit=100000L)
    stopifnot(fit$jerr==0L)
    v <- matrix(as.numeric(coef(fit))[-1],nrow=length(programs),ncol=rank)
    for(k in seq_len(rank)) {
      if(sum(abs(u[,k]))==0 || sum(abs(v[,k]))==0) {
        u[,k] <- 0
        v[,k] <- 0
      } else {
        balance <- optimize(function(logscale) {
          scale <- exp(logscale)
          penalty_value(u[,k]*scale,v[,k]/scale,lambda)
        },interval=c(-8,8))$minimum
        u[,k] <- u[,k]*exp(balance)
        v[,k] <- v[,k]/exp(balance)
      }
    }
    current <- evaluate(u,v)
    if(!is.finite(current) || current>previous+1e-7*(1+abs(previous))) {
      stop(sprintf('Block objective increased: rank=%d lambda=%g start=%d iteration=%d old=%.10g new=%.10g',
        rank,lambda,start,iter,previous,current))
    }
    trace <- c(trace,current)
    if(abs(previous-current)<convergence_tolerance*(1+abs(previous))) {
      converged <- TRUE
      break
    }
    previous <- current
  }
  dimnames(u) <- dimnames(v) <- list(programs,paste0('Pattern',seq_len(rank)))
  list(u=u*sqrt(target_rms),v=v*sqrt(target_rms),target_rms=target_rms,
    rank=rank,lambda=lambda,start=start,converged=converged,
    iterations=iter,penalized_objective=current,objective_trace=trace,
    effective_patterns=sum(colSums(abs(u))*colSums(abs(v))>1e-12))
}

# Cache only identical training subsets inside this execution. Two outer folds
# share an inner training set; recomputing it would change no evidence.
cache <- new.env(parent=emptyenv())
fit_audit <- list()
get_fit <- function(training_partners,rank,lambda) {
  key <- paste(paste(sort(training_partners),collapse='__'),rank,lambda,sep='|')
  if(exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
  idx <- which(obs$partner %chin% training_partners)
  fits <- lapply(1:2,function(s) fit_patterns(idx,rank,lambda,s))
  valid <- which(vapply(fits,function(f) f$converged,logical(1)))
  stopifnot(length(valid)>0L)
  best <- valid[which.min(vapply(fits[valid],function(f) f$penalized_objective,numeric(1)))]
  fit <- fits[[best]]
  fit$training_partners <- training_partners
  fit$training_observation_ids <- obs$observation_id[idx]
  fit_audit[[length(fit_audit)+1L]] <<- rbindlist(lapply(fits,function(f) data.table(
    training_set=paste(sort(training_partners),collapse='__'),rank=rank,lambda=lambda,start=f$start,
    converged=f$converged,iterations=f$iterations,penalized_objective=f$penalized_objective,
    selected_start=f$start==best,effective_patterns=f$effective_patterns)))
  assign(key,fit,envir=cache)
  fit
}
grid <- CJ(rank=ranks,lambda=penalties)
predictions <- audits <- selections <- loadings <- matrices <- contributions <- list()
for(held in partners) {
  tr_partners <- setdiff(partners,held)
  tr <- which(obs$partner %chin% tr_partners)
  te <- which(obs$partner==held)
  status(paste('exclude',held,'| fit sparse patterns with inner partner exclusions'))
  loss <- matrix(NA_real_,nrow(grid)+1L,length(tr_partners))
  for(j in seq_along(tr_partners)) {
    val <- which(obs$partner==tr_partners[j])
    training <- setdiff(tr_partners,tr_partners[j])
    for(k in seq_len(nrow(grid))) {
      fit <- get_fit(training,grid$rank[k],grid$lambda[k])
      pred <- rowSums(pattern_contributions(a[val,,drop=FALSE],b[val,,drop=FALSE],fit$u,fit$v))
      loss[k,j] <- objective(y[val],pred,obs[val])
    }
    loss[nrow(loss),j] <- objective(y[val],rep(0,length(val)),obs[val])
  }
  tuning <- rbind(grid,data.table(rank=0L,lambda=Inf))
  tuning[, `:=`(internal_loss=rowMeans(loss),fold_id=held)]
  selected <- which.min(tuning$internal_loss)
  tuning[,selected:=seq_len(.N)==selected]
  selections[[held]] <- tuning
  if(tuning$rank[selected]==0L) {
    fit <- list(u=matrix(0,40,1,dimnames=list(programs,'Pattern1')),
      v=matrix(0,40,1,dimnames=list(programs,'Pattern1')),rank=0L,lambda=Inf,
      effective_patterns=0L,converged=TRUE,training_partners=tr_partners,
      training_observation_ids=obs$observation_id[tr],objective_trace=numeric())
  } else fit <- get_fit(tr_partners,tuning$rank[selected],tuning$lambda[selected])
  model_file <- file.path(root,'models',paste0('3h_exclude_',held,'.qs'))
  qsave(fit,model_file)
  fit <- qread(model_file)
  parts <- pattern_contributions(a,b,fit$u,fit$v)
  pred <- rowSums(parts)
  swap <- rowSums(pattern_contributions(b,a,fit$u,fit$v))
  zero <- rowSums(pattern_contributions(a*0,b,fit$u,fit$v))
  w <- (fit$u%*%t(fit$v)+fit$v%*%t(fit$u))/2
  stopifnot(max(abs(pred-swap))<1e-10,max(abs(zero))==0,
    max(abs(pred-rowSums((a%*%w)*b)))<1e-10,all(is.finite(pred)))
  source_ids <- function(z) unique(unlist(strsplit(z,';',fixed=TRUE)))
  overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]),source_ids(obs$source_sample_ids_AB[te]))
  stopifnot(length(overlap)==0L)
  out <- obs[,.(observation_id,condition_id,anchor,partner,nominal_A_uM,nominal_B_uM,
    biological_replicate,profile_A,profile_B)]
  out[, `:=`(time_h=3L,fold_id=held,fold_type='single_partner',
    model='Sparse_patterns',split=fifelse(partner==held,'held','training'),
    predicted_single_sum_log2=offset,predicted_interaction_log2=pred,
    predicted_ratio_AB=2^(offset+pred))]
  write_table(out,paste0('prediction_3h_',held))
  audits[[held]] <- data.table(fold_id=held,rank=fit$rank,effective_patterns=fit$effective_patterns,
    lambda=fit$lambda,source_overlap=length(overlap),swap_error=max(abs(pred-swap)),zero_error=max(abs(zero)),
    replay_error=max(abs(pred-rowSums((a%*%w)*b))),
    prediction_sha256=digest(file.path(root,'data',paste0('prediction_3h_',held,'.tsv.gz')),file=TRUE,algo='sha256'))
  out[, `:=`(ratio_A=obs$ratio_A,ratio_B=obs$ratio_B,ratio_AB=obs$ratio_AB)]
  predictions[[held]] <- out
  loadings[[held]] <- rbindlist(lapply(seq_len(ncol(fit$u)),function(k) data.table(
    fold_id=held,pattern=k,program=programs,U=fit$u[,k],V=fit$v[,k])))
  matrices[[held]] <- data.table(fold_id=held,program_A=rep(programs,40),
    program_B=rep(programs,each=40),coefficient=as.vector(w))
  contributions[[held]] <- rbindlist(lapply(seq_len(ncol(parts)),function(k) data.table(
    observation_id=obs$observation_id,fold_id=held,pattern=k,contribution_log2=parts[,k])))
  write_table(rbindlist(predictions),'predictions')
  write_table(rbindlist(audits),'validation')
  write_table(rbindlist(selections),'selection')
  write_table(rbindlist(loadings),'pattern_loadings')
  write_table(rbindlist(matrices),'interaction_matrices')
  write_table(rbindlist(contributions),'pattern_contributions')
  write_table(rbindlist(fit_audit),'optimization_audit')
  status(paste('saved',held,'| selected',fit$effective_patterns,'patterns; lambda',fit$lambda))
}
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
status('Sparse interaction models complete; primary predictions saved.')
