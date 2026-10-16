#!/usr/bin/env Rscript
# Sparse low-rank relationships between separate program responses. Individual
# measured effects anchor the target, so g cannot repair single-decoder errors.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(digest)})
root <- '/work/results/program_pilot/66_bootstrapped_program_interaction_patterns'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
nonlinear <- '/work/results/program_pilot/64_nonlinear_response_composition'
ranks <- c(1L,3L,6L)
penalties <- c(.1,.03,.01,.003)
alpha <- .5
inner_bootstraps <- 6L
final_bootstraps <- 20L
max_iterations <- 750L
convergence_tolerance <- 1e-5
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
  inner_bootstraps=inner_bootstraps,final_bootstraps=final_bootstraps,
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
fit_patterns <- function(idx,rank,lambda,start,partner_weights) {
  aa <- a[idx,,drop=FALSE]
  bb <- b[idx,,drop=FALSE]
  md <- obs[idx]
  yy <- y[idx]
  w <- balanced_weights(md)
  w <- w*partner_weights[md$partner]
  w <- as.numeric(w/sum(w))
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
    .25*sum(w*(yy-pred)^2)+.25*sum(w*(center_context(yy,md)-center_context(pred,md))^2)+
      penalty_value(u,v,lambda)
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
  list(u=u*sqrt(target_rms),v=v*sqrt(target_rms),target_rms=target_rms,partner_weights=partner_weights,
    rank=rank,lambda=lambda,start=start,converged=converged,
    iterations=iter,penalized_objective=current,objective_trace=trace,
    effective_patterns=sum(colSums(abs(u))*colSums(abs(v))>1e-12))
}

# Bayesian bootstrap draws preserve whole partners, and all their conditions.
# Reuse identical draws across candidate settings for comparable inner errors.
bootstrap_weights <- function(training,b) {
  code <- sum(2^(match(sort(training),partners)-1L))
  set.seed(66011L+100L*code+b)
  w <- if(b==0L) rep(1,length(training)) else rexp(length(training))
  setNames(w/mean(w),training)
}
fit_audit <- list()
fit_member <- function(training,rank,lambda,b) {
  idx <- which(obs$partner %chin% training)
  pw <- bootstrap_weights(training,b)
  fits <- lapply(1:2,function(s) fit_patterns(idx,rank,lambda,s,pw))
  valid <- which(vapply(fits,function(f) f$converged,logical(1)))
  stopifnot(length(valid)>0L)
  selected <- valid[which.min(vapply(fits[valid],function(f) f$penalized_objective,numeric(1)))]
  best <- fits[[selected]]
  best$W <- (best$u%*%t(best$v)+best$v%*%t(best$u))/2
  fit_audit[[length(fit_audit)+1L]] <<- rbindlist(lapply(fits,function(f) data.table(
    training_set=paste(sort(training),collapse='__'),rank=rank,lambda=lambda,bootstrap=b,start=f$start,
    converged=f$converged,iterations=f$iterations,penalized_objective=f$penalized_objective,
    selected_start=f$start==selected,effective_patterns=f$effective_patterns)))
  best
}
cache <- new.env(parent=emptyenv())
get_ensemble <- function(training,rank,lambda,n_bootstraps) {
  key <- paste(paste(sort(training),collapse='__'),rank,lambda,n_bootstraps,sep='|')
  if(exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
  members <- lapply(seq_len(n_bootstraps),function(b) fit_member(training,rank,lambda,b))
  w <- Reduce('+',lapply(members,function(f) f$W))/n_bootstraps
  representative <- which.min(vapply(members,function(f) sum((f$W-w)^2),numeric(1)))
  fit <- list(W=w,members=members,rank=rank,lambda=lambda,
    effective_patterns=median(vapply(members,function(f) f$effective_patterns,numeric(1))),
    representative_member=representative,u=members[[representative]]$u,v=members[[representative]]$v,
    converged=TRUE,training_partners=training,training_observation_ids=obs[partner %chin% training,observation_id])
  assign(key,fit,envir=cache)
  fit
}
grid <- CJ(rank=ranks,lambda=penalties)
predictions <- audits <- selections <- edges <- member_predictions <- list()
for(held in partners) {
  tr_partners <- setdiff(partners,held)
  tr <- which(obs$partner %chin% tr_partners)
  te <- which(obs$partner==held)
  status(paste('exclude',held,'| grouped bootstrap inner selection'))
  loss <- matrix(NA_real_,nrow(grid)+1L,length(tr_partners))
  for(j in seq_along(tr_partners)) {
    val <- which(obs$partner==tr_partners[j])
    training <- setdiff(tr_partners,tr_partners[j])
    for(k in seq_len(nrow(grid))) {
      fit <- get_ensemble(training,grid$rank[k],grid$lambda[k],inner_bootstraps)
      pred <- rowSums((a[val,,drop=FALSE]%*%fit$W)*b[val,,drop=FALSE])
      loss[k,j] <- objective(y[val],pred,obs[val])
    }
    loss[nrow(loss),j] <- objective(y[val],rep(0,length(val)),obs[val])
    status(paste('exclude',held,'| inner partner',tr_partners[j],'complete'))
    write_table(rbindlist(fit_audit),'optimization_audit')
  }
  tuning <- rbind(grid,data.table(rank=0L,lambda=Inf))
  tuning[,`:=`(internal_loss=rowMeans(loss),fold_id=held)]
  selected <- which.min(tuning$internal_loss)
  tuning[,selected:=seq_len(.N)==selected]
  selections[[held]] <- tuning
  if(tuning$rank[selected]==0L) {
    zero <- list(u=matrix(0,40,1,dimnames=list(programs,'Pattern1')),
      v=matrix(0,40,1,dimnames=list(programs,'Pattern1')),
      W=matrix(0,40,40,dimnames=list(programs,programs)),effective_patterns=0L)
    fit <- c(zero,list(rank=0L,lambda=Inf,converged=TRUE,members=rep(list(zero),final_bootstraps),
      training_partners=tr_partners,training_observation_ids=obs$observation_id[tr],representative_member=1L))
    unbagged <- zero
  } else {
    fit <- get_ensemble(tr_partners,tuning$rank[selected],tuning$lambda[selected],final_bootstraps)
    unbagged <- fit_member(tr_partners,fit$rank,fit$lambda,0L)
  }
  fit$unbagged <- unbagged
  file <- file.path(root,'models',paste0('3h_exclude_',held,'.qs'))
  qsave(fit,file)
  fit <- qread(file)
  w <- fit$W
  pred <- rowSums((a%*%w)*b)
  swap <- rowSums((b%*%w)*a)
  member_pred <- vapply(fit$members,function(f) rowSums((a%*%f$W)*b),numeric(nrow(obs)))
  stopifnot(max(abs(pred-swap))<1e-10,max(abs(pred-rowMeans(member_pred)))<1e-10,
    max(abs(rowSums(((a*0)%*%w)*b)))==0,all(is.finite(pred)))
  source_ids <- function(z) unique(unlist(strsplit(z,';',fixed=TRUE)))
  overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]),source_ids(obs$source_sample_ids_AB[te]))
  stopifnot(length(overlap)==0L,!any(obs$partner[te] %chin% fit$training_partners))
  out <- obs[,.(observation_id,condition_id,anchor,partner,nominal_A_uM,nominal_B_uM,
    biological_replicate,profile_A,profile_B)]
  out[,`:=`(time_h=3L,fold_id=held,fold_type='single_partner',model='Bootstrap_patterns',
    split=fifelse(partner==held,'held','training'),predicted_single_sum_log2=offset,
    predicted_interaction_log2=pred,predicted_ratio_AB=2^(offset+pred),
    unbagged_interaction_log2=rowSums((a%*%fit$unbagged$W)*b),
    bootstrap_sd_log2=apply(member_pred,1,sd))]
  write_table(out,paste0('prediction_3h_',held))
  audits[[held]] <- data.table(fold_id=held,rank=fit$rank,effective_patterns=fit$effective_patterns,
    lambda=fit$lambda,source_overlap=length(overlap),swap_error=max(abs(pred-swap)),zero_error=0,
    replay_error=max(abs(pred-rowMeans(member_pred))),
    prediction_sha256=digest(file.path(root,'data',paste0('prediction_3h_',held,'.tsv.gz')),file=TRUE,algo='sha256'))
  out[,`:=`(ratio_A=obs$ratio_A,ratio_B=obs$ratio_B,ratio_AB=obs$ratio_AB)]
  predictions[[held]] <- out
  wc <- vapply(fit$members,function(f) as.vector(f$W),numeric(1600))
  edges[[held]] <- data.table(fold_id=held,program_A=rep(programs,40),program_B=rep(programs,each=40),
    mean_coefficient=as.vector(w),sd_coefficient=apply(wc,1,sd),
    nonzero_fraction=rowMeans(abs(wc)>1e-8),positive_fraction=rowMeans(wc>1e-8),negative_fraction=rowMeans(wc< -1e-8))
  member_predictions[[held]] <- rbindlist(lapply(seq_len(ncol(member_pred)),function(k) data.table(
    observation_id=obs$observation_id,fold_id=held,bootstrap=k,predicted_g=member_pred[,k])))
  write_table(rbindlist(predictions),'predictions')
  write_table(rbindlist(audits),'validation')
  write_table(rbindlist(selections),'selection')
  write_table(rbindlist(edges),'interaction_edge_stability')
  write_table(rbindlist(member_predictions),'bootstrap_predictions')
  write_table(rbindlist(fit_audit),'optimization_audit')
  status(paste('saved',held,'| capacity',fit$rank,'median active',fit$effective_patterns,'lambda',fit$lambda))
}
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
status('Bootstrapped interaction models complete; primary predictions saved.')

