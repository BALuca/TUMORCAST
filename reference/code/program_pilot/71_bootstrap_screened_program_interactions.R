#!/usr/bin/env Rscript
# Bootstrap simple interactions before joint elastic-net fitting; enforce heredity.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(digest)})
root <- '/work/results/program_pilot/71_bootstrap_screened_program_interactions'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data'
lambda_grid <- 10^seq(1.5,-3.5,length.out=32L)
alpha_grid <- c(.1,.5,1)
bootstrap_draws <- 500L
bootstrap_interval <- c(.05, .95)
bootstrap_seed <- 71011L
set.seed(bootstrap_seed)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn=1)
for(folder in c('data','models','figures')) dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
status <- function(x) {
  message(format(Sys.time(),tz='UTC',usetz=TRUE),' | ',x)
  writeLines(x,file.path(root,'STATUS.txt'))
}
inputs <- file.path(prior,c('matched_molecular_links.tsv.gz','representation_scores.qs'))
stopifnot(all(file.exists(inputs)))
obs <- fread(inputs[1])[time_h==3 & molecular_supported==TRUE]
setorder(obs,anchor,partner,nominal_A_uM,nominal_B_uM,biological_replicate)
bank <- qread(inputs[2])
stopifnot(nrow(obs)==496L,uniqueN(obs$condition_id)==198L,!anyDuplicated(obs$observation_id))
response <- bank$original40[bank$metadata[time_h==3,profile_id],,drop=FALSE]
rms <- sqrt(colMeans(response^2))
stopifnot(all(rms>0),ncol(response)==40L)
response <- sweep(response,2,rms,'/')
a <- response[obs$profile_A,,drop=FALSE]
b <- response[obs$profile_B,,drop=FALSE]
stopifnot(all(is.finite(a)),all(is.finite(b)),all(obs$ratio_A>0),all(obs$ratio_B>0),all(obs$ratio_AB>0))
cross_features <- function(a,b) {
  ij <- which(upper.tri(matrix(0,ncol(a),ncol(a)),diag=TRUE),arr.ind=TRUE)
  x <- vapply(seq_len(nrow(ij)),function(k) {
    i <- ij[k,1]
    j <- ij[k,2]
    (a[,i]*b[,j]+a[,j]*b[,i])/2
  },numeric(nrow(a)))
  colnames(x) <- paste(colnames(a)[ij[,1]],colnames(a)[ij[,2]],sep='__')
  x
}
smooth_basis <- function(x) {
  z <- tanh(x)
  colnames(z) <- paste0(colnames(x),'_sat')
  cbind(x,z)
}
designs <- list(Full=cross_features(a,b),Smooth=cross_features(smooth_basis(a),smooth_basis(b)))
stopifnot(ncol(designs$Full)==820L,ncol(designs$Smooth)==3240L,
  max(abs(designs$Smooth-cross_features(smooth_basis(b),smooth_basis(a))))<1e-12,
  max(abs(cross_features(smooth_basis(a*0),smooth_basis(b))))==0)
offset <- log2(obs$ratio_A*obs$ratio_B)
target <- log2(obs$ratio_AB)-offset
partners <- sort(unique(obs$partner))
write_table(data.table(path=inputs,sha256=vapply(inputs,digest,character(1),file=TRUE,algo='sha256')),'input_manifest')
qsave(list(lambda_grid=lambda_grid,alpha_grid=alpha_grid,rms=rms,primary_time_h=3,
  bootstrap_draws=bootstrap_draws,bootstrap_interval=bootstrap_interval,
  bootstrap_unit='whole partner including both anchors, all doses and repeats',
  original_weight_sha256=bank$original_weight_sha256,nonlinear_basis='x and tanh(x), all symmetric products'),
  file.path(root,'models/settings.qs'))

# Reuse62 weighting/contrast objective; dose trajectories remain intact.
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
  .5*colSums((pred-y)^2*w)+.5*colSums((center_context(pred,md)-center_context(y,md))^2*w)
}
# Screen a simple slope with an intercept. Bootstrap partner-level weighted
# sufficient statistics to reproduce ordinary weighted regression efficiently.
# The mathematical slope formula is checked against lm.wfit on real inputs.
screen_cache <- new.env(parent=emptyenv())
bootstrap_screen <- function(training) {
  training <- sort(training)
  key <- paste(training, collapse='__')
  if (exists(key, screen_cache, inherits=FALSE)) return(screen_cache[[key]])
  idx <- which(obs$partner %chin% training)
  xx <- designs$Full[idx,,drop=FALSE]
  yy <- target[idx]
  ww <- balanced_weights(obs[idx])
  group <- factor(obs$partner[idx], levels=training)
  group_w <- as.numeric(rowsum(ww, group, reorder=FALSE)[training,,drop=FALSE])
  group_x <- rowsum(xx*ww, group, reorder=FALSE)[training,,drop=FALSE]
  group_xx <- rowsum(xx^2*ww, group, reorder=FALSE)[training,,drop=FALSE]
  group_xy <- rowsum(xx*(ww*yy), group, reorder=FALSE)[training,,drop=FALSE]
  group_y <- as.numeric(rowsum(ww*yy, group, reorder=FALSE)[training,,drop=FALSE])
  set.seed(bootstrap_seed + sum(2^match(training, partners)))
  frequency <- t(replicate(bootstrap_draws,
    tabulate(sample.int(length(training), length(training), replace=TRUE), nbins=length(training))))
  total_w <- as.numeric(frequency %*% group_w)
  sum_x <- frequency %*% group_x
  sum_xx <- frequency %*% group_xx
  sum_xy <- frequency %*% group_xy
  mean_y <- as.numeric(frequency %*% group_y) / total_w
  denominator <- sum_xx - sum_x^2/total_w
  slopes <- (sum_xy - sum_x*mean_y)/denominator
  slopes[denominator <= pmax(1e-14, abs(sum_xx)*1e-12)] <- NA_real_
  for (column in c(1L, 410L, 820L)) {
    sample_w <- ww * frequency[1, as.integer(group)]
    checked <- lm.wfit(cbind(1, xx[,column]), yy, w=sample_w)$coefficients[2]
    if (is.finite(slopes[1,column])) stopifnot(abs(checked-slopes[1,column]) < 1e-8*(1+abs(checked)))
  }
  bounds <- apply(slopes, 2, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(c(NA_real_,NA_real_))
    quantile(z, bootstrap_interval, names=FALSE)
  })
  table <- data.table(term=colnames(xx), lower=as.numeric(bounds[1,]), upper=as.numeric(bounds[2,]),
    median_slope=apply(slopes,2,median,na.rm=TRUE),
    positive_fraction=colMeans(slopes>0,na.rm=TRUE), finite_fraction=colMeans(is.finite(slopes)))
  table[, retained := finite_fraction >= .95 & is.finite(lower) & is.finite(upper) & (lower>0 | upper<0)]
  table[, `:=`(training_partners=key, bootstrap_draws=bootstrap_draws)]
  result <- list(table=table, terms=table[retained==TRUE,term], training_partners=training)
  screen_cache[[key]] <- result
  write_table(table, paste0('bootstrap_screen_',key))
  result
}
parent_terms <- lapply(designs, function(x) vapply(strsplit(colnames(x),'__',fixed=TRUE),
  function(z) paste(sort(sub('_sat$','',z)),collapse='__'),character(1)))
cache <- new.env(parent=emptyenv())
fit_path <- function(training,variant,alpha) {
  key <- paste(paste(sort(training),collapse='__'),variant,alpha,sep='|')
  if(exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
  idx <- which(obs$partner %chin% training)
  x <- designs[[variant]][idx,,drop=FALSE]
  y <- target[idx]
  md <- obs[idx]
  weights <- balanced_weights(md)
  sx <- sqrt(colSums(x^2*weights))
  screen <- bootstrap_screen(training)
  keep <- sx>1e-10 & parent_terms[[variant]] %chin% screen$terms
  xx <- sweep(x[,keep,drop=FALSE],2,sx[keep],'/')
  yy <- c(y,center_context(y,md))
  ww <- c(weights,weights)/2
  sy <- sqrt(sum(ww*yy^2))
  stopifnot(sy>0)
  # glmnet requires at least two columns. An inert zero column supports a screen
  # retaining exactly one real interaction; no surviving terms means g=0.
  padded <- sum(keep)==1L
  fit <- NULL
  if (sum(keep)>0L) {
    if (padded) xx <- cbind(xx, inert_zero=0)
    fit <- glmnet(rbind(xx,center_context(xx,md)),yy/sy,weights=ww,alpha=alpha,
      lambda=lambda_grid,standardize=FALSE,intercept=FALSE,thresh=1e-8,maxit=100000L)
    stopifnot(fit$jerr==0L,length(fit$lambda)==length(lambda_grid))
  }
  result <- list(fit=fit,keep=keep,sx=sx[keep],sy=sy,columns=colnames(x),alpha=alpha,variant=variant,
    training_partners=training,training_observation_ids=obs$observation_id[idx],
    screen=screen$table, retained_simple_terms=screen$terms, padded=padded)
  assign(key,result,envir=cache)
  result
}
predict_path <- function(fit,x) {
  stopifnot(identical(colnames(x),fit$columns))
  if (is.null(fit$fit)) return(matrix(0,nrow(x),length(lambda_grid)))
  xx <- sweep(x[,fit$keep,drop=FALSE],2,fit$sx,'/')
  if (fit$padded) xx <- cbind(xx,inert_zero=0)
  as.matrix(predict(fit$fit,newx=xx,s=lambda_grid))*fit$sy
}
safe_cor <- function(x,y) if(sd(x)<1e-10 || sd(y)<1e-10) NA_real_ else cor(x,y)
predictions <- selections <- audits <- coefficients <- list()
for(held in partners) {
  training <- setdiff(partners,held)
  tr <- which(obs$partner %chin% training)
  te <- which(obs$partner==held)
  status(paste('exclude',held,'| bootstrap screen inside each training fold; elastic-net paths'))
  tuning_parts <- list()
  for(variant in names(designs)) for(alpha in alpha_grid) {
    loss <- matrix(NA_real_,length(training),length(lambda_grid))
    for(j in seq_along(training)) {
      val <- which(obs$partner==training[j])
      fit <- fit_path(setdiff(training,training[j]),variant,alpha)
      pred <- predict_path(fit,designs[[variant]][val,,drop=FALSE])
      loss[j,] <- objective(target[val],pred,obs[val])
    }
    tuning_parts[[paste(variant,alpha)]] <- data.table(variant=variant,alpha=alpha,
      lambda=lambda_grid,lambda_index=seq_along(lambda_grid),internal_loss=colMeans(loss))
  }
  tuning <- rbindlist(tuning_parts)
  zero_loss <- mean(vapply(training,function(pp) {
    idx <- which(obs$partner==pp)
    objective(target[idx],matrix(0,length(idx),1),obs[idx])
  },numeric(1)))
  tuning <- rbind(tuning,data.table(variant=names(designs),alpha=0,lambda=Inf,lambda_index=NA_integer_,internal_loss=zero_loss))
  tuning[,`:=`(fold_id=held,selected_within_variant=FALSE,selected_overall=FALSE)]
  winners <- tuning[, .I[which.min(internal_loss)],by=variant]$V1
  overall <- winners[which.min(tuning$internal_loss[winners])]
  tuning[winners,selected_within_variant:=TRUE]
  tuning[overall,selected_overall:=TRUE]
  selections[[held]] <- tuning
  for(chosen in winners) {
    variant <- tuning$variant[chosen]
    zero <- is.infinite(tuning$lambda[chosen])
    fit <- fit_path(training,variant,if(zero) alpha_grid[1] else tuning$alpha[chosen])
    fit$selected_lambda <- tuning$lambda[chosen]
    fit$zero <- zero
    fit$selected_overall <- chosen==overall
    pred <- if(zero) rep(0,nrow(obs)) else predict_path(fit,designs[[variant]])[,tuning$lambda_index[chosen]]
    raw <- setNames(rep(0,ncol(designs[[variant]])),colnames(designs[[variant]]))
    if(!zero && !is.null(fit$fit)) {
      beta <- as.numeric(coef(fit$fit,s=fit$selected_lambda))[-1]
      if(fit$padded) beta <- beta[-length(beta)]
      raw[fit$keep] <- beta*fit$sy/fit$sx
    }
    stopifnot(all(parent_terms[[variant]][raw!=0] %chin% fit$retained_simple_terms),
      !held %chin% fit$training_partners)
    fit$raw_coefficients <- raw
    file <- file.path(root,'models',paste0(held,'_',variant,'.qs'))
    qsave(fit,file)
    saved <- qread(file)
    replay <- as.numeric(designs[[variant]]%*%saved$raw_coefficients)
    stopifnot(max(abs(pred-replay))<1e-8,all(is.finite(2^(offset+pred))))
    source_ids <- function(x) unique(unlist(strsplit(x,';',fixed=TRUE)))
    overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]),source_ids(obs$source_sample_ids_AB[te]))
    stopifnot(length(overlap)==0L,!held %chin% fit$training_partners)
    out <- obs[,.(observation_id,condition_id,anchor,partner,nominal_A_uM,nominal_B_uM,biological_replicate)]
    out[,`:=`(time_h=3L,fold_id=held,model=variant,split=fifelse(partner==held,'held','training'),
      selected_overall=chosen==overall,predicted_single_sum_log2=offset,
      predicted_interaction_log2=pred,predicted_ratio_AB=2^(offset+pred))]
    output_name <- paste0('prediction_',held,'_',variant)
    write_table(out,output_name)
    key <- paste(held,variant)
    audits[[key]] <- data.table(fold_id=held,model=variant,source_overlap=length(overlap),
      replay_error=max(abs(pred-replay)),active_terms=sum(raw!=0),alpha=fit$alpha,lambda=fit$selected_lambda,
      retained_simple_terms=length(fit$retained_simple_terms),available_terms=sum(fit$keep),
      selected_overall=chosen==overall,prediction_sha256=digest(file.path(root,'data',paste0(output_name,'.tsv.gz')),file=TRUE,algo='sha256'))
    coefficients[[key]] <- data.table(fold_id=held,model=variant,term=names(raw),coefficient=raw)
    out[,`:=`(ratio_A=obs$ratio_A,ratio_B=obs$ratio_B,ratio_AB=obs$ratio_AB)]
    predictions[[key]] <- out
  }
  current <- rbindlist(predictions)
  write_table(current,'predictions')
  write_table(rbindlist(selections),'selection')
  write_table(rbindlist(audits),'validation')
  write_table(rbindlist(coefficients),'coefficients')
  means <- current[split=='held',.(observed=mean(log2(ratio_AB/(ratio_A*ratio_B))),
    predicted=mean(predicted_interaction_log2)),by=.(fold_id,model,condition_id)]
  quick <- means[,.(n=.N,held_r=safe_cor(observed,predicted)),by=.(fold_id,model)]
  write_table(quick,'held_correlations_by_partner')
  print(quick[fold_id==held])
  status(paste('saved',held,'|',tuning$variant[overall],'selected by training exclusions'))
}
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
status('Bootstrap-screened interaction fits complete; held correlations saved.')
