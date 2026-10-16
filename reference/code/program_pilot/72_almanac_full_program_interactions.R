#!/usr/bin/env Rscript
# Full symmetric original40 interactions, bootstrap screened within drug exclusions.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(digest)})
root <- '/work/results/program_pilot/72_almanac_full_program_interactions'
script_path <- '/work/src/program_pilot/72_almanac_full_program_interactions.R'
input_path <- file.path(root,'data/prepared_inputs.qs')
lambda_grid <- 10^seq(1.5,-3.5,length.out=32L)
alpha_grid <- c(.1,.5,1)
bootstrap_draws <- 500L
seed <- 72011L
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn=1)
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv.gz')),sep='\t')
status <- function(x) {
  message(format(Sys.time(),tz='UTC',usetz=TRUE),' | ',x)
  writeLines(x,file.path(root,'STATUS.txt'))
}
options(error=function() {
  status('FAILED; see model log')
  traceback(5)
  quit(save='no',status=1L)
})
input <- qread(input_path)
d <- input$metadata
refs <- input$reference_map
rms <- sqrt(colMeans(input$query_scores^2))
stopifnot(all(rms>0),ncol(input$A)==40L,!anyDuplicated(d$combination_id))
A <- sweep(input$A,2,rms,'/')
B <- sweep(input$B,2,rms,'/')
ij <- which(upper.tri(matrix(0,40,40),diag=TRUE),arr.ind=TRUE)
X <- vapply(seq_len(nrow(ij)),function(k) {
  i <- ij[k,1]
  j <- ij[k,2]
  (A[,i]*B[,j]+A[,j]*B[,i])/2
},numeric(nrow(d)))
colnames(X) <- paste(colnames(A)[ij[,1]],colnames(A)[ij[,2]],sep='__')
stopifnot(ncol(X)==820L,all(is.finite(X)),all(is.finite(d$observed_g)))
for (k in c(1L,410L,820L)) {
  i <- ij[k,1]
  j <- ij[k,2]
  stopifnot(max(abs(X[,k]-(B[,i]*A[,j]+B[,j]*A[,i])/2))<1e-12)
}
drugs <- sort(unique(c(d$drug_a,d$drug_b)))
run_signature <- digest(list(input=digest(input_path,file=TRUE,algo='sha256'),
  script=digest(script_path,file=TRUE,algo='sha256'),lambda_grid=lambda_grid,alpha_grid=alpha_grid,
  bootstrap_draws=bootstrap_draws,seed=seed),algo='sha256')
qsave(list(rms=rms,lambda_grid=lambda_grid,alpha_grid=alpha_grid,bootstrap_draws=bootstrap_draws,
  bootstrap_interval=c(.05,.95),screen='weighted single-term slope with intercept',
  objective='pair/recipient-balanced log interaction',run_signature=run_signature,
  weight_sha256=input$weight_sha256),file.path(root,'models/settings.qs'))
rm(input,A,B)

# Equal mass per pair, equal recipient mass within each pair, equal observations
# within pair/recipient. Dose trajectories are retained as individual conditions.
balanced_weights <- function(md) {
  group <- paste(md$pair,md$depmap_id,sep='|')
  n <- as.numeric(table(group)[group])
  cells <- md[,.(cells=uniqueN(depmap_id)),by=pair]
  w <- 1/n/cells$cells[match(md$pair,cells$pair)]
  w/sum(w)
}
getrefs <- function(rows) refs[row %in% rows,unique(ref)]
purge <- function(tr,te) setdiff(tr,refs[ref %chin% getrefs(te),unique(row)])
safe_cor <- function(x,y) {
  if (length(x)<3L || sd(x)<1e-10 || sd(y)<1e-10) return(NA_real_)
  cor(x,y)
}

# Reuse71's exact weighted univariate slope bootstrap with pair blocks instead
# of the single-anchor panel's partner blocks. Inner screening sees training only.
bootstrap_screen <- function(tr,screen_seed,label) {
  x <- X[tr,,drop=FALSE]
  y <- d$observed_g[tr]
  w <- balanced_weights(d[tr])
  units <- sort(unique(d$pair[tr]))
  group <- factor(d$pair[tr],levels=units)
  sw <- as.numeric(rowsum(w,group,reorder=FALSE)[units,,drop=FALSE])
  sx <- rowsum(x*w,group,reorder=FALSE)[units,,drop=FALSE]
  sxx <- rowsum(x^2*w,group,reorder=FALSE)[units,,drop=FALSE]
  sxy <- rowsum(x*(w*y),group,reorder=FALSE)[units,,drop=FALSE]
  sy <- as.numeric(rowsum(w*y,group,reorder=FALSE)[units,,drop=FALSE])
  set.seed(screen_seed)
  frequency <- t(replicate(bootstrap_draws,
    tabulate(sample.int(length(units),length(units),replace=TRUE),nbins=length(units))))
  total_w <- as.numeric(frequency%*%sw)
  sum_x <- frequency%*%sx
  sum_xx <- frequency%*%sxx
  denominator <- sum_xx-sum_x^2/total_w
  slopes <- (frequency%*%sxy-sum_x*as.numeric(frequency%*%sy)/total_w)/denominator
  slopes[denominator<=pmax(1e-14,abs(sum_xx)*1e-12)] <- NA_real_
  for (column in c(1L,410L,820L)) {
    checked <- lm.wfit(cbind(1,x[,column]),y,w=w*frequency[1,as.integer(group)])$coefficients[2]
    if (is.finite(slopes[1,column])) stopifnot(abs(checked-slopes[1,column])<1e-8*(1+abs(checked)))
  }
  bounds <- apply(slopes,2,function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(c(NA_real_,NA_real_))
    quantile(z,c(.05,.95),names=FALSE)
  })
  result <- data.table(term=colnames(X),lower=bounds[1,],upper=bounds[2,],
    median_slope=apply(slopes,2,median,na.rm=TRUE),positive_fraction=colMeans(slopes>0,na.rm=TRUE),
    finite_fraction=colMeans(is.finite(slopes)))
  result[, retained := finite_fraction>=.95 & is.finite(lower) & is.finite(upper) & (lower>0 | upper<0)]
  result[, `:=`(training_conditions=length(tr),training_pairs=length(units),screen=label)]
  write_table(result,paste0('screen_',label))
  result
}
prepare_fit <- function(tr,screen) {
  w <- balanced_weights(d[tr])
  sx <- sqrt(colSums(X[tr,,drop=FALSE]^2*w))
  keep <- screen$retained & sx>1e-10
  x <- sweep(X[tr,keep,drop=FALSE],2,sx[keep],'/')
  padded <- sum(keep)==1L
  if (padded) x <- cbind(x,inert_zero=0)
  sy <- sqrt(sum(w*d$observed_g[tr]^2))
  stopifnot(sy>0)
  list(x=x,y=d$observed_g[tr]/sy,w=w,sx=sx[keep],sy=sy,keep=keep,padded=padded)
}
fit_path <- function(prepared,alpha) {
  if (!any(prepared$keep)) return(NULL)
  fit <- glmnet(prepared$x,prepared$y,weights=prepared$w,alpha=alpha,lambda=lambda_grid,
    standardize=FALSE,intercept=FALSE,thresh=1e-8,maxit=100000L)
  stopifnot(fit$jerr==0L,length(fit$lambda)==length(lambda_grid))
  fit
}
predict_path <- function(fit,prepared,rows) {
  if (is.null(fit)) return(matrix(0,length(rows),length(lambda_grid)))
  x <- sweep(X[rows,prepared$keep,drop=FALSE],2,prepared$sx,'/')
  if (prepared$padded) x <- cbind(x,inert_zero=0)
  as.matrix(predict(fit,newx=x,s=lambda_grid))*prepared$sy
}

audits <- coefficients <- quick <- split_audits <- list()
for (ff in seq_along(drugs)) {
  fit <- fit_data <- NULL
  held <- drugs[ff]
  tag <- sprintf('%02d_%s',ff,gsub('[^A-Za-z0-9]+','_',held))
  model_path <- file.path(root,'models',paste0('exclude_',tag,'.qs'))
  te <- which(d$drug_a==held | d$drug_b==held)
  tr <- purge(setdiff(seq_len(nrow(d)),te),te)
  stopifnot(length(tr)>0L,length(te)>0L,!held %chin% c(d$drug_a[tr],d$drug_b[tr]),
    length(intersect(getrefs(tr),getrefs(te)))==0L)
  status(sprintf('%d/%d exclude %s | %d training conditions, %d pairs; %d held conditions',
    ff,length(drugs),held,length(tr),uniqueN(d$pair[tr]),length(te)))
  if (file.exists(model_path)) {
    saved <- qread(model_path)
    stopifnot(identical(saved$run_signature,run_signature),identical(saved$training_rows,tr),identical(saved$held_rows,te))
    raw_beta <- saved$raw_coefficients
    selection <- saved$selection
    screen <- saved$screen
    message('Reused verified completed fold: ',held)
  } else {
    # Five single-drug inner exclusions selected by coverage rank, outcome-free.
    available <- sort(unique(c(d$drug_a[tr],d$drug_b[tr])))
    counts <- vapply(available,function(u) sum(d$drug_a[tr]==u | d$drug_b[tr]==u),integer(1))
    ranked <- available[order(counts,available)]
    inner_drugs <- ranked[unique(round(seq(1,length(ranked),length.out=5)))]
    losses <- array(NA_real_,c(length(inner_drugs),length(lambda_grid),length(alpha_grid)))
    zero_losses <- numeric(length(inner_drugs))
    for (ii in seq_along(inner_drugs)) {
      inner_held <- inner_drugs[ii]
      iv <- tr[d$drug_a[tr]==inner_held | d$drug_b[tr]==inner_held]
      it <- purge(setdiff(tr,iv),iv)
      stopifnot(length(it)>0L,length(iv)>0L,
        !any(c(held,inner_held) %chin% c(d$drug_a[it],d$drug_b[it])),
        length(intersect(getrefs(it),getrefs(iv)))==0L)
      inner_screen <- bootstrap_screen(it,seed+100L*ff+ii,paste0(tag,'_inner',ii))
      fit_data <- prepare_fit(it,inner_screen)
      wv <- balanced_weights(d[iv])
      for (aa in seq_along(alpha_grid)) {
        fit <- fit_path(fit_data,alpha_grid[aa])
        predicted <- predict_path(fit,fit_data,iv)
        losses[ii,,aa] <- colSums((predicted-d$observed_g[iv])^2*wv)
      }
      zero_losses[ii] <- sum(d$observed_g[iv]^2*wv)
      split_audits[[paste(ff,ii)]] <- data.table(held_drug=held,inner_drug=inner_held,
        training_conditions=length(it),training_pairs=uniqueN(d$pair[it]),validation_conditions=length(iv),
        shared_reference_groups=0L,retained_terms=sum(inner_screen$retained))
    }
    score <- apply(losses,c(2,3),mean)
    choice <- which(score==min(score),arr.ind=TRUE)[1,]
    zero <- mean(zero_losses)<=min(score)
    selection <- list(alpha=alpha_grid[choice[2]],lambda=lambda_grid[choice[1]],
      lambda_index=choice[1],zero=zero,inner_drugs=inner_drugs)
    tuning <- data.table(lambda=rep(lambda_grid,length(alpha_grid)),
      alpha=rep(alpha_grid,each=length(lambda_grid)),internal_objective=as.vector(score),
      zero_objective=mean(zero_losses))
    write_table(tuning,paste0('tuning_',tag))
    screen <- bootstrap_screen(tr,seed+100L*ff,paste0(tag,'_outer'))
    fit_data <- prepare_fit(tr,screen)
    fit <- fit_path(fit_data,selection$alpha)
    raw_beta <- setNames(rep(0,ncol(X)),colnames(X))
    if (!zero && !is.null(fit)) {
      beta <- as.numeric(coef(fit,s=selection$lambda))[-1]
      if (fit_data$padded) beta <- beta[-length(beta)]
      raw_beta[fit_data$keep] <- beta*fit_data$sy/fit_data$sx
    }
    qsave(list(raw_coefficients=raw_beta,fit=fit,screen=screen,selection=selection,
      training_rows=tr,held_rows=te,held_drug=held,run_signature=run_signature,rms=rms),model_path)
    if (!zero && !is.null(fit)) {
      direct <- predict_path(fit,fit_data,te)[,selection$lambda_index]
      stopifnot(max(abs(direct-as.numeric(X[te,,drop=FALSE]%*%raw_beta)))<1e-8)
    }
  }
  use <- c(tr,te)
  prediction <- as.numeric(X[use,,drop=FALSE]%*%raw_beta)
  stopifnot(all(is.finite(prediction)),all(raw_beta[!screen$retained]==0))
  out <- d[use,.(combination_id,pair,recipient,depmap_id,SCREENER,drug_a,drug_b,dose_a,dose_b)]
  out[, `:=`(held_drug=held,split=c(rep('training',length(tr)),rep('held',length(te))),predicted_g=prediction)]
  write_table(out,paste0('prediction_',tag))
  prediction_hash <- digest(file.path(root,'data',paste0('prediction_',tag,'.tsv.gz')),file=TRUE,algo='sha256')
  audits[[tag]] <- data.table(held_drug=held,training_conditions=length(tr),training_pairs=uniqueN(d$pair[tr]),
    held_conditions=length(te),held_pairs=uniqueN(d$pair[te]),purged_conditions=nrow(d)-length(use),
    retained_terms=sum(screen$retained),active_terms=sum(raw_beta!=0),alpha=selection$alpha,
    lambda=selection$lambda,zero=selection$zero,shared_reference_groups=0L,prediction_sha256=prediction_hash)
  coefficients[[tag]] <- data.table(held_drug=held,term=names(raw_beta),coefficient=raw_beta)
  quick[[tag]] <- data.table(held_drug=held,held_conditions=length(te),
    training_interaction_r=safe_cor(d$observed_g[tr],prediction[seq_along(tr)]),
    held_interaction_r=safe_cor(d$observed_g[te],tail(prediction,length(te))))
  write_table(rbindlist(audits),'validation')
  write_table(rbindlist(coefficients),'coefficients')
  write_table(rbindlist(quick),'correlations_by_drug')
  if (length(split_audits)) write_table(rbindlist(split_audits),'inner_validation')
  print(quick[[tag]])
  rm(fit,fit_data)
  invisible(gc(FALSE))
}
files <- list.files(file.path(root,'data'),pattern='^prediction_[0-9].*tsv.gz$',full.names=TRUE)
all <- rbindlist(lapply(files,fread))
held <- all[split=='held']
stopifnot(length(files)==length(drugs),nrow(held)==2*nrow(d),
  held[, .N,by=combination_id][,all(N==2L)])
write_table(all,'predictions')
write_table(data.table(path=c(input_path,script_path),
  sha256=vapply(c(input_path,script_path),digest,character(1),file=TRUE,algo='sha256')),'run_manifest')
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
writeLines(system2('git',c('-C','/work','rev-parse','HEAD'),stdout=TRUE),file.path(root,'git_revision.txt'))
status('COMPLETE: all40 drug exclusions saved; full820 interaction model, bootstrap screening and exact reference purging')
