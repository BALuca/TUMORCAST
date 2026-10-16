#!/usr/bin/env Rscript
# Existing two signatures and training-only route assignments are fixed.
# Eight shared coefficients predict molecular-footprint level and positive
# amplitude. No per-drug viability coefficients or outcome-fitted program gains.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(glmnet);library(ggplot2);library(digest)})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

base_root <- '/work/results/program_pilot/104_erbb2_met_lapatinib_shared_decoder'
root <- file.path(base_root,'linear')
parent <- '/work/results/program_pilot/95_leave_one_cell_out_correlation_signatures'
source_root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data'
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
programs <- sprintf('P%02d',1:40)
signature_amplitude <- .25
calibration_penalty <- .005
assay_penalty <- .01
interaction_penalty <- .05
excluded_drugs <- character()
for (folder in c('data','models','figures')) dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
status <- function(message) {
  cat(format(Sys.time()),message,'\n')
  flush.console()
  writeLines(message,file.path(root,'STATUS.txt'))
}
pearson <- function(x,y) if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
key <- function(cell,drug) paste(cell,drug,sep='|')
selection <- fread(file.path(source_root,'selected_doses.tsv'))
single <- fread(file.path(source_root,'individual_training_outcomes.tsv'))
pairs <- fread(file.path(source_root,'combination_training_outcomes.tsv'))
parent_folds <- fread(file.path(parent,'data/recipient_folds.tsv'))
bank_selection <- copy(selection)
drugs <- bank_selection$drug
selection <- selection[!drug %in% excluded_drugs]
included_indices <- match(selection$drug,drugs)
single <- single[!drug %in% excluded_drugs]
pairs <- pairs[!drugA %in% excluded_drugs & !drugB %in% excluded_drugs]
short <- c('Azacytidine','Bendamustine','Bortezomib','Busulfan','Cytarabine',
           'Lapatinib','Topotecan','Vandetanib','Vismodegib','Crizotinib')
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
stopifnot(!anyDuplicated(receptors$depmap_id))
molecular <- merge(as.data.table(archive$prep$m),bank_selection,by=c('drug','dose'))[depmap_id!='ACH-000248']
prism <- single[assay=='secondary']
nci <- single[assay!='secondary' & depmap_id!='ACH-000248']
profiles <- merge(molecular,prism[,.(depmap_id,drug,dose,S=-y)],by=c('depmap_id','drug','dose'))
rm(archive)
stopifnot(nrow(profiles)>0L,nrow(pairs)>0L,
          !anyDuplicated(profiles[,.(depmap_id,drug)]),!anyDuplicated(baseline$depmap_id),
          all(is.finite(as.matrix(profiles[,..programs]))),all(is.finite(single$y)))
loo_cells <- sort(unique(c(profiles$depmap_id,nci$depmap_id,pairs$depmap_id)))
folds <- data.table(depmap_id=loo_cells,fold=seq_along(loo_cells))
write_table(selection,'selected_doses')
write_table(folds,'recipient_folds')

# Lapatinib uses ERBB2/HER2 and MET only, through predicted displacements.
# Other drugs retain the v95 RNA-only map. Nonmatched RNA recipients need
# the same RNA fit with their own cell removed; no functional labels enter it.
fit_rna <- function(train_rna) {
  ids <- sort(unique(train_rna$depmap_id))
  raw <- as.matrix(baseline[match(ids,depmap_id),..programs])
  bmean <- colMeans(raw)
  bsd <- apply(raw,2,sd)
  normalized <- sweep(sweep(raw,2,bmean),2,bsd,'/')
  basis <- svd(normalized,nu=0,nv=8)$v
  bscale <- sqrt(colMeans((normalized %*% basis)^2))
  B <- array(0,c(8,40,10))
  mu <- sigma <- matrix(0,10,40)
  bm <- matrix(0,10,8)
  for (g in 1:10) {
    rows <- which(train_rna$drug==drugs[g])
    D <- as.matrix(train_rna[rows,..programs])
    mu[g,] <- colMeans(D)
    sigma[g,] <- apply(D,2,sd)
    features <- sweep(normalized[match(train_rna$depmap_id[rows],ids),,drop=FALSE] %*% basis,2,bscale,'/')
    bm[g,] <- colMeans(features)
    features <- sweep(features,2,bm[g,])
    target <- sweep(sweep(D,2,mu[g,]),2,sigma[g,],'/')
    fit <- glmnet(features,target,family='mgaussian',alpha=0,lambda=1,intercept=FALSE,standardize=FALSE)
    B[,,g] <- vapply(coef(fit),function(x)as.numeric(x)[-1],numeric(8))
  }
  list(bmean=bmean,bsd=bsd,basis=basis,bscale=bscale,B=B,mu=mu,sigma=sigma,bm=bm,training_cells=ids)
}
fit_lapatinib_receptors <- function(train_rna) {
  train <- train_rna[drug=="Lapatinib ditosylate"]
  X <- as.matrix(receptors[match(train$depmap_id,depmap_id),.(ERBB2,MET)])
  D <- as.matrix(train[,..programs])
  xm <- colMeans(X); xs <- apply(X,2,sd)
  dm <- colMeans(D); ds <- pmax(apply(D,2,sd),1e-8)
  X <- sweep(sweep(X,2,xm),2,xs,"/")
  Y <- sweep(sweep(D,2,dm),2,ds,"/")
  fit <- glmnet(X,Y,family="mgaussian",alpha=0,lambda=1,intercept=FALSE,standardize=FALSE)
  B <- vapply(coef(fit),function(x)as.numeric(x)[-1],numeric(2))
  list(xmean=xm,xsd=xs,mu=dm,sigma=ds,B=B,training_cells=train$depmap_id,lambda=1)
}
predict_rna <- function(model,cells,drug) {
  raw <- as.matrix(baseline[match(cells,depmap_id),..programs])
  features <- sweep(sweep(sweep(raw,2,model$bmean),2,model$bsd,'/') %*% model$basis,2,model$bscale,'/')
  g <- match(drug,drugs)
  D <- matrix(0,length(cells),40)
  for (j in unique(g)) {
    rows <- which(g==j)
    Z <- sweep(features[rows,,drop=FALSE],2,model$bm[j,]) %*% model$B[,,j]
    D[rows,] <- sweep(sweep(Z,2,model$sigma[j,],'*'),2,model$mu[j,],'+')
  }
  rows <- which(drug=="Lapatinib ditosylate")
  if(length(rows)) {
    lap <- model$lap_receptors
    stopifnot(!is.null(lap))
    X <- as.matrix(receptors[match(cells[rows],depmap_id),.(ERBB2,MET)])
    Z <- sweep(sweep(X,2,lap$xmean),2,lap$xsd,"/") %*% lap$B
    D[rows,] <- sweep(sweep(Z,2,lap$sigma,"*"),2,lap$mu,"+")
  }
  stopifnot(all(is.finite(D)))
  D
}

# The same three RNA-only footprint coordinates used by v100/v102 supply all
# shared calibration functions. Route labels are frozen training correlations.
make_context <- function(template) {
  fm <- colMeans(template$mu)
  fs <- pmax(apply(template$mu,2,sd),1e-8)
  F <- sweep(sweep(template$mu,2,fm),2,fs,'/')
  fv <- svd(F,nu=0,nv=3)$v
  fscale <- sqrt(colMeans((F %*% fv)^2))
  Q <- cbind(1,sweep(F %*% fv,2,fscale,'/'))
  list(mu=template$mu,sigma=template$sigma,Q=Q,W=template$W*signature_amplitude,
       route=template$membership,fm=fm,fs=fs,fv=fv,fscale=fscale)
}
signature_node <- function(D,g,ctx) {
  X <- (D-ctx$mu[g,,drop=FALSE])/ctx$sigma[g,,drop=FALSE]
  both <- X %*% ctx$W
  route <- ctx$route[g]
  score <- both[cbind(seq_len(nrow(D)),route)]
  # Pair interactions reference zero displacement. Keep each constituent
  # separate and populate only its assigned route; no summed RNA displacement.
  absolute <- (D/ctx$sigma[g,,drop=FALSE]) %*% ctx$W
  activity <- matrix(0,nrow(D),2)
  activity[cbind(seq_len(nrow(D)),route)] <- absolute[cbind(seq_len(nrow(D)),route)]
  list(score=score,both=both,activity=activity,route=route)
}
pair_features <- function(A,B) cbind(A[,1]*B[,1],A[,1]*B[,2]+A[,2]*B[,1],A[,2]*B[,2])
geometry <- function(y,group,control=NULL) {
  g <- match(group,unique(group))
  w <- rep(1,length(y))
  if(!is.null(control)) w <- 1/as.numeric(table(control)[control])
  totals <- as.numeric(rowsum(w,g,reorder=FALSE))
  w <- w/totals[g]/length(totals)
  gw <- as.numeric(rowsum(w,g,reorder=FALSE))
  ym <- as.numeric(rowsum(y*w,g,reorder=FALSE))/gw
  list(y=y,w=w,within=max(sum(w*(y-ym[g])^2),1e-8))
}
loss <- function(pred,geo) {
  error <- pred-geo$y
  list(value=sum(geo$w*error^2)/geo$within,gradient=2*geo$w*error/geo$within)
}
unpack <- function(theta) list(level=theta[1:4],log_amplitude=theta[5:8],assay=theta[9:10],Gamma=theta[11:13])
pack <- function(par) unlist(par,use.names=FALSE)
prepare <- function(held,template,rna) {
  train <- profiles[!depmap_id %in% held]
  ns <- nci[!depmap_id %in% held]
  pp <- pairs[!depmap_id %in% held]
  query <- unique(rbind(ns[,.(depmap_id,drug)],pp[,.(depmap_id,drug=drugA)],pp[,.(depmap_id,drug=drugB)]))
  qkey <- key(query$depmap_id,query$drug)
  ctx <- make_context(template)
  gobs <- match(train$drug,drugs)
  gquery <- match(query$drug,drugs)
  obs <- signature_node(as.matrix(train[,..programs]),gobs,ctx)
  native <- signature_node(predict_rna(rna,query$depmap_id,query$drug),gquery,ctx)
  is <- match(key(ns$depmap_id,ns$drug),qkey)
  ia <- match(key(pp$depmap_id,pp$drugA),qkey)
  ib <- match(key(pp$depmap_id,pp$drugB),qkey)
  products <- pair_features(native$activity[ia,,drop=FALSE],native$activity[ib,,drop=FALSE])
  ctx$interaction_scale <- sqrt(colMeans(products^2))
  ctx$interaction_supported <- ctx$interaction_scale>1e-10
  # Unrepresented route pairs cannot identify an interaction coefficient.
  # Unit scaling for their zero training column permits an explicitly fixed 0.
  ctx$interaction_scale[!ctx$interaction_supported] <- 1
  stopifnot(any(ctx$interaction_supported),max(abs(tapply(obs$score,train$drug,mean)))<1e-10)
  list(ctx=ctx,train=train,ns=ns,pp=pp,obs=obs,native=native,
       Qobs=ctx$Q[gobs,,drop=FALSE],Qquery=ctx$Q[gquery,,drop=FALSE],is=is,ia=ia,ib=ib,
       products=sweep(products,2,ctx$interaction_scale,'/'),
       assay=match(ns$assay,c('nci_FF','nci_FG')),pair_assay=match(pp$SCREENER,c('FF','FG')),
       geoO=geometry(train$S,train$drug),geoN=geometry(-ns$y,ns$group,ns$control_key),
       geoP=geometry(-pp$y,pp$pair),geoI=geometry(-pp$I,pp$pair),
       training_cells=unique(c(train$depmap_id,ns$depmap_id,pp$depmap_id)))
}
fit_decoder <- function(p,template,job) {
  Qactive <- p$ctx$Q[included_indices,,drop=FALSE]
  initial_level <- as.numeric(solve(crossprod(Qactive)/nrow(Qactive)+diag(calibration_penalty,4),
                                   crossprod(Qactive,template$reference[included_indices])/nrow(Qactive)))
  initial <- c(initial_level,rep(0,4),-1,-1,rep(0,3))
  cache_theta <- cache_value <- NULL
  objective <- function(theta) {
    if(identical(theta,cache_theta)) return(cache_value)
    par <- unpack(theta)
    ao <- as.numeric(exp(p$Qobs %*% par$log_amplitude))
    aq <- as.numeric(exp(p$Qquery %*% par$log_amplitude))
    po <- as.numeric(p$Qobs %*% par$level)+ao*p$obs$score
    native <- as.numeric(p$Qquery %*% par$level)+aq*p$native$score
    scales <- exp(par$assay)
    pn <- native[p$is]*scales[p$assay]
    # Positive amplitude changes interaction magnitude, not its route assignment.
    F <- p$products*(aq[p$ia]*aq[p$ib])
    J <- as.numeric(F %*% par$Gamma)
    pp <- (native[p$ia]+native[p$ib]+J)*scales[p$pair_assay]
    pi <- J*scales[p$pair_assay]
    lo <- loss(po,p$geoO)
    ln <- loss(pn,p$geoN)
    lp <- loss(pp,p$geoP)
    li <- loss(pi,p$geoI)
    go <- lo$gradient
    gn <- .2*ln$gradient*scales[p$assay]
    gp <- .2*lp$gradient*scales[p$pair_assay]
    gj <- .2*(lp$gradient+li$gradient)*scales[p$pair_assay]
    QA <- p$Qquery[p$ia,,drop=FALSE]
    QB <- p$Qquery[p$ib,,drop=FALSE]
    QN <- p$Qquery[p$is,,drop=FALSE]
    level_gradient <- as.numeric(crossprod(p$Qobs,go)+crossprod(QN,gn)+crossprod(QA+QB,gp))
    amplitude_gradient <- as.numeric(crossprod(p$Qobs,go*ao*p$obs$score)+
      crossprod(QN,gn*aq[p$is]*p$native$score[p$is])+
      crossprod(QA,gp*aq[p$ia]*p$native$score[p$ia])+
      crossprod(QB,gp*aq[p$ib]*p$native$score[p$ib])+crossprod(QA+QB,gj*J))
    assay_gradient <- vapply(1:2,function(a)
      .2*sum(ln$gradient[p$assay==a]*pn[p$assay==a])+
      .2*sum(lp$gradient[p$pair_assay==a]*pp[p$pair_assay==a])+
      .2*sum(li$gradient[p$pair_assay==a]*pi[p$pair_assay==a]),numeric(1))
    penalty <- calibration_penalty*(sum(par$level^2)+sum(par$log_amplitude^2))+
      assay_penalty*sum(par$assay^2)+interaction_penalty*sum(par$Gamma^2)
    gradient <- pack(list(level=level_gradient+2*calibration_penalty*par$level,
      log_amplitude=amplitude_gradient+2*calibration_penalty*par$log_amplitude,
      assay=assay_gradient+2*assay_penalty*par$assay,
      Gamma=as.numeric(crossprod(F,gj))+2*interaction_penalty*par$Gamma))
    value <- list(value=lo$value+.2*ln$value+.2*lp$value+.2*li$value+penalty,gradient=gradient)
    stopifnot(is.finite(value$value),all(is.finite(value$gradient)))
    cache_theta <<- theta
    cache_value <<- value
    value
  }
  if(job==1L) {
    probe <- initial+c(.01,-.01,.02,.01,.03,-.02,.01,.02,.01,-.01,.01,-.02,.01)
    analytic <- objective(probe)$gradient
    numeric <- vapply(seq_along(probe),function(i) {
      plus <- minus <- probe
      plus[i] <- plus[i]+1e-5
      minus[i] <- minus[i]-1e-5
      (objective(plus)$value-objective(minus)$value)/2e-5
    },numeric(1))
    difference <- abs(analytic-numeric)/pmax(1,abs(analytic),abs(numeric))
    write_table(data.table(parameter=1:13,analytic,numeric,relative_difference=difference),'analytic_gradient_check')
    stopifnot(max(difference)<1e-5)
  }
  lower <- -c(rep(20,4),rep(2,4),rep(4,2),rep(10,3))
  upper <- -lower
  unsupported_indices <- 10L+which(!p$ctx$interaction_supported)
  lower[unsupported_indices] <- 0
  upper[unsupported_indices] <- 0
  opt <- optim(pmin(pmax(initial,lower),upper),function(t)objective(t)$value,function(t)objective(t)$gradient,
               method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=1000,factr=1e8,pgtol=1e-5))
  list(par=unpack(opt$par),context=p$ctx,training_cells=p$training_cells,
       convergence=opt$convergence,evaluations=opt$counts[1],max_gradient=max(abs(objective(opt$par)$gradient)))
}

predict_single <- function(model,meta,measured_D=NULL) {
  D <- if(is.null(measured_D)) predict_rna(model$rna,meta$depmap_id,meta$drug) else measured_D
  g <- match(meta$drug,drugs)
  node <- signature_node(D,g,model$context)
  Q <- model$context$Q[g,,drop=FALSE]
  level <- as.numeric(Q %*% model$par$level)
  amplitude <- as.numeric(exp(Q %*% model$par$log_amplitude))
  sensitivity <- level+amplitude*node$score
  assay <- if('assay' %in% names(meta)) match(meta$assay,c('secondary','nci_FF','nci_FG')) else rep(1L,nrow(meta))
  scale <- c(1,exp(model$par$assay))[assay]
  stopifnot(all(amplitude>0),!anyNA(scale))
  list(sensitivity=sensitivity*scale,score=node$score,route=node$route,level=level,amplitude=amplitude,
       activity=node$activity*amplitude,D=D)
}
predict_pair <- function(model,meta) {
  A <- predict_single(model,meta[,.(depmap_id,drug=drugA)])
  B <- predict_single(model,meta[,.(depmap_id,drug=drugB)])
  products <- sweep(pair_features(A$activity,B$activity),2,model$context$interaction_scale,'/')
  J <- as.numeric(products %*% model$par$Gamma)
  scale <- exp(model$par$assay)[match(meta$SCREENER,c('FF','FG'))]
  list(sensitivity=(A$sensitivity+B$sensitivity+J)*scale,departure=J*scale)
}
render_heatmap <- function(template,label) {
  ro <- hclust(dist(template$correlations),method='ward.D2')$order
  co <- template$drug_tree$order
  co <- co[co %in% included_indices]
  annotation <- data.frame(Route=paste0('Signature ',template$membership[co]),
                           Dose=paste0(bank_selection$dose[co],' µM'),row.names=drugs[co])
  stopifnot(identical(colnames(template$correlations[,co]),rownames(annotation)))
  left <- draw_heatmap(template$correlations[ro,co,drop=FALSE],scale_rows=FALSE,
    value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),hm_name='Training Pearson r',
    col_anno_df=annotation,col_anno_cols=c('Route','Dose'),
    palette_overrides=list(Route=c('Signature 1'='#D58B36','Signature 2'='#775DA6'),
                          Dose=c('0.05 µM'='#DCEAF2','0.5 µM'='#6BAED6','5 µM'='#08519C')),
    show_row_names=TRUE,row_names_side='left',row_names_gp=grid::gpar(fontsize=8),
    show_column_names=TRUE,column_labels=short[co],column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
    width=grid::unit(105,'mm'),title='Observed training correlations')
  right <- draw_heatmap(template$collapsed[ro,,drop=FALSE],scale_rows=FALSE,
    value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),hm_name='Collapsed mean r',
    show_row_names=FALSE,show_column_names=TRUE,column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
    width=grid::unit(24,'mm'),title='Two fixed signatures')
  ComplexHeatmap::draw(left+right,heatmap_legend_side='bottom',annotation_legend_side='right',
    column_title=paste('v104 |',label),padding=grid::unit(c(15,5,5,5),'mm'),ht_gap=grid::unit(8,'mm'))
  grid::grid.text('Route and program weights fixed from training cells. Eight shared footprint-calibration coefficients.\nPositive amplitude preserves the signature direction; calibration receives no drug identity or baseline features.',
    x=.5,y=grid::unit(4,'mm'),just=c('centre','bottom'),gp=grid::gpar(fontsize=8))
}

measured_predictions <- baseline_predictions <- combo_predictions <- audits <- calibrations <- list()
direction_checks <- coefficient_tables <- interaction_tables <- list()
final_parent <- qread(file.path(parent,'models/final.qs'))
book <- file.path(root,'figures/training_correlations_and_collapsed_signatures_LOO.pdf')
cairo_pdf(book,width=9,height=9,onefile=TRUE)
book_device <- dev.cur()
final_job <- length(loo_cells)+1L
status(paste('v104:',length(loo_cells),'whole-cell exclusions; fixed signatures/routes; shared positive calibration.'))
for (job in seq_len(final_job)) {
  held <- if(job==final_job) character() else loo_cells[job]
  label <- if(length(held)) sprintf('loo_%02d',job) else 'final'
  parent_index <- match(held,parent_folds$depmap_id)
  if(length(parent_index) && !is.na(parent_index)) {
    ancestor <- qread(file.path(parent,'models',sprintf('loo_%02d.qs',parent_index)))
    template <- ancestor$template
    rna <- ancestor$rna
  } else {
    template <- final_parent$template
    rna <- if(length(held) && held %in% molecular$depmap_id) fit_rna(molecular[depmap_id!=held]) else final_parent$rna
  }
  stopifnot(!any(held %in% template$training_cells),!any(held %in% rna$training_cells))
  rna$lap_receptors <- fit_lapatinib_receptors(molecular[!depmap_id %in% held])
  stopifnot(!any(held %in% rna$lap_receptors$training_cells))
  prepared <- prepare(held,template,rna)
  fit <- fit_decoder(prepared,template,job)
  model <- c(fit,list(template=template,rna=rna,held_cells=held,selected=selection))
  stopifnot(!any(held %in% model$training_cells),identical(model$context$W,template$W*signature_amplitude),
            identical(model$context$route,template$membership))
  qsave(model,file.path(root,'models',paste0(label,'.qs')))
  # Verify the agreed protection on every training treatment, not just parameter
  # signs: Pearson and centered score/prediction relationships must be preserved.
  training_result <- predict_single(model,prepared$train,as.matrix(prepared$train[,..programs]))
  training_checks <- data.table(drug=prepared$train$drug,observed=prepared$train$S,
                                score=training_result$score,predicted=training_result$sensitivity,
                                level=training_result$level,amplitude=training_result$amplitude)[,
    .(score_r=pearson(observed,score),prediction_r=pearson(observed,predicted),
      score_prediction_r=pearson(score,predicted),level=unique(level),amplitude=unique(amplitude)),by=drug]
  stopifnot(max(abs(training_checks$score_r-training_checks$prediction_r))<1e-10,
            max(abs(training_checks$score_prediction_r-1))<1e-10,all(training_checks$amplitude>0))
  direction_checks[[label]] <- training_checks[,fit:=label]
  evaluation <- if(length(held)) 'held_cell' else 'fitted'
  test <- if(length(held)) profiles[depmap_id==held] else profiles
  if(nrow(test)) {
    result <- predict_single(model,test,as.matrix(test[,..programs]))
    p <- test[,.(depmap_id,drug,dose,observed_log2=-S)]
    p[,`:=`(predicted_log2=-result$sensitivity,signature_score=result$score,route=result$route,
             calibration_level=result$level,positive_amplitude=result$amplitude,fold=job,evaluation=evaluation)]
    measured_predictions[[label]] <- p
  }
  test <- if(length(held)) prism[depmap_id==held] else prism
  if(nrow(test)) {
    result <- predict_single(model,test)
    p <- test[,.(depmap_id,drug,dose,observed_log2=y)]
    p[,`:=`(predicted_log2=-result$sensitivity,signature_score=result$score,route=result$route,
             calibration_level=result$level,positive_amplitude=result$amplitude,fold=job,evaluation=evaluation)]
    baseline_predictions[[label]] <- p
  }
  test <- if(length(held)) pairs[depmap_id==held] else pairs
  if(nrow(test)) {
    result <- predict_pair(model,test)
    p <- test[,.(depmap_id,pair,SCREENER,drugA,drugB,doseA,doseB,observed_log2=y,observed_departure_log2=I)]
    p[,`:=`(predicted_log2=-result$sensitivity,predicted_departure_log2=-result$departure,fold=job,evaluation=evaluation)]
    combo_predictions[[label]] <- p
  }
  calibration_level <- as.numeric(model$context$Q %*% model$par$level)
  calibration_amplitude <- as.numeric(exp(model$context$Q %*% model$par$log_amplitude))
  calibrations[[label]] <- data.table(fit=label,drug=drugs,route=model$context$route,
                                     predicted_level=calibration_level,positive_amplitude=calibration_amplitude)
  coefficient_tables[[label]] <- data.table(fit=label,term=rep(c('intercept','RNA_PC1','RNA_PC2','RNA_PC3'),2),
    function_name=rep(c('absolute_level','log_positive_amplitude'),each=4),
    coefficient=c(model$par$level,model$par$log_amplitude))
  interaction_tables[[label]] <- data.table(fit=label,term=c('route1_route1','cross_route','route2_route2'),
    training_supported=model$context$interaction_supported,
    product_RMS=model$context$interaction_scale,coefficient=model$par$Gamma)
  zeros <- signature_node(matrix(0,10,40),1:10,model$context)$activity
  stopifnot(max(abs(zeros))<1e-12)
  audits[[label]] <- data.table(fit=label,held_cell=if(length(held)) held else NA_character_,
    convergence=fit$convergence,evaluations=fit$evaluations,max_gradient=fit$max_gradient,
    functional_leakage=length(intersect(held,model$training_cells)),RNA_leakage=length(intersect(held,rna$training_cells)),
    minimum_positive_amplitude=min(calibration_amplitude),
    maximum_correlation_change=max(abs(training_checks$score_r-training_checks$prediction_r)))
  dev.set(book_device)
  render_heatmap(template,paste(label,if(length(held)) held else 'all-data fit'))
  if(job %in% c(1L,30L,final_job)) {
    png(file.path(root,'figures',paste0(label,'_correlations_and_signatures.png')),width=9,height=9,units='in',res=150,type='cairo')
    render_heatmap(template,label)
    dev.off()
  }
  write_table(rbindlist(audits),'fit_audit')
  if(job %% 10L==0L || job==final_job) status(paste('Fits completed:',job,'of',final_job))
}
dev.set(book_device)
dev.off()
measured <- rbindlist(measured_predictions)
replay <- rbindlist(baseline_predictions)
combo <- rbindlist(combo_predictions)
external <- copy(replay[evaluation=='fitted' & !depmap_id %in% loo_cells])
stopifnot(!any(external$depmap_id %in% model$training_cells),!any(external$depmap_id %in% model$rna$training_cells))
external[,evaluation:='held_cell']
replay <- rbind(replay,external)
stopifnot(nrow(measured[evaluation=='held_cell'])==nrow(profiles),
          nrow(replay[evaluation=='held_cell'])==nrow(prism),nrow(combo[evaluation=='held_cell'])==nrow(pairs),
          !anyDuplicated(replay[evaluation=='held_cell',.(depmap_id,drug)]))
write_table(measured,'measured_predictions_frozen')
write_table(replay,'baseline_predictions_frozen')
write_table(combo,'combination_predictions_frozen')
write_table(rbindlist(direction_checks),'signature_direction_protection_checks')
write_table(rbindlist(calibrations),'molecular_footprint_calibration_values')
write_table(rbindlist(coefficient_tables),'shared_calibration_coefficients')
write_table(rbindlist(interaction_tables),'interaction_feature_scales')

# HGF has no observed functional correlation pattern establishing its route.
# Interrogate both route assumptions, explicitly; do not choose a rescuing route
# from Wilson outcomes. Molecular footprint scale transport is unchanged.
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
hmean <- as.numeric(hgf[context=='Equal-context mean',..programs])
ctx <- model$context
Qh <- cbind(1,sweep(matrix((hmean-ctx$fm)/ctx$fs,nrow=1) %*% ctx$fv,2,ctx$fscale,'/'))
scale_map <- solve(crossprod(ctx$Q)/10+diag(c(0,1,1,1)),crossprod(ctx$Q,log(ctx$sigma))/10)
hscale_raw <- as.numeric(exp(Qh %*% scale_map))
hscale <- pmin(pmax(hscale_raw,apply(ctx$sigma,2,min)),apply(ctx$sigma,2,max))
lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate'))
hgf_queries <- list()
for (route in 1:2) {
  hctx <- ctx
  hctx$Q <- rbind(ctx$Q,Qh)
  hctx$mu <- rbind(ctx$mu,hmean)
  hctx$sigma <- rbind(ctx$sigma,hscale)
  hctx$route <- c(ctx$route,route)
  hnode <- signature_node(as.matrix(hgf[,..programs]),rep(11L,nrow(hgf)),hctx)
  hlevel <- as.numeric(Qh %*% model$par$level)
  hamplitude <- as.numeric(exp(Qh %*% model$par$log_amplitude))
  Sh <- hlevel+hamplitude*hnode$score
  LA <- lap$activity[rep(1L,nrow(hgf)),,drop=FALSE]
  HA <- hnode$activity*hamplitude
  F <- sweep(pair_features(LA,HA),2,ctx$interaction_scale,'/')
  J <- as.numeric(F %*% model$par$Gamma)
  hgf_queries[[route]] <- data.table(cell='AU565',dose=5,hgf_context=hgf$context,assumed_HGF_route=route,
    lapatinib_pct_vehicle=100*2^(-lap$sensitivity),hgf_alone_pct_vehicle=100*2^(-Sh),
    combination_pct_vehicle=100*2^(-lap$sensitivity-Sh-J),
    lapatinib_in_hgf_pct_hgf_control=100*2^(-lap$sensitivity-J),interaction_log2=-J,
    qualification='HGF route conditional, not identified; RNA-inferred SD; HGF/Wilson labels never fitted; no lapatinib-HGF training')
}
query <- rbindlist(hgf_queries)
write_table(query,'au565_hgf_5uM_predictions_frozen')
write_table(data.table(program=programs,inferred_SD=hscale_raw,used_SD=hscale),'hgf_scale_transfer')
paths <- file.path(root,'data',c('measured_predictions_frozen.tsv','baseline_predictions_frozen.tsv',
                               'combination_predictions_frozen.tsv','au565_hgf_5uM_predictions_frozen.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')

# Evaluation begins after predictions are frozen. Pearson compares cell lines
# receiving the same drug and dose, never pooled comparisons between drugs.
metrics <- measured[,.(cells=.N,r=pearson(observed_log2,predicted_log2),
  score_r=pearson(-observed_log2,signature_score),observed_sensitive=sum(observed_log2<= -1),
  predicted_sensitive=sum(predicted_log2<= -1),observed_median_pct=median(100*2^observed_log2),
  predicted_median_pct=median(100*2^predicted_log2)),by=.(drug,dose,evaluation)]
bridge_metrics <- replay[,.(cells=.N,r=pearson(observed_log2,predicted_log2),
  observed_sensitive=sum(observed_log2<= -1),predicted_sensitive=sum(predicted_log2<= -1),
  observed_median_pct=median(100*2^observed_log2),predicted_median_pct=median(100*2^predicted_log2)),by=.(drug,dose,evaluation)]
combo_metrics <- combo[,.(cells=.N,r_total=pearson(observed_log2,predicted_log2),
  r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_departure_SD=sd(observed_departure_log2),predicted_departure_SD=sd(predicted_departure_log2)),
  by=.(pair,SCREENER,evaluation)]
write_table(metrics,'measured_decoder_correlations')
write_table(bridge_metrics,'baseline_bridge_correlations')
write_table(combo_metrics,'combination_correlations')

plot_single <- function(table,heading,stem) {
  data <- copy(table[evaluation=='held_cell'])
  data[,panel:=factor(paste0(short[match(drug,drugs)],' | ',dose,' µM'),
                      levels=paste0(short[included_indices],' | ',selection$dose,' µM'))]
  stats <- data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed_log2,predicted_log2)),'\nn = ',.N)),by=panel]
  plot <- ggplot(data,aes(observed_log2,predicted_log2))+geom_point(color='#0072B2',alpha=.7,size=1.2)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
    geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=2.8)+
    facet_wrap(~panel,ncol=5,scales='free')+theme_bw(base_size=9)+theme(aspect.ratio=1,panel.grid.minor=element_blank())+
    labs(title=paste('v104 |',heading),x='Measured absolute log2 abundance',y='Predicted absolute log2 abundance',
      caption='Pearson across held-out cells at fixed drug/dose. 0 = 100% abundance; -1 = 50%.\nTwo fixed signatures/routes; eight shared calibration coefficients; positive amplitude.')
  ggsave(file.path(root,'figures',paste0(stem,'.pdf')),plot,width=13,height=6.8,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(stem,'.png')),plot,width=13,height=6.8,dpi=150)
}
plot_single(measured,'measured-displacement shared decoder','held_measured_absolute_decoder')
plot_single(replay,'baseline -> displacement -> shared decoder','held_baseline_absolute_bridge')
for (endpoint in c('total','departure')) {
  data <- copy(combo[evaluation=='held_cell'])
  data[,observed:=if(endpoint=='total') observed_log2 else observed_departure_log2]
  data[,predicted:=if(endpoint=='total') predicted_log2 else predicted_departure_log2]
  stats <- data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed,predicted)),'\nn = ',.N)),by=pair]
  plot <- ggplot(data,aes(observed,predicted))+geom_point(color='#0072B2',size=1.5)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
    geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
    facet_wrap(~pair,ncol=3,scales='free')+theme_bw(base_size=9)+theme(aspect.ratio=1)+
    labs(title=paste('v104 held-cell ALMANAC',endpoint),x=paste('Measured',endpoint,'log2 abundance'),
         y=paste('Predicted',endpoint,'log2 abundance'))
  ggsave(file.path(root,'figures',paste0('held_combination_',endpoint,'.pdf')),plot,width=10,height=6.4,device=cairo_pdf)
}

# Show calibration and direction separately on the principal lapatinib panel.
lap_data <- measured[evaluation=='held_cell' & drug=='Lapatinib ditosylate']
bridge_lap <- replay[evaluation=='held_cell' & drug=='Lapatinib ditosylate']
lap_plot_data <- rbind(lap_data[,.(observed_log2,predicted=signature_score,view='Held raw signature')],
                       lap_data[,.(observed_log2,predicted=-predicted_log2,view='Held calibrated sensitivity')],
                       bridge_lap[,.(observed_log2,predicted=-predicted_log2,view='Held full baseline bridge')])
lap_plot_data[,view:=factor(view,levels=c('Held raw signature','Held calibrated sensitivity','Held full baseline bridge'))]
lap_plot_data[,observed_S:=-observed_log2]
lap_stats <- lap_plot_data[,.(label=paste0('r = ',sprintf('%.3f',pearson(observed_S,predicted)),'\nn = ',.N)),by=view]
plot <- ggplot(lap_plot_data,aes(observed_S,predicted))+geom_point(color='#0072B2',alpha=.7,size=1.3)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
  geom_text(data=lap_stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
  facet_wrap(~view,nrow=1,scales='free')+theme_bw(base_size=10)+theme(aspect.ratio=1)+
  labs(title='v104 | Lapatinib 5 µM: fixed signature -> shared positive calibration -> baseline bridge',
       x='Measured sensitivity (-log2 abundance)',y='Score / predicted sensitivity')
ggsave(file.path(root,'figures/lapatinib_signature_calibration_bridge.pdf'),plot,width=10,height=4.1,device=cairo_pdf)
ggsave(file.path(root,'figures/lapatinib_signature_calibration_bridge.png'),plot,width=10,height=4.1,dpi=150)
write_table(data.table(assay=c('nci_FF','nci_FG'),shared_endpoint_scale=exp(model$par$assay)),'assay_scale')
writeLines(c('v104: fixed v95 correlation signatures and fold-training route assignments; shared positive affine calibration.',
  sprintf('%d drug/dose contexts; %d measured RNA/PRISM profiles, %d PRISM outcomes, %d ALMANAC pairs.',nrow(selection),nrow(profiles),nrow(prism),nrow(pairs)),
  sprintf('%d whole-cell exclusions across functional sources. AU565 excluded from RNA/functional fits.',length(loo_cells)),
  'Signature score: ((D - training RNA mean) / training RNA SD) times fixed W (v95 W x 0.25), selecting its training-derived route.',
  'Absolute sensitivity: Q beta + exp(Q eta) score. Q = intercept plus three RNA-only PCs of mean drug footprints.',
  'Beta and eta each have four shared coefficients; no drug identity, per-drug functional coefficient, or direct baseline input.',
  'No free gate, per-program outcome-fitted calibration gains/shifts, tanh, or quadratic readout.',
  'All treatment-level training Pearson correlations must agree with their uncalibrated score to 1e-10.',
  'Calibration values vary across whole-cell fits; identical within-fit ranking does not guarantee identical pooled LOO Pearson.',
  'Absolute endpoints fitted jointly: PRISM + .2 NCI single + .2 pair total + .2 pair departure, with original source weighting/geometry.',
  'L2 penalties: beta/eta .005, log NCI assay scales .01, interaction .05; fixed before evaluation.',
  'Two shared NCI assay scales and supported shared route-product interactions are also fitted. Unsupported route2-route2 coefficient fixed 0 and marked unlearned.',
  'Pair features use separate zero-referenced signature projections, multiplied by positive amplitude. No sum of displacement vectors.',
  'Raw zero displacement is not forced to 100% single-drug abundance in the mean-plus-centered-score model; zero-referenced pair interaction remains exactly 0 for an absent constituent.',
  'Lapatinib RNA predictor: ERBB2/HER2 + MET, standardized in RNA training cells; multivariate ridge lambda 1 predicts 40 standardized displacements, returned to measured units. Receptors never enter the functional decoder.',
  'Other baseline RNA maps unchanged and RNA-only. NCI single/pair training uses predicted RNA inputs.',
  paste('Included drugs:',paste(selection$drug,collapse='; ')),
  paste('Excluded functional contexts:',paste(excluded_drugs,collapse='; ')),
  'Fixed historical v95 signature bank and RNA footprint coordinates retained. Exclusion is an ablation of functional fitting, not a new signature-discovery analysis.',
  'HGF route is not identified by training data; both route assumptions are reported without selecting on Wilson labels. No HGF functional labels fitted.',
  'HGF receiver scales RNA-inferred/clamped as before. Query lapatinib 5 uM; Wilson reference experiment is 1 uM.',
  'Positive affine functional calibration is not demonstrated plate batch correction.',
  'Previously revealed labels/selected doses make all comparisons retrospective.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
input_paths <- c(archive_path,file.path(source_root,c('selected_doses.tsv','individual_training_outcomes.tsv',
                                                   'combination_training_outcomes.tsv')),
                file.path(parent,'models/final.qs'),'/work/src/program_pilot/104_erbb2_met_lapatinib_shared_decoder.R')
write_table(data.table(path=input_paths,sha256=vapply(input_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(metrics[evaluation=='held_cell'])
print(bridge_metrics[evaluation=='held_cell'])
print(combo_metrics[evaluation=='held_cell'])
print(query[hgf_context=='Equal-context mean'])
status(paste('Complete. Solver failures:',sum(rbindlist(audits)$convergence!=0L)))
