#!/usr/bin/env Rscript
# Frozen correlation signatures, learned input calibration and molecular mixtures,
# shared absolute functional readout, and a symmetric ALMANAC route interaction.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(ggplot2); library(digest)})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/96_fixed_signatures_calibrated_shared_absolute_decoder'
parent <- '/work/results/program_pilot/95_leave_one_cell_out_correlation_signatures'
source_root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
programs <- sprintf('P%02d', 1:40)
within_weight <- .8
max_iterations <- 300L
signature_amplitude <- .25
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive=TRUE, showWarnings=FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep='\t')
status <- function(message) {
  cat(format(Sys.time()), message, '\n')
  flush.console()
  writeLines(message, file.path(root, 'STATUS.txt'))
}
key <- function(cell, drug) paste(cell, drug, sep='|')
pearson <- function(x, y) if (length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
selection <- fread(file.path(source_root, 'data/selected_doses.tsv'))
single <- fread(file.path(source_root, 'data/individual_training_outcomes.tsv'))
pairs <- fread(file.path(source_root, 'data/combination_training_outcomes.tsv'))
parent_folds <- fread(file.path(parent, 'data/recipient_folds.tsv'))
drugs <- selection$drug
short <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
           'Lapatinib', 'Topotecan', 'Vandetanib', 'Vismodegib', 'Crizotinib')
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
molecular <- merge(as.data.table(archive$prep$m), selection, by=c('drug','dose'))[depmap_id!='ACH-000248']
prism <- single[assay=='secondary']
nci <- single[assay!='secondary' & depmap_id!='ACH-000248']
profiles <- merge(molecular, prism[,.(depmap_id, drug, dose, sensitivity=-y)], by=c('depmap_id','drug','dose'))
rm(archive)
stopifnot(nrow(profiles)==286L, nrow(pairs)==193L, all(is.finite(single$y)), !anyDuplicated(baseline$depmap_id))
loo_cells <- sort(unique(c(profiles$depmap_id, nci$depmap_id, pairs$depmap_id)))
folds <- data.table(depmap_id=loo_cells, fold=seq_along(loo_cells))
write_table(selection, 'selected_doses')
write_table(folds, 'recipient_folds')
status(paste(length(loo_cells), 'whole-cell LOO fits; fixed v95 correlation signatures; absolute outcomes.'))

# Cached RNA-only v95 models can be replayed; additional RNA recipients require
# the identical RNA-only fit with their own cell removed.
fit_rna <- function(train_rna) {
  ids <- sort(unique(train_rna$depmap_id))
  raw <- as.matrix(baseline[match(ids, depmap_id), ..programs])
  bmean <- colMeans(raw)
  bsd <- apply(raw, 2, sd)
  normalized <- sweep(sweep(raw, 2, bmean), 2, bsd, '/')
  basis <- svd(normalized, nu=0, nv=8)$v
  bscale <- sqrt(colMeans((normalized %*% basis)^2))
  B <- array(0, c(8,40,10))
  mu <- sigma <- matrix(0,10,40)
  bm <- matrix(0,10,8)
  for (g in 1:10) {
    rows <- which(train_rna$drug==drugs[g])
    D <- as.matrix(train_rna[rows, ..programs])
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
predict_rna <- function(model, cells, drug) {
  raw <- as.matrix(baseline[match(cells,depmap_id),..programs])
  features <- sweep(sweep(sweep(raw,2,model$bmean),2,model$bsd,'/') %*% model$basis,2,model$bscale,'/')
  g <- match(drug,drugs)
  D <- matrix(0,length(cells),40)
  for (j in unique(g)) {
    rows <- which(g==j)
    Z <- sweep(features[rows,,drop=FALSE],2,model$bm[j,]) %*% model$B[,,j]
    D[rows,] <- sweep(sweep(Z,2,model$sigma[j,],'*'),2,model$mu[j,],'+')
  }
  stopifnot(all(is.finite(D)))
  D
}

# Signed correlation templates remain frozen. Footprint coordinates are RNA-only.
make_context <- function(template) {
  fm <- colMeans(template$mu)
  fs <- pmax(apply(template$mu,2,sd),1e-8)
  F <- sweep(sweep(template$mu,2,fm),2,fs,'/')
  fv <- svd(F,nu=0,nv=3)$v
  fscale <- sqrt(colMeans((F %*% fv)^2))
  Q <- cbind(1,sweep(F %*% fv,2,fscale,'/'))
  W <- template$W*signature_amplitude
  # Only shifts in the two signature directions affect the decoder. Expressing
  # them before W avoids fitting 160 unidentifiable program-shift coefficients.
  shift_basis <- W %*% solve(crossprod(W))
  list(mu=template$mu,sigma=template$sigma,Q=Q,W=W,shift_basis=shift_basis,
       fm=fm,fs=fs,fv=fv,fscale=fscale)
}
sizes <- c(gain=160L,shift=8L,mix=4L,C=2L,assay=2L,Gamma=3L)
ends <- cumsum(sizes)
starts <- ends-sizes+1L
indices <- lapply(names(sizes),function(n)starts[n]:ends[n])
names(indices) <- names(sizes)
unpack <- function(theta) list(gain=matrix(theta[indices$gain],4,40),shift=matrix(theta[indices$shift],4,2),
  mix=matrix(theta[indices$mix],4,1),C=theta[indices$C],assay=theta[indices$assay],Gamma=theta[indices$Gamma])
pack <- function(par) unlist(par[names(sizes)],use.names=FALSE)

# Calibration really precedes the signatures: learned positive program gains
# and a molecular-footprint shift vector. The reference at D=0 is subtracted.
route_node <- function(D,g,par,ctx) {
  X <- (D-ctx$mu[g,,drop=FALSE])/ctx$sigma[g,,drop=FALSE]
  zero <- -ctx$mu/ctx$sigma
  gain <- exp(ctx$Q %*% par$gain)
  program_shift <- (ctx$Q %*% par$shift) %*% t(ctx$shift_basis)
  calibrated <- X*gain[g,,drop=FALSE]+program_shift[g,,drop=FALSE]
  calibrated0 <- zero*gain+program_shift
  h <- tanh(calibrated %*% ctx$W)
  h0 <- tanh(calibrated0 %*% ctx$W)
  a <- as.numeric(plogis(ctx$Q %*% par$mix))
  alpha <- cbind(a,1-a)
  activity <- (h-h0[g,,drop=FALSE])*alpha[g,,drop=FALSE]
  list(activity=activity,h=h,h0=h0,g=g,X=X,zero=zero,gain=gain,alpha=alpha)
}
route_back <- function(node,gradient,ctx) {
  g <- node$g
  dh <- gradient*node$alpha[g,,drop=FALSE]*(1-node$h^2)
  dh0 <- -gradient*node$alpha[g,,drop=FALSE]*(1-node$h0[g,,drop=FALSE]^2)
  grouped <- function(z) {
    out <- matrix(0,10,ncol(z))
    temp <- rowsum(z,g,reorder=FALSE)
    out[as.integer(rownames(temp)),] <- temp
    out
  }
  dg <- grouped((dh %*% t(ctx$W))*node$X+(dh0 %*% t(ctx$W))*node$zero[g,,drop=FALSE])*node$gain
  da <- grouped(gradient*(node$h-node$h0[g,,drop=FALSE]))
  dm <- node$alpha[,1]*node$alpha[,2]*(da[,1]-da[,2])
  list(gain=crossprod(ctx$Q,dg),shift=crossprod(ctx$Q,grouped(dh+dh0)),mix=crossprod(ctx$Q,dm))
}
pair_features <- function(A,B) cbind(A[,1]*B[,1],A[,1]*B[,2]+A[,2]*B[,1],A[,2]*B[,2])
gamma_matrix <- function(v) matrix(c(v[1],v[2],v[2],v[3]),2,2)

# Absolute targets remain intact. Within-treatment variation contributes 80%
# of each data objective so drug averages cannot dominate model selection.
geometry <- function(y,group,control=NULL) {
  g <- match(group,unique(group))
  w <- rep(1,length(y))
  if (!is.null(control)) w <- 1/as.numeric(table(control)[control])
  totals <- as.numeric(rowsum(w,g,reorder=FALSE))
  w <- w/totals[g]/length(totals)
  gw <- as.numeric(rowsum(w,g,reorder=FALSE))
  ym <- as.numeric(rowsum(y*w,g,reorder=FALSE))/gw
  list(y=y,g=g,w=w,gw=gw,within=max(sum(w*(y-ym[g])^2),1e-8),absolute=max(sum(w*y^2),1e-8))
}
loss <- function(pred,geo) {
  error <- pred-geo$y
  em <- as.numeric(rowsum(error*geo$w,geo$g,reorder=FALSE))/geo$gw
  centered <- error-em[geo$g]
  list(value=within_weight*sum(geo$w*centered^2)/geo$within+
    (1-within_weight)*sum(geo$w*error^2)/geo$absolute,
    gradient=2*geo$w*(within_weight*centered/geo$within+(1-within_weight)*error/geo$absolute))
}
prepare <- function(held, template, rna) {
  train <- profiles[!depmap_id %in% held]
  ns <- nci[!depmap_id %in% held]
  pp <- pairs[!depmap_id %in% held]
  query <- unique(rbind(ns[,.(depmap_id,drug)],pp[,.(depmap_id,drug=drugA)],pp[,.(depmap_id,drug=drugB)]))
  qkey <- key(query$depmap_id,query$drug)
  ctx <- make_context(template)
  list(ctx=ctx,train=train,ns=ns,pp=pp,query=query,
    Dobs=as.matrix(train[,..programs]),gobs=match(train$drug,drugs),
    Dquery=predict_rna(rna,query$depmap_id,query$drug),gquery=match(query$drug,drugs),
    is=match(key(ns$depmap_id,ns$drug),qkey),ia=match(key(pp$depmap_id,pp$drugA),qkey),
    ib=match(key(pp$depmap_id,pp$drugB),qkey),assay=match(ns$assay,c('nci_FF','nci_FG')),
    pair_assay=match(pp$SCREENER,c('FF','FG')),
    geoO=geometry(train$sensitivity,train$drug),
    geoN=geometry(-ns$y,ns$group,ns$control_key),
    geoP=geometry(-pp$y,pp$pair),geoI=geometry(-pp$I,pp$pair),
    training_cells=unique(c(train$depmap_id,ns$depmap_id,pp$depmap_id)))
}

fit_decoder <- function(p,template,job) {
  target <- ifelse(template$membership==1,1,-1)
  mix <- solve(crossprod(p$ctx$Q)+diag(1,4),crossprod(p$ctx$Q,target))
  initial_par <- list(gain=matrix(0,4,40),shift=matrix(0,4,2),mix=mix,
    C=c(1,1),assay=c(-1,-1),Gamma=c(0,0,0))
  initial <- pack(initial_par)
  cache_theta <- cache_value <- NULL
  objective <- function(theta) {
    if (identical(theta,cache_theta)) return(cache_value)
    par <- unpack(theta)
    obs <- route_node(p$Dobs,p$gobs,par,p$ctx)
    native <- route_node(p$Dquery,p$gquery,par,p$ctx)
    R <- native$activity
    scales <- exp(par$assay)
    N <- R[p$is,,drop=FALSE]
    A <- R[p$ia,,drop=FALSE]
    B <- R[p$ib,,drop=FALSE]
    F <- pair_features(A,B)
    po <- as.numeric(obs$activity %*% par$C)
    pn <- as.numeric(N %*% par$C)*scales[p$assay]
    rawI <- as.numeric(F %*% par$Gamma)
    pp <- (as.numeric((A+B) %*% par$C)+rawI)*scales[p$pair_assay]
    pi <- rawI*scales[p$pair_assay]
    lo <- loss(po,p$geoO)
    ln <- loss(pn,p$geoN)
    lp <- loss(pp,p$geoP)
    li <- loss(pi,p$geoI)
    go <- lo$gradient
    gn <- .2*ln$gradient
    gp <- .2*lp$gradient
    gi <- .2*li$gradient
    gj <- (gp+gi)*scales[p$pair_assay]
    gradN <- (gn*scales[p$assay]) %o% par$C
    gradA <- (gp*scales[p$pair_assay]) %o% par$C+gj*(B %*% gamma_matrix(par$Gamma))
    gradB <- (gp*scales[p$pair_assay]) %o% par$C+gj*(A %*% gamma_matrix(par$Gamma))
    gradient_native <- matrix(0,nrow(R),2)
    for (block in list(list(index=p$is,value=gradN),list(index=p$ia,value=gradA),list(index=p$ib,value=gradB))) {
      summed <- rowsum(block$value,block$index,reorder=FALSE)
      ii <- as.integer(rownames(summed))
      gradient_native[ii,] <- gradient_native[ii,,drop=FALSE]+summed
    }
    on <- route_back(obs,go %o% par$C,p$ctx)
    nn <- route_back(native,gradient_native,p$ctx)
    Cgrad <- as.numeric(crossprod(obs$activity,go)+crossprod(N,gn*scales[p$assay])+
                          crossprod(A+B,gp*scales[p$pair_assay]))
    assaygrad <- vapply(1:2,function(a)sum(gn[p$assay==a]*pn[p$assay==a])+
      sum(gp[p$pair_assay==a]*pp[p$pair_assay==a])+sum(gi[p$pair_assay==a]*pi[p$pair_assay==a]),numeric(1))
    penalty <- sum(par$gain^2)+sum(par$shift^2)+.05*sum(par$mix^2)+
      .005*sum(par$C^2)+.01*sum(par$assay^2)+.05*sum(par$Gamma^2)
    gradient <- list(gain=on$gain+nn$gain+2*par$gain,shift=on$shift+nn$shift+2*par$shift,
      mix=on$mix+nn$mix+.1*par$mix,C=Cgrad+.01*par$C,assay=assaygrad+.02*par$assay,
      Gamma=as.numeric(crossprod(F,gj))+.1*par$Gamma)
    value <- list(value=lo$value+.2*ln$value+.2*lp$value+.2*li$value+penalty,gradient=pack(gradient))
    stopifnot(is.finite(value$value),all(is.finite(value$gradient)))
    cache_theta <<- theta
    cache_value <<- value
    value
  }
  if (job==1L) {
    positions <- unique(unlist(lapply(indices,function(i)c(head(i,1),tail(i,1)))))
    probe <- initial
    probe[indices$Gamma] <- c(.01,-.01,.02)
    analytic <- objective(probe)$gradient[positions]
    numerical <- vapply(positions,function(i) {
      plus <- minus <- probe
      plus[i] <- plus[i]+1e-5
      minus[i] <- minus[i]-1e-5
      (objective(plus)$value-objective(minus)$value)/2e-5
    },numeric(1))
    difference <- abs(analytic-numerical)/pmax(1,abs(analytic),abs(numerical))
    write_table(data.table(position=positions,analytic,numerical,relative_difference=difference),'analytic_gradient_check')
    stopifnot(max(difference)<1e-5)
  }
  lower <- upper <- rep(0,length(initial))
  for (name in names(sizes)) {
    bound <- switch(name,gain=.5,shift=1,mix=3,C=20,assay=4,Gamma=10)
    lower[indices[[name]]] <- -bound
    upper[indices[[name]]] <- bound
  }
  opt <- optim(pmin(pmax(initial,lower),upper),function(t)objective(t)$value,function(t)objective(t)$gradient,
    method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=max_iterations,factr=1e8,pgtol=1e-4))
  if (opt$convergence==1L) opt <- optim(opt$par,function(t)objective(t)$value,function(t)objective(t)$gradient,
    method='L-BFGS-B',lower=lower,upper=upper,control=list(maxit=max_iterations,factr=1e8,pgtol=1e-4))
  list(par=unpack(opt$par),context=p$ctx,training_cells=p$training_cells,
       convergence=opt$convergence,evaluations=opt$counts,max_gradient=max(abs(objective(opt$par)$gradient)))
}

predict_single <- function(model, meta, measured_D=NULL) {
  D <- if (is.null(measured_D)) predict_rna(model$rna,meta$depmap_id,meta$drug) else measured_D
  node <- route_node(D,match(meta$drug,drugs),model$par,model$context)
  assay <- if ('assay' %in% names(meta)) match(meta$assay,c('secondary','nci_FF','nci_FG')) else rep(1L,nrow(meta))
  scales <- c(1,exp(model$par$assay))[assay]
  list(sensitivity=as.numeric(node$activity %*% model$par$C)*scales,activity=node$activity,D=D)
}
predict_pair <- function(model,meta) {
  A <- predict_single(model,meta[,.(depmap_id,drug=drugA)])
  B <- predict_single(model,meta[,.(depmap_id,drug=drugB)])
  J <- as.numeric(pair_features(A$activity,B$activity) %*% model$par$Gamma)
  scales <- exp(model$par$assay)[match(meta$SCREENER,c('FF','FG'))]
  list(sensitivity=(A$sensitivity+B$sensitivity+J)*scales,departure=J*scales)
}

render_heatmap <- function(template,label) {
  ro <- hclust(dist(template$correlations),method='ward.D2')$order
  co <- template$drug_tree$order
  annotation <- data.frame(Signature=paste0('Signature ',template$membership[co]),
    Dose=paste0(selection$dose[co],' µM'),row.names=drugs[co])
  left <- draw_heatmap(template$correlations[ro,co,drop=FALSE],scale_rows=FALSE,
    value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),hm_name='Training Pearson r',
    col_anno_df=annotation,col_anno_cols=c('Signature','Dose'),
    palette_overrides=list(Signature=c('Signature 1'='#D58B36','Signature 2'='#775DA6'),
      Dose=c('0.05 µM'='#DCEAF2','0.5 µM'='#6BAED6','5 µM'='#08519C')),
    show_row_names=TRUE,row_names_side='left',row_names_gp=grid::gpar(fontsize=8),
    show_column_names=TRUE,column_labels=short[co],column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
    cluster_rows=FALSE,cluster_columns=FALSE,width=grid::unit(105,'mm'),
    title='Observed training correlations',column_title_gp=grid::gpar(fontsize=11))
  right <- draw_heatmap(template$collapsed[ro,,drop=FALSE],scale_rows=FALSE,
    value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),hm_name='Collapsed mean r',
    show_row_names=FALSE,show_column_names=TRUE,column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
    cluster_rows=FALSE,cluster_columns=FALSE,width=grid::unit(24,'mm'),
    title='Two frozen signatures',column_title_gp=grid::gpar(fontsize=11))
  drawn <- ComplexHeatmap::draw(left+right,heatmap_legend_side='bottom',annotation_legend_side='right',
    column_title=paste0('v96 | ',label),column_title_gp=grid::gpar(fontsize=12),
    padding=grid::unit(c(15,5,5,5),'mm'),ht_gap=grid::unit(8,'mm'))
  stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn,name='Training Pearson r')),1:40))
  grid::grid.text('Signatures are fixed training-correlation directions. Learned gain/shift acts BEFORE them.\nMolecular-footprint mixtures and one shared readout predict absolute -log2 viability; zero displacement gives 100% viability.',
    x=.5,y=grid::unit(4,'mm'),just=c('centre','bottom'),gp=grid::gpar(fontsize=8))
}

measured_predictions <- baseline_predictions <- combo_predictions <- audits <- calibrations <- mixtures <- list()
final_parent <- qread(file.path(parent,'models/final.qs'))
book <- file.path(root,'figures/training_correlations_and_collapsed_signatures_LOO.pdf')
cairo_pdf(book,width=9,height=9,onefile=TRUE)
book_device <- dev.cur()
final_job <- length(loo_cells)+1L
for (job in seq_len(final_job)) {
  label <- if (job==final_job) 'final' else paste0('loo_',sprintf('%02d',job))
  held <- if (job==final_job) character() else loo_cells[job]
  parent_index <- match(held,parent_folds$depmap_id)
  if (length(parent_index) && !is.na(parent_index)) {
    ancestor <- qread(file.path(parent,'models',paste0('loo_',sprintf('%02d',parent_index),'.qs')))
    template <- ancestor$template
    rna <- ancestor$rna
  } else {
    template <- final_parent$template
    rna <- if (length(held) && held %in% molecular$depmap_id) fit_rna(molecular[depmap_id!=held]) else final_parent$rna
  }
  stopifnot(!any(held %in% template$training_cells),!any(held %in% rna$training_cells))
  status(paste(label,if (length(held)) held else '',': fitting learned calibration, molecular mixture, shared absolute readout and ALMANAC interaction.'))
  prepared <- prepare(held,template,rna)
  fit <- fit_decoder(prepared,template,job)
  model <- c(fit,list(template=template,rna=rna,held_cells=held,selected=selection))
  stopifnot(!any(held %in% model$training_cells))
  qsave(model,file.path(root,'models',paste0(label,'.qs')))
  evaluation <- if (job==final_job) 'fitted' else 'held_cell'
  test <- if (job==final_job) profiles else profiles[depmap_id %in% held]
  if (nrow(test)) {
    result <- predict_single(model,test,as.matrix(test[,..programs]))
    p <- test[,.(depmap_id,drug,dose,observed_log2=-sensitivity)]
    p[,`:=`(predicted_log2=-result$sensitivity,fold=job,evaluation=evaluation)]
    measured_predictions[[label]] <- p
  }
  test <- if (job==final_job) prism else prism[depmap_id %in% held]
  if (nrow(test)) {
    result <- predict_single(model,test)
    p <- test[,.(depmap_id,drug,dose,observed_log2=y)]
    p[,`:=`(predicted_log2=-result$sensitivity,fold=job,evaluation=evaluation)]
    baseline_predictions[[label]] <- p
  }
  test_pairs <- if (job==final_job) pairs else pairs[depmap_id %in% held]
  if (nrow(test_pairs)) {
    result <- predict_pair(model,test_pairs)
    p <- test_pairs[,.(depmap_id,pair,SCREENER,drugA,drugB,doseA,doseB,observed_log2=y,observed_departure_log2=I)]
    p[,`:=`(predicted_log2=-result$sensitivity,predicted_departure_log2=-result$departure,fold=job,evaluation=evaluation)]
    combo_predictions[[label]] <- p
  }
  zeros <- route_node(matrix(0,10,40),1:10,model$par,model$context)$activity
  zero_error <- max(abs(zeros))
  stopifnot(zero_error<1e-12,identical(model$context$W,template$W*signature_amplitude))
  audits[[label]] <- data.table(fit=label,held_cell=if (length(held)) held else NA_character_,
    convergence=fit$convergence,evaluations=fit$evaluations[1],max_gradient=fit$max_gradient,
    functional_leakage=length(intersect(held,model$training_cells)),
    RNA_leakage=length(intersect(held,rna$training_cells)),zero_displacement_effect=zero_error)
  gain <- exp(model$context$Q %*% model$par$gain)
  shift <- (model$context$Q %*% model$par$shift) %*% t(model$context$shift_basis)
  calibrations[[label]] <- data.table(fit=label,drug=rep(drugs,each=40),program=rep(programs,10),
    learned_gain=as.numeric(t(gain)),learned_shift=as.numeric(t(shift)),
    effective_input_gain=as.numeric(t(gain/model$context$sigma)))
  alpha1 <- as.numeric(plogis(model$context$Q %*% model$par$mix))
  mixtures[[label]] <- data.table(fit=label,drug=drugs,signature_1=alpha1,signature_2=1-alpha1)
  # Molecular templates differ in the 30 matched-cell LOO fits. Other folds use
  # the same frozen map, which is also shown with their own exclusion label.
  dev.set(book_device)
  render_heatmap(template,paste(label,if (length(held)) paste('held',held) else 'all-data fit'))
  if (job %in% c(1L,30L,final_job)) {
    png(file.path(root,'figures',paste0(label,'_correlations_and_signatures.png')),width=9,height=9,units='in',res=150,type='cairo')
    render_heatmap(template,label)
    dev.off()
  }
  write_table(rbindlist(audits),'fit_audit')
}
dev.set(book_device)
dev.off()
measured <- rbindlist(measured_predictions)
replay <- rbindlist(baseline_predictions)
combo <- rbindlist(combo_predictions)
# Truly untrained PRISM cells can use the final model without duplicate fits.
external <- copy(replay[evaluation=='fitted' & !depmap_id %in% loo_cells])
stopifnot(!any(external$depmap_id %in% model$training_cells),!any(external$depmap_id %in% model$rna$training_cells))
external[,evaluation:='held_cell']
replay <- rbind(replay,external)
stopifnot(nrow(replay[evaluation=='held_cell'])==nrow(prism),
  nrow(measured[evaluation=='held_cell'])==nrow(profiles),nrow(combo[evaluation=='held_cell'])==nrow(pairs),
  !anyDuplicated(replay[evaluation=='held_cell',.(depmap_id,drug)]))
write_table(measured,'measured_predictions_frozen')
write_table(replay,'baseline_predictions_frozen')
write_table(combo,'combination_predictions_frozen')
write_table(rbindlist(calibrations),'learned_program_calibration')
write_table(rbindlist(mixtures),'molecular_signature_mixtures')

# Transported HGF is an extrapolation, never a fitting label. Its receiver scale
# is RNA-inferred exactly as in v93; the shared footprint maps handle calibration.
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
hmean <- as.numeric(hgf[context=='Equal-context mean',..programs])
ctx <- model$context
Qh <- cbind(1,sweep(matrix((hmean-ctx$fm)/ctx$fs,nrow=1) %*% ctx$fv,2,ctx$fscale,'/'))
scale_map <- solve(crossprod(ctx$Q)/10+diag(c(0,1,1,1)),crossprod(ctx$Q,log(ctx$sigma))/10)
hscale_raw <- as.numeric(exp(Qh %*% scale_map))
hscale <- pmin(pmax(hscale_raw,apply(ctx$sigma,2,min)),apply(ctx$sigma,2,max))
ctx$Q <- rbind(ctx$Q,Qh)
ctx$mu <- rbind(ctx$mu,hmean)
ctx$sigma <- rbind(ctx$sigma,hscale)
hh <- route_node(as.matrix(hgf[,..programs]),rep(11L,nrow(hgf)),model$par,ctx)
lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate'))
Rlap <- lap$activity[rep(1L,nrow(hgf)),,drop=FALSE]
Sh <- as.numeric(hh$activity %*% model$par$C)
J <- as.numeric(pair_features(Rlap,hh$activity) %*% model$par$Gamma)
query <- data.table(cell='AU565',dose=5,hgf_context=hgf$context,
  lapatinib_pct_vehicle=100*2^(-lap$sensitivity),hgf_alone_pct_vehicle=100*2^(-Sh),
  combination_pct_vehicle=100*2^(-lap$sensitivity-Sh-J),
  lapatinib_in_hgf_pct_hgf_control=100*2^(-lap$sensitivity-J),interaction_log2=-J,
  qualification='HGF receiver scale RNA-inferred; HGF/Wilson labels never fitted; no lapatinib+HGF combination training')
write_table(query,'au565_hgf_5uM_predictions_frozen')
write_table(data.table(program=programs,inferred_SD=hscale_raw,used_SD=hscale),'hgf_scale_transfer')
paths <- file.path(root,'data',c('measured_predictions_frozen.tsv','baseline_predictions_frozen.tsv',
  'combination_predictions_frozen.tsv','au565_hgf_5uM_predictions_frozen.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')

# All predictions are frozen before evaluation. Endpoints are absolute log2
# abundance; r compares cells at the same treatment, never different drugs.
metrics <- measured[,.(cells=.N,r=pearson(observed_log2,predicted_log2)),by=.(drug,dose,evaluation)]
baseline_metrics <- replay[,.(cells=.N,r=pearson(observed_log2,predicted_log2),
  observed_sensitive=sum(observed_log2<= -1),predicted_sensitive=sum(predicted_log2<= -1),
  observed_median_pct=median(100*2^observed_log2),predicted_median_pct=median(100*2^predicted_log2)),by=.(drug,dose,evaluation)]
combo_metrics <- combo[,.(cells=.N,r_total=pearson(observed_log2,predicted_log2),
  r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_departure_SD=sd(observed_departure_log2),predicted_departure_SD=sd(predicted_departure_log2)),by=.(pair,SCREENER,evaluation)]
write_table(metrics,'measured_decoder_correlations')
write_table(baseline_metrics,'baseline_bridge_correlations')
write_table(combo_metrics,'combination_correlations')

plot_single <- function(table,heading,stem) {
  data <- copy(table[evaluation=='held_cell'])
  data[,panel:=factor(paste0(short[match(drug,drugs)],' | ',dose,' µM'),levels=paste0(short,' | ',selection$dose,' µM'))]
  stats <- data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed_log2,predicted_log2)),'\nn = ',.N)),by=panel]
  plot <- ggplot(data,aes(observed_log2,predicted_log2))+geom_point(color='#0072B2',alpha=.7,size=1.2)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
    geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.15,size=2.8)+
    facet_wrap(~panel,ncol=5,scales='free')+theme_bw(base_size=9)+
    theme(aspect.ratio=1,panel.grid.minor=element_blank())+
    labs(title=paste('v96 |',heading),x='Measured absolute log2 abundance',y='Predicted absolute log2 abundance',
      caption='Whole-cell exclusions; Pearson across cells at fixed drug/dose. 0 = 100% viability; -1 = 50%; -2 = 25%.\nNo free drug-specific functional intercept. Correlation-derived signature vectors remain fixed.')
  ggsave(file.path(root,'figures',paste0(stem,'.pdf')),plot,width=13,height=6.8,device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(stem,'.png')),plot,width=13,height=6.8,dpi=150)
}
plot_single(measured,'held-cell measured-displacement decoder','held_measured_absolute_decoder')
plot_single(replay,'held-cell baseline → displacement → shared absolute decoder','held_baseline_absolute_bridge')
for (endpoint in c('total','departure')) {
  data <- copy(combo[evaluation=='held_cell'])
  data[,observed:=if (endpoint=='total') observed_log2 else observed_departure_log2]
  data[,predicted:=if (endpoint=='total') predicted_log2 else predicted_departure_log2]
  stats <- data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed,predicted)),'\nn = ',.N)),by=pair]
  plot <- ggplot(data,aes(observed,predicted))+geom_point(color='#0072B2',size=1.5)+
    geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
    geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
    facet_wrap(~pair,ncol=3,scales='free')+theme_bw(base_size=9)+theme(aspect.ratio=1)+
    labs(title=paste('v96 held-cell ALMANAC',endpoint),x=paste('Measured',endpoint,'log2 abundance'),
      y=paste('Predicted',endpoint,'log2 abundance'))
  ggsave(file.path(root,'figures',paste0('held_combination_',endpoint,'.pdf')),plot,width=10,height=6.4,device=cairo_pdf)
}
final_gain <- exp(model$context$Q %*% model$par$gain)
final_mix <- as.numeric(plogis(model$context$Q %*% model$par$mix))
write_table(data.table(signature=1:2,shared_sensitivity_readout=model$par$C),'shared_readout')
write_table(data.table(assay=c('nci_FF','nci_FG'),shared_endpoint_scale=exp(model$par$assay)),'assay_scale')
writeLines(c('v96: fixed training-only correlation signatures -> learned program calibration -> molecular-footprint mixture -> shared absolute readout.',
  'Ten selected drug/dose contexts unchanged. 62 whole-cell LOO exclusions cover matched PRISM and ALMANAC recipients.',
  'The 30 matched-cell correlation templates and RNA-only baseline models are frozen v95 artifacts. Other RNA-only exclusions are refitted identically.',
  'Calibration gains: 4 RNA-footprint coefficients x 40 programs, exp(QG), shrunk toward 1; coefficient bounds +/-0.5.',
  'RNA means/SDs are calibration anchors, not the complete calibration. Learned effective input slope = exp(QG)/RNA_SD.',
  'Shifts live in the two-signature span (4 x 2 coefficients). They are applied as a program vector before W; reference D=0 is subtracted after tanh.',
  'Mixture = sigmoid(Qm), 4 shared molecular-footprint coefficients; no independent drug readout coefficients.',
  'W is fixed v95 L2-normalized correlation template x 0.25. Shared sensitivity readout has 2 coefficients and no intercept.',
  'A fixed tanh activation lets calibration shifts affect sensitivity without a free functional intercept. Zero D gives exactly 100% viability.',
  'Two global positive ALMANAC assay scales align endpoint amplitude with PRISM; they have no drug specificity or intercept.',
  'The symmetric combination interaction has 3 shared coefficients and retains separate route vectors for the two drugs.',
  'Joint fitting: measured PRISM functional loss + .2 NCI single + .2 ALMANAC total + .2 ALMANAC departure; each uses 80% within-treatment and 20% absolute outcomes.',
  'Penalties: gain 1, shift 1, mix .05, readout .005, log assay scale .01, Gamma .05; fixed before evaluation.',
  'The baseline -> displacement map remains RNA-only and frozen during functional fitting. NCI and combination functional training use its predicted displacements.',
  'HGF query at lapatinib 5 uM: receiver SD RNA-inferred/clamped, no HGF or Wilson labels fitted. Wilson reference dose is 1 uM.',
  'Previously exposed labels and selected doses make these retrospective comparisons.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
inputs <- c(archive_path,file.path(source_root,'data',c('selected_doses.tsv','individual_training_outcomes.tsv','combination_training_outcomes.tsv')),
  file.path(parent,'models/final.qs'),'/work/src/program_pilot/96_fixed_signatures_calibrated_shared_absolute_decoder.R')
write_table(data.table(path=inputs,sha256=vapply(inputs,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(metrics[evaluation=='held_cell'])
print(baseline_metrics[evaluation=='held_cell'])
print(combo_metrics[evaluation=='held_cell'])
print(query[hgf_context=='Equal-context mean'])
status(paste('Complete. Solver failures:',sum(rbindlist(audits)$convergence!=0L)))
