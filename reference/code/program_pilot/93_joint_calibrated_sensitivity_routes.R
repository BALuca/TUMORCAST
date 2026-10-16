#!/usr/bin/env Rscript
# Joint baseline -> measured programs -> calibrated shared nonlinear routes.
# Combinations retain separate inputs and share a symmetric route interaction.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(qs); library(ggplot2); library(digest)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive=TRUE, showWarnings=FALSE)
programs <- sprintf('P%02d', 1:40)
K <- 2L
seed <- 93012L
molecular_weight <- 10
within_weight <- .8
maximum_iterations <- 400L
assays <- c('secondary', 'nci_FF', 'nci_FG')
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep='\t')
status <- function(x) {
  cat(format(Sys.time()), x, '\n')
  flush.console()
  writeLines(x, file.path(root, 'STATUS.txt'))
}
pearson <- function(x,y) if(length(x)>=4L && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
key <- function(cell,drug) paste(cell,drug,sep='|')
macro_weights <- function(group) {
  counts <- table(group)
  1/as.numeric(counts[group])/length(counts)
}
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
archive <- readRDS(archive_path)
base <- as.data.table(archive$prep$raw)
v87 <- '/work/results/program_pilot/87_selected_dose_joint_program_bridge'
selection <- rbind(fread(file.path(v87,'data/selected_doses.tsv'))[,.(drug,dose)], data.table(drug='crizotinib',dose=5))
drugs <- selection$drug
molecular <- merge(as.data.table(archive$prep$m),selection,by=c('drug','dose'))[depmap_id!='ACH-000248']
single <- fread(file.path(v87,'data/individual_training_outcomes.tsv'))
criz <- fread('/work/results/program_pilot/91_lapatinib_crizotinib_5uM_joint_bridge/data/individual_training_outcomes.tsv')[drug=='crizotinib']
criz[,control_key:=paste('secondary',depmap_id,drug,dose,sep='|')]
keep <- c('depmap_id','drug','dose','y','assay','source','control_key')
single <- rbind(single[,..keep],criz[,..keep])
single[,group:=paste(assay,drug,dose,sep='|')]
single[,row_id:=.I]
pairs <- fread(file.path(v87,'data/combination_training_outcomes.tsv'))
folds <- fread('/work/results/program_pilot/91_lapatinib_crizotinib_5uM_joint_bridge/data/recipient_folds.tsv')
stopifnot(nrow(selection)==10L, nrow(pairs)==193L, !anyDuplicated(folds$depmap_id),
  !'Vemurafenib' %in% drugs, all(single$drug %in% drugs),
  all(abs(single$dose-selection$dose[match(single$drug,drugs)])<1e-10),
  all(is.finite(as.matrix(molecular[,..programs]))),all(is.finite(single$y)),
  all(single$depmap_id %in% base$depmap_id),all(pairs$depmap_id %in% base$depmap_id),
  !anyNA(folds$fold[match(unique(c(single$depmap_id,pairs$depmap_id)),folds$depmap_id)]))
write_table(selection,'selected_doses')
write_table(single,'individual_training_outcomes')
write_table(pairs,'combination_training_outcomes')
write_table(folds,'recipient_folds')
status(paste('Inputs:',nrow(single),'single outcomes;',nrow(molecular),'RNA profiles;',nrow(pairs),'combinations.'))
rm(archive)

# The calibration and mixing maps have four molecular-footprint coefficients,
# including an intercept, rather than arbitrary independent drug-only parameters.
sizes <- c(B=8*40*10,W=40*K,gain=4*40,shift=4*40,mix=4*K,C=3*K,Gamma=3)
ends <- cumsum(sizes)
starts <- ends-sizes+1
indices <- lapply(names(sizes),function(n) starts[n]:ends[n])
names(indices) <- names(sizes)
unpack <- function(theta) list(B=array(theta[indices$B],c(8,40,10)), W=matrix(theta[indices$W],40,K),
  gain=matrix(theta[indices$gain],4,40),shift=matrix(theta[indices$shift],4,40),
  mix=matrix(theta[indices$mix],4,K),C=matrix(theta[indices$C],3,K),Gamma=theta[indices$Gamma])
pack <- function(par) unlist(par[names(sizes)],use.names=FALSE)
getbaseline <- function(ids) as.matrix(base[match(ids,depmap_id),..programs])
geometry <- function(y,weight,group) {
  g <- match(group,unique(group))
  gw <- as.numeric(rowsum(weight,g,reorder=FALSE))
  ym <- as.numeric(rowsum(y*weight,g,reorder=FALSE))/gw
  list(y=y,weight=weight,group=g,gw=gw,within_scale=max(sum(weight*(y-ym[g])^2),1e-8),
    absolute_scale=max(sum(weight*y^2),1e-8))
}
loss <- function(pred,geo) {
  error <- pred-geo$y
  em <- as.numeric(rowsum(error*geo$weight,geo$group,reorder=FALSE))/geo$gw
  centered <- error-em[geo$group]
  list(value=within_weight*sum(geo$weight*centered^2)/geo$within_scale+
    (1-within_weight)*sum(geo$weight*error^2)/geo$absolute_scale,
    gradient=2*geo$weight*(within_weight*centered/geo$within_scale+(1-within_weight)*error/geo$absolute_scale))
}
prepare <- function(excluded) {
  m <- molecular[!depmap_id %in% excluded]
  ss <- single[!depmap_id %in% excluded]
  pp <- pairs[!depmap_id %in% excluded]
  ids <- sort(unique(m$depmap_id))
  rawbase <- getbaseline(ids)
  bmean <- colMeans(rawbase)
  bsd <- apply(rawbase,2,sd)
  bn <- sweep(sweep(rawbase,2,bmean),2,bsd,'/')
  bv <- svd(bn,nu=0,nv=8)$v
  bscale <- sqrt(colMeans((bn%*%bv)^2))
  query <- unique(rbind(ss[,.(depmap_id,drug)],pp[,.(depmap_id,drug=drugA)],
    pp[,.(depmap_id,drug=drugB)],m[,.(depmap_id,drug)]))
  query[,g:=match(drug,drugs)]
  qkey <- key(query$depmap_id,query$drug)
  s <- sweep(sweep(sweep(getbaseline(query$depmap_id),2,bmean),2,bsd,'/')%*%bv,2,bscale,'/')
  mi <- match(key(m$depmap_id,m$drug),qkey)
  mg <- match(m$drug,drugs)
  mu <- sigma <- matrix(0,10,40)
  bm <- matrix(0,10,8)
  B <- array(0,c(8,40,10))
  X <- matrix(0,nrow(query),40)
  observed <- matrix(0,nrow(m),40)
  group_rows <- group_mrows <- list()
  for(g in 1:10) {
    ii <- which(query$g==g)
    mm <- which(mg==g)
    stopifnot(length(mm)>=10)
    D <- as.matrix(m[mm,..programs])
    mu[g,] <- colMeans(D)
    sigma[g,] <- apply(D,2,sd)
    stopifnot(all(sigma[g,]>1e-10))
    bm[g,] <- colMeans(s[mi[mm],,drop=FALSE])
    s[ii,] <- sweep(s[ii,,drop=FALSE],2,bm[g,])
    observed[mm,] <- sweep(sweep(D,2,mu[g,]),2,sigma[g,],'/')
    fit <- glmnet(s[mi[mm],,drop=FALSE],observed[mm,,drop=FALSE],family='mgaussian',alpha=0,
      lambda=1,intercept=FALSE,standardize=FALSE,thresh=1e-9)
    B[,,g] <- vapply(coef(fit),function(z) as.numeric(z)[-1],numeric(8))
    X[ii,] <- s[ii,,drop=FALSE]%*%B[,,g]
    group_rows[[g]] <- ii
    group_mrows[[g]] <- mm
  }
  # Three RNA-only footprint descriptor directions condition calibration/routes.
  # The functional layer still receives all 40 programs, each on its drug's scale.
  fm <- colMeans(mu)
  fs <- pmax(apply(mu,2,sd),1e-8)
  F <- sweep(sweep(mu,2,fm),2,fs,'/')
  fv <- svd(F,nu=0,nv=3)$v
  fscale <- sqrt(colMeans((F%*%fv)^2))
  Q <- cbind(1,sweep(F%*%fv,2,fscale,'/'))
  scale_map <- solve(crossprod(Q)/10+diag(c(0,1,1,1)),crossprod(Q,log(sigma))/10)
  is <- match(key(ss$depmap_id,ss$drug),qkey)
  ia <- match(key(pp$depmap_id,pp$drugA),qkey)
  ib <- match(key(pp$depmap_id,pp$drugB),qkey)
  obs <- merge(ss[assay=='secondary'],m[,.(depmap_id,drug,molecular_row=.I)],by=c('depmap_id','drug'))
  source_weight <- c(prism_secondary=.8,nci_single=.2)
  sw <- numeric(nrow(ss))
  for(src in unique(ss$source)) {
    ii <- which(ss$source==src)
    units <- 1/as.numeric(table(ss$control_key[ii])[ss$control_key[ii]])
    totals <- tapply(units,ss$group[ii],sum)
    sw[ii] <- source_weight[src]*units/as.numeric(totals[ss$group[ii]])/uniqueN(ss$group[ii])
  }
  sw <- sw/sum(sw)
  pg <- paste(pp$SCREENER,pp$pair)
  list(query=query,s=s,B=B,X=X,mu=mu,sigma=sigma,Q=Q,observed=observed,m=m,mg=mg,mi=mi,
    group_rows=group_rows,group_mrows=group_mrows,ss=ss,pp=pp,is=is,ia=ia,ib=ib,obs=obs,
    assay=match(ss$assay,assays),pair_assay=match(paste0('nci_',pp$SCREENER),assays),
    geoS=geometry(ss$y,sw,ss$group),geoP=geometry(pp$y,macro_weights(pg),pg),
    geoI=geometry(pp$I,macro_weights(pg),pg),geoO=geometry(obs$y,macro_weights(obs$drug),obs$drug),
    context=list(bmean=bmean,bsd=bsd,bv=bv,bscale=bscale,bm=bm,mu=mu,sigma=sigma,Q=Q,
      fm=fm,fs=fs,fv=fv,fscale=fscale,scale_map=scale_map,training_RNA_cells=ids,excluded=excluded))
}
route_node <- function(X,g,par,Q,mu,sigma) {
  gain <- exp(Q%*%par$gain)
  shift <- Q%*%par$shift
  logits <- Q%*%par$mix
  alpha <- exp(logits-apply(logits,1,max))
  alpha <- alpha/rowSums(alpha)
  zero <- -mu/sigma
  calibrated <- X*gain[g,,drop=FALSE]+shift[g,,drop=FALSE]
  calibrated0 <- (zero*gain+shift)[g,,drop=FALSE]
  h <- tanh(calibrated%*%par$W)
  h0 <- tanh(calibrated0%*%par$W)
  activity <- (h-h0)*alpha[g,,drop=FALSE]
  list(activity=activity,X=X,g=g,gain=gain,zero=zero,calibrated=calibrated,
    calibrated0=calibrated0,h=h,h0=h0,alpha=alpha)
}
route_back <- function(node,gradient,par,Q) {
  g <- node$g
  dh <- gradient*node$alpha[g,,drop=FALSE]*(1-node$h^2)
  dh0 <- -gradient*node$alpha[g,,drop=FALSE]*(1-node$h0^2)
  dx <- dh%*%t(par$W)
  dx0 <- dh0%*%t(par$W)
  grouped <- function(z) {
    out <- matrix(0,nrow(Q),ncol(z))
    temp <- rowsum(z,g,reorder=FALSE)
    out[as.integer(rownames(temp)),] <- temp
    out
  }
  gg <- grouped((dx*node$X+dx0*node$zero[g,,drop=FALSE])*node$gain[g,,drop=FALSE])
  gs <- grouped(dx+dx0)
  ga <- grouped(gradient*(node$h-node$h0))
  gl <- node$alpha*(ga-rowSums(ga*node$alpha))
  list(X=dx*node$gain[g,,drop=FALSE],W=crossprod(node$calibrated,dh)+crossprod(node$calibrated0,dh0),
    gain=crossprod(Q,gg),shift=crossprod(Q,gs),mix=crossprod(Q,gl))
}
gamma_matrix <- function(v) matrix(c(v[1],v[2],v[2],v[3]),2,2)
pair_features <- function(A,B) cbind(A[,1]*B[,1],A[,1]*B[,2]+A[,2]*B[,1],A[,2]*B[,2])
# Training-only per-drug ridge signatures initialize the two shared directions;
# they do not remain independent drug decoders after joint optimization.
initialize <- function(p,job) {
  teacher <- matrix(0,10,40)
  for(g in 1:10) {
    rows <- which(p$obs$drug==drugs[g])
    ii <- p$obs$molecular_row[rows]
    y <- p$obs$y[rows]-mean(p$obs$y[rows])
    fit <- glmnet(p$observed[ii,,drop=FALSE],y,alpha=0,lambda=.1,
      intercept=FALSE,standardize=FALSE,thresh=1e-9)
    teacher[g,] <- as.numeric(coef(fit))[-1]
  }
  set.seed(seed+job)
  norm <- sqrt(rowSums(teacher^2))
  km <- kmeans(teacher/pmax(norm,1e-8),centers=2,nstart=20)
  target <- cbind(as.numeric(km$cluster==1),as.numeric(km$cluster==2))-.5
  mix <- solve(crossprod(p$Q)+diag(1,4),crossprod(p$Q,target))*2
  par <- list(B=p$B,W=t(km$centers)*.25,gain=matrix(0,4,40),shift=matrix(0,4,40),
    mix=mix,C=matrix(-1,3,K),Gamma=rep(0,3))
  node <- route_node(p$X,p$query$g,par,p$Q,p$mu,p$sigma)
  design <- matrix(0,nrow(p$ss),3*K)
  for(a in 1:3) {
    ii <- which(p$assay==a)
    design[ii,(a-1)*K+1:K] <- node$activity[p$is[ii],,drop=FALSE]
  }
  on <- route_node(p$observed[p$obs$molecular_row,,drop=FALSE],match(p$obs$drug,drugs),par,p$Q,p$mu,p$sigma)
  oc <- matrix(0,nrow(p$obs),3*K)
  oc[,1:K] <- on$activity
  views <- function(X,geo) {
    xm <- rowsum(X*geo$weight,geo$group,reorder=FALSE)/geo$gw
    ym <- as.numeric(rowsum(geo$y*geo$weight,geo$group,reorder=FALSE))/geo$gw
    list(X=rbind(X-xm[geo$group,,drop=FALSE],X),y=c(geo$y-ym[geo$group],geo$y),
      w=c(geo$weight*within_weight/geo$within_scale,geo$weight*(1-within_weight)/geo$absolute_scale))
  }
  s <- views(design,p$geoS)
  o <- views(oc,p$geoO)
  fit <- glmnet(rbind(s$X,o$X),c(s$y,o$y),weights=c(s$w,o$w),alpha=0,lambda=.01,
    intercept=FALSE,standardize=TRUE,thresh=1e-9)
  par$C <- t(matrix(as.numeric(coef(fit))[-1],K,3))
  par
}
fit_joint <- function(p,job) {
  initial <- pack(initialize(p,job))
  cached_theta <- cached <- NULL
  objective <- function(theta) {
    if(identical(theta,cached_theta)) return(cached)
    par <- unpack(theta)
    X <- matrix(0,nrow(p$query),40)
    for(g in 1:10) {
      ii <- p$group_rows[[g]]
      X[ii,] <- p$s[ii,,drop=FALSE]%*%par$B[,,g]
    }
    node <- route_node(X,p$query$g,par,p$Q,p$mu,p$sigma)
    act <- node$activity
    predS <- rowSums(act[p$is,,drop=FALSE]*par$C[p$assay,,drop=FALSE])
    A <- act[p$ia,,drop=FALSE]
    B <- act[p$ib,,drop=FALSE]
    feat <- pair_features(A,B)
    predI <- as.numeric(feat%*%par$Gamma)
    predP <- rowSums((A+B)*par$C[p$pair_assay,,drop=FALSE])+predI
    obsnode <- route_node(p$observed[p$obs$molecular_row,,drop=FALSE],match(p$obs$drug,drugs),par,p$Q,p$mu,p$sigma)
    predO <- as.numeric(obsnode$activity%*%par$C[1,])
    ls <- loss(predS,p$geoS)
    lp <- loss(predP,p$geoP)
    li <- loss(predI,p$geoI)
    lo <- loss(predO,p$geoO)
    gs <- ls$gradient
    gp <- .2*lp$gradient
    gi <- .2*li$gradient
    go <- lo$gradient
    gj <- gp+gi
    ga <- gp*par$C[p$pair_assay,,drop=FALSE]+gj*(B%*%gamma_matrix(par$Gamma))
    gb <- gp*par$C[p$pair_assay,,drop=FALSE]+gj*(A%*%gamma_matrix(par$Gamma))
    native_grad <- matrix(0,nrow(p$query),K)
    add <- function(index,values) {
      summed <- rowsum(values,index,reorder=FALSE)
      ii <- as.integer(rownames(summed))
      native_grad[ii,] <<- native_grad[ii,,drop=FALSE]+summed
    }
    add(p$is,gs*par$C[p$assay,,drop=FALSE])
    add(p$ia,ga)
    add(p$ib,gb)
    ng <- route_back(node,native_grad,par,p$Q)
    og <- route_back(obsnode,go*matrix(par$C[1,],nrow(p$obs),K,byrow=TRUE),par,p$Q)
    cgrad <- matrix(0,3,K)
    for(a in 1:3) {
      ii <- which(p$assay==a)
      pi <- which(p$pair_assay==a)
      cgrad[a,] <- colSums(act[p$is[ii],,drop=FALSE]*gs[ii])+
        colSums((A[pi,,drop=FALSE]+B[pi,,drop=FALSE])*gp[pi])
    }
    cgrad[1,] <- cgrad[1,]+colSums(obsnode$activity*go)
    penalty <- .01*sum(par$W^2)+sum(par$gain^2)+sum(par$shift^2)+.05*sum(par$mix^2)+
      .005*sum(par$C^2)+.05*sum(par$Gamma^2)
    grad <- list(B=par$B*(2*molecular_weight/(10*40)),
      W=ng$W+og$W+.02*par$W,gain=ng$gain+og$gain+2*par$gain,
      shift=ng$shift+og$shift+2*par$shift,mix=ng$mix+og$mix+.1*par$mix,
      C=cgrad+.01*par$C,Gamma=as.numeric(crossprod(feat,gj))+.1*par$Gamma)
    value <- ls$value+.2*lp$value+.2*li$value+lo$value+penalty
    dx <- ng$X
    for(g in 1:10) {
      ii <- p$group_rows[[g]]
      mm <- p$group_mrows[[g]]
      qi <- p$mi[mm]
      error <- X[qi,,drop=FALSE]-p$observed[mm,,drop=FALSE]
      value <- value+molecular_weight/10*(mean(error^2)+sum(par$B[,,g]^2)/40)
      dx[qi,] <- dx[qi,,drop=FALSE]+2*molecular_weight*error/(10*length(error))
      grad$B[,,g] <- grad$B[,,g]+crossprod(p$s[ii,,drop=FALSE],dx[ii,,drop=FALSE])
    }
    result <- list(value=value,gradient=pack(grad))
    stopifnot(is.finite(value),all(is.finite(result$gradient)))
    cached_theta <<- theta
    cached <<- result
    result
  }
  if(job==1) {
    positions <- unique(unlist(lapply(indices,function(i)c(head(i,1),tail(i,1)))))
    probe <- initial
    probe[indices$Gamma] <- c(.01,-.01,.02)
    analytic <- objective(probe)$gradient[positions]
    numerical <- vapply(positions,function(j) {
      plus <- minus <- probe
      step <- 1e-5*max(1,abs(probe[j]))
      plus[j] <- plus[j]+step
      minus[j] <- minus[j]-step
      (objective(plus)$value-objective(minus)$value)/(2*step)
    },numeric(1))
    relative <- abs(analytic-numerical)/pmax(1,abs(analytic),abs(numerical))
    write_table(data.table(position=positions,analytic,numerical,relative),'analytic_gradient_check')
    stopifnot(max(relative)<1e-5)
  }
  lower <- rep(-Inf,length(initial))
  upper <- rep(Inf,length(initial))
  for(n in c('gain','shift')) {
    lower[indices[[n]]] <- -.5
    upper[indices[[n]]] <- .5
  }
  for(n in c('W','mix')) {
    lower[indices[[n]]] <- -3
    upper[indices[[n]]] <- 3
  }
  lower[indices$C] <- -20
  upper[indices$C] <- 20
  initial <- pmin(pmax(initial,lower),upper)
  status(paste('Fit',job,'joint optimization'))
  opt <- optim(initial,function(t)objective(t)$value,function(t)objective(t)$gradient,
    method='L-BFGS-B',lower=lower,upper=upper,
    control=list(maxit=maximum_iterations,factr=1e8,pgtol=1e-4,lmm=10))
  if(opt$convergence==1) opt <- optim(opt$par,function(t)objective(t)$value,function(t)objective(t)$gradient,
    method='L-BFGS-B',lower=lower,upper=upper,
    control=list(maxit=maximum_iterations,factr=1e8,pgtol=1e-4,lmm=10))
  status(paste('Fit',job,'finished; solver',opt$convergence))
  list(par=unpack(opt$par),context=p$context,training_single_cells=unique(p$ss$depmap_id),
    training_pair_cells=unique(p$pp$depmap_id),selected=selection,
    optimization=list(convergence=opt$convergence,evaluations=opt$counts,max_gradient=max(abs(objective(opt$par)$gradient))))
}
predict_RNA <- function(model,cell,drug) {
  context <- model$context
  s <- sweep(sweep(sweep(getbaseline(cell),2,context$bmean),2,context$bsd,'/')%*%context$bv,2,context$bscale,'/')
  g <- match(drug,drugs)
  X <- matrix(0,length(cell),40)
  for(j in unique(g)) {
    ii <- which(g==j)
    X[ii,] <- sweep(s[ii,,drop=FALSE],2,context$bm[j,])%*%model$par$B[,,j]
  }
  list(X=X,D=X*context$sigma[g,,drop=FALSE]+context$mu[g,,drop=FALSE],g=g)
}
predict_single <- function(model,meta) {
  q <- predict_RNA(model,meta$depmap_id,meta$drug)
  context <- model$context
  node <- route_node(q$X,q$g,model$par,context$Q,context$mu,context$sigma)
  list(prediction=rowSums(node$activity*model$par$C[match(meta$assay,assays),,drop=FALSE]),
    activity=node$activity,D=q$D)
}
predict_pair <- function(model,meta) {
  aa <- predict_single(model,meta[,.(depmap_id,drug=drugA,assay=paste0('nci_',SCREENER))])
  bb <- predict_single(model,meta[,.(depmap_id,drug=drugB,assay=paste0('nci_',SCREENER))])
  J <- as.numeric(pair_features(aa$activity,bb$activity)%*%model$par$Gamma)
  list(total=aa$prediction+bb$prediction+J,departure=J)
}
sp <- op <- cp <- audits <- list()
for(job in 1:4) {
  label <- if(job==4) 'final' else paste0('fold_',job)
  held <- if(job==4) character() else folds[fold==job,depmap_id]
  excluded <- union('ACH-000248',held)
  status(paste(label,'preparing'))
  p <- prepare(excluded)
  model <- fit_joint(p,job)
  model$held_cells <- held
  stopifnot(!any(held %in% model$training_single_cells),!any(held %in% model$context$training_RNA_cells),
    !any(held %in% model$training_pair_cells))
  qsave(model,file.path(root,'models',paste0(label,'.qs')))
  test <- if(job==4) single[assay=='secondary'] else single[depmap_id %in% held & assay=='secondary']
  out <- test[,.(depmap_id,drug,dose,observed_log2=y)]
  out[,`:=`(fold=job,evaluation=if(job==4) 'fitted' else 'held_cell',predicted_log2=predict_single(model,test)$prediction)]
  sp[[label]] <- out
  measured <- if(job==4) molecular else molecular[depmap_id %in% held]
  obs <- merge(measured,test[,.(depmap_id,drug,observed_log2=y)],by=c('depmap_id','drug'))
  g <- match(obs$drug,drugs)
  context <- model$context
  X <- (as.matrix(obs[,..programs])-context$mu[g,,drop=FALSE])/context$sigma[g,,drop=FALSE]
  node <- route_node(X,g,model$par,context$Q,context$mu,context$sigma)
  truth <- obs[,.(depmap_id,drug,dose,observed_log2)]
  truth[,`:=`(fold=job,evaluation=if(job==4) 'fitted' else 'held_cell',
    predicted_log2=as.numeric(node$activity%*%model$par$C[1,]))]
  op[[label]] <- truth
  pt <- if(job==4) pairs else pairs[depmap_id %in% held]
  result <- predict_pair(model,pt)
  pout <- pt[,.(row_id,depmap_id,pair,SCREENER,drugA,drugB,doseA,doseB,observed_log2=y,observed_departure_log2=I)]
  pout[,`:=`(fold=job,evaluation=if(job==4) 'fitted' else 'held_cell',
    predicted_log2=result$total,predicted_departure_log2=result$departure)]
  cp[[label]] <- pout
  audits[[label]] <- data.table(job=label,convergence=model$optimization$convergence,
    function_evaluations=model$optimization$evaluations[1],max_gradient=model$optimization$max_gradient)
  write_table(rbindlist(sp),'single_predictions_frozen')
  write_table(rbindlist(op),'measured_displacement_predictions_frozen')
  write_table(rbindlist(cp),'combination_predictions_frozen')
  write_table(rbindlist(audits),'fit_audit')
}
qsave(list(single=rbindlist(sp),measured=rbindlist(op),combination=rbindlist(cp),audits=rbindlist(audits)),
  file.path(root,'data/completed_fits.qs'))
status('All four fits complete; predictions frozen')
