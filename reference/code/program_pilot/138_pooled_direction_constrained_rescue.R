#!/usr/bin/env Rscript
# v138: only RAF/MEK ligand outcomes; emphasize contrasts across ligands within each cell/drug.
.libPaths(c("/work/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest);library(glmnet)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
parent <- "/work/results/program_pilot/107_within_pair_departure_interaction/linear"
root <- "/work/results/program_pilot/138_pooled_direction_constrained_rescue/linear"
time_matched_root <- '/work/results/program_pilot/131_precision_weighted_HGF_bridge/linear'
external_root <- '/work/results/program_pilot/111_external_bypass_combination_transfer/linear'
programs <- sprintf("P%02d",1:40)
for(folder in c("data","models","figures"))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,"data",paste0(name,".tsv")),sep="\t")
pearson <- function(x,y)if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10)cor(x,y) else NA_real_
archive_path <- "/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds"
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
basis_reference <- readRDS("/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge/prepared.rds")
stopifnot(max(abs(archive$prep$W-basis_reference$W[rownames(archive$prep$W),programs]))<1e-12)
rm(archive,basis_reference)
bank_selection <- fread("/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv")
drugs <- bank_selection$drug
contexts <- fread(file.path(parent,"data/molecular_contexts.tsv"))
extra_contexts <- contexts[11:nrow(contexts)]
context_key <- function(drug,dose)paste(drug,sprintf("%.10g",dose),sep="|")
context_index <- function(drug,dose)match(context_key(drug,dose),context_key(contexts$drug,contexts$dose))

predict_rna <- function(model, cells, drug, dose = NULL) {
    if (is.null(dose)) 
        dose <- bank_selection$dose[match(drug, drugs)]
    raw <- as.matrix(baseline[match(cells, depmap_id), ..programs])
    features <- sweep(sweep(sweep(raw, 2, model$bmean), 2, model$bsd, "/") %*% model$basis, 2, model$bscale, "/")
    g <- match(drug, drugs)
    D <- matrix(0, length(cells), 40)
    for (j in unique(g)) {
        rows <- which(g == j)
        Z <- sweep(features[rows, , drop = FALSE], 2, model$bm[j, ]) %*% model$B[, , j]
        D[rows, ] <- sweep(sweep(Z, 2, model$sigma[j, ], "*"), 2, model$mu[j, ], "+")
    }
    rows <- which(drug == "Lapatinib ditosylate")
    if (length(rows)) {
        lap <- model$lap_receptors
        stopifnot(!is.null(lap))
        X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
        Z <- sweep(sweep(X, 2, lap$xmean), 2, lap$xsd, "/") %*% lap$B
        D[rows, ] <- sweep(sweep(Z, 2, lap$sigma, "*"), 2, lap$mu, "+")
    }
    for (i in seq_len(nrow(extra_contexts))) {
        dr <- extra_contexts$drug[i]
        ds <- extra_contexts$dose[i]
        rows <- which(drug == dr & abs(dose - ds) < 1e-10)
        if (!length(rows)) 
            next
        map <- model$extra_rna[[context_key(dr, ds)]]
        stopifnot(!is.null(map))
        if (map$kind == "receptors") {
            X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
            X <- sweep(sweep(X, 2, map$xmean), 2, map$xsd, "/")
        }
        else X <- sweep(features[rows, , drop = FALSE], 2, map$xmean)
        D[rows, ] <- sweep(sweep(X %*% map$B, 2, map$sigma, "*"), 2, map$mu, "+")
    }
    stopifnot(all(is.finite(D)), !anyNA(context_index(drug, dose)))
    D
}

signature_node <- function(D, g, ctx) {
    X <- (D - ctx$mu[g, , drop = FALSE])/ctx$sigma[g, , drop = FALSE]
    both <- X %*% ctx$W
    route <- ctx$route[g]
    score <- both[cbind(seq_len(nrow(D)), route)]
    absolute <- (D/ctx$sigma[g, , drop = FALSE]) %*% ctx$W
    activity <- matrix(0, nrow(D), 2)
    activity[cbind(seq_len(nrow(D)), route)] <- absolute[cbind(seq_len(nrow(D)), route)]
    list(score = score, both = both, activity = activity, interaction_activity = absolute, route = route)
}

pair_features <- function(A, B, W, basis) {
    stopifnot(ncol(A)==42L,ncol(B)==42L,nrow(A)==nrow(B),all(dim(W)==c(40L,2L)))
    rawA <- A[,1:40,drop=FALSE]
    rawB <- B[,1:40,drop=FALSE]
    # a and b are signed program contributions to each fixed sensitivity signature.
    # Reverse a sensitivity-directed move: R=sum(a+ b- + b+ a-).
    # Reinforce sensitivity-directed moves: K=sum(a+ b+).
    # Collapse AFTER programwise products, preserving magnitudes and zero reference.
    # Both-negative moves are not reinforcement of sensitivity. Pool routes BEFORE
    # functional fitting: neither route has its own functional coefficient.
    reversal <- reinforcement <- matrix(0, nrow(A), 2L)
    for (k in 1:2) {
      a <- sweep(rawA, 2, W[, k], '*')
      b <- sweep(rawB, 2, W[, k], '*')
      reversal[, k] <- rowSums(pmax(a, 0) * pmax(-b, 0) + pmax(b, 0) * pmax(-a, 0))
      reinforcement[, k] <- rowSums(pmax(a, 0) * pmax(b, 0))
    }
    features <- cbind(reversal = rowMeans(reversal), reinforcement = rowMeans(reinforcement))
    stopifnot(all(is.finite(features)))
    features
}

predict_single <- function(model, meta, measured_D = NULL) {
    dose <- if ("dose" %in% names(meta)) 
        meta$dose
    else bank_selection$dose[match(meta$drug, drugs)]
    D <- if (is.null(measured_D)) 
        predict_rna(model$rna, meta$depmap_id, meta$drug, dose)
    else measured_D
    g <- context_index(meta$drug, dose)
    node <- signature_node(D, g, model$context)
    Q <- model$context$Q[g, , drop = FALSE]
    level <- as.numeric(Q %*% model$par$level)
    amplitude <- as.numeric(exp(Q %*% model$par$log_amplitude))
    sensitivity <- level + amplitude * node$score
    assay <- if ("assay" %in% names(meta)) 
        match(meta$assay, c("secondary", "nci_FF", "nci_FG"))
    else rep(1L, nrow(meta))
    scale <- c(1, exp(model$par$assay))[assay]
    stopifnot(all(amplitude > 0), !anyNA(scale))
    list(sensitivity = sensitivity * scale, score = node$score, route = node$route, level = level, amplitude = amplitude, 
        activity = node$activity * amplitude, interaction_activity = node$interaction_activity * amplitude,
        program_activity = cbind(D/model$context$sigma[g,,drop=FALSE],node$both),D = D)
}

predict_pair <- function(model, meta) {
    A <- predict_single(model, meta[, .(depmap_id, drug = drugA, dose = doseA)])
    B <- predict_single(model, meta[, .(depmap_id, drug = drugB, dose = doseB)])
    products <- sweep(pair_features(A$program_activity,B$program_activity,model$context$W,model$interaction_basis), 2, model$context$interaction_scale,
        "/")
    J <- as.numeric(products %*% model$par$Gamma)
    scale <- exp(model$par$assay)[match(meta$SCREENER, c("FF", "FG"))]
    list(sensitivity = (A$sensitivity + B$sensitivity + J) * scale, departure = J * scale)
}


external_features <- function(original,held_cell=character()) {
  ctx <- original$context
  Aall <- Ball <- matrix(NA_real_,nrow(external),42)
  maps <- list()
  molecular_predictions <- list()
  for(dr in unique(external$drug)) {
    ii <- which(external$drug==dr)
    if(dr=='Lapatinib ditosylate') {
      native <- predict_single(original,external[ii,.(depmap_id,drug,dose)])
      Ad <- native$program_activity[,1:40,drop=FALSE]
      centred_A <- native$program_activity[,41:42,drop=FALSE]
      predicted <- native$D
      maps[[dr]] <- list(frozen_parent=parent,conditioning='unchanged frozen ERBB2/MET molecular map; native0.5/5uM')
    } else {
    train <- external_molecular[drug==dr & !depmap_id %in% held_cell]
    stopifnot(nrow(train)>=6L,!anyDuplicated(train$depmap_id))
    ids <- train$depmap_id
    D <- as.matrix(train[,..programs])
    rna <- original$rna
    raw <- as.matrix(baseline[match(ids,depmap_id),..programs])
    X <- sweep(sweep(sweep(raw,2,rna$bmean),2,rna$bsd,'/') %*% rna$basis,2,rna$bscale,'/')
    braf <- mutation$BRAF_hotspot[match(ids,mutation$depmap_id)]
    braf_mean <- mean(braf)
    braf_sd <- sd(braf)
    stopifnot(braf_sd>0)
    X <- cbind(X,(braf-braf_mean)/braf_sd)
    xm <- colMeans(X)
    X <- sweep(X,2,xm)
    mu <- colMeans(D)
    sigma <- pmax(apply(D,2,sd),1e-8)
    Y <- sweep(sweep(D,2,mu),2,sigma,'/')
    fit <- glmnet(X,Y,family='mgaussian',alpha=0,lambda=1,standardize=FALSE,intercept=FALSE)
    B <- vapply(coef(fit),function(x)as.numeric(x)[-1],numeric(9))
    ii <- which(external$drug==dr)
    raw <- as.matrix(baseline[match(external$depmap_id[ii],depmap_id),..programs])
    Xq <- sweep(sweep(sweep(raw,2,rna$bmean),2,rna$bsd,'/') %*% rna$basis,2,rna$bscale,'/')
    braf_query <- mutation$BRAF_hotspot[match(external$depmap_id[ii],mutation$depmap_id)]
    Xq <- cbind(Xq,(braf_query-braf_mean)/braf_sd)
    predicted <- sweep(sweep(sweep(Xq,2,xm) %*% B,2,sigma,'*'),2,mu,'+')
    Qd <- cbind(1,sweep(matrix((mu-ctx$fm)/ctx$fs,nrow=1) %*% ctx$fv,2,ctx$fscale,'/'))
    Ad <- sweep(predicted,2,sigma,'/')
    centred_A <- sweep(Ad,2,mu/sigma,'-') %*% ctx$W
    maps[[dr]] <- list(mu=mu,sigma=sigma,xmean=xm,B=B,training_cells=ids,
      BRAF_mean=braf_mean,BRAF_sd=braf_sd,lambda=1,
      conditioning='eight frozen baseline RNA directions plus BRAF-hotspot flag; functional labels excluded')
    }
    mi <- match(paste(external$depmap_id[ii],external$gene[ii]),
                paste(external_ligands$depmap_id,external_ligands$gene))
    stopifnot(!anyNA(mi))
    Dh <- as.matrix(external_ligands[mi,..programs])
    Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
    scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
                       crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
    minimum <- matrix(apply(ctx$sigma[1:10,],2,min),nrow(Dh),40,byrow=TRUE)
    maximum <- matrix(apply(ctx$sigma[1:10,],2,max),nrow(Dh),40,byrow=TRUE)
    scale <- pmin(pmax(exp(Qh %*% scale_map),minimum),maximum)
    Ah <- Dh/scale
    Aall[ii,] <- cbind(Ad,centred_A)
    Ball[ii,] <- cbind(Ah,Ah %*% ctx$W)
    colnames(predicted) <- programs
    molecular_predictions[[dr]] <- cbind(external[ii,.(depmap_id,gene,drug,dose)],as.data.table(predicted))
  }
  stopifnot(all(is.finite(Aall)),all(is.finite(Ball)))
  # Alternative source views express uncertainty about the same physical
  # experiment. Four HGF RNA studies share ONE total experiment weight.
  view_rows <- integer()
  view_labels <- character()
  Aviews <- Bviews <- matrix(numeric(),0,42)
  standalone_views <- numeric()
  for(i in seq_len(nrow(external))) {
    if(external$gene[i]=='HGF') {
      candidates <- source_views[source_views$CellLine==external$CellLine[i]]
      Dh <- as.matrix(candidates[,..programs])
      labels <- candidates$source_profile
    } else {
      mi <- match(paste(external$depmap_id[i],external$gene[i]),
        paste(external_ligands$depmap_id,external_ligands$gene))
      Dh <- as.matrix(external_ligands[mi,..programs])
      labels <- 'Characterized_mean'
    }
    Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
    scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
      crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
    minimum <- matrix(apply(ctx$sigma[1:10,],2,min),nrow(Dh),40,byrow=TRUE)
    maximum <- matrix(apply(ctx$sigma[1:10,],2,max),nrow(Dh),40,byrow=TRUE)
    scale <- pmin(pmax(exp(Qh %*% scale_map),minimum),maximum)
    H <- Dh/scale
    standalone <- -((H %*% ctx$W)*as.numeric(exp(Qh %*% original$par$log_amplitude)))[,1]
    Aviews <- rbind(Aviews,Aall[rep(i,nrow(Dh)),,drop=FALSE])
    Bviews <- rbind(Bviews,cbind(H,H %*% ctx$W))
    standalone_views <- c(standalone_views,standalone)
    view_rows <- c(view_rows,rep(i,nrow(Dh)))
    view_labels <- c(view_labels,labels)
  }
  list(A=Aall,B=Ball,Aviews=Aviews,Bviews=Bviews,view_rows=view_rows,view_labels=view_labels,
    standalone_views=standalone_views,rna=maps,molecular_predictions=rbindlist(molecular_predictions))
}


# Fixed configuration: no outcome-dependent source/context filtering.
within_weight <- .8
penalty_grid <- c(.003,.03,.3)
previous <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
original <- qread(file.path(previous,'models/final.qs'))
external <- fread(file.path(previous,'data/external_combination_outcomes.tsv'))[study=='Straussman2012']
external_ligands <- fread(file.path(previous,'data/external_ligand_programs.tsv'))
external_molecular <- fread(file.path(previous,'data/external_exact_drug_molecular_observations.tsv'))
mutation <- fread('/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv',
  select=c('V1','BRAF (673)'))
setnames(mutation,c('depmap_id','BRAF_hotspot'))
profiles <- fread(file.path(previous,'data/time_matched_HGF_program_inputs.tsv'))
source_views <- profiles[!grepl('^Leave_out_|^Primary_',source_profile)]
stopifnot(nrow(external)==156,uniqueN(external$gene)==13,uniqueN(external$depmap_id)==6,
  all(external$drug %in% c('PD-184352','PLX-4720')),
  !anyDuplicated(external[,.(depmap_id,drug,dose,gene)]))
external[,observation_id:=seq_len(.N)]
write_table(external,'interaction_training_outcomes')
write_table(source_views,'HGF_molecular_source_views')
for(name in c('baseline_predictions_frozen.tsv','baseline_bridge_correlations.tsv',
 'measured_predictions_frozen.tsv','measured_decoder_correlations.tsv','molecular_contexts.tsv')) {
  stopifnot(file.copy(file.path(previous,'data',name),file.path(root,'data',name),overwrite=TRUE))
}
terms <- c('reversal','reinforcement')
basis <- list(labels=terms,method='Equal route means before fitting; one nonnegative rescue strength and one nonnegative reinforcement strength')

# Baseline affects molecular displacements only. RNA fitting excludes validation cells;
# the entire single-drug decoder and all molecular calibrations remain frozen.
molecular_cache <- new.env(parent=emptyenv())
features_for <- function(held_cells) {
  key <- paste0('cells:',paste(sort(held_cells),collapse=','))
  if(!exists(key,molecular_cache,inherits=FALSE)) {
    features <- external_features(original,held_cells)
    features$products <- pair_features(features$Aviews,features$Bviews,original$context$W,basis)
    assign(key,features,molecular_cache)
  }
  get(key,molecular_cache,inherits=FALSE)
}

# Each physical cell/drug/ligand observation has equal weight. Four HGF RNA views
# divide that observation's weight; they are not four independent experiments.
geometry <- function(X,y,row_ids) {
  meta <- external[row_ids]
  group <- paste(meta$depmap_id,meta$drug,meta$dose,sep='|')
  physical <- unique(data.table(row=row_ids,group=group))
  n_per_group <- table(physical$group)
  views <- table(row_ids)
  w <- 1/as.numeric(views[as.character(row_ids)])/as.numeric(n_per_group[group])/length(n_per_group)
  stopifnot(abs(sum(w)-1)<1e-10)
  Xm <- matrix(0,nrow(X),ncol(X)); ym <- numeric(length(y))
  for(g in unique(group)) {
    ii <- which(group==g)
    wg <- w[ii]/sum(w[ii])
    Xm[ii,] <- matrix(colSums(X[ii,,drop=FALSE]*wg),length(ii),ncol(X),byrow=TRUE)
    ym[ii] <- sum(wg*y[ii])
  }
  list(X=X,Xc=X-Xm,Xm=Xm,y=y,yc=y-ym,ym=ym,w=w,group=group)
}

fit_g <- function(features,train_rows,penalty) {
  v <- which(features$view_rows %in% train_rows)
  er <- features$view_rows[v]
  stopifnot(setequal(unique(er),train_rows))
  z <- geometry(features$products[v,,drop=FALSE],-external$observed_departure_log2[er],er)
  variance_c <- sum(z$w*z$yc^2)
  variance_m <- sum(z$w*z$ym^2)
  stopifnot(variance_c>1e-10,variance_m>1e-10)
  scale <- sqrt(within_weight*colSums(z$Xc^2*z$w)+(1-within_weight)*colSums(z$Xm^2*z$w))
  stopifnot(all(scale>1e-10))
  Xc <- sweep(z$Xc,2,scale,'/')
  Xm <- sweep(z$Xm,2,scale,'/')
  gram <- within_weight*crossprod(Xc,Xc*z$w)/variance_c+
    (1-within_weight)*crossprod(Xm,Xm*z$w)/variance_m+diag(penalty,length(terms))
  rhs <- within_weight*as.numeric(crossprod(Xc,z$w*z$yc))/variance_c+
    (1-within_weight)*as.numeric(crossprod(Xm,z$w*z$ym))/variance_m
  # In abundance units, g = beta_R * R - beta_K * K, with beta_R,beta_K >= 0.
  # This is an explicit mechanistic restriction, not evidence that rescue occurs.
  # The existing prediction interface stores coefficients in sensitivity units.
  sign_map <- c(-1,1)
  G <- gram * outer(sign_map,sign_map)
  b <- rhs * sign_map
  solution <- optim(c(0,0),function(theta)as.numeric(crossprod(theta,G %*% theta)-2*crossprod(theta,b)),
    gr=function(theta)2*as.numeric(G %*% theta-b),method='L-BFGS-B',lower=c(0,0),
    control=list(factr=1e3,pgtol=1e-10,maxit=1000))
  theta <- solution$par
  gradient <- as.numeric(G %*% theta-b)
  kkt <- max(abs(gradient[theta>1e-8]),pmax(-gradient[theta<=1e-8],0),0)
  stopifnot(solution$convergence==0,all(theta>=0),kkt<1e-6)
  gamma <- sign_map * theta
  model <- original
  model$par$Gamma <- gamma
  model$context$interaction_scale <- scale
  model$interaction_basis <- basis
  model$external_rna <- features$rna
  model$interaction_training <- list(source='Straussman2012 only',
    physical_observation_ids=train_rows,physical_observations=length(train_rows),
    within_group='different ligands within fixed cell/drug/dose',within_weight=within_weight,
    penalty=penalty,unit='reported log2 net-growth departure',
    single_effects='frozen v136; no new ligand-alone calibration',
    direction_constraint='g_abundance = beta_R R - beta_K K; both beta nonnegative',
    kkt_residual=kkt)
  # Main-effect predictions are byte-identical; only g and its molecular fit metadata change.
  stopifnot(identical(model$rna,original$rna),identical(model$template,original$template),
    identical(model$par$level,original$par$level),
    identical(model$par$log_amplitude,original$par$log_amplitude),
    identical(model$ligand_single_conversion,original$ligand_single_conversion))
  model
}
predict_g <- function(model,features) {
  value <- -as.numeric(sweep(features$products,2,model$context$interaction_scale,'/') %*% model$par$Gamma)
  as.numeric(tapply(value,factor(features$view_rows,levels=seq_len(nrow(external))),mean))
}

# Inner whole-cell validation chooses the ridge penalty using the SAME across-ligand objective.
# Outer excluded cells/ligands are absent from every inner fit and penalty score.
audit <- list()
fit_selected <- function(held_cells=character(),held_genes=character(),label='final') {
  available <- which(!external$depmap_id %in% held_cells & !external$gene %in% held_genes)
  cells <- sort(unique(external$depmap_id[available]))
  groups <- rep(1:3,length.out=length(cells))
  predictions <- matrix(NA_real_,nrow(external),length(penalty_grid))
  for(k in unique(groups)) {
    validation <- cells[groups==k]
    excluded <- unique(c(held_cells,validation))
    features <- features_for(excluded)
    training <- available[!external$depmap_id[available] %in% validation]
    rows <- available[external$depmap_id[available] %in% validation]
    for(j in seq_along(penalty_grid)) {
      fit <- fit_g(features,training,penalty_grid[j])
      predictions[rows,j] <- predict_g(fit,features)[rows]
    }
  }
  stopifnot(all(is.finite(predictions[available,,drop=FALSE])))
  z <- geometry(predictions[available,,drop=FALSE],external$observed_departure_log2[available],available)
  score <- within_weight*colSums((z$Xc-z$yc)^2*z$w)/sum(z$w*z$yc^2)+
    (1-within_weight)*colSums((z$Xm-z$ym)^2*z$w)/sum(z$w*z$ym^2)
  chosen <- penalty_grid[which.min(score)]
  audit[[label]] <<- data.table(fit=label,penalty=penalty_grid,internal_selection_loss=score,
    selected=penalty_grid==chosen,held_cells=paste(held_cells,collapse=','),held_genes=paste(held_genes,collapse=','))
  features <- features_for(held_cells)
  model <- fit_g(features,available,chosen)
  list(model=model,predicted=predict_g(model,features),features=features)
}

# Predictions are frozen before evaluation; all outcomes have been examined in earlier pilots.
held <- list()
for(cell in sort(unique(external$depmap_id))) {
  label <- paste0('held_cell_',cell)
  fit <- fit_selected(held_cells=cell,label=label)
  rows <- which(external$depmap_id==cell)
  stopifnot(!any(rows %in% fit$model$interaction_training$physical_observation_ids))
  held[[label]] <- external[rows,.(observation_id,CellLine,depmap_id,drug,dose,gene)][,
    `:=`(predicted_departure_log2=fit$predicted[rows],evaluation='held_cell')]
  qsave(fit$model,file.path(root,'models',paste0(label,'.qs')))
  cat('Cell complete:',cell,'\n')
}
for(gene in sort(unique(external$gene))) {
  label <- paste0('held_ligand_',gene)
  fit <- fit_selected(held_genes=gene,label=label)
  rows <- which(external$gene==gene)
  stopifnot(!any(rows %in% fit$model$interaction_training$physical_observation_ids))
  held[[label]] <- external[rows,.(observation_id,CellLine,depmap_id,drug,dose,gene)][,
    `:=`(predicted_departure_log2=fit$predicted[rows],evaluation='held_ligand')]
  qsave(fit$model,file.path(root,'models',paste0(label,'.qs')))
  cat('Ligand complete:',gene,'\n')
}
fit <- fit_selected()
model <- fit$model
qsave(model,file.path(root,'models/final.qs'))
qsave(fit$features,file.path(root,'models/final_molecular_features.qs'))
write_table(rbindlist(audit),'training_only_penalty_selection')
predictions <- rbindlist(held)
write_table(predictions,'held_predictions_frozen')
write_table(data.table(term=terms,coefficient_log2_abundance=-model$par$Gamma,
  training_feature_scale=model$context$interaction_scale),'shared_interaction_coefficients')

# AU565 lapatinib 5 uM + HGF query: same vector, scale and single effects as v136.
ctx <- model$context
au <- profiles[CellLine=='AU565']
Dh <- as.matrix(au[,..programs])
Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
hscale <- pmin(pmax(exp(Qh %*% scale_map),matrix(apply(ctx$sigma[1:10,],2,min),nrow(Dh),40,byrow=TRUE)),
  matrix(apply(ctx$sigma[1:10,],2,max),nrow(Dh),40,byrow=TRUE))
H <- Dh/hscale
lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=5))
L <- lap$program_activity[rep(1,nrow(Dh)),,drop=FALSE]
P <- pair_features(L,cbind(H,H %*% ctx$W),ctx$W,basis)
stopifnot(max(abs(P-pair_features(cbind(H,H %*% ctx$W),L,ctx$W,basis)))<1e-12,
  all(pair_features(L,L*0,ctx$W,basis)==0))
contributions <- -sweep(sweep(P,2,model$context$interaction_scale,'/'),2,model$par$Gamma,'*')
interaction <- rowSums(contributions)
standalone <- -((H %*% ctx$W)*as.numeric(exp(Qh %*% model$par$log_amplitude)))[,1]*model$ligand_single_conversion
query <- data.table(source_profile=au$source_profile,cell='AU565',dose=5,
  lapatinib_pct_control=100*2^(-lap$sensitivity),HGF_alone_pct_control=100*2^standalone,
  interaction_log2=interaction,combination_pct_control=100*2^(-lap$sensitivity+standalone+interaction),
  endpoint_qualification='Conditional projection: RAF/MEK net-growth departure transferred to PRISM final-abundance units; not assay validated.')
write_table(query,'AU565_HGF_predictions_frozen')
write_table(data.table(source_profile=au$source_profile,contributions),'HGF_term_contributions')
paths <- file.path(root,'data',c('held_predictions_frozen.tsv','AU565_HGF_predictions_frozen.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')
scored <- merge(predictions,external[,.(observation_id,observed_departure_log2)],by='observation_id')
write_table(scored,'held_predictions_scored')
metrics <- scored[,.(n_ligands=.N,r=pearson(observed_departure_log2,predicted_departure_log2)),
  by=.(evaluation,CellLine,drug,dose)]
write_table(metrics,'across_ligand_held_correlations')
within <- scored[evaluation=='held_cell',.(n_cells=.N,r=pearson(observed_departure_log2,predicted_departure_log2)),by=.(gene,drug,dose)]
write_table(within,'across_cell_held_correlations')
hgf <- scored[gene=='HGF'][,.(evaluation,CellLine,drug,dose,observed_departure_log2,predicted_departure_log2)]
write_table(hgf,'HGF_training_source_held_predictions')
old <- fread(file.path(previous,'data/external_held_cell_predictions_frozen.tsv'))[study=='Straussman2012']
comparison <- rbindlist(list(
  old[,.(version='v136',r=pearson(observed_departure_log2,predicted_departure_log2)),by=.(CellLine,drug,dose)],
  scored[evaluation=='held_cell',.(version='v138',r=pearson(observed_departure_log2,predicted_departure_log2)),by=.(CellLine,drug,dose)]))
write_table(comparison,'comparison_across_ligand_held_cell_correlations')
plot_data <- scored[evaluation=='held_cell']
plot_data[,panel:=paste(CellLine,drug,sep=' | ')]
ann <- plot_data[,.(label=sprintf('r = %.2f',pearson(observed_departure_log2,predicted_departure_log2))),by=panel]
p <- ggplot(plot_data,aes(observed_departure_log2,predicted_departure_log2))+
  geom_point(aes(colour=gene=='HGF'),size=1.5)+scale_colour_manual(values=c('FALSE'='#666666','TRUE'='#D55E00'),guide='none')+
  geom_smooth(method='lm',formula=y~x,se=FALSE,colour='#0072B2',linewidth=.4)+
  geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=2.6)+
  facet_wrap(~panel,ncol=4,scales='free')+theme_bw(base_size=8)+theme(aspect.ratio=1)+
  labs(x='Measured growth departure (log2)',y='Predicted growth departure (log2)',
    title='v138 | held-out cells, across 13 ligands within each drug',
    caption='Orange = HGF. Previously examined data: retrospective evaluation. Different axes are independently scaled.')
ggsave(file.path(root,'figures/held_cell_across_ligand_scatters.pdf'),p,width=9,height=7.3,device=cairo_pdf)
ggsave(file.path(root,'figures/held_cell_across_ligand_scatters.png'),p,width=9,height=7.3,dpi=150)
writeLines(c('v138: RAF/MEK ligand data only; no ALMANAC, MEMA, HER2 or Wilson interaction outcomes.',
 '156 observations; 13 ligands x 6 cells x 2 drugs. Keep neutral/sensitizing ligands. No context selection by outcome.',
 'Two pooled g features: equal means of the two reversal and reinforcement scores. All single-drug components frozen.',
 'g_abundance = beta_R R - beta_K K, both coefficients nonnegative. This imposes direction, not rescue magnitude.',
 '80% within-cell/drug across-ligand contrasts, 20% cell/drug mean departure; ridge grid .003/.03/.3.',
 'Nested whole-cell selection; outer whole-cell and whole-ligand validation. No direct baseline-to-function terms.',
 'Four HGF molecular source views share one physical-observation weight. Weighted centering respects this.',
 'The AU565 query requires an unvalidated net-growth-to-final-abundance transfer assumption; no time extrapolation.',
 'Source selection follows retrospective diagnostics; repeated outcomes are not independent validation.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
input_paths <- c('/work/src/program_pilot/138_pooled_direction_constrained_rescue.R',file.path(previous,'models/final.qs'),
  file.path(previous,'data/external_combination_outcomes.tsv'),file.path(previous,'data/external_ligand_programs.tsv'),
  file.path(previous,'data/external_exact_drug_molecular_observations.tsv'),file.path(previous,'data/time_matched_HGF_program_inputs.tsv'))
write_table(data.table(path=input_paths,sha256=vapply(input_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(query[source_profile=='Primary_study_mean'])
print(comparison[,.(median_r=median(r)),by=version])
print(hgf[,.(median_observed=median(observed_departure_log2),median_predicted=median(predicted_departure_log2)),by=.(evaluation,drug)])
cat('v138 complete.\n')
