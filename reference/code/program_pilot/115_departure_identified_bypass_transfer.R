#!/usr/bin/env Rscript
# Same nonlinear g; MEMA informs identifiable ligand contrasts, not unknown no-ligand levels.
.libPaths(c("/work/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest);library(glmnet)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
parent <- "/work/results/program_pilot/107_within_pair_departure_interaction/linear"
root <- "/work/results/program_pilot/115_departure_identified_bypass_transfer/linear"
external_root <- '/work/results/program_pilot/111_external_bypass_combination_transfer/linear'
mema_root <- "/work/grant_rebuild_2026-10-04/transfer_solution"
challenge_genes <- c("HGF","NRG1")
primary_endpoint <- "normalized_cell_count"
source_weight <- .2
interaction_penalty <- .05
assay_penalty <- .01
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
    stopifnot(ncol(A)==40L,ncol(B)==40L,nrow(A)==nrow(B),all(dim(W)==c(40L,2L)))
    SA <- A %*% W
    SB <- B %*% W
    # RBF movement relative to zero: each modifier is exactly zero at no input.
    # Centers and length scale depend only on characterized molecular profiles.
    transform <- function(X) {
      distance <- pmax(outer(rowSums(X^2),rowSums(basis$centres^2),'+')-2*tcrossprod(X,basis$centres),0)
      zero <- exp(-rowSums(basis$centres^2)/(2*basis$bandwidth^2))
      sweep(exp(-distance/(2*basis$bandwidth^2)),2,zero,'-')
    }
    KA <- transform(A)
    KB <- transform(B)
    features <- do.call(cbind,lapply(1:2,function(k)KB*SA[,k]+KA*SB[,k]))
    colnames(features) <- paste0(rep(c('route1_','route2_'),each=nrow(basis$centres)),rep(basis$labels,2))
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
        program_activity = D/model$context$sigma[g,,drop=FALSE]*amplitude,D = D)
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


# Existing aggregation: average isoform log ratios within ECM, then ECMs equally.
labels <- fread(file.path(mema_root,'all_gene_endpoints.tsv'))
ligands <- fread(file.path(mema_root,'recipient_ligand_program_scores.tsv'))
stopifnot(nrow(ligands)==56L,!anyDuplicated(ligands[,.(CellLine,gene)]),
          all(is.finite(as.matrix(ligands[,..programs]))),
          all(labels$ecm_n==48L),all(labels$endpoint %in% c('normalized_cell_count','raw_cell_count')))
# Source HGF labels are not joined into training; challenge scoring occurs after freeze.
training_labels <- labels[!gene %in% challenge_genes]
stopifnot(uniqueN(training_labels$gene)==26L,nrow(training_labels)==104L)
write_table(training_labels,'mema_nonchallenge_training_outcomes')
write_table(ligands,'mema_molecular_inputs')
external <- fread(file.path(external_root,'data/external_combination_outcomes.tsv'))
external_ligands <- fread(file.path(external_root,'data/external_ligand_programs.tsv'))
external_molecular <- fread(file.path(external_root,'data/external_exact_drug_molecular_observations.tsv'))
mutation_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv'
mutation <- fread(mutation_path,select=c('V1','BRAF (673)'))
setnames(mutation,c('depmap_id','BRAF_hotspot'))
stopifnot(!anyDuplicated(mutation$depmap_id),all(mutation$BRAF_hotspot %in% c(0,1)),
  all(c(external$depmap_id,external_molecular$depmap_id) %in% mutation$depmap_id))
write_table(mutation[depmap_id %in% c(external$depmap_id,external_molecular$depmap_id)],'external_BRAF_molecular_covariates')
for(name in c('external_combination_outcomes','external_ligand_programs','external_exact_drug_molecular_observations',
 'external_source_hashes','external_ligand_mapping','external_lincs_observation_operator')) {
  stopifnot(file.copy(file.path(external_root,'data',paste0(name,'.tsv')),
    file.path(root,'data',paste0(name,'.tsv')),overwrite=TRUE))
}
stopifnot(nrow(external)==156L,uniqueN(external$gene)==13L,
          !any(external$drug=='Lapatinib ditosylate'),!any(external$depmap_id %in% ligands$depmap_id))
pairs <- fread(file.path(parent,'data/combination_training_outcomes.tsv'))
folds <- fread(file.path(parent,'data/recipient_folds.tsv'))
stopifnot(nrow(pairs)==307L,nrow(folds)==62L,
          !any(pairs$depmap_id %in% ligands$depmap_id))
write_table(pairs,'combination_training_outcomes')
write_table(folds,'recipient_folds')
for(name in c('measured_predictions_frozen.tsv','baseline_predictions_frozen.tsv',
 'measured_decoder_correlations.tsv','baseline_bridge_correlations.tsv','selected_doses.tsv','molecular_contexts.tsv')) {
 source <- file.path(parent,'data',name)
 target <- file.path(root,'data',name)
 stopifnot(file.copy(source,target,overwrite=TRUE),identical(digest(source,file=TRUE),digest(target,file=TRUE)))
}

# Frozen single-effect evaluations at 0.5 and 5 uM interpolate in log dose to 0.75 uM.
# This interpolates functional evaluations and interaction activities, not gene displacements.
mema_features <- function(model) {
 cells <- ligands$depmap_id
 low <- predict_single(model,data.table(depmap_id=cells,drug='Lapatinib ditosylate',dose=.5))
 high <- predict_single(model,data.table(depmap_id=cells,drug='Lapatinib ditosylate',dose=5))
 weight <- log(.75/.5)/log(5/.5)
 L <- (1-weight)*low$program_activity+weight*high$program_activity
 S <- (1-weight)*low$sensitivity+weight*high$sensitivity
 D <- as.matrix(ligands[,..programs])
 ctx <- model$context
 Q <- cbind(1,sweep(sweep(sweep(D,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
 scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
                    crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
 raw_scale <- exp(Q %*% scale_map)
 minimum <- matrix(apply(ctx$sigma[1:10,],2,min),nrow(D),40,byrow=TRUE)
 maximum <- matrix(apply(ctx$sigma[1:10,],2,max),nrow(D),40,byrow=TRUE)
 scale <- pmin(pmax(raw_scale,minimum),maximum)
 amplitude <- as.numeric(exp(Q %*% model$par$log_amplitude))
 H <- (D/scale)*amplitude
 list(A=L,B=H,S=S)
}

center_geometry <- function(X,y,group) {
 levels <- unique(group)
 index <- match(group,levels)
 weight <- 1/as.numeric(table(group)[group])/length(levels)
 means <- vapply(seq_len(ncol(X)),function(k)
   as.numeric(tapply(X[,k],factor(group,levels=levels),mean)),numeric(length(levels)))
 ymean <- as.numeric(tapply(y,factor(group,levels=levels),mean))
 Xcentered <- X-means[index,,drop=FALSE]
 ycentered <- y-ymean[index]
 list(X=X,Xc=Xcentered,y=y,yc=ycentered,w=weight,
      variance=sum(weight*ycentered^2),group=group)
}

# Infer ONLY molecular coordinates for the two exact added inhibitors. Their
# functional outcomes cannot change f or any molecular map. LINCS covers partial
# program measurements: preparation records the landmark-loading assumption.
external_features <- function(original,held_cell=character()) {
  ctx <- original$context
  Aall <- Ball <- matrix(NA_real_,nrow(external),40)
  maps <- list()
  molecular_predictions <- list()
  for(dr in unique(external$drug)) {
    train <- external_molecular[drug==dr & !depmap_id %in% held_cell]
    stopifnot(nrow(train)>=8L,!anyDuplicated(train$depmap_id))
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
    Ad <- sweep(predicted,2,sigma,'/')*as.numeric(exp(Qd %*% original$par$log_amplitude))
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
    Ah <- (Dh/scale)*as.numeric(exp(Qh %*% original$par$log_amplitude))
    Aall[ii,] <- Ad
    Ball[ii,] <- Ah
    maps[[dr]] <- list(mu=mu,sigma=sigma,xmean=xm,B=B,training_cells=ids,
      BRAF_mean=braf_mean,BRAF_sd=braf_sd,lambda=1,
      conditioning='eight frozen baseline RNA directions plus BRAF-hotspot flag; functional labels excluded')
    molecular_predictions[[dr]] <- cbind(external[ii,.(depmap_id,gene,drug,dose)],as.data.table(predicted))
  }
  stopifnot(all(is.finite(Aall)),all(is.finite(Ball)))
  list(A=Aall,B=Ball,rna=maps,molecular_predictions=rbindlist(molecular_predictions))
}

# Only Gamma and a positive MEMA assay scale change. g is one shared rule across both sources.
fit_interaction <- function(original,mfeatures,endpoint_name,excluded_gene=character(),held_cell=character()) {
 train_pairs <- pairs[!depmap_id %in% held_cell]
 A <- predict_single(original,train_pairs[,.(depmap_id,drug=drugA,dose=doseA)])
 B <- predict_single(original,train_pairs[,.(depmap_id,drug=drugB,dose=doseB)])
 ef <- external_features(original,held_cell)
 ei <- which(!external$gene %in% excluded_gene & !external$depmap_id %in% held_cell)
 # Mean prototypes from training-only molecular inputs. Intervention labels
 # identify prototypes here, but are never passed into a prediction function.
 mr <- which(!ligands$gene %in% excluded_gene)
 molecular_rows <- rbind(A$program_activity,B$program_activity,mfeatures$B[mr,,drop=FALSE],
   ef$A[ei,,drop=FALSE],ef$B[ei,,drop=FALSE])
 groups <- c(paste0('drug_',train_pairs$drugA),paste0('drug_',train_pairs$drugB),
   paste0('ligand_',ligands$gene[mr]),paste0('drug_',external$drug[ei]),paste0('ligand_',external$gene[ei]))
 group_levels <- sort(unique(groups))
 centres <- t(vapply(group_levels,function(g)colMeans(molecular_rows[groups==g,,drop=FALSE]),numeric(40)))
 distances <- as.matrix(dist(centres))
 diag(distances) <- Inf
 bandwidth <- 2*median(apply(distances,1,min))
 stopifnot(is.finite(bandwidth),bandwidth>1e-8)
 basis <- list(centres=centres,labels=group_levels,bandwidth=bandwidth,
   method='molecular mean prototypes; two times median nearest-prototype spacing; no outcomes')
 Pn <- pair_features(A$program_activity,B$program_activity,original$context$W,basis)
 nci_scale <- exp(original$par$assay)[match(train_pairs$SCREENER,c('FF','FG'))]
 nci_group <- paste(train_pairs$pair,train_pairs$SCREENER,train_pairs$doseA,train_pairs$doseB,sep='|')
 ng <- center_geometry(Pn*nci_scale,-train_pairs$I,nci_group)
 # Select endpoint explicitly to avoid ambiguity between a column and the argument.
 train_labels <- training_labels[training_labels$endpoint==endpoint_name & !training_labels$gene %in% excluded_gene]
 row <- match(paste(train_labels$CellLine,train_labels$gene),paste(ligands$CellLine,ligands$gene))
 stopifnot(!anyNA(row),!any(train_labels$gene %in% challenge_genes),
           !any(train_labels$gene %in% excluded_gene),!any(train_pairs$depmap_id %in% held_cell))
 Pm <- pair_features(mfeatures$A,mfeatures$B,original$context$W,basis)
 mg <- center_geometry(Pm[row,,drop=FALSE],-train_labels$y,train_labels$CellLine)
 base <- mfeatures$S[row]
 # These are exact departure contrasts: log[(D+H)/D]-log[H/vehicle].
 # They remove drug and ligand main effects without an inferred assay baseline.
 Pe <- pair_features(ef$A,ef$B,original$context$W,basis)
 eg <- center_geometry(Pe[ei,,drop=FALSE],-external$observed_departure_log2[ei],
   paste(external$drug[ei],external$gene[ei],sep='|'))
 # Equal sources; molecular feature scaling uses training within-group variation only.
 product_scale <- sqrt((colSums(ng$Xc^2*ng$w)+colSums(mg$Xc^2*mg$w)+
   colSums(eg$X^2*eg$w))/3)
 stopifnot(all(product_scale>1e-10),ng$variance>1e-10,mg$variance>1e-10)
 Xn <- sweep(ng$Xc,2,product_scale,'/')
 Xm <- sweep(mg$Xc,2,product_scale,'/')
 raw_m <- sweep(mg$X,2,product_scale,'/')
 Xe <- sweep(eg$X,2,product_scale,'/')
 # Profile the exact joint objective over ONE positive study conversion.
 # ALMANAC targets are measured departures; MEMA targets are absolute conditional
 # drug/vehicle responses. No receiver intercept or held-out outcome enters fitting.
 gram_n <- source_weight*crossprod(Xn,Xn*ng$w)/ng$variance
 rhs_n <- source_weight*as.numeric(crossprod(Xn,ng$w*ng$yc))/ng$variance
 # External log2-count departure is carried in the common g units. No freely
 # fitted external multiplier can absorb all rescue and destroy identifiability.
 gram_n <- gram_n+source_weight*crossprod(Xe,Xe*eg$w)/eg$variance
 rhs_n <- rhs_n+source_weight*as.numeric(crossprod(Xe,eg$w*eg$y))/eg$variance
 gram_m <- source_weight*crossprod(Xm,Xm*mg$w)/mg$variance
 solve_at_scale <- function(log_scale) {
   s <- exp(log_scale)
   gram <- gram_n+s^2*gram_m+diag(interaction_penalty,ncol(raw_m))
   target <- mg$yc
   rhs <- rhs_n+source_weight*s*as.numeric(crossprod(Xm,mg$w*target))/mg$variance
   gamma <- as.numeric(solve(gram,rhs))
   nci_residual <- ng$yc-as.numeric(Xn %*% gamma)
   mema_residual <- mg$yc-s*as.numeric(Xm %*% gamma)
   loss <- source_weight*sum(ng$w*nci_residual^2)/ng$variance+
     source_weight*sum(eg$w*(eg$y-as.numeric(Xe %*% gamma))^2)/eg$variance+
     source_weight*sum(mg$w*mema_residual^2)/mg$variance+
     interaction_penalty*sum(gamma^2)+assay_penalty*log_scale^2
   list(loss=loss,Gamma=gamma,gram=gram,rhs=rhs)
 }
 optimum <- optimize(function(z)solve_at_scale(z)$loss,interval=c(-4,4),tol=1e-8)
 log_scale <- optimum$minimum
 assay_scale <- exp(log_scale)
 solution <- solve_at_scale(log_scale)
 Gamma <- solution$Gamma
 gram <- solution$gram
 rhs <- solution$rhs
 converged <- is.finite(optimum$objective) && abs(log_scale)<3.999
 iteration <- 1L
 normal_residual <- max(abs(as.numeric(gram %*% Gamma)-rhs))
 stopifnot(normal_residual<1e-8,converged)
 model <- original
 model$par$Gamma <- Gamma
 model$context$interaction_scale <- product_scale
 model$context$interaction_supported <- rep(TRUE,ncol(raw_m))
 model$mema_assay_scale <- assay_scale
 model$interaction_basis <- basis
 model$external_rna <- ef$rna
 model$interaction_training <- list(method='ALMANAC departures + within-receiver MEMA contrasts + absolute RAF/MEK departures',
   endpoint=endpoint_name,challenge_genes=challenge_genes,training_genes=unique(train_labels$gene),
   excluded_gene=excluded_gene,held_cell=held_cell,source_weight=source_weight,
   interaction_penalty=interaction_penalty,assay_penalty=assay_penalty,
   iterations=iteration,converged=converged,normal_residual=normal_residual,
   external_rows=length(ei),external_HGF_rows=sum(external$gene[ei]=='HGF'),
   external_unit_assumption='one common log2-count departure unit; no study-specific external multiplier')
 stopifnot(identical(original$rna,model$rna),identical(original$template,model$template),
           identical(original$par$level,model$par$level),
           identical(original$par$log_amplitude,model$par$log_amplitude),
           identical(original$par$assay,model$par$assay))
 list(model=model,training=train_labels,iterations=iteration,normal_residual=normal_residual,
      external_features=ef,external_training_rows=ei)
}

predict_mema <- function(model,features,rows) {
 products <- pair_features(features$A[rows,,drop=FALSE],features$B[rows,,drop=FALSE],model$context$W,model$interaction_basis)
 interaction <- as.numeric(sweep(products,2,model$context$interaction_scale,'/') %*% model$par$Gamma)
 data.table(CellLine=ligands$CellLine[rows],gene=ligands$gene[rows],
   predicted_log2=-features$S[rows]-model$mema_assay_scale*interaction,
   uncalibrated_drug_log2=-features$S[rows],interaction_log2=-interaction,
   mema_assay_scale=model$mema_assay_scale)
}

original_final <- qread(file.path(parent,'models/final.qs'))
features_final <- mema_features(original_final)
all_predictions <- nci_predictions <- audits <- coefficient_tables <- list()
final_models <- list()
for(endpoint in c('normalized_cell_count','raw_cell_count')) {
 genes <- sort(unique(training_labels$gene))
 for(job in seq_len(length(genes)+1L)) {
   held_gene <- if(job>length(genes))character() else genes[job]
   label <- if(length(held_gene))paste0('ligand_',held_gene) else 'final'
   fit <- fit_interaction(original_final,features_final,endpoint,held_gene)
   model <- fit$model
   rows <- if(length(held_gene))which(ligands$gene==held_gene) else seq_len(nrow(ligands))
   output <- predict_mema(model,features_final,rows)
   output[,`:=`(endpoint=endpoint,evaluation=if(length(held_gene))'held_gene' else 'fitted',fit=label)]
   if(!length(held_gene))output[gene %in% challenge_genes,evaluation:='challenge']
   all_predictions[[paste(endpoint,label)]] <- output
   audits[[paste(endpoint,label)]] <- data.table(endpoint,fit=label,excluded_gene=if(length(held_gene))held_gene else NA_character_,
      training_ligands=uniqueN(fit$training$gene),training_contexts=nrow(fit$training),
      challenge_leakage=sum(fit$training$gene %in% challenge_genes),
      external_HGF_training_rows=sum(external$gene[fit$external_training_rows]=='HGF'),iterations=fit$iterations,
      normal_residual=fit$normal_residual,mema_assay_scale=model$mema_assay_scale)
   if(!length(held_gene)) {
      final_models[[endpoint]] <- model
      qsave(model,file.path(root,'models',paste0('final_',endpoint,'.qs')))
   }
 }
 cat(endpoint,': 26 complete ligand exclusions plus final fit complete\n')
}
# External whole-cell validation. A cell's ligand+drug outcomes never train its g.
external_predictions <- list()
for(cell in sort(unique(external$depmap_id))) {
  old_index <- match(cell,folds$depmap_id)
  original <- if(is.na(old_index))original_final else qread(file.path(parent,'models',sprintf('loo_%02d.qs',old_index)))
  fit <- fit_interaction(original,mema_features(original),primary_endpoint,held_cell=cell)
  rows <- which(external$depmap_id==cell)
  gamma <- fit$model$par$Gamma
  scale <- fit$model$context$interaction_scale
  products <- pair_features(fit$external_features$A[rows,,drop=FALSE],fit$external_features$B[rows,,drop=FALSE],
    fit$model$context$W,fit$model$interaction_basis)
  predicted <- -as.numeric(sweep(products,2,scale,'/') %*% gamma)
  external_predictions[[cell]] <- external[rows][,predicted_departure_log2:=predicted]
}
external_held <- rbindlist(external_predictions)
write_table(external_held,'external_held_cell_predictions_frozen')
external_metrics <- external_held[,.(ligands=.N,r_departure=pearson(observed_departure_log2,predicted_departure_log2),
  observed_HGF_departure=observed_departure_log2[gene=='HGF'],
  predicted_HGF_departure=predicted_departure_log2[gene=='HGF']),by=.(CellLine,drug,dose)]
write_table(external_metrics,'external_held_cell_correlations')
write_table(external_features(original_final)$molecular_predictions,'external_predicted_drug_programs')
# Retain primary models under the established filenames; recompute ALMANAC whole-cell predictions.
for(job in seq_len(nrow(folds)+1L)) {
 held <- if(job>nrow(folds))character() else folds$depmap_id[job]
 label <- if(length(held))sprintf('loo_%02d',job) else 'final'
 original <- qread(file.path(parent,'models',paste0(label,'.qs')))
 if(length(held)) {
   features <- mema_features(original)
   fit <- fit_interaction(original,features,primary_endpoint,held_cell=held)
   model <- fit$model
 } else model <- final_models[[primary_endpoint]]
 qsave(model,file.path(root,'models',paste0(label,'.qs')))
 test <- if(length(held))pairs[depmap_id==held] else pairs
 if(nrow(test)) {
   predicted <- predict_pair(model,test)
   output <- test[,.(depmap_id,pair,SCREENER,drugA,drugB,doseA,doseB,
      observed_log2=y,observed_departure_log2=I)]
   output[,`:=`(predicted_log2=-predicted$sensitivity,predicted_departure_log2=-predicted$departure,
                evaluation=if(length(held))'held_cell' else 'fitted',fit=label)]
   nci_predictions[[label]] <- output
 }
 coefficient_tables[[label]] <- data.table(fit=label,
    term=paste0(rep(c('route1_','route2_'),each=length(model$interaction_basis$labels)),rep(model$interaction_basis$labels,2)),
    coefficient_sensitivity=model$par$Gamma,product_scale=model$context$interaction_scale,
    mema_assay_scale=model$mema_assay_scale)
 if(job %% 20L==0L)cat('ALMANAC cell exclusions completed:',job,'of',nrow(folds)+1L,'\n')
}
mema_pred <- rbindlist(all_predictions)
combo <- rbindlist(nci_predictions)
write_table(mema_pred,'mema_predictions_frozen')
write_table(combo,'combination_predictions_frozen')
write_table(rbindlist(audits),'fit_audit')
write_table(rbindlist(coefficient_tables),'shared_interaction_coefficients')
# Original v107 HGF vectors and standalone effects remain unchanged for the 5 uM query.
model <- final_models[[primary_endpoint]]
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
canonical <- ligands[CellLine=='AU565' & gene=='HGF',..programs]
canonical[,context:='Study-balanced HGF (same recipe as training)']
setcolorder(canonical,c('context',programs))
hgf <- rbind(canonical,hgf)
hmean <- as.numeric(hgf[context=='Study-balanced HGF (same recipe as training)',..programs])
ctx <- model$context
Dh <- as.matrix(hgf[,..programs])
Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
minimum <- matrix(apply(ctx$sigma[1:10,],2,min),nrow(Dh),40,byrow=TRUE)
maximum <- matrix(apply(ctx$sigma[1:10,],2,max),nrow(Dh),40,byrow=TRUE)
hscale <- pmin(pmax(exp(Qh %*% scale_map),minimum),maximum)
lap_meta <- data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=5)
lap <- predict_single(model,lap_meta)
lap_pct <- 100*2^(-lap$sensitivity)
Hactivity <- ((Dh/hscale) %*% ctx$W)*as.numeric(exp(Qh %*% model$par$log_amplitude))
Lactivity <- lap$interaction_activity[rep(1L,nrow(hgf)),,drop=FALSE]
Lprogram <- lap$program_activity[rep(1L,nrow(hgf)),,drop=FALSE]
Hprogram <- (Dh/hscale)*as.numeric(exp(Qh %*% model$par$log_amplitude))
J <- as.numeric(sweep(pair_features(Lprogram,Hprogram,ctx$W,model$interaction_basis),2,ctx$interaction_scale,'/') %*% model$par$Gamma)
hgf_query <- rbindlist(lapply(1:2,function(route)data.table(cell='AU565',dose=5,hgf_context=hgf$context,assumed_HGF_route=route,
 lapatinib_pct_control=lap_pct,hgf_increment_pct=100*2^(-Hactivity[,route]),
 combination_pct_control=lap_pct*2^(-Hactivity[,route]-J),
 lapatinib_pct_hgf_control=lap_pct*2^(-J),interaction_log2=-J,
 interaction_relative_pct=100*(2^(-J)-1),
 qualification='Lapatinib+HGF and Wilson labels excluded; HGF+RAF/MEK trains g; standalone HGF increment uncalibrated.')))
write_table(hgf_query,'au565_hgf_5uM_predictions_frozen')
freeze_paths <- file.path(root,'data',c('mema_predictions_frozen.tsv','combination_predictions_frozen.tsv','au565_hgf_5uM_predictions_frozen.tsv'))
write_table(data.table(path=freeze_paths,sha256=vapply(freeze_paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')
# Challenge and ordinary validation labels are joined only after prediction freeze.
scored <- merge(mema_pred,labels[,.(CellLine,gene,endpoint,observed_log2=y,ecm_n)],
                by=c('CellLine','gene','endpoint'),all.x=TRUE)
stopifnot(!anyNA(scored$observed_log2))
write_table(scored,'mema_predictions_scored')
metrics <- scored[,.(ligands=.N,r=pearson(observed_log2,predicted_log2),
 observed_median_pct=median(100*2^observed_log2),predicted_median_pct=median(100*2^predicted_log2)),
 by=.(endpoint,CellLine,evaluation)]
write_table(metrics,'mema_correlations')
challenge <- scored[evaluation=='challenge'][,.(endpoint,CellLine,gene,
 measured_pct_ligand_control=100*2^observed_log2,predicted_pct_ligand_control=100*2^predicted_log2,
 interaction_relative_pct=100*(2^interaction_log2-1),mema_assay_scale)]
write_table(challenge,'held_HGF_NRG1_predictions')
combo_metrics <- combo[,.(cells=.N,r_total=pearson(observed_log2,predicted_log2),
 r_departure=pearson(observed_departure_log2,predicted_departure_log2),
 observed_departure_SD=sd(observed_departure_log2),predicted_departure_SD=sd(predicted_departure_log2)),
 by=.(pair,SCREENER,evaluation)]
write_table(combo_metrics,'combination_correlations')
# Square, independently scaled panels: comparisons across ligands in one recipient at fixed 0.75 uM.
plot_data <- scored[evaluation=='held_gene']
plot_data[,panel:=paste(CellLine,ifelse(endpoint==primary_endpoint,'Normalized count','Raw count'),sep=' | ')]
plot_stats <- plot_data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed_log2,predicted_log2)),'\nn = ',.N,' ligands')),by=panel]
p <- ggplot(plot_data,aes(observed_log2,predicted_log2))+geom_point(color='#0072B2',size=1.5)+
 geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
 geom_text(data=plot_stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
 facet_wrap(~panel,ncol=2,scales='free')+theme_bw(base_size=9)+theme(aspect.ratio=1)+
 labs(title='v115 | MEMA leave-one-ligand-out prediction',x='Measured log2 drug/vehicle under the same ligand',
 y='Predicted log2 drug/vehicle under the same ligand',
 caption='Lapatinib 0.75 µM, 72 h. Each point is one ligand averaged over ECM conditions.\nMEMA HGF/NRG1 excluded. HGF+RAF/MEK trains g. All Wilson outcomes excluded. Retrospective evaluation.')
ggsave(file.path(root,'figures/mema_held_ligand_scatters.pdf'),p,width=7,height=7.1,device=cairo_pdf)
ggsave(file.path(root,'figures/mema_held_ligand_scatters.png'),p,width=7,height=7.1,dpi=160)
writeLines(c('v115: same g and molecular predictors as v114; MEMA trains within-receiver ligand contrasts only.',
 'MEMA has no released zero-ligand control. Its receiver means cannot identify absolute g plus a study baseline.',
 'Independent RAF/MEK exact departure contrasts identify absolute g. MEMA contributes variation across ligands.',
 'Positive MEMA contrast scale is fitted jointly from ligand variation only. It does not recalibrate the frozen drug head.',
 'MEMA absolute drug baseline remains the unadjusted PRISM prediction: absolute study transfer is not validated.',
 'BRAF status is an RNA-stage covariate, never a direct functional predictor. Original lapatinib and all single-drug fits are unchanged.',
 'Prototypes and bandwidth use only characterized training molecular inputs. Shared ridge coefficients are outcome-learned.',
 'Existing molecular maps, signatures, single-effect calibrations and NCI assay scales are frozen from v107.',
 'PLX-4720 and PD-184352 are exact compounds with local 24h LINCS RNA; no vemurafenib substitute.',
 'Added molecular-only baseline maps predict external drug programs at 2 and 1 uM; landmark-span approximation is recorded.',
 'Straussman2012 Table8: fixed 25ng/mL, 13 supported ligands, six melanoma cells, 156 exact departure observations.',
 'HGF+RAF/MEK trains g. Lapatinib+HGF, MEMA HGF/NRG1, and all Wilson outcomes remain excluded. Unseen-combination transfer.',
 'External departure unit is common log2-count fold; this assay-transfer assumption is not independently calibrated.',
 'MEMA: lapatinib 0.75 uM, 72 h; normalized count primary, raw count mandatory sensitivity.',
 'Molecular and endpoint aggregation reused from prior preparation; 26 genes x two recipients train; MEMA HGF and NRG1 held out.',
 'MEMA endpoint is drug/vehicle under identical ligand and ECM. It is not an absolute no-ligand rescue endpoint.',
 'Drug functional evaluations/activities at 0.5 and 5 uM interpolate in log dose; no extrapolation or molecular-vector interpolation.',
 'g uses both existing signature activities for every constituent. No ligand identity or recipient intercept enters g.',
 'Center ALMANAC departures/products within pair; center MEMA features and outcomes within receiver; preserve absolute external departures.',
 'Equal source weights .2; ridge .05; source variance normalization. One positive MEMA study scale, log-scale penalty .01.',
 'Gamma and study conversion minimize one objective, profiled over a positive scale. No receiver/ligand intercepts. Fixed before scoring.',
 '26 leave-gene-out fits per endpoint remove each held gene from external outcomes too. ALMANAC62 and external6 whole-cell exclusions.',
 'Primary AU565 5uM HGF query uses the identical study-balanced geneFC-then-transport recipe used in g training.',
 'Legacy average-after-transport HGF and all source contexts remain diagnostic. Ligand-growth increment is uncalibrated. Wilson not fitted.',
 'MEMA and Wilson outcomes have previously been exposed; these are retrospective development results.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
input_paths <- c('/work/src/program_pilot/115_departure_identified_bypass_transfer.R',mutation_path,
 file.path(mema_root,c('all_gene_endpoints.tsv','recipient_ligand_program_scores.tsv')),
 file.path(parent,'models/final.qs'),archive_path)
write_table(data.table(path=input_paths,sha256=vapply(input_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(external_metrics)
print(metrics[evaluation=='held_gene'])
print(challenge)
print(hgf_query[hgf_context=='Study-balanced HGF (same recipe as training)'])
cat('Complete. Single-drug outputs byte-identical to v107; lapatinib+HGF and all Wilson outcomes excluded.\n')
