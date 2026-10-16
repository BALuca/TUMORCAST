#!/usr/bin/env Rscript
# Shared g trained across HGF molecular source views; ligand-single conversion measured across all ligands.
.libPaths(c("/work/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest);library(glmnet)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
parent <- "/work/results/program_pilot/107_within_pair_departure_interaction/linear"
root <- "/work/results/program_pilot/133_counts_and_growth_order_bridge/linear"
time_matched_root <- '/work/results/program_pilot/118_time_matched_HGF_bypass_transfer/linear'
external_root <- '/work/results/program_pilot/111_external_bypass_combination_transfer/linear'
mema_root <- "/work/grant_rebuild_2026-10-04/transfer_solution"
challenge_genes <- c("HGF","NRG1")
primary_endpoint <- "normalized_cell_count"
source_weight <- .2
interaction_penalty <- .05
assay_penalty <- .01
external_within_weight <- .8
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
    stopifnot(ncol(A)==43L,ncol(B)==43L,nrow(A)==nrow(B),all(dim(W)==c(40L,2L)))
    rawA <- A[,1:40,drop=FALSE]
    rawB <- B[,1:40,drop=FALSE]
    SA <- rawA %*% W
    SB <- rawB %*% W
    # RBF movement relative to zero: each modifier is exactly zero at no input.
    # Centers and length scale depend only on characterized molecular profiles.
    transform <- function(X) {
      distance <- pmax(outer(rowSums(X^2),rowSums(basis$centres^2),'+')-2*tcrossprod(X,basis$centres),0)
      zero <- exp(-rowSums(basis$centres^2)/(2*basis$bandwidth^2))
      sweep(exp(-distance/(2*basis$bandwidth^2)),2,zero,'-')
    }
    KA <- transform(rawA)
    KB <- transform(rawB)
    gateA <- -expm1(-rowSums(rawA^2)/(2*basis$bandwidth^2))
    gateB <- -expm1(-rowSums(rawB^2)/(2*basis$bandwidth^2))
    absolute <- do.call(cbind,lapply(1:2,function(k)KB*SA[,k]+KA*SB[,k]))
    modulation <- do.call(cbind,lapply(1:2,function(k)
      KB*(gateA*A[,40+k])+KA*(gateB*B[,40+k])))
    features <- cbind(absolute,modulation)
    colnames(features) <- paste0(rep(c('absolute1_','absolute2_','modulation1_','modulation2_'),each=nrow(basis$centres)),rep(basis$labels,4))
    # Four direct shared interactions retain the independently defined RNA
    # programme; its information cannot disappear inside a high-dimensional RBF.
    consensus_absolute <- do.call(cbind,lapply(1:2,function(k)
      SA[,k]*B[,43]+SB[,k]*A[,43]))
    consensus_modulation <- do.call(cbind,lapply(1:2,function(k)
      gateA*A[,40+k]*B[,43]+gateB*B[,40+k]*A[,43]))
    consensus <- cbind(consensus_absolute,consensus_modulation)
    colnames(consensus)<-c('consensus_absolute1','consensus_absolute2','consensus_modulation1','consensus_modulation2')
    features <- cbind(features,consensus)
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
        program_activity = cbind(D/model$context$sigma[g,,drop=FALSE],node$both,
          predict_auxiliary_coordinate(model,meta$depmap_id,meta$drug,dose)),D = D)
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
ligands <- fread(file.path(root,'data/time_matched_mema_molecular_inputs.tsv'))
stopifnot(nrow(ligands)==56L,!anyDuplicated(ligands[,.(CellLine,gene)]),
          all(is.finite(as.matrix(ligands[,..programs]))),
          all(labels$ecm_n==48L),all(labels$endpoint %in% c('normalized_cell_count','raw_cell_count')))
# Source HGF labels are not joined into training; challenge scoring occurs after freeze.
training_labels <- labels[!gene %in% challenge_genes | (CellLine=='HCC1954' & gene=='HGF')]
stopifnot(uniqueN(training_labels$gene)==27L,nrow(training_labels)==106L,
  !any(training_labels$CellLine=='AU565' & training_labels$gene=='HGF'),
  !any(training_labels$gene=='NRG1'))
write_table(training_labels,'mema_nonchallenge_training_outcomes')
write_table(ligands,'mema_molecular_inputs')
external <- fread('/work/results/program_pilot/123_HER2_ligand_combination_transfer/linear/data/expanded_external_combination_outcomes.tsv')
external_ligands <- fread(file.path(root,'data/expanded_external_ligand_programs.tsv'))
for(name in c('HGF_time_context_eligibility','time_matched_HGF_source_hashes')) {
  stopifnot(file.copy(file.path(time_matched_root,'data',paste0(name,'.tsv')),
    file.path(root,'data',paste0(name,'.tsv')),overwrite=TRUE))
}
source_views <- fread(file.path(root,'data/time_matched_HGF_program_inputs.tsv'))[
  !grepl('^Leave_out_|^Primary_',source_profile)]
stopifnot(uniqueN(source_views$source_profile)==4L)
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
write_table(external,'external_combination_outcomes')
write_table(external_ligands,'external_ligand_programs')
stopifnot(nrow(external)==172L,uniqueN(external$gene)==14L,
          !any(external$drug=='Lapatinib ditosylate' & external$gene=='HGF'),
          !any(external$depmap_id %in% ligands$depmap_id))
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
 H <- D/scale
 H <- cbind(H,H %*% ctx$W,ligand_coordinate(ligands$CellLine,ligands$gene))
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
  Aall <- Ball <- matrix(NA_real_,nrow(external),43)
  maps <- list()
  molecular_predictions <- list()
  for(dr in unique(external$drug)) {
    ii <- which(external$drug==dr)
    if(dr=='Lapatinib ditosylate') {
      native <- predict_single(original,external[ii,.(depmap_id,drug,dose)])
      Ad <- native$program_activity[,1:40,drop=FALSE]
      centred_A <- native$program_activity[,41:42,drop=FALSE]
      auxiliary_prediction <- native$program_activity[,43]
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
    am <- merge(train[,.(depmap_id,drug,dose)],auxiliary_external,by=c('depmap_id','drug','dose'),all.x=TRUE)
    am <- am[match(train$depmap_id,am$depmap_id)]
    stopifnot(!anyNA(am$HGF_RNA_coordinate))
    auxiliary_mean <- mean(am$HGF_RNA_coordinate)
    auxiliary_sd <- max(sd(am$HGF_RNA_coordinate),1e-8)
    auxiliary_B <- as.numeric(solve(crossprod(X)/nrow(X)+diag(1,ncol(X)),
      crossprod(X,(am$HGF_RNA_coordinate-auxiliary_mean)/auxiliary_sd)/nrow(X)))
    auxiliary_prediction <- (auxiliary_mean+as.numeric(sweep(Xq,2,xm) %*% auxiliary_B)*auxiliary_sd)/auxiliary_sd
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
    Aall[ii,] <- cbind(Ad,centred_A,auxiliary_prediction)
    Ball[ii,] <- cbind(Ah,Ah %*% ctx$W,ligand_coordinate(external$CellLine[ii],external$gene[ii]))
    colnames(predicted) <- programs
    molecular_predictions[[dr]] <- cbind(external[ii,.(depmap_id,gene,drug,dose)],as.data.table(predicted))
  }
  stopifnot(all(is.finite(Aall)),all(is.finite(Ball)))
  # Alternative source views express uncertainty about the same physical
  # experiment. Four HGF RNA studies share ONE total experiment weight.
  view_rows <- integer()
  view_labels <- character()
  Aviews <- Bviews <- matrix(numeric(),0,43)
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
    input_names <- if(external$gene[i]=='HGF')labels else external$gene[i]
    Bviews <- rbind(Bviews,cbind(H,H %*% ctx$W,
      ligand_coordinate(rep(external$CellLine[i],nrow(Dh)),input_names)))
    standalone_views <- c(standalone_views,standalone)
    view_rows <- c(view_rows,rep(i,nrow(Dh)))
    view_labels <- c(view_labels,labels)
  }
  list(A=Aall,B=Ball,Aviews=Aviews,Bviews=Bviews,view_rows=view_rows,view_labels=view_labels,
    standalone_views=standalone_views,rna=maps,molecular_predictions=rbindlist(molecular_predictions))
}


# Only Gamma and a positive MEMA assay scale change. g is one shared rule across both sources.

# Inner validation with whole external receiver groups. All lapatinib+HGF and
# Wilson labels remain excluded. Outcomes excluded by an outer job also remain
# excluded throughout its penalty choice. No rescue-threshold selection.
penalty_audit <- list()
fit_selected <- function(original,mfeatures,endpoint_name,excluded_gene=character(),held_cell=character()) {
  penalty <- select_penalty(original,mfeatures,endpoint_name,excluded_gene,held_cell)
  fit_interaction(original,mfeatures,endpoint_name,excluded_gene,held_cell,penalty)
}


auxiliary_root <- file.path(root,'data')
auxiliary_molecular <- fread(file.path(auxiliary_root,"tahoe_HGF_coordinate_molecular_observations.tsv"))
auxiliary_external <- fread(file.path(auxiliary_root,"external_exact_drug_HGF_coordinate.tsv"))
auxiliary_ligands <- fread(file.path(auxiliary_root,"ligand_HGF_coordinate_molecular_inputs.tsv"))
auxiliary_ligand_scale <- sd(auxiliary_ligands$HGF_RNA_coordinate[!grepl("^Leave_out_",auxiliary_ligands$input)])
auxiliary_key <- function(cell, input) paste(cell, input, sep = "|")
ligand_coordinate <- function(cells, inputs) {
    rows <- match(auxiliary_key(cells, inputs), auxiliary_key(auxiliary_ligands$CellLine, auxiliary_ligands$input))
    stopifnot(!anyNA(rows))
    auxiliary_ligands$HGF_RNA_coordinate[rows]/auxiliary_ligand_scale
}
fit_auxiliary_rna <- function(original, held_cells = character()) {
    maps <- list()
    for (i in seq_len(nrow(contexts))) {
        dr <- contexts$drug[i]
        ds <- contexts$dose[i]
        train <- auxiliary_molecular[drug == dr & abs(dose - ds) < 1e-10 & !depmap_id %in% c("ACH-000248", held_cells)]
        stopifnot(nrow(train) >= 3L, !anyDuplicated(train$depmap_id))
        if (dr == "Lapatinib ditosylate") {
            X <- as.matrix(receptors[match(train$depmap_id, depmap_id), .(ERBB2, MET)])
            kind <- "ERBB2_MET"
            xs <- apply(X, 2, sd)
        }
        else {
            rna <- original$rna
            raw <- as.matrix(baseline[match(train$depmap_id, depmap_id), ..programs])
            X <- sweep(sweep(sweep(raw, 2, rna$bmean), 2, rna$bsd, "/") %*% rna$basis, 2, rna$bscale, "/")
            kind <- "baseline_programs"
            xs <- rep(1, ncol(X))
        }
        xm <- colMeans(X)
        X <- sweep(sweep(X, 2, xm), 2, xs, "/")
        dm <- mean(train$HGF_RNA_coordinate)
        ds_response <- max(sd(train$HGF_RNA_coordinate), 1e-08)
        y <- (train$HGF_RNA_coordinate - dm)/ds_response
        B <- as.numeric(solve(crossprod(X)/nrow(X) + diag(1, ncol(X)), crossprod(X, y)/nrow(X)))
        maps[[context_key(dr, ds)]] <- list(kind = kind, xmean = xm, xsd = xs, mu = dm, sigma = ds_response, B = B, 
            training_cells = train$depmap_id, lambda = 1)
    }
    maps
}
predict_auxiliary_coordinate <- function(model, cells, drug, dose) {
    value <- numeric(length(cells))
    groups <- context_key(drug, dose)
    for (group in unique(groups)) {
        rows <- which(groups == group)
        map <- model$auxiliary_rna[[group]]
        stopifnot(!is.null(map))
        if (map$kind == "ERBB2_MET") 
            X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
        else {
            rna <- model$rna
            raw <- as.matrix(baseline[match(cells[rows], depmap_id), ..programs])
            X <- sweep(sweep(sweep(raw, 2, rna$bmean), 2, rna$bsd, "/") %*% rna$basis, 2, rna$bscale, "/")
        }
        X <- sweep(sweep(X, 2, map$xmean), 2, map$xsd, "/")
        value[rows] <- (map$mu + as.numeric(X %*% map$B) * map$sigma)/map$sigma
    }
    stopifnot(all(is.finite(value)))
    value
}
source("/work/src/program_pilot/133_assay_observation_fit.R")
original_final <- qread(file.path(parent,'models/final.qs'))
original_final$auxiliary_rna <- fit_auxiliary_rna(original_final)
features_final <- mema_features(original_final)
all_predictions <- nci_predictions <- audits <- coefficient_tables <- list()
final_models <- list()
for(endpoint in c('normalized_cell_count','raw_cell_count')) {
 genes <- sort(unique(training_labels$gene))
 for(job in seq_len(length(genes)+1L)) {
   held_gene <- if(job>length(genes))character() else genes[job]
   label <- if(length(held_gene))paste0('ligand_',held_gene) else 'final'
   fit <- fit_selected(original_final,features_final,endpoint,held_gene)
   model <- fit$model
   rows <- if(length(held_gene))which(ligands$gene==held_gene) else seq_len(nrow(ligands))
   output <- predict_mema(model,features_final,rows)
   output[,`:=`(endpoint=endpoint,evaluation=if(length(held_gene))'held_gene' else 'fitted',fit=label)]
   if(!length(held_gene))output[gene %in% challenge_genes & !(CellLine=='HCC1954' & gene=='HGF'),evaluation:='challenge']
   all_predictions[[paste(endpoint,label)]] <- output
   audits[[paste(endpoint,label)]] <- data.table(endpoint,fit=label,excluded_gene=if(length(held_gene))held_gene else NA_character_,
      training_ligands=uniqueN(fit$training$gene),training_contexts=nrow(fit$training),
      challenge_leakage=sum(fit$training$gene=='NRG1' | (fit$training$CellLine=='AU565' & fit$training$gene=='HGF')),
      HCC1954_HGF_training_rows=sum(fit$training$CellLine=='HCC1954' & fit$training$gene=='HGF'),
      external_HGF_training_rows=sum(external$gene[fit$external_training_rows]=='HGF'),iterations=fit$iterations,
      normal_residual=fit$normal_residual,mema_assay_scale=model$mema_assay_scale)
   if(!length(held_gene)) {
      final_models[[endpoint]] <- model
      qsave(model,file.path(root,'models',paste0('final_',endpoint,'.qs')))
   }
 }
 cat(endpoint,': 27 complete ligand exclusions plus final fit complete\n')
}
# External whole-cell validation. A cell's ligand+drug outcomes never train its g.
external_predictions <- list()
for(cell in sort(unique(external$depmap_id))) {
  old_index <- match(cell,folds$depmap_id)
  original <- if(is.na(old_index))original_final else qread(file.path(parent,'models',sprintf('loo_%02d.qs',old_index)))
  original$auxiliary_rna <- fit_auxiliary_rna(original,cell)
  fit <- fit_selected(original,mema_features(original),primary_endpoint,held_cell=cell)
  rows <- which(external$depmap_id==cell)
  gamma <- fit$model$par$Gamma
  scale <- fit$model$context$interaction_scale
  predicted <- predict_external_departure(fit$model,fit$external_features)[rows]
  external_predictions[[cell]] <- external[rows][,`:=`(predicted_departure_log2=predicted,
    observed_scoring_readout=ifelse(study=='elife32271',observed_departure_log2,observed_drug_arm_increment_log2),
    comparison_endpoint=ifelse(study=='elife32271','72h final-cell-abundance departure',
      'within-cell/drug ligand ordering: reported net-growth versus predicted count increment'))]
}
external_held <- rbindlist(external_predictions)
write_table(external_held,'external_held_cell_predictions_frozen')
external_metrics <- external_held[,.(ligands=.N,r=pearson(observed_scoring_readout,predicted_departure_log2),
  observed_HGF_readout=if(any(gene=='HGF'))observed_scoring_readout[gene=='HGF'] else NA_real_,
  predicted_HGF_count_score=if(any(gene=='HGF'))predicted_departure_log2[gene=='HGF'] else NA_real_),
  by=.(CellLine,drug,dose,comparison_endpoint)]
write_table(external_metrics,'external_held_cell_correlations')
write_table(external_held[study=='elife32271',.(cells=.N,
  r_departure_across_cells=pearson(observed_departure_log2,predicted_departure_log2),
  observed_median_departure=median(observed_departure_log2),predicted_median_departure=median(predicted_departure_log2)),
  by=.(drug,dose,gene)],'HER2_external_held_cell_correlations')
write_table(external_features(original_final)$molecular_predictions,'external_predicted_drug_programs')
# Retain primary models under the established filenames; recompute ALMANAC whole-cell predictions.
for(job in seq_len(nrow(folds)+1L)) {
 held <- if(job>nrow(folds))character() else folds$depmap_id[job]
 label <- if(length(held))sprintf('loo_%02d',job) else 'final'
 original <- qread(file.path(parent,'models',paste0(label,'.qs')))
 original$auxiliary_rna <- fit_auxiliary_rna(original,held)
 if(length(held)) {
   features <- mema_features(original)
   fit <- fit_selected(original,features,primary_endpoint,held_cell=held)
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
    term=c(paste0(rep(c('absolute1_','absolute2_','modulation1_','modulation2_'),each=length(model$interaction_basis$labels)),rep(model$interaction_basis$labels,4)),
      'consensus_absolute1','consensus_absolute2','consensus_modulation1','consensus_modulation2'),
    coefficient_sensitivity=model$par$Gamma,product_scale=model$context$interaction_scale,
    mema_assay_scale=model$mema_assay_scale)
 if(job %% 20L==0L)cat('ALMANAC cell exclusions completed:',job,'of',nrow(folds)+1L,'\n')
}
mema_pred <- rbindlist(all_predictions)
combo <- rbindlist(nci_predictions)
write_table(mema_pred,'mema_predictions_frozen')
write_table(combo,'combination_predictions_frozen')
write_table(rbindlist(audits),'fit_audit')
write_table(rbindlist(penalty_audit),'training_only_penalty_selection')
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
Hactivity <- ((Dh/hscale) %*% ctx$W)*as.numeric(exp(Qh %*% model$par$log_amplitude))*model$ligand_single_conversion
Lactivity <- lap$interaction_activity[rep(1L,nrow(hgf)),,drop=FALSE]
Lprogram <- lap$program_activity[rep(1L,nrow(hgf)),,drop=FALSE]
Hprogram <- Dh/hscale
Hprogram <- cbind(Hprogram,Hprogram %*% ctx$W)
# Only primary HGF has the new RNA coordinate in this query. The four source
# studies and source omissions are interrogated by the separate validation script.
hgf <- hgf[1]
Hactivity <- Hactivity[1,,drop=FALSE]
Lprogram <- lap$program_activity
Hprogram <- cbind(Hprogram[1,,drop=FALSE],ligand_coordinate('AU565','HGF'))
J <- as.numeric(sweep(pair_features(Lprogram,Hprogram,ctx$W,model$interaction_basis),2,ctx$interaction_scale,'/') %*% model$par$Gamma)
hgf_query <- rbindlist(lapply(1:2,function(route)data.table(cell='AU565',dose=5,hgf_context=hgf$context,assumed_HGF_route=route,
 lapatinib_pct_control=lap_pct,hgf_increment_pct=100*2^(-Hactivity[,route]),
 combination_pct_control=lap_pct*2^(-Hactivity[,route]-J),
 lapatinib_pct_hgf_control=lap_pct*2^(-J),interaction_log2=-J,
 interaction_relative_pct=100*(2^(-J)-1),
 qualification='HCC1954 lapatinib+HGF trains g; AU565+HGF and all Wilson labels excluded. Novel-recipient transfer, not an unseen pair.')))
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
 labs(title='v133 | MEMA leave-one-ligand-out prediction',x='Measured log2 drug/vehicle under the same ligand',
 y='Predicted log2 drug/vehicle under the same ligand',
 caption='Lapatinib 0.75 µM, 72 h. Each point is one ligand averaged over ECM conditions.\nHCC1954 HGF can train final g; whole-ligand exclusions remove it. AU565 HGF and all Wilson outcomes excluded. Retrospective.')
ggsave(file.path(root,'figures/mema_held_ligand_scatters.pdf'),p,width=7,height=7.1,device=cairo_pdf)
ggsave(file.path(root,'figures/mema_held_ligand_scatters.png'),p,width=7,height=7.1,dpi=160)
writeLines(c("v133: same molecular f/g architecture and training scope as132; correct functional observation units.",
"g uses five-day PRISM-equivalent log-count units.72h cell-count observations use3/5 exposure conversion in training and prediction.",
"Straussman ligand assay is day1-to-day7, baseline-subtracted growth. Its normalized drug-arm readout supplies within-cell/drug order only.",
"No claim that its growth-fold departure is an absolute cell-count interaction. Mixed count regression and growth-order logistic likelihood fit the same Gamma.",
"One positive ranking temperature is a readout nuisance; absolute magnitude is anchored only to proper count endpoints. No growth fold calibrates HGF rescue magnitude.",
"Time projection assumes constant relative net growth after the molecular response; long-term state adaptation remains unmodeled.",
"AU565 HGF, all MEMA NRG1 and all Wilson outcomes remain excluded. Test is novel-recipient transfer, not an unseen drug-ligand pair.",
"One prior held-out HCC1954 HGF aggregate per endpoint becomes training. The endpoint remains matched drug/vehicle under identical HGF and ECM.",
"Illumina: normexp background correction plus quantile log2 normalization. Affymetrix: released RMA log2 values.",
"Within the one study family, platforms weighted by inverse replicate-contrast-noise variance. Source-study means retain equal family weights.",
"RNA coordinate weights, native drug RNA measurements, all functional outcomes and all modeling choices remain unchanged from v129.",
"Original40, their two sensitivity signatures and all single-drug outputs are unchanged. New programme is defined from RNA only and orthogonal to the original40.",
"Auxiliary molecular maps never use functional outcomes. New programme enters only g through symmetric zero-anchored products, not an HGF-specific coefficient.",
"Historical Topotecan gene-cache plate pooling differs from old40; this limitation remains explicit.",
"Drug modulation scores are exactly the two RNA-centred signatures from frozen f. Ligand scores retain their zero reference; HGF mean response is never subtracted.",
"A molecular norm gate makes all modulation terms exactly zero for absent intervention. Same RNA-only RBF width, training objective and inner validation as v127.",
"Added NRG1 in6 other breast cancer lines and EGF in2, at native0.5/5uM. All AU565, lapatinib+HGF and Wilson outcomes excluded.", 
"Numeric Prism XML from eLife32271. Paired replicates; log-abundance/log-dose interpolation. No extrapolation.", 
"First0.0001 graph abscissa is the vehicle placeholder. New ligand dose10nM,72h,2D CellTiter-Glo.", 
"Frozen40 molecular programs, two signatures, ERBB2/MET lapatinib maps and all single-drug outputs unchanged from v107.", 
"Shared symmetric g uses fixed signature scores times zero-referenced molecular RBF modifiers. No direct baseline functional predictor.", 
"Training: centered ALMANAC307 departures; within-recipient26-ligand MEMA contrasts;156 RAF/MEK plus16 HER2 exact departures.", 
"Four independent HGF RNA studies nearest24h; four alternative views share one physical-observation weight.", 
"One nonnegative ligand-single conversion learns actual no-drug arms, deduplicated across inhibitor/dose columns.", 
"g penalties0.003/0.03/0.3 use inner external whole-cell groups; outer excluded genes/cells stay excluded.", 
"External leave-cell-out, MEMA leave-gene-out, ALMANAC leave-cell-out; source-study and whole-compound validation separate.", 
"Absolute PRISM-to-MEMA baseline remains unvalidated. RNA single-ligand transport is a prediction, not observed combination RNA.", 
"Primary AU565 query5uM. Wilson1uM,12days is a different endpoint. Previously seen labels make all evaluations retrospective."
),file.path(root,"RUN_DETAILS.txt"))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
input_paths <- c('/work/src/program_pilot/133_counts_and_growth_order_bridge.R',
 '/work/src/program_pilot/133_assay_observation_fit.R',mutation_path,
 file.path(root,'data/HGF_input_reconstruction.qs'),
 file.path(auxiliary_root,'HGF_consensus_coordinate_definition.qs'),
 file.path(mema_root,'all_gene_endpoints.tsv'),file.path(root,'data/time_matched_mema_molecular_inputs.tsv'),
 file.path(parent,'models/final.qs'),archive_path)
write_table(data.table(path=input_paths,sha256=vapply(input_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(external_metrics)
print(metrics[evaluation=='held_gene'])
print(challenge)
print(hgf_query[hgf_context=='Study-balanced HGF (same recipe as training)'])
cat('Complete. Single-drug outputs byte-identical to v107; AU565+HGF and all Wilson outcomes excluded. HCC1954+HGF trains g.\n')
