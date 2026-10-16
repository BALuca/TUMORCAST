.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/92_shared_molecular_pattern_sensitivity_decoder'
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
selected <- fread(file.path(root,'data/selected_doses.tsv'))
P <- sprintf('P%02d',1:40)
m <- merge(as.data.table(a$prep$m),selected,by=c('drug','dose'))
b <- as.data.table(a$prep$raw)
single <- fread(file.path(root,'data/single_predictions_frozen.tsv'))
pair_pred <- fread(file.path(root,'data/interaction_predictions_frozen.tsv'))
pairs <- fread(file.path(root,'data/combination_outcomes.tsv'))
key <- function(cell,drug)paste(cell,drug)
# Independently replay the serialized representation, membership weights and rules.
contributions <- function(model,D) {
  Z <- D%*%model$rep$transform
  common <- as.numeric(Z%*%model$beta[1:6])
  if(model$K==1) return(matrix(common,ncol=1))
  d <- vapply(1:model$K,function(k)rowSums(sweep(Z,2,model$rep$gates[[model$K]]$centers[k,])^2),numeric(nrow(D)))
  logg <- -d/(2*model$rep$gates[[model$K]]$radius)
  logg <- logg-apply(logg,1,max)
  g <- exp(logg)/rowSums(exp(logg))
  stopifnot(max(abs(rowSums(g)-1))<1e-12)
  vapply(1:model$K,function(k)g[,k]*(common+as.numeric(Z%*%model$beta[6+((k-1)*6+1:6)])),numeric(nrow(D)))
}
query <- function(model,cell,drug) {
  ix <- match(key(cell,drug),key(m$depmap_id,m$drug))
  D <- matrix(0,length(cell),40)
  have <- which(!is.na(ix))
  D[have,] <- as.matrix(m[ix[have],..P])
  miss <- which(is.na(ix))
  if(length(miss)) {
    raw <- as.matrix(b[match(cell[miss],depmap_id),..P])
    z <- sweep(sweep(sweep(raw,2,model$rep$bmean),2,model$rep$bsd,'/')%*%model$rep$BV,2,model$rep$bscale,'/')
    for(dr in unique(drug[miss])) {
      ii <- which(drug[miss]==dr)
      op <- model$rep$operators[[dr]]
      D[miss[ii],] <- sweep(sweep(z[ii,,drop=FALSE],2,op$xmean)%*%op$beta,2,op$ymean,'+')
    }
  }
  D
}
pair_score <- function(model,A,B) {
  ca <- tanh(contributions(model,A));cb <- tanh(contributions(model,B))
  ij <- which(upper.tri(matrix(0,model$K,model$K),diag=TRUE),arr.ind=TRUE)
  X <- matrix(0,nrow(A),nrow(ij))
  for(j in 1:nrow(ij)) {
    p <- ij[j,1];q <- ij[j,2]
    X[,j] <- if(p==q)ca[,p]*cb[,p] else (ca[,p]*cb[,q]+ca[,q]*cb[,p])/sqrt(2)
  }
  as.numeric(X%*%model$pair_beta[1:ncol(X)])
}
checks <- list()
for(fold in 1:3) for(K in 1:3) {
  model <- qread(file.path(root,'models',paste0('fold_',fold,'_K',K,'.qs')))
  # Use explicit vectors because data.table column names also name these loop variables.
  ii <- which(single$evaluation=='held_cell' & single$fold==fold & single$K==K)
  z <- single[ii]
  pi <- which(pair_pred$evaluation=='held_cell' & pair_pred$fold==fold & pair_pred$K==K)
  p <- pair_pred[pi]
  stopifnot(!any(z$depmap_id %in% model$training_cells),!any(z$depmap_id %in% model$rep$training_cells),
    !any(p$depmap_id %in% model$training_pair_cells),!any(p$depmap_id %in% model$rep$training_cells))
  D <- query(model,z$depmap_id,z$drug)
  ds <- max(abs(rowSums(contributions(model,D))-z$predicted_score))
  meta <- pairs[match(p$row_id,row_id)]
  A <- query(model,meta$depmap_id,meta$drugA)
  B <- query(model,meta$depmap_id,meta$drugB)
  replay <- pair_score(model,A,B)
  dp <- max(abs(replay-p$predicted_departure_score))
  symmetric <- max(abs(replay-pair_score(model,B,A)))
  zero <- max(abs(pair_score(model,A*0,B)))
  stopifnot(ds<1e-10,dp<1e-10,symmetric<1e-12,zero<1e-12)
  checks[[paste(fold,K)]] <- data.table(fold,K,single_replay_difference=ds,pair_replay_difference=dp,
    pair_symmetry_difference=symmetric,pair_zero_input_difference=zero)
}
fwrite(rbindlist(checks),file.path(root,'data/prediction_validation.tsv'),sep='\t')
# Assay-scale diagnostics expose near-flat predictions even when r is positive.
choice <- fread(file.path(root,'data/model_selection.tsv'))[selected & job!='all_data',.(fold=as.integer(sub('fold_','',job)),K)]
z <- merge(pair_pred[evaluation=='held_cell'],choice,by=c('fold','K'))
z[,`:=`(yc=observed_departure_log2-mean(observed_departure_log2),pc=predicted_departure_score-mean(predicted_departure_score)),by=.(pair,SCREENER,fold)]
fwrite(z[,.(cells=.N,observed_departure_SD_log2=sd(yc),predicted_departure_SD_log2=sd(pc),both_RNA_inputs_measured=sum(both_inputs_measured)),by=.(pair,SCREENER)],
  file.path(root,'data/interaction_assay_scale.tsv'),sep='\t')
cat('Validated all nine held-cell models: saved predictions replay, excluded cells, symmetric pair rule and zero-input behavior.\n')
