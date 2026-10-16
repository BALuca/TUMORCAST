.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(digest)})
root <- '/work/results/program_pilot/91_lapatinib_crizotinib_5uM_joint_bridge'
pred <- fread(file.path(root,'data/single_predictions_frozen.tsv'))
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
P <- sprintf('P%02d',1:40)
families <- unique(pred$family)
checks <- list()
for(k in 1:3) for(family in families) {
  model <- qread(file.path(root,'models',paste0('recipient_',k,'_',gsub(' ','_',family),'.qs')))
  wanted_family <- family
  z <- pred[evaluation=='held_cell' & fold==k & family==wanted_family]
  stopifnot(!any(z$depmap_id %in% model$training_cells),!any(z$depmap_id %in% model$training_molecular_cells))
  baseline <- as.matrix(a$prep$raw[match(z$depmap_id,a$prep$raw$depmap_id),P])
  s <- sweep(sweep(sweep(baseline,2,model$context$baseline_mean),2,model$context$baseline_sd,'/')%*%
    model$context$baseline_projection,2,model$context$baseline_projection_scale,'/')
  x <- matrix(0,nrow(z),40)
  for(g in 1:2) {
    ii <- which(z$drug==model$selected_drugs[g])
    x[ii,] <- sweep(sweep(s[ii,,drop=FALSE],2,model$context$group_baseline_mean[g,])%*%
      model$molecular_coefficient[,,g],2,model$group_response_mean[g,],'+')
  }
  latent <- tanh(x%*%model$decoder$u)
  pairs <- which(upper.tri(matrix(0,4,4),diag=TRUE),arr.ind=TRUE)
  h <- cbind(tanh(x),x^2/(1+x^2),latent[,pairs[,1],drop=FALSE]*latent[,pairs[,2],drop=FALSE])
  replay <- as.numeric(x%*%model$decoder$w+h%*%model$decoder$a)
  error <- max(abs(replay-z$predicted_log2))
  stopifnot(error<1e-10,model$combination_training_rows==0L,!model$functional_baseline_inputs)
  checks[[paste(k,family)]] <- data.table(fold=k,family,maximum_replay_difference=error)
}
fwrite(rbindlist(checks),file.path(root,'data/prediction_replay_validation.tsv'),sep='\t')
old <- fread('/work/results/program_pilot/75_hgf_pattern_partner_transfer/data/protected_441_monotherapy_predictions.tsv')
stopifnot(uniqueN(old$depmap_id)==441L)
common <- pred[evaluation=='held_cell' & drug=='Lapatinib ditosylate' & depmap_id %in% old$depmap_id]
metrics <- common[,.(cells=.N,pearson_log2=cor(observed_log2,predicted_log2)),by=family]
stopifnot(all(metrics$cells==441L))
fwrite(metrics,file.path(root,'data/original_441_cells_at_5uM.tsv'),sep='\t')
paths <- file.path(root,'data',c('single_predictions_frozen.tsv','measured_displacement_predictions_frozen.tsv'))
fwrite(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),
  file.path(root,'data/prediction_hashes.tsv'),sep='\t')
print(metrics)
