#!/usr/bin/env Rscript
# Canonical HGF query uses the identical ligand recipe used in g training.
# Legacy mean-of-transported contexts is shown as a provenance diagnostic only.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(digest)})
root <- '/work/results/program_pilot/113_nonlinear_molecular_bypass_interaction/linear/data'
P <- sprintf('P%02d',1:40)
lig <- fread('/work/grant_rebuild_2026-10-04/transfer_solution/recipient_ligand_program_scores.tsv')
legacy <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
canonical <- as.numeric(lig[CellLine=='AU565' & gene=='HGF',..P])
old <- as.numeric(legacy[context=='Equal-context mean',..P])
comparison <- data.table(program=P,canonical_HGF=canonical,legacy_HGF=old,difference=canonical-old)
fwrite(comparison,file.path(root,'HGF_query_recipe_program_comparison.tsv'),sep='\t')
cat('Canonical versus legacy HGF cosine:',sum(canonical*old)/sqrt(sum(canonical^2)*sum(old^2)),
    '; norm ratio:',sqrt(sum(canonical^2)/sum(old^2)),'\n')
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
receptors <- as.data.table(a$prep$receptors)
versions <- c('108_mema_calibrated_shared_interaction','109_joint_conditional_interaction',
 '110_signature_anchored_program_interaction','111_external_bypass_combination_transfer',
 '112_external_signature_anchored_interaction','113_nonlinear_molecular_bypass_interaction',
 '114_braf_conditioned_bypass_transfer')
output <- list()
for(version in versions) {
  model <- qread(file.path('/work/results/program_pilot',version,'linear/models/final.qs'))
  ctx <- model$context
  recipe_rows <- list(canonical_same_as_training=canonical,legacy_mean_after_transport=old)
  lap <- model$rna$lap_receptors
  X <- as.matrix(receptors[depmap_id=='ACH-000248',.(ERBB2,MET)])
  Dlap <- sweep(sweep(sweep(sweep(X,2,lap$xmean),2,lap$xsd,'/') %*% lap$B,2,lap$sigma,'*'),2,lap$mu,'+')
  j <- 6L
  score <- ((Dlap-ctx$mu[j,,drop=FALSE])/ctx$sigma[j,,drop=FALSE]) %*% ctx$W
  route <- ctx$route[j]
  amp <- as.numeric(exp(ctx$Q[j,,drop=FALSE] %*% model$par$log_amplitude))
  S <- as.numeric(ctx$Q[j,,drop=FALSE] %*% model$par$level)+amp*score[1,route]
  L <- Dlap/ctx$sigma[j,,drop=FALSE]*amp
  for(recipe in names(recipe_rows)) {
    Dh <- matrix(recipe_rows[[recipe]],nrow=1)
    Qh <- cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
    scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
      crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
    hscale <- pmin(pmax(as.numeric(exp(Qh %*% scale_map)),apply(ctx$sigma[1:10,],2,min)),
      apply(ctx$sigma[1:10,],2,max))
    H <- sweep(Dh,2,hscale,'/')*as.numeric(exp(Qh %*% model$par$log_amplitude))
    SA <- L %*% ctx$W
    SB <- H %*% ctx$W
    forward <- reciprocal <- NA_real_
    if(length(model$par$Gamma)==3L) {
      features <- cbind(SA[,1]*SB[,1],SA[,1]*SB[,2]+SA[,2]*SB[,1],SA[,2]*SB[,2])
    } else if(is.null(model$interaction_basis)) {
      features <- do.call(cbind,lapply(1:2,function(k)H*SA[,k]+L*SB[,k]))
    } else {
      basis <- model$interaction_basis
      transform <- function(X) {
        distance <- pmax(outer(rowSums(X^2),rowSums(basis$centres^2),'+')-2*tcrossprod(X,basis$centres),0)
        sweep(exp(-distance/(2*basis$bandwidth^2)),2,
          exp(-rowSums(basis$centres^2)/(2*basis$bandwidth^2)),'-')
      }
      Fforward <- do.call(cbind,lapply(1:2,function(k)transform(H)*SA[,k]))
      Freciprocal <- do.call(cbind,lapply(1:2,function(k)transform(L)*SB[,k]))
      features <- Fforward+Freciprocal
      forward <- -as.numeric(sweep(Fforward,2,ctx$interaction_scale,'/') %*% model$par$Gamma)
      reciprocal <- -as.numeric(sweep(Freciprocal,2,ctx$interaction_scale,'/') %*% model$par$Gamma)
    }
    g <- -as.numeric(sweep(features,2,ctx$interaction_scale,'/') %*% model$par$Gamma)
    h <- as.numeric(-SB[1,1])
    output[[paste(version,recipe)]] <- data.table(version,recipe,dose=5,cell='AU565',
      lapatinib_pct_vehicle=100*2^(-S),HGF_increment_pct=100*2^h,
      interaction_log2_abundance=g,interaction_relative_pct=100*(2^g-1),
      directional_component=forward,reciprocal_component=reciprocal,
      combination_pct_vehicle=100*2^(-S+h+g),conditional_drug_pct_HGF_control=100*2^(-S+g),
      selection_rule='Use canonical recipe because it matches training; never selected by rescue magnitude.')
  }
}
prediction <- rbindlist(output)
fwrite(prediction,file.path(root,'canonical_and_legacy_HGF_query_predictions_frozen.tsv'),sep='\t')
fwrite(data.table(path=file.path(root,'canonical_and_legacy_HGF_query_predictions_frozen.tsv'),
 sha256=digest(file.path(root,'canonical_and_legacy_HGF_query_predictions_frozen.tsv'),file=TRUE,algo='sha256')),
 file.path(root,'canonical_HGF_query_freeze_hash.tsv'),sep='\t')
print(prediction[,.(version,recipe,HGF_increment_pct,interaction_relative_pct,combination_pct_vehicle)])
