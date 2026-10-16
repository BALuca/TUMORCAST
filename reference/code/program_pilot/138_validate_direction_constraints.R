#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/138_pooled_direction_constrained_rescue/linear'
old <- qread('/work/results/program_pilot/136_programwise_reversal_reinforcement/linear/models/final.qs')
paths <- list.files(file.path(root,'models'),pattern='^(final|held_cell_.+|held_ligand_.+)\\.qs$',full.names=TRUE)
checks <- rbindlist(lapply(paths,function(path) {
  m <- qread(path)
  stopifnot(length(m$par$Gamma)==2, m$par$Gamma[1]<=0,m$par$Gamma[2]>=0,
    identical(m$rna,old$rna),identical(m$template,old$template),
    identical(m$par$level,old$par$level),identical(m$par$log_amplitude,old$par$log_amplitude))
  data.table(model=basename(path),rescue_strength=-m$par$Gamma[1],
    reinforcement_strength=m$par$Gamma[2],kkt_residual=m$interaction_training$kkt_residual)
}))
stopifnot(nrow(checks)==20)
fwrite(checks,file.path(root,'data/all_fit_direction_validation.tsv'),sep='\t')
for(e in parse('/work/src/program_pilot/138_pooled_direction_constrained_rescue.R')) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('pair_features'))) eval(e)
}
m <- qread(file.path(root,'models/final.qs'))
f <- qread(file.path(root,'models/final_molecular_features.qs'))
X <- pair_features(f$Aviews,f$Bviews,m$context$W,m$interaction_basis)
stopifnot(identical(X,f$products),max(abs(X-pair_features(f$Bviews,f$Aviews,m$context$W,m$interaction_basis)))<1e-12,
  all(pair_features(f$Aviews,f$Bviews*0,m$context$W,m$interaction_basis)==0))
# Show that opposite signature orientations cannot produce opposing reversal terms.
W <- matrix(0,40,2); W[1,] <- c(1,-1)
A <- B <- matrix(0,1,42); A[1,1] <- 2; B[1,1] <- -3
stopifnot(identical(as.numeric(pair_features(A,B,W,NULL)),c(6,0)))
fwrite(data.table(check=c('all_20_fits_obey_direction_constraints','single_drug_components_unchanged',
  'symmetry_on_real_inputs','zero_when_partner_absent','opposite_signatures_pool_without_cancellation'),passed=TRUE),
 file.path(root,'data/structural_validation.tsv'),sep='\t')
cat('Validated all 20 fits; molecular symmetry/zero checks and shared direction constraints pass.\n')
