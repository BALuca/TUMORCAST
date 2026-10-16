#!/usr/bin/env Rscript
# Inspect the frozen fit's training geometry; no new predictor or parameter selection.
prefix <- readLines('/work/src/program_pilot/136_validate_and_compare.R')
cut <- which(grepl('^mf <-', prefix))
eval(parse(text = prefix[seq_len(cut - 1L)]))
frozen <- original_final
fit_diagnostics <- fit_interaction
body(fit_diagnostics)[[length(body(fit_diagnostics))]] <- quote(list(
  model=model, ng=ng, mg=mg, eg=eg, Xn=Xn, Xm=Xm, Xec=Xec, Xem=Xem,
  scale=assay_scale, penalty=penalty))
d <- fit_diagnostics(frozen, mema_features(frozen), primary_endpoint,
  penalty=frozen$interaction_training$interaction_penalty)
print(cbind(reconstructed=d$model$par$Gamma,frozen=frozen$par$Gamma)); stopifnot(max(abs(d$model$par$Gamma-frozen$par$Gamma))<1e-6)
terms <- frozen$interaction_basis$labels
blocks <- list(ALMANAC=list(X=d$Xn,y=-d$ng$yc,w=d$ng$w),
  MEMA=list(X=d$Xm*d$scale,y=-d$mg$yc,w=d$mg$w),
  External_within=list(X=d$Xec,y=-d$eg$yc,w=d$eg$w),
  External_means=list(X=d$Xem,y=-d$eg$ym,w=d$eg$w))
associations <- rbindlist(lapply(names(blocks),function(n) {
  b <- blocks[[n]]
  data.table(source=n,term=terms,weighted_r=vapply(seq_along(terms),function(j)
    cov.wt(cbind(b$X[,j],b$y),wt=b$w,cor=TRUE)$cor[1,2],numeric(1)),
    feature_SD=sqrt(colSums(sweep(b$X,2,colSums(b$X*b$w),'-')^2*b$w)))
}))
write_table(associations,'g_training_feature_associations')
print(associations)
correlations <- rbindlist(lapply(names(blocks),function(n) {
  b <- blocks[[n]]
  C <- cov.wt(b$X,wt=b$w,cor=TRUE)$cor
  data.table(source=n,term=rep(terms,each=4),other=rep(terms,4),r=as.numeric(t(C)))
}))
write_table(correlations,'g_training_feature_collinearity')
print(correlations[term=='reversal_1' & other=='reversal_2'])
print(frozen$interaction_training[c('interaction_penalty','external_HGF_rows','ligand_single_conversion')])
print(fread(file.path(root,'data/external_held_cell_predictions_frozen.tsv'))[gene=='HGF',
 .(observed_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2)),by=drug])
