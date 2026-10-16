.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/91_lapatinib_crizotinib_5uM_joint_bridge'
a <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
P <- sprintf('P%02d',1:40)
pairs <- which(upper.tri(matrix(0,4,4),diag=TRUE),arr.ind=TRUE)
rows <- list()
for(family in c('Frozen molecular bridge','Linked molecular bridge')) {
  model <- qread(file.path(root,'models',paste0('final_',gsub(' ','_',family),'.qs')))
  stopifnot(!'ACH-000248' %in% model$training_cells,!'ACH-000248' %in% model$training_molecular_cells,
    max(abs(model$decoder$interaction))==0)
  b <- as.matrix(a$prep$raw[match('ACH-000248',a$prep$raw$depmap_id),P])
  s <- sweep(sweep(sweep(b,2,model$context$baseline_mean),2,model$context$baseline_sd,'/')%*%
    model$context$baseline_projection,2,model$context$baseline_projection_scale,'/')
  g <- match('Lapatinib ditosylate',model$selected_drugs)
  lap <- sweep(sweep(s,2,model$context$group_baseline_mean[g,])%*%model$molecular_coefficient[,,g],
    2,model$group_response_mean[g,],'+')
  h <- sweep(as.matrix(hgf[,..P]),2,model$context$response_sd,'/')
  functional <- function(x) {
    latent <- tanh(x%*%model$decoder$u)
    features <- cbind(tanh(x),x^2/(1+x^2),latent[,pairs[,1],drop=FALSE]*latent[,pairs[,2],drop=FALSE])
    as.numeric(x%*%model$decoder$w+features%*%model$decoder$a)
  }
  fl <- functional(lap)
  fh <- functional(h)
  rows[[family]] <- data.table(family,dose=5,cell='AU565',hgf_context=hgf$context,
    lapatinib_pct_vehicle=100*2^fl,hgf_alone_extrapolated_pct_vehicle=100*2^fh,
    additive_combination_extrapolated_pct_vehicle=100*2^(fl+fh),
    lapatinib_in_hgf_pct_hgf_control=100*2^fl,interaction_log2=0,
    qualification='No combination supervision; interaction is untrained zero, not evidence of absent rescue')
}
z <- rbindlist(rows)
fwrite(z,file.path(root,'data/au565_5uM_hgf_query.tsv'),sep='\t')
print(z[hgf_context=='Equal-context mean',.(family,lapatinib_pct_vehicle,hgf_alone_extrapolated_pct_vehicle,
  additive_combination_extrapolated_pct_vehicle,lapatinib_in_hgf_pct_hgf_control)])
