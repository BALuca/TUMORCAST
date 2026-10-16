#!/usr/bin/env Rscript
# Descriptive profile comparison only; this does not select rows or refit g.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
out <- file.path(root,'data/single_signatures_beside_combination_benefit')
m <- qread(file.path(root,'models/final.qs'))
p <- fread('/work/results/program_pilot/134_HGF_g_diagnostics/linear/data/AU565_program_level_opposition.tsv')[source_profile=='Primary_study_mean']
programs <- sprintf('P%02d',1:40)
p <- p[match(programs,program)]
D <- matrix(p$HGF_displacement,nrow=1)
ctx <- m$context
Q <- cbind(1,sweep(sweep(sweep(D,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
scale <- pmin(pmax(exp(Q %*% map),matrix(apply(ctx$sigma[1:10,],2,min),1)),matrix(apply(ctx$sigma[1:10,],2,max),1))
H <- as.numeric(D/scale)
corr <- fread(file.path(out,'program_combination_benefit_correlations.tsv'))[dataset!='HGF_focus']
ranked <- corr[,{
  ix <- match(program,programs)
  list(n_cells_or_ligands=first(n),profile_match_scaled_HGF=if(all(is.finite(r)))cor(r,H[ix]) else NA_real_,
    profile_match_raw_HGF=if(all(is.finite(r)))cor(r,p$HGF_displacement[ix]) else NA_real_,
    profile_match_opposite_lapatinib_sensitivity=if(all(is.finite(r)))cor(r,-m$template$correlations[program,'Lapatinib ditosylate']) else NA_real_)
},by=.(dataset,partner,context)]
setorder(ranked,-profile_match_scaled_HGF,na.last=TRUE)
fwrite(ranked,file.path(out,'HGF_profile_match_by_training_context.tsv'),sep='\t')
summary <- ranked[,.(contexts=.N,median_HGF_match=median(profile_match_scaled_HGF,na.rm=TRUE),
  minimum_HGF_match=min(profile_match_scaled_HGF,na.rm=TRUE),maximum_HGF_match=max(profile_match_scaled_HGF,na.rm=TRUE)),by=.(dataset,partner)]
fwrite(summary,file.path(out,'HGF_profile_match_by_source.tsv'),sep='\t')
print(summary)
cat('\nHGF training contexts:\n')
print(ranked[grepl('^HGF',context)])
cat('\nHighest ligand matches:\n')
print(head(ranked[partner=='B' & dataset!='ALMANAC'],12))
cat('\nHighest drug matches:\n')
print(head(ranked[partner=='A' | dataset=='ALMANAC'],8))
# Check whether the across-ligand pattern survives removing HGF itself.
# Reuse molecular reconstruction only; do not run the historical plotting workflow.
env <- new.env(parent=.GlobalEnv)
x <- readLines('/work/src/program_pilot/135_plot_combination_departure_signatures.R')
x <- gsub('envir = .GlobalEnv','envir = env',x,fixed=TRUE)
cut <- which(grepl('^# Each physical observation',x))
eval(parse(text=x[seq_len(cut-1L)]),envir=env)
rows <- which(env$external$study=='Straussman2012' & env$external$gene!='HGF')
meta <- copy(env$external[rows])
X <- env$EB[rows,,drop=FALSE]
meta[,hgf_direction_score:=as.numeric(X %*% (H/sqrt(sum(H^2))))]
checks <- meta[,{
  ii <- .I
  rr <- apply(X[ii,,drop=FALSE],2,function(z)cor(z,outcome))
  list(ligands=.N,program_profile_match_HGF=cor(rr,H),
    r_HGF_direction_vs_rescue=cor(hgf_direction_score,outcome))
},by=.(CellLine,drug,dose)]
fwrite(checks,file.path(out,'across_ligand_signal_excluding_HGF.tsv'),sep='\t')
print(checks)
print(checks[,.(median_profile_match=median(program_profile_match_HGF),
  range_min=min(program_profile_match_HGF),range_max=max(program_profile_match_HGF),
  positive_contexts=sum(program_profile_match_HGF>0),contexts=.N,
  median_score_r=median(r_HGF_direction_vs_rescue))])
