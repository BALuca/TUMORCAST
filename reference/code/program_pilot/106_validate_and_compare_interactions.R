.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
new <- '/work/results/program_pilot/106_lapatinib_supported_dose_interactions/linear'
old <- '/work/results/program_pilot/105_erbb2_met_without_crizotinib_vismodegib/linear'
for(label in c(sprintf('loo_%02d',1:62),'final')) {
 a <- qread(file.path(old,'models',paste0(label,'.qs')))
 b <- qread(file.path(new,'models',paste0(label,'.qs')))
 stopifnot(identical(a$template,b$template),identical(a$rna$lap_receptors,b$rna$lap_receptors),
           identical(a$rna$B,b$rna$B),
           isTRUE(all.equal(as.data.frame(a$selected),as.data.frame(b$selected),check.attributes=FALSE)),
           !any(b$held_cells %in% unlist(lapply(b$rna$extra_rna,`[[`,'training_cells'))))
}
audit <- fread(file.path(new,'data/fit_audit.tsv'))
stopifnot(nrow(audit)==63L,all(audit$convergence==0L),all(audit$functional_leakage==0L),all(audit$RNA_leakage==0L))
check <- fread(file.path(new,'data/analytic_gradient_check.tsv'))
stopifnot(max(check$relative_difference)<1e-5)
pair <- fread(file.path(new,'data/combination_training_outcomes.tsv'))
stopifnot(nrow(pair)==307L,uniqueN(pair$pair)==8L,!any(pair$depmap_id=='ACH-000248'),
          !any(c(pair$drugA,pair$drugB) %in% c('crizotinib','Vismodegib')),
          !anyDuplicated(pair[,.(depmap_id,pair,SCREENER,doseA,doseB)]))
q <- fread(file.path(new,'data/au565_hgf_5uM_zero_reference_query.tsv'))
q[,hgf_addition_relative_pct:=100*(combination_pct_zero_control/lapatinib_pct_zero_control-1)]
q[,interaction_relative_pct:=100*(2^interaction_log2-1)]
fwrite(q,file.path(new,'data/au565_hgf_5uM_zero_reference_summary.tsv'),sep='\t')
print(q[hgf_context=='Equal-context mean',.(assumed_HGF_route,lapatinib_pct_zero_control,
     hgf_alone_pct_zero_control,combination_pct_zero_control,hgf_addition_relative_pct,interaction_relative_pct)])
comparison <- rbindlist(lapply(c(old,new),function(root) {
  d <- fread(file.path(root,'data/baseline_predictions_frozen.tsv'))[evaluation=='held_cell' & drug=='Lapatinib ditosylate']
  m <- fread(file.path(root,'data/measured_predictions_frozen.tsv'))[evaluation=='held_cell' & drug=='Lapatinib ditosylate']
  data.table(version=if(root==old)105L else 106L,endpoint='Lapatinib 5 uM log2 abundance across held cells',
    all_cells=nrow(d),all_r=cor(d$observed_log2,d$predicted_log2),matched_cells=nrow(m),
    matched_predicted_bridge_r=cor(d[match(m$depmap_id,depmap_id),observed_log2],d[match(m$depmap_id,depmap_id),predicted_log2]),
    matched_measured_decoder_r=cor(m$observed_log2,m$predicted_log2),
    observed_sensitive=sum(d$observed_log2<= -1),predicted_sensitive=sum(d$predicted_log2<= -1))
}))
fwrite(comparison,file.path(new,'data/lapatinib_comparison_with_v105.tsv'),sep='\t')
print(comparison)
print(fread(file.path(new,'data/combination_correlations.tsv'))[evaluation=='held_cell'])
cat('Verified: all 63 fits preserve original signatures, selected doses and original molecular predictors; no cell leakage; 307 valid pair observations.\n')
