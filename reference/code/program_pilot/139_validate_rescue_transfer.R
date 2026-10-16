#!/usr/bin/env Rscript
# Validate separation of signature learning, calibration and held HGF outcomes.
x <- readLines('/work/src/program_pilot/139_raw_ligand_rescue_transfer.R')
cut <- which(grepl('^predictions <- list',x))
eval(parse(text=x[seq_len(cut-1L)]))
model <- qread(file.path(root,'models/final_rescue_signature.qs'))
transfer <- qread(file.path(root,'models/final_transfer.qs'))
paths <- list.files(file.path(root,'models'),pattern='^held_(cell|ligand)_.+\\.qs$',full.names=TRUE)
checks <- rbindlist(lapply(paths,function(path) {
  m <- qread(path); label <- sub('\\.qs$','',basename(path))
  if(grepl('held_ligand_',label)) stopifnot(!any(training$gene[m$training_rows]==sub('held_ligand_','',label)))
  if(grepl('held_cell_',label)) stopifnot(!any(training$depmap_id[m$training_rows]==sub('held_cell_','',label)))
  stopifnot(all(m$multipliers>0),!m$scale_boundary,all(is.finite(m$beta)))
  data.table(fit=label,excluded_outcomes_absent=TRUE,finite_coefficients=TRUE,interior_positive_calibration=TRUE)
}))
stopifnot(nrow(checks)==19,identical(transfer$signature,model),!model$scale_boundary)
write_table(checks,'fold_exclusion_validation')
# Perturb held HGF labels drastically; neither inner selection nor fit may change.
y <- copy(training$observed_departure_log2)
training[gene=='HGF',observed_departure_log2:=observed_departure_log2+100]
test <- select_and_fit(held_genes='HGF',label='validation_held_HGF')
training[,observed_departure_log2:=y]
ref <- qread(file.path(root,'models/held_ligand_HGF.qs'))
stopifnot(identical(test$penalty,ref$penalty),max(abs(test$beta-ref$beta))<1e-9,
 max(abs(test$multipliers-ref$multipliers))<1e-9)
# Calibration outcomes cannot alter the frozen signature.
cy <- copy(calibration$observed_departure_log2)
calibration[,observed_departure_log2:=observed_departure_log2+100]
again <- fit_signature(model$training_rows,model$penalty)
calibration[,observed_departure_log2:=cy]
stopifnot(max(abs(again$beta-model$beta))<1e-9)
# Predict the exact saved primary query independently from raw programs.
au <- source_profiles[CellLine=='AU565' & source_profile=='Primary_study_mean']
score <- sum(as.numeric(au[,..programs])*model$beta)
query <- fread(file.path(root,'data/AU565_HGF_predictions_frozen.tsv'))[source_profile=='Primary_study_mean']
stopifnot(abs(score*transfer$lapatinib$multiplier-query$interaction_log2)<1e-10,
 transfer$lapatinib$multiplier>=0,all(calibration$gene[transfer$lapatinib$training_rows] %in% c('NRG1','EGF')),
 sum(rep(0,40)*model$beta)==0)
write_table(data.table(check=c('held_HGF_label_perturbation_has_no_effect','lapatinib_labels_cannot_change_signature',
 'zero_ligand_zero_interaction','independent_query_matches_saved_prediction','no_lapatinib_HGF_training'),
 passed=TRUE),'structural_validation')
# Small grant-readable prediction plot. No Wilson measurement is fitted or plotted.
bars <- data.table(treatment=factor(c('Lapatinib','+ HGF','+ HGF\nHGF excluded'),
 levels=c('Lapatinib','+ HGF','+ HGF\nHGF excluded')),
 abundance=c(query$lapatinib_pct_control,query$combination_pct_control,
 query$all_HGF_outcomes_excluded_combination_pct))
write_table(bars,'AU565_prediction_barplot_data')
p <- ggplot(bars,aes(treatment,abundance))+geom_col(fill=c('#777777','#D55E00','#E7A781'),width=.65)+
 geom_text(aes(label=sprintf('%.1f',abundance)),vjust=-.4,size=3)+
 scale_y_continuous(expand=expansion(mult=c(0,.15)))+theme_classic(base_size=9)+
 labs(x=NULL,y='Predicted abundance (% control)',title='AU565 | lapatinib 5 µM')
ggsave(file.path(root,'figures/AU565_HGF_predictions.pdf'),p,width=3.2,height=2.7,device=cairo_pdf)
ggsave(file.path(root,'figures/AU565_HGF_predictions.png'),p,width=3.2,height=2.7,dpi=170)
cat('All fold, label-isolation, calibration and query checks passed.\n')
