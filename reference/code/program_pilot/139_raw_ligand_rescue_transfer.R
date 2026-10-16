#!/usr/bin/env Rscript
# Shared raw-ligand rescue signature, then frozen-signature lapatinib calibration.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest)})
RhpcBLASctl::blas_set_num_threads(1L)
setDTthreads(1L)
root <- '/work/results/program_pilot/139_raw_ligand_rescue_transfer/linear'
source_root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
programs <- sprintf('P%02d',1:40)
within_weight <- .8
penalty_grid <- c(.001,.01,.1,1,10)
for(folder in c('data','models','figures'))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
pearson <- function(x,y)if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10)cor(x,y) else NA_real_
all_outcomes <- fread(file.path(source_root,'data/external_combination_outcomes.tsv'))
ligands <- fread(file.path(source_root,'data/external_ligand_programs.tsv'))
source_profiles <- fread(file.path(source_root,'data/time_matched_HGF_program_inputs.tsv'))
training <- all_outcomes[study=='Straussman2012']
calibration <- all_outcomes[study=='elife32271']
stopifnot(nrow(training)==156,uniqueN(training$gene)==13,uniqueN(training$depmap_id)==6,
 nrow(calibration)==16,all(calibration$gene %in% c('NRG1','EGF')),
 all(calibration$drug=='Lapatinib ditosylate'),!any(calibration$gene=='HGF'))
key <- function(d)paste(d$depmap_id,d$gene,sep='|')
stopifnot(!anyDuplicated(key(ligands)))
ix <- match(key(training),key(ligands)); ci <- match(key(calibration),key(ligands))
stopifnot(!anyNA(ix),!anyNA(ci))
X <- as.matrix(ligands[ix,..programs]); CX <- as.matrix(ligands[ci,..programs])
stopifnot(all(is.finite(X)),all(is.finite(CX)))
training[,observation_id:=.I]
calibration[,observation_id:=.I]
write_table(training,'RAF_MEK_training_outcomes')
write_table(calibration,'lapatinib_calibration_outcomes')
write_table(cbind(training[,.(observation_id,CellLine,drug,dose,gene)],as.data.table(X)),'raw_ligand_inputs')
write_table(cbind(calibration[,.(observation_id,CellLine,drug,dose,gene)],as.data.table(CX)),'lapatinib_raw_ligand_inputs')

# Weighted contrasts across ligands within fixed cell/drug. Only training outcomes
# enter centering and fitting; uncentered predictions preserve zero ligand => zero g.
prepare_geometry <- function(meta,design) {
  group <- paste(meta$depmap_id,meta$drug,meta$dose,sep='|')
  w <- 1/as.numeric(table(group)[group])/length(unique(group))
  XM <- matrix(0,nrow(design),ncol(design)); ym <- numeric(nrow(design))
  y <- meta$observed_departure_log2
  for(g in unique(group)) {
    ii <- which(group==g)
    XM[ii,] <- matrix(colMeans(design[ii,,drop=FALSE]),length(ii),ncol(design),byrow=TRUE)
    ym[ii] <- mean(y[ii])
  }
  xc <- design-XM; yc <- y-ym
  vc <- sum(w*yc^2); vm <- sum(w*ym^2)
  stopifnot(vc>1e-10,vm>1e-10)
  parts <- lapply(c('PD-184352','PLX-4720'),function(dr) {
    i <- which(meta$drug==dr)
    list(G=within_weight*crossprod(xc[i,,drop=FALSE],xc[i,,drop=FALSE]*w[i])/vc+
      (1-within_weight)*crossprod(XM[i,,drop=FALSE],XM[i,,drop=FALSE]*w[i])/vm,
      b=within_weight*as.numeric(crossprod(xc[i,,drop=FALSE],w[i]*yc[i]))/vc+
      (1-within_weight)*as.numeric(crossprod(XM[i,,drop=FALSE],w[i]*ym[i]))/vm)
  })
  names(parts) <- c('PD-184352','PLX-4720')
  list(parts=parts,penalty_unit=sum(diag(parts[[1]]$G+parts[[2]]$G))/40)
}

# MEK multiplier is fixed to 1, resolving the beta/multiplier scale ambiguity.
# A single positive RAF multiplier is profiled; the raw program vector is never
# divided by a ligand-specific or program-specific SD. Ridge is isotropic in raw units.
fit_signature <- function(rows,penalty) {
  prep <- prepare_geometry(training[rows],X[rows,,drop=FALSE])
  one <- prep$parts[[1]]; two <- prep$parts[[2]]
  at_scale <- function(log_a) {
    a <- exp(log_a)
    # Penalize the effective program weights for BOTH drugs equally. Penalizing
    # beta alone would let an expanding RAF multiplier evade regularization.
    ridge <- penalty*prep$penalty_unit*(1+a^2)/2
    G <- one$G+a^2*two$G+diag(ridge,40)
    b <- one$b+a*two$b
    beta <- as.numeric(solve(G,b))
    list(objective=-sum(beta*b),beta=beta,normal_residual=max(abs(G %*% beta-b)))
  }
  # Coarse numerical bracketing prevents assuming the profiled curve is unimodal.
  brackets <- seq(-4,4,by=.5)
  values <- vapply(brackets,function(a)at_scale(a)$objective,numeric(1))
  j <- which.min(values)
  interval <- brackets[c(max(1,j-1),min(length(brackets),j+1))]
  opt <- optimize(function(a)at_scale(a)$objective,interval=interval,tol=1e-7)
  solution <- at_scale(opt$minimum)
  stopifnot(solution$normal_residual<1e-7)
  list(beta=setNames(solution$beta,programs),multipliers=c('PD-184352'=1,'PLX-4720'=exp(opt$minimum)),
    penalty=penalty,penalty_unit=prep$penalty_unit,training_rows=rows,
    scale_boundary=abs(opt$minimum)>3.99,normal_residual=solution$normal_residual,
    description='Raw signed ligand programs; shared ridge beta; positive relative RAF/MEK calibration; no intercept')
}
predict_signature <- function(model,design,drug) {
  multiplier <- model$multipliers[drug]
  stopifnot(!anyNA(multiplier))
  as.numeric(design %*% model$beta)*as.numeric(multiplier)
}
selection_audit <- list()
select_and_fit <- function(held_cells=character(),held_genes=character(),label='final') {
  rows <- which(!training$depmap_id %in% held_cells & !training$gene %in% held_genes)
  genes <- sort(unique(training$gene[rows]))
  fold <- setNames(rep(1:3,length.out=length(genes)),genes)
  pred <- matrix(NA_real_,length(rows),length(penalty_grid))
  for(k in 1:3) {
    valid <- which(fold[training$gene[rows]]==k)
    train_rows <- rows[-valid]
    for(j in seq_along(penalty_grid)) {
      model <- fit_signature(train_rows,penalty_grid[j])
      pred[valid,j] <- predict_signature(model,X[rows[valid],,drop=FALSE],training$drug[rows[valid]])
    }
  }
  stopifnot(all(is.finite(pred)))
  # Evaluate the same 80/20 contrast objective on inner held-out ligands, grouped
  # within each recipient/inhibitor; validation labels cannot enter the fitted beta.
  group <- paste(training$depmap_id[rows],training$drug[rows],sep='|')
  yc <- training$observed_departure_log2[rows]
  ym <- ave(yc,group,FUN=mean); yc <- yc-ym
  pm <- apply(pred,2,function(z)ave(z,group,FUN=mean))
  pc <- pred-pm
  w <- 1/as.numeric(table(group)[group])/length(unique(group))
  loss <- within_weight*colSums((pc-yc)^2*w)/sum(w*yc^2)+
    (1-within_weight)*colSums((pm-ym)^2*w)/sum(w*ym^2)
  chosen <- which.min(loss)
  selection_audit[[label]] <<- data.table(fit=label,penalty=penalty_grid,
    internal_selection_loss=loss,selected=seq_along(loss)==chosen,
    excluded_cells=paste(held_cells,collapse=','),excluded_ligands=paste(held_genes,collapse=','))
  fit_signature(rows,penalty_grid[chosen])
}

predictions <- list(); coefficients <- list(); fitted_models <- list()
for(gene in sort(unique(training$gene))) {
  label <- paste0('held_ligand_',gene)
  model <- select_and_fit(held_genes=gene,label=label)
  rows <- which(training$gene==gene)
  stopifnot(!any(rows %in% model$training_rows))
  predictions[[label]] <- training[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=predict_signature(model,X[rows,,drop=FALSE],training$drug[rows]),evaluation='held_ligand')]
  fitted_models[[label]] <- model
  coefficients[[label]] <- data.table(fit=label,program=programs,beta=model$beta)
  qsave(model,file.path(root,'models',paste0(label,'.qs')))
  cat('Ligand complete:',gene,'\n')
}
for(cell in sort(unique(training$depmap_id))) {
  label <- paste0('held_cell_',cell)
  model <- select_and_fit(held_cells=cell,label=label)
  rows <- which(training$depmap_id==cell)
  stopifnot(!any(rows %in% model$training_rows))
  predictions[[label]] <- training[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=predict_signature(model,X[rows,,drop=FALSE],training$drug[rows]),evaluation='held_cell')]
  fitted_models[[label]] <- model
  coefficients[[label]] <- data.table(fit=label,program=programs,beta=model$beta)
  qsave(model,file.path(root,'models',paste0(label,'.qs')))
  cat('Cell complete:',cell,'\n')
}
model <- select_and_fit()
fitted_models[['final']] <- model
coefficients[['final']] <- data.table(fit='final',program=programs,beta=model$beta)
qsave(model,file.path(root,'models/final_rescue_signature.qs'))
write_table(rbindlist(coefficients),'program_coefficients_all_fits')
write_table(rbindlist(selection_audit),'training_only_penalty_selection')
write_table(rbindlist(lapply(names(fitted_models),function(n) {
  m <- fitted_models[[n]]
  data.table(fit=n,RAF_multiplier=m$multipliers[['PLX-4720']],scale_boundary=m$scale_boundary,
    penalty=m$penalty,normal_residual=m$normal_residual)
})),'relative_drug_calibration_all_fits')

# Freeze beta before accessing the lapatinib calibration outcomes. This one
# nonnegative coefficient pools both native doses; no dose-specific/HGF-specific fit.
calibration_score <- as.numeric(CX %*% model$beta)
calibrate <- function(rows,score=calibration_score) {
  group <- paste(calibration$gene[rows],calibration$dose[rows],sep='|')
  w <- 1/as.numeric(table(group)[group])/length(unique(group))
  y <- calibration$observed_departure_log2[rows]
  denominator <- sum(w*score[rows]^2)
  stopifnot(denominator>1e-12)
  numerator <- sum(w*score[rows]*y)
  list(multiplier=max(0,numerator/denominator),numerator=numerator,denominator=denominator,
    training_rows=rows,raw_signed_multiplier=numerator/denominator)
}
lap_scale <- calibrate(seq_len(nrow(calibration)))
cal_pred <- list()
for(cell in unique(calibration$depmap_id)) {
  rows <- which(calibration$depmap_id==cell)
  scale <- calibrate(setdiff(seq_len(nrow(calibration)),rows))
  cal_pred[[cell]] <- calibration[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=scale$multiplier*calibration_score[rows],multiplier=scale$multiplier,
      evaluation='held_cell_calibration')]
}
for(gene in unique(calibration$gene)) {
  rows <- which(calibration$gene==gene)
  scale <- calibrate(setdiff(seq_len(nrow(calibration)),rows))
  cal_pred[[gene]] <- calibration[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=scale$multiplier*calibration_score[rows],multiplier=scale$multiplier,
      evaluation='held_ligand_calibration')]
}
qsave(list(signature=model,lapatinib=lap_scale),file.path(root,'models/final_transfer.qs'))
write_table(data.table(multiplier=lap_scale$multiplier,raw_signed_multiplier=lap_scale$raw_signed_multiplier,
 observations=16,ligands='NRG1,EGF',doses='0.5,5',HGF_outcomes=0),'lapatinib_calibration')
write_table(cbind(calibration,raw_rescue_score=calibration_score,
 fitted_departure_log2=lap_scale$multiplier*calibration_score),'lapatinib_calibration_fitted')

# Keep the existing single-drug predictions, including the native 5 uM AU565 f.
previous_query <- fread(file.path(source_root,'data/au565_hgf_5uM_predictions_frozen.tsv'))[
 hgf_context=='Study-balanced HGF (same recipe as training)' & assumed_HGF_route==1]
stopifnot(nrow(previous_query)==1)
au <- source_profiles[CellLine=='AU565']
QX <- as.matrix(au[,..programs])
query_score <- as.numeric(QX %*% model$beta)
query <- data.table(source_profile=au$source_profile,cell='AU565',dose=5,
 lapatinib_pct_control=previous_query$lapatinib_pct_control,
 HGF_alone_pct_control=previous_query$hgf_increment_pct,raw_rescue_score=query_score,
 lapatinib_multiplier=lap_scale$multiplier,interaction_log2=lap_scale$multiplier*query_score)
query[,combination_pct_control:=lapatinib_pct_control*HGF_alone_pct_control/100*2^interaction_log2]
# Diagnostic transfer using the signature trained with every RAF/MEK HGF outcome excluded.
held_hgf_model <- fitted_models[['held_ligand_HGF']]
held_hgf_score <- as.numeric(CX %*% held_hgf_model$beta)
held_hgf_scale <- calibrate(seq_len(nrow(calibration)),held_hgf_score)
query[,all_HGF_outcomes_excluded_interaction:=held_hgf_scale$multiplier*as.numeric(QX %*% held_hgf_model$beta)]
query[,all_HGF_outcomes_excluded_combination_pct:=lapatinib_pct_control*HGF_alone_pct_control/100*
  2^all_HGF_outcomes_excluded_interaction]
write_table(query,'AU565_HGF_predictions_frozen')
write_table(rbindlist(predictions),'RAF_MEK_predictions_frozen')
write_table(rbindlist(cal_pred),'lapatinib_predictions_frozen')
freeze <- file.path(root,'data',c('AU565_HGF_predictions_frozen.tsv','RAF_MEK_predictions_frozen.tsv','lapatinib_predictions_frozen.tsv'))
write_table(data.table(path=freeze,sha256=vapply(freeze,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')

# Previously inspected outcomes: all evaluations remain retrospective.
scored <- merge(rbindlist(predictions),training[,.(observation_id,observed_departure_log2)],by='observation_id')
cs <- merge(rbindlist(cal_pred),calibration[,.(observation_id,observed_departure_log2)],by='observation_id')
write_table(scored,'RAF_MEK_predictions_scored')
write_table(cs,'lapatinib_predictions_scored')
metrics <- scored[,.(ligands=.N,r=pearson(observed_departure_log2,predicted_departure_log2)),
 by=.(evaluation,CellLine,drug,dose)]
write_table(metrics,'across_ligand_held_correlations')
cm <- cs[,.(cells=.N,r=pearson(observed_departure_log2,predicted_departure_log2),
 measured_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2)),by=.(evaluation,gene,dose)]
write_table(cm,'lapatinib_held_calibration_correlations')
hgf <- scored[gene=='HGF']
write_table(hgf,'HGF_source_held_predictions')

# Compact evidence plots; each panel compares ligands within one fixed cell/drug.
plot_data <- scored[evaluation=='held_ligand']
plot_data[,panel:=paste(CellLine,drug,sep=' | ')]
ann <- plot_data[,.(label=sprintf('r = %.2f',pearson(observed_departure_log2,predicted_departure_log2))),by=panel]
p <- ggplot(plot_data,aes(observed_departure_log2,predicted_departure_log2))+
 geom_point(aes(colour=gene=='HGF'),size=1.6)+scale_colour_manual(values=c('FALSE'='#555555','TRUE'='#D55E00'),guide='none')+
 geom_smooth(method='lm',formula=y~x,se=FALSE,colour='#0072B2',linewidth=.4)+
 geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=2.5)+
 facet_wrap(~panel,ncol=4,scales='free')+theme_bw(base_size=8)+theme(aspect.ratio=1)+
 labs(x='Measured growth departure (log2)',y='Predicted growth departure (log2)',
 title='v139 | each ligand excluded from its own signature fit',
 caption='Orange = HGF. Across 13 ligands per cell/drug. Retrospective evaluation; independent axis scales.')
ggsave(file.path(root,'figures/held_ligand_scatters.pdf'),p,width=9,height=7.2,device=cairo_pdf)
ggsave(file.path(root,'figures/held_ligand_scatters.png'),p,width=9,height=7.2,dpi=150)
source('/work/src/R/lib_heatmap.R')
betas <- dcast(rbindlist(coefficients),program~fit,value.var='beta')
B <- as.matrix(betas[match(programs,program),-1]); rownames(B) <- programs
ord <- hclust(dist(B),method='ward.D2')$order
ht <- draw_heatmap(B[ord,,drop=FALSE],scale_rows=FALSE,
 value_range=c(-max(abs(B)),0,max(abs(B))),value_colors=c('#2166AC','white','#B2182B'),
 hm_name='Raw-program beta',show_row_names=TRUE,show_column_names=TRUE,
 row_names_gp=grid::gpar(fontsize=7),column_names_gp=grid::gpar(fontsize=6),column_names_rot=60,
 height=grid::unit(135,'mm'),width=grid::unit(140,'mm'),title='One rescue signature | final and held-out fits')
cairo_pdf(file.path(root,'figures/rescue_signature_by_fold.pdf'),width=8.8,height=8)
ComplexHeatmap::draw(ht,heatmap_legend_side='right')
dev.off()

# Explicit scope and provenance; no grant documents or prior models are edited.
writeLines(c('v139: 40 raw signed ligand displacements -> one shared ridge rescue signature.',
 'No ligand-specific/program-specific standardization, no opposition products, no single-drug signature weights.',
 'RAF/MEK outcomes only learn beta. MEK multiplier fixed to 1; one positive relative RAF multiplier.',
 '80% across-ligand contrasts within cell/drug, 20% group mean; no intercept. Raw-zero input gives zero interaction.',
 'Ridge penalizes mean squared effective program weights across both drugs, so calibration cannot evade shrinkage.',
 'Ridge selected in inner whole-ligand folds. Outer leave-ligand-out and leave-cell-out predictions saved.',
 'Then beta frozen; 16 lapatinib NRG1/EGF outcomes fit one nonnegative multiplier across both doses.',
 'All lapatinib+HGF and Wilson outcomes excluded. Additional all-HGF-outcomes-excluded transfer is reported.',
 'Baseline enters transported molecular inputs only. Existing drug-alone and HGF-alone predictions retained.',
 'Drug-specific scalar calibration is a restricted pilot, not a general molecular model of arbitrary drug partners.',
 'RAF/MEK net-growth and lapatinib 72h abundance differ; scalar transfer is assessed on held calibration cells/ligands.',
 'AU565 5uM query uses existing single-drug baseline. It is not a dose/time-matched Wilson validation.',
 'All previously seen outcomes remain retrospective. No model or source chosen to achieve an HGF threshold.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
paths <- c('/work/src/program_pilot/139_raw_ligand_rescue_transfer.R',
 file.path(source_root,'data/external_combination_outcomes.tsv'),file.path(source_root,'data/external_ligand_programs.tsv'),
 file.path(source_root,'data/time_matched_HGF_program_inputs.tsv'),file.path(source_root,'models/final.qs'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(query[source_profile=='Primary_study_mean'])
print(metrics[,.(median_r=median(r),min_r=min(r),max_r=max(r)),by=evaluation])
print(hgf[,.(observed_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2)),by=.(evaluation,drug)])
print(cm)
cat('v139 complete.\n')
