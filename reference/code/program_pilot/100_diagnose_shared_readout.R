#!/usr/bin/env Rscript
# Diagnostic only: replay saved v100 models, decompose their existing objective
# and trace measured/predicted RNA through every decoder stage. No models fit.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest)})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
base_root <- '/work/results/program_pilot/100_shared_signatures_scaled_interactions'
root <- file.path(base_root,'diagnostic')
for (folder in c('data','figures')) dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
source_root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data'
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
programs <- sprintf('P%02d',1:40)
selection <- fread(file.path(source_root,'selected_doses.tsv'))
single <- fread(file.path(source_root,'individual_training_outcomes.tsv'))
pairs <- fread(file.path(source_root,'combination_training_outcomes.tsv'))
drugs <- selection$drug
short <- c('Azacytidine','Bendamustine','Bortezomib','Busulfan','Cytarabine',
           'Lapatinib','Topotecan','Vandetanib','Vismodegib','Crizotinib')
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
molecular <- merge(as.data.table(archive$prep$m),selection,by=c('drug','dose'))[depmap_id!='ACH-000248']
prism <- single[assay=='secondary']
nci <- single[assay!='secondary' & depmap_id!='ACH-000248']
profiles <- merge(molecular,prism[,.(depmap_id,drug,dose,S=-y)],by=c('depmap_id','drug','dose'))
rm(archive)
stopifnot(nrow(profiles)==286L,nrow(pairs)==193L,!anyDuplicated(profiles[,.(depmap_id,drug)]))
write_table <- function(x,name) fwrite(x,file.path(root,'data',paste0(name,'.tsv')),sep='\t')
pearson <- function(x,y) if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_
key <- function(cell,drug) paste(cell,drug,sep='|')

# Load only named function definitions from each frozen script. Its top-level
# training loops, optimizer and output-writing statements are never executed.
functions_needed <- c('route_node','route_back','predict_rna','readout_features','readout_forward',
                      'readout_gradient','geometry','loss','pair_features','gamma_matrix')
read_functions <- function(path) {
  env <- new.env(parent=globalenv())
  for (expression in parse(path)) {
    if (is.call(expression) && identical(expression[[1]],as.name('<-')) &&
        is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% functions_needed &&
        is.call(expression[[3]]) && identical(expression[[3]][[1]],as.name('function'))) {
      eval(expression,envir=env)
    }
  }
  stopifnot(all(vapply(functions_needed,exists,logical(1),envir=env,inherits=FALSE)))
  env
}
old_functions <- read_functions(file.path(base_root,'before_molecular_gate/100_shared_signatures_scaled_interactions.R'))
gate_functions <- read_functions('/work/src/program_pilot/100_shared_signatures_scaled_interactions.R')
configs <- CJ(implementation=c('original','gate'),arm=c('linear','quadratic'))
configs[,folder:=ifelse(implementation=='original',file.path(base_root,'before_molecular_gate',arm),file.path(base_root,arm))]
raw_associations <- rbindlist(lapply(seq_along(drugs),function(g) {
  p <- profiles[drug==drugs[g]]
  data.table(drug=drugs[g],dose=selection$dose[g],cells=nrow(p),program=programs,
             r=vapply(programs,function(name)pearson(p[[name]],p$S),numeric(1)))
}))
write_table(raw_associations,'observed_program_sensitivity_correlations')
write_table(profiles[,.(profiles=.N,cells=uniqueN(depmap_id),sensitive=sum(S>=1),
                        observed_median_pct=median(100*2^(-S))),by=.(drug,dose)],'matched_context_counts')
write_table(profiles[,.(profiles=.N,cells=uniqueN(depmap_id)),by=.(drug,dose,contributing_plates)],
            'molecular_plate_provenance')
stopifnot(!anyNA(match(prism$depmap_id,baseline$depmap_id)),
          all(is.finite(as.matrix(baseline[,..programs]))),all(is.finite(profiles$S)),
          !anyDuplicated(prism[,.(depmap_id,drug,dose)]))
write_table(pairs[,.(profiles=.N,cells=uniqueN(depmap_id),doseA=paste(sort(unique(doseA)),collapse=','),
                    doseB=paste(sort(unique(doseB)),collapse=',')),by=.(pair,SCREENER)],'combination_coverage')

stages <- contributions <- calibration_checks <- gradient_tables <- source_summaries <- list()
rna_checks <- input_checks <- gate_checks <- optimizer_checks <- pair_checks <- list()
held_stages <- held_traces <- list()
for (configuration in seq_len(nrow(configs))) {
  config <- configs[configuration]
  tag <- paste(config$implementation,config$arm,sep='_')
  cat('Auditing frozen',tag,'\n')
  functions <- if(config$implementation=='original') old_functions else gate_functions
  model <- qread(file.path(config$folder,'models/final.qs'))
  ctx <- model$context
  par <- model$par
  g <- match(profiles$drug,drugs)
  D <- as.matrix(profiles[,..programs])
  node <- functions$route_node(D,g,par,ctx)
  ungated <- if(config$implementation=='gate') node$base_activity else node$activity
  features <- functions$readout_features(node$activity,length(par$C))
  term_names <- c('route1','route2','route1_squared','route1_route2','route2_squared')[seq_along(par$C)]
  parts <- sweep(features,2,par$C,'*')
  pred <- rowSums(parts)
  frozen <- fread(file.path(config$folder,'data/measured_predictions_frozen.tsv'))[evaluation=='fitted']
  positions <- match(key(profiles$depmap_id,profiles$drug),key(frozen$depmap_id,frozen$drug))
  stopifnot(!anyNA(positions),max(abs(pred+frozen$predicted_log2[positions]))<1e-10)
  # Track transformations. At a fixed treatment, positive scaling and subtracting
  # a constant cannot alter any individual program's Pearson correlation.
  X <- (D-ctx$mu[g,,drop=FALSE])/ctx$sigma[g,,drop=FALSE]
  gain <- exp(ctx$Q %*% par$gain)
  shift <- (ctx$Q %*% par$shift) %*% t(ctx$shift_basis)
  adjusted <- X*gain[g,,drop=FALSE]+shift[g,,drop=FALSE]
  pre_tanh <- adjusted %*% ctx$W
  raw_scores <- X %*% ctx$W
  probability <- if(config$implementation=='gate') node$probability else rep(NA_real_,nrow(D))
  for (j in 1:10) {
    rows <- which(g==j)
    values <- cbind(raw_scores[rows,,drop=FALSE],pre_tanh[rows,,drop=FALSE],
                    ungated[rows,,drop=FALSE],node$activity[rows,,drop=FALSE],pred[rows])
    colnames(values) <- c('raw_signature1','raw_signature2','calibrated_pre_tanh1','calibrated_pre_tanh2',
                         'ungated_signature1','ungated_signature2','gated_signature1','gated_signature2','decoder')
    stages[[length(stages)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
      drug=drugs[j],dose=selection$dose[j],stage=colnames(values),cells=length(rows),
      r=vapply(seq_len(ncol(values)),function(k)pearson(values[,k],profiles$S[rows]),numeric(1)),
      mean=colMeans(values),SD=apply(values,2,sd))
    contributions[[length(contributions)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
      drug=drugs[j],term=term_names,coefficient=par$C,mean_contribution=colMeans(parts[rows,,drop=FALSE]),
      covariance_with_observed_S=vapply(seq_len(ncol(parts)),function(k)cov(parts[rows,k],profiles$S[rows]),numeric(1)),
      r=vapply(seq_len(ncol(parts)),function(k)pearson(parts[rows,k],profiles$S[rows]),numeric(1)))
    gate_checks[[length(gate_checks)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
      drug=drugs[j],template_route=model$template$membership[j],gate_mean=mean(probability[rows]),gate_SD=sd(probability[rows]),
      observed_mean_S=mean(profiles$S[rows]),predicted_mean_S=mean(pred[rows]),
      observed_SD=sd(profiles$S[rows]),predicted_SD=sd(pred[rows]),
      tanh_saturation1=mean(abs(pre_tanh[rows,1])>2),tanh_saturation2=mean(abs(pre_tanh[rows,2])>2))
    for (k in 1:40) {
      calibration_checks[[length(calibration_checks)+1L]] <- data.table(implementation=config$implementation,
        arm=config$arm,drug=drugs[j],program=programs[k],raw_r=pearson(D[rows,k],profiles$S[rows]),
        calibrated_r=pearson(adjusted[rows,k],profiles$S[rows]),gain=gain[j,k],shift=shift[j,k])
    }
  }
  # Save all molecular activities, predictions and exact readout contributions.
  activity_table <- cbind(profiles[,.(depmap_id,drug,dose,S)],
    data.table(raw_signature1=raw_scores[,1],raw_signature2=raw_scores[,2],
               ungated_signature1=ungated[,1],ungated_signature2=ungated[,2],
               gated_signature1=node$activity[,1],gated_signature2=node$activity[,2],
               gate=probability,predicted_S=pred))
  write_table(activity_table,paste0(tag,'_molecular_trace'))
  features_table <- as.data.table(parts)
  setnames(features_table,term_names)
  write_table(cbind(profiles[,.(depmap_id,drug,dose,S)],features_table),paste0(tag,'_readout_contributions'))

  # Reconstruct every training source and its optimization geometry exactly.
  query <- unique(rbind(nci[,.(depmap_id,drug)],pairs[,.(depmap_id,drug=drugA)],pairs[,.(depmap_id,drug=drugB)]))
  query_key <- key(query$depmap_id,query$drug)
  is <- match(key(nci$depmap_id,nci$drug),query_key)
  ia <- match(key(pairs$depmap_id,pairs$drugA),query_key)
  ib <- match(key(pairs$depmap_id,pairs$drugB),query_key)
  Dquery <- functions$predict_rna(model$rna,query$depmap_id,query$drug)
  native <- functions$route_node(Dquery,match(query$drug,drugs),par,ctx)
  R <- native$activity
  N <- R[is,,drop=FALSE]
  A <- R[ia,,drop=FALSE]
  B <- R[ib,,drop=FALSE]
  products <- sweep(functions$pair_features(A,B),2,ctx$interaction_scale,'/')
  J <- as.numeric(products %*% par$Gamma)
  scales <- exp(par$assay)
  assay <- match(nci$assay,c('nci_FF','nci_FG'))
  pair_assay <- match(pairs$SCREENER,c('FF','FG'))
  pn <- functions$readout_forward(N,par$C)*scales[assay]
  pp <- (functions$readout_forward(A,par$C)+functions$readout_forward(B,par$C)+J)*scales[pair_assay]
  pi <- J*scales[pair_assay]
  sources <- list(PRISM=list(pred=pred,geo=functions$geometry(profiles$S,profiles$drug),weight=1,group=profiles$drug),
    NCI_single=list(pred=pn,geo=functions$geometry(-nci$y,nci$group,nci$control_key),weight=.2,group=nci$group),
    NCI_total=list(pred=pp,geo=functions$geometry(-pairs$y,pairs$pair),weight=.2,group=pairs$pair),
    NCI_departure=list(pred=pi,geo=functions$geometry(-pairs$I,pairs$pair),weight=.2,group=pairs$pair))
  for (name in names(sources)) {
    source <- sources[[name]]
    source_table <- data.table(group=source$group,observed_S=source$geo$y,predicted_S=source$pred)
    source_summaries[[length(source_summaries)+1L]] <- source_table[,.(profiles=.N,
      observed_mean_S=mean(observed_S),predicted_mean_S=mean(predicted_S),
      observed_SD=sd(observed_S),predicted_SD=sd(predicted_S),r=pearson(observed_S,predicted_S)),by=group][,
      `:=`(implementation=config$implementation,arm=config$arm,source=name)]
    error <- source$pred-source$geo$y
    group_error <- as.numeric(rowsum(error*source$geo$w,source$geo$g,reorder=FALSE))/source$geo$gw
    row_gradients <- list(within=2*source$geo$w*(error-group_error[source$geo$g])/source$geo$within*source$weight,
                         mean=2*source$geo$w*group_error[source$geo$g]/source$geo$within*source$weight)
    for (component in names(row_gradients)) {
      gradient <- row_gradients[[component]]
      if(name=='PRISM') {
        Cgradient <- crossprod(features,gradient)
        back <- functions$route_back(node,gradient*functions$readout_gradient(node$activity,par$C),ctx)
        Gamma_gradient <- rep(0,3)
        assay_gradient <- rep(0,2)
      } else {
        native_gradient <- matrix(0,nrow(R),2)
        if(name=='NCI_single') {
          multiplied <- gradient*scales[assay]
          blocks <- list(list(index=is,value=multiplied*functions$readout_gradient(N,par$C)))
          Cgradient <- crossprod(functions$readout_features(N,length(par$C)),multiplied)
          Gamma_gradient <- rep(0,3)
          assay_gradient <- vapply(1:2,function(a)sum(gradient[assay==a]*pn[assay==a]),numeric(1))
        } else {
          multiplied <- gradient*scales[pair_assay]
          Gamma_gradient <- as.numeric(crossprod(products,multiplied))
          interaction_matrix <- functions$gamma_matrix(par$Gamma/ctx$interaction_scale)
          ga <- multiplied*(B %*% interaction_matrix)
          gb <- multiplied*(A %*% interaction_matrix)
          if(name=='NCI_total') {
            ga <- ga+multiplied*functions$readout_gradient(A,par$C)
            gb <- gb+multiplied*functions$readout_gradient(B,par$C)
            Cgradient <- crossprod(functions$readout_features(A,length(par$C))+
                                     functions$readout_features(B,length(par$C)),multiplied)
          } else Cgradient <- rep(0,length(par$C))
          blocks <- list(list(index=ia,value=ga),list(index=ib,value=gb))
          values <- if(name=='NCI_total') pp else pi
          assay_gradient <- vapply(1:2,function(a)sum(gradient[pair_assay==a]*values[pair_assay==a]),numeric(1))
        }
        for (block in blocks) {
          summed <- rowsum(block$value,block$index,reorder=FALSE)
          ii <- as.integer(rownames(summed))
          native_gradient[ii,] <- native_gradient[ii,,drop=FALSE]+summed
        }
        back <- functions$route_back(native,native_gradient,ctx)
      }
      vector <- unlist(c(back,list(C=as.numeric(Cgradient),assay=assay_gradient,Gamma=Gamma_gradient)))
      gradient_tables[[length(gradient_tables)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
        source=name,component=component,parameter=names(vector),gradient=as.numeric(vector))
      if(name=='PRISM') {
        for (drug in drugs) {
          rows <- which(profiles$drug==drug)
          drug_gradient <- as.numeric(crossprod(features[rows,,drop=FALSE],gradient[rows]))
          gradient_tables[[length(gradient_tables)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
            source=paste0('PRISM_drug:',drug),component=component,parameter=paste0('C',seq_along(par$C)),gradient=drug_gradient)
        }
      }
    }
  }
  penalty <- list(gain=2*par$gain,shift=2*par$shift)
  if(config$implementation=='gate') penalty$gate <- .2*par$gate
  penalty <- c(penalty,list(C=.01*par$C,assay=.02*par$assay,Gamma=.1*par$Gamma))
  vector <- unlist(penalty)
  gradient_tables[[length(gradient_tables)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
    source='penalty',component='penalty',parameter=names(vector),gradient=as.numeric(vector))

  # Pair product information and endpoint units are evaluated without a refit.
  for (pair in unique(pairs$pair)) {
    rows <- which(pairs$pair==pair)
    F <- products[rows,,drop=FALSE]
    centered <- sweep(F,2,colMeans(F))
    singular <- svd(centered,nu=0,nv=0)$d
    pair_checks[[length(pair_checks)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
      pair=pair,term=c('route1A_route1B','cross_route','route2A_route2B'),
      feature_SD=apply(F,2,sd),feature_departure_r=vapply(1:3,function(k)pearson(F[,k],-pairs$I[rows]),numeric(1)),
      singular_values=singular,coefficient=par$Gamma,
      observed_departure_SD=sd(pairs$I[rows]),predicted_departure_SD=sd(pi[rows]))
  }
  audit <- fread(file.path(config$folder,'data/fit_audit.tsv'))
  optimizer_checks[[tag]] <- audit[,`:=`(implementation=config$implementation,arm=config$arm)]
  # Held RNA and held function are compared in exactly the same 286 profiles.
  folds <- fread(file.path(config$folder,'data/recipient_folds.tsv'))
  for (cell in unique(profiles$depmap_id)) {
    fold <- folds[depmap_id==cell,fold]
    held_model <- qread(file.path(config$folder,'models',sprintf('loo_%02d.qs',fold)))
    rows <- which(profiles$depmap_id==cell)
    measured_D <- D[rows,,drop=FALSE]
    predicted_D <- functions$predict_rna(held_model$rna,profiles$depmap_id[rows],profiles$drug[rows])
    measured_node <- functions$route_node(measured_D,match(profiles$drug[rows],drugs),held_model$par,held_model$context)
    predicted_node <- functions$route_node(predicted_D,match(profiles$drug[rows],drugs),held_model$par,held_model$context)
    held_X <- (measured_D-held_model$context$mu[match(profiles$drug[rows],drugs),,drop=FALSE])/
      held_model$context$sigma[match(profiles$drug[rows],drugs),,drop=FALSE]
    held_raw <- held_X %*% held_model$context$W
    held_ungated <- if(config$implementation=='gate') measured_node$base_activity else measured_node$activity
    held_traces[[length(held_traces)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
      depmap_id=cell,drug=profiles$drug[rows],S=profiles$S[rows],
      raw_signature1=held_raw[,1],raw_signature2=held_raw[,2],
      ungated_signature1=held_ungated[,1],ungated_signature2=held_ungated[,2],
      gated_signature1=measured_node$activity[,1],gated_signature2=measured_node$activity[,2],
      measured_decoder=functions$readout_forward(measured_node$activity,held_model$par$C),
      predicted_RNA_decoder=functions$readout_forward(predicted_node$activity,held_model$par$C))
    for (i in seq_along(rows)) {
      rna_checks[[length(rna_checks)+1L]] <- data.table(implementation=config$implementation,arm=config$arm,
        depmap_id=cell,drug=profiles$drug[rows[i]],program=programs,observed_D=measured_D[i,],predicted_D=predicted_D[i,])
    }
    stopifnot(!cell %in% held_model$training_cells,!cell %in% held_model$rna$training_cells)
  }
  # Whole-cell fold stability of signature sign and support.
  template_stability <- rbindlist(lapply(seq_len(nrow(folds)),function(index) {
    held_model <- qread(file.path(config$folder,'models',sprintf('loo_%02d.qs',index)))
    data.table(fold=index,program=rep(programs,2),route=rep(1:2,each=40),
               beta=as.numeric(held_model$template$W))
  }))
  write_table(template_stability,paste0(tag,'_template_fold_stability'))
}
write_table(rbindlist(stages),'transformation_stage_correlations')
write_table(rbindlist(contributions),'readout_covariance_decomposition')
write_table(rbindlist(calibration_checks),'calibration_correlation_preservation')
write_table(rbindlist(gate_checks),'route_selection_and_saturation')
write_table(rbindlist(gradient_tables),'objective_gradient_decomposition')
write_table(rbindlist(source_summaries),'training_source_endpoint_calibration')
write_table(rbindlist(optimizer_checks),'optimizer_and_leakage_checks')
write_table(rbindlist(pair_checks),'combination_feature_diagnosis')
rna_data <- rbindlist(rna_checks)
write_table(rna_data,'held_rna_predictions')
write_table(rna_data[,.(cells=.N,r=pearson(observed_D,predicted_D),observed_SD=sd(observed_D),
                       predicted_SD=sd(predicted_D),SD_ratio=sd(predicted_D)/sd(observed_D)),
                    by=.(implementation,arm,drug,program)],'held_rna_program_agreement')
held_trace <- rbindlist(held_traces)
write_table(held_trace,'held_cell_molecular_trace')
influence <- rbindlist(lapply(c('raw_signature1','raw_signature2','measured_decoder'),function(stage) {
  x <- held_trace[implementation=='original' & arm=='linear' & drug=='Lapatinib ditosylate']
  data.table(stage=stage,omitted_cell=x$depmap_id,omitted_S=x$S,
    all_cell_r=pearson(x$S,x[[stage]]),
    diagnostic_r_without_cell=vapply(seq_len(nrow(x)),function(i)pearson(x$S[-i],x[[stage]][-i]),numeric(1)))
}))
write_table(influence,'lapatinib_association_cell_influence_no_filtering')
held_long <- melt(held_trace,id.vars=c('implementation','arm','depmap_id','drug','S'),
                  variable.name='stage',value.name='value')
write_table(held_long[,.(cells=.N,r=pearson(S,value),mean=mean(value),SD=sd(value)),
                     by=.(implementation,arm,drug,stage)],'held_transformation_stage_correlations')
# A converged objective should have cancelling source/penalty gradients. This
# distinguishes an implementation error from a poor compromise under the loss.
gradients <- rbindlist(gradient_tables)
stationarity <- gradients[!startsWith(source,'PRISM_drug:'),.(net_gradient=sum(gradient)),
                          by=.(implementation,arm,parameter)]
write_table(stationarity,'objective_stationarity')
print(stationarity[,.(maximum_gradient=max(abs(net_gradient))),by=.(implementation,arm)])
stopifnot(max(abs(stationarity$net_gradient))<.02)
calibration <- rbindlist(calibration_checks)
stopifnot(max(abs(calibration$raw_r-calibration$calibrated_r),na.rm=TRUE)<1e-10)

# Descriptive variance partition, never a model performance criterion.
observed_groups <- profiles[,.(observed_mean=mean(S)),by=drug]
group_means <- observed_groups$observed_mean[match(profiles$drug,observed_groups$drug)]
write_table(data.table(profiles=nrow(profiles),cells=uniqueN(profiles$depmap_id),
                      between_treatment_fraction=sum((group_means-mean(profiles$S))^2)/
                        sum((profiles$S-mean(profiles$S))^2)), 'observed_endpoint_variance_partition')
summary <- rbindlist(lapply(seq_len(nrow(configs)),function(i) {
  config <- configs[i]
  measured <- fread(file.path(config$folder,'data/measured_decoder_correlations.tsv'))[evaluation=='held_cell']
  bridge <- fread(file.path(config$folder,'data/baseline_bridge_correlations.tsv'))[evaluation=='held_cell']
  setnames(measured,c('cells','r'),c('matched_cells','held_measured_decoder_r'))
  setnames(bridge,c('cells','r'),c('all_prism_cells','held_baseline_bridge_r'))
  table <- merge(measured[,.(drug,dose,matched_cells,held_measured_decoder_r)],bridge,by=c('drug','dose'))
  table[,`:=`(implementation=config$implementation,arm=config$arm)]
  table
}))
write_table(summary,'held_endpoint_summary')
held_stage_stats <- held_long[,.(cells=.N,r=pearson(S,value)),by=.(implementation,arm,drug,stage)]
association_matrix <- dcast(raw_associations,program~drug,value.var='r')
association_matrix <- as.matrix(association_matrix[match(programs,program),..drugs])
rownames(association_matrix) <- programs
colnames(association_matrix) <- short
ro <- hclust(dist(association_matrix),method='ward.D2')$order
co <- hclust(dist(t(association_matrix)),method='ward.D2')$order
stages_to_plot <- c('raw_signature1','raw_signature2','ungated_signature1','ungated_signature2',
                    'measured_decoder','predicted_RNA_decoder')
stage_table <- held_stage_stats[implementation=='original' & arm=='linear']
stage_matrix <- matrix(NA_real_,length(stages_to_plot),10,dimnames=list(stages_to_plot,drugs))
for (s in stages_to_plot) {
  z <- stage_table[stage==s]
  stage_matrix[s,] <- z$r[match(drugs,z$drug)]
}
colnames(stage_matrix) <- short
rownames(stage_matrix) <- c('Held raw signature 1','Held raw signature 2','Held calibrated signature 1',
                            'Held calibrated signature 2','Held measured decoder','Held baseline bridge (matched cells)')
render_audit_heatmap <- function() {
  left <- draw_heatmap(association_matrix[ro,co],scale_rows=FALSE,hm_name='Observed Pearson r',
    value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),
    show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=8),
    column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
    title='Observed program associations with sensitivity',width=grid::unit(105,'mm'))
  ComplexHeatmap::draw(left,heatmap_legend_side='bottom')
}
cairo_pdf(file.path(root,'figures/observed_signal_all_drugs.pdf'),width=7.5,height=9)
render_audit_heatmap()
dev.off()
cairo_pdf(file.path(root,'figures/held_signal_through_bridge.pdf'),width=9,height=4.5)
heat <- draw_heatmap(stage_matrix[,co],scale_rows=FALSE,hm_name='Held Pearson r',
  value_range=c(-1,0,1),value_colors=c('#2166AC','white','#B2182B'),
  show_row_names=TRUE,show_column_names=TRUE,row_names_gp=grid::gpar(fontsize=8),
  column_names_rot=45,column_names_gp=grid::gpar(fontsize=8),
  title='Original v100 linear: where held-cell sensitivity signal is lost')
ComplexHeatmap::draw(heat,heatmap_legend_side='bottom')
dev.off()

# Lapatinib panels all use whole-cell held-out quantities on the same 30 cells.
# The observed-program panel is explicitly a descriptive, previously exposed
# association; it is not presented as a trained predictor.
lap <- held_trace[implementation=='original' & arm=='linear' & drug=='Lapatinib ditosylate']
best <- raw_associations[drug=='Lapatinib ditosylate'][which.max(abs(r)),program]
raw_lap <- profiles[drug=='Lapatinib ditosylate']
lap[,observed_program:=raw_lap[[best]][match(depmap_id,raw_lap$depmap_id)]]
lap_plot <- melt(lap,id.vars=c('depmap_id','drug','S','implementation','arm'),
                 measure.vars=c('observed_program','raw_signature1','raw_signature2',
                                'ungated_signature1','measured_decoder','predicted_RNA_decoder'))
labels <- c(paste0('Observed ',best,' (descriptive)'), 'Held raw signature 1','Held raw signature 2',
            'Held calibrated signature 1','Held shared measured decoder','Held baseline bridge')
lap_plot[,panel:=factor(variable,levels=levels(variable),labels=labels)]
stats <- lap_plot[,.(label=paste0('r = ',sprintf('%.3f',pearson(S,value)),'\nn = ',.N)),by=panel]
plot <- ggplot(lap_plot,aes(S,value))+geom_point(aes(color=S>=1),size=1.5)+
  scale_color_manual(values=c('FALSE'='#0072B2','TRUE'='#D55E00'),labels=c('>=50% abundance','<50% abundance'))+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
  geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
  facet_wrap(~panel,ncol=3,scales='free')+theme_bw(base_size=9)+
  theme(aspect.ratio=1,panel.grid.minor=element_blank(),legend.position='bottom')+
  labs(title='Lapatinib 5 µM: signal survives cell holdout, then the shared decoder reverses it',
       x='Observed sensitivity (-log2 PRISM abundance)',y='Program / intermediate / predicted sensitivity',color=NULL,
       caption='Same 30 recipients. Orange: measured abundance below 50%. Higher observed sensitivity means lower abundance.')
ggsave(file.path(root,'figures/lapatinib_signal_loss.pdf'),plot,width=10,height=7,device=cairo_pdf)
ggsave(file.path(root,'figures/lapatinib_signal_loss.png'),plot,width=10,height=7,dpi=150)

gradient_plot_data <- gradients[implementation=='original' & arm=='linear' & parameter %in% c('C1','C2') &
                                 !startsWith(source,'PRISM_drug:')]
gradient_plot_data[,source_label:=paste(source,component)]
gradient_plot <- ggplot(gradient_plot_data,aes(source_label,gradient,fill=component))+
  geom_col()+geom_hline(yintercept=0,color='grey30')+facet_wrap(~parameter,ncol=2)+coord_flip()+
  theme_bw(base_size=10)+theme(legend.position='bottom')+
  labs(title='Original v100: opposing training pressures on shared coefficients',x=NULL,
       y='Training objective derivative',
       caption='Positive derivative pushes a coefficient down; negative pushes it up. These are optimization diagnostics.')
ggsave(file.path(root,'figures/source_gradient_conflict.pdf'),gradient_plot,width=9,height=5.5,device=cairo_pdf)

# Baseline map agreement for the strongest observed lapatinib programs.
top_programs <- raw_associations[drug=='Lapatinib ditosylate'][order(-abs(r))][1:6,program]
rna_plot_data <- rna_data[implementation=='original' & arm=='linear' &
                           drug=='Lapatinib ditosylate' & program %in% top_programs]
rna_stats <- rna_plot_data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed_D,predicted_D)),
                                          '\nSD ratio = ',sprintf('%.2f',sd(predicted_D)/sd(observed_D)))),by=program]
rna_plot <- ggplot(rna_plot_data,aes(observed_D,predicted_D))+geom_point(color='#0072B2',size=1.3)+
  geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
  geom_text(data=rna_stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=2.8)+
  facet_wrap(~program,ncol=3,scales='free')+theme_bw(base_size=9)+theme(aspect.ratio=1)+
  labs(title='Lapatinib: the RNA-only baseline map weakly recovers sensitivity-associated displacements',
       x='Measured displacement',y='Held-cell baseline-predicted displacement')
ggsave(file.path(root,'figures/lapatinib_baseline_to_program_agreement.pdf'),rna_plot,width=9,height=6.3,device=cairo_pdf)
q <- fread(file.path(base_root,'linear/data/au565_hgf_5uM_predictions_frozen.tsv'))
write_table(q,'hgf_query_components_and_scale_qualification')
inputs <- c('/work/src/program_pilot/100_diagnose_shared_readout.R',
            '/work/src/program_pilot/100_shared_signatures_scaled_interactions.R',
            file.path(base_root,'before_molecular_gate/100_shared_signatures_scaled_interactions.R'),
            configs[,file.path(folder,'models/final.qs')],archive_path,
            file.path(source_root,c('selected_doses.tsv','individual_training_outcomes.tsv','combination_training_outcomes.tsv')))
write_table(data.table(path=inputs,sha256=vapply(inputs,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
writeLines(c('Diagnostic only: no model fitting, tuning, optimization or modified functional training.',
             'Four existing v100 implementations audited: original/gated x linear/quadratic.',
             'Training predictions reconstructed exactly; raw/calibrated program correlations agree to numerical precision.',
             'Training gradient components split by source and within-treatment versus treatment means; net gradients match converged fits.',
             'Whole-cell signature scores, measured decoder, and predicted-RNA decoder replayed on the same paired recipients.',
             'All 40 held-RNA program predictions audited; pair feature associations/rank and departure amplitude audited.',
             'Clinical sensitivity/resistance and HGF rescue are not established by program correlations alone.',
             'Previously exposed labels/dose selection make these retrospective diagnostics.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
cat('All frozen-model reconstructions and gradient decompositions verified. No models fitted.\n')
