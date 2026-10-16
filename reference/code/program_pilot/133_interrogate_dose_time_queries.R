#!/usr/bin/env Rscript
# Fixed-model queries. No Wilson or AU565 combination outcome is used here.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
expressions<-parse('/work/src/program_pilot/133_counts_and_growth_order_bridge.R')
is_assignment<-function(x,name)is.call(x)&&identical(x[[1]],as.name('<-'))&&
  is.symbol(x[[2]])&&identical(as.character(x[[2]]),name)
begin_fit<-which(vapply(expressions,is_assignment,logical(1),name='original_final'))
stopifnot(length(begin_fit)==1L)
for(expression in expressions[seq_len(begin_fit-1L)])eval(expression,envir=.GlobalEnv)
model<-qread(file.path(root,'models/final_normalized_cell_count.qs'))
ctx<-model$context
profiles<-fread(file.path(root,'data/time_matched_HGF_program_inputs.tsv'))[CellLine=='AU565']
Dh<-as.matrix(profiles[,..programs])
Qh<-cbind(1,sweep(sweep(sweep(Dh,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
scale_map<-solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
  crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
minimum<-matrix(apply(ctx$sigma[1:10,],2,min),nrow(Dh),40,byrow=TRUE)
maximum<-matrix(apply(ctx$sigma[1:10,],2,max),nrow(Dh),40,byrow=TRUE)
hscale<-pmin(pmax(exp(Qh %*% scale_map),minimum),maximum)
H40<-Dh/hscale
input<-ifelse(profiles$source_profile=='Primary_study_mean','HGF',profiles$source_profile)
H<-cbind(H40,H40 %*% ctx$W,ligand_coordinate(rep('AU565',nrow(Dh)),input))
ligand_growth5<- -((H40 %*% ctx$W)*as.numeric(exp(Qh %*% model$par$log_amplitude)))[,1]*
  model$ligand_single_conversion
low<-predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=.5))
high<-predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=5))
queries<-list()
for(ds in c(.5,.75,1,2.5,5)) {
  weight<-log(ds/.5)/log(5/.5)
  A<-(1-weight)*low$program_activity+weight*high$program_activity
  sensitivity5<-(1-weight)*low$sensitivity+weight*high$sensitivity
  phi<-pair_features(A[rep(1L,nrow(H)),,drop=FALSE],H,ctx$W,model$interaction_basis)
  interaction5<- -as.numeric(sweep(phi,2,ctx$interaction_scale,'/') %*% model$par$Gamma)
  for(days in c(3,5,12)) {
    t<-days/model$interaction_reference_days
    alone<-100*2^(-sensitivity5*t)
    combined<-100*2^((-sensitivity5+ligand_growth5+interaction5)*t)
    queries[[paste(ds,days)]]<-data.table(cell='AU565',dose_uM=ds,days,
      source_profile=profiles$source_profile,lapatinib_pct_vehicle=alone,
      HGF_alone_pct_vehicle=100*2^(ligand_growth5*t),
      combination_pct_vehicle=combined,interaction_log2_native=interaction5*t,
      conditional_drug_pct_HGF_control=100*2^((-sensitivity5+interaction5)*t),
      restored_lost_viability_pct=100*(combined-alone)/(100-alone),
      qualification=if(days==5)'five-day reference prediction' else
        'duration projection assuming constant relative net growth; state adaptation unmodeled',
      evaluation='HCC1954 HGF training allowed; AU565 HGF and Wilson excluded')
  }
}
query<-rbindlist(queries)
write_table(query,'AU565_HGF_dose_time_predictions_frozen')
freeze_path<-file.path(root,'data/AU565_HGF_dose_time_predictions_frozen.tsv')
write_table(data.table(path=freeze_path,sha256=digest(freeze_path,file=TRUE,algo='sha256')),
  'AU565_dose_time_prediction_hash')
print(query[source_profile=='Primary_study_mean' & dose_uM %in% c(1,5)])
