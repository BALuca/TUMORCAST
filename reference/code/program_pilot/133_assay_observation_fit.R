# Same f/g molecular architecture as132. Counts determine magnitude; normalized
# net-growth measurements inform only within-drug/cell ordering across ligands.
reference_days <- 5
count_assay_days <- 3

fit_interaction <- function(original,mfeatures,endpoint_name,excluded_gene=character(),held_cell=character(),
                            penalty=interaction_penalty) {
  original$auxiliary_rna<-fit_auxiliary_rna(original,held_cell)
  train_pairs<-pairs[!depmap_id %in% held_cell]
  A<-predict_single(original,train_pairs[,.(depmap_id,drug=drugA,dose=doseA)])
  B<-predict_single(original,train_pairs[,.(depmap_id,drug=drugB,dose=doseB)])
  ef<-external_features(original,held_cell)
  ei<-which(!external$gene %in% excluded_gene & !external$depmap_id %in% held_cell)
  ev<-which(ef$view_rows %in% ei)
  mr<-which(!ligands$gene %in% excluded_gene)
  molecular_rows<-rbind(A$program_activity,B$program_activity,mfeatures$B[mr,,drop=FALSE],
    ef$Aviews[ev,,drop=FALSE],ef$Bviews[ev,,drop=FALSE])
  groups<-c(paste0('drug_',train_pairs$drugA),paste0('drug_',train_pairs$drugB),
    paste0('ligand_',ligands$gene[mr]),paste0('drug_',external$drug[ef$view_rows[ev]]),
    paste0('ligand_',external$gene[ef$view_rows[ev]],'_',ef$view_labels[ev]))
  group_levels<-sort(unique(groups))
  centres<-t(vapply(group_levels,function(g)colMeans(molecular_rows[groups==g,1:40,drop=FALSE]),numeric(40)))
  distances<-as.matrix(dist(centres));diag(distances)<-Inf
  bandwidth<-2*median(apply(distances,1,min))
  stopifnot(is.finite(bandwidth),bandwidth>1e-8)
  basis<-list(centres=centres,labels=group_levels,bandwidth=bandwidth,
    method='unchanged RNA-only molecular prototypes; no functional labels')
  Pn<-pair_features(A$program_activity,B$program_activity,original$context$W,basis)
  # Existing single-agent NCI conversions already carry their native observation
  # calibration. Gamma remains in the same five-day PRISM-equivalent count units.
  nci_scale<-exp(original$par$assay)[match(train_pairs$SCREENER,c('FF','FG'))]
  ng<-center_geometry(Pn*nci_scale,-train_pairs$I,
    paste(train_pairs$pair,train_pairs$SCREENER,train_pairs$doseA,train_pairs$doseB,sep='|'))
  train_labels<-training_labels[training_labels$endpoint==endpoint_name & !training_labels$gene %in% excluded_gene &
    !CellLine %in% ligands[depmap_id %in% held_cell,CellLine]]
  row<-match(paste(train_labels$CellLine,train_labels$gene),paste(ligands$CellLine,ligands$gene))
  stopifnot(!anyNA(row),!any(train_labels$gene=='NRG1'),
    !any(train_labels$CellLine=='AU565' & train_labels$gene=='HGF'))
  Pm<-pair_features(mfeatures$A,mfeatures$B,original$context$W,basis)*count_assay_days/reference_days
  mg<-center_geometry(Pm[row,,drop=FALSE],-train_labels$y,train_labels$CellLine)
  Pe<-pair_features(ef$Aviews,ef$Bviews,original$context$W,basis)
  er<-ef$view_rows[ev]

  # Estimate one zero-referenced ligand main-effect conversion using only
  # independently measured final-cell abundance controls, never net growth.
  count_views<-ev[external$study[er]=='elife32271']
  count_rows<-ef$view_rows[count_views]
  ids<-paste(external$depmap_id[count_rows],external$gene[count_rows],ef$view_labels[count_views])
  keep<-!duplicated(ids)
  if(any(keep)) {
    p<-ef$standalone_views[count_views][keep]*count_assay_days/reference_days
    target_growth<-external$observed_ligand_log2[count_rows][keep]
    stopifnot(all(is.finite(target_growth)))
    ligand_conversion<-max(0,sum(p*target_growth)/sum(p^2))
  } else ligand_conversion<-0

  # Breast cancer experiments measure final abundance after72h. They supply
  # an identifiable magnitude calibration for the common five-day g.
  if(length(count_views)) {
    eg<-center_geometry(Pe[count_views,,drop=FALSE]*count_assay_days/reference_days,
      -external$observed_departure_log2[count_rows],
      paste(external$drug[count_rows],external$dose[count_rows],external$gene[count_rows],sep='|'))
    eg$Xm<-eg$X-eg$Xc;eg$ym<-eg$y-eg$yc
    eg$variance<-max(sum(eg$w*eg$yc^2),1e-8)
    eg$mean_variance<-max(sum(eg$w*eg$ym^2),1e-8)
  } else {
    eg<-list(X=matrix(0,0,ncol(Pe)),Xc=matrix(0,0,ncol(Pe)),Xm=matrix(0,0,ncol(Pe)),
      y=numeric(),yc=numeric(),ym=numeric(),w=numeric(),variance=1,mean_variance=1)
  }

  # For the GFP study, drug-arm growth ordering is invariant to unknown initial
  # cell counts. Its baseline-subtracted fold values cannot calibrate final-cell
  # abundance. Each unordered ligand pair retains one weight across RNA views.
  rank_parts<-list();rank_targets<-rank_bases<-rank_groups<-rank_pairs<-list()
  growth_rows<-ei[external$study[ei]=='Straussman2012']
  recipient_groups<-unique(paste(external$depmap_id[growth_rows],external$drug[growth_rows],sep='|'))
  for(group in recipient_groups) {
    rows<-growth_rows[paste(external$depmap_id[growth_rows],external$drug[growth_rows],sep='|')==group]
    index<-combn(rows,2)
    for(j in seq_len(ncol(index))) {
      i1<-index[1,j];i2<-index[2,j]
      direction<-sign(external$observed_drug_arm_increment_log2[i1]-external$observed_drug_arm_increment_log2[i2])
      if(direction==0)next
      v1<-which(ef$view_rows==i1);v2<-which(ef$view_rows==i2)
      combinations<-CJ(v1=v1,v2=v2)
      key<-paste(group,i1,i2,sep='|')
      rank_parts[[key]]<-Pe[combinations$v1,,drop=FALSE]-Pe[combinations$v2,,drop=FALSE]
      rank_targets[[key]]<-rep(direction,nrow(combinations))
      rank_bases[[key]]<-(ef$standalone_views[combinations$v1]-ef$standalone_views[combinations$v2])*ligand_conversion
      rank_groups[[key]]<-rep(group,nrow(combinations));rank_pairs[[key]]<-rep(key,nrow(combinations))
    }
  }
  if(length(rank_parts)) {
    Xr_raw<-do.call(rbind,rank_parts);yr<-unlist(rank_targets,use.names=FALSE)
    br<-unlist(rank_bases,use.names=FALSE);rg<-unlist(rank_groups,use.names=FALSE);rp<-unlist(rank_pairs,use.names=FALSE)
    number_pairs<-unique(data.table(group=rg,pair=rp))[,.(pairs=.N),by=group]
    wr<-1/as.numeric(table(rp)[rp])/number_pairs$pairs[match(rg,number_pairs$group)]/uniqueN(rg)
    stopifnot(abs(sum(wr)-1)<1e-12)
  } else {
    Xr_raw<-matrix(0,0,ncol(Pe));yr<-br<-wr<-numeric()
  }
  product_scale<-sqrt((colSums(ng$Xc^2*ng$w)+colSums(mg$Xc^2*mg$w)+
    external_within_weight*colSums(eg$Xc^2*eg$w)+(1-external_within_weight)*colSums(eg$Xm^2*eg$w)+
    colSums(Xr_raw^2*wr))/4)
  supported<-product_scale>1e-10;product_scale[!supported]<-1
  scale_X<-function(X) {
    X<-sweep(X,2,product_scale,'/');X[,!supported]<-0;X
  }
  Xn<-scale_X(ng$Xc);Xm<-scale_X(mg$Xc);Xec<-scale_X(eg$Xc);Xem<-scale_X(eg$Xm);Xr<-scale_X(Xr_raw)
  ng$variance<-max(ng$variance,1e-8);mg$variance<-max(mg$variance,1e-8)
  p<-ncol(Xn)
  gram<-source_weight*(crossprod(Xn,Xn*ng$w)/ng$variance+crossprod(Xm,Xm*mg$w)/mg$variance+
    external_within_weight*crossprod(Xec,Xec*eg$w)/eg$variance+
    (1-external_within_weight)*crossprod(Xem,Xem*eg$w)/eg$mean_variance)+diag(penalty,p)
  rhs<-source_weight*(as.numeric(crossprod(Xn,ng$w*ng$yc))/ng$variance+
    as.numeric(crossprod(Xm,mg$w*mg$yc))/mg$variance+
    external_within_weight*as.numeric(crossprod(Xec,eg$w*eg$yc))/eg$variance+
    (1-external_within_weight)*as.numeric(crossprod(Xem,eg$w*eg$ym))/eg$mean_variance)
  initial<-c(as.numeric(solve(gram,rhs)),0,0)
  objective<-function(parameters,gradient=FALSE) {
    gamma<-parameters[seq_len(p)];lm<-parameters[p+1];lr<-parameters[p+2]
    sm<-exp(lm);sr<-exp(lr)
    rn<-ng$yc-as.numeric(Xn %*% gamma)
    vm<-as.numeric(Xm %*% gamma);rm<-mg$yc-sm*vm
    rc<-eg$yc-as.numeric(Xec %*% gamma);re<-eg$ym-as.numeric(Xem %*% gamma)
    difference<-br-as.numeric(Xr %*% gamma)
    margin<-yr*sr*difference
    rank_loss<-pmax(-margin,0)+log1p(exp(-abs(margin)))
    loss<-source_weight*(sum(ng$w*rn^2)/ng$variance+sum(mg$w*rm^2)/mg$variance+
      external_within_weight*sum(eg$w*rc^2)/eg$variance+
      (1-external_within_weight)*sum(eg$w*re^2)/eg$mean_variance+sum(wr*rank_loss))+
      penalty*sum(gamma^2)+assay_penalty*(lm^2+lr^2)
    if(!gradient)return(loss)
    error_probability<-plogis(-margin)
    gg<- -2*source_weight*(as.numeric(crossprod(Xn,ng$w*rn))/ng$variance+
      sm*as.numeric(crossprod(Xm,mg$w*rm))/mg$variance+
      external_within_weight*as.numeric(crossprod(Xec,eg$w*rc))/eg$variance+
      (1-external_within_weight)*as.numeric(crossprod(Xem,eg$w*re))/eg$mean_variance)+
      source_weight*sr*as.numeric(crossprod(Xr,wr*yr*error_probability))+2*penalty*gamma
    gm<- -2*source_weight*sm*sum(mg$w*rm*vm)/mg$variance+2*assay_penalty*lm
    gr<- -source_weight*sum(wr*margin*error_probability)+2*assay_penalty*lr
    c(gg,gm,gr)
  }
  lower<-c(rep(-Inf,p),-4,-4);upper<-c(rep(Inf,p),4,4)
  optimum<-optim(initial,objective,function(z)objective(z,TRUE),method='L-BFGS-B',
    lower=lower,upper=upper,control=list(maxit=400,factr=1e7,pgtol=1e-6))
  if(optimum$convergence!=0L)optimum<-optim(optimum$par,objective,function(z)objective(z,TRUE),method='L-BFGS-B',
    lower=lower,upper=upper,control=list(maxit=1000,factr=1e7,pgtol=1e-6))
  stopifnot(optimum$convergence==0L,all(is.finite(optimum$par)))
  gradient_difference<-NA_real_
  if(!length(held_cell)&&!length(excluded_gene)) {
    indices<-unique(c(1,p,p+1,p+2));eps<-1e-6
    numerical<-vapply(indices,function(i) {
      hi<-lo<-optimum$par;hi[i]<-hi[i]+eps;lo[i]<-lo[i]-eps
      (objective(hi)-objective(lo))/(2*eps)
    },numeric(1))
    gradient_difference<-max(abs(numerical-objective(optimum$par,TRUE)[indices]))
    stopifnot(gradient_difference<1e-4)
  }
  model<-original;model$par$Gamma<-optimum$par[seq_len(p)]
  model$context$interaction_scale<-product_scale;model$context$interaction_supported<-supported
  model$mema_assay_scale<-exp(optimum$par[p+1]);model$growth_rank_temperature<-exp(optimum$par[p+2])
  model$interaction_basis<-basis;model$external_rna<-ef$rna
  model$ligand_single_conversion<-ligand_conversion
  model$interaction_reference_days<-reference_days
  model$interaction_training<-list(method='actual-count calibration + within-drug/cell net-growth ordering',
    endpoint=endpoint_name,challenge_genes='NRG1',challenge_context='AU565+HGF; all Wilson outcomes',
    training_genes=unique(train_labels$gene),excluded_gene=excluded_gene,held_cell=held_cell,
    source_weight=source_weight,interaction_penalty=penalty,assay_penalty=assay_penalty,
    external_within_weight=external_within_weight,reference_days=reference_days,
    count_assay_days=count_assay_days,external_rows=length(ei),external_molecular_views=length(ev),
    exact_count_rows=length(count_rows),growth_ranking_comparisons=length(yr),
    external_HGF_rows=sum(external$gene[ei]=='HGF'),ligand_single_conversion=ligand_conversion,
    iterations=optimum$counts[['function']],converged=TRUE,
    normal_residual=max(abs(objective(optimum$par,TRUE)[seq_len(p)])),
    checked_gradient_difference=gradient_difference,
    growth_readout_qualification='normalized baseline-subtracted growth supplies order, not count-effect magnitude')
  stopifnot(identical(original$rna,model$rna),identical(original$template,model$template),
    identical(original$par$level,model$par$level),identical(original$par$log_amplitude,model$par$log_amplitude),
    identical(original$par$assay,model$par$assay))
  list(model=model,training=train_labels,iterations=model$interaction_training$iterations,
    normal_residual=model$interaction_training$normal_residual,external_features=ef,external_training_rows=ei)
}

select_penalty<-function(original,mfeatures,endpoint_name,excluded_gene=character(),held_cell=character()) {
  grid<-c(.003,.03,.3);cells<-sort(setdiff(unique(external$depmap_id),held_cell))
  groups<-rep(1:3,length.out=length(cells));predictions<-matrix(NA_real_,nrow(external),length(grid))
  for(fold in unique(groups)) {
    validation_cells<-cells[groups==fold]
    rows<-which(external$depmap_id %in% validation_cells & !external$gene %in% excluded_gene)
    for(k in seq_along(grid)) {
      fit<-fit_interaction(original,mfeatures,endpoint_name,excluded_gene,unique(c(held_cell,validation_cells)),grid[k])
      predictions[rows,k]<-predict_external_departure(fit$model,fit$external_features)[rows]
    }
  }
  rows<-which(!external$depmap_id %in% held_cell & !external$gene %in% excluded_gene & external$study=='elife32271')
  scores<-numeric(length(grid))
  if(length(rows)) {
    group<-paste(external$drug[rows],external$dose[rows],external$gene[rows],sep='|')
    geo<-center_geometry(predictions[rows,,drop=FALSE],external$observed_departure_log2[rows],group)
    ym<-geo$y-geo$yc;vc<-max(sum(geo$w*geo$yc^2),1e-8);vm<-max(sum(geo$w*ym^2),1e-8)
    for(k in seq_along(grid))scores[k]<-external_within_weight*sum(geo$w*(geo$yc-geo$Xc[,k])^2)/vc+
      (1-external_within_weight)*sum(geo$w*(ym-(geo$X[,k]-geo$Xc[,k]))^2)/vm
  }
  # Growth observations contribute their actual within-cell ligand order.
  growth<-which(!external$depmap_id %in% held_cell & !external$gene %in% excluded_gene & external$study=='Straussman2012')
  order_scores<-matrix(0,length(unique(paste(external$depmap_id[growth],external$drug[growth]))),length(grid))
  levels<-unique(paste(external$depmap_id[growth],external$drug[growth]))
  for(i in seq_along(levels)) {
    rr<-growth[paste(external$depmap_id[growth],external$drug[growth])==levels[i]]
    ix<-combn(rr,2);direction<-sign(external$observed_drug_arm_increment_log2[ix[1,]]-
      external$observed_drug_arm_increment_log2[ix[2,]])
    keep<-direction!=0
    for(k in seq_along(grid))order_scores[i,k]<-mean(direction[keep]*
      sign(predictions[ix[1,keep],k]-predictions[ix[2,keep],k])<0)
  }
  if(length(growth))scores<-scores+colMeans(order_scores)
  chosen<-grid[which.min(scores)]
  key<-paste(endpoint_name,paste(excluded_gene,collapse=','),paste(held_cell,collapse=','),sep='|')
  penalty_audit[[key]]<<-data.table(endpoint=endpoint_name,excluded_gene=paste(excluded_gene,collapse=','),
    held_cell=paste(held_cell,collapse=','),candidate_penalty=grid,internal_validation_loss=scores,
    selected=grid==chosen,validation='inner held-cell count variation/calibration plus within-cell growth order; no Wilson outcomes')
  chosen
}

predict_external_departure<-function(model,features) {
  products<-pair_features(features$Aviews,features$Bviews,model$context$W,model$interaction_basis)
  value<- -as.numeric(sweep(products,2,model$context$interaction_scale,'/') %*% model$par$Gamma)
  values<-as.numeric(tapply(value,factor(features$view_rows,levels=seq_len(nrow(external))),mean))
  main_effect<-features$standalone_views*model$ligand_single_conversion
  main_effect<-as.numeric(tapply(main_effect,factor(features$view_rows,levels=seq_len(nrow(external))),mean))
  # True counts at72h. For growth-ranking data the output is an order score,
  # not a claim of a calibrated growth-effect magnitude.
  ifelse(external$study=='elife32271',values*count_assay_days/reference_days,values+main_effect)
}

predict_mema<-function(model,features,rows) {
  products<-pair_features(features$A[rows,,drop=FALSE],features$B[rows,,drop=FALSE],model$context$W,model$interaction_basis)
  interaction<-as.numeric(sweep(products,2,model$context$interaction_scale,'/') %*% model$par$Gamma)
  data.table(CellLine=ligands$CellLine[rows],gene=ligands$gene[rows],
    predicted_log2= -features$S[rows]*count_assay_days/reference_days-
      model$mema_assay_scale*interaction*count_assay_days/reference_days,
    uncalibrated_drug_log2= -features$S[rows]*count_assay_days/reference_days,
    interaction_log2= -interaction*count_assay_days/reference_days,mema_assay_scale=model$mema_assay_scale)
}
