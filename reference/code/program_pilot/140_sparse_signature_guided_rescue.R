#!/usr/bin/env Rscript
# Sparse shared rescue weights on signed contributions to a fixed sensitivity route.
# g(A,B) = a_drug * max(S_route' A, 0) * sum_p theta_p (-S_route,p * B_p).
# theta >= 0; raw molecular inputs retain their zero and program scales.
# HGF outcomes never select weights, regularization, routes or calibration.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(ggplot2); library(digest)})
RhpcBLASctl::blas_set_num_threads(1L)
setDTthreads(1L)
root <- '/work/results/program_pilot/140_sparse_signature_guided_rescue/linear'
source_root <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
programs <- sprintf('P%02d', 1:40)
penalty_grid <- c(.005, .02, .08, .25, .6, 1)
elastic_alpha <- .8
within_weight <- .8
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive=TRUE, showWarnings=FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep='\t')
pearson <- function(x,y) if (length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10) cor(x,y) else NA_real_

parent <- qread(file.path(source_root, 'models/final.qs'))
outcomes <- fread(file.path(source_root, 'data/external_combination_outcomes.tsv'))
ligands <- fread(file.path(source_root, 'data/external_ligand_programs.tsv'))
drug_molecular <- fread(file.path(source_root, 'data/external_predicted_drug_programs.tsv'))
profiles <- fread(file.path(source_root, 'data/time_matched_HGF_program_inputs.tsv'))
S <- parent$template$collapsed
stopifnot(identical(rownames(S), programs), ncol(S)==2)
# One scalar per entire signature fixes its units without changing relative program weights.
S <- sweep(S, 2, sqrt(colSums(S^2)), '/')
outcomes[, observation_id := .I]
training <- outcomes[study=='Straussman2012' & gene!='HGF']
challenge <- outcomes[study=='Straussman2012' & gene=='HGF']
calibration <- outcomes[study=='elife32271']
stopifnot(nrow(training)==144, nrow(challenge)==12, nrow(calibration)==16,
          all(calibration$gene %in% c('NRG1','EGF')), uniqueN(training$gene)==12)
ligand_key <- function(x) paste(x$depmap_id,x$gene,sep='|')
drug_key <- function(x) paste(x$depmap_id,x$drug,x$dose,sep='|')
drug_molecular <- unique(drug_molecular[,c('depmap_id','drug','dose',programs),with=FALSE])
stopifnot(!anyDuplicated(ligand_key(ligands)), !anyDuplicated(drug_key(drug_molecular)))

# Native lapatinib route is frozen. RAF/MEK have no PRISM-derived route assignment;
# assign their strongest positive mean molecular projection, without functional labels.
assign_routes <- function(meta) {
  result <- c('Lapatinib ditosylate'=unname(parent$template$membership['Lapatinib ditosylate']))
  for (drug_name in c('PD-184352','PLX-4720')) {
    keys <- unique(drug_key(meta[drug==drug_name]))
    ix <- match(keys, drug_key(drug_molecular))
    stopifnot(length(ix)>0, !anyNA(ix))
    scores <- as.numeric(colMeans(as.matrix(drug_molecular[ix,..programs])) %*% S)
    stopifnot(max(scores)>0)
    result[drug_name] <- which.max(scores)
  }
  result
}
make_features <- function(meta, routes, ligand_matrix=NULL, drug_matrix=NULL) {
  li <- match(ligand_key(meta), ligand_key(ligands))
  di <- match(drug_key(meta), drug_key(drug_molecular))
  if (is.null(ligand_matrix)) {
    stopifnot(!anyNA(li))
    ligand_matrix <- as.matrix(ligands[li,..programs])
  }
  if (is.null(drug_matrix)) {
    stopifnot(!anyNA(di))
    drug_matrix <- as.matrix(drug_molecular[di,..programs])
  }
  route <- as.integer(routes[meta$drug])
  stopifnot(!anyNA(route), all(is.finite(ligand_matrix)), all(is.finite(drug_matrix)),
            nrow(ligand_matrix)==nrow(meta), nrow(drug_matrix)==nrow(meta))
  weights <- t(S[,route,drop=FALSE])
  activity <- pmax(rowSums(drug_matrix*weights),0)
  design <- -ligand_matrix*weights*activity
  colnames(design) <- programs
  list(X=design, activity=activity, route=route, drug=drug_matrix, ligand=ligand_matrix)
}

# Fit weighted within-cell/drug ligand contrasts (80%) and group means (20%).
# This is an optimization transformation only; prediction uses the uncentered inputs.
prepare_geometry <- function(meta, design) {
  group <- paste(meta$depmap_id,meta$drug,meta$dose,sep='|')
  w <- 1/as.numeric(table(group)[group])/length(unique(group))
  XM <- matrix(0,nrow(design),ncol(design))
  ym <- numeric(nrow(design))
  y <- meta$observed_departure_log2
  for (g in unique(group)) {
    ii <- which(group==g)
    XM[ii,] <- matrix(colMeans(design[ii,,drop=FALSE]),length(ii),ncol(design),byrow=TRUE)
    ym[ii] <- mean(y[ii])
  }
  xc <- design-XM
  yc <- y-ym
  vc <- sum(w*yc^2)
  vm <- sum(w*ym^2)
  stopifnot(vc>1e-10,vm>1e-10)
  A <- rbind(xc*sqrt(within_weight*w/vc),XM*sqrt((1-within_weight)*w/vm))
  ystar <- c(yc*sqrt(within_weight*w/vc),ym*sqrt((1-within_weight)*w/vm))
  lambda_unit <- max(as.numeric(crossprod(A,ystar)),0)/elastic_alpha
  stopifnot(lambda_unit>0)
  list(A=A,y=ystar,drug=rep(meta$drug,2),lambda_unit=lambda_unit)
}

fit_signature <- function(rows, penalty) {
  meta <- training[rows]
  routes <- assign_routes(meta)
  features <- make_features(meta,routes)
  geometry <- prepare_geometry(meta,features$X)
  lambda <- penalty*geometry$lambda_unit
  at_scale <- function(log_a) {
    a <- exp(log_a)
    multiplier <- ifelse(geometry$drug=='PLX-4720',a,1)
    A <- geometry$A*multiplier
    # Penalize effective weights on both drugs so their scalar cannot evade shrinkage.
    l1 <- lambda*elastic_alpha*(1+a)/2
    l2 <- lambda*(1-elastic_alpha)*(1+a^2)/2
    # Express the ridge term as extra rows so Gaussian response normalization in
    # glmnet cannot change the intended raw-unit L1/L2 penalty balance.
    augmented_A <- rbind(A,diag(sqrt(l2),40))
    augmented_y <- c(geometry$y,rep(0,40))
    fit <- glmnet(augmented_A,augmented_y,alpha=1,lambda=l1/nrow(augmented_A),
                  intercept=FALSE,standardize=FALSE,lower.limits=0,thresh=1e-13,maxit=100000)
    theta <- as.numeric(coef(fit))[-1]
    residual <- as.numeric(A %*% theta)-geometry$y
    objective <- sum(residual^2)/2+l1*sum(theta)+l2*sum(theta^2)/2
    gradient <- as.numeric(crossprod(A,residual))+l1+l2*theta
    active <- theta>1e-7
    kkt <- max(c(abs(gradient[active]),pmax(-gradient[!active],0),0))
    list(theta=theta,objective=objective,kkt=kkt)
  }
  brackets <- seq(-4,4,by=1)
  values <- vapply(brackets,function(a)at_scale(a)$objective,numeric(1))
  j <- which.min(values)
  interval <- brackets[c(max(1,j-1),min(length(brackets),j+1))]
  opt <- optimize(function(a)at_scale(a)$objective,interval=interval,tol=1e-5)
  solution <- at_scale(opt$minimum)
  # A null model has no identifiable relative calibration; retain the reference value.
  if (max(solution$theta)<1e-10) {
    opt$minimum <- 0
    solution <- at_scale(0)
  }
  stopifnot(all(solution$theta>=0),solution$kkt<1e-4)
  list(theta=setNames(solution$theta,programs),routes=routes,
       multipliers=c('PD-184352'=1,'PLX-4720'=exp(opt$minimum)),
       penalty=penalty,lambda=lambda,training_rows=rows,kkt=solution$kkt,
       scale_boundary=abs(opt$minimum)>3.99,signature=S)
}
predict_signature <- function(model, meta, ligand_matrix=NULL, drug_matrix=NULL) {
  features <- make_features(meta,model$routes,ligand_matrix,drug_matrix)
  multiplier <- model$multipliers[meta$drug]
  stopifnot(!anyNA(multiplier))
  as.numeric(features$X %*% model$theta)*as.numeric(multiplier)
}
selection_audit <- list()
select_and_fit <- function(held_cells=character(),held_genes=character(),label='final') {
  rows <- which(!training$depmap_id %in% held_cells & !training$gene %in% held_genes)
  genes <- sort(unique(training$gene[rows]))
  fold <- setNames(rep(1:3,length.out=length(genes)),genes)
  pred <- matrix(NA_real_,length(rows),length(penalty_grid))
  for (k in 1:3) {
    valid <- which(fold[training$gene[rows]]==k)
    for (j in seq_along(penalty_grid)) {
      model <- fit_signature(rows[-valid],penalty_grid[j])
      pred[valid,j] <- predict_signature(model,training[rows[valid]])
    }
  }
  stopifnot(all(is.finite(pred)))
  group <- paste(training$depmap_id[rows],training$drug[rows],sep='|')
  y <- training$observed_departure_log2[rows]
  ym <- ave(y,group,FUN=mean)
  yc <- y-ym
  pm <- apply(pred,2,function(z)ave(z,group,FUN=mean))
  pc <- pred-pm
  w <- 1/as.numeric(table(group)[group])/length(unique(group))
  loss <- within_weight*colSums((pc-yc)^2*w)/sum(w*yc^2)+
          (1-within_weight)*colSums((pm-ym)^2*w)/sum(w*ym^2)
  chosen <- which.min(loss)
  selection_audit[[label]] <<- data.table(fit=label,penalty=penalty_grid,
    internal_selection_loss=loss,selected=seq_along(loss)==chosen,
    excluded_cells=paste(held_cells,collapse=','),excluded_ligands=paste(c('HGF',held_genes),collapse=','))
  fit_signature(rows,penalty_grid[chosen])
}

predictions <- list()
models <- list()
coefficients <- list()
for (gene in sort(unique(training$gene))) {
  label <- paste0('held_ligand_',gene)
  model <- select_and_fit(held_genes=gene,label=label)
  rows <- which(training$gene==gene)
  predictions[[label]] <- training[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=predict_signature(model,training[rows]),evaluation='held_ligand')]
  models[[label]] <- model
  cat('Completed',label,'active programs:',sum(model$theta>1e-7),'\n')
}
for (cell in sort(unique(training$depmap_id))) {
  label <- paste0('held_cell_',cell)
  model <- select_and_fit(held_cells=cell,label=label)
  rows <- which(training$depmap_id==cell)
  predictions[[label]] <- training[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=predict_signature(model,training[rows]),evaluation='held_cell')]
  models[[label]] <- model
  cat('Completed',label,'active programs:',sum(model$theta>1e-7),'\n')
}
model <- select_and_fit()
models[['final']] <- model
for (label in names(models)) {
  m <- models[[label]]
  qsave(m,file.path(root,'models',paste0(label,'.qs')))
  coefficients[[label]] <- data.table(fit=label,program=programs,theta=m$theta)
}
coef_table <- rbindlist(coefficients)
write_table(coef_table,'program_weights_all_fits')
write_table(rbindlist(selection_audit),'training_only_penalty_selection')
write_table(rbindlist(lapply(names(models),function(label) {
  m <- models[[label]]
  data.table(fit=label,penalty=m$penalty,active_programs=sum(m$theta>1e-7),
    RAF_multiplier=m$multipliers[['PLX-4720']],scale_boundary=m$scale_boundary,kkt=m$kkt,
    MEK_route=m$routes[['PD-184352']],RAF_route=m$routes[['PLX-4720']])
})),'fit_diagnostics')

# Frozen program weights; one nonnegative lapatinib calibration from NRG1/EGF only.
cal_features <- make_features(calibration,model$routes)
cal_score <- as.numeric(cal_features$X %*% model$theta)
calibrate <- function(rows) {
  group <- paste(calibration$gene[rows],calibration$dose[rows],sep='|')
  w <- 1/as.numeric(table(group)[group])/length(unique(group))
  denominator <- sum(w*cal_score[rows]^2)
  multiplier <- if (denominator<=1e-12) 0 else
    max(0,sum(w*cal_score[rows]*calibration$observed_departure_log2[rows])/denominator)
  list(multiplier=multiplier,training_rows=rows,identified=denominator>1e-12)
}
lap_scale <- calibrate(seq_len(nrow(calibration)))
cal_predictions <- list()
for (cell in unique(calibration$depmap_id)) {
  rows <- which(calibration$depmap_id==cell)
  scale <- calibrate(setdiff(seq_len(nrow(calibration)),rows))
  cal_predictions[[cell]] <- calibration[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=scale$multiplier*cal_score[rows],evaluation='held_cell_calibration')]
}
for (gene in unique(calibration$gene)) {
  rows <- which(calibration$gene==gene)
  scale <- calibrate(setdiff(seq_len(nrow(calibration)),rows))
  cal_predictions[[gene]] <- calibration[rows,.(observation_id,CellLine,drug,dose,gene)][,
    `:=`(predicted_departure_log2=scale$multiplier*cal_score[rows],evaluation='held_ligand_calibration')]
}
qsave(list(signature=model,lapatinib=lap_scale),file.path(root,'models/final_transfer.qs'))

# Load only the existing, unchanged single-drug prediction definitions.
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
contexts <- fread(file.path(source_root,'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
required <- c('context_key','context_index','predict_rna','signature_node','predict_single')
for (e in parse('/work/src/program_pilot/138_pooled_direction_constrained_rescue.R')) {
  if (is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% required) eval(e)
}
# The historical single-drug helpers use these bank names for context lookup.
drugs <- bank_selection$drug
query_rows <- list()
query_contributions <- list()
for (dose in c(.5,5)) {
  meta <- data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=dose)
  lap <- predict_single(parent,meta)
  au <- profiles[CellLine=='AU565']
  D <- lap$D[rep(1,nrow(au)),,drop=FALSE]
  B <- as.matrix(au[,..programs])
  feature <- make_features(meta[rep(1,nrow(au))],model$routes,ligand_matrix=B,drug_matrix=D)
  contribution <- sweep(feature$X,2,model$theta,'*')*lap_scale$multiplier
  interaction <- rowSums(contribution)
  query_rows[[as.character(dose)]] <- data.table(cell='AU565',dose=dose,source_profile=au$source_profile,
    lapatinib_pct_control=100*2^(-lap$sensitivity),HGF_alone_pct_control=100,
    drug_activity=feature$activity,interaction_log2=interaction,
    combination_pct_control=100*2^(-lap$sensitivity+interaction))
  query_contributions[[as.character(dose)]] <- data.table(dose=dose,source_profile=au$source_profile,contribution)
}
query <- rbindlist(query_rows)
challenge_predictions <- challenge[,.(observation_id,CellLine,drug,dose,gene)][,
  `:=`(predicted_departure_log2=predict_signature(model,challenge),evaluation='all_HGF_outcomes_excluded')]
write_table(query,'AU565_HGF_predictions_frozen')
write_table(rbindlist(query_contributions),'AU565_HGF_program_contributions')
write_table(rbindlist(predictions),'RAF_MEK_predictions_frozen')
write_table(rbindlist(cal_predictions),'lapatinib_predictions_frozen')
write_table(challenge_predictions,'HGF_source_predictions_frozen')
write_table(training,'RAF_MEK_training_outcomes')
write_table(calibration,'lapatinib_calibration_outcomes')
write_table(data.table(multiplier=lap_scale$multiplier,identified=lap_scale$identified,
  observations=nrow(calibration),ligands='NRG1,EGF',HGF_outcomes=0),'lapatinib_calibration')
write_table(data.table(program=programs,S),'fixed_sensitivity_patterns')
write_table(cbind(training[,.(observation_id,CellLine,drug,dose,gene)],
  as.data.table(make_features(training,model$routes)$X)),'signature_guided_training_features')
freeze <- file.path(root,'data',c('AU565_HGF_predictions_frozen.tsv','RAF_MEK_predictions_frozen.tsv',
  'lapatinib_predictions_frozen.tsv','HGF_source_predictions_frozen.tsv'))
write_table(data.table(path=freeze,sha256=vapply(freeze,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')

scored <- merge(rbindlist(predictions),training[,.(observation_id,observed_departure_log2)],by='observation_id')
cs <- merge(rbindlist(cal_predictions),calibration[,.(observation_id,observed_departure_log2)],by='observation_id')
hs <- merge(challenge_predictions,challenge[,.(observation_id,observed_departure_log2)],by='observation_id')
metrics <- scored[,.(ligands=.N,r=pearson(observed_departure_log2,predicted_departure_log2)),
  by=.(evaluation,CellLine,drug,dose)]
cal_metrics <- cs[,.(cells=.N,r=pearson(observed_departure_log2,predicted_departure_log2),
  measured_median=median(observed_departure_log2),predicted_median=median(predicted_departure_log2)),by=.(evaluation,gene,dose)]
write_table(scored,'RAF_MEK_predictions_scored')
write_table(cs,'lapatinib_predictions_scored')
write_table(hs,'HGF_source_predictions_scored')
write_table(metrics,'across_ligand_held_correlations')
write_table(cal_metrics,'lapatinib_held_calibration_correlations')
write_table(coef_table[fit!='final',.(selected_folds=sum(theta>1e-7),folds=.N,median_theta=median(theta)),by=program],
  'program_selection_stability')

# Structural checks: zero inputs, assigned-route direction, actual fold exclusions.
sample_meta <- training[1:3]
zeros <- matrix(0,3,40)
stopifnot(all(predict_signature(model,sample_meta,ligand_matrix=zeros)==0),
          all(predict_signature(model,sample_meta,drug_matrix=zeros)==0))
test_ligands <- -t(S[,model$routes[sample_meta$drug],drop=FALSE])
stopifnot(all(predict_signature(model,sample_meta,ligand_matrix=test_ligands)>=0),
          all(predict_signature(model,sample_meta,ligand_matrix=-test_ligands)<=0))
for (label in names(models)) {
  fit <- models[[label]]
  stopifnot(!any(training$gene[fit$training_rows]=='HGF'))
  if (grepl('held_ligand_',label)) stopifnot(!any(training$gene[fit$training_rows]==sub('held_ligand_','',label)))
  if (grepl('held_cell_',label)) stopifnot(!any(training$depmap_id[fit$training_rows]==sub('held_cell_','',label)))
}
write_table(data.table(check=c('zero_drug_zero_g','zero_ligand_zero_g','opposition_has_rescue_sign',
  'reinforcement_has_harm_sign','all_HGF_outcomes_excluded','outer_fold_exclusions','nonnegative_sparse_weights'),passed=TRUE),
  'structural_validation')

# Compact result figures and per-fold weight heatmap.
bars <- query[source_profile=='Primary_study_mean']
bars <- rbind(bars[,.(dose,treatment='Lapatinib',abundance=lapatinib_pct_control)],
              bars[,.(dose,treatment='+ HGF',abundance=combination_pct_control)])
bars[,treatment:=factor(treatment,levels=c('Lapatinib','+ HGF'))]
p <- ggplot(bars,aes(treatment,abundance,fill=treatment))+geom_col(width=.65)+
  geom_text(aes(label=sprintf('%.1f',abundance)),vjust=-.4,size=3)+
  facet_wrap(~dose,labeller=label_both,scales='free_y')+
  scale_fill_manual(values=c('Lapatinib'='#777777','+ HGF'='#D55E00'),guide='none')+
  scale_y_continuous(expand=expansion(mult=c(0,.15)))+theme_classic(base_size=9)+
  labs(x=NULL,y='Predicted abundance (% control)',title='AU565 | lapatinib dose in µM')
ggsave(file.path(root,'figures/AU565_HGF_predictions.pdf'),p,width=4.4,height=2.6,device=cairo_pdf)
ggsave(file.path(root,'figures/AU565_HGF_predictions.png'),p,width=4.4,height=2.6,dpi=170)
plot_data <- scored[evaluation=='held_ligand']
plot_data[,panel:=paste(CellLine,drug,sep=' | ')]
ann <- plot_data[,.(label=sprintf('r = %.2f',pearson(observed_departure_log2,predicted_departure_log2))),by=panel]
p <- ggplot(plot_data,aes(observed_departure_log2,predicted_departure_log2))+
  geom_point(size=1.6,colour='#555555')+geom_smooth(method='lm',formula=y~x,se=FALSE,colour='#0072B2',linewidth=.4)+
  geom_text(data=ann,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=2.5)+
  facet_wrap(~panel,ncol=4,scales='free')+theme_bw(base_size=8)+theme(aspect.ratio=1)+
  labs(x='Measured growth departure (log2)',y='Predicted growth departure (log2)',
       title='v140 | held-out ligands, HGF excluded from every fit',
       caption='Across 12 ligands per fixed cell/drug. Frozen molecular inputs and single-drug patterns; retrospective evaluation.')
ggsave(file.path(root,'figures/held_ligand_scatters.pdf'),p,width=9,height=7.2,device=cairo_pdf)
ggsave(file.path(root,'figures/held_ligand_scatters.png'),p,width=9,height=7.2,dpi=150)
source('/work/src/R/lib_heatmap.R')
wide <- dcast(coef_table,program~fit,value.var='theta')
B <- as.matrix(wide[match(programs,program),-1])
rownames(B) <- programs
ro <- hclust(dist(B),method='ward.D2')$order
co <- hclust(dist(t(B)),method='ward.D2')$order
ht <- draw_heatmap(B[ro,co,drop=FALSE],scale_rows=FALSE,
  value_range=c(0,max(B)/2,max(B)),value_colors=c('white','#FDAE61','#A50026'),
  hm_name='Rescue weight',show_row_names=TRUE,show_column_names=TRUE,
  row_names_gp=grid::gpar(fontsize=7),column_names_gp=grid::gpar(fontsize=6),column_names_rot=60,
  height=grid::unit(130,'mm'),width=grid::unit(140,'mm'),title='Shared nonnegative program weights | final and held-out fits')
cairo_pdf(file.path(root,'figures/rescue_weights_by_fold.pdf'),width=8.8,height=8)
ComplexHeatmap::draw(ht,heatmap_legend_side='right')
dev.off()

writeLines(c('v140: fixed sensitivity patterns -> signed ligand program contributions -> shared sparse nonnegative rescue weights.',
  'g = drug scalar * positive raw drug projection on its assigned route * sum(theta * negative signed ligand contributions).',
  'The two routes are never added. No per-program or per-ligand molecular SD scaling; no intercept.',
  'Lapatinib uses its frozen route. RAF/MEK routes use positive mean molecular projection without functional outcomes.',
  'Primary signature fitting: 144 RAF/MEK combinations, 12 non-HGF ligands, 6 cells. HGF excluded throughout.',
  'Elastic-net alpha=0.8; penalty chosen in inner whole-ligand folds. Nonnegative theta shared across drugs.',
  'MEK scalar fixed at 1; positive RAF scalar profiled with penalties on both drugs effective weights.',
  'Frozen theta then one nonnegative lapatinib scalar calibrated on 16 NRG1/EGF observations at 0.5 and 5 uM.',
  'Zero drug or zero ligand molecular input gives zero interaction. Baseline enters molecular inputs only.',
  'Outer splits hold out functional cells/ligands; the existing molecular maps and single-drug patterns stay frozen.',
  'Existing single-drug predictions retained. HGF-alone abundance remains 100% under the parent zero standalone conversion.',
  'This is a typed drug-ligand pilot. ALMANAC is not used; no new general drug-drug model is claimed.',
  'All previously viewed outcomes are retrospective; AU565 queries are not dose/time-matched Wilson validation.'),
  file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
paths <- c('/work/src/program_pilot/140_sparse_signature_guided_rescue.R',
  file.path(source_root,'models/final.qs'),file.path(source_root,'data/external_combination_outcomes.tsv'),
  file.path(source_root,'data/external_ligand_programs.tsv'),file.path(source_root,'data/external_predicted_drug_programs.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(query[source_profile=='Primary_study_mean'])
print(coef_table[fit=='final' & theta>1e-7][order(-theta)])
print(metrics[,.(median_r=median(r,na.rm=TRUE),min_r=min(r,na.rm=TRUE),max_r=max(r,na.rm=TRUE)),by=evaluation])
print(cal_metrics)
cat('v140 complete.\n')
