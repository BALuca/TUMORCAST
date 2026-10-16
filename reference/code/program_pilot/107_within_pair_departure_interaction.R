#!/usr/bin/env Rscript
# Refit only shared g using within-combination departure variation. All f and RNA components frozen from v106.
.libPaths(c("/work/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
parent <- "/work/results/program_pilot/106_lapatinib_supported_dose_interactions/linear"
root <- "/work/results/program_pilot/107_within_pair_departure_interaction/linear"
interaction_penalty <- .05
departure_weight <- .2
programs <- sprintf("P%02d",1:40)
for(folder in c("data","models","figures"))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,"data",paste0(name,".tsv")),sep="\t")
pearson <- function(x,y)if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10)cor(x,y) else NA_real_
archive_path <- "/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds"
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
rm(archive)
bank_selection <- fread("/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv")
drugs <- bank_selection$drug
contexts <- fread(file.path(parent,"data/molecular_contexts.tsv"))
extra_contexts <- contexts[11:nrow(contexts)]
context_key <- function(drug,dose)paste(drug,sprintf("%.10g",dose),sep="|")
context_index <- function(drug,dose)match(context_key(drug,dose),context_key(contexts$drug,contexts$dose))
# Single-effect predictions remain frozen. Interaction inputs retain both signature projections.

predict_rna <- function(model, cells, drug, dose = NULL) {
    if (is.null(dose)) 
        dose <- bank_selection$dose[match(drug, drugs)]
    raw <- as.matrix(baseline[match(cells, depmap_id), ..programs])
    features <- sweep(sweep(sweep(raw, 2, model$bmean), 2, model$bsd, "/") %*% model$basis, 2, model$bscale, "/")
    g <- match(drug, drugs)
    D <- matrix(0, length(cells), 40)
    for (j in unique(g)) {
        rows <- which(g == j)
        Z <- sweep(features[rows, , drop = FALSE], 2, model$bm[j, ]) %*% model$B[, , j]
        D[rows, ] <- sweep(sweep(Z, 2, model$sigma[j, ], "*"), 2, model$mu[j, ], "+")
    }
    rows <- which(drug == "Lapatinib ditosylate")
    if (length(rows)) {
        lap <- model$lap_receptors
        stopifnot(!is.null(lap))
        X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
        Z <- sweep(sweep(X, 2, lap$xmean), 2, lap$xsd, "/") %*% lap$B
        D[rows, ] <- sweep(sweep(Z, 2, lap$sigma, "*"), 2, lap$mu, "+")
    }
    for (i in seq_len(nrow(extra_contexts))) {
        dr <- extra_contexts$drug[i]
        ds <- extra_contexts$dose[i]
        rows <- which(drug == dr & abs(dose - ds) < 1e-10)
        if (!length(rows)) 
            next
        map <- model$extra_rna[[context_key(dr, ds)]]
        stopifnot(!is.null(map))
        if (map$kind == "receptors") {
            X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
            X <- sweep(sweep(X, 2, map$xmean), 2, map$xsd, "/")
        }
        else X <- sweep(features[rows, , drop = FALSE], 2, map$xmean)
        D[rows, ] <- sweep(sweep(X %*% map$B, 2, map$sigma, "*"), 2, map$mu, "+")
    }
    stopifnot(all(is.finite(D)), !anyNA(context_index(drug, dose)))
    D
}

signature_node <- function(D, g, ctx) {
    X <- (D - ctx$mu[g, , drop = FALSE])/ctx$sigma[g, , drop = FALSE]
    both <- X %*% ctx$W
    route <- ctx$route[g]
    score <- both[cbind(seq_len(nrow(D)), route)]
    absolute <- (D/ctx$sigma[g, , drop = FALSE]) %*% ctx$W
    activity <- matrix(0, nrow(D), 2)
    activity[cbind(seq_len(nrow(D)), route)] <- absolute[cbind(seq_len(nrow(D)), route)]
    list(score = score, both = both, activity = activity,
         interaction_activity = absolute, route = route)
}

pair_features <- function(A, B) cbind(A[, 1] * B[, 1], A[, 1] * B[, 2] + A[, 2] * B[, 1], A[, 2] * B[, 2])

predict_single <- function(model, meta, measured_D = NULL) {
    dose <- if ("dose" %in% names(meta)) 
        meta$dose
    else bank_selection$dose[match(meta$drug, drugs)]
    D <- if (is.null(measured_D)) 
        predict_rna(model$rna, meta$depmap_id, meta$drug, dose)
    else measured_D
    g <- context_index(meta$drug, dose)
    node <- signature_node(D, g, model$context)
    Q <- model$context$Q[g, , drop = FALSE]
    level <- as.numeric(Q %*% model$par$level)
    amplitude <- as.numeric(exp(Q %*% model$par$log_amplitude))
    sensitivity <- level + amplitude * node$score
    assay <- if ("assay" %in% names(meta)) 
        match(meta$assay, c("secondary", "nci_FF", "nci_FG"))
    else rep(1L, nrow(meta))
    scale <- c(1, exp(model$par$assay))[assay]
    stopifnot(all(amplitude > 0), !anyNA(scale))
    list(sensitivity = sensitivity * scale, score = node$score, route = node$route, level = level, amplitude = amplitude, 
        activity = node$activity * amplitude,
        interaction_activity = node$interaction_activity * amplitude, D = D)
}

predict_pair <- function(model, meta) {
    A <- predict_single(model, meta[, .(depmap_id, drug = drugA, dose = doseA)])
    B <- predict_single(model, meta[, .(depmap_id, drug = drugB, dose = doseB)])
    products <- sweep(pair_features(A$interaction_activity, B$interaction_activity),
                      2, model$context$interaction_scale, "/")
    J <- as.numeric(products %*% model$par$Gamma)
    scale <- exp(model$par$assay)[match(meta$SCREENER, c("FF", "FG"))]
    list(sensitivity = (A$sensitivity + B$sensitivity + J) * scale, departure = J * scale)
}


pairs <- fread(file.path(parent,'data/combination_training_outcomes.tsv'))
folds <- fread(file.path(parent,'data/recipient_folds.tsv'))
stopifnot(nrow(pairs)==307L,uniqueN(pairs$pair)==8L,nrow(folds)==62L,
          !anyDuplicated(pairs[,.(depmap_id,pair,SCREENER,doseA,doseB)]),
          !'ACH-000248' %in% pairs$depmap_id)
# Preserve single-drug outputs exactly; no upstream, signature or single-effect fitting.
copy_names <- c('selected_doses.tsv','molecular_contexts.tsv','recipient_folds.tsv',
 'combination_training_outcomes.tsv','nci_single_training_outcomes.tsv',
 'measured_predictions_frozen.tsv','baseline_predictions_frozen.tsv',
 'measured_decoder_correlations.tsv','baseline_bridge_correlations.tsv',
 'lapatinib_comparison_with_v105.tsv','added_lapatinib_pair_coverage.tsv')
for(name in copy_names) {
 source <- file.path(parent,'data',name)
 destination <- file.path(root,'data',name)
 stopifnot(file.copy(source,destination,overwrite=TRUE),
           identical(digest(source,file=TRUE),digest(destination,file=TRUE)))
}
combination_predictions <- fit_tables <- audits <- list()
for(job in seq_len(nrow(folds)+1L)) {
 held <- if(job>nrow(folds))character() else folds$depmap_id[job]
 label <- if(length(held))sprintf('loo_%02d',job) else 'final'
 source_path <- file.path(parent,'models',paste0(label,'.qs'))
 original <- qread(source_path)
 model <- original
 stopifnot(identical(held,model$held_cells),!any(held %in% model$training_cells))
 train <- pairs[!depmap_id %in% held]
 A <- predict_single(model,train[,.(depmap_id,drug=drugA,dose=doseA)])
 B <- predict_single(model,train[,.(depmap_id,drug=drugB,dose=doseB)])
 # Fixed calibrated molecular activities; assay conversion also remains frozen.
 products <- pair_features(A$interaction_activity,B$interaction_activity)
 endpoint_scale <- exp(model$par$assay)[match(train$SCREENER,c('FF','FG'))]
 X <- products*endpoint_scale
 group <- paste(train$pair,train$SCREENER,train$doseA,train$doseB,sep='|')
 groups <- unique(group)
 group_index <- match(group,groups)
 weight <- 1/as.numeric(table(group)[group])/length(groups)
 Xmean <- vapply(1:3,function(k)as.numeric(tapply(X[,k],factor(group,levels=groups),mean)),numeric(length(groups)))
 Xcentered <- X-Xmean[group_index,,drop=FALSE]
 target <- -train$I
 target_mean <- as.numeric(tapply(target,factor(group,levels=groups),mean))
 target_centered <- target-target_mean[group_index]
 target_scale <- sqrt(sum(weight*target_centered^2))
 within_scale <- sqrt(colSums(Xcentered^2*weight))
 supported <- within_scale>1e-10
 stopifnot(target_scale>1e-10,all(supported))
 within_scale[!supported] <- 1
 Z <- sweep(Xcentered[,supported,drop=FALSE],2,within_scale[supported],'/')
 Y <- target_centered/target_scale
 # Original .2 departure weight and .05 L2 penalty, now on standardized within-pair variation.
 gram <- crossprod(Z,Z*weight)
 rhs <- as.numeric(crossprod(Z,weight*Y))
 ridge <- gram+diag(interaction_penalty/departure_weight,ncol(Z))
 theta <- as.numeric(solve(ridge,rhs))
 residual <- max(abs(as.numeric(ridge %*% theta)-rhs))
 stopifnot(residual<1e-10)
 Gamma <- rep(0,3)
 Gamma[supported] <- theta*target_scale
 model$par$Gamma <- Gamma
 # Prediction uses UNcentered products: no named-pair intercept or training-pair lookup.
 # Consequently every interaction remains exactly zero if either response is absent.
 model$context$interaction_scale <- within_scale
 model$context$interaction_supported <- supported
 model$interaction_training_cells <- unique(train$depmap_id)
 model$interaction_training <- list(method='within-pair departure-only ridge; both signature activities',
   feature_scale=within_scale,target_scale=target_scale,penalty=interaction_penalty,
   loss_weight=departure_weight,normal_equation_residual=residual,
   centering='training only; prediction remains zero-anchored uncentered product')
 stopifnot(identical(original$rna,model$rna),identical(original$template,model$template),
           identical(original$par$level,model$par$level),
           identical(original$par$log_amplitude,model$par$log_amplitude),
           identical(original$par$assay,model$par$assay),
           identical(A$sensitivity,predict_single(model,train[,.(depmap_id,drug=drugA,dose=doseA)])$sensitivity),
           !any(held %in% model$interaction_training_cells),
           max(abs(pair_features(matrix(0,nrow(train),2),B$interaction_activity)))==0)
 qsave(model,file.path(root,'models',paste0(label,'.qs')))
 test <- if(length(held))pairs[depmap_id==held] else pairs
 if(nrow(test)) {
   prediction <- predict_pair(model,test)
   output <- test[,.(depmap_id,pair,SCREENER,drugA,drugB,doseA,doseB,
     observed_log2=y,observed_departure_log2=I)]
   output[,`:=`(predicted_log2=-prediction$sensitivity,
     predicted_departure_log2=-prediction$departure,fold=job,
     evaluation=if(length(held))'held_cell' else 'fitted')]
   combination_predictions[[label]] <- output
 }
 fit_tables[[label]] <- data.table(fit=label,term=c('route1_route1','cross_route','route2_route2'),
   training_supported=supported,within_product_SD=within_scale,
   standardized_coefficient=replace(rep(0,3),which(supported),theta),
   coefficient_sensitivity=Gamma,target_within_SD=target_scale)
 audits[[label]] <- data.table(fit=label,held_cell=if(length(held))held else NA_character_,
   training_pairs=nrow(train),functional_leakage=length(intersect(held,model$interaction_training_cells)),
   unchanged_single_effects=TRUE,normal_equation_residual=residual)
 if(job %% 10L==0L || job>nrow(folds))cat('g fits completed:',job,'of',nrow(folds)+1L,'\n')
}
combo <- rbindlist(combination_predictions)
stopifnot(nrow(combo[evaluation=='held_cell'])==nrow(pairs))
write_table(combo,'combination_predictions_frozen')
write_table(rbindlist(fit_tables),'within_pair_interaction_coefficients')
write_table(rbindlist(audits),'fit_audit')
# HGF: same molecular vector, inferred scales, amplitudes and two conditional routes.
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
hmean <- as.numeric(hgf[context=='Equal-context mean',..programs])
ctx <- model$context
Qh <- cbind(1,sweep(matrix((hmean-ctx$fm)/ctx$fs,nrow=1) %*% ctx$fv,2,ctx$fscale,'/'))
scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
hscale_raw <- as.numeric(exp(Qh %*% scale_map))
hscale <- pmin(pmax(hscale_raw,apply(ctx$sigma[1:10,],2,min)),apply(ctx$sigma[1:10,],2,max))
lap <- predict_single(model,data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate'))
hgf_queries <- list()
for(route in 1:2) {
 hctx <- ctx
 hctx$mu <- rbind(ctx$mu,hmean)
 hctx$sigma <- rbind(ctx$sigma,hscale)
 hctx$route <- c(ctx$route,route)
 hnode <- signature_node(as.matrix(hgf[,..programs]),rep(nrow(ctx$Q)+1L,nrow(hgf)),hctx)
 hamplitude <- as.numeric(exp(Qh %*% model$par$log_amplitude))
 Sh <- hamplitude*hnode$activity[,route]
 lap_effect <- sum(lap$activity)
 HA <- hnode$interaction_activity*hamplitude
 LA <- lap$interaction_activity[rep(1L,nrow(hgf)),,drop=FALSE]
 J <- as.numeric(sweep(pair_features(LA,HA),2,ctx$interaction_scale,'/') %*% model$par$Gamma)
 hgf_queries[[route]] <- data.table(cell='AU565',dose=5,hgf_context=hgf$context,assumed_HGF_route=route,
   lapatinib_pct_zero_control=100*2^(-lap_effect),hgf_alone_pct_zero_control=100*2^(-Sh),
   combination_pct_zero_control=100*2^(-lap_effect-Sh-J),
   lapatinib_in_hgf_pct_hgf_control=100*2^(-lap_effect-J),interaction_log2=-J,
   interaction_relative_pct=100*(2^(-J)-1),hgf_addition_relative_pct=100*(2^(-Sh-J)-1),
   qualification='Frozen v106 single effects; g uses both signatures and needs no HGF route assignment. HGF standalone route unlearned; no HGF/Wilson labels fitted. Zero-reference query.')
}
query <- rbindlist(hgf_queries)
stopifnot(max(query[,.(spread=max(interaction_log2)-min(interaction_log2)),by=hgf_context]$spread)<1e-12)
write_table(query,'au565_hgf_5uM_zero_reference_query')
old_hgf <- fread(file.path(parent,'data/au565_hgf_5uM_zero_reference_summary.tsv'))
stopifnot(max(abs(query$lapatinib_pct_zero_control-old_hgf$lapatinib_pct_zero_control))<1e-10,
          max(abs(query$hgf_alone_pct_zero_control-old_hgf$hgf_alone_pct_zero_control))<1e-10)
paths <- file.path(root,'data',c('measured_predictions_frozen.tsv','baseline_predictions_frozen.tsv',
 'combination_predictions_frozen.tsv','au565_hgf_5uM_zero_reference_query.tsv'))
write_table(data.table(path=paths,sha256=vapply(paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')
# Evaluate only after predictions are frozen. All correlations are within fixed pair/dose.
metrics <- combo[,.(cells=.N,r_total=pearson(observed_log2,predicted_log2),
 r_departure=pearson(observed_departure_log2,predicted_departure_log2),
 observed_departure_SD=sd(observed_departure_log2),predicted_departure_SD=sd(predicted_departure_log2),
 observed_median_pct=median(100*2^observed_log2),predicted_median_pct=median(100*2^predicted_log2),
 observed_sensitive=sum(observed_log2<= -1),predicted_sensitive=sum(predicted_log2<= -1)),
 by=.(pair,SCREENER,evaluation)]
write_table(metrics,'combination_correlations')
old_metrics <- fread(file.path(parent,'data/combination_correlations.tsv'))[evaluation=='held_cell']
comparison <- merge(old_metrics,metrics[evaluation=='held_cell'],by=c('pair','SCREENER','evaluation'),suffixes=c('_v106','_v107'))
write_table(comparison,'combination_comparison_with_v106')
previous_metrics <- fread(file.path(root,'data/previous_assigned_route_107_correlations.tsv'))[evaluation=='held_cell']
previous_comparison <- merge(previous_metrics,metrics[evaluation=='held_cell'],
  by=c('pair','SCREENER','evaluation'),suffixes=c('_previous','_both_signatures'))
write_table(previous_comparison,'combination_comparison_with_previous_107')
previous_hgf <- fread(file.path(root,'data/previous_assigned_route_107_hgf.tsv'))
write_table(merge(previous_hgf,query,by=c('cell','dose','hgf_context','assumed_HGF_route'),
  suffixes=c('_previous','_both_signatures')),'hgf_comparison_with_previous_107')
old_predictions <- fread(file.path(parent,'data/combination_predictions_frozen.tsv'))[evaluation=='held_cell']
old_absolute <- old_predictions[,.(observed_median_pct=median(100*2^observed_log2),
  predicted_median_pct=median(100*2^predicted_log2),observed_sensitive=sum(observed_log2<= -1),
  predicted_sensitive=sum(predicted_log2<= -1),r=pearson(observed_log2,predicted_log2)),by=pair][,version:=106L]
new_absolute <- metrics[evaluation=='held_cell',.(pair,observed_median_pct,predicted_median_pct,
  observed_sensitive,predicted_sensitive,r=r_total)][,version:=107L]
write_table(rbind(old_absolute,new_absolute),'absolute_combination_calibration_comparison')
for(endpoint in c('total','departure')) {
 data <- copy(combo[evaluation=='held_cell'])
 data[,observed:=if(endpoint=='total')observed_log2 else observed_departure_log2]
 data[,predicted:=if(endpoint=='total')predicted_log2 else predicted_departure_log2]
 stats <- data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed,predicted)),'\nn = ',.N)),by=pair]
 plot <- ggplot(data,aes(observed,predicted))+geom_point(color='#0072B2',size=1.5)+
   geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
   geom_text(data=stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
   facet_wrap(~pair,ncol=4,scales='free')+theme_bw(base_size=9)+
   theme(aspect.ratio=1,panel.grid.minor=element_blank(),strip.text=element_text(size=8))+
   labs(title=paste('v107 | Held-cell combination',endpoint),
     subtitle='Both signature activities enter g; frozen single effects; within-pair departure training',
     x=paste('Measured',endpoint,'log2 abundance'),y=paste('Predicted',endpoint,'log2 abundance'),
     caption='Each panel compares cells at fixed combination doses. No new named-pair intercepts; zero-input interaction preserved.\nPreviously revealed outcomes and selected doses make evaluation retrospective.')
 ggsave(file.path(root,'figures',paste0('held_combination_',endpoint,'.pdf')),plot,width=12,height=7,device=cairo_pdf)
 ggsave(file.path(root,'figures',paste0('held_combination_',endpoint,'.png')),plot,width=12,height=7,dpi=150)
}
writeLines(c('v107: ONLY g refitted from frozen v106 whole-cell-excluded models; same 307 pair observations, eight pairs.',
 'All upstream predictors, two signatures, route assignments, single-effect calibrations and assay scales remain identical.',
 'Fit observed log2 departure I=yAB-yA-yB. Center response and product features by training pair/screen/doses.',
 'Equal pair weights; within-pair product SD scaling. Target standardized by pooled training within-pair SD.',
 'g uses both zero-referenced signature activities for every constituent: A1B1, A1B2+A2B1, A2B2. No assigned-route masking in g.',
 'Objective: .2 departure fit + .05 L2 penalty. Closed-form ridge; all three shared interaction coefficients supported.',
 'Prediction uses uncentered raw activity products. No pair-specific intercept; absent-input interaction remains exactly 0.',
 'Absolute interaction means are therefore not fitted: absolute combination calibration can worsen and is evaluated explicitly.',
 '63 fits: 62 whole-cell exclusions plus all-data fit. No held-cell interaction labels enter its fit.',
 'Original f was learned jointly with g in v106. Freezing it isolates this g training change; it is not newly trained independently.',
 'HGF molecule and both standalone route assumptions preserved. g is identical under both assignments. HGF/Wilson labels never fitted. Query 5 uM lapatinib.',
 'Source intervention/dose selection and previously revealed labels make comparisons retrospective.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
source_paths <- c('/work/src/program_pilot/107_within_pair_departure_interaction.R',
 file.path(parent,'data/combination_training_outcomes.tsv'),file.path(parent,'models/final.qs'),archive_path)
write_table(data.table(path=source_paths,sha256=vapply(source_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(comparison[,.(pair,r_departure_v106,r_departure_v107,r_total_v106,r_total_v107,
 predicted_departure_SD_v106,predicted_departure_SD_v107)])
print(query[hgf_context=='Equal-context mean'])
cat('Complete. Single-drug files byte-identical; molecular/f parameters unchanged in every fit.\n')
