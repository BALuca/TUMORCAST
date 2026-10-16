#!/usr/bin/env Rscript
# Joint conditional-response calibration; the same three g terms and frozen f/RNA.
.libPaths(c("/work/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(ggplot2);library(digest)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
parent <- "/work/results/program_pilot/107_within_pair_departure_interaction/linear"
root <- "/work/results/program_pilot/109_joint_conditional_interaction/linear"
mema_root <- "/work/grant_rebuild_2026-10-04/transfer_solution"
challenge_genes <- c("HGF","NRG1")
primary_endpoint <- "normalized_cell_count"
source_weight <- .2
interaction_penalty <- .05
assay_penalty <- .01
programs <- sprintf("P%02d",1:40)
for(folder in c("data","models","figures"))dir.create(file.path(root,folder),recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,name)fwrite(x,file.path(root,"data",paste0(name,".tsv")),sep="\t")
pearson <- function(x,y)if(length(x)>=4 && sd(x)>1e-10 && sd(y)>1e-10)cor(x,y) else NA_real_
archive_path <- "/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds"
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
basis_reference <- readRDS("/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge/prepared.rds")
stopifnot(max(abs(archive$prep$W-basis_reference$W[rownames(archive$prep$W),programs]))<1e-12)
rm(archive,basis_reference)
bank_selection <- fread("/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv")
drugs <- bank_selection$drug
contexts <- fread(file.path(parent,"data/molecular_contexts.tsv"))
extra_contexts <- contexts[11:nrow(contexts)]
context_key <- function(drug,dose)paste(drug,sprintf("%.10g",dose),sep="|")
context_index <- function(drug,dose)match(context_key(drug,dose),context_key(contexts$drug,contexts$dose))

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
    list(score = score, both = both, activity = activity, interaction_activity = absolute, route = route)
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
        activity = node$activity * amplitude, interaction_activity = node$interaction_activity * amplitude, D = D)
}

predict_pair <- function(model, meta) {
    A <- predict_single(model, meta[, .(depmap_id, drug = drugA, dose = doseA)])
    B <- predict_single(model, meta[, .(depmap_id, drug = drugB, dose = doseB)])
    products <- sweep(pair_features(A$interaction_activity, B$interaction_activity), 2, model$context$interaction_scale, 
        "/")
    J <- as.numeric(products %*% model$par$Gamma)
    scale <- exp(model$par$assay)[match(meta$SCREENER, c("FF", "FG"))]
    list(sensitivity = (A$sensitivity + B$sensitivity + J) * scale, departure = J * scale)
}


# Existing aggregation: average isoform log ratios within ECM, then ECMs equally.
labels <- fread(file.path(mema_root,'all_gene_endpoints.tsv'))
ligands <- fread(file.path(mema_root,'recipient_ligand_program_scores.tsv'))
stopifnot(nrow(ligands)==56L,!anyDuplicated(ligands[,.(CellLine,gene)]),
          all(is.finite(as.matrix(ligands[,..programs]))),
          all(labels$ecm_n==48L),all(labels$endpoint %in% c('normalized_cell_count','raw_cell_count')))
# Source HGF labels are not joined into training; challenge scoring occurs after freeze.
training_labels <- labels[!gene %in% challenge_genes]
stopifnot(uniqueN(training_labels$gene)==26L,nrow(training_labels)==104L)
write_table(training_labels,'mema_nonchallenge_training_outcomes')
write_table(ligands,'mema_molecular_inputs')
pairs <- fread(file.path(parent,'data/combination_training_outcomes.tsv'))
folds <- fread(file.path(parent,'data/recipient_folds.tsv'))
stopifnot(nrow(pairs)==307L,nrow(folds)==62L,
          !any(pairs$depmap_id %in% ligands$depmap_id))
write_table(pairs,'combination_training_outcomes')
write_table(folds,'recipient_folds')
for(name in c('measured_predictions_frozen.tsv','baseline_predictions_frozen.tsv',
 'measured_decoder_correlations.tsv','baseline_bridge_correlations.tsv','selected_doses.tsv','molecular_contexts.tsv')) {
 source <- file.path(parent,'data',name)
 target <- file.path(root,'data',name)
 stopifnot(file.copy(source,target,overwrite=TRUE),identical(digest(source,file=TRUE),digest(target,file=TRUE)))
}

# Frozen single-effect evaluations at 0.5 and 5 uM interpolate in log dose to 0.75 uM.
# This interpolates functional evaluations and interaction activities, not gene displacements.
mema_features <- function(model) {
 cells <- ligands$depmap_id
 low <- predict_single(model,data.table(depmap_id=cells,drug='Lapatinib ditosylate',dose=.5))
 high <- predict_single(model,data.table(depmap_id=cells,drug='Lapatinib ditosylate',dose=5))
 weight <- log(.75/.5)/log(5/.5)
 L <- (1-weight)*low$interaction_activity+weight*high$interaction_activity
 S <- (1-weight)*low$sensitivity+weight*high$sensitivity
 D <- as.matrix(ligands[,..programs])
 ctx <- model$context
 Q <- cbind(1,sweep(sweep(sweep(D,2,ctx$fm),2,ctx$fs,'/') %*% ctx$fv,2,ctx$fscale,'/'))
 scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),
                    crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
 raw_scale <- exp(Q %*% scale_map)
 minimum <- matrix(apply(ctx$sigma[1:10,],2,min),nrow(D),40,byrow=TRUE)
 maximum <- matrix(apply(ctx$sigma[1:10,],2,max),nrow(D),40,byrow=TRUE)
 scale <- pmin(pmax(raw_scale,minimum),maximum)
 amplitude <- as.numeric(exp(Q %*% model$par$log_amplitude))
 H <- (D/scale) %*% ctx$W*amplitude
 list(products=pair_features(L,H),S=S)
}

center_geometry <- function(X,y,group) {
 levels <- unique(group)
 index <- match(group,levels)
 weight <- 1/as.numeric(table(group)[group])/length(levels)
 means <- vapply(seq_len(ncol(X)),function(k)
   as.numeric(tapply(X[,k],factor(group,levels=levels),mean)),numeric(length(levels)))
 ymean <- as.numeric(tapply(y,factor(group,levels=levels),mean))
 Xcentered <- X-means[index,,drop=FALSE]
 ycentered <- y-ymean[index]
 list(X=X,Xc=Xcentered,y=y,yc=ycentered,w=weight,
      variance=sum(weight*ycentered^2),group=group)
}

# Only Gamma and a positive MEMA assay scale change. g is one shared rule across both sources.
fit_interaction <- function(original,mfeatures,endpoint_name,excluded_gene=character(),held_cell=character()) {
 train_pairs <- pairs[!depmap_id %in% held_cell]
 A <- predict_single(original,train_pairs[,.(depmap_id,drug=drugA,dose=doseA)])
 B <- predict_single(original,train_pairs[,.(depmap_id,drug=drugB,dose=doseB)])
 Pn <- pair_features(A$interaction_activity,B$interaction_activity)
 nci_scale <- exp(original$par$assay)[match(train_pairs$SCREENER,c('FF','FG'))]
 nci_group <- paste(train_pairs$pair,train_pairs$SCREENER,train_pairs$doseA,train_pairs$doseB,sep='|')
 ng <- center_geometry(Pn*nci_scale,-train_pairs$I,nci_group)
 # Select endpoint explicitly to avoid ambiguity between a column and the argument.
 train_labels <- training_labels[training_labels$endpoint==endpoint_name & !training_labels$gene %in% excluded_gene]
 row <- match(paste(train_labels$CellLine,train_labels$gene),paste(ligands$CellLine,ligands$gene))
 stopifnot(!anyNA(row),!any(train_labels$gene %in% challenge_genes),
           !any(train_labels$gene %in% excluded_gene),!any(train_pairs$depmap_id %in% held_cell))
 mg <- center_geometry(mfeatures$products[row,,drop=FALSE],-train_labels$y,train_labels$CellLine)
 base <- mfeatures$S[row]
 # Equal sources; molecular feature scaling uses training within-group variation only.
 product_scale <- sqrt(.5*colSums(ng$Xc^2*ng$w)+.5*colSums(mg$Xc^2*mg$w))
 stopifnot(all(product_scale>1e-10),ng$variance>1e-10,mg$variance>1e-10)
 Xn <- sweep(ng$Xc,2,product_scale,'/')
 Xm <- sweep(mg$Xc,2,product_scale,'/')
 raw_m <- sweep(mg$X,2,product_scale,'/')
 # Profile the exact joint objective over ONE positive study conversion.
 # ALMANAC targets are measured departures; MEMA targets are absolute conditional
 # drug/vehicle responses. No receiver intercept or held-out outcome enters fitting.
 gram_n <- source_weight*crossprod(Xn,Xn*ng$w)/ng$variance
 rhs_n <- source_weight*as.numeric(crossprod(Xn,ng$w*ng$yc))/ng$variance
 gram_m <- source_weight*crossprod(raw_m,raw_m*mg$w)/mg$variance
 solve_at_scale <- function(log_scale) {
   s <- exp(log_scale)
   gram <- gram_n+s^2*gram_m+diag(interaction_penalty,3)
   target <- mg$y-s*base
   rhs <- rhs_n+source_weight*s*as.numeric(crossprod(raw_m,mg$w*target))/mg$variance
   gamma <- as.numeric(solve(gram,rhs))
   nci_residual <- ng$yc-as.numeric(Xn %*% gamma)
   mema_residual <- mg$y-s*(base+as.numeric(raw_m %*% gamma))
   loss <- source_weight*sum(ng$w*nci_residual^2)/ng$variance+
     source_weight*sum(mg$w*mema_residual^2)/mg$variance+
     interaction_penalty*sum(gamma^2)+assay_penalty*log_scale^2
   list(loss=loss,Gamma=gamma,gram=gram,rhs=rhs)
 }
 optimum <- optimize(function(z)solve_at_scale(z)$loss,interval=c(-4,4),tol=1e-8)
 log_scale <- optimum$minimum
 assay_scale <- exp(log_scale)
 solution <- solve_at_scale(log_scale)
 Gamma <- solution$Gamma
 gram <- solution$gram
 rhs <- solution$rhs
 converged <- is.finite(optimum$objective) && abs(log_scale)<3.999
 iteration <- 1L
 normal_residual <- max(abs(as.numeric(gram %*% Gamma)-rhs))
 stopifnot(normal_residual<1e-8,converged)
 model <- original
 model$par$Gamma <- Gamma
 model$context$interaction_scale <- product_scale
 model$context$interaction_supported <- rep(TRUE,3)
 model$mema_assay_scale <- assay_scale
 model$interaction_training <- list(method='ALMANAC pair-centered + absolute MEMA jointly calibrated shared ridge',
   endpoint=endpoint_name,challenge_genes=challenge_genes,training_genes=unique(train_labels$gene),
   excluded_gene=excluded_gene,held_cell=held_cell,source_weight=source_weight,
   interaction_penalty=interaction_penalty,assay_penalty=assay_penalty,
   iterations=iteration,converged=converged,normal_residual=normal_residual)
 stopifnot(identical(original$rna,model$rna),identical(original$template,model$template),
           identical(original$par$level,model$par$level),
           identical(original$par$log_amplitude,model$par$log_amplitude),
           identical(original$par$assay,model$par$assay))
 list(model=model,training=train_labels,iterations=iteration,normal_residual=normal_residual)
}

predict_mema <- function(model,features,rows) {
 interaction <- as.numeric(sweep(features$products[rows,,drop=FALSE],2,model$context$interaction_scale,'/') %*% model$par$Gamma)
 data.table(CellLine=ligands$CellLine[rows],gene=ligands$gene[rows],
   predicted_log2=-model$mema_assay_scale*(features$S[rows]+interaction),
   uncalibrated_drug_log2=-features$S[rows],interaction_log2=-interaction,
   mema_assay_scale=model$mema_assay_scale)
}

original_final <- qread(file.path(parent,'models/final.qs'))
features_final <- mema_features(original_final)
all_predictions <- nci_predictions <- audits <- coefficient_tables <- list()
final_models <- list()
for(endpoint in c('normalized_cell_count','raw_cell_count')) {
 genes <- sort(unique(training_labels$gene))
 for(job in seq_len(length(genes)+1L)) {
   held_gene <- if(job>length(genes))character() else genes[job]
   label <- if(length(held_gene))paste0('ligand_',held_gene) else 'final'
   fit <- fit_interaction(original_final,features_final,endpoint,held_gene)
   model <- fit$model
   rows <- if(length(held_gene))which(ligands$gene==held_gene) else seq_len(nrow(ligands))
   output <- predict_mema(model,features_final,rows)
   output[,`:=`(endpoint=endpoint,evaluation=if(length(held_gene))'held_gene' else 'fitted',fit=label)]
   if(!length(held_gene))output[gene %in% challenge_genes,evaluation:='challenge']
   all_predictions[[paste(endpoint,label)]] <- output
   audits[[paste(endpoint,label)]] <- data.table(endpoint,fit=label,excluded_gene=if(length(held_gene))held_gene else NA_character_,
      training_ligands=uniqueN(fit$training$gene),training_contexts=nrow(fit$training),
      challenge_leakage=sum(fit$training$gene %in% challenge_genes),iterations=fit$iterations,
      normal_residual=fit$normal_residual,mema_assay_scale=model$mema_assay_scale)
   if(!length(held_gene)) {
      final_models[[endpoint]] <- model
      qsave(model,file.path(root,'models',paste0('final_',endpoint,'.qs')))
   }
 }
 cat(endpoint,': 26 complete ligand exclusions plus final fit complete\n')
}
# Retain primary models under the established filenames; recompute ALMANAC whole-cell predictions.
for(job in seq_len(nrow(folds)+1L)) {
 held <- if(job>nrow(folds))character() else folds$depmap_id[job]
 label <- if(length(held))sprintf('loo_%02d',job) else 'final'
 original <- qread(file.path(parent,'models',paste0(label,'.qs')))
 if(length(held)) {
   features <- mema_features(original)
   fit <- fit_interaction(original,features,primary_endpoint,held_cell=held)
   model <- fit$model
 } else model <- final_models[[primary_endpoint]]
 qsave(model,file.path(root,'models',paste0(label,'.qs')))
 test <- if(length(held))pairs[depmap_id==held] else pairs
 if(nrow(test)) {
   predicted <- predict_pair(model,test)
   output <- test[,.(depmap_id,pair,SCREENER,drugA,drugB,doseA,doseB,
      observed_log2=y,observed_departure_log2=I)]
   output[,`:=`(predicted_log2=-predicted$sensitivity,predicted_departure_log2=-predicted$departure,
                evaluation=if(length(held))'held_cell' else 'fitted',fit=label)]
   nci_predictions[[label]] <- output
 }
 coefficient_tables[[label]] <- data.table(fit=label,term=c('signature1_signature1','cross_signature','signature2_signature2'),
    coefficient_sensitivity=model$par$Gamma,product_scale=model$context$interaction_scale,
    mema_assay_scale=model$mema_assay_scale)
 if(job %% 20L==0L)cat('ALMANAC cell exclusions completed:',job,'of',nrow(folds)+1L,'\n')
}
mema_pred <- rbindlist(all_predictions)
combo <- rbindlist(nci_predictions)
write_table(mema_pred,'mema_predictions_frozen')
write_table(combo,'combination_predictions_frozen')
write_table(rbindlist(audits),'fit_audit')
write_table(rbindlist(coefficient_tables),'shared_interaction_coefficients')
# Original v107 HGF vectors and standalone effects remain unchanged for the 5 uM query.
model <- final_models[[primary_endpoint]]
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
hmean <- as.numeric(hgf[context=='Equal-context mean',..programs])
ctx <- model$context
Qh <- cbind(1,sweep(matrix((hmean-ctx$fm)/ctx$fs,nrow=1) %*% ctx$fv,2,ctx$fscale,'/'))
scale_map <- solve(crossprod(ctx$Q[1:10,])/10+diag(c(0,1,1,1)),crossprod(ctx$Q[1:10,],log(ctx$sigma[1:10,]))/10)
hscale <- pmin(pmax(as.numeric(exp(Qh %*% scale_map)),apply(ctx$sigma[1:10,],2,min)),apply(ctx$sigma[1:10,],2,max))
lap_meta <- data.table(depmap_id='ACH-000248',drug='Lapatinib ditosylate',dose=5)
lap <- predict_single(model,lap_meta)
lap_pct <- 100*2^(-lap$sensitivity)
Hactivity <- (sweep(as.matrix(hgf[,..programs]),2,hscale,'/') %*% ctx$W)*as.numeric(exp(Qh %*% model$par$log_amplitude))
Lactivity <- lap$interaction_activity[rep(1L,nrow(hgf)),,drop=FALSE]
J <- as.numeric(sweep(pair_features(Lactivity,Hactivity),2,ctx$interaction_scale,'/') %*% model$par$Gamma)
hgf_query <- rbindlist(lapply(1:2,function(route)data.table(cell='AU565',dose=5,hgf_context=hgf$context,assumed_HGF_route=route,
 lapatinib_pct_control=lap_pct,hgf_increment_pct=100*2^(-Hactivity[,route]),
 combination_pct_control=lap_pct*2^(-Hactivity[,route]-J),
 lapatinib_pct_hgf_control=lap_pct*2^(-J),interaction_log2=-J,
 interaction_relative_pct=100*(2^(-J)-1),
 qualification='HGF/Wilson functional labels excluded; frozen absolute drug fit plus uncalibrated standalone HGF increment; g independent of HGF route.')))
write_table(hgf_query,'au565_hgf_5uM_predictions_frozen')
freeze_paths <- file.path(root,'data',c('mema_predictions_frozen.tsv','combination_predictions_frozen.tsv','au565_hgf_5uM_predictions_frozen.tsv'))
write_table(data.table(path=freeze_paths,sha256=vapply(freeze_paths,digest,character(1),file=TRUE,algo='sha256')),'prediction_hashes')
# Challenge and ordinary validation labels are joined only after prediction freeze.
scored <- merge(mema_pred,labels[,.(CellLine,gene,endpoint,observed_log2=y,ecm_n)],
                by=c('CellLine','gene','endpoint'),all.x=TRUE)
stopifnot(!anyNA(scored$observed_log2))
write_table(scored,'mema_predictions_scored')
metrics <- scored[,.(ligands=.N,r=pearson(observed_log2,predicted_log2),
 observed_median_pct=median(100*2^observed_log2),predicted_median_pct=median(100*2^predicted_log2)),
 by=.(endpoint,CellLine,evaluation)]
write_table(metrics,'mema_correlations')
challenge <- scored[evaluation=='challenge'][,.(endpoint,CellLine,gene,
 measured_pct_ligand_control=100*2^observed_log2,predicted_pct_ligand_control=100*2^predicted_log2,
 interaction_relative_pct=100*(2^interaction_log2-1),mema_assay_scale)]
write_table(challenge,'held_HGF_NRG1_predictions')
combo_metrics <- combo[,.(cells=.N,r_total=pearson(observed_log2,predicted_log2),
 r_departure=pearson(observed_departure_log2,predicted_departure_log2),
 observed_departure_SD=sd(observed_departure_log2),predicted_departure_SD=sd(predicted_departure_log2)),
 by=.(pair,SCREENER,evaluation)]
write_table(combo_metrics,'combination_correlations')
# Square, independently scaled panels: comparisons across ligands in one recipient at fixed 0.75 uM.
plot_data <- scored[evaluation=='held_gene']
plot_data[,panel:=paste(CellLine,ifelse(endpoint==primary_endpoint,'Normalized count','Raw count'),sep=' | ')]
plot_stats <- plot_data[,.(label=paste0('r = ',sprintf('%.2f',pearson(observed_log2,predicted_log2)),'\nn = ',.N,' ligands')),by=panel]
p <- ggplot(plot_data,aes(observed_log2,predicted_log2))+geom_point(color='#0072B2',size=1.5)+
 geom_smooth(method='lm',formula=y~x,se=FALSE,color='grey30',linewidth=.5)+
 geom_text(data=plot_stats,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
 facet_wrap(~panel,ncol=2,scales='free')+theme_bw(base_size=9)+theme(aspect.ratio=1)+
 labs(title='v109 | MEMA leave-one-ligand-out prediction',x='Measured log2 drug/vehicle under the same ligand',
 y='Predicted log2 drug/vehicle under the same ligand',
 caption='Lapatinib 0.75 µM, 72 h. Each point is one ligand averaged over ECM conditions.\nHGF and NRG1 excluded from fitting; all their isoforms/ECM conditions kept together. Retrospective evaluation.')
ggsave(file.path(root,'figures/mema_held_ligand_scatters.pdf'),p,width=7,height=7.1,device=cairo_pdf)
ggsave(file.path(root,'figures/mema_held_ligand_scatters.png'),p,width=7,height=7.1,dpi=160)
writeLines(c('v109: unchanged three-coefficient g; exact joint fit to ALMANAC departures and absolute MEMA conditional responses.',
 'Molecular maps, signatures, single-effect calibrations and NCI assay scales are frozen from v107.',
 'MEMA: lapatinib 0.75 uM, 72 h; normalized count primary, raw count mandatory sensitivity.',
 'Molecular and endpoint aggregation reused from prior preparation; 26 genes x two recipients train; HGF and NRG1 globally held out.',
 'MEMA endpoint is drug/vehicle under identical ligand and ECM. It is not an absolute no-ligand rescue endpoint.',
 'Drug functional evaluations/activities at 0.5 and 5 uM interpolate in log dose; no extrapolation or molecular-vector interpolation.',
 'g uses both existing signature activities for every constituent. No ligand identity or recipient intercept enters g.',
 'Center ALMANAC departures/products within pair; preserve absolute conditional MEMA responses and features.',
 'Equal source weights .2; ridge .05; source variance normalization. One positive MEMA study scale, log-scale penalty .01.',
 'Gamma and study conversion minimize one objective, profiled over a positive scale. No receiver/ligand intercepts. Fixed before scoring.',
 '26 leave-gene-out fits per endpoint. ALMANAC 62 whole-cell exclusions recalibrate g on nonchallenge MEMA and unheld ALMANAC only.',
 'Original HGF 12-context query retained for AU565 5 uM; ligand-growth prediction remains uncalibrated. Wilson labels not fitted.',
 'MEMA and Wilson outcomes have previously been exposed; these are retrospective development results.'),file.path(root,'RUN_DETAILS.txt'))
capture.output(sessionInfo(),file=file.path(root,'SESSION_INFO.txt'))
input_paths <- c('/work/src/program_pilot/109_joint_conditional_interaction.R',
 file.path(mema_root,c('all_gene_endpoints.tsv','recipient_ligand_program_scores.tsv')),
 file.path(parent,'models/final.qs'),archive_path)
write_table(data.table(path=input_paths,sha256=vapply(input_paths,digest,character(1),file=TRUE,algo='sha256')),'source_hashes')
print(metrics[evaluation=='held_gene'])
print(challenge)
print(hgf_query[hgf_context=='Equal-context mean'])
cat('Complete. Single-drug outputs byte-identical to v107; HGF/NRG1/Wilson labels never fitted.\n')
