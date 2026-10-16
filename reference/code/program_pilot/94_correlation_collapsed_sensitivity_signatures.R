#!/usr/bin/env Rscript
# Training-only correlation prototypes -> two fixed scores -> shrunk drug readouts.
# This tests the measured-displacement decoder before joint/combination learning.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(qs); library(glmnet); library(ggplot2); library(digest)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/94_correlation_collapsed_sensitivity_signatures'
previous <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
programs <- sprintf('P%02d', 1:40)
seed <- 94012L
lambda_grid <- c(.01, .1, 1, 10)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive=TRUE, showWarnings=FALSE)
}
write_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep='\t')
}
status <- function(message) {
  cat(message, '\n')
  flush.console()
  writeLines(message, file.path(root, 'STATUS.txt'))
}
pearson <- function(x, y) {
  if (length(x) < 4 || sd(x) < 1e-10 || sd(y) < 1e-10) return(NA_real_)
  cor(x, y)
}

selection <- fread(file.path(previous, 'data/selected_doses.tsv'))
folds <- fread(file.path(previous, 'data/recipient_folds.tsv'))
outcomes <- fread(file.path(previous, 'data/individual_training_outcomes.tsv'))[assay=='secondary']
drugs <- selection$drug
short <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
           'Lapatinib', 'Topotecan', 'Vandetanib', 'Vismodegib', 'Crizotinib')
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
molecular <- merge(as.data.table(archive$prep$m), selection, by=c('drug', 'dose'))[
  depmap_id!='ACH-000248']
profiles <- merge(molecular, outcomes[, .(depmap_id, drug, dose, sensitivity=-y)],
                  by=c('depmap_id', 'drug', 'dose'))
profiles[, fold:=folds$fold[match(depmap_id, folds$depmap_id)]]
stopifnot(nrow(profiles)==286L, length(drugs)==10L,
          !anyDuplicated(profiles[, .(depmap_id, drug)]), !anyNA(profiles$fold),
          !anyDuplicated(baseline$depmap_id), all(is.finite(as.matrix(profiles[, ..programs]))))
rm(archive)
write_table(selection, 'selected_doses')
write_table(folds, 'recipient_folds')
status('v94: 286 measured profiles; ten selected contexts; same three whole-cell folds.')

# Build two prototypes from training correlations, keeping the program weights fixed.
# Route 1 is always the drug cluster containing lapatinib, solely for stable labels.
learn_templates <- function(train) {
  mu <- sigma <- correlations <- matrix(NA_real_, 10, 40, dimnames=list(drugs, programs))
  reference <- numeric(10)
  counts <- integer(10)
  for (g in seq_along(drugs)) {
    rows <- which(train$drug==drugs[g])
    stopifnot(length(rows)>=5L)
    D <- as.matrix(train[rows, ..programs])
    mu[g, ] <- colMeans(D)
    sigma[g, ] <- apply(D, 2, sd)
    correlations[g, ] <- as.numeric(cor(D, train$sensitivity[rows]))
    reference[g] <- mean(train$sensitivity[rows])
    counts[g] <- length(rows)
  }
  stopifnot(all(is.finite(correlations)), all(sigma>1e-10))
  tree <- hclust(dist(correlations), method='ward.D2')
  cluster <- cutree(tree, 2)
  lap_cluster <- cluster[match('Lapatinib ditosylate', drugs)]
  membership <- ifelse(cluster==lap_cluster, 1L, 2L)
  collapsed <- vapply(1:2, function(k) colMeans(correlations[membership==k, , drop=FALSE]), numeric(40))
  dimnames(collapsed) <- list(programs, c('Signature 1', 'Signature 2'))
  norms <- sqrt(colSums(collapsed^2))
  stopifnot(all(norms>1e-10))
  W <- sweep(collapsed, 2, norms, '/')
  g <- match(train$drug, drugs)
  Z <- (as.matrix(train[, ..programs])-mu[g, , drop=FALSE])/sigma[g, , drop=FALSE]
  raw_scores <- Z %*% W
  score_sd <- t(vapply(1:10, function(j) apply(raw_scores[g==j, , drop=FALSE], 2, sd), numeric(2)))
  stopifnot(all(score_sd>1e-10))
  list(mu=mu, sigma=sigma, correlations=t(correlations), collapsed=collapsed, W=W,
       score_sd=score_sd, reference=reference, counts=counts, membership=membership,
       drug_tree=tree, training_cells=unique(train$depmap_id))
}

# Signature score calibration uses training RNA variance, without functional shifts.
scores <- function(template, D, drug) {
  g <- match(drug, drugs)
  stopifnot(!anyNA(g), ncol(D)==40L)
  Z <- (D-template$mu[g, , drop=FALSE])/template$sigma[g, , drop=FALSE]
  (Z %*% template$W)/template$score_sd[g, , drop=FALSE]
}

# Global two-coefficient readout + shrunk two-coefficient deviations per drug.
# No intercept: the fitted target is sensitivity relative to its training mean.
fit_readout <- function(template, train, lambda) {
  S <- scores(template, as.matrix(train[, ..programs]), train$drug)
  g <- match(train$drug, drugs)
  target <- train$sensitivity-template$reference[g]
  design <- matrix(0, nrow(train), 22)
  design[, 1:2] <- S
  for (j in 1:10) design[g==j, 2+(2*j-1):(2*j)] <- S[g==j, , drop=FALSE]
  weights <- 1/(10*template$counts[g])
  penalty <- diag(c(rep(.001, 2), rep(lambda/10, 20)))
  beta <- solve(crossprod(design, design*weights)+penalty, crossprod(design, target*weights))
  coefficients <- t(matrix(beta[-(1:2)], 2, 10))
  coefficients <- sweep(coefficients, 2, beta[1:2], '+')
  list(coefficients=coefficients, shared=as.numeric(beta[1:2]), lambda=lambda)
}
predict_readout <- function(template, readout, D, drug) {
  S <- scores(template, D, drug)
  g <- match(drug, drugs)
  rowSums(S*readout$coefficients[g, , drop=FALSE])
}

# Inner whole-cell folds select only the shrinkage strength. Rebuild prototypes
# in every inner split: inner validation labels never enter its correlation map.
tune_readout <- function(train, job) {
  set.seed(seed+job)
  cells <- sample(sort(unique(train$depmap_id)))
  inner_fold <- setNames(rep(1:3, length.out=length(cells)), cells)
  losses <- list()
  for (k in 1:3) {
    held <- cells[inner_fold==k]
    inner_train <- train[!depmap_id %in% held]
    validation <- train[depmap_id %in% held]
    template <- learn_templates(inner_train)
    stopifnot(!any(held %in% template$training_cells))
    for (lambda in lambda_grid) {
      readout <- fit_readout(template, inner_train, lambda)
      predicted <- predict_readout(template, readout, as.matrix(validation[, ..programs]), validation$drug)
      error <- predicted-(validation$sensitivity-template$reference[match(validation$drug, drugs)])
      loss <- data.table(drug=validation$drug, error)[, .(loss=mean(error^2)), by=drug]
      losses[[length(losses)+1L]] <- data.table(inner_fold=k, lambda=lambda, loss=mean(loss$loss))
    }
  }
  audit <- rbindlist(losses)
  average <- audit[, .(loss=mean(loss)), by=lambda]
  best <- average[order(loss, -lambda), lambda][1]
  list(lambda=best, audit=audit)
}

# RNA-only baseline -> displacement mapping, fixed independently of viability.
# This secondary replay distinguishes decoder failure from displacement prediction.
fit_rna <- function(train_rna) {
  ids <- sort(unique(train_rna$depmap_id))
  raw <- as.matrix(baseline[match(ids, depmap_id), ..programs])
  bmean <- colMeans(raw)
  bsd <- apply(raw, 2, sd)
  normalized <- sweep(sweep(raw, 2, bmean), 2, bsd, '/')
  basis <- svd(normalized, nu=0, nv=8)$v
  bscale <- sqrt(colMeans((normalized %*% basis)^2))
  B <- array(0, c(8, 40, 10))
  mu <- sigma <- matrix(0, 10, 40)
  bm <- matrix(0, 10, 8)
  for (g in 1:10) {
    rows <- which(train_rna$drug==drugs[g])
    D <- as.matrix(train_rna[rows, ..programs])
    mu[g, ] <- colMeans(D)
    sigma[g, ] <- apply(D, 2, sd)
    features <- sweep(normalized[match(train_rna$depmap_id[rows], ids), , drop=FALSE] %*% basis, 2, bscale, '/')
    bm[g, ] <- colMeans(features)
    features <- sweep(features, 2, bm[g, ])
    target <- sweep(sweep(D, 2, mu[g, ]), 2, sigma[g, ], '/')
    fit <- glmnet(features, target, family='mgaussian', alpha=0, lambda=1,
                  intercept=FALSE, standardize=FALSE)
    B[, , g] <- vapply(coef(fit), function(x) as.numeric(x)[-1], numeric(8))
  }
  list(bmean=bmean, bsd=bsd, basis=basis, bscale=bscale, B=B,
       mu=mu, sigma=sigma, bm=bm, training_cells=ids)
}
predict_rna <- function(model, cells, drug) {
  raw <- as.matrix(baseline[match(cells, depmap_id), ..programs])
  stopifnot(all(is.finite(raw)))
  features <- sweep(sweep(sweep(raw, 2, model$bmean), 2, model$bsd, '/') %*% model$basis,
                    2, model$bscale, '/')
  g <- match(drug, drugs)
  D <- matrix(0, length(cells), 40)
  for (j in unique(g)) {
    rows <- which(g==j)
    Z <- sweep(features[rows, , drop=FALSE], 2, model$bm[j, ]) %*% model$B[, , j]
    D[rows, ] <- sweep(sweep(Z, 2, model$sigma[j, ], '*'), 2, model$mu[j, ], '+')
  }
  D
}

# Every fold shows its own observed training map and its two collapsed columns.
render_heatmap <- function(template, label) {
  observed <- template$correlations
  ro <- hclust(dist(observed), method='ward.D2')$order
  co <- template$drug_tree$order
  annotation <- data.frame(Signature=paste0('Signature ', template$membership[co]),
    Dose=paste0(selection$dose[co], ' µM'), Cells=as.character(template$counts[co]), row.names=drugs[co])
  left <- draw_heatmap(observed[ro, co, drop=FALSE], scale_rows=FALSE,
    value_range=c(-1, 0, 1), value_colors=c('#2166AC', 'white', '#B2182B'),
    hm_name='Training Pearson r', col_anno_df=annotation, col_anno_cols=c('Signature', 'Dose', 'Cells'),
    palette_overrides=list(Signature=c('Signature 1'='#D58B36', 'Signature 2'='#775DA6'),
      Dose=c('0.05 µM'='#DCEAF2', '0.5 µM'='#6BAED6', '5 µM'='#08519C')),
    show_row_names=TRUE, row_names_side='left', row_names_gp=grid::gpar(fontsize=8),
    show_column_names=TRUE, column_labels=short[co], column_names_rot=45,
    column_names_gp=grid::gpar(fontsize=8), cluster_rows=FALSE, cluster_columns=FALSE,
    width=grid::unit(105, 'mm'), title='Observed training-cell correlations',
    column_title_gp=grid::gpar(fontsize=11))
  right <- draw_heatmap(template$collapsed[ro, , drop=FALSE], scale_rows=FALSE,
    value_range=c(-1, 0, 1), value_colors=c('#2166AC', 'white', '#B2182B'),
    hm_name='Collapsed mean r', show_row_names=FALSE, show_column_names=TRUE,
    column_names_rot=45, column_names_gp=grid::gpar(fontsize=8),
    cluster_rows=FALSE, cluster_columns=FALSE, width=grid::unit(24, 'mm'),
    title='Two fixed signatures', column_title_gp=grid::gpar(fontsize=11))
  drawn <- ComplexHeatmap::draw(left+right, heatmap_legend_side='bottom', annotation_legend_side='right',
    ht_gap=grid::unit(8, 'mm'), column_title=paste0('v94 | ', label, ' | training data only'),
    column_title_gp=grid::gpar(fontsize=13), padding=grid::unit(c(15, 5, 5, 5), 'mm'))
  stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn, name='Training Pearson r')), 1:40))
  grid::grid.text('Positive = greater sensitivity (-log2 abundance). Right: equal-drug mean correlations within each cluster.\nWeights are fixed after L2 normalization; each drug/signature score uses its training RNA SD. No row z-scoring.',
    x=.5, y=grid::unit(4, 'mm'), just=c('centre', 'bottom'), gp=grid::gpar(fontsize=8))
}

measured_predictions <- baseline_predictions <- coefficient_tables <- membership_tables <- list()
template_tables <- tuning_tables <- audits <- list()
book_path <- file.path(root, 'figures/training_correlations_and_collapsed_signatures_per_fold.pdf')
cairo_pdf(book_path, width=9, height=9, onefile=TRUE)
book_device <- dev.cur()
for (job in 1:4) {
  label <- if (job==4) 'final' else paste0('fold_', job)
  held <- if (job==4) character() else folds[fold==job, depmap_id]
  train <- profiles[!depmap_id %in% held]
  test <- if (job==4) profiles else profiles[depmap_id %in% held]
  status(paste(label, ': deriving training-only templates and choosing shrinkage in inner cell folds.'))
  template <- learn_templates(train)
  tuning <- tune_readout(train, job)
  readout <- fit_readout(template, train, tuning$lambda)
  rna <- fit_rna(molecular[!depmap_id %in% held])
  stopifnot(!any(held %in% template$training_cells), !any(held %in% rna$training_cells))
  model <- list(template=template, readout=readout, rna=rna, held_cells=held, selected=selection)
  qsave(model, file.path(root, 'models', paste0(label, '.qs')))
  evaluation <- if (job==4) 'fitted' else 'held_cell'
  predicted <- predict_readout(template, readout, as.matrix(test[, ..programs]), test$drug)
  p <- test[, .(depmap_id, drug, dose, observed_sensitivity=sensitivity)]
  p[, observed_deviation:=observed_sensitivity-template$reference[match(drug, drugs)]]
  p[, `:=`(predicted_deviation=predicted, fold=job, evaluation=evaluation)]

  # Individual-drug ridge is a diagnostic comparator on exactly the same cells.
  comparator <- numeric(nrow(test))
  for (g in 1:10) {
    tr <- which(train$drug==drugs[g])
    te <- which(test$drug==drugs[g])
    Ztr <- sweep(sweep(as.matrix(train[tr, ..programs]), 2, template$mu[g, ]), 2, template$sigma[g, ], '/')
    Zte <- sweep(sweep(as.matrix(test[te, ..programs]), 2, template$mu[g, ]), 2, template$sigma[g, ], '/')
    fit <- glmnet(Ztr, train$sensitivity[tr]-template$reference[g], alpha=0, lambda=.1,
                  intercept=FALSE, standardize=FALSE)
    comparator[te] <- as.numeric(predict(fit, Zte, s=.1))
  }
  p[, individual_ridge_deviation:=comparator]
  measured_predictions[[label]] <- p
  all_test <- if (job==4) outcomes else outcomes[depmap_id %in% held]
  Dhat <- predict_rna(rna, all_test$depmap_id, all_test$drug)
  bp <- all_test[, .(depmap_id, drug, dose, observed_sensitivity=-y)]
  bp[, observed_deviation:=observed_sensitivity-template$reference[match(drug, drugs)]]
  bp[, `:=`(predicted_deviation=predict_readout(template, readout, Dhat, all_test$drug),
            fold=job, evaluation=evaluation)]
  baseline_predictions[[label]] <- bp
  coefficient_tables[[label]] <- data.table(fit=label, drug=drugs, dose=selection$dose,
    signature_1_beta=readout$coefficients[, 1], signature_2_beta=readout$coefficients[, 2],
    shared_1=readout$shared[1], shared_2=readout$shared[2], lambda=tuning$lambda)
  membership_tables[[label]] <- data.table(fit=label, drug=drugs, signature=template$membership,
                                          training_cells=template$counts)
  template_tables[[label]] <- data.table(fit=label, program=programs,
    collapsed_r_1=template$collapsed[, 1], collapsed_r_2=template$collapsed[, 2],
    weight_1=template$W[, 1], weight_2=template$W[, 2])
  tuning_tables[[label]] <- copy(tuning$audit)[, fit:=label]
  audits[[label]] <- data.table(fit=label, training_cells=uniqueN(train$depmap_id),
    held_measured_cells=if (job==4) 0L else uniqueN(test$depmap_id),
    held_cell_leakage=length(intersect(held, template$training_cells)), lambda=tuning$lambda)
  dev.set(book_device)
  render_heatmap(template, label)
  cairo_pdf(file.path(root, 'figures', paste0(label, '_training_correlations_and_signatures.pdf')), width=9, height=9)
  render_heatmap(template, label)
  dev.off()
  png(file.path(root, 'figures', paste0(label, '_training_correlations_and_signatures.png')),
      width=9, height=9, units='in', res=150, type='cairo')
  render_heatmap(template, label)
  dev.off()
}
dev.set(book_device)
dev.off()
measured <- rbindlist(measured_predictions)
replay <- rbindlist(baseline_predictions)
write_table(measured, 'measured_predictions_frozen')
write_table(replay, 'baseline_predictions_frozen')
write_table(rbindlist(coefficient_tables), 'drug_signature_coefficients')
write_table(rbindlist(membership_tables), 'drug_cluster_memberships')
write_table(rbindlist(template_tables), 'collapsed_signature_weights')
write_table(rbindlist(tuning_tables), 'inner_shrinkage_selection')
write_table(rbindlist(audits), 'fold_audit')
prediction_paths <- file.path(root, 'data', c('measured_predictions_frozen.tsv', 'baseline_predictions_frozen.tsv'))
write_table(data.table(path=prediction_paths, sha256=vapply(prediction_paths, digest, character(1), file=TRUE,
                                                          algo='sha256')), 'prediction_hashes')

# Evaluate only after writing predictions. Fold-centered correlations remove
# differences between fold-specific reference levels; both versions are retained.
for (table in list(measured, replay)) {
  table[, `:=`(observed_centered=observed_sensitivity-mean(observed_sensitivity),
               predicted_centered=predicted_deviation-mean(predicted_deviation)), by=.(drug, evaluation, fold)]
}
measured[, ridge_centered:=individual_ridge_deviation-mean(individual_ridge_deviation), by=.(drug, evaluation, fold)]
metrics <- measured[, .(cells=.N, r_within_fold=pearson(observed_centered, predicted_centered),
  r_pooled_deviation=pearson(observed_deviation, predicted_deviation),
  r_individual_ridge_within_fold=pearson(observed_centered, ridge_centered)), by=.(drug, dose, evaluation)]
baseline_metrics <- replay[, .(cells=.N, r_within_fold=pearson(observed_centered, predicted_centered),
  r_pooled_deviation=pearson(observed_deviation, predicted_deviation)), by=.(drug, dose, evaluation)]
write_table(metrics, 'measured_decoder_correlations')
write_table(baseline_metrics, 'baseline_bridge_correlations')
write_table(measured[evaluation=='held_cell', .(cells=.N, r=pearson(observed_sensitivity, predicted_deviation)),
                     by=.(drug, dose, fold)], 'per_fold_measured_correlations')

plot_predictions <- function(table, endpoint, stem) {
  data <- copy(table[evaluation=='held_cell'])
  data[, panel:=factor(paste0(short[match(drug, drugs)], ' | ', dose, ' µM'),
    levels=paste0(short, ' | ', selection$dose, ' µM'))]
  stats <- data[, .(label=paste0('r within fold = ', sprintf('%.2f', pearson(observed_centered, predicted_centered)),
                                '\nn = ', .N)), by=panel]
  plot <- ggplot(data, aes(observed_centered, predicted_centered))+
    geom_hline(yintercept=0, color='grey80', linewidth=.3)+
    geom_point(aes(color=factor(fold)), size=1.5, alpha=.8)+
    geom_smooth(method='lm', formula=y~x, se=FALSE, color='grey30', linewidth=.5)+
    geom_text(data=stats, aes(x=-Inf, y=Inf, label=label), inherit.aes=FALSE, hjust=-.05, vjust=1.15, size=2.8)+
    facet_wrap(~panel, ncol=5, scales='free')+
    scale_color_manual(values=c('#0072B2', '#D55E00', '#009E73'), name='Held fold')+
    theme_bw(base_size=9)+theme(aspect.ratio=1, legend.position='bottom', panel.grid.minor=element_blank())+
    labs(title=paste0('v94 | ', endpoint), x='Measured sensitivity (centered within held fold)',
         y='Predicted sensitivity deviation (centered within held fold)',
         caption='Sensitivity = -log2 PRISM abundance. No drug-average prediction or absolute-viability calibration.\nPearson r compares cells receiving the same drug/dose after centering within each held fold.')
  ggsave(file.path(root, 'figures', paste0(stem, '.pdf')), plot, width=13, height=6.8, device=cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(stem, '.png')), plot, width=13, height=6.8, dpi=150)
}
plot_predictions(measured, 'held-cell decoder supplied with measured displacements', 'held_measured_decoder')
plot_predictions(replay, 'held-cell baseline → RNA-only displacement → fixed signatures', 'held_baseline_bridge')
writeLines(c('v94 tests correlation-collapsed signatures on the ten v93 selected drug/dose contexts.',
  'All three outer whole-cell folds are inherited. AU565 is excluded from RNA and functional training.',
  'Training-only signed Pearson maps -> Ward.D2 two drug clusters -> equal-drug mean correlation templates.',
  'Template vectors are L2 normalized and frozen. Scores use measured training drug/program mean and SD, then training drug/score SD.',
  'Only two signed functional coefficients per drug are fitted, shrunk to shared coefficients. No functional intercept or program-wise learned shifts/gains.',
  'Target is within-drug sensitivity deviation, not absolute viability. Drug mean is a reference for observed deviations only, never added to predictions.',
  'Lambda .01/.1/1/10 is selected by inner whole-cell folds, rebuilding templates in each inner training split.',
  'Individual-drug 40-program ridge lambda .1 is a fixed diagnostic comparator, not part of the main predictor.',
  'Secondary baseline replay fits RNA-only 8-axis baseline -> 40 displacement ridge lambda 1, using outer-training RNA only.',
  'Only measured RNA/PRISM matched cells train the functional decoder. Other PRISM cells are secondary evaluation; no ALMANAC or combination refitting here.',
  'Primary statistic: pooled within-held-fold Pearson across cells within each fixed drug/dose. Pooled deviation correlations also saved.',
  'Dose selection and previous test labels were already exposed; these are retrospective comparisons, not independent confirmation.'),
  file.path(root, 'RUN_DETAILS.txt'))
capture.output(sessionInfo(), file=file.path(root, 'SESSION_INFO.txt'))
inputs <- c(archive_path, file.path(previous, 'data', c('selected_doses.tsv', 'recipient_folds.tsv',
                                                     'individual_training_outcomes.tsv')),
            '/work/src/program_pilot/94_correlation_collapsed_sensitivity_signatures.R')
write_table(data.table(path=inputs, sha256=vapply(inputs, digest, character(1), file=TRUE, algo='sha256')), 'source_hashes')
print(metrics[evaluation=='held_cell'])
print(baseline_metrics[evaluation=='held_cell'])
status('Complete: fold-specific correlation maps and collapsed signatures, held-cell measured decoder and baseline replay.')
