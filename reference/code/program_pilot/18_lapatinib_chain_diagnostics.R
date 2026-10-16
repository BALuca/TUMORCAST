#!/usr/bin/env Rscript
# Diagnose both links with the existing five frozen, training-only program sets.
# No clustering or gene screening is rerun. All model selection excludes outer test cells.
# Scores are unweighted means within the fixed hard gene clusters, at baseline and after treatment.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
})
data.table::setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(2L)
RhpcBLASctl::omp_set_num_threads(2L)
source_dir <- '/work/results/program_pilot/17_lapatinib_baseline_response_sensitivity'
out <- '/work/results/program_pilot/18_lapatinib_chain_diagnostics'
data_dir <- file.path(out, 'data')
figure_dir <- file.path(out, 'figures')
model_dir <- file.path(out, 'models')
for (path in c(out, data_dir, figure_dir, model_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(status = 'preparing hard cluster means', completed_folds = 0L),
  file.path(out, 'diagnostic_status.tsv'), sep = '\t')
# Preserve the previous soft-score summaries for the requested like-for-like comparison.
for (name in c('diagnostic_functional_metrics.tsv', 'diagnostic_molecular_metrics.tsv')) {
  old <- file.path(data_dir, name)
  saved <- file.path(data_dir, paste0('soft_', name))
  if (file.exists(old) && !file.exists(saved)) stopifnot(file.copy(old, saved))
}
lambda_grid <- c(10, 1, .1, .01)
labels_path <- '/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge/joint40/prism_screen_audit/preferred_screen_functional_labels.tsv'
expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
fold_paths <- file.path(source_dir, 'models', sprintf('fold_%02d.qs', 1:5))
stopifnot(all(file.exists(fold_paths)))
labels <- fread(labels_path)[drug == 'Lapatinib ditosylate' & dose == 2.5,
  .(depmap_id, outer_fold = fold, observed = corrected_y)]
folds <- fread(file.path(source_dir, 'data', 'recipient_folds.tsv'))
fold_map <- setNames(folds$fold, folds$depmap_id)
banks <- lapply(fold_paths, qread)
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
# Recover only the lapatinib gene displacement, retaining all existing normalization and recipients.
gene_response <- do.call(rbind, lapply(panel$input_paths, function(path) {
  molecular <- qread(path)
  row <- which(molecular$metadata$drug == 'lapatinib' & molecular$metadata$dose == 5)
  stopifnot(length(row) == 1L)
  molecular$response[row, , drop = FALSE]
}))
rownames(gene_response) <- panel$cell_map$depmap_id
genes <- sort(unique(c(unlist(lapply(banks, `[[`, 'genes')), 'ERBB2', 'MET')))
headers <- names(fread(expression_path, nrows = 0))
symbols <- sub(' \\(ENSG.*', '', headers[-1L])
ix <- which(symbols %chin% genes)
tab <- fread(expression_path, select = c(1L, ix + 1L))
expression <- as.matrix(tab[, -1L])
rownames(expression) <- tab[[1L]]
colnames(expression) <- symbols[ix]
stopifnot(!anyDuplicated(colnames(expression)), all(genes %chin% colnames(expression)), all(is.finite(expression)))
rm(tab)
paths <- c(fold_paths, labels_path, expression_path, panel_path, panel$input_paths,
  '/work/src/program_pilot/18_lapatinib_chain_diagnostics.R')
hashes <- data.table(path = paths, sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(hashes, file.path(data_dir, 'diagnostic_input_hashes.tsv'), sep = '\t')

# Penalty selection is internal to each outer training set. These inner scores
# condition on its frozen program basis; only the outer predictions are evidence.
fit_ridge <- function(x, y, ids, intercept = TRUE) {
  multi <- is.matrix(y)
  family <- if (multi) 'mgaussian' else 'gaussian'
  center <- if (multi) colMeans(y) else 0
  scale <- if (multi) apply(y, 2L, sd) else 1
  stopifnot(all(is.finite(scale)), all(scale > 0))
  target <- if (multi) sweep(sweep(y, 2L, center, '-'), 2L, scale, '/') else y
  inner <- as.integer(factor(fold_map[ids]))
  stopifnot(!anyNA(inner), length(unique(inner)) >= 3)
  cv <- cv.glmnet(x, target, family = family, alpha = 0, lambda = lambda_grid,
    foldid = inner, standardize = TRUE, intercept = intercept, thresh = 1e-8,
    parallel = FALSE, grouped = FALSE)
  list(fit = cv$glmnet.fit, lambda = cv$lambda.min, center = center, scale = scale,
    multi = multi, inner_mse = cv$cvm, training_ids = ids)
}
predict_ridge <- function(fit, x) {
  raw <- predict(fit$fit, newx = x, s = fit$lambda)
  if (!fit$multi) return(as.numeric(raw))
  z <- matrix(raw, nrow(x), length(fit$center))
  sweep(sweep(z, 2L, fit$scale, '*'), 2L, fit$center, '+')
}
cor_safe <- function(x, y) if (length(x) >= 3L && sd(x) > 0 && sd(y) > 0) cor(x, y) else NA_real_
design <- function(d, b, kind) {
  if (kind == 'linear') return(d)
  if (kind == 'quadratic') return(cbind(d, d^2))
  stopifnot(kind == 'conditioned')
  cbind(d, d^2, d * b)
}
molecular_predictions <- functional_predictions <- penalty_rows <- support_rows <- list()
score_checks <- list()
for (k in 1:5) {
  message('Stage diagnostics: fold ', k, '/5')
  bank <- banks[[k]]
  ids <- labels$depmap_id
  stopifnot(identical(names(bank$partition), bank$genes), all(bank$genes %chin% colnames(gene_response)))
  cluster_size <- tabulate(bank$partition, nbins = 40L)
  stopifnot(all(cluster_size > 0L))
  weights <- matrix(0, length(bank$genes), 40L,
    dimnames = list(bank$genes, sprintf('P%02d', 1:40)))
  weights[cbind(seq_along(bank$genes), bank$partition)] <- 1 / cluster_size[bank$partition]
  stopifnot(all(rowSums(weights > 0) == 1L), max(abs(colSums(weights) - 1)) < 1e-10)
  baseline <- expression[, bank$genes, drop = FALSE] %*% weights
  measured <- gene_response[, bank$genes, drop = FALSE] %*% weights
  # Verify the matrix projection against direct hard-cluster arithmetic means.
  direct_means <- sapply(1:40, function(j) rowMeans(gene_response[, bank$genes[bank$partition == j], drop = FALSE]))
  stopifnot(max(abs(measured - direct_means)) < 1e-10)
  response_correlation <- cor(measured[bank$discovery_ids, , drop = FALSE])
  eigenvalues <- pmax(eigen(response_correlation, symmetric = TRUE, only.values = TRUE)$values, 0)
  score_checks[[k]] <- data.table(fold = k, scoring = 'hard_cluster_mean',
    median_program_correlation = median(response_correlation[upper.tri(response_correlation)]),
    dimensions_for_90_percent = which(cumsum(eigenvalues) / sum(eigenvalues) >= .9)[1],
    maximum_direct_mean_error = max(abs(measured - direct_means)))
  molecular_ids <- intersect(rownames(measured), rownames(expression))
  train_m <- intersect(bank$discovery_ids, molecular_ids)
  test_m <- intersect(bank$excluded_ids, molecular_ids)
  train_f <- labels[outer_fold != k, depmap_id]
  test_f <- labels[outer_fold == k, depmap_id]
  stopifnot(!length(intersect(train_m, bank$excluded_ids)), !length(intersect(train_f, bank$excluded_ids)))
  y <- setNames(labels$observed, labels$depmap_id)
  bmu <- colMeans(baseline[train_f, , drop = FALSE])
  bsd <- apply(baseline[train_f, , drop = FALSE], 2L, sd)
  b <- sweep(sweep(baseline, 2L, bmu, '-'), 2L, bsd, '/')
  features <- list(receptors = expression[, c('ERBB2', 'MET'), drop = FALSE],
    baseline_programs = baseline,
    baseline_programs_receptors = cbind(baseline, expression[, c('ERBB2', 'MET'), drop = FALSE]))
  mean_response <- colMeans(measured[train_m, , drop = FALSE])
  sd_response <- apply(measured[train_m, , drop = FALSE], 2L, sd)
  fit_store <- list()
  add_function <- function(cohort, model, query_ids, prediction) {
    functional_predictions[[length(functional_predictions) + 1L]] <<- data.table(fold = k,
      cohort = cohort, model = model, depmap_id = query_ids, observed = y[query_ids], predicted = prediction)
  }
  for (feature_name in names(features)) {
    x <- features[[feature_name]]
    fit <- fit_ridge(x[train_m, , drop = FALSE], measured[train_m, , drop = FALSE], train_m)
    pred <- predict_ridge(fit, x[ids, , drop = FALSE])
    dimnames(pred) <- list(ids, colnames(measured))
    pred_m <- predict_ridge(fit, x[test_m, , drop = FALSE])
    molecular_predictions[[length(molecular_predictions) + 1L]] <- rbindlist(lapply(seq_along(test_m), function(i)
      data.table(fold = k, model = feature_name, depmap_id = test_m[i], program = colnames(measured),
        observed = measured[test_m[i], ], predicted = pred_m[i, ],
        training_mean = mean_response, training_sd = sd_response)))
    penalty_rows[[length(penalty_rows) + 1L]] <- data.table(fold = k, stage = 'molecular',
      model = feature_name, lambda = fit$lambda)
    # Cross-fit molecular training inputs by cell, conditional on the fixed outer bank.
    crossfit <- pred
    for (j in setdiff(1:5, k)) {
      target <- intersect(train_m[fold_map[train_m] == j], ids)
      if (!length(target)) next
      inner_train <- setdiff(train_m, target)
      inner_fit <- fit_ridge(x[inner_train, , drop = FALSE], measured[inner_train, , drop = FALSE], inner_train)
      crossfit[target, ] <- predict_ridge(inner_fit, x[target, , drop = FALSE])
    }
    for (kind in c('quadratic', 'conditioned')) {
      name <- paste(feature_name, kind, sep = '__')
      f <- fit_ridge(design(crossfit[train_f, , drop = FALSE], b[train_f, , drop = FALSE], kind),
        y[train_f], train_f, intercept = FALSE)
      p <- predict_ridge(f, design(pred[test_f, , drop = FALSE], b[test_f, , drop = FALSE], kind))
      add_function('all_functional', name, test_f, p)
      penalty_rows[[length(penalty_rows) + 1L]] <- data.table(fold = k, stage = 'chain', model = name, lambda = f$lambda)
      fit_store[[name]] <- f
    }
    fit_store[[paste0('upstream_', feature_name)]] <- fit
  }
  # Direct baseline prediction explicitly isolates information already present before treatment.
  for (feature_name in names(features)) {
    f <- fit_ridge(features[[feature_name]][train_f, , drop = FALSE], y[train_f], train_f)
    add_function('all_functional', paste0('direct_', feature_name), test_f,
      predict_ridge(f, features[[feature_name]][test_f, , drop = FALSE]))
    fit_store[[paste0('direct_', feature_name)]] <- f
  }
  # Fair measured-response test: TRAIN and TEST on measured responses, same paired cohort.
  paired_train <- intersect(train_f, molecular_ids)
  paired_test <- intersect(test_f, molecular_ids)
  support_rows[[k]] <- data.table(fold = k, molecular_train = length(train_m), molecular_test = length(test_m),
    paired_train = length(paired_train), paired_test = length(paired_test),
    paired_sensitive_train = sum(y[paired_train] <= -1), paired_sensitive_test = sum(y[paired_test] <= -1))
  for (kind in c('linear', 'quadratic', 'conditioned')) {
    f <- fit_ridge(design(measured[paired_train, , drop = FALSE], b[paired_train, , drop = FALSE], kind),
      y[paired_train], paired_train, intercept = FALSE)
    p <- predict_ridge(f, design(measured[paired_test, , drop = FALSE], b[paired_test, , drop = FALSE], kind))
    add_function('paired_fit_and_test', paste0('measured_', kind), paired_test, p)
    penalty_rows[[length(penalty_rows) + 1L]] <- data.table(fold = k, stage = 'measured_function', model = kind, lambda = f$lambda)
    fit_store[[paste0('measured_', kind)]] <- f
  }
  for (feature_name in names(features)) {
    f <- fit_ridge(features[[feature_name]][paired_train, , drop = FALSE], y[paired_train], paired_train)
    add_function('paired_fit_and_test', paste0('direct_', feature_name), paired_test,
      predict_ridge(f, features[[feature_name]][paired_test, , drop = FALSE]))
  }
  add_function('paired_fit_and_test', 'training_mean', paired_test, rep(mean(y[paired_train]), length(paired_test)))
  qsave(list(models = fit_store, source_bank_hash = hashes$sha256[k], fold = k, excluded_ids = bank$excluded_ids,
    scoring = 'hard_cluster_mean', weights = weights, measured_displacements = measured),
    file.path(model_dir, sprintf('diagnostic_fold_%02d.qs', k)))
  fwrite(rbindlist(molecular_predictions), file.path(data_dir, 'diagnostic_molecular_predictions.tsv'), sep = '\t')
  fwrite(rbindlist(functional_predictions), file.path(data_dir, 'diagnostic_functional_predictions.tsv'), sep = '\t')
  fwrite(rbindlist(penalty_rows), file.path(data_dir, 'diagnostic_selected_penalties.tsv'), sep = '\t')
  fwrite(rbindlist(support_rows), file.path(data_dir, 'diagnostic_training_support.tsv'), sep = '\t')
  fwrite(rbindlist(score_checks), file.path(data_dir, 'hard_cluster_score_checks.tsv'), sep = '\t')
  fwrite(data.table(status = 'running', completed_folds = k), file.path(out, 'diagnostic_status.tsv'), sep = '\t')
}
mp <- rbindlist(molecular_predictions)
fp <- rbindlist(functional_predictions)
stopifnot(all(is.finite(fp$predicted)), all(is.finite(mp$predicted)))
fm <- fp[, .(n = .N, pearson = cor_safe(observed, predicted), rmse = sqrt(mean((observed - predicted)^2))), by = .(cohort, model)]
mm <- mp[, .(cells = uniqueN(depmap_id), standardized_mse = mean(((observed - predicted) / training_sd)^2),
  mean_response_mse = mean(((observed - training_mean) / training_sd)^2)), by = model]
mm[, improvement_over_mean := 1 - standardized_mse / mean_response_mse]
fwrite(fm, file.path(data_dir, 'diagnostic_functional_metrics.tsv'), sep = '\t')
fwrite(mm, file.path(data_dir, 'diagnostic_molecular_metrics.tsv'), sep = '\t')
soft_metrics <- fread(file.path(data_dir, 'soft_diagnostic_functional_metrics.tsv'))
comparison <- merge(soft_metrics, fm, by = c('cohort', 'model'), suffixes = c('_soft', '_hard'))
stopifnot(all(comparison$n_soft == comparison$n_hard))
fwrite(comparison, file.path(data_dir, 'soft_vs_hard_functional_metrics.tsv'), sep = '\t')
fwrite(mp[, .(n = .N, r = cor_safe(observed, predicted), mse = mean((observed - predicted)^2),
  mean_mse = mean((observed - training_mean)^2)), by = .(fold, model, program)],
  file.path(data_dir, 'diagnostic_molecular_program_metrics.tsv'), sep = '\t')
for (cohort_name in unique(fp$cohort)) {
  d <- fp[cohort == cohort_name]
  annotation <- fm[cohort == cohort_name, .(model, label = sprintf('n=%d; r=%.2f; RMSE=%.2f', n, pearson, rmse))]
  d <- merge(d, annotation, by = 'model')
  d[, panel := paste(gsub('__', '\n', model, fixed = TRUE), label, sep = '\n')]
  plot <- ggplot(d, aes(observed, predicted)) + geom_point(size = .7, alpha = .6) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .4) +
    facet_wrap(~panel, scales = 'free', ncol = 3) + theme_classic(base_size = 9) +
    theme(aspect.ratio = 1, strip.text = element_text(size = 7)) +
    labs(x = 'Measured log2 abundance', y = 'Held-out prediction', title = gsub('_', ' ', cohort_name))
  ggsave(file.path(figure_dir, paste0('diagnostic_', cohort_name, '.pdf')), plot, width = 9, height = 9)
  ggsave(file.path(figure_dir, paste0('diagnostic_', cohort_name, '.png')), plot, width = 9, height = 9, dpi = 130)
}
report <- c('NEW PROGRAM CHAIN: STAGE DIAGNOSTICS', '',
  'HARD CLUSTER MEANS for both baseline expression and treatment displacement; no overlapping soft weights.',
  'Same five frozen training-only program sets; no reclustering. Retrospective evaluation.',
  '5 uM molecular response -> 2.5 uM corrected PRISM abundance; exposure proxy.',
  'Ridge penalties selected from 10, 1, .1, .01 using outer-training recipients only.',
  'Internal penalty selection conditions on the outer bank; only outer predictions are reported as validation.', '',
  'BASELINE -> DISPLACEMENT', capture.output(print(mm)), '',
  'FUNCTIONAL PREDICTION', capture.output(print(fm)), '',
  'Paired-fit tests train and test on measured responses; direct-baseline controls use the same paired cohort.',
  'All-functional chains use cross-fitted predicted responses for molecularly measured training recipients.',
  'The all-functional and paired-fit cohorts have different training sizes; compare within each cohort.', '',
  'SENSITIVE-RECIPIENT SUPPORT', capture.output(print(rbindlist(support_rows))), '',
  'SOFT VS HARD: SAME FOLDS AND COHORTS', capture.output(print(comparison)))
writeLines(report, file.path(out, 'STAGE_DIAGNOSTICS.txt'))
fwrite(data.table(status = 'complete', completed_folds = 5L), file.path(out, 'diagnostic_status.tsv'), sep = '\t')
print(mm)
print(fm)
