#!/usr/bin/env Rscript
# Lean joint baseline/response test on the frozen MEK bank; no program discovery.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

family <- 'MEK'
pca_fraction <- .90
grid <- data.table(candidate = 1:3, cost = c(10, 100, 100), gamma_multiplier = c(1, 4, 16))
parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
root <- '/work/results/program_pilot/29_joint_baseline_displacement_svr_MEK_holdout'
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
status <- function(stage) {
  message(format(Sys.time()), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading frozen inputs')
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
prepared <- qread(file.path(parent, 'data/prepared_scores.qs'))
settings <- qread(file.path(parent, 'settings.qs'))
old <- qread(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'))
folds <- fread(file.path(parent, 'data/nonlinear_training_folds.tsv'))
labels <- as.data.table(panel$sensitivity)
train <- labels[depmap_id %in% old$training_ids & !drug %in% settings$excluded_drugs &
  is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
prior_path <- file.path(parent, 'data/contextual_predictions.tsv.gz')
prior <- fread(prior_path)[model == 'baseline_reference']
query <- prior[, .(observation_id, depmap_id, recipient, drug, dose)]
stopifnot(!anyDuplicated(query$observation_id), !anyDuplicated(train$observation_id),
  uniqueN(train$depmap_id) == 26L, uniqueN(query$depmap_id) == 9L,
  !any(query$depmap_id %in% train$depmap_id), !any(query$drug %in% settings$excluded_drugs))
fold <- folds$fold[match(train$depmap_id, folds$depmap_id)]
stopifnot(!anyNA(fold), identical(sort(unique(fold)), 1:3))
y <- log2(train$prism_abundance_pct / 100)
d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
dt <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
stopifnot(ncol(d) == 40L, all(is.finite(d)), all(is.finite(dt)))
sources <- c(panel_path, prior_path, file.path(parent, c('settings.qs', 'data/prepared_scores.qs',
  'models/nonlinear_RBF_SVR_D.qs', 'data/nonlinear_training_folds.tsv')),
  '/work/src/program_pilot/29_joint_baseline_displacement_svr.R',
  '/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R')
manifest <- data.table(path = sources, sha256 = vapply(sources, digest::digest,
  character(1), file = TRUE, algo = 'sha256'))
save_table(manifest, 'input_hashes')
signature <- digest::digest(list(manifest, grid, pca_fraction), algo = 'sha256')
qsave(list(signature, family, grid, pca_fraction, training_ids = old$training_ids,
  held_ids = settings$held_ids, selection = 'three-view whole-cell validation',
  endpoint = 'log2 PRISM abundance relative to vehicle',
  anchoring = 'absolute joint readout; no post-fit zero subtraction'), file.path(root, 'settings.qs'))

# Context compression and molecular means use fitting recipients exclusively.
learn_baseline <- function(ids) {
  pca <- prcomp(prepared$baseline[ids, , drop = FALSE], center = TRUE, scale. = TRUE)
  components <- seq_len(which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= pca_fraction)[1L])
  list(pca = pca, components = components, training_ids = ids)
}
baseline_input <- function(pca, ids) {
  x <- predict(pca$pca, newdata = prepared$baseline[ids, , drop = FALSE])[, pca$components, drop = FALSE]
  colnames(x) <- paste0('Baseline_PC', pca$components)
  x
}
mean_response <- function(ids, md) {
  rows <- which(prepared$metadata$depmap_id %in% ids)
  keys <- paste(prepared$metadata$drug[rows], prepared$metadata$dose[rows])
  sums <- rowsum(prepared$measured[rows, , drop = FALSE], keys, reorder = FALSE)
  counts <- tabulate(match(keys, rownames(sums)), nrow(sums))
  requested <- match(paste(md$drug, md$dose), rownames(sums))
  stopifnot(!anyNA(requested), !any(ids %in% settings$held_ids))
  (sums / counts)[requested, , drop = FALSE]
}
validation_loss <- function(prediction, md, outcome) {
  views <- context_response_views(as.matrix(prediction), outcome, md)
  sum(abs(views$x[, 1L] - views$y) * views$weights) / sum(views$weights)
}

# Only two fitted variants: cell-specific D and training treatment-mean D.
validation <- list()
oof <- list(actual = matrix(NA_real_, nrow(train), nrow(grid)),
  average = matrix(NA_real_, nrow(train), nrow(grid)))
for (k in 1:3) {
  fitting <- which(fold != k)
  testing <- which(fold == k)
  ids <- unique(train$depmap_id[fitting])
  pca <- learn_baseline(ids)
  bf <- baseline_input(pca, train$depmap_id[fitting])
  bv <- baseline_input(pca, train$depmap_id[testing])
  mean_fitting <- mean_response(ids, train[fitting])
  mean_testing <- mean_response(ids, train[testing])
  for (variant in names(oof)) {
    xx <- cbind(bf, if (variant == 'actual') d[fitting, , drop = FALSE] else mean_fitting)
    vx <- cbind(bv, if (variant == 'actual') d[testing, , drop = FALSE] else mean_testing)
    for (j in seq_len(nrow(grid))) {
      status(paste('fold', k, variant, 'capacity', j, 'of 3'))
      path <- file.path(root, 'models', sprintf('%s_candidate%d_fold%d.qs', variant, j, k))
      if (file.exists(path)) {
        saved <- qread(path)
        stopifnot(identical(saved$signature, signature), identical(saved$validation_rows, testing))
        prediction <- saved$prediction
      } else {
        fit <- fit_program_svr(xx, y[fitting], 'radial', grid$cost[j], grid$gamma_multiplier[j] / ncol(xx))
        prediction <- predict_program_svr(fit, vx)
        qsave(list(signature = signature, model = fit, pca = pca, training_ids = ids,
          validation_rows = testing, prediction = prediction), path)
      }
      oof[[variant]][testing, j] <- prediction
      validation[[paste(k, variant, j)]] <- data.table(fold = k, variant = variant, candidate = j,
        loss = validation_loss(prediction, train[testing], y[testing]))
      save_table(rbindlist(validation), 'validation')
    }
  }
}
pca <- learn_baseline(old$training_ids)
bf <- baseline_input(pca, train$depmap_id)
bq <- baseline_input(pca, query$depmap_id)
mean_fitting <- mean_response(old$training_ids, train)
mean_query <- mean_response(old$training_ids, query)
predictions <- selections <- list()
for (variant in names(oof)) {
  loss <- apply(oof[[variant]], 2L, validation_loss, md = train, outcome = y)
  best <- which.min(loss)
  selections[[variant]] <- cbind(data.table(variant = variant, baseline_components = length(pca$components)), grid[best])
  save_table(rbindlist(selections), 'selected_settings')
  xx <- cbind(bf, if (variant == 'actual') d else mean_fitting)
  vx <- cbind(bq, if (variant == 'actual') dt else mean_query)
  status(paste('final', variant, 'selected capacity', best))
  fit <- fit_program_svr(xx, y, 'radial', grid$cost[best], grid$gamma_multiplier[best] / ncol(xx))
  qsave(list(signature = signature, model = fit, pca = pca, selected = grid[best],
    training_ids = old$training_ids, training_drugs = unique(train$drug),
    mean_fitting = mean_fitting, mean_query = mean_query), file.path(root, 'models', paste0(variant, '.qs')))
  prediction <- predict_program_svr(fit, vx)
  tab <- copy(query)
  tab[, `:=`(model = paste0('Joint_', variant), predicted = prediction)]
  predictions[[variant]] <- tab
  save_table(rbindlist(predictions), 'outcome_free_predictions')
}
reference <- copy(prior[, .(observation_id, depmap_id, recipient, drug, dose, predicted)])
reference[, model := 'Baseline_drug_dose']
predictions <- rbindlist(c(predictions, list(reference)), use.names = TRUE)
save_table(predictions, 'outcome_free_predictions')
prediction_path <- file.path(root, 'data/outcome_free_predictions.tsv')
save_table(data.table(path = prediction_path, sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')

# Native assay-scale reporting isolates variation across cells at fixed treatment.
predictions[, observed_pct := labels$prism_abundance_pct[match(observation_id, labels$observation_id)]]
predictions[, predicted_pct := 100 * 2^predicted]
stopifnot(all(is.finite(predictions$predicted_pct)), all(is.finite(predictions$observed_pct)))
predictions[, `:=`(observed_centered = observed_pct - mean(observed_pct),
  predicted_centered = predicted_pct - mean(predicted_pct)), by = .(model, drug, dose)]
save_table(predictions, 'predictions')
summary <- predictions[, .(observations = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug),
  within_treatment_r = cor_safe(observed_centered, predicted_centered),
  sensitive = sum(observed_pct <= 50), sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  resistant = sum(observed_pct > 50), false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)), by = model]
save_table(summary, 'summary')
per_drug <- predictions[, .(conditions = .N,
  within_treatment_r = cor_safe(observed_centered, predicted_centered)), by = .(model, drug)]
save_table(per_drug, 'per_drug')
save_table(predictions[, .(within_cell_r = cor_safe(observed_pct, predicted_pct)), by = .(model, recipient, dose)], 'per_cell_dose')
plot <- ggplot(predictions, aes(observed_centered, predicted_centered)) +
  geom_point(alpha = .12, size = .6) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#2166AC') +
  geom_text(data = summary, aes(x = -Inf, y = Inf, label = sprintf('r = %.2f', within_treatment_r)),
    inherit.aes = FALSE, hjust = -.1, vjust = 1.3, size = 3) +
  facet_wrap(vars(model), nrow = 1, scales = 'free') + theme_classic(base_size = 9) +
  theme(aspect.ratio = 1) + labs(x = 'Measured difference from drug/dose mean (percentage points)',
    y = 'Predicted difference (percentage points)', title = 'Does actual molecular response distinguish held cells?')
ggsave(file.path(root, 'figures/within_treatment_predictions.pdf'), plot, width = 10, height = 3.8)
ggsave(file.path(root, 'figures/within_treatment_predictions.png'), plot, width = 10, height = 3.8, dpi = 140)
print(summary)
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('COMPLETE')
