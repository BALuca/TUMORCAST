#!/usr/bin/env Rscript
# Fit the visible sensitive/resistant pattern, rather than treating its tails as outliers.
suppressPackageStartupMessages({library(data.table); library(qs); library(mgcv); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

root <- '/work/results/program_pilot/40_constrained_program_sensitivity_classifier'
basis_root <- '/work/results/program_pilot/39_constrained_additive_sensitivity_decoder'
score_path <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps/data/ERBB_all_training_inputs.qs'
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
input <- qread(score_path)
prepared_path <- file.path(parent, 'data/prepared_scores.qs')
prepared <- qread(prepared_path)
full_basis_path <- file.path(basis_root, 'models/constrained_decoder.qs')
full_source <- qread(full_basis_path)
settings <- qread(file.path(basis_root, 'settings.qs'))
md <- copy(input$annotations)
d <- input$displacement
programs <- colnames(d)
outcome <- as.numeric(md$measured_pct <= 50)
stopifnot(nrow(d) == 6260L, ncol(d) == 40L, identical(rownames(d), md$observation_id),
  all(is.finite(d)), !anyDuplicated(md$observation_id),
  identical(prepared$signature, full_source$source_signature),
  !any(md$depmap_id %in% settings$excluded_ids), !any(md$drug %in% settings$excluded_drugs))

# Reuse the exact frozen basis/scales from each corresponding fitting-cell set.
make_design <- function(ids, displacement, basis) {
  b <- sweep(sweep(prepared$baseline[ids, programs, drop = FALSE], 2L, basis$baseline_mean, '-'),
    2L, basis$baseline_sd, '/')
  curves <- do.call(cbind, lapply(seq_along(programs), function(j) {
    sweep(PredictMat(basis$smooths[[j]], data.frame(value = displacement[, j])), 2L, basis$curve_sd[[j]], '/')
  }))
  interaction <- sweep(displacement * b, 2L, basis$interaction_sd, '/')
  X <- cbind(curves, interaction)
  colnames(X) <- basis$feature_names
  stopifnot(ncol(X) == 160L, all(is.finite(X)))
  X
}
fit_classifier <- function(rows, basis_source) {
  stopifnot(identical(sort(unique(md$depmap_id[rows])), basis_source$training_ids),
    identical(basis_source$source_signature, prepared$signature))
  X <- make_design(md$depmap_id[rows], d[rows, , drop = FALSE], basis_source$basis)
  target <- outcome[rows]
  weights <- observation_weights(md[rows])
  # Equal total S/R evidence. Drug/recipient weights are retained within each class.
  for (value in 0:1) weights[target == value] <- weights[target == value] / sum(weights[target == value])
  weights <- weights / mean(weights)
  curvature <- matrix(0, 160L, 160L)
  for (j in seq_along(programs)) {
    sc <- basis_source$basis$curve_sd[[j]]
    penalty <- basis_source$basis$smooths[[j]]$S[[1L]] / outer(sc, sc)
    index <- (3L * j - 2L):(3L * j)
    curvature[index, index] <- penalty / sum(diag(penalty))
  }
  amplitude <- diag(c(rep(1, 120L), rep(20, 40L)))
  fit_warnings <- character()
  fit <- withCallingHandlers(gam(target ~ X, family = binomial(), weights = weights, method = 'REML',
    paraPen = list(X = list(amplitude, curvature)), control = gam.control(nthreads = 1L)),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      # These are fractional prior CASE weights, not binomial trial counts.
      # mgcv accepts them; record the initialization warning without repeating it.
      if (identical(conditionMessage(w), 'non-integer #successes in a binomial glm!')) invokeRestart('muffleWarning')
    })
  stopifnot(isTRUE(fit$converged), all(is.finite(coef(fit))))
  list(model = fit, basis = basis_source$basis, training_ids = basis_source$training_ids,
    context_keys = basis_source$context_keys, source_signature = prepared$signature,
    warnings = unique(fit_warnings))
}
predict_classifier <- function(fitted, ids, displacement) {
  X <- make_design(ids, displacement, fitted$basis)
  plogis(drop(cbind(1, X) %*% coef(fitted$model)))
}
measure <- function(metadata, score, keys, cohort) {
  z <- copy(metadata)
  z[, score := score]
  rbindlist(lapply(c(FALSE, TRUE), function(fixed) {
    zz <- if (fixed) z[condition %in% keys] else z
    data.table(cohort = cohort, axis = if (fixed) 'Informative fixed-drug/dose settings' else 'All observations',
      observations = nrow(zz), sensitive = sum(zz$measured_pct <= 50),
      sensitive_called = sum(zz$measured_pct <= 50 & zz$score >= .5),
      resistant = sum(zz$measured_pct > 50), false_sensitive = sum(zz$measured_pct > 50 & zz$score >= .5),
      sensitive_median_score = median(zz$score[zz$measured_pct <= 50]),
      resistant_median_score = median(zz$score[zz$measured_pct > 50]))
  }))
}

status('fitting one constrained class-balanced sensitivity classifier')
full <- fit_classifier(seq_len(nrow(md)), full_source)
qsave(full, file.path(root, 'models', 'constrained_classifier.qs'))
score <- predict_classifier(full, md$depmap_id, d)
stopifnot(max(abs(predict_classifier(qread(file.path(root, 'models', 'constrained_classifier.qs')), md$depmap_id, d) - score)) < 1e-10)
save_table(cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
  data.table(score = score)), 'training_predictions')
metrics <- list(measure(md, score, full$context_keys, 'Training fit'))
save_table(rbindlist(metrics), 'classification_counts')
save_table(data.table(effective_df = sum(full$model$edf), maximum_coefficients = 161L,
  amplitude_penalty = full$model$sp[1L], curvature_penalty = full$model$sp[2L],
  zero_displacement_score = plogis(coef(full$model)[1L])), 'fit_constraints')
status(paste('training sensitive called=', metrics[[1L]]$sensitive_called[1L], '/', metrics[[1L]]$sensitive[1L]))

oof <- rep(NA_real_, nrow(md))
eligible <- rep(FALSE, nrow(md))
fold_id <- rep(NA_integer_, nrow(md))
basis_paths <- file.path(basis_root, 'models', paste0('fold_', 1:5, '.qs'))
for (k in 1:5) {
  saved <- qread(basis_paths[k])
  tr <- which(md$depmap_id %in% saved$fitted$training_ids)
  va <- saved$validation_rows
  stopifnot(!length(intersect(md$depmap_id[tr], md$depmap_id[va])), length(tr) + length(va) == nrow(md))
  status(paste('whole-cell fold', k, '/5'))
  fitted <- fit_classifier(tr, saved$fitted)
  pp <- predict_classifier(fitted, md$depmap_id[va], d[va, , drop = FALSE])
  qsave(list(fitted = fitted, validation_rows = va, prediction = pp), file.path(root, 'models', paste0('fold_', k, '.qs')))
  oof[va] <- pp
  eligible[va] <- md$condition[va] %in% fitted$context_keys
  fold_id[va] <- k
  metrics[[k + 1L]] <- measure(md[va], pp, fitted$context_keys, paste0('Cell fold ', k))
  save_table(rbindlist(metrics), 'classification_counts')
  qsave(list(prediction = oof, eligible = eligible, fold = fold_id, observation_id = md$observation_id),
    file.path(root, 'data', 'cell_fold_predictions.qs'))
}
stopifnot(all(is.finite(oof)), !anyNA(fold_id))
status('freezing original held-cell scores')
query <- copy(as.data.table(prepared$metadata)[depmap_id %in% settings$excluded_ids & drug %in% md$drug])
query[, condition := paste(drug, dose, sep = '__')]
hd <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
query[, score := predict_classifier(full, depmap_id, hd)]
save_table(query, 'outcome_free_held_scores')
frozen_path <- file.path(root, 'data', 'outcome_free_held_scores.tsv')
save_table(data.table(path = frozen_path, sha256 = digest::digest(file = frozen_path, algo = 'sha256')), 'prediction_hash')
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
query[, measured_pct := labels$prism_abundance_pct[match(observation_id, labels$observation_id)]]
save_table(query[!is.finite(measured_pct) | measured_pct <= 0], 'unsupported_held_endpoints')
query <- query[is.finite(measured_pct) & measured_pct > 0]
query[, recipient := panel$cell_map$depmap_name[match(depmap_id, panel$cell_map$depmap_id)]]
stopifnot(uniqueN(query$depmap_id) == 9L, !anyNA(query$recipient))
metrics[[7L]] <- measure(query, query$score, full$context_keys, 'Original held cells')
save_table(rbindlist(metrics), 'classification_counts')
tab <- rbindlist(list(
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Training fit', score = score, cell_context = md$condition %in% full$context_keys, fold = 0L)),
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Development cell folds', score = oof, cell_context = eligible, fold = fold_id)),
  cbind(query[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct, score)],
    data.table(cohort = 'Original held cells', cell_context = query$condition %in% full$context_keys, fold = 0L))), use.names = TRUE)
stopifnot(is.numeric(tab$score), all(is.finite(tab$score)),
  setequal(unique(tab$cohort), c('Training fit', 'Development cell folds', 'Original held cells')))
tab[, actual_class := factor(ifelse(measured_pct <= 10, '<=10%', ifelse(measured_pct <= 50, '10-50%', '>50%')),
  levels = c('<=10%', '10-50%', '>50%'))]
save_table(tab, 'evaluated_scores')
pair_order <- tab[cell_context == TRUE, {
  si <- which(measured_pct <= 40)
  ri <- which(measured_pct >= 60)
  pairs <- CJ(sensitive_row = si, resistant_row = ri)
  pairs <- pairs[measured_pct[resistant_row] - measured_pct[sensitive_row] >= 25]
  gap <- score[pairs$sensitive_row] - score[pairs$resistant_row]
  list(clear_pairs = length(gap), correctly_ordered_fraction = if (length(gap)) mean((gap > 0) + .5 * (gap == 0)) else NA_real_)
}, by = .(cohort, fold, condition)]
save_table(pair_order, 'within_treatment_pair_order')
save_table(pair_order[, .(contexts = sum(is.finite(correctly_ordered_fraction)),
  clear_pairs = sum(clear_pairs), mean_context_correct_order_pct = 100 * mean(correctly_ordered_fraction, na.rm = TRUE)),
  by = cohort], 'within_treatment_pair_summary')
plot_data <- rbindlist(list(cbind(tab, axis = 'All observations'),
  cbind(tab[cell_context == TRUE], axis = 'Informative fixed-drug/dose settings')))
plot_data[, cohort := factor(cohort, levels = c('Training fit', 'Development cell folds', 'Original held cells'))]
plot_data[, axis := factor(axis, levels = c('All observations', 'Informative fixed-drug/dose settings'))]
p <- ggplot(plot_data, aes(actual_class, score, fill = actual_class)) +
  geom_boxplot(width = .65, outlier.shape = NA) + geom_hline(yintercept = .5, linetype = 2) +
  facet_grid(axis ~ cohort) + scale_fill_manual(values = c('#D55E00', '#E69F00', '#0072B2')) +
  theme_classic(base_size = 10) + theme(legend.position = 'none') +
  labs(x = 'Measured abundance (% vehicle)', y = 'Class-balanced sensitivity score',
    title = 'Can a constrained model learn the visible sensitivity pattern?',
    caption = 'Score >=0.5 calls sensitive. Scores are not prevalence-calibrated probabilities or predicted viability. Held results are retrospective.')
ggsave(file.path(root, 'figures', 'sensitivity_pattern_capture.pdf'), p, device = grDevices::cairo_pdf, width = 10, height = 6)
ggsave(file.path(root, 'figures', 'sensitivity_pattern_capture.png'), p, width = 10, height = 6, dpi = 180)
sources <- c(score_path, prepared_path, full_basis_path, basis_paths,
  '/work/src/program_pilot/40_constrained_program_sensitivity_classifier.R', '/work/src/program_pilot/response_models.R')
save_table(data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'input_hashes')
qsave(list(source_signature = prepared$signature, training_ids = full$training_ids,
  excluded_ids = settings$excluded_ids, excluded_drugs = settings$excluded_drugs,
  objective = 'class-balanced binomial deviance; <=50% abundance is sensitive',
  constraints = 'three curve coefficients/program; two shared REML penalties; baseline amplitude multiplier 20',
  interpretation = 'class-balanced score; not predicted viability or naturally calibrated probability'), file.path(root, 'settings.qs'))
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
print(rbindlist(metrics))
status('COMPLETE')
