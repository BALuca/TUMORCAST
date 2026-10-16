#!/usr/bin/env Rscript
# Version 22: PCA baseline and common hotspot conditioning of the frozen SVR.
# Response programs, held cells, excluded families and fitting views stay fixed.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
})
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(2L)
RhpcBLASctl::omp_set_num_threads(2L)
set.seed(20261010L)
options(warn = 1)

# Outcome-independent baseline compression and fixed mutation panel.
pca_variance <- .90
penalty_grid <- c(1000, 100, 10, 1, .1, .01)
driver_genes <- c('APC', 'BRAF', 'KRAS', 'NRAS', 'HRAS', 'EGFR', 'ERBB2', 'PIK3CA', 'PTEN',
  'TP53', 'RB1', 'STK11', 'CDKN2A', 'NF1', 'SMARCA4', 'KEAP1', 'CTNNB1', 'IDH1', 'IDH2',
  'KIT', 'PDGFRA', 'FGFR2', 'FGFR3', 'ALK', 'MET', 'MTOR', 'AKT1', 'MAP2K1', 'MAP2K2')
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, grepl('^--family-transfer=(MEK|ERBB)$', args))
family <- sub('--family-transfer=', '', args, fixed = TRUE)
parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
root <- paste0('/work/results/program_pilot/22_compressed_baseline_mutation_conditioning_', family, '_holdout')
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
commit <- function(x, name) {
  path <- file.path(root, 'models', paste0(name, '.qs'))
  tmp <- tempfile(tmpdir = dirname(path))
  qsave(x, tmp, preset = 'high')
  stopifnot(file.rename(tmp, path))
}
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', family, ' | ', stage)
  fwrite(data.table(stage = stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {
  status('FAILED: see run.log')
  traceback(8)
  quit(save = 'no', status = 1L)
})
status('loading frozen parent inputs')
settings <- qread(file.path(parent, 'settings.qs'))
prepared <- qread(file.path(parent, 'data/prepared_scores.qs'))
base <- qread(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'))
parent_context <- qread(file.path(parent, 'contextual_settings.qs'))
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
mutation_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv'
header <- names(fread(mutation_path, nrows = 0L))
symbols <- sub(' \\(.*', '', header[-1L])
stopifnot(all(driver_genes %in% symbols))
selected_columns <- match(driver_genes, symbols)
mutation_table <- fread(mutation_path, select = c(1L, selected_columns + 1L))
setnames(mutation_table, c('depmap_id', driver_genes))
stopifnot(!anyDuplicated(mutation_table$depmap_id), all(rownames(prepared$baseline) %in% mutation_table$depmap_id))
mutations <- as.matrix(mutation_table[match(rownames(prepared$baseline), depmap_id), ..driver_genes])
dimnames(mutations) <- list(rownames(prepared$baseline), paste0('MUT_', driver_genes))
stopifnot(!anyNA(mutations), all(mutations %in% 0:1))
labels <- as.data.table(panel$sensitivity)
train <- labels[depmap_id %in% base$training_ids & !drug %in% settings$excluded_drugs &
  is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
train[, observed := log2(prism_abundance_pct / 100)]
parent_predictions <- fread(file.path(parent, 'data/contextual_predictions.tsv.gz'))
query <- parent_predictions[model == 'anchored_actual', .(observation_id, depmap_id, recipient, drug, dose)]
d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
dt <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
programs <- colnames(d)
folds <- fread(file.path(parent, 'data/nonlinear_training_folds.tsv'))
fold <- folds$fold[match(train$depmap_id, folds$depmap_id)]
validation <- fread(file.path(parent, 'data/nonlinear_validation.tsv'))
candidate <- validation[model == 'RBF_SVR_D' & cost == base$selected$cost &
  gamma_multiplier == base$selected$gamma_multiplier, candidate]
stopifnot(length(candidate) == 1L, !anyNA(fold), all(is.finite(d)), all(is.finite(dt)),
  !anyDuplicated(query$observation_id), setequal(unique(train$depmap_id), base$training_ids),
  !length(intersect(base$training_ids, settings$held_ids)), identical(base$training_ids, parent_context$training_ids))

# Hash all reused fitting checkpoints and code. Parent files remain read-only.
checkpoint_paths <- c(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'),
  file.path(parent, 'models', sprintf('nonlinear_RBF_SVR_D_candidate%d_fold%d.qs', candidate, 1:3)),
  file.path(parent, 'models', sprintf('contextual_substrate_fold%d_inner%d.qs',
    rep(1:3, each = 2L), rep(1:2, 3L))))
paths <- c('/work/src/program_pilot/22_compressed_baseline_mutation_conditioning.R',
  '/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R',
  '/work/src/R/lib_heatmap.R', panel_path, mutation_path, file.path(parent, 'settings.qs'),
  file.path(parent, 'data/prepared_scores.qs'), checkpoint_paths)
manifest <- data.table(path = paths, sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
signature <- digest::digest(list(manifest, pca_variance, driver_genes, penalty_grid), algo = 'sha256')
settings_path <- file.path(root, 'settings.qs')
if (file.exists(settings_path)) stopifnot(identical(qread(settings_path)$signature, signature))
qsave(list(signature = signature, family = family, held_ids = settings$held_ids,
  training_ids = base$training_ids, excluded_drugs = settings$excluded_drugs,
  pca_variance = pca_variance, driver_genes = driver_genes, penalty_grid = penalty_grid,
  parent = parent, evaluation = 'Retrospective; previously inspected held outcomes'), settings_path)
write_table(manifest, 'input_hashes')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))

# PCA is fitted to distinct training baselines. Mutation flags stay on their
# binary scale, preventing rare mutations being amplified by variance scaling.
learn_baseline <- function(ids) {
  stopifnot(!length(intersect(ids, settings$held_ids)))
  pca <- prcomp(prepared$baseline[ids, , drop = FALSE], center = TRUE, scale. = TRUE)
  explained <- pca$sdev^2 / sum(pca$sdev^2)
  components <- seq_len(which(cumsum(explained) >= pca_variance)[1L])
  mutation_mean <- colMeans(mutations[ids, , drop = FALSE])
  active_mutations <- names(mutation_mean)[mutation_mean > 0 & mutation_mean < 1]
  list(pca = pca, components = components, mutation_mean = mutation_mean,
    active_mutations = active_mutations, training_ids = ids, explained = explained)
}
baseline_input <- function(representation, ids, kind) {
  pc <- predict(representation$pca, newdata = prepared$baseline[ids, , drop = FALSE])[, representation$components, drop = FALSE]
  pc <- sweep(pc, 2L, representation$pca$sdev[representation$components], '/')
  colnames(pc) <- sprintf('PC%02d', representation$components)
  if (kind == 'pca') return(pc)
  mm <- mutations[ids, representation$active_mutations, drop = FALSE]
  mm <- sweep(mm, 2L, representation$mutation_mean[representation$active_mutations], '-')
  cbind(pc, mm)
}
mean_response <- function(ids, metadata) {
  rows <- which(prepared$metadata$depmap_id %in% ids)
  keys <- paste(prepared$metadata$drug[rows], prepared$metadata$dose[rows])
  totals <- rowsum(prepared$measured[rows, , drop = FALSE], keys, reorder = FALSE)
  counts <- tabulate(match(keys, rownames(totals)), nbins = nrow(totals))
  idx <- match(paste(metadata$drug, metadata$dose), rownames(totals))
  (totals / counts)[idx, , drop = FALSE]
}
fit_context <- function(x, md) {
  views <- context_response_views(x, md$observed, md)
  glmnet(views$x, views$y, alpha = 0, lambda = penalty_grid, intercept = FALSE,
    standardize = FALSE, weights = views$weights, penalty.factor = c(0, rep(1, ncol(x) - 1L)), thresh = 1e-6)
}
score_context <- function(prediction, md) {
  views <- context_response_views(as.matrix(prediction), md$observed, md)
  colSums(abs(views$x - views$y) * views$weights) / sum(views$weights)
}
predict_compressed <- function(model, ids, response) {
  b <- baseline_input(model$representation, ids, model$kind)
  predict_context_response(model, b, response)
}
variants <- CJ(kind = c('pca', 'pca_mutation'), response = c('actual', 'average'))
variants[, name := paste(kind, response, sep = '_')]
losses <- setNames(lapply(variants$name, function(x) rep(0, length(penalty_grid))), variants$name)
oof <- list(actual = rep(NA_real_, nrow(train)), average = rep(NA_real_, nrow(train)))
mean_oof <- matrix(NA_real_, nrow(train), ncol(d), dimnames = dimnames(d))
cv_predictions <- cv_losses <- representations <- mutation_coverage <- list()
for (k in 1:3) {
  status(paste('context fold', k, '/3'))
  tr <- which(fold != k)
  va <- which(fold == k)
  ids <- unique(train$depmap_id[tr])
  representation <- learn_baseline(ids)
  stopifnot(!length(intersect(ids, unique(train$depmap_id[va]))))
  representations[[as.character(k)]] <- data.table(fold = k, fitting_recipients = length(ids),
    components = length(representation$components), retained_variance = sum(representation$explained[representation$components]))
  mutation_coverage[[as.character(k)]] <- data.table(fold = k, gene = driver_genes,
    positive = colSums(mutations[ids, , drop = FALSE]), recipients = length(ids),
    included = colnames(mutations) %in% representation$active_mutations)
  outer <- qread(checkpoint_paths[1L + k])
  stopifnot(identical(outer$validation_rows, va), identical(outer$signature, base$signature))
  validation_zero <- predict_program_svr(outer$model, matrix(0, 1L, ncol(d)))
  validation_s <- predict_program_svr(outer$model, d[va, , drop = FALSE]) - validation_zero
  stopifnot(max(abs(validation_s + validation_zero - base$oof[va])) < 1e-10)
  validation_mean <- mean_response(ids, train[va])
  oof$actual[va] <- validation_s
  supported <- complete.cases(validation_mean)
  oof$average[va[supported]] <- predict_program_svr(outer$model, validation_mean[supported, , drop = FALSE]) - validation_zero
  mean_oof[va, ] <- validation_mean
  inner_s <- list(actual = rep(NA_real_, length(tr)), average = rep(NA_real_, length(tr)))
  inner_mean <- matrix(NA_real_, length(tr), ncol(d), dimnames = list(NULL, programs))
  for (j in 1:2) {
    inner <- qread(file.path(parent, 'models', sprintf('contextual_substrate_fold%d_inner%d.qs', k, j)))
    iv <- inner$validation_rows
    stopifnot(!length(intersect(inner$training_rows, va)), !length(intersect(iv, va)))
    ii <- match(iv, tr)
    zero <- predict_program_svr(inner$model, matrix(0, 1L, ncol(d)))
    inner_s$actual[ii] <- predict_program_svr(inner$model, d[iv, , drop = FALSE]) - zero
    mm <- mean_response(unique(train$depmap_id[inner$training_rows]), train[iv])
    ok <- complete.cases(mm)
    inner_s$average[ii[ok]] <- predict_program_svr(inner$model, mm[ok, , drop = FALSE]) - zero
    inner_mean[ii, ] <- mm
  }
  stopifnot(all(is.finite(inner_s$actual)))
  for (v in seq_len(nrow(variants))) {
    variant <- variants[v]
    btr <- baseline_input(representation, train$depmap_id[tr], variant$kind)
    bva <- baseline_input(representation, train$depmap_id[va], variant$kind)
    scaling <- list(b_mean = setNames(rep(0, ncol(btr)), colnames(btr)),
      b_sd = setNames(rep(1, ncol(btr)), colnames(btr)), d_sd = apply(d[tr, ], 2L, sd))
    dd <- if (variant$response == 'actual') d[tr, ] else inner_mean
    vd <- if (variant$response == 'actual') d[va, ] else validation_mean
    fr <- which(complete.cases(dd) & is.finite(inner_s[[variant$response]]))
    vr <- which(complete.cases(vd) & is.finite(oof[[variant$response]][va]))
    xx <- context_response_design(btr[fr, , drop = FALSE], dd[fr, , drop = FALSE], inner_s[[variant$response]][fr], scaling)
    vx <- context_response_design(bva[vr, , drop = FALSE], vd[vr, , drop = FALSE], oof[[variant$response]][va[vr]], scaling)
    fit <- fit_context(xx, train[tr[fr]])
    pp <- as.matrix(predict(fit, newx = vx, s = penalty_grid))
    score <- score_context(pp, train[va[vr]])
    losses[[variant$name]] <- losses[[variant$name]] + score / 3
    commit(list(signature = signature, fit = fit, representation = representation, scaling = scaling,
      validation_rows = va[vr], prediction = pp), paste0(variant$name, '_fold', k))
    cv_losses[[paste(k, v)]] <- data.table(fold = k, variant = variant$name, lambda = penalty_grid,
      internal_balanced_absolute_log2_error = score, fitting_rows = length(fr), excluded_fitting_rows = length(tr) - length(fr))
    cv_predictions[[paste(k, v)]] <- data.table(variant = variant$name, fold = k,
      lambda = rep(penalty_grid, each = length(vr)), observation_id = rep(train$observation_id[va[vr]], length(penalty_grid)),
      predicted = as.vector(pp))
    write_table(rbindlist(cv_losses), 'development_selection')
  }
}
write_table(rbindlist(representations), 'fold_baseline_compression')
write_table(rbindlist(mutation_coverage), 'fold_mutation_coverage')
cv <- rbindlist(cv_predictions)
fwrite(cv, file.path(root, 'data/development_predictions.tsv.gz'), sep = '\t')
representation <- learn_baseline(base$training_ids)
query_mean <- mean_response(base$training_ids, query)
stopifnot(all(is.finite(query_mean)))
final_predictions <- list()
components <- coefficient_tables <- selection <- verification <- list()
models <- list()
for (v in seq_len(nrow(variants))) {
  variant <- variants[v]
  status(paste('final model', variant$name))
  b <- baseline_input(representation, train$depmap_id, variant$kind)
  scaling <- list(b_mean = setNames(rep(0, ncol(b)), colnames(b)),
    b_sd = setNames(rep(1, ncol(b)), colnames(b)), d_sd = apply(d, 2L, sd))
  input <- if (variant$response == 'actual') d else mean_oof
  rows <- which(complete.cases(input) & is.finite(oof[[variant$response]]))
  xx <- context_response_design(b[rows, , drop = FALSE], input[rows, , drop = FALSE], oof[[variant$response]][rows], scaling)
  fit <- fit_context(xx, train[rows])
  lambda <- penalty_grid[which.min(losses[[variant$name]])]
  model <- list(signature = signature, kind = variant$kind, representation = representation,
    fit = fit, lambda = lambda, scaling = scaling, svr = base$model,
    svr_zero = predict_program_svr(base$model, matrix(0, 1L, ncol(d))), training_ids = base$training_ids)
  commit(model, variant$name)
  input_test <- if (variant$response == 'actual') dt else query_mean
  result <- predict_compressed(model, query$depmap_id, input_test)
  zero <- predict_compressed(model, query$depmap_id, input_test * 0)$predicted
  replay <- predict_compressed(qread(file.path(root, 'models', paste0(variant$name, '.qs'))), query$depmap_id, input_test)$predicted
  stopifnot(all(is.finite(result$predicted)), max(abs(zero)) == 0, max(abs(replay - result$predicted)) < 1e-12)
  models[[variant$name]] <- model
  final_predictions[[variant$name]] <- cbind(copy(query), data.table(model = variant$name, predicted = result$predicted))
  components[[variant$name]] <- cbind(copy(query), data.table(model = variant$name,
    shared = result$shared, contextual = result$contextual, predicted = result$predicted))
  coefficients <- as.numeric(coef(fit, s = lambda))[-1L]
  coefficient_tables[[variant$name]] <- data.table(model = variant$name,
    baseline_feature = rep(colnames(b), times = ncol(d)), displaced_program = rep(programs, each = ncol(b)),
    coefficient = coefficients[-1L])
  selection[[variant$name]] <- data.table(model = variant$name, lambda = lambda,
    pca_components = length(representation$components), mutation_features = if (variant$kind == 'pca') 0L else length(representation$active_mutations),
    interaction_features = ncol(xx) - 1L, svr_gain = coefficients[1L])
  verification[[variant$name]] <- data.table(model = variant$name, predictions = nrow(query),
    zero_input_max_abs = max(abs(zero)), saved_prediction_replay_max_abs = max(abs(replay - result$predicted)))
}
controls <- parent_predictions[model %in% c('anchored_actual', 'anchored_average', 'context_actual', 'context_average')]
controls[, model := paste0('previous_', model)]
predictions <- rbindlist(c(final_predictions, list(controls)))
stopifnot(!anyDuplicated(predictions[, .(observation_id, model)]), all(is.finite(predictions$predicted)))
prediction_path <- file.path(root, 'data/predictions.tsv.gz')
fwrite(predictions, prediction_path, sep = '\t')
write_table(data.table(path = prediction_path, sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')
write_table(rbindlist(components), 'components')
write_table(rbindlist(coefficient_tables), 'interaction_coefficients')
write_table(rbindlist(selection), 'selected_models')
write_table(rbindlist(verification), 'verification')
write_table(data.table(component = seq_along(representation$explained), variance_fraction = representation$explained,
  cumulative_variance = cumsum(representation$explained), retained = seq_along(representation$explained) %in% representation$components), 'pca_variance')
rotation <- representation$pca$rotation[, representation$components, drop = FALSE]
write_table(cbind(data.table(program = rownames(rotation)), as.data.table(rotation)), 'pca_loadings')
mutation_long <- melt(mutation_table[depmap_id %in% rownames(prepared$baseline)],
  id.vars = 'depmap_id', variable.name = 'gene', value.name = 'hotspot_flag')
mutation_long[, `:=`(cell = panel$cell_map$depmap_name[match(depmap_id, panel$cell_map$depmap_id)],
  cohort = ifelse(depmap_id %in% base$training_ids, 'Training cells',
    ifelse(depmap_id %in% settings$held_ids, 'Held cells', 'No functional training data')))]
write_table(mutation_long, 'mutation_features')
write_table(mutation_long[, .(recipients = .N, positive = sum(hotspot_flag)), by = .(gene, cohort)], 'mutation_coverage')
inactive <- setdiff(colnames(mutations), representation$active_mutations)
unsupported <- melt(as.data.table(mutations[settings$held_ids, inactive, drop = FALSE], keep.rownames = 'depmap_id'),
  id.vars = 'depmap_id', variable.name = 'mutation', value.name = 'flag')[flag == 1]
write_table(unsupported, 'unsupported_held_mutations')

# Reporting begins after outcome-free predictions are frozen and hashed.
status('reporting held-cell prediction and mutation-associated displacement')
evaluation <- copy(predictions)
ii <- match(evaluation$observation_id, labels$observation_id)
evaluation[, `:=`(observed_pct = labels$prism_abundance_pct[ii], predicted_pct = 100 * 2^predicted,
  cell = recipient, drug_set = ifelse(drug %in% settings$excluded_drugs, 'Excluded family', 'Seen drugs'))]
evaluation[, observed := log2(observed_pct / 100)]
evaluation[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted),
  yp = observed_pct - mean(observed_pct), pp = predicted_pct - mean(predicted_pct)), by = .(model, drug, dose)]
write_table(evaluation, 'evaluation')
by_cell <- evaluation[, .(r = cor_safe(observed, predicted)), by = .(model, drug_set, cell, dose)]
summary <- evaluation[, .(across_cells_r_log2 = cor_safe(yc, pc), across_cells_r_pct = cor_safe(yp, pp),
  sensitive = sum(observed_pct <= 50), sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  resistant = sum(observed_pct > 50), false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)), by = .(model, drug_set)]
summary <- merge(summary, by_cell[, .(median_across_drugs_fixed_dose_r = median(r, na.rm = TRUE)),
  by = .(model, drug_set)], by = c('model', 'drug_set'))
write_table(summary, 'summary')
cv[, `:=`(drug = train$drug[match(observation_id, train$observation_id)],
  dose = train$dose[match(observation_id, train$observation_id)],
  cell = train$recipient[match(observation_id, train$observation_id)], observed = train$observed[match(observation_id, train$observation_id)])]
cv[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)), by = .(variant, lambda, drug, dose)]
cv_axes <- cv[, .(across_cells_r = cor_safe(yc, pc)), by = .(variant, lambda)]
cv_cell <- cv[, .(r = cor_safe(observed, predicted)), by = .(variant, lambda, cell, dose)]
cv_axes <- merge(cv_axes, cv_cell[, .(across_drugs_r = median(r, na.rm = TRUE)), by = .(variant, lambda)],
  by = c('variant', 'lambda'))
write_table(cv_axes, 'development_axes')

# Descriptive BRAF group differences use all measured cells, with no outcome
# screen. The permutation compares the centroid gap to random equally sized
# cell groups; lineage confounding prevents a causal mutation interpretation.
braf_profiles <- braf_summary <- braf_differences <- list()
for (dose_value in c(.05, .5, 5)) {
  md <- prepared$metadata[drug == 'vemurafenib' & dose == dose_value]
  stopifnot(nrow(md) > 0L, !anyDuplicated(md$depmap_id))
  dd <- prepared$measured[match(md$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  flag <- mutation_table$BRAF[match(md$depmap_id, mutation_table$depmap_id)]
  stopifnot(!anyNA(flag))
  nm <- sum(flag == 1)
  nw <- sum(flag == 0)
  stopifnot(nm >= 2L, nw >= 3L)
  dz <- sweep(dd, 2L, apply(d, 2L, sd), '/')
  delta <- colMeans(dd[flag == 1, , drop = FALSE]) - colMeans(dd[flag == 0, , drop = FALSE])
  statistic <- sqrt(mean((colMeans(dz[flag == 1, , drop = FALSE]) - colMeans(dz[flag == 0, , drop = FALSE]))^2))
  null <- replicate(499L, {
    selected <- sample.int(nrow(dd), nm)
    sqrt(mean((colMeans(dz[selected, , drop = FALSE]) - colMeans(dz[-selected, , drop = FALSE]))^2))
  })
  braf_summary[[as.character(dose_value)]] <- data.table(dose = dose_value, mutant_cells = nm, other_cells = nw,
    standardized_centroid_distance = statistic, random_group_median = median(null),
    permutation_p = (1 + sum(null >= statistic)) / (length(null) + 1),
    raw_centroid_difference = sqrt(mean(delta^2)), cosine_centroids =
      sum(colMeans(dd[flag == 1, ]) * colMeans(dd[flag == 0, ])) /
      sqrt(sum(colMeans(dd[flag == 1, ])^2) * sum(colMeans(dd[flag == 0, ])^2)))
  braf_differences[[as.character(dose_value)]] <- data.table(dose = dose_value, program = programs,
    mutant_mean = colMeans(dd[flag == 1, ]), other_mean = colMeans(dd[flag == 0, ]), difference = delta)
  braf_profiles[[as.character(dose_value)]] <- data.table(dose = dose_value,
    cell = rep(panel$cell_map$depmap_name[match(md$depmap_id, panel$cell_map$depmap_id)], ncol(dd)),
    BRAF = rep(flag, ncol(dd)), program = rep(programs, each = nrow(dd)), displacement = as.vector(dd))
}
write_table(rbindlist(braf_summary), 'braf_displacement_summary')
write_table(rbindlist(braf_differences), 'braf_displacement_differences')
write_table(rbindlist(braf_profiles), 'braf_displacement_profiles')

# Reuse exact prior failure-case IDs, comparing the two new actual-response
# decoders with the prior forty-baseline decoder in native assay units.
cases <- fread(file.path(parent, 'data/contextual_audit_cases.tsv'))[case %in% c('C01', 'C17', 'C19', 'C25')]
original_cases <- fread('/work/results/program_pilot/20_multidrug_baseline_response_function_MEK_holdout/data/diagnostic_geometry_example_resolution.tsv')
original_cases[, case := sprintf('C%02d', seq_len(.N))]
case_values <- rbindlist(lapply(seq_len(nrow(cases)), function(k) {
  pair <- cases[k]
  pair_ids <- original_cases[case == pair$case, c(observation_i, observation_j)]
  rows <- evaluation[observation_id %in% pair_ids & model %in% c('previous_context_actual', 'pca_actual', 'pca_mutation_actual')]
  rows[, `:=`(case = pair$case, endpoint = match(observation_id, pair_ids))]
  rows
}))
stopifnot(nrow(case_values) == 24L)
write_table(case_values, 'failure_cases')

# Compact figures preserve raw ranges and the two distinct prediction axes.
plot_save <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 140)
}
model_labels <- c(previous_anchored_actual = 'Response SVR', previous_context_actual = 'Previous 40 baseline',
  pca_actual = 'Compressed baseline', pca_mutation_actual = 'Compressed + mutations')
main <- evaluation[model %in% names(model_labels)]
main[, model_label := factor(model_labels[model], levels = model_labels)]
focal <- c('lapatinib', 'crizotinib', 'alpelisib', 'trametinib', 'gefitinib', 'everolimus', 'bortezomib', 'vemurafenib')
pdf(file.path(root, 'figures/focal_drugs.pdf'), width = 11.5, height = 8)
for (drug_name in focal) {
  rows <- main[drug == drug_name]
  annotation <- rows[, .(label = sprintf('r=%.2f; sensitive %d/%d', cor_safe(observed_pct, predicted_pct),
    sum(observed_pct <= 50 & predicted_pct <= 50), sum(observed_pct <= 50))), by = .(dose, model_label)]
  p <- ggplot(rows, aes(observed_pct, predicted_pct, color = cell)) + geom_point(size = 1.4) +
    geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
    geom_hline(yintercept = 50, linetype = 2, color = 'grey65') + geom_vline(xintercept = 50, linetype = 2, color = 'grey65') +
    facet_grid(dose ~ model_label, scales = 'free') + scale_color_brewer(palette = 'Set1') +
    geom_text(data = annotation, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = 1.03, vjust = 1.2, size = 2) + theme_classic(base_size = 9) +
    theme(aspect.ratio = 1, legend.position = 'bottom') +
    labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)', color = 'Held cell',
      title = paste('Version 22 |', family, '|', drug_name),
      subtitle = 'Rows: dose (uM); all transformations fitted in training recipients',
      caption = 'Previously inspected outcomes: retrospective evaluation. Dashed lines: 50% vehicle.')
  print(p)
  if (drug_name %in% c('lapatinib', 'vemurafenib')) plot_save(p, paste0('focal_', drug_name), 11.5, 8)
}
dev.off()
p <- ggplot(main, aes(yp, pp, color = cell)) + geom_point(size = .6, alpha = .3) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .6) +
  facet_grid(drug_set ~ model_label, scales = 'free') + scale_color_brewer(palette = 'Set1') +
  theme_classic(base_size = 9) + theme(aspect.ratio = 1, legend.position = 'bottom') +
  labs(x = 'Observed difference from drug/dose mean (percentage points)',
    y = 'Predicted difference from drug/dose mean (percentage points)', color = 'Held cell',
    title = 'Across cells receiving the same drug and dose')
plot_save(p, 'within_treatment', 11.5, 6.5)
cv_long <- melt(cv_axes, id.vars = c('variant', 'lambda'), measure.vars = c('across_cells_r', 'across_drugs_r'))
p <- ggplot(cv_long, aes(lambda, value, color = variant)) + geom_line() + geom_point(size = 1.5) +
  facet_wrap(~variable, scales = 'free_y') + scale_x_log10() + theme_classic(base_size = 10) +
  theme(legend.position = 'bottom') + labs(x = 'Penalty (larger = more shrinkage)', y = 'Development Pearson r',
    color = NULL, title = 'Whole-cell development validation; no held-cell penalty screening')
plot_save(p, 'development_penalty', 8.5, 4.5)
pc_table <- fread(file.path(root, 'data/pca_variance.tsv'))
p <- ggplot(pc_table[component <= 15], aes(component, 100 * cumulative_variance)) + geom_line() +
  geom_point(aes(color = retained), size = 2) + geom_hline(yintercept = 90, linetype = 2) +
  theme_classic(base_size = 10) + labs(x = 'Baseline components', y = 'Cumulative baseline variance (%)',
    color = 'Retained', title = paste('Version 22 |', family, '| baseline compression'))
plot_save(p, 'baseline_compression', 5.2, 4)
case_values[, model_label := factor(model_labels[model], levels = model_labels)]
case_plot <- copy(case_values)
case_plot[, measurement := as.character(model_label)]
measured <- unique(case_plot[, .(case, endpoint, cell, observed_pct)])
measured[, `:=`(measurement = 'Measured PRISM', predicted_pct = observed_pct)]
case_plot <- rbindlist(list(case_plot, measured), fill = TRUE)
case_plot[, measurement := factor(measurement, levels = c('Measured PRISM', 'Previous 40 baseline',
  'Compressed baseline', 'Compressed + mutations'))]
case_labels <- c(C01 = 'CFPAC-1: anastrozole / ixazomib (0.5 uM)',
  C17 = 'Hexylresorcinol: HT-29 / Hs 766T',
  C19 = 'Raltitrexed: HT-29 / BT-474 (0.05 uM)',
  C25 = 'Vemurafenib: LOX IMVI / NCI-H460 (0.5 uM)')
case_plot[, sample_label := ifelse(endpoint == 1L, 'First named sample', 'Second named sample')]
p <- ggplot(case_plot, aes(measurement, predicted_pct, fill = sample_label)) +
  geom_col(position = 'dodge', width = .75) + geom_hline(yintercept = 50, linetype = 2, color = 'grey65') +
  facet_wrap(~case, ncol = 2, scales = 'free_y', labeller = labeller(case = case_labels)) +
  scale_fill_manual(values = c('First named sample' = '#2166AC', 'Second named sample' = '#B2182B')) +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = 'bottom') +
  labs(x = NULL, y = 'Abundance (% vehicle)', fill = NULL,
    title = 'Exact prior failure cases: measured and predicted',
    subtitle = 'Colors identify the first and second sample named in each panel')
plot_save(p, 'failure_cases', 10, 7)
braf <- rbindlist(braf_profiles)
braf[, genotype := ifelse(BRAF == 1, cell, 'Other cells')]
p <- ggplot(braf, aes(match(program, programs), displacement, group = interaction(cell, dose), color = genotype)) +
  geom_line(aes(alpha = BRAF == 1), linewidth = .6) + facet_wrap(~dose, nrow = 1, scales = 'free_y') +
  scale_alpha_manual(values = c('FALSE' = .12, 'TRUE' = 1), guide = 'none') +
  theme_classic(base_size = 9) + theme(legend.position = 'bottom') +
  labs(x = 'Frozen program (P01-P40)', y = 'Raw displacement', color = 'BRAF-positive cell',
    title = 'Vemurafenib displacement: BRAF-positive cells versus other cells',
    subtitle = 'All measured cells; no profile normalization; genotype comparison is retrospective and lineage-confounded')
plot_save(p, 'braf_displacement', 11, 4.2)
pdf(file.path(root, 'figures/interaction_weights.pdf'), width = 9, height = 6)
for (name in variants$name) {
  tab <- coefficient_tables[[name]]
  features <- unique(tab$baseline_feature)
  mm <- matrix(NA_real_, length(features), length(programs), dimnames = list(features, programs))
  mm[cbind(match(tab$baseline_feature, features), match(tab$displaced_program, programs))] <- tab$coefficient
  ro <- hclust(dist(mm), method = 'ward.D2')$order
  co <- hclust(dist(t(mm)), method = 'ward.D2')$order
  limit <- max(abs(mm))
  if (limit == 0) limit <- 1
  hm <- draw_heatmap(mm[ro, co, drop = FALSE], scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE, value_range = c(-limit, 0, limit),
    value_colors = c('#2166AC', '#FFFFFF', '#B2182B'), hm_name = 'Log2 effect',
    title = paste(name, '| baseline feature by displaced program'), row_names_gp = grid::gpar(fontsize = 8))
  ComplexHeatmap::draw(hm)
}
dev.off()
status('COMPLETE')
print(rbindlist(selection))
print(summary)
print(rbindlist(braf_summary))
