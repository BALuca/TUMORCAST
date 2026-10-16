#!/usr/bin/env Rscript
# A low-complexity, robust displacement-to-function map using the frozen 40 programs.
suppressPackageStartupMessages({library(data.table); library(qs); library(mgcv); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

root <- '/work/results/program_pilot/39_constrained_additive_sensitivity_decoder'
score_path <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps/data/ERBB_all_training_inputs.qs'
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout'
decoder_root <- '/work/results/program_pilot/33_balanced_cell_and_dose_contrast_decoder'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
curve_basis_size <- 4L # One zero-point constraint leaves three coefficients per program.
baseline_penalty_multiplier <- 20
absolute_mass <- .5
for (path in file.path(root, c('data', 'models', 'figures'))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(root, 'figures', paste0(stem, '.pdf')), plot,
    device = grDevices::cairo_pdf, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(stem, '.png')), plot, width = width, height = height, dpi = 180)
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})

status('loading unchanged ERBB training programs')
input <- qread(score_path)
prepared_path <- file.path(parent, 'data/prepared_scores.qs')
prepared <- qread(prepared_path)
ridge_path <- file.path(decoder_root, 'models/ERBB_balanced_decoder.qs')
ridge <- qread(ridge_path)
fold_paths <- file.path(decoder_root, 'models', paste0('ERBB_fold_', 1:5, '.qs'))
folds <- lapply(fold_paths, qread)
md <- copy(input$annotations)
d <- input$displacement
programs <- colnames(d)
y <- log2(md$measured_pct / 100)
stopifnot(nrow(md) == 6260L, ncol(d) == 40L, !anyDuplicated(md$observation_id),
  identical(rownames(d), md$observation_id), identical(prepared$signature, ridge$source_signature),
  all(is.finite(d)), all(is.finite(y)), uniqueN(md$depmap_id) == 26L,
  !any(md$depmap_id %in% ridge$excluded_ids), !any(md$drug %in% ridge$excluded_drugs))
sources <- c(score_path, prepared_path, ridge_path, fold_paths,
  '/work/src/program_pilot/39_constrained_additive_sensitivity_decoder.R',
  '/work/src/program_pilot/response_models.R')
manifest <- data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256'))
signature <- digest::digest(manifest, algo = 'sha256')
save_table(manifest, 'input_hashes')
qsave(list(signature = signature, source_signature = prepared$signature, programs = programs,
  training_ids = ridge$training_ids, excluded_ids = ridge$excluded_ids, excluded_drugs = ridge$excluded_drugs,
  curve_basis_size = curve_basis_size, baseline_penalty_multiplier = baseline_penalty_multiplier,
  absolute_mass = absolute_mass, family = 'scaled-t', smoothing = 'REML',
  objective = 'absolute abundance and centered same-treatment abundance; equal target-energy-normalized masses',
  anchoring = 'each spline and baseline interaction equals zero at zero displacement'), file.path(root, 'settings.qs'))

# Reused for full/fold fits and arbitrary future molecular queries. Knots and all
# scales are learned from fitting recipients only. The baseline term has no intercept.
make_design <- function(ids, displacement, basis) {
  b <- sweep(sweep(prepared$baseline[ids, programs, drop = FALSE], 2L, basis$baseline_mean, '-'),
    2L, basis$baseline_sd, '/')
  main <- do.call(cbind, lapply(seq_along(programs), function(j) {
    raw <- PredictMat(basis$smooths[[j]], data.frame(value = displacement[, j]))
    sweep(raw, 2L, basis$curve_sd[[j]], '/')
  }))
  interaction <- sweep(displacement * b, 2L, basis$interaction_sd, '/')
  x <- cbind(main, interaction)
  colnames(x) <- basis$feature_names
  stopifnot(ncol(x) == 160L, all(is.finite(x)))
  x
}

fit_additive <- function(rows, context_keys) {
  dd <- d[rows, , drop = FALSE]
  metadata <- md[rows]
  ids <- sort(unique(metadata$depmap_id))
  bb <- prepared$baseline[ids, programs, drop = FALSE]
  basis <- list(baseline_mean = colMeans(bb), baseline_sd = apply(bb, 2L, sd),
    smooths = vector('list', 40L), curve_sd = vector('list', 40L))
  stopifnot(all(basis$baseline_sd > 0))
  curvature <- matrix(0, 160L, 160L)
  for (j in seq_along(programs)) {
    sm <- smoothCon(s(value, bs = 'cr', k = curve_basis_size, pc = list(value = 0)),
      data.frame(value = dd[, j]), absorb.cons = TRUE)[[1L]]
    sc <- apply(sm$X, 2L, sd)
    stopifnot(ncol(sm$X) == 3L, all(sc > 0))
    penalty <- sm$S[[1L]] / outer(sc, sc)
    # Equalize the curvature penalty units across the 40 distinct bases.
    penalty <- penalty / sum(diag(penalty))
    index <- (3L * j - 2L):(3L * j)
    curvature[index, index] <- penalty
    basis$smooths[[j]] <- sm
    basis$curve_sd[[j]] <- sc
  }
  bz <- sweep(sweep(prepared$baseline[metadata$depmap_id, programs, drop = FALSE],
    2L, basis$baseline_mean, '-'), 2L, basis$baseline_sd, '/')
  basis$interaction_sd <- apply(dd * bz, 2L, sd)
  basis$feature_names <- c(unlist(lapply(programs, function(p) paste0(p, '_curve_', 1:3))), paste0(programs, '_DB'))
  stopifnot(all(basis$interaction_sd > 0))
  x_absolute <- make_design(metadata$depmap_id, dd, basis)
  yy <- y[rows]
  informative <- which(metadata$condition %in% context_keys)
  stopifnot(length(informative) > 0L)
  cell_metadata <- metadata[informative]
  x_cell <- x_absolute[informative, , drop = FALSE]
  y_cell <- yy[informative]
  # Center the evaluated curves, not their displacement arguments. This is the
  # same f(B,D) on both axes and permits nonlinear thresholds without changing units.
  for (key in unique(cell_metadata$condition)) {
    ii <- which(cell_metadata$condition == key)
    x_cell[ii, ] <- sweep(x_cell[ii, , drop = FALSE], 2L, colMeans(x_cell[ii, , drop = FALSE]), '-')
    y_cell[ii] <- y_cell[ii] - mean(y_cell[ii])
  }
  wa <- observation_weights(metadata)
  wa <- wa / sum(wa)
  cell_metadata[, context_n := .N, by = condition]
  wc <- 1 / cell_metadata$context_n
  wc <- wc / sum(wc)
  energy_absolute <- sum(wa * yy^2)
  energy_cell <- sum(wc * y_cell^2)
  stopifnot(energy_absolute > 0, energy_cell > 0)
  X <- rbind(x_absolute, x_cell)
  target <- c(yy, y_cell)
  weights <- c(absolute_mass * wa / energy_absolute, (1 - absolute_mass) * wc / energy_cell)
  weights <- weights / mean(weights)
  # Two shared regularizers: amplitude and curvature. Baseline deviations share
  # the amplitude penalty but pay 20 times more, rather than choosing 40 free fits.
  amplitude <- diag(c(rep(1, 120L), rep(baseline_penalty_multiplier, 40L)))
  fit <- gam(target ~ X - 1, weights = weights, family = scat(), method = 'REML',
    paraPen = list(X = list(amplitude, curvature)), control = gam.control(nthreads = 1L))
  stopifnot(all(is.finite(coef(fit))), isTRUE(fit$converged), all(fit$sp > 0))
  list(model = fit, basis = basis, training_ids = ids, context_keys = context_keys,
    objective_energy = c(absolute = energy_absolute, cell = energy_cell),
    effective_df = sum(fit$edf), signature = signature, source_signature = prepared$signature)
}
predict_additive <- function(fitted, ids, displacement) {
  drop(make_design(ids, displacement, fitted$basis) %*% coef(fitted$model))
}
summarize <- function(metadata, prediction, context_keys, cohort) {
  z <- copy(metadata)
  z[, `:=`(observed = log2(measured_pct / 100), predicted = prediction, predicted_pct = 100 * 2^prediction)]
  cell <- z[condition %in% context_keys]
  cell[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)), by = condition]
  data.table(cohort = cohort, observations = nrow(z), cells = uniqueN(z$depmap_id),
    global_r = cor_safe(z$observed, z$predicted), fixed_treatment_observations = nrow(cell),
    fixed_treatment_r = cor_safe(cell$yc, cell$pc), sensitive = sum(z$measured_pct <= 50),
    sensitive_called = sum(z$measured_pct <= 50 & z$predicted_pct <= 50),
    resistant = sum(z$measured_pct > 50), false_sensitive = sum(z$measured_pct > 50 & z$predicted_pct <= 50),
    predicted_sensitive_median_pct = median(z$predicted_pct[z$measured_pct <= 50]),
    predicted_resistant_median_pct = median(z$predicted_pct[z$measured_pct > 50]))
}

status('fitting one constrained robust additive decoder')
full <- fit_additive(seq_len(nrow(md)), ridge$screen$cell_keys)
qsave(full, file.path(root, 'models', 'constrained_decoder.qs'))
prediction <- predict_additive(full, md$depmap_id, d)
zero <- predict_additive(full, full$training_ids, matrix(0, length(full$training_ids), 40L,
  dimnames = list(full$training_ids, programs)))
stopifnot(max(abs(zero)) < 1e-10,
  max(abs(predict_additive(qread(file.path(root, 'models', 'constrained_decoder.qs')), md$depmap_id, d) - prediction)) < 1e-10)
save_table(data.table(depmap_id = full$training_ids, zero_displacement_pct_vehicle = 100 * 2^zero), 'zero_input_audit')
save_table(data.table(feature = names(coef(full$model)), coefficient = coef(full$model), effective_df = full$model$edf),
  'coefficients_and_effective_df')
save_table(data.table(total_effective_df = full$effective_df, maximum_coefficients = 160L,
  amplitude_penalty = full$model$sp[1L], curvature_penalty = full$model$sp[2L],
  baseline_penalty_multiplier = baseline_penalty_multiplier,
  robust_t_df = full$model$family$getTheta(TRUE)[1L], robust_scale = full$model$family$getTheta(TRUE)[2L]), 'fit_constraints')
training <- copy(md)
training[, `:=`(predicted = prediction, predicted_pct = 100 * 2^prediction, observed = y)]
save_table(training, 'training_predictions')
metrics <- list(summarize(md, prediction, ridge$screen$cell_keys, 'Training fit'))
save_table(rbindlist(metrics), 'readout_summary')
status(paste('training fit saved | effective df=', round(full$effective_df, 1),
  '| fixed drug/dose r=', round(metrics[[1L]]$fixed_treatment_r, 3)))

# Expose each learned response curve and its baseline-dependent adjustment.
curves <- rbindlist(lapply(seq_along(programs), function(j) {
  values <- seq(quantile(d[, j], .01), quantile(d[, j], .99), length.out = 101L)
  main <- sweep(PredictMat(full$basis$smooths[[j]], data.frame(value = values)), 2L, full$basis$curve_sd[[j]], '/')
  index <- (3L * j - 2L):(3L * j)
  g <- drop(main %*% coef(full$model)[index])
  adjustment <- values / full$basis$interaction_sd[j] * coef(full$model)[120L + j]
  rbindlist(lapply(c(-1, 0, 1), function(b) data.table(program = programs[j], displacement = values,
    baseline_sd = factor(b, levels = c(-1, 0, 1)), log2_abundance_contribution = g + b * adjustment)))
}))
save_table(curves, 'program_response_curves')
p <- ggplot(curves, aes(displacement, log2_abundance_contribution, color = baseline_sd)) +
  geom_hline(yintercept = 0, color = 'grey85') + geom_vline(xintercept = 0, color = 'grey85') +
  geom_line(linewidth = .5) + facet_wrap(~program, ncol = 8, scales = 'free') +
  scale_x_continuous(n.breaks = 3) +
  scale_color_manual(values = c('#0072B2', '#333333', '#D55E00'), labels = c('-1 SD', 'Mean baseline', '+1 SD')) +
  theme_classic(base_size = 9) + theme(legend.position = 'bottom') +
  labs(x = 'Program displacement (hard gene mean)', y = 'Contribution to log2 abundance', color = 'Baseline',
    title = paste0('40 constrained program curves; total effective df = ', round(full$effective_df, 1)),
    caption = 'Curves span fitting 1st–99th displacement percentiles; every curve passes through zero. Baseline terms are strongly penalized.')
save_plot(p, 'program_response_curves', 11, 8)

status('checking the same frozen decoder specification in five whole-cell folds')
oof <- rep(NA_real_, nrow(md))
eligible <- rep(FALSE, nrow(md))
fold_id <- rep(NA_integer_, nrow(md))
for (k in 1:5) {
  saved <- folds[[k]]
  tr <- which(md$depmap_id %in% saved$training_ids)
  va <- which(md$depmap_id %in% saved$validation_ids)
  stopifnot(!length(intersect(md$depmap_id[tr], md$depmap_id[va])), length(tr) + length(va) == nrow(md))
  status(paste('whole-cell fold', k, '/5'))
  fitted <- fit_additive(tr, saved$screen$cell_keys)
  pp <- predict_additive(fitted, md$depmap_id[va], d[va, , drop = FALSE])
  qsave(list(fitted = fitted, validation_ids = saved$validation_ids, validation_rows = va, prediction = pp),
    file.path(root, 'models', paste0('fold_', k, '.qs')))
  oof[va] <- pp
  eligible[va] <- md$condition[va] %in% saved$screen$cell_keys
  fold_id[va] <- k
  metrics[[k + 1L]] <- summarize(md[va], pp, saved$screen$cell_keys, paste0('Cell fold ', k))
  save_table(rbindlist(metrics), 'readout_summary')
  qsave(list(prediction = oof, eligible = eligible, fold = fold_id, observation_id = md$observation_id),
    file.path(root, 'data', 'cell_fold_predictions.qs'))
}
stopifnot(all(is.finite(oof)), !anyNA(fold_id))

# Original excluded recipients never enter the fit. Freeze their predictions
# before accessing their previously revealed functional outcomes.
status('freezing original held-cell predictions')
query <- copy(as.data.table(prepared$metadata)[depmap_id %in% ridge$excluded_ids & drug %in% md$drug])
query[, condition := paste(drug, dose, sep = '__')]
hd <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
query[, predicted := predict_additive(full, depmap_id, hd)]
stopifnot(all(is.finite(query$predicted)))
save_table(query, 'outcome_free_held_predictions')
prediction_path <- file.path(root, 'data', 'outcome_free_held_predictions.tsv')
save_table(data.table(path = prediction_path, sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')
panel <- qread(panel_path)
lookup <- as.data.table(panel$cell_map)
query[, recipient := lookup$depmap_name[match(depmap_id, lookup$depmap_id)]]
stopifnot(!anyNA(query$recipient))
labels <- as.data.table(panel$sensitivity)
query[, measured_pct := labels$prism_abundance_pct[match(observation_id, labels$observation_id)]]
save_table(query[!is.finite(measured_pct) | measured_pct <= 0], 'unsupported_held_endpoints')
query <- query[is.finite(measured_pct) & measured_pct > 0]
stopifnot(uniqueN(query$depmap_id) == 9L, !any(query$depmap_id %in% md$depmap_id))
metrics[[7L]] <- summarize(query, query$predicted, ridge$screen$cell_keys, 'Original held cells')
save_table(rbindlist(metrics), 'readout_summary')

tab <- rbindlist(list(
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Training fit', predicted = prediction, cell_context = md$condition %in% ridge$screen$cell_keys, fold = 0L)),
  cbind(md[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Development cell folds', predicted = oof, cell_context = eligible, fold = fold_id)),
  cbind(query[, .(observation_id, depmap_id, recipient, drug, dose, condition, measured_pct)],
    data.table(cohort = 'Original held cells', predicted = query$predicted,
      cell_context = query$condition %in% ridge$screen$cell_keys, fold = 0L))))
tab[, `:=`(observed = log2(measured_pct / 100), predicted_pct = 100 * 2^predicted)]
tab[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)), by = .(cohort, fold, condition)]
save_table(tab, 'evaluated_predictions')
save_table(tab[cell_context == TRUE, .(observations = .N, cells = uniqueN(depmap_id),
  r = cor_safe(observed, predicted)), by = .(cohort, fold, drug, dose)], 'per_treatment_correlations')
plot_data <- rbindlist(list(
  tab[, .(cohort, drug, x = observed, y = predicted, axis = 'All treatments combined')],
  tab[cell_context == TRUE, .(cohort, drug, x = yc, y = pc, axis = 'Cells at fixed drug/dose')]))
plot_data[, cohort := factor(cohort, levels = c('Training fit', 'Development cell folds', 'Original held cells'))]
plot_data[, axis := factor(axis, levels = c('All treatments combined', 'Cells at fixed drug/dose'))]
cor_labels <- plot_data[, .(label = sprintf('r = %.2f; n = %s', cor_safe(x, y), format(.N, big.mark = ','))), by = .(cohort, axis)]
p <- ggplot(plot_data, aes(x, y)) + geom_point(alpha = .13, size = .45, color = '#0072B2') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#333333', linewidth = .5) +
  geom_text(data = cor_labels, aes(x = -Inf, y = Inf, label = label), hjust = -.03, vjust = 1.2,
    inherit.aes = FALSE, size = 3) + facet_wrap(vars(axis, cohort), ncol = 3, scales = 'free') +
  theme_classic(base_size = 10) + theme(aspect.ratio = 1) +
  labs(x = 'Measured log2 abundance (centered for fixed-treatment panels)',
    y = 'Predicted log2 abundance (centered for fixed-treatment panels)', title = 'Constrained additive molecular decoder',
    caption = 'Development folds condition on the frozen outcome-informed bank; original nine held cells are retrospective. No drug/cell identity inputs.')
save_plot(p, 'measured_predicted_readout', 10, 7.3)
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(c('CONSTRAINED ADDITIVE DECODER', '',
  'Same 40 hard-mean displacement programs and exact 6,260 fitting observations as script 35.',
  '40 anchored cubic curves, three coefficients each; 40 linear matched-baseline interactions.',
  'Two shared REML penalties; standardized baseline interaction coefficients pay 20 times more amplitude penalty.',
  'Robust scaled-t likelihood. Half absolute and half same-treatment centered losses, normalized by fitting target energy.',
  'Contexts with weak/no clear cross-cell sensitivity variation do not enter the centered objective.',
  'No free baseline main effect, identity predictor, unrestricted interactions, parameter grid or zero subtraction.',
  'Zero displacement gives exactly 100% vehicle response at every baseline.',
  'Save each model immediately. Same five whole-cell folds; all scales/knots learned from fitting recipients.',
  'Discovery was already frozen using development outcomes; cell folds evaluate the readout conditional on that bank.',
  'Original held outcomes are retrospective and unused for selection; predictions hashed before outcome matching.',
  'This iteration isolates measured molecular displacement-to-function fitting. Baseline-to-displacement fitting is unchanged.'),
  file.path(root, 'READING_GUIDE.txt'))
print(rbindlist(metrics))
status('COMPLETE')
