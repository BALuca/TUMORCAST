#!/usr/bin/env Rscript
# Diagnose the current decoder's training signal, objective and linear fitting capacity.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet)})
source('/work/src/program_pilot/response_models.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/37_training_program_signal_capture_audit'
score_root <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps'
decoder_root <- '/work/results/program_pilot/33_balanced_cell_and_dose_contrast_decoder'
families <- c('ERBB', 'MEK')
dir.create(file.path(root, 'data'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, 'models'), showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})

# Subtract context means for each column; preserve observation order exactly.
center_context <- function(mat, groups) {
  mat <- as.matrix(mat)
  groups <- factor(groups)
  means <- rowsum(mat, groups, reorder = TRUE) / as.numeric(table(groups))
  mat - means[as.integer(groups), , drop = FALSE]
}
summaries <- program_tables <- coefficient_tables <- list()
sources <- c('/work/src/program_pilot/37_training_program_signal_capture_audit.R',
  '/work/src/program_pilot/response_models.R')
for (family_name in families) {
  status(paste(family_name, 'auditing the exact training matrix'))
  source_path <- file.path(score_root, 'data', paste0(family_name, '_all_training_inputs.qs'))
  model_path <- file.path(decoder_root, 'models', paste0(family_name, '_balanced_decoder.qs'))
  prepared_path <- file.path(paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_',
    family_name, '_holdout'), 'data/prepared_scores.qs')
  sources <- c(sources, source_path, model_path, prepared_path)
  input <- qread(source_path)
  decoder <- qread(model_path)
  prepared <- qread(prepared_path)
  md <- copy(input$annotations)
  d <- input$displacement
  programs <- colnames(d)
  y <- log2(md$measured_pct / 100)
  b <- sweep(sweep(prepared$baseline[md$depmap_id, , drop = FALSE], 2L, decoder$baseline_mean, '-'),
    2L, decoder$baseline_sd, '/')
  x <- cbind(d, d * b)
  colnames(x) <- c(paste0('D_', programs), paste0('DB_', programs))
  stopifnot(ncol(d) == 40L, identical(rownames(d), md$observation_id), all(is.finite(x)),
    !any(md$depmap_id %in% decoder$excluded_ids), !any(md$drug %in% decoder$excluded_drugs))
  current <- predict_ridge(decoder$model, x)
  stopifnot(max(abs(current - log2(md$predicted_pct / 100))) < 1e-10)
  weak <- as.numeric(pred_path(decoder$model, x, s = 1e-6))
  eligible <- which(md$condition %in% decoder$screen$cell_keys)
  cy <- as.numeric(center_context(matrix(y[eligible], ncol = 1), md$condition[eligible]))
  cd <- center_context(d[eligible, , drop = FALSE], md$condition[eligible])
  cx <- center_context(x[eligible, , drop = FALSE], md$condition[eligible])
  within_weights <- 1 / as.numeric(table(md$condition[eligible]))[match(md$condition[eligible],
    names(table(md$condition[eligible])))]
  global_weights <- observation_weights(md)
  predictions <- list('Current selected ridge' = current, 'Same objective, weak ridge' = weak)
  fitting_rank <- c('Current selected ridge' = NA_integer_, 'Same objective, weak ridge' = NA_integer_)

  # Ordinary QR linear fits are training-only diagnostics, not selected predictors.
  full_designs <- list('Absolute linear D' = cbind(Intercept = 1, d),
    'Absolute linear D + DB' = cbind(Intercept = 1, x))
  for (name in names(full_designs)) {
    fit <- lm.wfit(full_designs[[name]], y, w = global_weights)
    predictions[[name]] <- fit$fitted.values
    fitting_rank[name] <- fit$rank
    coefficient_tables[[paste(family_name, name)]] <- data.table(family = family_name, model = name,
      feature = colnames(full_designs[[name]]), coefficient = fit$coefficients)
    qsave(fit, file.path(root, 'models', paste0(family_name, '_', gsub('[^a-z]+', '_', tolower(name)), '.qs')))
  }
  fixed_designs <- list('Cell-only centered D' = cd, 'Cell-only centered D + DB' = cx)
  for (name in names(fixed_designs)) {
    fit <- lm.wfit(fixed_designs[[name]], cy, w = within_weights)
    coefficient <- fit$coefficients
    coefficient[is.na(coefficient)] <- 0
    predictions[[name]] <- rep(NA_real_, nrow(md))
    predictions[[name]][eligible] <- as.numeric(fixed_designs[[name]] %*% coefficient)
    fitting_rank[name] <- fit$rank
    coefficient_tables[[paste(family_name, name)]] <- data.table(family = family_name, model = name,
      feature = colnames(fixed_designs[[name]]), coefficient = fit$coefficients)
    qsave(fit, file.path(root, 'models', paste0(family_name, '_', gsub('[^a-z]+', '_', tolower(name)), '.qs')))
  }
  for (name in names(predictions)) {
    p <- predictions[[name]]
    stopifnot(all(is.finite(p[eligible])))
    centered <- as.numeric(center_context(matrix(p[eligible], ncol = 1), md$condition[eligible]))
    global <- all(is.finite(p))
    summaries[[paste(family_name, name)]] <- data.table(family = family_name, model = name,
      global_observations = if (global) nrow(md) else 0L, global_r = if (global) cor_safe(y, p) else NA_real_,
      fixed_context_observations = length(eligible), fixed_contexts = uniqueN(md$condition[eligible]),
      fixed_context_cell_r = cor_safe(cy, centered), fitted_rank = fitting_rank[name],
      measured_sensitive_median_pct = median(md$measured_pct[md$measured_pct <= 50]),
      measured_resistant_median_pct = median(md$measured_pct[md$measured_pct > 50]),
      predicted_sensitive_median_pct = if (global) median(100 * 2^p[md$measured_pct <= 50]) else NA_real_,
      predicted_resistant_median_pct = if (global) median(100 * 2^p[md$measured_pct > 50]) else NA_real_)
  }
  program_tables[[family_name]] <- rbindlist(lapply(seq_along(programs), function(j) {
    z <- d[, j] / input$program_sd[j]
    context_r <- vapply(split(eligible, md$condition[eligible]), function(ii) cor_safe(d[ii, j], y[ii]), numeric(1))
    data.table(family = family_name, program = programs[j], global_r = cor_safe(d[, j], y),
      fixed_context_cell_r = cor_safe(cd[, j], cy),
      resistant_minus_sensitive_training_sd = mean(z[md$measured_pct > 50]) - mean(z[md$measured_pct <= 50]),
      positive_contexts = sum(context_r > 0, na.rm = TRUE), negative_contexts = sum(context_r < 0, na.rm = TRUE),
      evaluable_contexts = sum(is.finite(context_r)))
  }))
  qsave(list(observation_id = md$observation_id, eligible_rows = eligible, predictions = predictions),
    file.path(root, 'data', paste0(family_name, '_diagnostic_predictions.qs')))
}
summary <- rbindlist(summaries)
save_table(summary, 'training_readout_comparison')
save_table(rbindlist(program_tables), 'individual_program_signal')
save_table(rbindlist(coefficient_tables), 'diagnostic_linear_coefficients')
save_table(data.table(path = unique(sources), sha256 = vapply(unique(sources), digest::digest,
  character(1), file = TRUE, algo = 'sha256')), 'input_hashes')
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(c('TRAINING SIGNAL CAPTURE AUDIT', '',
  'Global r: all training observations, combining drug/dose/cell differences.',
  'Fixed-context cell r: same 34/29 informative contexts used for script-33 r=0.165/0.283.',
  'Within each fixed context, both outcome and prediction are centered before pooled correlation.',
  'Simple linear models are unpenalized training capacity/objective diagnostics, not independent predictors.',
  'Cell-only centered models use observed training context means and provide no absolute predictions.',
  'Changing losses, intercepts and penalties is explicitly identified; this is not a selected new architecture.',
  'No outer recipients or family-excluded drugs were read for evaluation.'), file.path(root, 'READING_GUIDE.txt'))
print(summary)
status('COMPLETE')
