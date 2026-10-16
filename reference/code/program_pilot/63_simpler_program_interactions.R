#!/usr/bin/env Rscript
# A constrained grant model:40 same-program interaction terms. The single-effect
# decoder is frozen from62; all other fitting and exclusion rules are unchanged.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(digest)
})

root <- '/work/results/program_pilot/63_simpler_program_interactions'
previous <- '/work/results/program_pilot/62_separate_response_interaction_bridge'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data'
lambda_grid <- 10^seq(3, -4, length.out = 24L)
rna_times <- c(3, 24)
set.seed(63011L)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}

inputs <- c(file.path(prior, c('matched_molecular_links.tsv.gz', 'representation_scores.qs')),
  file.path(previous, 'models', paste0(rna_times, 'h_Original40_single_decoder.qs')),
  file.path(previous, 'data/predictions.tsv.gz'))
stopifnot(all(file.exists(inputs)))
links <- fread(inputs[1])
bank <- qread(inputs[2])
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1),
  algo = 'sha256', file = TRUE)), 'input_manifest')
stopifnot(ncol(bank$original40) == 40L, all(is.finite(bank$original40)))

# Every feature transformation preserves zero displacement. Use the same bank.
scores <- list(Original40 = bank$original40)
qsave(list(lambda_grid = lambda_grid, primary_time_h = 3, secondary_time_h = 24,
  contrast_objective_fraction = .5, outer = 'one partner primary; two partners matched-contrast diagnostic',
  single_measurements_available = TRUE, original_weight_sha256 = bank$original_weight_sha256),
  file.path(root, 'models/settings.qs'))

# The symmetric product map is reused for molecular and functional controls.
cross_features <- function(a, b) {
  stopifnot(identical(dim(a), dim(b)), identical(colnames(a), colnames(b)))
  ij <- which(upper.tri(matrix(0, ncol(a), ncol(a)), diag = TRUE), arr.ind = TRUE)
  x <- vapply(seq_len(nrow(ij)), function(k) {
    i <- ij[k, 1]
    j <- ij[k, 2]
    (a[, i] * b[, j] + a[, j] * b[, i]) / 2
  }, numeric(nrow(a)))
  colnames(x) <- paste(colnames(a)[ij[, 1]], colnames(a)[ij[, 2]], sep = '__')
  x
}
balanced_weights <- function(md) {
  group <- interaction(md$anchor, md$partner, md$biological_replicate, drop = TRUE)
  w <- 1 / as.numeric(table(group)[as.character(group)])
  w / sum(w)
}
center_context <- function(x, md) {
  group <- interaction(md$anchor, md$nominal_A_uM, md$biological_replicate, drop = TRUE)
  result <- x
  for (idx in split(seq_len(nrow(md)), group)) {
    if (is.matrix(x)) result[idx, ] <- sweep(x[idx, , drop = FALSE], 2, colMeans(x[idx, , drop = FALSE]), '-')
    else result[idx] <- x[idx] - mean(x[idx])
  }
  result
}

# Ridge scaling uses RMS around biological zero, with no intercept. Explicit
# training contrasts change the objective, not the prediction transformation.
fit_path <- function(x, y, weights, md = NULL, block = rep('response', ncol(x))) {
  rms <- sqrt(colSums(x^2 * weights) / sum(weights))
  keep <- is.finite(rms) & rms > 1e-10
  stopifnot(sum(keep) >= 2L)
  xx <- sweep(x[, keep, drop = FALSE], 2, rms[keep], '/')
  penalty <- as.numeric(table(block[keep])[block[keep]])
  if (!is.null(md)) {
    xx <- rbind(xx, center_context(xx, md))
    y <- c(y, center_context(y, md))
    weights <- c(weights * .5, weights * .5)
  }
  fit <- glmnet(xx, y, weights = weights, alpha = 0, lambda = lambda_grid,
    standardize = FALSE, intercept = FALSE, penalty.factor = penalty,
    thresh = 1e-8, maxit = 100000L)
  list(fit = fit, keep = keep, rms = rms[keep], columns = colnames(x), block = block)
}
predict_path <- function(fit, x) {
  stopifnot(identical(colnames(x), fit$columns))
  xx <- sweep(x[, fit$keep, drop = FALSE], 2, fit$rms, '/')
  as.matrix(predict(fit$fit, newx = xx, s = lambda_grid))
}
objective <- function(y, prediction, md) {
  w <- balanced_weights(md)
  absolute <- colSums((prediction - y)^2 * w)
  centered <- colSums((center_context(prediction, md) - center_context(y, md))^2 * w)
  .5 * absolute + .5 * centered
}
safe_cor <- function(a, b) {
  if (length(a) < 4L || sd(a) < 1e-10 || sd(b) < 1e-10) return(NA_real_)
  cor(a, b)
}

all_predictions <- list()
all_tuning <- list()
all_coefficients <- list()
all_audits <- list()
prediction_counter <- 0L
previous_predictions <- fread(file.path(previous, 'data/predictions.tsv.gz'))
for (tt in rna_times) {
  obs <- copy(links[time_h == tt & molecular_supported == TRUE])
  setorder(obs, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)
  stopifnot(nrow(obs) == 496L, uniqueN(obs$condition_id) == 198L,
    !anyDuplicated(obs$observation_id), all(obs$ratio_AB > 0), all(obs$ratio_A > 0), all(obs$ratio_B > 0))
  partners <- sort(unique(obs$partner))
  held_sets <- c(lapply(partners, function(x) x), combn(partners, 2, simplify = FALSE))

  # Exact reuse isolates interaction complexity from single-effect calibration.
  frozen <- qread(file.path(previous, 'models', paste0(tt, 'h_Original40_single_decoder.qs')))
  f_predictions <- list(Original40 = setNames(
    predict_path(frozen$fit, bank$original40)[, frozen$selected], rownames(bank$original40)))

  # All control features are single-input quantities. They are symmetric under
  # swapping drugs and vanish for an absent (dose0, relative count1) constituent.
  control_a <- cbind(effect = log2(obs$ratio_A), exposure = log1p(obs$nominal_A_uM),
    effect_squared = log2(obs$ratio_A)^2)
  control_b <- cbind(effect = log2(obs$ratio_B), exposure = log1p(obs$nominal_B_uM),
    effect_squared = log2(obs$ratio_B)^2)
  control_x <- cross_features(control_a, control_b)
  colnames(control_x) <- paste0('functional_', colnames(control_x))
  additive_observed <- log2(obs$ratio_A) + log2(obs$ratio_B)
  designs <- offsets <- blocks <- list()
  for (representation in names(scores)) {
    a <- scores[[representation]][obs$profile_A, , drop = FALSE]
    b <- scores[[representation]][obs$profile_B, , drop = FALSE]
    xx <- a * b
    colnames(xx) <- paste(colnames(a), colnames(a), sep = '__')
    symmetry_error <- max(abs(xx - b * a))
    zero_error <- max(abs((a * 0) * b))
    stopifnot(ncol(xx) == 40L)
    stopifnot(symmetry_error < 1e-12, zero_error == 0)
    conditioned <- paste0(representation, '_interaction')
    bridge <- paste0(representation, '_bridge')
    designs[[conditioned]] <- cbind(control_x, xx)
    offsets[[conditioned]] <- additive_observed
    blocks[[conditioned]] <- c(rep('functional', ncol(control_x)), rep('molecular', ncol(xx)))
    designs[[bridge]] <- xx
    offsets[[bridge]] <- f_predictions[[representation]][obs$profile_A] +
      f_predictions[[representation]][obs$profile_B]
    blocks[[bridge]] <- rep('molecular', ncol(xx))
  }
  old_offset <- unique(previous_predictions[time_h == tt & model == 'Original40_bridge',
    .(observation_id, predicted_single_sum_log2)])
  stopifnot(!anyDuplicated(old_offset$observation_id))
  offset_error <- max(abs(offsets$Original40_bridge -
    old_offset$predicted_single_sum_log2[match(obs$observation_id, old_offset$observation_id)]))
  stopifnot(is.finite(offset_error), offset_error < 1e-12)
  # Only y[training] is supplied to interaction fitting and nested selection.
  y <- log2(obs$ratio_AB)
  for (held in held_sets) {
    fold_id <- paste(held, collapse = '__')
    fold_type <- if (length(held) == 1L) 'single_partner' else 'two_partners'
    tr <- which(!obs$partner %chin% held)
    te <- which(obs$partner %chin% held)
    training_partners <- setdiff(partners, held)
    source_ids <- function(x) unique(unlist(strsplit(x, ';', fixed = TRUE)))
    overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]), source_ids(obs$source_sample_ids_AB[te]))
    stopifnot(length(overlap) == 0L, !any(obs$partner[tr] %chin% held))
    for (model_name in names(designs)) {
      status(paste(tt, 'h | exclude', paste(held, collapse = ' + '), '|', model_name))
      x <- designs[[model_name]]
      offset <- offsets[[model_name]]
      target <- y - offset
      inner_loss <- matrix(NA_real_, length(training_partners), length(lambda_grid) + 1L)
      for (j in seq_along(training_partners)) {
        ii <- tr[obs$partner[tr] != training_partners[j]]
        vv <- tr[obs$partner[tr] == training_partners[j]]
        fit <- fit_path(x[ii, , drop = FALSE], target[ii], balanced_weights(obs[ii]),
          md = obs[ii], block = blocks[[model_name]])
        pred <- cbind(predict_path(fit, x[vv, , drop = FALSE]), zero = 0)
        inner_loss[j, ] <- objective(target[vv], pred, obs[vv])
      }
      mean_loss <- colMeans(inner_loss)
      selected <- which.min(mean_loss)
      fit <- fit_path(x[tr, , drop = FALSE], target[tr], balanced_weights(obs[tr]),
        md = obs[tr], block = blocks[[model_name]])
      prediction_matrix <- cbind(predict_path(fit, x), zero = 0)
      predicted_g <- prediction_matrix[, selected]
      predicted_log2 <- offset + predicted_g
      stopifnot(all(is.finite(predicted_log2)), all(is.finite(2^predicted_log2)))
      predicted_zero <- cbind(predict_path(fit, x[1:2, , drop = FALSE] * 0), zero = 0)[, selected]
      stopifnot(max(abs(predicted_zero)) < 1e-12)
      prediction_counter <- prediction_counter + 1L
      tuning <- data.table(time_h = tt, fold_id = fold_id, fold_type = fold_type, model = model_name,
        lambda = c(lambda_grid, Inf), internal_loss = mean_loss,
        selected = seq_along(mean_loss) == selected)
      all_tuning[[prediction_counter]] <- tuning
      model_file <- file.path(root, 'models', paste0(tt, 'h_', fold_id, '_', model_name, '.qs'))
      qsave(list(fit = fit, selected = selected, lambda = c(lambda_grid, Inf)[selected],
        held_partners = held, fold_type = fold_type, training_partners = training_partners,
        training_observation_ids = obs$observation_id[tr], inference_inputs = model_name,
        single_effect_offset = if (grepl('_bridge$', model_name)) 'f(xA)+f(xB)' else 'measured log2 singles'), model_file)
      saved <- qread(model_file)
      replay <- cbind(predict_path(saved$fit, x[te, , drop = FALSE]), zero = 0)[, saved$selected]
      stopifnot(max(abs(replay - predicted_g[te])) < 1e-10)
      pred <- obs[, .(observation_id, condition_id, cell, anchor, partner,
        nominal_A_uM, nominal_B_uM, biological_replicate, profile_A, profile_B)]
      pred[, `:=`(time_h = tt, fold_id = fold_id, fold_type = fold_type, model = model_name,
        split = fifelse(partner %chin% held, 'held', 'training'),
        predicted_single_sum_log2 = as.numeric(offset), predicted_interaction_log2 = predicted_g,
        predicted_ratio_AB = as.numeric(2^predicted_log2))]
      # Predictions are flushed and hashed before observed combination values join.
      pred_path <- paste0('predictions_', tt, 'h_', fold_id, '_', model_name)
      write_table(pred, pred_path)
      all_audits[[prediction_counter]] <- data.table(time_h = tt, fold_id = fold_id, fold_type = fold_type, model = model_name,
        training_partners = paste(training_partners, collapse = ';'), held_partners = paste(held, collapse = ';'),
        training_observations = length(tr), held_observations = length(te), source_overlap = length(overlap),
        molecular_terms = 40L, frozen_single_offset_error = offset_error, symmetry_error = symmetry_error,
        zero_input_error = max(abs(predicted_zero)), replay_error = max(abs(replay - predicted_g[te])),
        prediction_sha256 = digest(file.path(root, 'data', paste0(pred_path, '.tsv.gz')), file = TRUE, algo = 'sha256'))
      pred[, `:=`(ratio_A = obs$ratio_A, ratio_B = obs$ratio_B, ratio_AB = obs$ratio_AB)]
      all_predictions[[prediction_counter]] <- pred
      coefficient <- if (selected <= length(lambda_grid)) as.matrix(coef(fit$fit, s = lambda_grid[selected]))[-1, 1] else
        setNames(rep(0, sum(fit$keep)), colnames(x)[fit$keep])
      all_coefficients[[prediction_counter]] <- data.table(time_h = tt, fold_id = fold_id, fold_type = fold_type, model = model_name,
        term = colnames(x)[fit$keep], scaled_coefficient = as.numeric(coefficient),
        raw_coefficient = as.numeric(coefficient) / fit$rms)
    }
    write_table(rbindlist(all_predictions), 'predictions')
    write_table(rbindlist(all_tuning), 'interaction_selection')
    write_table(rbindlist(all_audits), 'validation')
    write_table(rbindlist(all_coefficients), 'interaction_coefficients')
  }
}
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('All models complete; predictions and validation saved.')
