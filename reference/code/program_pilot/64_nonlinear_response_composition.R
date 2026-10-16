#!/usr/bin/env Rscript
# Smooth response-to-function and symmetric response-pair interactions. Kernel
# features preserve biological zero and use only available individual inputs.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(kernlab)
  library(digest)
})
root <- '/work/results/program_pilot/64_nonlinear_response_composition'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data'
linear_root <- '/work/results/program_pilot/62_separate_response_interaction_bridge'
lambda_grid <- 10^seq(3, -4, length.out = 24L)
width_factors <- c(.25, 1, 4)
rna_times <- c(3, 24)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(64011L)
options(warn = 1)
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
inputs <- c(file.path(prior, c('matched_molecular_links.tsv.gz', 'representation_scores.qs')),
  file.path(linear_root, 'data/single_predictions.tsv.gz'),
  file.path(linear_root, 'models', paste0(rna_times, 'h_Original40_single_decoder.qs')))
stopifnot(all(file.exists(inputs)))
links <- fread(inputs[1])
bank <- qread(inputs[2])
single_all <- fread(inputs[3])[representation == 'Original40']
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')),
  'input_manifest')
qsave(list(lambda_grid = lambda_grid, width_factors = width_factors, primary_time = 3,
  interaction_kernel = 'symmetric tensor product of origin-centered Gaussian kernels',
  objective = 'equal raw and within-anchor-dose residual losses, balanced by partner/anchor/experiment',
  bank_sha256 = bank$original_weight_sha256), file.path(root, 'models/settings.qs'))

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
fit_path <- function(x, y, weights, md = NULL) {
  # Eigenvalue magnitudes must be retained: unit-variance standardization would
  # destroy Gaussian-kernel regularization and exaggerate tiny eigen-directions.
  if (!is.null(md)) {
    x <- rbind(x, center_context(x, md))
    y <- c(y, center_context(y, md))
    weights <- c(weights * .5, weights * .5)
  }
  fit <- glmnet(x, y, alpha = 0, lambda = lambda_grid, weights = weights,
    standardize = FALSE, intercept = FALSE, thresh = 1e-8, maxit = 100000L)
  stopifnot(fit$jerr == 0L)
  fit
}
predict_path <- function(fit, x) as.matrix(predict(fit, newx = x, s = lambda_grid))
objective <- function(y, prediction, md) {
  w <- balanced_weights(md)
  .5 * colSums((prediction - y)^2 * w) +
    .5 * colSums((center_context(prediction, md) - center_context(y, md))^2 * w)
}
origin_kernel <- function(x, sigma) {
  ordinary <- as.matrix(kernelMatrix(rbfdot(sigma = sigma), x))
  to_zero <- exp(-sigma * rowSums(x^2))
  ordinary - outer(to_zero, rep(1, nrow(x))) - outer(rep(1, nrow(x)), to_zero) + 1
}
pair_kernel <- function(k, a, b) {
  .5 * (k[a, a, drop = FALSE] * k[b, b, drop = FALSE] +
    k[a, b, drop = FALSE] * k[b, a, drop = FALSE])
}
normalize_kernel <- function(k) k / mean(diag(k))
kernel_basis <- function(k) {
  k <- (k + t(k)) / 2
  eig <- eigen(k, symmetric = TRUE)
  stopifnot(min(eig$values) > -1e-7 * max(eig$values))
  keep <- eig$values > max(eig$values) * 1e-9
  values <- eig$values[keep]
  vectors <- eig$vectors[, keep, drop = FALSE]
  x <- sweep(vectors, 2, sqrt(values), '*')
  colnames(x) <- sprintf('K%03d', seq_len(ncol(x)))
  rownames(x) <- rownames(k)
  list(features = x, vectors = vectors, values = values,
    discarded_numerical_rank = sum(!keep), minimum_eigenvalue = min(eig$values))
}

predictions <- selections <- audits <- single_outputs <- single_selections <- list()
counter <- 0L
for (tt in rna_times) {
  status(paste(tt, 'h | prepare outcome-free kernel bases'))
  obs <- copy(links[time_h == tt & molecular_supported == TRUE])
  setorder(obs, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)
  stopifnot(nrow(obs) == 496L, !anyDuplicated(obs$observation_id), all(obs$ratio_AB > 0))
  single <- copy(single_all[time_h == tt])
  profiles <- bank$metadata[time_h == tt, profile_id]
  response <- bank$original40[profiles, , drop = FALSE]
  response_rms <- sqrt(colMeans(response^2))
  stopifnot(all(response_rms > 0), all(single$profile_id %in% profiles))
  x <- sweep(response, 2, response_rms, '/')
  sigma_base <- 1 / median(as.numeric(dist(x))^2)
  single_index <- match(single$profile_id, profiles)
  a_index <- match(obs$profile_A, profiles)
  b_index <- match(obs$profile_B, profiles)
  stopifnot(!anyNA(a_index), !anyNA(b_index))

  controls_a <- cbind(effect = log2(obs$ratio_A), exposure = log1p(obs$nominal_A_uM),
    effect_squared = log2(obs$ratio_A)^2)
  controls_b <- cbind(effect = log2(obs$ratio_B), exposure = log1p(obs$nominal_B_uM),
    effect_squared = log2(obs$ratio_B)^2)
  control <- rbind(controls_a, controls_b)
  control_rms <- sqrt(colMeans(control^2))
  control <- sweep(control, 2, control_rms, '/')
  distinct_control <- unique(control)
  control_sigma_base <- 1 / median(as.numeric(dist(distinct_control))^2)
  n <- nrow(obs)
  single_bases <- rna_bases <- functional_bases <- combined_bases <- list()
  kernel_audits <- list()
  for (j in seq_along(width_factors)) {
    k <- origin_kernel(x, sigma_base * width_factors[j])
    dimnames(k) <- list(profiles, profiles)
    single_bases[[j]] <- kernel_basis(normalize_kernel(k))
    kr <- pair_kernel(k, a_index, b_index)
    symmetry_error <- max(abs(kr - pair_kernel(k, b_index, a_index)))
    k_zero <- origin_kernel(rbind(rep(0, ncol(x)), x), sigma_base * width_factors[j])
    zero_error <- max(abs(k_zero[1, ]))
    stopifnot(symmetry_error < 1e-10, zero_error < 1e-10)
    kr <- normalize_kernel(kr)
    kc <- origin_kernel(control, control_sigma_base * width_factors[j])
    kf <- normalize_kernel(pair_kernel(kc, seq_len(n), n + seq_len(n)))
    dimnames(kr) <- dimnames(kf) <- list(obs$observation_id, obs$observation_id)
    rna_bases[[j]] <- kernel_basis(kr)
    functional_bases[[j]] <- kernel_basis(kf)
    combined_bases[[j]] <- kernel_basis((kr + kf) / 2)
    kernel_audits[[j]] <- data.table(time_h = tt, width_factor = width_factors[j],
      rna_sigma = sigma_base * width_factors[j], functional_sigma = control_sigma_base * width_factors[j],
      symmetry_error = symmetry_error, origin_kernel_error = zero_error,
      single_rank = ncol(single_bases[[j]]$features), rna_pair_rank = ncol(rna_bases[[j]]$features),
      functional_pair_rank = ncol(functional_bases[[j]]$features), combined_rank = ncol(combined_bases[[j]]$features))
  }
  write_table(rbindlist(kernel_audits), paste0('kernel_validation_', tt, 'h'))
  qsave(list(profiles = profiles, response_rms = response_rms, control_rms = control_rms,
    sigma_base = sigma_base, control_sigma_base = control_sigma_base, width_factors = width_factors,
    single = single_bases, RNA = rna_bases, functional = functional_bases, combined = combined_bases,
    observation_ids = obs$observation_id), file.path(root, 'models', paste0(tt, 'h_kernel_bases.qs')))

  status(paste(tt, 'h | nonlinear single decoder | same grouped drug folds as62'))
  y_single <- log2(single$ratio)
  loss <- matrix(NA_real_, length(width_factors), length(lambda_grid))
  oof_single <- list()
  for (j in seq_along(width_factors)) {
    sx <- single_bases[[j]]$features[single_index, , drop = FALSE]
    oof <- matrix(NA_real_, nrow(single), length(lambda_grid))
    for (ff in sort(unique(single$fold))) {
      tr <- which(single$fold != ff)
      te <- which(single$fold == ff)
      fit <- fit_path(sx[tr, , drop = FALSE], y_single[tr], single$weight[tr])
      oof[te, ] <- predict_path(fit, sx[te, , drop = FALSE])
    }
    loss[j, ] <- colSums((oof - y_single)^2 * single$weight)
    oof_single[[j]] <- oof
  }
  best <- which(loss == min(loss), arr.ind = TRUE)[1, ]
  sx <- single_bases[[best[1]]]$features[single_index, , drop = FALSE]
  f_model <- fit_path(sx, y_single, single$weight)
  nonlinear_f <- setNames(predict_path(f_model, single_bases[[best[1]]]$features)[, best[2]], profiles)
  single[, `:=`(predicted_nonlinear_fitted_log2 = nonlinear_f[profile_id],
    predicted_nonlinear_excluded_drug_log2 = oof_single[[best[1]]][, best[2]])]
  single_outputs[[as.character(tt)]] <- single
  write_table(rbindlist(single_outputs), 'single_predictions')
  single_grid <- data.table(band = rep(seq_along(width_factors), length(lambda_grid)),
    lambda_index = rep(seq_along(lambda_grid), each = length(width_factors)), internal_loss = as.vector(loss))
  single_grid[, `:=`(width_factor = width_factors[as.integer(band)], lambda = lambda_grid[as.integer(lambda_index)], time_h = tt)]
  single_selections[[as.character(tt)]] <- single_grid
  write_table(rbindlist(single_selections), 'single_selection')
  qsave(list(fit = f_model, width_index = best[1], lambda_index = best[2],
    width_factor = width_factors[best[1]], lambda = lambda_grid[best[2]],
    training_individual_ids = single$sample_id), file.path(root, 'models', paste0(tt, 'h_single_decoder.qs')))

  frozen <- qread(file.path(linear_root, 'models', paste0(tt, 'h_Original40_single_decoder.qs')))
  frozen_x <- sweep(response[, frozen$fit$keep, drop = FALSE], 2, frozen$fit$rms, '/')
  linear_f <- setNames(as.numeric(predict(frozen$fit$fit, newx = frozen_x, s = frozen$lambda)), profiles)
  measured_offset <- log2(obs$ratio_A) + log2(obs$ratio_B)
  offsets <- list(Nonlinear_control = measured_offset, Nonlinear_conditioned = measured_offset,
    Nonlinear_g = linear_f[obs$profile_A] + linear_f[obs$profile_B],
    Nonlinear_bridge = nonlinear_f[obs$profile_A] + nonlinear_f[obs$profile_B])
  bases <- list(Nonlinear_control = functional_bases, Nonlinear_conditioned = combined_bases,
    Nonlinear_g = rna_bases, Nonlinear_bridge = rna_bases)
  partners <- sort(unique(obs$partner))
  held_sets <- c(lapply(partners, function(x) x), combn(partners, 2, simplify = FALSE))
  y <- log2(obs$ratio_AB)
  for (held in held_sets) {
    fold_id <- paste(held, collapse = '__')
    fold_type <- if (length(held) == 1) 'single_partner' else 'two_partners'
    tr <- which(!obs$partner %chin% held)
    te <- which(obs$partner %chin% held)
    training_partners <- setdiff(partners, held)
    source_ids <- function(x) unique(unlist(strsplit(x, ';', fixed = TRUE)))
    overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]), source_ids(obs$source_sample_ids_AB[te]))
    stopifnot(length(overlap) == 0)
    for (model in names(offsets)) {
      status(paste(tt, 'h | exclude', paste(held, collapse = '+'), '|', model))
      target <- y - offsets[[model]]
      losses <- array(NA_real_, c(length(training_partners), length(lambda_grid), length(width_factors)))
      for (j in seq_along(width_factors)) {
        xx <- bases[[model]][[j]]$features
        stopifnot(identical(rownames(xx), obs$observation_id))
        for (ii in seq_along(training_partners)) {
          fit_idx <- tr[obs$partner[tr] != training_partners[ii]]
          val_idx <- tr[obs$partner[tr] == training_partners[ii]]
          fit <- fit_path(xx[fit_idx, , drop = FALSE], target[fit_idx], balanced_weights(obs[fit_idx]), obs[fit_idx])
          candidate <- predict_path(fit, xx[val_idx, , drop = FALSE])
          losses[ii, , j] <- objective(target[val_idx], candidate, obs[val_idx])
        }
      }
      mean_loss <- apply(losses, c(2, 3), mean)
      choice <- which(mean_loss == min(mean_loss), arr.ind = TRUE)[1, ]
      zero_loss <- mean(vapply(training_partners, function(pp) {
        idx <- tr[obs$partner[tr] == pp]
        objective(target[idx], matrix(0, nrow = length(idx), ncol = 1), obs[idx])
      }, numeric(1)))
      zero <- zero_loss <= min(mean_loss)
      xx <- bases[[model]][[choice[2]]]$features
      fit <- fit_path(xx[tr, , drop = FALSE], target[tr], balanced_weights(obs[tr]), obs[tr])
      predicted_g <- if (zero) rep(0, nrow(obs)) else predict_path(fit, xx)[, choice[1]]
      ratio <- 2^(offsets[[model]] + predicted_g)
      stopifnot(all(is.finite(ratio)), all(ratio > 0))
      file <- file.path(root, 'models', paste0(tt, 'h_', fold_id, '_', model, '.qs'))
      qsave(list(fit = fit, lambda_index = choice[1], width_index = choice[2],
        lambda = if (zero) Inf else lambda_grid[choice[1]], width_factor = width_factors[choice[2]],
        zero_interaction = zero, model = model, held_partners = held, training_partners = training_partners,
        training_observation_ids = obs$observation_id[tr], basis_file = paste0(tt, 'h_kernel_bases.qs')), file)
      saved <- qread(file)
      replay <- if (saved$zero_interaction) rep(0, length(te)) else
        predict_path(saved$fit, xx[te, , drop = FALSE])[, saved$lambda_index]
      zero_prediction <- predict_path(fit, xx[1:2, , drop = FALSE] * 0)[, choice[1]]
      stopifnot(max(abs(replay - predicted_g[te])) < 1e-10, max(abs(zero_prediction)) < 1e-10)
      out <- obs[, .(observation_id, condition_id, cell, anchor, partner, nominal_A_uM,
        nominal_B_uM, biological_replicate, profile_A, profile_B)]
      out[, `:=`(time_h = tt, fold_id = fold_id, fold_type = fold_type, model = model,
        split = fifelse(partner %chin% held, 'held', 'training'),
        predicted_single_sum_log2 = as.numeric(offsets[[model]]), predicted_interaction_log2 = predicted_g,
        predicted_ratio_AB = as.numeric(ratio))]
      output_name <- paste0('prediction_', tt, 'h_', fold_id, '_', model)
      write_table(out, output_name)
      counter <- counter + 1L
      audits[[counter]] <- data.table(time_h = tt, fold_id = fold_id, fold_type = fold_type, model = model,
        source_overlap = length(overlap), replay_error = max(abs(replay - predicted_g[te])),
        zero_input_error = max(abs(zero_prediction)), basis_dimension = ncol(xx),
        training_observations = length(tr), held_observations = length(te),
        prediction_sha256 = digest(file.path(root, 'data', paste0(output_name, '.tsv.gz')), file = TRUE, algo = 'sha256'))
      selections[[counter]] <- data.table(time_h = tt, fold_id = fold_id, fold_type = fold_type,
        model = model, width_factor = width_factors[choice[2]], lambda = if (zero) Inf else lambda_grid[choice[1]],
        zero_interaction = zero, selected_internal_loss = if (zero) zero_loss else min(mean_loss))
      grid <- data.table(lambda_index = rep(seq_along(lambda_grid), length(width_factors)),
        width_index = rep(seq_along(width_factors), each = length(lambda_grid)), internal_loss = as.vector(mean_loss))
      grid[, `:=`(lambda = lambda_grid[as.integer(lambda_index)], width_factor = width_factors[as.integer(width_index)])]
      write_table(grid, paste0('tuning_', tt, 'h_', fold_id, '_', model))
      out[, `:=`(ratio_A = obs$ratio_A, ratio_B = obs$ratio_B, ratio_AB = obs$ratio_AB)]
      predictions[[counter]] <- out
    }
    write_table(rbindlist(predictions), 'predictions')
    write_table(rbindlist(audits), 'validation')
    write_table(rbindlist(selections), 'selection')
  }
}
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('All nonlinear models complete; predictions and checks saved.')
