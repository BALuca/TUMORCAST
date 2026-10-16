#!/usr/bin/env Rscript
# Independent replay of the two-input combination model, with checks that joint
# learning changes shared molecular rules and preserves measured RNA means.
suppressPackageStartupMessages({library(data.table); library(qs)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/82_joint_program_rules_and_pair_interactions'
original <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original'
setwd(original)
source('modular_pair_pilot/bank_operator.R')
archive <- readRDS('modular_pair_pilot/inputs.rds')
PREP <- archive$prep
settings <- qread(file.path(root, 'models/settings.qs'))
native <- settings$native_doses
PREP$m <- PREP$m[PREP$m$drug %in% settings$selected_drugs, ]
single <- fread(file.path(root, 'data/held_single_predictions_frozen.tsv'))
pair <- fread(file.path(root, 'data/held_combination_predictions_frozen.tsv'))
oracle <- fread(file.path(root, 'data/held_measured_displacement_function_frozen.tsv'))
latent_pairs <- which(upper.tri(matrix(0, 4L, 4L), diag = TRUE), arr.ind = TRUE)
checks <- mean_checks <- list()

displacement <- function(cell, drug, dose, model) {
  c <- model$context
  baseline <- getmatrix(PREP$raw, cell, P)
  coordinates <- sweep(sweep(sweep(baseline, 2L, c$baseline_mean), 2L, c$baseline_sd, '/') %*%
    c$baseline_projection, 2L, c$baseline_projection_scale, '/')
  g <- (match(drug, model$selected_drugs) - 1L) * 3L + match(dose, native)
  sweep(sweep(coordinates, 2L, c$group_baseline_mean[g, ]) %*% model$molecular_coefficient[, , g],
    2L, model$group_response_mean[g, ], '+')
}
function_value <- function(response, model, assay = 'secondary') {
  par <- model$decoder
  value <- as.numeric(response %*% par$w) + as.numeric(tanh(response) %*% par$a[1:40]) +
    as.numeric((response^2 / (1 + response^2)) %*% par$a[40L + 1:40])
  latent <- as.numeric(tanh(response %*% par$u))
  value <- value + sum(latent[latent_pairs[, 1]] * latent[latent_pairs[, 2]] * par$a[80L + 1:10])
  if (assay != 'secondary') value <- value + as.numeric(response %*%
    par$deviation[, match(assay, c('primary', 'nci_FF', 'nci_FG'))])
  value
}
interaction_value <- function(a, b, model) {
  par <- model$decoder
  value <- sum(tanh(a) * tanh(b) * par$interaction[1:40])
  ha <- as.numeric(tanh(a %*% par$u))
  hb <- as.numeric(tanh(b %*% par$u))
  for (j in seq_len(nrow(latent_pairs))) {
    p <- latent_pairs[j, 1]
    q <- latent_pairs[j, 2]
    h <- if (p == q) ha[p] * hb[p] else (ha[p] * hb[q] + ha[q] * hb[p]) / sqrt(2)
    value <- value + h * par$interaction[40L + j]
  }
  value
}
single_point <- function(cell, drug, dose, model, assay = 'secondary') {
  values <- vapply(native, function(d) function_value(displacement(cell, drug, d, model), model, assay), numeric(1))
  approx(log(native), values, xout = log(dose), rule = 1)$y
}
pair_point <- function(row, model) {
  assay <- paste0('nci_', row$SCREENER)
  interaction <- matrix(0, 3L, 3L)
  for (i in 1:3) for (j in 1:3) interaction[i, j] <- interaction_value(
    displacement(row$depmap_id, row$drugA, native[i], model),
    displacement(row$depmap_id, row$drugB, native[j], model), model)
  at_a <- apply(interaction, 2L, function(v) approx(log(native), v, xout = log(row$doseA))$y)
  value <- approx(log(native), at_a, xout = log(row$doseB))$y
  c(total = single_point(row$depmap_id, row$drugA, row$doseA, model, assay) +
    single_point(row$depmap_id, row$drugB, row$doseB, model, assay) + value, departure = value)
}
paths <- list.files(file.path(root, 'models'), pattern = '__.*bridge.qs$', full.names = TRUE)
stopifnot(length(paths) == 26L)
for (path in paths) {
  model <- qread(path)
  job <- sub('__.*', '', basename(path))
  stopifnot(length(model$selected_drugs) == 10L, model$functional_baseline_inputs == FALSE,
    identical(model$functional_features, P), identical(dim(model$molecular_coefficient), c(8L, 40L, 30L)),
    !any(model$training_cells %in% model$held_cells), !any(model$training_molecular_cells %in% model$held_cells),
    !AU %in% c(model$training_cells, model$training_molecular_cells))
  rows <- single[family == model$family & depmap_id %in% model$held_cells]
  rows <- rows[, head(.SD, 2L), by = .(drug, dose)]
  single_error <- pair_error <- departure_error <- 0
  if (nrow(rows)) {
    prediction <- vapply(seq_len(nrow(rows)), function(i)
      single_point(rows$depmap_id[i], rows$drug[i], rows$dose[i], model), numeric(1))
    single_error <- max(abs(prediction - rows$predicted_log2))
  }
  combos <- if (grepl('^recipient_', job)) pair[family == model$family & regime == 'recipient' &
    depmap_id %in% model$held_cells] else if (!is.na(model$held_partner)) pair[
      family == model$family & regime == 'partner' & excluded_partner == model$held_partner] else pair[0]
  combos <- combos[, head(.SD, 1L), by = .(pair, doseA, doseB, SCREENER)]
  if (nrow(combos)) {
    predictions <- vapply(seq_len(nrow(combos)), function(i) pair_point(combos[i], model), numeric(2))
    pair_error <- max(abs(predictions[1L, ] - combos$predicted_log2))
    departure_error <- max(abs(predictions[2L, ] - combos$predicted_departure_log2))
  }
  measured <- as.data.table(PREP$m)[depmap_id %in% model$training_molecular_cells]
  group_check <- measured[, {
    dr <- .BY$drug
    dd <- .BY$dose
    values <- vapply(depmap_id, function(cell) as.numeric(displacement(cell, dr, dd, model)), numeric(40))
    physical <- sweep(t(values), 2L, model$context$response_sd, '*')
    .(cells = .N, largest_mean_difference = max(abs(colMeans(physical) - colMeans(as.matrix(.SD)))))
  }, by = .(drug, dose), .SDcols = P]
  group_check[, `:=`(job = job, family = model$family)]
  mean_checks[[path]] <- group_check
  a <- displacement(measured$depmap_id[1], measured$drug[1], measured$dose[1], model)
  b <- displacement(measured$depmap_id[2], measured$drug[2], measured$dose[2], model)
  symmetry_error <- abs(interaction_value(a, b, model) - interaction_value(b, a, model))
  zero_error <- abs(interaction_value(a, b * 0, model))
  stopifnot(single_error < 1e-8, pair_error < 1e-8, departure_error < 1e-8,
    max(group_check$largest_mean_difference) < 1e-8, symmetry_error < 1e-12, zero_error == 0)
  checks[[path]] <- data.table(job, family = model$family, single_conditions = nrow(rows),
    pair_conditions = nrow(combos), maximum_single_difference = single_error,
    maximum_pair_difference = pair_error, maximum_departure_difference = departure_error,
    symmetry_error, zero_error, solver_status = model$optimization$convergence,
    scaled_gradient = model$optimization$maximum_scaled_gradient)
}
correlation <- function(x, y) if (length(x) >= 4 && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
projection <- oracle[, .(cells = uniqueN(depmap_id),
  program_projection_agreement = correlation(measured_displacement_prediction, predicted_displacement_prediction),
  measured_displacement_to_viability = correlation(observed_log2, measured_displacement_prediction),
  predicted_displacement_to_viability = correlation(observed_log2, predicted_displacement_prediction)),
  by = .(family, drug, dose)]
fwrite(projection, file.path(root, 'data/viability_relevant_displacement_validation.tsv'), sep = '\t')
fwrite(rbindlist(checks), file.path(root, 'data/independent_prediction_replay.tsv'), sep = '\t')
fwrite(rbindlist(mean_checks), file.path(root, 'data/measured_RNA_mean_preservation.tsv'), sep = '\t')
print(rbindlist(checks)[, .(models = .N, maximum_single_difference = max(maximum_single_difference),
  maximum_pair_difference = max(maximum_pair_difference), unconverged = sum(solver_status != 0))])
print(projection[drug == 'Lapatinib ditosylate'])
cat('Independent replay, joint rule inputs, fixed RNA means, and pair interaction checks passed.\n')
