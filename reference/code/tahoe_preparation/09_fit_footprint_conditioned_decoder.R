#!/usr/bin/env Rscript

## Stage 2A: add one globally shared 9 x 9 displacement-by-footprint ridge
## correction to the exact nested v2 program decoder.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
source(file.path(WORK_ROOT, "src", "R", "lib_heatmap.R"))
RESULT <- stage_result("09_functional_modes")
ARTIFACT <- stage_artifact("09_functional_modes")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(stage_result("08_summary"),
    "stage1_reproduction_complete.flag")))
  fail("Stage 2A is forbidden until the Stage 1 reproduction gate passes.")

make_assignment <- function(lines, salt, folds = N_INNER_FOLDS) {
  set.seed(SEED + salt)
  shuffled <- sample(sort(unique(lines)), length(unique(lines)), FALSE)
  setNames(rep(seq_len(min(folds, length(shuffled))),
    length.out = length(shuffled)), shuffled)
}

make_footprint <- function(lines) {
  use <- score_meta$depmap_id %chin% lines
  means <- data.table(perturbation_id = score_meta$perturbation_id[use],
    D[use, , drop = FALSE])[, lapply(.SD, mean), by = perturbation_id,
      .SDcols = programs]
  matrix <- as.matrix(means[, ..programs])
  rownames(matrix) <- means$perturbation_id
  matrix
}

choose_nested_base <- function(lines, salt) {
  rows <- function_rows$depmap_id %chin% lines & chain
  assignment <- make_assignment(lines, 500L + salt)
  losses <- matrix(NA_real_, max(assignment), length(DECODER_LAMBDA_GRID))
  for (subfold in seq_len(max(assignment))) {
    validation_lines <- names(assignment)[assignment == subfold]
    learn <- rows & !function_rows$depmap_id %chin% validation_lines
    valid <- rows & function_rows$depmap_id %chin% validation_lines
    path <- ridge_path(D[function_position[learn], , drop = FALSE],
      function_rows$prism_sensitization[learn],
      D[function_position[valid], , drop = FALSE], DECODER_LAMBDA_GRID)
    losses[subfold, ] <- colMeans(sweep(path, 1L,
      function_rows$prism_sensitization[valid], "-") ^ 2)
  }
  index <- which.min(colMeans(losses))
  fit <- fit_ridge(D[function_position[rows], , drop = FALSE],
    function_rows$prism_sensitization[rows], DECODER_LAMBDA_GRID[index])
  list(raw = raw_ridge_coefficients(fit), lambda = DECODER_LAMBDA_GRID[index],
    mean_validation_mse = colMeans(losses)[index])
}

make_interactions <- function(function_indices, footprint) {
  positions <- function_position[function_indices]
  reference <- footprint[function_rows$perturbation_id[function_indices],
    programs, drop = FALSE]
  displacement <- D[positions, programs, drop = FALSE] - reference
  pieces <- lapply(seq_along(programs), function(index)
    reference * displacement[, index])
  design <- do.call(cbind, pieces)
  colnames(design) <- unlist(lapply(programs, function(context_program)
    paste0(context_program, "_by_", programs)), use.names = FALSE)
  if (any(!is.finite(design))) fail("Interaction design is nonfinite.")
  list(x = design, displacement = displacement, reference = reference)
}

scale_interactions <- function(x, weights, scale = NULL, keep = NULL) {
  if (is.null(scale)) {
    scale <- sqrt(colSums(x ^ 2 * weights) / sum(weights))
    keep <- is.finite(scale) & scale > 0
  }
  if (!any(keep)) fail("No variable interaction feature remains.")
  list(x = sweep(x[, keep, drop = FALSE], 2L, scale[keep], "/"),
    scale = scale, keep = keep)
}

interaction_path <- function(x_train, y_train, perturbations, x_test) {
  weights <- balanced_weights(perturbations)
  scaled_train <- scale_interactions(x_train, weights)
  scaled_test <- scale_interactions(x_test, rep(1, nrow(x_test)),
    scaled_train$scale, scaled_train$keep)
  root_weight <- sqrt(weights)
  decomposition <- svd(scaled_train$x * root_weight,
    nu = min(dim(scaled_train$x)), nv = min(dim(scaled_train$x)))
  uy <- as.numeric(crossprod(decomposition$u, y_train * root_weight))
  cross <- scaled_test$x %*% decomposition$v
  predictions <- vapply(INTERACTION_LAMBDA_GRID, function(lambda)
    as.numeric(cross %*% (uy * decomposition$d /
      (decomposition$d ^ 2 + lambda))), numeric(nrow(x_test)))
  list(predictions = predictions, scale = scaled_train$scale,
    keep = scaled_train$keep)
}

fit_interaction <- function(x, y, perturbations, lambda) {
  weights <- balanced_weights(perturbations)
  scaled <- scale_interactions(x, weights)
  root_weight <- sqrt(weights)
  decomposition <- svd(scaled$x * root_weight,
    nu = min(dim(scaled$x)), nv = min(dim(scaled$x)))
  beta <- decomposition$v %*% (as.numeric(crossprod(decomposition$u,
    y * root_weight)) * decomposition$d /
      (decomposition$d ^ 2 + lambda))
  list(lambda = lambda, scale = scaled$scale, keep = scaled$keep,
    feature_names = colnames(x), beta_standardized = as.numeric(beta))
}

predict_interaction <- function(fit, x) {
  if (is.null(fit)) return(rep(0, nrow(x)))
  scaled <- scale_interactions(x, rep(1, nrow(x)), fit$scale, fit$keep)
  as.numeric(scaled$x %*% fit$beta_standardized)
}

raw_interaction_beta <- function(fit, names) {
  answer <- setNames(rep(0, length(names)), names)
  if (!is.null(fit)) answer[fit$feature_names[fit$keep]] <-
    fit$beta_standardized / fit$scale[fit$keep]
  answer
}

evaluate_function <- function(table, prediction_column, cohort, model) {
  current <- table[evaluation_included == TRUE]
  if (cohort == "activity_supported") current <- current[
    activity_supported == TRUE]
  rbindlist(lapply(c("pooled", "within_perturbation",
    "strict_cell_line_by_perturbation"), function(metric) {
      residual <- function_residual_table(current, metric, prediction_column)
      weights <- balanced_weights(residual$drug)
      data.table(cohort = cohort, prediction_model = model, metric = metric,
        n_rows = nrow(residual), n_conditions = uniqueN(residual$perturbation_id),
        n_drugs = uniqueN(residual$drug), n_lines = uniqueN(residual$depmap_id),
        pearson = safe_cor(residual$observed_residual,
          residual$predicted_residual),
        drug_balanced_pearson = weighted_cor(residual$observed_residual,
          residual$predicted_residual, weights),
        mse = mean((residual$observed_residual -
          residual$predicted_residual) ^ 2),
        predicted_residual_sd = sd(residual$predicted_residual),
        observed_residual_sd = sd(residual$observed_residual),
        correlation_method = "pearson")
    }))
}

## Regenerated v2 inputs only.
stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- stability[outer_fold == OUTER_FOLD & stable == TRUE,
  unique(program)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS)
  fail("Stage 2A requires the frozen nine-program response bank.")
score_long <- fread(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"))[program %chin% programs]
score_meta <- dcast(score_long, observation_id + depmap_id + cell_line +
  perturbation_id + drug + dose + dose_panel + role ~ program,
  value.var = "score")
D <- as.matrix(score_meta[, ..programs])
rownames(D) <- score_meta$observation_id
function_rows <- fread(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"))
function_rows[, `:=`(observation_id = paste(depmap_id, perturbation_id,
  sep = "__"), drug = tahoe_drug, dose = tahoe_dose_uM)]
function_position <- match(function_rows$observation_id,
  score_meta$observation_id)
chain <- function_rows$chain_observed == 1L & !is.na(function_position)
folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- folds[outer_fold != OUTER_FOLD, depmap_id]
held_lines <- folds[outer_fold == OUTER_FOLD, depmap_id]
inner_assignment <- make_inner_assignment(training_lines)
activity <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
active_ids <- activity[activity_supported == TRUE, perturbation_id]

outer_coefficients <- fread(file.path(stage_result(
  "05_program_function_bridge"), "function_decoder_coefficients.tsv"))
outer_decoder <- list(intercept = outer_coefficients[
  term == "(Intercept)", coefficient], beta = setNames(outer_coefficients[
    term != "(Intercept)", coefficient], outer_coefficients[
      term != "(Intercept)", term]))
outer_footprint_dt <- fread(file.path(stage_result(
  "05_program_function_bridge"), "reference_perturbation_footprints.tsv.gz"))
outer_footprint <- as.matrix(outer_footprint_dt[, ..programs])
rownames(outer_footprint) <- outer_footprint_dt$perturbation_id

## Five training-only validation folds. Each fold rebuilds reference means and
## refits the exact nested base decoder before testing the correction grid.
selection_parts <- vector("list", N_INNER_FOLDS)
for (inner_fold in seq_len(N_INNER_FOLDS)) {
  validation_lines <- names(inner_assignment)[inner_assignment == inner_fold]
  learn_lines <- names(inner_assignment)[inner_assignment != inner_fold]
  train <- which(function_rows$depmap_id %chin% learn_lines & chain)
  valid <- which(function_rows$depmap_id %chin% validation_lines & chain)
  local_decoder <- choose_nested_base(learn_lines, inner_fold)
  local_footprint <- make_footprint(learn_lines)
  valid_before_footprint <- length(valid)
  valid <- valid[function_rows$perturbation_id[valid] %chin%
    rownames(local_footprint)]
  train_design <- make_interactions(train, local_footprint)
  valid_design <- make_interactions(valid, local_footprint)
  base_train <- local_decoder$raw$intercept + as.numeric(
    D[function_position[train], programs, drop = FALSE] %*%
      local_decoder$raw$beta[programs])
  base_valid <- local_decoder$raw$intercept + as.numeric(
    D[function_position[valid], programs, drop = FALSE] %*%
      local_decoder$raw$beta[programs])
  path <- interaction_path(train_design$x,
    function_rows$prism_sensitization[train] - base_train,
    function_rows$perturbation_id[train], valid_design$x)
  predictions <- cbind(base_valid, sweep(path$predictions, 1L,
    base_valid, "+"))
  candidates <- c(Inf, INTERACTION_LAMBDA_GRID)
  weights <- balanced_weights(function_rows$perturbation_id[valid])
  selection_parts[[inner_fold]] <- data.table(inner_fold = inner_fold,
    candidate = c("no_interaction", rep("interaction_ridge",
      length(INTERACTION_LAMBDA_GRID))), lambda = candidates,
    validation_mse = vapply(seq_along(candidates), function(index)
      sum(weights * (function_rows$prism_sensitization[valid] -
        predictions[, index]) ^ 2) / sum(weights), numeric(1L)),
    validation_pearson = vapply(seq_along(candidates), function(index)
      safe_cor(function_rows$prism_sensitization[valid],
        predictions[, index]), numeric(1L)),
    nested_base_lambda = local_decoder$lambda,
    nested_base_validation_mse = local_decoder$mean_validation_mse,
    n_training_lines = length(learn_lines),
    n_validation_lines = length(validation_lines),
    n_training_rows = length(train), n_validation_rows = length(valid),
    n_validation_rows_excluded_no_reference =
      valid_before_footprint - length(valid),
    n_features = ncol(train_design$x), condition_equal_weighting = TRUE,
    interaction_scaling = "training-only weighted RMS; no centering",
    correction_intercept = FALSE, perturbation_identity_feature = FALSE)
}
selection_folds <- rbindlist(selection_parts)
selection <- selection_folds[, .(
  mean_validation_mse = mean(validation_mse),
  mean_validation_pearson = mean(validation_pearson),
  n_inner_folds = .N), by = .(candidate, lambda)]
selected_index <- which.min(selection$mean_validation_mse)
selection[, `:=`(selected = seq_len(.N) == selected_index,
  selected_lambda = lambda[selected_index],
  no_interaction_selected = candidate[selected_index] == "no_interaction")]
fwrite(rbindlist(list(selection_folds[, record_type := "inner_fold"],
  selection[, `:=`(record_type = "aggregate", inner_fold = NA_integer_,
    validation_mse = NA_real_, validation_pearson = NA_real_,
    nested_base_lambda = NA_real_, nested_base_validation_mse = NA_real_,
    n_training_lines = NA_integer_, n_validation_lines = NA_integer_,
    n_training_rows = NA_integer_, n_validation_rows = NA_integer_,
    n_validation_rows_excluded_no_reference = NA_integer_,
    n_features = 81L, condition_equal_weighting = TRUE,
    interaction_scaling = "training-only weighted RMS; no centering",
    correction_intercept = FALSE, perturbation_identity_feature = FALSE)]),
  fill = TRUE), file.path(RESULT, "interaction_inner_selection.tsv"),
  sep = "\t")

## Fit the selected correction on all 26 training lines.
outer_train <- which(function_rows$depmap_id %chin% training_lines & chain)
outer_train_design <- make_interactions(outer_train, outer_footprint)
outer_base_train <- outer_decoder$intercept + as.numeric(
  D[function_position[outer_train], programs, drop = FALSE] %*%
    outer_decoder$beta[programs])
selected_lambda <- selection[selected == TRUE, lambda]
interaction_fit <- if (is.infinite(selected_lambda)) NULL else
  fit_interaction(outer_train_design$x,
    function_rows$prism_sensitization[outer_train] - outer_base_train,
    function_rows$perturbation_id[outer_train], selected_lambda)
saveRDS(interaction_fit, file.path(ARTIFACT,
  "footprint_conditioned_decoder.rds"))
raw_beta <- raw_interaction_beta(interaction_fit,
  colnames(outer_train_design$x))
coefficient_matrix <- matrix(raw_beta, nrow = length(programs), byrow = TRUE,
  dimnames = list(context_program = programs, reference_program = programs))
coefficient_long <- as.data.table(as.table(coefficient_matrix))
setnames(coefficient_long, c("context_program", "reference_program",
  "coefficient"))
coefficient_long[, `:=`(lambda = selected_lambda,
  no_interaction_selected = is.null(interaction_fit),
  no_intercept = TRUE, n_training_lines = length(training_lines),
  n_training_rows = length(outer_train))]
fwrite(coefficient_long, file.path(RESULT,
  "interaction_coefficient_matrix.tsv"), sep = "\t")

## Frozen seven-line test predictions.
held <- which(function_rows$outer_fold == OUTER_FOLD & chain)
held_design <- make_interactions(held, outer_footprint)
base_prediction <- outer_decoder$intercept + as.numeric(
  D[function_position[held], programs, drop = FALSE] %*%
    outer_decoder$beta[programs])
mean_prediction <- outer_decoder$intercept + as.numeric(
  outer_footprint[function_rows$perturbation_id[held], programs,
    drop = FALSE] %*% outer_decoder$beta[programs])
correction <- predict_interaction(interaction_fit, held_design$x)
held_predictions <- data.table(
  observation_id = function_rows$observation_id[held],
  depmap_id = function_rows$depmap_id[held],
  cell_line = function_rows$cell_line[held],
  perturbation_id = function_rows$perturbation_id[held],
  drug = function_rows$drug[held], dose = function_rows$dose[held],
  observed_function = function_rows$prism_sensitization[held],
  evaluation_included = function_rows$chain_observed[held] == 1L,
  activity_supported = function_rows$perturbation_id[held] %chin% active_ids,
  mean_prediction = mean_prediction, base_prediction = base_prediction,
  interaction_correction = correction,
  footprint_conditioned_prediction = base_prediction + correction,
  base_context_contribution = base_prediction - mean_prediction,
  footprint_conditioned_context_contribution =
    base_prediction + correction - mean_prediction)
fwrite(held_predictions, file.path(ARTIFACT,
  "held_functional_predictions.tsv.gz"), sep = "\t", compress = "gzip")

performance <- rbindlist(lapply(c("complete_matched_panel",
  "activity_supported"), function(cohort) rbindlist(list(
    evaluate_function(held_predictions, "base_prediction", cohort,
      "nested_base_decoder"),
    evaluate_function(held_predictions, "footprint_conditioned_prediction",
      cohort, "footprint_conditioned_decoder")))))
comparison <- merge(
  performance[prediction_model == "nested_base_decoder"],
  performance[prediction_model == "footprint_conditioned_decoder"],
  by = c("cohort", "metric"), suffixes = c("_base", "_interaction"))
comparison[, `:=`(
  pearson_change = pearson_interaction - pearson_base,
  drug_balanced_pearson_change = drug_balanced_pearson_interaction -
    drug_balanced_pearson_base,
  mse_change = mse_interaction - mse_base,
  predicted_sd_change = predicted_residual_sd_interaction -
    predicted_residual_sd_base)]
fwrite(performance, file.path(RESULT,
  "measured_functional_performance.tsv"), sep = "\t")
fwrite(comparison, file.path(RESULT,
  "measured_functional_comparison.tsv"), sep = "\t")

activity_comparison <- comparison[cohort == "activity_supported"]
improves_context <- activity_comparison[metric %chin% c(
  "within_perturbation", "strict_cell_line_by_perturbation"),
  any(pearson_change > 0)]
propagate_baseline <- !is.null(interaction_fit) && improves_context
gate <- data.table(
  training_selected_nonzero_correction = !is.null(interaction_fit),
  activity_within_pearson_change = activity_comparison[
    metric == "within_perturbation", pearson_change],
  activity_strict_pearson_change = activity_comparison[
    metric == "strict_cell_line_by_perturbation", pearson_change],
  held_within_or_strict_improved = improves_context,
  propagate_to_baseline = propagate_baseline,
  rule = "nonzero training selection and positive held within and/or strict Pearson change")
fwrite(gate, file.path(RESULT, "baseline_propagation_gate.tsv"), sep = "\t")

## Paired drug-block bootstrap of held activity-panel changes.
bootstrap_parts <- list()
activity_held <- held_predictions[evaluation_included == TRUE &
  activity_supported == TRUE]
for (metric_index in seq_along(c("within_perturbation",
    "strict_cell_line_by_perturbation"))) {
  metric <- c("within_perturbation",
    "strict_cell_line_by_perturbation")[metric_index]
  base_residual <- function_residual_table(activity_held, metric,
    "base_prediction")
  interaction_residual <- function_residual_table(activity_held, metric,
    "footprint_conditioned_prediction")
  stopifnot(identical(base_residual$observation_id,
    interaction_residual$observation_id))
  drugs <- sort(unique(base_residual$drug))
  counts <- as.numeric(table(factor(base_residual$drug, levels = drugs))[
    match(base_residual$drug, drugs)])
  set.seed(SEED + 9000L + metric_index)
  bootstrap_parts[[metric]] <- rbindlist(lapply(seq_len(N_BOOTSTRAPS),
    function(draw) {
      sampled <- sample(drugs, length(drugs), replace = TRUE)
      multiplicity <- as.numeric(table(factor(sampled, levels = drugs))[
        match(base_residual$drug, drugs)])
      balanced <- multiplicity / counts
      base_r <- weighted_cor(base_residual$observed_residual,
        base_residual$predicted_residual, multiplicity)
      interaction_r <- weighted_cor(interaction_residual$observed_residual,
        interaction_residual$predicted_residual, multiplicity)
      base_balanced <- weighted_cor(base_residual$observed_residual,
        base_residual$predicted_residual, balanced)
      interaction_balanced <- weighted_cor(
        interaction_residual$observed_residual,
        interaction_residual$predicted_residual, balanced)
      data.table(bootstrap = draw, metric = metric,
        base_pearson = base_r, interaction_pearson = interaction_r,
        pearson_change = interaction_r - base_r,
        base_drug_balanced_pearson = base_balanced,
        interaction_drug_balanced_pearson = interaction_balanced,
        drug_balanced_pearson_change = interaction_balanced - base_balanced)
    }))
}
bootstrap <- rbindlist(bootstrap_parts)
fwrite(bootstrap, file.path(ARTIFACT,
  "measured_functional_change_drug_bootstrap.tsv.gz"), sep = "\t",
  compress = "gzip")
bootstrap_intervals <- bootstrap[, .(
  pearson_change_ci_lower = quantile(pearson_change, 0.025, na.rm = TRUE),
  pearson_change_ci_upper = quantile(pearson_change, 0.975, na.rm = TRUE),
  probability_pearson_change_positive = mean(pearson_change > 0,
    na.rm = TRUE),
  drug_balanced_change_ci_lower = quantile(drug_balanced_pearson_change,
    0.025, na.rm = TRUE),
  drug_balanced_change_ci_upper = quantile(drug_balanced_pearson_change,
    0.975, na.rm = TRUE)), by = metric]
fwrite(bootstrap_intervals, file.path(RESULT,
  "measured_functional_change_bootstrap_intervals.tsv"), sep = "\t")

## Per-drug performance and prediction amplitude.
per_drug <- rbindlist(lapply(sort(unique(activity_held$drug)),
  function(drug_name) {
    current <- activity_held[drug == drug_name]
    metrics <- rbindlist(lapply(c("within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric) {
        base <- function_residual_table(current, metric, "base_prediction")
        conditioned <- function_residual_table(current, metric,
          "footprint_conditioned_prediction")
        data.table(metric = metric,
          base_pearson = safe_cor(base$observed_residual,
            base$predicted_residual),
          footprint_conditioned_pearson = safe_cor(
            conditioned$observed_residual, conditioned$predicted_residual))
      }))
    metrics[, `:=`(drug = drug_name, n_rows = nrow(current),
      n_conditions = uniqueN(current$perturbation_id),
      observed_context_sd = sd(current$observed_function -
        current$mean_prediction),
      base_context_sd = sd(current$base_context_contribution),
      footprint_conditioned_context_sd = sd(
        current$footprint_conditioned_context_contribution),
      correction_sd = sd(current$interaction_correction))]
    metrics[, pearson_change := footprint_conditioned_pearson - base_pearson]
    metrics
  }))
fwrite(per_drug, file.path(RESULT,
  "measured_functional_per_drug.tsv"), sep = "\t")

amplitude <- activity_held[, .(
  n_rows = .N, n_conditions = uniqueN(perturbation_id),
  observed_context_sd = sd(observed_function - mean_prediction),
  base_context_sd = sd(base_context_contribution),
  footprint_conditioned_context_sd = sd(
    footprint_conditioned_context_contribution),
  correction_sd = sd(interaction_correction),
  base_to_observed_sd_ratio = sd(base_context_contribution) /
    sd(observed_function - mean_prediction),
  footprint_conditioned_to_observed_sd_ratio = sd(
    footprint_conditioned_context_contribution) /
    sd(observed_function - mean_prediction))]
fwrite(amplitude, file.path(RESULT, "prediction_amplitude.tsv"), sep = "\t")

## SVD is interpretation only: these are readout modes over frozen programs.
decomposition <- svd(coefficient_matrix)
for (mode in seq_len(ncol(decomposition$u))) {
  anchor <- which.max(abs(decomposition$u[, mode]))
  if (decomposition$u[anchor, mode] < 0) {
    decomposition$u[, mode] <- -decomposition$u[, mode]
    decomposition$v[, mode] <- -decomposition$v[, mode]
  }
}
singular_values <- data.table(mode = paste0("M", seq_along(decomposition$d)),
  singular_value = decomposition$d,
  fraction = if (sum(decomposition$d) > 0) decomposition$d /
    sum(decomposition$d) else 0,
  cumulative_fraction = if (sum(decomposition$d) > 0)
    cumsum(decomposition$d) / sum(decomposition$d) else 0,
  interpretation_only = TRUE,
  terminology = "function-specific readout mode over frozen response programs")
weights <- rbindlist(lapply(seq_along(decomposition$d), function(mode)
  rbindlist(list(data.table(mode = paste0("M", mode), side = "context",
    program = programs, weight = decomposition$u[, mode]),
  data.table(mode = paste0("M", mode), side = "reference_footprint",
    program = programs, weight = decomposition$v[, mode])))))
fwrite(singular_values, file.path(RESULT,
  "functional_mode_singular_values.tsv"), sep = "\t")
fwrite(weights, file.path(RESULT, "functional_mode_weights.tsv"), sep = "\t")

mode_contributions <- rbindlist(lapply(seq_along(decomposition$d),
  function(mode) data.table(
    observation_id = held_predictions$observation_id,
    depmap_id = held_predictions$depmap_id,
    perturbation_id = held_predictions$perturbation_id,
    drug = held_predictions$drug, evaluation_included =
      held_predictions$evaluation_included,
    activity_supported = held_predictions$activity_supported,
    mode = paste0("M", mode),
    context_mode_score = as.numeric(held_design$displacement %*%
      decomposition$u[, mode]),
    reference_mode_score = as.numeric(held_design$reference %*%
      decomposition$v[, mode]),
    mode_contribution = decomposition$d[mode] *
      as.numeric(held_design$displacement %*% decomposition$u[, mode]) *
      as.numeric(held_design$reference %*% decomposition$v[, mode]))))
reconstruction <- mode_contributions[, .(
  reconstructed_correction = sum(mode_contribution)), by = observation_id]
if (max(abs(reconstruction$reconstructed_correction - correction[
    match(reconstruction$observation_id,
      held_predictions$observation_id)])) > 1e-8)
  fail("Functional modes do not reconstruct the interaction correction.")
fwrite(mode_contributions, file.path(ARTIFACT,
  "functional_mode_contributions_by_response.tsv.gz"), sep = "\t",
  compress = "gzip")
mode_by_drug <- mode_contributions[evaluation_included == TRUE &
  activity_supported == TRUE, .(
    n_responses = .N, signed_mean_contribution = mean(mode_contribution),
    mean_absolute_contribution = mean(abs(mode_contribution)),
    contribution_sd = sd(mode_contribution)), by = .(drug, mode)]
fwrite(mode_by_drug, file.path(RESULT,
  "functional_mode_contributions_by_drug.tsv"), sep = "\t")

## Drug-bootstrap stability of the fitted readout modes.
N_MODE_BOOTSTRAPS <- 200L
if (!is.null(interaction_fit)) {
  training_drugs <- sort(unique(function_rows$drug[outer_train]))
  mode_stability_parts <- vector("list", N_MODE_BOOTSTRAPS)
  set.seed(SEED + 12000L)
  for (draw in seq_len(N_MODE_BOOTSTRAPS)) {
    sampled_drugs <- sample(training_drugs, length(training_drugs), TRUE)
    sampled_indices <- unlist(lapply(seq_along(sampled_drugs), function(index)
      outer_train[function_rows$drug[outer_train] == sampled_drugs[index]]),
      use.names = FALSE)
    sampled_condition <- unlist(lapply(seq_along(sampled_drugs), function(index)
      paste(index, function_rows$perturbation_id[outer_train[
        function_rows$drug[outer_train] == sampled_drugs[index]]], sep = "__")),
      use.names = FALSE)
    fit <- fit_interaction(outer_train_design$x[
      match(sampled_indices, outer_train), , drop = FALSE],
      function_rows$prism_sensitization[sampled_indices] -
        outer_base_train[match(sampled_indices, outer_train)],
      sampled_condition, selected_lambda)
    draw_beta <- raw_interaction_beta(fit,
      colnames(outer_train_design$x))
    draw_matrix <- matrix(draw_beta, nrow = length(programs), byrow = TRUE)
    draw_svd <- svd(draw_matrix)
    remaining <- seq_along(draw_svd$d)
    rows <- vector("list", length(programs))
    for (mode in seq_along(decomposition$d)) {
      similarity <- vapply(remaining, function(candidate)
        mean(c(abs(sum(decomposition$u[, mode] * draw_svd$u[, candidate])),
          abs(sum(decomposition$v[, mode] * draw_svd$v[, candidate])))),
        numeric(1L))
      candidate <- remaining[which.max(similarity)]
      sign <- ifelse(sum(decomposition$u[, mode] *
        draw_svd$u[, candidate]) < 0, -1, 1)
      rows[[mode]] <- data.table(bootstrap = draw,
        mode = paste0("M", mode), matched_bootstrap_mode = candidate,
        singular_value = draw_svd$d[candidate],
        left_similarity = sum(decomposition$u[, mode] *
          (sign * draw_svd$u[, candidate])),
        right_similarity = sum(decomposition$v[, mode] *
          (sign * draw_svd$v[, candidate])),
        combined_absolute_similarity = similarity[which.max(similarity)])
      remaining <- setdiff(remaining, candidate)
    }
    mode_stability_parts[[draw]] <- rbindlist(rows)
  }
  mode_stability <- rbindlist(mode_stability_parts)
} else {
  mode_stability <- data.table(bootstrap = integer(), mode = character(),
    matched_bootstrap_mode = integer(), singular_value = numeric(),
    left_similarity = numeric(), right_similarity = numeric(),
    combined_absolute_similarity = numeric())
}
fwrite(mode_stability, file.path(ARTIFACT,
  "functional_mode_training_drug_bootstrap.tsv.gz"), sep = "\t",
  compress = "gzip")
mode_stability_summary <- if (nrow(mode_stability)) mode_stability[, .(
  median_singular_value = median(singular_value),
  singular_value_q10 = quantile(singular_value, 0.1),
  median_left_similarity = median(left_similarity),
  left_similarity_q10 = quantile(left_similarity, 0.1),
  median_right_similarity = median(right_similarity),
  right_similarity_q10 = quantile(right_similarity, 0.1),
  median_combined_absolute_similarity = median(
    combined_absolute_similarity)), by = mode] else
      data.table(mode = singular_values$mode,
        median_singular_value = NA_real_, singular_value_q10 = NA_real_,
        median_left_similarity = NA_real_, left_similarity_q10 = NA_real_,
        median_right_similarity = NA_real_, right_similarity_q10 = NA_real_,
        median_combined_absolute_similarity = NA_real_)
fwrite(mode_stability_summary, file.path(RESULT,
  "functional_mode_stability.tsv"), sep = "\t")

## One plot per PDF.
theme_v2 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1,
    plot.title = element_text(face = "bold"), legend.position = "bottom")
plot_rows <- held_predictions[evaluation_included == TRUE &
  activity_supported == TRUE]
scatter <- ggplot(plot_rows, aes(x = footprint_conditioned_prediction,
    y = observed_function)) +
  geom_point(alpha = 0.4, size = 0.8, color = "#276FBF") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.6, color = "#C44536") +
  annotate("text", x = Inf, y = -Inf,
    label = sprintf("Pearson r = %.3f", safe_cor(
      plot_rows$footprint_conditioned_prediction,
      plot_rows$observed_function)), hjust = 1.05, vjust = -0.6,
    size = 2.8) +
  labs(title = "Footprint-conditioned measured function",
    x = "Predicted PRISM sensitization", y = "Observed PRISM sensitization") +
  theme_v2()
ggsave(file.path(RESULT,
  "09a_footprint_conditioned_observed_vs_predicted.pdf"), scatter,
  width = 4.25, height = 4.15, device = cairo_pdf, bg = "white")

coefficient_heatmap <- draw_heatmap(coefficient_matrix,
  scale_rows = FALSE, show_row_names = TRUE, show_column_names = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  title = "Footprint-conditioned functional readout",
  column_title_gp = grid::gpar(fontsize = 10))
cairo_pdf(file.path(RESULT, "09b_interaction_coefficient_heatmap.pdf"),
  width = 4.8, height = 4.6)
draw(coefficient_heatmap)
dev.off()

leading_modes <- if (sum(singular_values$singular_value) > 0) {
  which(c(0, head(singular_values$cumulative_fraction, -1L)) < 0.9)
} else integer()
leading_modes <- head(leading_modes, 5L)
mode_plot_index <- data.table()
for (mode_index in leading_modes) {
  current <- weights[mode == paste0("M", mode_index)]
  current[, program := factor(program, levels = programs)]
  plot <- ggplot(current, aes(x = program, y = weight, fill = side)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(context = "#276FBF",
      reference_footprint = "#E08E45")) +
    labs(title = sprintf("Functional readout mode M%d", mode_index),
      x = "Frozen response program", y = "SVD weight", fill = "Mode side") +
    theme_v2()
  filename <- sprintf("09c_functional_mode_M%d_weights.pdf", mode_index)
  ggsave(file.path(RESULT, filename), plot, width = 4.6, height = 4.15,
    device = cairo_pdf, bg = "white")
  mode_plot_index <- rbind(mode_plot_index, data.table(
    mode = paste0("M", mode_index), figure = filename,
    singular_value = decomposition$d[mode_index]))
}

p_singular <- ggplot(singular_values,
  aes(x = factor(mode, levels = mode), y = fraction)) +
  geom_col(width = 0.68, fill = "#6A51A3") +
  geom_point(aes(y = cumulative_fraction), color = "#C44536", size = 1.7) +
  geom_line(aes(y = cumulative_fraction, group = 1), color = "#C44536",
    linewidth = 0.45) +
  labs(title = "Functional-mode singular-value spectrum",
    subtitle = "Bars: individual fraction; red points: cumulative fraction",
    x = "Interpretation-only functional mode", y = "Fraction") +
  theme_v2()
ggsave(file.path(RESULT, "09d_functional_mode_singular_values.pdf"),
  p_singular, width = 4.8, height = 4.15, device = cairo_pdf, bg = "white")

stability_plot <- melt(mode_stability_summary,
  id.vars = "mode", measure.vars = c("left_similarity_q10",
    "right_similarity_q10"), variable.name = "side", value.name = "q10")
stability_plot[, side := factor(side,
  levels = c("left_similarity_q10", "right_similarity_q10"),
  labels = c("Context side", "Reference-footprint side"))]
p_stability <- ggplot(stability_plot,
  aes(x = factor(mode, levels = singular_values$mode), y = q10,
    color = side, group = side)) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35,
    color = "grey50") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.8) +
  scale_color_manual(values = c("Context side" = "#276FBF",
    "Reference-footprint side" = "#E08E45")) +
  labs(title = "Drug-bootstrap stability of functional modes",
    subtitle = "Tenth-percentile matched loading similarity across 200 bootstraps",
    x = "Functional mode", y = "q10 loading similarity", color = NULL) +
  theme_v2()
ggsave(file.path(RESULT, "09e_functional_mode_stability.pdf"),
  p_stability, width = 5.2, height = 4.2, device = cairo_pdf, bg = "white")

leading_labels <- paste0("M", leading_modes)
p_contribution <- ggplot(mode_by_drug[mode %chin% leading_labels],
  aes(x = factor(mode, levels = leading_labels),
    y = mean_absolute_contribution, fill = mode)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.14, alpha = 0.24, size = 0.5) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Functional-mode contributions across held-out drugs",
    x = "Leading functional mode",
    y = "Mean absolute contribution per drug") + theme_v2()
ggsave(file.path(RESULT, "09f_functional_mode_contributions_by_drug.pdf"),
  p_contribution, width = 5.0, height = 4.2, device = cairo_pdf, bg = "white")

comparison_plot <- melt(activity_comparison,
  id.vars = "metric", measure.vars = c("pearson_base",
    "pearson_interaction"), variable.name = "model", value.name = "pearson")
comparison_plot[, `:=`(
  metric = factor(metric, levels = c("pooled", "within_perturbation",
    "strict_cell_line_by_perturbation"), labels = c("Pooled",
      "Within perturbation", "Strict")),
  model = factor(model, levels = c("pearson_base", "pearson_interaction"),
    labels = c("Nested base", "Footprint-conditioned")))]
p_comparison <- ggplot(comparison_plot,
  aes(x = metric, y = pearson, fill = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.3f", pearson)),
    position = position_dodge(width = 0.72), vjust = -0.25, size = 2.3) +
  scale_fill_manual(values = c("Nested base" = "#4D4D4D",
    "Footprint-conditioned" = "#276FBF")) +
  labs(title = "Held-out functional performance",
    x = NULL, y = "Pearson correlation", fill = NULL) + theme_v2()
ggsave(file.path(RESULT, "09g_measured_functional_comparison.pdf"),
  p_comparison, width = 5.5, height = 4.2, device = cairo_pdf, bg = "white")

amplitude_plot <- data.table(
  component = factor(c("Observed context", "Nested base context",
    "Footprint-conditioned context", "Interaction correction"),
    levels = c("Observed context", "Nested base context",
      "Footprint-conditioned context", "Interaction correction")),
  sd = c(amplitude$observed_context_sd, amplitude$base_context_sd,
    amplitude$footprint_conditioned_context_sd, amplitude$correction_sd))
p_amplitude <- ggplot(amplitude_plot, aes(component, sd, fill = component)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.3f", sd)), vjust = -0.25, size = 2.4) +
  scale_fill_manual(values = c("Observed context" = "#4D4D4D",
    "Nested base context" = "#72B7B2",
    "Footprint-conditioned context" = "#276FBF",
    "Interaction correction" = "#E08E45"), guide = "none") +
  labs(title = "Held-out functional prediction amplitude",
    x = NULL, y = "Standard deviation") + theme_v2() +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))
ggsave(file.path(RESULT, "09h_prediction_amplitude.pdf"), p_amplitude,
  width = 5.3, height = 4.2, device = cairo_pdf, bg = "white")

per_drug_plot <- per_drug[is.finite(base_pearson) &
  is.finite(footprint_conditioned_pearson)]
per_drug_plot[, metric_label := factor(metric,
  levels = c("within_perturbation", "strict_cell_line_by_perturbation"),
  labels = c("Within perturbation", "Strict"))]
p_per_drug <- ggplot(per_drug_plot,
  aes(base_pearson, footprint_conditioned_pearson, color = metric_label)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55) +
  facet_wrap(~metric_label, scales = "free") +
  scale_color_manual(values = c("Within perturbation" = "#276FBF",
    "Strict" = "#E08E45"), guide = "none") +
  labs(title = "Per-drug functional performance",
    x = "Nested-base Pearson r",
    y = "Footprint-conditioned Pearson r") + theme_v2()
ggsave(file.path(RESULT, "09i_per_drug_performance.pdf"), p_per_drug,
  width = 7.6, height = 4.1, device = cairo_pdf, bg = "white")

bootstrap_plot <- merge(bootstrap_intervals,
  activity_comparison[, .(metric, pearson_change)], by = "metric")
bootstrap_plot[, metric_label := factor(metric,
  levels = c("within_perturbation", "strict_cell_line_by_perturbation"),
  labels = c("Within perturbation", "Strict"))]
p_bootstrap <- ggplot(bootstrap_plot,
  aes(x = pearson_change, y = metric_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey50") +
  geom_segment(aes(x = pearson_change_ci_lower,
    xend = pearson_change_ci_upper, yend = metric_label), linewidth = 0.6,
    color = "#276FBF") +
  geom_point(size = 2, color = "#276FBF") +
  labs(title = "Drug-block uncertainty in functional change",
    x = "Footprint-conditioned minus nested-base Pearson r", y = NULL) +
  theme_v2()
ggsave(file.path(RESULT, "09j_bootstrap_change_intervals.pdf"), p_bootstrap,
  width = 5.4, height = 3.6, device = cairo_pdf, bg = "white")

fwrite(mode_plot_index, file.path(RESULT,
  "functional_mode_figure_index.tsv"), sep = "\t")
figure_index <- rbindlist(list(data.table(
  figure = c("09a_footprint_conditioned_observed_vs_predicted.pdf",
    "09b_interaction_coefficient_heatmap.pdf",
    "09d_functional_mode_singular_values.pdf",
    "09e_functional_mode_stability.pdf",
    "09f_functional_mode_contributions_by_drug.pdf",
    "09g_measured_functional_comparison.pdf",
    "09h_prediction_amplitude.pdf", "09i_per_drug_performance.pdf",
    "09j_bootstrap_change_intervals.pdf"),
  description = c(
    "Held activity-panel observed versus footprint-conditioned function",
    "Nine-by-nine footprint-conditioned interaction coefficients",
    "Singular-value spectrum for all function-specific readout modes",
    "Drug-bootstrap stability of context and reference mode weights",
    "Leading-mode contribution distributions across held-out drugs",
    "Nested-base versus footprint-conditioned held-out performance",
    "Observed and predicted functional-context amplitude",
    "Per-drug nested-base versus footprint-conditioned performance",
    "Drug-block bootstrap intervals for held-out correlation changes")),
  mode_plot_index[, .(figure, description = paste(
    "Context and reference weights for", mode))]))
pdf_page_count <- function(path) {
  content <- suppressWarnings(system2("strings", path, stdout = TRUE,
    stderr = TRUE))
  as.integer(sum(grepl("/Type /Page([^s]|$)", content)))
}
figure_index[, pages := vapply(file.path(RESULT, figure), pdf_page_count,
  integer(1L))]
figure_index[, passed_one_page := pages == 1L]
if (!all(figure_index$passed_one_page))
  fail("A Stage 09 diagnostic PDF does not contain exactly one page.")
fwrite(figure_index, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")

selected_validation <- selection[selected == TRUE]
zero_validation <- selection[candidate == "no_interaction"]
stage_summary <- data.table(
  training_selected_nonzero = !is.null(interaction_fit),
  selected_lambda = selected_lambda,
  no_interaction_validation_mse = zero_validation$mean_validation_mse,
  selected_validation_mse = selected_validation$mean_validation_mse,
  validation_mse_change = selected_validation$mean_validation_mse -
    zero_validation$mean_validation_mse,
  held_activity_within_base = activity_comparison[
    metric == "within_perturbation", pearson_base],
  held_activity_within_conditioned = activity_comparison[
    metric == "within_perturbation", pearson_interaction],
  held_activity_within_change = gate$activity_within_pearson_change,
  held_activity_strict_base = activity_comparison[
    metric == "strict_cell_line_by_perturbation", pearson_base],
  held_activity_strict_conditioned = activity_comparison[
    metric == "strict_cell_line_by_perturbation", pearson_interaction],
  held_activity_strict_change = gate$activity_strict_pearson_change,
  propagated_to_baseline = propagate_baseline)
fwrite(stage_summary, file.path(RESULT, "stage2a_summary.tsv"), sep = "\t")
writeLines(c("# Footprint-conditioned functional decoder",
  "",
  sprintf(paste0("Training-only validation selected a nonzero correction ",
    "(lambda %.6g; MSE change %.6f versus the exact no-interaction model)."),
    selected_lambda, stage_summary$validation_mse_change),
  sprintf(paste0("On the held activity panel, within-perturbation Pearson ",
    "changed by %.4f and strict Pearson changed by %.4f."),
    gate$activity_within_pearson_change,
    gate$activity_strict_pearson_change),
  "The prespecified propagation gate failed, so no baseline conditioner was refit.",
  "SVD outputs are interpretation-only function-specific readout modes over the frozen response programs; they are not new response programs or evidence of distinct mechanisms."),
  file.path(RESULT, "README.md"))

audit <- data.table(
  check = c("stage1_gate_present", "nine_by_nine_interactions",
    "no_interaction_candidate_present", "no_correction_intercept",
    "no_drug_identity_predictor", "no_per_drug_models",
    "training_only_reference_means", "training_only_interaction_scaling",
    "condition_equal_selection", "held_out_not_used_for_selection"),
  passed = c(file.exists(file.path(stage_result("08_summary"),
    "stage1_reproduction_complete.flag")), nrow(coefficient_long) == 81L,
    any(selection$candidate == "no_interaction"), TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE))
fwrite(audit, file.path(RESULT, "functional_modes_required_checks.tsv"),
  sep = "\t")
if (!all(audit$passed)) fail("A Stage 2A required check failed.")
write_artifact_manifest(ARTIFACT,
  "scripts/09_fit_footprint_conditioned_decoder.R")
writeLines("complete", file.path(RESULT, "stage2a_complete.flag"))
message(sprintf(paste0("Stage 2A complete: selected %s; activity within ",
  "change %.4f; strict change %.4f; baseline propagation %s."),
  ifelse(is.null(interaction_fit), "no interaction", sprintf("lambda %.6g",
    selected_lambda)), gate$activity_within_pearson_change,
  gate$activity_strict_pearson_change, gate$propagate_to_baseline))
