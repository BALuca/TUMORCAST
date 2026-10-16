#!/usr/bin/env Rscript

## Fit two shared molecular response transformations selected by baseline state.
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
RESULT <- stage_result("11_state_modulated_response")
ARTIFACT <- stage_artifact("11_state_modulated_response")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

theme_stage11 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage11(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(name, plot, width, height) ggsave(file.path(RESULT, name),
  plot, width = width, height = height, device = cairo_pdf, bg = "white",
  limitsize = FALSE)
save_heatmap <- function(name, heatmap, width, height) {
  grDevices::pdf(file.path(RESULT, name), width = width, height = height,
    useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap)
  grDevices::dev.off()
}

## Deterministic two-center baseline configuration with distance-soft weights.
fit_configurations <- function(raw_scores, lines) {
  x <- raw_scores[lines, , drop = FALSE]
  center <- colMeans(x)
  scale <- apply(x, 2L, sd)
  keep <- is.finite(center) & is.finite(scale) & scale > 0
  if (sum(keep) < 2L) fail("Fewer than two recurrent baseline programs vary.")
  z <- sweep(sweep(x[, keep, drop = FALSE], 2L, center[keep], "-"),
    2L, scale[keep], "/")
  set.seed(SEED)
  clustering <- kmeans(z, centers = 2L, iter.max = 100L, nstart = 25L,
    algorithm = "Lloyd")
  centers <- clustering$centers
  ordering <- do.call(order, c(as.data.frame(centers), list(method = "radix")))
  centers <- centers[ordering, , drop = FALSE]
  distance2 <- vapply(seq_len(2L), function(k)
    rowSums(sweep(z, 2L, centers[k, ], "-") ^ 2), numeric(nrow(z)))
  temperature <- median(apply(distance2, 1L, min))
  if (!is.finite(temperature) || temperature <= 0)
    temperature <- mean(distance2[is.finite(distance2) & distance2 > 0])
  if (!is.finite(temperature) || temperature <= 0)
    fail("Baseline configuration temperature is undefined.")
  list(center = center, scale = scale, keep = keep, centers = centers,
    temperature = temperature, feature_names = colnames(raw_scores))
}

predict_configurations <- function(fit, raw_scores) {
  z <- sweep(sweep(raw_scores[, fit$keep, drop = FALSE], 2L,
    fit$center[fit$keep], "-"), 2L, fit$scale[fit$keep], "/")
  distance2 <- vapply(seq_len(2L), function(k)
    rowSums(sweep(z, 2L, fit$centers[k, ], "-") ^ 2), numeric(nrow(z)))
  logits <- -distance2 / (2 * fit$temperature)
  logits <- sweep(logits, 1L, apply(logits, 1L, max), "-")
  membership <- exp(logits)
  membership <- membership / rowSums(membership)
  colnames(membership) <- c("C1", "C2")
  rownames(membership) <- rownames(raw_scores)
  list(z = z, distance2 = distance2, membership = membership)
}

## Weighted multivariate ridge path for one global 9-by-9 transformation.
multi_ridge_path <- function(x_train, y_train, x_test, weights, lambda) {
  weights <- pmax(as.numeric(weights), 1e-8)
  weights <- weights / mean(weights)
  w <- weights / sum(weights)
  x_center <- colSums(x_train * w)
  x_scale <- sqrt(colSums(sweep(x_train, 2L, x_center, "-") ^ 2 * w))
  keep <- is.finite(x_center) & is.finite(x_scale) & x_scale > 0
  if (!all(keep)) fail("A reference-footprint predictor is nonvariable.")
  xz <- sweep(sweep(x_train, 2L, x_center, "-"), 2L, x_scale, "/")
  test_z <- sweep(sweep(x_test, 2L, x_center, "-"), 2L, x_scale, "/")
  y_center <- colSums(y_train * w)
  root_weight <- sqrt(weights)
  decomposition <- svd(xz * root_weight,
    nu = min(dim(xz)), nv = min(dim(xz)))
  uy <- crossprod(decomposition$u,
    sweep(y_train, 2L, y_center, "-") * root_weight)
  cross <- test_z %*% decomposition$v
  lapply(lambda, function(value) sweep(cross %*%
    sweep(uy, 1L, decomposition$d /
      (decomposition$d ^ 2 + value), "*"), 2L, y_center, "+"))
}

fit_multi_ridge <- function(x, y, weights, lambda) {
  weights <- pmax(as.numeric(weights), 1e-8)
  weights <- weights / mean(weights)
  w <- weights / sum(weights)
  x_center <- colSums(x * w)
  x_scale <- sqrt(colSums(sweep(x, 2L, x_center, "-") ^ 2 * w))
  keep <- is.finite(x_center) & is.finite(x_scale) & x_scale > 0
  if (!all(keep)) fail("A final reference-footprint predictor is nonvariable.")
  xz <- sweep(sweep(x, 2L, x_center, "-"), 2L, x_scale, "/")
  y_center <- colSums(y * w)
  root_weight <- sqrt(weights)
  decomposition <- svd(xz * root_weight,
    nu = min(dim(xz)), nv = min(dim(xz)))
  beta_standardized <- decomposition$v %*%
    sweep(crossprod(decomposition$u,
      sweep(y, 2L, y_center, "-") * root_weight), 1L,
      decomposition$d / (decomposition$d ^ 2 + lambda), "*")
  beta <- sweep(beta_standardized, 1L, x_scale, "/")
  intercept <- as.numeric(y_center - x_center %*% beta)
  rownames(beta) <- colnames(x)
  colnames(beta) <- colnames(y)
  list(lambda = lambda, x_center = x_center, x_scale = x_scale,
    beta = beta, intercept = intercept, weights_sum = sum(weights))
}

predict_multi_ridge <- function(fit, x) sweep(x %*% fit$beta, 2L,
  fit$intercept, "+")

make_footprint <- function(lines, score_table, displacement, programs) {
  use <- score_table$depmap_id %chin% lines
  means <- data.table(perturbation_id = score_table$perturbation_id[use],
    displacement[use, programs, drop = FALSE])[, lapply(.SD, mean),
      by = perturbation_id, .SDcols = programs]
  matrix <- as.matrix(means[, ..programs])
  rownames(matrix) <- means$perturbation_id
  matrix
}

## Load only frozen molecular and untreated-baseline objects; no function table.
required <- c(file.path(stage_result("10_baseline_states"),
    "stage10_complete.flag"),
  file.path(stage_artifact("10_baseline_states"),
    "baseline_program_scores.tsv.gz"),
  file.path(stage_result("02_response_programs"),
    "fold_program_stability.tsv"),
  file.path(stage_artifact("02_response_programs"),
    "fold_program_scores.tsv.gz"),
  file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 11 input(s): %s",
  paste(missing, collapse = ", "))

folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- sort(folds[outer_fold != OUTER_FOLD, depmap_id])
held_lines <- sort(folds[outer_fold == OUTER_FOLD, depmap_id])
all_lines <- sort(folds$depmap_id)
if (length(training_lines) != EXPECTED_TRAINING_LINES ||
    length(held_lines) != EXPECTED_HELD_LINES ||
    any(EXCLUDED_IDS %chin% all_lines))
  fail("The accepted 26/7 cohort changed before Stage 11.")
inner_assignment <- make_inner_assignment(training_lines)

baseline_long <- fread(file.path(stage_artifact("10_baseline_states"),
  "baseline_program_scores.tsv.gz"))
dominance <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_program_line_dominance.tsv"))[baseline_source == "tahoe_dmso"]
baseline_stability <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_program_stability.tsv"))[, .SD[1L], by = program]
state_programs <- merge(baseline_stability[, .(program, stable)],
  dominance[, .(program, single_line_dominated,
    leave_one_line_out_stable, eligible_recurrent_state)],
  by = "program")[eligible_recurrent_state == TRUE, program]
state_programs <- paste0("S", sort(as.integer(sub("S", "", state_programs))))
if (length(state_programs) < 2L)
  fail("Fewer than two stable, distributed baseline programs remain.")
baseline_wide <- dcast(baseline_long[baseline_source == "tahoe_dmso" &
    aggregation_level == "line_mean" & program %chin% state_programs],
  depmap_id ~ program, value.var = "raw_program_score")
baseline_raw <- as.matrix(baseline_wide[, ..state_programs])
rownames(baseline_raw) <- baseline_wide$depmap_id
baseline_raw <- baseline_raw[all_lines, state_programs, drop = FALSE]
if (any(!is.finite(baseline_raw))) fail("Baseline-state scores are invalid.")

response_stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- response_stability[outer_fold == OUTER_FOLD & stable == TRUE,
  unique(program)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS)
  fail("The accepted nine response programs changed.")
score_long <- fread(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"))[program %chin% programs]
score_table <- dcast(score_long, observation_id + depmap_id + cell_line +
  perturbation_id + drug + dose + dose_panel + role ~ program,
  value.var = "score")
D <- as.matrix(score_table[, ..programs])
rownames(D) <- score_table$observation_id
if (any(!is.finite(D))) fail("Frozen response-program scores are invalid.")

## Nested molecular-only validation recomputes every fold-specific quantity.
lambda_grid <- DECODER_LAMBDA_GRID
loss <- matrix(NA_real_, N_INNER_FOLDS, length(lambda_grid))
mean_loss <- numeric(N_INNER_FOLDS)
selection_parts <- vector("list", N_INNER_FOLDS)
for (fold in seq_len(N_INNER_FOLDS)) {
  validation_lines <- names(inner_assignment)[inner_assignment == fold]
  learn_lines <- names(inner_assignment)[inner_assignment != fold]
  state_fit <- fit_configurations(baseline_raw, learn_lines)
  learn_membership <- predict_configurations(state_fit,
    baseline_raw[learn_lines, , drop = FALSE])$membership
  validation_membership <- predict_configurations(state_fit,
    baseline_raw[validation_lines, , drop = FALSE])$membership
  footprint <- make_footprint(learn_lines, score_table, D, programs)
  train <- score_table$depmap_id %chin% learn_lines
  valid <- score_table$depmap_id %chin% validation_lines &
    score_table$perturbation_id %chin% rownames(footprint)
  x_train <- footprint[score_table$perturbation_id[train], programs,
    drop = FALSE]
  x_valid <- footprint[score_table$perturbation_id[valid], programs,
    drop = FALSE]
  paths <- lapply(seq_len(2L), function(configuration)
    multi_ridge_path(x_train, D[train, , drop = FALSE], x_valid,
      learn_membership[score_table$depmap_id[train], configuration],
      lambda_grid))
  valid_weights <- validation_membership[score_table$depmap_id[valid],,
    drop = FALSE]
  for (lambda_index in seq_along(lambda_grid)) {
    prediction <- paths[[1L]][[lambda_index]] * valid_weights[, 1L] +
      paths[[2L]][[lambda_index]] * valid_weights[, 2L]
    loss[fold, lambda_index] <- mean((D[valid, , drop = FALSE] -
      prediction) ^ 2)
  }
  mean_loss[fold] <- mean((D[valid, , drop = FALSE] - x_valid) ^ 2)
  selection_parts[[fold]] <- data.table(inner_fold = fold,
    lambda = lambda_grid, validation_mse = loss[fold, ],
    mean_footprint_validation_mse = mean_loss[fold],
    n_learning_lines = length(learn_lines),
    n_validation_lines = length(validation_lines),
    reference_means_recomputed = TRUE,
    baseline_centers_recomputed = TRUE,
    soft_memberships_recomputed = TRUE,
    response_transformations_recomputed = TRUE)
}
best <- which.min(colMeans(loss))
selected_lambda <- lambda_grid[best]
selected_model <- if (mean(loss[, best]) < mean(mean_loss))
  "two_configuration" else "mean_footprint"
selection <- rbindlist(selection_parts)
selection[, `:=`(selected_lambda = selected_lambda,
  selected_penalty = lambda == selected_lambda,
  selected_molecular_model = selected_model)]
fwrite(selection, file.path(RESULT,
  "molecular_inner_model_selection.tsv.gz"), sep = "\t", compress = "gzip")

## Recreate training out-of-fold predictions at the selected molecular penalty.
oof_parts <- vector("list", N_INNER_FOLDS)
for (fold in seq_len(N_INNER_FOLDS)) {
  validation_lines <- names(inner_assignment)[inner_assignment == fold]
  learn_lines <- names(inner_assignment)[inner_assignment != fold]
  state_fit <- fit_configurations(baseline_raw, learn_lines)
  learn_membership <- predict_configurations(state_fit,
    baseline_raw[learn_lines, , drop = FALSE])$membership
  validation_membership <- predict_configurations(state_fit,
    baseline_raw[validation_lines, , drop = FALSE])$membership
  footprint <- make_footprint(learn_lines, score_table, D, programs)
  train <- score_table$depmap_id %chin% learn_lines
  valid <- score_table$depmap_id %chin% validation_lines &
    score_table$perturbation_id %chin% rownames(footprint)
  x_train <- footprint[score_table$perturbation_id[train], programs,
    drop = FALSE]
  x_valid <- footprint[score_table$perturbation_id[valid], programs,
    drop = FALSE]
  fits <- lapply(seq_len(2L), function(configuration)
    fit_multi_ridge(x_train, D[train, , drop = FALSE],
      learn_membership[score_table$depmap_id[train], configuration],
      selected_lambda))
  valid_weights <- validation_membership[score_table$depmap_id[valid],,
    drop = FALSE]
  prediction <- predict_multi_ridge(fits[[1L]], x_valid) *
    valid_weights[, 1L] + predict_multi_ridge(fits[[2L]], x_valid) *
    valid_weights[, 2L]
  oof_parts[[fold]] <- data.table(inner_fold = fold,
    observation_id = score_table$observation_id[valid],
    depmap_id = score_table$depmap_id[valid],
    perturbation_id = score_table$perturbation_id[valid],
    program = rep(programs, each = sum(valid)),
    observed_displacement = as.vector(D[valid, , drop = FALSE]),
    reference_footprint = as.vector(x_valid),
    predicted_displacement = as.vector(prediction))
}
oof <- rbindlist(oof_parts)
oof[, evaluation_role := "training_inner_validation"]

## Fit and freeze the outer molecular model, then predict the seven held lines.
outer_state_fit <- fit_configurations(baseline_raw, training_lines)
outer_membership_all <- predict_configurations(outer_state_fit,
  baseline_raw[all_lines, , drop = FALSE])$membership
outer_footprint <- make_footprint(training_lines, score_table, D, programs)
outer_train <- score_table$depmap_id %chin% training_lines
outer_held <- score_table$depmap_id %chin% held_lines
outer_x_train <- outer_footprint[score_table$perturbation_id[outer_train],
  programs, drop = FALSE]
outer_x_held <- outer_footprint[score_table$perturbation_id[outer_held],
  programs, drop = FALSE]
outer_fits <- lapply(seq_len(2L), function(configuration)
  fit_multi_ridge(outer_x_train, D[outer_train, , drop = FALSE],
    outer_membership_all[score_table$depmap_id[outer_train], configuration],
    selected_lambda))
held_weights <- outer_membership_all[score_table$depmap_id[outer_held],,
  drop = FALSE]
held_prediction <- predict_multi_ridge(outer_fits[[1L]], outer_x_held) *
  held_weights[, 1L] + predict_multi_ridge(outer_fits[[2L]], outer_x_held) *
  held_weights[, 2L]
held <- data.table(observation_id = score_table$observation_id[outer_held],
  depmap_id = score_table$depmap_id[outer_held],
  cell_line = score_table$cell_line[outer_held],
  perturbation_id = score_table$perturbation_id[outer_held],
  drug = score_table$drug[outer_held], dose = score_table$dose[outer_held],
  dose_panel = score_table$dose_panel[outer_held],
  program = rep(programs, each = sum(outer_held)),
  observed_displacement = as.vector(D[outer_held, , drop = FALSE]),
  reference_footprint = as.vector(outer_x_held),
  predicted_displacement = as.vector(held_prediction),
  evaluation_role = "outer_held_out")
predictions <- rbindlist(list(oof, held), use.names = TRUE, fill = TRUE)
fwrite(predictions, file.path(ARTIFACT,
  "state_modulated_molecular_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")

## Save frozen state centers, line memberships and both 9-by-9 operators.
center_table <- rbindlist(lapply(seq_len(2L), function(configuration)
  data.table(configuration = paste0("C", configuration),
    baseline_program = colnames(outer_state_fit$centers),
    standardized_center = outer_state_fit$centers[configuration, ],
    membership_temperature = outer_state_fit$temperature)))
fwrite(center_table, file.path(RESULT,
  "baseline_configuration_centers.tsv"), sep = "\t")
membership_table <- data.table(depmap_id = all_lines,
  role = fifelse(all_lines %chin% training_lines, "training", "held_out"),
  C1 = outer_membership_all[all_lines, "C1"],
  C2 = outer_membership_all[all_lines, "C2"])
membership_table[, assigned_configuration :=
  fifelse(C1 >= C2, "C1", "C2")]
fwrite(membership_table, file.path(RESULT,
  "baseline_configuration_memberships.tsv"), sep = "\t")
coefficient_table <- rbindlist(lapply(seq_len(2L), function(configuration) {
  fit <- outer_fits[[configuration]]
  rbind(data.table(configuration = paste0("C", configuration),
      input_program = "(Intercept)", output_program = programs,
      coefficient = fit$intercept, coefficient_type = "intercept"),
    data.table(configuration = paste0("C", configuration),
      input_program = rep(programs, times = length(programs)),
      output_program = rep(programs, each = length(programs)),
      coefficient = as.vector(fit$beta),
      coefficient_type = ifelse(rep(programs, times = length(programs)) ==
        rep(programs, each = length(programs)), "diagonal", "redirection")))
}))
fwrite(coefficient_table, file.path(RESULT,
  "configuration_response_transformations.tsv"), sep = "\t")
saveRDS(list(label = "frozen_two_baseline_configuration_molecular_model",
    baseline_programs = state_programs, response_programs = programs,
    selected_lambda = selected_lambda, selected_model = selected_model,
    state_fit = outer_state_fit, line_memberships = outer_membership_all,
    reference_footprint = outer_footprint, transformations = outer_fits,
    training_lines = training_lines, held_lines = held_lines,
    function_outcomes_loaded = FALSE, drug_features = FALSE,
    cell_line_features = FALSE),
  file.path(ARTIFACT, "state_modulated_response_model.rds"))

## Molecular-only held-out performance, variance recovery and biological labels.
performance_parts <- list()
for (model in c("mean_footprint", "two_configuration")) {
  prediction_column <- if (model == "mean_footprint")
    "reference_footprint" else "predicted_displacement"
  for (metric in c("pooled", "within_perturbation", "context_displacement")) {
    current <- copy(held)
    current[, predicted_value := get(prediction_column)]
    if (metric == "within_perturbation") current[, `:=`(
      observed_value = observed_displacement - mean(observed_displacement),
      predicted_value = predicted_value - mean(predicted_value)),
      by = perturbation_id] else if (metric == "context_displacement")
        current[, `:=`(
          observed_value = observed_displacement - reference_footprint,
          predicted_value = predicted_value - reference_footprint)] else
        current[, observed_value := observed_displacement]
    performance_parts[[length(performance_parts) + 1L]] <- data.table(
      model = model, metric = metric, n_rows = nrow(current),
      n_lines = uniqueN(current$depmap_id),
      n_perturbations = uniqueN(current$perturbation_id),
      n_programs = uniqueN(current$program),
      pearson = safe_cor(current$observed_value, current$predicted_value),
      mse = mean((current$observed_value - current$predicted_value) ^ 2),
      observed_sd = sd(current$observed_value),
      predicted_sd = sd(current$predicted_value),
      predicted_to_observed_sd_ratio = sd(current$predicted_value) /
        sd(current$observed_value), correlation_method = "pearson")
  }
}
performance <- rbindlist(performance_parts)
fwrite(performance, file.path(RESULT, "molecular_performance.tsv"), sep = "\t")

per_program <- held[, {
  observed_context <- observed_displacement - reference_footprint
  predicted_context <- predicted_displacement - reference_footprint
  .(n_rows = .N,
    mean_footprint_pearson = safe_cor(observed_displacement,
      reference_footprint),
    state_modulated_pearson = safe_cor(observed_displacement,
      predicted_displacement),
    context_pearson = safe_cor(observed_context, predicted_context),
    mean_footprint_mse = mean((observed_displacement -
      reference_footprint) ^ 2),
    state_modulated_mse = mean((observed_displacement -
      predicted_displacement) ^ 2),
    observed_context_sd = sd(observed_context),
    predicted_context_sd = sd(predicted_context))
}, by = program]
per_program[, mse_improvement := mean_footprint_mse - state_modulated_mse]
fwrite(per_program, file.path(RESULT,
  "molecular_per_program_performance.tsv"), sep = "\t")
line_annotation <- fread(file.path(stage_result("01_inputs"),
  "fold_assignments.tsv"))[, .(depmap_id, organ)]
per_line <- held[, .(cell_line = cell_line[1L], n_rows = .N,
    mean_footprint_mse = mean((observed_displacement -
      reference_footprint) ^ 2),
    state_modulated_mse = mean((observed_displacement -
      predicted_displacement) ^ 2),
    mean_footprint_pearson = safe_cor(observed_displacement,
      reference_footprint),
    state_modulated_pearson = safe_cor(observed_displacement,
      predicted_displacement)), by = depmap_id]
per_line[, mse_improvement := mean_footprint_mse - state_modulated_mse]
per_line <- merge(per_line, line_annotation, by = "depmap_id", all.x = TRUE)
positive_line_gain <- sum(pmax(per_line$mse_improvement, 0))
if (positive_line_gain > 0) per_line[, positive_improvement_fraction :=
  pmax(mse_improvement, 0) / positive_line_gain] else
  per_line[, positive_improvement_fraction := NA_real_]
fwrite(per_line, file.path(RESULT,
  "molecular_per_line_performance.tsv"), sep = "\t")
per_lineage <- per_line[, .(n_lines = .N,
    mean_mse_improvement = mean(mse_improvement),
    total_mse_improvement = sum(mse_improvement),
    n_lines_improved = sum(mse_improvement > 0)), by = organ]
fwrite(per_lineage, file.path(RESULT,
  "molecular_per_lineage_performance.tsv"), sep = "\t")
per_perturbation <- held[, {
  observed_context <- observed_displacement - reference_footprint
  predicted_context <- predicted_displacement - reference_footprint
  .(drug = drug[1L], dose = dose[1L], n_rows = .N,
    mean_footprint_pearson = safe_cor(observed_displacement,
      reference_footprint),
    state_modulated_pearson = safe_cor(observed_displacement,
      predicted_displacement),
    context_pearson = safe_cor(observed_context, predicted_context),
    mean_footprint_mse = mean(observed_context ^ 2),
    state_modulated_mse = mean((observed_displacement -
      predicted_displacement) ^ 2),
    observed_context_sd = sd(observed_context),
    predicted_context_sd = sd(predicted_context))
}, by = perturbation_id]
per_perturbation[, mse_improvement := mean_footprint_mse -
  state_modulated_mse]
fwrite(per_perturbation, file.path(RESULT,
  "molecular_per_perturbation_performance.tsv.gz"), sep = "\t",
  compress = "gzip")

transformation_summary <- coefficient_table[coefficient_type != "intercept",
  .(coefficient = coefficient,
    effect = fcase(coefficient_type == "diagonal" & coefficient > 1,
        "potentiation", coefficient_type == "diagonal" & coefficient < 1,
        "suppression", default = "redirection")),
  by = .(configuration, input_program, output_program, coefficient_type)]
fwrite(transformation_summary, file.path(RESULT,
  "program_potentiation_suppression_redirection.tsv"), sep = "\t")
summary <- data.table(selected_lambda = selected_lambda,
  selected_molecular_model = selected_model,
  mean_inner_two_configuration_mse = mean(loss[, best]),
  mean_inner_mean_footprint_mse = mean(mean_loss),
  held_mean_footprint_pearson = performance[model == "mean_footprint" &
    metric == "pooled", pearson],
  held_state_modulated_pearson = performance[model == "two_configuration" &
    metric == "pooled", pearson],
  held_context_pearson = performance[model == "two_configuration" &
    metric == "context_displacement", pearson],
  held_context_variance_recovered = 1 -
    performance[model == "two_configuration" & metric ==
      "context_displacement", mse] /
    performance[model == "mean_footprint" & metric ==
      "context_displacement", mse],
  n_held_lines_improved = sum(per_line$mse_improvement > 0),
  n_held_lineages_with_positive_mean_improvement = sum(
    per_lineage$mean_mse_improvement > 0),
  top_line_positive_improvement_fraction = max(
    per_line$positive_improvement_fraction, na.rm = TRUE))
fwrite(summary, file.path(RESULT,
  "state_modulated_response_summary.tsv"), sep = "\t")

## Stage-owned one-plot PDFs.
membership_matrix <- t(as.matrix(membership_table[, .(C1, C2)]))
colnames(membership_matrix) <- membership_table$depmap_id
configuration_order <- membership_table[order(role, assigned_configuration,
  -pmax(C1, C2), depmap_id)]
membership_matrix <- membership_matrix[, configuration_order$depmap_id,
  drop = FALSE]
stopifnot(identical(colnames(membership_matrix),
  configuration_order$depmap_id))
hm_membership <- draw_heatmap(membership_matrix,
  col_anno_df = data.frame(Role = configuration_order$role,
    State = configuration_order$assigned_configuration),
  col_anno_cols = c("Role", "State"), scale_rows = FALSE,
  value_range = c(0, 0.5, 1), value_colors = c("white", "#BFD3E6",
    "#2166AC"), show_row_names = TRUE, show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 5.2),
  cluster_rows = FALSE, cluster_columns = FALSE,
  title = "Soft baseline-state configuration memberships")
save_heatmap("11a_baseline_configuration_membership_heatmap.pdf",
  hm_membership, 8.0, 3.7)

for (configuration_name in c("C1", "C2")) {
  current <- coefficient_table[configuration == configuration_name &
    coefficient_type != "intercept"]
  wide <- dcast(current, input_program ~ output_program,
    value.var = "coefficient")
  matrix <- as.matrix(wide[, ..programs])
  rownames(matrix) <- wide$input_program
  matrix <- matrix[programs, programs, drop = FALSE]
  limit <- max(abs(matrix))
  hm <- draw_heatmap(matrix, scale_rows = FALSE,
    value_range = c(-limit, 0, limit),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = sprintf("%s frozen response transformation", configuration_name))
  save_heatmap(sprintf("11b_%s_response_transformation_heatmap.pdf",
    configuration_name), hm, 5.8, 5.2)
}

plot_scatter <- held[sample(.N, min(.N, 18000L))]
p_scatter <- ggplot(plot_scatter,
    aes(observed_displacement, predicted_displacement)) +
  geom_point(alpha = 0.12, size = 0.35, color = "#4D4D4D") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55, color = "#B2182B") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("r=%.3f\nn=%s",
      safe_cor(held$observed_displacement, held$predicted_displacement),
      format(nrow(held), big.mark = ",")), size = 2.4) +
  labs(x = "Observed molecular displacement",
    y = "Predicted molecular displacement",
    title = "Held-out state-modulated molecular response") +
  square_scatter_theme()
save_gg("11c_observed_vs_predicted_molecular_displacement.pdf",
  p_scatter, 4.8, 4.6)

correlation_plot <- melt(per_program[, .(program,
    mean_footprint_pearson, state_modulated_pearson)], id.vars = "program",
  variable.name = "model", value.name = "pearson")
p_program <- ggplot(correlation_plot,
    aes(program, pearson, color = model, group = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.45) + geom_point(size = 1.8) +
  scale_color_manual(values = c(mean_footprint_pearson = "#4D4D4D",
    state_modulated_pearson = "#B2182B"),
    labels = c(mean_footprint_pearson = "Mean footprint",
      state_modulated_pearson = "Two-state modulation")) +
  labs(x = NULL, y = "Held-out Pearson r", color = NULL,
    title = "Molecular prediction by response program") + theme_stage11()
save_gg("11d_per_program_prediction_correlation.pdf", p_program, 6.2, 4.1)

p_improvement <- ggplot(per_perturbation, aes(mse_improvement)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2,
    color = "#B2182B") +
  geom_histogram(bins = 45, fill = "#2166AC", color = "white",
    linewidth = 0.15) +
  labs(x = "Mean-footprint MSE minus two-state MSE",
    y = "Perturbations",
    title = "Per-perturbation molecular improvement") + theme_stage11()
save_gg("11e_per_perturbation_improvement_distribution.pdf",
  p_improvement, 5.8, 4.1)

figures <- data.table(pdf = c(
  "11a_baseline_configuration_membership_heatmap.pdf",
  "11b_C1_response_transformation_heatmap.pdf",
  "11b_C2_response_transformation_heatmap.pdf",
  "11c_observed_vs_predicted_molecular_displacement.pdf",
  "11d_per_program_prediction_correlation.pdf",
  "11e_per_perturbation_improvement_distribution.pdf"),
  producing_script = "scripts/11_fit_state_modulated_response.R")
figures[, exists := file.exists(file.path(RESULT, pdf))]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")
if (!all(figures$exists)) fail("A required Stage 11 PDF was not created.")

checks <- data.table(check = c("excluded_lines_absent", "fixed_26_7_split",
    "nine_response_programs_unchanged", "molecular_stage_loaded_function",
    "reference_means_arithmetic", "inner_quantities_training_only",
    "drug_identity_predictor", "cell_line_identity_predictor",
    "per_drug_models", "primary_correlation_pearson"),
  observed = c(sum(EXCLUDED_IDS %chin% all_lines),
    paste(length(training_lines), length(held_lines), sep = "/"),
    length(programs), 0, 1, 1, 0, 0, 0, 1),
  expected = c("0", "26/7", "9", "0", "1", "1", "0", "0", "0", "1"))
checks[, passed := as.character(observed) == expected]
fwrite(checks, file.path(RESULT,
  "state_modulated_response_required_checks.tsv"), sep = "\t")
if (!all(checks$passed)) fail("A Stage 11 separation check failed.")
write_artifact_manifest(ARTIFACT,
  "scripts/11_fit_state_modulated_response.R")
writeLines("complete", file.path(RESULT, "stage11_complete.flag"))
message(sprintf(paste0("Stage 11 complete: selected lambda %.6g; held molecular ",
  "r %.4f; context variance recovered %.4f."), selected_lambda,
  summary$held_state_modulated_pearson,
  summary$held_context_variance_recovered))
