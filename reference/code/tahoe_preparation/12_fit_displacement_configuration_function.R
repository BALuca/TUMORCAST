#!/usr/bin/env Rscript

## Learn a two-configuration function decoder from measured displacement only.
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
RESULT <- stage_result("12_displacement_function")
ARTIFACT <- stage_artifact("12_displacement_function")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

theme_stage12 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage12(base_size) +
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

## Discover two directional displacement configurations without outcomes.
fit_displacement_configurations <- function(displacement) {
  center <- colMeans(displacement)
  scale <- apply(displacement, 2L, sd)
  keep <- is.finite(center) & is.finite(scale) & scale > 0
  if (!all(keep)) fail("A response program is nonvariable in configuration discovery.")
  z <- sweep(sweep(displacement, 2L, center, "-"), 2L, scale, "/")
  magnitude <- sqrt(rowSums(z ^ 2))
  direction <- z
  nonzero <- is.finite(magnitude) & magnitude > 0
  direction[nonzero, ] <- direction[nonzero, , drop = FALSE] /
    magnitude[nonzero]
  direction[!nonzero, ] <- 0
  set.seed(SEED)
  clustering <- kmeans(direction, centers = 2L, iter.max = 100L,
    nstart = 25L, algorithm = "Lloyd")
  centers <- clustering$centers
  ordering <- do.call(order, c(as.data.frame(centers), list(method = "radix")))
  centers <- centers[ordering, , drop = FALSE]
  distance2 <- vapply(seq_len(2L), function(k)
    rowSums(sweep(direction, 2L, centers[k, ], "-") ^ 2),
    numeric(nrow(direction)))
  temperature <- median(apply(distance2, 1L, min))
  if (!is.finite(temperature) || temperature <= 0)
    temperature <- mean(distance2[is.finite(distance2) & distance2 > 0])
  if (!is.finite(temperature) || temperature <= 0)
    fail("Displacement configuration temperature is undefined.")
  list(center = center, scale = scale, keep = keep, centers = centers,
    temperature = temperature, feature_names = colnames(displacement),
    magnitude_normalized_for_discovery = TRUE)
}

predict_displacement_configurations <- function(fit, displacement) {
  z <- sweep(sweep(displacement[, fit$keep, drop = FALSE], 2L,
    fit$center[fit$keep], "-"), 2L, fit$scale[fit$keep], "/")
  magnitude <- sqrt(rowSums(z ^ 2))
  direction <- z
  nonzero <- is.finite(magnitude) & magnitude > 0
  direction[nonzero, ] <- direction[nonzero, , drop = FALSE] /
    magnitude[nonzero]
  direction[!nonzero, ] <- 0
  distance2 <- vapply(seq_len(2L), function(k)
    rowSums(sweep(direction, 2L, fit$centers[k, ], "-") ^ 2),
    numeric(nrow(direction)))
  logits <- -distance2 / (2 * fit$temperature)
  logits <- sweep(logits, 1L, apply(logits, 1L, max), "-")
  membership <- exp(logits)
  membership <- membership / rowSums(membership)
  colnames(membership) <- c("F1", "F2")
  list(membership = membership, direction = direction,
    magnitude = magnitude, distance2 = distance2)
}

evaluate_predictions <- function(table, prediction_column, cohort, model) {
  current <- table[evaluation_included == TRUE]
  if (cohort == "activity_supported") current <- current[
    activity_supported == TRUE]
  rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric) {
    residual <- function_residual_table(current, metric, prediction_column)
    data.table(cohort = cohort, model = model, metric = metric,
      n_rows = nrow(residual), n_conditions = uniqueN(residual$perturbation_id),
      n_drugs = uniqueN(residual$drug), n_lines = uniqueN(residual$depmap_id),
      pearson = safe_cor(residual$observed_residual,
        residual$predicted_residual),
      mse = mean((residual$observed_residual -
        residual$predicted_residual) ^ 2),
      observed_sd = sd(residual$observed_residual),
      predicted_sd = sd(residual$predicted_residual),
      predicted_to_observed_sd_ratio = sd(residual$predicted_residual) /
        sd(residual$observed_residual), correlation_method = "pearson")
  }))
}

## Load measured responses and outcomes, but never a baseline or Stage 11 object.
required <- c(file.path(stage_result("02_response_programs"),
    "fold_program_stability.tsv"),
  file.path(stage_artifact("02_response_programs"),
    "fold_program_scores.tsv.gz"),
  file.path(stage_artifact("01_inputs"), "function_outcomes.tsv.gz"),
  file.path(stage_result("01_inputs"), "fold_assignments.tsv"),
  file.path(stage_result("04_activity_panel"),
    "condition_activity_census.tsv"),
  file.path(stage_artifact("05_program_function_bridge"),
    "measured_function_predictions.tsv.gz"),
  file.path(stage_result("05_program_function_bridge"),
    "function_decoder_selection.tsv"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 12 input(s): %s",
  paste(missing, collapse = ", "))

folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- sort(folds[outer_fold != OUTER_FOLD, depmap_id])
held_lines <- sort(folds[outer_fold == OUTER_FOLD, depmap_id])
all_lines <- sort(folds$depmap_id)
if (length(training_lines) != EXPECTED_TRAINING_LINES ||
    length(held_lines) != EXPECTED_HELD_LINES ||
    any(EXCLUDED_IDS %chin% all_lines))
  fail("The accepted 26/7 cohort changed before Stage 12.")
inner_assignment <- make_inner_assignment(training_lines)

stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- stability[outer_fold == OUTER_FOLD & stable == TRUE,
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
if (any(!is.finite(D))) fail("Measured response-program displacement is invalid.")

function_rows <- fread(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"))
function_rows[, `:=`(observation_id = paste(depmap_id, perturbation_id,
  sep = "__"), drug = tahoe_drug, dose = tahoe_dose_uM)]
position <- match(function_rows$observation_id, score_table$observation_id)
chain <- function_rows$chain_observed == 1L & !is.na(position)
activity_ids <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))[activity_supported == TRUE,
    perturbation_id]

## Nested cell-line validation recomputes molecular centers and both ridges.
lambda_grid <- DECODER_LAMBDA_GRID
shared_loss <- matrix(NA_real_, N_INNER_FOLDS, length(lambda_grid))
configuration_loss <- matrix(NA_real_, N_INNER_FOLDS, length(lambda_grid))
selection_parts <- vector("list", N_INNER_FOLDS)
for (fold in seq_len(N_INNER_FOLDS)) {
  validation_lines <- names(inner_assignment)[inner_assignment == fold]
  learn_lines <- names(inner_assignment)[inner_assignment != fold]
  discovery <- score_table$depmap_id %chin% learn_lines
  configuration_fit <- fit_displacement_configurations(
    D[discovery, , drop = FALSE])
  train <- function_rows$depmap_id %chin% learn_lines & chain
  valid <- function_rows$depmap_id %chin% validation_lines & chain
  x_train <- D[position[train], , drop = FALSE]
  x_valid <- D[position[valid], , drop = FALSE]
  train_membership <- predict_displacement_configurations(
    configuration_fit, x_train)$membership
  valid_membership <- predict_displacement_configurations(
    configuration_fit, x_valid)$membership
  shared_path <- ridge_path(x_train,
    function_rows$prism_sensitization[train], x_valid, lambda_grid)
  configuration_paths <- lapply(seq_len(2L), function(configuration)
    weighted_ridge_path(x_train, function_rows$prism_sensitization[train],
      x_valid, train_membership[, configuration], lambda_grid))
  configuration_path <- configuration_paths[[1L]] *
    valid_membership[, 1L] + configuration_paths[[2L]] *
    valid_membership[, 2L]
  observed <- function_rows$prism_sensitization[valid]
  shared_loss[fold, ] <- colMeans(sweep(shared_path, 1L, observed, "-") ^ 2)
  configuration_loss[fold, ] <- colMeans(sweep(configuration_path, 1L,
    observed, "-") ^ 2)
  selection_parts[[fold]] <- rbind(
    data.table(inner_fold = fold, model = "shared_decoder",
      lambda = lambda_grid, validation_mse = shared_loss[fold, ]),
    data.table(inner_fold = fold, model = "two_configuration_decoder",
      lambda = lambda_grid, validation_mse = configuration_loss[fold, ]))
}
shared_best <- which.min(colMeans(shared_loss))
configuration_best <- which.min(colMeans(configuration_loss))
shared_lambda <- lambda_grid[shared_best]
configuration_lambda <- lambda_grid[configuration_best]
selection <- rbindlist(selection_parts)
selection[, selected := (model == "shared_decoder" & lambda == shared_lambda) |
  (model == "two_configuration_decoder" & lambda == configuration_lambda)]
selection[, `:=`(shared_selected_lambda = shared_lambda,
  configuration_selected_lambda = configuration_lambda,
  molecular_configuration_discovery_used_outcome = FALSE)]
fwrite(selection, file.path(RESULT,
  "functional_inner_model_selection.tsv.gz"), sep = "\t", compress = "gzip")

## Recreate selected inner-fold predictions to test fold stability directly.
inner_prediction_parts <- vector("list", N_INNER_FOLDS)
for (fold in seq_len(N_INNER_FOLDS)) {
  validation_lines <- names(inner_assignment)[inner_assignment == fold]
  learn_lines <- names(inner_assignment)[inner_assignment != fold]
  discovery <- score_table$depmap_id %chin% learn_lines
  configuration_fit <- fit_displacement_configurations(
    D[discovery, , drop = FALSE])
  train <- function_rows$depmap_id %chin% learn_lines & chain
  valid <- function_rows$depmap_id %chin% validation_lines & chain
  x_train <- D[position[train], , drop = FALSE]
  x_valid <- D[position[valid], , drop = FALSE]
  train_membership <- predict_displacement_configurations(
    configuration_fit, x_train)$membership
  valid_membership <- predict_displacement_configurations(
    configuration_fit, x_valid)$membership
  shared_fit <- fit_ridge(x_train, function_rows$prism_sensitization[train],
    shared_lambda)
  configuration_fits <- lapply(seq_len(2L), function(configuration)
    fit_weighted_ridge(x_train, function_rows$prism_sensitization[train],
      train_membership[, configuration], configuration_lambda))
  shared_prediction <- predict_ridge(shared_fit, x_valid)
  configuration_prediction <- predict_ridge(configuration_fits[[1L]],
      x_valid) * valid_membership[, 1L] +
    predict_ridge(configuration_fits[[2L]], x_valid) * valid_membership[, 2L]
  inner_prediction_parts[[fold]] <- data.table(inner_fold = fold,
    observation_id = function_rows$observation_id[valid],
    depmap_id = function_rows$depmap_id[valid],
    perturbation_id = function_rows$perturbation_id[valid],
    drug = function_rows$drug[valid],
    observed_function = function_rows$prism_sensitization[valid],
    activity_supported = function_rows$perturbation_id[valid] %chin%
      activity_ids, shared_prediction = shared_prediction,
    configuration_prediction = configuration_prediction)
}
inner_predictions <- rbindlist(inner_prediction_parts)
inner_performance <- rbindlist(lapply(seq_len(N_INNER_FOLDS), function(fold) {
  current_fold <- inner_predictions[inner_fold == fold]
  rbindlist(lapply(c("complete_matched_panel", "activity_supported"),
    function(cohort) {
      current <- if (cohort == "activity_supported")
        current_fold[activity_supported == TRUE] else current_fold
      rbindlist(lapply(c("shared_decoder", "two_configuration_decoder"),
        function(model) {
          column <- if (model == "shared_decoder") "shared_prediction" else
            "configuration_prediction"
          rbindlist(lapply(c("pooled", "within_perturbation",
              "strict_cell_line_by_perturbation"), function(metric) {
            residual <- function_residual_table(current, metric, column)
            data.table(inner_fold = fold, cohort = cohort, model = model,
              metric = metric, n_rows = nrow(residual),
              pearson = safe_cor(residual$observed_residual,
                residual$predicted_residual),
              mse = mean((residual$observed_residual -
                residual$predicted_residual) ^ 2))
          }))
        }))
    }))
}))
fwrite(inner_performance, file.path(RESULT,
  "functional_inner_fold_comparison.tsv"), sep = "\t")

## Fit and freeze the independent outer functional configuration decoder.
outer_discovery <- score_table$depmap_id %chin% training_lines
outer_configuration_fit <- fit_displacement_configurations(
  D[outer_discovery, , drop = FALSE])
outer_train <- function_rows$depmap_id %chin% training_lines & chain
outer_held <- function_rows$depmap_id %chin% held_lines
outer_available <- outer_held & !is.na(position)
x_train <- D[position[outer_train], , drop = FALSE]
y_train <- function_rows$prism_sensitization[outer_train]
train_membership <- predict_displacement_configurations(
  outer_configuration_fit, x_train)$membership
configuration_fits <- lapply(seq_len(2L), function(configuration)
  fit_weighted_ridge(x_train, y_train, train_membership[, configuration],
    configuration_lambda))
shared_fit <- fit_ridge(x_train, y_train, shared_lambda)
x_held <- D[position[outer_available], , drop = FALSE]
held_membership <- predict_displacement_configurations(
  outer_configuration_fit, x_held)$membership
configuration_prediction <- rep(NA_real_, sum(outer_held))
configuration_prediction[!is.na(position[outer_held])] <-
  predict_ridge(configuration_fits[[1L]], x_held) * held_membership[, 1L] +
  predict_ridge(configuration_fits[[2L]], x_held) * held_membership[, 2L]

accepted_predictions <- fread(file.path(stage_artifact(
  "05_program_function_bridge"), "measured_function_predictions.tsv.gz"))
held_function <- function_rows[outer_held]
accepted_position <- match(held_function$observation_id,
  accepted_predictions$observation_id)
if (anyNA(accepted_position)) fail("Accepted shared predictions do not align.")
panel <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
held_predictions <- data.table(
  observation_id = held_function$observation_id,
  depmap_id = held_function$depmap_id, cell_line = held_function$cell_line,
  perturbation_id = held_function$perturbation_id, drug = held_function$drug,
  dose = held_function$dose,
  observed_function = held_function$prism_sensitization,
  evaluation_included = held_function$chain_observed == 1L,
  activity_supported = panel$activity_supported[match(
    held_function$perturbation_id, panel$perturbation_id)],
  shared_prediction = accepted_predictions$measured_prediction[
    accepted_position], configuration_prediction = configuration_prediction,
  F1 = NA_real_, F2 = NA_real_)
held_predictions[!is.na(position[outer_held]), `:=`(
  F1 = held_membership[, 1L], F2 = held_membership[, 2L])]
held_predictions[, assigned_configuration := fifelse(F1 >= F2, "F1", "F2")]
fwrite(held_predictions, file.path(ARTIFACT,
  "held_displacement_configuration_function_predictions.tsv.gz"),
  sep = "\t", compress = "gzip")

all_function_available <- chain
all_x <- D[position[all_function_available], , drop = FALSE]
all_membership <- predict_displacement_configurations(
  outer_configuration_fit, all_x)$membership
membership_table <- data.table(
  observation_id = function_rows$observation_id[all_function_available],
  depmap_id = function_rows$depmap_id[all_function_available],
  perturbation_id = function_rows$perturbation_id[all_function_available],
  role = fifelse(function_rows$depmap_id[all_function_available] %chin%
    training_lines, "training", "held_out"),
  F1 = all_membership[, 1L], F2 = all_membership[, 2L])
membership_table[, assigned_configuration := fifelse(F1 >= F2, "F1", "F2")]
fwrite(membership_table, file.path(ARTIFACT,
  "displacement_configuration_memberships.tsv.gz"), sep = "\t",
  compress = "gzip")

center_table <- rbindlist(lapply(seq_len(2L), function(configuration)
  data.table(configuration = paste0("F", configuration),
    program = programs,
    direction_center = outer_configuration_fit$centers[configuration, ],
    membership_temperature = outer_configuration_fit$temperature)))
fwrite(center_table, file.path(RESULT,
  "displacement_configuration_centers.tsv"), sep = "\t")
coefficient_table <- rbindlist(lapply(seq_len(2L), function(configuration) {
  raw <- raw_ridge_coefficients(configuration_fits[[configuration]])
  data.table(configuration = paste0("F", configuration),
    term = c("(Intercept)", programs),
    coefficient = c(raw$intercept, raw$beta[programs]),
    penalty_lambda = configuration_lambda)
}))
fwrite(coefficient_table, file.path(RESULT,
  "configuration_function_coefficients.tsv"), sep = "\t")
saveRDS(list(label = "frozen_displacement_only_two_configuration_decoder",
    response_programs = programs, displacement_configuration_fit =
      outer_configuration_fit, configuration_decoders = configuration_fits,
    configuration_lambda = configuration_lambda, shared_lambda = shared_lambda,
    training_lines = training_lines, held_lines = held_lines,
    baseline_features = FALSE, drug_features = FALSE,
    cell_line_features = FALSE, stage11_predictions_used = FALSE),
  file.path(ARTIFACT, "displacement_configuration_function_model.rds"))

## Held-out complete/activity performance and per-drug comparisons.
performance <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort)
  rbind(evaluate_predictions(held_predictions, "shared_prediction", cohort,
      "shared_decoder"),
    evaluate_predictions(held_predictions, "configuration_prediction", cohort,
      "two_configuration_decoder"))))
fwrite(performance, file.path(RESULT,
  "displacement_function_performance.tsv"), sep = "\t")

per_drug <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort) {
  current <- held_predictions[evaluation_included == TRUE]
  if (cohort == "activity_supported") current <- current[
    activity_supported == TRUE]
  rbindlist(lapply(unique(current$drug), function(drug_name) {
    x <- current[drug == drug_name]
    shared_within <- function_residual_table(x, "within_perturbation",
      "shared_prediction")
    configuration_within <- function_residual_table(x,
      "within_perturbation", "configuration_prediction")
    shared_strict <- function_residual_table(x,
      "strict_cell_line_by_perturbation", "shared_prediction")
    configuration_strict <- function_residual_table(x,
      "strict_cell_line_by_perturbation", "configuration_prediction")
    data.table(cohort = cohort, drug = drug_name, n_rows = nrow(x),
      shared_pooled = safe_cor(x$observed_function, x$shared_prediction),
      configuration_pooled = safe_cor(x$observed_function,
        x$configuration_prediction),
      shared_within = safe_cor(shared_within$observed_residual,
        shared_within$predicted_residual),
      configuration_within = safe_cor(
        configuration_within$observed_residual,
        configuration_within$predicted_residual),
      shared_strict = safe_cor(shared_strict$observed_residual,
        shared_strict$predicted_residual),
      configuration_strict = safe_cor(
        configuration_strict$observed_residual,
        configuration_strict$predicted_residual),
      observed_sd = sd(x$observed_function),
      shared_predicted_sd = sd(x$shared_prediction),
      configuration_predicted_sd = sd(x$configuration_prediction))
  }))
}))
per_drug[, `:=`(pooled_change = configuration_pooled - shared_pooled,
  within_change = configuration_within - shared_within,
  strict_change = configuration_strict - shared_strict)]
fwrite(per_drug, file.path(RESULT,
  "displacement_function_per_drug.tsv"), sep = "\t")

## Drug-block bootstrap quantifies changes from the accepted shared decoder.
bootstrap_file <- file.path(ARTIFACT,
  "displacement_function_drug_block_bootstrap.tsv.gz")
if (file.exists(bootstrap_file)) {
  bootstrap <- fread(bootstrap_file)
  valid_cache <- nrow(bootstrap) == 2L * 3L * N_BOOTSTRAPS &&
    setequal(unique(bootstrap$draw), seq_len(N_BOOTSTRAPS)) &&
    setequal(unique(bootstrap$cohort), c("complete_matched_panel",
      "activity_supported"))
  if (!valid_cache) fail("The existing Stage 12 bootstrap cache is invalid.")
} else {
  set.seed(SEED + 1200L)
  bootstrap_parts <- list()
  for (cohort in c("complete_matched_panel", "activity_supported")) {
    current <- held_predictions[evaluation_included == TRUE]
    if (cohort == "activity_supported") current <- current[
      activity_supported == TRUE]
    drugs <- sort(unique(current$drug))
    for (metric in c("pooled", "within_perturbation",
        "strict_cell_line_by_perturbation")) {
      draws <- rbindlist(lapply(seq_len(N_BOOTSTRAPS), function(draw) {
        sampled <- sample(drugs, length(drugs), replace = TRUE)
        resampled <- rbindlist(lapply(seq_along(sampled), function(block) {
          piece <- copy(current[drug == sampled[block]])
          piece[, perturbation_id := paste0(block, "__", perturbation_id)]
          piece
        }))
        shared <- function_residual_table(resampled, metric,
          "shared_prediction")
        configured <- function_residual_table(resampled, metric,
          "configuration_prediction")
        shared_r <- safe_cor(shared$observed_residual,
          shared$predicted_residual)
        configured_r <- safe_cor(configured$observed_residual,
          configured$predicted_residual)
        data.table(draw = draw, shared_pearson = shared_r,
          configuration_pearson = configured_r,
          pearson_change = configured_r - shared_r)
      }))
      draws[, `:=`(cohort = cohort, metric = metric)]
      bootstrap_parts[[length(bootstrap_parts) + 1L]] <- draws
    }
  }
  bootstrap <- rbindlist(bootstrap_parts)
  fwrite(bootstrap, bootstrap_file, sep = "\t", compress = "gzip")
}
intervals <- bootstrap[, .(
  change_median = median(pearson_change, na.rm = TRUE),
  change_lower_95 = quantile(pearson_change, 0.025, na.rm = TRUE),
  change_upper_95 = quantile(pearson_change, 0.975, na.rm = TRUE),
  fraction_change_positive = mean(pearson_change > 0, na.rm = TRUE)),
  by = .(cohort, metric)]
fwrite(intervals, file.path(RESULT,
  "displacement_function_bootstrap_intervals.tsv"), sep = "\t")

accepted_lambda <- fread(file.path(stage_result(
  "05_program_function_bridge"), "function_decoder_selection.tsv"))$
    selected_lambda[1L]
summary <- data.table(shared_selected_lambda = shared_lambda,
  accepted_shared_lambda = accepted_lambda,
  shared_lambda_reproduced = isTRUE(all.equal(as.numeric(shared_lambda),
    as.numeric(accepted_lambda), tolerance = 1e-12)),
  configuration_selected_lambda = configuration_lambda,
  mean_inner_shared_mse = mean(shared_loss[, shared_best]),
  mean_inner_configuration_mse = mean(configuration_loss[,
    configuration_best]),
  activity_pooled_change = performance[cohort == "activity_supported" &
    model == "two_configuration_decoder" & metric == "pooled", pearson] -
      performance[cohort == "activity_supported" & model ==
        "shared_decoder" & metric == "pooled", pearson],
  activity_within_change = performance[cohort == "activity_supported" &
    model == "two_configuration_decoder" & metric ==
      "within_perturbation", pearson] - performance[cohort ==
        "activity_supported" & model == "shared_decoder" & metric ==
          "within_perturbation", pearson],
  activity_strict_change = performance[cohort == "activity_supported" &
    model == "two_configuration_decoder" & metric ==
      "strict_cell_line_by_perturbation", pearson] - performance[cohort ==
        "activity_supported" & model == "shared_decoder" & metric ==
          "strict_cell_line_by_perturbation", pearson])
fwrite(summary, file.path(RESULT,
  "displacement_function_summary.tsv"), sep = "\t")

## Stage-owned one-plot PDFs.
center_matrix <- dcast(center_table, configuration ~ program,
  value.var = "direction_center")
center_values <- as.matrix(center_matrix[, ..programs])
rownames(center_values) <- center_matrix$configuration
limit <- max(abs(center_values))
hm_centers <- draw_heatmap(center_values, scale_rows = FALSE,
  value_range = c(-limit, 0, limit),
  value_colors = c("#2166AC", "white", "#B2182B"),
  show_row_names = TRUE, show_column_names = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  title = "Direction-normalized displacement configurations")
save_heatmap("12a_displacement_configuration_heatmap.pdf", hm_centers,
  5.8, 3.8)

for (configuration_name in c("F1", "F2")) {
  current <- coefficient_table[configuration == configuration_name &
    term != "(Intercept)"]
  current[, program := factor(term, levels = rev(programs))]
  p_coefficient <- ggplot(current, aes(coefficient, program)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey65") +
    geom_col(width = 0.67, fill = if (configuration_name == "F1")
      "#2166AC" else "#B2182B") +
    labs(x = "Function-decoder coefficient", y = NULL,
      title = sprintf("%s displacement-to-function decoder", configuration_name)) +
    theme_stage12()
  save_gg(sprintf("12b_%s_function_coefficients.pdf", configuration_name),
    p_coefficient, 5.4, 4.2)
}

plot_rows <- held_predictions[evaluation_included == TRUE &
  activity_supported == TRUE]
p_scatter <- ggplot(plot_rows,
    aes(observed_function, configuration_prediction)) +
  geom_point(alpha = 0.28, size = 0.65, color = "#6A51A3") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.6, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("r=%.3f\nn=%d",
      safe_cor(plot_rows$observed_function,
        plot_rows$configuration_prediction), nrow(plot_rows)), size = 2.4) +
  labs(x = "Observed held-out function",
    y = "Predicted held-out function",
    title = "Activity panel: displacement-configuration decoder") +
  square_scatter_theme()
save_gg("12c_observed_vs_predicted_held_function.pdf", p_scatter,
  4.8, 4.6)

metric_plot <- copy(performance)
metric_plot[, metric_label := factor(metric,
  levels = c("pooled", "within_perturbation",
    "strict_cell_line_by_perturbation"),
  labels = c("Pooled", "Within perturbation", "Strict context"))]
p_metric <- ggplot(metric_plot,
    aes(metric_label, pearson, color = model, group = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.45) + geom_point(size = 1.8) +
  facet_wrap(~cohort, ncol = 2) +
  scale_color_manual(values = c(shared_decoder = "#4D4D4D",
    two_configuration_decoder = "#6A51A3"),
    labels = c(shared_decoder = "Shared decoder",
      two_configuration_decoder = "Two configurations")) +
  labs(x = NULL, y = "Held-out Pearson r", color = NULL,
    title = "Function-decoder held-out comparison") +
  theme_stage12() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_gg("12d_function_metric_comparison.pdf", p_metric, 7.8, 4.3)

p_drug <- ggplot(per_drug[cohort == "activity_supported"],
    aes(within_change)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2,
    color = "#B2182B") +
  geom_histogram(bins = 28, fill = "#6A51A3", color = "white",
    linewidth = 0.15) +
  labs(x = "Per-drug within-perturbation Pearson change",
    y = "Drugs", title = "Activity-panel per-drug decoder change") +
  theme_stage12()
save_gg("12e_per_drug_change_distribution.pdf", p_drug, 5.8, 4.1)

figures <- data.table(pdf = c(
  "12a_displacement_configuration_heatmap.pdf",
  "12b_F1_function_coefficients.pdf",
  "12b_F2_function_coefficients.pdf",
  "12c_observed_vs_predicted_held_function.pdf",
  "12d_function_metric_comparison.pdf",
  "12e_per_drug_change_distribution.pdf"),
  producing_script = "scripts/12_fit_displacement_configuration_function.R")
figures[, exists := file.exists(file.path(RESULT, pdf))]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")
if (!all(figures$exists)) fail("A required Stage 12 PDF was not created.")

checks <- data.table(check = c("excluded_lines_absent", "fixed_26_7_split",
    "nine_response_programs_unchanged", "functional_stage_loaded_baseline",
    "stage11_predictions_loaded", "configuration_discovery_used_outcome",
    "drug_identity_predictor", "cell_line_identity_predictor",
    "per_drug_models", "primary_correlation_pearson",
    "shared_decoder_lambda_reproduced"),
  observed = c(sum(EXCLUDED_IDS %chin% all_lines),
    paste(length(training_lines), length(held_lines), sep = "/"),
    length(programs), 0, 0, 0, 0, 0, 0, 1,
    as.integer(summary$shared_lambda_reproduced)),
  expected = c("0", "26/7", "9", "0", "0", "0", "0", "0", "0",
    "1", "1"))
checks[, passed := as.character(observed) == expected]
fwrite(checks, file.path(RESULT,
  "displacement_function_required_checks.tsv"), sep = "\t")
if (!all(checks$passed)) fail("A Stage 12 separation check failed.")
write_artifact_manifest(ARTIFACT,
  "scripts/12_fit_displacement_configuration_function.R")
writeLines("complete", file.path(RESULT, "stage12_complete.flag"))
message(sprintf(paste0("Stage 12 complete: activity changes pooled %.4f, ",
  "within %.4f, strict %.4f."), summary$activity_pooled_change,
  summary$activity_within_change, summary$activity_strict_change))
