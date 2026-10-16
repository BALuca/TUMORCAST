#!/usr/bin/env Rscript

## Fit one global continuous baseline-by-reference molecular susceptibility
## model.  This stage uses frozen molecular programs, frozen recurrent baseline
## programs and fold-specific frozen functional-axis coefficients.  It never
## loads functional outcomes and produces only molecular predictions.
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
RESULT <- stage_result("14_continuous_baseline_susceptibility")
ARTIFACT <- stage_artifact("14_continuous_baseline_susceptibility")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

theme_stage14 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage14(base_size) +
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
fmt <- function(x, digits = 3L) ifelse(is.finite(x),
  formatC(x, digits = digits, format = "f"), "NA")

make_footprint <- function(lines, score_table, displacement, programs) {
  use <- score_table$depmap_id %chin% lines
  means <- data.table(perturbation_id = score_table$perturbation_id[use],
    displacement[use, programs, drop = FALSE])[, lapply(.SD, mean),
      by = perturbation_id, .SDcols = programs]
  matrix <- as.matrix(means[, ..programs])
  rownames(matrix) <- means$perturbation_id
  matrix
}

standardize_baseline <- function(raw_scores, learning_lines, label) {
  center <- colMeans(raw_scores[learning_lines, , drop = FALSE])
  scale <- apply(raw_scores[learning_lines, , drop = FALSE], 2L, sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0))
    fail("%s baseline programs cannot be standardized.", label)
  z <- sweep(sweep(raw_scores, 2L, center, "-"), 2L, scale, "/")
  list(scores = z, center = center, scale = scale,
    learning_lines = learning_lines)
}

make_interactions <- function(line_ids, perturbation_ids, baseline_z,
                              footprint, baseline_programs, programs) {
  if (anyNA(match(line_ids, rownames(baseline_z))) ||
      anyNA(match(perturbation_ids, rownames(footprint))))
    fail("An interaction row lacks a frozen baseline score or footprint.")
  S <- baseline_z[line_ids, baseline_programs, drop = FALSE]
  R <- footprint[perturbation_ids, programs, drop = FALSE]
  pieces <- lapply(baseline_programs, function(baseline_program)
    sweep(R, 1L, S[, baseline_program], "*"))
  X <- do.call(cbind, pieces)
  colnames(X) <- unlist(lapply(baseline_programs, function(b)
    paste0(b, "__", programs)), use.names = FALSE)
  X
}

axis_vector <- function(axis, fold_label, programs) {
  current <- axis[fold == fold_label]
  if (nrow(current) != length(programs) ||
      !setequal(current$program, programs) ||
      uniqueN(current$ridge_penalty) != 1L)
    fail("Frozen functional axis is incomplete for %s.", fold_label)
  setNames(current$coefficient[match(programs, current$program)], programs)
}

molecular_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]; predicted <- predicted[keep]
  data.table(n_rows = length(observed), pearson = safe_cor(observed, predicted),
    observed_sd = sd(observed), predicted_sd = sd(predicted),
    predicted_to_observed_sd_ratio = sd(predicted) / sd(observed),
    mse = mean((observed - predicted) ^ 2),
    rmse = sqrt(mean((observed - predicted) ^ 2)),
    correlation_method = "pearson")
}

required <- c(
  file.path(stage_result("01_inputs"), "fold_assignments.tsv"),
  file.path(stage_result("02_response_programs"),
    "fold_program_stability.tsv"),
  file.path(stage_artifact("02_response_programs"),
    "fold_program_scores.tsv.gz"),
  file.path(stage_result("04_activity_panel"),
    "condition_activity_census.tsv"),
  file.path(stage_result("04_activity_panel"),
    "molecular_context_variance_census.tsv.gz"),
  file.path(stage_result("05_program_function_bridge"),
    "functional_molecular_axis_coefficients.tsv"),
  file.path(stage_artifact("10_baseline_states"),
    "baseline_program_scores.tsv.gz"),
  file.path(stage_result("10_baseline_states"),
    "baseline_program_stability.tsv"),
  file.path(stage_result("10_baseline_states"),
    "baseline_program_line_dominance.tsv"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 14 input(s): %s",
  paste(missing, collapse = ", "))

folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- sort(folds[outer_fold != OUTER_FOLD, depmap_id])
held_lines <- sort(folds[outer_fold == OUTER_FOLD, depmap_id])
all_lines <- sort(folds$depmap_id)
if (length(training_lines) != EXPECTED_TRAINING_LINES ||
    length(held_lines) != EXPECTED_HELD_LINES ||
    !identical(held_lines, sort(EXPECTED_HELD_IDS)) ||
    any(EXCLUDED_IDS %chin% all_lines))
  fail("The accepted 26/7 cohort changed before Stage 14.")
inner_assignment <- make_inner_assignment(training_lines)

baseline_stability <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_program_stability.tsv"))[, .SD[1L], by = program]
dominance <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_program_line_dominance.tsv"))[baseline_source == "tahoe_dmso"]
state_programs <- merge(baseline_stability[, .(program, stable)],
  dominance[, .(program, single_line_dominated,
    leave_one_line_out_stable, eligible_recurrent_state)], by = "program")[
      eligible_recurrent_state == TRUE, program]
state_programs <- paste0("S", sort(as.integer(sub("S", "", state_programs))))
if (!identical(state_programs, paste0("S", 2:8)))
  fail("Stage 14 requires exactly recurrent-eligible programs S2--S8.")
baseline_long <- fread(file.path(stage_artifact("10_baseline_states"),
  "baseline_program_scores.tsv.gz"))
baseline_wide <- dcast(baseline_long[baseline_source == "tahoe_dmso" &
    aggregation_level == "line_mean" & program %chin% state_programs],
  depmap_id ~ program, value.var = "raw_program_score")
baseline_raw <- as.matrix(baseline_wide[, ..state_programs])
rownames(baseline_raw) <- baseline_wide$depmap_id
baseline_raw <- baseline_raw[all_lines, state_programs, drop = FALSE]
if (any(!is.finite(baseline_raw)))
  fail("Frozen recurrent baseline-program scores are invalid.")

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

axis <- fread(file.path(stage_result("05_program_function_bridge"),
  "functional_molecular_axis_coefficients.tsv"))
axis <- axis[program %chin% programs]
if (!setequal(unique(axis$fold), c(paste0("inner_", 1:5),
    "outer_training")))
  fail("Stage 05 did not export every frozen fold-specific functional axis.")
fwrite(axis, file.path(RESULT,
  "functional_molecular_axis_coefficients.tsv"), sep = "\t")

eligibility <- fread(file.path(stage_result("04_activity_panel"),
  "molecular_context_variance_census.tsv.gz"))
forbidden_gate_fields <- grep("function|prism", names(eligibility),
  value = TRUE, ignore.case = TRUE)
allowed_audit_field <- "prism_used_for_molecular_classification"
if (length(setdiff(forbidden_gate_fields, allowed_audit_field)))
  fail("The Stage 14 molecular gate contains a functional field.")
if (any(eligibility[[allowed_audit_field]] != FALSE) ||
    any(eligibility$performance_used_for_selection != FALSE) ||
    any(eligibility$validation_or_held_values_used_for_selection != FALSE))
  fail("The Stage 14 eligibility table is not molecular-training-only.")
if (nrow(eligibility) != EXPECTED_DISCOVERY_CONDITIONS *
    (N_INNER_FOLDS + 1L) ||
    anyDuplicated(eligibility[, .(fold, perturbation_id)]))
  fail("Fold-specific molecular eligibility is incomplete.")
fwrite(eligibility, file.path(RESULT, "molecular_condition_eligibility.tsv"),
  sep = "\t")

## Nested validation: every footprint, baseline standardization, interaction
## standardization and molecular target is restricted to inner-learning lines.
lambda_grid <- DECODER_LAMBDA_GRID
selection_parts <- vector("list", N_INNER_FOLDS)
validation_store <- vector("list", N_INNER_FOLDS)
for (fold in seq_len(N_INNER_FOLDS)) {
  label <- paste0("inner_", fold)
  validation_lines <- names(inner_assignment)[inner_assignment == fold]
  learning_lines <- names(inner_assignment)[inner_assignment != fold]
  gate <- eligibility[fold == label & molecular_context_variable == TRUE,
    perturbation_id]
  footprint <- make_footprint(learning_lines, score_table, D, programs)
  baseline_fit <- standardize_baseline(baseline_raw, learning_lines, label)
  beta <- axis_vector(axis, label, programs)
  learning_rows <- score_table$depmap_id %chin% learning_lines &
    score_table$perturbation_id %chin% gate &
    score_table$perturbation_id %chin% rownames(footprint)
  validation_rows <- score_table$depmap_id %chin% validation_lines &
    score_table$perturbation_id %chin% gate &
    score_table$perturbation_id %chin% rownames(footprint)
  if (sum(learning_rows) < length(state_programs) * length(programs) ||
      sum(validation_rows) < 3L)
    fail("%s has too few molecularly eligible rows.", label)
  R_learning <- footprint[score_table$perturbation_id[learning_rows],
    programs, drop = FALSE]
  R_validation <- footprint[score_table$perturbation_id[validation_rows],
    programs, drop = FALSE]
  z_learning <- as.numeric((D[learning_rows, programs, drop = FALSE] -
    R_learning) %*% beta)
  z_validation <- as.numeric((D[validation_rows, programs, drop = FALSE] -
    R_validation) %*% beta)
  X_learning <- make_interactions(score_table$depmap_id[learning_rows],
    score_table$perturbation_id[learning_rows], baseline_fit$scores,
    footprint, state_programs, programs)
  X_validation <- make_interactions(score_table$depmap_id[validation_rows],
    score_table$perturbation_id[validation_rows], baseline_fit$scores,
    footprint, state_programs, programs)
  path <- ridge_path(X_learning, z_learning, X_validation, lambda_grid)
  fold_metrics <- rbindlist(lapply(seq_along(lambda_grid), function(index) {
    predicted <- path[, index]
    data.table(inner_fold = fold, fold = label,
      lambda = lambda_grid[index],
      validation_mse = mean((z_validation - predicted) ^ 2),
      validation_pearson = safe_cor(z_validation, predicted),
      validation_observed_sd = sd(z_validation),
      validation_predicted_sd = sd(predicted),
      validation_sd_ratio = sd(predicted) / sd(z_validation),
      n_learning_lines = length(learning_lines),
      n_validation_lines = length(validation_lines),
      n_learning_rows = sum(learning_rows),
      n_validation_rows = sum(validation_rows),
      n_eligible_conditions = length(gate),
      footprint_learning_only = TRUE,
      functional_axis_fold_specific_and_frozen = TRUE,
      baseline_scaling_learning_only = TRUE,
      interaction_scaling_learning_only = TRUE,
      selection_target = "molecular_z")
  }))
  selection_parts[[fold]] <- fold_metrics
  validation_store[[fold]] <- list(observed = z_validation, path = path,
    observation_id = score_table$observation_id[validation_rows],
    depmap_id = score_table$depmap_id[validation_rows],
    perturbation_id = score_table$perturbation_id[validation_rows])
}
selection <- rbindlist(selection_parts)
mean_loss <- selection[, .(mean_validation_mse = mean(validation_mse)),
  by = lambda]
selected_lambda <- mean_loss[which.min(mean_validation_mse), lambda]
selection[, `:=`(selected_lambda = selected_lambda,
  selected_penalty = lambda == selected_lambda)]
fwrite(selection, file.path(RESULT,
  "susceptibility_inner_model_selection.tsv"), sep = "\t")

selected_index <- match(selected_lambda, lambda_grid)
oof <- rbindlist(lapply(seq_len(N_INNER_FOLDS), function(fold)
  data.table(inner_fold = fold,
    observation_id = validation_store[[fold]]$observation_id,
    depmap_id = validation_store[[fold]]$depmap_id,
    perturbation_id = validation_store[[fold]]$perturbation_id,
    observed_z = validation_store[[fold]]$observed,
    predicted_z = validation_store[[fold]]$path[, selected_index])))

## Freeze one outer model using all 26 training lines and predict every held
## molecular response.  Failed-gate conditions receive exactly zero correction.
outer_label <- "outer_1"
outer_gate <- eligibility[fold == outer_label &
  molecular_context_variable == TRUE, perturbation_id]
outer_footprint <- make_footprint(training_lines, score_table, D, programs)
outer_baseline <- standardize_baseline(baseline_raw, training_lines,
  "outer_training")
outer_beta <- axis_vector(axis, "outer_training", programs)
outer_train <- score_table$depmap_id %chin% training_lines &
  score_table$perturbation_id %chin% outer_gate
outer_held <- score_table$depmap_id %chin% held_lines
X_train <- make_interactions(score_table$depmap_id[outer_train],
  score_table$perturbation_id[outer_train], outer_baseline$scores,
  outer_footprint, state_programs, programs)
R_train <- outer_footprint[score_table$perturbation_id[outer_train],
  programs, drop = FALSE]
z_train <- as.numeric((D[outer_train, programs, drop = FALSE] - R_train) %*%
  outer_beta)
outer_fit <- fit_ridge(X_train, z_train, selected_lambda)
outer_raw <- raw_ridge_coefficients(outer_fit)

X_held <- make_interactions(score_table$depmap_id[outer_held],
  score_table$perturbation_id[outer_held], outer_baseline$scores,
  outer_footprint, state_programs, programs)
R_held <- outer_footprint[score_table$perturbation_id[outer_held], programs,
  drop = FALSE]
C_held <- D[outer_held, programs, drop = FALSE] - R_held
observed_z <- as.numeric(C_held %*% outer_beta)
eligible_held <- score_table$perturbation_id[outer_held] %chin% outer_gate
ungated_predicted_z <- predict_ridge(outer_fit, X_held)
predicted_z <- fifelse(eligible_held, ungated_predicted_z, 0)
beta_norm2 <- sum(outer_beta ^ 2)
if (!is.finite(beta_norm2) || beta_norm2 <= 0)
  fail("The frozen functional molecular axis has zero norm.")
predicted_context <- outer(predicted_z / beta_norm2, outer_beta)
colnames(predicted_context) <- programs
D_hat <- R_held + predicted_context
projection_error <- max(abs(as.numeric((D_hat - R_held) %*% outer_beta) -
  predicted_z))
if (!is.finite(projection_error) || projection_error > 1e-10)
  fail("Reconstructed displacements do not recover predicted molecular z.")

feature_map <- tstrsplit(colnames(X_train), "__", fixed = TRUE)
coefficients <- data.table(baseline_program = feature_map[[1L]],
  response_program = feature_map[[2L]],
  interaction = colnames(X_train),
  coefficient = outer_raw$beta[colnames(X_train)],
  standardized_coefficient = as.numeric(outer_fit$beta_standardized)[
    match(colnames(X_train), outer_fit$feature_names[outer_fit$keep])],
  intercept = outer_raw$intercept, ridge_penalty = selected_lambda,
  n_training_lines = length(training_lines),
  n_training_rows = sum(outer_train))
coefficients[is.na(standardized_coefficient), standardized_coefficient := 0]
fwrite(coefficients, file.path(RESULT,
  "baseline_response_susceptibility_coefficients.tsv"), sep = "\t")

contribution_matrix <- sweep(X_held, 2L,
  outer_raw$beta[colnames(X_held)], "*")
contribution_matrix[!eligible_held, ] <- 0
applied_intercept <- fifelse(eligible_held, outer_raw$intercept, 0)
contribution_check <- max(abs(rowSums(contribution_matrix) +
  applied_intercept - predicted_z))
if (!is.finite(contribution_check) || contribution_check > 1e-10)
  fail("Interaction contributions do not reconstruct predicted molecular z.")
held_index <- which(outer_held)
contributions <- data.table(
  observation_id = rep(score_table$observation_id[held_index],
    each = ncol(X_held)),
  depmap_id = rep(score_table$depmap_id[held_index], each = ncol(X_held)),
  perturbation_id = rep(score_table$perturbation_id[held_index],
    each = ncol(X_held)),
  molecular_eligible = rep(eligible_held, each = ncol(X_held)),
  interaction = rep(colnames(X_held), times = nrow(X_held)),
  baseline_program = rep(feature_map[[1L]], times = nrow(X_held)),
  response_program = rep(feature_map[[2L]], times = nrow(X_held)),
  baseline_score_standardized = as.vector(t(do.call(cbind,
    lapply(state_programs, function(b) matrix(outer_baseline$scores[
      score_table$depmap_id[held_index], b], nrow = nrow(X_held),
      ncol = length(programs)))))),
  reference_program_displacement = as.vector(t(do.call(cbind,
    lapply(seq_along(state_programs), function(i) R_held)))),
  coefficient = rep(outer_raw$beta[colnames(X_held)], times = nrow(X_held)),
  contribution = as.vector(t(contribution_matrix)))
fwrite(contributions, file.path(RESULT,
  "baseline_response_susceptibility_contributions.tsv.gz"), sep = "\t",
  compress = "gzip")

activity <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
held <- data.table(observation_id = score_table$observation_id[held_index],
  depmap_id = score_table$depmap_id[held_index],
  cell_line = score_table$cell_line[held_index],
  perturbation_id = score_table$perturbation_id[held_index],
  drug = score_table$drug[held_index], dose = score_table$dose[held_index],
  dose_panel = score_table$dose_panel[held_index],
  molecular_eligible = eligible_held,
  molecular_gate_reason = eligibility[fold == outer_label,
    gate_failure_reason[match(score_table$perturbation_id[held_index],
      perturbation_id)]],
  activity_supported = activity$activity_supported[match(
    score_table$perturbation_id[held_index], activity$perturbation_id)],
  observed_z = observed_z, predicted_z = predicted_z,
  ungated_predicted_z = ungated_predicted_z,
  applied_intercept = applied_intercept,
  frozen_axis_norm_squared = beta_norm2)
for (program in programs) {
  set(held, j = paste0("reference_", program), value = R_held[, program])
  set(held, j = paste0("observed_context_", program),
    value = C_held[, program])
  set(held, j = paste0("observed_displacement_", program),
    value = D[outer_held, program])
  set(held, j = paste0("reconstructed_predicted_displacement_", program),
    value = D_hat[, program])
}
for (feature in colnames(X_held)) set(held,
  j = paste0("contribution_", feature),
  value = contribution_matrix[, feature])
fwrite(held, file.path(RESULT,
  "held_function_relevant_molecular_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")

reconstructed <- data.table(
  observation_id = rep(held$observation_id, each = length(programs)),
  depmap_id = rep(held$depmap_id, each = length(programs)),
  cell_line = rep(held$cell_line, each = length(programs)),
  perturbation_id = rep(held$perturbation_id, each = length(programs)),
  drug = rep(held$drug, each = length(programs)),
  dose = rep(held$dose, each = length(programs)),
  molecular_eligible = rep(held$molecular_eligible, each = length(programs)),
  program = rep(programs, times = nrow(held)),
  observed_displacement = as.vector(t(D[outer_held, programs, drop = FALSE])),
  reference_footprint = as.vector(t(R_held)),
  observed_context_displacement = as.vector(t(C_held)),
  predicted_functional_context_displacement = as.vector(t(predicted_context)),
  reconstructed_predicted_displacement = as.vector(t(D_hat)),
  functional_axis_coefficient = rep(outer_beta, times = nrow(held)))
fwrite(reconstructed, file.path(RESULT,
  "held_reconstructed_program_displacements.tsv.gz"), sep = "\t",
  compress = "gzip")

per_line <- held[, {
  eligible <- molecular_eligible == TRUE
  list(n_rows = .N, n_eligible_rows = sum(eligible),
    pearson_all = safe_cor(observed_z, predicted_z),
    pearson_eligible = safe_cor(observed_z[eligible], predicted_z[eligible]),
    observed_sd_eligible = sd(observed_z[eligible]),
    predicted_sd_eligible = sd(predicted_z[eligible]),
    predicted_to_observed_sd_ratio_eligible =
      sd(predicted_z[eligible]) / sd(observed_z[eligible]),
    rmse_eligible = sqrt(mean((observed_z[eligible] -
      predicted_z[eligible]) ^ 2)))
}, by = .(depmap_id, cell_line)]
fwrite(per_line, file.path(RESULT, "susceptibility_per_line.tsv"), sep = "\t")

per_perturbation <- held[, {
  list(n_rows = .N, molecular_eligible = molecular_eligible[1L],
    activity_supported = activity_supported[1L],
    pearson = safe_cor(observed_z, predicted_z),
    observed_sd = sd(observed_z), predicted_sd = sd(predicted_z),
    predicted_to_observed_sd_ratio = sd(predicted_z) / sd(observed_z),
    rmse = sqrt(mean((observed_z - predicted_z) ^ 2)))
}, by = .(perturbation_id, drug, dose)]
fwrite(per_perturbation, file.path(RESULT,
  "susceptibility_per_perturbation.tsv"), sep = "\t")

baseline_contribution <- contributions[, .(
  baseline_program_contribution = sum(contribution)),
  by = .(observation_id, baseline_program)]
baseline_contribution <- merge(baseline_contribution,
  held[, .(observation_id, molecular_eligible, observed_z)],
  by = "observation_id", all.x = TRUE)
per_baseline <- baseline_contribution[molecular_eligible == TRUE, .(
  n_rows = .N,
  contribution_vs_observed_z_pearson = safe_cor(observed_z,
    baseline_program_contribution),
  mean_absolute_contribution = mean(abs(baseline_program_contribution)),
  signed_mean_contribution = mean(baseline_program_contribution),
  contribution_sd = sd(baseline_program_contribution)),
  by = baseline_program]
per_baseline <- merge(per_baseline,
  coefficients[, .(coefficient_l1 = sum(abs(coefficient)),
    coefficient_l2 = sqrt(sum(coefficient ^ 2))), by = baseline_program],
  by = "baseline_program", all = TRUE)
fwrite(per_baseline, file.path(RESULT,
  "susceptibility_per_baseline_program.tsv"), sep = "\t")

summary <- rbindlist(list(
  cbind(data.table(cohort = "inner_validation_molecular_eligible"),
    molecular_metrics(oof$observed_z, oof$predicted_z)),
  cbind(data.table(cohort = "held_all_molecular_rows"),
    molecular_metrics(held$observed_z, held$predicted_z)),
  cbind(data.table(cohort = "held_molecular_eligible"),
    molecular_metrics(held[molecular_eligible == TRUE, observed_z],
      held[molecular_eligible == TRUE, predicted_z])),
  cbind(data.table(cohort = "held_activity_supported"),
    molecular_metrics(held[activity_supported == TRUE, observed_z],
      held[activity_supported == TRUE, predicted_z])),
  cbind(data.table(cohort = "held_eligible_and_activity_supported"),
    molecular_metrics(held[molecular_eligible == TRUE &
      activity_supported == TRUE, observed_z],
      held[molecular_eligible == TRUE &
        activity_supported == TRUE, predicted_z]))), fill = TRUE)
summary[, `:=`(selected_lambda = selected_lambda,
  n_outer_eligible_conditions = length(outer_gate),
  susceptibility_predictors = "S2-S8 x frozen P1-P10 (stable nine only)",
  target = "function-relevant molecular displacement z",
  functional_outcomes_loaded = FALSE)]
fwrite(summary, file.path(RESULT, "susceptibility_summary.tsv"), sep = "\t")

saveRDS(list(label = "continuous_baseline_by_reference_susceptibility",
    baseline_programs = state_programs, response_programs = programs,
    selected_lambda = selected_lambda, fit = outer_fit,
    raw_coefficients = outer_raw, baseline_center = outer_baseline$center,
    baseline_scale = outer_baseline$scale, functional_axis = outer_beta,
    functional_axis_norm_squared = beta_norm2,
    molecularly_eligible_conditions = outer_gate),
  file.path(ARTIFACT, "continuous_susceptibility_model.rds"))

## Five owned, one-page molecular diagnostic plots.
coefficient_wide <- dcast(coefficients,
  baseline_program ~ response_program, value.var = "coefficient")
coefficient_matrix <- as.matrix(coefficient_wide[, ..programs])
rownames(coefficient_matrix) <- coefficient_wide$baseline_program
coefficient_matrix <- coefficient_matrix[state_programs, programs,
  drop = FALSE]
coefficient_limit <- max(abs(coefficient_matrix))
coefficient_heatmap <- draw_heatmap(coefficient_matrix,
  scale_rows = FALSE, value_range = c(-coefficient_limit, 0,
    coefficient_limit), value_colors = c("#2166AC", "white", "#B2182B"),
  show_row_names = TRUE, show_column_names = TRUE,
  cluster_rows = FALSE, cluster_columns = FALSE,
  title = "Continuous baseline-response susceptibility coefficients")
save_heatmap("14a_susceptibility_coefficient_heatmap.pdf",
  coefficient_heatmap, 6.2, 4.6)

plot_held <- held[molecular_eligible == TRUE]
p_z <- ggplot(plot_held, aes(observed_z, predicted_z)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.35, size = 0.8, color = "#2166AC") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.65, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("r=%s; SD ratio=%s\nn=%s; %d conditions",
      fmt(safe_cor(plot_held$observed_z, plot_held$predicted_z)),
      fmt(sd(plot_held$predicted_z) / sd(plot_held$observed_z)),
      format(nrow(plot_held), big.mark = ","),
      uniqueN(plot_held$perturbation_id)), size = 2.35) +
  labs(x = "Observed molecular z", y = "Predicted molecular z",
    title = "Held-out function-relevant molecular displacement") +
  square_scatter_theme()
save_gg("14b_observed_vs_predicted_held_molecular_z.pdf", p_z, 5.0, 4.8)

p_line <- ggplot(per_line, aes(reorder(cell_line, pearson_eligible),
    pearson_eligible)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_segment(aes(xend = cell_line, y = 0, yend = pearson_eligible),
    linewidth = 0.6, color = "grey65") +
  geom_point(size = 2.2, color = "#2166AC") +
  coord_flip() +
  labs(x = NULL, y = "Held molecular-z Pearson r",
    title = "Continuous susceptibility performance by held-out line") +
  theme_stage14()
save_gg("14c_per_line_molecular_z_performance.pdf", p_line, 5.8, 4.5)

plot_perturbation <- per_perturbation[molecular_eligible == TRUE &
  is.finite(pearson)][order(pearson)]
plot_perturbation[, rank := seq_len(.N)]
p_perturbation <- ggplot(plot_perturbation, aes(rank, pearson)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_point(alpha = 0.75, size = 1.15, color = "#6A51A3") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.5, color = "black") +
  labs(x = "Eligible perturbation rank", y = "Held molecular-z Pearson r",
    title = "Continuous susceptibility performance by perturbation") +
  theme_stage14()
save_gg("14d_per_perturbation_molecular_z_performance.pdf",
  p_perturbation, 6.0, 4.3)

plot_amplitude <- per_perturbation[molecular_eligible == TRUE &
  is.finite(observed_sd) & observed_sd > 0 & is.finite(predicted_sd)]
p_amplitude <- ggplot(plot_amplitude, aes(observed_sd, predicted_sd)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3,
    linewidth = 0.45, color = "grey40") +
  geom_point(alpha = 0.75, size = 1.2, color = "#B2182B") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55, color = "black") +
  labs(x = "Observed molecular-z SD", y = "Predicted molecular-z SD",
    title = "Predicted-to-observed molecular-z amplitude") +
  square_scatter_theme()
save_gg("14e_predicted_to_observed_molecular_z_amplitude.pdf",
  p_amplitude, 4.8, 4.6)

figures <- data.table(pdf = c(
  "14a_susceptibility_coefficient_heatmap.pdf",
  "14b_observed_vs_predicted_held_molecular_z.pdf",
  "14c_per_line_molecular_z_performance.pdf",
  "14d_per_perturbation_molecular_z_performance.pdf",
  "14e_predicted_to_observed_molecular_z_amplitude.pdf"),
  producing_script = "scripts/14_fit_continuous_baseline_susceptibility.R")
figures[, exists := file.exists(file.path(RESULT, pdf))]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")

checks <- data.table(check = c(
  "accepted_26_training_7_held_split",
  "excluded_lines_absent", "exactly_S2_to_S8",
  "accepted_nine_response_programs", "sixty_three_interactions_only",
  "fold_specific_frozen_functional_axis", "molecular_training_only_gate",
  "no_functional_outcome_loaded", "inner_footprints_recomputed",
  "inner_baseline_and_interaction_scaling_learning_only",
  "held_gate_failure_means_zero_correction",
  "program_reconstruction_projection_identity",
  "interaction_contribution_reconstruction", "all_required_pdfs_created"),
  value = c(length(training_lines) == 26L && length(held_lines) == 7L,
    !any(EXCLUDED_IDS %chin% all_lines),
    identical(state_programs, paste0("S", 2:8)), length(programs) == 9L,
    nrow(coefficients) == 63L,
    setequal(unique(axis$fold), c(paste0("inner_", 1:5),
      "outer_training")),
    all(eligibility$selection_training_derived) &&
      !any(eligibility$performance_used_for_selection) &&
      !any(eligibility$validation_or_held_values_used_for_selection),
    TRUE, all(selection$footprint_learning_only),
    all(selection$baseline_scaling_learning_only) &&
      all(selection$interaction_scaling_learning_only),
    all(held[molecular_eligible == FALSE, predicted_z] == 0),
    projection_error <= 1e-10, contribution_check <= 1e-10,
    all(figures$exists)),
  detail = c("fixed outer split", "H596 and H2122 excluded",
    paste(state_programs, collapse = ","), paste(programs, collapse = ","),
    "7 baseline x 9 reference-response interactions; no main or identity features",
    "Stage 05 coefficients at the accepted shared-decoder penalty",
    "activity and dispersion use fold-training Tahoe molecular data only",
    "Stage 14 input manifest contains no functional-outcome table",
    "reference footprints use inner-learning lines",
    "standardization is fitted on inner-learning rows only",
    "ineligible conditions retain the arithmetic reference footprint",
    sprintf("maximum absolute error %.3g", projection_error),
    sprintf("maximum absolute error %.3g", contribution_check),
    "five one-page-owned molecular diagnostics"))
fwrite(checks, file.path(RESULT, "susceptibility_required_checks.tsv"),
  sep = "\t")
if (!all(checks$value)) fail("A required Stage 14 check failed.")
write_artifact_manifest(ARTIFACT,
  "scripts/14_fit_continuous_baseline_susceptibility.R")
writeLines("complete", file.path(RESULT, "stage14_complete.flag"))
message(sprintf(paste0("Stage 14 complete: lambda %.6g; held eligible ",
  "molecular z r=%s; SD ratio=%s."), selected_lambda,
  fmt(summary[cohort == "held_molecular_eligible", pearson]),
  fmt(summary[cohort == "held_molecular_eligible",
    predicted_to_observed_sd_ratio])))
