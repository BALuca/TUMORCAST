#!/usr/bin/env Rscript

## Fit the accepted shared nine-program ridge decoder, then calculate the
## measured function-relevant molecular context score.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("05_program_function_bridge")
ARTIFACT <- stage_artifact("05_program_function_bridge")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

plot_stage05 <- function() {
  theme_pilot <- function(base_size = 8.5) theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
      strip.background = element_rect(fill = "grey95", color = "grey75"),
      plot.title = element_text(face = "bold"))
  square_scatter_theme <- function(base_size = 8.2) theme_pilot(base_size) +
    theme(aspect.ratio = 1)
  save_gg <- function(filename, plot, width, height) ggsave(
    file.path(RESULT, filename), plot, width = width, height = height,
    device = cairo_pdf, bg = "white", limitsize = FALSE)
  fmt <- function(x, digits = 3L) ifelse(is.finite(x),
    formatC(x, digits = digits, format = "f"), "NA")

  decoder_path <- fread(file.path(RESULT,
    "function_decoder_inner_selection.tsv"))
  decoder_mean <- decoder_path[, .(
    mean_validation_mse = mean(validation_mse),
    standard_error = sd(validation_mse) / sqrt(.N)), by = lambda]
  selected_lambda <- decoder_path$selected_lambda[1L]
  p_selection <- ggplot(decoder_mean, aes(log10(lambda),
      mean_validation_mse)) +
    geom_ribbon(aes(ymin = mean_validation_mse - standard_error,
      ymax = mean_validation_mse + standard_error), alpha = 0.12,
      fill = "#3B7A57", color = NA) +
    geom_line(linewidth = 0.65, color = "#3B7A57") +
    geom_point(size = 0.9, color = "#3B7A57") +
    geom_point(data = decoder_mean[lambda == selected_lambda], size = 3,
      shape = 21, fill = "gold", color = "black") +
    labs(x = expression(log[10](lambda)),
      y = "Mean validation MSE",
      title = "Shared program-to-function decoder penalty path",
      subtitle = "Gold point: selected penalty") + theme_pilot()
  save_gg("05a_function_decoder_model_selection.pdf", p_selection, 5.8, 4.2)

  stability <- fread(file.path(stage_result("02_response_programs"),
    "fold_program_stability.tsv"))
  programs <- stability[outer_fold == OUTER_FOLD & stable == TRUE,
    unique(program)]
  programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
  score_long <- fread(file.path(stage_artifact("02_response_programs"),
    "fold_program_scores.tsv.gz"))[role == "held_out" &
      program %chin% programs]
  footprint <- fread(file.path(RESULT,
    "reference_perturbation_footprints.tsv.gz"))
  footprint_long <- melt(footprint, id.vars = "perturbation_id",
    measure.vars = programs, variable.name = "program",
    value.name = "predicted_displacement")
  reference_long <- merge(score_long[, .(observation_id, depmap_id,
    perturbation_id, program, measured_displacement = score)],
    footprint_long, by = c("perturbation_id", "program"), all.x = TRUE,
    sort = FALSE)
  reference_plot <- ggplot(reference_long,
    aes(measured_displacement, predicted_displacement)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = 3,
      linewidth = 0.35, color = "grey45") +
    geom_point(alpha = 0.12, size = 0.35, color = "#4D4D4D") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 0.55, color = "black") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
      label = sprintf("r=%s\n%d lines; %d perturbations; n=%s",
        fmt(safe_cor(reference_long$measured_displacement,
          reference_long$predicted_displacement)),
        uniqueN(reference_long$depmap_id),
        uniqueN(reference_long$perturbation_id),
        format(nrow(reference_long), big.mark = ",")), size = 2.35) +
    labs(x = "Measured total displacement",
      y = "Predicted total displacement",
      title = "Reference perturbation footprint") + square_scatter_theme()
  save_gg("05b_molecular_reference_prediction.pdf", reference_plot, 5.0, 4.8)

  contribution <- fread(file.path(RESULT,
    "program_contribution_summary.tsv"))
  p_contribution <- ggplot(contribution,
    aes(program, mean_absolute_contribution)) +
    geom_col(width = 0.68, fill = "#6A51A3") +
    labs(x = NULL, y = "Mean absolute functional contribution",
      title = "Response-program contribution to function") + theme_pilot()
  save_gg("05c_program_function_contributions.pdf", p_contribution, 5.6, 4.2)

  function_rows <- fread(file.path(ARTIFACT,
    "measured_function_predictions.tsv.gz"))
  plot_current <- function_rows[evaluation_included == TRUE &
    activity_supported == TRUE]
  plot_current[, observed_function_context := observed_function -
    training_reference_function]
  p_measured <- ggplot(plot_current,
    aes(measured_z, observed_function_context)) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 0.7) +
    labs(title = sprintf("Nine programs: measured functional score (r = %s)",
      fmt(safe_cor(plot_current$measured_z,
        plot_current$observed_function_context))),
      x = "Measured molecular functional score z",
      y = "Observed PRISM context") + square_scatter_theme()
  save_gg("05d_measured_functional_score_current_programs.pdf",
    p_measured, 3.8, 3.8)

  model_columns <- c(perturbation_average = "mean_prediction",
    complete_measured = "measured_prediction")
  model_labels <- c(perturbation_average = "Perturbation average",
    complete_measured = "Complete measured response")
  model_colors <- c(perturbation_average = "#4D4D4D",
    complete_measured = "#3B7A57")
  metric_labels <- c(pooled = "Pooled", within_perturbation =
    "Within perturbation", strict_cell_line_by_perturbation =
    "Strict cell line by perturbation")
  for (model in names(model_columns)) for (metric in names(metric_labels)) {
    current <- function_residual_table(copy(plot_current), metric,
      model_columns[[model]])
    label <- sprintf("r=%s\nobs SD=%s; pred SD=%s\nn=%s",
      fmt(safe_cor(current$observed_residual, current$predicted_residual)),
      fmt(sd(current$observed_residual)),
      fmt(sd(current$predicted_residual)),
      format(nrow(current), big.mark = ","))
    plot <- ggplot(current, aes(observed_residual, predicted_residual)) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
      geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
      geom_point(alpha = 0.22, size = 0.45, color = model_colors[[model]]) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
        linewidth = 0.55, color = "black") +
      annotate("text", x = -Inf, y = Inf, label = label,
        hjust = -0.05, vjust = 1.1, size = 2.15) +
      labs(x = "Observed function", y = "Predicted function",
        title = paste(model_labels[[model]], metric_labels[[metric]],
          sep = " — ")) + square_scatter_theme(7.5)
    save_gg(sprintf("05e_%s_%s.pdf", model, metric), plot, 4.4, 4.2)
  }
}

if ("--plots-only" %in% commandArgs(trailingOnly = TRUE)) {
  plot_stage05()
  message("Stage 05 plots refreshed from existing bridge artifacts.")
  quit(save = "no", status = 0L)
}

stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- stability[outer_fold == OUTER_FOLD & stable == TRUE,
  unique(program)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS)
  fail("The reproduced response bank is not nine stable programs.")
scores_long <- fread(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"))[program %chin% programs]
scores <- dcast(scores_long, observation_id + depmap_id + cell_line +
  perturbation_id + drug + dose + dose_panel + role ~ program,
  value.var = "score")
setkey(scores, observation_id)
D <- as.matrix(scores[, ..programs]); rownames(D) <- scores$observation_id
if (any(!is.finite(D))) fail("Reproduced response-program scores are invalid.")

function_rows <- fread(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"))
function_rows[, `:=`(observation_id = paste(depmap_id, perturbation_id,
  sep = "__"), drug = tahoe_drug, dose = tahoe_dose_uM)]
position <- match(function_rows$observation_id, scores$observation_id)
chain <- function_rows$chain_observed == 1L & !is.na(position)
folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- folds[outer_fold != OUTER_FOLD, depmap_id]
held_lines <- folds[outer_fold == OUTER_FOLD, depmap_id]
inner <- make_inner_assignment(training_lines)

## Preserve the accepted five-fold inner selection for the shared decoder.
loss <- matrix(NA_real_, N_INNER_FOLDS, length(DECODER_LAMBDA_GRID))
correlation <- loss
selection_parts <- vector("list", N_INNER_FOLDS)
for (fold in seq_len(N_INNER_FOLDS)) {
  validation_lines <- names(inner)[inner == fold]
  learn_lines <- names(inner)[inner != fold]
  train <- function_rows$depmap_id %chin% learn_lines & chain
  valid <- function_rows$depmap_id %chin% validation_lines & chain
  path <- ridge_path(D[position[train], , drop = FALSE],
    function_rows$prism_sensitization[train],
    D[position[valid], , drop = FALSE], DECODER_LAMBDA_GRID)
  loss[fold, ] <- colMeans(sweep(path, 1L,
    function_rows$prism_sensitization[valid], "-") ^ 2)
  correlation[fold, ] <- apply(path, 2L, function(predicted)
    safe_cor(function_rows$prism_sensitization[valid], predicted))
  selection_parts[[fold]] <- data.table(inner_fold = fold,
    lambda = DECODER_LAMBDA_GRID, validation_mse = loss[fold, ],
    validation_pearson = correlation[fold, ],
    n_training_lines = length(learn_lines),
    n_validation_lines = length(validation_lines),
    n_training_rows = sum(train), n_validation_rows = sum(valid))
}
best <- which.min(colMeans(loss)); selected_lambda <- DECODER_LAMBDA_GRID[best]
decoder_train <- function_rows$depmap_id %chin% training_lines & chain
decoder_fit <- fit_ridge(D[position[decoder_train], , drop = FALSE],
  function_rows$prism_sensitization[decoder_train], selected_lambda)
decoder <- raw_ridge_coefficients(decoder_fit)

## Export the already accepted functional direction for outer training and for
## each existing inner-learning split.  These fold-specific coefficient
## vectors are frozen inputs for later molecular-only susceptibility fitting;
## this export does not change decoder selection or any prediction below.
axis_parts <- lapply(seq_len(N_INNER_FOLDS), function(fold) {
  learn_lines <- names(inner)[inner != fold]
  learn <- function_rows$depmap_id %chin% learn_lines & chain
  fit <- fit_ridge(D[position[learn], , drop = FALSE],
    function_rows$prism_sensitization[learn], selected_lambda)
  raw <- raw_ridge_coefficients(fit)
  data.table(fold = paste0("inner_", fold), program = programs,
    coefficient = raw$beta[programs], intercept = raw$intercept,
    ridge_penalty = selected_lambda,
    n_training_lines = length(learn_lines))
})
axis_parts[[length(axis_parts) + 1L]] <- data.table(
  fold = "outer_training", program = programs,
  coefficient = decoder$beta[programs], intercept = decoder$intercept,
  ridge_penalty = selected_lambda,
  n_training_lines = length(training_lines))
functional_axis <- rbindlist(axis_parts)
fwrite(functional_axis, file.path(RESULT,
  "functional_molecular_axis_coefficients.tsv"), sep = "\t")
selection <- rbindlist(selection_parts)
selection[, `:=`(selected = lambda == selected_lambda,
  selected_lambda = selected_lambda)]
selection_summary <- data.table(selected_lambda = selected_lambda,
  mean_inner_mse = colMeans(loss)[best],
  mean_inner_pearson = mean(correlation[, best]),
  standard_error = sd(loss[, best]) / sqrt(N_INNER_FOLDS),
  n_training_chain_rows = sum(decoder_train),
  n_inner_folds = N_INNER_FOLDS, predictors = paste(programs,
    collapse = ","), global_intercept = TRUE, drug_features = FALSE,
  baseline_features = FALSE, direct_transcriptome_features = FALSE)
coefficients <- data.table(term = c("(Intercept)", programs),
  coefficient = c(decoder$intercept, decoder$beta[programs]),
  penalty_lambda = selected_lambda,
  permitted_design = "nine measured response programs + one shared intercept")
fwrite(selection, file.path(RESULT, "function_decoder_inner_selection.tsv"),
  sep = "\t")
fwrite(selection_summary,
  file.path(RESULT, "function_decoder_selection.tsv"), sep = "\t")
fwrite(coefficients, file.path(RESULT, "function_decoder_coefficients.tsv"),
  sep = "\t")
saveRDS(decoder_fit, file.path(ARTIFACT, "function_decoder_model.rds"))

## Arithmetic training-line reference footprints for all 1,137 conditions.
train_score <- scores$depmap_id %chin% training_lines
footprint_dt <- data.table(perturbation_id = scores$perturbation_id[train_score],
  D[train_score, , drop = FALSE])[, lapply(.SD, mean),
    by = perturbation_id, .SDcols = programs]
footprint <- as.matrix(footprint_dt[, ..programs])
rownames(footprint) <- footprint_dt$perturbation_id
if (nrow(footprint) != EXPECTED_DISCOVERY_CONDITIONS)
  fail("Reference footprints do not cover all Tahoe conditions.")
fwrite(footprint_dt, file.path(RESULT, "reference_perturbation_footprints.tsv.gz"),
  sep = "\t", compress = "gzip")

footprint_position <- match(scores$perturbation_id, rownames(footprint))
delta <- D - footprint[footprint_position, programs, drop = FALSE]
measured_z <- as.numeric(delta %*% decoder$beta[programs])
measured <- data.table(observation_id = scores$observation_id,
  depmap_id = scores$depmap_id, cell_line = scores$cell_line,
  perturbation_id = scores$perturbation_id, drug = scores$drug,
  dose = scores$dose, dose_panel = scores$dose_panel, role = scores$role,
  measured_z = measured_z)
contribution <- data.table(
  observation_id = rep(scores$observation_id, each = length(programs)),
  depmap_id = rep(scores$depmap_id, each = length(programs)),
  perturbation_id = rep(scores$perturbation_id, each = length(programs)),
  drug = rep(scores$drug, each = length(programs)),
  program = rep(programs, times = nrow(scores)),
  context_program_displacement = as.vector(t(delta)),
  decoder_coefficient = rep(decoder$beta[programs], times = nrow(scores)))
contribution[, function_relevant_contribution :=
  context_program_displacement * decoder_coefficient]
if (max(abs(contribution[, sum(function_relevant_contribution),
    by = observation_id]$V1 - measured_z)) > 1e-9)
  fail("Program contributions do not reconstruct the measured target.")
fwrite(measured, file.path(ARTIFACT,
  "measured_function_relevant_scores.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(contribution, file.path(ARTIFACT,
  "measured_function_relevant_contributions.tsv.gz"), sep = "\t",
  compress = "gzip")
fwrite(contribution[depmap_id %chin% training_lines,
  .(mean_absolute_contribution = mean(abs(function_relevant_contribution)),
    signed_mean_contribution = mean(function_relevant_contribution)),
  by = program], file.path(RESULT, "program_contribution_summary.tsv"),
  sep = "\t")

## Held functional predictions from the mean footprint and complete response.
held_function <- function_rows[outer_fold == OUTER_FOLD]
held_function_position <- match(held_function$observation_id,
  scores$observation_id)
held_footprint_position <- match(held_function$perturbation_id,
  rownames(footprint))
mean_D <- footprint[held_footprint_position, programs, drop = FALSE]
mean_prediction <- decoder$intercept +
  as.numeric(mean_D %*% decoder$beta[programs])
measured_prediction <- rep(NA_real_, nrow(held_function))
available <- !is.na(held_function_position)
measured_prediction[available] <- decoder$intercept + as.numeric(
  D[held_function_position[available], programs, drop = FALSE] %*%
    decoder$beta[programs])
measured_z_function <- measured_prediction - mean_prediction
training_function_mean <- function_rows[outer_fold != OUTER_FOLD &
  chain_observed == 1L, .(training_reference_function = mean(
    prism_sensitization)), by = perturbation_id]
panel <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
function_predictions <- data.table(
  observation_id = held_function$observation_id,
  depmap_id = held_function$depmap_id, cell_line = held_function$cell_line,
  perturbation_id = held_function$perturbation_id, drug = held_function$drug,
  dose = held_function$dose, observed_function = held_function$prism_sensitization,
  training_reference_function = training_function_mean$
    training_reference_function[match(held_function$perturbation_id,
      training_function_mean$perturbation_id)],
  evaluation_included = held_function$chain_observed == 1L,
  activity_supported = panel$activity_supported[match(
    held_function$perturbation_id, panel$perturbation_id)],
  mean_prediction = mean_prediction,
  measured_prediction = measured_prediction,
  measured_z = measured_z_function)
fwrite(function_predictions, file.path(ARTIFACT,
  "measured_function_predictions.tsv.gz"), sep = "\t", compress = "gzip")

## Molecular-reference metrics use the accepted program-long definitions.
held_score <- scores$depmap_id %chin% held_lines
reference_long <- data.table(
  observation_id = rep(scores$observation_id[held_score], each = length(programs)),
  depmap_id = rep(scores$depmap_id[held_score], each = length(programs)),
  perturbation_id = rep(scores$perturbation_id[held_score], each = length(programs)),
  program = rep(programs, times = sum(held_score)),
  measured_displacement = as.vector(t(D[held_score, , drop = FALSE])),
  predicted_displacement = as.vector(t(footprint[
    footprint_position[held_score], programs, drop = FALSE])))
reference_long[, `:=`(
  observed_within = measured_displacement - mean(measured_displacement),
  predicted_within = predicted_displacement - mean(predicted_displacement)),
  by = perturbation_id]
molecular_reference_performance <- data.table(
  n_program_observations = nrow(reference_long),
  n_response_rows = uniqueN(reference_long$observation_id),
  n_cell_lines = uniqueN(reference_long$depmap_id),
  n_perturbations = uniqueN(reference_long$perturbation_id),
  held_out_pearson = safe_cor(reference_long$measured_displacement,
    reference_long$predicted_displacement),
  within_perturbation_pearson = safe_cor(reference_long$observed_within,
    reference_long$predicted_within), correlation_method = "pearson")
fwrite(molecular_reference_performance,
  file.path(RESULT, "molecular_reference_performance.tsv"), sep = "\t")

## Shared bridge and measured-target performance on complete and active panels.
performance_parts <- list()
for (cohort in c("complete_matched_panel", "activity_supported")) {
  current <- function_predictions[evaluation_included == TRUE]
  if (cohort == "activity_supported") current <- current[
    activity_supported == TRUE]
  for (model in c("perturbation_average", "complete_measured")) {
    column <- if (model == "perturbation_average") "mean_prediction" else
      "measured_prediction"
    for (metric in c("pooled", "within_perturbation",
        "strict_cell_line_by_perturbation")) {
      residual <- function_residual_table(current, metric, column)
      w <- balanced_weights(residual$drug)
      performance_parts[[length(performance_parts) + 1L]] <- data.table(
        cohort = cohort, prediction_model = model, metric = metric,
        n_rows = nrow(residual), n_conditions = uniqueN(residual$perturbation_id),
        n_drugs = uniqueN(residual$drug), n_lines = uniqueN(residual$depmap_id),
        pearson = safe_cor(residual$observed_residual,
          residual$predicted_residual),
        drug_balanced_pearson = weighted_cor(residual$observed_residual,
          residual$predicted_residual, w),
        mse = mean((residual$observed_residual -
          residual$predicted_residual) ^ 2), correlation_method = "pearson")
    }
  }
}
bridge_performance <- rbindlist(performance_parts)
fwrite(bridge_performance,
  file.path(RESULT, "program_function_performance.tsv"), sep = "\t")
active_context <- function_predictions[evaluation_included == TRUE &
  activity_supported == TRUE]
active_context[, observed_function_context := observed_function -
  training_reference_function]
measured_performance <- endpoint_metrics(active_context,
  "observed_function_context", "measured_z")
fwrite(measured_performance,
  file.path(RESULT, "measured_function_relevant_performance.tsv"), sep = "\t")

write_artifact_manifest(ARTIFACT,
  "scripts/05_fit_program_function_bridge.R")
writeLines("complete", file.path(RESULT, "stage05_complete.flag"))
plot_stage05()
message(sprintf("Stage 05 complete: decoder lambda %.6g; measured score r=%.3f.",
  selected_lambda, measured_performance$pearson))
