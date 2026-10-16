#!/usr/bin/env Rscript

## Compose the Stage 14 reconstructed molecular displacement with the already
## frozen shared response-to-function decoder.  Functional outcomes enter only
## here for held-out evaluation; no model is fitted or recalibrated.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("15_continuous_susceptibility_chain")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)

theme_stage15 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage15(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(name, plot, width, height) ggsave(file.path(RESULT, name),
  plot, width = width, height = height, device = cairo_pdf, bg = "white",
  limitsize = FALSE)
fmt <- function(x, digits = 3L) ifelse(is.finite(x),
  formatC(x, digits = digits, format = "f"), "NA")

evaluate_function <- function(table, model, cohort) {
  current <- table[prediction_model == model & evaluation_included == TRUE]
  if (cohort == "activity_supported")
    current <- current[activity_supported == TRUE]
  rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric) {
    residual <- function_residual_table(current, metric,
      "predicted_function")
    w <- balanced_weights(residual$drug)
    error <- residual$observed_residual - residual$predicted_residual
    data.table(cohort = cohort, prediction_model = model, metric = metric,
      n_rows = nrow(residual), n_conditions = uniqueN(residual$perturbation_id),
      n_drugs = uniqueN(residual$drug), n_lines = uniqueN(residual$depmap_id),
      pearson = safe_cor(residual$observed_residual,
        residual$predicted_residual),
      drug_balanced_pearson = weighted_cor(residual$observed_residual,
        residual$predicted_residual, w),
      mse = mean(error ^ 2), rmse = sqrt(mean(error ^ 2)),
      observed_sd = sd(residual$observed_residual),
      predicted_sd = sd(residual$predicted_residual),
      predicted_to_observed_sd_ratio = sd(residual$predicted_residual) /
        sd(residual$observed_residual), correlation_method = "pearson")
  }))
}

molecular_metrics <- function(table, prediction, cohort) {
  current <- copy(table)
  if (cohort == "held_activity_supported")
    current <- current[activity_supported == TRUE]
  if (cohort == "held_molecular_eligible")
    current <- current[molecular_eligible == TRUE]
  if (cohort == "held_eligible_and_activity_supported")
    current <- current[molecular_eligible == TRUE &
      activity_supported == TRUE]
  if (cohort == "held_matched_activity_panel")
    current <- current[matched_activity_panel == TRUE]
  keep <- is.finite(current$observed_z) & is.finite(current[[prediction]])
  current <- current[keep]
  observed <- current$observed_z; predicted <- current[[prediction]]
  data.table(endpoint = "molecular_z", cohort = cohort,
    metric = "pooled_pearson", prediction_model = prediction,
    n_rows = nrow(current), n_conditions = uniqueN(current$perturbation_id),
    n_drugs = uniqueN(current$drug), n_lines = uniqueN(current$depmap_id),
    pearson = safe_cor(observed, predicted),
    observed_sd = sd(observed), predicted_sd = sd(predicted),
    predicted_to_observed_sd_ratio = sd(predicted) / sd(observed),
    mse = mean((observed - predicted) ^ 2),
    rmse = sqrt(mean((observed - predicted) ^ 2, na.rm = TRUE)),
    correlation_method = "pearson")
}

required <- c(
  file.path(stage_result("14_continuous_baseline_susceptibility"),
    "stage14_complete.flag"),
  file.path(stage_result("14_continuous_baseline_susceptibility"),
    "held_function_relevant_molecular_predictions.tsv.gz"),
  file.path(stage_result("14_continuous_baseline_susceptibility"),
    "held_reconstructed_program_displacements.tsv.gz"),
  file.path(stage_result("05_program_function_bridge"),
    "function_decoder_coefficients.tsv"),
  file.path(stage_artifact("05_program_function_bridge"),
    "measured_function_predictions.tsv.gz"),
  file.path(stage_artifact("06_baseline_bridge"),
    "baseline_function_predictions.tsv.gz"),
  file.path(stage_artifact("06_baseline_bridge"),
    "baseline_molecular_predictions.tsv.gz"),
  file.path(stage_artifact("11_state_modulated_response"),
    "state_modulated_molecular_predictions.tsv.gz"),
  file.path(stage_artifact("13_separated_chain"),
    "separated_chain_held_predictions.tsv.gz"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 15 input(s): %s",
  paste(missing, collapse = ", "))

reconstructed <- fread(file.path(stage_result(
  "14_continuous_baseline_susceptibility"),
  "held_reconstructed_program_displacements.tsv.gz"))
programs <- paste0("P", sort(as.integer(sub("P", "",
  unique(reconstructed$program)))))
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS)
  fail("Stage 14 reconstructed displacement is not nine-dimensional.")
Dhat_wide <- dcast(reconstructed, observation_id ~ program,
  value.var = "reconstructed_predicted_displacement")
R_wide <- dcast(reconstructed, observation_id ~ program,
  value.var = "reference_footprint")
Dhat <- as.matrix(Dhat_wide[, ..programs]); rownames(Dhat) <-
  Dhat_wide$observation_id
R <- as.matrix(R_wide[, ..programs]); rownames(R) <- R_wide$observation_id
R <- R[rownames(Dhat), programs, drop = FALSE]

decoder <- fread(file.path(stage_result("05_program_function_bridge"),
  "function_decoder_coefficients.tsv"))
intercept <- decoder[term == "(Intercept)", coefficient]
beta <- setNames(decoder[term %chin% programs, coefficient],
  decoder[term %chin% programs, term])
if (length(intercept) != 1L || length(beta) != length(programs))
  fail("The accepted shared response-to-function decoder is incomplete.")
new_function <- intercept + as.numeric(Dhat %*% beta[programs])
reference_function <- intercept + as.numeric(R %*% beta[programs])

held_molecular <- fread(file.path(stage_result(
  "14_continuous_baseline_susceptibility"),
  "held_function_relevant_molecular_predictions.tsv.gz"))
held_position <- match(rownames(Dhat), held_molecular$observation_id)
if (anyNA(held_position)) fail("Stage 14 molecular tables do not align.")
formula_prediction <- reference_function +
  held_molecular$predicted_z[held_position]
decoder_composition_error <- max(abs(new_function - formula_prediction))
if (!is.finite(decoder_composition_error) ||
    decoder_composition_error > 1e-10)
  fail("Frozen decoder composition is not f(R) + predicted molecular z.")

## Align the six prespecified functional comparisons.  The fifth prediction is
## the explicit f(R)+measured-z composition and the sixth is direct application
## of the same shared decoder to measured displacement; equality is audited.
accepted <- fread(file.path(stage_artifact("05_program_function_bridge"),
  "measured_function_predictions.tsv.gz"))
baseline_function <- fread(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_function_predictions.tsv.gz"))[baseline_source == "tahoe_dmso"]
state_function <- fread(file.path(stage_artifact("13_separated_chain"),
  "separated_chain_held_predictions.tsv.gz"))[
    prediction_model == "state_predicted_displacement_shared_decoder"]
index_baseline <- match(accepted$observation_id,
  baseline_function$observation_id)
index_state <- match(accepted$observation_id, state_function$observation_id)
index_new <- match(accepted$observation_id, rownames(Dhat))
if (anyNA(index_baseline) || anyNA(index_state))
  fail("An existing frozen functional comparator did not align.")
new_values <- fifelse(is.na(index_new), NA_real_, new_function[index_new])
new_molecular_z <- fifelse(is.na(index_new), NA_real_,
  held_molecular$predicted_z[match(accepted$observation_id,
    held_molecular$observation_id)])
new_molecular_eligible <- fifelse(is.na(index_new), FALSE,
  held_molecular$molecular_eligible[match(accepted$observation_id,
    held_molecular$observation_id)])
new_reference_function <- fifelse(is.na(index_new), NA_real_,
  reference_function[index_new])
explicit_complete <- accepted$mean_prediction + accepted$measured_z
complete_identity_error <- max(abs(explicit_complete -
  accepted$measured_prediction), na.rm = TRUE)
if (!is.finite(complete_identity_error) || complete_identity_error > 1e-10)
  fail("Measured-z composition does not equal the shared-decoder ceiling.")

base <- accepted[, .(observation_id, depmap_id, cell_line,
  perturbation_id, drug, dose, observed_function,
  training_reference_function, evaluation_included, activity_supported,
  mean_prediction, measured_z)]
base[, `:=`(continuous_molecular_eligible = new_molecular_eligible,
  continuous_predicted_molecular_z = new_molecular_z,
  continuous_reference_function = new_reference_function)]
model_values <- list(
  mean_perturbation_footprint = accepted$mean_prediction,
  existing_baseline_conditioned =
    baseline_function$predicted_function[index_baseline],
  failed_C1_C2_state_modulated =
    state_function$predicted_function[index_state],
  continuous_susceptibility = new_values,
  complete_measured_displacement = explicit_complete,
  measured_displacement_shared_decoder_ceiling =
    accepted$measured_prediction)
predictions <- rbindlist(lapply(names(model_values), function(model)
  cbind(copy(base), prediction_model = model,
    predicted_function = model_values[[model]])), use.names = TRUE)
fwrite(predictions, file.path(RESULT,
  "continuous_susceptibility_chain_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")

models <- names(model_values)
performance <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort)
  rbindlist(lapply(models, function(model)
    evaluate_function(predictions, model, cohort)))))
fwrite(performance, file.path(RESULT,
  "continuous_susceptibility_chain_performance.tsv"), sep = "\t")

per_drug <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort) {
  rbindlist(lapply(models, function(model) {
    current <- predictions[prediction_model == model &
      evaluation_included == TRUE]
    if (cohort == "activity_supported")
      current <- current[activity_supported == TRUE]
    current[, {
      within <- function_residual_table(.SD, "within_perturbation",
        "predicted_function")
      strict <- function_residual_table(.SD,
        "strict_cell_line_by_perturbation", "predicted_function")
      list(n_rows = .N, n_conditions = uniqueN(perturbation_id),
        pooled_pearson = safe_cor(observed_function, predicted_function),
        within_pearson = safe_cor(within$observed_residual,
          within$predicted_residual),
        strict_pearson = safe_cor(strict$observed_residual,
          strict$predicted_residual),
        observed_sd = sd(observed_function),
        predicted_sd = sd(predicted_function),
        predicted_to_observed_sd_ratio = sd(predicted_function) /
          sd(observed_function))
    }, by = drug][, `:=`(cohort = cohort, prediction_model = model)]
  }))
}))
setcolorder(per_drug, c("cohort", "prediction_model", "drug",
  setdiff(names(per_drug), c("cohort", "prediction_model", "drug"))))
fwrite(per_drug, file.path(RESULT,
  "continuous_susceptibility_chain_per_drug.tsv.gz"), sep = "\t",
  compress = "gzip")

per_line <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort) {
  rbindlist(lapply(models, function(model) {
    current <- predictions[prediction_model == model &
      evaluation_included == TRUE]
    if (cohort == "activity_supported")
      current <- current[activity_supported == TRUE]
    current[, .(n_rows = .N,
      pooled_pearson = safe_cor(observed_function, predicted_function),
      observed_context_pearson = safe_cor(
        observed_function - training_reference_function,
        predicted_function - mean_prediction),
      observed_sd = sd(observed_function),
      predicted_sd = sd(predicted_function),
      predicted_to_observed_sd_ratio = sd(predicted_function) /
        sd(observed_function)), by = .(depmap_id, cell_line)][,
          `:=`(cohort = cohort, prediction_model = model)]
  }))
}))
setcolorder(per_line, c("cohort", "prediction_model", "depmap_id",
  "cell_line", setdiff(names(per_line), c("cohort", "prediction_model",
    "depmap_id", "cell_line"))))
fwrite(per_line, file.path(RESULT,
  "continuous_susceptibility_chain_per_line.tsv"), sep = "\t")

## Metric-specific fraction of the measured functional-context gain recovered.
gain <- rbindlist(lapply(c("complete_matched_panel", "activity_supported"),
  function(cohort_name) rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric_name) {
    current <- performance[cohort == cohort_name & metric == metric_name]
    mean_sse <- current[prediction_model == "mean_perturbation_footprint",
      mse * n_rows]
    ceiling_sse <- current[prediction_model ==
      "measured_displacement_shared_decoder_ceiling", mse * n_rows]
    current[, .(cohort, prediction_model, metric, n_rows,
      sse = mse * n_rows,
      measured_functional_context_gain_recovered =
        (mean_sse - mse * n_rows) / (mean_sse - ceiling_sse))]
  }))))
fwrite(gain, file.path(RESULT,
  "functional_context_gain_recovered.tsv"), sep = "\t")

## Molecular comparison on exactly shared held rows and panels.
baseline_molecular <- fread(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_molecular_predictions.tsv.gz"))[baseline_source == "tahoe_dmso"]
state_molecular <- fread(file.path(stage_artifact("11_state_modulated_response"),
  "state_modulated_molecular_predictions.tsv.gz"))[
    evaluation_role == "outer_held_out"]
state_molecular[, context_contribution :=
  (predicted_displacement - reference_footprint) * beta[program]]
state_z <- state_molecular[, .(failed_C1_C2_predicted_z =
  sum(context_contribution)), by = observation_id]
state_molecular[, observed_contribution :=
  (observed_displacement - reference_footprint) * beta[program]]
state_observed <- state_molecular[, .(state_observed_z =
  sum(observed_contribution)), by = observation_id]
molecular <- held_molecular[, .(observation_id, depmap_id, cell_line,
  perturbation_id, drug, dose, molecular_eligible, activity_supported,
  observed_z, continuous_predicted_z = predicted_z)]
molecular <- merge(molecular, baseline_molecular[, .(observation_id,
  existing_baseline_predicted_z = predicted_z,
  existing_baseline_observed_z = measured_z)], by = "observation_id",
  all.x = TRUE, sort = FALSE)
molecular <- merge(molecular, state_z, by = "observation_id",
  all.x = TRUE, sort = FALSE)
molecular <- merge(molecular, state_observed, by = "observation_id",
  all.x = TRUE, sort = FALSE)
matched_activity_ids <- accepted[evaluation_included == TRUE &
  activity_supported == TRUE, observation_id]
molecular[, matched_activity_panel := observation_id %chin%
  matched_activity_ids]
observed_alignment_error <- max(c(abs(molecular$observed_z -
  molecular$existing_baseline_observed_z), abs(molecular$observed_z -
  molecular$state_observed_z)), na.rm = TRUE)
if (!is.finite(observed_alignment_error) || observed_alignment_error > 1e-9)
  fail("Existing molecular comparators use a different measured-z target.")

molecular_models <- c("existing_baseline_predicted_z",
  "failed_C1_C2_predicted_z", "continuous_predicted_z")
molecular_cohorts <- c("held_all_molecular_rows",
  "held_activity_supported", "held_molecular_eligible",
  "held_eligible_and_activity_supported", "held_matched_activity_panel")
molecular_comparison <- rbindlist(lapply(molecular_cohorts,
  function(cohort) rbindlist(lapply(molecular_models, function(model)
    molecular_metrics(molecular, model, cohort)))))

function_comparison <- performance[prediction_model %chin%
    c("existing_baseline_conditioned", "failed_C1_C2_state_modulated",
      "continuous_susceptibility"),
  .(endpoint = "end_to_end_function", cohort, metric, prediction_model,
    n_rows, n_conditions, n_drugs, n_lines, pearson, observed_sd,
    predicted_sd, predicted_to_observed_sd_ratio, mse, rmse,
    correlation_method)]
comparison <- rbindlist(list(molecular_comparison, function_comparison),
  use.names = TRUE, fill = TRUE)
comparison[, existing_bridge_pearson := {
  existing_name <- if (endpoint[1L] == "molecular_z")
    "existing_baseline_predicted_z" else "existing_baseline_conditioned"
  pearson[prediction_model == existing_name][1L]
}, by = .(endpoint, cohort, metric)]
comparison[, change_from_existing_baseline_bridge := pearson -
  existing_bridge_pearson]
fwrite(comparison, file.path(RESULT,
  "comparison_with_existing_baseline_bridge.tsv"), sep = "\t")

matched_molecular <- molecular_comparison[
  cohort == "held_matched_activity_panel"]
eligible_molecular <- molecular_comparison[
  cohort == "held_molecular_eligible"]
activity_function <- performance[cohort == "activity_supported"]
value_at <- function(table, model, metric = NULL, column = "pearson") {
  current <- table[prediction_model == model]
  if (!is.null(metric)) {
    metric_name <- metric
    current <- current[metric == metric_name]
  }
  current[[column]][1L]
}

## Distribution checks use the same matched-activity molecular rows.
line_molecular <- molecular[matched_activity_panel == TRUE, .(
  continuous_r = safe_cor(observed_z, continuous_predicted_z),
  C1_C2_r = safe_cor(observed_z, failed_C1_C2_predicted_z),
  existing_bridge_r = safe_cor(observed_z, existing_baseline_predicted_z)),
  by = .(depmap_id, cell_line)]
perturbation_molecular <- molecular[molecular_eligible == TRUE, .(
  continuous_r = safe_cor(observed_z, continuous_predicted_z),
  C1_C2_r = safe_cor(observed_z, failed_C1_C2_predicted_z),
  existing_bridge_r = safe_cor(observed_z, existing_baseline_predicted_z)),
  by = .(perturbation_id, drug, dose)]

matched_continuous_r <- value_at(matched_molecular,
  "continuous_predicted_z")
matched_state_r <- value_at(matched_molecular,
  "failed_C1_C2_predicted_z")
matched_existing_r <- value_at(matched_molecular,
  "existing_baseline_predicted_z")
eligible_continuous_r <- value_at(eligible_molecular,
  "continuous_predicted_z")
eligible_amplitude <- value_at(eligible_molecular,
  "continuous_predicted_z", column = "predicted_to_observed_sd_ratio")
new_within <- value_at(activity_function, "continuous_susceptibility",
  "within_perturbation")
existing_within <- value_at(activity_function,
  "existing_baseline_conditioned", "within_perturbation")
new_strict <- value_at(activity_function, "continuous_susceptibility",
  "strict_cell_line_by_perturbation")
existing_strict <- value_at(activity_function,
  "existing_baseline_conditioned", "strict_cell_line_by_perturbation")
new_pooled <- value_at(activity_function, "continuous_susceptibility",
  "pooled")
existing_pooled <- value_at(activity_function,
  "existing_baseline_conditioned", "pooled")
n_positive_lines <- line_molecular[is.finite(continuous_r),
  sum(continuous_r > 0)]
n_better_lines <- line_molecular[is.finite(continuous_r) &
  is.finite(C1_C2_r), sum(continuous_r > C1_C2_r)]
n_positive_perturbations <- perturbation_molecular[
  is.finite(continuous_r), sum(continuous_r > 0)]

susceptibility_supported <- is.finite(eligible_continuous_r) &&
  eligible_continuous_r > 0.05 &&
  is.finite(matched_continuous_r) && matched_continuous_r > 0 &&
  matched_continuous_r > matched_state_r &&
  matched_continuous_r > matched_existing_r &&
  is.finite(eligible_amplitude) && eligible_amplitude > 0.05 &&
  n_positive_lines >= 3L && n_positive_perturbations >= 5L
composition_supported <- is.finite(new_within) && is.finite(new_strict) &&
  new_within > existing_within && new_strict > existing_strict
chain_supported <- susceptibility_supported && composition_supported
failure_location <- if (!susceptibility_supported)
  "function-relevant molecular susceptibility" else if
    (!composition_supported) "downstream composition" else "none"
stopping_rule <- if (!susceptibility_supported)
  "stop baseline-conditioning modifications" else if (!composition_supported)
    "retain molecular susceptibility result; downstream chain unsupported" else
      "continuous susceptibility chain supported"

summary <- data.table(
  eligible_molecular_z_pearson = eligible_continuous_r,
  matched_activity_molecular_z_pearson = matched_continuous_r,
  matched_activity_C1_C2_molecular_z_pearson = matched_state_r,
  matched_activity_existing_bridge_molecular_z_pearson = matched_existing_r,
  eligible_molecular_z_sd_ratio = eligible_amplitude,
  n_held_lines_with_positive_matched_molecular_z = n_positive_lines,
  n_held_lines_beating_C1_C2 = n_better_lines,
  n_eligible_perturbations_with_positive_molecular_z =
    n_positive_perturbations,
  activity_function_pooled_continuous = new_pooled,
  activity_function_pooled_existing_baseline = existing_pooled,
  activity_function_within_continuous = new_within,
  activity_function_within_existing_baseline = existing_within,
  activity_function_strict_continuous = new_strict,
  activity_function_strict_existing_baseline = existing_strict,
  susceptibility_supported = susceptibility_supported,
  downstream_composition_supported = composition_supported,
  continuous_chain_supported = chain_supported,
  failure_location = failure_location, stopping_rule = stopping_rule,
  interpretation = paste("Baseline program Sb potentiates or suppresses",
    "the function-relevant molecular consequence of perturbations that",
    "engage response program Pj; no direct baseline-to-function path."))
fwrite(summary, file.path(RESULT,
  "continuous_susceptibility_chain_summary.tsv"), sep = "\t")

## Four owned, one-page diagnostics.
model_labels <- c(mean_perturbation_footprint = "Mean footprint",
  existing_baseline_conditioned = "Existing baseline-conditioned",
  failed_C1_C2_state_modulated = "Failed C1/C2",
  continuous_susceptibility = "Continuous susceptibility",
  complete_measured_displacement = "Complete measured (f(R)+z)",
  measured_displacement_shared_decoder_ceiling =
    "Measured shared-decoder ceiling")
metric_labels <- c(pooled = "Pooled", within_perturbation = "Within",
  strict_cell_line_by_perturbation = "Strict")
plot_metric <- copy(activity_function[is.finite(pearson)])
plot_metric[, `:=`(model_label = model_labels[prediction_model],
  metric_label = metric_labels[metric])]
p_metric <- ggplot(plot_metric,
    aes(pearson, factor(model_label, levels = rev(model_labels)),
      color = metric_label, shape = metric_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(size = 2.0, position = position_dodge(width = 0.45)) +
  labs(x = "Held activity-panel Pearson r", y = NULL,
    color = NULL, shape = NULL,
    title = "Frozen continuous-susceptibility chain metrics") +
  theme_stage15()
save_gg("15a_final_metric_comparison.pdf", p_metric, 7.2, 4.9)

plot_function <- predictions[prediction_model ==
  "continuous_susceptibility" & evaluation_included == TRUE &
  activity_supported == TRUE]
p_function <- ggplot(plot_function,
    aes(observed_function, predicted_function)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.32, size = 0.75, color = "#2166AC") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.65, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("pooled r=%s\nwithin r=%s; strict r=%s",
      fmt(new_pooled), fmt(new_within), fmt(new_strict)), size = 2.35) +
  labs(x = "Observed function", y = "Predicted function",
    title = "Continuous susceptibility on the activity panel") +
  square_scatter_theme()
save_gg("15b_observed_vs_predicted_activity_function.pdf", p_function,
  5.0, 4.8)

drug_change <- dcast(per_drug[cohort == "activity_supported" &
    prediction_model %chin% c("existing_baseline_conditioned",
      "continuous_susceptibility")], drug ~ prediction_model,
  value.var = "within_pearson")
drug_change[, change := continuous_susceptibility -
  existing_baseline_conditioned]
drug_change <- drug_change[is.finite(change)]
setorder(drug_change, change)
drug_change[, rank := seq_len(.N)]
p_drug <- ggplot(drug_change, aes(rank, change)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
  geom_point(alpha = 0.72, size = 1.15, color = "#B2182B") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.5, color = "black") +
  labs(x = "Activity-panel drug rank",
    y = "Change in within-drug Pearson r",
    title = "Per-drug change from the existing baseline bridge") +
  theme_stage15()
save_gg("15c_per_drug_change_from_existing_baseline.pdf", p_drug,
  6.0, 4.3)

molecular_by_drug <- molecular[matched_activity_panel == TRUE, .(
  molecular_z_pearson = safe_cor(observed_z, continuous_predicted_z)),
  by = drug]
functional_by_drug <- dcast(per_drug[cohort == "activity_supported" &
    prediction_model %chin% c("existing_baseline_conditioned",
      "continuous_susceptibility")], drug ~ prediction_model,
  value.var = "within_pearson")
functional_by_drug[, function_within_change :=
  continuous_susceptibility - existing_baseline_conditioned]
link_plot <- merge(molecular_by_drug, functional_by_drug, by = "drug")
link_plot <- link_plot[is.finite(molecular_z_pearson) &
  is.finite(function_within_change)]
p_link <- ggplot(link_plot,
    aes(molecular_z_pearson, function_within_change)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_point(alpha = 0.7, size = 1.2, color = "#6A51A3") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("Across-drug r=%s",
      fmt(safe_cor(link_plot$molecular_z_pearson,
        link_plot$function_within_change))), size = 2.35) +
  labs(x = "Predicted vs observed molecular-z Pearson r",
    y = "Within-function r change vs existing baseline",
    title = "Molecular susceptibility and end-to-end performance") +
  square_scatter_theme()
save_gg("15d_molecular_z_vs_end_to_end_function.pdf", p_link, 5.0, 4.8)

figures <- data.table(pdf = c("15a_final_metric_comparison.pdf",
  "15b_observed_vs_predicted_activity_function.pdf",
  "15c_per_drug_change_from_existing_baseline.pdf",
  "15d_molecular_z_vs_end_to_end_function.pdf"),
  producing_script =
    "scripts/15_evaluate_continuous_susceptibility_chain.R")
figures[, exists := file.exists(file.path(RESULT, pdf))]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")

checks <- data.table(check = c("six_requested_comparators",
  "frozen_shared_decoder_only", "decoder_equals_reference_plus_predicted_z",
  "complete_measured_composition_identity", "molecular_targets_align",
  "both_evaluation_panels", "all_three_function_metrics",
  "all_function_passes_through_program_displacement",
  "all_required_pdfs_created"),
  value = c(setequal(models, c("mean_perturbation_footprint",
      "existing_baseline_conditioned", "failed_C1_C2_state_modulated",
      "continuous_susceptibility", "complete_measured_displacement",
      "measured_displacement_shared_decoder_ceiling")),
    TRUE, decoder_composition_error <= 1e-10,
    complete_identity_error <= 1e-10, observed_alignment_error <= 1e-9,
    setequal(unique(performance$cohort), c("complete_matched_panel",
      "activity_supported")),
    setequal(unique(performance$metric), c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation")), TRUE, all(figures$exists)),
  detail = c("mean, existing baseline, failed C1/C2, continuous, explicit measured composition, measured ceiling",
    "accepted Stage 05 ridge coefficients; no fitting in Stage 15",
    sprintf("maximum absolute error %.3g", decoder_composition_error),
    sprintf("maximum absolute error %.3g", complete_identity_error),
    sprintf("maximum absolute error %.3g", observed_alignment_error),
    "complete matched and activity-supported panels",
    "pooled, within-perturbation and strict context",
    "continuous function = shared decoder(reconstructed nine-program displacement)",
    "four one-page-owned frozen-chain diagnostics"))
fwrite(checks, file.path(RESULT,
  "continuous_susceptibility_chain_required_checks.tsv"), sep = "\t")
if (!all(checks$value)) fail("A required Stage 15 check failed.")
writeLines("complete", file.path(RESULT, "stage15_complete.flag"))
message(sprintf(paste0("Stage 15 complete: matched activity molecular z r=%s; ",
  "activity within %.3f vs %.3f; strict %.3f vs %.3f; failure location: %s."),
  fmt(matched_continuous_r), new_within, existing_within, new_strict,
  existing_strict, failure_location))
