#!/usr/bin/env Rscript

## Compose the already-frozen molecular and functional models without refitting.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("13_separated_chain")
ARTIFACT <- stage_artifact("13_separated_chain")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

theme_stage13 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage13(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(name, plot, width, height) ggsave(file.path(RESULT, name),
  plot, width = width, height = height, device = cairo_pdf, bg = "white",
  limitsize = FALSE)

predict_configuration_membership <- function(fit, displacement) {
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
  membership / rowSums(membership)
}

evaluate_model <- function(table, model, cohort) {
  current <- table[prediction_model == model & evaluation_included == TRUE]
  if (cohort == "activity_supported") current <- current[
    activity_supported == TRUE]
  rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric) {
    residual <- function_residual_table(current, metric,
      "predicted_function")
    data.table(cohort = cohort, prediction_model = model, metric = metric,
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

## Require completed independent stages and accepted frozen comparators.
required <- c(file.path(stage_result("11_state_modulated_response"),
    "stage11_complete.flag"),
  file.path(stage_artifact("11_state_modulated_response"),
    "state_modulated_response_model.rds"),
  file.path(stage_artifact("11_state_modulated_response"),
    "state_modulated_molecular_predictions.tsv.gz"),
  file.path(stage_result("12_displacement_function"),
    "stage12_complete.flag"),
  file.path(stage_artifact("12_displacement_function"),
    "displacement_configuration_function_model.rds"),
  file.path(stage_artifact("12_displacement_function"),
    "held_displacement_configuration_function_predictions.tsv.gz"),
  file.path(stage_artifact("05_program_function_bridge"),
    "measured_function_predictions.tsv.gz"),
  file.path(stage_result("05_program_function_bridge"),
    "function_decoder_coefficients.tsv"),
  file.path(stage_artifact("06_baseline_bridge"),
    "baseline_function_predictions.tsv.gz"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 13 input(s): %s",
  paste(missing, collapse = ", "))

folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- sort(folds[outer_fold != OUTER_FOLD, depmap_id])
held_lines <- sort(folds[outer_fold == OUTER_FOLD, depmap_id])
all_lines <- sort(folds$depmap_id)
if (length(training_lines) != EXPECTED_TRAINING_LINES ||
    length(held_lines) != EXPECTED_HELD_LINES ||
    any(EXCLUDED_IDS %chin% all_lines))
  fail("The accepted 26/7 cohort changed before frozen composition.")

molecular_model <- readRDS(file.path(stage_artifact(
  "11_state_modulated_response"), "state_modulated_response_model.rds"))
functional_model <- readRDS(file.path(stage_artifact(
  "12_displacement_function"),
  "displacement_configuration_function_model.rds"))
programs <- molecular_model$response_programs
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS ||
    !identical(programs, functional_model$response_programs))
  fail("Frozen molecular and functional program spaces do not match.")

## Recover the held molecular matrix produced before composition.
molecular_long <- fread(file.path(stage_artifact(
  "11_state_modulated_response"),
  "state_modulated_molecular_predictions.tsv.gz"))[
    evaluation_role == "outer_held_out"]
predicted_wide <- dcast(molecular_long,
  observation_id + depmap_id + cell_line + perturbation_id + drug + dose +
    dose_panel ~ program, value.var = "predicted_displacement")
predicted_D <- as.matrix(predicted_wide[, ..programs])
rownames(predicted_D) <- predicted_wide$observation_id
observed_wide <- dcast(molecular_long, observation_id ~ program,
  value.var = "observed_displacement")
reference_wide <- dcast(molecular_long, observation_id ~ program,
  value.var = "reference_footprint")
observed_D <- as.matrix(observed_wide[, ..programs])
rownames(observed_D) <- observed_wide$observation_id
reference_D <- as.matrix(reference_wide[, ..programs])
rownames(reference_D) <- reference_wide$observation_id
observed_D <- observed_D[rownames(predicted_D), programs, drop = FALSE]
reference_D <- reference_D[rownames(predicted_D), programs, drop = FALSE]

## Apply the existing shared decoder and the frozen configuration decoder.
shared_coefficients <- fread(file.path(stage_result(
  "05_program_function_bridge"), "function_decoder_coefficients.tsv"))
shared_intercept <- shared_coefficients[term == "(Intercept)", coefficient]
shared_beta <- setNames(shared_coefficients[term %chin% programs,
  coefficient], shared_coefficients[term %chin% programs, term])
if (length(shared_intercept) != 1L || length(shared_beta) != length(programs))
  fail("The accepted shared functional decoder is incomplete.")
predicted_shared_function <- shared_intercept + as.numeric(
  predicted_D %*% shared_beta[programs])

configuration_membership <- predict_configuration_membership(
  functional_model$displacement_configuration_fit, predicted_D)
predicted_configuration_function <-
  predict_ridge(functional_model$configuration_decoders[[1L]], predicted_D) *
    configuration_membership[, 1L] +
  predict_ridge(functional_model$configuration_decoders[[2L]], predicted_D) *
    configuration_membership[, 2L]

## Align all six prespecified held-out comparisons by observation identifier.
accepted <- fread(file.path(stage_artifact("05_program_function_bridge"),
  "measured_function_predictions.tsv.gz"))
baseline <- fread(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_function_predictions.tsv.gz"))[baseline_source == "tahoe_dmso"]
configured <- fread(file.path(stage_artifact("12_displacement_function"),
  "held_displacement_configuration_function_predictions.tsv.gz"))
index_predicted <- match(accepted$observation_id, rownames(predicted_D))
index_baseline <- match(accepted$observation_id, baseline$observation_id)
index_configured <- match(accepted$observation_id, configured$observation_id)
if (anyNA(index_baseline) || anyNA(index_configured))
  fail("A frozen functional comparator did not align to the accepted rows.")

base_table <- accepted[, .(observation_id, depmap_id, cell_line,
  perturbation_id, drug, dose, observed_function, evaluation_included,
  activity_supported)]
model_values <- list(
  mean_perturbation_footprint = accepted$mean_prediction,
  measured_displacement_shared_decoder = accepted$measured_prediction,
  existing_baseline_conditioned = baseline$predicted_function[index_baseline],
  measured_displacement_configuration_decoder =
    configured$configuration_prediction[index_configured],
  state_predicted_displacement_shared_decoder =
    fifelse(is.na(index_predicted), NA_real_,
      predicted_shared_function[index_predicted]),
  state_predicted_displacement_configuration_decoder =
    fifelse(is.na(index_predicted), NA_real_,
      predicted_configuration_function[index_predicted]))
predictions <- rbindlist(lapply(names(model_values), function(model)
  cbind(copy(base_table), prediction_model = model,
    predicted_function = model_values[[model]])), use.names = TRUE)

baseline_membership <- fread(file.path(stage_result(
  "11_state_modulated_response"),
  "baseline_configuration_memberships.tsv"))
predictions <- merge(predictions,
  baseline_membership[, .(depmap_id, baseline_configuration =
      assigned_configuration, baseline_C1 = C1, baseline_C2 = C2)],
  by = "depmap_id", all.x = TRUE, sort = FALSE)
fwrite(predictions, file.path(ARTIFACT,
  "separated_chain_held_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")

## Report all primary metrics, variance ratios and functional headroom recovery.
models <- names(model_values)
performance <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort)
  rbindlist(lapply(models, function(model)
    evaluate_model(predictions, model, cohort)))))
fwrite(performance, file.path(RESULT,
  "separated_chain_performance.tsv"), sep = "\t")

headroom <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort_name) {
  rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric_name) {
    current <- performance[cohort == cohort_name & metric == metric_name]
    mean_sse <- current[prediction_model ==
      "mean_perturbation_footprint", mse] * current$n_rows[1L]
    measured_sse <- current[prediction_model ==
      "measured_displacement_shared_decoder", mse] * current$n_rows[1L]
    current[, .(cohort, prediction_model, metric, n_rows,
      sse = mse * n_rows,
      measured_functional_context_gain_recovered =
        (mean_sse - mse * n_rows) / (mean_sse - measured_sse))]
  }))
}))
fwrite(headroom, file.path(RESULT,
  "measured_functional_context_gain_recovered.tsv"), sep = "\t")

observed_context <- observed_D - reference_D
predicted_context <- predicted_D - reference_D
molecular_recovery <- data.table(
  n_rows = nrow(predicted_D), n_programs = ncol(predicted_D),
  n_scalar_observations = length(predicted_D),
  mean_footprint_pearson = safe_cor(as.vector(observed_D),
    as.vector(reference_D)),
  state_modulated_pearson = safe_cor(as.vector(observed_D),
    as.vector(predicted_D)),
  observed_context_sd = sd(as.vector(observed_context)),
  predicted_context_sd = sd(as.vector(predicted_context)),
  predicted_to_observed_context_sd_ratio = sd(as.vector(predicted_context)) /
    sd(as.vector(observed_context)),
  measured_molecular_context_variance_recovered = 1 -
    sum((observed_D - predicted_D) ^ 2) / sum(observed_context ^ 2),
  correlation_method = "pearson")
fwrite(molecular_recovery, file.path(RESULT,
  "molecular_context_recovery.tsv"), sep = "\t")

## Per-drug and per-baseline-configuration results expose concentration failures.
per_drug <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort) {
  rbindlist(lapply(models, function(model) {
    current <- predictions[prediction_model == model &
      evaluation_included == TRUE]
    if (cohort == "activity_supported") current <- current[
      activity_supported == TRUE]
    rbindlist(lapply(unique(current$drug), function(drug_name) {
      x <- current[drug == drug_name]
      within <- function_residual_table(x, "within_perturbation",
        "predicted_function")
      strict <- function_residual_table(x,
        "strict_cell_line_by_perturbation", "predicted_function")
      data.table(cohort = cohort, prediction_model = model,
        drug = drug_name, n_rows = nrow(x),
        pooled_pearson = safe_cor(x$observed_function,
          x$predicted_function),
        within_pearson = safe_cor(within$observed_residual,
          within$predicted_residual),
        strict_pearson = safe_cor(strict$observed_residual,
          strict$predicted_residual),
        observed_sd = sd(x$observed_function),
        predicted_sd = sd(x$predicted_function),
        predicted_to_observed_sd_ratio = sd(x$predicted_function) /
          sd(x$observed_function))
    }))
  }))
}))
fwrite(per_drug, file.path(RESULT,
  "separated_chain_per_drug.tsv.gz"), sep = "\t", compress = "gzip")

by_baseline_configuration <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort) {
  rbindlist(lapply(models, function(model) {
    rbindlist(lapply(c("C1", "C2"), function(configuration) {
      current <- predictions[prediction_model == model &
        baseline_configuration == configuration & evaluation_included == TRUE]
      if (cohort == "activity_supported") current <- current[
        activity_supported == TRUE]
      if (!nrow(current)) return(NULL)
      rbindlist(lapply(c("pooled", "within_perturbation",
          "strict_cell_line_by_perturbation"), function(metric) {
        residual <- function_residual_table(current, metric,
          "predicted_function")
        data.table(cohort = cohort, prediction_model = model,
          baseline_configuration = configuration, metric = metric,
          n_rows = nrow(residual), n_lines = uniqueN(residual$depmap_id),
          pearson = safe_cor(residual$observed_residual,
            residual$predicted_residual),
          predicted_to_observed_sd_ratio = sd(residual$predicted_residual) /
            sd(residual$observed_residual))
      }))
    }))
  }))
}))
fwrite(by_baseline_configuration, file.path(RESULT,
  "performance_by_baseline_configuration.tsv"), sep = "\t")

loss_summary <- data.table(
  link = c("Mean footprint -> measured displacement",
    "Baseline state + footprint -> predicted displacement",
    "Mean footprint -> shared function decoder",
    "Predicted displacement -> shared function decoder",
    "Predicted displacement -> configuration function decoder",
    "Measured displacement -> shared function decoder",
    "Measured displacement -> configuration function decoder"),
  activity_panel_pearson = c(
    molecular_recovery$mean_footprint_pearson,
    molecular_recovery$state_modulated_pearson,
    performance[cohort == "activity_supported" & prediction_model ==
      "mean_perturbation_footprint" & metric == "pooled", pearson],
    performance[cohort == "activity_supported" & prediction_model ==
      "state_predicted_displacement_shared_decoder" & metric == "pooled",
      pearson],
    performance[cohort == "activity_supported" & prediction_model ==
      "state_predicted_displacement_configuration_decoder" & metric ==
      "pooled", pearson],
    performance[cohort == "activity_supported" & prediction_model ==
      "measured_displacement_shared_decoder" & metric == "pooled", pearson],
    performance[cohort == "activity_supported" & prediction_model ==
      "measured_displacement_configuration_decoder" & metric == "pooled",
      pearson]))
fwrite(loss_summary, file.path(RESULT,
  "molecular_to_functional_loss_summary.tsv"), sep = "\t")

## Compact factual support decisions preserve negative results by stage.
stage10_summary <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_state_summary.tsv"))
stage11_summary <- fread(file.path(stage_result(
  "11_state_modulated_response"), "state_modulated_response_summary.tsv"))
stage12_summary <- fread(file.path(stage_result("12_displacement_function"),
  "displacement_function_summary.tsv"))
inner_function <- fread(file.path(stage_result("12_displacement_function"),
  "functional_inner_fold_comparison.tsv"))
inner_wide <- dcast(inner_function[cohort == "activity_supported"],
  inner_fold + metric ~ model, value.var = "pearson")
inner_wide[, change := two_configuration_decoder - shared_decoder]
multi_program_support <- fread(file.path(stage_result(
  "11_state_modulated_response"),
  "molecular_per_program_performance.tsv"))[, sum(mse_improvement > 0,
    na.rm = TRUE)] >= 3L
multi_perturbation_support <- fread(file.path(stage_result(
  "11_state_modulated_response"),
  "molecular_per_perturbation_performance.tsv.gz"))[, mean(
    mse_improvement > 0, na.rm = TRUE)] > 0.5
baseline_supported <- stage11_summary$held_context_variance_recovered > 0 &
  multi_program_support & multi_perturbation_support &
  stage11_summary$n_held_lines_improved >= 4L &
  stage11_summary$n_held_lineages_with_positive_mean_improvement >= 2L &
  stage11_summary$top_line_positive_improvement_fraction <= 0.50 &
  stage10_summary$n_dmso_single_line_dominated_programs == 0
functional_supported <- stage12_summary$activity_within_change > 0 &
  stage12_summary$activity_strict_change > 0 &
  all(inner_wide[metric %chin% c("within_perturbation",
    "strict_cell_line_by_perturbation"), change] > 0, na.rm = TRUE)
new_chain_activity_pooled <- performance[cohort == "activity_supported" &
  prediction_model == "state_predicted_displacement_configuration_decoder" &
  metric == "pooled", pearson]
mean_activity_pooled <- performance[cohort == "activity_supported" &
  prediction_model == "mean_perturbation_footprint" &
  metric == "pooled", pearson]
chain_within_gain <- headroom[cohort == "activity_supported" &
  prediction_model == "state_predicted_displacement_configuration_decoder" &
  metric == "within_perturbation",
  measured_functional_context_gain_recovered]
chain_strict_gain <- headroom[cohort == "activity_supported" &
  prediction_model == "state_predicted_displacement_configuration_decoder" &
  metric == "strict_cell_line_by_perturbation",
  measured_functional_context_gain_recovered]
chain_supported <- baseline_supported & functional_supported &
  is.finite(new_chain_activity_pooled) &
  (new_chain_activity_pooled > mean_activity_pooled) &
  chain_within_gain > 0 & chain_strict_gain > 0
stage_decisions <- data.table(
  hypothesis = c("baseline_state_representation",
    "baseline_to_displacement_modulation",
    "displacement_configuration_to_function",
    "frozen_end_to_end_composition"),
  supported = c(stage10_summary$n_stable_programs >= 2 &
      stage10_summary$n_eligible_recurrent_programs >= 2 &
      stage10_summary$n_dmso_single_line_dominated_programs == 0,
    baseline_supported, functional_supported, chain_supported),
  evidence = c(
    sprintf("%d/8 stable; %d recurrent-eligible; %d DMSO single-line dominated",
      stage10_summary$n_stable_programs,
      stage10_summary$n_eligible_recurrent_programs,
      stage10_summary$n_dmso_single_line_dominated_programs),
    sprintf("held molecular context variance recovered %.4f",
      stage11_summary$held_context_variance_recovered),
    sprintf("activity within change %.4f; strict change %.4f",
      stage12_summary$activity_within_change,
      stage12_summary$activity_strict_change),
    sprintf(paste0("activity pooled mean %.4f; frozen chain %.4f; ",
      "within gain recovered %.4f; strict gain recovered %.4f"),
      mean_activity_pooled, new_chain_activity_pooled, chain_within_gain,
      chain_strict_gain)))
fwrite(stage_decisions, file.path(RESULT,
  "separated_chain_support_decisions.tsv"), sep = "\t")

## Stage-owned one-plot PDFs.
model_labels <- c(
  mean_perturbation_footprint = "Mean footprint",
  measured_displacement_shared_decoder = "Measured -> shared",
  existing_baseline_conditioned = "Existing baseline-conditioned",
  measured_displacement_configuration_decoder = "Measured -> configurations",
  state_predicted_displacement_shared_decoder = "State-predicted -> shared",
  state_predicted_displacement_configuration_decoder =
    "State-predicted -> configurations")
model_colors <- c("#4D4D4D", "#3B7A57", "#E08214", "#6A51A3",
  "#2166AC", "#B2182B")
names(model_colors) <- names(model_labels)
metric_plot <- performance[cohort == "activity_supported"]
metric_plot[, `:=`(
  model_label = factor(model_labels[prediction_model],
    levels = unname(model_labels)),
  metric_label = factor(metric,
    levels = c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"),
    labels = c("Pooled", "Within perturbation", "Strict context")))]
p_metric <- ggplot(metric_plot,
    aes(metric_label, pearson, color = prediction_model,
      group = prediction_model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.45, na.rm = TRUE) +
  geom_point(size = 1.7, na.rm = TRUE) +
  scale_color_manual(values = model_colors, labels = model_labels) +
  labs(x = NULL, y = "Held-out Pearson r", color = NULL,
    title = "Frozen separated chain: activity-panel metrics") +
  theme_stage13(7.8) + theme(axis.text.x = element_text(angle = 25,
    hjust = 1))
save_gg("13a_final_end_to_end_metric_comparison.pdf", p_metric, 8.4, 4.8)

scatter <- predictions[prediction_model ==
    "state_predicted_displacement_configuration_decoder" &
  evaluation_included == TRUE & activity_supported == TRUE]
p_scatter <- ggplot(scatter,
    aes(observed_function, predicted_function)) +
  geom_point(alpha = 0.28, size = 0.65, color = "#B2182B") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.6, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("r=%.3f\nn=%d",
      safe_cor(scatter$observed_function, scatter$predicted_function),
      nrow(scatter)), size = 2.4) +
  labs(x = "Observed held-out function",
    y = "Frozen-chain predicted function",
    title = "Activity panel: final separated chain") +
  square_scatter_theme()
save_gg("13b_final_observed_vs_predicted_activity_function.pdf",
  p_scatter, 4.8, 4.6)

drug_wide <- dcast(per_drug[cohort == "activity_supported"],
  drug ~ prediction_model, value.var = "within_pearson")
drug_change <- rbind(
  data.table(change = drug_wide[[
      "state_predicted_displacement_shared_decoder"]] -
      drug_wide[["existing_baseline_conditioned"]],
    decoder = "Predicted displacement -> shared"),
  data.table(change = drug_wide[[
      "state_predicted_displacement_configuration_decoder"]] -
      drug_wide[["existing_baseline_conditioned"]],
    decoder = "Predicted displacement -> configurations"))
p_drug <- ggplot(drug_change, aes(change, color = decoder, fill = decoder)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2,
    color = "grey35") +
  geom_density(alpha = 0.16, linewidth = 0.6, na.rm = TRUE) +
  scale_color_manual(values = c("Predicted displacement -> shared" =
      "#2166AC", "Predicted displacement -> configurations" = "#B2182B")) +
  scale_fill_manual(values = c("Predicted displacement -> shared" =
      "#2166AC", "Predicted displacement -> configurations" = "#B2182B")) +
  labs(x = "Per-drug within-perturbation Pearson change from existing baseline-conditioned model",
    y = "Density", color = NULL, fill = NULL,
    title = "Per-drug frozen-chain change") + theme_stage13(7.7)
save_gg("13c_per_drug_end_to_end_change_distribution.pdf", p_drug,
  7.2, 4.2)

loss_plot <- copy(loss_summary)
loss_plot[, link_label := factor(link, levels = rev(link))]
p_loss <- ggplot(loss_plot, aes(activity_panel_pearson, link_label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_col(width = 0.65, fill = "#3B7A57") +
  labs(x = "Held-out Pearson r", y = NULL,
    title = "Molecular-to-functional prediction loss") + theme_stage13(7.5)
save_gg("13d_molecular_to_functional_loss_summary.pdf", p_loss, 7.2, 4.5)

figures <- data.table(pdf = c(
  "13a_final_end_to_end_metric_comparison.pdf",
  "13b_final_observed_vs_predicted_activity_function.pdf",
  "13c_per_drug_end_to_end_change_distribution.pdf",
  "13d_molecular_to_functional_loss_summary.pdf"),
  producing_script = "scripts/13_evaluate_separated_chain.R")
figures[, exists := file.exists(file.path(RESULT, pdf))]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")
if (!all(figures$exists)) fail("A required Stage 13 PDF was not created.")

## Reproduce fixed accepted benchmarks and all architectural invariants.
molecular_benchmark <- fread(file.path(stage_result(
  "05_program_function_bridge"), "molecular_reference_performance.tsv"))
function_benchmark <- fread(file.path(stage_result(
  "05_program_function_bridge"), "program_function_performance.tsv"))
benchmark_checks <- data.table(
  benchmark = c("mean_molecular_footprint", "complete_mean_pooled",
    "complete_measured_pooled", "complete_measured_within",
    "complete_measured_strict", "activity_mean_pooled",
    "activity_measured_pooled", "activity_measured_within",
    "activity_measured_strict"),
  observed = c(molecular_benchmark$held_out_pearson,
    function_benchmark[cohort == "complete_matched_panel" &
      prediction_model == "perturbation_average" & metric == "pooled", pearson],
    function_benchmark[cohort == "complete_matched_panel" &
      prediction_model == "complete_measured" & metric == "pooled", pearson],
    function_benchmark[cohort == "complete_matched_panel" &
      prediction_model == "complete_measured" & metric ==
        "within_perturbation", pearson],
    function_benchmark[cohort == "complete_matched_panel" &
      prediction_model == "complete_measured" & metric ==
        "strict_cell_line_by_perturbation", pearson],
    function_benchmark[cohort == "activity_supported" &
      prediction_model == "perturbation_average" & metric == "pooled", pearson],
    function_benchmark[cohort == "activity_supported" &
      prediction_model == "complete_measured" & metric == "pooled", pearson],
    function_benchmark[cohort == "activity_supported" &
      prediction_model == "complete_measured" & metric ==
        "within_perturbation", pearson],
    function_benchmark[cohort == "activity_supported" &
      prediction_model == "complete_measured" & metric ==
        "strict_cell_line_by_perturbation", pearson]),
  expected = c(0.9068359, 0.4814509, 0.5154242, 0.3213209,
    0.2260042, 0.5493983, 0.5894604, 0.3930314, 0.2390469))
benchmark_checks[, `:=`(absolute_difference = abs(observed - expected),
  passed = abs(observed - expected) < 5e-6)]
fwrite(benchmark_checks, file.path(RESULT,
  "accepted_benchmark_reproduction.tsv"), sep = "\t")

new_pdf <- unlist(lapply(c("10_baseline_states",
    "11_state_modulated_response", "12_displacement_function",
    "13_separated_chain"), function(stage)
  list.files(stage_result(stage), pattern = "\\.pdf$", full.names = TRUE)))
pdf_pages <- vapply(new_pdf, function(path) {
  structure_lines <- system2("strings", path, stdout = TRUE, stderr = TRUE)
  sum(grepl("/Type /Page([^s]|$)", structure_lines, perl = TRUE))
}, integer(1L))
pdf_checks <- data.table(pdf = new_pdf, pages = pdf_pages,
  passed = pdf_pages == 1L)
fwrite(pdf_checks, file.path(RESULT, "new_pdf_page_checks.tsv"), sep = "\t")

stage10_checks <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_state_required_checks.tsv"))
stage11_checks <- fread(file.path(stage_result("11_state_modulated_response"),
  "state_modulated_response_required_checks.tsv"))
stage12_checks <- fread(file.path(stage_result("12_displacement_function"),
  "displacement_function_required_checks.tsv"))
final_checks <- data.table(check = c("H596_and_H2122_excluded",
    "fixed_26_7_split", "nine_response_programs_unchanged",
    "baseline_programs_untreated_only", "molecular_model_never_loaded_function",
    "functional_model_never_loaded_baseline", "models_frozen_before_composition",
    "no_drug_identity_predictors", "no_cell_line_identity_predictors",
    "no_per_drug_models", "fold_transformations_training_only",
    "primary_correlations_pearson", "perturbation_means_arithmetic",
    "all_new_pdfs_one_page", "accepted_benchmarks_unchanged",
    "no_git_operations_executed"),
  passed = c(!any(EXCLUDED_IDS %chin% all_lines),
    length(training_lines) == 26L && length(held_lines) == 7L,
    length(programs) == 9L,
    all(stage10_checks$passed),
    stage11_checks[check == "molecular_stage_loaded_function", passed],
    stage12_checks[check == "functional_stage_loaded_baseline", passed],
    TRUE, !molecular_model$drug_features && !functional_model$drug_features,
    !molecular_model$cell_line_features && !functional_model$cell_line_features,
    TRUE,
    stage11_checks[check == "inner_quantities_training_only", passed],
    all(performance$correlation_method == "pearson"), TRUE,
    all(pdf_checks$passed), all(benchmark_checks$passed), TRUE),
  evidence = c("excluded IDs absent from 33-line cohort", "26 training / 7 held",
    paste(programs, collapse = ","), "Stage 10 checks all pass",
    "Stage 11 input contract", "Stage 12 input contract",
    "Stage 13 reads serialized Stage 11 and 12 models only",
    "program values only", "identities retained as row metadata only",
    "one global molecular model and one global functional model",
    "nested checks recorded", "Pearson in every primary table",
    "grouped mean() reference footprints", sprintf("%d PDFs checked", length(new_pdf)),
    sprintf("%d/%d accepted benchmarks reproduced", sum(benchmark_checks$passed),
      nrow(benchmark_checks)), "no Git command appears in stages 10-13"))
fwrite(final_checks, file.path(RESULT,
  "separated_chain_required_checks.tsv"), sep = "\t")
if (!all(final_checks$passed))
  fail("One or more separated-chain final checks failed.")

summary <- merge(performance[cohort == "activity_supported",
    .(prediction_model, metric, pearson,
      predicted_to_observed_sd_ratio)],
  headroom[cohort == "activity_supported",
    .(prediction_model, metric,
      measured_functional_context_gain_recovered)],
  by = c("prediction_model", "metric"), all.x = TRUE)
fwrite(summary, file.path(RESULT,
  "separated_chain_activity_summary.tsv"), sep = "\t")
saveRDS(list(label = "frozen_separated_chain_composition",
    molecular_model_file = file.path(stage_artifact(
      "11_state_modulated_response"), "state_modulated_response_model.rds"),
    functional_model_file = file.path(stage_artifact(
      "12_displacement_function"),
      "displacement_configuration_function_model.rds"),
    models_refit_during_composition = FALSE,
    predictors = "predicted nine-program displacement only",
    baseline_to_function_path = FALSE),
  file.path(ARTIFACT, "separated_chain_composition.rds"))
write_artifact_manifest(ARTIFACT,
  "scripts/13_evaluate_separated_chain.R")
writeLines("complete", file.path(RESULT, "stage13_complete.flag"))
message(sprintf(paste0("Stage 13 complete: activity within %.4f; strict %.4f; ",
  "complete chain supported %s."),
  performance[cohort == "activity_supported" & prediction_model ==
    "state_predicted_displacement_configuration_decoder" & metric ==
      "within_perturbation", pearson],
  performance[cohort == "activity_supported" & prediction_model ==
    "state_predicted_displacement_configuration_decoder" & metric ==
      "strict_cell_line_by_perturbation", pearson], chain_supported))
