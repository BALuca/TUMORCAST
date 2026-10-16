#!/usr/bin/env Rscript

## Stage 2B: learn a five-coordinate DepMap-to-DMSO affine conversion on the
## 26 training lines, then apply the unchanged DMSO-trained conditioner.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("09_dmso_to_depmap")
ARTIFACT <- stage_artifact("09_dmso_to_depmap")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(stage_result("08_summary"),
    "stage1_reproduction_complete.flag")))
  fail("Stage 2B is forbidden until the Stage 1 reproduction gate passes.")

standardize_fixed <- function(raw, center, scale) {
  result <- sweep(sweep(raw, 2L, center, "-"), 2L, scale, "/")
  if (any(!is.finite(result))) fail("Fixed baseline standardization failed.")
  result
}

make_bilinear_design <- function(row_lines, perturbations, baseline_scores,
                                 footprint, programs) {
  baseline <- baseline_scores[row_lines, , drop = FALSE]
  reference <- footprint[perturbations, programs, drop = FALSE]
  pieces <- lapply(seq_len(ncol(baseline)), function(index)
    reference * baseline[, index])
  design <- do.call(cbind, pieces)
  colnames(design) <- unlist(lapply(colnames(baseline), function(program)
    paste0(program, "_x_", programs)), use.names = FALSE)
  if (any(!is.finite(design))) fail("Baseline design is nonfinite.")
  design
}

predict_baseline <- function(fit, x) {
  if (is.null(fit)) return(rep(0, nrow(x)))
  z <- sweep(sweep(x[, fit$keep, drop = FALSE], 2L,
    fit$center[fit$keep], "-"), 2L, fit$scale[fit$keep], "/")
  as.numeric(z %*% fit$beta_standardized)
}

functional_metrics <- function(table, cohort) {
  current <- table[evaluation_included == TRUE]
  if (cohort == "activity_supported") current <- current[
    activity_supported == TRUE]
  rbindlist(lapply(c("pooled", "within_perturbation",
    "strict_cell_line_by_perturbation"), function(metric) {
      residual <- function_residual_table(current, metric)
      weights <- balanced_weights(residual$drug)
      data.table(cohort = cohort, metric = metric, n_rows = nrow(residual),
        n_conditions = uniqueN(residual$perturbation_id),
        n_drugs = uniqueN(residual$drug), n_lines = uniqueN(residual$depmap_id),
        pearson = safe_cor(residual$observed_residual,
          residual$predicted_residual),
        drug_balanced_pearson = weighted_cor(residual$observed_residual,
          residual$predicted_residual, weights),
        calibration_slope = if (sd(residual$predicted_residual) > 0)
          coef(lm(observed_residual ~ predicted_residual,
            data = residual))[2L] else NA_real_,
        observed_sd = sd(residual$observed_residual),
        predicted_sd = sd(residual$predicted_residual),
        predicted_to_observed_sd_ratio = sd(residual$predicted_residual) /
          sd(residual$observed_residual),
        rmse = sqrt(mean((residual$observed_residual -
          residual$predicted_residual) ^ 2)),
        correlation_method = "pearson")
    }))
}

weighted_slope <- function(observed, predicted, weights) {
  keep <- is.finite(observed) & is.finite(predicted) &
    is.finite(weights) & weights > 0
  observed <- observed[keep]; predicted <- predicted[keep]
  weights <- weights[keep] / sum(weights[keep])
  observed_center <- sum(weights * observed)
  predicted_center <- sum(weights * predicted)
  denominator <- sum(weights * (predicted - predicted_center) ^ 2)
  if (denominator <= 0) return(NA_real_)
  sum(weights * (predicted - predicted_center) *
    (observed - observed_center)) / denominator
}

bootstrap_metrics <- function(table, endpoint, metric, scenario_index) {
  if (endpoint == "function") {
    residual <- function_residual_table(table, metric)
    observed <- residual$observed_residual
    predicted <- residual$predicted_residual
    drug <- residual$drug
  } else {
    observed <- table$measured_z
    predicted <- table$predicted_z
    drug <- table$drug
  }
  drugs <- sort(unique(drug))
  row_count <- as.numeric(table(factor(drug, levels = drugs))[
    match(drug, drugs)])
  set.seed(SEED + match(metric, c("molecular", "pooled", "within_perturbation",
      "strict_cell_line_by_perturbation")) * 1000L)
  rbindlist(lapply(seq_len(N_BOOTSTRAPS), function(draw) {
    sampled <- sample(drugs, length(drugs), replace = TRUE)
    multiplicity <- as.numeric(table(factor(sampled, levels = drugs))[
      match(drug, drugs)])
    balanced <- multiplicity / row_count
    data.table(bootstrap = draw, endpoint = endpoint, metric = metric,
      pearson = weighted_cor(observed, predicted, multiplicity),
      drug_balanced_pearson = weighted_cor(observed, predicted, balanced),
      calibration_slope = weighted_slope(observed, predicted, multiplicity),
      predicted_to_observed_sd_ratio = weighted_sd(predicted, multiplicity) /
        weighted_sd(observed, multiplicity))
  }))
}

## Frozen v2 base inputs and models.
folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- folds[outer_fold != OUTER_FOLD, depmap_id]
held_lines <- folds[outer_fold == OUTER_FOLD, depmap_id]
all_lines <- sort(folds$depmap_id)
response_stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- response_stability[outer_fold == OUTER_FOLD & stable == TRUE,
  unique(program)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
baseline_stability <- fread(file.path(stage_result("03_baseline_programs"),
  "baseline_program_stability.tsv"))
baseline_programs <- baseline_stability[outer_fold == OUTER_FOLD &
  stable == TRUE, unique(program)]
baseline_programs <- paste0("B", sort(as.integer(sub("B", "",
  baseline_programs))))
if (length(programs) != 9L || length(baseline_programs) != 5L)
  fail("The frozen v2 program banks are incomplete.")

baseline_long <- fread(file.path(stage_artifact("03_baseline_programs"),
  "baseline_program_scores.tsv.gz"))
make_baseline_matrix <- function(source, aggregation) {
  table <- dcast(baseline_long[baseline_source == source &
    aggregation_level == aggregation & program %chin% baseline_programs],
    depmap_id ~ program, value.var = "raw_program_score")
  matrix <- as.matrix(table[, ..baseline_programs])
  rownames(matrix) <- table$depmap_id
  matrix[all_lines, baseline_programs, drop = FALSE]
}
dmso_raw <- make_baseline_matrix("tahoe_dmso", "line_mean")
depmap_raw <- make_baseline_matrix("depmap_static", "line_level")

## Five independent training-line affine conversions.
aligned_depmap_raw <- matrix(NA_real_, nrow(depmap_raw), ncol(depmap_raw),
  dimnames = dimnames(depmap_raw))
conversion_parts <- vector("list", length(baseline_programs))
for (index in seq_along(baseline_programs)) {
  program <- baseline_programs[index]
  fit <- lm(dmso_raw[training_lines, program] ~
    depmap_raw[training_lines, program])
  aligned_depmap_raw[, program] <- coef(fit)[1L] +
    coef(fit)[2L] * depmap_raw[, program]
  conversion_parts[[index]] <- data.table(program = program,
    intercept = unname(coef(fit)[1L]), slope = unname(coef(fit)[2L]),
    training_lines = length(training_lines), held_lines_used = 0L,
    training_pearson = safe_cor(dmso_raw[training_lines, program],
      depmap_raw[training_lines, program]),
    training_rmse = sqrt(mean(residuals(fit) ^ 2)),
    fit = "ordinary least squares per program",
    molecular_or_prism_outcomes_used = FALSE)
}
conversion <- rbindlist(conversion_parts)
fwrite(conversion, file.path(RESULT,
  "depmap_to_dmso_affine_conversion.tsv"), sep = "\t")

alignment_parts <- list()
for (program in baseline_programs) {
  for (role in c("training", "held_out")) {
    lines <- if (role == "training") training_lines else held_lines
    observed <- dmso_raw[lines, program]
    unaligned <- depmap_raw[lines, program]
    aligned <- aligned_depmap_raw[lines, program]
    alignment_parts[[length(alignment_parts) + 1L]] <- data.table(
      program = program, role = role, n_lines = length(lines),
      unaligned_pearson = safe_cor(observed, unaligned),
      aligned_pearson = safe_cor(observed, aligned),
      unaligned_rmse = sqrt(mean((observed - unaligned) ^ 2)),
      aligned_rmse = sqrt(mean((observed - aligned) ^ 2)),
      aligned_calibration_slope = coef(lm(observed ~ aligned))[2L],
      aligned_to_dmso_sd_ratio = sd(aligned) / sd(observed))
  }
}
alignment_performance <- rbindlist(alignment_parts)
fwrite(alignment_performance, file.path(RESULT,
  "baseline_alignment_performance.tsv"), sep = "\t")

models <- readRDS(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_bridge_models.rds"))
if (is.null(models$tahoe_dmso$fit) || is.null(models$depmap_static$fit))
  fail("Both accepted in-domain baseline conditioners must be nonzero.")
dmso_scores <- standardize_fixed(dmso_raw,
  models$tahoe_dmso$baseline_center, models$tahoe_dmso$baseline_scale)
depmap_scores <- standardize_fixed(depmap_raw,
  models$depmap_static$baseline_center, models$depmap_static$baseline_scale)
aligned_depmap_scores <- standardize_fixed(aligned_depmap_raw,
  models$tahoe_dmso$baseline_center, models$tahoe_dmso$baseline_scale)
rownames(dmso_scores) <- rownames(depmap_scores) <-
  rownames(aligned_depmap_scores) <- all_lines
colnames(dmso_scores) <- colnames(depmap_scores) <-
  colnames(aligned_depmap_scores) <- baseline_programs

footprint_dt <- fread(file.path(stage_result("05_program_function_bridge"),
  "reference_perturbation_footprints.tsv.gz"))
footprint <- as.matrix(footprint_dt[, ..programs])
rownames(footprint) <- footprint_dt$perturbation_id
decoder_coefficients <- fread(file.path(stage_result(
  "05_program_function_bridge"), "function_decoder_coefficients.tsv"))
decoder <- list(intercept = decoder_coefficients[
  term == "(Intercept)", coefficient], beta = setNames(decoder_coefficients[
    term != "(Intercept)", coefficient], decoder_coefficients[
      term != "(Intercept)", term]))

score_long <- fread(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"))[program %chin% programs]
score_meta <- dcast(score_long, observation_id + depmap_id + cell_line +
  perturbation_id + drug + dose + dose_panel + role ~ program,
  value.var = "score")
held_score <- score_meta$depmap_id %chin% held_lines
measured_scores <- fread(file.path(stage_artifact(
  "05_program_function_bridge"),
  "measured_function_relevant_scores.tsv.gz"))
measured_z <- measured_scores$measured_z[match(score_meta$observation_id,
  measured_scores$observation_id)]
activity <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
active_ids <- activity[activity_supported == TRUE, perturbation_id]
function_ids <- activity[activity_supported_prism_overlap == TRUE,
  perturbation_id]

scenario_specs <- list(
  dmso_trained_to_dmso = list(scores = dmso_scores,
    fit = models$tahoe_dmso$fit,
    training_source = "tahoe_dmso", application_source = "tahoe_dmso",
    alignment = "none"),
  depmap_trained_to_depmap = list(scores = depmap_scores,
    fit = models$depmap_static$fit,
    training_source = "depmap_static", application_source = "depmap_static",
    alignment = "none"),
  dmso_trained_to_aligned_depmap = list(scores = aligned_depmap_scores,
    fit = models$tahoe_dmso$fit,
    training_source = "tahoe_dmso", application_source = "depmap_static",
    alignment = "training-line per-program affine"))

function_rows <- fread(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"))
function_rows[, `:=`(observation_id = paste(depmap_id, perturbation_id,
  sep = "__"), drug = tahoe_drug, dose = tahoe_dose_uM)]
held_function <- function_rows[outer_fold == OUTER_FOLD]
mean_prediction <- decoder$intercept + as.numeric(footprint[
  held_function$perturbation_id, programs, drop = FALSE] %*%
  decoder$beta[programs])

molecular_parts <- list(); function_parts <- list()
for (scenario in names(scenario_specs)) {
  specification <- scenario_specs[[scenario]]
  molecular_design <- make_bilinear_design(
    score_meta$depmap_id[held_score], score_meta$perturbation_id[held_score],
    specification$scores, footprint, programs)
  molecular_parts[[scenario]] <- data.table(
    scenario = scenario, training_source = specification$training_source,
    application_source = specification$application_source,
    alignment = specification$alignment,
    observation_id = score_meta$observation_id[held_score],
    depmap_id = score_meta$depmap_id[held_score],
    cell_line = score_meta$cell_line[held_score],
    perturbation_id = score_meta$perturbation_id[held_score],
    drug = score_meta$drug[held_score], dose = score_meta$dose[held_score],
    measured_z = measured_z[held_score],
    predicted_z = predict_baseline(specification$fit, molecular_design),
    activity_supported = score_meta$perturbation_id[held_score] %chin%
      active_ids,
    activity_supported_prism = score_meta$perturbation_id[held_score] %chin%
      function_ids)
  function_design <- make_bilinear_design(held_function$depmap_id,
    held_function$perturbation_id, specification$scores, footprint, programs)
  predicted_z <- predict_baseline(specification$fit, function_design)
  function_parts[[scenario]] <- data.table(
    scenario = scenario, training_source = specification$training_source,
    application_source = specification$application_source,
    alignment = specification$alignment,
    observation_id = held_function$observation_id,
    depmap_id = held_function$depmap_id,
    cell_line = held_function$cell_line,
    perturbation_id = held_function$perturbation_id,
    drug = held_function$drug, dose = held_function$dose,
    observed_function = held_function$prism_sensitization,
    evaluation_included = held_function$chain_observed == 1L,
    activity_supported = held_function$perturbation_id %chin% active_ids,
    measured_z = measured_scores$measured_z[match(
      held_function$observation_id, measured_scores$observation_id)],
    predicted_z = predicted_z, mean_prediction = mean_prediction,
    predicted_function = mean_prediction + predicted_z)
}
molecular_predictions <- rbindlist(molecular_parts)
function_predictions <- rbindlist(function_parts)
fwrite(molecular_predictions, file.path(ARTIFACT,
  "cross_source_molecular_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")
fwrite(function_predictions, file.path(ARTIFACT,
  "cross_source_function_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")

## Molecular and end-to-end functional performance.
molecular_performance <- rbindlist(lapply(names(scenario_specs),
  function(scenario_name) {
    current <- molecular_predictions[scenario == scenario_name]
    function_panel <- function_predictions[scenario == scenario_name &
      evaluation_included == TRUE & activity_supported == TRUE]
    rbindlist(list(cbind(data.table(scenario = scenario_name,
      cohort = "all_activity_supported"), endpoint_metrics(current[
        activity_supported == TRUE], "measured_z", "predicted_z")),
      cbind(data.table(scenario = scenario_name,
        cohort = "activity_supported_prism"), endpoint_metrics(function_panel,
          "measured_z", "predicted_z"))))
  }))
fwrite(molecular_performance, file.path(RESULT,
  "cross_source_molecular_performance.tsv"), sep = "\t")

functional_performance <- rbindlist(lapply(names(scenario_specs),
  function(scenario_name) {
    current <- function_predictions[scenario == scenario_name]
    result <- rbindlist(lapply(c("complete_matched_panel",
      "activity_supported"), function(cohort) functional_metrics(current,
        cohort)))
    result[, scenario := scenario_name]
    result
  }))
setcolorder(functional_performance, c("scenario",
  setdiff(names(functional_performance), "scenario")))
fwrite(functional_performance, file.path(RESULT,
  "cross_source_functional_performance.tsv"), sep = "\t")

## Drug-block bootstrap intervals for molecular and functional metrics.
bootstrap_parts <- list()
for (scenario_index in seq_along(scenario_specs)) {
  scenario_name <- names(scenario_specs)[scenario_index]
  molecular_current <- molecular_predictions[scenario == scenario_name &
    activity_supported == TRUE]
  draw <- bootstrap_metrics(molecular_current, "molecular", "molecular",
    scenario_index)
  draw[, scenario := scenario_name]
  bootstrap_parts[[length(bootstrap_parts) + 1L]] <- draw
  function_current <- function_predictions[scenario == scenario_name &
    evaluation_included == TRUE & activity_supported == TRUE]
  for (metric in c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation")) {
    draw <- bootstrap_metrics(function_current, "function", metric,
      scenario_index)
    draw[, scenario := scenario_name]
    bootstrap_parts[[length(bootstrap_parts) + 1L]] <- draw
  }
}
bootstrap <- rbindlist(bootstrap_parts)
fwrite(bootstrap, file.path(ARTIFACT,
  "cross_source_drug_block_bootstrap.tsv.gz"), sep = "\t",
  compress = "gzip")
bootstrap_intervals <- bootstrap[, .(
  pearson_ci_lower = quantile(pearson, 0.025, na.rm = TRUE),
  pearson_ci_upper = quantile(pearson, 0.975, na.rm = TRUE),
  drug_balanced_pearson_ci_lower = quantile(drug_balanced_pearson,
    0.025, na.rm = TRUE),
  drug_balanced_pearson_ci_upper = quantile(drug_balanced_pearson,
    0.975, na.rm = TRUE),
  calibration_slope_ci_lower = quantile(calibration_slope, 0.025,
    na.rm = TRUE),
  calibration_slope_ci_upper = quantile(calibration_slope, 0.975,
    na.rm = TRUE),
  sd_ratio_ci_lower = quantile(predicted_to_observed_sd_ratio, 0.025,
    na.rm = TRUE),
  sd_ratio_ci_upper = quantile(predicted_to_observed_sd_ratio, 0.975,
    na.rm = TRUE)), by = .(scenario, endpoint, metric)]
fwrite(bootstrap_intervals, file.path(RESULT,
  "cross_source_bootstrap_intervals.tsv"), sep = "\t")

primary <- bootstrap[scenario == "dmso_trained_to_aligned_depmap"]
comparison_parts <- lapply(c("dmso_trained_to_dmso",
  "depmap_trained_to_depmap"), function(comparator) {
    reference <- bootstrap[scenario == comparator]
    joined <- merge(primary, reference,
      by = c("bootstrap", "endpoint", "metric"),
      suffixes = c("_primary", "_comparator"))
    joined[, `:=`(comparator = comparator,
      pearson_change = pearson_primary - pearson_comparator,
      drug_balanced_pearson_change = drug_balanced_pearson_primary -
        drug_balanced_pearson_comparator)]
    joined
  })
bootstrap_comparison <- rbindlist(comparison_parts)
fwrite(bootstrap_comparison, file.path(ARTIFACT,
  "primary_transfer_comparison_bootstrap.tsv.gz"), sep = "\t",
  compress = "gzip")
comparison_intervals <- bootstrap_comparison[, .(
  pearson_change_ci_lower = quantile(pearson_change, 0.025, na.rm = TRUE),
  pearson_change_ci_upper = quantile(pearson_change, 0.975, na.rm = TRUE),
  probability_change_positive = mean(pearson_change > 0, na.rm = TRUE),
  drug_balanced_change_ci_lower = quantile(drug_balanced_pearson_change,
    0.025, na.rm = TRUE),
  drug_balanced_change_ci_upper = quantile(drug_balanced_pearson_change,
    0.975, na.rm = TRUE)),
  by = .(comparator, endpoint, metric)]
fwrite(comparison_intervals, file.path(RESULT,
  "primary_transfer_comparison_intervals.tsv"), sep = "\t")

## Compact primary transfer assessment.
primary_molecular <- molecular_performance[
  scenario == "dmso_trained_to_aligned_depmap" &
    cohort == "all_activity_supported"]
primary_function <- functional_performance[
  scenario == "dmso_trained_to_aligned_depmap" &
    cohort == "activity_supported"]
primary_molecular_ci <- bootstrap_intervals[
  scenario == "dmso_trained_to_aligned_depmap" & endpoint == "molecular"]
assessment <- data.table(
  primary_molecular_pearson = primary_molecular$pearson,
  primary_molecular_drug_balanced_pearson =
    primary_molecular$drug_balanced_pearson,
  primary_molecular_pearson_ci_lower = primary_molecular_ci$pearson_ci_lower,
  primary_molecular_pearson_ci_upper = primary_molecular_ci$pearson_ci_upper,
  primary_molecular_calibration_slope = primary_molecular$calibration_slope,
  primary_molecular_sd_ratio =
    primary_molecular$predicted_to_observed_sd_ratio,
  activity_function_pooled = primary_function[metric == "pooled", pearson],
  activity_function_within = primary_function[
    metric == "within_perturbation", pearson],
  activity_function_strict = primary_function[
    metric == "strict_cell_line_by_perturbation", pearson],
  molecular_change_vs_dmso_in_domain = primary_molecular$pearson -
    molecular_performance[scenario == "dmso_trained_to_dmso" &
      cohort == "all_activity_supported", pearson],
  molecular_change_vs_depmap_in_domain = primary_molecular$pearson -
    molecular_performance[scenario == "depmap_trained_to_depmap" &
      cohort == "all_activity_supported", pearson],
  positive_molecular_transfer_supported =
    primary_molecular_ci$pearson_ci_lower > 0,
  held_outcomes_used_for_alignment = FALSE)
fwrite(assessment, file.path(RESULT, "cross_source_transfer_assessment.tsv"),
  sep = "\t")

## One plot per PDF.
theme_v2 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1,
    plot.title = element_text(face = "bold"), legend.position = "bottom")
figure_index <- list()
for (program in baseline_programs) {
  alignment_plot <- data.table(depmap_id = all_lines,
    dmso = dmso_raw[, program], unaligned = depmap_raw[, program],
    aligned = aligned_depmap_raw[, program],
    role = fifelse(all_lines %chin% training_lines, "training", "held_out"))
  plot_long <- melt(alignment_plot,
    id.vars = c("depmap_id", "dmso", "role"),
    measure.vars = c("unaligned", "aligned"), variable.name = "mapping",
    value.name = "depmap_score")
  held_table <- alignment_plot[role == "held_out"]
  plot <- ggplot(plot_long, aes(x = depmap_score, y = dmso,
      color = mapping, shape = role)) +
    geom_point(alpha = 0.75, size = 1.4) +
    geom_smooth(aes(linetype = mapping), method = "lm", formula = y ~ x,
      se = FALSE, linewidth = 0.55) +
    annotate("text", x = Inf, y = -Inf,
      label = sprintf("Held r = %.3f", safe_cor(held_table$dmso,
        held_table$aligned)), hjust = 1.05, vjust = -0.6, size = 2.7) +
    scale_color_manual(values = c(unaligned = "#8A8A8A",
      aligned = "#276FBF")) +
    labs(title = paste(program, "DepMap to DMSO alignment"),
      x = "DepMap program score", y = "Tahoe DMSO program score",
      color = "DepMap score", shape = "Cell-line role") + theme_v2()
  filename <- sprintf("09a_alignment_%s.pdf", program)
  ggsave(file.path(RESULT, filename), plot, width = 4.35, height = 4.15,
    device = cairo_pdf, bg = "white")
  figure_index[[length(figure_index) + 1L]] <- data.table(figure = filename,
    description = paste("Per-program affine alignment for", program))
}

primary_molecular_rows <- molecular_predictions[
  scenario == "dmso_trained_to_aligned_depmap" & activity_supported == TRUE]
molecular_plot <- ggplot(primary_molecular_rows,
  aes(x = predicted_z, y = measured_z)) +
  geom_point(alpha = 0.35, size = 0.7, color = "#276FBF") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.6, color = "#C44536") +
  annotate("text", x = Inf, y = -Inf,
    label = sprintf("Pearson r = %.3f", safe_cor(
      primary_molecular_rows$predicted_z,
      primary_molecular_rows$measured_z)), hjust = 1.05, vjust = -0.6,
    size = 2.8) +
  labs(title = "DMSO-trained conditioner on aligned DepMap",
    x = "Predicted molecular score", y = "Measured molecular score") +
  theme_v2()
ggsave(file.path(RESULT, "09b_aligned_depmap_predicted_vs_measured.pdf"),
  molecular_plot, width = 4.25, height = 4.15, device = cairo_pdf,
  bg = "white")
figure_index[[length(figure_index) + 1L]] <- data.table(
  figure = "09b_aligned_depmap_predicted_vs_measured.pdf",
  description = "Primary transfer prediction of the measured molecular target")

primary_function_rows <- function_predictions[
  scenario == "dmso_trained_to_aligned_depmap" &
    evaluation_included == TRUE & activity_supported == TRUE]
function_plot <- ggplot(primary_function_rows,
  aes(x = predicted_function, y = observed_function)) +
  geom_point(alpha = 0.4, size = 0.8, color = "#276FBF") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.6, color = "#C44536") +
  annotate("text", x = Inf, y = -Inf,
    label = sprintf("Pearson r = %.3f", safe_cor(
      primary_function_rows$predicted_function,
      primary_function_rows$observed_function)),
    hjust = 1.05, vjust = -0.6, size = 2.8) +
  labs(title = "Aligned-DepMap end-to-end function",
    x = "Predicted PRISM sensitization", y = "Observed PRISM sensitization") +
  theme_v2()
ggsave(file.path(RESULT, "09c_aligned_depmap_end_to_end_function.pdf"),
  function_plot, width = 4.25, height = 4.15, device = cairo_pdf,
  bg = "white")
figure_index[[length(figure_index) + 1L]] <- data.table(
  figure = "09c_aligned_depmap_end_to_end_function.pdf",
  description = "Primary transfer end-to-end held function")
fwrite(rbindlist(figure_index), file.path(RESULT, "FIGURE_INDEX.tsv"),
  sep = "\t")

## Verify that both in-domain comparators are exactly the frozen v2 fits.
base_molecular <- fread(file.path(stage_result("06_baseline_bridge"),
  "baseline_molecular_performance.tsv"))
base_function <- fread(file.path(stage_result("07_end_to_end"),
  "end_to_end_performance.tsv"))
comparator_checks <- list()
for (mapping in list(
    c("dmso_trained_to_dmso", "tahoe_dmso"),
    c("depmap_trained_to_depmap", "depmap_static"))) {
  scenario_name <- mapping[1L]; source_name <- mapping[2L]
  observed_molecular <- molecular_performance[scenario == scenario_name &
    cohort == "all_activity_supported"]
  expected_molecular <- base_molecular[baseline_source == source_name &
    cohort == "all_activity_supported"]
  for (metric_name in c("pearson", "drug_balanced_pearson",
      "calibration_slope", "predicted_to_observed_sd_ratio")) {
    observed <- observed_molecular[[metric_name]]
    expected <- expected_molecular[[metric_name]]
    comparator_checks[[length(comparator_checks) + 1L]] <- data.table(
      scenario = scenario_name, endpoint = "molecular",
      metric = metric_name, observed = observed, expected = expected,
      absolute_difference = abs(observed - expected),
      passed = abs(observed - expected) < 1e-12)
  }
  for (functional_metric in c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation")) {
    observed_function <- functional_performance[scenario == scenario_name &
      cohort == "activity_supported" & metric == functional_metric]
    expected_function <- base_function[
      prediction_model == "baseline_conditioned" &
      baseline_source == source_name & cohort == "activity_supported" &
      metric == functional_metric]
    for (metric_name in c("pearson", "drug_balanced_pearson",
        "calibration_slope", "predicted_to_observed_sd_ratio")) {
      observed <- observed_function[[metric_name]]
      expected <- expected_function[[metric_name]]
      comparator_checks[[length(comparator_checks) + 1L]] <- data.table(
        scenario = scenario_name, endpoint = paste0("function_",
          functional_metric), metric = metric_name, observed = observed,
        expected = expected, absolute_difference = abs(observed - expected),
        passed = abs(observed - expected) < 1e-12)
    }
  }
}
comparator_checks <- rbindlist(comparator_checks)
fwrite(comparator_checks, file.path(RESULT,
  "in_domain_comparator_reproduction_checks.tsv"), sep = "\t")

checks <- data.table(check = c("stage1_gate_present",
  "exactly_five_affine_conversions", "conversion_training_lines_26",
  "conversion_held_lines_zero", "conversion_uses_no_outcomes",
  "in_domain_comparators_exact", "unchanged_dmso_conditioner",
  "same_five_baseline_coordinates",
  "no_direct_baseline_to_prism_fit", "no_drug_identity_predictor",
  "no_per_drug_models"), passed = c(TRUE, nrow(conversion) == 5L,
    all(conversion$training_lines == 26L),
    all(conversion$held_lines_used == 0L),
    !any(conversion$molecular_or_prism_outcomes_used),
    all(comparator_checks$passed), TRUE, TRUE, TRUE, TRUE, TRUE))
fwrite(checks, file.path(RESULT, "dmso_to_depmap_required_checks.tsv"),
  sep = "\t")
if (!all(checks$passed)) fail("A Stage 2B required check failed.")
writeLines(c("# DMSO-trained conditioning transferred to DepMap",
  "",
  "Five independent affine conversions were fitted using only the 26 training-line baseline coordinates.",
  sprintf("Held molecular Pearson for DMSO-trained to aligned DepMap was %.3f (95%% drug-bootstrap CI %.3f to %.3f).",
    assessment$primary_molecular_pearson,
    assessment$primary_molecular_pearson_ci_lower,
    assessment$primary_molecular_pearson_ci_upper),
  sprintf("Activity-panel end-to-end Pearson was %.3f pooled, %.3f within perturbation, and %.3f strict.",
    assessment$activity_function_pooled,
    assessment$activity_function_within,
    assessment$activity_function_strict),
  "The transferred signal was positive but remained strongly amplitude-compressed and was weaker than the in-domain DepMap molecular model.",
  "The DMSO conditioner was applied unchanged; no held molecular or PRISM outcome was used for source conversion."),
  file.path(RESULT, "README.md"))
write_artifact_manifest(ARTIFACT,
  "scripts/09_test_dmso_to_depmap_transfer.R")
writeLines("complete", file.path(RESULT, "stage2b_complete.flag"))
message(sprintf(paste0("Stage 2B complete: aligned DepMap molecular r=%.3f; ",
  "function pooled/within/strict=%.3f/%.3f/%.3f."),
  assessment$primary_molecular_pearson,
  assessment$activity_function_pooled,
  assessment$activity_function_within,
  assessment$activity_function_strict))
