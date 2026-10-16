#!/usr/bin/env Rscript

## Frozen end-to-end evaluation of the response-informed baseline
## susceptibility programs.  This stage applies the accepted Stage 05 decoder
## without fitting, recalibration, or model selection.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("17_baseline_susceptibility_chain")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(RESULT, c("stage17_complete.flag",
  "17e_per_drug_improvement_distribution.pdf",
  "17f_ntc_vs_depmap_end_to_end_concordance.pdf",
  "baseline_susceptibility_chain_molecular_comparison.tsv",
  "baseline_susceptibility_chain_performance.tsv",
  "baseline_susceptibility_chain_per_line.tsv",
  "baseline_susceptibility_chain_per_perturbation.tsv")), force = TRUE)
theme_stage17 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage17(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(name, plot, width, height) ggsave(
  file.path(RESULT, name), plot, width = width, height = height,
  device = cairo_pdf, bg = "white", limitsize = FALSE)
fmt <- function(x, digits = 3L) ifelse(is.finite(x),
  formatC(x, digits = digits, format = "f"), "NA")

required <- c(
  file.path(stage_result("16_baseline_susceptibility_programs"),
    "stage16_complete.flag"),
  file.path(stage_result("16_baseline_susceptibility_programs"),
    "held_susceptibility_molecular_predictions.tsv.gz"),
  file.path(stage_result("16_baseline_susceptibility_programs"),
    "stage16_required_checks.tsv"),
  file.path(stage_result("14_continuous_baseline_susceptibility"),
    "held_function_relevant_molecular_predictions.tsv.gz"),
  file.path(stage_result("05_program_function_bridge"),
    "function_decoder_coefficients.tsv"),
  file.path(stage_artifact("05_program_function_bridge"),
    "measured_function_predictions.tsv.gz"),
  file.path(stage_artifact("06_baseline_bridge"),
    "baseline_function_predictions.tsv.gz"),
  file.path(stage_artifact("06_baseline_bridge"),
    "baseline_molecular_predictions.tsv.gz"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 17 input(s): %s",
  paste(missing, collapse = ", "))

model_labels <- c(
  existing_baseline_bridge = "Existing baseline bridge",
  current_S2_S8_continuous = "Current S2-S8 continuous",
  previous_coexpression_susceptibility_programs =
    "Previous coexpression susceptibility programs",
  new_NTC_susceptibility_programs =
    "Direct-axis O1-O3 NTC modules",
  new_DepMap_transfer_susceptibility_programs =
    "Direct-axis O1-O3 DepMap transfer",
  zero_context_reference_footprint = "Zero context / reference",
  complete_measured_molecular_z_ceiling = "Complete measured ceiling")
models <- names(model_labels)
row_models <- setdiff(models,
  "previous_coexpression_susceptibility_programs")

decoder <- fread(file.path(stage_result("05_program_function_bridge"),
  "function_decoder_coefficients.tsv"))
programs <- decoder[term != "(Intercept)", term]
programs <- programs[grepl("^P[0-9]+$", programs)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS)
  fail("The frozen decoder does not contain exactly nine response programs.")
intercept <- decoder[term == "(Intercept)", coefficient]
beta <- setNames(decoder[term %chin% programs, coefficient],
  decoder[term %chin% programs, term])
if (length(intercept) != 1L || length(beta) != length(programs))
  fail("The frozen Stage 05 decoder is incomplete.")

## Align the executable molecular comparators on identical held rows; the
## failed coexpression method is appended from its frozen compact summary.
stage16 <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"),
  "held_susceptibility_molecular_predictions.tsv.gz"))
ntc <- stage16[baseline_source == "tahoe_ntc"]
depmap <- stage16[baseline_source == "depmap_static"]
ntc[, predicted_z := as.numeric(predicted_z)]
depmap[, predicted_z := as.numeric(predicted_z)]
if (!setequal(ntc$observation_id, depmap$observation_id))
  fail("NTC and DepMap Stage 16 predictions do not cover identical rows.")
setkey(ntc, observation_id); setkey(depmap, observation_id)

continuous <- fread(file.path(stage_result(
  "14_continuous_baseline_susceptibility"),
  "held_function_relevant_molecular_predictions.tsv.gz"))
baseline_molecular <- fread(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_molecular_predictions.tsv.gz"))[baseline_source == "tahoe_dmso"]
setkey(continuous, observation_id); setkey(baseline_molecular, observation_id)

molecular <- ntc[, c("observation_id", "depmap_id", "cell_line",
  "perturbation_id", "drug", "dose", "dose_panel", "molecular_eligible",
  "activity_supported", "observed_z", paste0("reference_", programs),
  paste0("observed_context_", programs)), with = FALSE]
molecular[, existing_baseline_bridge :=
  baseline_molecular[.SD, on = "observation_id", x.predicted_z]]
molecular[, current_S2_S8_continuous :=
  continuous[.SD, on = "observation_id", x.predicted_z]]
molecular[, new_NTC_susceptibility_programs :=
  ntc[.SD, on = "observation_id", x.predicted_z]]
molecular[, new_DepMap_transfer_susceptibility_programs :=
  depmap[.SD, on = "observation_id", x.predicted_z]]
molecular[, zero_context_reference_footprint := 0]
molecular[, complete_measured_molecular_z_ceiling := observed_z]

accepted <- fread(file.path(stage_artifact("05_program_function_bridge"),
  "measured_function_predictions.tsv.gz"))
matched_activity_ids <- accepted[evaluation_included == TRUE &
  activity_supported == TRUE, observation_id]
molecular[, matched_activity_panel := observation_id %chin%
  matched_activity_ids]

alignment_errors <- c(
  baseline = max(abs(molecular$observed_z - baseline_molecular[
    molecular, on = "observation_id", x.measured_z]), na.rm = TRUE),
  continuous = max(abs(molecular$observed_z - continuous[
    molecular, on = "observation_id", x.observed_z]), na.rm = TRUE),
  depmap = max(abs(molecular$observed_z - depmap[
    molecular, on = "observation_id", x.observed_z]), na.rm = TRUE))
if (any(!is.finite(alignment_errors)) || max(alignment_errors) > 1e-9)
  fail("The molecular comparator targets are not aligned.")

molecular_long <- melt(molecular,
  id.vars = setdiff(names(molecular), row_models), measure.vars = row_models,
  variable.name = "prediction_model", value.name = "predicted_z",
  variable.factor = FALSE)

metric_values <- function(table, observed, predicted) {
  keep <- is.finite(table[[observed]]) & is.finite(table[[predicted]])
  current <- table[keep]
  if (!nrow(current)) return(data.table(n_rows = 0L, n_lines = 0L,
    n_conditions = 0L, n_drugs = 0L, pearson = NA_real_,
    drug_balanced_pearson = NA_real_, observed_sd = NA_real_,
    predicted_sd = NA_real_, predicted_to_observed_sd_ratio = NA_real_,
    mse = NA_real_, rmse = NA_real_, correlation_method = "pearson",
    reconstructed_molecular_displacement_pearson = NA_real_,
    reconstructed_context_displacement_pearson = NA_real_))
  y <- current[[observed]]; p <- current[[predicted]]
  w <- balanced_weights(current$drug)
  error <- y - p
  result <- data.table(n_rows = nrow(current),
    n_lines = uniqueN(current$depmap_id),
    n_conditions = uniqueN(current$perturbation_id),
    n_drugs = uniqueN(current$drug), pearson = safe_cor(y, p),
    drug_balanced_pearson = weighted_cor(y, p, w),
    observed_sd = sd(y), predicted_sd = sd(p),
    predicted_to_observed_sd_ratio = sd(p) / sd(y),
    mse = mean(error ^ 2), rmse = sqrt(mean(error ^ 2)),
    correlation_method = "pearson")
  reference_columns <- paste0("reference_", programs)
  context_columns <- paste0("observed_context_", programs)
  if (all(c(reference_columns, context_columns) %chin% names(current))) {
    reference <- as.matrix(current[, ..reference_columns])
    observed_context <- as.matrix(current[, ..context_columns])
    predicted_context <- outer(p / sum(beta ^ 2), beta[programs])
    result[, `:=`(reconstructed_molecular_displacement_pearson = safe_cor(
        as.numeric(reference + observed_context),
        as.numeric(reference + predicted_context)),
      reconstructed_context_displacement_pearson = safe_cor(
        as.numeric(observed_context), as.numeric(predicted_context)))]
  } else result[, `:=`(
    reconstructed_molecular_displacement_pearson = NA_real_,
    reconstructed_context_displacement_pearson = NA_real_)]
  result
}

molecular_subset <- function(table, cohort) {
  if (cohort == "held_all_molecular_rows") return(table)
  if (cohort == "held_activity_supported")
    return(table[activity_supported == TRUE])
  if (cohort == "held_molecularly_eligible")
    return(table[molecular_eligible == TRUE])
  if (cohort == "held_eligible_and_activity_supported")
    return(table[molecular_eligible == TRUE & activity_supported == TRUE])
  if (cohort == "held_matched_PRISM_activity_panel")
    return(table[matched_activity_panel == TRUE])
  panel <- as.integer(sub("held_molecularly_eligible_dose_panel_", "", cohort))
  if (is.finite(panel)) return(table[molecular_eligible == TRUE &
    dose_panel == panel])
  fail("Unknown molecular cohort: %s", cohort)
}

molecular_cohorts <- c("held_all_molecular_rows",
  "held_activity_supported", "held_molecularly_eligible",
  "held_eligible_and_activity_supported",
  "held_matched_PRISM_activity_panel",
  paste0("held_molecularly_eligible_dose_panel_", 1:3))
molecular_comparison <- rbindlist(lapply(molecular_cohorts, function(cohort)
  rbindlist(lapply(row_models, function(model) {
    current <- molecular_subset(
      molecular_long[prediction_model == model], cohort)
    cbind(data.table(endpoint = "molecular_z", cohort = cohort,
      metric = "pooled", prediction_model = model),
      metric_values(current, "observed_z", "predicted_z"))
  }))))

## Versioned compact record of the failed coexpression-grouped Stage 16 model.
## It remains a comparator, not an executable branch of the pipeline.
previous_molecular <- data.table(endpoint = "molecular_z",
  cohort = c("held_molecularly_eligible",
    "held_matched_PRISM_activity_panel"), metric = "pooled",
  prediction_model = "previous_coexpression_susceptibility_programs",
  n_rows = c(367L, 646L), n_lines = 7L, n_conditions = c(54L, 95L),
  n_drugs = c(42L, 54L), pearson = c(0.291059577068679,
    0.167835558717202), drug_balanced_pearson = c(0.347209910018229,
    0.165827879732250), observed_sd = c(0.442917204146064,
    0.344417357229911), predicted_sd = c(0.151689745244577,
    0.107585755374168), predicted_to_observed_sd_ratio = c(
    0.342478783449002, 0.312370306303553), mse = c(0.180673411358165,
    0.117807347720165), rmse = c(0.425056950723271,
    0.343230749963003), correlation_method = "pearson",
  reconstructed_molecular_displacement_pearson = NA_real_,
  reconstructed_context_displacement_pearson = NA_real_)
molecular_comparison <- rbindlist(list(molecular_comparison,
  previous_molecular), use.names = TRUE, fill = TRUE)
fwrite(molecular_comparison, file.path(RESULT,
  "molecular_performance.tsv"), sep = "\t")

## Decode Stage 14 and Stage 16 reconstructed displacements using the frozen
## Stage 05 coefficients.  No function outcome is used in this composition.
decode_wide <- function(table, prefix) {
  columns <- paste0(prefix, programs)
  if (!all(columns %chin% names(table)))
    fail("Reconstructed displacement columns are incomplete: %s", prefix)
  matrix <- as.matrix(table[, ..columns]); storage.mode(matrix) <- "double"
  data.table(observation_id = table$observation_id,
    predicted_function = intercept + as.numeric(matrix %*% beta[programs]))
}
decode_reference <- function(table) {
  columns <- paste0("reference_", programs)
  matrix <- as.matrix(table[, ..columns]); storage.mode(matrix) <- "double"
  intercept + as.numeric(matrix %*% beta[programs])
}

ntc_function <- decode_wide(ntc, "reconstructed_displacement_")
depmap_function <- decode_wide(depmap, "reconstructed_displacement_")
continuous_function <- decode_wide(continuous,
  "reconstructed_predicted_displacement_")
ntc_reference <- decode_reference(ntc)
depmap_reference <- decode_reference(depmap)
continuous_reference <- decode_reference(continuous)
ntc_composition_error <- max(abs(ntc_function$predicted_function -
  (ntc_reference + ntc$predicted_z)))
depmap_composition_error <- max(abs(depmap_function$predicted_function -
  (depmap_reference + depmap$predicted_z)))
continuous_composition_error <- max(abs(
  continuous_function$predicted_function -
  (continuous_reference + continuous$predicted_z)))
baseline_function <- fread(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_function_predictions.tsv.gz"))[baseline_source == "tahoe_dmso"]
setkey(baseline_function, observation_id)
setkey(ntc_function, observation_id); setkey(depmap_function, observation_id)
setkey(continuous_function, observation_id)

function_base <- accepted[, .(observation_id, depmap_id, cell_line,
  perturbation_id, drug, dose, observed_function,
  training_reference_function, evaluation_included, activity_supported,
  mean_prediction, measured_prediction, measured_z)]
function_base[, existing_baseline_bridge := baseline_function[
  .SD, on = "observation_id", x.predicted_function]]
function_base[, current_S2_S8_continuous := continuous_function[
  .SD, on = "observation_id", x.predicted_function]]
function_base[, new_NTC_susceptibility_programs := ntc_function[
  .SD, on = "observation_id", x.predicted_function]]
function_base[, new_DepMap_transfer_susceptibility_programs := depmap_function[
  .SD, on = "observation_id", x.predicted_function]]
function_base[, zero_context_reference_footprint := mean_prediction]
function_base[, complete_measured_molecular_z_ceiling := measured_prediction]
function_long <- melt(function_base,
  id.vars = setdiff(names(function_base), row_models),
  measure.vars = row_models,
  variable.name = "prediction_model", value.name = "predicted_function",
  variable.factor = FALSE)

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
    data.table(endpoint = "end_to_end_function", cohort = cohort,
      metric = metric, prediction_model = model, n_rows = nrow(residual),
      n_lines = uniqueN(residual$depmap_id),
      n_conditions = uniqueN(residual$perturbation_id),
      n_drugs = uniqueN(residual$drug),
      pearson = safe_cor(residual$observed_residual,
        residual$predicted_residual),
      drug_balanced_pearson = weighted_cor(residual$observed_residual,
        residual$predicted_residual, w),
      observed_sd = sd(residual$observed_residual),
      predicted_sd = sd(residual$predicted_residual),
      predicted_to_observed_sd_ratio = sd(residual$predicted_residual) /
        sd(residual$observed_residual), mse = mean(error ^ 2),
      rmse = sqrt(mean(error ^ 2)), correlation_method = "pearson")
  }))
}

performance <- rbindlist(lapply(c("complete_matched_panel",
    "activity_supported"), function(cohort)
  rbindlist(lapply(row_models, function(model)
    evaluate_function(function_long, model, cohort)))))
previous_function <- data.table(endpoint = "end_to_end_function",
  cohort = rep(c("complete_matched_panel", "activity_supported"), each = 3L),
  metric = rep(c("pooled", "within_perturbation",
    "strict_cell_line_by_perturbation"), 2L),
  prediction_model = "previous_coexpression_susceptibility_programs",
  n_rows = rep(c(1752L, 646L), each = 3L), n_lines = 7L,
  n_conditions = rep(c(262L, 95L), each = 3L),
  n_drugs = rep(c(88L, 54L), each = 3L),
  pearson = c(0.481627076598682, 0.0995105170139224,
    0.0639526125272674, 0.550031131248659, 0.141792510378784,
    0.0586937136332477),
  drug_balanced_pearson = c(0.474593106554096, 0.0991845181503362,
    0.0650292954728914, 0.511011424497038, 0.127195195234123,
    0.0341753185399496),
  observed_sd = c(1.3539797277128, 0.662849615669154,
    0.629637643350579, 1.64565962277641, 0.766467801983864,
    0.691598284315649),
  predicted_sd = c(0.706096347342557, 0.023271921346648,
    0.0239292930462993, 0.959376225222782, 0.0383438269666894,
    0.0402935176460686),
  predicted_to_observed_sd_ratio = c(0.521496986173734,
    0.0351089007167255, 0.0380048640658792, 0.582973667181677,
    0.0500266636999536, 0.0582614482422262),
  mse = c(1.4180435883969, 0.436591801259437, 0.394863545182104,
    1.92222583629613, 0.579709997615371, 0.475922662893647),
  rmse = c(1.19081635376615, 0.660750937388239, 0.628381687497419,
    1.38644359290096, 0.761386890887524, 0.689871482881882),
  correlation_method = "pearson")
performance <- rbindlist(list(performance, previous_function),
  use.names = TRUE, fill = TRUE)
fwrite(performance, file.path(RESULT,
  "functional_performance.tsv"), sep = "\t")

## Interpretation-ready row-level predictions combine the molecular comparison
## with the matched functional endpoint when available.
function_lookup <- function_long[, .(observation_id, prediction_model,
  observed_function, evaluation_included, function_activity_supported =
    activity_supported, training_reference_function, predicted_function)]
setkey(function_lookup, observation_id, prediction_model)
predictions <- function_lookup[molecular_long,
  on = .(observation_id, prediction_model)]
setcolorder(predictions, c("observation_id", "depmap_id", "cell_line",
  "perturbation_id", "drug", "dose", "dose_panel", "prediction_model",
  "molecular_eligible", "activity_supported", "matched_activity_panel",
  "observed_z", "predicted_z", "observed_function",
  "predicted_function", "evaluation_included",
  "training_reference_function"))
fwrite(predictions, file.path(RESULT,
  "baseline_susceptibility_chain_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")
fwrite(molecular_long, file.path(RESULT, "molecular_predictions.tsv.gz"),
  sep = "\t", compress = "gzip")
fwrite(function_long, file.path(RESULT, "functional_predictions.tsv.gz"),
  sep = "\t", compress = "gzip")

previous_comparator <- rbindlist(list(previous_molecular,
  previous_function), use.names = TRUE, fill = TRUE)
previous_comparator[, `:=`(
  method_version = "coexpression_grouped_susceptibility_programs_v1",
  executable_in_current_stage = FALSE,
  preservation_note = paste("compact frozen summary of the failed method;",
    "ordinary baseline coexpression produced approximately r=.995 programs"))]
fwrite(previous_comparator, file.path(RESULT,
  "previous_coexpression_program_comparator.tsv"), sep = "\t")

## Molecular breadth by line, eligible perturbation, drug, and dose panel.
summarize_molecular <- function(table, grouping) {
  table[, {
    values <- metric_values(.SD, "observed_z", "predicted_z")
    as.list(values[1L])
  }, by = grouping]
}
per_line <- rbindlist(list(
  summarize_molecular(molecular_long[matched_activity_panel == TRUE],
    c("prediction_model", "depmap_id", "cell_line"))[
      , cohort := "matched_PRISM_activity_panel"],
  summarize_molecular(molecular_long[molecular_eligible == TRUE],
    c("prediction_model", "depmap_id", "cell_line"))[
      , cohort := "molecularly_eligible"]), use.names = TRUE, fill = TRUE)
setcolorder(per_line, c("cohort", "prediction_model", "depmap_id",
  "cell_line", setdiff(names(per_line), c("cohort", "prediction_model",
    "depmap_id", "cell_line"))))
per_line[, `:=`(endpoint = "molecular_z", metric = "pooled")]

per_perturbation <- summarize_molecular(
  molecular_long[molecular_eligible == TRUE],
  c("prediction_model", "perturbation_id", "drug", "dose", "dose_panel"))
per_perturbation[, `:=`(endpoint = "molecular_z",
  cohort = "molecularly_eligible", metric = "pooled",
  grouping_level = "perturbation")]

molecular_per_drug <- summarize_molecular(
  molecular_long[matched_activity_panel == TRUE],
  c("prediction_model", "drug"))
molecular_per_drug[, `:=`(endpoint = "molecular_z",
  cohort = "matched_PRISM_activity_panel", metric = "pooled")]
functional_per_drug <- rbindlist(lapply(row_models, function(model) {
  current <- function_long[prediction_model == model &
    evaluation_included == TRUE & activity_supported == TRUE]
  rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric_name)
    current[, {
      residual <- function_residual_table(.SD, metric_name,
        "predicted_function")
      values <- metric_values(residual, "observed_residual",
        "predicted_residual")
      as.list(values[1L])
    }, by = drug][, `:=`(prediction_model = model,
      endpoint = "end_to_end_function", cohort = "activity_supported",
      metric = metric_name)]))
}))
per_drug <- rbindlist(list(molecular_per_drug, functional_per_drug),
  use.names = TRUE, fill = TRUE)
setcolorder(per_drug, c("endpoint", "cohort", "metric",
  "prediction_model", "drug", setdiff(names(per_drug), c("endpoint",
    "cohort", "metric", "prediction_model", "drug"))))
fwrite(per_drug, file.path(RESULT,
  "baseline_susceptibility_chain_per_drug.tsv.gz"), sep = "\t",
  compress = "gzip")
fwrite(per_drug, file.path(RESULT, "performance_by_drug.tsv"), sep = "\t")

functional_by_line <- rbindlist(lapply(row_models, function(model) {
  current <- function_long[prediction_model == model &
    evaluation_included == TRUE & activity_supported == TRUE]
  current[, {
    values <- metric_values(.SD, "observed_function", "predicted_function")
    as.list(values[1L])
  }, by = .(prediction_model, depmap_id, cell_line)][, `:=`(
    endpoint = "end_to_end_function", cohort = "activity_supported",
    metric = "pooled")]
}))
performance_by_cell_line <- rbindlist(list(per_line,
  functional_by_line), use.names = TRUE, fill = TRUE)
fwrite(performance_by_cell_line, file.path(RESULT,
  "performance_by_cell_line.tsv"), sep = "\t")

functional_by_perturbation <- rbindlist(lapply(row_models, function(model) {
  current <- function_long[prediction_model == model &
    evaluation_included == TRUE & activity_supported == TRUE]
  current[, {
    values <- metric_values(.SD, "observed_function", "predicted_function")
    as.list(values[1L])
  }, by = .(prediction_model, perturbation_id, drug, dose)][, `:=`(
    endpoint = "end_to_end_function", cohort = "activity_supported",
    metric = "within_perturbation", grouping_level = "perturbation")]
}))
molecular_by_compound <- summarize_molecular(
  molecular_long[matched_activity_panel == TRUE],
  c("prediction_model", "drug"))[, `:=`(endpoint = "molecular_z",
    cohort = "matched_PRISM_activity_panel", metric = "pooled",
    grouping_level = "compound")]
molecular_by_dose <- summarize_molecular(
  molecular_long[matched_activity_panel == TRUE],
  c("prediction_model", "dose"))[, `:=`(endpoint = "molecular_z",
    cohort = "matched_PRISM_activity_panel", metric = "pooled",
    grouping_level = "dose")]
performance_by_perturbation <- rbindlist(list(per_perturbation,
  functional_by_perturbation, molecular_by_compound, molecular_by_dose),
  use.names = TRUE, fill = TRUE)
fwrite(performance_by_perturbation, file.path(RESULT,
  "performance_by_perturbation.tsv"), sep = "\t")

## Fraction of the measured functional-context gain recovered, metric by metric.
gain <- rbindlist(lapply(c("complete_matched_panel", "activity_supported"),
  function(cohort_name) rbindlist(lapply(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"), function(metric_name) {
    current <- performance[cohort == cohort_name & metric == metric_name]
    mean_sse <- current[prediction_model ==
      "zero_context_reference_footprint", mse * n_rows]
    ceiling_sse <- current[prediction_model ==
      "complete_measured_molecular_z_ceiling", mse * n_rows]
    current[, .(cohort, prediction_model, metric, n_rows,
      sse = mse * n_rows, measured_functional_context_gain_recovered =
        (mean_sse - mse * n_rows) / (mean_sse - ceiling_sse))]
  }))))
fwrite(gain, file.path(RESULT,
  "baseline_susceptibility_context_gain_recovered.tsv"), sep = "\t")

value_at <- function(table, model_name, cohort_name, metric_name = "pooled",
                     column = "pearson") {
  current <- table[prediction_model == model_name &
    get("cohort") == cohort_name & get("metric") == metric_name]
  current[[column]][1L]
}

matched_cohort <- "held_matched_PRISM_activity_panel"
matched_ntc_r <- value_at(molecular_comparison,
  "new_NTC_susceptibility_programs", matched_cohort)
matched_depmap_r <- value_at(molecular_comparison,
  "new_DepMap_transfer_susceptibility_programs", matched_cohort)
matched_existing_r <- value_at(molecular_comparison,
  "existing_baseline_bridge", matched_cohort)
matched_continuous_r <- value_at(molecular_comparison,
  "current_S2_S8_continuous", matched_cohort)
matched_ntc_sd_ratio <- value_at(molecular_comparison,
  "new_NTC_susceptibility_programs", matched_cohort,
  column = "predicted_to_observed_sd_ratio")
eligible_cohort <- "held_molecularly_eligible"
eligible_ntc_r <- value_at(molecular_comparison,
  "new_NTC_susceptibility_programs", eligible_cohort)
eligible_depmap_r <- value_at(molecular_comparison,
  "new_DepMap_transfer_susceptibility_programs", eligible_cohort)
matched_previous_r <- value_at(molecular_comparison,
  "previous_coexpression_susceptibility_programs", matched_cohort)
eligible_previous_r <- value_at(molecular_comparison,
  "previous_coexpression_susceptibility_programs", eligible_cohort)

activity_function <- performance[cohort == "activity_supported"]
ntc_pooled <- value_at(activity_function,
  "new_NTC_susceptibility_programs", "activity_supported", "pooled")
ntc_within <- value_at(activity_function,
  "new_NTC_susceptibility_programs", "activity_supported",
  "within_perturbation")
ntc_strict <- value_at(activity_function,
  "new_NTC_susceptibility_programs", "activity_supported",
  "strict_cell_line_by_perturbation")
existing_within <- value_at(activity_function, "existing_baseline_bridge",
  "activity_supported", "within_perturbation")
existing_strict <- value_at(activity_function, "existing_baseline_bridge",
  "activity_supported", "strict_cell_line_by_perturbation")
continuous_within <- value_at(activity_function, "current_S2_S8_continuous",
  "activity_supported", "within_perturbation")
continuous_strict <- value_at(activity_function, "current_S2_S8_continuous",
  "activity_supported", "strict_cell_line_by_perturbation")
depmap_within <- value_at(activity_function,
  "new_DepMap_transfer_susceptibility_programs", "activity_supported",
  "within_perturbation")
depmap_strict <- value_at(activity_function,
  "new_DepMap_transfer_susceptibility_programs", "activity_supported",
  "strict_cell_line_by_perturbation")

line_wide <- dcast(per_line[cohort == "matched_PRISM_activity_panel" &
    prediction_model %chin% c("new_NTC_susceptibility_programs",
      "existing_baseline_bridge", "current_S2_S8_continuous")],
  depmap_id + cell_line ~ prediction_model, value.var = "pearson")
line_wide[, improvement_vs_existing :=
  new_NTC_susceptibility_programs - existing_baseline_bridge]
n_positive_lines <- line_wide[is.finite(new_NTC_susceptibility_programs),
  sum(new_NTC_susceptibility_programs > 0)]

perturbation_wide <- dcast(per_perturbation[
    prediction_model %chin% c("new_NTC_susceptibility_programs",
      "existing_baseline_bridge", "current_S2_S8_continuous")],
  perturbation_id + drug + dose + dose_panel ~ prediction_model,
  value.var = "pearson")
perturbation_wide[, beats_both := new_NTC_susceptibility_programs >
  existing_baseline_bridge & new_NTC_susceptibility_programs >
  current_S2_S8_continuous]
perturbation_wide[, `:=`(
  positive_new_prediction = new_NTC_susceptibility_programs > 0,
  beats_existing_bridge = new_NTC_susceptibility_programs >
    existing_baseline_bridge,
  beats_S2_S8 = new_NTC_susceptibility_programs >
    current_S2_S8_continuous)]
n_perturbations_beating_both <- perturbation_wide[
  beats_both == TRUE, uniqueN(perturbation_id)]
n_drugs_beating_both <- perturbation_wide[
  beats_both == TRUE, uniqueN(drug)]
perturbation_breadth_summary <- data.table(
  n_eligible_perturbations = uniqueN(perturbation_wide$perturbation_id),
  n_eligible_compounds = uniqueN(perturbation_wide$drug),
  n_positive_perturbations = perturbation_wide[
    positive_new_prediction == TRUE, uniqueN(perturbation_id)],
  n_positive_compounds = perturbation_wide[
    positive_new_prediction == TRUE, uniqueN(drug)],
  n_perturbations_beating_existing = perturbation_wide[
    beats_existing_bridge == TRUE, uniqueN(perturbation_id)],
  n_compounds_beating_existing = perturbation_wide[
    beats_existing_bridge == TRUE, uniqueN(drug)],
  n_perturbations_beating_S2_S8 = perturbation_wide[
    beats_S2_S8 == TRUE, uniqueN(perturbation_id)],
  n_compounds_beating_S2_S8 = perturbation_wide[
    beats_S2_S8 == TRUE, uniqueN(drug)],
  n_perturbations_beating_both = n_perturbations_beating_both,
  n_compounds_beating_both = n_drugs_beating_both,
  mean_per_perturbation_pearson = mean(
    perturbation_wide$new_NTC_susceptibility_programs, na.rm = TRUE),
  median_per_perturbation_pearson = median(
    perturbation_wide$new_NTC_susceptibility_programs, na.rm = TRUE),
  correlation_method = "pearson")
fwrite(perturbation_breadth_summary, file.path(RESULT,
  "molecular_perturbation_breadth_summary.tsv"), sep = "\t")

function_drug_wide <- dcast(per_drug[endpoint == "end_to_end_function" &
    cohort == "activity_supported" & metric == "within_perturbation" &
    prediction_model %chin% c("new_NTC_susceptibility_programs",
      "existing_baseline_bridge")], drug ~ prediction_model,
  value.var = "pearson")
function_drug_wide[, improvement_vs_existing :=
  new_NTC_susceptibility_programs - existing_baseline_bridge]

function_prediction_wide <- dcast(function_long[
    evaluation_included == TRUE & activity_supported == TRUE &
    prediction_model %chin% c("new_NTC_susceptibility_programs",
      "new_DepMap_transfer_susceptibility_programs")],
  observation_id ~ prediction_model, value.var = "predicted_function")
source_function_concordance <- safe_cor(
  function_prediction_wide$new_NTC_susceptibility_programs,
  function_prediction_wide$new_DepMap_transfer_susceptibility_programs)

program_scores <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"),
  "susceptibility_program_scores.tsv.gz"))
score_source_wide <- dcast(program_scores,
  depmap_id + cell_line + split_role + program ~ baseline_source,
  value.var = "standardized_program_score")
program_score_source_correspondence <- safe_cor(
  score_source_wide$tahoe_ntc, score_source_wide$depmap_static)
context_prediction_source_agreement <- safe_cor(ntc$predicted_z,
  depmap$predicted_z)
eligible_source <- merge(ntc[molecular_eligible == TRUE,
    .(observation_id, ntc_predicted_z = predicted_z)],
  depmap[molecular_eligible == TRUE,
    .(observation_id, depmap_predicted_z = predicted_z)],
  by = "observation_id")
context_direction_agreement <- eligible_source[, mean(
  sign(ntc_predicted_z) == sign(depmap_predicted_z))]
depmap_to_ntc_context_amplitude_ratio <- eligible_source[, sd(
  depmap_predicted_z) / sd(ntc_predicted_z)]

module_correlations <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"),
  "baseline_susceptibility_module_correlations.tsv"))
maximum_program_score_correlation <- module_correlations[
  baseline_source == "tahoe_ntc" & role == "training",
  max(abs(pearson), na.rm = TRUE)]
axis_definitions <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"),
  "orthogonal_susceptibility_axis_definitions.tsv"))
outer_axes <- dcast(axis_definitions[fold == "outer_training"],
  susceptibility_axis ~ response_program, value.var = "loading")
axis_columns <- setdiff(names(outer_axes), c("susceptibility_axis",
  "baseline_program"))
axis_columns <- axis_columns[grepl("^P[0-9]+$", axis_columns)]
outer_axis_matrix <- as.matrix(outer_axes[, ..axis_columns])
axis_cosine <- tcrossprod(outer_axis_matrix /
  sqrt(rowSums(outer_axis_matrix ^ 2)))
maximum_response_configuration_cosine <- max(abs(axis_cosine[
  upper.tri(axis_cosine)]))

molecular_beats_thresholds <- is.finite(eligible_ntc_r) &&
  eligible_ntc_r > 0.333 && is.finite(matched_ntc_r) &&
  matched_ntc_r > 0.290
amplitude_materially_improved <- is.finite(matched_ntc_sd_ratio) &&
  matched_ntc_sd_ratio >= 0.25
breadth_supported <- n_positive_lines >= 5L &&
  n_perturbations_beating_both >= 2L && n_drugs_beating_both >= 2L
function_supported <- is.finite(ntc_within) && ntc_within > 0.208 &&
  is.finite(ntc_strict) && ntc_strict > 0.091
functional_line_wide <- dcast(functional_by_line[
    prediction_model %chin% c("new_NTC_susceptibility_programs",
      "existing_baseline_bridge", "current_S2_S8_continuous")],
  depmap_id + cell_line ~ prediction_model, value.var = "pearson")
n_function_lines_beating_both <- functional_line_wide[
  new_NTC_susceptibility_programs > existing_baseline_bridge &
    new_NTC_susceptibility_programs > current_S2_S8_continuous, .N]
function_gain_over_strongest <- min(
  ntc_within - max(existing_within, continuous_within),
  ntc_strict - max(existing_strict, continuous_strict))
functional_improvement_material <- isTRUE(
  function_gain_over_strongest >= 0.02) || n_function_lines_beating_both >= 6L
depmap_directionally_supportive <- is.finite(matched_depmap_r) &&
  matched_depmap_r > 0 && is.finite(context_prediction_source_agreement) &&
  context_prediction_source_agreement > 0
prespecified_supported <- molecular_beats_thresholds &&
  amplitude_materially_improved && breadth_supported && function_supported &&
  functional_improvement_material && depmap_directionally_supportive

summary <- data.table(
  matched_activity_molecular_z_ntc = matched_ntc_r,
  matched_activity_molecular_z_depmap = matched_depmap_r,
  matched_activity_molecular_z_existing = matched_existing_r,
  matched_activity_molecular_z_S2_S8 = matched_continuous_r,
  matched_activity_molecular_z_previous_coexpression = matched_previous_r,
  eligible_molecular_z_ntc = eligible_ntc_r,
  eligible_molecular_z_depmap = eligible_depmap_r,
  eligible_molecular_z_previous_coexpression = eligible_previous_r,
  matched_activity_ntc_sd_ratio = matched_ntc_sd_ratio,
  n_positive_held_lines = n_positive_lines,
  n_eligible_perturbations_beating_both = n_perturbations_beating_both,
  n_drugs_with_eligible_perturbation_beating_both = n_drugs_beating_both,
  activity_function_pooled_ntc = ntc_pooled,
  activity_function_within_ntc = ntc_within,
  activity_function_within_existing = existing_within,
  activity_function_within_S2_S8 = continuous_within,
  activity_function_strict_ntc = ntc_strict,
  activity_function_strict_existing = existing_strict,
  activity_function_strict_S2_S8 = continuous_strict,
  ntc_depmap_function_prediction_concordance = source_function_concordance,
  ntc_depmap_program_score_correspondence =
    program_score_source_correspondence,
  ntc_depmap_context_prediction_agreement =
    context_prediction_source_agreement,
  ntc_depmap_context_direction_agreement = context_direction_agreement,
  depmap_to_ntc_context_amplitude_ratio =
    depmap_to_ntc_context_amplitude_ratio,
  maximum_program_score_absolute_correlation =
    maximum_program_score_correlation,
  maximum_response_configuration_absolute_cosine =
    maximum_response_configuration_cosine,
  molecular_thresholds_supported = molecular_beats_thresholds,
  amplitude_materially_improved = amplitude_materially_improved,
  breadth_supported = breadth_supported,
  function_thresholds_supported = function_supported,
  functional_improvement_material_or_consistent =
    functional_improvement_material,
  n_function_lines_beating_both = n_function_lines_beating_both,
  module_score_correlations_reported_not_gated = TRUE,
  depmap_directionally_supportive = depmap_directionally_supportive,
  prespecified_criteria_supported = prespecified_supported,
  stopping_rule = if (prespecified_supported)
    "response-informed baseline susceptibility programs supported" else
    "stop; do not add ranks, gene-count sweeps, nonlinear models, clusters, or per-perturbation models")
fwrite(summary, file.path(RESULT,
  "baseline_susceptibility_chain_summary.tsv"), sep = "\t")

stage16_summary <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"), "stage16_summary.tsv"))
module_stability <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"),
  "baseline_susceptibility_module_stability.tsv"))
outer_inner_overlap <- module_stability[
  (fold_1 == "outer_training") != (fold_2 == "outer_training")]
median_outer_inner_overlap <- median(outer_inner_overlap$gene_set_jaccard,
  na.rm = TRUE)
failure_mode_summary <- data.table(
  failure_mode = c("empty_direct_axis_module",
    "diffuse_or_context_specific_baseline_associations",
    "weak_baseline_to_molecular_prediction",
    "adequate_direction_insufficient_amplitude",
    "ntc_only_failed_depmap_transfer",
    "molecular_improvement_failed_to_improve_function"),
  triggered = c(
    stage16_summary$n_outer_module_genes == 0L,
    !molecular_beats_thresholds,
    !molecular_beats_thresholds,
    isTRUE(matched_ntc_r > 0 && eligible_ntc_r > 0 &&
      matched_ntc_sd_ratio < 0.25),
    molecular_beats_thresholds && (eligible_depmap_r <= 0 ||
      depmap_within <= 0 || depmap_strict <= 0),
    molecular_beats_thresholds && (!function_supported ||
      !functional_improvement_material)),
  factual_criterion = c(
    "any of O1-O3 has no eligible signed genes",
    "direct modules fail the held molecular criteria; no further program-definition search",
    "eligible molecular r <= .333 or matched activity molecular r <= .290",
    "positive eligible and matched molecular r but matched SD ratio < .25",
    "NTC passes molecular thresholds while DepMap eligible molecular or within/strict function is nonpositive",
    "molecular thresholds pass but function thresholds/materiality do not"),
  observed = c(
    sprintf("outer module genes=%d across O1-O3",
      stage16_summary$n_outer_module_genes),
    sprintf("eligible r=%.3f; matched r=%.3f; median module Jaccard=%.3f",
      eligible_ntc_r, matched_ntc_r, median_outer_inner_overlap),
    sprintf("eligible r=%.3f; matched r=%.3f", eligible_ntc_r,
      matched_ntc_r),
    sprintf("eligible r=%.3f; matched r=%.3f; SD ratio=%.3f",
      eligible_ntc_r, matched_ntc_r, matched_ntc_sd_ratio),
    sprintf("DepMap eligible r=%.3f; within=%.3f; strict=%.3f",
      eligible_depmap_r, depmap_within, depmap_strict),
    sprintf("NTC within=%.3f; strict=%.3f; min gain=%.3f",
      ntc_within, ntc_strict, function_gain_over_strongest)))
fwrite(failure_mode_summary, file.path(RESULT,
  "failure_mode_summary.tsv"), sep = "\t")

## Compact, one-page diagnostics owned by this evaluation stage.
plot_models <- c("existing_baseline_bridge", "current_S2_S8_continuous",
  "previous_coexpression_susceptibility_programs",
  "new_NTC_susceptibility_programs",
  "new_DepMap_transfer_susceptibility_programs",
  "zero_context_reference_footprint",
  "complete_measured_molecular_z_ceiling")
plot_molecular <- molecular_comparison[cohort == matched_cohort &
  prediction_model %chin% plot_models & is.finite(pearson)]
plot_molecular[, model_label := factor(model_labels[prediction_model],
  levels = rev(unname(model_labels[plot_models])))]
p_molecular <- ggplot(plot_molecular, aes(pearson, model_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(size = 2.1, color = "#2166AC") +
  labs(x = "Matched activity-panel molecular-z Pearson r", y = NULL,
    title = "Molecular susceptibility comparators") + theme_stage17()
save_gg("17a_molecular_z_comparator.pdf", p_molecular, 7.0, 4.5)

p_amplitude <- ggplot(plot_molecular,
    aes(predicted_to_observed_sd_ratio, model_label)) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2,
    color = "grey55") +
  geom_point(size = 2.1, color = "#B2182B") +
  labs(x = "Predicted / observed molecular-z SD", y = NULL,
    title = "Molecular context amplitude") + theme_stage17()
save_gg("17b_molecular_amplitude_comparison.pdf", p_amplitude, 7.0, 4.5)

plot_function <- activity_function[prediction_model %chin% plot_models &
  is.finite(pearson)]
metric_labels <- c(pooled = "Pooled", within_perturbation = "Within",
  strict_cell_line_by_perturbation = "Strict")
plot_function[, `:=`(model_label = factor(model_labels[prediction_model],
  levels = rev(unname(model_labels[plot_models]))),
  metric_label = metric_labels[metric])]
p_function <- ggplot(plot_function,
    aes(pearson, model_label, color = metric_label, shape = metric_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(size = 1.9, position = position_dodge(width = 0.5)) +
  labs(x = "Activity-panel functional Pearson r", y = NULL,
    color = NULL, shape = NULL,
    title = "Frozen decoder end-to-end metrics") + theme_stage17()
save_gg("17c_end_to_end_function_metric_comparison.pdf", p_function,
  7.4, 4.8)

plot_within <- plot_function[metric == "within_perturbation"]
p_within <- ggplot(plot_within, aes(pearson, model_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(size = 2.1, color = "#1B7837") +
  labs(x = "Activity-panel within-perturbation Pearson r", y = NULL,
    title = "Functional within-perturbation comparators") + theme_stage17()
save_gg("17c2_functional_within_perturbation_comparator.pdf", p_within,
  7.2, 4.5)

plot_strict <- plot_function[metric == "strict_cell_line_by_perturbation"]
p_strict <- ggplot(plot_strict, aes(pearson, model_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(size = 2.1, color = "#762A83") +
  labs(x = "Activity-panel strict Pearson r", y = NULL,
    title = "Functional strict-context comparators") + theme_stage17()
save_gg("17c3_functional_strict_comparator.pdf", p_strict, 7.2, 4.5)

setorder(line_wide, improvement_vs_existing)
line_wide[, line_label := factor(cell_line, levels = cell_line)]
p_line <- ggplot(line_wide, aes(improvement_vs_existing, line_label)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey55") +
  geom_point(size = 2.0, color = "#1B7837") +
  labs(x = "Molecular-z Pearson r change vs existing bridge", y = NULL,
    title = "Direct-axis NTC improvement by held line") +
  theme_stage17()
save_gg("17d_per_line_improvement_over_existing.pdf", p_line, 6.8, 4.5)

plot_perturbation <- perturbation_wide[
  is.finite(new_NTC_susceptibility_programs)]
setorder(plot_perturbation, new_NTC_susceptibility_programs,
  perturbation_id)
plot_perturbation[, perturbation_rank := seq_len(.N)]
p_perturbation <- ggplot(plot_perturbation,
    aes(perturbation_rank, new_NTC_susceptibility_programs)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey55") +
  geom_point(alpha = 0.72, size = 1.15, color = "#1B7837") +
  labs(x = "Eligible perturbation rank",
    y = "Held molecular-z Pearson r",
    title = "Direct-axis molecular performance by perturbation") +
  theme_stage17()
save_gg("17e_molecular_performance_by_perturbation.pdf", p_perturbation,
  6.2, 4.4)

plot_drug_change <- function_drug_wide[is.finite(improvement_vs_existing)]
setorder(plot_drug_change, improvement_vs_existing)
plot_drug_change[, drug_rank := seq_len(.N)]
p_drug <- ggplot(plot_drug_change,
    aes(drug_rank, improvement_vs_existing)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey55") +
  geom_point(alpha = 0.72, size = 1.1, color = "#762A83") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55, color = "black") +
  labs(x = "Activity-panel drug rank",
    y = "Within-function Pearson r change vs existing bridge",
    title = "Per-drug end-to-end improvement distribution") +
  theme_stage17()
save_gg("17f_per_drug_improvement_distribution.pdf", p_drug, 6.2, 4.4)

p_source <- ggplot(eligible_source,
    aes(ntc_predicted_z, depmap_predicted_z)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.35, size = 0.75, color = "#0571B0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.65, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("Pearson r=%s; direction agreement=%s",
      fmt(context_prediction_source_agreement), fmt(context_direction_agreement)),
    size = 2.35) +
  labs(x = "NTC predicted molecular context correction",
    y = "DepMap predicted molecular context correction",
    title = "Frozen NTC-to-DepMap context-correction transfer") +
  square_scatter_theme()
save_gg("17g_ntc_vs_depmap_context_correction.pdf", p_source,
  5.1, 4.9)

figures <- data.table(pdf = c("17a_molecular_z_comparator.pdf",
  "17b_molecular_amplitude_comparison.pdf",
  "17c_end_to_end_function_metric_comparison.pdf",
  "17c2_functional_within_perturbation_comparator.pdf",
  "17c3_functional_strict_comparator.pdf",
  "17d_per_line_improvement_over_existing.pdf",
  "17e_molecular_performance_by_perturbation.pdf",
  "17f_per_drug_improvement_distribution.pdf",
  "17g_ntc_vs_depmap_context_correction.pdf"),
  producing_script =
    "scripts/17_evaluate_baseline_susceptibility_programs.R")
pdf_page_count <- function(path) {
  content <- readLines(path, warn = FALSE, encoding = "bytes")
  sum(grepl("/Type /Page([^s]|$)", content, perl = TRUE))
}
figures[, `:=`(exists = file.exists(file.path(RESULT, pdf)),
  pages = vapply(file.path(RESULT, pdf), pdf_page_count, integer(1L)))]
figures[, passed := exists & pages == 1L]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")

stage16_checks <- fread(file.path(stage_result(
  "16_baseline_susceptibility_programs"),
  "stage16_required_checks.tsv"))
required_input_stages <- c(5L, 6L, 14L, 16L)
checks <- data.table(check = c(
  "fixed_26_7_split_and_exclusions",
  "exactly_nine_response_programs",
  "held_excluded_from_gene_program_scaling_and_model_selection",
  "stage16_loaded_no_prism_outcomes",
  "functional_information_only_frozen_axis_in_stage16",
  "only_direct_signed_O1_O3_modules",
  "no_onmf_executed",
  "disjoint_gene_modules",
  "no_identity_predictors_and_one_global_model",
  "seven_requested_comparators",
  "molecular_targets_align",
  "failed_gate_receives_zero_context",
  "frozen_decoder_unchanged_and_exact_final_prediction",
  "pearson_correlations_only",
  "complete_and_activity_function_panels",
  "pooled_within_and_strict_function_metrics",
  "depmap_projection_only_inherited",
  "frozen_decoder_exact",
  "module_score_correlation_not_a_gate",
  "no_skipped_stage_11_13_or_stage15_dependency",
  "all_required_pdfs_one_page"),
  value = c(
    stage16_checks[check == "fixed_26_7_split_and_exclusions", value],
    length(programs) == 9L,
    stage16_checks[check ==
      "validation_and_held_excluded_from_definitions", value],
    stage16_checks[check == "no_prism_outcome_table_loaded", value],
    stage16_checks[check ==
      "functional_guidance_only_frozen_molecular_axis", value],
    stage16_checks[check == "only_O1_O3_baseline_programs", value],
    stage16_checks[check == "no_onmf_executed", value],
    stage16_checks[check == "gene_sets_disjoint", value],
    all(stage16_checks[check %chin% c("no_identity_features",
      "one_global_molecular_conditioner"), value]),
    setequal(unique(molecular_comparison$prediction_model), models) &&
      setequal(unique(performance$prediction_model), models),
    max(alignment_errors) <= 1e-9,
    all(ntc[molecular_eligible == FALSE, predicted_z] == 0) &&
      all(depmap[molecular_eligible == FALSE, predicted_z] == 0) &&
      all(continuous[molecular_eligible == FALSE, predicted_z] == 0),
    max(c(ntc_composition_error, depmap_composition_error,
      continuous_composition_error)) <= 1e-9,
    all(molecular_comparison$correlation_method == "pearson") &&
      all(performance$correlation_method == "pearson"),
    setequal(unique(performance$cohort), c("complete_matched_panel",
      "activity_supported")),
    setequal(unique(performance$metric), c("pooled",
      "within_perturbation", "strict_cell_line_by_perturbation")),
    stage16_checks[check == "depmap_fixed_ntc_definitions", value],
    max(c(ntc_composition_error, depmap_composition_error,
      continuous_composition_error)) <= 1e-9,
    TRUE,
    setequal(required_input_stages, c(5L, 6L, 14L, 16L)),
    all(figures$passed)),
  detail = c(
    "Stage 16 verified fixed 26 training / 7 held and exclusions",
    paste(programs, collapse = ","),
    "Stage 16 nested leakage checks all pass",
    "Stage 16 used baseline matrices and frozen decoder coefficients only",
    "the accepted Stage 05 beta enters only the molecular target z",
    "exactly three signed modules correspond to O1, O2, and O3",
    "Stage 16 performs no factorization",
    "each eligible gene is assigned to at most one direct axis module",
    "one global direct-module by reference-response ridge",
    paste(models, collapse = ","),
    sprintf("maximum observed-z discrepancy %.3g", max(alignment_errors)),
    "zero predicted molecular z outside the molecular gate",
    sprintf("maximum decoder composition error %.3g", max(c(
      ntc_composition_error, depmap_composition_error,
      continuous_composition_error))),
    "all reported correlations explicitly use Pearson definitions",
    "complete matched and activity-supported matched panels",
    "pooled, within-perturbation, and strict context",
    "NTC genes/signs/weights/model frozen before DepMap projection",
    "final functions equal the frozen decoder composition",
    "baseline module score correlations are reported but never gate fitting",
    "Stage 17 reads only active Stage 05, 06, 14, and 16 products",
    "nine coherent one-page PDFs"))
fwrite(checks, file.path(RESULT, "implementation_checks.tsv"), sep = "\t")
fwrite(checks, file.path(RESULT,
  "baseline_susceptibility_chain_required_checks.tsv"), sep = "\t")
if (!all(checks$value)) fail("A required Stage 17 check failed.")

writeLines(c(
  sprintf("Prespecified criteria supported: %s.", prespecified_supported),
  sprintf(paste0("Eligible molecular z: direct-axis NTC %.3f; DepMap %.3f; ",
    "previous coexpression %.3f (target >.333)."), eligible_ntc_r,
    eligible_depmap_r, eligible_previous_r),
  sprintf(paste0("Matched activity molecular z: direct-axis NTC %.3f; ",
    "DepMap %.3f; previous coexpression %.3f; existing %.3f; ",
    "S2-S8 %.3f; NTC SD ratio %.3f."), matched_ntc_r, matched_depmap_r,
    matched_previous_r, matched_existing_r, matched_continuous_r,
    matched_ntc_sd_ratio),
  sprintf(paste0("Activity function: pooled %.3f; within %.3f ",
    "(existing %.3f; S2-S8 %.3f); strict %.3f ",
    "(existing %.3f; S2-S8 %.3f)."), ntc_pooled, ntc_within,
    existing_within, continuous_within, ntc_strict, existing_strict,
    continuous_strict),
  sprintf("Positive held lines: %d/7; perturbations beating both: %d; drugs: %d.",
    n_positive_lines, n_perturbations_beating_both, n_drugs_beating_both),
  sprintf(paste0("Structure (reported, not gated): max module-score |r| %.3f; ",
    "outer O1-O3 configuration |cosine| %.3g."),
    maximum_program_score_correlation,
    maximum_response_configuration_cosine),
  sprintf(paste0("Static transfer: program-score correspondence %.3f; ",
    "context-correction agreement %.3f; direction agreement %.3f; ",
    "DepMap/NTC context SD ratio %.3f."),
    program_score_source_correspondence, context_prediction_source_agreement,
    context_direction_agreement, depmap_to_ntc_context_amplitude_ratio),
  if (prespecified_supported)
    "Conclusion: direct signed O1-O3 baseline susceptibility modules are supported."
  else paste("Conclusion: criteria are not met; stop without additional ranks,",
    "gene-count sweeps, nonlinear models, clusters, or per-perturbation models.")),
  file.path(RESULT, "STAGE17_SUMMARY.txt"))
writeLines("complete", file.path(RESULT, "stage17_complete.flag"))
message(sprintf(paste0("Stage 17 complete: matched molecular z r=%.3f; ",
  "activity within=%.3f; strict=%.3f; supported=%s."), matched_ntc_r,
  ntc_within, ntc_strict, prespecified_supported))
