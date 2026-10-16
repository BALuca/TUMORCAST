#!/usr/bin/env Rscript

## Reproduce the accepted DMSO-q99 molecular activity panel using training
## molecular responses only.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("04_activity_panel")
ARTIFACT <- stage_artifact("04_activity_panel")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

plot_stage04 <- function() {
  census <- fread(file.path(RESULT, "condition_activity_census.tsv"))
  magnitude <- fread(file.path(stage_artifact("02_response_programs"),
    "dmso_vs_treated_response_magnitude.tsv.gz"))
  theme_activity <- function(base_size = 8.5) theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
      plot.title = element_text(face = "bold"))
  save_activity <- function(name, plot, width, height) ggsave(
    file.path(RESULT, name), plot, width = width, height = height,
    device = cairo_pdf, bg = "white", limitsize = FALSE)

  p_activity <- ggplot(census, aes(median_training_magnitude,
      n_training_lines_active, color = activity_supported)) +
    geom_hline(yintercept = MIN_ACTIVE_LINES, linetype = 2,
      linewidth = 0.4, color = "grey45") +
    geom_point(alpha = 0.55, size = 0.9) +
    scale_color_manual(values = c(`TRUE` = "#2E7D32",
      `FALSE` = "grey70"), name = "Activity supported") +
    labs(x = "Median training-line molecular magnitude",
      y = "Training lines above DMSO q99",
      title = "Training-only molecular activity support") +
    theme_activity() + theme(aspect.ratio = 1)
  save_activity("04a_activity_support_map.pdf", p_activity, 4.8, 4.6)

  magnitude[, response_label := factor(fcase(
    response_type == "training_DMSO_pseudo",
      "Training DMSO pseudo-response",
    response_type == "treated" & role == "training", "Training treated",
    response_type == "treated" & role == "held_out", "Held-out treated"),
    levels = c("Training DMSO pseudo-response", "Training treated",
      "Held-out treated"))]
  q99 <- unique(magnitude$empirical_dmso_q99)[1L]
  p_magnitude <- ggplot(magnitude, aes(magnitude, color = response_label,
      fill = response_label)) +
    geom_density(alpha = 0.12, linewidth = 0.55, adjust = 1.1) +
    geom_vline(xintercept = q99, linetype = 2, linewidth = 0.45,
      color = "black") +
    labs(x = "Noise-standardized response magnitude", y = "Density",
      color = NULL, fill = NULL,
      title = "DMSO-null and treated-response magnitudes",
      subtitle = "Dashed line: training DMSO empirical q99") +
    theme_activity()
  save_activity("04b_dmso_vs_treated_response_magnitude.pdf", p_magnitude,
    6.2, 4.3)
}

if ("--plots-only" %in% commandArgs(trailingOnly = TRUE)) {
  plot_stage04()
  message("Stage 04 plots refreshed from existing activity tables.")
  quit(save = "no", status = 0L)
}

observed <- fread(file.path(stage_artifact("02_response_programs"),
  "observed_response_metadata.tsv.gz"))
dmso <- fread(file.path(stage_artifact("02_response_programs"),
  "dmso_response_magnitude.tsv"))
function_rows <- fread(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"))
folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- folds[outer_fold != OUTER_FOLD, depmap_id]
held_lines <- folds[outer_fold == OUTER_FOLD, depmap_id]
threshold <- unname(quantile(dmso[depmap_id %chin% training_lines,
  magnitude], 0.99, type = 7L))
if (!is.finite(threshold)) fail("The training DMSO q99 is nonfinite.")

training <- observed[depmap_id %chin% training_lines]
census <- training[, .(
  drug = drug[1L], dose = dose[1L], dose_panel = dose_panel[1L],
  n_training_lines_observed = uniqueN(depmap_id),
  n_training_lines_active = sum(
    selected_gene_noise_standardized_magnitude > threshold),
  median_training_magnitude = median(
    selected_gene_noise_standardized_magnitude)), by = perturbation_id]
census[, `:=`(
  minimum_coverage_lines = ceiling(length(training_lines) / 2),
  individual_activity_dmso_q99 = threshold,
  coverage_pass = n_training_lines_observed >= ceiling(length(training_lines) / 2),
  activity_pass = n_training_lines_active >= MIN_ACTIVE_LINES)]
census[, activity_supported := coverage_pass & activity_pass]
census[, complete_prism_overlap := perturbation_id %chin%
  unique(function_rows[outer_fold == OUTER_FOLD & chain_observed == 1L,
    perturbation_id])]
census[, activity_supported_prism_overlap :=
  activity_supported & complete_prism_overlap]
census[, `:=`(
  n_held_molecular_rows = vapply(perturbation_id, function(id)
    observed[depmap_id %chin% held_lines & perturbation_id == id, .N],
    integer(1L)),
  n_held_function_rows = vapply(perturbation_id, function(id)
    function_rows[depmap_id %chin% held_lines & perturbation_id == id &
      chain_observed == 1L, .N], integer(1L)),
  selection_uses_training_molecular_only = TRUE,
  prism_used_for_activity_selection = FALSE,
  held_out_values_used_for_selection = FALSE)]
setorder(census, drug, dose)

membership <- merge(observed,
  census[, .(perturbation_id, activity_supported, complete_prism_overlap,
    activity_supported_prism_overlap)], by = "perturbation_id", all.x = TRUE,
  sort = FALSE)
setorder(membership, selected_response_row)
fwrite(census, file.path(RESULT, "condition_activity_census.tsv"), sep = "\t")
fwrite(membership, file.path(ARTIFACT, "molecular_activity_membership.tsv.gz"),
  sep = "\t", compress = "gzip")

summary <- data.table(
  cohort = c("all_molecular", "activity_supported",
    "complete_matched_prism", "activity_supported_prism"),
  n_conditions = c(uniqueN(census$perturbation_id),
    census[activity_supported == TRUE, uniqueN(perturbation_id)],
    census[complete_prism_overlap == TRUE, uniqueN(perturbation_id)],
    census[activity_supported_prism_overlap == TRUE,
      uniqueN(perturbation_id)]),
  n_drugs = c(uniqueN(census$drug),
    census[activity_supported == TRUE, uniqueN(drug)],
    census[complete_prism_overlap == TRUE, uniqueN(drug)],
    census[activity_supported_prism_overlap == TRUE, uniqueN(drug)]),
  n_held_molecular_rows = c(sum(census$n_held_molecular_rows),
    census[activity_supported == TRUE, sum(n_held_molecular_rows)],
    census[complete_prism_overlap == TRUE, sum(n_held_molecular_rows)],
    census[activity_supported_prism_overlap == TRUE,
      sum(n_held_molecular_rows)]),
  n_held_function_rows = c(sum(census$n_held_function_rows),
    census[activity_supported == TRUE, sum(n_held_function_rows)],
    census[complete_prism_overlap == TRUE, sum(n_held_function_rows)],
    census[activity_supported_prism_overlap == TRUE,
      sum(n_held_function_rows)]))
fwrite(summary, file.path(RESULT, "activity_panel_summary.tsv"), sep = "\t")
checks <- data.table(check = c("complete_prism_conditions",
  "complete_prism_drugs", "complete_prism_held_rows",
  "activity_supported_conditions",
  "activity_supported_drugs", "activity_supported_held_molecular_rows",
  "activity_prism_conditions", "activity_prism_drugs",
  "activity_prism_held_rows", "selection_outcome_independent"),
  observed = c(summary[cohort == "complete_matched_prism", n_conditions],
    summary[cohort == "complete_matched_prism", n_drugs],
    summary[cohort == "complete_matched_prism", n_held_function_rows],
    summary[cohort == "activity_supported", n_conditions],
    summary[cohort == "activity_supported", n_drugs],
    summary[cohort == "activity_supported", n_held_molecular_rows],
    summary[cohort == "activity_supported_prism", n_conditions],
    summary[cohort == "activity_supported_prism", n_drugs],
    summary[cohort == "activity_supported_prism", n_held_function_rows], 1L),
  expected = c(262L, 88L, 1752L, 356L, 216L, 2473L, 95L, 54L,
    646L, 1L))
checks[, passed := observed == expected]
fwrite(checks, file.path(RESULT, "activity_panel_checks.tsv"), sep = "\t")
if (!all(checks$passed)) fail("The accepted activity panel was not reproduced.")
write_artifact_manifest(ARTIFACT, "scripts/04_define_activity_panel.R")
writeLines("complete", file.path(RESULT, "stage04_complete.flag"))
plot_stage04()
message("Stage 04 complete: 356 active conditions; 95 overlap PRISM.")
