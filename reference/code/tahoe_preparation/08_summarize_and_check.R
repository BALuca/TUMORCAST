#!/usr/bin/env Rscript

## Enforce the Stage-1 reproduction gate and index stage-owned diagnostics.
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("08_summary")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)

targets <- fread(file.path(root, "reproduction_targets.tsv"))
molecular_reference <- fread(file.path(stage_result(
  "05_program_function_bridge"), "molecular_reference_performance.tsv"))
bridge <- fread(file.path(stage_result("05_program_function_bridge"),
  "program_function_performance.tsv"))
measured <- fread(file.path(stage_result("05_program_function_bridge"),
  "measured_function_relevant_performance.tsv"))
baseline <- fread(file.path(stage_result("06_baseline_bridge"),
  "baseline_molecular_performance.tsv"))
end_to_end <- fread(file.path(stage_result("07_end_to_end"),
  "end_to_end_performance.tsv"))
headroom <- fread(file.path(stage_result("07_end_to_end"),
  "functional_headroom_recovered.tsv"))

lookup <- function(group, metric, source) {
  if (group == "molecular_reference") return(switch(metric,
    held_out_pearson = molecular_reference$held_out_pearson,
    within_perturbation_pearson = molecular_reference$within_perturbation_pearson))
  if (group %chin% c("complete_matched_prism",
      "activity_supported_prism")) {
    target_cohort <- if (group == "complete_matched_prism")
      "complete_matched_panel" else "activity_supported"
    model <- if (grepl("perturbation_average", metric))
      "perturbation_average" else "complete_measured"
    metric_name <- if (grepl("pooled", metric)) "pooled" else
      if (grepl("within", metric)) "within_perturbation" else
        "strict_cell_line_by_perturbation"
    return(bridge[cohort == target_cohort & prediction_model == model &
      get("metric") == metric_name, pearson][1L])
  }
  if (group == "measured_function_relevant") return(measured[[metric]][1L])
  if (group == "baseline_all_active") {
    column <- switch(metric, pearson = "pearson", partial_r2 = "partial_r2",
      sd_ratio = "predicted_to_observed_sd_ratio")
    return(baseline[baseline_source == source &
      cohort == "all_activity_supported", get(column)][1L])
  }
  if (group == "baseline_function_panel")
    return(baseline[baseline_source == source &
      cohort == "activity_supported_prism", pearson][1L])
  if (group == "end_to_end_activity") {
    if (metric == "headroom") return(headroom[baseline_source == source &
      cohort == "activity_supported", headroom_recovered][1L])
    metric_name <- switch(metric, pooled = "pooled",
      within = "within_perturbation",
      strict = "strict_cell_line_by_perturbation")
    return(end_to_end[baseline_source == source &
      prediction_model == "baseline_conditioned" &
      cohort == "activity_supported" & get("metric") == metric_name,
      pearson][1L])
  }
  NA_real_
}
targets[, observed := mapply(lookup, target_group, metric,
  fifelse(is.na(baseline_source), "", baseline_source))]
targets[, `:=`(absolute_difference = abs(observed - target),
  passed = is.finite(observed) & abs(observed - target) <= tolerance)]
fwrite(targets, file.path(RESULT, "reproduction_checks.tsv"), sep = "\t")

cohort_checks <- rbindlist(list(
  fread(file.path(stage_result("04_activity_panel"),
    "activity_panel_checks.tsv"))[, .(check, observed, expected, passed)],
  data.table(check = "all_numeric_reproduction_targets",
    observed = sum(targets$passed), expected = nrow(targets),
    passed = all(targets$passed))))
fwrite(cohort_checks, file.path(RESULT, "required_checks.tsv"), sep = "\t")
if (!all(targets$passed) || !all(cohort_checks$passed)) {
  first <- targets[passed == FALSE][1L]
  fail("Stage-1 reproduction gate failed first at %s/%s/%s: observed %.10f; target %.10f; tolerance %.4g.",
    first$target_group, first$metric, first$baseline_source,
    first$observed, first$target, first$tolerance)
}

refactor_manifest <- data.table(
  original_script = c(
    "00_audit_overlap.R", "02_prepare_tahoe_prism.R",
    "03_run_response_onmf_single_fold.py", "03_run_response_onmf.py",
    "03b_run_ntc_baseline_onmf_single_fold.py",
    "07_context_signal_failure_audit.R", "04_fit_selected_chain.R",
    "09_fit_function_relevant_context_chain.R",
    "09_fit_function_relevant_context_chain.R",
    "09_fit_function_relevant_context_chain.R",
    "05_plot_and_summarize.R",
    "03_run_response_onmf_single_fold.py",
    "09a_build_rich_response_programs.py",
    "04_fit_selected_chain.R and discarded downstream branches"),
  original_section_or_function = c(
    "exact-name PRISM overlap and dose matching",
    "corrected cohort, baselines and function assembly",
    "paired response normalization, gene selection, fold-1 oNMF and scoring",
    "posneg correlation, repeated oNMF, consensus and orientation helpers",
    "DMSO line-mean baseline discovery and DepMap scoring",
    "make_fold_census activity q99/coverage rule",
    "shared ridge decoder and held function prediction",
    "measured_z and reference-footprint construction",
    "single-inner-fold global bilinear baseline bridge",
    "end-to-end metrics, headroom and block bootstrap",
    "all retained diagnostic plotting blocks",
    "PC1 confounding diagnostic", "rank-18 response-program experiment",
    "PCA/full-gene/context-residual/two-head/shrinkage/nearest-neighbor/history branches"),
  new_v2_script = c("01_prepare_inputs.R", "01_prepare_inputs.R",
    "02_run_response_onmf.py", "02_onmf_helpers.py",
    "03_run_baseline_onmf.py", "04_define_activity_panel.R",
    "05_fit_program_function_bridge.R", "05_fit_program_function_bridge.R",
    "06_fit_function_relevant_baseline_bridge.R",
    "07_evaluate_end_to_end.R",
    "01_prepare_inputs.R through 07_evaluate_end_to_end.R",
    NA_character_, NA_character_, NA_character_),
  disposition = c(rep("minimally adapted", 11L), rep("dropped", 3L)),
  reason = c(
    "rebuild joins from canonical local files and redirect v2 outputs",
    "remove old-output dependencies, reticulate and five-fold scaffolding",
    "redirect outputs and remove PCA-only raw matrix/diagnostic",
    "genuinely reused by both native oNMF stages",
    "redirect score matrices to artifacts and summaries to results",
    "retain only the accepted training-only activity definition",
    "retain the exact shared decoder selection and fit",
    "retain the accepted function-relevant molecular target",
    "retain intentionally configured one active inner fold",
    "retain requested metrics without historical duplicate predictions",
    "move each retained plot unchanged to the script producing its underlying data",
    "discarded PCA analysis", "discarded higher-rank experiment",
    "explicitly excluded from v2 active chain"))
fwrite(refactor_manifest, file.path(RESULT, "refactor_manifest.tsv"), sep = "\t")

## Stage 08 owns no scientific plot. It verifies and indexes the plots owned
## by the scripts that produced each underlying table or artifact.
plot_paths <- list.files(RESULTS_ROOT, pattern = "\\.pdf$", recursive = TRUE,
  full.names = TRUE)
plot_paths <- plot_paths[grepl("/0[1-7]_[^/]+/", plot_paths)]
pdf_pages <- function(path) {
  content <- suppressWarnings(system2("strings", path, stdout = TRUE,
    stderr = TRUE))
  as.integer(sum(grepl("/Type /Page([^s]|$)", content)))
}
owner_map <- c("01_inputs" = "scripts/01_prepare_inputs.R",
  "02_response_programs" = "scripts/02_learn_response_programs.R",
  "03_baseline_programs" = "scripts/03_learn_baseline_programs.R",
  "04_activity_panel" = "scripts/04_define_activity_panel.R",
  "05_program_function_bridge" = "scripts/05_fit_program_function_bridge.R",
  "06_baseline_bridge" = "scripts/06_fit_function_relevant_baseline_bridge.R",
  "07_end_to_end" = "scripts/07_evaluate_end_to_end.R")
figure_index <- data.table(
  stage = basename(dirname(plot_paths)),
  figure = basename(plot_paths),
  relative_path = substring(plot_paths, nchar(RESULTS_ROOT) + 2L),
  pages = vapply(plot_paths, pdf_pages, integer(1L)))
figure_index[, `:=`(
  producing_script = unname(owner_map[stage]),
  description = gsub("_", " ", sub("\\.pdf$", "", figure)),
  passed_one_page = pages == 1L)]
setorder(figure_index, stage, figure)
if (!nrow(figure_index) || anyNA(figure_index$pages) ||
    !all(figure_index$passed_one_page))
  fail("The stage-owned diagnostic suite is missing or contains a non-one-page PDF.")
fwrite(figure_index, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")

summary_lines <- c("# Validated Tahoe--PRISM v2 base", "",
  "All accepted Stage-1 benchmarks reproduced within their prespecified tolerances.",
  "Every surviving diagnostic is generated and stored by the stage that produces its underlying data; Stage 08 only indexes and verifies those figures.",
  "The response bank contains nine stable outcome-independent signed oNMF programs.",
  paste("The five DMSO/NTC baseline programs are retained as the benchmark;",
    "their stability and lineage/individual-line concentration limitations remain visible."),
  "", "The active molecular panel contains 356 conditions from 216 drugs; the activity-supported PRISM panel contains 95 conditions from 54 drugs and 646 held-out rows.",
  "", "Run: `Rscript tahoe_prism_v2/run_all.R`", "",
  "Large matrices, observation-level predictions, model objects and bootstrap draws are in `artifacts/`; directly shareable tables and figures are in `results/`.")
writeLines(summary_lines, file.path(RESULT, "README.md"))
writeLines("complete", file.path(RESULT, "stage1_reproduction_complete.flag"))
message(sprintf("Stage 08 complete: all reproduction gates passed; %d stage-owned plots indexed.",
  nrow(figure_index)))
