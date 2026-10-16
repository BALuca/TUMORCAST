#!/usr/bin/env Rscript
# Outcome-informed descriptive pilot. Gene admission uses retrospective PRISM
# and ALMANAC correlations plus held-out regression gain and loading variation.
# NMF fits full prediction-improvement gain times measured response sign;
# full-data gene screening makes this retrospective, not independent validation.
# Run from /work. The molecular-only 04 script/results remain a separate experiment.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(cluster)
  library(patchwork)
  library(RcppML)
  library(glmnet)
})
source("/work/src/program_pilot/config.R")
arguments <- commandArgs(trailingOnly = TRUE)
check(length(arguments) == 0L || identical(arguments, "--render-only"),
  "Supported invocation: Rscript 13_response_signed_functional_support_programs_nmf.R [--render-only]")
render_only <- identical(arguments, "--render-only")
nmf_tolerance <- 1e-4
nmf_max_iterations <- 100L
nmf_workers <- 8L # Distinct ranks run independently; each rank still has exactly one fit.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
local_resamples <- 0L # Single NMF fit per rank, without bootstraps or repeated starts.
local_minimum_rank <- 2L # This version selects the silhouette peak across all tested ranks.
local_rank_grid <- 2:20 # User-requested upper limit for this version.
functional_folds <- 5L # Hold whole drugs together; every dose remains a separate row.
functional_ridge_lambda <- 0.1 # Fixed mild ridge penalty on training-standardized predictor and outcome.
functional_input_transform <- "full cross-fitted gain times sign(50 - assay abundance); separate positive/negative channels; no clipping"
regression_gain_p_threshold <- 0.05 # One-sided nominal test: positive held-out gain across whole drugs.
minimum_loading_iqr <- 0.05 # Signed gain spread, in units of training outcome variance; robust to isolated outliers.
# This separately requested experiment keeps its outputs and gene rule local.
# Test each drug/dose row equally; pool doses without dose adjustment or averaging.
outcome_p_threshold <- 0.05 # Two-sided nominal threshold; BH values remain diagnostic only.
outcome_minimum_abs_correlation <- 0 # Earlier lenient rule has no effect-size floor.
local_gene_limit <- Inf # Retain every significant gene in the existing expressed coding universe.
output_dir <- file.path(work_dir, "results/program_pilot/13_response_signed_functional_support_programs_nmf")
data_dir <- file.path(output_dir, "data")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
for (path in c(output_dir, data_dir, figure_dir, model_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
set.seed(seed)
options(warn = 1)
theme_set(theme_classic(base_size = 10))
message("Outcome-selected displacement pilot; seed ", seed, "; results: ", output_dir)

source(file.path(code_dir, "functional_programs.R"))
# All eligible cells come from the explicit panel preparation, including assay-only lines.
prep_path <- file.path(work_dir, "results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs")
check(file.exists(prep_path), "Run prepare_functional_panel.R before the full panel")
panel <- qread(prep_path)
cell_map <- panel$cell_map
metadata <- copy(panel$metadata)
sensitivity <- copy(panel$sensitivity)
genes <- panel$genes
input_files <- c(gene_audit_path, panel$input_paths)
responses <- null_list <- baseline_list <- checks <- list()
for (id in cell_map$depmap_id) {
  molecular <- qread(panel$input_paths[match(id, cell_map$depmap_id)])
  rows <- which(metadata$depmap_id == id)
  response <- molecular$response[match(metadata$observation_id[rows], rownames(molecular$response)), , drop = FALSE]
  check(identical(colnames(response), genes) && all(is.finite(response)), "Invalid prepared panel response")
  responses[[id]] <- response
  null_list[[id]] <- molecular$null[[id]]
  baseline_list[[id]] <- molecular$baseline[[id]]
  checks[[id]] <- data.table(depmap_id=id, profiles=nrow(response), null_profiles=nrow(null_list[[id]]), finite=TRUE)
}
prepared <- list(response=do.call(rbind,responses), null=null_list, baseline=baseline_list, checks=rbindlist(checks))
metadata <- metadata[match(rownames(prepared$response), observation_id)]
rm(responses, molecular, response, null_list, baseline_list)
invisible(gc())
fwrite(prepared$checks, file.path(data_dir, "normalization_checks.tsv"), sep="\t")
fwrite(cell_map, file.path(data_dir, "panel_cell_mapping.tsv"), sep="\t")
fwrite(metadata, file.path(data_dir, "molecular_metadata.tsv.gz"), sep="\t")
fwrite(sensitivity, file.path(data_dir, "sensitivity_alignment.tsv"), sep="\t")
almanac_path <- file.path(work_dir, "almanac_conditionals/raw/ComboDrugGrowth_Nov2017.csv")
drug_mapping_path <- file.path(work_dir, "almanac_conditionals/repaired/results/drug_mapping.tsv")
prism_info_path <- file.path(cache_dir, "prism/secondary-screen-replicate-collapsed-treatment-info.csv")
prism_matrix_path <- file.path(cache_dir, "prism/secondary-screen-replicate-collapsed-logfold-change.csv")
function_mapping_path <- file.path(work_dir, "tahoe_prism_v2/artifacts/01_inputs/function_outcomes.tsv.gz")
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm = TRUE) / 25) * 25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm = TRUE) / 25) * 25)
message("Sensitivity bar maxima (% vehicle): ALMANAC ", almanac_bar_max, "; PRISM ", prism_bar_max)

# Select genes by two separate assay-specific outcome correlations. Each matched
# drug/dose observation has equal weight, with all doses pooled. Retain either-sign
# genes with nominal p <=0.05 in either assay, without a numerical cap.
# BH adjusted values are diagnostics only. Molecular variation/control metrics
# remain diagnostics; regression gain and signed-loading variation gate admission below.
line_genes <- list()
gene_selection <- list()
outcome_test_tables <- list()
outcome_test_summary <- list()
gene_selection_signatures <- list()
condition_checks <- list()
for (id in cell_map$depmap_id) {
  rows <- which(metadata$depmap_id == id)
  x <- prepared$response[rows, , drop = FALSE]
  null <- prepared$null[[id]]
  baseline <- prepared$baseline[[id]]
  noise_variance <- matrixStats::colVars(null)
  abundance_bin <- pmin(10L, ceiling(rank(baseline, ties.method = "first") / length(baseline) * 10))
  noise_floor <- ave(noise_variance, abundance_bin, FUN = median)
  noise <- pmax(noise_variance, noise_floor)
  check(all(noise > 0), "Nonpositive abundance-stratified control variance")
  null_threshold <- as.numeric(matrixStats::colQuantiles(abs(null), probs = 0.95))
  exceeds <- sweep(abs(x), 2L, null_threshold, ">")
  supported_drugs <- colSums(rowsum(exceeds * 1L, metadata$drug[rows], reorder = FALSE) > 0)
  variance <- matrixStats::colVars(x)
  focus <- which(metadata$focal[rows])
  focal_variance <- matrixStats::colVars(x[focus, , drop = FALSE])
  table <- data.table(depmap_id = id, recipient = metadata$recipient[rows[1]], gene = genes,
    response_variance = variance, control_variance = noise,
    excess_variance = pmax(variance - noise, 0), focal_excess_variance = pmax(focal_variance - noise, 0),
    supported_drugs = supported_drugs, focal_supported_conditions = colSums(exceeds[focus, , drop = FALSE]),
    abundance_bin = abundance_bin)
  # Keep original-scale variability visible without imposing the old SD/noise gates.
  table[, `:=`(response_sd = sqrt(response_variance), variance_ratio = response_variance / control_variance)]
  matched_outcomes <- sensitivity[match(metadata$observation_id[rows], observation_id)]
  check(identical(matched_outcomes$observation_id, metadata$observation_id[rows]),
    "Outcome rows do not align with gene displacement")
  tests <- rbindlist(lapply(c("PRISM", "ALMANAC"), function(assay_name) {
    outcome_column <- if (assay_name == "PRISM") "prism_abundance_pct" else "almanac_abundance_pct"
    result <- gene_outcome_tests(x, matched_outcomes[[outcome_column]], metadata$drug[rows])
    result[, `:=`(depmap_id = id, recipient = metadata$recipient[rows[1]], assay = assay_name)]
    result
  }))
  tests[, q_value := p.adjust(p_value, method = "BH", n = sum(!is.na(p_value))), by = assay]
  tests[, `:=`(significant = !is.na(p_value) & p_value <= outcome_p_threshold &
      abs(correlation) >= outcome_minimum_abs_correlation,
    outcome_scale = "abundance (% vehicle)", dose_handling = "all matched doses pooled; equal condition weight",
    multiplicity_family = "diagnostic BH across estimable genes within this cell and assay",
    selection_threshold = sprintf("nominal two-sided p <=%g and |r| >=%g in the same assay; BH diagnostic only",
      outcome_p_threshold, outcome_minimum_abs_correlation))]
  outcome_test_tables[[id]] <- tests
  outcome_test_summary[[id]] <- tests[, .(candidate_genes = .N, tested_genes = sum(test_status == "tested"),
    significant_genes = sum(significant), matched_conditions = matched_conditions[1],
    matched_drugs = matched_drugs[1], non_estimable_genes = sum(test_status != "tested")),
    by = .(depmap_id, recipient, assay)]
  prism_tests <- tests[assay == "PRISM", .(gene, prism_correlation = correlation,
    prism_p_value = p_value, prism_q_value = q_value, prism_outcome_eligible = significant)]
  almanac_tests <- tests[assay == "ALMANAC", .(gene, almanac_correlation = correlation,
    almanac_p_value = p_value, almanac_q_value = q_value, almanac_outcome_eligible = significant)]
  table <- merge(table, prism_tests, by = "gene", all.x = TRUE, sort = FALSE)
  table <- merge(table, almanac_tests, by = "gene", all.x = TRUE, sort = FALSE)
  table[, best_outcome_p_value := pmin(prism_p_value, almanac_p_value, na.rm = TRUE)]
  table[is.na(prism_p_value) & is.na(almanac_p_value), best_outcome_p_value := NA_real_]
  table[, best_outcome_q_value := pmin(prism_q_value, almanac_q_value, na.rm = TRUE)]
  table[is.na(prism_q_value) & is.na(almanac_q_value), best_outcome_q_value := NA_real_]
  table[, strongest_outcome_correlation := pmax(abs(prism_correlation), abs(almanac_correlation), na.rm = TRUE)]
  table[is.na(prism_correlation) & is.na(almanac_correlation), strongest_outcome_correlation := NA_real_]
  table[, outcome_eligible := prism_outcome_eligible | almanac_outcome_eligible]
  # Preserve qualifying-assay diagnostics; admission has no numerical gene cap.
  table[, best_qualifying_p_value := pmin(fifelse(prism_outcome_eligible, prism_p_value, NA_real_),
    fifelse(almanac_outcome_eligible, almanac_p_value, NA_real_), na.rm = TRUE)]
  table[, strongest_qualifying_correlation := pmax(fifelse(prism_outcome_eligible, abs(prism_correlation), NA_real_),
    fifelse(almanac_outcome_eligible, abs(almanac_correlation), NA_real_), na.rm = TRUE)]
  table[outcome_eligible == FALSE, `:=`(best_qualifying_p_value = NA_real_, strongest_qualifying_correlation = NA_real_)]
  chosen <- sort(table[outcome_eligible == TRUE, gene])
  check(length(chosen) >= 3L, paste("Fewer than three outcome-significant genes for", id,
    "at nominal p", outcome_p_threshold, "and |r| >=", outcome_minimum_abs_correlation,
    "-- no variance-selected substitute"))
  table[, selected := gene %chin% chosen]
  table[, selection_role := "outcome-significant in PRISM or ALMANAC; assays tested separately"]
  table[selected == FALSE, selection_role := "not selected"]
  table[, exclusion_reason := fcase(is.na(best_outcome_p_value), "no estimable outcome correlation",
    !outcome_eligible, "fails nominal p or absolute-correlation gate in both assays",
    selected, "selected", default = "invalid selection state")]
  check(identical(table$selected, table$outcome_eligible), "A significant gene was omitted from selection")
  line_genes[[id]] <- chosen
  gene_selection[[id]] <- table
  gene_selection_signatures[[id]] <- digest::digest(list(rule = "separate Pearson assays; either-sign nominal significance",
    p_threshold = outcome_p_threshold, minimum_abs_correlation = outcome_minimum_abs_correlation,
    cap = local_gene_limit, observations = metadata$observation_id[rows],
    outcomes = matched_outcomes[, .(observation_id, prism_abundance_pct, almanac_abundance_pct)],
    tests = tests[, .(gene, assay, correlation, p_value, q_value)], selected_genes = chosen), algo = "sha256")
  message(metadata$recipient[rows[1]], ": ", sum(table$outcome_eligible), " outcome-significant genes; retained ",
    length(chosen), "; PRISM/ALMANAC matched conditions ",
    paste(outcome_test_summary[[id]]$matched_conditions, collapse = "/"))
  null_energy <- mean(rowMeans(null^2))
  energy <- rowMeans(x^2)
  support <- energy / (energy + null_energy)
  metadata[rows, `:=`(condition_support = support, global_displacement = rowMeans(x),
    raw_response_rms = sqrt(energy), control_rms = sqrt(null_energy),
    response_control_rms_ratio = sqrt(energy / null_energy))]
  condition_checks[[id]] <- metadata[rows, .(observation_id, depmap_id, recipient, drug, dose, focal,
    n_plates, condition_support, raw_response_rms, control_rms, response_control_rms_ratio,
    global_displacement, contributing_plates)]
}
fwrite(rbindlist(gene_selection), file.path(data_dir, "gene_outcome_screen.tsv.gz"), sep = "\t")
selection_audit <- rbindlist(gene_selection)
fwrite(selection_audit[, .(genes = .N), by = .(recipient, exclusion_reason)],
  file.path(data_dir, "gene_filter_counts.tsv"), sep = "\t")
fwrite(rbindlist(outcome_test_tables), file.path(data_dir, "gene_outcome_correlations.tsv.gz"), sep = "\t")
fwrite(rbindlist(outcome_test_summary), file.path(data_dir, "gene_outcome_test_summary.tsv"), sep = "\t")
# Use a finite display floor for log10(0); every untransformed p-value is exported.
selection_audit[, displayed_minus_log10_p := -log10(pmax(best_outcome_p_value, .Machine$double.xmin))]
p <- ggplot(selection_audit, aes(displayed_minus_log10_p, fill = selected)) +
  geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.6, na.rm = TRUE) +
  geom_vline(xintercept = -log10(outcome_p_threshold), lty = 2) +
  facet_wrap(~recipient, ncol = 4) + scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#27818A")) +
  labs(x = "-log10(best outcome nominal p-value)", y = "Genes", fill = "Candidate",
    title = "Initial outcome screen: PRISM and ALMANAC tested separately",
    subtitle = sprintf("All doses pooled equally; nominal p <=%g AND |r| >=%g; regression admission follows",
      outcome_p_threshold, outcome_minimum_abs_correlation))
save_plot(p, "gene_outcome_filter", 11.5, 5.7)
fwrite(rbindlist(condition_checks), file.path(data_dir, "condition_signal_vs_control.tsv"), sep = "\t")

# Reuse validated saved banks on restart; never repeat a completed NMF rank grid.
# Discover within cell lines across drugs and doses. These are candidate
# banks; biological reproducibility is unavailable from aggregated PDEX profiles.
discoveries <- list()
regression_test_tables <- list()
for (id in cell_map$depmap_id) {
  rows <- which(metadata$depmap_id == id)
  scope <- paste0("cell_", gsub("[^A-Za-z0-9]+", "_", metadata$recipient[rows[1]]))
  scope_model_dir <- file.path(output_dir, scope, "models")
  dir.create(scope_model_dir, recursive = TRUE, showWarnings = FALSE)
  bank_path <- file.path(scope_model_dir, "program_bank.qs")
  scope_data_dir <- file.path(output_dir, scope, "data")
  dir.create(scope_data_dir, recursive = TRUE, showWarnings = FALSE)

  # Fit all Pearson-screen candidates once. Keep the pre-admission support so
  # restarting or adjusting this filter never requires repeating regressions.
  candidates <- line_genes[[id]]
  current_outcomes <- sensitivity[match(metadata$observation_id[rows], observation_id)]
  regression_signature <- digest::digest(list(
    candidate_screen = gene_selection_signatures[[id]],
    expression = digest::digest(prepared$response[rows, candidates, drop = FALSE], algo = "sha256"),
    folds = functional_folds, ridge_lambda = functional_ridge_lambda,
    glmnet_version = as.character(packageVersion("glmnet")), sensitivity_sign = TRUE), algo = "sha256")
  regression_path <- file.path(scope_model_dir, "functional_regression.qs")
  if (file.exists(regression_path)) {
    regression <- qread(regression_path)
    check(identical(regression$input_signature, regression_signature),
      paste("Cached regression inputs disagree for", scope))
    support <- regression$support
  } else {
    check(!render_only, paste("No saved pre-admission regression for", scope))
    message("Fitting regression admission for ", scope, "; ", length(candidates), " candidate genes")
    support <- functional_support(prepared$response[rows, candidates, drop = FALSE], metadata[rows],
      current_outcomes, allow_zero = TRUE, sensitivity_sign = TRUE)
    support$channels <- NULL # Exactly reconstructed after gene admission below.
    regression <- list(input_signature = regression_signature, support = support)
    qsave(regression, regression_path, preset = "high")
  }
  check(identical(rownames(support$loadings), candidates) &&
    identical(colnames(support$loadings), paste0(support$observations$assay, ":", support$observations$observation_id)) &&
    all(is.finite(support$loadings)), paste("Invalid pre-admission regression support for", scope))

  # Evaluate usefulness on the original, unsigned held-out prediction gain.
  # Each drug is one test unit; doses remain separate in fitting and clustering.
  # Require predictive significance and robust loading variation in the SAME assay.
  regression_tests <- rbindlist(lapply(unique(support$observations$assay), function(assay_name) {
    columns <- which(support$observations$assay == assay_name)
    gain <- support$raw_prediction_gain[, columns, drop = FALSE]
    loading <- support$loadings[, columns, drop = FALSE]
    drug <- support$observations$drug[columns]
    drug_sums <- rowsum(t(gain), drug, reorder = FALSE)
    drug_counts <- as.numeric(table(factor(drug, levels = rownames(drug_sums))))
    drug_gain <- t(sweep(drug_sums, 1L, drug_counts, "/"))
    mean_gain <- rowMeans(drug_gain)
    gain_se <- matrixStats::rowSds(drug_gain) / sqrt(ncol(drug_gain))
    estimable <- is.finite(gain_se) & gain_se > 0
    statistic <- p_value <- rep(NA_real_, length(candidates))
    statistic[estimable] <- mean_gain[estimable] / gain_se[estimable]
    p_value[estimable] <- stats::pt(statistic[estimable], df = ncol(drug_gain) - 1L, lower.tail = FALSE)
    loading_iqr <- matrixStats::rowIQRs(loading)
    result <- data.table(depmap_id = id, recipient = metadata$recipient[rows[1]],
      gene = candidates, assay = assay_name, matched_conditions = ncol(gain), matched_drugs = ncol(drug_gain),
      mean_gain_equal_condition = rowMeans(gain), mean_gain_equal_drug = mean_gain,
      drug_gain_standard_error = gain_se, statistic = statistic, p_value = p_value,
      test_status = fifelse(estimable, "tested", "constant drug-level gain: not estimable"),
      loading_iqr = loading_iqr, minimum_loading_iqr = minimum_loading_iqr)
    result[, q_value := p.adjust(p_value, method = "BH")]
    result[, `:=`(regression_significant = !is.na(p_value) & p_value <= regression_gain_p_threshold,
      loading_variation_eligible = loading_iqr >= minimum_loading_iqr)]
    result[, selected_in_assay := regression_significant & loading_variation_eligible]
    result
  }))
  admission <- regression_tests[, .(regression_significant = any(regression_significant),
    loading_variation_eligible = any(loading_variation_eligible),
    functional_eligible = any(selected_in_assay)), by = gene]
  chosen <- candidates[candidates %chin% admission[functional_eligible == TRUE, gene]]
  check(length(chosen) >= 3L, paste("Fewer than three genes pass held-out regression and loading variation for", scope,
    "-- no weak-gene substitute"))
  line_genes[[id]] <- chosen
  table <- merge(gene_selection[[id]], admission, by = "gene", all.x = TRUE, sort = FALSE)
  table[, selected := gene %chin% chosen]
  table[, selection_role := fifelse(selected,
    "outcome candidate; significant held-out regression gain and variable signed loading in the same assay", "not selected")]
  table[, exclusion_reason := fcase(is.na(best_outcome_p_value), "no estimable outcome correlation",
    !outcome_eligible, "fails nominal p or absolute-correlation gate in both assays",
    !regression_significant, "no significant positive held-out regression gain",
    !loading_variation_eligible, "signed loading IQR below threshold",
    !functional_eligible, "prediction significance and loading variation do not qualify in the same assay",
    selected, "selected", default = "invalid selection state")]
  gene_selection[[id]] <- table
  regression_test_tables[[id]] <- regression_tests
  fwrite(regression_tests, file.path(scope_data_dir, "gene_regression_admission.tsv.gz"), sep = "\t")
  dir.create(file.path(output_dir, scope, "figures"), recursive = TRUE, showWarnings = FALSE)
  gain_plot <- ggplot(regression_tests, aes(-log10(pmax(p_value, .Machine$double.xmin)), fill = selected_in_assay)) +
    geom_histogram(bins = 35L, position = "identity", alpha = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = -log10(regression_gain_p_threshold), lty = 2) + facet_wrap(~assay, ncol = 1L) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#27818A")) +
    labs(x = "-log10(held-out gain nominal p)", y = "Candidate genes", fill = "Pass both gates")
  variation_plot <- ggplot(regression_tests, aes(loading_iqr, fill = selected_in_assay)) +
    geom_histogram(bins = 35L, position = "identity", alpha = 0.6) +
    geom_vline(xintercept = minimum_loading_iqr, lty = 2) + facet_wrap(~assay, ncol = 1L) +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#27818A")) +
    labs(x = "Signed-loading IQR (normalized gain units)", y = "Candidate genes", fill = "Pass both gates")
  admission_plot <- gain_plot + variation_plot + plot_layout(guides = "collect") +
    plot_annotation(title = paste(gsub("_", " ", scope), "| regression gene admission"),
      subtitle = sprintf("Retain %d / %d genes; same-assay gain p <=%g and signed-loading IQR >=%g; either assay",
        length(chosen), length(candidates), regression_gain_p_threshold, minimum_loading_iqr))
  ggsave(file.path(output_dir, scope, "figures/regression_gene_filter.pdf"), admission_plot,
    width = 10, height = if (uniqueN(regression_tests$assay) == 2L) 6.2 else 4.2, device = cairo_pdf)
  ggsave(file.path(output_dir, scope, "figures/regression_gene_filter.png"), admission_plot,
    width = 10, height = if (uniqueN(regression_tests$assay) == 2L) 6.2 else 4.2, dpi = 160, bg = "white")
  message(scope, ": regression/loading admission retained ", length(chosen), " / ", length(candidates), " genes")
  gene_selection_signatures[[id]] <- digest::digest(list(regression_signature = regression_signature,
    p_threshold = regression_gain_p_threshold, minimum_loading_iqr = minimum_loading_iqr,
    rule = "positive drug-level held-out gain and variable signed loading; same assay; either assay",
    tests = regression_tests, selected_genes = chosen), algo = "sha256")

  # Preserve every retained value and coefficient. Only gene rows are removed.
  support$loadings <- support$loadings[chosen, , drop = FALSE]
  support$raw_prediction_gain <- support$raw_prediction_gain[chosen, , drop = FALSE]
  support$coefficients <- support$coefficients[gene %chin% chosen]
  support$channels <- cbind(pmax(support$loadings, 0), pmax(-support$loadings, 0))
  colnames(support$channels) <- c(paste0("positive:", colnames(support$loadings)),
    paste0("negative:", colnames(support$loadings)))
  support$settings$gene_admission <- list(p_threshold = regression_gain_p_threshold,
    minimum_loading_iqr = minimum_loading_iqr, test = "one-sided held-out drug-mean gain t test; same-assay IQR")

  bank <- if (file.exists(bank_path)) qread(bank_path) else NULL
  reuse_bank <- !is.null(bank) && identical(bank$gene_selection_signature, gene_selection_signatures[[id]])
  check(!render_only || reuse_bank, paste("Saved clustering predates current regression admission for", scope))
  if (reuse_bank) {
    bank <- restore_functional_bank(bank, prepared$response)
    check(identical(bank$genes, line_genes[[id]]) &&
      identical(bank$metadata$observation_id, metadata$observation_id[rows]) &&
      identical(bank$ranks$rank, as.integer(local_rank_grid[local_rank_grid <=
        min(length(bank$genes) - 1L, length(rows) - 1L)])) &&
      nrow(bank$resampling) == nrow(bank$ranks) * local_resamples &&
      identical(bank$functional_labels_used, TRUE) &&
      identical(bank$functional_labels_used_for_gene_selection, TRUE) &&
      identical(bank$rank_selection_uses_outcome_objective, TRUE) &&
      identical(bank$gene_selection_signature, gene_selection_signatures[[id]]),
      paste("Saved bank settings, outcomes or gene-admission rule disagree for", scope))
    # Check the actual saved molecular structure against current primary data,
    # not just filenames. Render-only never changes weights or selected ranks.
    current_correlation <- gene_correlation(prepared$response[rows, bank$genes, drop = FALSE],
      group_weights(metadata$drug[rows]) * metadata$condition_support[rows])
    check(max(abs(current_correlation - bank$correlation)) < 1e-10,
      paste("Saved molecular correlation differs from current inputs for", scope))
    current_selection <- select_silhouette(bank$ranks)
    check(identical(bank$ranks$selected, current_selection$selected) &&
      identical(bank$ranks$selection_status, current_selection$selection_status) &&
      identical(bank$ranks$selection_rule, current_selection$selection_rule) &&
      identical(bank$ranks$minimum_size_qc_pass, current_selection$minimum_size_qc_pass) &&
      identical(bank$ranks$silhouette_candidate, current_selection$silhouette_candidate),
      paste("Saved cut does not follow the maximum mean gene silhouette rule for", scope))
    channels <- support$channels
    check(identical(digest::digest(channels, algo = "sha256"), bank$nmf_input_signature) &&
      bank$nmf_settings$fits_per_rank == 1L && bank$nmf_settings$resamples == 0L &&
      bank$nmf_settings$tolerance == nmf_tolerance && bank$nmf_settings$max_iterations == nmf_max_iterations &&
      identical(bank$nmf_settings$package_version, as.character(packageVersion("RcppML"))) &&
      identical(bank$nmf_settings$input_transform, functional_input_transform) &&
      identical(bank$functional_support$settings, support$settings) &&
      max(abs(bank$clustering_correlation - gene_correlation(t(support$loadings), rep(1, ncol(support$loadings))))) < 1e-10,
      paste("Saved single-fit NMF settings or functional loading inputs differ for", scope))
    check(all(bank$ranks$nmf_fits == 1L) &&
      identical(as.integer(bank$ranks$nmf_seed),
        as.integer(seed + 100000L * match(id, cell_map$depmap_id) + bank$ranks$rank)),
      paste("Saved NMF seeds or initialization counts differ for", scope))
    message("Rendering saved molecular analysis: ", scope)
  } else {
    if (!is.null(bank)) message("Replacing clustering fitted before the current regression-admission rule: ", scope)
    message("Discovering ", scope, " over ", length(rows), " conditions")
    bank <- discover_local(prepared$response[rows, , drop = FALSE], prepared$null[[id]],
      metadata[rows], line_genes[[id]], "drug", scope, "cell", support = support, sensitivity_sign = TRUE)
    bank$gene_selection_signature <- gene_selection_signatures[[id]]
    bank$gene_selection_rule <- sprintf(
      paste0("Pearson candidate p <=%g; held-out positive drug-mean gain p <=%g and signed-loading IQR >=%g",
        " in the same assay; either assay; no gene cap"),
      outcome_p_threshold, regression_gain_p_threshold, minimum_loading_iqr)
    save_functional_bank(bank, bank_path)
  }
  if (is.null(bank$derived_matrix_signatures)) save_functional_bank(bank, bank_path)
  discoveries[[scope]] <- bank[c("scope", "orientation", "metadata", "genes", "weights", "scores", "ranks", "resampling", "diagnostics", "programs")]
  render_diagnostics(bank)
  if (is.null(bank$weights)) {
    message("Program map withheld for ", bank$scope, ": ", bank$diagnostics$selection_status)
    render_unresolved(bank)
  } else {
    message("Rendering ", bank$scope)
    render_local(bank)
  }
  render_functional_loadings(bank)
  writeLines(c(readLines(file.path(output_dir, scope, "README.txt")),
    sprintf("Additional regression admission: one-sided held-out drug-mean gain p <=%g and signed-loading IQR >=%g in the same assay.",
      regression_gain_p_threshold, minimum_loading_iqr),
    "IQR units: squared-error gain divided by training outcome variance. BH values are diagnostic only.",
    "data/gene_regression_admission.tsv.gz audits every regression candidate and assay."),
    file.path(output_dir, scope, "README.txt"))
  message("Figures ready: ", file.path(output_dir, bank$scope, "figures"))
  rm(support, regression_tests, admission, regression, bank)
  invisible(gc())
}

fwrite(rbindlist(gene_selection), file.path(data_dir, "gene_selection.tsv.gz"), sep = "\t")
fwrite(rbindlist(gene_selection)[, .(genes = .N), by = .(recipient, exclusion_reason)],
  file.path(data_dir, "gene_filter_counts.tsv"), sep = "\t")
fwrite(rbindlist(regression_test_tables), file.path(data_dir, "gene_regression_admission.tsv.gz"), sep = "\t")
fwrite(rbindlist(lapply(discoveries, `[[`, "ranks")), file.path(data_dir, "rank_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(lapply(discoveries, `[[`, "resampling")), file.path(data_dir, "resampling_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(lapply(discoveries, `[[`, "diagnostics")), file.path(data_dir, "discovery_diagnostics.tsv"), sep = "\t")
program_summary <- rbindlist(lapply(discoveries, `[[`, "programs"))
if (!ncol(program_summary)) program_summary <- data.table(scope = character(), orientation = character(),
  program = character(), seed_genes = integer(), effective_weighted_genes = numeric(),
  mean_gene_correlation = numeric(), mean_control_gene_correlation = numeric(), response_rms = numeric(),
  control_rms = numeric(), response_control_rms_ratio = numeric(), marker_genes = character())
fwrite(program_summary, file.path(data_dir, "program_diagnostics.tsv"), sep = "\t")
freeze <- data.table(scope = names(discoveries), functional_labels_used = TRUE,
  functional_labels_used_for_gene_selection = TRUE, rank_selection_uses_outcome_objective = TRUE,
  bank_path = file.path(output_dir, names(discoveries), "models/program_bank.qs"))
freeze[, sha256 := vapply(bank_path, digest::digest, character(1), file = TRUE, algo = "sha256")]
fwrite(freeze, file.path(data_dir, "molecular_bank_freeze.tsv"), sep = "\t")
message("All ", nrow(cell_map), " cell-line fits and their figures are ready")

# Raw RNA signal and functional abundance are shown on separate descriptive
# scales. The panel is a complete 7 x 7 x 3 view, with no outcome-based sorting.
signal <- metadata[focal == TRUE, .(observation_id, raw_response_rms, response_control_rms_ratio)]
signal <- merge(sensitivity[focal == TRUE], signal, by = "observation_id", all.x = TRUE, sort = FALSE)
signal[, drug := factor(drug, levels = local_drugs)]
p <- ggplot(signal, aes(dose, response_control_rms_ratio, color = recipient, group = recipient)) +
  geom_hline(yintercept = 1, lty = 2, color = "grey60") + geom_line() + geom_point() +
  scale_x_log10(breaks = c(0.05, 0.5, 5)) + facet_wrap(~drug, ncol = 4) +
  labs(x = "Tahoe dose (uM)", y = "Whole-transcriptome RNA RMS / plate-control RMS",
    title = "Molecular support across the complete focal panel", color = "Recipient")
save_plot(p, "panel_molecular_support", 11.5, 5.7)
p <- ggplot(signal, aes(dose, prism_abundance_pct, color = recipient, group = recipient)) +
  geom_hline(yintercept = 100, lty = 2, color = "grey60") + geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  scale_x_log10(breaks = c(0.05, 0.5, 5)) + facet_wrap(~drug, ncol = 4) +
  labs(x = "Tahoe query dose (uM)", y = "PRISM abundance (% vehicle)",
    title = "Sensitivity annotations at Tahoe doses", subtitle = "Bracketed log-dose interpolation; no extrapolation",
    color = "Recipient")
save_plot(p, "panel_sensitivity", 11.5, 5.7)

# Program/null coherence exposes shared-control artifacts rather than claiming
# that a clean diagonal establishes biological validity.
programs <- rbindlist(lapply(discoveries, `[[`, "programs"))
if (nrow(programs) >= 3L) {
  association <- cor(programs$mean_control_gene_correlation, programs$mean_gene_correlation)
  p <- ggplot(programs, aes(mean_control_gene_correlation, mean_gene_correlation, color = scope)) +
    geom_point(alpha = 0.75) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    theme(aspect.ratio = 1) + labs(x = "Within-module correlation in vehicle pseudo-responses",
      y = "Within-module correlation in treated responses", color = "Cell-line bank",
      title = sprintf("Control structure versus program coherence | Pearson r=%.2f", association))
  save_plot(p, "program_control_structure", 5.8, 4.8)
} else {
  message("Program/control association figure unavailable: fewer than three selected programs")
}

# Record exploratory dose trends and plate confounding for original program
# scores. With three doses these correlations are descriptive, not tests.
dose_rows <- list()
plate_rows <- list()
for (bank in discoveries) {
  if (is.null(bank$weights)) next
  for (j in seq_len(ncol(bank$scores))) {
    md <- copy(bank$metadata)
    md[, `:=`(value = bank$scores[, j], scope = bank$scope, program = colnames(bank$scores)[j])]
    dose_rows[[length(dose_rows) + 1L]] <- md[, .(
      n_doses = uniqueN(dose), log_dose_correlation = if (uniqueN(dose) >= 3 && sd(value) > 0)
        cor(log(dose), value) else NA_real_, minimum_score = min(value), maximum_score = max(value)),
      by = .(scope, program, recipient, drug)]
    md[, primary_plate := sub(",.*", "", contributing_plates)]
    fit <- lm(value ~ factor(primary_plate), data = md)
    plate_rows[[length(plate_rows) + 1L]] <- data.table(scope = bank$scope,
      program = colnames(bank$scores)[j], primary_plate_r_squared = summary(fit)$r.squared,
      interpretation = "Drug/dose and plate are confounded; association is not a causal batch estimate")
  }
}
dose_summary <- if (length(dose_rows)) rbindlist(dose_rows) else data.table(scope = character(), program = character(),
  recipient = character(), drug = character(), n_doses = integer(), log_dose_correlation = numeric(),
  minimum_score = numeric(), maximum_score = numeric())
plate_summary <- if (length(plate_rows)) rbindlist(plate_rows) else data.table(scope = character(), program = character(),
  primary_plate_r_squared = numeric(), interpretation = character())
fwrite(dose_summary, file.path(data_dir, "program_dose_trends.tsv"), sep = "\t")
fwrite(plate_summary, file.path(data_dir, "program_plate_associations.tsv"), sep = "\t")

# The index points to each cell's own outputs; shared tables summarize the run.
artifact_index <- rbindlist(lapply(discoveries, function(bank) {
  artifacts <- c("silhouette.pdf", "cophenetic.pdf", "diagnostics.pdf", "gene_correlation.pdf",
    "displacement_overview.pdf", "displacement_pages.pdf", "program_scores.pdf", "focal_panel.pdf")
  if (is.null(bank$weights)) artifacts <- c("silhouette.pdf", "cophenetic.pdf", "diagnostics.pdf",
    "gene_correlation.pdf", "displacement_unresolved.pdf", "focal_panel.pdf")
  artifacts <- c(artifacts, "functional_loadings.pdf", "functional_loadings_focal.pdf", "regression_gene_filter.pdf")
  data.table(recipient = bank$metadata$recipient[1], scope = bank$scope,
    selected_rank = bank$diagnostics$selected_rank, selection_status = bank$diagnostics$selection_status,
    genes = length(bank$genes),
    artifact = artifacts, path = file.path(output_dir, bank$scope, "figures", artifacts))
}))
fwrite(artifact_index, file.path(output_dir, "artifact_index.tsv"), sep = "\t")
writeLines(c(sprintf("Cell-line-first Tahoe displacement pilot: %d drugs, %d profiles, %d recipients.", uniqueN(metadata$drug), nrow(metadata), nrow(cell_map)),
  "Creating script: /work/src/program_pilot/13_response_signed_functional_support_programs_nmf.R",
  "Each cell_* folder owns its models, data and figures. artifact_index.tsv lists figure paths.",
  "Each cell-line bank is saved and its figures rendered before the next cell line; validated banks are reused on restart.",
  sprintf("Gene admission: separate PRISM/ALMANAC tests; p <=%g AND |r| >=%g in the same assay, either sign.",
    outcome_p_threshold, outcome_minimum_abs_correlation),
  "All matched drug/dose rows pooled with equal weight; no dose averaging, splitting or adjustment.",
  "BH adjusted values are diagnostic only; retain every gene passing both admission stages without a numerical cap.",
  sprintf("Regression admission: held-out drug-mean gain one-sided p <=%g AND signed-loading IQR >=%g in the same assay; either assay.",
    regression_gain_p_threshold, minimum_loading_iqr),
  "Gain tests use whole drugs as units; regression and clustering retain every dose separately.",
  "data/gene_regression_admission.tsv.gz exports gains, nominal p, diagnostic BH q, loading IQR and same-assay decisions.",
  "No displacement SD, variance-ratio, drug-recurrence or abundance-bin gate is used.",
  "Rank: maximum positive mean gene silhouette among cuts with all seeds >=12 genes; exact ties choose smaller rank.",
  "All ranks 2--20 are eligible, subject to size QC; retain the complete diagnostic curve.",
  "silhouette.pdf/png: selection curve; cophenetic.pdf/png: diagnostic only; diagnostics.pdf/png: four panels.",
  "displacement_pages.pdf: readable pages with sensitivity bars where observed/bracketed.",
  "Within each program: PRISM first, ALMANAC-only second, no data last; ascending % vehicle (most sensitive first).",
  "Sensitivity: purple <=50% vehicle; green >50%; grey cross means no value. Measured zero stays purple.",
  "gene_correlation.pdf/png: all retained genes in block order; off-block values are preserved.",
  "data/gene_selection.tsv.gz includes every candidate gene and its exclusion reason.",
  "One NMF fit per rank on gene-by-condition cross-fitted prediction-improvement profiles; separate signed assay channels.",
  "One ridge regression per gene/assay; five deterministic whole-drug folds; training-only standardization; lambda=0.1.",
  "Loading=original variance-scaled squared-error gain * sign(50 - measured abundance); all negative gains retained.",
  "Missing functional outcomes are absent from the NMF input and grey in the loading maps; no imputation.",
  "functional_loadings.pdf/png and functional_loadings_focal.pdf/png match the corresponding displacement order exactly.",
  "Silhouette uses Pearson distance on functional loadings; select its peak; signed displacement scores are unchanged.",
  "Consensus/cophenetic, PAC, resampling ARI and the second-fit centering diagnostic are unavailable; no repeats or subsampling.",
  "Gene seed means, overlapping weights, signed scores and sensitivity matching/sorting retain the original pipeline.",
  "Outcomes inform gene admission and functional loadings. Full-data gene screening makes these retrospective maps.",
  "Unsupported selections have raw annotated maps but no program dictionary.",
  "The original 04 and 05 experiments remain unchanged in ../04_mockup_local_programs and ../05_mockup_outcome_programs."),
  file.path(output_dir, "README.txt"))

# License provenance is an inventory, not a commercial-use clearance. No new
# resource was downloaded. R package copyleft obligations are recorded.
packages <- c("data.table", "qs", "ComplexHeatmap", "circlize", "ggplot2", "cluster", "clue", "patchwork", "matrixStats", "RcppML", "glmnet")
package_inventory <- rbindlist(lapply(packages, function(package) {
  description <- packageDescription(package)
  data.table(resource = package, version = description$Version, license = description$License,
    provenance = find.package(package), commercial_review = if (grepl("GPL", description$License))
      "Copyleft: review distribution obligations" else "Review license terms before distribution")
}))
fwrite(package_inventory, file.path(data_dir, "package_licenses.tsv"), sep = "\t")
writeLines(c("Tahoe PDEX, PRISM secondary screen, ALMANAC: existing cached research resources.",
  "No commercial-use clearance was established by this exploratory run.",
  "Resource origins and normalization are retained in the input manifests and source README files."),
  file.path(output_dir, "resource_provenance.txt"))
write_manifest("13_response_signed_functional_support_programs_nmf", c(input_files, prep_path, almanac_path, drug_mapping_path,
  prism_info_path, prism_matrix_path, function_mapping_path), list(drugs = sort(unique(metadata$drug)),
    focal_supplement = local_drugs, recipients = cell_map$depmap_name,
    gene_limit_per_line = "none; retain all genes passing outcome and regression/loading admission",
    gene_selection = "Pearson candidates then positive held-out drug-mean gain and signed-loading IQR in the same assay; either assay; BH diagnostic only",
    outcome_p_threshold = outcome_p_threshold, outcome_minimum_abs_correlation = outcome_minimum_abs_correlation,
    regression_gain_p_threshold = regression_gain_p_threshold, minimum_loading_iqr = minimum_loading_iqr,
    regression_admission_test = "one-sided one-sample t test of whole-drug held-out mean gain against zero",
    loading_variation_scale = "IQR of signed prediction-improvement gain normalized by training outcome variance",
    gene_cap_order = "none; every qualifying gene retained, ordered by gene symbol",
    outcome_correlation_dose_handling = "all matched drug/dose rows pooled with equal condition weight; no dose adjustment",
    resamples = local_resamples, ranks = local_rank_grid, programs = "single-fit gene/signature NMF clusters; unchanged overlapping gene-seed weights",
    nmf_fits_per_rank = 1L, nmf_tolerance = nmf_tolerance, nmf_max_iterations = nmf_max_iterations,
    nmf_package_version = as.character(packageVersion("RcppML")),
    nmf_input = functional_input_transform, functional_folds = functional_folds,
    functional_ridge_lambda = functional_ridge_lambda, glmnet_version = as.character(packageVersion("glmnet")),
    functional_loss_scaling = "training outcome variance; equal dose-row weight; separate assays; multiply full gain by sign(50 - measured abundance)",
    rank_distance = "same Pearson distance, evaluated on functional loading profiles",
    rank_rule = "maximum positive mean gene silhouette on functional loadings; size-qualified cut; smaller rank breaks exact ties",
    minimum_selectable_rank = local_minimum_rank, minimum_seed_genes = local_minimum_genes,
    sensitivity = "exact ALMANAC; PRISM bracketing log-dose/log-abundance interpolation; no extrapolation",
    within_program_display_order = "PRISM ascending % vehicle, ALMANAC-only ascending % vehicle, then no data; ALMANAC breaks PRISM ties",
    render_only = render_only,
    frozen_banks = freeze[, .(scope, sha256)],
    scientific_status = "retrospective outcome-informed exploratory programs; nominal screening; no independent evaluation"))
message("Complete: ", nrow(cell_map), " cell-line fits, full-panel maps, silhouette diagnostics and gene-outcome tables")
fwrite(metadata, file.path(data_dir, "molecular_metadata.tsv.gz"), sep="\t")
fwrite(data.table(depmap_id=cell_map$depmap_id, recipient=cell_map$depmap_name, status="complete"),
  file.path(data_dir, "panel_completed.tsv"), sep="\t")
