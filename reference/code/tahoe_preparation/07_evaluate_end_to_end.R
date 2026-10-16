#!/usr/bin/env Rscript

## Reconstruct held-out function and report pooled, within, strict,
## drug-balanced, amplitude, headroom and drug-block bootstrap metrics.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
RESULT <- stage_result("07_end_to_end")
ARTIFACT <- stage_artifact("07_end_to_end")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

plot_stage07 <- function() {
  theme_stage07 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
      plot.title = element_text(face = "bold"))
  square_scatter_theme <- function(base_size = 8.2)
    theme_stage07(base_size) + theme(aspect.ratio = 1)
  save_gg <- function(filename, plot, width, height) ggsave(
    file.path(RESULT, filename), plot, width = width, height = height,
    device = cairo_pdf, bg = "white", limitsize = FALSE)
  fmt <- function(x, digits = 3L) ifelse(is.finite(x),
    formatC(x, digits = digits, format = "f"), "NA")
  source_labels <- c(tahoe_dmso = "Tahoe DMSO", depmap_static = "DepMap")
  source_colors <- c(tahoe_dmso = "#B2182B", depmap_static = "#2166AC")

  performance <- fread(file.path(RESULT, "end_to_end_performance.tsv"))
  activity_performance <- performance[cohort == "activity_supported" &
    metric %chin% c("within_perturbation",
      "strict_cell_line_by_perturbation")]
  activity_performance[, display_model := fcase(
    prediction_model == "perturbation_average", "Perturbation average",
    prediction_model == "complete_measured", "Complete measured response",
    prediction_model == "baseline_conditioned" &
      baseline_source == "tahoe_dmso", "Tahoe DMSO baseline",
    prediction_model == "baseline_conditioned" &
      baseline_source == "depmap_static", "DepMap baseline")]
  activity_performance[, metric_label := fifelse(metric ==
    "within_perturbation", "Within perturbation",
    "Strict line + perturbation")]
  p_end <- ggplot(activity_performance, aes(x = pearson,
    y = reorder(display_model, pearson, FUN = max, na.rm = TRUE),
    color = metric_label, shape = metric_label)) +
    geom_vline(xintercept = 0, color = "grey75") +
    geom_point(size = 2) +
    labs(title = "End-to-end functional context correlation",
      x = "Pearson correlation", y = NULL, color = NULL, shape = NULL) +
    theme_stage07()
  save_gg("07a_end_to_end_context_correlation.pdf", p_end, 5.2, 4.3)

  headroom <- fread(file.path(RESULT,
    "functional_headroom_recovered.tsv"))[cohort == "activity_supported"]
  headroom[, display := source_labels[baseline_source]]
  p_head <- ggplot(headroom, aes(x = headroom_recovered,
    y = reorder(display, headroom_recovered), color = baseline_source)) +
    geom_vline(xintercept = c(0, 0.10, 0.25, 0.50),
      color = c("grey50", "grey85", "grey85", "grey85"),
      linetype = c("solid", "dashed", "dashed", "dashed")) +
    geom_point(size = 2.2) +
    scale_color_manual(values = source_colors, guide = "none") +
    labs(title = "Measured functional headroom recovered",
      x = "Error-reduction fraction", y = NULL) + theme_stage07()
  save_gg("07b_functional_headroom_recovered.pdf", p_head, 4.8, 3.3)

  function_predictions <- fread(file.path(stage_artifact(
    "06_baseline_bridge"), "baseline_function_predictions.tsv.gz"))
  activity_ids <- fread(file.path(stage_result("04_activity_panel"),
    "condition_activity_census.tsv"))[activity_supported == TRUE,
      perturbation_id]
  metric_labels <- c(pooled = "Pooled", within_perturbation =
    "Within perturbation", strict_cell_line_by_perturbation =
    "Strict cell line by perturbation")
  intervals <- fread(file.path(RESULT,
    "drug_block_bootstrap_intervals.tsv"))
  for (source_name in names(source_labels)) for (metric in names(metric_labels)) {
    metric_name <- metric
    current <- function_predictions[baseline_source == source_name &
      evaluation_included == TRUE & perturbation_id %chin% activity_ids]
    residual <- function_residual_table(current, metric_name)
    interval <- intervals[prediction_model == "baseline_conditioned" &
      baseline_source == source_name & metric == metric_name]
    label <- sprintf("r=%s [%s, %s]\nobs SD=%s; pred SD=%s\nn=%s",
      fmt(safe_cor(residual$observed_residual,
        residual$predicted_residual)), fmt(interval$pearson_ci_lower[1L]),
      fmt(interval$pearson_ci_upper[1L]),
      fmt(sd(residual$observed_residual)),
      fmt(sd(residual$predicted_residual)),
      format(nrow(residual), big.mark = ","))
    plot <- ggplot(residual,
      aes(observed_residual, predicted_residual)) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
      geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
      geom_point(alpha = 0.22, size = 0.45,
        color = source_colors[[source_name]]) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
        linewidth = 0.55, color = "black") +
      annotate("text", x = -Inf, y = Inf, label = label,
        hjust = -0.05, vjust = 1.1, size = 2.15) +
      labs(x = "Observed function", y = "Predicted function",
        title = paste(source_labels[[source_name]],
          metric_labels[[metric_name]], sep = " — ")) + square_scatter_theme(7.5)
    save_gg(sprintf("07c_%s_%s_observed_vs_predicted.pdf",
      source_name, metric_name), plot, 4.4, 4.2)
  }

  per_drug <- fread(file.path(RESULT, "end_to_end_per_drug.tsv"))
  per_drug[, source_label := source_labels[baseline_source]]
  per_drug_long <- melt(per_drug,
    id.vars = c("baseline_source", "source_label", "drug"),
    measure.vars = c("pooled_pearson", "within_pearson", "strict_pearson"),
    variable.name = "metric", value.name = "pearson")
  per_drug_long[, metric := factor(metric,
    levels = c("pooled_pearson", "within_pearson", "strict_pearson"),
    labels = c("Pooled", "Within perturbation", "Strict"))]
  p_per_drug <- ggplot(per_drug_long, aes(metric, pearson,
      color = source_label)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey60") +
    geom_boxplot(aes(group = interaction(metric, source_label)),
      position = position_dodge(width = 0.65), width = 0.55,
      outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.16,
      dodge.width = 0.65), alpha = 0.28, size = 0.55) +
    scale_color_manual(values = c("Tahoe DMSO" = "#B2182B",
      "DepMap" = "#2166AC")) +
    labs(x = NULL, y = "Per-drug Pearson correlation", color = NULL,
      title = "End-to-end performance across drugs") + theme_stage07()
  save_gg("07d_end_to_end_per_drug_performance.pdf", p_per_drug, 6.2, 4.4)

  p_amplitude <- ggplot(per_drug,
    aes(source_label, predicted_to_observed_sd_ratio,
      color = source_label)) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35,
      color = "grey45") +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35) +
    geom_jitter(width = 0.14, alpha = 0.32, size = 0.6) +
    scale_color_manual(values = c("Tahoe DMSO" = "#B2182B",
      "DepMap" = "#2166AC"), guide = "none") +
    labs(x = NULL, y = "Predicted / observed within-drug SD",
      title = "End-to-end prediction amplitude by drug") + theme_stage07()
  save_gg("07e_end_to_end_per_drug_prediction_amplitude.pdf", p_amplitude,
    5.2, 4.3)

  interval_plot <- intervals[metric %chin% c("within_perturbation",
    "strict_cell_line_by_perturbation")]
  interval_plot[, display_model := fcase(
    prediction_model == "perturbation_average", "Perturbation average",
    prediction_model == "complete_measured", "Complete measured response",
    prediction_model == "baseline_conditioned" &
      baseline_source == "tahoe_dmso", "Tahoe DMSO baseline",
    prediction_model == "baseline_conditioned" &
      baseline_source == "depmap_static", "DepMap baseline")]
  interval_plot <- merge(interval_plot,
    performance[cohort == "activity_supported",
      .(prediction_model, baseline_source, metric, pearson)],
    by = c("prediction_model", "baseline_source", "metric"), all.x = TRUE)
  interval_plot[, metric_label := fifelse(metric == "within_perturbation",
    "Within perturbation", "Strict")]
  p_intervals <- ggplot(interval_plot,
    aes(pearson, reorder(display_model, pearson), color = metric_label)) +
    geom_vline(xintercept = 0, color = "grey75") +
    geom_errorbarh(aes(xmin = pearson_ci_lower, xmax = pearson_ci_upper),
      height = 0.18, linewidth = 0.45) +
    geom_point(size = 1.8) +
    labs(x = "Pearson correlation with 95% drug-block interval", y = NULL,
      color = NULL, title = "Held-out functional uncertainty") +
    theme_stage07()
  save_gg("07f_end_to_end_bootstrap_intervals.pdf", p_intervals, 6.4, 4.6)
}

if ("--plots-only" %in% commandArgs(trailingOnly = TRUE)) {
  plot_stage07()
  message("Stage 07 plots refreshed from existing end-to-end results.")
  quit(save = "no", status = 0L)
}

measured <- fread(file.path(stage_artifact("05_program_function_bridge"),
  "measured_function_predictions.tsv.gz"))
baseline <- fread(file.path(stage_artifact("06_baseline_bridge"),
  "baseline_function_predictions.tsv.gz"))
activity_ids <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))[activity_supported == TRUE,
  perturbation_id]

common <- measured[, .(observation_id, depmap_id, cell_line,
  perturbation_id, drug, dose, observed_function, evaluation_included,
  mean_prediction, measured_prediction)]
predictions <- rbindlist(list(
  common[, .(observation_id, depmap_id, cell_line, perturbation_id, drug,
    dose, observed_function, evaluation_included,
    prediction_model = "perturbation_average",
    baseline_source = NA_character_, predicted_function = mean_prediction)],
  common[, .(observation_id, depmap_id, cell_line, perturbation_id, drug,
    dose, observed_function, evaluation_included,
    prediction_model = "complete_measured",
    baseline_source = NA_character_, predicted_function = measured_prediction)],
  baseline[, .(observation_id, depmap_id, cell_line, perturbation_id, drug,
    dose, observed_function, evaluation_included,
    prediction_model = "baseline_conditioned", baseline_source,
    predicted_function)]), use.names = TRUE)
predictions[, activity_supported := perturbation_id %chin% activity_ids]

performance_parts <- list()
for (cohort in c("complete_matched_panel", "activity_supported")) {
  cohort_rows <- predictions[evaluation_included == TRUE]
  if (cohort == "activity_supported") cohort_rows <- cohort_rows[
    activity_supported == TRUE]
  for (model in unique(cohort_rows[, paste(prediction_model,
      baseline_source, sep = "__")])) {
    current <- cohort_rows[paste(prediction_model, baseline_source,
      sep = "__") == model]
    for (metric in c("pooled", "within_perturbation",
        "strict_cell_line_by_perturbation")) {
      residual <- function_residual_table(current, metric)
      w <- balanced_weights(residual$drug)
      error <- residual$observed_residual - residual$predicted_residual
      performance_parts[[length(performance_parts) + 1L]] <- data.table(
        cohort = cohort, prediction_model = current$prediction_model[1L],
        baseline_source = current$baseline_source[1L], metric = metric,
        n_rows = nrow(residual), n_conditions = uniqueN(residual$perturbation_id),
        n_drugs = uniqueN(residual$drug), n_lines = uniqueN(residual$depmap_id),
        pearson = safe_cor(residual$observed_residual,
          residual$predicted_residual),
        drug_balanced_pearson = weighted_cor(residual$observed_residual,
          residual$predicted_residual, w),
        rmse = sqrt(mean(error ^ 2)), mae = mean(abs(error)),
        drug_balanced_rmse = sqrt(sum(w * error ^ 2) / sum(w)),
        observed_residual_sd = sd(residual$observed_residual),
        predicted_residual_sd = sd(residual$predicted_residual),
        predicted_to_observed_sd_ratio = sd(residual$predicted_residual) /
          sd(residual$observed_residual),
        calibration_slope = if (sd(residual$predicted_residual) > 0)
          coef(lm(observed_residual ~ predicted_residual,
            data = residual))[2L] else NA_real_, correlation_method = "pearson")
    }
  }
}
performance <- rbindlist(performance_parts, fill = TRUE)
fwrite(performance, file.path(RESULT, "end_to_end_performance.tsv"), sep = "\t")

## Fraction of the complete-measured functional headroom recovered.
headroom_parts <- list()
for (cohort in c("complete_matched_panel", "activity_supported")) {
  keep_ids <- if (cohort == "activity_supported") activity_ids else
    unique(common$perturbation_id)
  base_rows <- common[evaluation_included == TRUE &
    perturbation_id %chin% keep_ids]
  for (source in unique(baseline$baseline_source)) {
    conditioned <- baseline[baseline_source == source &
      evaluation_included == TRUE & perturbation_id %chin% keep_ids]
    joined <- merge(base_rows, conditioned[, .(observation_id,
      baseline_prediction = predicted_function)], by = "observation_id",
      all = FALSE)
    sse_mean <- sum((joined$observed_function - joined$mean_prediction) ^ 2)
    sse_baseline <- sum((joined$observed_function -
      joined$baseline_prediction) ^ 2)
    sse_measured <- sum((joined$observed_function -
      joined$measured_prediction) ^ 2)
    w <- balanced_weights(joined$drug)
    wsse_mean <- sum(w * (joined$observed_function -
      joined$mean_prediction) ^ 2)
    wsse_baseline <- sum(w * (joined$observed_function -
      joined$baseline_prediction) ^ 2)
    wsse_measured <- sum(w * (joined$observed_function -
      joined$measured_prediction) ^ 2)
    headroom_parts[[length(headroom_parts) + 1L]] <- data.table(
      cohort = cohort, baseline_source = source, n_rows = nrow(joined),
      sse_mean = sse_mean, sse_baseline = sse_baseline,
      sse_measured = sse_measured,
      headroom_recovered = (sse_mean - sse_baseline) /
        (sse_mean - sse_measured),
      drug_balanced_sse_mean = wsse_mean,
      drug_balanced_sse_baseline = wsse_baseline,
      drug_balanced_sse_measured = wsse_measured,
      drug_balanced_headroom_recovered = (wsse_mean - wsse_baseline) /
        (wsse_mean - wsse_measured))
  }
}
headroom <- rbindlist(headroom_parts)
fwrite(headroom, file.path(RESULT, "functional_headroom_recovered.tsv"),
  sep = "\t")

## Per-drug context performance and prediction amplitude.
per_drug <- rbindlist(lapply(unique(baseline$baseline_source), function(source) {
  current <- baseline[baseline_source == source & evaluation_included == TRUE &
    perturbation_id %chin% activity_ids]
  rbindlist(lapply(unique(current$drug), function(drug_name) {
    x <- current[drug == drug_name]
    within <- function_residual_table(x, "within_perturbation")
    strict <- function_residual_table(x, "strict_cell_line_by_perturbation")
    data.table(baseline_source = source, drug = drug_name,
      n_rows = nrow(x), n_conditions = uniqueN(x$perturbation_id),
      pooled_pearson = safe_cor(x$observed_function, x$predicted_function),
      within_pearson = safe_cor(within$observed_residual,
        within$predicted_residual),
      strict_pearson = safe_cor(strict$observed_residual,
        strict$predicted_residual),
      predicted_sd = sd(x$predicted_function),
      observed_sd = sd(x$observed_function),
      predicted_to_observed_sd_ratio = sd(x$predicted_function) /
        sd(x$observed_function))
  }))
}))
fwrite(per_drug, file.path(RESULT, "end_to_end_per_drug.tsv"), sep = "\t")

## Drug-block bootstrap for all activity-panel models and metrics.
bootstrap_parts <- list(); activity_predictions <- predictions[
  evaluation_included == TRUE & activity_supported == TRUE]
groups <- unique(activity_predictions[, .(prediction_model, baseline_source)])
for (g in seq_len(nrow(groups))) {
  current <- activity_predictions[prediction_model == groups$prediction_model[g] &
    fifelse(is.na(baseline_source), "", baseline_source) ==
      fifelse(is.na(groups$baseline_source[g]), "", groups$baseline_source[g])]
  for (metric_index in seq_along(c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation"))) {
    metric <- c("pooled", "within_perturbation",
      "strict_cell_line_by_perturbation")[metric_index]
    residual <- function_residual_table(current, metric)
    drugs <- sort(unique(residual$drug))
    original_counts <- table(factor(residual$drug, levels = drugs))
    set.seed(SEED + g * 10000L + metric_index * 1000L)
    draws <- rbindlist(lapply(seq_len(N_BOOTSTRAPS), function(draw) {
      sampled <- sample(drugs, length(drugs), replace = TRUE)
      multiplicity <- as.numeric(table(factor(sampled,
        levels = drugs))[match(residual$drug, drugs)])
      drug_weight <- multiplicity / as.numeric(original_counts[
        match(residual$drug, drugs)])
      data.table(bootstrap = draw,
        pearson = weighted_cor(residual$observed_residual,
          residual$predicted_residual, multiplicity),
        drug_balanced_pearson = weighted_cor(residual$observed_residual,
          residual$predicted_residual, drug_weight))
    }))
    draws[, `:=`(prediction_model = groups$prediction_model[g],
      baseline_source = groups$baseline_source[g], metric = metric)]
    bootstrap_parts[[length(bootstrap_parts) + 1L]] <- draws
  }
}
bootstrap <- rbindlist(bootstrap_parts, fill = TRUE)
fwrite(bootstrap, file.path(ARTIFACT, "drug_block_bootstrap.tsv.gz"),
  sep = "\t", compress = "gzip")
intervals <- bootstrap[, .(
  pearson_ci_lower = quantile(pearson, 0.025, na.rm = TRUE),
  pearson_ci_upper = quantile(pearson, 0.975, na.rm = TRUE),
  drug_balanced_pearson_ci_lower = quantile(drug_balanced_pearson,
    0.025, na.rm = TRUE),
  drug_balanced_pearson_ci_upper = quantile(drug_balanced_pearson,
    0.975, na.rm = TRUE)),
  by = .(prediction_model, baseline_source, metric)]
fwrite(intervals, file.path(RESULT, "drug_block_bootstrap_intervals.tsv"),
  sep = "\t")

write_artifact_manifest(ARTIFACT, "scripts/07_evaluate_end_to_end.R")
writeLines("complete", file.path(RESULT, "stage07_complete.flag"))
plot_stage07()
message("Stage 07 complete: end-to-end metrics and 1,000 drug bootstraps per model/metric.")
