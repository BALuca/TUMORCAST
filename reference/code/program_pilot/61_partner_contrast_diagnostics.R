# Matched partner contrasts from frozen excluded-partner predictions; no model fitting.
# Run: Rscript /work/src/program_pilot/61_partner_contrast_diagnostics.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

result_root <- "/work/results/program_pilot/61_individual_responses_combination_transfer"
input_path <- file.path(result_root, "data/ablation_predictions.tsv.gz")
script_path <- "/work/src/program_pilot/61_partner_contrast_diagnostics.R"
models <- c("ContextOnly_SVR", "Original40_RNAonly_SVR", "Biology_RNAonly_SVR",
            "GenePC12_RNAonly_SVR", "Singles_logcounts_SVR")
model_labels <- c("Context only", "Original 40: RNA only", "Biology: RNA only",
                  "Gene PCs (12): RNA only", "Measured singles")
large_difference_pp <- 10
nominal_rounding_tolerance_pct <- 3.5

# All five models and both predeclared RNA times are retained. Differences are
# relative cell counts at 72 h, in percentage points of vehicle; not synergy.
predictions <- fread(input_path)
required <- c("observation_id", "time_h", "split", "held_partner", "model", "cell",
              "partner", "anchor", "nominal_A_uM", "nominal_B_uM", "dose_key_A",
              "dose_key_B", "biological_replicate", "dose_A_uM", "dose_B_uM",
              "ratio_AB", "predicted_relative_AB", "source_actual_doses_A",
              "source_actual_doses_B", "source_sample_ids_AB")
stopifnot(all(required %in% names(predictions)))
held <- predictions[split == "held_partner" & model %in% models, ..required]
stopifnot(setequal(unique(held$model), models), setequal(unique(held$time_h), c(3, 24)))
stopifnot(held[, anyDuplicated(paste(model, time_h, observation_id))] == 0L)
stopifnot(all(held$held_partner == held$partner))
stopifnot(all(is.finite(held$ratio_AB)), all(is.finite(held$predicted_relative_AB)))
stopifnot(all(held$nominal_A_uM > 0), all(held$nominal_B_uM > 0))
stopifnot(held[, .N, by = .(time_h, observation_id)][, all(N == length(models))])
stopifnot(held[, uniqueN(ratio_AB), by = .(time_h, observation_id)][, all(V1 == 1L)])

# A contrast requires the same anchor, both nominal concentrations and the same
# biological experiment. Alphabetical partner 1 minus partner 2 fixes orientation
# without using outcomes. Technical wells are already averaged in the input.
match_keys <- c("model", "time_h", "cell", "anchor", "nominal_A_uM", "nominal_B_uM",
                "dose_key_A", "dose_key_B", "biological_replicate")
stopifnot(held[, anyDuplicated(.SD), .SDcols = c(match_keys, "partner")] == 0L)
contrasts <- merge(held, held, by = match_keys, allow.cartesian = TRUE,
                   suffixes = c("_1", "_2"))
contrasts <- contrasts[partner_1 < partner_2]
stopifnot(nrow(contrasts) > 0L)
contrasts[, `:=`(
  partner_pair = paste(partner_1, "minus", partner_2),
  measured_difference_pp = 100 * (ratio_AB_1 - ratio_AB_2),
  predicted_difference_pp = 100 * (predicted_relative_AB_1 - predicted_relative_AB_2),
  actual_anchor_mismatch_pct = 100 * abs(dose_A_uM_1 - dose_A_uM_2) / nominal_A_uM,
  actual_partner_mismatch_pct = 100 * abs(dose_B_uM_1 - dose_B_uM_2) / nominal_B_uM,
  anchor_1_rounding_pct = 100 * abs(dose_A_uM_1 / nominal_A_uM - 1),
  anchor_2_rounding_pct = 100 * abs(dose_A_uM_2 / nominal_A_uM - 1),
  partner_1_rounding_pct = 100 * abs(dose_B_uM_1 / nominal_B_uM - 1),
  partner_2_rounding_pct = 100 * abs(dose_B_uM_2 / nominal_B_uM - 1)
)]
rounding_columns <- c("anchor_1_rounding_pct", "anchor_2_rounding_pct",
                      "partner_1_rounding_pct", "partner_2_rounding_pct")
stopifnot(max(as.matrix(contrasts[, ..rounding_columns])) <= nominal_rounding_tolerance_pct)
fwrite(contrasts, file.path(result_root, "data/partner_contrast_replicates.tsv.gz"))

# Average only contrasts formed inside each biological replicate. An unmatched
# replicate is never compared with another experiment to manufacture a contrast.
contrast_keys <- c("model", "time_h", "cell", "anchor", "nominal_A_uM", "nominal_B_uM",
                   "dose_key_A", "dose_key_B", "partner_1", "partner_2", "partner_pair")
averaged <- contrasts[, .(
  n_biological_replicates = uniqueN(biological_replicate),
  biological_replicates = paste(sort(unique(biological_replicate)), collapse = ";"),
  measured_difference_pp = mean(measured_difference_pp),
  measured_min_pp = min(measured_difference_pp),
  measured_max_pp = max(measured_difference_pp),
  measured_sd_pp = sd(measured_difference_pp),
  predicted_difference_pp = mean(predicted_difference_pp),
  predicted_min_pp = min(predicted_difference_pp),
  predicted_max_pp = max(predicted_difference_pp),
  observed_replicate_direction_consistent = all(sign(measured_difference_pp) == sign(mean(measured_difference_pp))),
  max_actual_anchor_mismatch_pct = max(actual_anchor_mismatch_pct),
  max_actual_partner_mismatch_pct = max(actual_partner_mismatch_pct)
), by = contrast_keys]
averaged[, `:=`(
  measured_large = abs(measured_difference_pp) >= large_difference_pp,
  predicted_large = abs(predicted_difference_pp) >= large_difference_pp,
  direction_correct = sign(measured_difference_pp) == sign(predicted_difference_pp)
)]
fwrite(averaged, file.path(result_root, "data/partner_contrast_means.tsv.gz"))

# Coverage lists all possible shared dose contexts, including any context having
# only one partner. Same micromolar concentration does not mean equal potency.
observations <- held[model == models[1]]
context_keys <- setdiff(match_keys, "model")
coverage <- observations[, .(
  n_partners = uniqueN(partner),
  partners = paste(sort(unique(partner)), collapse = ";"),
  n_available_observations = .N,
  n_pairwise_contrasts = choose(uniqueN(partner), 2)
), by = context_keys]
coverage_summary <- coverage[, .(
  n_observations = sum(n_available_observations),
  n_shared_contexts = sum(n_partners >= 2),
  n_single_partner_contexts = sum(n_partners == 1),
  n_observations_with_comparator = sum(n_available_observations[n_partners >= 2]),
  n_replicate_contrasts = sum(n_pairwise_contrasts)
), by = .(time_h)]
fwrite(coverage, file.path(result_root, "data/partner_contrast_context_coverage.tsv"))
fwrite(coverage_summary, file.path(result_root, "data/partner_contrast_coverage_summary.tsv"))
pair_coverage <- averaged[model == models[1], .(
  n_dose_contexts = .N,
  n_replicate_contrasts = sum(n_biological_replicates),
  n_two_or_more_replicates = sum(n_biological_replicates >= 2),
  partner_doses_uM = paste(sort(unique(nominal_B_uM)), collapse = ";"),
  anchor_doses_uM = paste(sort(unique(nominal_A_uM)), collapse = ";"),
  n_measured_large = sum(measured_large),
  n_large_replicated = sum(measured_large & n_biological_replicates >= 2),
  n_large_replicate_consistent = sum(measured_large & n_biological_replicates >= 2 &
                                     observed_replicate_direction_consistent)
), by = .(time_h, anchor, partner_pair)]
fwrite(pair_coverage, file.path(result_root, "data/partner_contrast_pair_coverage.tsv"))

# Descriptive association and assay-scale amplitudes; no independence-based
# inference. Pairwise contrasts share observations and outer training fits differ.
summarize_contrasts <- function(values) {
  measured <- values$measured_difference_pp
  predicted <- values$predicted_difference_pp
  large <- values$measured_large
  list(
    n_dose_contexts = nrow(values),
    pearson_r = if (sd(measured) > 0 && sd(predicted) > 0) cor(measured, predicted) else NA_real_,
    predicted_per_measured_slope = if (sd(measured) > 0) cov(measured, predicted) / var(measured) else NA_real_,
    median_abs_measured_pp = median(abs(measured)),
    median_abs_predicted_pp = median(abs(predicted)),
    n_measured_large = sum(large),
    n_large_direction_correct = sum(values$direction_correct & large),
    n_large_also_predicted_large_correct = sum(values$direction_correct & large & values$predicted_large),
    n_large_wrong_direction = sum(!values$direction_correct & large)
  )
}
metrics <- averaged[, summarize_contrasts(.SD), by = .(time_h, model)]
anchor_metrics <- averaged[, summarize_contrasts(.SD), by = .(time_h, model, anchor)]
pair_metrics <- averaged[, summarize_contrasts(.SD), by = .(time_h, model, anchor, partner_pair)]
fwrite(metrics, file.path(result_root, "data/partner_contrast_metrics.tsv"))
fwrite(anchor_metrics, file.path(result_root, "data/partner_contrast_anchor_metrics.tsv"))
fwrite(pair_metrics, file.path(result_root, "data/partner_contrast_pair_metrics.tsv"))

# Readable five-model overview: biological-replicate means, observed replicate
# range as horizontal strokes, independently scaled x/y axes and Pearson trends.
averaged[, model_label := factor(model, levels = models, labels = model_labels)]
metrics[, model_label := factor(model, levels = models, labels = model_labels)]
metrics[, label := sprintf("r = %.2f; n = %d\nMedian |predicted| = %.1f pp", pearson_r,
                            n_dose_contexts, median_abs_predicted_pp)]
for (rna_time in c(3, 24)) {
  plot_data <- averaged[time_h == rna_time]
  labels <- metrics[time_h == rna_time]
  title_status <- if (rna_time == 3) "Primary: 3 h RNA" else "Secondary: 24 h RNA"
  figure <- ggplot(plot_data, aes(measured_difference_pp, predicted_difference_pp)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_segment(aes(x = measured_min_pp, xend = measured_max_pp,
                     yend = predicted_difference_pp, colour = anchor),
                 alpha = 0.12, linewidth = 0.25) +
    geom_point(aes(colour = anchor), size = 1.3, alpha = 0.60) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.55) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.04, vjust = 1.10, size = 3.0, lineheight = 1.0) +
    facet_wrap(~model_label, ncol = 3, scales = "free_y") +
    scale_colour_manual(values = c(alpelisib = "#087E8B", dactolisib = "#C16622")) +
    labs(title = paste(title_status, "— changing only the partner"),
         subtitle = "Same anchor, nominal anchor dose, nominal partner dose and biological replicate",
         x = "Measured combination cell-count difference (percentage points of vehicle)",
         y = "Predicted combination cell-count difference (pp)", colour = "Anchor",
         caption = paste0("Alphabetical partner 1 minus partner 2; positive = more remaining cells with partner 1. ",
                          "Points: matched-replicate means; strokes: observed replicate range.\n",
                          "Contrasts share observations; two partners use separate excluded-partner fits. ",
                          "No inferential tests. Same concentration does not imply equal potency.")) +
    theme_bw(base_size = 10) +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey93"),
          strip.background = element_rect(fill = "grey95"), legend.position = "bottom",
          axis.title = element_text(size = 10), plot.caption = element_text(size = 8, hjust = 0),
          plot.title = element_text(face = "bold", size = 13), plot.margin = margin(6, 8, 6, 6))
  figure_stem <- file.path(result_root, "figures", paste0("partner_contrasts_", rna_time, "h"))
  ggsave(paste0(figure_stem, ".pdf"), figure, width = 10.2, height = 7.5, device = cairo_pdf)
  ggsave(paste0(figure_stem, ".png"), figure, width = 10.2, height = 7.5, dpi = 160)
}

# Explicit source and limitation record accompanies the plotted source tables.
primary_original <- metrics[time_h == 3 & model == "Original40_RNAonly_SVR"]
secondary_original <- metrics[time_h == 24 & model == "Original40_RNAonly_SVR"]
primary_singles <- metrics[time_h == 3 & model == "Singles_logcounts_SVR"]
large_measured <- averaged[time_h == 3 & model == models[1] & measured_large]
notes <- c(
  "MATCHED PARTNER CONTRAST DIAGNOSTIC — WORKFLOW 61",
  "Purpose: determine whether aggregate combination-cell-count prediction contains partner-specific information.",
  "RESULT AND INTERPRETATION",
  sprintf(paste0("Original40 primary RNA-only predicts matched partner differences with r=%.3f; ",
                 "%d/%d observed differences >=10 pp have the correct predicted direction, ",
                 "and %d also have a correctly directed prediction >=10 pp."),
          primary_original$pearson_r, primary_original$n_large_direction_correct,
          primary_original$n_measured_large, primary_original$n_large_also_predicted_large_correct),
  sprintf("Measured-single input predicts these differences more consistently (r=%.3f; direction %d/%d).",
          primary_singles$pearson_r, primary_singles$n_large_direction_correct, primary_singles$n_measured_large),
  sprintf("At the secondary24 h RNA time, Original40 r=%.3f; direction %d/%d.",
          secondary_original$pearson_r, secondary_original$n_large_direction_correct,
          secondary_original$n_measured_large),
  sprintf("All dose-context means have2-3 biological replicates; %d/%d large measured contrasts agree in sign in every replicate.",
          sum(large_measured$observed_replicate_direction_consistent), nrow(large_measured)),
  "There is reproducible measured partner variation and partial recovery by original40 RNA; aggregate predictability overstates that recovery.",
  "This tests total count differences between partners, not prediction of combination-specific departure from single-drug effects.",
  "It does not establish synergy prediction, cross-cell transport, or a general compositional mechanism.",
  "DESIGN AND LIMITATIONS",
  "No new models, representation selection, tuning, outcome censoring or inferential tests were performed.",
  "Input: frozen held_partner predictions from ablation_predictions.tsv.gz; all five models retained.",
  "Primary RNA time 3 h; secondary 24 h. Functional endpoint: imaged cell count at 72 h / vehicle count.",
  "Contrast = 100 * (relative cell count for alphabetical partner 1 - partner 2). This is not a synergy contrast.",
  "Match: same cell, anchor, nominal anchor/partner doses and biological replicate; technical wells already averaged.",
  "First form each within-replicate contrast, then average available replicates. Replicate range is descriptive.",
  "Every feasible partner pair is included; grids vary by partner/replicate, so contrasts are not a balanced panel.",
  "Pairwise contrasts share observations and are dependent; different excluded-partner predictions use different fits.",
  "Context-only can therefore give nonzero partner contrasts solely because the two outer training fits differ.",
  "Its contrasts do not demonstrate a partner-sensing mechanism; they quantify this separate-fit background.",
  "Same nominal micromolar dose does not mean equal target engagement or potency across different compounds.",
  sprintf("Maximum actual anchor-dose mismatch within a nominal match: %.5f%%.",
          max(contrasts$actual_anchor_mismatch_pct)),
  sprintf("Maximum actual partner-dose mismatch within a nominal match: %.5f%%.",
          max(contrasts$actual_partner_mismatch_pct)),
  sprintf("Largest actual-versus-nominal rounding: %.5f%% (allowed %.1f%%).",
          max(as.matrix(contrasts[, ..rounding_columns])), nominal_rounding_tolerance_pct),
  "10 percentage points is a predeclared descriptive magnitude threshold, not a significance or reliability cutoff.",
  "Metrics are over equally weighted matched dose-context means, not independent biological units.",
  "Positive Pearson association may coexist with heavily attenuated amplitudes; consult pp and direction counts.",
  paste("Input MD5:", unname(tools::md5sum(input_path))),
  paste("Script MD5:", unname(tools::md5sum(script_path))),
  paste("Run UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "", "COVERAGE", capture.output(print(coverage_summary)), "", "METRICS", capture.output(print(metrics)),
  "", "BIOLOGICAL REPLICATE SUPPORT (primary; all pairs dependent)",
  capture.output(print(averaged[model == models[1] & time_h == 3,
    .(n_contexts = .N, n_measured_large = sum(measured_large),
      n_large_sign_consistent = sum(measured_large & observed_replicate_direction_consistent)),
    by = n_biological_replicates]))
)
writeLines(notes, file.path(result_root, "PARTNER_CONTRAST_DIAGNOSTIC.txt"))
writeLines(capture.output(sessionInfo()), file.path(result_root, "data/partner_contrast_session_info.txt"))
print(coverage_summary)
print(metrics)
