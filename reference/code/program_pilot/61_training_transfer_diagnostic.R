#!/usr/bin/env Rscript
# No fitting: distinguish training capacity from transfer using saved predictions.
# Outer folds remain separate because their training observations overlap.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
root <- "/work/results/program_pilot/61_individual_responses_combination_transfer"
data_dir <- file.path(root, "data")
figure_dir <- file.path(root, "figures")
input <- file.path(data_dir, "ablation_predictions.tsv.gz")
models <- c("Original40_RNAonly_SVR", "ContextOnly_SVR", "Singles_logcounts_SVR")
model_labels <- c(Original40_RNAonly_SVR = "Original 40 programs: individual RNA + context",
                  ContextOnly_SVR = "Context only: doses, anchor and experiment",
                  Singles_logcounts_SVR = "Individual function + context")
stages <- c("training_fit", "inner_partner_OOF", "held_partner")
stage_labels <- c(training_fit = "Training fit", inner_partner_OOF = "Inner excluded partner\n(selection-used)",
                  held_partner = "Outer excluded partner")
write_table <- function(x, name) fwrite(x, file.path(data_dir, name), sep = "\t", na = "NA")
safe_cor <- function(a, b) {
  if (length(a) < 4L || sd(a) < 1e-10 || sd(b) < 1e-10) return(NA_real_)
  cor(a, b)
}
pred <- fread(input)[model %chin% models]
stopifnot(setequal(pred$split, stages), setequal(pred$time_h, c(3L, 24L)),
          !anyDuplicated(pred[, .(model, time_h, held_partner, split, observation_id)]),
          all(is.finite(pred$relative_AB)), all(is.finite(pred$predicted_relative_AB)))
stopifnot(all(pred[split == "held_partner", partner == held_partner]),
          all(pred[split != "held_partner", partner != held_partner]))
pred[, `:=`(observed_pct = 100 * relative_AB, predicted_pct = 100 * predicted_relative_AB)]
# One point per nominal dose condition within each model/time/outer-fold/stage.
# Averaging biological experiments is only a plotting view, not extra replication.
condition <- pred[, .(observed_pct = mean(observed_pct), predicted_pct = mean(predicted_pct),
                     biological_experiments = uniqueN(biological_replicate),
                     observed_min_pct = min(observed_pct), observed_max_pct = max(observed_pct),
                     observed_GR_mean = mean(GR_AB), predicted_GR_mean = mean(predicted_GR_AB)),
                 by = .(model, time_h, held_partner, split, condition_id, anchor, partner,
                        nominal_A_uM, nominal_B_uM)]
write_table(condition, "training_inner_held_condition_points.tsv.gz")
metrics <- condition[, .(dose_conditions = .N,
                         constituent_biological_observations = sum(biological_experiments),
                         Pearson_r = safe_cor(observed_pct, predicted_pct),
                         observed_median_pct = median(observed_pct), predicted_median_pct = median(predicted_pct),
                         median_prediction_minus_observed_pp = median(predicted_pct - observed_pct),
                         median_absolute_difference_pp = median(abs(predicted_pct - observed_pct)),
                         conditions_within10pp = sum(abs(predicted_pct - observed_pct) <= 10),
                         observed_at_most50pct = sum(observed_pct <= 50),
                         predicted_at_most50pct = sum(predicted_pct <= 50),
                         both_at_most50pct = sum(observed_pct <= 50 & predicted_pct <= 50),
                         false_at_most50pct = sum(observed_pct > 50 & predicted_pct <= 50)),
                     by = .(model, time_h, held_partner, split)]
write_table(metrics, "training_inner_held_metrics.tsv")
observation_metrics <- pred[, .(biological_observations = .N, unique_dose_conditions = uniqueN(condition_id),
                               Pearson_r = safe_cor(observed_pct, predicted_pct),
                               median_absolute_difference_pp = median(abs(predicted_pct - observed_pct)),
                               observed_median_pct = median(observed_pct), predicted_median_pct = median(predicted_pct)),
                           by = .(model, time_h, held_partner, split)]
write_table(observation_metrics, "training_inner_held_biological_observation_metrics.tsv")
anchor_metrics <- condition[, .(dose_conditions = .N, Pearson_r = safe_cor(observed_pct, predicted_pct),
                                 observed_median_pct = median(observed_pct), predicted_median_pct = median(predicted_pct),
                                 conditions_within10pp = sum(abs(predicted_pct - observed_pct) <= 10)),
                             by = .(model, time_h, held_partner, split, anchor)]
write_table(anchor_metrics, "training_inner_held_anchor_metrics.tsv")
# Summarize six model-fold diagnostics, never concatenate repeated training rows.
fold_summary <- metrics[, .(outer_folds = .N, median_fold_r = median(Pearson_r),
                            minimum_fold_r = min(Pearson_r), maximum_fold_r = max(Pearson_r),
                            median_fold_absolute_difference_pp = median(median_absolute_difference_pp),
                            median_fold_fraction_within10pp = median(conditions_within10pp / dose_conditions)),
                        by = .(model, time_h, split)]
write_table(fold_summary, "training_inner_held_fold_summary.tsv")

condition[, `:=`(held_partner = factor(held_partner, levels = sort(unique(as.character(held_partner)))),
                   split = factor(split, levels = stages))]
metrics[, `:=`(held_partner = factor(held_partner, levels = levels(condition$held_partner)),
               split = factor(split, levels = stages))]
metrics[, panel_label := sprintf("r = %.3f; n = %d\nwithin 10 pp: %d/%d", Pearson_r, dose_conditions,
                                  conditions_within10pp, dose_conditions)]
colors <- c(alpelisib = "#1B9E77", dactolisib = "#D95F02")
plots <- list()
for (m in models) {
  for (tt in c(3, 24)) {
    z <- condition[model == m & time_h == tt]
    ann <- metrics[model == m & time_h == tt]
    p <- ggplot(z, aes(observed_pct, predicted_pct)) +
      geom_point(aes(color = anchor), alpha = .55, size = .85) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#333333", linewidth = .45) +
      geom_label(data = ann, aes(x = -Inf, y = Inf, label = panel_label),
                 inherit.aes = FALSE, hjust = -.02, vjust = 1.03, size = 2.5,
                 fill = "white", alpha = .86, linewidth = 0, label.padding = grid::unit(.10, "lines")) +
      facet_grid(held_partner ~ split, scales = "free", labeller = labeller(split = stage_labels)) +
      scale_color_manual(values = colors, name = "Anchor") +
      scale_x_continuous(expand = expansion(mult = c(.03, .04))) +
      scale_y_continuous(expand = expansion(mult = c(.03, .16))) +
      labs(title = paste0(model_labels[m], " | ", tt, "-hour RNA"),
           subtitle = "Each row is a separately fitted outer exclusion; rows share training data. Points average biological repeats.",
           x = "Observed combined cell count (% of vehicle)", y = "Predicted combined cell count (% of vehicle)",
           caption = "n = dose conditions, not independent experiments. Inner OOF selected tuning and is not independent validation.\nAxes vary independently; gray lines are fitted trends. Each held partner is absent from all joint training outcomes.") +
      theme_bw(base_size = 9) +
      theme(aspect.ratio = 1, panel.grid.minor = element_blank(), legend.position = "top",
            strip.text = element_text(size = 8.2), strip.background = element_rect(fill = "#F1F1F1"),
            axis.text = element_text(size = 7.2), axis.title = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 8), plot.caption = element_text(size = 8, hjust = 0),
            panel.spacing = grid::unit(.35, "lines"), plot.margin = margin(7, 7, 7, 7))
    plots[[paste(m, tt, sep = "_")]] <- p
  }
}
pdf_path <- file.path(figure_dir, "training_inner_held_diagnostics.pdf")
pdf(pdf_path, width = 8.7, height = 15.3, onefile = TRUE, useDingbats = FALSE)
for (p in plots) print(p)
dev.off()
ggsave(file.path(figure_dir, "training_inner_held_diagnostics.png"), plots[[1L]], width = 8.7, height = 15.3, dpi = 150)
# Small companion plot shows all model/time/fold trajectories without pooling outcomes.
summary_metrics <- copy(metrics)
summary_metrics[, model := factor(model, levels = models)]
short_labels <- c(Original40_RNAonly_SVR = "Original 40 RNA + context", ContextOnly_SVR = "Context only",
                   Singles_logcounts_SVR = "Individual function + context")
summary_plot <- ggplot(summary_metrics, aes(split, Pearson_r, color = held_partner, group = held_partner)) +
  geom_line(linewidth = .45, alpha = .75) + geom_point(size = 1.5) +
  facet_grid(model ~ time_h, labeller = labeller(model = short_labels, time_h = function(x) paste0(x, "h RNA"))) +
  scale_x_discrete(labels = c(training_fit = "Training", inner_partner_OOF = "Inner OOF*", held_partner = "Held partner")) +
  labs(title = "Does the fitted relationship survive partner exclusion?", x = NULL,
       y = "Pearson r for combined cell counts", color = "Outer excluded partner",
       caption = "Each line tracks one outer fold; folds overlap in training observations. *Inner OOF was used to select tuning.\nCorrelations use dose-condition means across biological repeats. This is a diagnostic of saved models; no refitting.") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text.y = element_text(size = 8),
        plot.caption = element_text(hjust = 0, size = 8), axis.text.x = element_text(size = 8))
ggsave(file.path(figure_dir, "training_inner_held_fold_summary.pdf"), summary_plot, width = 8.5, height = 8.4)
ggsave(file.path(figure_dir, "training_inner_held_fold_summary.png"), summary_plot, width = 8.5, height = 8.4, dpi = 150)
writeLines(c(
  "Training-capacity versus transfer diagnostic, using existing predictions only. No fitting, tuning or model search.",
  sprintf("Original40 RNA-only fits its own training data: median fold correlation %.3f at3h and %.3f at24h.",
          fold_summary[model == "Original40_RNAonly_SVR" & time_h == 3 & split == "training_fit", median_fold_r],
          fold_summary[model == "Original40_RNAonly_SVR" & time_h == 24 & split == "training_fit", median_fold_r]),
  sprintf("For selected inner-partner OOF predictions those correlations fall to %.3f and %.3f; median fractions within10pp fall from %.1f%%/%.1f%% in training to %.1f%%/%.1f%%.",
          fold_summary[model == "Original40_RNAonly_SVR" & time_h == 3 & split == "inner_partner_OOF", median_fold_r],
          fold_summary[model == "Original40_RNAonly_SVR" & time_h == 24 & split == "inner_partner_OOF", median_fold_r],
          100 * fold_summary[model == "Original40_RNAonly_SVR" & time_h == 3 & split == "training_fit", median_fold_fraction_within10pp],
          100 * fold_summary[model == "Original40_RNAonly_SVR" & time_h == 24 & split == "training_fit", median_fold_fraction_within10pp],
          100 * fold_summary[model == "Original40_RNAonly_SVR" & time_h == 3 & split == "inner_partner_OOF", median_fold_fraction_within10pp],
          100 * fold_summary[model == "Original40_RNAonly_SVR" & time_h == 24 & split == "inner_partner_OOF", median_fold_fraction_within10pp]),
  "The saved model has fitting capacity. Its relationships lose accuracy when a partner is excluded; near-perfect training fits do not establish transfer.",
  "The simpler context-only model fits training data less tightly but does not incur the same transfer penalty. Training fit alone therefore cannot attribute useful information to RNA.",
  "Primary plot PDF pages1-2: original40 RNA-only3h/24h; pages3-4: context-only; pages5-6: measured-single-function comparator.",
  "Every outer fold is kept separate. Training rows recur across folds and are not counted as independent measurements.",
  "Inner partner OOF predictions used the selected tuning candidate and thus are selection-used diagnostics, not independent validation.",
  "Points average2-3 biological repeats within one nominal dose condition. The separate observation-level table retains individual biological experiments.",
  "A 50%-of-vehicle cell count threshold describes growth suppression, not necessarily net cell loss. A10-percentage-point tolerance is descriptive, not significance.",
  "Compare training-fit correlation and amplitude with inner and outer exclusions; direct combined counts avoid the shared-reference issue of posthoc excess-effect correlations.",
  "Six-fold summaries report distributions of fold diagnostics. They are not a pooled training performance estimate or confidence interval.",
  capture.output(print(fold_summary))
), file.path(root, "TRAINING_TRANSFER_DIAGNOSTIC.txt"))
write_table(data.table(path = c(input, "/work/src/program_pilot/61_training_transfer_diagnostic.R"),
                       sha256 = vapply(c(input, "/work/src/program_pilot/61_training_transfer_diagnostic.R"),
                                       digest, character(1), file = TRUE, algo = "sha256")),
             "training_inner_held_source_manifest.tsv")
print(fold_summary)
