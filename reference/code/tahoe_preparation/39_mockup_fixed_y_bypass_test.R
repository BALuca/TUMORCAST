#!/usr/bin/env Rscript
## Short Aim 2 mockup using the fixed Y-shaped BRAF/ERBB -> P3 -> function model.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

root <- "/work/tahoe_prism_v2"
out <- file.path(root, "results/39_fixed_y_bypass_test")
unlink(out, recursive = TRUE, force = TRUE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, name),
  sep = "\t", na = "NA")

stage38 <- file.path(root,
  "results/38_mockup_baseline_conditioned_braf_transfer")
stage37 <- file.path(root, "results/37_source_aware_bypass_validation")
prediction <- fread(file.path(stage38,
  "tahoe_to_gse68453_function_predictions.tsv"))[
    normalization == "program_scaled_P3"]
factorial <- fread(file.path(stage37, "gse68453_program_factorial.tsv"))
functional_stats <- fread(file.path(stage37,
  "gse68453_functional_prediction_stats.tsv"))

## The topology is fixed.  The two intervention arms meet at P3, so their
## conditional molecular effect is the factorial difference-in-differences.
p3 <- factorial[program == "P3", .(cell_line, control_P3 = control,
  lapatinib_P3 = lapatinib, plx_P3 = plx, combination_P3 = plx_lap,
  P3_interaction = interaction,
  extra_P3_reduction = -interaction)]
d <- merge(prediction[, .(cell_line, baseline_predicted_bypass =
  predicted_residual_viability, bliss)], p3, by = "cell_line")
write_table(d, "fixed_y_line_data.tsv")

edge_test <- function(x, y, edge, role) {
  ok <- is.finite(x) & is.finite(y)
  test <- cor.test(x[ok], y[ok])
  data.table(edge, role, n = sum(ok), pearson = unname(test$estimate),
    p_value = test$p.value, ci_low = test$conf.int[1L],
    ci_high = test$conf.int[2L])
}
edges <- rbindlist(list(
  edge_test(d$baseline_predicted_bypass, d$plx_P3,
    "Baseline prediction -> residual P3 during PLX4720", "Aim 1 to Aim 2"),
  edge_test(d$baseline_predicted_bypass, d$extra_P3_reduction,
    "Baseline prediction -> lapatinib-specific P3 reduction",
    "Conditional molecular effect"),
  edge_test(d$extra_P3_reduction, d$bliss,
    "Lapatinib-specific P3 reduction -> Bliss benefit",
    "Molecular interaction to functional interaction"),
  edge_test(d$baseline_predicted_bypass, d$bliss,
    "Baseline prediction -> Bliss benefit", "End-to-end prediction")))
write_table(edges, "fixed_y_edge_tests.tsv")

program_comparison <- copy(functional_stats)
setnames(program_comparison, "pearson", "extra_program_reduction_to_bliss")
program_comparison[, selected_by_fixed_graph := program == "P3"]
write_table(program_comparison, "downstream_program_comparison.tsv")

## Predict explicit EGF bypass using the fixed P2 -> P3 edge.  Learn the edge
## from variation among GSE68453 lines during PLX4720, then transfer it without
## refitting to the independent GSE127988 EGF time course during vemurafenib.
p23_training <- dcast(factorial[program %chin% c("P2", "P3"),
  .(cell_line, program, plx)], cell_line ~ program, value.var = "plx")
p23_fit <- lm(P3 ~ 0 + P2, data = p23_training)
p23_training[, predicted_P3 := predict(p23_fit, newdata = p23_training)]
p23_test <- cor.test(p23_training$P2, p23_training$P3)
p23_training_stats <- data.table(n_lines = nrow(p23_training),
  slope_through_origin = unname(coef(p23_fit)["P2"]),
  pearson = unname(p23_test$estimate), p_value = p23_test$p.value,
  rmse = sqrt(mean((p23_training$predicted_P3 - p23_training$P3)^2)))

egf <- fread(file.path(stage37, "gse127988_program_scores.tsv"))
egf_wide <- dcast(egf[program %chin% c("P2", "P3"),
  .(arm, egf_h, biological_set, program, score)],
  arm + egf_h + biological_set ~ program, value.var = "score",
  fun.aggregate = mean)
egf_wide[, predicted_P3 := predict(p23_fit, newdata = egf_wide)]
egf_stats <- egf_wide[, {
  test <- cor.test(predicted_P3, P3)
  .(n_conditions = .N, pearson = unname(test$estimate),
    p_value = test$p.value,
    rmse = sqrt(mean((predicted_P3 - P3)^2)))
}, by = arm]
egf_day_stats <- egf_wide[, .(n_time_points = .N,
  pearson = cor(predicted_P3, P3),
  rmse = sqrt(mean((predicted_P3 - P3)^2))),
  by = .(arm, biological_set)]
write_table(p23_training, "gse68453_learned_P2_to_P3.tsv")
write_table(p23_training_stats, "gse68453_learned_P2_to_P3_stats.tsv")
write_table(egf_wide, "gse127988_predicted_egf_P3_response.tsv")
write_table(egf_stats, "gse127988_predicted_egf_P3_stats.tsv")
write_table(egf_day_stats, "gse127988_predicted_egf_P3_by_day.tsv")

theme_panel <- theme_classic(base_size = 9.5) + theme(aspect.ratio = 1,
  plot.title = element_text(size = 10.3, face = "bold"))
edge_value <- function(pattern) edges[grepl(pattern, edge), pearson]
pa <- ggplot(d, aes(baseline_predicted_bypass, plx_P3)) +
  geom_point(size = 2.2, color = "#24745E") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) + theme_panel +
  labs(x = "Baseline-predicted bypass", y = "P3 response during PLX4720",
    title = sprintf("Bypass state, n=12; r=%.2f",
      edge_value("residual P3")))
pb <- ggplot(d[is.finite(extra_P3_reduction)],
  aes(baseline_predicted_bypass, extra_P3_reduction)) +
  geom_point(size = 2.2, color = "#276C9B") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) + theme_panel +
  labs(x = "Baseline-predicted bypass",
    y = "Extra P3 reduction with lapatinib",
    title = sprintf("Conditional effect, n=5; r=%.2f",
      edge_value("lapatinib-specific P3 reduction$")))
pc <- ggplot(d, aes(baseline_predicted_bypass, bliss)) +
  geom_point(size = 2.2, color = "#5B4775") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) + theme_panel +
  labs(x = "Baseline-predicted bypass",
    y = "PLX4720 + lapatinib Bliss benefit",
    title = sprintf("Bliss benefit, n=12; r=%.2f",
      edge_value("Baseline prediction -> Bliss")))
p <- ggpubr::ggarrange(pa, pb, pc, nrow = 1, align = "hv")
p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(
  "Fixed Y model: pretreatment state predicts BRAF-ERBB bypass",
  face = "bold", size = 13))
ggsave(file.path(out, "figures/01_fixed_y_bypass_test.pdf"), p,
  width = 8.4, height = 3.45, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/01_fixed_y_bypass_test.png"), p,
  width = 8.4, height = 3.45, dpi = 240, bg = "white")

program_colors <- c("P2" = "#377EB8", "P3" = "#4DAF4A")
pa <- ggplot(p23_training, aes(P2, P3)) +
  geom_hline(yintercept = 0, color = "#D8DDDF", linewidth = .4) +
  geom_vline(xintercept = 0, color = "#D8DDDF", linewidth = .4) +
  geom_point(size = 2.4, color = program_colors["P3"]) +
  geom_abline(slope = coef(p23_fit)["P2"], intercept = 0,
    color = "#263840", linewidth = .8) + theme_panel +
  labs(x = "Measured P2 response during PLX4720",
    y = "Measured P3 response during PLX4720",
    title = sprintf("Learn P2 -> P3, n=12; r=%.2f",
      p23_training_stats$pearson))
egf_vem <- egf_wide[arm == "Vem"]
vem_stat <- egf_stats[arm == "Vem"]
pb <- ggplot(egf_vem, aes(predicted_P3, P3,
  color = factor(biological_set))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
    color = "#9AA4A9", linewidth = .6) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) +
  scale_color_manual(values = c("3" = "#D88736", "4" = "#5B4775"),
    name = "Biological day") + theme_panel +
  theme(legend.position = "bottom") +
  labs(x = "P3 predicted from EGF-induced P2",
    y = "Measured P3 response after EGF",
    title = sprintf("Predict EGF bypass, n=14; r=%.2f", vem_stat$pearson))
egf_means <- egf_wide[, .(P3 = mean(P3), predicted_P3 = mean(predicted_P3)),
  by = .(arm, egf_h)]
pc <- ggplot(egf_means, aes(egf_h, P3, color = arm)) +
  geom_hline(yintercept = 0, color = "#D8DDDF", linewidth = .4) +
  geom_line(linewidth = .8) + geom_point(size = 2.1) +
  geom_line(data = egf_means[arm == "Vem"], aes(y = predicted_P3),
    inherit.aes = TRUE, color = "#263840", linewidth = .7,
    linetype = "dashed") +
  scale_color_manual(values = c("Vem" = "#276C9B",
    "Vem + Cobi" = "#A33870"), name = "Background") + theme_panel +
  theme(legend.position = "bottom") +
  labs(x = "Hours after adding EGF", y = "Measured P3 response",
    title = "Cobimetinib blocks the response")
p <- ggpubr::ggarrange(pa, pb, pc, nrow = 1, align = "hv")
p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(
  "A program relationship learned elsewhere predicts EGF-driven bypass",
  face = "bold", size = 13))
ggsave(file.path(out, "figures/02_predict_egf_bypass.pdf"), p,
  width = 8.6, height = 3.65, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/02_predict_egf_bypass.png"), p,
  width = 8.6, height = 3.65, dpi = 240, bg = "white")

writeLines(c("# Fixed Y-shaped bypass test", "",
  "The topology was specified before this analysis: PLX4720 and lapatinib are the two intervention arms, P3 is their shared downstream response, and Bliss is the functional interaction. The molecular interaction is the factorial difference-in-differences: -(combination - PLX4720 - lapatinib + control).", "",
  sprintf(paste0("The Tahoe-derived baseline prediction identifies residual P3 ",
    "during PLX4720 across 12 lines (r = %.2f, p = %.3f) and predicts Bliss ",
    "combination benefit (r = %.2f, p = %.3f)."),
    edges[role == "Aim 1 to Aim 2", pearson],
    edges[role == "Aim 1 to Aim 2", p_value],
    edges[role == "End-to-end prediction", pearson],
    edges[role == "End-to-end prediction", p_value]),
  sprintf(paste0("In the five lines with complete factorial RNA, baseline ",
    "prediction is directionally associated with the lapatinib-specific P3 ",
    "reduction (r = %.2f). That molecular interaction is associated with Bliss ",
    "benefit at r = %.2f. P3 gives the strongest molecular-to-functional ",
    "correlation among the four frozen programs."),
    edges[role == "Conditional molecular effect", pearson],
    edges[role == "Molecular interaction to functional interaction", pearson]),
  "", "## Explicit EGF prediction", "",
  sprintf(paste0("A P2-to-P3 response coefficient learned from 12 GSE68453 ",
    "melanoma lines during PLX4720 was frozen and transferred to GSE127988. ",
    "When EGF is added during vemurafenib, its measured P2 response predicts ",
    "the downstream P3 response at r = %.2f across 14 day-by-time conditions. ",
    "The correlations are %.2f and %.2f on the two biological days. ",
    "Cobimetinib suppresses both the EGF-induced P2 input and P3 output."),
    egf_stats[arm == "Vem", pearson],
    egf_day_stats[arm == "Vem" & biological_set == 3, pearson],
    egf_day_stats[arm == "Vem" & biological_set == 4, pearson]),
  "", "Interpretation: current data support the fixed pathway direction and the end-to-end baseline prediction. The direct prediction of the molecular combination interaction remains a five-line pilot."),
  file.path(out, "README.md"))

message("Stage 39 fixed-Y mockup complete: ", out)
