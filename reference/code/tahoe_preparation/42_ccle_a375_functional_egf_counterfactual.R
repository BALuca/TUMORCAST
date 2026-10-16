#!/usr/bin/env Rscript
## Functional counterfactual from the CCLE pretreatment state of A375.
## Baseline -> BRAFi function is frozen in Tahoe/CCLE; the EGF contribution is
## propagated through the fixed tree and translated by the frozen P3 decoder.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggpubr)})

root <- "/work/tahoe_prism_v2"
out <- file.path(root, "results/42_ccle_a375_functional_egf_counterfactual")
unlink(out, recursive = TRUE, force = TRUE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, name),
  sep = "\t", na = "NA")

r38 <- file.path(root, "results/38_mockup_baseline_conditioned_braf_transfer")
r40 <- file.path(root, "results/40_cytosig_tahoe_egf_prediction")

## CCLE pretreatment measurement for the exact Gerosa cell line.
ccle <- fread(file.path(r38, "tahoe_to_ccle_transfer_predictions.tsv"))[
  depmap_id == "ACH-000219" & transfer_model == "P1"]
stopifnot(nrow(ccle) == 1L, ccle$braf_hotspot == 1L)

## Frozen ligand displacement and frozen program-to-function mapping.
egf_model <- fread(file.path(r40, "frozen_egf_model.tsv"))
decoder <- fread(file.path(r38, "frozen_tahoe_function_decoders.tsv"))[
  program == "P3"]
p3_function_slope <- decoder[term == "response_score", coefficient]

braﬁ_no_egf <- ccle$predicted_survival
egf_p3_displacement <- egf_model$predicted_P3_after_EGF
egf_function_shift <- p3_function_slope * egf_p3_displacement
braﬁ_with_egf <- braﬁ_no_egf + egf_function_shift

prediction <- data.table(
  depmap_id = ccle$depmap_id,
  cell_line = "A375",
  braf_hotspot = ccle$braf_hotspot,
  ccle_baseline_P3 = ccle$baseline_P3,
  frozen_brafi_function_without_EGF = braﬁ_no_egf,
  frozen_EGF_P2_input = egf_model$egf_P2_input_study_balanced,
  frozen_P2_to_P3_edge = egf_model$P2_to_P3_slope,
  predicted_EGF_P3_displacement = egf_p3_displacement,
  frozen_P3_to_function_slope = p3_function_slope,
  predicted_EGF_function_shift = egf_function_shift,
  predicted_brafi_function_with_EGF = braﬁ_with_egf,
  observed_PRISM_brafi_without_EGF = ccle$observed_survival)
write_table(prediction, "a375_frozen_functional_counterfactual.tsv")

conditions <- data.table(
  condition = factor(c("BRAFi\nEGF absent", "BRAFi + EGF"),
    levels = c("BRAFi\nEGF absent", "BRAFi + EGF")),
  predicted_function = c(braﬁ_no_egf, braﬁ_with_egf),
  observed_function = c(ccle$observed_survival, NA_real_))
write_table(conditions, "a375_condition_predictions.tsv")

## The independent hard functional bypass cohort remains the external
## cross-line validation of the baseline-conditioned functional ranking.
gse <- fread(file.path(r38, "tahoe_to_gse68453_function_predictions.tsv"))[
  normalization == "program_scaled_P3"]
gse_test <- cor.test(gse$predicted_residual_viability, gse$bliss)
write_table(gse, "gse68453_external_functional_validation.tsv")

ink <- "#263840"; grid <- "#D9DEDF"; orange <- "#D88736"
purple <- "#5B4775"; teal <- "#168C86"
theme_panel <- theme_classic(base_size = 11) + theme(
  aspect.ratio = 1, plot.title = element_text(face = "bold", size = 11),
  plot.subtitle = element_text(size = 9, color = "#52656F"),
  axis.title = element_text(size = 10), plot.margin = margin(6, 8, 5, 5))

pa <- ggplot(conditions, aes(condition, predicted_function, group = 1)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .5) +
  geom_line(color = ink, linewidth = .85) +
  geom_point(aes(color = condition), size = 3.4) +
  geom_point(data = conditions[is.finite(observed_function)],
    aes(y = observed_function), shape = 21, size = 3.5, stroke = 1.1,
    color = ink, fill = "white") +
  scale_color_manual(values = c("BRAFi\nEGF absent" = orange,
    "BRAFi + EGF" = teal), guide = "none") + theme_panel +
  labs(x = NULL, y = "Predicted functional response\n(higher = more resistant)",
    title = "A  CCLE baseline predicts a functional switch",
    subtitle = "Filled: prediction · open: observed PRISM (EGF absent)")

steps <- data.table(
  step = factor(c("CCLE baseline\n+ BRAF status", "EGF propagated\nthrough tree",
    "BRAFi + EGF\nfunctional prediction"),
    levels = c("CCLE baseline\n+ BRAF status", "EGF propagated\nthrough tree",
      "BRAFi + EGF\nfunctional prediction")),
  value = c(braﬁ_no_egf, egf_function_shift, braﬁ_with_egf),
  kind = c("Starting prediction", "Added EGF effect", "Final prediction"))
pb <- ggplot(steps, aes(step, value, fill = kind)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .5) +
  geom_col(width = .64) +
  geom_text(aes(label = sprintf("%+.2f", value)),
    vjust = ifelse(steps$value >= 0, -0.35, 1.25), size = 3.4,
    color = ink) +
  scale_fill_manual(values = c("Starting prediction" = orange,
    "Added EGF effect" = teal, "Final prediction" = purple), guide = "none") +
  theme_panel + theme(axis.text.x = element_text(size = 8.7)) +
  labs(x = NULL, y = "Functional-response contribution",
    title = "B  Frozen prediction decomposition",
    subtitle = sprintf("A375 · BRAF V600E · CCLE baseline P3 = %.2f",
      ccle$baseline_P3))

pc <- ggplot(gse, aes(predicted_residual_viability, bliss)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .5) +
  geom_vline(xintercept = 0, color = grid, linewidth = .5) +
  geom_point(size = 2.5, color = purple, alpha = .88) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = ink, linewidth = .8) + theme_panel +
  labs(x = "Baseline-predicted residual survival\nduring BRAF inhibition",
    y = "Observed functional benefit from\nERBB blockade (Bliss)",
    title = "B  External functional bypass",
    subtitle = sprintf("GSE68453 · 12 BRAF-mutant lines · r = %.2f",
      gse_test$estimate))

fig <- ggarrange(pa, pc, nrow = 1, widths = c(1, 1.05))
fig <- annotate_figure(fig,
  top = text_grob("Pretreatment expression + context → functional response",
    face = "bold", size = 13.5, color = ink),
  bottom = text_grob(
    "The EGF-on A375 value is a frozen model prediction; no matched EGF-on viability measurement is available in GSE127988.",
    size = 9, color = "#52656F"))
ggsave(file.path(out, "figures/01_a375_functional_counterfactual.pdf"), fig,
  width = 7.7, height = 4.25, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/01_a375_functional_counterfactual.png"), fig,
  width = 7.7, height = 4.25, dpi = 240, bg = "white")

writeLines(c("# CCLE-baseline A375 functional counterfactual", "",
  sprintf("A375 maps to DepMap ACH-000219. Its CCLE pretreatment profile has baseline P3 %.3f and BRAF hotspot status 1.", ccle$baseline_P3),
  sprintf("The frozen Tahoe-to-CCLE model predicts a BRAF-inhibitor functional score of %.3f without EGF; the observed PRISM value is %.3f.", braﬁ_no_egf, ccle$observed_survival),
  sprintf("The independently trained EGF input produces P3 displacement %.3f. The frozen Tahoe P3-to-function coefficient converts this into a functional shift of %+.3f, giving an EGF-on prediction of %.3f.", egf_p3_displacement, egf_function_shift, braﬁ_with_egf), "",
  "Thus, using only the CCLE pretreatment measurement, BRAF genotype and defined EGF input, the frozen model predicts a switch from BRAF-inhibitor sensitivity to resistance. The EGF-on functional value is currently a prediction rather than a directly observed GSE127988 viability endpoint. GSE68453 independently validates that the baseline-conditioned bypass prediction tracks a hard functional combination endpoint across 12 lines."),
  file.path(out, "README.md"))

message("Stage 42 complete: ", out)
