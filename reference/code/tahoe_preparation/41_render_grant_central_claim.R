#!/usr/bin/env Rscript
## Grant-facing synthesis of the evidence for baseline-conditioned response
## and context-dependent bypass. All panels use existing frozen analyses.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggpubr)})

root <- "/work/tahoe_prism_v2"
out <- file.path(root, "results/41_grant_central_claim")
unlink(out, recursive = TRUE, force = TRUE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, name),
  sep = "\t", na = "NA")

r38 <- file.path(root, "results/38_mockup_baseline_conditioned_braf_transfer")
r40 <- file.path(root, "results/40_cytosig_tahoe_egf_prediction")

## A. External baseline prediction of BRAF-inhibitor viability.
ccle <- fread(file.path(r38, "tahoe_to_ccle_transfer_predictions.tsv"))[
  primary_tissue == "skin" & transfer_model == "P1"]
stopifnot(nrow(ccle) == 38L)
ccle_test <- cor.test(ccle$predicted_survival, ccle$observed_survival)
ccle[, BRAF := factor(braf_hotspot, levels = c(0, 1),
  labels = c("BRAF wild type", "BRAF hotspot"))]
write_table(ccle, "panel_A_ccle_melanoma_predictions.tsv")

## B. Independent training of the fixed EGF/ERBB -> MEK/ERK edge.
ligand <- fread(file.path(r40, "cytosig_ligand_training_scores.tsv"))
edge <- fread(file.path(r40, "frozen_egf_model.tsv"))
stopifnot(nrow(ligand) == 67L)
write_table(ligand, "panel_B_cytosig_ligand_training.tsv")

## C. Held EGF challenge; predictions were written before outcomes were read.
auc <- fread(file.path(r40, "held_gse127988_bypass_auc.tsv"))
auc[, arm_label := factor(arm, levels = c("Vem", "Vem + Cobi"),
  labels = c("BRAF inhibitor", "BRAF + MEK\ninhibitors"))]
auc[, day := factor(biological_set)]
write_table(auc, "panel_C_gse127988_held_bypass.tsv")

## D. External hard functional outcome: combination benefit in 12 lines.
functional <- fread(file.path(r38, "tahoe_to_gse68453_function_predictions.tsv"))[
  normalization == "program_scaled_P3"]
stopifnot(nrow(functional) == 12L)
functional_test <- cor.test(functional$predicted_residual_viability,
  functional$bliss)
write_table(functional, "panel_D_gse68453_functional_predictions.tsv")

ink <- "#263840"
grid <- "#D9DEDF"
orange <- "#D88736"
purple <- "#5B4775"
teal <- "#168C86"
blue <- "#377EB8"
theme_grant <- theme_classic(base_size = 10.5) +
  theme(aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 10.8, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 8.8, color = "#52656F",
      margin = margin(b = 6)),
    axis.title = element_text(size = 9.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    plot.margin = margin(6, 7, 5, 5))

pa <- ggplot(ccle, aes(predicted_survival, observed_survival, color = BRAF)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .4) +
  geom_vline(xintercept = 0, color = grid, linewidth = .4) +
  geom_point(size = 2.15, alpha = .86) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = ink, linewidth = .8) +
  scale_color_manual(values = c("BRAF wild type" = "#9AA7AD",
    "BRAF hotspot" = orange)) + theme_grant +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Prediction from baseline expression\nand BRAF status",
    y = "Observed BRAF-inhibitor viability",
    title = "A  No-EGF drug response",
    subtitle = sprintf("CCLE/PRISM · 38 held melanoma lines · r = %.2f",
      ccle_test$estimate))

pb <- ggplot(ligand, aes(P2, P3, color = treatment)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .4) +
  geom_vline(xintercept = 0, color = grid, linewidth = .4) +
  geom_point(size = 1.95, alpha = .76) +
  geom_abline(slope = edge$P2_to_P3_slope, intercept = 0,
    color = ink, linewidth = .8) +
  scale_color_manual(values = c(EGF = teal, NRG1 = blue)) + theme_grant +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "P2: EGF–ERBB input", y = "P3: MEK–ERK output",
    title = "B  Train the EGF→MEK–ERK edge",
    subtitle = sprintf("CytoSig · 67 signatures, 25 studies · r = %.2f",
      edge$P2_P3_pearson))

pc <- ggplot(auc, aes(arm_label, observed_auc, group = day, color = day)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .4) +
  geom_line(color = "#879398", linewidth = .75) +
  geom_point(size = 2.7) +
  geom_point(aes(y = tree_predicted_auc), shape = 4, stroke = 1.35,
    size = 3.7, color = ink) + theme_grant +
  scale_color_manual(values = c("3" = orange, "4" = purple),
    labels = c("Day 1", "Day 2")) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Integrated held P3 response",
    title = "C  Held EGF-bypass prediction",
    subtitle = "A375 EGF challenge · × frozen prediction")

pd <- ggplot(functional,
  aes(predicted_residual_viability, bliss)) +
  geom_hline(yintercept = 0, color = grid, linewidth = .4) +
  geom_vline(xintercept = 0, color = grid, linewidth = .4) +
  geom_point(size = 2.35, color = purple, alpha = .88) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = ink, linewidth = .8) + theme_grant +
  labs(x = "Baseline-predicted residual survival\nduring BRAF inhibition",
    y = "Observed benefit from\nERBB-pathway blockade (Bliss)",
    title = "D  Functional bypass prediction",
    subtitle = sprintf("GSE68453 · 12 BRAF-mutant lines · r = %.2f",
      functional_test$estimate))

fig <- ggarrange(pa, pb, pc, pd, nrow = 2, ncol = 2,
  common.legend = FALSE, align = "hv")
fig <- annotate_figure(fig,
  top = text_grob("Pretreatment state and defined inputs predict BRAF response and bypass",
    face = "bold", size = 14, color = ink),
  bottom = text_grob(
    "Expression + genotype define the starting state; treatment and EGF define the inputs; the fixed tree propagates their effects.",
    size = 9.3, color = "#52656F"))
ggsave(file.path(out, "figures/01_central_claim_evidence.pdf"), fig,
  width = 7.8, height = 7.5, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/01_central_claim_evidence.png"), fig,
  width = 7.8, height = 7.5, dpi = 240, bg = "white")

stats <- data.table(
  panel = c("A", "B", "C", "D"),
  role = c("No-EGF functional response", "Independent edge training",
    "Held exogenous-EGF molecular bypass", "Hard functional bypass"),
  dataset = c("CCLE/PRISM", "CytoSig", "GSE127988", "GSE68453"),
  n = c(nrow(ccle), nrow(ligand), uniqueN(auc$biological_set), nrow(functional)),
  primary_result = c(
    sprintf("Pearson r = %.3f; p = %.3g", ccle_test$estimate, ccle_test$p.value),
    sprintf("Pearson r = %.3f; p = %.3g", edge$P2_P3_pearson,
      edge$P2_P3_p_value),
    sprintf("predicted AUC %.3f; observed %.3f and %.3f; MEK-blocked %.3f and %.3f",
      unique(auc[arm == "Vem", tree_predicted_auc]),
      auc[arm == "Vem", observed_auc][1], auc[arm == "Vem", observed_auc][2],
      auc[arm == "Vem + Cobi", observed_auc][1],
      auc[arm == "Vem + Cobi", observed_auc][2]),
    sprintf("Pearson r = %.3f; p = %.3g", functional_test$estimate,
      functional_test$p.value)))
write_table(stats, "panel_statistics.tsv")

writeLines(c("# Evidence for the grant's central prediction claim", "",
  "The figure separates training from validation and shows four linked results.", "",
  sprintf("**A. No added EGF.** Baseline expression and BRAF status predict held-out BRAF-inhibitor viability in 38 CCLE/PRISM melanoma lines (r = %.2f, p = %.3g).", ccle_test$estimate, ccle_test$p.value),
  sprintf("**B. Independent graph training.** Across 67 CytoSig EGF/NRG1 signatures, the frozen P2 ligand input is coupled to P3 MEK-ERK output (r = %.2f, p = %.3g).", edge$P2_P3_pearson, edge$P2_P3_p_value),
  sprintf("**C. Added EGF.** The frozen model predicts an integrated P3 bypass response of %.2f in held A375; the two observed biological days are %.2f and %.2f. Closing the edge with cobimetinib reduces the observed response to %.2f and %.2f.", unique(auc[arm == "Vem", tree_predicted_auc]), auc[arm == "Vem", observed_auc][1], auc[arm == "Vem", observed_auc][2], auc[arm == "Vem + Cobi", observed_auc][1], auc[arm == "Vem + Cobi", observed_auc][2]),
  sprintf("**D. Functional bypass.** The baseline-conditioned prediction identifies benefit from ERBB-pathway blockade across 12 BRAF-mutant GSE68453 melanoma lines (Bliss r = %.2f, p = %.3g).", functional_test$estimate, functional_test$p.value), "",
  "Panels A and D are external functional tests. Panel C is a held molecular bypass test with an independently trained input and edge. GSE127988 contains one cell line, so it tests context transfer rather than cross-line variation in EGF-bypass magnitude."),
  file.path(out, "README.md"))

message("Stage 41 complete: ", out)
