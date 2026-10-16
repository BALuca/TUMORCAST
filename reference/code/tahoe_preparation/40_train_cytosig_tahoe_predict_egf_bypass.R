#!/usr/bin/env Rscript
## Freeze an EGF-bypass model in CytoSig + Tahoe, then test it in GSE127988.
## GSE127988 post-EGF measurements are used only after predictions are frozen.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggpubr)})

root <- "/work/tahoe_prism_v2"
work <- dirname(root)
out <- file.path(root, "results/40_cytosig_tahoe_egf_prediction")
unlink(out, recursive = TRUE, force = TRUE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, name),
  sep = "\t", na = "NA")

## Frozen topology and signed programs from Stage 36.
frozen <- readRDS(file.path(root,
  "artifacts/36_source_aware_braf_erbb_programs/frozen_source_aware_programs.rds"))
membership <- as.data.table(frozen$membership)
signs <- frozen$gene_sign
programs <- c("P1", "P2", "P3", "P4")

## TRAINING ONLY: score all CytoSig profiles, retaining the 51 EGF and 16 NRG1
## experiments for the ligand input and the fixed P2 -> P3 edge magnitude.
cyto <- suppressWarnings(fread(file.path(work,
  "data/01_CytoSig_data/diff.merge.gz"), showProgress = FALSE))
setnames(cyto, 1L, "gene")
cyto <- cyto[gene %chin% membership$gene]
cx <- as.matrix(cyto[, -1]); rownames(cx) <- cyto$gene
score <- sapply(programs, function(p) {
  g <- membership[program == p, gene]
  colMeans(sweep(cx[g, , drop = FALSE], 1, signs[g], "*"), na.rm = TRUE)
})
cm0 <- fread(file.path(work, "data/01_CytoSig_data/meta_info.txt"), skip = 1L,
  header = FALSE,
  col.names = c("profile_id", "treatment", "condition", "duration", "dose"))
cm <- merge(data.table(profile_id = colnames(cx)), cm0, by = "profile_id",
  all.x = TRUE, sort = FALSE)[match(colnames(cx), profile_id)]
cm[, study := sub(".*@((?:GSE[0-9]+|E-MTAB-[0-9]+)).*", "\\1",
  profile_id, perl = TRUE)]
cm[!grepl("^(GSE[0-9]+|E-MTAB-[0-9]+)$", study),
  study := paste0("independent_", .I)]
ligand <- cbind(cm, as.data.table(score))[treatment %chin% c("EGF", "NRG1")]
ligand[, study_weight := 1 / .N, by = study]
ligand[, study_weight := study_weight / sum(study_weight)]

## The graph supplies the direction P2 -> P3. Only its numerical slope and
## the EGF input magnitude are estimated. Study weights stop large series from
## dominating the estimates.
edge_fit <- lm(P3 ~ 0 + P2, data = ligand, weights = study_weight)
edge <- unname(coef(edge_fit)["P2"])
egf <- ligand[treatment == "EGF"]
egf[, egf_weight := (1 / .N), by = study]
egf[, egf_weight := egf_weight / sum(egf_weight)]
egf_input <- egf[, sum(egf_weight * P2)]
expected_p3 <- edge * egf_input
edge_cor <- cor.test(ligand$P2, ligand$P3)

training_stats <- data.table(
  training_source = "CytoSig control-normalized perturbations",
  n_egf = egf[, .N], n_nrg1 = ligand[treatment == "NRG1", .N],
  n_studies = ligand[, uniqueN(study)],
  egf_P2_input_study_balanced = egf_input,
  P2_to_P3_slope = edge,
  P2_P3_pearson = unname(edge_cor$estimate),
  P2_P3_p_value = edge_cor$p.value,
  predicted_P3_after_EGF = expected_p3)
write_table(ligand, "cytosig_ligand_training_scores.tsv")
write_table(training_stats, "frozen_egf_model.tsv")

## Tahoe supplies the independently learned baseline/drug context. These
## quantities are audited here but are not re-estimated in the held-out data.
tahoe <- fread(file.path(root,
  "results/38_mockup_baseline_conditioned_braf_transfer/tahoe_baseline_response_training_data.tsv"))
tahoe_context <- rbindlist(list(
  tahoe[compound == "vemurafenib" & program == "P1",
    .(training_source = "Tahoe", context = "BRAF inhibitor / P1",
      n_lines = .N, mean_response = mean(response_score),
      braf_hotspot_difference = mean(response_score[braf_hotspot == 1]) -
        mean(response_score[braf_hotspot == 0]))],
  tahoe[compound == "cobimetinib" & program == "P3",
    .(training_source = "Tahoe", context = "MEK inhibitor / P3",
      n_lines = .N, mean_response = mean(response_score),
      braf_hotspot_difference = mean(response_score[braf_hotspot == 1]) -
        mean(response_score[braf_hotspot == 0]))]))
write_table(tahoe_context, "tahoe_context_training_audit.tsv")

## FREEZE PREDICTIONS. The only test-set facts used here are available at
## baseline/design time: A375 is BRAF V600E; EGF is added; cobimetinib either
## is or is not present. The no-tree model gets the same inputs but has no
## allowed P2 -> P3 route, hence predicts no post-EGF P3 displacement.
design <- CJ(arm = c("Vem", "Vem + Cobi"), egf_h = c(.5, 1, 2, 3, 4, 8),
  biological_set = c(3L, 4L))
design[, `:=`(braf_v600 = 1L, egf_added = 1L,
  mek_edge_open = as.integer(arm == "Vem"))]
design[, tree_prediction := braf_v600 * egf_added * mek_edge_open * expected_p3]
design[, no_tree_prediction := 0]
write_table(design, "frozen_predictions_before_outcome_access.tsv")

## HELD TEST: only now attach measured post-EGF P3 scores from GSE127988.
observed <- fread(file.path(root,
  "results/37_source_aware_bypass_validation/gse127988_time_course_means.tsv"))[
    program == "P3" & egf_h > 0,
    .(arm = fifelse(as.character(drug_background) == "Vemurafenib",
      "Vem", "Vem + Cobi"), egf_h, biological_set, observed_P3 = score)]
test <- merge(design, observed, by = c("arm", "egf_h", "biological_set"),
  all.x = TRUE)
stopifnot(!anyNA(test$observed_P3))
write_table(test, "held_gse127988_predictions.tsv")

rmse <- function(y, p) sqrt(mean((y - p)^2))
metrics <- rbindlist(list(
  test[, .(model = "No tree", n = .N, rmse = rmse(observed_P3,
    no_tree_prediction), mean_prediction = mean(no_tree_prediction),
    mean_observed = mean(observed_P3)), by = arm],
  test[, .(model = "Fixed tree", n = .N, rmse = rmse(observed_P3,
    tree_prediction), mean_prediction = mean(tree_prediction),
    mean_observed = mean(observed_P3)), by = arm]))
metrics[, rmse_improvement_over_no_tree :=
  metrics[model == "No tree", rmse][match(arm, metrics[model == "No tree", arm])] - rmse]
write_table(metrics, "held_gse127988_model_comparison.tsv")

## Integrated response is the primary bypass endpoint.
trapz <- function(t, y) sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
auc <- test[order(egf_h), .(
  observed_auc = trapz(c(0, egf_h), c(0, observed_P3)),
  tree_predicted_auc = trapz(c(0, egf_h), c(0, tree_prediction)),
  no_tree_predicted_auc = 0), by = .(arm, biological_set)]
write_table(auc, "held_gse127988_bypass_auc.tsv")

## Compact audit figure.
theme_panel <- theme_classic(base_size = 10) +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold", size = 10.5))
p1 <- ggplot(ligand, aes(P2, P3, size = study_weight, color = treatment)) +
  geom_hline(yintercept = 0, color = "#D9DEDF") +
  geom_vline(xintercept = 0, color = "#D9DEDF") +
  geom_point(alpha = .72) +
  geom_abline(slope = edge, intercept = 0, color = "#263840", linewidth = .8) +
  scale_color_manual(values = c(EGF = "#168C86", NRG1 = "#377EB8")) +
  scale_size(range = c(1.4, 3.3), guide = "none") + theme_panel +
  theme(legend.position = "bottom") +
  labs(x = "CytoSig P2 ligand input", y = "CytoSig P3 output",
    title = sprintf("Train fixed edge, n=67; r=%.2f", edge_cor$estimate),
    color = NULL)

plot_test <- melt(test,
  id.vars = c("arm", "egf_h", "biological_set", "observed_P3"),
  measure.vars = c("tree_prediction", "no_tree_prediction"),
  variable.name = "model", value.name = "prediction")
plot_test[, model := factor(model,
  levels = c("tree_prediction", "no_tree_prediction"),
  labels = c("Fixed tree", "No tree"))]
p2 <- ggplot(plot_test, aes(egf_h, observed_P3,
  group = biological_set, color = factor(biological_set))) +
  geom_hline(yintercept = 0, color = "#D9DEDF") +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  geom_line(aes(y = prediction, group = interaction(model, biological_set),
    linetype = model), color = "#263840", linewidth = .8) +
  facet_wrap(~arm) +
  scale_color_manual(values = c("3" = "#D88736", "4" = "#5B4775"),
    name = "Biological day") +
  scale_linetype_manual(values = c("Fixed tree" = "solid", "No tree" = "dotted")) +
  theme_panel + theme(legend.position = "bottom", aspect.ratio = .82) +
  labs(x = "Hours after EGF", y = "Held P3 response",
    title = "Predict EGF bypass and MEK blockade", linetype = NULL)

p3 <- ggplot(auc, aes(arm, observed_auc, group = biological_set,
  color = factor(biological_set))) +
  geom_hline(yintercept = 0, color = "#D9DEDF") +
  geom_line(color = "#7D898E", linewidth = .7) + geom_point(size = 2.6) +
  geom_point(aes(y = tree_predicted_auc), shape = 4, stroke = 1.2,
    size = 3.2, color = "#263840") + theme_panel +
  scale_color_manual(values = c("3" = "#D88736", "4" = "#5B4775"),
    guide = "none") +
  labs(x = NULL, y = "Integrated P3 response",
    title = "Held bypass endpoint\n× = frozen prediction")

fig <- ggarrange(p1, p2, p3, nrow = 1, widths = c(1, 1.65, .85))
fig <- annotate_figure(fig, top = text_grob(
  "CytoSig + Tahoe training predicts EGF bypass in held-out A375",
  face = "bold", size = 13))
ggsave(file.path(out, "figures/01_frozen_training_held_egf_test.pdf"), fig,
  width = 10.4, height = 3.8, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/01_frozen_training_held_egf_test.png"), fig,
  width = 10.4, height = 3.8, dpi = 220, bg = "white")

vem <- metrics[arm == "Vem"]
writeLines(c("# CytoSig + Tahoe trained EGF bypass prediction", "",
  paste0("The topology was fixed before fitting. CytoSig estimates the EGF P2 input and the P2-to-P3 edge from 51 EGF and 16 NRG1 signatures. Tahoe supplies the BRAF-mutant inhibitor context and the MEK-blocked P3 branch. No GSE127988 post-EGF value enters training."), "",
  sprintf(paste0("In held-out A375, the fixed tree predicts a positive P3 bypass response under vemurafenib (%.3f) and no routed response when cobimetinib closes the MEK-ERK edge. The observed mean responses are %.3f and %.3f, respectively."), expected_p3, test[arm == "Vem", mean(observed_P3)], test[arm == "Vem + Cobi", mean(observed_P3)]),
  sprintf("Under vemurafenib, RMSE is %.3f for the fixed tree and %.3f for the same-input no-tree model.", vem[model == "Fixed tree", rmse], vem[model == "No tree", rmse]), "",
  "This is a held-condition pilot in one BRAF V600E line with two biological days. It tests whether the fixed tree routes an independently learned EGF input to bypass and predicts MEK blockade. It does not by itself test cross-line variation in the size of EGF bypass."),
  file.path(out, "README.md"))

message("Stage 40 complete: ", out)
