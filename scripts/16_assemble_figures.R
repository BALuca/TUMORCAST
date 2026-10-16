#!/usr/bin/env Rscript
# Recreate the supplied four-panel joined figure using the regenerated vector
# panels. Geometry follows the original assembly: 3-pt gaps, 12-pt header, 17 cm.
# SVG composition and Cairo export replace the historical Python/PyMuPDF step.

suppressPackageStartupMessages({
  library(xml2)
  library(rsvg)
})

panel_dir <- "results/15_plot_polished_panels/panels"
output_dir <- "results/16_assemble_figures"
panels <- c("B2_2A_shared_decoder_bars", "B2_3A_measured_input_scatter",
  "B2_3A_predicted_input_scatter", "B2_2B_lapatinib_HGF_program_scatter")
gap_pt <- 3
header_pt <- 12
width_cm <- 17
resolution <- 600
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)

# Keep vector paths intact. Prefix the glyph/clip identifiers so definitions in
# different source SVGs cannot collide when placed on the same page.
sources <- lapply(seq_along(panels), function(i) {
  path <- file.path(panel_dir, paste0(panels[i], ".svg"))
  stopifnot(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  prefix <- paste0("panel", i, "_")
  text <- gsub('id="', paste0('id="', prefix), text, fixed = TRUE)
  text <- gsub('href="#', paste0('href="#', prefix), text, fixed = TRUE)
  text <- gsub("url(#", paste0("url(#", prefix), text, fixed = TRUE)
  document <- read_xml(text)
  viewport <- as.numeric(strsplit(xml_attr(document, "viewBox"), "[ ,]+")[[1]])
  stopifnot(length(viewport) == 4, all(is.finite(viewport)), all(viewport[3:4] > 0))
  list(document = document, width = viewport[3], height = viewport[4])
})
native_width <- sum(vapply(sources, `[[`, numeric(1), "width")) + gap_pt * (length(panels) - 1)
native_height <- max(vapply(sources, `[[`, numeric(1), "height")) + header_pt
page_width_pt <- width_cm / 2.54 * 72
scale <- page_width_pt / native_width
page_height_pt <- native_height * scale
figure <- xml_new_root("svg", xmlns = "http://www.w3.org/2000/svg",
  width = paste0(page_width_pt, "pt"), height = paste0(page_height_pt, "pt"),
  viewBox = paste(0, 0, native_width, native_height))
xml_add_child(figure, "rect", x = "0", y = "0", width = as.character(native_width),
  height = as.character(native_height), fill = "white")
x <- 0
layout <- vector("list", length(panels))
for (i in seq_along(panels)) {
  source <- sources[[i]]
  panel <- xml_root(source$document)
  xml_set_attrs(panel, c(x = as.character(x), y = as.character(header_pt),
    width = as.character(source$width), height = as.character(source$height)))
  xml_add_child(figure, panel)
  letter <- xml_add_child(figure, "text", x = as.character(x), y = "9",
    `font-family` = "Helvetica", `font-weight` = "bold", `font-size` = "9", fill = "black")
  xml_text(letter) <- LETTERS[i]
  layout[[i]] <- data.frame(panel = panels[i], x_pt = x, y_pt = header_pt,
    width_pt = source$width, height_pt = source$height, whole_figure_scale = scale)
  x <- x + source$width + gap_pt
}
svg_path <- file.path(output_dir, "B2_Figure_2.svg")
write_xml(figure, svg_path)
rsvg_pdf(svg_path, file.path(output_dir, "B2_Figure_2.pdf"))
rsvg_png(svg_path, file.path(output_dir, "B2_Figure_2.png"),
  width = round(page_width_pt / 72 * resolution), height = round(page_height_pt / 72 * resolution))
write.table(do.call(rbind, layout), file.path(output_dir, "joined_panel_layout.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
# Reproduce the molecular-bridge and ligand-rescue support composite.
# !/usr/bin/env Rscript
# Assemble existing evidence only. No model fitting, selection, or grant edits.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
root <- "results/12_fit_ligand_rescue"
figure_root <- file.path(root, "figures/grant_support")
data_root <- file.path(root, "data/grant_support")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
baseline_path <- "results/08_fit_single_drug_bridge/data/baseline_predictions_frozen.tsv"
measured_path <- "results/08_fit_single_drug_bridge/data/measured_predictions_frozen.tsv"
opposition_path <- "results/14_plot_model_diagnostics/AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv"
comparison_path <- file.path(root, "data/prediction_beside_Wilson.tsv")
baseline <- fread(baseline_path)[drug == "Lapatinib ditosylate" & dose == 5 & evaluation == "held_cell"]
measured <- fread(measured_path)[drug == "Lapatinib ditosylate" & dose == 5 & evaluation == "held_cell"]
stopifnot(nrow(baseline) == 441, nrow(measured) == 30, !anyDuplicated(baseline$depmap_id))
baseline[, RNA_overlap := depmap_id %in% measured$depmap_id]
stopifnot(sum(baseline$RNA_overlap) == 30)
r_all <- cor(baseline$observed_log2, baseline$predicted_log2)
r_overlap <- baseline[RNA_overlap == TRUE, cor(observed_log2, predicted_log2)]
writeLines(sprintf(
  "All 441: r=%.6f; same RNA-overlap 30: r=%.6f; endpoint=log2 PRISM abundance; dose=5 uM",
  r_all, r_overlap
), file.path(data_root, "single_drug_endpoint.txt"))
fwrite(baseline, file.path(data_root, "panel_A_single_drug_predictions.tsv"), sep = "\t")
theme_base <- theme_classic(base_size = 8) + theme(
  plot.title = element_text(size = 9, face = "bold"),
  plot.subtitle = element_text(size = 7), plot.margin = margin(4, 4, 4, 4), axis.title = element_text(size = 8)
)
a <- ggplot(baseline, aes(observed_log2, predicted_log2)) +
  geom_point(data = baseline[RNA_overlap == FALSE], size = .65, color = "#A8AFB5", alpha = .6) +
  geom_point(data = baseline[RNA_overlap == TRUE], size = 1.1, color = "#2769A5") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .35, color = "black") +
  labs(
    title = "Lapatinib bridge", subtitle = sprintf("441 cells: r = %.2f\nRNA-overlap 30 (blue): r = %.2f", r_all, r_overlap),
    x = "Observed log2 PRISM abundance", y = "Predicted log2 abundance"
  ) +
  theme_base +
  theme(aspect.ratio = 1)

# Selection is inherited from single-drug associations, with no new HGF selection.
opposition <- fread(opposition_path)[strong_lapatinib_association == TRUE]
stopifnot(nrow(opposition) == 7, all(abs(opposition$lapatinib_sensitivity_r) >= .5))
opposition[, opposed_direction := lapatinib_displacement * HGF_displacement < 0]
stopifnot(sum(opposition$opposed_direction) == 6)
setorder(opposition, lapatinib_sensitivity_r)
opposition[, program := factor(program, levels = program)]
fwrite(opposition, file.path(data_root, "panel_B_program_opposition.tsv"), sep = "\t")
plot_b <- melt(opposition,
  id.vars = "program", measure.vars = c("lapatinib_displacement", "HGF_displacement"),
  variable.name = "condition", value.name = "displacement"
)
plot_b[, condition := factor(condition,
  levels = c("lapatinib_displacement", "HGF_displacement"),
  labels = c("Lapatinib", "HGF")
)]
b <- ggplot(plot_b, aes(program, displacement, fill = condition)) +
  geom_hline(yintercept = 0, linewidth = .25, color = "#555555") +
  geom_col(position = position_dodge(width = .75), width = .65) +
  scale_fill_manual(values = c("Lapatinib" = "#2769A5", "HGF" = "#D86A25")) +
  labs(
    title = "Program opposition", subtitle = "Opposite shifts in 6 of 7 programs\nSingle-drug association |r| >= 0.5",
    x = NULL, y = "Program displacement", fill = NULL
  ) +
  theme_base +
  theme(
    legend.position = "bottom", legend.text = element_text(size = 7),
    legend.key.size = grid::unit(.3, "cm"), legend.margin = margin(0, 0, 0, 0), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )

comparison <- fread(comparison_path)
comparison[, panel := factor(panel,
  levels = c("Prediction\n5 uM, five-day baseline", "Wilson measurement\n1 uM, 12 days"),
  labels = c("Prediction\n5 uM, five-day baseline", "Wilson\n1 uM, 12 days")
)]
comparison[, condition := factor(condition, levels = c("Lapatinib", "Lapatinib + HGF"), labels = c("Lap", "Lap + HGF"))]
stopifnot(nrow(comparison) == 4, all(is.finite(comparison$abundance)))
fwrite(comparison, file.path(data_root, "panel_C_prediction_and_measurement.tsv"), sep = "\t")
c <- ggplot(comparison, aes(condition, abundance, fill = condition)) +
  geom_col(width = .65) +
  geom_errorbar(
    data = comparison[is.finite(sd)],
    aes(ymin = abundance - sd, ymax = abundance + sd), width = .15, linewidth = .35
  ) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c("Lap" = "#2769A5", "Lap + HGF" = "#D86A25")) +
  labs(
    title = "HGF rescue hypothesis", subtitle = "HGF and Wilson outcomes excluded from fitting",
    x = NULL, y = "Abundance (% control)"
  ) +
  theme_base +
  theme(
    legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 7), strip.background = element_blank()
  )
figure <- (a + b + c) + plot_layout(widths = c(1, 1, 1.35)) + plot_annotation(
  tag_levels = "A",
  theme = theme(plot.tag = element_text(face = "bold", size = 10))
)
ggsave(file.path(figure_root, "v106_v144_grant_support.pdf"), figure, width = 7.5, height = 2.9, device = cairo_pdf)
ggsave(file.path(figure_root, "v106_v144_grant_support.png"), figure, width = 7.5, height = 2.9, dpi = 260)
for (extension in c("pdf", "png")) {
  copied <- file.copy(file.path(figure_root, paste0("v106_v144_grant_support.", extension)),
    file.path(figure_root, paste0("molecular_bridge_and_ligand_rescue_support.", extension)),
    overwrite = TRUE)
  stopifnot(copied)
}
claims <- data.table(
  panel = c("A", "B", "C"), supported_claim = c(
    "Frozen single-drug molecular bridge predicts lapatinib log2 abundance across held-out recipients at 5 uM",
    "Transported HGF response opposes lapatinib shifts in six of seven programs selected by strong single-drug sensitivity association",
    "Bounded rescue pilot trained on NRG1 gives an HGF reversal prediction without HGF functional labels"
  ),
  limits = c(
    "Retrospective development; r=.34 across 441, r=.59 in the specified 30 RNA-overlap recipients",
    "Descriptive association; AU565 HGF RNA response is reconstructed, not measured; amplitude comparability remains uncertain",
    "Hypothesis, not quantitative validation: dose/duration differ from Wilson; neutral HGF-alone assumed; constant-recovery comparator matches or exceeds held interaction correlation"
  )
)
fwrite(claims, file.path(data_root, "supported_claims_and_limits.tsv"), sep = "\t")
paths <- c("scripts/16_assemble_figures.R", baseline_path, measured_path, opposition_path, comparison_path)
fwrite(data.table(path = paths, sha256 = vapply(paths, digest, character(1), file = TRUE, algo = "sha256")),
  file.path(data_root, "source_hashes.tsv"),
  sep = "\t"
)
cat("Joined figure and grant support figure saved.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
sink()
