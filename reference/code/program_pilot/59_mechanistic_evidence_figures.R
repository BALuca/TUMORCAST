# Display observed genetic-intervention and drug-combination measurements for dataset collection 59.
# No fitting, program discovery or outcome selection; source SD remains source SD.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("/work/src/R/lib_heatmap.R")
root <- "/work/results/program_pilot/59_mechanistic_evidence_datasets"
data_root <- file.path(root, "data")
fig_root <- file.path(root, "figures")
dir.create(fig_root, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
  legend.position = "bottom", plot.title = element_text(face = "bold"),
  plot.caption = element_text(size = 8, hjust = 0)))

# Genetic control: all four source experiments at every dose, not selected responsive doses.
myc <- fread(file.path(data_root, "mek_measurements.tsv.gz"))[
  source_panel %in% c("Fig8F", "Fig8G", "Fig8H")]
stopifnot(all(myc$unit == "percent genotype-matched vehicle"), !anyNA(myc$value))
myc[, state := factor(biological_state, levels = c("LacZ_control", "MYC_overexpression"),
                       labels = c("LacZ control", "MYC overexpression"))]
myc_summary <- myc[, .(value = mean(value), sd = sd(value), n_experiments = uniqueN(replicate)),
                   by = .(source_panel, cell_line, treatment, state, dose_A_uM)]
stopifnot(all(myc_summary$n_experiments == 4))
plots <- list()
for (panel in c("Fig8F", "Fig8G", "Fig8H")) {
  d <- copy(myc[source_panel == panel])
  s <- copy(myc_summary[source_panel == panel])
  positive_doses <- sort(unique(d[dose_A_uM > 0, dose_A_uM]))
  vehicle_x <- log10(min(positive_doses)) - log10(4)
  d[, dose_position := fifelse(dose_A_uM == 0, vehicle_x, log10(pmax(dose_A_uM, 1e-12)))]
  s[, dose_position := fifelse(dose_A_uM == 0, vehicle_x, log10(pmax(dose_A_uM, 1e-12)))]
  selected <- unique(round(seq(1, length(positive_doses), length.out = 4)))
  breaks <- c(vehicle_x, log10(positive_doses[selected]))
  labels <- c("Veh", sprintf("%g", positive_doses[selected]))
  plots[[panel]] <- ggplot(d, aes(dose_position, value, color = state)) +
    geom_hline(yintercept = 100, color = "grey70", linetype = 2) +
    geom_line(aes(group = interaction(state, replicate)), linewidth = 0.25, alpha = 0.25) +
    geom_point(size = 1.1, alpha = 0.5) +
    geom_errorbar(data = s, aes(ymin = value - sd, ymax = value + sd), width = 0.06, linewidth = 0.4) +
    geom_line(data = s, aes(group = state), linewidth = 0.85) +
    geom_point(data = s, size = 2) +
    scale_color_manual(values = c("LacZ control" = "#3B6FB6", "MYC overexpression" = "#C64B32")) +
    scale_x_continuous(breaks = breaks, labels = labels) +
    labs(title = paste(unique(d$cell_line), "—", unique(d$treatment)),
         x = "Drug concentration (µM; log scale)",
         y = "Viability (% matched genotype vehicle)", color = NULL) +
    theme(axis.text.x = element_text(size = 8))
}
myc_plot <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = "Does maintaining MYC change drug sensitivity?",
    subtitle = "MYC overexpression versus LacZ control; 10-day crystal-violet outgrowth assay",
    caption = "Points/thin lines: four source experiments. Thick lines: means; bars: ±1 SD.\nVehicle is shown separately to the left of the log-dose axis. Each genotype is normalized to its own vehicle.") &
  theme(legend.position = "bottom")
ggsave(file.path(fig_root, "MYC_genetic_control_dose_curves.pdf"), myc_plot, width = 10.8, height = 4.3, device = cairo_pdf)
ggsave(file.path(fig_root, "MYC_genetic_control_dose_curves.png"), myc_plot, width = 10.8, height = 4.3, dpi = 170)
fwrite(myc_summary, file.path(data_root, "figure_MYC_genetic_control_summary.tsv"), sep = "\t")

# Orthogonal functional readouts: cell number and ATP remain distinct measurements.
diaz <- fread(file.path(data_root, "diaz_functional_condition_summaries.tsv.gz"))[
  assay %in% c("InCell", "CTG")]
stopifnot(!anyDuplicated(diaz, by = c("cell_line", "arm", "time_h", "assay")))
arm_levels <- c("T", "M", "W", "TM", "TW", "MW", "TMW")
arm_colors <- c(T = "#377EB8", M = "#4DAF4A", W = "#E69F00", TM = "#984EA3",
                TW = "#E41A1C", MW = "#A65628", TMW = "#202020")
diaz[, arm := factor(arm, levels = arm_levels)]
functional_plots <- list()
for (assay_name in c("InCell", "CTG")) {
  d <- diaz[assay == assay_name]
  label <- if (assay_name == "InCell") "Counted nuclei (% same-time vehicle)" else "ATP luminescence (% same-time vehicle)"
  title <- if (assay_name == "InCell") "Direct cell-number measurements after single and combined drugs" else
    "ATP measurements from the same treatment design"
  p <- ggplot(d, aes(time_h, value, color = arm, group = arm)) +
    geom_hline(yintercept = 100, color = "grey65", linetype = 2) +
    geom_ribbon(data = d[is.finite(sd)], aes(ymin = value - sd, ymax = value + sd, fill = arm),
                color = NA, alpha = 0.07, show.legend = FALSE) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.7) +
    facet_wrap(~cell_line, nrow = 1) +
    scale_color_manual(values = arm_colors, drop = FALSE) +
    scale_fill_manual(values = arm_colors, drop = FALSE) +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 24, 48)) +
    labs(title = title, x = "Hours after treatment", y = label, color = "Treatment",
      caption = paste0("T: tamoxifen 20 µM; M: mefloquine 10 µM; W: withaferin A 5 µM. Letters denote combinations.\n",
        "Points/lines: published condition means; ribbons: source SD where available. Missing SD is not treated as zero.\n",
        "The dashed 100% line is the normalization reference; it is not an independently plotted vehicle trajectory."))
  functional_plots[[assay_name]] <- p
  file_name <- if (assay_name == "InCell") "diaz_counted_nuclei_timecourses" else "diaz_ATP_timecourses"
  ggsave(file.path(fig_root, paste0(file_name, ".pdf")), p, width = 8.6, height = 4.8, device = cairo_pdf)
  ggsave(file.path(fig_root, paste0(file_name, ".png")), p, width = 8.6, height = 4.8, dpi = 170)
}

# Fixed bank43 projection: preserve program order and actual signed displacement, without row scaling.
scores <- fread(file.path(data_root, "diaz_bank43_scores.tsv.gz"))
meta <- fread(file.path(data_root, "diaz_samples.tsv.gz"))
stopifnot(!anyDuplicated(meta$sample_id))
scores <- merge(scores, meta[, .(sample_id, cell_line, arm, time_h)], by = "sample_id", all.x = TRUE)
stopifnot(!anyNA(scores$cell_line))
scores <- scores[time_h %in% c(3, 6, 9, 12) & arm %in% arm_levels]
means <- scores[, .(displacement = mean(displacement), n_RNA_samples = .N), by = .(cell_line, arm, time_h, program)]
program_order <- sort(unique(means$program))
color_limit <- max(abs(means$displacement))
for (cell in c("MCF7", "LNCaP")) {
  d <- means[cell_line == cell]
  columns <- unique(d[, .(arm, time_h)])
  columns[, arm := factor(arm, levels = arm_levels)]
  setorder(columns, arm, time_h)
  columns[, label := paste(arm, paste0(time_h, "h"), sep = "_")]
  d[, label := paste(arm, paste0(time_h, "h"), sep = "_")]
  wide <- dcast(d, program ~ label, value.var = "displacement")
  matrix_columns <- columns$label
  mat <- as.matrix(wide[match(program_order, program), ..matrix_columns])
  rownames(mat) <- program_order
  stopifnot(identical(colnames(mat), columns$label), all(is.finite(mat)))
  annotation <- data.frame(Arm = columns$arm, Time = paste0(columns$time_h, "h"), row.names = columns$label)
  hm <- draw_heatmap(mat, col_anno_df = annotation, col_anno_cols = c("Arm", "Time"),
    palette_overrides = list(Arm = arm_colors), scale_rows = FALSE,
    value_range = c(-color_limit, 0, color_limit), value_colors = c("#2166AC", "white", "#B2182B"),
    hm_name = "Mean log2\ndisplacement", show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE, column_split = columns$arm,
    row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 90, title = paste0(cell, ": fixed bank43; measured RNA, mean of source replicates"))
  pdf(file.path(fig_root, paste0("diaz_bank43_", cell, ".pdf")), width = 9, height = 7)
  ComplexHeatmap::draw(hm, padding = grid::unit(c(4, 5, 5, 4), "mm"))
  dev.off()
  png(file.path(fig_root, paste0("diaz_bank43_", cell, ".png")), width = 1530, height = 1190, res = 170)
  ComplexHeatmap::draw(hm, padding = grid::unit(c(4, 5, 5, 4), "mm"))
  dev.off()
}
fwrite(means, file.path(data_root, "figure_diaz_program_means.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(fig_root, "figure_session_info.txt"))
