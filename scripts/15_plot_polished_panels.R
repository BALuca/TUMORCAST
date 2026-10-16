#!/usr/bin/env Rscript
# Exact accepted panel styling is retained for reproduction. Full source data and
# correlation tests accompany the plots; new diagnostic figures follow AGENTS.md.
# Reproduce separate grant panels using ggplot2.
# Each file is a single panel: PDF, SVG and 600-dpi PNG, with a white background.

required <- c("readr", "dplyr", "tidyr", "ggplot2", "ggrepel", "stringr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})
root <- "results/15_plot_polished_panels"
output <- file.path(root, "panels")
audit <- file.path(root, "audit")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
sink(file.path(root, "run.log"), split = TRUE)
programs <- sprintf("P%02d", 1:40)
annotations <- read_tsv("data/raw/annotations/program_labels.tsv", show_col_types = FALSE)
stopifnot(!anyDuplicated(annotations$program))

# All plotted values come from the raw-data rebuild in stages 01-14. Molecular ridge reconstruction
# below is checked against the rebuilt recipient-excluded program predictions.
baseline <- read_tsv("results/08_fit_single_drug_bridge/data/baseline_predictions_frozen.tsv",
  show_col_types = FALSE) |> filter(evaluation == "held_cell")
measured <- read_tsv("results/08_fit_single_drug_bridge/data/measured_predictions_frozen.tsv",
  show_col_types = FALSE) |> filter(evaluation == "held_cell")
molecular <- read_tsv("results/05_score_molecular_inputs/program_displacements.tsv", show_col_types = FALSE)
receptors <- read_tsv("results/05_score_molecular_inputs/receptors.tsv", show_col_types = FALSE)
rna <- molecular |>
  filter(drug == "Lapatinib ditosylate", dose == 5) |>
  inner_join(receptors, by = "depmap_id", relationship = "one-to-one")
x <- as.matrix(rna[, c("ERBB2", "MET")])
y <- as.matrix(rna[, programs])
prediction <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
for (i in seq_len(nrow(y))) {
  train <- seq_len(nrow(y)) != i
  xm <- colMeans(x[train, , drop = FALSE])
  xs <- apply(x[train, , drop = FALSE], 2, sd)
  ym <- colMeans(y[train, , drop = FALSE])
  ys <- pmax(apply(y[train, , drop = FALSE], 2, sd), 1e-8)
  z <- sweep(sweep(x[train, , drop = FALSE], 2, xm, "-"), 2, xs, "/")
  yz <- sweep(sweep(y[train, , drop = FALSE], 2, ym, "-"), 2, ys, "/")
  beta <- solve(crossprod(z) + sum(train) * diag(ncol(x)), crossprod(z, yz))
  prediction[i, ] <- drop(matrix((x[i, ] - xm) / xs, nrow = 1) %*% beta) * ys + ym
}
held <- read_tsv("results/14_plot_model_diagnostics/held_predicted_displacements.tsv",
  show_col_types = FALSE) |> filter(drug == "Lapatinib ditosylate")
check <- bind_cols(tibble(depmap_id = rna$depmap_id), as_tibble(prediction)) |>
  inner_join(held, by = "depmap_id", suffix = c("_new", "_saved"), relationship = "one-to-one")
delta <- max(abs(as.matrix(check[, paste0(programs, "_new")]) - as.matrix(check[, paste0(programs, "_saved")])))
stopifnot(is.finite(delta), delta < 1e-8)
lap_m <- measured |> filter(drug == "Lapatinib ditosylate")
lap_b <- baseline |> filter(drug == "Lapatinib ditosylate", depmap_id %in% lap_m$depmap_id)
stopifnot(nrow(lap_m) == 30L, nrow(lap_b) == 30L)
functional <- lap_m |> inner_join(rna, by = c("depmap_id", "drug", "dose"), relationship = "one-to-one")
program_stats <- bind_rows(lapply(seq_along(programs), function(j) {
  molecular_test <- cor.test(prediction[, j], y[, j], method = "pearson")
  functional_test <- cor.test(functional[[programs[j]]], functional$observed_log2, method = "pearson")
  tibble(
    program = programs[j], prediction_r = unname(molecular_test$estimate),
    prediction_p = molecular_test$p.value, viability_r = unname(functional_test$estimate),
    viability_p = functional_test$p.value
  )
})) |>
  mutate(prediction_q = p.adjust(prediction_p, "BH"), viability_q = p.adjust(viability_p, "BH"))
selected <- program_stats |>
  filter(prediction_r > 0, prediction_q < .05, viability_q < .10) |>
  left_join(annotations, by = "program", relationship = "one-to-one") |>
  arrange(desc(prediction_r))
stopifnot(nrow(selected) == 7L, !anyNA(selected$plot_label))
write_tsv(program_stats, file.path(audit, "program_prediction_and_viability_statistics.tsv"))
write_tsv(selected, file.path(audit, "selected_programs_with_annotations.tsv"))

# Editable plotting style for standalone panels and lettered composite figures.
blue <- "#426E9D"
orange <- "#C36F38"
teal <- "#287C75"
red <- "#A73C3C"
font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))
update_geom_defaults("label", list(family = font_family))
update_geom_defaults("text_repel", list(family = font_family))

# A self-contained Pearson correlation layer, using stats::cor.test for r and two-sided P.
# Both statistics are computed from the plotted observations, without rounding the inputs.
StatCor <- ggproto("StatCor", Stat,
  required_aes = c("x", "y"),
  compute_group = function(data, scales, corner = "top_left") {
    if (nrow(data) < 3L || any(!is.finite(data$x)) || any(!is.finite(data$y))) stop("Invalid correlation data")
    test <- cor.test(data$x, data$y, method = "pearson")
    p_label <- format.pval(test$p.value, digits = 2, eps = 1e-4)
    p_text <- if (startsWith(p_label, "<")) paste("P", p_label) else paste("P =", p_label)
    data.frame(x = if (corner == "top_left") -Inf else Inf, y = if (corner == "top_left") Inf else -Inf, label = sprintf("r = %.2f\n%s", unname(test$estimate), p_text))
  }
)
stat_cor <- function(size = 2.5, show_p = TRUE, corner = "top_left") {
  layer(
    stat = StatCor, geom = "text", position = "identity", inherit.aes = TRUE,
    mapping = if (show_p) aes(label = after_stat(label)) else aes(label = sub("\n.*", "", after_stat(label))),
    params = list(
      size = size, family = font_family, corner = corner,
      hjust = if (corner == "top_left") -.08 else 1.05,
      vjust = if (corner == "top_left") 1.15 else -.15, colour = "black"
    )
  )
}
for (tool in c("gs", "pdftocairo")) {
  if (!nzchar(Sys.which(tool))) stop("Install Ghostscript/Poppler and add ", tool, " to PATH")
}
# Short titles are set explicitly so standalone panels and composites use identical wording.
panel_titles <- c(
  B1_2A_lapatinib_scatter = "Lapatinib response predicted\nfrom untreated cells",
  B1_2B_AU565_prediction_bars = "Predicted HGF protection\nTUMORCAST model",
  B1_2C_Wilson_bars = "Measured HGF protection\nWilson et al.",
  B2_1A_program_predictability_bars = "Response prediction\nfrom untreated expression",
  B2_1B_program_viability_bars = "Program changes\nassociated with viability",
  B2_2A_shared_decoder_bars = "Molecular changes predict\nviability across eight drugs",
  B2_2B_lapatinib_HGF_program_scatter = "Lapatinib and HGF\nprogram changes",
  B2_3A_measured_input_scatter = "Viability from measured\nmolecular changes",
  B2_3A_predicted_input_scatter = "Viability from untreated\ncell expression",
  B2_3B_NRG1_interaction_scatter = "NRG1 interaction prediction"
)
for (program in as.character(selected$program)) {
  panel_titles[paste0("B2_1C_", program, "_displacement_scatter")] <- paste0(program, " response prediction\nfrom untreated cells")
}
paper_theme <- function(square = FALSE) {
  theme_classic(base_size = 8.5, base_family = font_family) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "black", linewidth = .3),
      axis.ticks = element_line(colour = "black", linewidth = .3),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.text = element_text(colour = "black", size = 7.5, lineheight = .9),
      axis.title = element_text(colour = "black", size = 8),
      axis.title.x = element_text(margin = margin(t = 3)),
      axis.title.y = element_text(margin = margin(r = 3)),
      plot.title = element_text(face = "plain", size = 8.5, hjust = .5, lineheight = .95, margin = margin(b = 2)),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 7, hjust = .5, lineheight = 1, margin = margin(b = 3)),
      legend.title = element_blank(), legend.text = element_text(size = 7.5),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 2, .7, .7, unit = "mm"),
      aspect.ratio = if (square) 1 else NULL
    )
}
export_figure <- function(plot, name, width, height, output) {
  destination <- output
  output <- tempfile("figure_export_")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  target <- file.path(output, paste0(name, ".pdf"))
  ggsave(target, plot,
    width = width, height = height, units = "cm", bg = "white",
    device = grDevices::pdf, family = font_family, useDingbats = FALSE
  )
  # PDF declares Helvetica. Ghostscript resolves standard PDF font metrics for raster previews.
  png_target <- file.path(output, paste0(name, ".png"))
  status <- system2("gs", c(
    "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r600",
    "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", shQuote(paste0("-sOutputFile=", png_target)),
    shQuote(target)
  ), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status"))) stop("PNG conversion failed: ", paste(status, collapse = "\n"))
  svg_source <- tempfile(fileext = ".pdf")
  status <- system2("gs", c(
    "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite", "-dEmbedAllFonts=true",
    "-dPDFSETTINGS=/prepress", shQuote(paste0("-sOutputFile=", svg_source)), shQuote(target)
  ),
  stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(status, "status"))) stop("SVG font preparation failed: ", paste(status, collapse = "\n"))
  status <- system2("pdftocairo", c("-svg", shQuote(svg_source), shQuote(file.path(output, paste0(name, ".svg")))),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(status, "status"))) stop("SVG conversion failed: ", paste(status, collapse = "\n"))
  unlink(svg_source)
  for (extension in c("pdf", "svg", "png")) {
    path <- file.path(output, paste0(name, ".", extension))
    if (!file.exists(path) || file.info(path)$size == 0) stop("Empty panel export: ", path)
    if (!file.copy(path, file.path(destination, basename(path)), overwrite = TRUE)) stop("Could not save export: ", path)
  }
}

manifest <- list()
square_checks <- list()
panel_plots <- list()
save_panel <- function(plot, name, width, height, square = FALSE) {
  if (!name %in% names(panel_titles)) stop("Missing panel title: ", name)
  panel_subtitles <- c(
    B1_2B_AU565_prediction_bars = "AU565 cells, 5 days",
    B1_2C_Wilson_bars = "AU565 cells, 12 days"
  )
  subtitle <- if (name %in% names(panel_subtitles)) unname(panel_subtitles[name]) else NULL
  plot <- plot + labs(title = unname(panel_titles[name]), subtitle = subtitle)
  height <- height + .5
  panel_plots[[name]] <<- plot
  export_figure(plot, name, width, height, output)
  manifest[[length(manifest) + 1L]] <<- tibble(
    panel = name, width_cm = width, height_cm = height,
    square_plot_box = square
  )
  if (square) {
    tmp <- tempfile(fileext = ".pdf")
    pdf(tmp, width = width / 2.54, height = height / 2.54, family = font_family)
    print(plot)
    grid::grid.force()
    paths <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)$name
    viewport <- grep("^panel[.][0-9]", paths, value = TRUE)[1]
    if (is.na(viewport)) stop("Cannot verify square plot box: ", name)
    grid::seekViewport(viewport)
    w <- grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
    h <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
    stopifnot(is.finite(w), is.finite(h), abs(w - h) < .05)
    square_checks[[length(square_checks) + 1L]] <<- tibble(panel = name, plot_width_mm = w, plot_height_mm = h)
    dev.off()
    unlink(tmp)
  }
}

# B1, standalone end-to-end scatter. Only the 30 molecular-overlap recipients.
lap_r <- cor(lap_b$observed_log2, lap_b$predicted_log2)
stopifnot(abs(lap_r - .5864346872851314) < 1e-4)
p <- ggplot(lap_b, aes(observed_log2, predicted_log2)) +
  geom_point(colour = blue, size = 1.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = red, linewidth = .45) +
  stat_cor(corner = "bottom_right") +
  scale_x_continuous(breaks = c(-2, -1, 0)) +
  labs(
    x = expression("Measured viability (log"[2] * ")"),
    y = expression("Predicted viability (log"[2] * ")")
  ) +
  paper_theme(TRUE)
save_panel(p, "B1_2A_lapatinib_scatter", 4.4, 4.0, TRUE)

# B1, model and Wilson barplots are separate panels with their own endpoints.
hgf <- read_tsv("results/12_fit_ligand_rescue/data/AU565_HGF_predictions_frozen.tsv", show_col_types = FALSE) |>
  filter(CellLine == "AU565", dose == 5, source_profile == "Primary_study_mean")
stopifnot(nrow(hgf) == 1L)
model_bars <- tibble(
  condition = factor(c("Lapatinib", "Lapatinib\n+ HGF"),
    levels = c("Lapatinib", "Lapatinib\n+ HGF")
  ),
  value = c(hgf$lapatinib_pct_control, hgf$predicted_combination_pct_control), sd = 0
)
wilson_source <- read_tsv("data/raw/wilson/wilson_observed_source.tsv", show_col_types = FALSE)
wilson_bars <- wilson_source |>
  filter(series == "Wilson et al.", condition %in% c("Lapatinib", "Lapatinib + HGF")) |>
  mutate(condition = factor(if_else(condition == "Lapatinib", "Lapatinib", "Lapatinib\n+ HGF"),
    levels = levels(model_bars$condition)
  )) |>
  rename(value = abundance)
stopifnot(nrow(wilson_bars) == 2L)
make_abundance_bars <- function(data, title, subtitle, ylabel, errors = FALSE) {
  p <- ggplot(data, aes(condition, value, fill = condition)) +
    geom_col(width = .57, colour = "black", linewidth = .25) +
    geom_text(aes(y = value + sd + 3, label = sprintf("%.1f", value)), size = 2.6, vjust = 0) +
    scale_fill_manual(values = c(blue, orange), guide = "none") +
    scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = c("Lapatinib" = "Alone", "Lapatinib\n+ HGF" = "+ HGF")) +
    labs(title = title, subtitle = subtitle, x = "Lapatinib", y = ylabel) +
    paper_theme()
  if (errors) p <- p + geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .13, linewidth = .35)
  p
}
save_panel(
  make_abundance_bars(model_bars, "AU565, model prediction", "Lapatinib (5 µM)", "Relative cell abundance"),
  "B1_2B_AU565_prediction_bars", 3.8, 3.6
)
save_panel(
  make_abundance_bars(wilson_bars, "AU565, Wilson et al.", "Lapatinib (1 µM, 12 days)", "Relative cell abundance", TRUE),
  "B1_2C_Wilson_bars", 3.8, 3.6
)
write_tsv(model_bars, file.path(audit, "B1_model_bar_values.tsv"))
write_tsv(wilson_bars, file.path(audit, "B1_Wilson_bar_values.tsv"))

# B2 figure 1, single-value summaries use bars anchored at zero, not isolated dots.
order <- selected$program
bar_labels <- setNames(selected$plot_label, selected$program)
selected <- selected |> mutate(program = factor(program, levels = rev(order)))
make_program_bars <- function(column, title, xlabel, colour, signed = FALSE, row_labels = bar_labels) {
  d <- selected |> mutate(
    value = .data[[column]], label_x = value + if_else(value < 0, -.04, .04),
    label_hjust = if_else(value < 0, 1, 0)
  )
  ggplot(d, aes(value, program)) +
    geom_col(fill = colour, width = .65) +
    geom_vline(xintercept = 0, linewidth = .25, colour = "black") +
    geom_text(aes(x = label_x, label = sprintf("%.2f", value), hjust = label_hjust), size = 2.5) +
    scale_y_discrete(labels = row_labels) +
    scale_x_continuous(
      limits = if (signed) c(-1.1, 1.05) else c(0, .86),
      breaks = if (signed) c(-.6, 0, .6) else c(0, .3, .6), expand = expansion(mult = 0)
    ) +
    labs(
      title = title, subtitle = if (signed) "Independent viability measurements (n = 30)" else "Recipient-excluded molecular predictions (n = 42)",
      x = xlabel, y = NULL
    ) +
    paper_theme() +
    theme(axis.text.y = element_text(size = 7.5))
}
save_panel(
  make_program_bars("prediction_r", "Predictable response programs", "Pearson r", blue),
  "B2_1A_program_predictability_bars", 7.3, 3.4
)
save_panel(make_program_bars(
  "viability_r", "Program association with viability", "Pearson r", teal, TRUE,
  setNames(order, order)
), "B2_1B_program_viability_bars", 3.8, 3.4)

# B2 figure 2A, measured-input decoder performance across all eight drug contexts.
shared <- measured |> inner_join(baseline,
  by = c("depmap_id", "drug", "dose"),
  suffix = c("_measured", "_predicted"), relationship = "one-to-one"
)
short <- c(
  "5-Azacytidine" = "Azacytidine", "Bendamustine" = "Bendamustine", "Bortezomib" = "Bortezomib",
  "Busulfan" = "Busulfan", "Cytarabine (hydrochloride)" = "Cytarabine",
  "Lapatinib ditosylate" = "Lapatinib", "Topotecan (hydrochloride)" = "Topotecan", "Vandetanib" = "Vandetanib"
)
drug_stats <- shared |>
  group_by(drug, dose) |>
  summarise(n = n(), r = cor(observed_log2_measured, predicted_log2_measured), .groups = "drop") |>
  mutate(name = paste0(unname(short[drug]), " (", sprintf("%g", dose), " µM)")) |>
  arrange(r)
stopifnot(sum(drug_stats$r > 0) == 7L, sum(drug_stats$r > .3) == 5L)
drug_stats <- drug_stats |> mutate(
  name = factor(name, levels = name), label_x = r + if_else(r < 0, -.025, .025),
  label_hjust = if_else(r < 0, 1, 0)
)
p <- ggplot(drug_stats, aes(r, name)) +
  geom_col(fill = blue, width = .65) +
  geom_vline(xintercept = 0, linewidth = .25) +
  geom_text(aes(x = label_x, label = sprintf("%.2f", r), hjust = label_hjust), size = 2.5) +
  scale_x_continuous(limits = c(-.28, .72), breaks = c(0, .3, .6), expand = expansion(mult = 0)) +
  labs(
    title = "Viability prediction\nacross eight drugs", subtitle = "Measured response inputs, excluded recipients",
    x = "Pearson r", y = NULL
  ) +
  paper_theme() +
  theme(axis.text.y = element_text(size = 7.5))
save_panel(p, "B2_2A_shared_decoder_bars", 6.8, 3.6)
write_tsv(drug_stats, file.path(audit, "shared_decoder_drug_correlations.tsv"))

# B2 figure 2B: the seven programs selected by single-drug sensitivity association.
# Selection is independent of the direction or correlation of the HGF response.
opposition <- read_tsv("results/14_plot_model_diagnostics/AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv", show_col_types = FALSE) |>
  mutate(composition_program = abs(lapatinib_sensitivity_r) >= .5)
stopifnot(
  nrow(opposition) == 40L, sum(opposition$composition_program) == 7L,
  sum(with(filter(opposition, composition_program), lapatinib_displacement * HGF_displacement < 0)) == 6L
)
sensitive_opposition <- opposition |> filter(composition_program)
opposition_test <- cor.test(sensitive_opposition$lapatinib_displacement, sensitive_opposition$HGF_displacement)
write_tsv(
  tibble(n = nrow(sensitive_opposition), r = unname(opposition_test$estimate), p = opposition_test$p.value),
  file.path(audit, "sensitivity_program_opposition_correlation.tsv")
)
p <- ggplot(sensitive_opposition, aes(lapatinib_displacement, HGF_displacement)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = red, linewidth = .4) +
  stat_cor(show_p = FALSE) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .25, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = .25, linetype = "dashed") +
  geom_point(colour = blue, size = 1.4) +
  geom_text_repel(aes(label = program),
    seed = 42, size = 2.5, colour = "black", min.segment.length = 0,
    nudge_y = ifelse(sensitive_opposition$program == "P29", -.012, 0),
    segment.colour = "grey55", segment.size = .2, box.padding = .3,
    point.padding = .2, max.overlaps = Inf
  ) +
  scale_x_continuous(
    breaks = c(-.5, 0, .25), labels = c("-0.5", "0", "0.25"),
    expand = expansion(mult = c(.08, .18))
  ) +
  scale_y_continuous(expand = expansion(mult = c(.07, .40))) +
  labs(
    title = "Lapatinib and HGF response programs", subtitle = "Opposite directions in 6 of 7 composition programs",
    x = "Lapatinib displacement", y = "HGF displacement"
  ) +
  paper_theme(TRUE) +
  theme(plot.subtitle = element_text(size = 7.5))
save_panel(p, "B2_2B_lapatinib_HGF_program_scatter", 4.4, 4.0, TRUE)
write_tsv(opposition, file.path(audit, "drug_HGF_program_displacements.tsv"))

# B2 figure 3A, measured and predicted inputs on identical recipients.
chain <- tibble(
  input = factor(c("Measured\nresponses", "Predicted\nresponses"),
    levels = c("Measured\nresponses", "Predicted\nresponses")
  ),
  r = c(cor(lap_m$observed_log2, lap_m$predicted_log2), lap_r)
)
make_chain_scatter <- function(data, colour) {
  ggplot(data, aes(observed_log2, predicted_log2)) +
    geom_point(colour = colour, size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = red, linewidth = .4) +
    stat_cor(corner = "bottom_right") +
    scale_x_continuous(breaks = c(-2, -1, 0)) +
    scale_y_continuous(expand = expansion(mult = c(.2, .15))) +
    labs(x = "Measured viability\n(log2)", y = "Predicted viability\n(log2)") +
    paper_theme(TRUE)
}
save_panel(make_chain_scatter(lap_m, teal), "B2_3A_measured_input_scatter", 4.4, 4.0, TRUE)
save_panel(make_chain_scatter(lap_b, blue), "B2_3A_predicted_input_scatter", 4.4, 4.0, TRUE)
write_tsv(lap_m, file.path(audit, "B2_3A_measured_input_points.tsv"))
write_tsv(lap_b, file.path(audit, "B2_3A_predicted_input_points.tsv"))
write_tsv(chain, file.path(audit, "prediction_chain_correlations.tsv"))

# B2 figure 3B, square conditional interaction scatter. Equal limits and equal scales.
combo <- read_tsv("results/12_fit_ligand_rescue/data/combination_predictions_frozen.tsv",
  show_col_types = FALSE) |> filter(evaluation == "held_cell", gene == "NRG1")
dose_stats <- combo |>
  group_by(dose) |>
  summarise(n = n(), r = cor(observed_departure_log2, predicted_departure_log2), .groups = "drop")
stopifnot(all(dose_stats$n == 6L))
dose_labels <- setNames(sprintf("%g µM, r = %.2f", dose_stats$dose, dose_stats$r), dose_stats$dose)
combo <- combo |> mutate(dose_label = factor(as.character(dose), levels = c("0.5", "5")))
p <- ggplot(combo, aes(observed_departure_log2, predicted_departure_log2, colour = dose_label, shape = dose_label)) +
  geom_abline(intercept = 0, slope = 1, linewidth = .3, colour = "grey65", linetype = "dashed") +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c("0.5" = orange, "5" = blue), labels = dose_labels) +
  scale_shape_manual(values = c("0.5" = 17, "5" = 16), labels = dose_labels) +
  scale_x_continuous(breaks = c(-1, 0, 1, 2)) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2)) +
  coord_fixed(ratio = 1, xlim = c(-1.5, 2.05), ylim = c(-1.5, 2.05), expand = FALSE) +
  annotate("text", x = 1.9, y = -1.35, label = "n = 6 per dose", hjust = 1, size = 2.5) +
  labs(
    title = "Lapatinib and NRG1", subtitle = "Measured single-agent controls supplied",
    x = expression("Measured interaction (log"[2] * ")"),
    y = expression("Predicted interaction (log"[2] * ")")
  ) +
  paper_theme(TRUE) +
  theme(
    legend.position = c(.02, .98), legend.justification = c(0, 1),
    legend.key.size = grid::unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = grid::unit(0, "mm"), legend.text = element_text(size = 7.1)
  )
save_panel(p, "B2_3B_NRG1_interaction_scatter", 4.4, 4.0, TRUE)
write_tsv(dose_stats, file.path(audit, "NRG1_interaction_correlations.tsv"))

# B2 molecular prediction: separate held-recipient scatters for the seven joint-criterion programs.
# Each point is one of the same 42 recipients used in the molecular-predictability summary.
for (program in as.character(selected$program)) {
  displacement <- tibble(depmap_id = rna$depmap_id, measured = y[, program], predicted = prediction[, program])
  stopifnot(nrow(displacement) == 42L, all(is.finite(as.matrix(displacement[, c("measured", "predicted")]))))
  p <- ggplot(displacement, aes(measured, predicted)) +
    geom_point(colour = blue, size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = red, linewidth = .4) +
    stat_cor() +
    scale_y_continuous(expand = expansion(mult = c(.05, .30))) +
    labs(x = "Measured program change", y = "Predicted program change") +
    paper_theme(TRUE)
  save_panel(p, paste0("B2_1C_", program, "_displacement_scatter"), 4.4, 4.0, TRUE)
  write_tsv(displacement, file.path(audit, paste0(program, "_displacement_scatter_points.tsv")))
}

# Composite figures are assembled from the same ggplot objects as the individual exports.
composite_output <- file.path(root, "composites")
dir.create(composite_output, recursive = TRUE, showWarnings = FALSE)
composite_manifest <- list()
composite_square_checks <- list()
save_composite <- function(names, name, widths, body_height = 4.3, gap = .05) {
  stopifnot(length(names) == length(widths), all(names %in% names(panel_plots)))
  width <- sum(widths) + gap * (length(names) - 1)
  height <- body_height
  offsets <- c(0, cumsum(widths + gap))[seq_along(names)]
  grobs <- list()
  plot_grobs <- lapply(seq_along(names), function(i) {
    plot <- panel_plots[[names[i]]] + labs(tag = LETTERS[i]) +
      theme(
        plot.tag = element_text(family = font_family, face = "bold", size = 9, hjust = 0, vjust = 1),
        plot.tag.position = "topleft"
      )
    ggplotGrob(plot)
  })
  if (name == "B2_Figure_1") {
    common_heights <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$heights))
    for (i in seq_along(plot_grobs)) plot_grobs[[i]]$heights <- common_heights
  }
  if (name == "B2_Figure_3") {
    common_widths <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$widths))
    for (i in seq_along(plot_grobs)) plot_grobs[[i]]$widths <- common_widths
  }
  for (i in seq_along(names)) {
    vp <- grid::viewport(
      x = grid::unit(offsets[i], "cm"), y = grid::unit(0, "cm"),
      width = grid::unit(widths[i], "cm"), height = grid::unit(body_height, "cm"),
      just = c("left", "bottom"), name = paste0(name, "_", LETTERS[i])
    )
    grobs[[length(grobs) + 1L]] <- grid::grobTree(plot_grobs[[i]], vp = vp)
  }
  figure <- do.call(grid::grobTree, grobs)
  export_figure(figure, name, width, height, composite_output)
  composite_manifest[[length(composite_manifest) + 1L]] <<-
    tibble(figure = name, width_cm = width, height_cm = height, panels = paste(names, collapse = "; "))
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp, width = width / 2.54, height = height / 2.54, family = font_family)
  grid::grid.newpage()
  grid::grid.draw(figure)
  grid::grid.force()
  for (i in seq_along(names)) {
    if (!isTRUE(panel_plots[[names[i]]]$theme$aspect.ratio == 1)) next
    grid::seekViewport(paste0(name, "_", LETTERS[i]))
    viewports <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)$name
    target <- grep("^panel[.][0-9]", viewports, value = TRUE)[1]
    if (is.na(target)) stop("Cannot locate composite scatter viewport")
    grid::downViewport(target)
    w <- grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
    h <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
    if (abs(w - h) > .05) stop("Non-square composite scatter: ", names[i])
    composite_square_checks[[length(composite_square_checks) + 1L]] <<-
      tibble(figure = name, panel = LETTERS[i], width_mm = w, height_mm = h)
    grid::upViewport(0)
  }
  dev.off()
  unlink(tmp)
}
save_composite(c("B1_2A_lapatinib_scatter", "B1_2B_AU565_prediction_bars", "B1_2C_Wilson_bars"),
  "B1_Figure_2", c(4.4, 3.8, 3.8),
  body_height = 4.6
)
save_composite(
  c("B2_1A_program_predictability_bars", "B2_1B_program_viability_bars"),
  "B2_Figure_1", c(7.3, 3.8)
)
save_composite(
  c("B2_2A_shared_decoder_bars", "B2_2B_lapatinib_HGF_program_scatter"),
  "B2_Figure_2", c(6.8, 4.4)
)
save_composite(
  c("B2_3A_measured_input_scatter", "B2_3A_predicted_input_scatter"),
  "B2_Figure_3", c(4.4, 4.4)
)
write_tsv(bind_rows(composite_manifest), file.path(root, "composite_manifest.tsv"))
write_tsv(bind_rows(composite_square_checks), file.path(audit, "composite_square_checks.tsv"))

write_tsv(bind_rows(manifest), file.path(root, "panel_manifest.tsv"))
write_tsv(bind_rows(square_checks), file.path(audit, "square_plot_box_checks.tsv"))
writeLines(capture.output(sessionInfo()), file.path(audit, "plot_session_info.txt"))
sink()
message("Saved ", length(manifest), " independent panels (PDF, SVG, PNG).")
message("All ", length(square_checks), " scatter plotting boxes were checked to be square.")
