#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
source(file.path(WORK_ROOT, "src", "R", "lib_heatmap.R"))
RESULT <- stage_result("03_baseline_programs")
ARTIFACT <- stage_artifact("03_baseline_programs")

theme_pilot <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.2) theme_pilot(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(filename, plot, width, height) ggsave(
  file.path(RESULT, filename), plot, width = width, height = height,
  device = cairo_pdf, bg = "white", limitsize = FALSE)
save_heatmap_pdf <- function(filename, heatmap, width, height) {
  grDevices::pdf(file.path(RESULT, filename), width = width, height = height,
    useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap)
  grDevices::dev.off()
}

plot_stage03 <- function() {
  baseline_scores <- fread(file.path(ARTIFACT,
    "baseline_program_scores.tsv.gz"))
  baseline_stability <- fread(file.path(RESULT,
    "baseline_program_stability.tsv"))
  baseline_diagnostics <- fread(file.path(RESULT,
    "baseline_program_diagnostics.tsv"))
  baseline_weights <- fread(file.path(RESULT,
    "baseline_program_weights.tsv.gz"))
  baseline_assignment <- fread(file.path(RESULT,
    "baseline_gene_program_assignment.tsv.gz"))

  ## Retained plotting blocks below are transferred from the accepted
  ## single-fold pipeline; paths and the native NPY loader are the only
  ## mechanical adaptations.
  stability_once <- baseline_stability[, .SD[1L], by = program]
  p_b01 <- ggplot(baseline_stability,
    aes(program, matched_pearson, color = stable)) +
    geom_hline(yintercept = 0.50, linetype = 2, linewidth = 0.4,
      color = "grey45") +
    geom_jitter(width = 0.10, height = 0, alpha = 0.70, size = 1.6) +
    geom_point(data = stability_once, aes(y = q10_matched_pearson),
      shape = 23, fill = "gold", color = "black", size = 3) +
    scale_color_manual(values = c(`TRUE` = "#2E7D32",
      `FALSE` = "#B2182B"), guide = "none") +
    labs(x = NULL, y = "Matched repeat membership Pearson r",
      title = "NTC baseline-program stability",
      subtitle = "Points: ten oNMF repeats; diamonds: q10; dashed: 0.50 threshold") +
    theme_pilot()
  save_gg("baseline_01_program_stability.pdf", p_b01, 5.8, 4.0)

  baseline_line <- baseline_scores[baseline_source == "tahoe_dmso" &
    aggregation_level == "line_mean",
    .(depmap_id, program, role,
      line_score = standardized_program_score)]
  baseline_plate <- merge(baseline_scores[baseline_source == "tahoe_dmso" &
    aggregation_level == "natural_plate"], baseline_line,
    by = c("depmap_id", "program", "role"), all.x = TRUE)
  p_b02 <- ggplot(baseline_plate, aes(line_score,
      standardized_program_score, color = role)) +
    geom_point(alpha = 0.42, size = 0.65) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 0.5) +
    facet_wrap(~program, scales = "free", ncol = 3) +
    scale_color_manual(values = c(training = "#2166AC",
      held_out = "#B2182B")) +
    labs(x = "Complete line-mean score", y = "Natural-plate score",
      title = "Natural DMSO plate reproducibility",
      subtitle = "Discovery and ICC use training lines; held lines are shown only as diagnostics") +
    square_scatter_theme()
  save_gg("baseline_02_plate_reproducibility.pdf", p_b02, 8.4,
    max(5.6, 2.6 * ceiling(uniqueN(baseline_plate$program) / 3)))

  baseline_correlation <- read_npy_float32(file.path(ARTIFACT,
    "baseline_gene_gene_correlation.npy"))
  assignment_original <- copy(baseline_assignment)
  setorder(assignment_original, selected_gene_index)
  if (nrow(baseline_correlation) != nrow(assignment_original))
    fail("Baseline correlation and assignment dimensions disagree.")
  rownames(baseline_correlation) <- assignment_original$gene
  colnames(baseline_correlation) <- assignment_original$gene
  block_index <- order(assignment_original$plotting_order)
  block_assignment <- assignment_original[block_index]
  baseline_correlation_block <- baseline_correlation[block_index,
    block_index, drop = FALSE]
  stopifnot(identical(rownames(baseline_correlation_block),
    block_assignment$gene), identical(colnames(baseline_correlation_block),
    block_assignment$gene))
  block_annotation <- data.frame(Program = block_assignment$assigned_program,
    Direction = block_assignment$direction)
  hm_b03a <- draw_heatmap(baseline_correlation_block,
    col_anno_df = block_annotation, col_anno_cols = "Program",
    row_anno_df = block_annotation, row_anno_cols = "Program",
    scale_rows = FALSE, value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = FALSE, show_column_names = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Training-line baseline gene correlation: program-block order")
  save_heatmap_pdf("baseline_03a_gene_gene_correlation_full.pdf", hm_b03a,
    7.2, 7.0)

  top_labeled <- baseline_assignment[order(assigned_program, -absolute_weight),
    head(.SD, 12L), by = assigned_program]
  top_labeled <- top_labeled[order(plotting_order)]
  labeled_genes <- top_labeled$gene
  hm_b03b <- draw_heatmap(baseline_correlation[labeled_genes,
    labeled_genes, drop = FALSE],
    col_anno_df = data.frame(Program = top_labeled$assigned_program),
    col_anno_cols = "Program",
    row_anno_df = data.frame(Program = top_labeled$assigned_program),
    row_anno_cols = "Program", scale_rows = FALSE,
    value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 5),
    column_names_gp = grid::gpar(fontsize = 5),
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Top-weight genes in frozen program-block order")
  save_heatmap_pdf("baseline_03b_gene_gene_correlation_labeled.pdf", hm_b03b,
    8.6, 8.0)

  retained_programs <- baseline_diagnostics$program
  weight_wide <- dcast(baseline_weights[program %chin% retained_programs],
    gene ~ program, value.var = "signed_weight")
  weight_matrix <- as.matrix(weight_wide[, ..retained_programs])
  rownames(weight_matrix) <- weight_wide$gene
  top_weight_rows <- baseline_assignment[
    assigned_program %chin% retained_programs][order(assigned_program,
      -absolute_weight), head(.SD, 30L), by = assigned_program]
  top_weight_rows <- top_weight_rows[order(assigned_program, -absolute_weight)]
  weight_matrix <- weight_matrix[top_weight_rows$gene, retained_programs,
    drop = FALSE]
  stopifnot(identical(rownames(weight_matrix), top_weight_rows$gene))
  hm_b04 <- draw_heatmap(weight_matrix,
    row_anno_df = data.frame(Program = top_weight_rows$assigned_program,
      Direction = top_weight_rows$direction),
    row_anno_cols = c("Program", "Direction"), scale_rows = FALSE,
    value_range = c(-max(abs(weight_matrix)), 0, max(abs(weight_matrix))),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 5),
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Frozen signed NTC baseline-program weights")
  save_heatmap_pdf("baseline_04_program_gene_weight_heatmap.pdf", hm_b04,
    7.4, 10.0)

  control_meta <- fread(file.path(CANONICAL_CACHE,
    "dmso_control_metadata.tsv"))
  control_genes <- fread(file.path(CANONICAL_CACHE,
    "dmso_control_genes.tsv"))
  gene_columns <- match(assignment_original$gene, control_genes$feature)
  if (anyNA(gene_columns)) fail("A selected baseline gene is absent from DMSO.")
  control_selected <- read_npy_float32_columns(file.path(CANONICAL_CACHE,
    "dmso_control_expression.npy"), gene_columns)
  retained_control <- control_meta$depmap_id %chin% unique(baseline_line$depmap_id)
  control_selected <- control_selected[control_meta$array_row[
    retained_control] + 1L, , drop = FALSE]
  control_meta <- control_meta[retained_control]
  control_log <- log2(1 + sweep(control_selected, 1L,
    rowSums(control_selected), "/") * 10000)
  line_ids <- sort(unique(control_meta$depmap_id))
  line_expression <- t(vapply(line_ids, function(id)
    colMeans(control_log[control_meta$depmap_id == id, , drop = FALSE]),
    numeric(ncol(control_log))))
  rownames(line_expression) <- line_ids
  line_roles <- unique(baseline_line[, .(depmap_id, role)])
  training_ids <- line_roles[role == "training", depmap_id]
  gene_center <- colMeans(line_expression[training_ids, , drop = FALSE])
  gene_scale <- apply(line_expression[training_ids, , drop = FALSE], 2L, sd)
  gene_scale[!is.finite(gene_scale) | gene_scale <= 0] <- 1
  expression_z <- sweep(sweep(line_expression, 2L, gene_center, "-"), 2L,
    gene_scale, "/")
  expression_heatmap <- t(expression_z[training_ids, , drop = FALSE])
  rownames(expression_heatmap) <- assignment_original$gene
  expression_heatmap <- expression_heatmap[block_index, , drop = FALSE]
  line_order <- hclust(dist(t(expression_heatmap)), method = "average")$order
  expression_heatmap <- expression_heatmap[, line_order, drop = FALSE]
  stopifnot(identical(rownames(expression_heatmap), block_assignment$gene))
  hm_b05 <- draw_heatmap(expression_heatmap,
    row_anno_df = block_annotation, row_anno_cols = "Program",
    scale_rows = FALSE, value_range = c(-2.5, 0, 2.5),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = FALSE, show_column_names = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "26 training-line DMSO means: gene z-scores by program block")
  save_heatmap_pdf("baseline_05_ntc_expression_by_program.pdf", hm_b05,
    8.0, 8.6)

  line_order_score <- baseline_line[program == retained_programs[1L]][
    order(role, line_score), depmap_id]
  baseline_plate[, depmap_id := factor(depmap_id, levels = line_order_score)]
  p_b06 <- ggplot(baseline_plate,
    aes(depmap_id, standardized_program_score, color = role)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_point(alpha = 0.42, size = 0.52,
      position = position_jitter(width = 0.12, height = 0)) +
    geom_point(data = baseline_line, aes(depmap_id, line_score),
      inherit.aes = FALSE, shape = 21, fill = "white", color = "black",
      size = 1.25) +
    facet_wrap(~program, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c(training = "#2166AC",
      held_out = "#B2182B")) +
    labs(x = "Cell line (ordered by B1 line mean)",
      y = "Standardized baseline-program score",
      title = "Natural plate scores grouped by cell line",
      subtitle = "Open circles: complete line means") +
    theme_pilot(7.5) + theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  save_gg("baseline_06_program_scores_by_line_and_plate.pdf", p_b06,
    8.6, max(5.7, 2.6 * ceiling(length(retained_programs) / 2)))

  line_score_wide <- dcast(baseline_line, depmap_id + role ~ program,
    value.var = "line_score")
  redundancy <- cor(as.matrix(line_score_wide[, ..retained_programs]),
    method = "pearson")
  hm_b07 <- draw_heatmap(redundancy[retained_programs, retained_programs,
    drop = FALSE], scale_rows = FALSE, value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Complete-line NTC baseline-program correlations")
  save_heatmap_pdf("baseline_07_program_redundancy_heatmap.pdf", hm_b07,
    5.2, 4.8)

  baseline_line_all <- baseline_scores[aggregation_level %chin%
    c("line_mean", "line_level")]
  baseline_concordance <- dcast(baseline_line_all,
    depmap_id + cell_line + program + role ~ baseline_source,
    value.var = "standardized_program_score")
  p_b09 <- ggplot(baseline_concordance,
    aes(tahoe_dmso, depmap_static, color = role)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
    geom_point(alpha = 0.75, size = 1.0) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 0.5) +
    facet_wrap(~program, scales = "free", ncol = 3) +
    scale_color_manual(values = c(training = "#2166AC",
      held_out = "#B2182B")) +
    labs(x = "DMSO baseline-program score",
      y = "DepMap frozen-program score",
      title = "DMSO–DepMap baseline-program concordance",
      subtitle = "DepMap uses frozen DMSO oNMF weights with independent training-line normalization") +
    square_scatter_theme()
  save_gg("baseline_09_dmso_depmap_concordance.pdf", p_b09, 8.4,
    max(5.6, 2.6 * ceiling(length(retained_programs) / 3)))

  concentration <- fread(file.path(RESULT,
    "baseline_program_line_concentration.tsv"))
  concentration[, source_label := fifelse(baseline_source == "tahoe_dmso",
    "Tahoe DMSO", "DepMap static")]
  p_concentration <- ggplot(concentration,
    aes(program, top_1_squared_score_fraction, fill = source_label)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.66) +
    geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35,
      color = "grey45") +
    scale_fill_manual(values = c("Tahoe DMSO" = "#2166AC",
      "DepMap static" = "#B2182B")) +
    labs(x = NULL, y = "Top-line fraction of squared score",
      fill = "Baseline source",
      title = "Baseline-program individual-line concentration") +
    theme_pilot()
  save_gg("baseline_10_program_line_concentration.pdf", p_concentration,
    5.8, 4.2)
}

plots_only <- "--plots-only" %in% commandArgs(trailingOnly = TRUE)
if (!plots_only) {
  status <- system2("/usr/bin/python3",
    file.path(root, "scripts", "03_run_baseline_onmf.py"),
    env = sprintf("TAHOE_V2_ROOT=%s", root))
  if (!identical(status, 0L))
    stop(sprintf("Native baseline-oNMF stage failed with status %d.", status),
      call. = FALSE)
}
plot_stage03()
message(if (plots_only) "Stage 03 plots refreshed from existing artifacts."
  else "Stage 03 complete, including baseline-program diagnostics.")
