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
RESULT <- stage_result("02_response_programs")
ARTIFACT <- stage_artifact("02_response_programs")

theme_pilot <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
save_gg <- function(filename, plot, width, height) ggsave(
  file.path(RESULT, filename), plot, width = width, height = height,
  device = cairo_pdf, bg = "white", limitsize = FALSE)
save_heatmap_pdf <- function(filename, heatmap, width, height) {
  grDevices::pdf(file.path(RESULT, filename), width = width, height = height,
    useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap)
  grDevices::dev.off()
}
zscore_columns <- function(matrix) {
  center <- colMeans(matrix)
  scale <- apply(matrix, 2L, sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  sweep(sweep(matrix, 2L, center, "-"), 2L, scale, "/")
}

plot_stage02 <- function() {
  response_stability <- fread(file.path(RESULT,
    "fold_program_stability.tsv"))
  response_weights <- fread(file.path(RESULT,
    "fold_program_weights.tsv.gz"))
  response_scores <- fread(file.path(ARTIFACT,
    "fold_program_scores.tsv.gz"))
  response_assignment <- fread(file.path(RESULT,
    "fold_gene_program_assignment.tsv.gz"))[stable == TRUE]
  response_bottleneck <- fread(file.path(RESULT,
    "program_bottleneck_audit.tsv"))
  response_observed <- fread(file.path(ARTIFACT,
    "observed_response_metadata.tsv.gz"))
  response_signal <- fread(file.path(ARTIFACT,
    "gene_signal_selection.tsv.gz"))[selected == TRUE][order(selection_rank)]
  response_programs <- response_stability[stable == TRUE, unique(program)]
  response_programs <- paste0("P", sort(as.integer(sub("P", "",
    response_programs))))
  if (length(response_programs) != 9L || "P6" %chin% response_programs)
    fail(paste("The corrected response basis must contain the nine stable",
      "programs and exclude P6."))

  ## Retained plotting blocks below are transferred from the accepted
  ## single-fold pipeline; paths and the native NPY loader are the only
  ## mechanical adaptations.
  p_r01 <- ggplot(response_stability,
    aes(program, matched_pearson, color = stable)) +
    geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.4) +
    geom_jitter(width = 0.10, height = 0, alpha = 0.45, size = 0.8) +
    geom_point(data = response_stability[, .SD[1L], by = program],
      aes(y = q10_matched_pearson), shape = 23, fill = "gold",
      color = "black", size = 2.5) +
    scale_color_manual(values = c(`TRUE` = "#2166AC",
      `FALSE` = "grey70"), guide = "none") +
    labs(x = NULL, y = "Matched repeat membership Pearson r",
      title = "Response-program oNMF stability",
      subtitle = "Ten fixed-seed repeats; diamonds show q10 stability") +
    theme_pilot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_gg("01_program_stability.pdf", p_r01, 6.0, 4.2)

  response_correlation <- read_npy_float32(file.path(ARTIFACT,
    "reference_gene_gene_cross_correlation.npy"))
  setorder(response_assignment, plotting_order)
  response_index <- response_assignment$selected_gene_index + 1L
  response_correlation <- response_correlation[response_index,
    response_index, drop = FALSE]
  rownames(response_correlation) <- response_assignment$gene
  colnames(response_correlation) <- response_assignment$gene
  response_split <- factor(response_assignment$assigned_program,
    levels = response_programs)
  response_annotation <- data.frame(Program = response_split)
  hm_r03a <- draw_heatmap(response_correlation,
    col_anno_df = response_annotation, col_anno_cols = "Program",
    row_anno_df = response_annotation, row_anno_cols = "Program",
    scale_rows = FALSE, value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    title = "Signed response gene correlation: frozen program-block order")
  save_heatmap_pdf("03a_gene_gene_cross_correlation_full.pdf", hm_r03a,
    7.4, 7.2)

  response_labeled <- response_assignment[order(assigned_program,
    -absolute_weight), head(.SD, 40L), by = assigned_program][
      order(plotting_order)]
  response_labeled_split <- factor(response_labeled$assigned_program,
    levels = response_programs)
  hm_r03b <- draw_heatmap(response_correlation[response_labeled$gene,
    response_labeled$gene, drop = FALSE],
    col_anno_df = data.frame(Program = response_labeled_split),
    col_anno_cols = "Program",
    row_anno_df = data.frame(Program = response_labeled_split),
    row_anno_cols = "Program", scale_rows = FALSE,
    value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 4),
    column_names_gp = grid::gpar(fontsize = 4),
    title = "Top response genes in frozen program-block order")
  save_heatmap_pdf("03b_gene_gene_cross_correlation_labeled.pdf", hm_r03b,
    18.0, 17.0)

  response_signed_weights <- response_weights[
    program %chin% response_programs,
    .(signed_weight = sum(fifelse(direction == "UP", weight, -weight))),
    by = .(gene, program)]
  response_weight_wide <- dcast(response_signed_weights, gene ~ program,
    value.var = "signed_weight", fill = 0)
  response_weight_matrix <- as.matrix(
    response_weight_wide[, ..response_programs])
  rownames(response_weight_matrix) <- response_weight_wide$gene
  response_feature_manifest <- response_assignment[order(assigned_program,
    -absolute_weight), head(.SD, 40L), by = assigned_program]
  setorder(response_feature_manifest, assigned_program, -absolute_weight, gene)
  response_weight_display <- response_weight_matrix[
    response_feature_manifest$gene, response_programs, drop = FALSE]
  hm_r03 <- draw_heatmap(response_weight_display,
    row_anno_df = data.frame(
      Program = factor(response_feature_manifest$assigned_program,
        levels = response_programs),
      Direction = response_feature_manifest$direction),
    row_anno_cols = c("Program", "Direction"), scale_rows = FALSE,
    value_range = c(-max(abs(response_weight_display)), 0,
      max(abs(response_weight_display))),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    title = "Outcome-blind signed response-program weights")
  save_heatmap_pdf("03_program_gene_weight_heatmap.pdf", hm_r03,
    7.2, 11.2)

  selected_gene_position <- match(response_feature_manifest$gene,
    response_signal$gene)
  if (anyNA(selected_gene_position))
    fail("A displayed response gene is absent from the selected matrix.")
  response_feature_values <- read_npy_float32_columns(file.path(ARTIFACT,
    "normalized_selected_responses.npy"), selected_gene_position)
  score_wide <- dcast(response_scores[program %chin% response_programs],
    array_row + observation_id + depmap_id + cell_line + perturbation_id +
      drug + dose + dose_panel ~ program, value.var = "score")
  score_row <- match(score_wide$observation_id,
    response_observed$observation_id)
  if (anyNA(score_row)) fail("Response score rows do not align.")
  response_feature_values <- response_feature_values[score_row, , drop = FALSE]
  direction_sign <- ifelse(response_feature_manifest$direction == "UP", 1, -1)
  response_feature_values <- sweep(response_feature_values, 2L,
    direction_sign, "*")
  score_matrix <- as.matrix(score_wide[, ..response_programs])
  score_z <- zscore_columns(score_matrix)
  perturbation_levels_response <- unique(score_wide$perturbation_id)
  perturbation_rows_response <- split(seq_len(nrow(score_wide)), factor(
    score_wide$perturbation_id, levels = perturbation_levels_response))
  program_means <- t(vapply(perturbation_rows_response, function(index)
    colMeans(score_matrix[index, , drop = FALSE]),
    numeric(length(response_programs))))
  feature_means <- t(vapply(perturbation_rows_response, function(index)
    colMeans(response_feature_values[index, , drop = FALSE]),
    numeric(ncol(response_feature_values))))
  rownames(program_means) <- perturbation_levels_response
  rownames(feature_means) <- perturbation_levels_response
  colnames(feature_means) <- response_feature_manifest$gene
  program_means_z <- zscore_columns(program_means)
  feature_means_z <- zscore_columns(feature_means)
  absolute_program <- abs(program_means_z)
  dominant_program <- max.col(absolute_program, ties.method = "first")
  dominant_strength <- absolute_program[cbind(seq_len(nrow(absolute_program)),
    dominant_program)]
  second_strength <- apply(absolute_program, 1L, function(value)
    sort(value, decreasing = TRUE)[2L])
  associated <- dominant_strength >= 1 &
    dominant_strength / pmax(second_strength, .Machine$double.eps) >= 1.10
  for (program_index in seq_along(response_programs)) {
    if (!any(associated & dominant_program == program_index)) {
      candidate <- which.max(absolute_program[, program_index])
      associated[candidate] <- TRUE
      dominant_program[candidate] <- program_index
    }
  }
  associated_index <- which(associated)
  block_order <- unlist(lapply(seq_along(response_programs), function(k) {
    current <- associated_index[dominant_program[associated_index] == k]
    if (length(current) <= 1L) return(current)
    local <- hclust(as.dist(1 - cor(t(feature_means_z[current, ,
      drop = FALSE]))), method = "average")$order
    current[local]
  }), use.names = FALSE)
  program_block <- factor(response_programs[dominant_program[block_order]],
    levels = response_programs)
  feature_block <- factor(response_feature_manifest$assigned_program,
    levels = response_programs)
  discovery_matrix <- t(feature_means_z[block_order, , drop = FALSE])
  rownames(discovery_matrix) <- response_feature_manifest$gene
  hm_r04 <- draw_heatmap(discovery_matrix,
    row_anno_df = data.frame(Program = feature_block),
    row_anno_cols = "Program",
    col_anno_df = data.frame(Program = program_block),
    col_anno_cols = "Program", scale_rows = FALSE,
    value_range = c(-2, 0, 2),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    title = "Response features by within-program clustered perturbation blocks")
  save_heatmap_pdf("04_program_discovery_input_heatmap.pdf", hm_r04,
    12.0, 10.8)

  program_block_matrix <- t(program_means_z[block_order, , drop = FALSE])
  hm_r05 <- draw_heatmap(program_block_matrix,
    col_anno_df = data.frame(Program = program_block),
    col_anno_cols = "Program", scale_rows = FALSE,
    value_range = c(-2, 0, 2),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = FALSE,
    title = "Program displacement in mandatory clustered program blocks")
  save_heatmap_pdf("05_program_displacement_heatmap.pdf", hm_r05,
    12.0, 5.8)

  feature_program_cor <- cor(feature_means_z[block_order, , drop = FALSE],
    program_means_z[block_order, , drop = FALSE])
  hm_rfc <- draw_heatmap(feature_program_cor,
    row_anno_df = data.frame(Program = feature_block),
    row_anno_cols = "Program", scale_rows = FALSE,
    value_range = c(-0.7, 0, 0.7),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    title = "Directional response-feature correlation with programs")
  save_heatmap_pdf("program_directional_feature_correlation.pdf", hm_rfc,
    7.2, 10.8)

  aligned_mean <- vapply(seq_along(response_programs), function(k) {
    current <- which(dominant_program[block_order] == k)
    direction <- sign(program_means_z[block_order[current], k])
    direction[direction == 0] <- 1
    colMeans(feature_means_z[block_order[current], , drop = FALSE] * direction)
  }, numeric(ncol(feature_means_z)))
  rownames(aligned_mean) <- response_feature_manifest$gene
  colnames(aligned_mean) <- response_programs
  hm_ram <- draw_heatmap(aligned_mean,
    row_anno_df = data.frame(Program = feature_block),
    row_anno_cols = "Program", scale_rows = FALSE,
    value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    title = "Direction-aligned mean response by program")
  save_heatmap_pdf("program_direction_aligned_mean_response.pdf", hm_ram,
    7.2, 10.8)

  perturbation_info_response <- unique(score_wide[, .(perturbation_id, drug,
    dose_panel)])[match(rownames(program_means_z), perturbation_id)]
  enrichment_order <- order(perturbation_info_response$dose_panel,
    tolower(perturbation_info_response$drug))
  enrichment_split <- factor(perturbation_info_response$dose_panel[
    enrichment_order], levels = 1:3,
    labels = c("0.05 uM", "0.5 uM", "5 uM"))
  hm_r06 <- draw_heatmap(program_means_z[enrichment_order, , drop = FALSE],
    row_anno_df = data.frame(Dose = enrichment_split),
    row_anno_cols = "Dose", scale_rows = FALSE,
    value_range = c(-2, 0, 2),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    title = "Drug-dose mean response-program displacement")
  save_heatmap_pdf("06_drug_program_enrichment_heatmap.pdf", hm_r06,
    7.2, 10.0)

  hm_r07 <- draw_heatmap(cor(score_matrix), scale_rows = FALSE,
    value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    title = "Response-program score redundancy")
  save_heatmap_pdf("07_program_redundancy_heatmap.pdf", hm_r07, 5.2, 4.8)

  p_rb <- ggplot(unique(response_bottleneck[, .(outer_fold, program,
    effective_drug_count, barcode_like)]),
    aes(program, effective_drug_count, fill = barcode_like)) +
    geom_col(width = 0.72) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#B2182B",
      `FALSE` = "#4D8C6F")) +
    labs(x = NULL, y = "Effective number of activating drugs",
      fill = "Barcode-like", title = "Response-program bottleneck audit") +
    theme_pilot()
  save_gg("program_bottleneck_audit.pdf", p_rb, 5.6, 4.8)
}

plots_only <- "--plots-only" %in% commandArgs(trailingOnly = TRUE)
if (!plots_only) {
  status <- system2("/usr/bin/python3",
    file.path(root, "scripts", "02_run_response_onmf.py"),
    env = sprintf("TAHOE_V2_ROOT=%s", root))
  if (!identical(status, 0L))
    stop(sprintf("Native response-oNMF stage failed with status %d.", status),
      call. = FALSE)
}
plot_stage02()
message(if (plots_only) "Stage 02 plots refreshed from existing artifacts."
  else "Stage 02 complete, including response-program diagnostics.")
