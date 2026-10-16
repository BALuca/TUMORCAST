#!/usr/bin/env Rscript

## Learn and diagnose rank-eight recurrent baseline states from untreated data.
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
RESULT <- stage_result("10_baseline_states")
ARTIFACT <- stage_artifact("10_baseline_states")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

theme_stage10 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage10(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(name, plot, width, height) ggsave(file.path(RESULT, name),
  plot, width = width, height = height, device = cairo_pdf, bg = "white",
  limitsize = FALSE)
save_heatmap <- function(name, heatmap, width, height) {
  grDevices::pdf(file.path(RESULT, name), width = width, height = height,
    useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap)
  grDevices::dev.off()
}

## Fail before factorization if any accepted untreated input is absent.
required <- c(file.path(CANONICAL_CACHE, "dmso_control_expression.npy"),
  file.path(CANONICAL_CACHE, "dmso_control_metadata.tsv"),
  file.path(CANONICAL_CACHE, "dmso_control_genes.tsv"),
  file.path(stage_artifact("01_inputs"), "prepared_inputs.rds"),
  file.path(stage_result("01_inputs"), "fold_assignments.tsv"),
  file.path(stage_artifact("03_baseline_programs"),
    "baseline_program_scores.tsv.gz"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 10 input(s): %s",
  paste(missing, collapse = ", "))

## Invoke only the existing native signed-oNMF machinery; R owns the stage.
status <- system2("/usr/bin/python3", file.path(SCRIPT_DIR,
  "10_run_recurrent_baseline_onmf.py"),
  env = sprintf("TAHOE_V2_ROOT=%s", root))
if (!identical(status, 0L))
  fail("Stage 10 native oNMF failed with status %d.", status)

scores <- fread(file.path(ARTIFACT, "baseline_program_scores.tsv.gz"))
loadings <- fread(file.path(RESULT, "baseline_program_loadings.tsv.gz"))
stability <- fread(file.path(RESULT, "baseline_program_stability.tsv"))
dominance <- fread(file.path(RESULT, "baseline_program_line_dominance.tsv"))
programs <- paste0("S", seq_len(8L))
if (!setequal(unique(scores$program), programs) ||
    !setequal(unique(loadings$program), programs))
  fail("Stage 10 did not produce the fixed rank-eight program bank.")

## Plot frozen baseline activities with training and held-out lines in blocks.
line_scores <- dcast(scores[baseline_source == "tahoe_dmso" &
    aggregation_level == "line_mean"], program ~ depmap_id,
  value.var = "standardized_program_score")
activity <- as.matrix(line_scores[, -"program"])
rownames(activity) <- line_scores$program
role_table <- unique(scores[baseline_source == "tahoe_dmso" &
    aggregation_level == "line_mean", .(depmap_id, role, cell_line)])
order_block <- function(ids) {
  if (length(ids) < 2L) return(ids)
  ids[hclust(dist(t(activity[, ids, drop = FALSE])),
    method = "average")$order]
}
training_ids <- sort(role_table[role == "training", depmap_id])
held_ids <- sort(role_table[role == "held_out", depmap_id])
column_order <- c(order_block(training_ids), order_block(held_ids))
row_order <- rownames(activity)[hclust(dist(activity), method = "average")$order]
activity <- activity[row_order, column_order, drop = FALSE]
annotation <- role_table[match(column_order, depmap_id)]
stopifnot(identical(colnames(activity), annotation$depmap_id))
hm_activity <- draw_heatmap(activity,
  col_anno_df = data.frame(Role = annotation$role), col_anno_cols = "Role",
  scale_rows = FALSE, value_range = c(-3, 0, 3),
  value_colors = c("#2166AC", "white", "#B2182B"),
  show_row_names = TRUE, show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 5.5),
  row_names_gp = grid::gpar(fontsize = 7),
  cluster_rows = FALSE, cluster_columns = FALSE,
  title = "Rank-eight untreated baseline-state activities")
save_heatmap("10a_baseline_program_activity_heatmap.pdf", hm_activity,
  8.2, 4.8)

## Plot top signed genes in deterministic program blocks.
top_genes <- loadings[order(program, -absolute_weight, gene),
  head(.SD, 18L), by = program]
top_genes <- top_genes[!duplicated(gene)]
setorder(top_genes, program, -absolute_weight, gene)
wide_weights <- dcast(loadings[gene %chin% top_genes$gene], gene ~ program,
  value.var = "signed_weight", fill = 0)
weight_matrix <- as.matrix(wide_weights[, ..programs])
rownames(weight_matrix) <- wide_weights$gene
weight_matrix <- weight_matrix[top_genes$gene, programs, drop = FALSE]
stopifnot(identical(rownames(weight_matrix), top_genes$gene))
limit <- max(abs(weight_matrix))
hm_genes <- draw_heatmap(weight_matrix,
  row_anno_df = data.frame(Program = top_genes$program,
    Direction = top_genes$direction),
  row_anno_cols = c("Program", "Direction"), scale_rows = FALSE,
  value_range = c(-limit, 0, limit),
  value_colors = c("#2166AC", "white", "#B2182B"),
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 5.5),
  cluster_rows = FALSE, cluster_columns = FALSE,
  title = "Frozen recurrent baseline-state top-gene loadings")
save_heatmap("10b_baseline_program_top_gene_heatmap.pdf", hm_genes,
  7.0, 9.2)

## Show repeat stability and the independent leave-one-line-out audit.
loo <- fread(file.path(RESULT,
  "baseline_program_leave_one_line_out_stability.tsv"))
stability_once <- stability[, .SD[1L], by = program]
loo_summary <- loo[, .(loo_q10 = quantile(matched_membership_pearson, 0.1),
  loo_median = median(matched_membership_pearson)), by = program]
p_stability <- ggplot(stability,
    aes(program, matched_pearson, color = "Repeated oNMF")) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35,
    color = "grey45") +
  geom_jitter(width = 0.09, height = 0, alpha = 0.5, size = 1.2) +
  geom_point(data = stability_once, aes(y = q10_matched_pearson),
    shape = 23, fill = "gold", color = "black", size = 2.7) +
  geom_point(data = loo_summary, aes(y = loo_q10,
    color = "Leave one line out"), shape = 18, size = 2.5) +
  scale_color_manual(values = c("Repeated oNMF" = "#2166AC",
    "Leave one line out" = "#B2182B")) +
  labs(x = NULL, y = "Matched loading-membership Pearson r",
    color = NULL, title = "Recurrent baseline-program stability",
    subtitle = "Gold diamonds: repeat q10; red diamonds: leave-one-line-out q10") +
  theme_stage10()
save_gg("10c_baseline_program_stability.pdf", p_stability, 6.0, 4.2)

## Flag stable programs if one training line still dominates their activity.
p_dominance <- ggplot(dominance,
    aes(program, top_1_squared_score_fraction, fill = baseline_source)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.68) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.4,
    color = "#B2182B") +
  scale_fill_manual(values = c(tahoe_dmso = "#2166AC",
    depmap_static = "#B2182B"), labels = c(tahoe_dmso = "Tahoe DMSO",
      depmap_static = "DepMap")) +
  labs(x = NULL, y = "Top-line fraction of squared score",
    fill = "Baseline source", title = "Baseline-program cell-line dominance",
    subtitle = "Programs above 0.50 are flagged, not interpreted as recurrent states") +
  theme_stage10()
save_gg("10d_baseline_program_cell_line_dominance.pdf", p_dominance,
  6.1, 4.2)

## Use DepMap only as an orthogonal frozen-loading reproducibility check.
source_wide <- dcast(scores[aggregation_level %chin% c("line_mean",
    "line_level")], depmap_id + cell_line + role + program ~ baseline_source,
  value.var = "standardized_program_score")
reproducibility <- fread(file.path(RESULT,
  "baseline_dmso_depmap_reproducibility.tsv"))[role == "all"]
p_reproducibility <- ggplot(source_wide,
    aes(tahoe_dmso, depmap_static, color = role)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.75, size = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.45) +
  geom_text(data = reproducibility, aes(x = -Inf, y = Inf,
    label = sprintf("r=%.2f", pearson)), inherit.aes = FALSE,
    hjust = -0.08, vjust = 1.2, size = 2.3) +
  facet_wrap(~program, scales = "free", ncol = 4) +
  scale_color_manual(values = c(training = "#2166AC",
    held_out = "#B2182B")) +
  labs(x = "Tahoe DMSO frozen program score",
    y = "DepMap frozen program score", color = NULL,
    title = "Cross-source reproducibility of untreated baseline states") +
  square_scatter_theme(7.2)
save_gg("10e_baseline_dmso_depmap_reproducibility.pdf",
  p_reproducibility, 8.0, 5.9)

## Record plot ownership and finish only after every required output exists.
figures <- data.table(
  pdf = c("10a_baseline_program_activity_heatmap.pdf",
    "10b_baseline_program_top_gene_heatmap.pdf",
    "10c_baseline_program_stability.pdf",
    "10d_baseline_program_cell_line_dominance.pdf",
    "10e_baseline_dmso_depmap_reproducibility.pdf"),
  producing_script = "scripts/10_learn_recurrent_baseline_states.R")
figures[, exists := file.exists(file.path(RESULT, pdf))]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")
if (!all(figures$exists)) fail("A required Stage 10 PDF was not created.")
writeLines("complete", file.path(RESULT, "stage10_complete.flag"))
message("Stage 10 complete: rank-eight untreated recurrent baseline states fitted and diagnosed.")
