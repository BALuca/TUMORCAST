#!/usr/bin/env Rscript
# Inspect the complete molecular bank and the lapatinib baseline -> displacement
# -> function chain. These calculations do not select or refit the functional model.

suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(ggplot2)
})
source("R/heatmaps.R")
source("R/prediction.R")

output_dir <- "results/14_plot_model_diagnostics"
programs <- sprintf("P%02d", 1:40)
genes_per_program <- 8L
program_colors <- setNames(pal_get(length(programs), "set_base"), programs)
RhpcBLASctl::blas_set_num_threads(2L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
prepared <- qread("data/processed/molecular_inputs.qs")$result
model <- qread("results/08_fit_single_drug_bridge/models/final.qs")
contexts <- read.delim("results/08_fit_single_drug_bridge/data/molecular_contexts.tsv")
ligands <- qread("data/processed/hgf_program_inputs.qs")$result
weights <- read.delim("data/processed/program_bank/frozen_program_weights.tsv")
relationships <- qread("data/processed/gene_relationships.qs")$result

# Independently reconstruct the two-receptor ridge prediction for each held RNA
# recipient. The closed-form solution uses the same standardized design and
# lambda = 1 as the fitted multivariate ridge map.
rna <- prepared$m %>% filter(drug == "Lapatinib ditosylate", dose == 5) %>%
  inner_join(prepared$receptors, by = "depmap_id", relationship = "one-to-one")
x <- as.matrix(rna[, c("ERBB2", "MET")])
y <- as.matrix(rna[, programs])
prediction <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
for (i in seq_len(nrow(y))) {
  train <- seq_len(nrow(y)) != i
  x_mean <- colMeans(x[train, , drop = FALSE])
  x_sd <- apply(x[train, , drop = FALSE], 2, sd)
  y_mean <- colMeans(y[train, , drop = FALSE])
  y_sd <- pmax(apply(y[train, , drop = FALSE], 2, sd), 1e-8)
  z <- sweep(sweep(x[train, , drop = FALSE], 2, x_mean), 2, x_sd, "/")
  target <- sweep(sweep(y[train, , drop = FALSE], 2, y_mean), 2, y_sd, "/")
  beta <- solve(crossprod(z) + sum(train) * diag(ncol(x)), crossprod(z, target))
  prediction[i, ] <- drop(matrix((x[i, ] - x_mean) / x_sd, nrow = 1) %*% beta) * y_sd + y_mean
}
held_predictions <- cbind(rna[, c("depmap_id", "drug", "dose")], as.data.frame(prediction))
write.table(held_predictions, file.path(output_dir, "held_predicted_displacements.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Save the raw HGF-versus-drug vectors underlying the polished opposition panel.
query <- data.frame(depmap_id = "ACH-000248", drug = "Lapatinib ditosylate", dose = 5)
lapatinib <- predict_lapatinib(model, query, prepared$receptors, contexts)$D
hgf <- ligands$program_inputs %>%
  filter(depmap_id == "ACH-000248", source_profile == "Primary_study_mean")
stopifnot(nrow(hgf) == 1L)
association <- model$template$correlations[, "Lapatinib ditosylate"]
opposition <- data.frame(program = programs, lapatinib_displacement = as.numeric(lapatinib),
  HGF_displacement = as.numeric(hgf[, programs]), lapatinib_sensitivity_r = association,
  strong_lapatinib_association = abs(association) >= 0.5)
write.table(opposition, file.path(output_dir, "AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Retain the earlier route-assumption diagnostic. HGF has no observed functional
# correlation pattern identifying a route, so both assumptions are displayed.
# Wilson outcomes do not choose the route, scale, or molecular footprint.
original_hgf <- ligands$original_hgf_contexts
hgf_mean <- as.numeric(original_hgf[original_hgf$context == "Equal-context mean", programs])
context <- model$context
hgf_standardized <- matrix((hgf_mean - context$fm) / context$fs, nrow = 1)
hgf_coordinates <- cbind(1, sweep(hgf_standardized %*% context$fv, 2, context$fscale, "/"))
scale_map <- solve(crossprod(context$Q[1:10, ]) / 10 + diag(c(0, 1, 1, 1)),
  crossprod(context$Q[1:10, ], log(context$sigma[1:10, ])) / 10)
inferred_scale <- as.numeric(exp(hgf_coordinates %*% scale_map))
used_scale <- pmin(pmax(inferred_scale, apply(context$sigma[1:10, ], 2, min)),
  apply(context$sigma[1:10, ], 2, max))
lap_prediction <- predict_lapatinib(model, query, prepared$receptors, contexts)
lap_context <- match("Lapatinib ditosylate|5", paste(contexts$drug, contexts$dose, sep = "|"))
lap_activity <- matrix(0, nrow(original_hgf), 2)
lap_absolute <- as.numeric((lapatinib / context$sigma[lap_context, ]) %*% context$W)
lap_activity[, lap_prediction$route] <- lap_absolute[lap_prediction$route] * lap_prediction$amplitude
hgf_amplitude <- as.numeric(exp(hgf_coordinates %*% model$par$log_amplitude))
hgf_absolute <- sweep(as.matrix(original_hgf[, programs]), 2, used_scale, "/") %*% context$W
route_queries <- vector("list", 2)
for (route in 1:2) {
  hgf_activity <- matrix(0, nrow(original_hgf), 2)
  hgf_activity[, route] <- hgf_absolute[, route] * hgf_amplitude
  products <- cbind(lap_activity[, 1] * hgf_activity[, 1],
    lap_activity[, 1] * hgf_activity[, 2] + lap_activity[, 2] * hgf_activity[, 1],
    lap_activity[, 2] * hgf_activity[, 2])
  interaction <- as.numeric(sweep(products, 2, context$interaction_scale, "/") %*% model$par$Gamma)
  lap_effect <- rowSums(lap_activity)
  hgf_effect <- rowSums(hgf_activity)
  route_queries[[route]] <- data.frame(cell = "AU565", dose = 5,
    hgf_context = original_hgf$context, assumed_HGF_route = route,
    lapatinib_pct_zero_control = 100 * 2^(-lap_effect),
    hgf_alone_pct_zero_control = 100 * 2^(-hgf_effect),
    combination_pct_zero_control = 100 * 2^(-lap_effect - hgf_effect - interaction),
    lapatinib_in_hgf_pct_hgf_control = 100 * 2^(-lap_effect - interaction),
    interaction_log2 = -interaction,
    qualification = paste("Zero-displacement contrasts; HGF route conditional, not identified;",
      "RNA-inferred SD; HGF/Wilson labels never fitted; no lapatinib-HGF training"))
}
write.table(bind_rows(route_queries),
  file.path(output_dir, "au565_hgf_5uM_zero_reference_query.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(program = programs, inferred_SD = inferred_scale, used_SD = used_scale),
  file.path(output_dir, "hgf_scale_transfer.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Show source-derived program relationships without using outcome labels.
W <- prepared$W[relationships$genes, programs]
program_covariance <- crossprod(W, relationships$joint %*% W)
program_scale <- sqrt(diag(program_covariance))
program_correlation <- program_covariance / outer(program_scale, program_scale)
groups <- cutree(hclust(dist(program_correlation), method = "ward.D2"), k = 6)
blocks <- block_row_order(program_correlation, groups)
order <- blocks$order
annotation <- data.frame(Block = paste("Block", blocks$assignment[order]), row.names = programs[order])
heatmap <- draw_heatmap(program_correlation[order, order], scale_rows = FALSE,
  value_range = seq(-1, 1, length.out = length(hm_value_colors)),
  value_colors = hm_value_colors, hm_name = "Program correlation",
  col_anno_df = annotation, col_anno_cols = "Block", row_anno_df = annotation, row_anno_cols = "Block",
  ann_row_args = list(show_legend = FALSE),
  show_row_names = TRUE, show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 10), width = grid::unit(150, "mm"), height = grid::unit(150, "mm"),
  title = "Source-derived relationships between all 40 programs")
pdf(file.path(output_dir, "program_relationships.pdf"), width = 10, height = 9)
ComplexHeatmap::draw(heatmap, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()
write.table(data.frame(program = rownames(program_correlation), program_correlation, check.names = FALSE),
  file.path(output_dir, "program_relationships.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Detail pages retain large gene-symbol labels: four programs per page, with
# eight strongest hard-assigned genes each. Signs orient genes to their program.
pdf(file.path(output_dir, "program_gene_relationships.pdf"), width = 10, height = 9, onefile = TRUE)
for (page in split(programs, ceiling(seq_along(programs) / 4))) {
  selected <- weights %>% filter(program %in% page) %>%
    group_by(program) %>% arrange(desc(abs(weight)), gene, .by_group = TRUE) %>%
    slice_head(n = genes_per_program) %>% ungroup()
  indices <- match(selected$gene, relationships$genes)
  signs <- sign(selected$weight)
  matrix <- relationships$joint[indices, indices] * outer(signs, signs)
  # Membership is fixed by the discovered bank; the display only orders genes
  # within their true program. It does not reassign genes to improve separation.
  within <- vapply(seq_len(nrow(matrix)), function(i) {
    mean(matrix[i, selected$program == selected$program[i]]) -
      mean(matrix[i, selected$program != selected$program[i]])
  }, numeric(1))
  order <- order(selected$program, -within, selected$gene)
  matrix <- matrix[order, order]
  dimnames(matrix) <- list(selected$gene[order], selected$gene[order])
  annotation <- data.frame(Program = selected$program[order], row.names = rownames(matrix))
  heatmap <- draw_heatmap(matrix, scale_rows = FALSE,
    value_range = seq(-1, 1, length.out = length(hm_value_colors)),
    value_colors = hm_value_colors, hm_name = "Signed gene r",
    col_anno_df = annotation, col_anno_cols = "Program", row_anno_df = annotation, row_anno_cols = "Program",
    ann_row_args = list(show_legend = FALSE),
    palette_overrides = list(Program = program_colors[unique(selected$program)]),
    show_row_names = TRUE, show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 11),
    column_names_gp = grid::gpar(fontsize = 11), column_names_rot = 60,
    width = grid::unit(145, "mm"), height = grid::unit(145, "mm"),
    title = paste("Gene relationships:", paste(page, collapse = ", ")))
  ComplexHeatmap::draw(heatmap, annotation_legend_side = "right", heatmap_legend_side = "right")
}
dev.off()
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
