#!/usr/bin/env Rscript

## Interpret the existing EIGHT joint programs. Never reselect or recluster genes.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
source("/work/src/R/lib_heatmap.R")
root <- "/work/tahoe_prism_v2"
input <- file.path(root, "results/20_tahoe_cytosig_programs")
out <- file.path(root, "results/21_program_drug_hierarchy")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
paths <- file.path(input, c("joint_program_membership.tsv", "selected_genes_by_drug.tsv",
  "program_source_coherence.tsv"))
stopifnot(all(file.exists(paths)))
input_md5 <- unname(tools::md5sum(paths))
membership <- as_tibble(fread(paths[1])) %>% arrange(plotting_order)
selected <- as_tibble(fread(paths[2]))
coherence <- as_tibble(fread(paths[3]))
programs <- paste0("P", 1:8)
drugs <- c("erlotinib", "crizotinib", "trametinib")
stopifnot(!anyDuplicated(membership$gene), setequal(membership$program, programs),
  !anyDuplicated(selected[, c("compound", "gene")]),
  all((selected %>% dplyr::count(compound))$n <= 100L),
  setequal(selected$compound, drugs), setequal(coherence$program, programs))
universe <- membership$gene
N <- length(universe)

## Condition on the genes that entered the programs. Using all measured genes
## here would make enrichment nearly inevitable because selection already used
## these drug lists. The original nominal-p gene rule and cap remain unchanged.
drug_genes <- lapply(drugs, function(d) intersect(selected$gene[selected$compound == d], universe))
names(drug_genes) <- drugs
enrichment <- bind_rows(lapply(programs, function(p) {
  pg <- membership$gene[membership$program == p]
  n <- length(pg)
  bind_rows(lapply(drugs, function(d) {
    dg <- drug_genes[[d]]; K <- length(dg); overlap <- intersect(pg, dg); k <- length(overlap)
    contingency <- matrix(c(k, n - k, K - k, N - n - K + k), nrow = 2, byrow = TRUE)
    stopifnot(all(contingency >= 0), sum(contingency) == N)
    ft <- fisher.test(contingency, alternative = "greater")
    tibble(program = p, compound = d, program_genes = n, background_genes = N,
      drug_genes_in_background = K, overlapping_genes = k, expected_overlap = n * K / N,
      fold_enrichment = (k / n) / (K / N), odds_ratio = unname(ft$estimate),
      nominal_p = ft$p.value,
      ## Jeffreys-smoothed proportions permit finite log2 enrichment at zero.
      ## This 0.5 count is ONLY for the descriptive hierarchy/color, never tests.
      hierarchy_log2_enrichment = log2(((k + .5) / (n + 1)) / ((K + .5) / (N + 1))),
      overlap_gene_symbols = paste(overlap, collapse = ";"))
  }))
})) %>% mutate(bh_q = p.adjust(nominal_p, method = "BH"),
  nominal_enrichment = nominal_p < .05, fdr_enrichment = bh_q < .05,
  test = "one-sided Fisher; BH across all 8 x 3 program-drug tests")
fwrite(enrichment, file.path(out, "program_drug_enrichment.tsv"), sep = "\t")
fwrite(membership %>% left_join(select(selected, compound, gene, selection_rank,
  best_training_p, strongest_training_abs_r), by = "gene"),
  file.path(out, "program_drug_gene_membership.tsv"), sep = "\t")

profile <- matrix(enrichment$hierarchy_log2_enrichment, nrow = 8, byrow = TRUE,
  dimnames = list(programs, drugs))
stopifnot(identical(enrichment$program, rep(programs, each = 3L)), all(is.finite(profile)))
## This clusters EIGHT already fixed programs by their three drug enrichments.
## It does not fit another gene model, choose k, or establish causal ancestry.
tree <- hclust(dist(profile, method = "euclidean"), method = "average")
program_order <- programs[tree$order]
sizes <- table(factor(membership$program, levels = programs))
descendants <- vector("list", nrow(tree$merge))
edges <- list()
for (i in seq_len(nrow(tree$merge))) {
  children <- tree$merge[i, ]
  leaves <- lapply(children, function(j) if (j < 0) programs[-j] else descendants[[j]])
  descendants[[i]] <- unlist(leaves, use.names = FALSE)
  for (j in 1:2) edges[[length(edges) + 1L]] <- tibble(parent = paste0("node", i),
    child = if (children[j] < 0) programs[-children[j]] else paste0("node", children[j]),
    parent_height = tree$height[i],
    child_height = if (children[j] < 0) 0 else tree$height[children[j]],
    child_programs = paste(leaves[[j]], collapse = ";"))
}
fwrite(bind_rows(edges), file.path(out, "program_drug_hierarchy_edges.tsv"), sep = "\t")
annotation <- bind_rows(lapply(programs, function(p) {
  e <- filter(enrichment, program == p)
  tibble(program = p, n_genes = unname(sizes[p]),
    nominal_enriched_drugs = paste(e$compound[e$nominal_enrichment], collapse = ";"),
    fdr_enriched_drugs = paste(e$compound[e$fdr_enrichment], collapse = ";"),
    n_nominal_enriched_drugs = sum(e$nominal_enrichment),
    n_fdr_enriched_drugs = sum(e$fdr_enrichment))
})) %>% left_join(coherence %>% select(-n_genes), by = "program") %>%
  mutate(hierarchy_order = match(program, program_order)) %>% arrange(hierarchy_order)
fwrite(annotation, file.path(out, "program_hierarchy_annotation.tsv"), sep = "\t")
fwrite(tibble(input = basename(paths), path = paths, md5 = input_md5,
  role = c("Frozen eight programs and original gene order", "Frozen capped training drug gene lists",
    "Existing source support; annotation only")), file.path(out, "input_manifest.tsv"), sep = "\t")

palette <- setNames(grDevices::hcl.colors(8, "Dark 3"), programs)
support <- coherence$sharing[match(programs, coherence$program)]
support <- c("Coherent in both sources" = "Both", "Coherent in Tahoe only" = "Tahoe",
  "Coherent in CytoSig only" = "CytoSig", "Not coherent in either source" = "Neither")[support]
ann <- data.frame(Program = factor(programs, levels = programs),
  Support = factor(unname(support), levels = c("Both", "Tahoe", "CytoSig", "Neither")), row.names = programs)
stopifnot(identical(rownames(profile), rownames(ann)))
text_cells <- matrix(paste0(enrichment$overlapping_genes, "/", enrichment$program_genes,
  ifelse(enrichment$fdr_enrichment, " *", ifelse(enrichment$nominal_enrichment, " +", ""))),
  nrow = 8, byrow = TRUE)
grDevices::cairo_pdf(file.path(out, "program_drug_hierarchy.pdf"), width = 9.1, height = 6.7, onefile = TRUE)
hm <- draw_heatmap(profile, row_anno_df = ann, row_anno_cols = c("Program", "Support"),
  palette_overrides = list(Program = palette,
    Support = c(Both = "#537E55", Tahoe = "#4179A8", CytoSig = "#B47534", Neither = "#BCBCBC")),
  ann_row_args = list(show_legend = c(Program = FALSE, Support = TRUE), show_annotation_name = FALSE,
    annotation_name_gp = grid::gpar(fontsize = 11),
    annotation_legend_param = list(Support = list(title = "Source support",
      title_gp = grid::gpar(fontsize = 12, fontface = "bold"), labels_gp = grid::gpar(fontsize = 12)))),
  scale_rows = FALSE, value_range = c(-3, 0, 3),
  value_colors = c("#2166AC", "white", "#B2182B"), hm_name = "log2 enrichment",
  cluster_rows = tree, row_dend_reorder = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_labels = paste0(programs, "  (", as.integer(sizes), " genes)"),
  row_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  column_labels = c("Erlotinib", "Crizotinib", "Trametinib"), column_names_rot = 0,
  column_names_gp = grid::gpar(fontsize = 13),
  row_dend_width = grid::unit(22, "mm"),
  column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
  title = "Eight fixed programs | Drug-gene enrichment hierarchy",
  heatmap_legend_param = list(at = -3:3, title = "log2 enrichment",
    title_gp = grid::gpar(fontsize = 12, fontface = "bold"), labels_gp = grid::gpar(fontsize = 11)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(text_cells[i, j], x, y, gp = grid::gpar(fontsize = 13,
      col = if (abs(profile[i, j]) >= 2) "white" else "#1F2937"))
  })
ComplexHeatmap::draw(hm, merge_legends = TRUE, heatmap_legend_side = "right",
  annotation_legend_side = "right", padding = grid::unit(c(23, 5, 6, 5), "mm"))
grid::grid.text(paste0("Cells: drug-associated genes / program genes. * BH q < 0.05; + nominal p < 0.05 only.\n",
  "Background: all ", N, " program genes. Colors and tree use log2 enrichment with a 0.5-count correction.\n",
  "Average-linkage Euclidean hierarchy of enrichment profiles; branch lengths describe similarity, not causality."),
  x = .5, y = grid::unit(10, "mm"), gp = grid::gpar(fontsize = 10))
dev.off()

## A separate, readable tree retains the same program and support colors.
leaf_y <- setNames(rev(seq_along(program_order)), program_order)
node_y <- vapply(descendants, function(ps) mean(range(leaf_y[ps])), numeric(1))
segments <- list()
for (i in seq_len(nrow(tree$merge))) {
  children <- tree$merge[i, ]
  ys <- vapply(children, function(j) if (j < 0) leaf_y[programs[-j]] else node_y[j], numeric(1))
  xs <- vapply(children, function(j) if (j < 0) 0 else tree$height[j], numeric(1))
  segments[[length(segments) + 1L]] <- tibble(x = -tree$height[i], xend = -tree$height[i], y = min(ys), yend = max(ys))
  segments[[length(segments) + 1L]] <- tibble(x = -tree$height[i], xend = -xs, y = ys, yend = ys)
}
leaf <- annotation %>% mutate(y = unname(leaf_y[program]),
  enriched_label = ifelse(n_fdr_enriched_drugs > 0, paste0(fdr_enriched_drugs, " *"),
    ifelse(n_nominal_enriched_drugs > 0, paste0(nominal_enriched_drugs, " +"), "No enrichment at p < 0.05")))
g <- ggplot() + geom_segment(data = bind_rows(segments), aes(x, y, xend = xend, yend = yend),
    linewidth = .75, color = "#64748B") +
  geom_point(data = leaf, aes(0, y, color = program), size = 3.4) +
  geom_text(data = leaf, aes(.10, y, label = paste0(program, "  |  ", n_genes, " genes")), hjust = 0, size = 4.3) +
  geom_text(data = leaf, aes(1.05, y, label = enriched_label), hjust = 0, size = 3.9) +
  scale_color_manual(values = palette, guide = "none") +
  scale_x_continuous(limits = c(-max(tree$height) - .05, 3), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(.5, 8.5), expand = expansion(mult = 0)) +
  labs(title = "Program hierarchy from drug-gene enrichment",
    subtitle = "The eight joint gene programs remain unchanged",
    caption = "* BH q < 0.05 across 24 tests; + nominal p < 0.05 only.\nBranches group similar drug-enrichment profiles; they do not establish upstream/downstream order.") +
  theme_void(base_size = 13) + theme(plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10, hjust = 0),
    plot.margin = margin(10, 10, 10, 10))
ggsave(file.path(out, "program_drug_tree.pdf"), g, width = 10.3, height = 5.6, device = cairo_pdf)

stopifnot(identical(unname(tools::md5sum(paths)), input_md5))
writeLines(c("# Drug enrichment of the eight frozen joint programs", "",
  "This step reads the existing program membership and top-100-per-drug selections. It never reselects genes, changes gene order, changes cluster membership, or fits programs in Wilson data.",
  "Each of 24 program-drug comparisons uses a one-sided Fisher exact test. The background is the 286 genes actually assigned to programs, so enrichment is conditional on the selection step. A gene selected by two drugs contributes to both drug sets. Gene lists are never topped up.",
  "The enrichment table reports raw overlap counts, expected counts, uncorrected fold enrichment, exact-test p, and BH q across all 24 comparisons. These are descriptive, selection-conditioned enrichment tests, not independent validation of drug specificity.",
  "The hierarchy clusters the eight programs by Euclidean distance between their three log2 enrichment values, with average linkage. A 0.5-count correction in each binary proportion gives finite log2 values when overlap is zero; it is used only for color and hierarchy, never Fisher tests. The explicit formula is log2(((k+0.5)/(n+1))/((K+0.5)/(N+1))).",
  "The program-drug hierarchy PDF shows this tree beside enrichment and overlap counts; the separate tree adds enrichment labels. Stars denote BH q<0.05; plus signs denote nominal p<0.05 only. Both levels are reported; no new threshold selects or removes programs.",
  "Source support is copied from the preceding analysis and is only an annotation: positive mean silhouette of the same fixed program within Tahoe and/or CytoSig. It is not used to build the drug-enrichment hierarchy.",
  "No ancestral or causal direction is inferred from these enrichment profiles. No multi-drug enriched trunk is imposed when the enrichment tests do not support one.",
  "", "Run: Rscript /work/tahoe_prism_v2/scripts/21_program_drug_hierarchy.R"), file.path(out, "README.md"))
print(enrichment %>% filter(nominal_enrichment) %>% select(program, compound,
  overlapping_genes, program_genes, fold_enrichment, nominal_p, bh_q))
message("Completed drug-enrichment hierarchy for the SAME eight programs.")
