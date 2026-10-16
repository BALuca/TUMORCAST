#!/usr/bin/env Rscript

## RAF/EGF pilot: discovery in Tahoe + PRISM + CytoSig EGF,
## followed by frozen-program molecular validation in GSE127988.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cluster)
})
source("/work/src/R/lib_heatmap.R")
RhpcBLASctl::blas_set_num_threads(8L)
root <- "/work/tahoe_prism_v2"
cache <- "/work/cache/tahoe_prism_full_chain"
pdex <- file.path(cache, "tahoe_pdex_selected")
out <- file.path(root, "results/22_raf_egf_bypass")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
drugs <- c("vemurafenib", "cobimetinib")
paths <- c(
  effect_metadata = file.path(root, "artifacts/01_inputs/response_metadata.tsv.gz"),
  all_genes = file.path(pdex, "genes.tsv"),
  normalization = file.path(root, "artifacts/02_response_programs/normalization_records.tsv.gz"),
  control_effects = file.path(pdex, "pdex_plate_reference.npy"),
  control_metadata = file.path(pdex, "pdex_plate_reference_metadata.tsv"),
  prism = file.path(root, "artifacts/01_inputs/function_outcomes.tsv.gz"),
  split = file.path(root, "results/01_inputs/fold_assignments.tsv"),
  cytosig_metadata = "/work/data/01_CytoSig_data/meta_info.txt",
  cytosig_effects = "/work/data/01_CytoSig_data/diff.merge.gz")
stopifnot(all(file.exists(paths)))
paths <- vapply(paths, normalizePath, character(1), mustWork = TRUE)
meta <- as_tibble(fread(paths[["effect_metadata"]])) %>%
  mutate(compound = tolower(trimws(drug))) %>% filter(compound %in% drugs, available == 1L)
genes <- fread(paths[["all_genes"]]) %>% as_tibble() %>% arrange(gene_index)
stopifnot(nrow(genes) == 62710L, !anyDuplicated(genes$feature),
  identical(as.integer(genes$gene_index), seq_len(nrow(genes)) - 1L),
  "vemurafenib" %in% meta$compound)
folds <- as_tibble(fread(paths[["split"]])) %>%
  transmute(depmap_id, tahoe_cell_name, split = if_else(outer_fold == 1L, "Holdout", "Training"))
stopifnot(!anyDuplicated(folds$depmap_id), sum(folds$split == "Training") == 26L,
  sum(folds$split == "Holdout") == 7L)
prism <- as_tibble(fread(paths[["prism"]])) %>% filter(chain_observed == 1L) %>%
  mutate(observation_id = paste(depmap_id, perturbation_id, sep = "__"))
norm <- as_tibble(fread(paths[["normalization"]]))
stopifnot(!anyDuplicated(prism$observation_id), !anyDuplicated(norm$observation_id))
meta <- meta %>% left_join(folds, by = "depmap_id") %>%
  left_join(select(prism, observation_id, prism_sensitization, prism_dose_uM), by = "observation_id") %>%
  left_join(select(norm, observation_id, treated_total, reference_total,
    treated_correction_factor, contributing_plates, negative_tolerance), by = "observation_id") %>%
  arrange(compound, dose, depmap_id)
stopifnot(!anyNA(meta$split), !anyNA(meta$treated_total), !anyDuplicated(meta$observation_id))
## Adequacy is assessed from overlap counts only, before expression selection.
coverage <- meta %>% group_by(compound, dose, split) %>% summarise(
  n_molecular_lines = n_distinct(depmap_id),
  n_prism_lines = n_distinct(depmap_id[is.finite(prism_sensitization)]), .groups = "drop")
adequate <- coverage %>% filter(split == "Training") %>% group_by(compound) %>%
  summarise(include_in_discovery = any(n_prism_lines >= 8L), .groups = "drop")
stopifnot(adequate$include_in_discovery[adequate$compound == "vemurafenib"])
drugs <- adequate$compound[adequate$include_in_discovery]
meta <- meta %>% filter(compound %in% drugs)
fwrite(left_join(coverage, adequate, by = "compound"), file.path(out, "drug_overlap.tsv"), sep = "\t")
meta$matrix_row <- seq_len(nrow(meta))

## Read rows of saved processed pseudobulk arrays in place, using their original
## gene order. The same saved sample-level normalization factors apply to ALL
## genes; no raw single cells or upstream feature-selection code is rerun.
read_npy_rows <- function(path, rows_zero_based) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  stopifnot(identical(magic, "\x93NUMPY"), version[1] %in% 1:3)
  len <- readBin(con, "integer", 1L, size = if (version[1] == 1L) 2L else 4L,
    signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", len))
  stopifnot(grepl("'descr': '<f4'", header, fixed = TRUE),
    grepl("'fortran_order': False", header, fixed = TRUE))
  shape <- as.integer(trimws(strsplit(sub(".*'shape': \\(([^)]*)\\).*", "\\1", header), ",")[[1L]]))
  shape <- shape[!is.na(shape)]
  stopifnot(length(shape) == 2L, all(rows_zero_based >= 0 & rows_zero_based < shape[1]))
  offset <- seek(con)
  answer <- matrix(NA_real_, length(rows_zero_based), shape[2])
  for (i in seq_along(rows_zero_based)) {
    seek(con, offset + 4 * shape[2] * rows_zero_based[i], origin = "start")
    v <- readBin(con, "numeric", n = shape[2], size = 4L, endian = "little")
    stopifnot(length(v) == shape[2]); answer[i, ] <- v
  }
  answer
}
control_meta <- as_tibble(fread(paths[["control_metadata"]]))
plate_keys <- unlist(lapply(seq_len(nrow(meta)), function(i)
  paste(meta$depmap_id[i], strsplit(as.character(meta$contributing_plates[i]), ",", fixed = TRUE)[[1]], sep = "__")))
control_meta$key <- paste(control_meta$depmap_id, control_meta$plate, sep = "__")
control_meta <- control_meta %>% filter(key %in% plate_keys)
stopifnot(!anyDuplicated(control_meta$key), setequal(control_meta$key, plate_keys))
controls <- read_npy_rows(paths[["control_effects"]], control_meta$array_row)
rownames(controls) <- control_meta$key
effects <- matrix(NA_real_, nrow(meta), nrow(genes), dimnames = list(meta$observation_id, genes$feature))
for (path in unique(meta$shard_path)) {
  ids <- which(meta$shard_path == path)
  raw <- read_npy_rows(file.path(cache, path), meta$shard_row[ids])
  for (j in seq_along(ids)) {
    i <- ids[j]
    keys <- paste(meta$depmap_id[i], strsplit(as.character(meta$contributing_plates[i]), ",", fixed = TRUE)[[1]], sep = "__")
    reference <- colMeans(controls[keys, , drop = FALSE])
    treated <- raw[j, ] + reference
    stopifnot(all(is.finite(treated)), min(treated) >= -meta$negative_tolerance[i])
    treated <- pmax(treated, 0)
    effects[i, ] <- log2(1 + treated / meta$treated_total[i] * 10000 /
      meta$treated_correction_factor[i]) - log2(1 + reference / meta$reference_total[i] * 10000)
  }
}
stopifnot(all(is.finite(effects)))
rm(controls, raw)
message("Reconstructed existing normalized pseudobulk effects for all ", ncol(effects),
  " genes, ", nrow(effects), " observations; no activity or old gene gate.")

## Discovery only: exact EGF profiles; never use Gerosa/GSE127988 here.
cyto_meta <- as_tibble(fread(paths[["cytosig_metadata"]], skip = 1L, header = FALSE,
  col.names = c("profile_id", "treatment", "condition", "duration", "dose"))) %>%
  filter(toupper(trimws(treatment)) == "EGF") %>%
  mutate(study = sub(".*@(GSE[0-9]+|E-MTAB-[0-9]+).*", "\\1", profile_id))
stopifnot(nrow(cyto_meta) >= 3L, !any(grepl("GSE127988", cyto_meta$profile_id)),
  all(grepl("^(GSE|E-MTAB-)", cyto_meta$study)))
con <- gzfile(paths[["cytosig_effects"]], "rt")
profile_header <- strsplit(readLines(con, 1L), "\t", fixed = TRUE)[[1]]; close(con)
stopifnot(all(cyto_meta$profile_id %in% profile_header))
cyto_dt <- suppressWarnings(fread(paths[["cytosig_effects"]], header = TRUE,
  col.names = c("gene", cyto_meta$profile_id), check.names = FALSE,
  select = c(1L, match(cyto_meta$profile_id, profile_header) + 1L), showProgress = FALSE))
stopifnot(!anyDuplicated(cyto_dt$gene))
universe <- intersect(genes$feature, cyto_dt$gene)
cyto <- t(as.matrix(cyto_dt[match(universe, gene), -1L, with = FALSE]))
colnames(cyto) <- universe; rownames(cyto) <- cyto_meta$profile_id
cyto[!is.finite(cyto)] <- NA_real_
studies <- sort(unique(cyto_meta$study))
study_effect <- t(vapply(studies, function(s)
  apply(cyto[cyto_meta$study == s, , drop = FALSE], 2L, median, na.rm = TRUE), numeric(length(universe))))
colnames(study_effect) <- universe; rownames(study_effect) <- studies
egf <- tibble(gene = universe, n_egf_profiles = colSums(is.finite(cyto)),
  n_egf_studies = colSums(is.finite(study_effect)),
  egf_median_effect = apply(study_effect, 2L, median, na.rm = TRUE))
## EGF feature ranking uses study-balanced median signed response. It does not
## assert significance across dependent time points or select on validation.
egf_selected <- egf %>% filter(n_egf_studies >= 3L, is.finite(egf_median_effect)) %>%
  arrange(desc(abs(egf_median_effect)), gene) %>% slice_head(n = 100L) %>%
  mutate(selection_rank = row_number())

pearson_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(tibble(n = n, pearson_r = NA_real_, nominal_p = NA_real_))
  r <- cor(x[ok], y[ok])
  tibble(n = n, pearson_r = r, nominal_p = 2 * pt(abs(r) * sqrt((n - 2) /
    pmax(0, 1 - r^2)), df = n - 2, lower.tail = FALSE))
}
gene_tests <- list()
for (drug in drugs) for (dose in sort(unique(meta$dose[meta$compound == drug]))) {
  ids <- which(meta$compound == drug & meta$dose == dose & meta$split == "Training" &
    is.finite(meta$prism_sensitization))
  n <- length(ids); y <- meta$prism_sensitization[ids]
  r <- rep(NA_real_, length(universe))
  if (n >= 8L && sd(y) > 0) r <- suppressWarnings(as.numeric(cor(effects[ids, universe], y)))
  p <- 2 * pt(abs(r) * sqrt((n - 2) / pmax(0, 1 - r^2)), df = n - 2, lower.tail = FALSE)
  gene_tests[[length(gene_tests) + 1L]] <- tibble(compound = drug, dose_uM = dose,
    gene = universe, n_training_lines = n, pearson_r = r, nominal_p = p,
    selected_at_dose = is.finite(p) & p < .05)
}
gene_tests <- bind_rows(gene_tests)
drug_selected <- gene_tests %>% filter(selected_at_dose) %>% group_by(compound, gene) %>%
  summarise(max_abs_r = max(abs(pearson_r)), minimum_p = min(nominal_p),
    significant_doses_uM = paste(sort(dose_uM), collapse = ";"), .groups = "drop") %>%
  arrange(compound, desc(max_abs_r), minimum_p, gene) %>% group_by(compound) %>%
  mutate(selection_rank = row_number()) %>% slice_head(n = 100L) %>% ungroup()
stopifnot("vemurafenib" %in% drug_selected$compound,
  all((drug_selected %>% dplyr::count(compound))$n <= 100L))
selected <- sort(union(drug_selected$gene, egf_selected$gene))
audit <- tibble(gene = selected, used_for_clustering = TRUE, exclusion_reason = NA_character_)
training <- which(meta$split == "Training")
r_tahoe <- suppressWarnings(cor(effects[training, selected, drop = FALSE]))
r_cyto <- suppressWarnings(cor(cyto[, selected, drop = FALSE], use = "pairwise.complete.obs"))
pair_n <- crossprod(is.finite(cyto[, selected, drop = FALSE]))
r_cyto[pair_n < 3L] <- NA_real_
invalid <- !is.finite(r_tahoe) | !is.finite(r_cyto)
## Remove only mathematically unestimable correlations, with an explicit audit.
## Do not refill a capped list or turn missing correlations into zero.
while (any(invalid)) {
  bad <- tibble(gene = selected, bad_pairs = rowSums(invalid),
    observed = colSums(is.finite(cyto[, selected, drop = FALSE]))) %>%
    arrange(desc(bad_pairs), observed, gene)
  remove <- bad$gene[1]; keep <- selected != remove
  audit$used_for_clustering[audit$gene == remove] <- FALSE
  audit$exclusion_reason[audit$gene == remove] <- "Gene-pair correlation not estimable in one or both sources"
  selected <- selected[keep]
  r_tahoe <- r_tahoe[keep, keep, drop = FALSE]; r_cyto <- r_cyto[keep, keep, drop = FALSE]
  invalid <- !is.finite(r_tahoe) | !is.finite(r_cyto)
}
stopifnot(length(selected) > 3L, length(selected) <= 300L)
matrices <- list(Joint = (r_tahoe + r_cyto) / 2, Tahoe = r_tahoe, CytoSig_EGF = r_cyto)
## One joint clustering. The new architecture fixes THREE programs. Average
## linkage is valid for correlation dissimilarity without Euclidean assumptions.
tree <- hclust(as.dist(1 - matrices$Joint), method = "average")
raw_groups <- cutree(tree, k = 3L)
gene_order <- selected[tree$order]
program_levels <- paste0("P", 1:3)
programs <- setNames(paste0("P", match(raw_groups[gene_order], unique(raw_groups[gene_order]))), gene_order)
palette <- setNames(c("#247A96", "#B95B32", "#6E59A5"), program_levels)
membership <- tibble(gene = gene_order, program = unname(programs), gene_order = seq_along(gene_order),
  selected_vemurafenib = gene %in% drug_selected$gene[drug_selected$compound == "vemurafenib"],
  selected_cobimetinib = gene %in% drug_selected$gene[drug_selected$compound == "cobimetinib"],
  selected_egf = gene %in% egf_selected$gene) %>% left_join(egf, by = "gene")
stopifnot(all(membership$selected_vemurafenib | membership$selected_cobimetinib | membership$selected_egf))
annotation <- data.frame(Program = factor(programs, levels = program_levels), row.names = gene_order)
coherence <- tibble(program = program_levels,
  n_genes = as.integer(table(factor(programs, levels = program_levels))))
for (s in names(matrices)) {
  distance <- 1 - matrices[[s]][gene_order, gene_order]
  distance[distance < 0] <- 0
  sil <- cluster::silhouette(as.integer(factor(programs, levels = program_levels)), as.dist(distance))
  coherence[[paste0(tolower(s), "_silhouette")]] <- vapply(program_levels,
    function(p) mean(sil[programs == p, "sil_width"]), numeric(1))
}
coherence <- coherence %>% mutate(candidate_role = case_when(
  tahoe_silhouette > 0 & cytosig_egf_silhouette > 0 ~ "Shared downstream candidate",
  cytosig_egf_silhouette > 0 ~ "EGF/receptor candidate",
  tahoe_silhouette > 0 ~ "Drug-specific candidate",
  TRUE ~ "Unresolved"))
## Freeze an equal-weight gene score. All genes within a program retain the
## same sign; a single program-level sign establishes a discovery orientation.
## Receptor-only programs face EGF induction. Others face PRISM resistance.
orientation <- list()
for (p in program_levels) {
  g <- gene_order[programs == p]
  dose_r <- vapply(sort(unique(meta$dose[meta$compound == "vemurafenib"])), function(d) {
    idx <- which(meta$split == "Training" & meta$compound == "vemurafenib" & meta$dose == d)
    pearson_test(rowMeans(effects[idx, g, drop = FALSE]), -meta$prism_sensitization[idx])$pearson_r
  }, numeric(1))
  egf_response <- mean(egf$egf_median_effect[match(g, egf$gene)], na.rm = TRUE)
  direction <- if (coherence$candidate_role[coherence$program == p] == "EGF/receptor candidate")
    sign(egf_response) else sign(median(dose_r, na.rm = TRUE))
  if (!is.finite(direction) || direction == 0) direction <- 1
  orientation[[p]] <- tibble(program = p, frozen_sign = direction,
    median_training_r_resistance_unoriented = median(dose_r, na.rm = TRUE),
    oriented_egf_effect = direction * egf_response)
}
coherence <- left_join(coherence, bind_rows(orientation), by = "program")
membership <- membership %>% left_join(select(coherence, program, candidate_role, frozen_sign), by = "program") %>%
  group_by(program) %>% mutate(frozen_weight = frozen_sign / n()) %>% ungroup() %>% arrange(gene_order)
weights <- matrix(0, nrow(membership), 3L, dimnames = list(gene_order, program_levels))
weights[cbind(seq_len(nrow(membership)), match(membership$program, program_levels))] <- membership$frozen_weight
stopifnot(all(abs(colSums(abs(weights)) - 1) < 1e-12))
predictions <- coherence %>% transmute(program, candidate_role, frozen_sign,
  expected_egf_direction = if_else(candidate_role %in% c("Shared downstream candidate", "EGF/receptor candidate"),
    "increase", "exploratory"),
  egf_discovery_agrees = oriented_egf_effect > 0,
  expected_cobimetinib_direction = if_else(candidate_role == "Shared downstream candidate", "decrease", "exploratory"),
  primary_window = "Trapezoidal mean displacement over the full 0-8 h EGF time course",
  statistical_unit = "Biological set/day (GEO set id 3 and 4); technical wells averaged",
  score_definition = "Program sign times mean gene log2(1+expression per 10000); differences from matched reference",
  missing_gene_rule = "Exact symbols only; frozen relative weights renormalized over measured genes; report coverage; no replacements",
  caution = "Candidate roles describe source coherence, not proven receptor specificity or causal ancestry")
fwrite(membership, file.path(out, "frozen_program_genes.tsv"), sep = "\t", na = "NA")
fwrite(coherence, file.path(out, "frozen_program_annotation.tsv"), sep = "\t", na = "NA")
fwrite(predictions, file.path(out, "frozen_predictions.tsv"), sep = "\t")
freeze_paths <- file.path(out, c("frozen_program_genes.tsv", "frozen_program_annotation.tsv", "frozen_predictions.tsv"))
freeze_hashes <- tools::md5sum(freeze_paths)
fwrite(tibble(path = freeze_paths, md5 = unname(freeze_hashes),
  frozen_before_validation_expression_load = TRUE), file.path(out, "freeze_manifest.tsv"), sep = "\t")
fwrite(gene_tests, file.path(out, "gene_prism_training_tests.tsv.gz"), sep = "\t", na = "NA")
fwrite(drug_selected, file.path(out, "selected_genes_by_drug.tsv"), sep = "\t")
fwrite(egf_selected, file.path(out, "selected_egf_genes.tsv"), sep = "\t")
fwrite(cyto_meta, file.path(out, "cytosig_egf_profiles.tsv"), sep = "\t")
fwrite(audit, file.path(out, "gene_selection_audit.tsv"), sep = "\t", na = "NA")
fwrite(membership %>% group_by(program) %>% summarise(n_genes = n(),
  vemurafenib_selected = sum(selected_vemurafenib), cobimetinib_selected = sum(selected_cobimetinib),
  egf_selected = sum(selected_egf), .groups = "drop"),
  file.path(out, "program_gene_origin_summary.tsv"), sep = "\t")
fwrite(tibble(object = names(paths), path = unname(paths),
  mapping = c("observation_id; available=1; drug; dose uM; shard_path; shard_row zero based",
    "All feature values, original gene_index order", "Saved pseudobulk normalization factors and contributing plates",
    "Saved DMSO pseudobulk controls", "array_row zero based; depmap_id; plate",
    "chain_observed=1; prism_sensitization higher means more sensitive",
    "outer_fold 1 holdout, otherwise training; original split",
    "Exact EGF treatment; study extracted from profile_id", "Signed gene-level treatment-control profiles")),
  file.path(out, "input_manifest.tsv"), sep = "\t")
message("FROZEN before validation expression: ", length(gene_order), " genes, three programs; ",
  length(universe), " shared measured genes; ", nrow(cyto), " EGF profiles from ", length(studies), " studies.")
print(coherence)

## Shared publication heatmap wrapper. No plot-specific gene clustering.
save_heatmap <- function(mat, filename, title, range, legend, column_annotation = NULL,
    row_annotation = annotation, labelled = FALSE, column_labels = NULL, width = NULL, height = NULL,
    row_labels = rownames(mat)) {
  stopifnot(is.matrix(mat))
  if (!is.null(row_annotation)) stopifnot(identical(rownames(mat), rownames(row_annotation)))
  if (!is.null(column_annotation)) stopifnot(identical(colnames(mat), rownames(column_annotation)))
  if (is.null(width)) width <- if (labelled && nrow(mat) == ncol(mat)) max(9, ncol(mat) * .155 + 3) else 12
  if (is.null(height)) height <- if (labelled) max(6, nrow(mat) * .155 + 2.4) else 8
  grDevices::cairo_pdf(file.path(out, filename), width = width, height = height)
  ht <- draw_heatmap(mat, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    row_anno_df = row_annotation, row_anno_cols = if (!is.null(row_annotation)) names(row_annotation) else NULL,
    col_anno_df = column_annotation, col_anno_cols = if (!is.null(column_annotation)) names(column_annotation) else NULL,
    palette_overrides = list(Program = palette, Arm = c(Vem = "#247A96", `Vem + Cobi` = "#B95B32")),
    ann_row_args = list(show_legend = FALSE),
    hm_name = legend, value_range = range, value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = labelled, row_labels = row_labels,
    show_column_names = !is.null(column_labels), column_labels = column_labels,
    row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 11),
    column_title_gp = grid::gpar(fontsize = 15, fontface = "bold"), title = title, na_col = "#BDBDBD",
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 11)))
  ComplexHeatmap::draw(ht, merge_legends = TRUE, padding = grid::unit(c(5, 5, 5, 5), "mm"))
  if (identical(rownames(mat), gene_order)) ComplexHeatmap::decorate_heatmap_body(legend, {
    for (b in head(cumsum(as.integer(table(factor(programs, levels = program_levels)))), -1L))
      grid::grid.lines(c(0, 1), rep(1 - b/nrow(mat), 2), gp = grid::gpar(col = "#777777", lwd = .45))
    if (identical(colnames(mat), gene_order)) for (b in head(cumsum(as.integer(table(factor(programs, levels = program_levels)))), -1L))
      grid::grid.lines(rep(b/ncol(mat), 2), c(0, 1), gp = grid::gpar(col = "#777777", lwd = .45))
  })
  dev.off()
}
write_matrix <- function(m, filename) fwrite(data.table(gene = rownames(m), m, check.names = FALSE),
  file.path(out, filename), sep = "\t", na = "NA")
for (s in names(matrices)) {
  m <- matrices[[s]][gene_order, gene_order]
  write_matrix(m, paste0(tolower(s), "_correlations.tsv.gz"))
  for (labelled in c(FALSE, TRUE)) save_heatmap(m,
    paste0(tolower(s), "_correlation", if (labelled) "_labelled" else "", ".pdf"),
    paste0(s, " | ", length(gene_order), " genes | Three fixed joint programs"), c(-1, 0, 1), "Pearson r",
    column_annotation = annotation, labelled = labelled, column_labels = if (labelled) gene_order else NULL)
}
curve <- bind_rows(lapply(2:min(8, length(gene_order) - 1L), function(k) tibble(rank = k,
  mean_silhouette = mean(cluster::silhouette(cutree(tree, k), as.dist(1 - matrices$Joint))[, "sil_width"]))))
fwrite(curve, file.path(out, "rank_silhouette.tsv"), sep = "\t")
p <- ggplot(curve, aes(rank, mean_silhouette)) + geom_line(color = "#247A96") + geom_point(size = 2) +
  geom_point(data = filter(curve, rank == 3L), color = "#B95B32", size = 4) +
  labs(x = "Number of clusters", y = "Mean silhouette", title = "Three programs fixed by the architecture",
    subtitle = "Ranks 2-8 are a diagnostic; they do not select the model") + theme_classic(base_size = 13)
ggsave(file.path(out, "rank_silhouette.pdf"), p, width = 6.4, height = 4.5, device = cairo_pdf)
sample_orders <- list()
for (s in c("Tahoe", "CytoSig_EGF")) {
  m <- if (s == "Tahoe") t(effects[training, gene_order]) else t(cyto[, gene_order])
  z <- sweep(sweep(m, 1, rowMeans(m, na.rm = TRUE), "-"), 1, apply(m, 1, sd, na.rm = TRUE), "/")
  activity <- vapply(program_levels, function(p) colMeans(z[programs == p, , drop = FALSE], na.rm = TRUE), numeric(ncol(z)))
  activity[!is.finite(activity)] <- 0
  dominant <- max.col(activity, ties.method = "first")
  ord <- unlist(lapply(1:3, function(k) {
    ids <- which(dominant == k)
    if (length(ids) < 3L) ids else ids[hclust(dist(activity[ids, , drop = FALSE]), method = "average")$order]
  }))
  ca <- data.frame(Program = factor(program_levels[dominant[ord]], levels = program_levels), row.names = colnames(z)[ord])
  sample_orders[[s]] <- tibble(source = s, sample_id = colnames(z)[ord], order = seq_along(ord),
    dominant_program = program_levels[dominant[ord]])
  if (s == "Tahoe") {
    mm <- meta[match(colnames(z)[ord], meta$observation_id), ]
    ca$Drug <- mm$compound; ca$Dose <- factor(mm$dose)
  }
  write_matrix(m[, ord], paste0(tolower(s), "_gene_sample_effects.tsv.gz"))
  for (labelled in c(FALSE, TRUE)) save_heatmap(z[, ord],
    paste0(tolower(s), "_gene_samples", if (labelled) "_labelled" else "", ".pdf"),
    paste0(s, " | Fixed genes; samples grouped by program"), c(-2, 0, 2), "Gene z-score", ca, labelled = labelled)
}
fwrite(bind_rows(sample_orders), file.path(out, "discovery_sample_order.tsv"), sep = "\t")
score_matrix <- effects[, gene_order] %*% weights
scores <- bind_rows(lapply(program_levels, function(p) meta %>% transmute(observation_id, depmap_id,
  cell_line = tahoe_cell_name, compound, dose_uM = dose, split, prism_sensitization,
  program = p, score = score_matrix[, p])))
associations <- scores %>% group_by(compound, dose_uM, split, program) %>%
  group_modify(~pearson_test(.x$score, -.x$prism_sensitization)) %>% ungroup() %>%
  group_by(split) %>% mutate(bh_q = p.adjust(nominal_p, "BH")) %>% ungroup()
fwrite(scores, file.path(out, "tahoe_program_scores.tsv"), sep = "\t")
fwrite(associations, file.path(out, "program_prism_associations.tsv"), sep = "\t", na = "NA")
grDevices::cairo_pdf(file.path(out, "program_prism_scatterplots.pdf"), width = 11.6, height = 8.4, onefile = TRUE)
for (p in program_levels) {
  d <- scores %>% filter(program == p, is.finite(prism_sensitization))
  a <- associations %>% filter(program == p) %>% mutate(label = sprintf("%s n=%d r=%.2f p=%.2g", split, n, pearson_r, nominal_p),
    text_vjust = if_else(split == "Training", 1.3, 2.7))
  print(ggplot(d, aes(score, -prism_sensitization, color = split)) + geom_point(size = 1.7, alpha = .8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .7) +
    geom_text(data = a, aes(x = -Inf, y = Inf, label = label, vjust = text_vjust), hjust = -.02,
      size = 3.1, inherit.aes = FALSE) + facet_grid(compound ~ dose_uM, scales = "free") +
    scale_color_manual(values = c(Training = "#247A96", Holdout = "#B95B32")) +
    labs(title = paste0(p, " | Frozen program versus PRISM"), x = "Signed mean gene-expression displacement",
      y = "PRISM resistance (-sensitization)", subtitle = "All doses in uM; training associations reuse gene-selection data") +
    theme_classic(base_size = 12) + theme(aspect.ratio = 1, legend.position = "bottom"))
}
dev.off()
if ("--discovery-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no", status = 0)

## HELD-OUT EXPRESSION IS FIRST READ BELOW THIS LINE.
geo_dir <- file.path(root, "artifacts/22_raf_egf_bypass/gse127988")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
geo_urls <- c(
  GSE127988_family.soft.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127988/soft/GSE127988_family.soft.gz",
  GSE127988_RPM_Gerosa_et_al.csv.gz = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127988/suppl/GSE127988_RPM_Gerosa_et_al.csv.gz")
for (n in names(geo_urls)) if (!file.exists(file.path(geo_dir, n)))
  download.file(geo_urls[[n]], file.path(geo_dir, n), mode = "wb", quiet = TRUE)
fwrite(tibble(file = names(geo_urls), url = unname(geo_urls),
  md5 = unname(tools::md5sum(file.path(geo_dir, names(geo_urls))))),
  file.path(out, "external_input_manifest.tsv"), sep = "\t")
soft <- readLines(gzfile(file.path(geo_dir, names(geo_urls)[1])), warn = FALSE)
starts <- which(startsWith(soft, "^SAMPLE = "))
geo_meta <- bind_rows(lapply(seq_along(starts), function(i) {
  block <- soft[starts[i]:(if (i == length(starts)) length(soft) else starts[i+1] - 1L)]
  field <- function(prefix) {
    hit <- block[startsWith(block, prefix)]
    stopifnot(length(hit) == 1L)
    substring(hit, nchar(prefix) + 1L)
  }
  title <- field("!Sample_title = ")
  tibble(gsm = field("^SAMPLE = "), title = title,
    sample_number = as.integer(sub(".*\\(S([0-9]+)\\)$", "\\1", title)),
    vem_uM = as.numeric(field("!Sample_characteristics_ch1 = vemurafenib for 24 hrs (um): ")),
    cobi_uM = as.numeric(field("!Sample_characteristics_ch1 = cobimetinib for 24 hrs (um): ")),
    egf_h = as.numeric(field("!Sample_characteristics_ch1 = egf at 100 ng/ml (h): ")),
    biological_set = as.integer(field("!Sample_characteristics_ch1 = set id: ")))
})) %>% arrange(sample_number) %>% mutate(sample_id = paste0("S", sample_number),
  experiment = if_else(biological_set %in% 1:2, "Vemurafenib dose", "EGF time course"),
  arm = case_when(vem_uM == 0 ~ "DMSO", cobi_uM > 0 ~ "Vem + Cobi", TRUE ~ "Vem"))
stopifnot(nrow(geo_meta) == 96L, identical(geo_meta$sample_number, 1:96),
  setequal(geo_meta$biological_set, 1:4))
csv_path <- file.path(geo_dir, names(geo_urls)[2])
headers <- fread(csv_path, nrows = 6L, header = FALSE, colClasses = "character", check.names = FALSE)
individual <- which(as.character(headers[4, ]) == "Individual")
sample_numbers <- as.integer(as.character(headers[5, ..individual]))
stopifnot(identical(sample_numbers, 1:96),
  all(as.numeric(as.character(headers[1, ..individual])) == geo_meta$vem_uM),
  all(as.numeric(as.character(headers[2, ..individual])) == geo_meta$cobi_uM),
  all(as.numeric(as.character(headers[3, ..individual])) == geo_meta$egf_h))
geo_meta$geo_reported_set <- geo_meta$biological_set
geo_meta$csv_reported_set <- as.integer(as.character(headers[6, ..individual]))
discrepancies <- geo_meta %>% filter(geo_reported_set != csv_reported_set)
## The publisher matrix assigns S15/S16 DMSO wells to dose-experiment set 1;
## GEO SOFT assigns them to set 2. Use the matrix's explicit set labels, audit
## both, and require complete agreement for the held-out time-course experiment.
stopifnot(all(discrepancies$sample_number %in% c(15L, 16L)),
  all(discrepancies$geo_reported_set == 2L), all(discrepancies$csv_reported_set == 1L))
fwrite(discrepancies %>% mutate(resolution = "Use published expression-matrix set ID; dose-series DMSO only"),
  file.path(out, "gse127988_metadata_discrepancies.tsv"), sep = "\t")
geo_meta$biological_set <- geo_meta$csv_reported_set
rpm <- fread(csv_path, skip = 7L, header = FALSE, select = c(2L, individual),
  col.names = c("gene", geo_meta$sample_id), check.names = FALSE)
stopifnot(nrow(rpm) == 18261L, all(as.matrix(rpm[, -1L, with = FALSE]) >= 0))
mapping <- rpm[, .(ensembl_rows = .N), by = gene]
## Exact published gene symbols; sum RPM over duplicate symbol rows before log.
## No aliases, validation-driven gene selection, reorientation or reclustering.
rpm <- rpm[, lapply(.SD, sum), by = gene]
available_genes <- gene_order[gene_order %in% rpm$gene]
external_weights <- weights[available_genes, , drop = FALSE]
coverage_weight <- colSums(abs(external_weights))
external_weights <- sweep(external_weights, 2L, ifelse(coverage_weight > 0, coverage_weight, 1), "/")
gene_coverage <- membership %>% mutate(measured_gse127988 = gene %in% available_genes) %>%
  left_join(mapping, by = "gene")
program_coverage <- gene_coverage %>% group_by(program) %>% summarise(
  n_frozen = n(), n_measured = sum(measured_gse127988), fraction_measured = mean(measured_gse127988), .groups = "drop")
fwrite(gene_coverage, file.path(out, "gse127988_gene_coverage.tsv"), sep = "\t", na = "NA")
fwrite(program_coverage, file.path(out, "gse127988_program_coverage.tsv"), sep = "\t")
gene_expression <- log2(1 + as.matrix(rpm[match(available_genes, gene), -1L, with = FALSE]) / 100)
rownames(gene_expression) <- available_genes
colnames(gene_expression) <- geo_meta$sample_id
absolute_score <- t(gene_expression) %*% external_weights
absolute_score[, coverage_weight == 0] <- NA_real_
sample_scores <- bind_rows(lapply(program_levels, function(p) geo_meta %>% mutate(
  program = p, absolute_score = absolute_score[, p])))
## Set ID identifies the independent experiment day. All technical wells of a
## condition (including the four DMSO wells in each time-course set) are averaged.
bio <- sample_scores %>% group_by(program, experiment, biological_set, arm, vem_uM, cobi_uM, egf_h) %>%
  summarise(score = mean(absolute_score), n_technical_wells = n(), .groups = "drop")
baseline <- bio %>% filter(arm == "DMSO") %>% select(program, experiment, biological_set, dmso_score = score)
stopifnot(!anyDuplicated(baseline[c("program", "experiment", "biological_set")]))
bio <- bio %>% left_join(baseline, by = c("program", "experiment", "biological_set")) %>%
  mutate(dmso_displacement = score - dmso_score)
timecourse <- bio %>% filter(experiment == "EGF time course", arm != "DMSO")
time_zero <- timecourse %>% filter(egf_h == 0) %>% select(program, biological_set, arm, time_zero_score = score)
timecourse <- timecourse %>% left_join(time_zero, by = c("program", "biological_set", "arm")) %>%
  mutate(egf_displacement = score - time_zero_score)
stopifnot(nrow(timecourse) == 3L * 2L * 2L * 7L,
  all(timecourse$n_technical_wells == 2L), !anyNA(timecourse$time_zero_score))
fwrite(geo_meta, file.path(out, "gse127988_sample_metadata.tsv"), sep = "\t")
fwrite(sample_scores, file.path(out, "gse127988_sample_scores.tsv"), sep = "\t", na = "NA")
fwrite(bio, file.path(out, "gse127988_biological_scores.tsv"), sep = "\t", na = "NA")
fwrite(timecourse, file.path(out, "gse127988_timecourse_scores.tsv"), sep = "\t", na = "NA")

## Time-integrated contrasts, one value per biological day. No genes, wells or
## time points are counted as independent biological replicates.
auc <- function(time, value) {
  o <- order(time); time <- time[o]; value <- value[o]
  stopifnot(identical(as.numeric(time), c(0, .5, 1, 2, 3, 4, 8)))
  sum(diff(time) * (head(value, -1) + tail(value, -1)) / 2) / 8
}
area <- timecourse %>% group_by(program, biological_set, arm) %>% summarise(
  egf_auc_mean = auc(egf_h, egf_displacement), total_auc_mean = auc(egf_h, dmso_displacement), .groups = "drop")
v <- area %>% filter(arm == "Vem")
c <- area %>% filter(arm == "Vem + Cobi")
paired <- inner_join(v, c, by = c("program", "biological_set"), suffix = c("_vem", "_combo"))
contrasts <- bind_rows(
  transmute(v, program, biological_set, contrast = "EGF during vemurafenib", effect = egf_auc_mean),
  transmute(c, program, biological_set, contrast = "EGF during Vem + Cobi", effect = egf_auc_mean),
  transmute(paired, program, biological_set, contrast = "Cobi change in EGF response",
    effect = egf_auc_mean_combo - egf_auc_mean_vem),
  transmute(paired, program, biological_set, contrast = "Cobi change in total program level",
    effect = total_auc_mean_combo - total_auc_mean_vem))
test_mean <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 2L || sd(x) == 0) return(tibble(n_biological = n, mean_effect = mean(x),
    ci_low = NA_real_, ci_high = NA_real_, nominal_p = NA_real_, sign_test_p = NA_real_))
  a <- t.test(x, mu = 0)
  tibble(n_biological = n, mean_effect = mean(x), ci_low = a$conf.int[1], ci_high = a$conf.int[2],
    nominal_p = a$p.value, sign_test_p = if (any(x != 0)) binom.test(sum(x > 0), sum(x != 0), p = .5)$p.value else 1)
}
tests <- contrasts %>% group_by(program, contrast) %>% group_modify(~test_mean(.x$effect)) %>% ungroup() %>%
  mutate(bh_q = p.adjust(nominal_p, method = "BH")) %>% left_join(predictions, by = "program") %>%
  left_join(program_coverage, by = "program") %>%
  mutate(expected_direction = case_when(contrast == "EGF during vemurafenib" ~ expected_egf_direction,
    contrast == "Cobi change in total program level" ~ expected_cobimetinib_direction, TRUE ~ "exploratory"),
    direction_agrees = case_when(expected_direction == "increase" ~ mean_effect > 0,
      expected_direction == "decrease" ~ mean_effect < 0, TRUE ~ NA),
    interpretation = "Exploratory biological-day t test (n=2); exact sign test also reported; no pseudoreplication")
fwrite(contrasts, file.path(out, "gse127988_day_level_contrasts.tsv"), sep = "\t")
fwrite(tests, file.path(out, "gse127988_validation_tests.tsv"), sep = "\t", na = "NA")

arm_colors <- c(Vem = "#247A96", `Vem + Cobi` = "#B95B32")
for (measure in c("egf_displacement", "dmso_displacement")) {
  d <- timecourse %>% mutate(displacement = .data[[measure]],
    panel = paste0(program, " | ", coherence$candidate_role[match(program, coherence$program)],
      " | ", program_coverage$n_measured[match(program, program_coverage$program)], "/",
      program_coverage$n_frozen[match(program, program_coverage$program)], " genes matched"))
  means <- d %>% group_by(panel, arm, egf_h) %>% summarise(displacement = mean(displacement), .groups = "drop")
  p <- ggplot(d, aes(egf_h, displacement, color = arm)) + geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = .4) +
    geom_line(aes(group = interaction(arm, biological_set)), linewidth = .5, alpha = .4) +
    geom_point(aes(shape = factor(biological_set)), size = 2.4) +
    geom_line(data = means, aes(group = arm), linewidth = .9) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") + scale_color_manual(values = arm_colors) +
    scale_x_continuous(breaks = c(0, .5, 1, 2, 3, 4, 8)) +
    labs(x = "Hours after EGF addition", y = "Frozen program displacement", color = NULL, shape = "Biological set",
      title = if (measure == "egf_displacement") "EGF response under RAF / MEK inhibition" else "Program state relative to vehicle",
      subtitle = if (measure == "egf_displacement") "Relative to the same drug arm at time zero" else "Relative to same-day DMSO",
      caption = "Points: two biological days, each averaging duplicate wells. Thick lines: day means.\nP1/P2 face PRISM resistance; P3 faces CytoSig EGF induction. Signs fixed in discovery.\nVemurafenib 1 uM; cobimetinib 1 uM; EGF 100 ng/mL.") +
    theme_classic(base_size = 12) + theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, face = "bold"),
      legend.position = "bottom", plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 9))
  ggsave(file.path(out, paste0("gse127988_", measure, ".pdf")), p, width = 8.2, height = 9.2, device = cairo_pdf)
}

## Validation heatmaps keep discovery gene order; chronological columns are
## essential for the perturbation test. No validation-driven sample clustering.
conditions <- geo_meta %>% filter(experiment == "EGF time course", arm != "DMSO") %>%
  distinct(arm, egf_h) %>% arrange(factor(arm, levels = c("Vem", "Vem + Cobi")), egf_h) %>%
  mutate(condition_id = paste0(arm, " | ", egf_h, " h"))
gene_displacement <- matrix(NA_real_, length(gene_order), nrow(conditions), dimnames = list(gene_order, conditions$condition_id))
for (j in seq_len(nrow(conditions))) {
  day_effects <- vapply(3:4, function(day) {
    target <- geo_meta$sample_id[geo_meta$biological_set == day & geo_meta$arm == conditions$arm[j] & geo_meta$egf_h == conditions$egf_h[j]]
    reference <- geo_meta$sample_id[geo_meta$biological_set == day & geo_meta$arm == conditions$arm[j] & geo_meta$egf_h == 0]
    rowMeans(gene_expression[, target, drop = FALSE]) - rowMeans(gene_expression[, reference, drop = FALSE])
  }, numeric(length(available_genes)))
  gene_displacement[available_genes, j] <- rowMeans(day_effects)
}
write_matrix(gene_displacement, "gse127988_gene_displacements.tsv.gz")
ca <- data.frame(Arm = conditions$arm, row.names = conditions$condition_id)
for (labelled in c(FALSE, TRUE)) save_heatmap(gene_displacement,
  paste0("gse127988_gene_displacement", if (labelled) "_labelled" else "", ".pdf"),
  "GSE127988 | EGF displacement | Discovery gene order", c(-1, 0, 1), "Log2 change", ca,
  labelled = labelled, column_labels = paste0(conditions$egf_h, " h"), width = 10.5)
program_displacement <- t(gene_displacement[available_genes, , drop = FALSE]) %*% external_weights
program_displacement <- t(program_displacement)
pa <- data.frame(Program = factor(program_levels, levels = program_levels), row.names = program_levels)
save_heatmap(program_displacement, "gse127988_program_heatmap.pdf", "Signed frozen programs | EGF response in A375",
  c(-max(abs(program_displacement), na.rm = TRUE), 0, max(abs(program_displacement), na.rm = TRUE)),
  "Program change", ca, row_annotation = pa, labelled = TRUE,
  column_labels = paste0(conditions$egf_h, " h"), width = 9.5, height = 3.4,
  row_labels = paste0(program_levels, " (", program_coverage$n_measured, "/", program_coverage$n_frozen, " genes)"))
stopifnot(identical(unname(tools::md5sum(freeze_paths)), unname(freeze_hashes)))
writeLines(c("# Frozen RAF–EGF bypass pilot", "",
  "Run: `Rscript /work/tahoe_prism_v2/scripts/pilot_raf_egf_bypass.R`. `--discovery-only` stops before loading external expression.", "",
  "Discovery: reuse processed Tahoe pseudobulk displacements, original normalization, original PRISM matches and original 26/7 line split. Use all available doses. A drug needs at least eight matched training lines at one dose; each gene test also needs eight at its dose. Cobimetinib is optional at this overlap check, not selected by its results.",
  "First intersect the FULL Tahoe and CytoSig gene universes. No old 3,000-gene, activity, or response-strength filter. Within each drug/dose select nominal two-sided Pearson p<0.05 genes, then cap at 100 per drug ranked by maximum absolute significant correlation (ties: minimum p, gene name). Training associations reuse selection data; p values are nominal screening evidence, not selection-adjusted inference.",
  "CytoSig: exact EGF profiles only, excluding the external accession. To represent EGF-only biology add up to 100 genes ranked by absolute median effect across study medians, with at least three measured studies. This EGF ranking is not a significance claim. Keep all individual EGF profiles for gene correlations. The profile count includes correlated time points and doses, so coherence is descriptive.",
  "Cluster the union once: arithmetic mean of Tahoe training and CytoSig EGF Pearson matrices, distance 1-r, average linkage, exactly three clusters. Only genes with unestimable correlations are removed; the audit lists them and capped lists are not replenished. Ranks 2-8 diagnose silhouette without selecting rank. All source heatmaps preserve the identical gene order.",
  "Candidate roles follow positive mean silhouette of each fixed cluster in the two sources: both = shared downstream candidate; EGF only = receptor candidate; Tahoe only = drug candidate; neither = unresolved. These labels are descriptive evidence, not proof of receptor exclusivity, survival causality or a forced one-of-each Y. Missing or repeated roles are a negative/incomplete architectural result.",
  "Freeze equal gene weights within each program. Receptor-only candidates face EGF induction; other programs face the median training vemurafenib PRISM-resistance correlation across doses. A single sign applies to every gene in a program. The shared-survival interpretation also needs positive EGF displacement and independent functional support. Program genes, signs, weights, role annotations and predictions are written and hashed BEFORE loading external expression, then checked unchanged at the end.",
  "GSE127988 is external molecular validation. Read individual RPM columns only and verify all 96 sample IDs, doses, times and set IDs against GEO SOFT metadata. Sum RPM for duplicate exact gene symbols, use log2(1+RPM/100) to match Tahoe's 10,000-count normalization units, and average fixed weighted genes. Available genes retain their frozen relative weights, renormalized once by measured coverage; no replacement genes, external SD fitting or sign changes. Gene/program coverage is reported explicitly.",
  "Coverage is based on exact symbol matches. Unmatched genes may include historical symbol differences as well as genes absent from the published matrix; no unverified alias mapping is applied. Gene heatmaps show raw expression displacement, while program plots apply the frozen program sign. A positive program score need not mean increased raw gene expression.",
  "GEO sets 1/2 are the dose experiment; sets 3/4 are two biological days of the EGF time course. Average duplicate wells within each day/condition before inference. Pool four DMSO wells per time-course day as the same-day vehicle reference. Score all 96 samples; infer EGF effects from the time-course experiment only.",
  "Metadata discrepancy: GEO labels dose-series DMSO samples S15/S16 as set 2, while the expression matrix explicitly labels them set 1. Use matrix set IDs and preserve both values in the discrepancy audit. All EGF time-course set IDs agree between the sources.",
  "Primary EGF contrast is the trapezoidal 0–8 h mean score relative to that arm's time zero. Cobimetinib's total effect is the paired difference between Vem+Cobi and Vem time-integrated scores; also report its difference in EGF-specific induction separately. Test at the biological-day level (n=2), reporting exploratory two-sided t tests, 95% intervals, BH across the 12 contrasts, and exact two-sided sign tests. With two days, sign tests cannot establish conventional significance. Time points and genes are not replicates.",
  "No untreated time-matched EGF-free time course is available; the EGF contrast also contains post-time-zero temporal change. There is no EGF-only arm in drug-naive cells. These are molecular perturbation contrasts, not a transcriptome-derived viability assay. No claim of complete cobimetinib blockade is required: receptor-driven signaling can persist under RAF/MEK inhibition.",
  "Evidence roles: PRISM provides scalable functional association for selection; Gerosa et al. independently demonstrates receptor-driven ERK pulses and persistence; GSE127988 tests frozen molecular displacement under joint perturbations. https://pubmed.ncbi.nlm.nih.gov/33113355/ ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988", "",
  paste0("Observed discovery: ", length(gene_order), " frozen genes; ", nrow(cyto), " EGF profiles; ", length(studies), " studies."),
  paste0(coherence$program, ": ", coherence$n_genes, " genes; ", coherence$candidate_role),
  "", "Measured validation coverage (partial projections, not complete-program measurements):",
  paste0(program_coverage$program, ": ", program_coverage$n_measured, "/", program_coverage$n_frozen,
    " genes (", round(100 * program_coverage$fraction_measured, 1), "%)."),
  "The source-coherence labels do not imply that a program predominantly contains EGF-selected genes. program_gene_origin_summary.tsv reports selection provenance separately. No clean one-of-each three-branch architecture was recovered in this fit.",
  "", "Key plots: gse127988_egf_displacement.pdf, gse127988_dmso_displacement.pdf, gse127988_program_heatmap.pdf, joint_correlation.pdf, and program_prism_scatterplots.pdf."),
  file.path(out, "README.md"))
print(program_coverage)
print(select(tests, program, contrast, mean_effect, nominal_p, bh_q, expected_direction, direction_agrees))
message("Completed frozen RAF/EGF discovery and GSE127988 molecular validation.")
