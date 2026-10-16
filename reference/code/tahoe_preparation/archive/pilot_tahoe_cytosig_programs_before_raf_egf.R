#!/usr/bin/env Rscript

## Replacement analysis requested by the user: all genes, all doses of the
## three drugs, and all local CytoSig profiles. No prior activity/gene gate.
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
out <- file.path(root, "results/20_tahoe_cytosig_programs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
drugs <- c("erlotinib", "crizotinib", "trametinib")
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
  setequal(unique(meta$compound), drugs))
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

## All local CytoSig treatment/control profiles, with no cytokine restriction.
## fread's decompression scratch file is temporary; the dataset is not copied
## into the project and no new CytoSig data are downloaded.
cyto_meta <- as_tibble(fread(paths[["cytosig_metadata"]], skip = 1L, header = FALSE,
  col.names = c("profile_id", "treatment", "condition", "duration", "dose")))
con <- gzfile(paths[["cytosig_effects"]], "rt")
profile_header <- strsplit(readLines(con, 1L), "\t", fixed = TRUE)[[1]]
close(con)
stopifnot(!anyDuplicated(cyto_meta$profile_id), !anyDuplicated(profile_header))
cyto_meta <- tibble(profile_id = profile_header) %>%
  left_join(cyto_meta, by = "profile_id") %>% mutate(metadata_available = !is.na(treatment))
cyto_dt <- suppressWarnings(fread(paths[["cytosig_effects"]], header = TRUE,
  col.names = c("gene", profile_header), check.names = FALSE, showProgress = FALSE))
stopifnot(!anyDuplicated(cyto_dt$gene), ncol(cyto_dt) == length(profile_header) + 1L)
cyto_genes <- cyto_dt$gene

## First define the full cross-source measured-gene universe. Then apply the
## requested training association rule and 100-per-drug cap, ranked by |r|.
shared_universe <- intersect(genes$feature, cyto_genes)
stopifnot(length(shared_universe) > 0)
message("Shared measured-gene universe: ", length(shared_universe))
## Two-sided Pearson tests, separately within each drug/dose and training set.
## This is the only statistical gene-selection rule: nominal p < 0.05.
gene_tests <- list()
for (drug in drugs) for (dose in sort(unique(meta$dose[meta$compound == drug]))) {
  rows <- which(meta$compound == drug & meta$dose == dose & meta$split == "Training" &
    is.finite(meta$prism_sensitization))
  y <- meta$prism_sensitization[rows]
  x <- effects[rows, shared_universe, drop = FALSE]
  r <- suppressWarnings(as.numeric(cor(x, y)))
  n <- length(rows)
  if (n < 3L || sd(y) == 0) r[] <- NA_real_
  statistic <- abs(r) * sqrt((n - 2) / pmax(0, 1 - r ^ 2))
  p <- 2 * pt(statistic, df = n - 2, lower.tail = FALSE)
  gene_tests[[length(gene_tests) + 1L]] <- tibble(compound = drug, dose_uM = dose,
    gene = colnames(x), n_training_lines = n, pearson_r = r, nominal_p = p,
    selected_at_dose = is.finite(p) & p < .05,
    status = if_else(is.finite(r), "tested", "untestable_zero_variance_or_insufficient_pairs"))
}
gene_tests <- bind_rows(gene_tests)
fwrite(gene_tests, file.path(out, "gene_prism_training_tests.tsv.gz"), sep = "\t", na = "NA")
candidate_genes <- gene_tests %>% filter(selected_at_dose) %>% group_by(compound, gene) %>%
  summarise(best_training_p = min(nominal_p), strongest_training_abs_r = max(abs(pearson_r)),
    significant_doses_uM = paste(sort(dose_uM), collapse = ";"), .groups = "drop") %>%
  arrange(compound, desc(strongest_training_abs_r), best_training_p, gene)

candidate_genes <- candidate_genes %>% mutate(measured_in_cytosig = gene %in% cyto_genes)
fwrite(candidate_genes, file.path(out, "significant_gene_candidates.tsv"), sep = "\t")
candidate_genes <- candidate_genes %>% filter(measured_in_cytosig) %>%
  group_by(compound) %>% mutate(selection_rank = row_number()) %>% slice_head(n = 100L) %>% ungroup()
stopifnot(all((candidate_genes %>% dplyr::count(compound))$n <= 100L))
print(candidate_genes %>% dplyr::count(compound, name = "genes_after_training_p_rank_cap"))
cyto_dt <- cyto_dt[gene %chin% candidate_genes$gene]
cyto <- t(as.matrix(cyto_dt[, -1L, with = FALSE]))
colnames(cyto) <- cyto_dt$gene; rownames(cyto) <- profile_header
cyto[!is.finite(cyto)] <- NA_real_
rm(cyto_dt); invisible(gc())
fwrite(cyto_meta, file.path(out, "cytosig_profiles.tsv"), sep = "\t", na = "NA")
fwrite(meta %>% group_by(compound, dose, split) %>% summarise(
  molecular_lines = n_distinct(depmap_id), prism_lines = sum(is.finite(prism_sensitization)),
  cell_line_ids = paste(depmap_id, collapse = ";"), .groups = "drop") %>%
  mutate(full_tahoe_genes = nrow(genes), cytosig_genes = length(cyto_genes),
    cytosig_profiles = nrow(cyto), shared_available_genes = sum(genes$feature %in% cyto_genes)),
  file.path(out, "input_coverage.tsv"), sep = "\t")
manifest <- tibble(object = names(paths), path = unname(paths), mapping = c(
  "depmap_id; compound=drug; dose=dose uM; shard_path; shard_row zero-based; available=1",
  "ALL feature values in gene_index order, no selected-gene or activity filter",
  "observation_id; treated_total; reference_total; treated_correction_factor; contributing_plates",
  "NPY float32 contributing-plate DMSO pseudobulk references; columns=all_genes.gene_index",
  "array_row zero-based; depmap_id; plate",
  "existing chain_observed=1 matches; response=prism_sensitization (higher=more sensitivity)",
  "depmap_id; outer_fold=1 holdout, otherwise training; original 26/7 split",
  "unnamed first column=profile_id; Treatment; Condition; Duration; Dose; ALL treatments",
  "unnamed first column=gene; all profile columns; signed processed treatment/control log2 changes"))
manifest <- bind_rows(manifest, tibble(object = "processed_tahoe_shard",
  path = normalizePath(file.path(cache, unique(meta$shard_path)), mustWork = TRUE),
  mapping = "row=shard_row; columns=all_genes.gene_index; processed arithmetic treatment-control displacement"))
fwrite(manifest, file.path(out, "input_manifest.tsv"), sep = "\t")
message("Loaded all ", nrow(cyto), " CytoSig profiles; ", ncol(cyto), " of the capped candidate genes are measured.")

## Ward.D2 hierarchical clustering of signed-correlation distances. The user fixed
## the final rank at eight. Silhouette remains a diagnostic over valid ranks.
fit_clusters <- function(r) {
  stopifnot(nrow(r) >= 3L, all(is.finite(r)), max(abs(r - t(r))) < 1e-10)
  r[r > 1] <- 1; r[r < -1] <- -1; diag(r) <- 1
  distance <- as.dist(1 - r)
  tree <- hclust(distance, method = "ward.D2")
  ranks <- 2:(nrow(r) - 1L)
  curve <- bind_rows(lapply(ranks, function(k) {
    groups <- cutree(tree, k)
    s <- cluster::silhouette(groups, distance)
    tibble(rank = k, mean_silhouette = mean(s[, "sil_width"]),
      minimum_cluster_size = min(table(groups)))
  }))
  stopifnot(nrow(r) > 8L)
  k <- 8L
  raw <- cutree(tree, k)
  cluster_order <- unique(raw[tree$order])
  groups <- setNames(match(raw, cluster_order), names(raw))
  list(tree = tree, groups = groups, k = k, curve = mutate(curve, selected = rank == k,
    selection_method = "user_fixed_rank_8", silhouette_maximum = mean_silhouette == max(mean_silhouette)),
    order = tree$order)
}
test_correlation <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(tibble(n = n, pearson_r = NA_real_, nominal_p = NA_real_))
  r <- cor(x[ok], y[ok])
  tibble(n = n, pearson_r = r, nominal_p = 2 * pt(abs(r) * sqrt((n - 2) /
    pmax(0, 1 - r ^ 2)), df = n - 2, lower.tail = FALSE))
}

## ONE union of the per-drug selected genes; ONE joint clustering.
selected_genes <- sort(unique(candidate_genes$gene))
selection_audit <- tibble(gene = selected_genes,
  selected_by = vapply(selected_genes, function(g) paste(sort(candidate_genes$compound[candidate_genes$gene == g]), collapse = ";"), character(1)),
  used_for_clustering = TRUE, exclusion_reason = NA_character_)
cyto_ok <- vapply(selected_genes, function(g) {
  v <- cyto[, g]; v <- v[is.finite(v)]
  length(v) >= 3 && sd(v) > 0
}, logical(1))
selection_audit$used_for_clustering[!cyto_ok] <- FALSE
selection_audit$exclusion_reason[!cyto_ok] <- "CytoSig correlation not estimable: fewer than three values or zero variance"
selected_genes <- selected_genes[cyto_ok]
training <- which(meta$split == "Training")
r_tahoe <- cor(effects[training, selected_genes, drop = FALSE])
r_cyto <- cor(cyto[, selected_genes, drop = FALSE], use = "pairwise.complete.obs")
cyto_pair_n <- crossprod(is.finite(cyto[, selected_genes, drop = FALSE]))
r_cyto[cyto_pair_n < 3L] <- NA_real_
invalid <- !is.finite(r_tahoe) | !is.finite(r_cyto)
bad <- which(invalid & upper.tri(invalid), arr.ind = TRUE)
fwrite(tibble(gene = selected_genes[bad[, 1]], other_gene = selected_genes[bad[, 2]],
  cytosig_joint_profiles = cyto_pair_n[bad]), file.path(out, "unestimable_gene_pairs.tsv"), sep = "\t")
## No correlation is imputed. Resolve only mathematically untestable pairs,
## using missingness counts, then observed-profile counts and gene name.
while (any(invalid)) {
  problem <- tibble(gene = selected_genes, bad_pairs = rowSums(invalid),
    observed = colSums(is.finite(cyto[, selected_genes, drop = FALSE]))) %>%
    arrange(desc(bad_pairs), observed, gene)
  remove <- problem$gene[1]
  selection_audit$used_for_clustering[selection_audit$gene == remove] <- FALSE
  selection_audit$exclusion_reason[selection_audit$gene == remove] <-
    paste0("Unestimable cross-source gene correlations: ", problem$bad_pairs[1], " pairs at removal")
  keep <- selected_genes != remove
  selected_genes <- selected_genes[keep]
  r_tahoe <- r_tahoe[keep, keep, drop = FALSE]
  r_cyto <- r_cyto[keep, keep, drop = FALSE]
  cyto_pair_n <- cyto_pair_n[keep, keep, drop = FALSE]
  invalid <- !is.finite(r_tahoe) | !is.finite(r_cyto)
}
stopifnot(length(selected_genes) >= 3L, length(selected_genes) <= 300L)
matrices <- list(Joint = (r_tahoe + r_cyto) / 2, Tahoe = r_tahoe, CytoSig = r_cyto)
fit <- fit_clusters(matrices$Joint)
gene_order <- selected_genes[fit$order]
program_levels <- paste0("P", seq_len(fit$k))
programs <- setNames(paste0("P", fit$groups[gene_order]), gene_order)
membership <- tibble(gene = gene_order, program = unname(programs), plotting_order = seq_along(gene_order)) %>%
  left_join(select(selection_audit, gene, selected_by), by = "gene")
stopifnot(!anyNA(membership$program), identical(names(programs), gene_order))
palette <- setNames(grDevices::hcl.colors(fit$k, "Dark 3"), program_levels)
split <- factor(programs, levels = program_levels)
annotation <- data.frame(Program = split, row.names = gene_order)
fwrite(candidate_genes, file.path(out, "selected_genes_by_drug.tsv"), sep = "\t")
fwrite(selection_audit, file.path(out, "gene_selection_audit.tsv"), sep = "\t", na = "NA")
fwrite(membership, file.path(out, "joint_program_membership.tsv"), sep = "\t")
fwrite(fit$curve, file.path(out, "rank_silhouette.tsv"), sep = "\t")
p_rank <- ggplot(fit$curve, aes(rank, mean_silhouette)) + geom_line(linewidth = .8, color = "#255F85") +
  geom_point(size = 1.3, color = "#255F85") +
  geom_point(data = filter(fit$curve, selected), size = 3.7, shape = 21, fill = "#E9A23B", stroke = .9) +
  labs(x = "Number of joint programs", y = "Mean silhouette width",
    title = "One joint clustering of the combined drug-associated genes",
    subtitle = paste0(length(gene_order), " genes | rank fixed manually at ", fit$k),
    caption = "Distance = 1 - (Tahoe r + CytoSig r)/2. Eight programs chosen by the user; silhouette is diagnostic.") +
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 10), plot.margin = margin(9, 9, 9, 9))
ggsave(file.path(out, "joint_rank_silhouette.pdf"), p_rank, width = 8.5, height = 5.5, device = cairo_pdf)

## All correlation figures use EXACTLY the same genes, order and program split.
## A compact overview is followed by a fully labelled page when many genes are
## present; both pages use the same matrix, with no representative-gene subset.
for (source_name in names(matrices)) {
  r <- matrices[[source_name]][gene_order, gene_order, drop = FALSE]
  stopifnot(identical(rownames(r), gene_order), identical(colnames(r), gene_order))
  fwrite(data.table(gene = gene_order, r, check.names = FALSE),
    file.path(out, paste0(tolower(source_name), "_gene_correlations.tsv.gz")), sep = "\t")
  for (labelled in c(FALSE, TRUE)) {
    suffix <- if (labelled) "_labelled" else ""
    size <- if (labelled) max(9, length(gene_order) * .16 + 3.4) else 9
    file <- file.path(out, paste0(tolower(source_name), "_correlation_heatmap", suffix, ".pdf"))
    grDevices::cairo_pdf(file, width = size + .5, height = size, onefile = TRUE)
    ht <- draw_heatmap(r, col_anno_df = annotation, col_anno_cols = "Program",
      row_anno_df = annotation, row_anno_cols = "Program", palette_overrides = list(Program = palette),
      ann_top_args = list(show_legend = TRUE, annotation_legend_param = list(Program = list(
        title_gp = grid::gpar(fontsize = 12, fontface = "bold"), labels_gp = grid::gpar(fontsize = 12)))),
      ann_row_args = list(show_legend = FALSE),
      scale_rows = FALSE, value_range = c(-1, 0, 1), value_colors = c("#2166AC", "white", "#B2182B"),
      hm_name = "Pearson r", show_row_names = labelled, show_column_names = labelled,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_labels = paste0(gene_order, " [", programs, "]"),
      row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 11),
      row_title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
      column_title_gp = grid::gpar(fontsize = 17, fontface = "bold"),
      row_gap = grid::unit(2, "mm"), column_gap = grid::unit(2, "mm"),
      title = paste0(source_name, " | Combined drug-associated genes\n",
        length(gene_order), " genes | ", fit$k, " joint programs | fixed gene order"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 11), at = c(-1, -.5, 0, .5, 1)))
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
      merge_legends = TRUE, padding = grid::unit(c(9, 7, 7, 7), "mm"))
    ComplexHeatmap::decorate_heatmap_body("Pearson r", {
      for (p in program_levels) {
        block <- which(programs == p)
        position <- (min(block) + max(block) - 1) / (2 * length(gene_order))
        grid::grid.rect(x = position, y = 1 - position, width = length(block) / length(gene_order),
          height = length(block) / length(gene_order), gp = grid::gpar(fill = NA, col = "#71717A", lwd = .35))
      }
    })
    grid::grid.text("Tahoe: all three drugs and doses, training observations only. CytoSig: all 2,056 profiles.",
      x = .5, y = grid::unit(3, "mm"), gp = grid::gpar(fontsize = 10))
    dev.off()
  }
}

## Genes x samples: gene order is frozen. ONLY sample columns are reordered.
## Row z-scores are for display; original signed responses and NA values are
## preserved in the exported matrices. Dominant positive program response
## groups samples into blocks, then samples are ordered within those blocks.
sample_orders <- list()
for (source_name in c("Tahoe", "CytoSig")) {
  expression <- if (source_name == "Tahoe") t(effects[training, gene_order, drop = FALSE]) else
    t(cyto[, gene_order, drop = FALSE])
  center <- rowMeans(expression, na.rm = TRUE)
  scale <- apply(expression, 1L, sd, na.rm = TRUE)
  z <- sweep(sweep(expression, 1L, center, "-"), 1L, scale, "/")
  stopifnot(identical(rownames(z), gene_order))
  activity_score <- vapply(program_levels, function(p) colMeans(z[programs == p, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(z)))
  rownames(activity_score) <- colnames(z)
  observed_score <- activity_score
  activity_score[!is.finite(activity_score)] <- -Inf
  has_score <- rowSums(is.finite(observed_score)) > 0L
  dominant <- max.col(activity_score, ties.method = "first")
  sample_group <- ifelse(has_score, program_levels[dominant], "Unmeasured")
  column_order <- unlist(lapply(c(program_levels, "Unmeasured"), function(p) {
    ids <- which(sample_group == p)
    if (length(ids) < 3L) return(ids)
    ## Euclidean distance on program-response profiles. For samples lacking
    ## a program score, order by the observed dominant score and ID instead.
    good <- apply(observed_score[ids, , drop = FALSE], 1L, function(v) all(is.finite(v)) && sd(v) > 0)
    if (all(good)) return(ids[hclust(dist(observed_score[ids, , drop = FALSE]), method = "ward.D2")$order])
    strength <- apply(activity_score[ids, , drop = FALSE], 1L, max)
    ids[order(-strength, colnames(z)[ids])]
  }), use.names = FALSE)
  stopifnot(setequal(column_order, seq_len(ncol(z))))
  z <- z[, column_order, drop = FALSE]
  expression <- expression[, column_order, drop = FALSE]
  col_split <- factor(sample_group[column_order], levels = c(program_levels, "Unmeasured"))
  col_annotation <- data.frame(Program = col_split, row.names = colnames(z))
  sample_orders[[source_name]] <- tibble(source = source_name, sample_id = colnames(z),
    plotting_order = seq_len(ncol(z)), dominant_program = as.character(col_split))
  if (source_name == "Tahoe") {
    m <- meta[match(colnames(z), meta$observation_id), ]
    col_annotation$Drug <- factor(m$compound, levels = drugs)
    col_annotation$Dose <- factor(m$dose)
    sample_orders[[source_name]] <- sample_orders[[source_name]] %>%
      left_join(select(meta, sample_id = observation_id, compound, dose, depmap_id, tahoe_cell_name), by = "sample_id")
  } else sample_orders[[source_name]] <- sample_orders[[source_name]] %>%
    left_join(rename(cyto_meta, cytosig_dose = dose), by = c("sample_id" = "profile_id"))
  fwrite(data.table(gene = gene_order, expression, check.names = FALSE),
    file.path(out, paste0(tolower(source_name), "_gene_sample_effects.tsv.gz")), sep = "\t", na = "NA")
  fwrite(tibble(source = source_name, gene = gene_order, display_center = center, display_scale = scale),
    file.path(out, paste0(tolower(source_name), "_heatmap_scaling.tsv")), sep = "\t")
  for (labelled in c(FALSE, TRUE)) {
    height <- if (labelled) max(8, length(gene_order) * .16 + 2.8) else 9
    suffix <- if (labelled) "_labelled" else ""
    grDevices::cairo_pdf(file.path(out, paste0(tolower(source_name), "_gene_sample_heatmap", suffix, ".pdf")),
      width = 14, height = height, onefile = TRUE)
    ht <- draw_heatmap(z, col_anno_df = col_annotation, col_anno_cols = names(col_annotation),
      row_anno_df = annotation, row_anno_cols = "Program",
      ann_top_args = list(show_legend = TRUE, annotation_legend_param = list(Program = list(
        title_gp = grid::gpar(fontsize = 12, fontface = "bold"), labels_gp = grid::gpar(fontsize = 12)))),
      ann_row_args = list(show_legend = FALSE),
      palette_overrides = list(Program = c(palette, Unmeasured = "#AAAAAA")),
      scale_rows = FALSE, value_range = c(-2, 0, 2), value_colors = c("#2166AC", "white", "#B2182B"),
      hm_name = "Gene z-score", show_row_names = labelled, show_column_names = FALSE,
      cluster_rows = FALSE, cluster_columns = FALSE, row_labels = paste0(gene_order, " [", programs, "]"),
      row_names_gp = grid::gpar(fontsize = 11), row_title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
      column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      row_gap = grid::unit(2, "mm"), column_gap = grid::unit(2, "mm"), na_col = "#B7B7B7",
      title = paste0(source_name, " | Genes x samples\n", length(gene_order), " genes | ", ncol(z),
        " samples | fixed genes; samples grouped by program"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 11), at = c(-2, -1, 0, 1, 2)))
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
      merge_legends = TRUE, padding = grid::unit(c(10, 7, 7, 7), "mm"))
    ComplexHeatmap::decorate_heatmap_body("Gene z-score", {
      for (boundary in cumsum(as.integer(table(split)))[-length(program_levels)])
        grid::grid.lines(x = c(0, 1), y = rep(1 - boundary / length(gene_order), 2), gp = grid::gpar(col = "#FFFFFF", lwd = .4))
      column_counts <- as.integer(table(col_split)); column_counts <- column_counts[column_counts > 0]
      for (boundary in head(cumsum(column_counts), -1))
        grid::grid.lines(x = rep(boundary / ncol(z), 2), y = c(0, 1), gp = grid::gpar(col = "#FFFFFF", lwd = .4))
    })
    grid::grid.text("Rows: fixed joint-program gene order. Columns: individual response profiles. Gray: missing. Color saturates at +/-2 SD.",
      x = .5, y = grid::unit(3, "mm"), gp = grid::gpar(fontsize = 10))
    dev.off()
  }
}
fwrite(bind_rows(sample_orders), file.path(out, "sample_plotting_order.tsv"), sep = "\t", na = "NA")

## Describe source support WITHOUT reclustering either source. Positive mean
## silhouette means the fixed program is more coherent internally than with
## its nearest other fixed program, in that source's correlation matrix.
coherence <- tibble(program = program_levels,
  n_genes = vapply(program_levels, function(p) sum(programs == p), integer(1)))
for (s in c("Tahoe", "CytoSig")) {
  distance_matrix <- 1 - matrices[[s]][gene_order, gene_order]
  distance_matrix[distance_matrix < 0] <- 0
  silhouette <- cluster::silhouette(as.integer(split), as.dist(distance_matrix))
  coherence[[paste0(tolower(s), "_mean_silhouette")]] <- vapply(program_levels,
    function(p) mean(silhouette[programs == p, "sil_width"]), numeric(1))
}
coherence <- coherence %>% mutate(sharing = case_when(
  tahoe_mean_silhouette > 0 & cytosig_mean_silhouette > 0 ~ "Coherent in both sources",
  tahoe_mean_silhouette > 0 ~ "Coherent in Tahoe only",
  cytosig_mean_silhouette > 0 ~ "Coherent in CytoSig only", TRUE ~ "Not coherent in either source"))
fwrite(coherence, file.path(out, "program_source_coherence.tsv"), sep = "\t")
tree_rows <- coherence %>% arrange(sharing, match(program, program_levels)) %>% mutate(y = rev(seq_len(n())))
tree_groups <- tree_rows %>% group_by(sharing) %>% summarise(low = min(y), high = max(y), y = mean(y), .groups = "drop")
segments <- bind_rows(tibble(x = 0, xend = 0, y = min(tree_groups$y), yend = max(tree_groups$y)),
  transmute(tree_groups, x = 0, xend = .8, y, yend = y),
  transmute(tree_groups, x = .8, xend = .8, y = low, yend = high),
  transmute(tree_rows, x = .8, xend = 1.7, y, yend = y))
p_tree <- ggplot() + geom_segment(data = segments, aes(x, y, xend = xend, yend = yend), linewidth = .7, color = "#586574") +
  geom_point(data = tree_rows, aes(1.7, y, color = program), size = 3) +
  geom_text(data = tree_groups, aes(.85, high + .34, label = sharing), hjust = 0, size = 3.4, fontface = "bold") +
  geom_text(data = tree_rows, aes(1.8, y, label = paste0(program, " (", n_genes, " genes)")), hjust = 0, size = 4.2) +
  geom_text(data = tree_rows, aes(3.05, y, label = sprintf("%.2f       %.2f", tahoe_mean_silhouette, cytosig_mean_silhouette)), hjust = 0, size = 3.8) +
  annotate("text", x = 3.05, y = max(tree_rows$y) + .7, label = "Mean silhouette\nTahoe   CytoSig", hjust = 0, size = 3.6) +
  scale_color_manual(values = palette, guide = "none") + xlim(-.1, 4.2) + ylim(0, max(tree_rows$y) + 1.1) +
  labs(title = "Joint programs | Descriptive source-sharing tree",
    subtitle = "The same fixed programs are evaluated in each source; no source is reclustered",
    caption = "Coherent = positive mean silhouette for the fixed program. Groups remain unresolved.\nThis describes co-response structure, not evolutionary ancestry or causal direction.") +
  theme_void(base_size = 13) + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 11),
    plot.caption = element_text(size = 10, hjust = 0), plot.margin = margin(10, 10, 10, 10))
ggsave(file.path(out, "joint_program_sharing_tree.pdf"), p_tree, width = 11, height = max(4.6, fit$k * .48 + 2.4), device = cairo_pdf)

## Score these same programs for EVERY drug/dose, with training-only gene
## scaling. No PRISM-based gene sign changes alter the learned co-expression
## programs: positive scores mean increased standardized program expression.
all_scores <- all_parameters <- list()
for (drug in drugs) for (dose in sort(unique(meta$dose[meta$compound == drug]))) {
  idx <- which(meta$compound == drug & meta$dose == dose)
  train_idx <- idx[meta$split[idx] == "Training"]
  center <- colMeans(effects[train_idx, gene_order, drop = FALSE])
  scale <- apply(effects[train_idx, gene_order, drop = FALSE], 2, sd)
  usable <- is.finite(scale) & scale > 0
  z <- sweep(sweep(effects[idx, gene_order, drop = FALSE], 2, center, "-"), 2, ifelse(usable, scale, 1), "/")
  all_parameters[[length(all_parameters) + 1L]] <- tibble(compound = drug, dose_uM = dose,
    gene = gene_order, program = unname(programs), training_center = center, training_scale = scale, used_for_score = usable)
  for (p in program_levels) {
    columns <- which(programs == p & usable)
    score <- if (length(columns)) rowMeans(z[, columns, drop = FALSE]) else rep(NA_real_, length(idx))
    all_scores[[length(all_scores) + 1L]] <- meta[idx, ] %>% transmute(compound, dose_uM = dose,
      observation_id, depmap_id, cell_line = tahoe_cell_name, split, prism_sensitization,
      program = p, n_scored_genes = length(columns), score = score)
  }
}
scores <- bind_rows(all_scores)
associations <- scores %>% group_by(compound, dose_uM, program, split) %>%
  group_modify(~test_correlation(.x$score, .x$prism_sensitization)) %>% ungroup()
fwrite(bind_rows(all_parameters), file.path(out, "program_gene_scaling.tsv"), sep = "\t", na = "NA")
fwrite(scores, file.path(out, "program_scores.tsv"), sep = "\t", na = "NA")
fwrite(associations, file.path(out, "program_prism_associations.tsv"), sep = "\t", na = "NA")
grDevices::cairo_pdf(file.path(out, "joint_program_prism_scatterplots.pdf"), width = 13.5, height = 13.6, onefile = TRUE)
for (p in program_levels) {
  doses <- sort(unique(meta$dose))
  d <- scores %>% filter(program == p, is.finite(score), is.finite(prism_sensitization)) %>%
    mutate(split = factor(split, levels = c("Training", "Holdout")),
      compound = factor(compound, levels = drugs),
      dose_label = factor(paste0(dose_uM, " uM"), levels = paste0(doses, " uM")))
  a <- associations %>% filter(program == p) %>% mutate(
    split = factor(split, levels = c("Training", "Holdout")),
    compound = factor(compound, levels = drugs),
    dose_label = factor(paste0(dose_uM, " uM"), levels = paste0(doses, " uM")),
    label = if_else(is.finite(pearson_r), sprintf("%s: n=%d, r=%.2f, p=%.2g", split, n, pearson_r, nominal_p),
      paste0(split, ": n=", n, ", r not estimable")),
    text_vjust = if_else(split == "Training", 1.3, 2.7))
  smooth <- semi_join(d, filter(a, is.finite(pearson_r)), by = c("compound", "dose_label", "split"))
  plot <- ggplot(d, aes(score, prism_sensitization, color = split)) +
    geom_point(size = 1.8, alpha = .8) +
    geom_smooth(data = smooth, method = "lm", formula = y ~ x, se = FALSE, linewidth = .7) +
    geom_text(data = a, aes(x = -Inf, y = Inf, label = label, color = split, vjust = text_vjust),
      inherit.aes = FALSE, hjust = -.02, size = 2.95, show.legend = FALSE) +
    facet_grid(compound ~ dose_label, scales = "free") +
    scale_color_manual(values = c(Training = "#2563A6", Holdout = "#D16D25"), name = NULL) +
    labs(x = "Frozen joint-program expression score", y = "PRISM sensitization",
      title = paste0(p, " | ", sum(programs == p), " genes | Same program across all drugs and doses"),
      subtitle = "Training selects genes and scaling; holdout applies the frozen program",
      caption = "Training associations reuse selection data. Trend lines and Pearson statistics are shown separately by split.\nHigher PRISM sensitization means greater drug sensitivity; positive program score means increased expression.") +
    theme_classic(base_size = 12) + theme(aspect.ratio = 1,
      strip.background = element_rect(fill = "#F1F4F6", color = NA),
      strip.text = element_text(size = 12, face = "bold"), legend.position = "bottom",
      plot.title = element_text(size = 16, face = "bold"), plot.caption = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10))
  print(plot)
}
dev.off()
writeLines(c("# Joint Tahoe-CytoSig programs", "",
  "1. Intersect the full Tahoe and CytoSig measured-gene universes. No old 3,000-gene or activity filter is applied.",
  "2. Within each drug/dose, test training Pearson association with PRISM at nominal two-sided p < 0.05. For each drug, rank significant genes by maximum absolute training correlation across doses, then minimum p, then gene name. Keep at most 100 per drug and combine the lists.",
  "3. Calculate one Tahoe correlation matrix over all training observations from the three drugs and all doses, and one CytoSig matrix over all 2,056 profiles. Average the two source matrices with equal weight.",
  "4. Cluster the joint matrix ONCE by hclust with method ward.D2, distance 1-r, and cut at EIGHT clusters as explicitly selected by the user. Silhouette over every valid rank is shown for diagnosis only; it does not choose the final rank.",
  "5. Freeze the joint gene order. The joint, Tahoe and CytoSig correlation matrices all use this exact order on both axes. No source-specific clustering is performed.",
  "6. Genes-by-samples heatmaps use this same row order. Sample columns are grouped by their strongest mean program z-score, then ordered within each group. Row scaling is display-only; missing values stay missing, and colors saturate at +/-2 SD.",
  "Compact and fully labelled PDFs show exactly the same genes and values. Large labelled PDF dimensions preserve readable labels; no gene subset is used for display.",
  "Unestimable gene-pair correlations are never imputed. The availability-only resolution removes the gene with most unestimable pairs, breaking ties by fewer observed profiles and then gene name; exclusions and original missing pairs are recorded. The per-drug cap is not topped up.",
  "The descriptive sharing tree groups the fixed joint programs by positive mean silhouette in each source. It does not infer ancestry or causal direction.",
  "The PRISM scatter PDF has one page per joint program, with all nine drug-dose combinations and separate training/holdout points, trends and statistics. Program score is mean gene-standardized expression change; gene scaling is learned within each drug/dose using training lines only. Genes with zero training SD cannot contribute to that dose's score.",
  "Training associations reuse selection data. The existing holdout has been inspected in prior project iterations, so these are held-out-from-this-fit results, not new external validation.",
  "All 2,056 CytoSig profiles are included; 54 lack metadata rows and are explicitly flagged. No new CytoSig data were downloaded.",
  "", "Run: Rscript /work/tahoe_prism_v2/scripts/pilot_tahoe_cytosig_programs.R"), file.path(out, "README.md"))
print(membership %>% dplyr::count(program))
print(coherence)
print(associations %>% filter(split == "Holdout"))
message("Completed ONE joint model: ", length(gene_order), " genes, ", fit$k, " programs; ", length(training), " Tahoe training profiles, ", nrow(cyto), " CytoSig profiles.")
