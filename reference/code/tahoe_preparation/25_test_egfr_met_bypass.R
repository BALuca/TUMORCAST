#!/usr/bin/env Rscript

## Joint EGFR--MET response programs: all doses, all CytoSig profiles.
## Restore the agreed small gene set, shared clustering and fixed block order.
if(file.exists("/work/renv/activate.R") && !requireNamespace("RhpcBLASctl",quietly=TRUE)) source("/work/renv/activate.R")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cluster)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"), mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
source(file.path(WORK_ROOT, "src/R/lib_heatmap.R"))
RhpcBLASctl::blas_set_num_threads(8L)
cache <- CANONICAL_CACHE
pdex <- TAHOE_PACKAGE
RESULT <- stage_result("25_egfr_met_bypass")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
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
  all(drugs %in% meta$compound))
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
stopifnot(all(adequate$include_in_discovery))
drugs <- adequate$compound[adequate$include_in_discovery]
meta <- meta %>% filter(compound %in% drugs)
fwrite(left_join(coverage, adequate, by = "compound"), file.path(RESULT, "drug_overlap.tsv"), sep = "\t")
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
abundance <- effects
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
    abundance[i, ] <- pmax(treated / meta$treated_total[i] * 10000 / meta$treated_correction_factor[i],
      reference / meta$reference_total[i] * 10000)
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

universe <- shared_universe
training <- which(meta$split == "Training")
gene_tests <- list()
for (drug in drugs) for (dose in sort(unique(meta$dose[meta$compound == drug]))) {
  ids <- which(meta$compound == drug & meta$dose == dose & meta$split == "Training" & is.finite(meta$prism_sensitization))
  x <- effects[ids, universe, drop = FALSE]; y <- meta$prism_sensitization[ids]; n <- length(y)
  if(n < 8L || sd(y)==0) next
  xc <- sweep(x, 2, colMeans(x)); yc <- y - mean(y)
  ssx <- colSums(xc^2); ssy <- sum(yc^2); sxy <- colSums(xc * yc)
  r <- sxy / sqrt(ssx * ssy)
  p <- 2 * pt(abs(r) * sqrt((n-2)/pmax(0,1-r^2)), n-2, lower.tail = FALSE)
  loo <- matrix(NA_real_, n, ncol(x))
  for (i in seq_len(n)) {
    vx <- ssx - n/(n-1)*xc[i, ]^2; vy <- ssy - n/(n-1)*yc[i]^2
    loo[i, ] <- (sxy - n/(n-1)*xc[i, ]*yc[i]) / sqrt(pmax(vx,0)*pmax(vy,0))
    loo[i, vx <= ssx*1e-10 | !is.finite(loo[i, ])] <- NA_real_
  }
  loo[is.finite(loo)] <- pmax(-1,pmin(1,loo[is.finite(loo)]))
  stable <- colSums(!is.finite(loo) | sweep(sign(loo), 2, sign(r), "!=") , na.rm = TRUE) == 0
  worst_r <- apply(abs(loo),2,min); worst_r[!is.finite(worst_r)] <- 0
  worst_p <- 2 * pt(worst_r * sqrt((n-3)/pmax(0,1-worst_r^2)),n-3,lower.tail=FALSE)
  detected <- colSums(abundance[ids, universe, drop = FALSE] >= .1)
  dominance <- apply(xc^2,2,max)/ssx
  gene_tests[[length(gene_tests)+1L]] <- tibble(compound=drug,dose_uM=dose,gene=universe,
    n_training_lines=n,pearson_r=r,nominal_p=p,minimum_loo_abs_r=worst_r,
    maximum_loo_p=worst_p,loo_sign_stable=stable,n_detected_lines=detected,
    median_abundance_cpm=apply(abundance[ids,universe,drop=FALSE],2,median)*100,
    max_single_line_variance_fraction=dominance,
    eligible=is.finite(p)&p<.05&stable&worst_p<.05&detected>=max(8L,ceiling(n/2)))
}
gene_tests <- bind_rows(gene_tests)
drug_selected <- gene_tests %>% filter(eligible) %>% arrange(compound,desc(abs(pearson_r)),nominal_p,gene,dose_uM) %>%
  distinct(compound,gene,.keep_all=TRUE) %>% group_by(compound) %>% slice_head(n=20L) %>%
  mutate(selection_rank=row_number()) %>% ungroup()

## At most 20 genes per drug keeps the union <=60, within the user's 30-60
## total target. Nominal p<.05 is retained. Detectability and single-line
## influence checks prevent the previously diagnosed sparse-gene artifacts.
## Ranking is by absolute training r among eligible genes, not held-RESULT fit.
candidate_genes <- drug_selected
stopifnot(all((candidate_genes %>% dplyr::count(compound))$n <= 20L))
fwrite(gene_tests, file.path(RESULT, "gene_prism_training_tests.tsv.gz"), sep="\t", na="NA")
cyto_dt <- cyto_dt[gene %chin% candidate_genes$gene]
cyto <- t(as.matrix(cyto_dt[, -1L, with = FALSE]))
colnames(cyto) <- cyto_dt$gene; rownames(cyto) <- profile_header
cyto[!is.finite(cyto)] <- NA_real_
rm(cyto_dt); invisible(gc())
fwrite(cyto_meta, file.path(RESULT, "cytosig_profiles.tsv"), sep = "\t", na = "NA")
fwrite(meta %>% group_by(compound, dose, split) %>% summarise(
  molecular_lines = n_distinct(depmap_id), prism_lines = sum(is.finite(prism_sensitization)),
  cell_line_ids = paste(depmap_id, collapse = ";"), .groups = "drop") %>%
  mutate(full_tahoe_genes = nrow(genes), cytosig_genes = length(cyto_genes),
    cytosig_profiles = nrow(cyto), shared_available_genes = sum(genes$feature %in% cyto_genes)),
  file.path(RESULT, "input_coverage.tsv"), sep = "\t")
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
fwrite(manifest, file.path(RESULT, "input_manifest.tsv"), sep = "\t")
message("Loaded all ", nrow(cyto), " CytoSig profiles; ", ncol(cyto), " of the capped candidate genes are measured.")

fit_clusters <- function(r) {
  stopifnot(nrow(r) >= 3L, all(is.finite(r)), max(abs(r - t(r))) < 1e-10)
  r[r > 1] <- 1; r[r < -1] <- -1; diag(r) <- 1
  distance <- as.dist(1 - r)
  tree <- hclust(distance, method = "average")
  ranks <- 2:min(8L, nrow(r) - 1L)
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
  cytosig_joint_profiles = cyto_pair_n[bad]), file.path(RESULT, "unestimable_gene_pairs.tsv"), sep = "\t")
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
stopifnot(length(selected_genes) >= 30L, length(selected_genes) <= 60L)
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
fwrite(candidate_genes, file.path(RESULT, "selected_genes_by_drug.tsv"), sep = "\t")
fwrite(selection_audit, file.path(RESULT, "gene_selection_audit.tsv"), sep = "\t", na = "NA")
membership$frozen_weight <- 1 / as.numeric(table(programs)[membership$program])
fwrite(membership, file.path(RESULT, "frozen_program_genes.tsv"), sep = "\t")
frozen_hash <- unname(tools::md5sum(file.path(RESULT, "frozen_program_genes.tsv")))
message("Frozen joint set: ",nrow(membership)," genes in ",length(program_levels)," programs.")
print(membership %>% dplyr::count(program))
fwrite(fit$curve, file.path(RESULT, "rank_silhouette.tsv"), sep = "\t")
p_rank <- ggplot(fit$curve, aes(rank, mean_silhouette)) + geom_line(linewidth = .8, color = "#255F85") +
  geom_point(size = 1.3, color = "#255F85") +
  geom_point(data = filter(fit$curve, selected), size = 3.7, shape = 21, fill = "#E9A23B", stroke = .9) +
  labs(x = "Number of joint programs", y = "Mean silhouette width",
    title = "One joint clustering of the combined drug-associated genes",
    subtitle = paste0(length(gene_order), " genes | rank fixed manually at ", fit$k),
    caption = "Distance = 1 - (Tahoe r + CytoSig r)/2. Eight programs chosen by the user; silhouette is diagnostic.") +
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 10), plot.margin = margin(9, 9, 9, 9))
ggsave(file.path(RESULT, "25d_rank_silhouette.pdf"), p_rank, width = 8.5, height = 5.5, device = cairo_pdf)

## Use observed off-diagonal magnitudes for the display scale, excluding the
## self-correlation diagonal. One symmetric scale is shared across sources.
off_diagonal <- unlist(lapply(matrices,function(m) m[upper.tri(m)]))
cor_limit <- max(.1,ceiling(10*as.numeric(quantile(abs(off_diagonal),.98)))/10)
cor_breaks <- seq(-cor_limit,cor_limit,length.out=5L)
cor_labels <- format(cor_breaks,trim=TRUE,digits=2)
cor_labels[c(1L,5L)] <- paste0(c("<=",">="),cor_labels[c(1L,5L)])
scale_audit <- bind_rows(lapply(names(matrices),function(name) {
  m <- matrices[[name]];v <- m[upper.tri(m)]
  tibble(source=name,off_diagonal_min=min(v),off_diagonal_max=max(v),display_min=-cor_limit,
    display_max=cor_limit,n_off_diagonal_pairs=length(v),n_saturated_pairs=sum(abs(v)>cor_limit),
    diagonal_display="Self-correlations masked; unchanged in exported matrices",
    scale_rule="98th percentile of pooled absolute off-diagonal r, rounded up to 0.1")
}))
fwrite(scale_audit,file.path(RESULT,"correlation_display_scale.tsv"),sep="\t")

## All correlation figures use EXACTLY the same genes, order and program split.
## A compact overview is followed by a fully labelled page when many genes are
## present; both pages use the same matrix, with no representative-gene subset.
for (source_name in names(matrices)) {
  r <- matrices[[source_name]][gene_order, gene_order, drop = FALSE]
  stopifnot(identical(rownames(r), gene_order), identical(colnames(r), gene_order))
  fwrite(data.table(gene = gene_order, r, check.names = FALSE),
    file.path(RESULT, paste0(tolower(source_name), "_gene_correlations.tsv.gz")), sep = "\t")
  for (labelled in TRUE) {
    suffix <- if (labelled) "_labelled" else ""
    size <- if (labelled) max(9, length(gene_order) * .16 + 3.4) else 9
    file <- file.path(RESULT, paste0(c(Joint="25a",Tahoe="25b",CytoSig="25c")[[source_name]], "_",tolower(source_name),"_correlation.pdf"))
    grDevices::cairo_pdf(file, width = size + .5, height = size, onefile = TRUE)
    display_r <- r; diag(display_r) <- NA_real_
    ht <- draw_heatmap(display_r, na_col="#ECECEC", col_anno_df = annotation, col_anno_cols = "Program",
      row_anno_df = annotation, row_anno_cols = "Program", palette_overrides = list(Program = palette),
      ann_top_args = list(show_legend = TRUE, annotation_legend_param = list(Program = list(
        title_gp = grid::gpar(fontsize = 12, fontface = "bold"), labels_gp = grid::gpar(fontsize = 12)))),
      ann_row_args = list(show_legend = FALSE),
      scale_rows = FALSE, value_range = c(-cor_limit, 0, cor_limit), value_colors = c("#2166AC", "white", "#B2182B"),
      hm_name = "Pearson r", show_row_names = labelled, show_column_names = labelled,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_order=seq_len(nrow(r)),column_order=seq_len(ncol(r)),
      row_labels = paste0(gene_order, " [", programs, "]"),
      row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 11),
      row_title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
      column_title_gp = grid::gpar(fontsize = 17, fontface = "bold"),
      row_gap = grid::unit(2, "mm"), column_gap = grid::unit(2, "mm"),
      title = paste0(source_name, " | Combined drug-associated genes\n",
        length(gene_order), " genes | ", fit$k, " joint programs | fixed gene order"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 11), at = cor_breaks, labels = cor_labels))
    drawn <- ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
      merge_legends = TRUE, padding = grid::unit(c(9, 7, 7, 7), "mm"))
    stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn)),seq_len(nrow(r))),
      identical(as.integer(ComplexHeatmap::column_order(drawn)),seq_len(ncol(r))))
    ComplexHeatmap::decorate_heatmap_body("Pearson r", {
      for (p in program_levels) {
        block <- which(programs == p)
        position <- (min(block) + max(block) - 1) / (2 * length(gene_order))
        grid::grid.rect(x = position, y = 1 - position, width = length(block) / length(gene_order),
          height = length(block) / length(gene_order), gp = grid::gpar(fill = NA, col = "#71717A", lwd = .35))
      }
    })
    grid::grid.text(paste0("Common scale ±",cor_limit,"; tail values saturated. Gray diagonal: self-correlations omitted. Both sources use the identical gene order."),
      x = .5, y = grid::unit(4, "mm"), gp = grid::gpar(fontsize = 10))
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
    file.path(RESULT, paste0(tolower(source_name), "_gene_sample_effects.tsv.gz")), sep = "\t", na = "NA")
  fwrite(tibble(source = source_name, gene = gene_order, display_center = center, display_scale = scale),
    file.path(RESULT, paste0(tolower(source_name), "_heatmap_scaling.tsv")), sep = "\t")
  for (labelled in TRUE) {
    height <- if (labelled) max(8, length(gene_order) * .16 + 2.8) else 9
    suffix <- if (labelled) "_labelled" else ""
    grDevices::cairo_pdf(file.path(RESULT, paste0(c(Tahoe="25e",CytoSig="25f")[[source_name]], "_",tolower(source_name),"_gene_samples.pdf")),
      width = 14, height = height, onefile = TRUE)
    ht <- draw_heatmap(z, col_anno_df = col_annotation, col_anno_cols = names(col_annotation),
      row_anno_df = annotation, row_anno_cols = "Program",
      ann_top_args = list(show_legend = TRUE, annotation_legend_param = list(Program = list(
        title_gp = grid::gpar(fontsize = 12, fontface = "bold"), labels_gp = grid::gpar(fontsize = 12)))),
      ann_row_args = list(show_legend = FALSE),
      palette_overrides = list(Program = c(palette, Unmeasured = "#AAAAAA")),
      scale_rows = FALSE, value_range = c(-2, 0, 2), value_colors = c("#2166AC", "white", "#B2182B"),
      hm_name = "Gene z-score", show_row_names = labelled, show_column_names = FALSE,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_order=seq_len(nrow(z)),column_order=seq_len(ncol(z)), row_labels = paste0(gene_order, " [", programs, "]"),
      row_names_gp = grid::gpar(fontsize = 11), row_title_gp = grid::gpar(fontsize = 15, fontface = "bold"),
      column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      row_gap = grid::unit(2, "mm"), column_gap = grid::unit(2, "mm"), na_col = "#B7B7B7",
      title = paste0(source_name, " | Genes x samples\n", length(gene_order), " genes | ", ncol(z),
        " samples | fixed genes; samples grouped by program"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 11), at = c(-2, -1, 0, 1, 2)))
    drawn <- ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
      merge_legends = TRUE, padding = grid::unit(c(10, 7, 7, 7), "mm"))
    stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn)),seq_len(nrow(z))),
      identical(as.integer(ComplexHeatmap::column_order(drawn)),seq_len(ncol(z))))
    ComplexHeatmap::decorate_heatmap_body("Gene z-score", {
      for (boundary in cumsum(as.integer(table(split)))[-length(program_levels)])
        grid::grid.lines(x = c(0, 1), y = rep(1 - boundary / length(gene_order), 2), gp = grid::gpar(col = "#FFFFFF", lwd = .4))
      column_counts <- as.integer(table(col_split)); column_counts <- column_counts[column_counts > 0]
      for (boundary in head(cumsum(column_counts), -1))
        grid::grid.lines(x = rep(boundary / ncol(z), 2), y = c(0, 1), gp = grid::gpar(col = "#FFFFFF", lwd = .4))
    })
    grid::grid.text("Rows: fixed joint-program gene order. Columns: individual response profiles. Gray: missing. Color saturates at +/-2 SD.",
      x = .5, y = grid::unit(4, "mm"), gp = grid::gpar(fontsize = 10))
    dev.off()
  }
}
fwrite(bind_rows(sample_orders), file.path(RESULT, "sample_plotting_order.tsv"), sep = "\t", na = "NA")


## Quantify each fixed cluster in both sources; do not recluster a source to
## improve its apparent agreement. Two-gene groups have only one correlation.
coherence <- bind_rows(lapply(program_levels, function(p) {
  g <- gene_order[programs == p]
  within <- vapply(matrices, function(m) if (length(g) >= 2L) mean(m[g,g][upper.tri(m[g,g])]) else NA_real_, numeric(1))
  tibble(program = p, n_genes = length(g), joint_mean_r = within[["Joint"]],
    tahoe_mean_r = within[["Tahoe"]], cytosig_mean_r = within[["CytoSig"]], genes = paste(g, collapse = ", "))
}))
for (s in c("Tahoe", "CytoSig")) {
  dd <- 1 - matrices[[s]][gene_order, gene_order]; dd[dd < 0] <- 0
  ss <- cluster::silhouette(as.integer(split), as.dist(dd))
  coherence[[paste0(tolower(s), "mean_silhouette")]] <- vapply(program_levels,
    function(p) mean(ss[programs == p, "sil_width"]), numeric(1))
}
coherence <- coherence %>% mutate(source_support = case_when(
  n_genes < 3L ~ "Two-gene cluster; limited evidence",
  tahoe_mean_r >= .4 & cytosig_mean_r >= .4 ~ "Coherent in both sources",
  tahoe_mean_r >= .4 ~ "Tahoe-supported", cytosig_mean_r >= .4 ~ "CytoSig-supported", TRUE ~ "Weak source coherence"))
fwrite(coherence, file.path(RESULT, "program_source_coherence.tsv"), sep = "\t", na = "NA")

## All three sources are shown together with order assertions after drawing.
comparison <- lapply(seq_along(matrices), function(i) {
  m <- matrices[[i]][gene_order, gene_order]; diag(m) <- NA_real_
  draw_heatmap(m, na_col="#ECECEC", scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = seq_along(gene_order), column_order = seq_along(gene_order),
    row_anno_df = if (i == 1L) annotation else NULL, row_anno_cols = if (i == 1L) "Program" else NULL,
    palette_overrides = list(Program = palette), ann_row_args = list(show_legend = TRUE),
    show_row_names = i == 1L, row_names_side = "left", show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 9),
    hm_name = paste0("r", i), title = names(matrices)[i], column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
    value_range = c(-cor_limit,0,cor_limit), value_colors = c("#2166AC", "white", "#B2182B"),
    show_heatmap_legend = i == 3L, heatmap_legend_param = list(title = "Pearson r", at = cor_breaks, labels = cor_labels))
})
comparison <- comparison[[1]] + comparison[[2]] + comparison[[3]]
draw_comparison <- function() {
  drawn <- ComplexHeatmap::draw(comparison, merge_legends = TRUE, padding = grid::unit(c(11,6,6,6), "mm"))
  for (i in 1:3) {
    name <- paste0("r", i)
    stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn, name = name)), seq_along(gene_order)),
      identical(as.integer(ComplexHeatmap::column_order(drawn, name = name)), seq_along(gene_order)))
    ComplexHeatmap::decorate_heatmap_body(name, {
      for (b in head(cumsum(as.integer(table(split))), -1)) {
        grid::grid.lines(c(0,1), rep(1 - b/length(gene_order), 2), gp = grid::gpar(col = "white", lwd = 1.7))
        grid::grid.lines(rep(b/length(gene_order), 2), c(0,1), gp = grid::gpar(col = "white", lwd = 1.7))
      }
    })
  }
  grid::grid.text(paste0("Shared scale: -",cor_limit," to +",cor_limit,". Extreme values saturate; gray self-correlation diagonal omitted. Raw correlations are retained in the tables."),
    x=.5,y=grid::unit(4,"mm"),gp=grid::gpar(fontsize=12))
}
cairo_pdf(file.path(RESULT, "25g_correlation_comparison.pdf"), width = 22, height = 9.4)
draw_comparison(); dev.off()
png(file.path(RESULT, "25g_correlation_comparison.png"), width = 3300, height = 1410, res = 150, type = "cairo")
draw_comparison(); dev.off()

## Render the actual joint gene tree used to define the eight programs.
hc <- fit$tree; yy <- setNames(rev(seq_along(gene_order)), gene_order); nodes <- segments <- list()
for (i in seq_len(nrow(hc$merge))) {
  child <- lapply(hc$merge[i,], function(j) if (j < 0) c(0, yy[hc$labels[-j]]) else nodes[[j]])
  xs <- vapply(child, `[`, numeric(1), 1L); ys <- vapply(child, `[`, numeric(1), 2L)
  nodes[[i]] <- c(hc$height[i], mean(ys))
  segments[[i]] <- tibble(x = c(xs,hc$height[i]), xend = rep(hc$height[i],3), y = c(ys,min(ys)), yend = c(ys,max(ys)))
}
leaves <- membership %>% mutate(y = unname(yy[gene]))
gene_tree <- ggplot() + geom_segment(data = bind_rows(segments), aes(x,y,xend=xend,yend=yend), color = "#506070", linewidth = .5) +
  geom_point(data = leaves, aes(0,y,color=program), size = 2) +
  geom_text(data = leaves, aes(-.03,y,label=paste0(gene," [",program,"]"),color=program), hjust = 0, size = 3.6, show.legend = FALSE) +
  scale_color_manual(values=palette) + scale_x_reverse(limits=c(max(hc$height)+.03,-.58),breaks=c(0,.25,.5,.75,1,1.25)) +
  scale_y_continuous(NULL,breaks=NULL,expand=expansion(add=.8)) +
  labs(x="Average-linkage distance: 1 - mean(Tahoe r, CytoSig r)",color=NULL,
    title="The joint gene tree that defines the programs",subtitle=paste(length(gene_order),"genes; eight clusters fixed by the user")) +
  theme_classic(base_size=12) + theme(axis.line.y=element_blank(),legend.position="bottom",plot.title=element_text(face="bold"))
ggsave(file.path(RESULT,"25h_gene_hierarchy.pdf"),gene_tree,width=9.8,height=11.3,device=cairo_pdf)

## One training center/scale and positive weights, shared across drugs/doses.
## All prediction/scoring parameters are frozen before associations are tested.
scaling <- tibble(gene=gene_order,program=unname(programs),training_center=colMeans(effects[training,gene_order]),
  training_scale=apply(effects[training,gene_order],2L,sd),frozen_weight=membership$frozen_weight)
stopifnot(all(is.finite(scaling$training_scale)),all(scaling$training_scale>0))
fwrite(scaling,file.path(RESULT,"frozen_gene_scaling.tsv"),sep="\t")
z <- sweep(sweep(effects[,gene_order],2L,scaling$training_center,"-"),2L,scaling$training_scale,"/")
scores <- bind_rows(lapply(program_levels,function(p) meta %>% mutate(program=p,score=rowMeans(z[,programs==p,drop=FALSE]))))
associations <- scores %>% group_by(compound,dose,program,split) %>%
  group_modify(~test_correlation(.x$score,.x$prism_sensitization)) %>% ungroup() %>%
  group_by(split) %>% mutate(bh_q=p.adjust(nominal_p,"BH")) %>% ungroup()
fwrite(scores %>% select(observation_id,compound,dose,depmap_id,tahoe_cell_name,split,program,score,prism_sensitization),
  file.path(RESULT,"program_scores.tsv"),sep="\t",na="NA")
fwrite(associations,file.path(RESULT,"program_prism_associations.tsv"),sep="\t",na="NA")
cairo_pdf(file.path(RESULT,"25i_program_prism_scatterplots.pdf"),width=11.5,height=11.5,onefile=TRUE)
for(p in program_levels) {
  d <- scores %>% filter(program==p,is.finite(prism_sensitization))
  a <- associations %>% filter(program==p) %>% mutate(label=sprintf("%s n=%d r=%.2f p=%.2g",split,n,pearson_r,nominal_p),
    v=if_else(split=="Training",1.3,2.7))
  smooth <- semi_join(d,a %>% filter(is.finite(pearson_r)),by=c("compound","dose","split"))
  print(ggplot(d,aes(score,prism_sensitization,color=split))+geom_point(size=1.7)+
    geom_smooth(data=smooth,method="lm",formula=y~x,se=FALSE,linewidth=.6)+
    geom_text(data=a,aes(x=-Inf,y=Inf,label=label,vjust=v),inherit.aes=FALSE,hjust=-.03,size=2.7)+
    facet_grid(compound~dose,scales="free")+scale_color_manual(values=c(Training="#287A91",Holdout="#B86C38"))+
    labs(x="Frozen program expression score",y="PRISM sensitization",title=paste0(p," | ",sum(programs==p)," genes"),
      subtitle="All doses in uM; the same training scaling and program weights in every panel",
      caption="Training associations reuse selection data. Holdout was inspected in earlier iterations.\nPositive program score means increased expression; higher sensitization means greater drug sensitivity.")+
    theme_classic(base_size=11)+theme(aspect.ratio=1,legend.position="bottom",plot.caption=element_text(hjust=0,size=9),plot.margin=margin(12,8,8,8)))
}
dev.off()

## Annotate programs by drug provenance. Counts are not causal directions.
enrichment <- bind_rows(lapply(program_levels,function(p) bind_rows(lapply(drugs,function(drug) {
  in_p <- gene_order %in% names(programs)[programs==p]; in_drug <- gene_order %in% candidate_genes$gene[candidate_genes$compound==drug]
  tab <- matrix(c(sum(in_p&in_drug),sum(in_p&!in_drug),sum(!in_p&in_drug),sum(!in_p&!in_drug)),2,byrow=TRUE)
  a <- fisher.test(tab,alternative="greater")
  tibble(program=p,compound=drug,n_selected_genes=tab[1,1],n_program_genes=sum(in_p),n_background=length(gene_order),
    nominal_p=a$p.value,odds_ratio=unname(a$estimate))
})))) %>% mutate(bh_q=p.adjust(nominal_p,"BH"))
fwrite(enrichment,file.path(RESULT,"program_drug_enrichment.tsv"),sep="\t",na="NA")

## Descriptive hierarchy of the frozen programs, using their between-program
## joint correlations. This tree contains no causal upstream/downstream arrows.
between <- matrix(0,length(program_levels),length(program_levels),dimnames=list(program_levels,program_levels))
for(p in program_levels) for(q in program_levels) between[p,q] <- if(p==q) 1 else
  mean(matrices$Joint[names(programs)[programs==p],names(programs)[programs==q]])
ptree <- hclust(as.dist(1-between),method="average")
py <- setNames(rev(seq_along(program_levels)),ptree$labels[ptree$order]);pnodes <- psegs <- list()
for(i in seq_len(nrow(ptree$merge))) {
  child <- lapply(ptree$merge[i,],function(j) if(j<0) c(0,py[ptree$labels[-j]]) else pnodes[[j]])
  xs <- vapply(child,`[`,numeric(1),1L);ys <- vapply(child,`[`,numeric(1),2L)
  pnodes[[i]] <- c(ptree$height[i],mean(ys))
  psegs[[i]] <- tibble(x=c(xs,ptree$height[i]),xend=rep(ptree$height[i],3),y=c(ys,min(ys)),yend=c(ys,max(ys)))
}
leaves <- coherence %>% mutate(y=unname(py[program]))
pplot <- ggplot()+geom_segment(data=bind_rows(psegs),aes(x,y,xend=xend,yend=yend),linewidth=.7,color="#506070")+
  geom_point(data=leaves,aes(0,y,color=program),size=3)+
  geom_text(data=leaves,aes(-.04,y,label=paste0(program," | ",n_genes," genes | ",source_support,"\n",
    sprintf("Within r: Tahoe %.2f; CytoSig %.2f",tahoe_mean_r,cytosig_mean_r))),hjust=0,size=3.7)+
  scale_color_manual(values=palette,guide="none")+scale_x_reverse(limits=c(max(ptree$height)+.08,-1.8),breaks=c(0,.25,.5,.75,1,1.25))+
  scale_y_continuous(NULL,breaks=NULL,expand=expansion(add=.7))+
  labs(x="Between-program joint correlation distance",title="Descriptive program hierarchy",
    subtitle="The eight frozen clusters; hierarchy uses the mean correlation between their genes",
    caption="Branch lengths describe co-response, not causal ordering. Source-support labels use mean within-program r >= 0.4.\nSmall clusters and weak source coherence remain explicit.")+theme_classic(base_size=12)+
  theme(axis.line.y=element_blank(),plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(RESULT,"25j_program_hierarchy.pdf"),pplot,width=12,height=7.5,device=cairo_pdf)

## HGF remains the bypass perturbation: use its individual profiles to form
## the median HGF effect. The correlation fit above still uses ALL CytoSig.
hgf_ids <- cyto_meta$profile_id[tolower(trimws(cyto_meta$treatment)) %in% c("hgf","hepatocyte growth factor")]
stopifnot(length(hgf_ids)>0L)
hgf_median <- apply(cyto[hgf_ids,gene_order,drop=FALSE],2L,median,na.rm=TRUE)
hgf_projection <- vapply(program_levels,function(p) mean(hgf_median[programs==p]/scaling$training_scale[programs==p]),numeric(1))
source_effects <- bind_rows(lapply(program_levels,function(p) meta %>% filter(split=="Training") %>%
  mutate(program=p,displacement=rowMeans(sweep(effects[training,gene_order[programs==p],drop=FALSE],
    2L,scaling$training_scale[programs==p],"/"))) %>% group_by(program,compound,dose) %>%
  summarise(mean_displacement=mean(displacement),.groups="drop")))
checks <- source_effects %>% mutate(hgf_displacement=hgf_projection[program],hgf_opposes=mean_displacement*hgf_displacement<0,
  n_hgf_profiles=length(hgf_ids),interpretation="Molecular direction relative to zero effect; no causal program ordering inferred")
fwrite(checks,file.path(RESULT,"program_hgf_directional_checks.tsv"),sep="\t",na="NA")
fwrite(tibble(gene=gene_order,program=unname(programs),hgf_median_effect=hgf_median),file.path(RESULT,"program_hgf_gene_effects.tsv"),sep="\t",na="NA")
wilson <- tibble(condition=c("Control","Erlotinib","Erlotinib + HGF","Erlotinib + HGF + crizotinib"),
  pERK=c(1,.180,.979,.085),pMET=c(1,1.098,11.660,1.095),
  source="Wilson Figure 2C, values reproduced in https://elifesciences.org/articles/04037",
  interpretation="Known molecular reference; Wilson does not supply transcriptomes for these conditions")
fwrite(wilson,file.path(RESULT,"wilson_khm3s_reference.tsv"),sep="\t")
stopifnot(identical(unname(tools::md5sum(file.path(RESULT,"frozen_program_genes.tsv"))),frozen_hash))
writeLines(c("# Joint EGFR–MET response programs", "",
  "Stage 25 now carries forward the agreed analysis, replacing the uncapped logical-branch implementation archived in superseded_rule_branches/.",
  paste0("Final model: ",length(gene_order)," genes; eight joint clusters. All three drugs, all available doses, all 2,056 CytoSig profiles."),"",
  "## How the programs are produced", "",
  "1. Intersect all measured Tahoe/CytoSig genes (19,339) before selection. Do not load the old 3,000-gene list or molecular-activity gate.",
  "2. Test gene-to-PRISM Pearson correlations separately by drug/dose in the original training lines. Require nominal p<0.05 and at least eight matched lines. Retain the previous repair's gene detectability requirement (>=10 CPM in at least half the lines and at least eight) and same-sign, finite nominal p<0.05 after every single-line deletion. These checks address previously observed sparse and single-line-driven genes.",
  "3. Rank eligible genes by absolute training r across doses (ties: nominal p, gene, dose). Take at most 20 genes per drug to meet the requested 30–60-total scale; combine and deduplicate. Genes failing correlation availability are audited, never imputed or topped up.",
  "4. Compute gene–gene Pearson matrices over all Tahoe training observations and all CytoSig profiles. Joint R = (R_Tahoe + R_CytoSig)/2. Fit average-linkage hclust to distance 1-R, once, then cut into eight clusters as previously requested. Rank 2–8 silhouette is diagnostic only.",
  "5. Freeze the resulting gene order. Every correlation matrix uses it on BOTH axes. White/outlined block boundaries identify programs; no source is independently reordered. Each gene-by-sample heatmap uses the identical rows; columns are grouped by their largest program z-score and ordered within that group. Thus the ordering exposes blocks without changing the underlying correlations or forcing off-block values to zero.","",
  "## Figures", "",
  "- 25a/25b/25c: joint, Tahoe and CytoSig gene–gene heatmaps, all genes labelled.",
  "- 25d: rank versus silhouette, eight marked as user-selected.",
  "- 25e/25f: genes versus individual Tahoe/CytoSig response samples, with the same gene order.",
  "- 25g: all three correlation matrices side by side (PDF and PNG).",
  "- 25h: the actual joint gene hierarchy.",
  "- 25i: program–PRISM scatterplots, one page per program, all drug/dose panels with split-specific trends.",
  "- 25j: descriptive hierarchy of the frozen programs, with source-coherence statistics.","",
  paste0("Correlation displays share ±",cor_limit," limits derived from the pooled off-diagonal values. Self-correlations are masked. Endpoint colors saturate; exact values are exported unchanged. See correlation_display_scale.tsv."),
  "The source-coherence table reports mean within-program r and silhouette. Eight clusters are a user-fixed partition, not a claim of eight independently supported biological programs; small or weak clusters remain visibly labelled. A diagonal block layout cannot make weak correlations strong. Drug enrichment describes selection provenance, not causal ordering.",
  "Scoring freezes one gene center/scale from all Tahoe training observations and positive equal within-program weights, used for every drug/dose. Positive scores mean increased expression; PRISM retains sensitization orientation. The current holdout was inspected previously, so these are held out from this fit, not a new independent validation set.",
  "HGF is used for the bypass-direction check via its median individual-profile effect; all CytoSig profiles are used for correlation discovery. The published Wilson pERK/pMET table remains a known mechanistic reference, not a transcriptomic projection or a fitted viability test. No Wilson observation selects genes or programs.","",
  "Run from /work: `Rscript /work/tahoe_prism_v2/scripts/25_test_egfr_met_bypass.R`."),file.path(RESULT,"README.md"))
fwrite(tibble(file=sort(list.files(RESULT,pattern="^25.*\\.pdf$"))),file.path(RESULT,"FIGURE_INDEX.tsv"),sep="\t")
message("Completed corrected Stage 25: ",length(gene_order)," genes, eight joint clusters, all source heatmaps restored.")
print(coherence)
