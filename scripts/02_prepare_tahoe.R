#!/usr/bin/env Rscript
# Rebuild sample-local normalization from released Tahoe mean-response arrays.
# All 62,710 genes determine composition scaling; only then retain discovery genes.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(qs)
  library(digest)
})
source("R/data_io.R")

raw_dir <- "data/raw/tahoe/pdex_selected"
cache_path <- "data/processed/tahoe_discovery.qs"
output_dir <- "results/02_prepare_tahoe"
workers <- 4L
minimum_cp10k <- 0.1
fold_change_trim <- c(0.3, 0.7)
abundance_trim <- c(0.05, 0.95)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
expression <- qread("data/processed/expression.qs")$result
genes <- expression$genes
discovery_cells <- expression$eligibility_cells
rm(expression)

raw_paths <- list.files(raw_dir, pattern = "\\.(npy|tsv|tsv.gz)$", recursive = TRUE, full.names = TRUE)
fingerprint <- list(
  inputs = vapply(c(raw_paths, "data/processed/expression.qs"), digest,
    character(1), file = TRUE, algo = "sha256"),
  code = vapply(c("scripts/02_prepare_tahoe.R", "R/data_io.R"), digest,
    character(1), file = TRUE, algo = "sha256"),
  minimum_cp10k = minimum_cp10k,
  fold_change_trim = fold_change_trim,
  abundance_trim = abundance_trim,
  R = R.version.string
)
cached <- if (file.exists(cache_path)) qread(cache_path) else NULL
if (!is.null(cached) && identical(cached$fingerprint, fingerprint)) {
  cat("Reusing matching Tahoe normalization cache.\n")
  result <- cached$result
} else {
  cat("Rebuilding Tahoe normalization: absent or stale cache.\n")
  gene_table <- fread(file.path(raw_dir, "genes.tsv")) %>% arrange(gene_index)
  gene_indices <- match(genes, gene_table$feature)
  stopifnot(!anyNA(gene_indices), !anyDuplicated(gene_table$feature))
  perturbations <- fread(file.path(raw_dir, "perturbations.tsv"))
  availability <- fread(file.path(raw_dir, "response_availability.tsv"))
  metadata <- availability %>%
    filter(available == 1, depmap_id %in% discovery_cells, depmap_id != "ACH-000248") %>%
    left_join(select(perturbations, condition_index, drug, dose, dose_unit, dose_panel),
      by = "condition_index") %>%
    arrange(depmap_id, condition_index) %>%
    mutate(observation_id = paste(depmap_id, perturbation_id, sep = "__"))
  stopifnot(nrow(metadata) == 36050L, !anyDuplicated(metadata$observation_id))
  plates <- fread(file.path(raw_dir, "pdex_plate_audit.tsv.gz")) %>%
    mutate(observation_id = paste(depmap_id, perturbation_id, sep = "__"))
  plate_lists <- split(plates$plate, plates$observation_id)
  reference_metadata <- fread(file.path(raw_dir, "pdex_plate_reference_metadata.tsv"))
  references <- read_npy_matrix(file.path(raw_dir, "pdex_plate_reference.npy"))
  reference_keys <- paste(reference_metadata$depmap_id, reference_metadata$plate, sep = "__")
  stopifnot(!anyDuplicated(reference_keys), ncol(references) == nrow(gene_table))
  rownames(references) <- reference_keys[order(reference_metadata$array_row)]

  # Each independent cell is a sizeable block. Forked workers share read-only
  # reference arrays. There is one calculation and no separate workflow layer.
  normalize_cell <- function(cell_id) {
    cell_metadata <- metadata[metadata$depmap_id == cell_id, ]
    shard_path <- file.path(raw_dir, "response_shards",
      paste0(cell_id, "__", cell_metadata$cell_line[1], ".npy"))
    shard <- read_npy_matrix(shard_path)
    effects <- matrix(NA_real_, nrow(cell_metadata), length(genes),
      dimnames = list(cell_metadata$observation_id, genes))
    records <- vector("list", nrow(cell_metadata))
    for (i in seq_len(nrow(cell_metadata))) {
      observation <- cell_metadata$observation_id[i]
      contributing <- plate_lists[[observation]]
      stopifnot(length(contributing) == cell_metadata$plate_count[i], !anyDuplicated(contributing))
      keys <- paste(cell_id, contributing, sep = "__")
      stopifnot(all(keys %in% rownames(references)))
      reference <- colMeans(references[keys, , drop = FALSE])
      treated <- shard[cell_metadata$condition_index[i] + 1L, ] + reference
      tolerance <- max(1e-6, 32 * 2^-23 * max(abs(c(treated, reference)), 1))
      stopifnot(min(c(treated, reference)) >= -tolerance, all(is.finite(treated)))
      minimum_treated <- min(treated)
      treated <- pmax(treated, 0)
      reference <- pmax(reference, 0)
      treated_total <- sum(treated)
      reference_total <- sum(reference)
      stopifnot(treated_total > 0, reference_total > 0)
      treated_proportion <- treated / treated_total
      reference_proportion <- reference / reference_total
      eligible <- treated_proportion * 10000 >= minimum_cp10k &
        reference_proportion * 10000 >= minimum_cp10k
      # Keep the original order of logarithms. Algebraically rearranging this
      # changes ties at trimming boundaries in sparse abundance profiles.
      log_treated <- log2(treated_proportion[eligible])
      log_reference <- log2(reference_proportion[eligible])
      log_ratio <- log_treated - log_reference
      log_abundance <- 0.5 * (log_treated + log_reference)
      ratio_bounds <- quantile(log_ratio, fold_change_trim)
      abundance_bounds <- quantile(log_abundance, abundance_trim)
      retained <- log_ratio >= ratio_bounds[1] & log_ratio <= ratio_bounds[2] &
        log_abundance >= abundance_bounds[1] & log_abundance <= abundance_bounds[2]
      stopifnot(sum(retained) > 0)
      offset <- median(log_ratio[retained])
      correction <- 2^offset
      effects[i, ] <- log2(1 + treated_proportion[gene_indices] * 10000 / correction) -
        log2(1 + reference_proportion[gene_indices] * 10000)
      records[[i]] <- data.frame(
        observation_id = observation, treated_total = treated_total,
        reference_total = reference_total, composition_offset_log2 = offset,
        treated_correction_factor = correction, n_scaling_eligible = sum(eligible),
        n_scaling_retained = sum(retained), minimum_reconstructed_treated = minimum_treated,
        negative_tolerance = tolerance, contributing_plates = paste(contributing, collapse = ",")
      )
    }
    message(cell_id, ": normalized ", nrow(cell_metadata), " profiles")
    list(effects = round_float32(effects), normalization = bind_rows(records))
  }
  cells <- unique(metadata$depmap_id)
  blocks <- parallel::mclapply(cells, normalize_cell, mc.cores = workers)
  if (any(vapply(blocks, inherits, logical(1), "try-error"))) stop("Tahoe normalization worker failed")
  effects <- do.call(rbind, lapply(blocks, `[[`, "effects"))
  normalization <- bind_rows(lapply(blocks, `[[`, "normalization"))
  effects <- effects[metadata$observation_id, , drop = FALSE]
  stopifnot(all(is.finite(effects)), identical(rownames(effects), metadata$observation_id))
  result <- list(metadata = metadata, normalization = normalization, effects = effects, genes = genes)
  save_cache(result, fingerprint, cache_path)
}
write.table(result$normalization, gzfile(file.path(output_dir, "normalization.tsv.gz")),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(result$metadata, gzfile(file.path(output_dir, "observed_profiles.tsv.gz")),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Original normalization records are comparison targets, never inputs to scaling.
reference_path <- "reference/results/tahoe_normalization/normalization_records.tsv.gz"
if (file.exists(reference_path)) {
  original <- fread(reference_path)
  index <- match(result$normalization$observation_id, original$observation_id)
  comparison <- data.frame(
    field = c("treated_total", "reference_total", "treated_correction_factor"),
    maximum_absolute_difference = vapply(
      c("treated_total", "reference_total", "treated_correction_factor"), function(field) {
        max(abs(result$normalization[[field]] - original[[field]][index]))
      }, numeric(1))
  )
  write.table(comparison, file.path(output_dir, "reference_comparison.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  print(comparison)
  # The historical array implementation and R differ at floating-point ties in
  # sparse trimming boundaries. Preserve the computed R result and report it;
  # do not substitute the saved historical scaling factors. The downstream
  # correlation and bank comparisons separately determine preservation there.
  stopifnot(comparison$maximum_absolute_difference[3] < 1e-4)
}
cat(nrow(result$effects), "profiles normalized from released means.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
