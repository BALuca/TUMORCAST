#!/usr/bin/env Rscript
# Rebuild source-balanced, signed gene relationships for molecular discovery.
# Profile centering removes a response-wide shift before gene co-response ranking.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(qs)
  library(digest)
})
source("R/data_io.R")

output_dir <- "results/03_build_gene_relationships"
cache_path <- "data/processed/gene_relationships.qs"
cytosig_path <- "data/raw/cytosig/diff.merge.gz"
metadata_path <- "data/raw/cytosig/meta_info.txt"
RhpcBLASctl::blas_set_num_threads(8)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
input_paths <- c("data/processed/tahoe_discovery.qs", cytosig_path, metadata_path)
fingerprint <- list(
  inputs = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256"),
  code = digest("scripts/03_build_gene_relationships.R", file = TRUE, algo = "sha256"),
  R = R.version.string, matrixStats = as.character(packageVersion("matrixStats"))
)
# Downstream staleness depends on the scientific data, not a changed timestamp,
# a different compression stream or updated upstream comparison diagnostics.
fingerprint$inputs[1] <- digest(qread(input_paths[1])$result, algo = "sha256")
cached <- if (file.exists(cache_path)) qread(cache_path) else NULL
if (!is.null(cached) && identical(cached$fingerprint, fingerprint)) {
  cat("Reusing matching gene-relationship cache.\n")
  result <- cached$result
} else {
  cat("Rebuilding gene relationships: absent or stale cache.\n")
  tahoe <- qread(input_paths[1])$result
  genes <- tahoe$genes
  relative <- tahoe$effects - rowMeans(tahoe$effects)
  drug_counts <- table(tahoe$metadata$drug)
  weights <- 1 / as.numeric(drug_counts[tahoe$metadata$drug])
  weights <- weights / sum(weights)
  profile_metadata <- tahoe$metadata
  rm(tahoe)
  ranks <- matrixStats::colRanks(relative, ties.method = "average", preserveShape = TRUE)
  rm(relative)
  centered <- sweep(ranks, 2, colSums(ranks * weights), "-")
  rm(ranks)
  covariance <- crossprod(centered * sqrt(weights))
  scale <- sqrt(pmax(diag(covariance), 1e-30))
  tahoe_correlation <- pmax(pmin(covariance / outer(scale, scale), 1), -1)
  diag(tahoe_correlation) <- 1
  dimnames(tahoe_correlation) <- list(genes, genes)
  rm(centered, covariance)
  gc()
  cat("Tahoe correlation complete. Reading CytoSig.\n")

  # The release omits a heading for the first (gene) column.
  connection <- gzfile(cytosig_path, "rt")
  profile_names <- strsplit(readLines(connection, n = 1), "\t", fixed = TRUE)[[1]]
  close(connection)
  cytosig <- fread(cytosig_path, skip = 1, header = FALSE,
    col.names = c("gene", profile_names), showProgress = FALSE)
  index <- match(genes, cytosig$gene)
  stopifnot(!anyNA(index))
  responses <- as.matrix(cytosig[index, -1, with = FALSE])
  rownames(responses) <- genes
  rm(cytosig)
  annotations <- fread(metadata_path, skip = 1, header = FALSE,
    col.names = c("profile_id", "treatment", "condition", "duration", "dose"))
  stopifnot(!anyDuplicated(annotations$profile_id))
  annotations <- data.frame(profile_id = colnames(responses)) %>%
    left_join(annotations, by = "profile_id")
  excluded <- apply(annotations, 1, function(values) {
    any(grepl("AU[ _-]?565|ACH-000248", values, ignore.case = TRUE))
  })
  write.table(annotations[excluded, ], file.path(output_dir, "excluded_profiles.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  annotations <- annotations[!excluded, ]
  responses <- responses[, !excluded, drop = FALSE]
  missing_treatment <- is.na(annotations$treatment) | !nzchar(annotations$treatment)
  annotations$treatment[missing_treatment] <-
    sub("[&@].*", "", annotations$profile_id[missing_treatment])
  annotations$study <- sub(".*@((?:GSE[0-9]+|E-(?:MTAB|MEXP|TABM|SGRP)-[0-9]+)).*",
    "\\1", annotations$profile_id, perl = TRUE)
  unknown_study <- !grepl("^(GSE[0-9]+|E-(?:MTAB|MEXP|TABM|SGRP)-[0-9]+)$", annotations$study)
  annotations$study[unknown_study] <- paste0("independent_", which(unknown_study))
  study_counts <- table(annotations$study)
  weights <- 1 / as.numeric(study_counts[annotations$study])
  weights <- weights / sum(weights)
  responses <- sweep(responses, 2, colMeans(responses, na.rm = TRUE), "-")
  ranks <- t(apply(responses, 1, rank, na.last = "keep", ties.method = "average"))
  ranks <- t(ranks)
  valid <- is.finite(ranks)
  cat(sum(!valid), "missing CytoSig gene-profile measurements; pairwise weighting retained.\n")
  ranks[!valid] <- 0
  if (all(valid)) {
    centered <- sweep(ranks, 2, colSums(ranks * weights), "-")
    covariance <- crossprod(centered * sqrt(weights))
    scale <- sqrt(pmax(diag(covariance), 1e-30))
    cytosig_correlation <- covariance / outer(scale, scale)
  } else {
    denominator <- crossprod(valid * sqrt(weights))
    sums <- crossprod(ranks * weights, valid)
    squares <- crossprod(ranks^2 * weights, valid)
    covariance <- crossprod(ranks * sqrt(weights)) - sums * t(sums) / denominator
    variance <- squares - sums^2 / denominator
    cytosig_correlation <- covariance / sqrt(pmax(variance, 0) * pmax(t(variance), 0))
    cytosig_correlation[!is.finite(cytosig_correlation)] <- 0
  }
  cytosig_correlation <- pmax(pmin(cytosig_correlation, 1), -1)
  diag(cytosig_correlation) <- 1
  dimnames(cytosig_correlation) <- list(genes, genes)
  tahoe_fisher <- atanh(pmax(pmin(tahoe_correlation, 0.999), -0.999))
  cytosig_fisher <- atanh(pmax(pmin(cytosig_correlation, 0.999), -0.999))
  joint <- tanh((tahoe_fisher + cytosig_fisher) / 2)
  diag(joint) <- 1
  stopifnot(all(is.finite(joint)), identical(rownames(joint), genes))
  result <- list(genes = genes, tahoe = tahoe_correlation,
    cytosig = cytosig_correlation, joint = joint, cytosig_profiles = annotations)
  save_cache(result, fingerprint, cache_path)
}
write.table(result$cytosig_profiles, file.path(output_dir, "cytosig_profiles.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Comparing each square independently avoids keeping two full banks in memory.
comparisons <- lapply(c("tahoe", "cytosig", "joint"), function(source) {
  path <- file.path("reference/results/program_discovery",
    paste0("discovery_", source, "_correlation.npy"))
  if (!file.exists(path)) return(NULL)
  original <- read_npy_matrix(path)
  difference <- max(abs(result[[source]] - original))
  data.frame(source = source, maximum_absolute_difference = difference)
})
comparisons <- bind_rows(comparisons)
write.table(comparisons, file.path(output_dir, "reference_comparison.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
print(comparisons)
if (nrow(comparisons)) stopifnot(all(comparisons$maximum_absolute_difference < 1e-6))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
