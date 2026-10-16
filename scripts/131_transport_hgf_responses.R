#!/usr/bin/env Rscript
# Transport reconstructed HGF responses onto original recipient RNA and the frozen bank.
# This rebuilds the ligand molecular inputs; rebuilding bank discovery is a separate stage.

suppressPackageStartupMessages({
  library(qs)
  library(digest)
  library(dplyr)
})

expression_path <- "data/raw/depmap/CCLE_expression_full.csv"
bank_path <- "data/processed/program_bank/frozen_program_weights.tsv"
source_path <- "data/processed/hgf_sources.qs"
cache_path <- "data/processed/hgf_program_inputs.qs"
output_dir <- "results/131_transport_hgf_responses"
programs <- sprintf("P%02d", 1:40)
target_cells <- c(
  "A2058", "C32", "COLO829", "G361", "SKMEL28", "SKMEL5",
  "SKBR3", "MCF7", "BT474", "ZR751", "HCC1569", "HCC1419", "AU565", "HCC1954"
)
minimum_cp10k <- 0.1
fold_change_trim <- c(0.3, 0.7)
abundance_trim <- c(0.05, 0.95)
comparison_tolerance <- 1e-10

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
sources <- qread(source_path)
source_code_hash <- digest("scripts/130_reconstruct_hgf_sources.R", file = TRUE, algo = "sha256")
if (sources$fingerprint$code_hash != source_code_hash) stop("HGF source cache is stale; rerun reconstruction")
current_source_hashes <- vapply(
  names(sources$fingerprint$input_hashes), digest, character(1), file = TRUE, algo = "sha256"
)
if (!identical(current_source_hashes, sources$fingerprint$input_hashes)) {
  stop("HGF source data changed; rerun source reconstruction before transport")
}

input_paths <- c(expression_path, bank_path, source_path, "data/raw/depmap/Model.csv")
if (!all(file.exists(input_paths))) stop("Missing required HGF transport input")
fingerprint <- list(
  inputs = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256"),
  code = digest("scripts/131_transport_hgf_responses.R", file = TRUE, algo = "sha256"),
  target_cells = target_cells,
  minimum_cp10k = minimum_cp10k,
  fold_change_trim = fold_change_trim,
  abundance_trim = abundance_trim,
  R_version = R.version.string,
  packages = vapply(c("dplyr", "qs", "digest"), function(package) {
    as.character(packageVersion(package))
  }, character(1))
)
reuse_cache <- FALSE
if (file.exists(cache_path)) {
  cached <- qread(cache_path)
  reuse_cache <- identical(cached$fingerprint, fingerprint)
  if (!reuse_cache) cat("HGF transport cache is stale: inputs, code, settings or package versions changed.\n")
}

if (reuse_cache) {
  cat("Reusing matching HGF transport cache.\n")
  result <- cached$result
} else {
  # Recover the same signed, hard-assigned 40-program weight matrix from its table.
  weights <- read.delim(bank_path)
  stopifnot(all(c("gene", "program", "weight") %in% names(weights)))
  stopifnot(!anyDuplicated(weights$gene), setequal(unique(weights$program), programs))
  genes <- sort(unique(weights$gene))
  W <- matrix(0, length(genes), length(programs), dimnames = list(genes, programs))
  W[cbind(match(weights$gene, genes), match(weights$program, programs))] <- weights$weight
  stopifnot(nrow(W) == 8842L, all(is.finite(W)))

  # Original recipient choices include calibration cells and explicit query cells.
  # AU565/HCC1954 are projected here; this script does not train or select a model.
  models <- read.csv("data/raw/depmap/Model.csv", check.names = FALSE)
  cells <- data.frame(
    CellLine = target_cells,
    depmap_id = models$ModelID[match(target_cells, models$StrippedCellLineName)]
  )
  stopifnot(!anyNA(cells$depmap_id), !anyDuplicated(cells$depmap_id))

  # Stream the original full CCLE release and retain only the 14 required rows.
  connection <- file(expression_path, "rt")
  header <- readLines(connection, n = 1)
  expression_rows <- list()
  repeat {
    lines <- readLines(connection, n = 100, warn = FALSE)
    if (!length(lines)) break
    ids <- gsub('"', "", sub(",.*", "", lines))
    for (i in which(ids %in% cells$depmap_id)) expression_rows[[ids[i]]] <- lines[i]
    if (length(expression_rows) == nrow(cells)) break
  }
  close(connection)
  stopifnot(setequal(names(expression_rows), cells$depmap_id))
  expression <- read.csv(
    text = paste(c(header, unlist(expression_rows, use.names = FALSE)), collapse = "\n"),
    check.names = FALSE
  )
  names(expression)[1] <- "depmap_id"
  stopifnot(!anyDuplicated(expression$depmap_id))

  # Preserve study-specific and leave-one-study-out profiles, with equal study means.
  gene_profiles <- sources$result$corrected_gene_profiles
  studies <- colnames(gene_profiles)
  profiles <- cbind(
    Primary_study_mean = rowMeans(gene_profiles, na.rm = TRUE),
    gene_profiles,
    vapply(studies, function(study) {
      rowMeans(gene_profiles[, studies != study, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(gene_profiles)))
  )
  colnames(profiles) <- c("Primary_study_mean", studies, paste0("Leave_out_", studies))
  # Historical source-only projections assign no change where a source lacks a gene.
  # Keep that assumption explicit; these zeros are not measured gene responses.
  profiles[!is.finite(profiles)] <- 0

  output_parts <- list()
  normalization_parts <- list()
  for (i in seq_len(nrow(cells))) {
    cat("Transporting HGF into", cells$CellLine[i], "\n")
    row <- expression[expression$depmap_id == cells$depmap_id[i], , drop = FALSE]
    stopifnot(nrow(row) == 1L)
    symbols <- sub(" \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$", "", names(row)[-1], perl = TRUE)
    baseline <- data.frame(
      gene = symbols,
      abundance = pmax(2^as.numeric(row[1, -1]) - 1, 0)
    ) |>
      summarise(abundance = mean(abundance), .by = gene)
    abundance <- baseline$abundance
    reference <- abundance / sum(abundance) * 1e4
    profile_index <- match(baseline$gene, rownames(profiles))
    program_index <- match(rownames(W), baseline$gene)
    stopifnot(!anyNA(program_index), all(is.finite(reference)), sum(abundance) > 0)

    projected <- matrix(NA_real_, ncol(profiles), ncol(W), dimnames = list(colnames(profiles), programs))
    for (j in seq_len(ncol(profiles))) {
      fold_change <- profiles[profile_index, j]
      missing_source_genes <- sum(!is.finite(fold_change))
      fold_change[!is.finite(fold_change)] <- 0
      treated <- abundance * 2^fold_change
      treated <- treated / sum(treated) * 1e4
      eligible <- treated >= minimum_cp10k & reference >= minimum_cp10k
      log_ratio <- log2(treated[eligible] / reference[eligible])
      mean_log_abundance <- (log2(treated[eligible]) + log2(reference[eligible])) / 2
      ratio_limits <- quantile(log_ratio, fold_change_trim)
      abundance_limits <- quantile(mean_log_abundance, abundance_trim)
      central <- log_ratio >= ratio_limits[1] & log_ratio <= ratio_limits[2] &
        mean_log_abundance >= abundance_limits[1] & mean_log_abundance <= abundance_limits[2]
      stopifnot(any(central))
      correction_factor <- 2^median(log_ratio[central])
      displacement <- log2(1 + treated / correction_factor) - log2(1 + reference)
      projected[j, ] <- as.numeric(displacement[program_index] %*% W)
      normalization_parts[[paste(i, j)]] <- data.frame(
        CellLine = cells$CellLine[i],
        source_profile = colnames(profiles)[j],
        baseline_genes = nrow(baseline),
        genes_without_source_response = missing_source_genes,
        normalization_genes = sum(central),
        correction_factor = correction_factor
      )
    }
    stopifnot(all(is.finite(projected)))
    output_parts[[i]] <- cbind(
      cells[rep(i, nrow(projected)), ],
      data.frame(gene = "HGF", source_profile = rownames(projected)),
      as.data.frame(projected)
    )
  }
  result <- list(
    program_inputs = bind_rows(output_parts),
    normalization = bind_rows(normalization_parts),
    program_weights = W
  )
}

# The old object is read only for comparison, never to generate these values.
reference_bank <- readRDS("reference/inputs/original_molecular_functional_inputs.rds")$prep$W
bank_difference <- max(abs(result$program_weights[rownames(reference_bank), programs] - reference_bank))
stopifnot(bank_difference < comparison_tolerance)
accepted <- read.delim("reference/results/131_time_matched_hgf_program_inputs.tsv")
rebuilt <- result$program_inputs
reference_keys <- paste(accepted$depmap_id, accepted$source_profile, sep = "|")
rebuilt_keys <- paste(rebuilt$depmap_id, rebuilt$source_profile, sep = "|")
stopifnot(!anyDuplicated(reference_keys), !anyDuplicated(rebuilt_keys), setequal(reference_keys, rebuilt_keys))
index <- match(reference_keys, rebuilt_keys)
difference <- max(abs(as.matrix(rebuilt[index, programs]) - as.matrix(accepted[, programs])))
stopifnot(difference < comparison_tolerance, nrow(rebuilt) == 126L)

if (!reuse_cache) {
  temporary_cache <- tempfile(pattern = "hgf_transport_", tmpdir = dirname(cache_path))
  qsave(list(fingerprint = fingerprint, result = result), temporary_cache)
  if (!file.rename(temporary_cache, cache_path)) stop("Cannot publish completed HGF transport cache")
}
write.table(
  rebuilt, file.path(output_dir, "hgf_program_inputs.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  result$normalization, file.path(output_dir, "normalization.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  data.frame(
    component = c("frozen_program_weights", "transported_hgf_program_inputs"),
    rows = c(nrow(result$program_weights), nrow(rebuilt)),
    maximum_absolute_difference = c(bank_difference, difference)
  ),
  file.path(output_dir, "comparison_with_original.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("All 126 HGF profiles across 14 recipients match the original 40-program inputs.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
sink()
