#!/usr/bin/env Rscript
# Transport reconstructed HGF responses onto original recipient RNA and the frozen bank.
# This rebuilds the ligand molecular inputs; rebuilding bank discovery is a separate stage.

suppressPackageStartupMessages({
  library(qs)
  library(digest)
  library(dplyr)
})
source("R/prediction.R")

expression_path <- "data/raw/depmap/CCLE_expression_full.csv"
bank_path <- "data/processed/program_bank/frozen_program_weights.tsv"
source_path <- "data/processed/hgf_sources.qs"
cache_path <- "data/processed/hgf_program_inputs.qs"
output_dir <- "results/11_transport_hgf_responses"
programs <- sprintf("P%02d", 1:40)
target_cells <- c(
  "A2058", "C32", "COLO829", "G361", "SKMEL28", "SKMEL5",
  "SKBR3", "MCF7", "BT474", "ZR751", "HCC1569", "HCC1419", "AU565", "HCC1954"
)
minimum_cp10k <- 0.1
fold_change_trim <- c(0.3, 0.7)
abundance_trim <- c(0.05, 0.95)
comparison_tolerance <- 1e-10
RhpcBLASctl::blas_set_num_threads(1L)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
sources <- qread(source_path)
source_code_hash <- digest("scripts/09_reconstruct_hgf_sources.R", file = TRUE, algo = "sha256")
if (sources$fingerprint$code_hash != source_code_hash) stop("HGF source cache is stale; rerun reconstruction")
current_source_hashes <- vapply(
  names(sources$fingerprint$input_hashes), digest, character(1), file = TRUE, algo = "sha256"
)
if (!identical(current_source_hashes, sources$fingerprint$input_hashes)) {
  stop("HGF source data changed; rerun source reconstruction before transport")
}

input_paths <- c(expression_path, bank_path, source_path, "data/raw/depmap/Model.csv",
  "data/raw/cytosig/diff.merge.gz", "data/raw/cytosig/meta_info.txt",
  "data/processed/her2_ligand_outcomes.tsv")
if (!all(file.exists(input_paths))) stop("Missing required HGF transport input")
fingerprint <- list(
  inputs = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256"),
  code = digest("scripts/11_transport_hgf_responses.R", file = TRUE, algo = "sha256"),
  prediction_code = digest("R/prediction.R", file = TRUE, algo = "sha256"),
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
  ligand_baselines <- list()
  transported_genes <- list()
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
    ligand_baselines[[cells$depmap_id[i]]] <- baseline
    reference <- abundance / sum(abundance) * 1e4
    profile_index <- match(baseline$gene, rownames(profiles))
    program_index <- match(rownames(W), baseline$gene)
    stopifnot(!anyNA(program_index), all(is.finite(reference)), sum(abundance) > 0)

    projected <- matrix(NA_real_, ncol(profiles), ncol(W), dimnames = list(colnames(profiles), programs))
    for (j in seq_len(ncol(profiles))) {
      transformed <- transport_gene_response(baseline,
        setNames(profiles[, j], rownames(profiles)), W,
        minimum_cp10k, fold_change_trim, abundance_trim)
      projected[j, ] <- transformed$programs
      transported_genes[[paste(cells$depmap_id[i], colnames(profiles)[j], sep = "|")]] <- transformed$genes
      normalization_parts[[paste(i, j)]] <- data.frame(
        CellLine = cells$CellLine[i],
        source_profile = colnames(profiles)[j],
        baseline_genes = nrow(baseline),
        genes_without_source_response = transformed$missing_source_genes,
        normalization_genes = transformed$normalization_genes,
        correction_factor = transformed$correction
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
    program_weights = W,
    recipient_baselines = ligand_baselines,
    transported_hgf_gene_responses = transported_genes
  )

  # Rebuild NRG1 and EGF consensus RNA responses from the same CytoSig release.
  # Average replicates within context, contexts within study, then studies equally.
  connection <- gzfile("data/raw/cytosig/diff.merge.gz", "rt")
  profile_names <- strsplit(readLines(connection, n = 1), "\t", fixed = TRUE)[[1]]
  close(connection)
  profile_annotations <- read.delim("data/raw/cytosig/meta_info.txt", skip = 1, header = FALSE,
    col.names = c("profile_id", "treatment", "condition", "duration", "dose"),
    quote = "\"", check.names = FALSE)
  all_ligand_metadata <- data.frame(profile_id = profile_names) %>%
    left_join(profile_annotations, by = "profile_id", relationship = "one-to-one") %>%
    mutate(treatment = ifelse(is.na(treatment) | !nzchar(treatment),
      sub("[&@].*", "", profile_id), treatment)) %>%
    mutate(treatment = toupper(trimws(treatment)),
      study = sub(".*@((?:GSE[0-9]+|E-(?:MTAB|MEXP|TABM|SGRP)-[0-9]+)).*",
        "\\1", profile_id, perl = TRUE)) %>%
    filter(treatment %in% c("NRG1", "EGF", "HGF"),
      !grepl("AU[ _-]?565", paste(profile_id, condition), ignore.case = TRUE)) %>%
    mutate(condition = ifelse(is.na(condition) | condition == "", "Unspecified context", condition),
      context_key = paste(treatment, study, condition, sep = "::"))
  ligand_metadata <- all_ligand_metadata %>%
    filter(treatment %in% c("NRG1", "EGF"),
      !grepl("HCC[ _-]?1954", paste(profile_id, condition), ignore.case = TRUE))
  column_names <- c("gene", profile_names)
  column_classes <- rep("NULL", length(column_names))
  column_classes[1] <- "character"
  column_classes[match(all_ligand_metadata$profile_id, column_names)] <- "numeric"
  ligand_source <- read.delim("data/raw/cytosig/diff.merge.gz", skip = 1, header = FALSE,
    col.names = column_names, colClasses = column_classes, quote = "", check.names = FALSE)
  ligand_fold_changes <- vapply(c("NRG1", "EGF"), function(ligand) {
    contexts <- unique(ligand_metadata$context_key[ligand_metadata$treatment == ligand])
    context_responses <- vapply(contexts, function(context) {
      columns <- ligand_metadata$profile_id[ligand_metadata$context_key == context]
      response <- rowMeans(as.matrix(ligand_source[, columns, drop = FALSE]), na.rm = TRUE)
      response[!is.finite(response)] <- 0
      response
    }, numeric(nrow(ligand_source)))
    context_studies <- ligand_metadata$study[match(contexts, ligand_metadata$context_key)]
    study_responses <- vapply(unique(context_studies), function(study) {
      rowMeans(context_responses[, context_studies == study, drop = FALSE])
    }, numeric(nrow(ligand_source)))
    rowMeans(study_responses)
  }, numeric(nrow(ligand_source)))
  rownames(ligand_fold_changes) <- ligand_source$gene
  assay_outcomes <- read.delim("data/processed/her2_ligand_outcomes.tsv")
  recipients <- unique(assay_outcomes[, c("CellLine", "depmap_id", "gene")])
  ligand_results <- vector("list", nrow(recipients))
  ligand_gene_results <- vector("list", nrow(recipients))
  for (i in seq_len(nrow(recipients))) {
    baseline <- ligand_baselines[[recipients$depmap_id[i]]]
    transformed <- transport_gene_response(baseline,
      setNames(ligand_fold_changes[, recipients$gene[i]], rownames(ligand_fold_changes)), W,
      minimum_cp10k, fold_change_trim, abundance_trim)
    projected <- matrix(transformed$programs, nrow = 1, dimnames = list(NULL, programs))
    ligand_results[[i]] <- cbind(recipients[i, ], as.data.frame(projected))
    ligand_gene_results[[i]] <- transformed$genes
  }
  names(ligand_gene_results) <- paste(recipients$depmap_id, recipients$gene, sep = "|")
  result$ligand_program_inputs <- bind_rows(ligand_results)
  result$ligand_fold_changes <- ligand_fold_changes
  result$ligand_source_metadata <- ligand_metadata
  result$ligand_gene_responses <- ligand_gene_results

  # Also reproduce the original 12-context HGF conversion used by the earlier
  # route-assumption diagnostic. It is distinct from the corrected study mean.
  original_metadata <- all_ligand_metadata %>% filter(treatment == "HGF")
  original_contexts <- unique(original_metadata$condition)
  original_hgf_genes <- vapply(original_contexts, function(context) {
    columns <- original_metadata$profile_id[original_metadata$condition == context]
    rowMeans(as.matrix(ligand_source[, columns, drop = FALSE]), na.rm = TRUE)
  }, numeric(nrow(ligand_source)))
  rownames(original_hgf_genes) <- ligand_source$gene
  original_hgf_programs <- t(vapply(original_contexts, function(context) {
    transport_gene_response(ligand_baselines[["ACH-000248"]],
      setNames(original_hgf_genes[, context], rownames(original_hgf_genes)), W,
      minimum_cp10k, fold_change_trim, abundance_trim)$programs
  }, numeric(length(programs))))
  colnames(original_hgf_programs) <- programs
  original_hgf_programs <- rbind(original_hgf_programs,
    `Equal-context mean` = colMeans(original_hgf_programs))
  result$original_hgf_contexts <- data.frame(context = rownames(original_hgf_programs),
    original_hgf_programs, check.names = FALSE)
  result$original_hgf_gene_profiles <- original_hgf_genes
}

# The old object is read only for comparison, never to generate these values.
reference_bank <- readRDS("reference/inputs/original_molecular_functional_inputs.rds")$prep$W
bank_difference <- max(abs(result$program_weights[rownames(reference_bank), programs] - reference_bank))
stopifnot(bank_difference < 1e-5)
accepted <- read.delim("reference/results/131_time_matched_hgf_program_inputs.tsv")
rebuilt <- result$program_inputs
reference_keys <- paste(accepted$depmap_id, accepted$source_profile, sep = "|")
rebuilt_keys <- paste(rebuilt$depmap_id, rebuilt$source_profile, sep = "|")
stopifnot(!anyDuplicated(reference_keys), !anyDuplicated(rebuilt_keys), setequal(reference_keys, rebuilt_keys))
index <- match(reference_keys, rebuilt_keys)
difference <- max(abs(as.matrix(rebuilt[index, programs]) - as.matrix(accepted[, programs])))
stopifnot(difference < 1e-5, nrow(rebuilt) == 126L)

# Isolate transport from small numerical changes in the independently rebuilt
# bank: the original coordinates must still reproduce the original profiles.
reference_scored <- t(vapply(result$transported_hgf_gene_responses, function(displacement) {
  as.numeric(displacement[rownames(reference_bank)] %*% reference_bank)
}, numeric(length(programs))))
gene_keys <- names(result$transported_hgf_gene_responses)
reference_coordinate_difference <- max(abs(reference_scored[match(reference_keys, gene_keys), ] -
  as.matrix(accepted[, programs])))
stopifnot(reference_coordinate_difference < comparison_tolerance)

# Check every required NRG1/EGF profile against the historical ligand analysis.
accepted_ligands <- read.delim("reference/results/ligand_inputs/external_ligand_programs.tsv")
ligand_keys <- paste(result$ligand_program_inputs$depmap_id,
  result$ligand_program_inputs$gene, sep = "|")
accepted_ligand_keys <- paste(accepted_ligands$depmap_id, accepted_ligands$gene, sep = "|")
ligand_index <- match(ligand_keys, accepted_ligand_keys)
stopifnot(!anyNA(ligand_index), !anyDuplicated(ligand_keys))
ligand_difference <- max(abs(as.matrix(result$ligand_program_inputs[, programs]) -
  as.matrix(accepted_ligands[ligand_index, programs])))
stopifnot(ligand_difference < 1e-5)
ligand_reference_scored <- t(vapply(result$ligand_gene_responses, function(displacement) {
  as.numeric(displacement[rownames(reference_bank)] %*% reference_bank)
}, numeric(length(programs))))
ligand_transport_difference <- max(abs(ligand_reference_scored -
  as.matrix(accepted_ligands[ligand_index, programs])))
stopifnot(ligand_transport_difference < comparison_tolerance)

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
write.table(result$ligand_program_inputs, file.path(output_dir, "her2_ligand_program_inputs.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(result$ligand_source_metadata, file.path(output_dir, "her2_ligand_source_profiles.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(result$original_hgf_contexts, "data/processed/hgf_original_context_programs.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(
    component = c("generated_program_weights", "transported_hgf_program_inputs",
      "hgf_transport_in_original_coordinates", "nrg1_egf_program_inputs",
      "nrg1_egf_transport_in_original_coordinates"),
    rows = c(nrow(result$program_weights), nrow(rebuilt), nrow(rebuilt),
      length(ligand_keys), length(ligand_keys)),
    maximum_absolute_difference = c(bank_difference, difference, reference_coordinate_difference,
      ligand_difference, ligand_transport_difference)
  ),
  file.path(output_dir, "comparison_with_original.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("All 126 HGF profiles across 14 recipients match the original 40-program inputs.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
sink()
