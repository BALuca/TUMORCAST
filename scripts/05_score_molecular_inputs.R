#!/usr/bin/env Rscript
# Reconstruct the selected 42-recipient, 35-drug molecular panel from Tahoe means,
# then project it and all CCLE baselines onto the independently regenerated bank.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(qs)
  library(digest)
  library(arrow)
  library(jsonlite)
})
source("R/data_io.R")

raw_dir <- "data/raw/tahoe"
selected_dir <- file.path(raw_dir, "pdex_selected")
effect_cache_path <- "data/processed/molecular_gene_effects.qs"
output_dir <- "results/05_score_molecular_inputs"
bank_path <- "data/processed/program_bank/frozen_program_weights.tsv"
native_doses <- c(0.05, 0.5, 5)
excluded_panel_cells <- c("ACH-000389", "ACH-000437", "ACH-000558")
# The historical ALMANAC/Tahoe intersection defines the intervention panel.
drug_names <- c(
  methotrexate = "Methotrexate", busulfan = "Busulfan", altretamine = "Altretamine",
  hydroxyurea = "Hydroxyurea", vinblastine = "Vinblastine (sulfate)",
  cytarabine = "Cytarabine (hydrochloride)", vincristine = "vincristine",
  azacitidine = "5-Azacytidine", doxorubicin = "Doxorubicin (hydrochloride)",
  paclitaxel = "Paclitaxel", decitabine = "Decitabine", bendamustine = "Bendamustine",
  sirolimus = "Rapamycin", idarubicin = "Idarubicin (hydrochloride)",
  oxaliplatin = "Oxaliplatin", mitoxantrone = "Mitoxantrone (dihydrochloride)",
  imiquimod = "Imiquimod (hydrochloride)", clofarabine = "Clofarabine",
  topotecan = "Topotecan (hydrochloride)", gemcitabine = "Gemcitabine",
  docetaxel = "Docetaxel", bortezomib = "Bortezomib", gefitinib = "Gefitinib",
  erlotinib = "Erlotinib", fulvestrant = "Fulvestrant", anastrozole = "Anastrozole",
  everolimus = "Everolimus", lapatinib = "Lapatinib ditosylate",
  lenalidomide = "Lenalidomide (hemihydrate)", sunitinib = "Sunitinib",
  vemurafenib = "Vemurafenib", vismodegib = "Vismodegib", crizotinib = "crizotinib",
  vandetanib = "Vandetanib", plicamycin = "Plicamycin"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
expression <- qread("data/processed/expression.qs")$result
genes <- expression$genes
cells <- setdiff(expression$eligibility_cells, excluded_panel_cells)
gene_table <- fread(file.path(selected_dir, "genes.tsv")) %>% arrange(gene_index)
gene_indices <- match(genes, gene_table$feature)
stopifnot(length(cells) == 42L, length(drug_names) == 35L, !anyNA(gene_indices))

# Input hashes include only allocated raw files/ranges, without expanding sparse
# Parquet holes. Each downloaded range has its original checksum in its JSON.
source_paths <- c("data/processed/expression.qs",
  list.files(selected_dir, pattern = "\\.(npy|tsv|tsv.gz)$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(raw_dir, "pdex_projection_chunks"), full.names = TRUE),
  file.path(raw_dir, "pdex_fetched_ranges.tsv"), file.path(raw_dir, "pdex_source_version.json"))
fingerprint <- list(
  inputs = vapply(source_paths, digest, character(1), file = TRUE, algo = "sha256"),
  code = vapply(c("scripts/05_score_molecular_inputs.R", "R/data_io.R"), digest,
    character(1), file = TRUE, algo = "sha256"),
  cells = cells, drugs = drug_names, doses = native_doses, R = R.version.string,
  sparse_files = file.info(file.path(raw_dir, c("pdex_selected_sparse.parquet",
    "pdex_missing_rowgroups.parquet")))[, c("size", "mtime")]
)
cached <- if (file.exists(effect_cache_path)) qread(effect_cache_path) else NULL
if (!is.null(cached) && identical(cached$fingerprint, fingerprint)) {
  cat("Reusing matching raw molecular-effect cache.\n")
  molecular <- cached$result
} else {
  cat("Rebuilding molecular effects from source means.\n")
  normalize <- function(treated, reference) {
    tolerance <- max(1e-6, 32 * 2^-23 * max(abs(c(treated, reference)), 1))
    stopifnot(all(is.finite(treated)), all(is.finite(reference)), min(treated, reference) >= -tolerance)
    treated <- pmax(treated, 0)
    reference <- pmax(reference, 0)
    treated_total <- sum(treated)
    reference_total <- sum(reference)
    stopifnot(treated_total > 0, reference_total > 0)
    treated_cp10k <- treated / treated_total * 10000
    reference_cp10k <- reference / reference_total * 10000
    eligible <- treated_cp10k >= 0.1 & reference_cp10k >= 0.1
    log_ratio <- log2(treated_cp10k[eligible]) - log2(reference_cp10k[eligible])
    log_abundance <- 0.5 * (log2(treated_cp10k[eligible]) + log2(reference_cp10k[eligible]))
    ratio_bounds <- quantile(log_ratio, c(0.3, 0.7))
    abundance_bounds <- quantile(log_abundance, c(0.05, 0.95))
    kept <- log_ratio >= ratio_bounds[1] & log_ratio <= ratio_bounds[2] &
      log_abundance >= abundance_bounds[1] & log_abundance <= abundance_bounds[2]
    stopifnot(sum(kept) > 0)
    correction <- 2^median(log_ratio[kept])
    list(effect = log2(1 + treated_cp10k[gene_indices] / correction) -
      log2(1 + reference_cp10k[gene_indices]), treated_total = treated_total,
      reference_total = reference_total, correction = correction)
  }
  availability <- fread(file.path(selected_dir, "response_availability.tsv"))
  conditions <- fread(file.path(selected_dir, "perturbations.tsv"))
  metadata <- availability %>%
    filter(depmap_id %in% cells, available == 1) %>%
    left_join(select(conditions, condition_index, drug, dose), by = "condition_index") %>%
    filter(drug %in% drug_names, dose %in% native_doses) %>%
    arrange(depmap_id, drug, dose)
  plates <- fread(file.path(selected_dir, "pdex_plate_audit.tsv.gz")) %>%
    mutate(key = paste(depmap_id, perturbation_id, sep = "__"))
  plate_lists <- split(plates$plate, plates$key)
  reference_metadata <- fread(file.path(selected_dir, "pdex_plate_reference_metadata.tsv")) %>%
    arrange(array_row)
  references <- read_npy_matrix(file.path(selected_dir, "pdex_plate_reference.npy"))
  rownames(references) <- paste(reference_metadata$depmap_id, reference_metadata$plate, sep = "__")
  identities <- fread(file.path(raw_dir, "public_tahoe_cell_metadata.tsv")) %>%
    distinct(cell_line = Cell_ID_Cellosaur, depmap_id = Cell_ID_DepMap)
  chunk_files <- list.files(file.path(raw_dir, "pdex_projection_chunks"), full.names = TRUE)
  source_index <- bind_rows(lapply(chunk_files, function(path) {
    chunk <- fromJSON(path)
    if (!length(chunk$profiles)) return(NULL)
    as.data.frame(chunk$profiles) %>%
      select(cell_line, plate, drug, dose) %>%
      distinct() %>% mutate(row_group = chunk$rowgroup)
  })) %>%
    left_join(identities, by = "cell_line") %>%
    filter(depmap_id %in% cells, drug %in% drug_names, dose %in% native_doses) %>%
    mutate(key = paste(depmap_id, drug, dose, sep = "|"))
  cached_plate_keys <- paste(plates$depmap_id, plates$drug, plates$dose, plates$plate, sep = "|")
  source_index <- source_index %>%
    filter(!paste(depmap_id, drug, dose, plate, sep = "|") %in% cached_plate_keys)
  pooled <- list()

  effects <- records <- list()
  for (id in unique(metadata$depmap_id)) {
    rows <- metadata[metadata$depmap_id == id, ]
    shard <- read_npy_matrix(file.path(selected_dir, "response_shards",
      paste0(id, "__", rows$cell_line[1], ".npy")))
    for (i in seq_len(nrow(rows))) {
      selected_plates <- plate_lists[[paste(id, rows$perturbation_id[i], sep = "__")]]
      reference <- colMeans(references[paste(id, selected_plates, sep = "__"), , drop = FALSE])
      treated <- shard[rows$condition_index[i] + 1L, ] + reference
      normalized <- normalize(treated, reference)
      key <- paste(id, rows$drug[i], rows$dose[i], sep = "|")
      if (key %in% source_index$key) {
        # Reconstruct the sum of means for the cached plates, then add any
        # additional source plates before estimating the composition offset.
        pooled[[key]] <- list(treated = treated * length(selected_plates),
          reference = reference * length(selected_plates), plates = selected_plates,
          metadata = data.frame(depmap_id = id, drug = rows$drug[i], dose = rows$dose[i]))
      }
      effects[[key]] <- normalized$effect
      records[[key]] <- data.frame(depmap_id = id, drug = rows$drug[i], dose = rows$dose[i],
        contributing_plates = paste(sort(selected_plates), collapse = ","),
        source = "released selected Tahoe means", treated_total = normalized$treated_total,
        reference_total = normalized$reference_total, correction = normalized$correction)
    }
    cat(id, ":", nrow(rows), "selected-release profiles\n")
  }
  rm(shard, references)

  # Remaining recipients/profiles are read directly from validated downloaded
  # Parquet row groups. Columns are zero-based in Arrow's row-group interface.
  reader <- ParquetFileReader$create(file.path(raw_dir, "pdex_selected_sparse.parquet"))
  groups <- unique(source_index$row_group)
  for (g in seq_along(groups)) {
    group_id <- groups[g]
    requested <- source_index[source_index$row_group == group_id, ]
    table <- as.data.frame(reader$ReadRowGroup(as.integer(group_id), c(0L, 1L, 2L, 3L, 11L, 12L)))
    table$drug <- sub("^\\[\\('(.+)', [0-9.eE+-]+, 'uM'\\)\\]$", "\\1", table$target)
    table$dose <- as.numeric(sub(".*', ([0-9.eE+-]+),.*", "\\1", table$target))
    for (i in seq_len(nrow(requested))) {
      r <- requested[i, ]
      profile <- table[table$cell_line == r$cell_line & table$plate == r$plate &
        table$drug == r$drug & table$dose == r$dose, ]
      stopifnot(nrow(profile) == nrow(gene_table), !anyDuplicated(profile$feature))
      index <- match(gene_table$feature, profile$feature)
      stopifnot(!anyNA(index))
      treated <- profile$target_mean[index]
      reference <- profile$ref_mean[index]
      key <- r$key
      if (is.null(pooled[[key]])) {
        pooled[[key]] <- list(treated = treated, reference = reference, plates = r$plate,
          metadata = r[, c("depmap_id", "drug", "dose")])
      } else {
        stopifnot(!r$plate %in% pooled[[key]]$plates)
        pooled[[key]]$treated <- pooled[[key]]$treated + treated
        pooled[[key]]$reference <- pooled[[key]]$reference + reference
        pooled[[key]]$plates <- c(pooled[[key]]$plates, r$plate)
      }
    }
    if (g %% 20L == 0L) cat(g, "/", length(groups), "source row groups read\n")
  }
  for (key in names(pooled)) {
    profile <- pooled[[key]]
    normalized <- normalize(profile$treated, profile$reference)
    effects[[key]] <- normalized$effect
    records[[key]] <- cbind(profile$metadata,
      data.frame(contributing_plates = paste(sort(profile$plates), collapse = ","),
        source = "public Tahoe means pooled across plates", treated_total = normalized$treated_total,
        reference_total = normalized$reference_total, correction = normalized$correction))
  }
  metadata <- bind_rows(records)
  metadata$sample_id <- paste0(metadata$depmap_id, "__", metadata$drug, "__", metadata$dose, "uM")
  effect_matrix <- do.call(rbind, effects)
  colnames(effect_matrix) <- genes
  stopifnot(identical(names(records), rownames(effect_matrix)), !anyDuplicated(metadata$sample_id))
  molecular <- list(metadata = metadata, effects = effect_matrix, genes = genes, drugs = drug_names)
  save_cache(molecular, fingerprint, effect_cache_path)
}
coverage <- molecular$metadata %>% count(drug, dose, name = "recipients")
write.table(coverage, file.path(output_dir, "molecular_coverage.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
if (nrow(molecular$metadata) != 4410L || any(coverage$recipients != 42L)) {
  stop("Source molecular panel is incomplete; inspect coverage before fitting")
}

# Use the accepted coordinates only to isolate the raw-preparation comparison.
# They do not generate the new bank, molecular effects or downstream model fits.
reference_path <- "reference/inputs/original_molecular_functional_inputs.rds"
if (file.exists(reference_path)) {
  original <- readRDS(reference_path)$prep
  programs <- colnames(original$W)
  reference_projection <- molecular$effects %*% original$W[molecular$genes, programs]
  original_key <- paste(original$m$depmap_id, original$m$drug, original$m$dose, sep = "|")
  rebuilt_key <- paste(molecular$metadata$depmap_id, molecular$metadata$drug, molecular$metadata$dose, sep = "|")
  index <- match(original_key, rebuilt_key)
  stopifnot(!anyNA(index), !anyDuplicated(rebuilt_key))
  difference <- apply(abs(reference_projection[index, ] - as.matrix(original$m[, programs])), 1, max)
  comparison <- data.frame(original$m[, c("depmap_id", "drug", "dose", "source")],
    maximum_program_difference = difference)
  write.table(comparison, file.path(output_dir, "raw_preparation_comparison.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Maximum raw-preparation difference in accepted coordinates:", max(difference), "\n")
  stopifnot(max(difference) < 1e-5)
  rm(original, reference_projection)
}
if (!file.exists(bank_path)) stop("Run 04_generate_programs.R before molecular projection")
weights <- read.delim(bank_path)
programs <- sort(unique(weights$program))
W <- matrix(0, length(genes), length(programs), dimnames = list(genes, programs))
W[cbind(match(weights$gene, genes), match(weights$program, programs))] <- weights$weight
RhpcBLASctl::blas_set_num_threads(4)
scores <- molecular$effects %*% W
baseline_scores <- expression$baseline[, genes, drop = FALSE] %*% W
metadata_columns <- c("depmap_id", "drug", "dose", "contributing_plates", "source", "sample_id")
prepared <- list(
  m = cbind(molecular$metadata[, metadata_columns], as.data.frame(scores)),
  raw = data.frame(depmap_id = rownames(baseline_scores), baseline_scores, check.names = FALSE),
  receptors = data.frame(depmap_id = rownames(expression$baseline),
    expression$baseline[, c("ERBB2", "MET"), drop = FALSE], check.names = FALSE),
  W = W, program_version = "all_eligible_program_first"
)
save_cache(prepared, list(inputs = c(fingerprint$inputs,
  bank = digest(bank_path, file = TRUE, algo = "sha256")), code = fingerprint$code),
  "data/processed/molecular_inputs.qs")
write.table(prepared$m, file.path(output_dir, "program_displacements.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(prepared$raw, file.path(output_dir, "baseline_program_scores.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(prepared$receptors, file.path(output_dir, "receptors.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
