#!/usr/bin/env Rscript
# Reconstruct the expression-only gene eligibility gate and CCLE baselines.
# No functional measurements or drug-response labels enter this stage.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(qs)
  library(digest)
  library(ensembldb)
})
source("R/data_io.R")

expression_path <- "data/raw/depmap/CCLE_expression_full.csv"
tahoe_genes_path <- "data/raw/tahoe/pdex_selected/genes.tsv"
cytosig_path <- "data/raw/cytosig/diff.merge.gz"
annotation_path <- "data/raw/annotations/EnsDb.Hsapiens.v86.sqlite"
cache_path <- "data/processed/expression.qs"
output_dir <- "results/01_prepare_expression"
minimum_median_tpm <- 5
minimum_tpm <- 1
minimum_fraction <- 0.8

# Original molecular cohort, fixed before functional selection. AU565 is excluded.
eligibility_cells <- c(
  "ACH-000120", "ACH-000138", "ACH-000139", "ACH-000148", "ACH-000164",
  "ACH-000178", "ACH-000222", "ACH-000236", "ACH-000348", "ACH-000389",
  "ACH-000396", "ACH-000437", "ACH-000463", "ACH-000493", "ACH-000496",
  "ACH-000521", "ACH-000552", "ACH-000555", "ACH-000558", "ACH-000580",
  "ACH-000601", "ACH-000669", "ACH-000681", "ACH-000750", "ACH-000757",
  "ACH-000790", "ACH-000793", "ACH-000842", "ACH-000853", "ACH-000861",
  "ACH-000875", "ACH-000890", "ACH-000900", "ACH-000916", "ACH-000927",
  "ACH-000932", "ACH-000940", "ACH-000943", "ACH-000950", "ACH-000954",
  "ACH-000957", "ACH-000958", "ACH-000997", "ACH-001190", "ACH-001333"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
input_paths <- c(expression_path, tahoe_genes_path, cytosig_path, annotation_path)
stopifnot(all(file.exists(input_paths)), !"ACH-000248" %in% eligibility_cells)
fingerprint <- list(
  inputs = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256"),
  code = vapply(c("scripts/01_prepare_expression.R", "R/data_io.R"), digest,
    character(1), file = TRUE, algo = "sha256"),
  cells = eligibility_cells,
  thresholds = c(minimum_median_tpm, minimum_tpm, minimum_fraction),
  R = R.version.string,
  annotation_reader = as.character(packageVersion("ensembldb"))
)
cached <- if (file.exists(cache_path)) qread(cache_path) else NULL
if (!is.null(cached) && identical(cached$fingerprint, fingerprint)) {
  cat("Reusing matching expression cache.\n")
  result <- cached$result
} else {
  cat("Rebuilding expression cache: absent or stale inputs, code or settings.\n")
  expression <- fread(expression_path, showProgress = FALSE, nThread = 4)
  cell_ids <- expression[[1]]
  source_names <- names(expression)[-1]
  symbols <- sub(" .*", "", source_names)
  ensembl_ids <- sub(".*\\((ENSG[0-9]+).*", "\\1", source_names)
  stopifnot(!anyDuplicated(cell_ids), all(eligibility_cells %in% cell_ids))
  expression_matrix <- as.matrix(expression[, -1, with = FALSE])
  rm(expression)

  # Duplicate symbols are averaged on the released log2(TPM + 1) scale,
  # matching the original baseline preprocessing, before converting to TPM.
  symbol_counts <- table(symbols)
  summed_expression <- rowsum(t(expression_matrix), symbols, reorder = TRUE)
  baseline <- t(summed_expression / as.numeric(symbol_counts[rownames(summed_expression)]))
  rownames(baseline) <- cell_ids
  rm(expression_matrix, summed_expression)
  stopifnot(ncol(baseline) == 53960L, all(is.finite(baseline)))

  annotation_database <- ensembldb::EnsDb(annotation_path)
  annotation <- as.data.frame(ensembldb::genes(
    annotation_database,
    columns = c("gene_id", "gene_name", "gene_biotype"), return.type = "DataFrame"
  ))
  mapping <- data.frame(
    gene = symbols,
    ensembl_id = ensembl_ids,
    biotype = annotation$gene_biotype[match(ensembl_ids, annotation$gene_id)]
  )
  biotypes <- mapping %>%
    group_by(gene) %>%
    summarise(
      protein_coding = all(!is.na(biotype) & biotype == "protein_coding"),
      .groups = "drop"
    )
  tpm <- pmax(2^baseline[eligibility_cells, , drop = FALSE] - 1, 0)
  eligibility <- data.frame(
    gene = colnames(tpm),
    median_TPM = matrixStats::colMedians(tpm),
    fraction_ge_1_TPM = colMeans(tpm >= minimum_tpm)
  ) %>%
    left_join(biotypes, by = "gene") %>%
    mutate(eligible = protein_coding & median_TPM >= minimum_median_tpm &
      fraction_ge_1_TPM >= minimum_fraction)

  # Require exact symbol coverage in both molecular discovery sources.
  tahoe_genes <- fread(tahoe_genes_path)$feature
  connection <- gzfile(cytosig_path, "rt")
  cytosig_header <- strsplit(readLines(connection, n = 1), "\t", fixed = TRUE)[[1]]
  cytosig_lines <- readLines(connection)
  close(connection)
  cytosig_genes <- sub("\t.*", "", cytosig_lines)
  eligibility <- eligibility %>%
    mutate(in_tahoe = gene %in% tahoe_genes, in_cytosig = gene %in% cytosig_genes,
      included = eligible & in_tahoe & in_cytosig)
  genes <- sort(eligibility$gene[eligibility$included])
  stopifnot(sum(eligibility$eligible) == 9063L, length(genes) == 8842L)
  result <- list(baseline = baseline, eligibility = eligibility, genes = genes,
    eligibility_cells = eligibility_cells, source_gene_mapping = mapping)
  save_cache(result, fingerprint, cache_path)
}
write.table(result$eligibility, file.path(output_dir, "gene_eligibility.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = result$genes), file.path(output_dir, "program_gene_universe.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
cat(nrow(result$baseline), "CCLE recipients;", length(result$genes), "admitted discovery genes.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
