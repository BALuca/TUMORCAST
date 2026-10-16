#!/usr/bin/env Rscript
# Rebuild controlled molecular contrasts from cached primary Tahoe matrices and
# released CytoSig signatures. No functional outcomes or learned bank are read.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
})
source("/work/src/program_pilot/config.R")
set.seed(seed)

# Reapply the documented expression/biotype eligibility rules. The annotation
# table is an input resource, not a sensitivity-selected list of genes.
audit <- fread(gene_audit_path)
check(all(c("gene", "protein_coding", "median_TPM", "fraction_ge_1_TPM") %in% names(audit)),
  "Gene eligibility resource has the wrong schema")
check(!anyDuplicated(audit$gene), "Duplicate gene symbols in eligibility resource")
tahoe_genes <- fread(file.path(tahoe_dir, "genes.tsv"))[order(gene_index)]
cyto <- fread(cytosig_path)
setnames(cyto, 1L, "gene")
audit[, admitted := protein_coding & median_TPM >= 5 & fraction_ge_1_TPM >= 0.8 &
  gene %chin% tahoe_genes$feature & gene %chin% cyto$gene]
genes <- sort(audit[admitted == TRUE, gene])
check(length(genes) > 1000L, "Unexpectedly small expressed shared coding universe")
fwrite(audit, file.path(data_dir, "gene_admission.tsv.gz"), sep = "\t")
gene_index <- match(genes, tahoe_genes$feature)

# The metadata includes plate-matched count normalization records. Rebuild the
# actual contrasts from primary shard counts; never use the previous response fit.
metadata <- fread(tahoe_metadata_path)
metadata[, drug := canonical_drug(drug)]
metadata[, excluded := depmap_id == excluded_cell | drug %chin% held_drugs]
fwrite(metadata[excluded == TRUE], file.path(data_dir, "excluded_tahoe_profiles.tsv.gz"), sep = "\t")
metadata <- metadata[excluded == FALSE]
check(!anyDuplicated(metadata$observation_id), "Duplicate Tahoe observations")
check(all(is.finite(metadata$treated_total) & metadata$treated_total > 0 &
  is.finite(metadata$reference_total) & metadata$reference_total > 0 &
  is.finite(metadata$treated_correction_factor) & metadata$treated_correction_factor > 0),
  "Missing or invalid plate normalization records")

reference_metadata <- fread(file.path(tahoe_dir, "pdex_plate_reference_metadata.tsv"))
reference_metadata[, key := paste(depmap_id, plate, sep = "__")]
check(!anyDuplicated(reference_metadata$key), "Duplicate Tahoe reference plate")
reference_all <- read_npy_matrix(file.path(tahoe_dir, "pdex_plate_reference.npy"))
reference <- reference_all[, gene_index, drop = FALSE]
rownames(reference) <- reference_metadata$key
check(nrow(reference_all) == nrow(reference_metadata), "Reference matrix and metadata disagree")
rm(reference_all)
invisible(gc())

tahoe <- matrix(NA_real_, nrow(metadata), length(genes),
  dimnames = list(metadata$observation_id, genes))
normalization_checks <- list()
for (shard in unique(metadata$shard_path)) {
  rows <- which(metadata$shard_path == shard)
  raw <- read_npy_matrix(file.path(cache_dir, shard))
  check(ncol(raw) == nrow(tahoe_genes), "Tahoe shard gene dimension changed")
  selected <- raw[metadata$shard_row[rows] + 1L, gene_index, drop = FALSE]
  for (j in seq_along(rows)) {
    i <- rows[j]
    plate_keys <- paste(metadata$depmap_id[i],
      trimws(strsplit(metadata$contributing_plates[i], ",", fixed = TRUE)[[1]]), sep = "__")
    check(all(plate_keys %in% rownames(reference)), "Unmatched Tahoe control plate")
    ref <- colMeans(reference[plate_keys, , drop = FALSE])
    treated <- pmax(selected[j, ] + ref, 0)
    treated_cp <- treated / metadata$treated_total[i] * 10000 / metadata$treated_correction_factor[i]
    reference_cp <- ref / metadata$reference_total[i] * 10000
    tahoe[i, ] <- log2(1 + treated_cp) - log2(1 + reference_cp)
  }
  normalization_checks[[shard]] <- data.table(shard = shard, profiles = length(rows),
    finite = all(is.finite(tahoe[rows, ])))
  message("Prepared ", length(rows), " profiles from ", basename(shard))
  rm(raw, selected)
  invisible(gc())
}
check(all(is.finite(tahoe)), "Nonfinite Tahoe molecular contrast")
fwrite(rbindlist(normalization_checks), file.path(data_dir, "tahoe_normalization_checks.tsv"), sep = "\t")
fwrite(metadata, file.path(data_dir, "tahoe_profile_metadata.tsv.gz"), sep = "\t")

# Independent ligand signatures are balanced by study, not by number of profiles.
cytosig <- t(as.matrix(cyto[match(genes, gene), -1L]))
colnames(cytosig) <- genes
cyto_metadata <- fread(cytosig_metadata_path, skip = 1L, header = FALSE,
  col.names = c("profile_id", "treatment", "condition", "duration", "dose"))
metadata_index <- match(rownames(cytosig), cyto_metadata$profile_id)
# The release lacks metadata for 54 matrix profiles. Exclude and report them;
# do not silently invent treatment/context annotations from an alternate source.
missing_metadata <- is.na(metadata_index)
fwrite(data.table(profile_id = rownames(cytosig)[missing_metadata],
  exclusion_reason = "profile absent from released meta_info.txt"),
  file.path(data_dir, "cytosig_missing_metadata.tsv"), sep = "\t")
cytosig <- cytosig[!missing_metadata, , drop = FALSE]
cyto_metadata <- cyto_metadata[metadata_index[!missing_metadata]]
check(identical(rownames(cytosig), cyto_metadata$profile_id), "CytoSig metadata join failed")
exclude <- apply(cyto_metadata, 1L, function(row) {
  any(grepl("AU[ _-]?565|ACH-000248|lapatinib|neratinib", row, ignore.case = TRUE))
})
fwrite(cyto_metadata[exclude], file.path(data_dir, "excluded_cytosig_profiles.tsv"), sep = "\t")
cytosig <- cytosig[!exclude, , drop = FALSE]
cyto_metadata <- cyto_metadata[!exclude]
cyto_metadata[, study := sub(".*@((?:GSE[0-9]+|E-(?:MTAB|MEXP|TABM|SGRP)-[0-9]+)).*",
  "\\1", profile_id, perl = TRUE)]
cyto_metadata[!grepl("^(GSE[0-9]+|E-)", study), study := paste0("independent_", .I)]
fwrite(cyto_metadata, file.path(data_dir, "cytosig_profile_metadata.tsv.gz"), sep = "\t")
fwrite(data.table(source = c("Tahoe", "CytoSig"),
  profiles = c(nrow(tahoe), nrow(cytosig)), genes = length(genes),
  missing_values = c(sum(!is.finite(tahoe)), sum(!is.finite(cytosig))),
  balancing = c("equal total weight per drug", "equal total weight per study")),
  file.path(data_dir, "molecular_source_coverage.tsv"), sep = "\t")

qsave(list(tahoe = tahoe, cytosig = cytosig, tahoe_metadata = metadata,
  cytosig_metadata = cyto_metadata, genes = genes, reference = reference),
  file.path(data_dir, "molecular_inputs.qs"), preset = "fast")
write_manifest("01_prepare_molecular",
  c(gene_audit_path, tahoe_metadata_path, cytosig_path, cytosig_metadata_path,
    file.path(tahoe_dir, "genes.tsv"), file.path(tahoe_dir, "pdex_plate_reference.npy"),
    file.path(cache_dir, unique(metadata$shard_path))),
  list(genes = length(genes), tahoe_profiles = nrow(tahoe), cytosig_profiles = nrow(cytosig),
    excluded_drugs = held_drugs, functional_labels_used = FALSE,
    contrast = "log2(1 + corrected counts per 10000) minus matched plate reference"))
message("Molecular preparation complete: ", nrow(tahoe), " Tahoe profiles and ", length(genes), " genes")
