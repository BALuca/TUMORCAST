#!/usr/bin/env Rscript
# Build source-specific and balanced gene correlations. Whole development cell
# lines and CytoSig studies are reserved for molecular rank selection.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
})
source("/work/src/program_pilot/config.R")
input_path <- file.path(data_dir, "molecular_inputs.qs")
check(file.exists(input_path), "Run 01_prepare_molecular.R first")
input <- qread(input_path)

# Fixed splits use identifiers only. The external functional test stays separate.
set.seed(seed)
line_ids <- sample(sort(unique(input$tahoe_metadata$depmap_id)))
line_fold <- setNames(rep(seq_len(n_molecular_folds), length.out = length(line_ids)), line_ids)
study_ids <- sample(sort(unique(input$cytosig_metadata$study)))
study_fold <- setNames(rep(seq_len(n_molecular_folds), length.out = length(study_ids)), study_ids)
input$tahoe_metadata[, molecular_fold := unname(line_fold[depmap_id])]
input$cytosig_metadata[, molecular_fold := unname(study_fold[study])]
fwrite(input$tahoe_metadata, file.path(data_dir, "tahoe_molecular_split.tsv.gz"), sep = "\t")
fwrite(input$cytosig_metadata, file.path(data_dir, "cytosig_molecular_split.tsv.gz"), sep = "\t")

# Remove the global displacement within each profile only for correlation
# discovery. Raw signed contrasts remain intact for maps, scoring, and validation.
# This matches the initial bank's relative co-response definition explicitly.
tahoe_relative <- sweep(input$tahoe, 1L, rowMeans(input$tahoe), "-")
cyto_relative <- sweep(input$cytosig, 1L, rowMeans(input$cytosig, na.rm = TRUE), "-")
tahoe_train <- input$tahoe_metadata$molecular_fold != 1L
cyto_train <- input$cytosig_metadata$molecular_fold != 1L
message("Building Tahoe discovery correlation from ", sum(tahoe_train), " profiles")
tahoe_correlation <- weighted_gene_correlation(tahoe_relative[tahoe_train, , drop = FALSE],
  group_weights(input$tahoe_metadata$drug[tahoe_train]))
message("Building independent CytoSig discovery correlation")
cytosig_correlation <- weighted_gene_correlation(cyto_relative[cyto_train, , drop = FALSE],
  group_weights(input$cytosig_metadata$study[cyto_train]))

# Equal source contribution in Fisher-z space. Clipping is only for the atanh
# singularity; source weighting is independent of the number of signatures.
clip <- function(x) pmax(-0.999999, pmin(0.999999, x))
joint_correlation <- tanh((atanh(clip(tahoe_correlation)) +
  atanh(clip(cytosig_correlation))) / 2)
diag(joint_correlation) <- 1
dimnames(joint_correlation) <- list(input$genes, input$genes)
check(all(is.finite(joint_correlation)), "Nonfinite joint correlation")

qsave(list(tahoe = tahoe_correlation, cytosig = cytosig_correlation,
  joint = joint_correlation, line_fold = line_fold, study_fold = study_fold,
  genes = input$genes), file.path(data_dir, "discovery_correlations.qs"), preset = "fast")
fwrite(data.table(source = c("Tahoe", "CytoSig"),
  discovery_profiles = c(sum(tahoe_train), sum(cyto_train)),
  excluded_development_profiles = c(sum(!tahoe_train), sum(!cyto_train)),
  discovery_units = c(uniqueN(input$tahoe_metadata$depmap_id[tahoe_train]),
    uniqueN(input$cytosig_metadata$study[cyto_train])),
  excluded_development_units = c(uniqueN(input$tahoe_metadata$depmap_id[!tahoe_train]),
    uniqueN(input$cytosig_metadata$study[!cyto_train]))),
  file.path(data_dir, "molecular_discovery_split_counts.tsv"), sep = "\t")
write_manifest("02_build_correlations", input_path,
  list(source_weight = c(Tahoe = 0.5, CytoSig = 0.5),
    correlation = "weighted Spearman; profiles centered over admitted genes; Fisher-z source average",
    selection_split = "whole Tahoe cell lines and whole CytoSig studies; fold 1 validation",
    functional_labels_used = FALSE))
message("Discovery correlations complete")
