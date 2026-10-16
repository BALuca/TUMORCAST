#!/usr/bin/env Rscript
# Verify the accepted bundles and recover their original figure files unchanged.
# Run after stages 01-16. Preserve originals byte-for-byte, then compare the
# independently rebuilt measurements, fitted predictions and figure coverage.

library(digest)

bundle_dir <- "reference/bundles"
figure_dir <- "reference/figures"
expected_hashes <- c(
  "TUMORCAST_molecular_bridge_and_ligand_rescue.tar.gz" =
    "3219d1cc2e2d7261ab93f39049b3cc4843c6191b6fdf76c9db1dcaaf1793f2d7",
  "TUMORCAST_ggplot2_panels_v2 (8).zip" =
    "1163077f6f31cca75da55aa74a537b5ef805866e881aa8c4035a6c4c6931d4c3",
  "TUMORCAST_original_panels_joined.zip" =
    "7e99d708761873babdc307db0c5681ebbc9810f54b6955e3f40a305acd642f2d"
)

# Reject changed archives before extracting any reference files.
for (name in names(expected_hashes)) {
  path <- file.path(bundle_dir, name)
  if (!file.exists(path)) stop("Missing original bundle: ", path)
  observed_hash <- digest(path, file = TRUE, algo = "sha256")
  if (observed_hash != expected_hashes[[name]]) stop("Original bundle has changed: ", path)
  message("Verified bundle: ", name)
}

# Preserve the different figure arrangements in separate, descriptive locations.
destinations <- c("molecular_bridge_and_ligand_rescue", "ggplot2", "joined")
figure_counts <- integer(length(expected_hashes))
for (i in seq_along(expected_hashes)) {
  archive <- file.path(bundle_dir, names(expected_hashes)[i])
  scratch <- tempfile(pattern = "reference_figures_")
  dir.create(scratch)

  if (i == 1L) {
    members <- untar(archive, list = TRUE)
    selected <- members[grepl("[.](pdf|png|svg)$", members)]
    untar(archive, files = selected, exdir = scratch)
  } else {
    members <- unzip(archive, list = TRUE)$Name
    selected <- members[grepl("[.](pdf|png|svg)$", members)]
    unzip(archive, files = selected, exdir = scratch)
  }

  if (!length(selected)) stop("No figure files found in ", archive)
  if (any(grepl("(^/|(^|/)[.][.](/|$))", selected))) stop("Unexpected archive member path")
  for (member in selected) {
    source <- file.path(scratch, member)
    relative <- sub("^[^/]+/", "", member)
    target <- file.path(figure_dir, destinations[i], relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    source_hash <- digest(source, file = TRUE, algo = "sha256")
    if (file.exists(target)) {
      if (digest(target, file = TRUE, algo = "sha256") != source_hash) {
        stop("Preserved reference figure differs from the original: ", target)
      }
    } else {
      if (!file.copy(source, target, copy.date = TRUE)) stop("Cannot preserve figure: ", target)
    }
    stopifnot(digest(target, file = TRUE, algo = "sha256") == source_hash)
  }
  figure_counts[i] <- length(selected)
  unlink(scratch, recursive = TRUE)
  message("Verified ", length(selected), " unchanged figure files in ", destinations[i])
}

message("Reference preservation complete: ", sum(figure_counts), " figure files, zero differing bytes.")

output_dir <- "results/17_verify_reproduction"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(bundle = names(expected_hashes), preserved_figures = figure_counts,
  differing_bytes = 0L), file.path(output_dir, "original_figure_preservation.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Compare measurements and predictions using their biological identifiers.
# This is deliberately separate from recovery of the accepted image bytes.
comparisons <- list(
  signature_decoder = list(
    rebuilt = "results/07_learn_response_signatures/data/measured_predictions_frozen.tsv",
    original = paste0("reference/results/95_leave_one_cell_out_correlation_signatures/",
      "measured_predictions_frozen.tsv"),
    keys = c("depmap_id", "drug", "dose", "evaluation"),
    values = c("observed_sensitivity", "predicted_deviation"), tolerance = 0.001),
  signature_baseline_bridge = list(
    rebuilt = "results/07_learn_response_signatures/data/baseline_predictions_frozen.tsv",
    original = paste0("reference/results/95_leave_one_cell_out_correlation_signatures/",
      "baseline_predictions_frozen.tsv"),
    keys = c("depmap_id", "drug", "dose", "evaluation"),
    values = c("observed_sensitivity", "predicted_deviation"), tolerance = 0.001),
  measured_single_drug = list(
    rebuilt = "results/08_fit_single_drug_bridge/data/measured_predictions_frozen.tsv",
    original = "reference/results/106_single_drug_bridge/measured_predictions_frozen.tsv",
    keys = c("depmap_id", "drug", "dose", "evaluation"),
    values = c("observed_log2", "predicted_log2"), tolerance = 0.001),
  baseline_single_drug = list(
    rebuilt = "results/08_fit_single_drug_bridge/data/baseline_predictions_frozen.tsv",
    original = "reference/results/106_single_drug_bridge/baseline_predictions_frozen.tsv",
    keys = c("depmap_id", "drug", "dose", "evaluation"),
    values = c("observed_log2", "predicted_log2"), tolerance = 0.001),
  drug_pairs = list(
    rebuilt = "results/08_fit_single_drug_bridge/data/combination_predictions_frozen.tsv",
    original = "reference/results/106_single_drug_bridge/combination_predictions_frozen.tsv",
    keys = c("depmap_id", "pair", "SCREENER", "doseA", "doseB", "evaluation"),
    values = c("observed_log2", "predicted_log2", "predicted_departure_log2"), tolerance = 0.001),
  ligand_rescue = list(
    rebuilt = "results/12_fit_ligand_rescue/data/combination_predictions_frozen.tsv",
    original = "reference/results/144_ligand_rescue/combination_predictions_frozen.tsv",
    keys = c("depmap_id", "gene", "dose", "evaluation"),
    values = c("observed_departure_log2", "predicted_departure_log2"), tolerance = 0.001),
  hgf_query = list(
    rebuilt = "results/12_fit_ligand_rescue/data/AU565_HGF_predictions_frozen.tsv",
    original = "reference/results/144_ligand_rescue/AU565_HGF_predictions_frozen.tsv",
    keys = c("CellLine", "dose", "source_profile"),
    values = c("lapatinib_pct_control", "predicted_combination_pct_control"), tolerance = 0.1),
  original_hgf_route_query = list(
    rebuilt = "results/14_plot_model_diagnostics/au565_hgf_5uM_zero_reference_query.tsv",
    original = "reference/results/106_single_drug_bridge/au565_hgf_5uM_zero_reference_query.tsv",
    keys = c("hgf_context", "assumed_HGF_route"),
    values = c("lapatinib_pct_zero_control", "hgf_alone_pct_zero_control",
      "combination_pct_zero_control", "lapatinib_in_hgf_pct_hgf_control"), tolerance = 0.1)
)
results <- lapply(names(comparisons), function(name) {
  specification <- comparisons[[name]]
  if (!file.exists(specification$rebuilt)) stop("Missing rebuilt result: ", specification$rebuilt)
  rebuilt <- read.delim(specification$rebuilt, check.names = FALSE)
  original <- read.delim(specification$original, check.names = FALSE)
  key <- function(table) do.call(paste, c(table[specification$keys], sep = "|"))
  rebuilt_key <- key(rebuilt)
  original_key <- key(original)
  stopifnot(!anyDuplicated(rebuilt_key), !anyDuplicated(original_key),
    setequal(rebuilt_key, original_key))
  index <- match(original_key, rebuilt_key)
  do.call(rbind, lapply(specification$values, function(column) {
    difference <- max(abs(rebuilt[[column]][index] - original[[column]]))
    data.frame(analysis = name, column = column, observations = nrow(original),
      maximum_absolute_difference = difference, tolerance = specification$tolerance,
      passed = is.finite(difference) && difference <= specification$tolerance)
  }))
})
results <- do.call(rbind, results)
write.table(results, file.path(output_dir, "numerical_reproduction.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
print(results)

# Every accepted figure has a freshly generated counterpart. Archive duplicates
# remain identifiable, including the two different B2 Figure 2 arrangements.
accepted_figures <- list.files(figure_dir, pattern = "[.](pdf|png|svg)$", recursive = TRUE, full.names = TRUE)
counterparts <- vapply(accepted_figures, function(path) {
  name <- basename(path)
  if (grepl("/ggplot2/panels/|/joined/original_panels/", path)) {
    return(file.path("results/15_plot_polished_panels/panels", name))
  }
  if (grepl("/ggplot2/composites/", path)) {
    return(file.path("results/15_plot_polished_panels/composites", name))
  }
  if (grepl("/joined/", path)) return(file.path("results/16_assemble_figures", name))
  if (grepl("/single_drug_bridge/", path)) return(file.path("results/08_fit_single_drug_bridge/figures", name))
  if (grepl("grant_support/", path)) {
    return(file.path("results/12_fit_ligand_rescue/figures/grant_support", name))
  }
  if (grepl("/ligand_rescue/", path)) return(file.path("results/12_fit_ligand_rescue/figures", name))
  stop("Unmapped accepted figure: ", path)
}, character(1))
coverage <- data.frame(accepted_figure = accepted_figures, regenerated_figure = counterparts,
  regenerated = file.exists(counterparts))
write.table(coverage, file.path(output_dir, "figure_reproduction.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
if (!all(results$passed) || !all(coverage$regenerated)) {
  stop("Reproduction is incomplete; inspect numerical_reproduction.tsv and figure_reproduction.tsv")
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
message("Raw-to-figure reconstruction verified; original figure files remain byte-for-byte unchanged.")
