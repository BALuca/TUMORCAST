#!/usr/bin/env Rscript
# Verify the accepted bundles and recover their original figure files unchanged.
# Run from the TUMORCAST project root. This script does not refit or redraw figures.

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
