#!/usr/bin/env Rscript
# Rebuild HGF gene responses from original GEO matrices and the CytoSig release.
# Preserve the native contrasts, normalization and precision weighting of analyses 130/131.

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(digest)
  library(qs)
  library(dplyr)
})

# Source files are original releases. Responses are HGF-treated minus untreated
# log2 expression. No functional outcomes enter this molecular reconstruction.
raw_dir <- "data/raw"
output_dir <- "results/09_reconstruct_hgf_sources"
cache_path <- "data/processed/hgf_sources.qs"
studies <- c("GSE26732", "GSE26735")
platforms <- c(GSE26732 = "GPL96", GSE26735 = "GPL2725")
control_columns <- 1:2
hgf_columns <- 9:10
comparison_tolerance <- 1e-10

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
cat("HGF source reconstruction\n")

# Store a small freshness record with the expensive result itself.
input_paths <- c(
  file.path(raw_dir, "cytosig", c("diff.merge.gz", "meta_info.txt")),
  file.path(raw_dir, "hgf_geo", paste0(studies, "_series_matrix.txt.gz")),
  file.path(raw_dir, "hgf_geo", paste0(unname(platforms), ".annot.gz")),
  file.path(raw_dir, "annotations/org.Hs.eg.sqlite")
)
if (!all(file.exists(input_paths))) stop("Missing HGF source files: ", paste(input_paths[!file.exists(input_paths)]))
packages <- c("dplyr", "limma", "AnnotationDbi", "qs")
fingerprint <- list(
  input_hashes = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256"),
  code_hash = digest("scripts/09_reconstruct_hgf_sources.R", file = TRUE, algo = "sha256"),
  packages = vapply(packages, function(package) as.character(packageVersion(package)), character(1)),
  R_version = R.version.string,
  control_columns = control_columns,
  hgf_columns = hgf_columns,
  illumina_background_offset = 16
)
reuse_cache <- FALSE
if (file.exists(cache_path)) {
  cached <- qread(cache_path)
  reuse_cache <- identical(cached$fingerprint, fingerprint)
  if (!reuse_cache) cat("Cache is stale: source data, code, settings or package versions changed.\n")
}

if (reuse_cache) {
  cat("Reusing matching HGF source cache.\n")
  result <- cached$result
} else {
  cat("Reconstructing HGF source responses from original data.\n")

  # Recover the original HGF profile selection from the complete CytoSig metadata.
  metadata <- read.delim(
    file.path(raw_dir, "cytosig/meta_info.txt"),
    skip = 1,
    header = FALSE,
    col.names = c("profile_id", "treatment", "condition", "duration", "dose"),
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  metadata <- metadata |>
    filter(toupper(treatment) == "HGF") |>
    mutate(
      study = sub(".*@((?:GSE[0-9]+|E-(?:MTAB|MEXP|TABM|SGRP)-[0-9]+)).*", "\\1", profile_id, perl = TRUE),
      hours = suppressWarnings(as.numeric(sub("[^0-9.].*$", "", duration))),
      # The original study methods specify 16 h for GSE135552 despite its 48 h CytoSig label.
      hours = if_else(study == "GSE135552", 16, hours),
      study_family = if_else(study %in% studies, "GSE26736", study),
      cancer_cell = grepl("cancer|H1975|Detroit|Fadu|Scc9|Pc3m|Mda.Mb435", condition, ignore.case = TRUE),
      eligible = cancer_cell & is.finite(hours) & hours >= 12 & hours <= 48 &
        !grepl("AU[ _-]?565|HCC[ _-]?1954", paste(profile_id, condition), ignore.case = TRUE),
      context_key = paste(study_family, condition, sep = "|")
    )
  selected <- metadata |>
    filter(eligible) |>
    mutate(distance_to_24h = abs(log(hours / 24))) |>
    filter(distance_to_24h <= min(distance_to_24h) + 1e-10, .by = context_key)
  metadata$selected <- metadata$profile_id %in% selected$profile_id
  stopifnot(n_distinct(selected$study_family) == 4L, n_distinct(selected$context_key) == 6L)
  stopifnot(nrow(selected) == 7L, !anyDuplicated(selected$profile_id))

  # Base R reads only the seven needed numeric columns from the full release.
  connection <- gzfile(file.path(raw_dir, "cytosig/diff.merge.gz"), "rt")
  # The release header omits the gene-column name; data rows contain that extra field.
  header <- c("gene", strsplit(readLines(connection, n = 1), "\t", fixed = TRUE)[[1]])
  close(connection)
  selected_columns <- match(selected$profile_id, header)
  stopifnot(!anyNA(selected_columns))
  column_types <- rep("NULL", length(header))
  column_types[1] <- "character"
  column_types[selected_columns] <- "numeric"
  connection <- gzfile(file.path(raw_dir, "cytosig/diff.merge.gz"), "rt")
  cytosig <- read.delim(
    connection,
    header = FALSE,
    skip = 1,
    col.names = header,
    row.names = NULL,
    colClasses = column_types,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    na.strings = c("NA", "NaN", ""),
    stringsAsFactors = FALSE
  )
  close(connection)
  names(cytosig)[1] <- "gene"
  stopifnot(!anyDuplicated(cytosig$gene), nrow(cytosig) == 19918L)

  # Retain context averaging followed by equal weighting of independent studies.
  context_names <- unique(selected$context_key)
  context_responses <- vapply(context_names, function(context) {
    profiles <- selected$profile_id[selected$context_key == context]
    values <- rowMeans(as.matrix(cytosig[, profiles, drop = FALSE]), na.rm = TRUE)
    # Preserve the historical rule for an entirely unmeasured context/gene.
    values[!is.finite(values)] <- 0
    values
  }, numeric(nrow(cytosig)))
  context_studies <- selected$study_family[match(context_names, selected$context_key)]
  study_names <- sort(unique(context_studies))
  original_profiles <- vapply(study_names, function(study) {
    rowMeans(context_responses[, context_studies == study, drop = FALSE])
  }, numeric(nrow(cytosig)))
  rownames(original_profiles) <- cytosig$gene

  # Reconstruct native 24 h HGF versus untreated wild-type contrasts on both platforms.
  response_parts <- list()
  diagnostic_parts <- list()
  sample_titles <- list()
  for (study in studies) {
    cat("Reading native source:", study, "\n")
    connection <- gzfile(file.path(raw_dir, "hgf_geo", paste0(study, "_series_matrix.txt.gz")), "rt")
    lines <- readLines(connection, warn = FALSE)
    close(connection)
    sample_titles[[study]] <- lines[grepl("^!Sample_title|^!Sample_geo_accession", lines)]
    begin <- grep("^!series_matrix_table_begin", lines)
    end <- grep("^!series_matrix_table_end", lines)
    stopifnot(length(begin) == 1L, length(end) == 1L, begin < end)
    native <- read.delim(text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"), check.names = FALSE)
    expression <- as.matrix(native[, -1, drop = FALSE])
    stopifnot(ncol(expression) >= max(hgf_columns), all(is.finite(expression)))

    connection <- gzfile(file.path(raw_dir, "hgf_geo", paste0(platforms[[study]], ".annot.gz")), "rt")
    lines <- readLines(connection, warn = FALSE)
    close(connection)
    begin <- grep("^!platform_table_begin", lines, ignore.case = TRUE)
    end <- grep("^!platform_table_end", lines, ignore.case = TRUE)
    stopifnot(length(begin) == 1L, length(end) == 1L)
    annotation <- read.delim(
      text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"),
      check.names = FALSE,
      colClasses = "character"
    )
    entrez_ids <- sub(" .*", "", sub("///.*", "", annotation[["Gene ID"]]))
    valid_ids <- unique(entrez_ids[grepl("^[0-9]+$", entrez_ids)])
    gene_symbols <- unname(mapIds(
      AnnotationDbi::loadDb(file.path(raw_dir, "annotations/org.Hs.eg.sqlite")),
      keys = valid_ids,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    ))
    annotation$gene <- gene_symbols[match(entrez_ids, valid_ids)]
    genes <- annotation$gene[match(native$ID_REF, annotation$ID)]

    # Affymetrix is already RMA log2. Illumina contains linear, background-corrected
    # signals including negatives. Keep both its historical floor-log diagnostic
    # and the accepted normexp(offset=16) plus quantile-normalized log2 response.
    if (study == "GSE26732") {
      expressions <- list(native_RMA_log2 = expression)
    } else {
      corrected <- backgroundCorrect(new("EListRaw", list(E = expression)), method = "normexp", offset = 16)$E
      expressions <- list(
        released_signal_floor_log2 = log2(pmax(expression, 0) + 1),
        normexp_quantile_log2 = normalizeBetweenArrays(log2(corrected), method = "quantile")
      )
    }

    for (method in names(expressions)) {
      values <- expressions[[method]]
      responses <- data.frame(
        gene = genes,
        logFC = rowMeans(values[, hgf_columns]) - rowMeans(values[, control_columns]),
        WT_difference = values[, 1] - values[, 2],
        HGF_difference = values[, 9] - values[, 10],
        replicate1_logFC = values[, 9] - values[, 1],
        replicate2_logFC = values[, 10] - values[, 2],
        control_mean = rowMeans(values[, control_columns])
      ) |>
        filter(!is.na(gene), gene != "") |>
        summarise(across(where(is.numeric), mean), .by = gene)
      profile <- selected$profile_id[selected$study == study]
      stopifnot(length(profile) == 1L)
      responses$cached_CytoSig_logFC <- cytosig[[profile]][match(responses$gene, cytosig$gene)]
      responses$study <- study
      responses$preprocessing <- method
      response_parts[[paste(study, method)]] <- responses
      diagnostic_parts[[paste(study, method)]] <- data.frame(
        study = study,
        preprocessing = method,
        genes = nrow(responses),
        response_SD = sd(responses$logFC),
        WT_difference_SD = sd(responses$WT_difference),
        HGF_difference_SD = sd(responses$HGF_difference),
        expected_mean_contrast_noise_SD = sqrt((var(responses$WT_difference) + var(responses$HGF_difference)) / 4),
        replicate_response_r = cor(responses$replicate1_logFC, responses$replicate2_logFC),
        cached_CytoSig_r = cor(responses$logFC, responses$cached_CytoSig_logFC, use = "complete.obs")
      )
    }
  }
  native_responses <- bind_rows(response_parts)
  diagnostics <- bind_rows(diagnostic_parts)

  # Weight the two native platforms by inverse replicate-contrast noise variance.
  precision <- diagnostics |>
    filter(
      (study == "GSE26732" & preprocessing == "native_RMA_log2") |
        (study == "GSE26735" & preprocessing == "normexp_quantile_log2")
    ) |>
    mutate(
      inverse_variance = 1 / expected_mean_contrast_noise_SD^2,
      platform_weight = inverse_variance / sum(inverse_variance)
    )
  measurements <- native_responses |>
    inner_join(dplyr::select(precision, study, preprocessing, inverse_variance), by = c("study", "preprocessing")) |>
    summarise(
      logFC = sum(inverse_variance * logFC) / sum(inverse_variance),
      native_platforms = n(),
      .by = gene
    )
  corrected_profiles <- original_profiles
  corrected_profiles[, "GSE26736"] <- measurements$logFC[match(rownames(original_profiles), measurements$gene)]
  stopifnot(all(is.finite(precision$platform_weight)), abs(sum(precision$platform_weight) - 1) < 1e-12)

  result <- list(
    original_gene_profiles = original_profiles,
    corrected_gene_profiles = corrected_profiles,
    native_gene_responses = native_responses,
    platform_diagnostics = diagnostics,
    precision = precision,
    eligibility = metadata,
    sample_titles = sample_titles,
    missing_annotation_policy = "Unmapped native annotations remain NA; study means use available observations."
  )
  temporary_cache <- tempfile(pattern = "hgf_sources_", tmpdir = dirname(cache_path))
  qsave(list(fingerprint = fingerprint, result = result), temporary_cache)
  if (!file.rename(temporary_cache, cache_path)) stop("Cannot publish completed HGF source cache")
}

# Compare every native gene contrast against accepted historical values.
reference <- read.delim("reference/results/130_hgf_sources/reconstructed_platform_gene_responses.tsv")
comparison <- inner_join(
  result$native_gene_responses,
  reference,
  by = c("gene", "study", "preprocessing"),
  suffix = c("_rebuilt", "_reference"),
  relationship = "one-to-one"
)
stopifnot(nrow(comparison) == nrow(reference), nrow(comparison) == nrow(result$native_gene_responses))
quantities <- c(
  "logFC", "WT_difference", "HGF_difference", "replicate1_logFC", "replicate2_logFC",
  "control_mean", "cached_CytoSig_logFC"
)
differences <- vapply(quantities, function(quantity) {
  rebuilt <- comparison[[paste0(quantity, "_rebuilt")]]
  accepted <- comparison[[paste0(quantity, "_reference")]]
  stopifnot(identical(is.na(rebuilt), is.na(accepted)))
  max(abs(rebuilt - accepted), na.rm = TRUE)
}, numeric(1))
stopifnot(all(differences < comparison_tolerance))

reference_profiles <- qread("reference/results/131_hgf_source_reconstruction.qs")
stopifnot(identical(dimnames(result$corrected_gene_profiles), dimnames(reference_profiles$corrected_gene_profiles)))
stopifnot(identical(is.na(result$corrected_gene_profiles), is.na(reference_profiles$corrected_gene_profiles)))
source_difference <- max(abs(result$corrected_gene_profiles - reference_profiles$corrected_gene_profiles), na.rm = TRUE)
stopifnot(source_difference < comparison_tolerance)

# Keep the scientifically useful summaries and one compact numerical preservation check.
write.table(
  result$platform_diagnostics,
  file.path(output_dir, "platform_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  result$precision,
  file.path(output_dir, "platform_precision_weights.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  result$eligibility,
  file.path(output_dir, "hgf_source_eligibility.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  data.frame(
    quantity = c(quantities, "corrected_source_profiles"),
    maximum_absolute_difference = c(differences, source_difference)
  ),
  file.path(output_dir, "comparison_with_original.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
print(result$precision[, c("study", "preprocessing", "genes", "replicate_response_r", "platform_weight")])
cat("Native gene contrasts and corrected HGF profiles reproduce the accepted analysis.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
sink()

# Inspect the matched native HGF replicates after successful reconstruction.
 #!/usr/bin/env Rscript
# Inspect agreement of the two matched native HGF-control replicate contrasts.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(digest)
})

output_dir <- "results/09_reconstruct_hgf_sources"
cache <- qread("data/processed/hgf_sources.qs")
source_hash <- digest("scripts/09_reconstruct_hgf_sources.R", file = TRUE, algo = "sha256")
if (cache$fingerprint$code_hash != source_hash) stop("HGF cache is stale; rerun its reconstruction first")
current_input_hashes <- vapply(
  names(cache$fingerprint$input_hashes), digest, character(1), file = TRUE, algo = "sha256"
)
if (!identical(current_input_hashes, cache$fingerprint$input_hashes)) {
  stop("HGF source data changed; rerun reconstruction before plotting")
}

# Show the two accepted preprocessing choices. Each point is one mapped gene.
points <- cache$result$native_gene_responses |>
  filter(
    (study == "GSE26732" & preprocessing == "native_RMA_log2") |
      (study == "GSE26735" & preprocessing == "normexp_quantile_log2")
  ) |>
  mutate(platform = if_else(study == "GSE26732", "Affymetrix: native RMA", "Illumina: normexp + quantile"))
stopifnot(all(is.finite(points$replicate1_logFC)), all(is.finite(points$replicate2_logFC)))

# Compute a stat_cor equivalent on exactly the observations in each panel.
platforms <- unique(points$platform)
statistics <- lapply(platforms, function(platform) {
  rows <- points[points$platform == platform, ]
  test <- cor.test(rows$replicate1_logFC, rows$replicate2_logFC, method = "pearson")
  p_label <- if (test$p.value < 2.2e-16) "P < 2.2e-16" else paste0("P = ", format.pval(test$p.value, digits = 2))
  data.frame(
    platform = platform,
    label = paste0("r = ", sprintf("%.2f", unname(test$estimate)), "\n", p_label, "\nn = ", nrow(rows))
  )
}) |>
  bind_rows()

plot <- ggplot(points, aes(replicate1_logFC, replicate2_logFC)) +
  geom_point(size = 0.35, alpha = 0.2, color = "#2568A0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7, color = "#BF3434") +
  geom_text(
    data = statistics,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.15,
    size = 3.8
  ) +
  facet_wrap(vars(platform), scales = "free", nrow = 1) +
  labs(
    title = "Native HGF response: replicate agreement",
    x = "Replicate 1: HGF - control (log2)",
    y = "Replicate 2: HGF - control (log2)",
    caption = "Each point is a gene. Pearson tests are descriptive across correlated genes."
  ) +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "plain"),
    plot.caption = element_text(size = 10, hjust = 0),
    panel.spacing = grid::unit(7, "mm")
  )

ggsave(file.path(output_dir, "hgf_replicate_agreement.pdf"), plot, width = 7.8, height = 4.4)
ggsave("tmp/hgf_replicate_agreement.png", plot, width = 7.8, height = 4.4, dpi = 150)
message("Saved the native HGF replicate diagnostic with square panels and independent axes.")
