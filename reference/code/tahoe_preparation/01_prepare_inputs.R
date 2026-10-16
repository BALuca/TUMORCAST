#!/usr/bin/env Rscript

## Rebuild corrected Tahoe, PRISM and DepMap analysis inputs from canonical
## local files. No accepted result table or model is read here.
plots_only <- "--plots-only" %in% commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  if (!plots_only) {
    library(nanoparquet)
    library(jsonlite)
  }
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
started <- proc.time()[["elapsed"]]
RESULT <- stage_result("01_inputs")
ARTIFACT <- stage_artifact("01_inputs")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

plot_stage01 <- function() {
  fold_file <- file.path(RESULT, "fold_assignments.tsv")
  if (!file.exists(fold_file)) fail("Stage 01 plotting requires %s.", fold_file)
  selected_lines <- fread(fold_file)
  control_metadata <- fread(file.path(CANONICAL_CACHE,
    "dmso_control_metadata.tsv"))

  ## Exact retained preparation-stage plotting block. The native column
  ## extractor only avoids loading the complete 62,710-gene matrix into R.
  gene_sample <- unique(as.integer(round(seq(1, EXPECTED_GENES,
    length.out = 200L))))
  control_expression <- read_npy_float32_columns(file.path(CANONICAL_CACHE,
    "dmso_control_expression.npy"), gene_sample)
  plot_rows <- rbindlist(lapply(sort(selected_lines$depmap_id), function(id) {
    rows <- which(control_metadata$depmap_id == id)
    rows <- rows[order(control_metadata$plate[rows])][1:2]
    source_rows <- control_metadata$array_row[rows] + 1L
    data.table(depmap_id = id,
      plate_1_expression = log1p(control_expression[source_rows[1L], ]),
      plate_2_expression = log1p(control_expression[source_rows[2L], ]))
  }))
  pair_r <- cor(plot_rows$plate_1_expression, plot_rows$plate_2_expression,
    use = "complete.obs")
  p <- ggplot(plot_rows, aes(plate_1_expression, plate_2_expression)) +
    geom_point(alpha = 0.16, size = 0.45, color = "#255F85") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 0.55, color = "#D1495B") +
    annotate("text", x = Inf, y = -Inf,
      label = sprintf("Pooled Pearson r = %.3f", pair_r),
      hjust = 1.05, vjust = -0.7, size = 3.0) +
    labs(x = "Plate 1 log1p DMSO expression",
      y = "Plate 2 log1p DMSO expression",
      title = "Tahoe plate-level DMSO reproducibility",
      subtitle = "Deterministic 200-gene sample in each of 33 retained cell lines") +
    theme_classic(base_size = 9) +
    theme(aspect.ratio = 1, plot.title = element_text(face = "bold"))
  ggsave(file.path(RESULT, "01_control_reproducibility.pdf"), p,
    width = 4.6, height = 4.6, units = "in", device = cairo_pdf)
}

if (plots_only) {
  plot_stage01()
  message("Stage 01 plots refreshed from existing inputs.")
  quit(save = "no", status = 0L)
}

clean_text <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  trimws(gsub("[[:space:]]+", " ", x))
}
name_key <- function(x) tolower(clean_text(x))
format_number <- function(x) {
  answer <- sprintf("%.12f", as.numeric(x))
  answer <- sub("0+$", "", answer); sub("\\.$", "", answer)
}

package_files <- c("manifest.json", "cells.tsv", "perturbations.tsv",
  "response_availability.tsv", "genes.tsv", "pdex_plate_reference.npy",
  "pdex_plate_reference_metadata.tsv", "pdex_plate_audit.tsv.gz")
package_paths <- file.path(TAHOE_PACKAGE, package_files)
control_paths <- file.path(CANONICAL_CACHE, c("dmso_control_expression.npy",
  "dmso_control_leave_one_plate_out.npy", "dmso_control_metadata.tsv",
  "dmso_control_genes.tsv", "response_line_sufficient_stats.npz"))
prism_paths <- file.path(PRISM_CACHE, c(
  "secondary-screen-replicate-treatment-info.csv",
  "secondary-screen-cell-line-info.csv",
  "secondary-screen-logfold-change.csv"))
depmap_path <- file.path(DATASET_CACHE, "depmap", "CCLE_expression_full.csv")
tahoe_cell_path <- file.path(DATASET_CACHE, "tahoe",
  "cell_line_metadata.parquet")
required <- c(package_paths, control_paths, prism_paths, depmap_path,
  tahoe_cell_path)
if (any(!file.exists(required)))
  fail("Missing canonical input(s): %s",
    paste(required[!file.exists(required)], collapse = ", "))

manifest <- fromJSON(package_paths[1L], simplifyVector = FALSE)
manifest_files <- rbindlist(lapply(manifest$files, as.data.table), fill = TRUE)
manifest_files[, full_path := file.path(TAHOE_PACKAGE, path)]
if (any(!file.exists(manifest_files$full_path)) ||
    any(file.info(manifest_files$full_path)$size != manifest_files$size))
  fail("The canonical Tahoe package does not match its local manifest.")

cells35 <- fread(file.path(TAHOE_PACKAGE, "cells.tsv"))
conditions <- fread(file.path(TAHOE_PACKAGE, "perturbations.tsv"))
availability <- fread(file.path(TAHOE_PACKAGE, "response_availability.tsv"))
genes <- fread(file.path(TAHOE_PACKAGE, "genes.tsv"))
if (nrow(cells35) != 35L || nrow(conditions) != EXPECTED_DISCOVERY_CONDITIONS ||
    uniqueN(conditions$drug) != EXPECTED_DISCOVERY_DRUGS ||
    nrow(genes) != EXPECTED_GENES || anyDuplicated(genes$feature))
  fail("The canonical Tahoe dimensions changed.")
if (nrow(availability) != 35L * EXPECTED_DISCOVERY_CONDITIONS)
  fail("Tahoe response availability is not a 35-line Cartesian grid.")

raw_cells <- as.data.table(read_parquet(tahoe_cell_path))
cell_annotation <- unique(raw_cells[, .(
  cell_line = clean_text(Cell_ID_Cellosaur),
  depmap_id = clean_text(Cell_ID_DepMap),
  tahoe_cell_name = clean_text(cell_name), organ = clean_text(Organ))])
cells35 <- merge(cells35[, .(depmap_id, cell_line,
  tahoe_cell_name = cell_name)],
  cell_annotation[, .(cell_line, organ)], by = "cell_line", all.x = TRUE,
  sort = FALSE)
if (anyNA(cells35$organ) || anyDuplicated(cells35$depmap_id) ||
    anyDuplicated(cells35$cell_line)) fail("Tahoe cell annotation join failed.")

## Rebuild response metadata from the uploaded primary package.
cell_order <- copy(cells35)
cell_order[, cell_index := .I - 1L]
response_metadata <- merge(availability,
  conditions[, .(condition_index, perturbation_id, drug, dose, dose_unit,
    dose_panel)], by = c("condition_index", "perturbation_id"),
  all.x = TRUE, sort = FALSE)
response_metadata <- merge(response_metadata,
  cell_order[, .(depmap_id, cell_line, cell_index)],
  by = c("depmap_id", "cell_line"), all.x = TRUE, sort = FALSE)
response_metadata[, `:=`(
  array_row = cell_index * EXPECTED_DISCOVERY_CONDITIONS + condition_index,
  n_plates = plate_count,
  shard_path = file.path("tahoe_pdex_selected", "response_shards",
    paste0(depmap_id, "__", cell_line, ".npy")),
  shard_row = condition_index)]
setorder(response_metadata, array_row)
response_metadata <- response_metadata[, .(array_row, cell_line, depmap_id,
  drug, dose, dose_unit, dose_panel, perturbation_id, available, n_plates,
  shard_path, shard_row)]
if (anyNA(response_metadata) || anyDuplicated(response_metadata$array_row))
  fail("Reconstructed response metadata are invalid.")

## Rebuild the exact-name, coverage-locked PRISM intersection.
prism_treatments <- fread(prism_paths[1L], showProgress = FALSE)
prism_cells <- fread(prism_paths[2L], showProgress = FALSE)
matrix_header <- names(fread(prism_paths[3L], nrows = 0L,
  check.names = FALSE, showProgress = FALSE))
matrix_rows <- fread(prism_paths[3L], select = 1L, showProgress = FALSE)
setnames(matrix_rows, 1L, "depmap_id")
matrix_rows[, depmap_id := clean_text(depmap_id)]
prism_treatments <- prism_treatments[
  perturbation_type == "experimental_treatment" & is.finite(dose) &
    dose > 0 & column_name %chin% matrix_header]
prism_treatments[, `:=`(prism_drug = clean_text(name),
  exact_name_key = name_key(name), prism_dose_uM = as.numeric(dose))]
conditions[, exact_name_key := name_key(drug)]
exact_keys <- intersect(unique(conditions$exact_name_key),
  unique(prism_treatments$exact_name_key))
exact_conditions <- conditions[exact_name_key %chin% exact_keys]
selected_cells <- cells35[depmap_id %chin% matrix_rows$depmap_id]
selected_cells <- merge(selected_cells, unique(prism_cells[, .(
  depmap_id = clean_text(depmap_id), prism_ccle_name = clean_text(ccle_name),
  prism_passed_str = as.character(passed_str_profiling))]),
  by = "depmap_id", all.x = TRUE, sort = FALSE)
setorder(selected_cells, depmap_id)
if (nrow(selected_cells) != 35L) fail("The raw Tahoe/PRISM line join changed.")

dose_matches <- exact_conditions[, {
  candidates <- sort(unique(prism_treatments[
    exact_name_key == .BY$exact_name_key, prism_dose_uM]))
  distance <- abs(log10(candidates) - log10(.BY$tahoe_dose_uM))
  minimum <- min(distance)
  selected <- min(candidates[distance <= minimum + 1e-12])
  .(prism_dose_uM = selected,
    absolute_log10_dose_difference = abs(log10(selected) -
      log10(.BY$tahoe_dose_uM)),
    dose_ratio = max(selected / .BY$tahoe_dose_uM,
      .BY$tahoe_dose_uM / selected))
}, by = .(exact_name_key, tahoe_drug = drug, tahoe_dose_uM = dose,
  dose_unit, dose_panel, perturbation_id)]
selected_columns <- merge(prism_treatments, dose_matches,
  by = c("exact_name_key", "prism_dose_uM"), all = FALSE,
  allow.cartesian = TRUE)
if (anyDuplicated(selected_columns$column_name))
  fail("A PRISM column maps to more than one Tahoe condition.")
selected_indices <- match(selected_columns$column_name, matrix_header)
prism_values <- fread(prism_paths[3L],
  select = unique(c(1L, selected_indices)), check.names = FALSE,
  showProgress = FALSE)
setnames(prism_values, 1L, "depmap_id")
prism_values[, depmap_id := clean_text(depmap_id)]
prism_values <- prism_values[depmap_id %chin% selected_cells$depmap_id]
long_function <- melt(prism_values, id.vars = "depmap_id",
  variable.name = "column_name", value.name = "prism_log2_fold_change",
  variable.factor = FALSE)
long_function <- merge(long_function,
  selected_columns[, .(column_name, exact_name_key,
    perturbation_id, prism_dose_uM, screen_id, broad_id,
    tahoe_drug, tahoe_dose_uM)],
  by = "column_name", all.x = TRUE)
function_values <- long_function[, .(
  n_prism_replicates = sum(is.finite(prism_log2_fold_change)),
  prism_log2_fold_change = if (all(!is.finite(prism_log2_fold_change)))
    NA_real_ else mean(prism_log2_fold_change[is.finite(
      prism_log2_fold_change)]),
  prism_sensitization = if (all(!is.finite(prism_log2_fold_change)))
    NA_real_ else -mean(prism_log2_fold_change[is.finite(
      prism_log2_fold_change)]),
  prism_screens = paste(sort(unique(screen_id)), collapse = ";"),
  prism_broad_ids = paste(sort(unique(broad_id)), collapse = ";")),
  by = .(depmap_id, exact_name_key, tahoe_drug,
    tahoe_dose_uM, perturbation_id, prism_dose_uM)]
minimum_complete_lines <- floor(0.75 * 35L)
drug_coverage <- function_values[, .(
  complete = .N == 3L && all(is.finite(prism_sensitization))),
  by = .(exact_name_key, tahoe_drug, depmap_id)][, .(
    n_lines_with_complete_three_dose_curve = sum(complete)),
  by = .(exact_name_key, tahoe_drug)]
selected_keys <- drug_coverage[
  n_lines_with_complete_three_dose_curve >= minimum_complete_lines,
  exact_name_key]
functional <- function_values[exact_name_key %chin% selected_keys]
if (uniqueN(functional$tahoe_drug) != EXPECTED_MATCHED_DRUGS ||
    uniqueN(functional$perturbation_id) != EXPECTED_MATCHED_CONDITIONS ||
    nrow(functional) != 35L * EXPECTED_MATCHED_CONDITIONS)
  fail("The exact-name PRISM intersection changed.")

## Fix the accepted corrected cohort and outer split.
selected_lines <- selected_cells[!depmap_id %chin% EXCLUDED_IDS]
response_metadata <- response_metadata[!depmap_id %chin% EXCLUDED_IDS]
functional <- functional[!depmap_id %chin% EXCLUDED_IDS]
fold_assignments <- selected_lines[, .(depmap_id, cell_line,
  tahoe_cell_name, organ)]
fold_assignments[, outer_fold := fifelse(depmap_id %chin% EXPECTED_HELD_IDS,
  OUTER_FOLD, 0L)]
setorder(fold_assignments, depmap_id)
if (nrow(selected_lines) != EXPECTED_LINES ||
    fold_assignments[outer_fold == OUTER_FOLD, uniqueN(depmap_id)] !=
      EXPECTED_HELD_LINES ||
    fold_assignments[outer_fold != OUTER_FOLD, uniqueN(depmap_id)] !=
      EXPECTED_TRAINING_LINES ||
    !identical(sort(fold_assignments[outer_fold == OUTER_FOLD, depmap_id]),
      sort(EXPECTED_HELD_IDS))) fail("The fixed 26/7 split changed.")
response_metadata <- merge(response_metadata,
  fold_assignments[, .(depmap_id, outer_fold)], by = "depmap_id",
  all.x = TRUE, sort = FALSE)
setorder(response_metadata, array_row)
response_metadata[, observation_id := paste(depmap_id, perturbation_id,
  sep = "__")]

## Prepare contemporaneous DMSO line means from the canonical float32 matrix.
control_meta <- fread(file.path(CANONICAL_CACHE, "dmso_control_metadata.tsv"))
control_genes <- fread(file.path(CANONICAL_CACHE, "dmso_control_genes.tsv"))
if (!identical(control_genes$feature, genes$feature))
  fail("DMSO and response gene orders differ.")
control <- read_npy_float32(file.path(CANONICAL_CACHE,
  "dmso_control_expression.npy"))
if (!identical(dim(control), c(35L * 14L, EXPECTED_GENES)))
  fail("The canonical DMSO matrix dimensions changed.")
line_ids <- sort(selected_lines$depmap_id)
tahoe_expression <- t(vapply(line_ids, function(id) {
  rows <- control_meta[depmap_id == id, array_row] + 1L
  if (length(rows) != 14L) fail("DMSO plate count changed for %s.", id)
  colMeans(control[rows, , drop = FALSE])
}, numeric(EXPECTED_GENES)))
rownames(tahoe_expression) <- line_ids
colnames(tahoe_expression) <- genes$feature
rm(control); gc(FALSE)

## Stream only the 33 selected rows from the canonical DepMap CSV.
expected_depmap_md5 <- "aadbbf68c587db63eba025b8bf7a792a"
if (!identical(unname(tools::md5sum(depmap_path)), expected_depmap_md5))
  fail("The canonical DepMap expression checksum changed.")
con <- file(depmap_path, open = "r"); on.exit(close(con), add = TRUE)
header <- strsplit(readLines(con, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1L]]
original_genes <- header[-1L]
depmap_raw <- matrix(NA_real_, nrow = EXPECTED_LINES,
  ncol = length(original_genes), dimnames = list(line_ids, original_genes))
wanted <- setNames(seq_along(line_ids), line_ids); matched <- logical(EXPECTED_LINES)
repeat {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) break
  comma <- regexpr(",", line, fixed = TRUE)[1L]
  id <- substr(line, 1L, comma - 1L); position <- unname(wanted[id])
  if (!is.na(position)) {
    values <- scan(text = substring(line, comma + 1L), sep = ",",
      quiet = TRUE, what = double())
    if (length(values) != ncol(depmap_raw) || any(!is.finite(values)))
      fail("Invalid DepMap row: %s", id)
    depmap_raw[position, ] <- values; matched[position] <- TRUE
  }
}
close(con); on.exit(NULL, add = FALSE)
if (!all(matched)) fail("DepMap lacks at least one selected line.")
parsed <- sub(" \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$", "", original_genes,
  perl = TRUE)
keep <- !is.na(parsed) & nzchar(parsed)
collapsed <- rowsum(t(depmap_raw[, keep, drop = FALSE]),
  group = parsed[keep], reorder = TRUE)
collapsed <- collapsed / as.numeric(table(parsed[keep])[rownames(collapsed)])
depmap_expression <- t(collapsed)
rm(depmap_raw, collapsed); gc(FALSE)

functional <- merge(functional,
  response_metadata[, .(depmap_id, perturbation_id, cell_line,
    molecular_observed = available, outer_fold)],
  by = c("depmap_id", "perturbation_id"), all.x = TRUE, sort = FALSE)
functional[, `:=`(function_observed = as.integer(is.finite(
  prism_sensitization)), chain_observed = as.integer(
  is.finite(prism_sensitization) & molecular_observed == 1L))]
setorder(functional, depmap_id, tahoe_drug, tahoe_dose_uM)

fwrite(fold_assignments, file.path(RESULT, "fold_assignments.tsv"), sep = "\t")
fwrite(selected_lines, file.path(RESULT, "cohort.tsv"), sep = "\t")
fwrite(drug_coverage, file.path(RESULT, "prism_drug_coverage.tsv"), sep = "\t")
fwrite(response_metadata, file.path(ARTIFACT, "response_metadata.tsv.gz"),
  sep = "\t", compress = "gzip")
fwrite(functional, file.path(ARTIFACT, "function_outcomes.tsv.gz"),
  sep = "\t", compress = "gzip")

prepared <- list(
  response = list(genes = genes$feature,
    metadata = as.data.frame(response_metadata),
    shard_root = TAHOE_PACKAGE),
  tahoe_baseline = list(expression = tahoe_expression,
    expression_scale = "released full-gene DMSO ref_mean; arithmetic 14-plate line mean"),
  depmap_baseline = list(expression = depmap_expression,
    expression_scale = "DepMap Public 22Q2 log2(TPM+1)"),
  fold_assignments = as.data.frame(fold_assignments),
  function_outcomes = as.data.frame(functional),
  parameters = list(seed = SEED, outer_fold = OUTER_FOLD,
    active_baseline_inner_folds = ACTIVE_BASELINE_INNER_FOLDS))
saveRDS(prepared, file.path(ARTIFACT, "prepared_inputs.rds"), compress = TRUE)
saveRDS(prepared$depmap_baseline,
  file.path(ARTIFACT, "depmap_baseline_expression.rds"), compress = TRUE)

## Source checksums. Tahoe package checksums come from the validated upload
## manifest; other canonical files are hashed directly.
source_manifest <- manifest_files[, .(relative_path = file.path(
  "cache/tahoe_prism_full_chain/tahoe_pdex_selected", path),
  file_size_bytes = size, checksum_type = "sha256_package_manifest",
  checksum = sha256, dimensions = fifelse(grepl("response_shards", path),
    "1137 x 62710 float32", fifelse(grepl("plate_counts", path),
      "1137", "see package manifest")))]
direct <- c(control_paths, prism_paths, depmap_path, tahoe_cell_path)
direct_manifest <- data.table(relative_path = substring(direct,
  nchar(WORK_ROOT) + 2L), file_size_bytes = file.info(direct)$size,
  checksum_type = "md5", checksum = unname(tools::md5sum(direct)),
  dimensions = "recorded by producing-stage contract")
source_manifest <- rbindlist(list(source_manifest, direct_manifest),
  use.names = TRUE)
fwrite(source_manifest, file.path(RESULT, "source_input_manifest.tsv"),
  sep = "\t")

census <- data.table(item = c("retained_cell_lines", "training_lines",
  "held_out_lines", "discovery_drugs", "discovery_conditions",
  "matched_prism_drugs", "matched_prism_conditions", "function_rows",
  "finite_function_rows", "response_genes", "tahoe_baseline_genes",
  "depmap_baseline_genes"), observed = c(nrow(selected_lines),
  sum(fold_assignments$outer_fold != OUTER_FOLD),
  sum(fold_assignments$outer_fold == OUTER_FOLD), uniqueN(conditions$drug),
  uniqueN(conditions$perturbation_id), uniqueN(functional$tahoe_drug),
  uniqueN(functional$perturbation_id), nrow(functional),
  sum(functional$function_observed), length(genes$feature),
  ncol(tahoe_expression), ncol(depmap_expression)))
fwrite(census, file.path(RESULT, "input_census.tsv"), sep = "\t")
write_artifact_manifest(ARTIFACT, "scripts/01_prepare_inputs.R")
writeLines(sprintf("complete\nruntime_seconds=%.3f\n",
  proc.time()[["elapsed"]] - started), file.path(RESULT,
  "stage01_complete.flag"))
plot_stage01()
message(sprintf("Stage 01 complete: %d lines; %d response conditions; %d PRISM conditions.",
  nrow(selected_lines), uniqueN(response_metadata$perturbation_id),
  uniqueN(functional$perturbation_id)))
