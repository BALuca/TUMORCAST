#!/usr/bin/env Rscript
# Reconstruct the 72 h NRG1/EGF calibration measurements from original Prism XML.
# Preserve the exact replicate pairing and log-dose interpolation of analysis 123.

suppressPackageStartupMessages({
  library(xml2)
  library(dplyr)
})

raw_dir <- "data/raw/ligand_assays"
output_dir <- "results/10_prepare_ligand_assays"
outcome_path <- "data/processed/her2_ligand_outcomes.tsv"
target_doses <- c(0.5, 5)
comparison_tolerance <- 1e-10
source_files <- c("elife-32271-fig1-data1-v1.pzfx", "elife-32271-fig1-figsupp1-data1-v1.pzfx")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(outcome_path), recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
paths <- file.path(raw_dir, source_files)
if (!all(file.exists(paths))) stop("Missing original eLife Prism source files")
documents <- lapply(paths, read_xml)

# The fixed table choices retain parental 2D assays. They exclude AU565,
# transfected cells, 3D assays and additional inhibitors. MCF7 is read once.
selected <- data.frame(
  source_number = c(1, 1, 2, 2, 2, 2, 2, 2),
  table_id = c("Table11", "Table11", "Table13", "Table17", "Table27", "Table29", "Table19", "Table21"),
  CellLine = c("SKBR3", "MCF7", "BT474", "ZR751", "HCC1569", "HCC1419", "BT474", "SKBR3"),
  gene = c(rep("NRG1", 6), "EGF", "EGF"),
  drug_column = c("SKBR3 LAP", "MCF7 LAP", "LAP", "LAP", "HCC1569 Lapatinib", "Lapatinib", "LAP", "LAP"),
  combination_column = c(
    "SKBR3 LAP NRG", "MCF7 LAP NRG", "LAP NRG", "LAP NRG",
    "HCC1569 LAP NRG", "Lapatinib + Heregulin", "LAP EGF", "LAP EGF"
  )
)
annotation <- read.csv("data/raw/depmap/Model.csv", check.names = FALSE)
selected$depmap_id <- annotation$ModelID[match(selected$CellLine, annotation$StrippedCellLineName)]
stopifnot(!anyNA(selected$depmap_id), !any(selected$depmap_id == "ACH-000248"))

# Return a dose-by-replicate matrix for one named Prism column. The replicate
# columns must remain in their original order for paired drug/ligand contrasts.
read_replicates <- function(table_node, column_title) {
  column <- xml_find_first(table_node, paste0('./YColumn[Title="', column_title, '"]'))
  if (inherits(column, "xml_missing")) stop("Missing Prism column: ", column_title)
  values <- lapply(xml_find_all(column, "./Subcolumn"), function(replicate) {
    as.numeric(xml_text(xml_find_all(replicate, "./d")))
  })
  stopifnot(length(values) >= 3L, all(lengths(values) == 8L))
  result <- do.call(cbind, values)
  stopifnot(all(is.finite(result)), all(result > 0))
  result
}

outcome_parts <- list()
for (i in seq_len(nrow(selected))) {
  cat("Reading", selected$CellLine[i], selected$gene[i], selected$table_id[i], "\n")
  table_node <- xml_find_first(
    documents[[selected$source_number[i]]],
    paste0('//Table[@ID="', selected$table_id[i], '"]')
  )
  if (inherits(table_node, "xml_missing")) stop("Missing Prism table: ", selected$table_id[i])
  doses <- as.numeric(xml_text(xml_find_all(table_node, "./XColumn/Subcolumn/d")))
  # The plotted 0.0001 abscissa is a graph-only placeholder for the zero-dose controls.
  # Some X columns have an extra plotting abscissa with no paired Y measurements.
  # The historical analysis uses the first eight rows, which have replicate counts.
  stopifnot(length(doses) >= 8L, doses[1] == 0.0001)
  doses <- doses[1:8]
  stopifnot(all(diff(doses) > 0))
  drug <- read_replicates(table_node, selected$drug_column[i])
  combination <- read_replicates(table_node, selected$combination_column[i])
  stopifnot(identical(dim(drug), dim(combination)))

  for (dose in target_doses) {
    stopifnot(dose >= min(doses[-1]), dose <= max(doses[-1]))
    # Interpolate log2 abundance against log dose within each original replicate.
    # This matches measured curves to the native Tahoe 0.5/5 uM doses without extrapolation.
    drug_log2 <- vapply(seq_len(ncol(drug)), function(replicate) {
      approx(log(doses[-1]), log2(drug[-1, replicate]), xout = log(dose))$y
    }, numeric(1))
    combination_log2 <- vapply(seq_len(ncol(combination)), function(replicate) {
      approx(log(doses[-1]), log2(combination[-1, replicate]), xout = log(dose))$y
    }, numeric(1))
    ligand_log2 <- log2(combination[1, ] / drug[1, ])
    departure_log2 <- combination_log2 - drug_log2 - ligand_log2

    outcome_parts[[paste(i, dose)]] <- data.frame(
      CellLine = selected$CellLine[i],
      depmap_id = selected$depmap_id[i],
      gene = selected$gene[i],
      drug = "Lapatinib ditosylate",
      dose = dose,
      observed_departure_log2 = mean(departure_log2),
      replicate_departure_SD = sd(departure_log2),
      observed_ligand_log2 = mean(ligand_log2),
      observed_drug_arm_increment_log2 = mean(combination_log2 - drug_log2),
      observed_drug_pct_vehicle = mean(100 * 2^drug_log2 / drug[1, ]),
      observed_combination_pct_vehicle = mean(100 * 2^combination_log2 / drug[1, ]),
      bio_replicates = ncol(drug),
      study = "elife32271",
      ligand_nM = 10,
      time_h = 72,
      source_table = selected$table_id[i],
      source_file = paths[selected$source_number[i]],
      outcome_match = "within-replicate log-abundance/log-dose interpolation; nativeTahoe0.5/5uM"
    )
  }
}
outcomes <- bind_rows(outcome_parts)
stopifnot(
  nrow(outcomes) == 16L,
  sum(outcomes$gene == "NRG1") == 12L,
  sum(outcomes$gene == "EGF") == 4L,
  n_distinct(outcomes$depmap_id) == 6L,
  !anyDuplicated(outcomes[, c("depmap_id", "gene", "drug", "dose")]),
  all(outcomes$time_h == 72)
)

# Check all numerical measurements against the accepted original preparation.
# The source_file field changes only because the original files now live here.
reference <- read.delim("reference/results/123_her2_ligand_independent_outcomes.tsv")
comparison <- inner_join(
  outcomes,
  reference,
  by = c("depmap_id", "gene", "drug", "dose"),
  suffix = c("_rebuilt", "_reference"),
  relationship = "one-to-one"
)
stopifnot(nrow(comparison) == nrow(reference), nrow(comparison) == nrow(outcomes))
quantities <- names(outcomes)[vapply(outcomes, is.numeric, logical(1))]
quantities <- setdiff(quantities, "dose")
differences <- vapply(quantities, function(quantity) {
  max(abs(comparison[[paste0(quantity, "_rebuilt")]] - comparison[[paste0(quantity, "_reference")]]))
}, numeric(1))
stopifnot(all(is.finite(differences)), all(differences < comparison_tolerance))

write.table(outcomes, outcome_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(
  data.frame(quantity = quantities, maximum_absolute_difference = differences),
  file.path(output_dir, "comparison_with_original.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("Reconstructed 12 NRG1 and 4 EGF observations across six recipients at 72 h.\n")
cat("All numerical measurements agree with the original preparation.\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
sink()
