# Prepare measured RNA, independent growth and baseline protein data from Niepel et al. 2017.
# This is a data adapter for analysis 59: no model fitting, outcome-based filtering or program learning.
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(rhdf5)
  library(qs)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
source_root <- "/work/v20_pilot_design_2026-10-09/source"
output_root <- "/work/results/program_pilot/59_mechanistic_evidence_datasets"
data_root <- file.path(output_root, "data")
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "source"), recursive = TRUE, showWarnings = FALSE)
study <- "Niepel_2017_LJP"
study_url <- "https://www.nature.com/articles/s41467-017-01383-w"
cells <- c("BT20", "HS578T", "MCF10A", "MCF7", "MDAMB231", "SKBR3")
gctx <- file.path(source_root, "GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx")

# Exact source aliases, and rounding-equivalent concentrations only: no interpolation.
canonical_drug <- function(x) {
  x <- tolower(trimws(x))
  aliases <- c("bez235" = "dactolisib", "nvp-bez235" = "dactolisib",
               "azd6244" = "selumetinib", "byl-719" = "alpelisib",
               "byl719" = "alpelisib", "brd-k54997624" = "alpelisib")
  keep <- x %in% names(aliases)
  x[keep] <- aliases[x[keep]]
  x
}
nominal_dose <- function(x) {
  grid <- c(0.04, 0.12, 0.37, 1.11, 3.33, 10)
  vapply(x, function(v) {
    if (is.na(v) || v < 0) return(NA_real_)
    if (v == 0) return(0)
    i <- which.min(abs(log(grid / v)))
    if (abs(log(grid[i] / v)) <= 0.035) grid[i] else NA_real_
  }, numeric(1))
}
write_data <- function(x, name) fwrite(x, file.path(data_root, paste0("matched_", name)), sep = "\t", na = "NA")

# Primary molecular samples: preserve every measured sample in the six study cells, both times, and vehicle.
rna <- fread(file.path(source_root, "LJP005_006_sample_metadata.tsv"))[cell_id %in% cells]
stopifnot(!anyDuplicated(rna$inst_id), all(rna$pert_time %in% c(3, 24)))
rna[, `:=`(study_id = study, sample_id = trimws(inst_id), cell_line = cell_id,
           drug = canonical_drug(pert_iname), dose_uM = fifelse(pert_type == "ctl_vehicle", 0, pert_dose),
           time_h = pert_time, assay = "L1000_measured_landmarks_Level2",
           replicate_id = sub(".*_X([^_]+)_.*", "\\1", det_plate),
           replicate_type = "L1000_assay_replicate; paper_methods_call_technical",
           control_group_id = det_plate,
           evaluation_status = "retrospective_resource_previously_inspected; no_new_holdout_claim")]
rna[, condition_id := paste(study_id, cell_line, pert_id, dose_uM, sep = "|")]
rna[, source_path := gctx]
column_ids <- trimws(h5read(gctx, "0/META/COL/id"))
rna[, gctx_column_index := match(sample_id, column_ids)]
stopifnot(!anyNA(rna$gctx_column_index), !anyDuplicated(rna$gctx_column_index))
setorder(rna, gctx_column_index)
controls <- rna[, .(n_samples = .N, n_vehicle = sum(pert_type == "ctl_vehicle")), by = det_plate]
controls[, plate_vehicle_normalization_eligible := n_vehicle >= 5]
rna <- merge(rna, controls[, .(det_plate, n_vehicle, plate_vehicle_normalization_eligible)],
             by = "det_plate", all.x = TRUE, sort = FALSE)
setorder(rna, gctx_column_index)

# Numeric phenotype data retain individual biological repeats, not a fitted curve or sensitivity label.
single <- fread(file.path(source_root, "41467_2017_1383_MOESM5_ESM.txt"))
stopifnot(all(single$CellLine %in% cells), all(is.finite(single$GRvalue)), all(single$RelCellCnt > 0))
drug_names <- unique(rna[pert_type == "trt_cp", .(BRDid = pert_id, drug)])
stopifnot(!anyDuplicated(drug_names$BRDid))
single <- merge(single, drug_names, by = "BRDid", all.x = TRUE, sort = FALSE)
single[is.na(drug), drug := canonical_drug(DrugName)]
single[, `:=`(study_id = study, sample_id = sprintf("LJP_growth_single_%05d", seq_len(.N)),
              cell_line = CellLine, dose_uM = Conc, time_h = 72,
              biological_replicate = Biological_Replicate,
              assay = "image_count_derived_growth", source_table = "Supplementary_Data_2",
              source_path = file.path(source_root, "41467_2017_1383_MOESM5_ESM.txt"))]
single[, condition_id := paste(study_id, cell_line, BRDid, dose_uM, sep = "|")]
write_data(single, "single_agent_growth.tsv.gz")
function_measurements <- melt(single[, .(study_id, sample_id, condition_id, cell_line, drug, dose_uM,
                                          time_h, biological_replicate, source_path,
                                          GR = GRvalue, relative_cell_count = RelCellCnt)],
                              id.vars = c("study_id", "sample_id", "condition_id", "cell_line", "drug", "dose_uM",
                                          "time_h", "biological_replicate", "source_path"),
                              variable.name = "feature", value.name = "value")
function_measurements[, unit := fifelse(feature == "GR", "GR: 1=control growth; 0=cytostasis; negative=net cell loss",
                                        "72h treated cells / 72h vehicle cells")]
function_measurements[, assay := "single_agent_function"]

# Preserve day0/72h counts from the deposited screen as a separate source representation,
# not additional independent replicates of the supplementary table above.
raw_counts <- as.data.table(read_excel(file.path(source_root, "Screen20237_LJP_processedData1.xls"), sheet = "Data"))
counts <- raw_counts[, .(study_id = study, sample_id = sprintf("LJP_count_source_%05d", seq_len(.N)),
                         cell_line = CellLine, drug = canonical_drug(DrugName), BRDid,
                         dose_uM = `Small Mol Concentration`, time_h = Timepoint,
                         biological_replicate = Replicate, library_plate = `LJP Library Plate`,
                         well = `Assay Well`, assay_barcode = Barcode,
                         n0 = `Total Cell Count Before Treatment`, n72 = `Total Cell Count After Treatment`,
                         ncontrol72 = `Total Control Cell Count`,
                         reported_relative_count = `Relative Cell Count`,
                         source_normalized_net_growth = `Relative Growth Rate After Treatment`)]
counts[, `:=`(recomputed_GR = 2^(log2(n72 / n0) / log2(ncontrol72 / n0)) - 1,
              recomputed_normalized_net_growth = (n72 - n0) / (ncontrol72 - n0),
              source_path = file.path(source_root, "Screen20237_LJP_processedData1.xls"),
              overlap_note = "underlying_source_representation; do_not_stack_with_supplement_as_extra_replicates")]
stopifnot(all(is.finite(counts$recomputed_GR)))
stopifnot(max(abs(counts$source_normalized_net_growth - counts$recomputed_normalized_net_growth)) < 1e-4)
write_data(counts, "single_agent_source_counts.tsv.gz")

# Combination screen includes all raw wells, separate constituent effects, vehicle and initial counts.
combo_raw <- as.data.table(read_excel(file.path(source_root, "Screen20259_LJP_drugcombo_all.xls"), sheet = "Data"))
combo <- combo_raw[, .(study_id = study, sample_id = sprintf("LJP_combo_%05d", seq_len(.N)),
                       cell_line = CellLine, drug_a = canonical_drug(DrugName), drug_b = canonical_drug(DrugName2),
                       dose_a_uM = as.numeric(`Small Mol 1 Concentration`),
                       dose_b_uM = as.numeric(`Small Mol 2 Concentration`),
                       biological_replicate = `Biological Replicate`, technical_replicate = `Technical Replicate`,
                       n0 = `Total Cell Count Before Treatment`, n72 = `Total Cell Count After Treatment`,
                       ncontrol72 = `Total Control Cell Count`, GR = `Normalized Growth Rate Inhibition Value`)]
combo[, `:=`(time_h = 72, nominal_dose_a_uM = nominal_dose(dose_a_uM),
              nominal_dose_b_uM = nominal_dose(dose_b_uM),
              source_path = file.path(source_root, "Screen20259_LJP_drugcombo_all.xls"))]
combo[, `:=`(recomputed_GR = 2^(log2(n72 / n0) / log2(ncontrol72 / n0)) - 1,
              relative_cell_count = n72 / ncontrol72, fold_vs_day0 = n72 / n0,
              net_cell_loss = n72 < n0,
              arm = fcase(drug_a == "-" & drug_b == "-", "vehicle", drug_b == "-", "single", default = "combination"))]
stopifnot(max(abs(combo$GR - combo$recomputed_GR)) < 1e-4)
combo[, condition_id := paste(study_id, cell_line, drug_a, dose_a_uM, drug_b, dose_b_uM, sep = "|")]
write_data(combo, "combination_counts.tsv.gz")
combo_measurements <- melt(combo[, .(study_id, sample_id, condition_id, cell_line, drug_a, drug_b,
                                      dose_a_uM, dose_b_uM, time_h, biological_replicate, technical_replicate,
                                      source_path, n0, n72, ncontrol72, GR, relative_cell_count, fold_vs_day0)],
                           id.vars = c("study_id", "sample_id", "condition_id", "cell_line", "drug_a", "drug_b",
                                       "dose_a_uM", "dose_b_uM", "time_h", "biological_replicate",
                                       "technical_replicate", "source_path"),
                           variable.name = "feature", value.name = "value")
combo_measurements[, `:=`(assay = "combination_screen_function",
                          unit = fcase(feature %in% c("n0", "n72", "ncontrol72"), "imaged_cell_count",
                                       feature == "GR", "GR: 1=control growth; 0=cytostasis; negative=net cell loss",
                                       default = "cell_count_ratio"))]

# Independent baseline proteins and phosphopeptides are actual measured abundances, not RNA inference.
protein_raw <- fread(file.path(source_root, "41467_2017_1383_MOESM8_ESM.txt"))
phospho_raw <- fread(file.path(source_root, "41467_2017_1383_MOESM9_ESM.txt"))
protein <- melt(protein_raw, id.vars = c("Gene name", "Replicate"), variable.name = "cell_line", value.name = "value")
setnames(protein, c("Gene name", "Replicate"), c("feature", "biological_replicate"))
protein[, `:=`(study_id = study, assay = "baseline_protein_mass_spectrometry", time_h = 0,
                unit = "source_reported_normalized_protein_abundance",
                sample_id = paste("LJP_baseline_protein", cell_line, biological_replicate, sep = "|"),
                source_path = file.path(source_root, "41467_2017_1383_MOESM8_ESM.txt"))]
phospho <- melt(phospho_raw, id.vars = c("Gene name", "Sequence", "Replicate"),
                variable.name = "cell_line", value.name = "value")
setnames(phospho, c("Gene name", "Replicate", "Sequence"), c("feature", "biological_replicate", "peptide_sequence"))
phospho[, `:=`(study_id = study, assay = "baseline_phosphopeptide_mass_spectrometry", time_h = 0,
                unit = "source_reported_normalized_phosphopeptide_abundance",
                sample_id = paste("LJP_baseline_phosphopeptide", cell_line, biological_replicate, sep = "|"),
                source_path = file.path(source_root, "41467_2017_1383_MOESM9_ESM.txt"))]
measurements <- rbindlist(list(function_measurements, combo_measurements, protein, phospho), fill = TRUE)
write_data(measurements, "measurements.tsv.gz")

# Explicit molecular/functional join: cell + compound identity + nominal dose; times remain separate.
rna_coverage <- rna[pert_type == "trt_cp", .(rna_samples = .N, rna_plates = uniqueN(det_plate),
                                             rna_sample_ids = paste(sample_id, collapse = ";")),
                    by = .(condition_id, cell_line, drug, dose_uM, time_h)]
growth_coverage <- single[, .(growth_samples = .N, growth_biological_replicates = uniqueN(biological_replicate),
                              growth_sample_ids = paste(sample_id, collapse = ";")), by = condition_id]
matched <- merge(rna_coverage, growth_coverage, by = "condition_id", all.x = TRUE, sort = FALSE)
matched[, `:=`(rna_time_h = time_h, function_time_h = 72,
                matching_rule = "same_cell_same_BRD_compound_same_nominal_dose; separate_assay_wells_not_same_cells",
                valid_function_match = !is.na(growth_samples),
                independent_evaluation_unit = cell_line)]
write_data(matched, "rna_growth_conditions.tsv.gz")
coverage <- matched[, .(RNA_conditions = .N, matched_conditions = sum(valid_function_match),
                         RNA_samples = sum(rna_samples), matched_drugs = uniqueN(drug[valid_function_match])),
                     by = .(cell_line, rna_time_h)]
write_data(coverage, "coverage.tsv")

# Compound combinations map to measured constituent RNA, not invented joint RNA.
rna_compound <- rna[pert_type == "trt_cp", .(rna_samples = .N), by = .(cell_line, drug, dose_uM, time_h)]
grids <- unique(combo[arm == "combination", .(cell_line, drug_a, drug_b, dose_a_uM, dose_b_uM,
                                              nominal_dose_a_uM, nominal_dose_b_uM)])
grids <- rbindlist(list(copy(grids)[, time_h := 3], copy(grids)[, time_h := 24]))
grids <- merge(grids, rna_compound, by.x = c("cell_line", "drug_a", "nominal_dose_a_uM", "time_h"),
               by.y = c("cell_line", "drug", "dose_uM", "time_h"), all.x = TRUE)
setnames(grids, "rna_samples", "rna_a_samples")
grids <- merge(grids, rna_compound, by.x = c("cell_line", "drug_b", "nominal_dose_b_uM", "time_h"),
               by.y = c("cell_line", "drug", "dose_uM", "time_h"), all.x = TRUE)
setnames(grids, "rna_samples", "rna_b_samples")
grids[, `:=`(joint_RNA_measured = FALSE, constituent_RNA_supported = !is.na(rna_a_samples) & !is.na(rna_b_samples),
              whole_partner_exclusion_unit = drug_b,
              evaluation_status = "retrospective: all existing outcomes previously available; no blind test claim")]
write_data(grids, "combination_rna_links.tsv")

# Store the actual measured 978-gene intensities compactly. No inferred genes are included.
# Raw level2 is intentionally preserved: plate-aware log2 quantile normalization and vehicle subtraction
# must be made explicit in the downstream biological comparison, not confused with source measurements.
gene_ids <- as.character(h5read(gctx, "0/META/ROW/id"))
gene_symbols <- unname(mapIds(org.Hs.eg.db, keys = gene_ids, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))
stopifnot(length(gene_ids) == 978L, !anyNA(gene_symbols), !anyDuplicated(gene_symbols))
rna_matrix <- matrix(NA_real_, nrow = length(gene_ids), ncol = nrow(rna),
                     dimnames = list(gene_symbols, rna$sample_id))
# Read contiguous chunks to avoid a very large HDF5 index selection and avoid duplicating the full source.
for (start in seq(1L, nrow(rna), by = 256L)) {
  keep <- start:min(start + 255L, nrow(rna))
  source_indices <- rna$gctx_column_index[keep]
  block <- h5read(gctx, "0/DATA/0/matrix", start = c(1L, min(source_indices)),
                  count = c(978L, max(source_indices) - min(source_indices) + 1L))
  rna_matrix[, keep] <- block[, source_indices - min(source_indices) + 1L, drop = FALSE]
}
stopifnot(all(is.finite(rna_matrix)), all(rna_matrix >= 0), identical(colnames(rna_matrix), rna$sample_id))
qsave(rna_matrix, file.path(data_root, "matched_measured_RNA_Level2.qs"), preset = "high")
gene_annotation <- data.table(gene_symbol = gene_symbols, source_entrez_id = gene_ids,
                              measurement_status = "directly_measured_L1000_landmark; not_inferred",
                              annotation_package = paste0("org.Hs.eg.db_", packageVersion("org.Hs.eg.db")))
write_data(gene_annotation, "gene_annotation.tsv")
write_data(rna, "samples.tsv.gz")
write_data(controls, "RNA_plate_controls.tsv")

# Dataset-level contrasts are specified by experimental design before any model fitting.
contrasts <- data.table(
  study_id = study,
  contrast_id = c("same_drug_dose_different_cells", "within_cell_dose_and_time", "partner_substitution", "adaptive_signaling_timecourse"),
  question = c("Does RNA distinguish different growth outcomes at identical exposure?",
               "Does reproducible RNA onset precede arrest or net cell loss across the whole dose trajectory?",
               "Do separate drug responses identify effective substitute combination partners beyond constituent strength?",
               "Does alpelisib-induced adaptation rebound at pAKT/pERK and is rebound prevented by partners?"),
  measured_comparator = c("978 measured RNA genes vs independent 72h GR and cell counts in 6 cells",
                          "3h/24h RNA and plate vehicles; 6-dose 72h growth curves in biological triplicates",
                          "BT20 constituent RNA and raw single/combination/day0/vehicle cell counts",
                          "published quantitative IF pERK/pAKT distributions for vehicle,A,B,AB at3/6/12/24h"),
  falsification = c("No transferable phenotype separation despite RNA differences larger than assay replicate variation",
                   "Apparent response is unreproducible or dose-only and does not separate cell-specific growth",
                   "Molecular similarity fails to recover substitute-partner benefit beyond single-agent effects",
                   "No rebound or no selective prevention; cannot substitute RNA gene-set changes for this measurement"),
  holdout_unit = c("entire_cell_all_drugs_doses_times", "entire_drug_family_or_cell_not_dose_points",
                   "entire_partner_all_doses_and_pairs", "independent_experimental_repeat"),
  availability = c("numeric_RNA_and_function_prepared", "numeric_RNA_function_controls_prepared",
                   "numeric_constituent_RNA_and_combo_function_prepared; no_joint_RNA",
                   "figure_only_in_local_supplement; underlying_numeric_IF_not_verified"))
write_data(contrasts, "contrasts.tsv")
studies <- data.table(study_id = study, title = "Common and cell-type specific responses to anti-cancer drugs",
                      source_url = study_url, cells = paste(cells, collapse = ";"), directly_measured_genes = 978L,
                      RNA_samples = nrow(rna), RNA_times_h = "3;24", function_time_h = 72,
                      single_function_rows = nrow(single), unique_single_conditions = uniqueN(single$condition_id),
                      combination_wells = nrow(combo[arm == "combination"]),
                      combination_pairs = uniqueN(combo[arm == "combination", .(drug_a, drug_b)]),
                      matched_RNA_function_conditions = nrow(matched[valid_function_match == TRUE]),
                      baseline_protein_measurements = nrow(protein), baseline_phosphopeptide_measurements = nrow(phospho),
                      joint_combination_RNA = FALSE,
                      license_note = "article_CC_BY4.0; source_HMS_LINCS_terms_and_GEO_repository_terms_retained",
                      caveat = "RNA landmark assay replicates are not declared independent biological experiments; studies already inspected" )
write_data(studies, "studies.tsv")

# Reused downloads are referenced in place with checksums; no large source-file copying.
inputs <- c(gctx, file.path(source_root, c("LJP005_006_sample_metadata.tsv", "41467_2017_1383_MOESM5_ESM.txt",
                "Screen20237_LJP_processedData1.xls", "Screen20259_LJP_drugcombo_all.xls",
                "41467_2017_1383_MOESM8_ESM.txt", "41467_2017_1383_MOESM9_ESM.txt", "supp_info.pdf")))
manifest <- data.table(study_id = study, local_path = inputs, bytes = file.info(inputs)$size,
                       retrieval_status = "existing_download_verified_readable_reused_in_place",
                       source_url = c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138",
                                      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138",
                                      study_url,
                                      "https://lincs.hms.harvard.edu/_static/db/prod-20200624/Screen20237_LJP_processedData1.xls",
                                      "https://lincs.hms.harvard.edu/_static/db/prod-20200624/Screen20259_LJP_drugcombo_all.xls",
                                      study_url, study_url, study_url))
manifest[, md5 := unname(tools::md5sum(local_path))]
write_data(manifest, "download_manifest.tsv")
capture.output(sessionInfo(), file = file.path(data_root, "matched_session_info.txt"))
writeLines(c("Actual numeric observations prepared; no model, gene selection or outcome-based reliability filtering fitted.",
             "RNA matrix rows are gene symbols; columns exactly match sample_id in matched_samples.tsv.gz.",
             "RNA values are raw Level2 deconvolved fluorescence, not log fold change or inferred expression.",
             "Growth uses independent same-condition assays; do not equate RNA replicate IDs with growth biological replicates.",
             "matched_single_agent_source_counts duplicates the underlying assay represented by supplementary growth; do not stack.",
             "Legacy count-table Relative Growth Rate After Treatment is (N72-N0)/(Ncontrol72-N0), NOT the published GR formula.",
             "Function joins use cell, BRD compound and dose. Combo nominal dose matching only tolerates rounding within3.5%.",
             "Biological repeats distinguish growth from day0 (cytostasis vs net cell loss), unlike vehicle-relative viability alone.",
             "pERK/pAKT intervention timecourse is present as figure, not yet numeric data; baseline phosphoproteomics is numeric.",
             "All study outcomes previously inspected. Defined exclusion units support retrospective controlled evaluations only."),
           file.path(data_root, "matched_README.txt"))
print(studies)
print(coverage)
cat("Recomputed combination GR max error:", max(abs(combo$GR - combo$recomputed_GR)), "\n")
cat("Source normalized net-growth max error:",
    max(abs(counts$source_normalized_net_growth - counts$recomputed_normalized_net_growth)), "\n")
