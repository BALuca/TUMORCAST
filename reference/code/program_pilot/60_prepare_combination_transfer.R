#!/usr/bin/env Rscript
# Prepare observed Jaaks combination wells and measured LINCS constituent RNA.
# No fitting or outcome-based selection. All MEK+PI3K family pairs are reserved.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(rhdf5)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(limma)
  library(digest)
})
options(datatable.print.nrows = 12)
setDTthreads(4)
root <- "/work/results/program_pilot/60_defensible_response_composition"
out <- file.path(root, "data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source_root <- "/work/jaaks_2022/raw"
lincs_root <- "/work/v20_pilot_design_2026-10-09/source"
cells <- c("BT20", "MCF7", "MDAMB231", "HS578T")
key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
cell_key <- function(x) toupper(key(x))
drug_key <- function(x) {
  ans <- key(x)
  ans[ans == "gdc0941"] <- "pictilisib"
  ans
}
write_table <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
source_paths <- file.path(source_root, paste0(c("Original", "Validation"), " screen_All tissues_fitted data.csv"))
fit <- rbindlist(lapply(source_paths, function(f) {
  x <- fread(f, select = c("BARCODE", "COMBI_ID", "CELL_LINE_NAME", "ANCHOR_ID", "ANCHOR_NAME",
                           "ANCHOR_TARGET", "ANCHOR_CONC", "LIBRARY_ID", "LIBRARY_NAME", "LIBRARY_TARGET", "LIBRARY_CONC"))
  x[, source_release := sub(" .*", "", basename(f))]
  x
}))
fit[, `:=`(cell = cell_key(CELL_LINE_NAME), drug_A = drug_key(ANCHOR_NAME), drug_B = drug_key(LIBRARY_NAME))]
fit[, `:=`(family_A_MEK = grepl("MEK", ANCHOR_TARGET, fixed = TRUE),
           family_B_MEK = grepl("MEK", LIBRARY_TARGET, fixed = TRUE),
           family_A_PI3K = grepl("PI3K", ANCHOR_TARGET, fixed = TRUE),
           family_B_PI3K = grepl("PI3K", LIBRARY_TARGET, fixed = TRUE),
           multi_drug_backbone = grepl("|", ANCHOR_NAME, fixed = TRUE) | grepl("|", LIBRARY_NAME, fixed = TRUE))]
fit[, heldout := (family_A_MEK & family_B_PI3K) | (family_B_MEK & family_A_PI3K)]
fit[, pair := fifelse(drug_A < drug_B, paste(drug_A, drug_B, sep = " + "), paste(drug_B, drug_A, sep = " + "))]
# This exclusion manifest covers every cell and both releases, before four-cell restriction.
family_manifest <- unique(fit[, .(source_release, cell, pair, drug_A, drug_B, ANCHOR_TARGET, LIBRARY_TARGET,
                                  heldout, multi_drug_backbone)])
write_table(family_manifest, "transfer_family_exclusion_manifest.tsv.gz")
family_drugs <- unique(rbind(fit[, .(drug = drug_A, name = ANCHOR_NAME, target = ANCHOR_TARGET,
                                     MEK = family_A_MEK, PI3K = family_A_PI3K)],
                            fit[, .(drug = drug_B, name = LIBRARY_NAME, target = LIBRARY_TARGET,
                                     MEK = family_B_MEK, PI3K = family_B_PI3K)]))
write_table(family_drugs, "transfer_drug_family_registry.tsv")
fit <- unique(fit[cell %chin% cells & !multi_drug_backbone])
fit[, `:=`(ANCHOR_CONC = as.numeric(ANCHOR_CONC), LIBRARY_CONC = as.numeric(LIBRARY_CONC))]
stopifnot(all(is.finite(fit$ANCHOR_CONC)), all(is.finite(fit$LIBRARY_CONC)))
stopifnot(!anyDuplicated(fit[, .(source_release, BARCODE, COMBI_ID, ANCHOR_CONC)]))

# Stream a ZIP through the ordinary unzip file utility, filtering in R before CSV parsing.
# A complete extraction would consume 1.2GB of scarce disk. Cache only the four study cells.
raw_list <- list()
for (release in c("Original", "Validation")) {
  cache <- file.path(out, paste0("transfer_", tolower(release), "_source_wells.tsv.gz"))
  zipfile <- file.path(source_root, paste0(release, "_screen_All_tissues_raw_data.csv.zip"))
  if (!file.exists(cache)) {
    connection <- pipe(paste("unzip -p", shQuote(zipfile)), open = "r")
    header <- readLines(connection, n = 1L)
    chunks <- list()
    nread <- 0L
    repeat {
      lines <- readLines(connection, n = 100000L)
      if (!length(lines)) break
      nread <- nread + length(lines)
      selected <- grep(",MCF7,|,BT-20,|,MDA-MB-231,|,Hs-578-T,", lines, value = TRUE)
      if (length(selected)) chunks[[length(chunks) + 1L]] <- fread(text = paste(c(header, selected), collapse = "\n"))
      if (nread %% 1000000L == 0L) message(release, ": scanned ", nread, " source records")
    }
    close(connection)
    raw <- rbindlist(chunks)
    raw[, source_release := release]
    write_table(raw, basename(cache))
  } else raw <- fread(cache)
  stopifnot(all(cell_key(raw$CELL_LINE_NAME) %chin% cells))
  raw_list[[release]] <- raw
}
raw <- rbindlist(raw_list)
raw[, `:=`(cell = cell_key(CELL_LINE_NAME), drug = drug_key(DRUG_NAME),
            plate_id = paste(source_release, BARCODE, sep = "|"),
            assay_day = substr(DATE_CREATED, 1, 10))]
stopifnot(all(raw$ASSAY == "Glo"), all(raw$DURATION == 4))
# Two drug records share one physical combination well; only one intensity is measured.
well <- unique(raw[, .(source_release, BARCODE, plate_id, cell, assay_day, SCAN_DATE, POSITION, INTENSITY)])
stopifnot(!anyDuplicated(well[, .(plate_id, POSITION)]))
controls <- raw[TAG %chin% c("B", "NC-1"),
                .(mean = mean(INTENSITY), sd = sd(INTENSITY), n = .N), by = .(plate_id, TAG)]
controls <- dcast(controls, plate_id ~ TAG, value.var = c("mean", "sd", "n"))
setnames(controls, c("mean_B", "mean_NC-1", "sd_B", "sd_NC-1", "n_B", "n_NC-1"),
         c("blank_mean", "vehicle_mean", "blank_sd", "vehicle_sd", "n_blank", "n_vehicle"))
controls[, `:=`(valid_denominator = is.finite(vehicle_mean - blank_mean) & vehicle_mean > blank_mean,
                vehicle_cv = vehicle_sd / vehicle_mean)]
write_table(controls, "transfer_plate_controls.tsv")
stopifnot(all(controls$valid_denominator))
raw <- merge(raw, controls, by = "plate_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(raw$blank_mean), !anyNA(raw$vehicle_mean))
raw[, viability := 100 * (INTENSITY - blank_mean) / (vehicle_mean - blank_mean)]
single <- raw[grepl("-S$", TAG) & is.finite(CONC),
              .(viability = mean(viability), sd = sd(viability), n_wells = .N,
                positions = paste(POSITION, collapse = ";")), by = .(plate_id, DRUG_ID, CONC)]
stopifnot(!anyDuplicated(single[, .(plate_id, DRUG_ID, CONC)]))
anchor <- raw[grepl("^A[0-9]+-C$", TAG),
              .(plate_id, POSITION, ANCHOR_ID = DRUG_ID, dose_A_uM = CONC)]
library <- raw[grepl("^L[0-9]+-D[0-9]+-C$", TAG),
               .(plate_id, POSITION, LIBRARY_ID = DRUG_ID, dose_B_uM = CONC,
                 viability_AB = viability, source_release, BARCODE, cell, assay_day)]
stopifnot(!anyDuplicated(anchor[, .(plate_id, POSITION)]), !anyDuplicated(library[, .(plate_id, POSITION)]))
combo <- merge(anchor, library, by = c("plate_id", "POSITION"))
unpaired <- rbindlist(list(anchor[!library, on = .(plate_id, POSITION)][, missing_record := "library"],
                           library[!anchor, on = .(plate_id, POSITION)][, missing_record := "anchor"]), fill = TRUE)
write_table(unpaired, "transfer_excluded_incomplete_well_annotations.tsv.gz")
message("Incomplete drug-pair annotation wells excluded: ", nrow(unpaired))
stopifnot(nrow(combo) > 0.9 * nrow(anchor), nrow(combo) > 0.9 * nrow(library))
combo[, COMBI_ID := paste(ANCHOR_ID, LIBRARY_ID, sep = ":")]
# Restrict to curves retained in the deposited fitted-release metadata (source QC), never to their fitted outcomes.
fit_join <- fit[, .(source_release, BARCODE, COMBI_ID, dose_A_uM = ANCHOR_CONC, drug_A, drug_B, pair,
                    heldout, target_A = ANCHOR_TARGET, target_B = LIBRARY_TARGET)]
combo <- merge(combo, fit_join, by = c("source_release", "BARCODE", "COMBI_ID", "dose_A_uM"))
single_A <- copy(single)
setnames(single_A, c("DRUG_ID", "CONC", "viability", "sd", "n_wells", "positions"),
         c("ANCHOR_ID", "dose_A_uM", "viability_A", "sd_A", "n_wells_A", "positions_A"))
single_B <- copy(single)
setnames(single_B, c("DRUG_ID", "CONC", "viability", "sd", "n_wells", "positions"),
         c("LIBRARY_ID", "dose_B_uM", "viability_B", "sd_B", "n_wells_B", "positions_B"))
combo <- merge(combo, single_A, by = c("plate_id", "ANCHOR_ID", "dose_A_uM"), all.x = TRUE)
combo <- merge(combo, single_B, by = c("plate_id", "LIBRARY_ID", "dose_B_uM"), all.x = TRUE)
combo[, `:=`(observation_id = paste(plate_id, POSITION, sep = "|"),
              curve_id = paste(plate_id, COMBI_ID, dose_A_uM, sep = "|"),
              function_time_h = 72, function_assay = "CellTiterGlo_ATP_relative_viability",
              dose_unit = "uM", explicit_biological_replicate_id_available = FALSE,
              replicate_note = "assay_day is date-created grouping; plate and well IDs do not establish biological independence",
              complete_measured_constituents = is.finite(viability_A) & is.finite(viability_B))]
stopifnot(!anyDuplicated(combo$observation_id), all(is.finite(combo$viability_AB)))
combo[, `:=`(bliss_expected_viability = viability_A * viability_B / 100,
              excess_inhibition = viability_A * viability_B / 100 - viability_AB)]
write_table(combo, "transfer_combination_observations.tsv.gz")
write_table(single, "transfer_measured_single_agent_wells.tsv.gz")
message("Observed combination wells: ", nrow(combo), "; complete measured constituents: ", sum(combo$complete_measured_constituents))

# Only directly measured 978-gene RNA, 24h, and constituent drugs present in the four-cell study.
metadata_path <- file.path(lincs_root, "GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz")
meta <- fread(metadata_path)
meta[, `:=`(cell = cell_key(cell_id), drug = drug_key(pert_iname))]
meta <- meta[cell %chin% cells & pert_time == 24]
drug_cells <- unique(rbind(combo[, .(cell, drug = drug_A)], combo[, .(cell, drug = drug_B)]))
keep_keys <- paste(drug_cells$cell, drug_cells$drug, sep = "|")
treatment_plates <- unique(meta[pert_type == "trt_cp" & paste(cell, drug, sep = "|") %chin% keep_keys, det_plate])
# Plate normalization uses every measured sample on selected plates, not an outcome-selected subset.
meta <- meta[det_plate %chin% treatment_plates]
stopifnot(all(meta$pert_type %chin% c("trt_cp", "ctl_vehicle")), all(meta$pert_time_unit == "h"))
plate <- meta[, .(n_vehicle = sum(pert_type == "ctl_vehicle"), n_samples = .N), by = .(cell, det_plate)]
plate[, eligible := n_vehicle >= 5]
write_table(plate, "transfer_RNA_plate_coverage.tsv")
meta <- meta[det_plate %chin% plate[eligible == TRUE, det_plate]]
stopifnot(all(meta[pert_type == "trt_cp", pert_dose_unit] == "um"))
meta[, `:=`(sample_id = trimws(inst_id), time_h = 24,
             source_series = sub("_.*", "", det_plate), dose_uM = pert_dose)]
# Source reporting sometimes rounds threefold-series nominal doses. This 3.5% rule
# merges only the known 10/3^k six-dose grid; all other measured doses remain exact.
grid <- c(0.04, 0.12, 0.37, 1.11, 3.33, 10)
for (g in grid) {
  eligible <- which(meta$pert_type == "trt_cp" & meta$pert_dose > 0)
  aligned <- eligible[abs(log(meta$pert_dose[eligible] / g)) <= 0.035]
  meta[aligned, dose_uM := g]
}
meta[pert_type == "ctl_vehicle", dose_uM := 0]
gctx <- file.path(lincs_root, "GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx")
ids <- trimws(h5read(gctx, "0/META/COL/id"))
meta[, gctx_column_index := match(sample_id, ids)]
stopifnot(!anyNA(meta$gctx_column_index), !anyDuplicated(meta$sample_id))
setorder(meta, gctx_column_index)
gene_ids <- as.character(h5read(gctx, "0/META/ROW/id"))
symbols <- unname(mapIds(org.Hs.eg.db, keys = gene_ids, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))
stopifnot(length(symbols) == 978L, !anyNA(symbols), !anyDuplicated(symbols))
expression <- matrix(NA_real_, nrow = 978L, ncol = nrow(meta), dimnames = list(symbols, meta$sample_id))
for (start in seq(1L, nrow(meta), by = 256L)) {
  keep <- start:min(start + 255L, nrow(meta))
  ix <- meta$gctx_column_index[keep]
  block <- h5read(gctx, "0/DATA/0/matrix", start = c(1L, min(ix)), count = c(978L, max(ix) - min(ix) + 1L))
  expression[, keep] <- block[, ix - min(ix) + 1L, drop = FALSE]
}
stopifnot(all(is.finite(expression)), all(expression >= 0))
# LINCS Level2 deconvolved fluorescence -> log2(1+x), then maintained limma
# quantile normalization across all included samples. The common distribution
# uses RNA only, never outcomes. The installed preprocessCore target-normalization
# binary cannot create its pthread even with R_THREADS=1; limma is serial here.
expression <- log2(1 + expression)
expression <- limma::normalizeQuantiles(expression)
dimnames(expression) <- list(symbols, meta$sample_id)
baseline_plate <- vapply(unique(meta$det_plate), function(p) rowMeans(expression[, meta$det_plate == p & meta$pert_type == "ctl_vehicle", drop = FALSE]), numeric(978L))
rownames(baseline_plate) <- symbols
displacement <- expression - baseline_plate[, match(meta$det_plate, colnames(baseline_plate)), drop = FALSE]
dimnames(displacement) <- dimnames(expression)
vehicle_residual <- vapply(unique(meta$det_plate), function(p) {
  max(abs(rowMeans(displacement[, meta$det_plate == p & meta$pert_type == "ctl_vehicle", drop = FALSE])))
}, numeric(1))
stopifnot(all(vehicle_residual < 1e-10))
baseline <- vapply(cells, function(c) rowMeans(baseline_plate[, unique(meta[cell == c, det_plate]), drop = FALSE]), numeric(978L))
rownames(baseline) <- symbols
meta[, selected_constituent := pert_type == "trt_cp" & paste(cell, drug, sep = "|") %chin% keep_keys]
meta[, profile_id := paste(cell, drug, dose_uM, "24h", sep = "|")]
profiles <- meta[selected_constituent == TRUE,
                 .(n_samples = .N, n_plates = uniqueN(det_plate), source_series = paste(sort(unique(source_series)), collapse = ";"),
                   sample_ids = paste(sample_id, collapse = ";"), replicate_type = "source_assay_replicates_not_assumed_biological"),
                 by = .(profile_id, cell, drug, dose_uM, time_h)]
profile_expr <- vapply(profiles$profile_id, function(p) rowMeans(expression[, meta$profile_id == p, drop = FALSE]), numeric(978L))
profile_delta <- vapply(profiles$profile_id, function(p) rowMeans(displacement[, meta$profile_id == p, drop = FALSE]), numeric(978L))
profile_se <- vapply(profiles$profile_id, function(p) {
  z <- displacement[, meta$profile_id == p, drop = FALSE]
  if (ncol(z) > 1L) apply(z, 1L, sd) / sqrt(ncol(z)) else rep(NA_real_, 978L)
}, numeric(978L))
rownames(profile_expr) <- rownames(profile_delta) <- rownames(profile_se) <- symbols
stopifnot(identical(colnames(profile_delta), profiles$profile_id), !anyDuplicated(profiles$profile_id),
          all(is.finite(profile_delta)), all(is.finite(baseline)))
qsave(list(expression = profile_expr, displacement = profile_delta, baseline = baseline,
           displacement_assay_replicate_se = profile_se, metadata = profiles,
           normalization = "log2(1+measured fluorescence); limma common quantile target across included RNA; subtract own-plate vehicle gene means"),
      file.path(out, "transfer_molecular.qs"), preset = "high")
write_table(profiles, "transfer_RNA_profiles.tsv.gz")
write_table(meta, "transfer_RNA_samples.tsv.gz")
write_table(data.table(gene_symbol = symbols, source_entrez_id = gene_ids, measured = TRUE), "transfer_measured_genes.tsv")

# Identify exact measured molecular doses or two bracketing doses. Never extrapolate.
queries <- unique(rbind(combo[, .(cell, drug = drug_A, dose_uM = dose_A_uM)],
                        combo[, .(cell, drug = drug_B, dose_uM = dose_B_uM)]))
links <- rbindlist(lapply(seq_len(nrow(queries)), function(i) {
  q <- queries[i]
  p <- profiles[cell == q$cell & drug == q$drug][order(dose_uM)]
  status <- "no_measured_RNA_drug_cell"
  low <- high <- NA_character_
  weight <- NA_real_
  if (nrow(p)) {
    exact <- which(abs(log(p$dose_uM / q$dose_uM)) <= 0.035)
    if (length(exact)) {
      j <- exact[which.min(abs(log(p$dose_uM[exact] / q$dose_uM)))]
      low <- high <- p$profile_id[j]
      weight <- 0
      status <- "measured_dose_within_3.5pct_rounding"
    } else if (q$dose_uM < min(p$dose_uM) || q$dose_uM > max(p$dose_uM)) status <- "outside_measured_RNA_dose_range"
    else {
      a <- max(which(p$dose_uM < q$dose_uM))
      b <- min(which(p$dose_uM > q$dose_uM))
      low <- p$profile_id[a]
      high <- p$profile_id[b]
      weight <- log(q$dose_uM / p$dose_uM[a]) / log(p$dose_uM[b] / p$dose_uM[a])
      status <- "linear_interpolation_in_logdose_between_measured_RNA_profiles"
    }
  }
  data.table(cell = q$cell, drug = q$drug, dose_uM = q$dose_uM, profile_low = low, profile_high = high,
             high_weight = weight, matching_status = status, supported = is.finite(weight))
}))
stopifnot(all(links$high_weight[links$supported] >= 0), all(links$high_weight[links$supported] <= 1))
write_table(links, "transfer_constituent_dose_links.tsv")
for (arm in c("A", "B")) {
  z <- copy(links)
  setnames(z, c("drug", "dose_uM"), c(paste0("drug_", arm), paste0("dose_", arm, "_uM")))
  other <- setdiff(names(z), c("cell", paste0("drug_", arm), paste0("dose_", arm, "_uM")))
  setnames(z, other, paste0(arm, "_", other))
  combo <- merge(combo, z, by = c("cell", paste0("drug_", arm), paste0("dose_", arm, "_uM")), all.x = TRUE)
}
combo[, molecular_supported := A_supported & B_supported & complete_measured_constituents]
write_table(combo, "transfer_molecular_links.tsv.gz")
coverage <- combo[, .(wells = .N, dose_contexts = uniqueN(paste(drug_A, dose_A_uM, drug_B, dose_B_uM)),
                      curves = uniqueN(curve_id), pairs = uniqueN(pair), assay_days = uniqueN(assay_day), plates = uniqueN(plate_id)),
                  by = .(cell, heldout, molecular_supported)]
write_table(coverage, "transfer_coverage.tsv")
repeat_coverage <- combo[, .(wells = .N, plates = uniqueN(plate_id), assay_days = uniqueN(assay_day),
                             releases = uniqueN(source_release)), by = .(cell, pair, drug_A, dose_A_uM, drug_B, dose_B_uM)]
write_table(repeat_coverage, "transfer_repeat_coverage.tsv.gz")
focal_coverage <- combo[pair == "pictilisib + trametinib",
                       .(functional_wells = .N, RNA_supported_wells = sum(molecular_supported),
                         assay_days = uniqueN(assay_day), RNA_supported_assay_days = uniqueN(assay_day[molecular_supported]),
                         functional_doses = paste(sort(unique(dose_B_uM)), collapse = ";")),
                       by = .(cell, source_release, drug_A, dose_A_uM, drug_B)]
write_table(focal_coverage, "transfer_focal_coverage.tsv")
# Independent reconstruction check only, after all matching/selection is fixed.
# The fitted-release anchor outcomes are never used as targets or selection rules.
reported_anchor <- rbindlist(lapply(source_paths, function(f) {
  z <- fread(f, select = c("BARCODE", "COMBI_ID", "ANCHOR_CONC", "ANCHOR_VIABILITY"))
  z[, `:=`(source_release = sub(" .*", "", basename(f)), dose_A_uM = suppressWarnings(as.numeric(ANCHOR_CONC)))]
  z
}))
anchor_check <- merge(unique(combo[, .(source_release, BARCODE, COMBI_ID, dose_A_uM, viability_A)]),
                      reported_anchor, by = c("source_release", "BARCODE", "COMBI_ID", "dose_A_uM"))
anchor_check <- anchor_check[, .(source_release, BARCODE, COMBI_ID, dose_A_uM,
                                 recomputed_viability_A = viability_A, source_reported_viability_A = 100 * ANCHOR_VIABILITY)]
write_table(anchor_check, "transfer_anchor_normalization_check.tsv.gz")
inputs <- c(source_paths, file.path(source_root, paste0(c("Original", "Validation"), "_screen_All_tissues_raw_data.csv.zip")),
            file.path(source_root, "Glossary.xlsx"), metadata_path, gctx,
            "/work/src/program_pilot/60_prepare_combination_transfer.R")
manifest <- data.table(path = inputs, bytes = file.info(inputs)$size, modified = as.character(file.info(inputs)$mtime))
manifest[, sha256 := vapply(path, function(p) if (file.info(p)$size < 3e8) digest(p, file = TRUE, algo = "sha256") else NA_character_, character(1))]
write_table(manifest, "transfer_input_manifest.tsv")
writeLines(c(
  "Primary source: Jaaks et al., Nature 2022, doi:10.1038/s41586-022-04437-2. Dose units uM from source Glossary.xlsx.",
  "Function is 72h ATP CellTiterGlo relative viability, not cell counts. Raw DURATION=4 denotes day4 after day1 treatment.",
  "Normalization: 100*(well intensity - mean blank B)/(mean vehicle NC-1 - mean blank B); no clipping.",
  "Combination wells contain separate anchor/library drug records sharing one physical intensity. These are merged once per well.",
  paste0("Source wells missing one drug annotation are excluded and listed: ", nrow(unpaired), " records across four cells."),
  "Constituent activities are actual same-plate same-dose single-drug wells, averaged where duplicated; no Hill-curve substitution.",
  "Only curves present in the source fitted-release metadata are retained. No fitted response/synergy value selects a curve.",
  "All MEK+PI3K target-family pairs, including dual PI3K/mTOR, are held out in every cell, dose, release and orientation.",
  "Multi-drug backbone records are excluded. Four-cell primary test is intervention transfer, not unseen-cell generalization.",
  "RNA is 24h directly measured 978-gene LINCS Level2; log2(1+x), limma common quantile target across included RNA, own-plate vehicle subtraction.",
  "Measured RNA samples are assay replicates, not assumed independent biological experiments. Condition means preserve sample IDs and series.",
  "Known threefold dose reporting differences within 3.5% are nominally aligned. Otherwise log-dose interpolation only inside the observed RNA range.",
  "All focal functional doses remain in raw tables, including molecular-unsupported doses; no extrapolated prediction is justified there.",
  "Jaaks source has plate, well and assay creation dates but no explicit biological replicate ID. Date groups are proxies, not proof of independence.",
  "Public historical outcomes have previously been available. Reserved-family evaluation is retrospective; no newly blinded claim.",
  paste0("Completed ", Sys.time(), "; observations ", nrow(combo), "; measured molecular profiles ", nrow(profiles))
), file.path(root, "TRANSFER_DATA_DESIGN.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "transfer_session_info.txt"))
print(coverage)
