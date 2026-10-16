# Assemble independently measured molecular and functional MEK-response data.
# Public source measurements are preserved; no outcome is inferred from RNA.
library(data.table)
library(readxl)
library(xml2)
library(rvest)
library(digest)
library(qs)
library(edgeR)

root <- "/work/results/program_pilot/59_mechanistic_evidence_datasets"
source_dir <- file.path(root, "source")
data_dir <- file.path(root, "data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 300)

# Source-data packages provide the numeric phosphoprotein and growth measurements.
downloads <- data.table(
  study_id = c(rep("vanBentum2025", 2), rep("Macaya2023", 3)),
  source_id = c("mek_vanBentum_source", "mek_vanBentum_supplement_description",
                "mek_Macaya_RNA", "mek_Macaya_proteome", "mek_Macaya_source"),
  url = c(
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-59404-y/MediaObjects/41467_2025_59404_MOESM21_ESM.zip",
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-59404-y/MediaObjects/41467_2025_59404_MOESM3_ESM.docx",
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-41828-z/MediaObjects/41467_2023_41828_MOESM3_ESM.xlsx",
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-41828-z/MediaObjects/41467_2023_41828_MOESM4_ESM.xlsx",
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-41828-z/MediaObjects/41467_2023_41828_MOESM7_ESM.xlsx"
  ),
  extension = c("zip", "docx", "xlsx", "xlsx", "xlsx")
)
downloads[, path := file.path(source_dir, paste0(source_id, ".", extension))]
# GEO includes per-sample counts, while the article supplement contains full
# gene-symbol annotations and published differential-expression statistics.
geo_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161218/soft/GSE161218_family.soft.gz"
geo_path <- file.path(source_dir, "mek_GSE161218_family.soft.gz")
if (!file.exists(geo_path)) download.file(geo_url, geo_path, mode = "wb")
soft <- readLines(gzfile(geo_path))
count_urls <- sub("^!Sample_supplementary_file_1 = ftp:", "https:",
                  soft[grepl("^!Sample_supplementary_file_1", soft)])
additional <- data.table(
  study_id = c(rep("Macaya2023", 10), "vanBentum2025"),
  source_id = c("mek_GSE161218_family.soft", paste0("mek_", sub("\\.txt\\.gz$", "", basename(count_urls))),
                "mek_vanBentum_targeted_phospho"),
  url = c(geo_url, count_urls,
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-59404-y/MediaObjects/41467_2025_59404_MOESM7_ESM.csv"),
  extension = c("gz", rep("txt.gz", 9), "csv")
)
additional[, path := file.path(source_dir, paste0(source_id, ".", extension))]
downloads <- rbind(downloads, additional, fill = TRUE)
for (i in seq_len(nrow(downloads))) {
  if (!file.exists(downloads$path[i])) {
    download.file(downloads$url[i], downloads$path[i], mode = "wb", method = "libcurl")
  }
}
downloads[, bytes := file.info(path)$size]
downloads[, sha256 := vapply(path, digest, character(1), algo = "sha256", file = TRUE)]
downloads[, downloaded_utc := format(Sys.time(), tz = "UTC", usetz = TRUE)]
downloads[, license := "CC-BY-4.0 article; source files distributed with article"]
fwrite(downloads, file.path(data_dir, "mek_download_manifest.tsv"), sep = "\t")

zip_path <- downloads[source_id == "mek_vanBentum_source", path]
zip_contents <- as.data.table(unzip(zip_path, list = TRUE))
fwrite(zip_contents, file.path(data_dir, "mek_vanBentum_zip_inventory.tsv"), sep = "\t")
print(zip_contents)

# Read RNA counts and annotate using the same-study supplementary gene mapping.
rna_results <- rbindlist(lapply(c("Tram vs Ctrl", "Lest vs Ctrl"), function(sheet) {
  tab <- as.data.table(read_excel(downloads[source_id == "mek_Macaya_RNA", path], sheet = sheet))
  tab[, contrast := sheet]
  tab
}), fill = TRUE)
fwrite(rna_results, file.path(data_dir, "mek_Macaya_RNA_published_contrasts.tsv.gz"), sep = "\t")
gene_map <- unique(rna_results[, .(gene_id = sub("\\..*$", "", gene_id), gene = gene_name)])
stopifnot(!anyDuplicated(gene_map$gene_id))
count_paths <- downloads[grepl("^mek_GSM", source_id), path]
count_tables <- lapply(count_paths, fread)
stopifnot(length(count_tables) == 9L)
stopifnot(all(vapply(count_tables[-1], function(x) identical(x$Geneid, count_tables[[1]]$Geneid), logical(1))))
sample_names <- sub("^mek_GSM[0-9]+_(.*)\\.counts\\.txt\\.gz$", "\\1", basename(count_paths))
counts <- do.call(cbind, lapply(count_tables, function(x) x[[ncol(x)]]))
colnames(counts) <- sample_names
count_gene_ids <- sub("\\..*$", "", count_tables[[1]]$Geneid)
symbols <- gene_map$gene[match(count_gene_ids, gene_map$gene_id)]
mapped <- !is.na(symbols) & nzchar(symbols)
# Only same-study unambiguous gene-symbol mappings enter the scoring matrix.
# Counts from entries with the same symbol are summed; unmapped entries remain
# preserved in source counts and an explicit mapping audit.
mapping_audit <- data.table(source_gene_id = count_tables[[1]]$Geneid,
                           gene = symbols, mapped = mapped,
                           total_count = rowSums(counts))
fwrite(mapping_audit, file.path(data_dir, "mek_Macaya_RNA_gene_mapping.tsv.gz"), sep = "\t")
symbol_counts <- rowsum(counts[mapped, , drop = FALSE], symbols[mapped], reorder = TRUE)
stopifnot(all(is.finite(symbol_counts)), all(symbol_counts >= 0), !anyDuplicated(rownames(symbol_counts)))
dge <- calcNormFactors(DGEList(symbol_counts), method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 0.5)
rna_samples <- data.table(
  study_id = "Macaya2023", sample_id = sample_names, cell_line = "H1792",
  treatment = c(rep("DMSO", 3), rep("lestaurtinib", 3), rep("trametinib", 3)),
  trametinib_uM = c(rep(0, 6), rep(0.5, 3)),
  lestaurtinib_uM = c(rep(0, 3), rep(0.625, 3), rep(0, 3)),
  time_h = 16, replicate = rep(1:3, 3), assay = "RNA-seq", source_file = count_paths,
  matched_function = "H1792 72h viability; same drug doses; separate wells and later endpoint",
  missing_arm = "No combination RNA-seq and no resistant or MYC-overexpressing RNA-seq in GSE161218"
)
stopifnot(identical(colnames(log_cpm), rna_samples$sample_id))
displacement <- sweep(log_cpm, 1, rowMeans(log_cpm[, rna_samples$treatment == "DMSO", drop = FALSE]))
qsave(list(counts = symbol_counts, log_cpm = log_cpm, displacement = displacement,
           metadata = rna_samples, normalization = "TMM log2 CPM, prior count 0.5; displacement vs mean DMSO"),
      file.path(data_dir, "mek_Macaya_RNA_gene_matrices.qs"))
fwrite(rna_samples, file.path(data_dir, "mek_Macaya_RNA_samples.tsv"), sep = "\t")

# A consistent long schema keeps expression, protein and growth values separate.
measurement_parts <- list()
add_measurements <- function(tab) {
  stopifnot(all(c("study_id", "sample_id", "cell_line", "assay", "feature", "value", "unit",
                  "source_file", "source_panel") %in% names(tab)))
  stopifnot(all(is.finite(tab$value)))
  measurement_parts[[length(measurement_parts) + 1L]] <<- tab
}
source_xlsx <- downloads[source_id == "mek_Macaya_source", path]
f8 <- as.matrix(read_excel(source_xlsx, sheet = "Fig 8", col_names = FALSE, .name_repair = "minimal"))
# Author-provided independent experiment values; every dose including controls
# is retained, with no fitted dose-response or retrospectively selected cutoff.
dose_blocks <- data.table(panel = c("8F", "8G", "8H"), cell = c("H2009", "H358", "H358"),
                         drug = c("trametinib", "trametinib", "sotorasib"),
                         row_start = c(11L, 25L, 42L), row_end = c(18L, 35L, 52L))
for (b in seq_len(nrow(dose_blocks))) {
  block <- dose_blocks[b]
  rr <- block$row_start:block$row_end
  stopifnot(all(is.finite(as.numeric(f8[rr, 2]))))
  for (j in 3:10) {
    state <- if (j <= 6L) "LacZ_control" else "MYC_overexpression"
    rep_id <- if (j <= 6L) j - 2L else j - 6L
    tab <- data.table(study_id = "Macaya2023", cell_line = block$cell,
      biological_state = state, treatment = block$drug, dose_A_uM = as.numeric(f8[rr, 2]) / 1000,
      dose_B_uM = 0, time_h = 240, replicate = as.character(rep_id), assay = "cell_viability",
      feature = "viability", value = as.numeric(f8[rr, j]), unit = "percent genotype-matched vehicle",
      source_file = source_xlsx, source_panel = paste0("Fig", block$panel),
      source_row = rr, source_col = j, measurement_role = "functional consequence after genetic intervention")
    tab[, sample_id := paste(cell_line, biological_state, treatment, dose_A_uM, replicate, sep = "_")]
    add_measurements(tab)
  }
}
# Quantitative PCR directly measures MYC RNA in matched mono/combination arms;
# these are not a score derived from MYC targets or function.
fs6 <- as.matrix(read_excel(source_xlsx, sheet = "Fig Suppl 6", col_names = FALSE, .name_repair = "minimal"))
qpcr_blocks <- data.table(cell = c("H1792", "H2009", "H1792", "H2009"),
                         time_h = c(24, 24, 48, 48), row_start = c(29L, 29L, 37L, 37L),
                         col_start = c(3L, 9L, 3L, 11L), n_conditions = c(4L, 4L, 6L, 6L))
for (b in seq_len(nrow(qpcr_blocks))) {
  block <- qpcr_blocks[b]
  treatments <- if (block$n_conditions == 4) c("DMSO", "trametinib", "lestaurtinib", "trametinib+lestaurtinib") else
    c("DMSO", "trametinib", "lestaurtinib", "midostaurin", "trametinib+lestaurtinib", "trametinib+midostaurin")
  for (j in seq_along(treatments)) {
    rr <- block$row_start + 0:2
    cc <- block$col_start + j - 1L
    tr <- treatments[j]
    tab <- data.table(study_id = "Macaya2023", cell_line = block$cell, biological_state = "parental",
      treatment = tr, dose_A_uM = if (grepl("trametinib", tr)) 0.5 else 0,
      dose_B_uM = if (grepl("lestaurtinib|midostaurin", tr)) 0.625 else 0,
      time_h = block$time_h, replicate = as.character(1:3), assay = "qPCR", feature = "MYC",
      value = as.numeric(fs6[rr, cc]), unit = "relative MYC mRNA expression",
      source_file = source_xlsx, source_panel = if (block$time_h == 24) "SupFig6E" else "SupFig6F",
      source_row = rr, source_col = cc, measurement_role = "direct transcript measurement")
    tab[, sample_id := paste(cell_line, treatment, time_h, replicate, sep = "_")]
    add_measurements(tab)
  }
}
# Parental versus acquired-resistance MYC qPCR provides a second independently
# manipulated context; missing source cell is retained in an exclusions table.
qpcr_missing <- list()
for (j in 4:15) {
  cl <- if (j <= 9) "H1792" else "H2009"
  state <- if (j %in% c(4:6, 10:12)) "parental" else "trametinib_resistant"
  rep_id <- ((j - 4) %% 3) + 1L
  vals <- suppressWarnings(as.numeric(f8[3:4, j]))
  valid <- is.finite(vals)
  tab <- data.table(study_id = "Macaya2023", cell_line = cl, biological_state = state,
    treatment = c("DMSO", "trametinib")[valid], dose_A_uM = c(0, 0.5)[valid], dose_B_uM = 0,
    time_h = 24, replicate = as.character(rep_id), assay = "qPCR", feature = "MYC", value = vals[valid],
    unit = "relative MYC mRNA expression", source_file = source_xlsx, source_panel = "Fig8D",
    source_row = (3:4)[valid], source_col = j, measurement_role = "acquired resistance state")
  tab[, sample_id := paste(cell_line, biological_state, treatment, time_h, replicate, sep = "_")]
  add_measurements(tab)
  if (any(!valid)) qpcr_missing[[length(qpcr_missing) + 1L]] <- data.table(
    source_panel = "Fig8D", row = (3:4)[!valid], column = j, reason = "Blank cell in source workbook")
}
fwrite(rbindlist(qpcr_missing), file.path(data_dir, "mek_source_missing_values.tsv"), sep = "\t")

# Published 72-hour viability summaries pair with the H1792 molecular arms.
# They are summaries (one number per arm), not independent replicate values.
f2 <- as.matrix(read_excel(source_xlsx, sheet = "Fig 2", col_names = FALSE, .name_repair = "minimal"))
for (bb in list(list(rows = 3:7, cellcol = 2L, cols = 3:5, status = "KRAS_wild_type"),
                list(rows = 3:8, cellcol = 7L, cols = 8:10, status = "KRAS_mutant"))) {
  for (j in seq_along(bb$cols)) {
    tr <- c("trametinib", "lestaurtinib", "trametinib+lestaurtinib")[j]
    tab <- data.table(study_id = "Macaya2023", cell_line = sub(" cells", "", f2[bb$rows, bb$cellcol]),
      biological_state = bb$status, treatment = tr, dose_A_uM = if (grepl("trametinib", tr)) 0.5 else 0,
      dose_B_uM = if (grepl("lestaurtinib", tr)) 0.625 else 0, time_h = 72,
      replicate = "published_summary", assay = "cell_viability", feature = "viability",
      value = as.numeric(f2[bb$rows, bb$cols[j]]), unit = "percent vehicle",
      source_file = source_xlsx, source_panel = "Fig2A", source_row = bb$rows, source_col = bb$cols[j],
      measurement_role = "functional consequence; published summary, no replicate uncertainty")
    tab[, sample_id := paste(cell_line, treatment, time_h, sep = "_")]
    add_measurements(tab)
  }
}

# Preserve experimentally measured protein fold changes at 48 h. The first
# five contrasts have explicit conditions in the methods. No inferred pathway
# activation or functional label is substituted for these measurements.
px <- as.matrix(read_excel(downloads[source_id == "mek_Macaya_proteome", path],
                           col_names = FALSE, .name_repair = "minimal"))
protein_parts <- list()
for (j in seq(6L, 22L, by = 4L)) {
  tab <- data.table(protein_id = px[-(1:2), 1], gene = px[-(1:2), 4], contrast = px[1, j],
                    fold_change = as.numeric(sub(",", ".", px[-(1:2), j + 1L], fixed = TRUE)),
                    p_value = as.numeric(sub(",", ".", px[-(1:2), j + 2L], fixed = TRUE)))
  tab[, source_row := seq_len(.N) + 2L]
  protein_parts[[length(protein_parts) + 1L]] <- tab
}
protein_results <- rbindlist(protein_parts)
fwrite(protein_results, file.path(data_dir, "mek_Macaya_protein_contrasts.tsv.gz"), sep = "\t")
valid_proteins <- protein_results[!is.na(gene) & gene != "0" & is.finite(fold_change) & fold_change > 0]
add_measurements(valid_proteins[, .(study_id = "Macaya2023", sample_id = paste0("H1792_", contrast),
  cell_line = "H1792", biological_state = ifelse(contrast %in% c("TrvsC", "TrLvsC"), "trametinib_resistant", "parental"),
  treatment = contrast, time_h = 48, replicate = "published_contrast",
  assay = "TMT_proteomics", feature = gene, source_feature_id = protein_id,
  value = fold_change, unit = "protein abundance fold change vs parental DMSO",
  source_file = downloads[source_id == "mek_Macaya_proteome", path], source_panel = "Supplementary Data 2",
  source_row, measurement_role = "independent protein abundance")])

# van Bentum: complete factorial ligand x MEKi phosphoprotein measurements.
# Display-source files preserve numeric measurements rather than extracting
# values from images; source row numbers remain for every observation.
akt <- as.data.table(read.csv(unz(zip_path, "SourceData_Figure2D.csv")))
akt[, source_row := .I + 1L]
# Baseline BSA observations are repeated for display across ligand panels.
# Preserve them with the explicit ligand-panel key; do not count these as extra
# independent replicates in pooled analyses.
akt[, replicate := as.character(seq_len(.N)), by = .(celline, stimulus, inhibitor, ligandGroup)]
add_measurements(akt[, .(study_id = "vanBentum2025",
  sample_id = paste(celline, stimulus, inhibitor, ligandGroup, replicate, sep = "_"),
  cell_line = celline, biological_state = "parental", treatment = paste(inhibitor, stimulus, sep = "+"),
  dose_A_uM = ifelse(inhibitor == "MEKi", 1, 0), time_h = 4,
  ligand = stimulus, ligand_panel = ligandGroup, ligand_duration_min = 20,
  replicate, assay = "Luminex_phosphoprotein", feature = "AKT_phosphorylation", value,
  unit = "log2 phospho-AKT fold change vs untreated control",
  source_file = zip_path, source_panel = "Fig2D", source_row,
  measurement_role = "direct kinase signaling; BSA controls repeated between display panels")])
ph <- as.data.table(read.csv(unz(zip_path, "SourceData_Figure4E.csv")))
ph[, source_row := .I + 1L]
ph[, replicate := as.character(seq_len(.N)), by = .(Modified.sequence, CellLine, Ligandmix, MEKi)]
add_measurements(ph[, .(study_id = "vanBentum2025",
  sample_id = paste(CellLine, Ligandmix, MEKi, replicate, sep = "_"), cell_line = CellLine,
  biological_state = "parental", treatment = paste(MEKi, Ligandmix, sep = "+"),
  dose_A_uM = ifelse(MEKi == "MEKi", 10, 0), time_h = 4, ligand = Ligandmix, ligand_duration_min = 30,
  replicate, assay = "SPIED_DIA_phosphoproteomics", feature = HGNC_id,
  source_feature_id = Modified.sequence, value = norm_intensity,
  unit = "log10 cyclic-loess normalized precursor intensity", source_file = zip_path,
  source_panel = "Fig4E", source_row, measurement_role = "direct ERK/JNK phosphosite measurement")])
# Independent full dose grids directly test whether pharmacological blockade
# of the measured alternate route changes growth; single-agent controls remain.
for (panel in c("SourceData_Figure6Bleft.csv", "SourceData_Figure6Dleft.csv")) {
  gr <- as.data.table(read.csv(unz(zip_path, panel)))
  partner <- if (grepl("6B", panel)) "pictilisib" else "JNK-IN-8"
  partner_col <- if (grepl("6B", panel)) "PI3Ki" else "JNKi"
  gr[, partner_uM := get(partner_col)]
  gr[, source_row := .I + 1L]
  gr[, sample_id := paste(CellLine, partner, MEKi, partner_uM, Well, Replicate, sep = "_")]
  for (endpoint in c("GR_post", "norm_GR")) {
    tab <- gr[is.finite(get(endpoint)), .(study_id = "vanBentum2025", sample_id,
      cell_line = CellLine, biological_state = "parental", treatment = paste0("selumetinib+", partner),
      dose_A_uM = MEKi, dose_B_uM = partner_uM, time_h = 48, replicate = Replicate, well = Well,
      assay = "live_cell_imaging", feature = if (endpoint == "GR_post") "doubling_time" else "normalized_doubling_time",
      value = get(endpoint), unit = if (endpoint == "GR_post") "hours" else "fold matched no-inhibitor doubling time",
      source_file = zip_path, source_panel = panel, source_row,
      measurement_role = "functional consequence of inhibitor factorial")]
    add_measurements(tab)
  }
}
# The article displays exact author-computed MYC densitometry numbers above
# immunoblot lanes. Transcribe those annotations without image densitometry.
# These representative-blot summaries do not establish replicate uncertainty.
blot <- data.table(
  cell_line = rep(c("H1792", "H2009", "H1792", "HCC44"), each = 4),
  treatment = c(rep(c("DMSO", "trametinib", "lestaurtinib", "trametinib+lestaurtinib"), 2),
                rep(c("DMSO", "trametinib", "midostaurin", "trametinib+midostaurin"), 2)),
  value = c(1, 1, 1, 0.7, 1, 0.8, 0.9, 0.6, 1, 0.7, 0.7, 0.5, 1, 0.5, 0.6, 0.1),
  source_panel = rep(c("Fig7B", "Fig7B", "Fig7C", "Fig7C"), each = 4)
)
blot[, `:=`(study_id = "Macaya2023", biological_state = "parental", time_h = 48,
             assay = "immunoblot_densitometry", feature = "MYC", replicate = "representative_blot",
             unit = "relative MYC densitometry printed by authors", source_file = file.path(source_dir, "mek_Macaya_Fig7.png"),
             measurement_role = "direct protein; author annotations transcribed, no replicate uncertainty")]
blot[, sample_id := paste(cell_line, treatment, time_h, source_panel, sep = "_")]
blot[, dose_A_uM := ifelse(grepl("trametinib", treatment), 0.5, 0)]
blot[, dose_B_uM := ifelse(grepl("lestaurtinib|midostaurin", treatment), 0.625, 0)]
add_measurements(blot)
measurements <- rbindlist(measurement_parts, fill = TRUE)
measurements[, measurement_id := sprintf("MEK%07d", .I)]
stopifnot(!anyDuplicated(measurements$measurement_id), all(is.finite(measurements$value)))
fwrite(measurements, file.path(data_dir, "mek_measurements.tsv.gz"), sep = "\t")
myc_dose_summary <- measurements[assay == "cell_viability" & time_h == 240,
  .(mean_percent = mean(value), sd_percent = sd(value), n_experiments = uniqueN(replicate)),
  by = .(study_id, source_panel, cell_line, biological_state, treatment, dose_A_uM)]
fwrite(myc_dose_summary, file.path(data_dir, "mek_MYC_intervention_dose_summary.tsv"), sep = "\t")
# Multiple wells in a biological experiment are technical observations. Average
# them within the recorded replicate before any uncertainty calculation.
growth_replicates <- measurements[assay == "live_cell_imaging",
  .(value = mean(value), n_wells = uniqueN(well)),
  by = .(study_id, cell_line, treatment, dose_A_uM, dose_B_uM, replicate, feature, unit, time_h)]
fwrite(growth_replicates, file.path(data_dir, "mek_growth_biological_replicates.tsv"), sep = "\t")
sample_cols <- c("study_id", "sample_id", "cell_line", "biological_state", "treatment", "dose_A_uM",
                 "dose_B_uM", "time_h", "replicate", "assay", "source_file", "source_panel")
samples <- unique(measurements[, ..sample_cols])
fwrite(rbindlist(list(samples, rna_samples), fill = TRUE), file.path(data_dir, "mek_samples.tsv"), sep = "\t")
fwrite(measurements[, .(n_measurements = .N, n_cells = uniqueN(cell_line), n_features = uniqueN(feature)),
                     by = .(study_id, assay, source_panel)], file.path(data_dir, "mek_measurement_inventory.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(data_dir, "mek_sessionInfo.txt"))

# The study and contrast registry states which hypothesis each intervention
# distinguishes, and explicitly prevents non-matched assays being conflated.
studies <- data.table(
  study_id = c("Macaya2023", "vanBentum2025"),
  doi = c("10.1038/s41467-023-41828-z", "10.1038/s41467-025-59404-y"),
  primary_url = c("https://www.nature.com/articles/s41467-023-41828-z",
                  "https://www.nature.com/articles/s41467-025-59404-y"),
  repository = c("GSE161218; PXD024023", "PXD050961"),
  biology = c("MYC-dependent MEK/KRAS resistance; pharmacological combinations and MYC overexpression",
              "MEK feedback to AKT/JNK under growth-factor input; direct PI3K/JNK blockade"),
  downloaded_molecular = c("H1792 RNA-seq9samples16h; protein48h; MYC qPCR24/48h; printed MYC immunoblot values",
                          "pAKT Luminex factorial3cells4ligands; ERK/JNK phosphosites factorial3cells; complete targeted intensity source"),
  downloaded_functional = c("Tram/Lest72hviability11cells; MYC vsLacZ 10daydosecurves3cell-drugcontexts4experiments",
                            "MEKi x PI3Ki and MEKi x JNKi dose grids in HCT116/DLD1; 48h doubling-time summaries with replicate IDs"),
  causal_control = c("LacZ vs forced MYC; single agents and combinations; parental vs acquired resistance",
                     "DMSO/BSA, ligand alone, MEKi alone, ligand+MEKi; alternate-route inhibitor alone and with MEKi"),
  limitation = c("No combo RNA-seq or RNA-seq in genetic MYC arms; 16hRNA/48hprotein/72hfunction are separate assays; genotype panel not isogenic",
                 "No matched RNA; molecularMEKi10uM vs growth0.2/1/5uM; phosphoGFmix vs standard growth media; not exact-dose molecular-function pairs"),
  license = "CC-BY-4.0"
)
fwrite(studies, file.path(data_dir, "mek_studies.tsv"), sep = "\t")
contrasts <- data.table(
  contrast_id = c("MYC_gain_of_function", "RNA_vs_protein_under_combination", "MEK_ligand_feedback", "Alternate_route_blockade"),
  study_id = c("Macaya2023", "Macaya2023", "vanBentum2025", "vanBentum2025"),
  question = c("Does increasing MYC cause resistance to MEK/KRAS inhibition?",
               "Does MYC RNA follow the directly measured MYC protein change under combinations?",
               "Does MEK inhibition enhance alternate signaling despite ERK target engagement?",
               "Does blocking the measured alternate route increase functional MEKi effect?"),
  arms = c("LacZ vs MYC expression; each with complete matched trametinib/sotorasib dose curve including vehicle",
           "DMSO,trametinib,lestaurtinib,tram+lest; parallel midostaurin arms at48h",
           "BSA/DMSO,BSA/MEKi,ligand/DMSO,ligand/MEKi",
           "vehicle,MEKi,PI3Ki,MEKi+PI3Ki; separately vehicle,MEKi,JNKi,MEKi+JNKi full dose grids"),
  independent_molecular_measurement = c("Forced MYC construct confirmed by source WB; no RNA-seq in these genetic arms",
                                       "Direct MYC qPCR24/48h and independent immunoblot48h; proteome48h",
                                       "Quantitative pAKT Luminex; mass-spectrometry ERK/JNK phosphosites",
                                       "Independent ligand/MEKi phosphoprotein factorial in the same parental cells"),
  independent_functional_measurement = c("10dayviability four experiments in H2009/H358; sotorasib H358",
                                        "72hviability; independent from RNA and protein measurements",
                                        "Separate PI3K/JNK combination growth assays",
                                        "Direct time-lapse confluence and doubling time at each inhibitor combination"),
  planned_statistic = c("Within-experiment MYC-minus-LacZ viability at each dose; compare whole dose curves, do not select best dose",
                       "Within-cell combination-minus-single MYC RNA and protein contrasts; compare with independent viability; preserve timing",
                       "Interaction:(ligand+MEKi - BSA+MEKi)-(ligand+DMSO - BSA+DMSO) per cell/ligand/site",
                       "Partner increment in doubling time conditional on MEKi, minus partner-alone increment; retain all doses"),
  falsifier = c("No reproducible sensitivity shift after forced MYC at matched drug exposures",
                "RNA matches protein across matched arms, opposing a direct MYC RNA/protein mismatch; otherwise test downstream footprint separately",
                "No resolved alternate-route activation after MEKi, or ERK not suppressed",
                "Alternate-route blocker alone explains all growth change; no context-specific combination increment"),
  matching_and_missing_controls = studies$limitation[c(1,1,2,2)],
  evaluation_status = "Retrospective mechanistic reference; no model fitted and no held-out performance claim"
)
fwrite(contrasts, file.path(data_dir, "mek_contrasts.tsv"), sep = "\t")
# Register downloaded original figures as provenance for transcribed annotation.
fig_urls <- c(
  "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-023-41828-z/MediaObjects/41467_2023_41828_Fig7_HTML.png",
  "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-025-59404-y/MediaObjects/41467_2025_59404_Fig2_HTML.png")
fig_paths <- file.path(source_dir, c("mek_Macaya_Fig7.png", "mek_vanBentum_Fig2.png"))
for (i in seq_along(fig_paths)) {
  if (!file.exists(fig_paths[i])) download.file(fig_urls[i], fig_paths[i], mode = "wb")
}
fig_manifest <- data.table(study_id = c("Macaya2023", "vanBentum2025"), source_id = basename(fig_paths),
                           url = fig_urls, extension = "png", path = fig_paths,
                           bytes = file.info(fig_paths)$size,
                           sha256 = vapply(fig_paths, digest, character(1), algo = "sha256", file = TRUE),
                           downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), license = "CC-BY-4.0")
fwrite(rbind(downloads, fig_manifest, fill = TRUE), file.path(data_dir, "mek_download_manifest.tsv"), sep = "\t")
writeLines(c(
  "Mechanistic source data, not a fitted model.",
  "Macaya2023: GSE161218 count matrices were downloaded for 9 samples. Same-study gene symbol mapping, TMM log2 CPM, and vehicle displacement are saved in qs.",
  "RNA 16h, TMT protein 48h, and viability72h share H1792 and drug doses but are independent assays, not paired aliquots.",
  "Combination transcriptome absent. Matched combination MYC qPCR24/48h exists and was extracted directly from source workbook.",
  "MYC gain-of-function response curves have four independent experiments and both LacZ and MYC vehicle controls. No outcome was inferred from expression.",
  "16 immunoblot values are exact numbers printed by authors above Fig7B/C lanes, manually transcribed, not image-estimated. They are representative-blot summaries, not replicate observations.",
  "vanBentum2025: ligand x MEKi pAKT and ERK/JNK are direct molecular assays. PI3K/JNK inhibitor growth grids are independent functional assays in the same cell lines.",
  "Molecular Selumetinib10uM differs from functional0.2/1/5uM; cultures and ligand exposures also differ. Do not treat those as exact molecular-function matched samples.",
  "pAKT source rows include BSA controls duplicated across ligand display panels: retain panel key, never count those duplicated controls as independent evidence.",
  "No claim that these datasets prove the Tahoe transcriptional representation carries the functional mechanism; these give independently controlled measurements against which to test that representation."
), file.path(root, "MEK_DATA_README.txt"))
cat("Saved", nrow(measurements), "numeric measurements; RNA matrix", nrow(log_cpm), "genes x", ncol(log_cpm), "samples\n")
