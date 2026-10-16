# 59: public ERBB ligand rescue data with explicit measurement and control provenance.
# This adapter prepares evidence; it does not fit a model or fill censored results with zero.
suppressPackageStartupMessages({
  library(data.table)
  library(curl)
  library(xml2)
  library(readxl)
  library(qs)
})

root <- "/work/results/program_pilot/59_mechanistic_evidence_datasets"
source_dir <- file.path(root, "source")
data_dir <- file.path(root, "data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 120)

# Download only compact processed tables and metadata, not the raw sequencing or MS files.
resources <- data.table(
  file = c("hgf_PMC11472020.html", "hgf_Shaban2024.html", "hgf_Shaban2024_S1.pdf",
           "hgf_Shaban2024_S2.xlsx", "hgf_Shaban2024_S3.xlsx", "hgf_PXD054180_project.json",
           "hgf_PXD054180_proteinGroups.txt", "hgf_PXD054180_supplement.zip", "hgf_PRJNA1029741_runs.tsv"),
  url = c(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC11472020/",
    "https://link.springer.com/article/10.1134/S000629792403009X",
    paste0("https://media.springernature.com/original/springer-static/esm/art%3A10.1134%2FS000629792403009X/",
           "MediaObjects/10541_2024_2658_MOESM1_ESM.pdf"),
    paste0("https://media.springernature.com/original/springer-static/esm/art%3A10.1134%2FS000629792403009X/",
           "MediaObjects/10541_2024_2658_MOESM2_ESM.xlsx"),
    paste0("https://media.springernature.com/original/springer-static/esm/art%3A10.1134%2FS000629792403009X/",
           "MediaObjects/10541_2024_2658_MOESM3_ESM.xlsx"),
    "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD054180",
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/09/PXD054180/proteinGroups.txt",
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11472020/supplementaryFiles",
    paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1029741&result=read_run&",
           "fields=run_accession,sample_accession,sample_title,experiment_title,fastq_ftp,fastq_bytes,fastq_md5&format=tsv")
  )
)
for (i in seq_len(nrow(resources))) {
  target <- file.path(source_dir, resources$file[i])
  if (!file.exists(target)) curl_download(resources$url[i], target, quiet = TRUE)
  stopifnot(file.info(target)$size > 0)
}
supp_dir <- file.path(source_dir, "hgf_PXD054180_supplement")
dir.create(supp_dir, showWarnings = FALSE)
unzip(file.path(source_dir, "hgf_PXD054180_supplement.zip"), files = "DataSheet1.ZIP", exdir = supp_dir)
inner <- file.path(supp_dir, "DataSheet1.ZIP")
entries <- unzip(inner, list = TRUE)$Name
unzip(inner, files = entries[grepl("\\.(xlsx|docx)$", entries)], exdir = supp_dir)

# RNA run aliases are the only deposited treatment-level sample identifiers found.
# Preserve donor identity; do not falsely call these six equally replicated arms.
runs <- fread(file.path(source_dir, "hgf_PRJNA1029741_runs.tsv"))
stopifnot(nrow(runs) == 12, !anyDuplicated(runs$sample_accession))
sample_records <- vector("list", nrow(runs))
for (i in seq_len(nrow(runs))) {
  id <- runs$sample_accession[i]
  target <- file.path(source_dir, paste0("hgf_", id, ".xml"))
  url <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", id, "?download=false")
  if (!file.exists(target)) curl_download(url, target, quiet = TRUE)
  alias <- xml_attr(xml_find_first(read_xml(target), "//SAMPLE"), "alias")
  token <- sub("^[0-9]+-", "", sub("_S.*", "", alias))
  arm <- switch(token, No = "control", Lap = "lapatinib", Egf = "EGF", Leg = "lapatinib+EGF",
                Sd = "serum_pooled", Lsd = "lapatinib+serum_pooled", Sd19 = "serum_DS19",
                Sd19L = "lapatinib+serum_DS19", Sd12 = "serum_DS12", stop("Unrecognized RNA alias: ", alias))
  sample_records[[i]] <- data.table(
    study_id = "Shaban2024_PRJNA1029741", sample_id = id, run_id = runs$run_accession[i],
    source_alias = alias, cell_line = "SK-BR-3", arm = arm, assay = "RNA-seq",
    lapatinib_nM = ifelse(grepl("lapatinib", arm), 150, 0),
    egf_ng_ml = ifelse(grepl("EGF", arm), 2, 0), serum_percent = ifelse(grepl("serum", arm), 5, 0),
    exposure_hours = 48, biological_replicate = NA_character_, technical_replicate = NA_character_,
    mapping_status = "Treatment inferred from deposited alias and corroborated by paper; replicate identity unspecified",
    numeric_data_status = "Raw FASTQ available; full processed count matrix not deposited in inspected supplements"
  )
  resources <- rbind(resources, data.table(file = basename(target), url = url))
}
samples <- rbindlist(sample_records)

# Protein LFQ values are numerical and unfiltered, but 01..12 do not have a verified arm key.
# Save those values without guessing that their order follows the paper's condition order.
protein <- fread(file.path(source_dir, "hgf_PXD054180_proteinGroups.txt"))
lfq_cols <- grep("^LFQ intensity ", names(protein), value = TRUE)
stopifnot(length(lfq_cols) == 36, all(c("Gene names", "Protein IDs", "Reverse", "Potential contaminant") %in% names(protein)))
protein[, source_row := .I]
protein_long <- melt(protein, id.vars = c("source_row", "Gene names", "Protein IDs", "Reverse", "Potential contaminant"),
                     measure.vars = lfq_cols, variable.name = "sample_id", value.name = "value", variable.factor = FALSE)
setnames(protein_long, c("Gene names", "Protein IDs"), c("feature", "source_protein_id"))
protein_long[, `:=`(study_id = "Kamashev2024_PXD054180", assay = "protein_LFQ",
                   unit = "MaxQuant LFQ intensity; zero denotes no quantified intensity",
                   value_status = fifelse(value == 0, "not_quantified_zero", "measured"),
                   arm = NA_character_, sample_id = sub("^LFQ intensity ", "", sample_id))]
protein_samples <- unique(protein_long[, .(study_id, sample_id)])
protein_samples[, `:=`(cell_line = "SK-BR-3", assay = "protein_LFQ", arm = NA_character_, exposure_hours = 48,
                      source_alias = sample_id, biological_replicate = sub("_.*", "", sample_id),
                      technical_replicate = sample_id,
                      mapping_status = "No verified source run-to-treatment key; do not assign by numerical order",
                      numeric_data_status = "Full LFQ values downloaded; arm assignment blocked by missing key")]
samples <- rbindlist(list(samples, protein_samples), fill = TRUE)
fwrite(protein_long, file.path(data_dir, "hgf_protein_LFQ_unmapped.tsv.gz"), sep = "\t", na = "NA")

# Published RNA coefficients are significance-censored. Keep every missing/censored cell explicitly.
rna <- as.data.table(read_excel(file.path(source_dir, "hgf_Shaban2024_S2.xlsx"), col_types = "text"))
rna <- rna[!is.na(SYMBOL)]
rna_cols <- names(rna)[3:9]
rna_long <- melt(rna[, c("SYMBOL", rna_cols), with = FALSE], id.vars = "SYMBOL", variable.name = "contrast_id",
                 value.name = "source_value", variable.factor = FALSE)
setnames(rna_long, "SYMBOL", "feature")
rna_long[, value := suppressWarnings(as.numeric(source_value))]
rna_long[, `:=`(study_id = "Shaban2024_PRJNA1029741", assay = "RNA_published_contrast",
                sample_id = contrast_id, unit = "author signed differential-expression coefficient; scale not explicit in XLSX",
                value_status = fifelse(is.finite(value), "reported_numeric",
                                      fifelse(grepl("padj", source_value), "censored_padj_above_0.05", "not_reported")))]
stopifnot(!anyDuplicated(rna_long[, .(feature, contrast_id)]))
rna_wide <- dcast(rna_long, feature ~ contrast_id, value.var = "value")
rna_matrix <- as.matrix(rna_wide[, -1])
rownames(rna_matrix) <- rna_wide$feature
qsave(list(values = rna_matrix, table = rna_long,
           warning = "This is a significance-censored contrast matrix, NOT a gene-by-sample expression matrix. NA is not zero."),
      file.path(data_dir, "hgf_RNA_censored_contrasts.qs"))

# Protein supplement uses explicit log2 fold changes and adjusted P values, also censored by selection.
protein_supp <- as.data.table(read_excel(file.path(supp_dir, "SKBR_suppementary_2024/Table S1.xlsx"), col_types = "text"))
protein_supp <- protein_supp[!is.na(SYMBOL)]
protein_contrasts <- c("lapatinib_vs_control", "lapatinib_EGF_vs_control", "lapatinib_serum_vs_control",
                       "EGF_vs_control", "serum_vs_control")
protein_effects <- rbindlist(lapply(seq_along(protein_contrasts), function(i) {
  col <- 12 + 2 * i
  data.table(study_id = "Kamashev2024_PXD054180", sample_id = protein_contrasts[i], contrast_id = protein_contrasts[i],
             feature = protein_supp$SYMBOL, assay = "protein_published_contrast", unit = "log2 fold change vs control",
             source_value = protein_supp[[col]], value = suppressWarnings(as.numeric(protein_supp[[col]])),
             adjusted_p = suppressWarnings(as.numeric(protein_supp[[col + 1]])))
}))
protein_effects[, value_status := fifelse(is.finite(value), "reported_numeric", "not_reported_after_DE_filter")]

# These are explicit numeric outcomes printed in the primary paper, NOT reconstructed replicates or digitized bars.
# Growth is increase in cell count relative to the untreated increase, not raw viable-cell abundance.
function_means <- data.table(
  study_id = "Kamashev2024_PXD054180", sample_id = c("lapatinib", "lapatinib+EGF", "lapatinib+serum"),
  feature = "growth_increment_relative_to_control", assay = "manual_cell_count_growth",
  value = c(6, 72, 63), unit = "percent untreated net cell-number increase",
  value_status = "paper_reported_group_summary_no_numeric_replicates", exposure_hours = 144,
  source_value = c("6%", "72%", "63%"),
  source_locator = "PMC11472020 section3.3/Fig2A; methods2.5 six-day growth assay; omics48h",
  caveat = "Results text references48h proteome treatment; separate growth methods specify6days; no per-replicate source values"
)
measurements <- rbindlist(list(rna_long, protein_effects, function_means), fill = TRUE)
fwrite(measurements, file.path(data_dir, "hgf_measurements.tsv.gz"), sep = "\t", na = "NA")
fwrite(samples, file.path(data_dir, "hgf_samples.tsv"), sep = "\t", na = "NA")

# Contrasts define what the experiment can discriminate before any program projection or outcome fitting.
contrasts <- data.table(
  contrast_id = c("drug_effect", "ligand_effect", "rescue_in_drug", "context_specific_ligand_effect",
                  "RNA_vs_protein_retention", "HGF_MET_specific_reversal"),
  test = c("L-C", "E-C", "LE-L", "(LE-L)-(E-C)",
           "Compare the same gene's RNA and protein drug/rescue contrasts", "((LHC-LC)-(LH-L))"),
  discriminates = c("Drug engagement/growth arrest", "Generic EGF stimulation", "EGF molecular/function restoration under drug",
                    "Specific modification of drug effect beyond ordinary ligand stimulation",
                    "RNA omission versus program compression for measured protein response",
                    "MET-dependent HGF protection versus direct MET-inhibitor toxicity"),
  measurement_requirement = c("C,L RNA and growth", "C,E RNA and growth", "L,LE RNA and growth",
                              "C,L,E,LE matched measurements", "Same conditions/time in RNA and protein with controls",
                              "L,LH,LC,LHC plus environment controls; HGF and MET-blocker molecular readouts"),
  status = c("Named published RNA/protein contrasts available; sample-level RNA needs quantification",
             "Named published contrasts available; significance-censored",
             "RNA direct contrast available; protein derived contrast needs unmapped LFQ arm key",
             "Design includes all four arms; full uncensored RNA needed for unbiased test",
             "Matched48h conditions published; RNA raw available; protein arm mapping unresolved",
             "Absent from this dataset; EGF is not HGF and no MET blocker arm exists")
)
fwrite(contrasts, file.path(data_dir, "hgf_contrasts.tsv"), sep = "\t", na = "NA")
studies <- data.table(
  study_id = c("Shaban2024_PRJNA1029741", "Kamashev2024_PXD054180", "Wilson2012", "Watson2018"),
  url = c("https://doi.org/10.1134/S000629792403009X", "https://doi.org/10.3389/fmolb.2024.1470496",
          "https://pmc.ncbi.nlm.nih.gov/articles/PMC3724525/", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5927625/"),
  setting = c("SK-BR-3 lapatinib +/- EGF or serum", "Same six conditions, protein48h and cell growth",
              "Matched drug+HGF+MET inhibitor signaling and survival examples", "HER2 subtype-specific HGF/NRG1 protection and blocker controls"),
  numeric_available = c("Gene DE coefficients with censoring;12raw RNA runs; no complete processed count matrix",
                        "Full protein LFQ; named filtered log2FC;3growth summary values transcribed from paper",
                        "Previously inspected source primarily blots/images; no quantitative phosphoprotein table verified",
                        "Existing MEMA cell-count/EdU screen; targeted crizotinib experiments not in MEMA matrix"),
  decisive_limit = c("Deposited RNA arms are unequally represented; paper replicate claim not matched by run manifest",
                     "Protein sample-to-condition key absent; no HGF/MET inhibitor; no growth replicate table",
                     "Do not use historical hand-coded +/-1 states as measured protein abundance",
                     "Baseline RNA GSE48213 is not drug+ligand RNA; MEMA alone does not test MET reversal"),
  license = c("Publisher rights reserved; public ENA raw sequencing", "PRIDE CC0; article/supplements CC BY",
              "Paper copyright Macmillan; no new numeric extraction", "Check release terms for existing MEMA files")
)
fwrite(studies, file.path(data_dir, "hgf_studies.tsv"), sep = "\t", na = "NA")
resources[, local_path := file.path(source_dir, file)]
resources[, `:=`(size_bytes = file.info(local_path)$size, md5 = unname(tools::md5sum(local_path)),
                 downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))]
fwrite(resources, file.path(data_dir, "hgf_download_manifest.tsv"), sep = "\t", na = "NA")
fwrite(samples[assay == "RNA-seq", .N, by = arm], file.path(data_dir, "hgf_RNA_arm_counts.tsv"), sep = "\t")
fwrite(measurements[, .(rows = .N, numerical_values = sum(is.finite(value))), by = .(study_id, assay, value_status)],
       file.path(data_dir, "hgf_measurement_coverage.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(data_dir, "hgf_session_info.txt"))
cat("Prepared", nrow(samples), "sample/run records;", sum(is.finite(measurements$value)),
    "published numerical measurements;", nrow(protein_long), "unmapped LFQ values.\n")
