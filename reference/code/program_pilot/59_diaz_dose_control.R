# Diaz2020 dose controls: prepare measured RNA and independent nuclei-count outcomes.
# Same analysis59 workflow. No model fitting or significance-based gene selection.
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(qs)
  library(Matrix)
})
root <- "/work/results/program_pilot/59_mechanistic_evidence_datasets"
archive <- file.path(root, "source/diaz_elife-52707-supp8-v2.zip")
source_dir <- file.path(root, "source/diaz_dose")
data_dir <- file.path(root, "data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
files <- c("self_synergy_fcounts_ucsc.mtrx", "pheno_ss.txt", "dose_ss.txt")
entries <- unzip(archive, list = TRUE)$Name
selected <- entries[basename(entries) %in% files]
stopifnot(length(selected) == 3)
unzip(archive, files = selected, exdir = source_dir)
paths <- setNames(file.path(source_dir, selected), basename(selected))

# Source featureCounts gene IDs are already gene symbols. Keep original symbols,
# including historical LOC names; do not ambiguously remap or discard them.
tab <- fread(paths[["self_synergy_fcounts_ucsc.mtrx"]], skip = 1)
stopifnot(names(tab)[1] == "Geneid", !anyDuplicated(tab$Geneid), all(nzchar(tab$Geneid)))
source_columns <- names(tab)[-(1:6)]
sample_ids <- sub(".*star_align/(RR[AB][0-9]+)_.*", "\\1", source_columns)
stopifnot(length(sample_ids) == 18, !anyDuplicated(sample_ids), all(grepl("^RR[AB][0-9]+$", sample_ids)))
counts <- as.matrix(tab[, ..source_columns])
rownames(counts) <- tab$Geneid
colnames(counts) <- sample_ids
stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == round(counts)))
pheno <- fread(paths[["pheno_ss.txt"]])
stopifnot(!anyDuplicated(pheno$Sample), setequal(sample_ids, pheno$Sample))
pheno <- pheno[match(sample_ids, Sample)]
samples <- data.table(
  study_id = "Diaz2020_MCF7_dose_control", sample_id = pheno$Sample, cell_line = "MCF7",
  drug = pheno$Drug, concentration_source = pheno[["Concentration (uM)"]],
  well = pheno$Well, RNA_yield_ng_ul = pheno[["RNA Yield (nG/uL)"]], RIN = pheno$RIN,
  exposure_hours = 24, DMSO_percent = 0.4,
  measurement = "RNA-seq", batch = "dose_experiment_about_two_years_after_combination_experiment"
)
samples[, dose_uM := suppressWarnings(as.numeric(concentration_source))]
samples[drug == "DMSO", dose_uM := 0]
stopifnot(all(is.finite(samples$dose_uM)), all(table(samples$drug, samples$concentration_source) %in% c(0, 2)))
samples[, condition_id := paste("Diaz_dose24h", drug, concentration_source, sep = "|")]
samples[, RNA_replicate := seq_len(.N), by = condition_id]
samples[, control_group_id := "Diaz_dose24h|DMSO|-"]
samples[, functional_measurement_id := paste0(condition_id, "|InCell_summary")]
samples[, functional_mean_repeated_in_source := TRUE]
samples[, sample_nature := "RNA replicate; biological independence not resolved; attached function is shared condition mean"]
controls <- which(samples$drug == "DMSO")
stopifnot(length(controls) == 2, identical(colnames(counts), samples$sample_id))

# TMM normalization is within this dose experiment only. Log2CPM minus the mean
# log2CPM of the two same-experiment DMSO controls yields signed displacement.
# No per-gene z-scaling, threshold clipping, outcome filtering, or joint-batch correction.
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
baseline <- rowMeans(log_cpm[, controls, drop = FALSE])
displacement <- sweep(log_cpm, 1, baseline, "-")
stopifnot(all(is.finite(displacement)), max(abs(rowMeans(displacement[, controls, drop = FALSE]))) < 1e-10)
log_expression_cp10k <- log2(1 + cpm(dge) / 100)
baseline_cp10k <- rowMeans(log_expression_cp10k[, controls, drop = FALSE])
displacement_cp10k <- sweep(log_expression_cp10k, 1, baseline_cp10k, "-")
stopifnot(max(abs(rowMeans(displacement_cp10k[, controls, drop = FALSE]))) < 1e-10)
samples[, `:=`(library_size = colSums(counts), TMM_factor = dge$samples$norm.factors)]
gene_manifest <- data.table(gene = rownames(counts), source_gene_id = rownames(counts),
                            total_counts = rowSums(counts), measured_samples = rowSums(counts > 0),
                            gene_mapping = "Source featureCounts gene symbol preserved")
qsave(list(counts = counts, genes = gene_manifest, samples = samples), file.path(data_dir, "diaz_dose_counts.qs"))
qsave(list(log2CPM = log_cpm, displacement = displacement, baseline_log2CPM = baseline,
           log2_one_plus_TMM_CP10K = log_expression_cp10k, displacement_CP10K = displacement_cp10k,
           baseline_log2_one_plus_TMM_CP10K = baseline_cp10k,
           samples = samples, genes = gene_manifest,
           normalization = paste("edgeR TMM within dose batch; original log2CPM prior.count1 retained;",
                                 "bank43 uses log2(1+TMM_CP10K); each subtracts own two DMSO controls")),
      file.path(data_dir, "diaz_dose_expression.qs"))
fwrite(samples, file.path(data_dir, "diaz_dose_samples.tsv"), sep = "\t", na = "NA")
fwrite(gene_manifest, file.path(data_dir, "diaz_dose_genes.tsv.gz"), sep = "\t", na = "NA")

# Fixed bank43 hard-cluster means use the identical CP10K scale as the combination
# adapter. Each batch keeps its own controls; this does not correct the two-year batch gap.
bank_path <- "/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs"
bank <- qread(bank_path)
coverage <- score_tables <- vector("list", ncol(bank$weights))
for (i in seq_len(ncol(bank$weights))) {
  program <- colnames(bank$weights)[i]
  members <- bank$genes[as.numeric(bank$weights[, program]) > 0]
  present <- intersect(members, rownames(displacement_cp10k))
  fraction <- length(present) / length(members)
  coverage[[i]] <- data.table(study_id = "Diaz2020_MCF7_dose_control", cell_line = "MCF7", program = program,
                              program_genes = length(members), measured_genes = length(present), fraction_measured = fraction)
  score_tables[[i]] <- data.table(sample_id = samples$sample_id, program = program,
    displacement = if (fraction >= 0.8) colMeans(displacement_cp10k[present, , drop = FALSE]) else NA_real_,
    baseline_reference = "own dose-experiment24h DMSO mean",
    score_status = "Fixed bank43 projection; partial membership renormalized if coverage >=80%",
    expression_scale = "log2(1+TMM_CP10K)")
}
fwrite(rbindlist(coverage), file.path(data_dir, "diaz_dose_bank43_coverage.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(score_tables), file.path(data_dir, "diaz_dose_bank43_scores.tsv"), sep = "\t", na = "NA")

# pheno_ss attaches the same summary to both RNA replicates. Deduplicate to one
# independently measured functional summary per drug/dose, with its reported SD.
attached <- copy(samples[, .(study_id, condition_id, functional_measurement_id, cell_line, drug,
                             concentration_source, dose_uM, exposure_hours)])
attached[, `:=`(mean = pheno[["InCell Viability"]], sd = pheno[["InCell Viability Stdev"]])]
stopifnot(attached[, all(vapply(.SD, uniqueN, integer(1)) == 1L),
                   by = condition_id, .SDcols = c("mean", "sd")]$V1)
functional <- unique(attached)
functional[, `:=`(assay = "InCell nuclei count", unit = "percent DMSO nuclei count", source_file = "pheno_ss.txt",
                  replicate_n = NA_integer_, statistic = "condition mean and SD; no individual function replicates supplied",
                  linkage = "Measured alongside RNA in the same dose experiment; paper Results/Fig1 supplement1B")]
stopifnot(nrow(functional) == 9, !anyDuplicated(functional$functional_measurement_id))
fwrite(functional, file.path(data_dir, "diaz_dose_matched_function.tsv"), sep = "\t", na = "NA")

# dose_ss is an additional functional dose experiment. Its values differ from
# pheno_ss and must never overwrite the matched RNA-associated measurements.
other <- fread(paths[["dose_ss.txt"]])
other_base <- data.table(study_id = "Diaz2020_dose_ss_function_separate", drug = other$Drug,
                        concentration_source = as.character(other[["Concentration (uM)"]]),
                        dose_uM = other[["Concentration (uM)"]], cell_line = "MCF7")
other_long <- rbindlist(list(
  cbind(copy(other_base), data.table(assay = "InCell nuclei count", mean = other[["InCell Viability"]], sd = other[["InCell Stdev"]])),
  cbind(copy(other_base), data.table(assay = "CellTiter-Glo ATP", mean = other[["CTG Viability"]], sd = other[["CTG Stdev"]]))
))
other_long[, `:=`(unit = "percent respective DMSO assay signal", source_file = "dose_ss.txt",
                  replicate_n = NA_integer_, exposure_hours = NA_real_,
                  linkage = "Separate functional table; source file omits time and experiment key; do not substitute for pheno_ss")]
other_long[, functional_measurement_id := paste(study_id, drug, concentration_source, assay, sep = "|")]
fwrite(other_long, file.path(data_dir, "diaz_dose_other_function.tsv"), sep = "\t", na = "NA")

# In the paper, 'sham combinations' mean one drug at a higher total concentration,
# not an observed combination of two different agents. Preserve that distinction.
sham <- data.table(
  drug = c("Tamoxifen", "Tamoxifen", "Mefloquine"),
  source_interpretation = c("5uM + 20uM interpreted as25uM total", "10uM + 20uM interpreted as30uM total",
                            "5uM + 10uM interpreted as15uM total"),
  total_dose_uM = c(25, 30, 15), constituent1_uM = c(5, 10, 5), constituent2_uM = c(20, 20, 10),
  RNA_available = c(TRUE, FALSE, TRUE),
  note = c("Single-drug dose-control; not two measured treatment components",
           "Paper reports too much cell death for RNA collection; do not invent transcriptome",
           "Single-drug dose-control; not two measured treatment components")
)
fwrite(sham, file.path(data_dir, "diaz_dose_sham_interpretations.tsv"), sep = "\t")
manifest <- data.table(file = c(archive, unname(paths), bank_path))
manifest[, `:=`(size_bytes = file.info(file)$size, md5 = unname(tools::md5sum(file)),
                 primary_url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7546737/",
                 supplement = "Supplementary file8 raw expression and dose annotation")]
fwrite(manifest, file.path(data_dir, "diaz_dose_source_manifest.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(data_dir, "diaz_dose_session_info.txt"))
writeLines(c(
  "Diaz2020 MCF7 dose controls:18RNA samples,9conditions, two RNA replicates per condition.",
  "RNA and matched InCell measurements24h. Tamoxifen5,10,20,25uM; mefloquine2.5,5,10,15uM; own DMSO.",
  "Independent function is9condition means with SD, NOT18independent viability observations.",
  "dose_ss supplies18additional dose conditions and two functional readouts; kept separate because values differ.",
  "Dose experiment ran about two years after combination experiment; separate batch and0.4%vs0.3%DMSO.",
  "Compare changes against their own controls. Do not claim a same-batch dose-matched contrast against combinations.",
  "All gene symbols retained; TMM log2CPM prior.count1, displacement subtracts mean ownDMSOlog2CPM.",
  "Bank43 projection separately uses log2(1+TMM_CP10K), exactly matching the combination adapter, own DMSO reference.",
  "Fixed hard-cluster scores require >=80% gene coverage; no batch combination or correction performed.",
  "RNA replicate biological independence is not resolved from sample metadata.",
  "30uM tamoxifen transcriptome is missing because excessive cell death prevented RNA collection.",
  paste("Counts dimensions:", nrow(counts), "genes x", ncol(counts), "samples"),
  "No functional fitting, differential-expression screening, gene selection, or plots performed."
), file.path(root, "DIAZ_DOSE_CONTROL_README.txt"))
cat("Prepared", nrow(counts), "genes x", ncol(counts), "RNA samples;", nrow(functional), "matched function summaries.\n")
