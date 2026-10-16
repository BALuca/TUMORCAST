# Shared settings for the sequential R program pilot.
# Run scripts from /work. Change this file for supported analysis settings;
# do not copy the pipeline into a new version directory.
work_dir <- "/work"
code_dir <- file.path(work_dir, "src/program_pilot")
output_dir <- file.path(work_dir, "results/program_pilot")
data_dir <- file.path(output_dir, "data")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
cache_dir <- file.path(work_dir, "cache/tahoe_prism_full_chain")
tahoe_dir <- file.path(cache_dir, "tahoe_pdex_selected")
source_dir <- file.path(work_dir, "v20_pilot_design_2026-10-09/source")
gene_audit_path <- file.path(work_dir,
  "almanac_conditionals/state_conditioned_bridge/program_first/results/gene_universe_audit.tsv")
tahoe_metadata_path <- file.path(work_dir,
  "almanac_conditionals/state_conditioned_bridge/program_first/results/full_tahoe_profile_metadata.tsv.gz")
cytosig_path <- file.path(work_dir, "data/01_CytoSig_data/diff.merge.gz")
cytosig_metadata_path <- file.path(work_dir, "data/01_CytoSig_data/meta_info.txt")

# No functional labels inform gene admission, molecular rank, or orientation.
held_drugs <- c("lapatinib", "neratinib")
excluded_cell <- "ACH-000248"
rank_grid <- c(2L, 4L, 8L, 12L, 20L, 40L, 80L)
n_starts <- 10L
onmf_iterations <- 500L
seed <- 20261009L
blas_threads <- 8L
minimum_program_genes <- 20L
minimum_start_agreement <- 0.50
n_molecular_folds <- 5L

# BT20 has a fixed baseline. Dose matching tolerates source rounding only.
cell_lines <- c("BT20", "HS578T", "MCF7", "MDAMB231", "SKBR3")
nominal_doses <- c(0.04, 0.12, 0.37, 1.11, 3.33, 10)
dose_log_tolerance <- 0.035
development_partners <- c("trametinib", "selumetinib", "saracatinib", "nvp-tae684")
penalty_grid <- c(0.01, 0.1, 1, 10, 100)

# Descriptive seven-drug continuation: local discovery uses the wider drug
# compendium in each of the seven measured recipients. This branch neither
# changes the earlier BT20 split nor fits functional outcomes.
local_drugs <- c("everolimus", "sirolimus", "erlotinib", "gefitinib",
  "lapatinib", "vemurafenib", "bortezomib")
local_cell_mapping_path <- file.path(work_dir,
  "almanac_tahoe_prism_core_2026-10-08/cell_mapping.tsv")
local_output_dir <- file.path(output_dir, "04_mockup_local_programs")
local_gene_limit <- 500L
# Hard variation admission precedes correlation and abundance-bin selection.
# SD is in signed log2(1+corrected CP10K) treatment-minus-vehicle units.
local_minimum_gene_sd <- 0.15
local_minimum_variance_ratio <- 2
local_rank_grid <- 2:80
# Ranks 2--5 remain visible but cannot be selected. Use the rightmost blue dot
# in the mean gene silhouette curve: positive silhouette and all seeds >=12 genes.
local_minimum_rank <- 6L
local_resamples <- 24L
local_minimum_genes <- 12L # Minimum seed size defining blue points eligible for selection.
local_display_genes_per_program <- 12L

# Source the existing heatmap library and the genuinely shared numerical helpers.
source(file.path(work_dir, "src/R/lib_heatmap.R"))
source(file.path(code_dir, "utilities.R"))
data.table::setDTthreads(blas_threads)
RhpcBLASctl::blas_set_num_threads(blas_threads)
for (path in c(output_dir, data_dir, figure_dir, model_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
