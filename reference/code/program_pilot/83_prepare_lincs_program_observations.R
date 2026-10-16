#!/usr/bin/env Rscript
# Local measured-landmark LINCS evidence for the existing program bridge.
# The partial observation operator assumes within-program gene changes follow
# the frozen signed loading pattern. This assumption is audited, not hidden.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(rhdf5)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(limma)
  library(digest)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/83_lincs_and_tahoe_joint_program_bridge'
source_root <- '/work/v20_pilot_design_2026-10-09/source'
reference <- '/work/results/program_pilot/82_joint_program_rules_and_pair_interactions'
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
metadata_path <- file.path(source_root, 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz')
expression_path <- file.path(source_root, 'GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx')
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
archive <- readRDS(archive_path)
W <- archive$prep$W
P <- colnames(W)
native <- c(.05, .5, 5)
instances <- fread(file.path(reference, 'data/lincs_selected_panel_instances.tsv'))[
  in_CCLE == TRUE & pert_time == 24 & pert_type == 'trt_cp' & depmap_id != 'ACH-000248']
panel <- sort(unique(instances$drug))
stopifnot(length(panel) == 9L, !'Topotecan (hydrochloride)' %in% panel)
raw_metadata <- fread(metadata_path)
plates <- unique(instances$det_plate)
metadata <- raw_metadata[det_plate %in% plates]
stopifnot(all(metadata$pert_time == 24), all(metadata$pert_time_unit == 'h'))
plate_qc <- metadata[, .(samples = .N, vehicles = sum(pert_type == 'ctl_vehicle')), by = det_plate]
plate_qc[, eligible := vehicles >= 5L]
fwrite(plate_qc, file.path(root, 'data/lincs_plate_coverage.tsv'), sep = '\t')
instances <- instances[det_plate %in% plate_qc[eligible == TRUE, det_plate]]
metadata <- metadata[det_plate %in% instances$det_plate]
metadata[, column_index := match(inst_id, h5read(expression_path, '/0/META/COL/id'))]
stopifnot(!anyNA(metadata$column_index), !anyDuplicated(metadata$inst_id))
setorder(metadata, column_index)
gene_ids <- as.character(h5read(expression_path, '/0/META/ROW/id'))
symbols <- unname(mapIds(org.Hs.eg.db, keys = gene_ids, keytype = 'ENTREZID', column = 'SYMBOL', multiVals = 'first'))
stopifnot(length(symbols) == 978L, !anyNA(symbols), !anyDuplicated(symbols))
expression <- matrix(NA_real_, 978L, nrow(metadata), dimnames = list(symbols, metadata$inst_id))
message('Reading local LINCS: ', nrow(metadata), ' wells across ', uniqueN(metadata$det_plate), ' plates')
for (start in seq(1L, nrow(metadata), by = 512L)) {
  keep <- start:min(start + 511L, nrow(metadata))
  ix <- metadata$column_index[keep]
  block <- h5read(expression_path, '/0/DATA/0/matrix', start = c(1L, min(ix)),
    count = c(978L, max(ix) - min(ix) + 1L))
  expression[, keep] <- block[, ix - min(ix) + 1L, drop = FALSE]
}
stopifnot(all(is.finite(expression)), all(expression >= 0))
expression <- log2(1 + expression)
# Each plate is processed independently: RNA from held cells does not alter
# normalization of training cells. Unselected compounds serve only plate QC.
for (plate in unique(metadata$det_plate)) {
  ii <- which(metadata$det_plate == plate)
  expression[, ii] <- limma::normalizeQuantiles(expression[, ii, drop = FALSE])
  vehicle <- ii[metadata$pert_type[ii] == 'ctl_vehicle']
  expression[, ii] <- sweep(expression[, ii, drop = FALSE], 1L,
    rowMeans(expression[, vehicle, drop = FALSE]), '-')
  stopifnot(max(abs(rowMeans(expression[, vehicle, drop = FALSE]))) < 1e-10)
}
genes <- intersect(rownames(W), symbols)
partial <- t(expression[genes, , drop = FALSE]) %*% W[genes, , drop = FALSE]
gram <- crossprod(W)
stopifnot(max(abs(gram - diag(diag(gram)))) < 1e-12)
fraction <- colSums(W[genes, , drop = FALSE]^2) / diag(gram)
stopifnot(all(fraction > 0), all(fraction < 1))
# For gene movements in the frozen program span: observed partial score =
# loading-energy fraction * full score. No weights are silently renormalized.
observer <- diag(as.numeric(fraction), 40L)
dimnames(observer) <- list(P, P)
selected <- match(instances$inst_id, metadata$inst_id)
stopifnot(!anyNA(selected))
observations <- cbind(instances[, .(inst_id, depmap_id, cell_id, drug, dose = pert_dose, det_plate)],
  as.data.table(partial[selected, , drop = FALSE]))
grid <- c(.04, .12, .37, 1.11, 3.33, 10, 20)
for (dose in grid) observations[abs(log(observations$dose / dose)) <= .035, dose := dose]
plate_scores <- observations[, lapply(.SD, mean), by = .(depmap_id, cell_id, drug, dose, det_plate), .SDcols = P]
curves <- plate_scores[, c(list(plates = .N), lapply(.SD, mean)),
  by = .(depmap_id, cell_id, drug, dose), .SDcols = P]
native_profiles <- curves[, {
  ord <- order(dose)
  values <- vapply(P, function(p) approx(log(dose[ord]), get(p)[ord], xout = log(native), rule = 1)$y, numeric(3))
  stopifnot(all(is.finite(values)))
  cbind(data.table(dose = native), as.data.table(values))
}, by = .(depmap_id, cell_id, drug)]
native_profiles[, molecular_source := 'LINCS']
qsave(list(metadata = native_profiles, observer = observer, partial_scores = native_profiles[, P, with = FALSE],
  source = expression_path, gene_symbols = genes, original_weights = W,
  observation_model = 'partial = diag(landmark loading-energy fractions) * full program score; loading-pattern assumption',
  normalization = 'log2(1+Level2 fluorescence); within-plate limma quantile normalization; subtract own-plate vehicle mean',
  times_hours = 24, native_doses = native, panel = panel), file.path(root, 'data/lincs_program_observations.qs'))
fwrite(native_profiles, file.path(root, 'data/lincs_native_partial_program_observations.tsv'), sep = '\t')
fwrite(curves, file.path(root, 'data/lincs_measured_dose_partial_program_observations.tsv'), sep = '\t')
fwrite(data.table(program = P, genes_observed = sum(rownames(W) %in% genes),
  loading_energy_fraction = as.numeric(fraction)), file.path(root, 'data/lincs_program_observation_operator.tsv'), sep = '\t')

# Audit the approximation on actual historical gene movements, without any
# viability labels. This is a molecular diagnostic, not an operator fitted to Y.
gene_path <- '/work/almanac_conditionals/state_conditioned_bridge/sensitivity_selected/artifacts/eligible_gene_measurements.rds'
historical <- readRDS(gene_path)
gene_matrix <- historical$response[, rownames(W), drop = FALSE]
full_score <- gene_matrix %*% W
partial_score <- gene_matrix[, genes, drop = FALSE] %*% W[genes, , drop = FALSE]
predicted_partial <- full_score %*% t(observer)
audit <- data.table(program = P, measured_gene_fraction = length(genes) / nrow(W),
  approximation_pearson = vapply(1:40, function(p) cor(partial_score[, p], predicted_partial[, p]), numeric(1)))
fwrite(audit, file.path(root, 'data/partial_observation_molecular_audit.tsv'), sep = '\t')
fwrite(data.table(path = c(metadata_path, expression_path, archive_path, gene_path),
  sha256 = vapply(c(metadata_path, expression_path, archive_path, gene_path), digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/lincs_input_manifest.tsv'), sep = '\t')
print(native_profiles[, .(cells = uniqueN(depmap_id), native_profiles = .N), by = drug])
print(audit[, .(median_partial_observation_pearson = median(approximation_pearson), minimum = min(approximation_pearson))])
capture.output(sessionInfo(), file = file.path(root, 'lincs_preparation_sessionInfo.txt'))
message('LINCS molecular preparation complete; no functional outcomes loaded.')
