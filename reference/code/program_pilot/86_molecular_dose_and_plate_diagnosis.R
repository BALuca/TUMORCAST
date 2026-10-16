#!/usr/bin/env Rscript
# Diagnose physical molecular dose trajectories, source arithmetic, DMSO plate
# background and independent Busulfan replicates. No functional model fitting.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
  library(arrow)
  library(digest)
})
setDTthreads(1L)
source('/work/src/program_pilot/utilities.R')
root <- '/work/results/program_pilot/86_molecular_dose_and_plate_diagnosis'
prior <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
cache <- '/work/cache/tahoe_prism_full_chain'
pdex <- file.path(cache, 'tahoe_pdex_selected')
bank <- '/work/almanac_conditionals/state_conditioned_bridge/program_first'
archive_path <- paste0('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/',
  'inputs/original/modular_pair_pilot/inputs.rds')
for (folder in c('data', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
P <- sprintf('P%02d', 1:40)
native_doses <- c(.05, .5, 5)
colors <- c('0.05' = '#1879a4', '0.5' = '#df9b32', '5' = '#ae3769')
cor_safe <- function(x, y) {
  if (length(x) < 3L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y)
}
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')

# Read only requested rows of historical NPY caches. Full control matrices use
# the existing R reader. All array indexing and gene orientation are explicit.
read_npy_rows <- function(path, rows) {
  con <- file(path, 'rb')
  on.exit(close(con))
  stopifnot(identical(readBin(con, 'raw', 6L), as.raw(c(147, 78, 85, 77, 80, 89))))
  version <- readBin(con, 'integer', 2L, size = 1L, signed = FALSE)
  size <- readBin(con, 'integer', 1L, size = if (version[1] == 1L) 2L else 4L,
    signed = FALSE, endian = 'little')
  header <- rawToChar(readBin(con, 'raw', size))
  dims <- as.integer(trimws(strsplit(sub(".*'shape': \\(([^)]*)\\).*", '\\1', header), ',')[[1]]))
  dims <- dims[!is.na(dims)]
  stopifnot(length(dims) == 2L, grepl("'fortran_order': False", header, fixed = TRUE),
    grepl("'<f[48]'", header), all(rows >= 0 & rows < dims[1]))
  bytes <- if (grepl("'<f8'", header, fixed = TRUE)) 8L else 4L
  offset <- seek(con)
  ans <- matrix(NA_real_, length(rows), dims[2])
  for (i in seq_along(rows)) {
    seek(con, offset + bytes * dims[2] * rows[i], origin = 'start')
    ans[i, ] <- readBin(con, 'numeric', dims[2], size = bytes, endian = 'little')
  }
  stopifnot(all(is.finite(ans)))
  ans
}

archive <- readRDS(archive_path)
observed <- fread(file.path(prior, 'data/observed_program_viability_points.tsv.gz'))
drugs <- unique(observed$drug)
molecular <- as.data.table(archive$prep$m)[drug %in% drugs]
W <- Matrix(as.matrix(archive$prep$W), sparse = TRUE)
stopifnot(identical(colnames(W), P), !anyDuplicated(molecular[, .(depmap_id, drug, dose)]))
full_gene_table <- fread(file.path(pdex, 'genes.tsv'))[order(gene_index)]
gene_index <- match(rownames(W), full_gene_table$feature)
stopifnot(!anyNA(gene_index), identical(full_gene_table$gene_index, 0:(nrow(full_gene_table) - 1L)))

# Replay independently cached full-compendium molecular coordinates. This is a
# source-calculation check, not an accuracy comparison against functional labels.
full_meta <- fread(file.path(bank, 'results/full_tahoe_profile_metadata.tsv.gz'))
full_meta[, matrix_row := .I - 1L]
source_rows <- merge(molecular[, .(depmap_id, drug, dose, source)], full_meta,
  by = c('depmap_id', 'drug', 'dose'))
gene_order <- fread(file.path(bank, 'results/full_tahoe_gene_order.tsv'))[order(matrix_column), gene]
gene_displacement <- read_npy_rows(file.path(bank, 'artifacts/full_tahoe_candidate_responses.npy'),
  source_rows$matrix_row)
stopifnot(ncol(gene_displacement) == length(gene_order), all(rownames(W) %in% gene_order))
replayed <- as.matrix(gene_displacement[, match(rownames(W), gene_order), drop = FALSE] %*% W)
old_keys <- paste(molecular$depmap_id, molecular$drug, molecular$dose)
source_keys <- paste(source_rows$depmap_id, source_rows$drug, source_rows$dose)
stored <- as.matrix(molecular[match(source_keys, old_keys), ..P])
replay_difference <- replayed - stored
source_rows[, maximum_program_difference := apply(abs(replay_difference), 1L, max)]
write_table(source_rows[, .(depmap_id, drug, dose, source, contributing_plates,
  dose_unit, n_plates, maximum_program_difference)], 'independent_cache_replay')
replay_summary <- source_rows[, .(profiles = .N,
  maximum_program_difference = max(maximum_program_difference),
  median_program_difference = median(maximum_program_difference)), by = drug]
write_table(replay_summary, 'independent_cache_replay_summary')
print(replay_summary)
cat('Source replay complete.\n')

# Reproduce the original trimmed paired normalizer in R, including the full
# 62,710-gene library totals and scaling-gene eligibility, before frozen scoring.
# This exact routine is reused for DMSO background and separate treated replicates.
paired_scores <- function(treated, reference) {
  stopifnot(length(treated) == nrow(full_gene_table), length(reference) == length(treated),
    all(is.finite(treated)), all(is.finite(reference)), min(treated) >= 0, min(reference) >= 0)
  tcp <- treated / sum(treated) * 10000
  rcp <- reference / sum(reference) * 10000
  eligible <- tcp >= .1 & rcp >= .1
  stopifnot(sum(eligible) >= 100L)
  log_t <- log2(tcp[eligible])
  log_r <- log2(rcp[eligible])
  m <- log_t - log_r
  a <- .5 * (log_t + log_r)
  mlim <- quantile(m, c(.30, .70), names = FALSE)
  alim <- quantile(a, c(.05, .95), names = FALSE)
  retained <- m >= mlim[1] & m <= mlim[2] & a >= alim[1] & a <= alim[2]
  stopifnot(sum(retained) >= 100L)
  factor <- 2^median(m[retained])
  delta <- log2(1 + tcp[gene_index] / factor) - log2(1 + rcp[gene_index])
  list(score = as.numeric(delta %*% W), factor = factor,
    treated_total = sum(treated), reference_total = sum(reference))
}

reference_meta <- fread(file.path(pdex, 'pdex_plate_reference_metadata.tsv'))
controls_full <- read_npy_matrix(file.path(pdex, 'pdex_plate_reference.npy'))
stopifnot(nrow(controls_full) == nrow(reference_meta), ncol(controls_full) == nrow(full_gene_table))
reference_meta[, array_index := .I]
ids <- sort(unique(observed$depmap_id))
reference_meta <- reference_meta[depmap_id %in% ids]
stopifnot(all(reference_meta[, .N, by = depmap_id]$N == 14L))

# Each control plate is compared against the other thirteen plates of its own
# cell line. This is inter-plate background variability, not within-well sampling
# error. It cannot alone establish noise in same-plate treatment/control contrasts.
background <- list()
for (i in seq_len(nrow(reference_meta))) {
  current <- reference_meta[i]
  other <- reference_meta[depmap_id == current$depmap_id & plate != current$plate, array_index]
  result <- paired_scores(controls_full[current$array_index, ], colMeans(controls_full[other, , drop = FALSE]))
  background[[i]] <- cbind(current[, .(depmap_id, plate)], as.data.table(as.list(setNames(result$score, P))))
}
background_wide <- rbindlist(background)
background_long <- melt(background_wide, id.vars = c('depmap_id', 'plate'), variable.name = 'program',
  value.name = 'displacement', variable.factor = FALSE)
background_summary <- background_long[, .(plates = .N, mean_background = mean(displacement),
  sd_background = sd(displacement), background_low = quantile(displacement, .025),
  background_high = quantile(displacement, .975)), by = .(depmap_id, program)]
write_table(background_long, 'dmso_background_identical_normalizer')
write_table(background_summary, 'dmso_background_per_line_program')
with_background <- merge(observed, background_summary, by = c('depmap_id', 'program'))
stopifnot(nrow(with_background) == nrow(observed))
diagnostic <- with_background[, .(cells = .N, mean_displacement = mean(observed_displacement),
  sd_displacement = sd(observed_displacement), median_background_sd = median(sd_background),
  outside_background_range = sum(observed_displacement < background_low | observed_displacement > background_high)),
  by = .(drug, dose, program)]
write_table(diagnostic, 'drug_dose_program_vs_interplate_background')
print(diagnostic[drug == 'Lapatinib ditosylate' & program %in% c('P11', 'P21', 'P40')])
print(diagnostic[drug == 'Busulfan' & program == 'P12'])
cat('Identically normalized DMSO background complete.\n')

# Whole dose paths, including all forty frozen programs. Per-cell matching avoids
# interpreting different recipient sets at high and low dose as dose progression.
mean_rows <- observed[, .(mean_displacement = mean(observed_displacement),
  sd_displacement = sd(observed_displacement), cells = .N), by = .(drug, dose, program)]
mean_paths <- dcast(mean_rows, drug + program ~ dose, value.var = 'mean_displacement')
setnames(mean_paths, c('0.05', '0.5', '5'), c('low', 'middle', 'high'))
mean_paths[, middle_outside_endpoints := middle < pmin(low, high) | middle > pmax(low, high)]
write_table(mean_paths, 'all_program_mean_dose_paths')
paired_paths <- dcast(observed, depmap_id + drug + program ~ dose, value.var = 'observed_displacement')
setnames(paired_paths, c('0.05', '0.5', '5'), c('low', 'middle', 'high'))
complete <- paired_paths[is.finite(low) & is.finite(middle) & is.finite(high)]
path_summary <- complete[, .(cells = .N,
  middle_outside_endpoints = sum(middle < pmin(low, high) | middle > pmax(low, high)),
  low_to_middle_r = cor_safe(low, middle), middle_to_high_r = cor_safe(middle, high),
  low_to_high_r = cor_safe(low, high)), by = .(drug, program)]
write_table(path_summary, 'same_cell_program_dose_paths')

# Correlations of forty-program dose-change vectors for the existing ten drugs.
# Shared-plate similarities are diagnostic associations, not causal batch estimates.
plate_groups <- molecular[dose == .05, .(primary_plate = min(as.integer(
  unlist(strsplit(contributing_plates, ',', fixed = TRUE))))), by = drug]
plate_groups[, primary_plate_group := (primary_plate - 1L) %/% 3L + 1L]
mean_paths <- merge(mean_paths, plate_groups, by = 'drug')
mean_paths[, `:=`(lower_change = middle - low, upper_change = high - middle)]
drug_pairs <- combn(drugs, 2L)
trajectory_correlations <- list()
for (i in seq_len(ncol(drug_pairs))) {
  aa <- mean_paths[drug == drug_pairs[1, i]][order(program)]
  bb <- mean_paths[drug == drug_pairs[2, i]][order(program)]
  trajectory_correlations[[i]] <- data.table(drugA = drug_pairs[1, i], drugB = drug_pairs[2, i],
    same_primary_plates = aa$primary_plate_group[1] == bb$primary_plate_group[1],
    low_to_middle_r = cor_safe(aa$lower_change, bb$lower_change),
    middle_to_high_r = cor_safe(aa$upper_change, bb$upper_change))
}
trajectory_correlations <- rbindlist(trajectory_correlations)
write_table(trajectory_correlations, 'cross_drug_dose_trajectory_correlations')
print(trajectory_correlations[same_primary_plates == TRUE])
cat('Dose trajectory comparisons complete.\n')
