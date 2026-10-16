#!/usr/bin/env Rscript
# Outcome-informed descriptive pilot. Gene admission uses retrospective PRISM
# and ALMANAC correlations; this is not independent predictive validation.
# Same pipeline as stage 05, with Ward clustering on sample/sample correlation.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(cluster)
  library(patchwork)
})
source("/work/src/program_pilot/config.R")
arguments <- commandArgs(trailingOnly = TRUE)
check(length(arguments) == 0L || identical(arguments, "--render-only"),
  "Supported invocation: Rscript 06_mockup_sample_programs.R [--render-only]")
render_only <- identical(arguments, "--render-only")
# This separately requested experiment keeps its outputs and gene rule local.
# Test each drug/dose row equally; pool doses without dose adjustment or averaging.
outcome_p_threshold <- 0.05 # Restore the earlier nominal significance threshold.
outcome_minimum_abs_correlation <- 0 # No effect-size floor in the earlier lenient screen.
output_dir <- file.path(work_dir, "results/program_pilot/06_mockup_sample_programs")
data_dir <- file.path(output_dir, "data")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
for (path in c(output_dir, data_dir, figure_dir, model_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
set.seed(seed)
options(warn = 1)
theme_set(theme_classic(base_size = 10))
message("Outcome-selected displacement pilot; seed ", seed, "; results: ", output_dir)

# The panel's aliases are explicit; no salt or synonym inference is performed.
panel_name <- function(x) {
  out <- canonical_drug(x)
  out[out == "rapamycin"] <- "sirolimus"
  out
}

# Weighted Pearson gene correlation. Equal drug mass and molecular/control
# support weights preserve weak conditions without inflating them by normalization.
# x: conditions x genes; weights: nonnegative vector of condition mass.
gene_correlation <- function(x, weights) {
  check(nrow(x) == length(weights) && all(is.finite(x)) &&
    all(is.finite(weights)) && all(weights >= 0) && sum(weights) > 0,
    "Invalid gene-correlation inputs")
  weights <- weights / sum(weights)
  centered <- sweep(x, 2L, colSums(x * weights), "-")
  covariance <- crossprod(centered * sqrt(weights))
  scale <- sqrt(diag(covariance))
  check(all(scale > 0), "Constant gene in local discovery")
  out <- covariance / outer(scale, scale)
  out[out < -1] <- -1
  out[out > 1] <- 1
  diag(out) <- 1
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

# Deterministic block numbering follows the original gene dendrogram. Renumber
# hard seed partitions solely for presentation; soft weights are saved separately.
ordered_partition <- function(tree, k) {
  partition <- cutree(tree, k)
  block_order <- unique(partition[tree$order])
  setNames(match(partition, block_order), names(partition))
}

# Reused by discovery and cached-bank validation. Select the rightmost blue
# point in the mean gene silhouette curve: the largest allowed rank with a
# positive mean silhouette and every seed at least local_minimum_genes large.
# Cophenetic, PAC and resampling agreement remain descriptive diagnostics.
select_silhouette <- function(rank_table) {
  required <- c("rank", "minimum_size", "mean_silhouette")
  check(all(required %in% names(rank_table)), "Missing silhouette rank-selection columns")
  check(nrow(rank_table) > 0L && !anyDuplicated(rank_table$rank) &&
    all(is.finite(rank_table$rank)) && all(rank_table$rank >= 2L) &&
    all(rank_table$rank == as.integer(rank_table$rank)) &&
    all(is.finite(rank_table$minimum_size)) && all(rank_table$minimum_size >= 0L) &&
    all(is.na(rank_table$mean_silhouette) | abs(rank_table$mean_silhouette) <= 1),
    "Invalid silhouette rank-selection diagnostics")
  out <- copy(rank_table)
  out[, `:=`(rank_allowed = rank >= local_minimum_rank,
    minimum_size_qc_pass = minimum_size >= local_minimum_genes)]
  out[, silhouette_candidate := rank_allowed & minimum_size_qc_pass &
    is.finite(mean_silhouette) & mean_silhouette > 0]
  candidates <- out[silhouette_candidate == TRUE, rank]
  selected_rank <- if (length(candidates)) max(candidates) else NA_integer_
  status <- if (length(candidates)) "selected" else
    "unresolved: no allowed size-qualified cut with positive mean gene silhouette"
  out[, `:=`(selected = !is.na(selected_rank) & rank == selected_rank,
    selection_status = status,
    selection_rule = "last blue dot in mean gene silhouette: largest allowed rank with positive silhouette and minimum seed size",
    cophenetic_change_to_next_rank = c(diff(consensus_cophenetic), NA_real_))]
  out
}

# Reused for the separate PRISM and ALMANAC tests in every cell line. Complete
# matched drug/dose rows have equal weight; missing assay values are excluded.
# Return two-sided Pearson tests on original signed gene displacement versus
# vehicle-relative abundance. BH correction is reported within each assay/cell
# for reference, while admission uses nominal significance without an effect-size floor.
gene_outcome_tests <- function(x, outcome, drugs) {
  check(nrow(x) == length(outcome) && length(drugs) == nrow(x) &&
    all(is.finite(x)) && !anyNA(drugs), "Invalid gene-outcome correlation inputs")
  check(all(is.na(outcome) | is.finite(outcome)), "Nonfinite measured outcome")
  matched <- which(!is.na(outcome))
  y <- outcome[matched]
  results <- lapply(seq_len(ncol(x)), function(j) {
    values <- x[matched, j]
    status <- if (length(matched) < 3L) "fewer than three matched conditions" else
      if (sd(y) == 0) "constant outcome" else if (sd(values) == 0) "constant gene displacement" else "tested"
    test <- if (status == "tested") cor.test(values, y, method = "pearson", alternative = "two.sided") else NULL
    data.table(gene = colnames(x)[j], matched_conditions = length(matched),
      matched_drugs = uniqueN(drugs[matched]), correlation = if (is.null(test)) NA_real_ else unname(test$estimate),
      p_value = if (is.null(test)) NA_real_ else test$p.value, test_status = status)
  })
  rbindlist(results)
}

# Map sample clusters to gene groups for the existing gene silhouette, weights
# and plots. Assign each gene to its strongest positive weighted Pearson
# association with a sample-cluster indicator. Empty groups fail the same size QC.
# Reused for full discovery, whole-drug resamples and the centering diagnostic.
genes_from_sample_clusters <- function(x, sample_partition, weights) {
  weights <- weights / sum(weights)
  weighted_x <- x * weights
  mass <- as.numeric(rowsum(weights, sample_partition, reorder = TRUE))
  cross <- t(rowsum(weighted_x, sample_partition, reorder = TRUE)) - outer(colSums(weighted_x), mass)
  # The positive gene SD is common to all groups for a gene and cancels in
  # argmax. Divide by indicator SD only; this gives the same Pearson assignment.
  association <- sweep(cross, 2L, sqrt(mass * (1 - mass)), "/")
  check(all(is.finite(association)), "Undefined gene-to-sample-cluster association")
  partition <- max.col(association, ties.method = "first")
  check(all(apply(association, 1L, max) > 0), "Gene has no positive sample-cluster association")
  setNames(partition, colnames(x))
}

# Freeze a program dictionary after outcome-informed gene admission and assess
# collection dependence. Resampling entire drugs retains related doses together;
# the gene filter is fixed, so this is not a resampling validation of supervision.
discover_local <- function(x, null, metadata, genes, unit_column, scope, orientation) {
  x <- x[, genes, drop = FALSE]
  null <- null[, genes, drop = FALSE]
  units <- metadata[[unit_column]]
  weights <- group_weights(units) * metadata$condition_support
  correlation <- gene_correlation(x, weights)
  # Signed-correlation Euclidean distance permits Ward's variance partition.
  # Average linkage produced only inadmissible tiny fragments in HT-29.
  distance <- as.dist(sqrt(2 * (1 - correlation)))
  tree <- hclust(distance, method = "ward.D2")
  # Only the discovery axis changes: same Pearson distance and Ward method,
  # now between drug/dose samples across the selected genes.
  condition_correlation <- cor(t(x))
  check(all(is.finite(condition_correlation)), "Undefined sample/sample Pearson correlation")
  condition_correlation[condition_correlation > 1] <- 1
  condition_correlation[condition_correlation < -1] <- -1
  diag(condition_correlation) <- 1
  condition_tree <- hclust(as.dist(sqrt(2 * (1 - condition_correlation))), method = "ward.D2")
  # Explore high cuts even when their seed sizes are unsupported. Size is an
  # explicit eligibility check, not an upper truncation hiding the rank curve.
  grid <- local_rank_grid[local_rank_grid < min(ncol(x), nrow(x))]
  check(length(grid) > 0, paste("No feasible ranks in", scope))
  condition_partitions <- lapply(grid, function(k) ordered_partition(condition_tree, k))
  names(condition_partitions) <- as.character(grid)
  partitions <- lapply(condition_partitions, function(partition) genes_from_sample_clusters(x, partition, weights))
  names(partitions) <- as.character(grid)
  rank_table <- rbindlist(lapply(seq_along(grid), function(j) {
    partition <- partitions[[j]]
    sizes <- tabulate(partition, nbins = grid[j])
    data.table(scope = scope, orientation = orientation, rank = grid[j],
      minimum_size = min(sizes), maximum_size = max(sizes),
      mean_silhouette = if (uniqueN(partition) >= 2L) mean(silhouette(partition, distance)[, "sil_width"]) else NA_real_,
      minimum_size_qc_pass = min(sizes) >= local_minimum_genes)
  }))
  consensus <- lapply(grid, function(k) matrix(0, ncol(x), ncol(x)))
  agreement_rows <- list()
  unit_ids <- sort(unique(units))
  # clue::cl_agreement consumes RNG even for cRand. Draw every bootstrap first
  # with a fixed per-cell seed so changing the rank grid cannot change samples.
  bootstrap_seed <- seed + match(metadata$depmap_id[1], cell_map$depmap_id)
  set.seed(bootstrap_seed)
  bootstrap_units <- replicate(local_resamples,
    sample(unit_ids, length(unit_ids), replace = TRUE), simplify = FALSE)
  for (resample in seq_len(local_resamples)) {
    sampled <- bootstrap_units[[resample]]
    multiplicity <- tabulate(match(sampled, unit_ids), nbins = length(unit_ids))
    indices <- unlist(lapply(sampled, function(unit) which(units == unit)))
    current_x <- x[indices, , drop = FALSE]
    current_weights <- weights[indices]
    current_tree <- hclust(as.dist(sqrt(2 * (1 - condition_correlation[indices, indices]))), method = "ward.D2")
    for (j in seq_along(grid)) {
      current_partition <- genes_from_sample_clusters(current_x, cutree(current_tree, grid[j]), current_weights)
      consensus[[j]] <- consensus[[j]] + outer(current_partition, current_partition, "==")
      ari <- as.numeric(clue::cl_agreement(clue::as.cl_partition(partitions[[j]]),
        clue::as.cl_partition(current_partition), method = "cRand"))
      agreement_rows[[length(agreement_rows) + 1L]] <- data.table(scope = scope,
        orientation = orientation, rank = grid[j], resample = resample,
        adjusted_rand = ari, sampled_units = sum(multiplicity > 0),
        omitted_units = sum(multiplicity == 0))
    }
  }
  agreement <- rbindlist(agreement_rows)
  for (j in seq_along(grid)) {
    consensus[[j]] <- consensus[[j]] / local_resamples
    d <- as.dist(1 - consensus[[j]])
    ct <- hclust(d, method = "average")
    probability <- consensus[[j]][upper.tri(consensus[[j]])]
    rank_table[rank == grid[j], `:=`(
      consensus_cophenetic = if (sd(as.numeric(d)) > 0) cor(as.numeric(d), as.numeric(cophenetic(ct))) else NA_real_,
      pac = mean(probability > 0.1 & probability < 0.9),
      median_resample_ari = median(agreement[rank == grid[j], adjusted_rand]),
      q10_resample_ari = as.numeric(quantile(agreement[rank == grid[j], adjusted_rand], 0.1)))]
  }
  rank_table <- select_silhouette(rank_table)
  selected_rank <- if (any(rank_table$selected)) rank_table[selected == TRUE, rank] else NA_integer_
  selection_status <- rank_table$selection_status[1]
  # An unresolved scientific criterion is a reported abstention, not an input
  # error or permission to select the largest rank. Other scopes still complete.
  if (is.na(selected_rank)) {
    null_correlation <- cor(null)
    diagnostics <- data.table(scope = scope, orientation = orientation, conditions = nrow(x),
      discovery_units = length(unit_ids), genes = length(genes), selected_rank = NA_integer_,
      selection_status = selection_status, linkage = "sample/sample Ward.D2 on sqrt(2*(1-r))",
      raw_tree_cophenetic = cor(as.numeric(distance), as.numeric(cophenetic(tree))),
      correlation_with_control_structure = cor(correlation[upper.tri(correlation)],
        null_correlation[upper.tri(null_correlation)]), centering_adjusted_rand = NA_real_,
      mean_reconstruction_relative_error = NA_real_, null_profiles = nrow(null),
      experimental_replicate_agreement_available = FALSE)
    return(list(scope = scope, orientation = orientation, genes = genes, correlation = correlation,
      tree = tree, condition_tree = condition_tree, condition_correlation = condition_correlation,
      all_condition_partitions = condition_partitions, metadata = metadata, ranks = rank_table, resampling = agreement,
      diagnostics = diagnostics, programs = data.table(), weights = NULL,
      all_partitions = partitions, bootstrap_seed = bootstrap_seed, bootstrap_units = bootstrap_units,
      functional_labels_used = TRUE, functional_labels_used_for_gene_selection = TRUE,
      rank_selection_uses_outcome_objective = FALSE))
  }
  partition <- partitions[[as.character(selected_rank)]]
  program_ids <- sprintf("P%02d", seq_len(selected_rank))
  # Seed modules define profiles; every admitted gene can have positive weight
  # in several programs. Positive correlation preserves co-displacement direction.
  seeds <- sapply(seq_len(selected_rank), function(j) rowMeans(x[, partition == j, drop = FALSE]))
  centered_x <- sweep(x, 2L, colSums(x * weights) / sum(weights), "-")
  centered_seeds <- sweep(seeds, 2L, colSums(seeds * weights) / sum(weights), "-")
  cross <- crossprod(centered_x * sqrt(weights), centered_seeds * sqrt(weights))
  scales <- outer(sqrt(colSums(centered_x^2 * weights)), sqrt(colSums(centered_seeds^2 * weights)))
  membership <- cross / scales
  membership[membership < 0] <- 0
  membership <- membership^2
  membership <- sweep(membership, 2L, colSums(membership), "/")
  dimnames(membership) <- list(genes, program_ids)
  check(all(is.finite(membership)) && max(abs(colSums(membership) - 1)) < 1e-10,
    "Invalid overlapping local weights")
  scores <- x %*% membership
  null_scores <- null %*% membership
  score_sd <- matrixStats::colSds(scores)
  check(all(score_sd > 0), "Constant local program score")
  null_correlation <- cor(null)
  centered_correlation <- gene_correlation(sweep(x, 1L, metadata$global_displacement, "-"), weights)
  centered_partition <- genes_from_sample_clusters(sweep(x, 1L, metadata$global_displacement, "-"),
    condition_partitions[[as.character(selected_rank)]], weights)
  q <- qr(membership)
  check(q$rank == selected_rank, "Dependent local program weight columns")
  projected <- x %*% qr.Q(q)
  response_energy <- rowSums(x^2)
  reconstruction_error <- pmax(response_energy - rowSums(projected^2), 0) / response_energy
  # Control comparisons are descriptive; vehicle contrasts share reference plates.
  diagnostics <- data.table(scope = scope, orientation = orientation,
    conditions = nrow(x), discovery_units = length(unit_ids), genes = length(genes),
    selected_rank = selected_rank, selection_status = selection_status, linkage = "sample/sample Ward.D2 on sqrt(2*(1-r))",
    raw_tree_cophenetic = cor(as.numeric(distance), as.numeric(cophenetic(tree))),
    correlation_with_control_structure = cor(correlation[upper.tri(correlation)],
      null_correlation[upper.tri(null_correlation)]),
    centering_adjusted_rand = as.numeric(clue::cl_agreement(clue::as.cl_partition(partition),
      clue::as.cl_partition(centered_partition), method = "cRand")),
    mean_reconstruction_relative_error = mean(reconstruction_error),
    null_profiles = nrow(null), experimental_replicate_agreement_available = FALSE)
  program_table <- rbindlist(lapply(seq_len(selected_rank), function(j) {
    idx <- which(partition == j)
    pair <- correlation[idx, idx, drop = FALSE]
    cp <- null_correlation[idx, idx, drop = FALSE]
    null_rms <- sqrt(mean(null_scores[, j]^2))
    data.table(scope = scope, orientation = orientation, program = program_ids[j],
      seed_genes = length(idx), effective_weighted_genes = 1 / sum(membership[, j]^2),
      mean_gene_correlation = mean(pair[upper.tri(pair)]),
      mean_control_gene_correlation = mean(cp[upper.tri(cp)]),
      response_rms = sqrt(mean(scores[, j]^2)), control_rms = null_rms,
      response_control_rms_ratio = sqrt(mean(scores[, j]^2)) / null_rms,
      marker_genes = paste(genes[order(membership[, j], decreasing = TRUE)[1:8]], collapse = ", "))
  }))
  assignments <- data.table(gene = genes, seed_program = program_ids[partition],
    tree_order = match(seq_along(genes), tree$order))
  list(scope = scope, orientation = orientation, genes = genes, correlation = correlation,
    consensus = consensus[[match(selected_rank, grid)]], tree = tree, partition = partition,
    condition_tree = condition_tree, condition_correlation = condition_correlation,
    condition_partition = condition_partitions[[as.character(selected_rank)]],
    all_condition_partitions = condition_partitions,
    weights = membership, scores = scores, null_scores = null_scores, score_sd = score_sd,
    metadata = metadata, ranks = rank_table, resampling = agreement, diagnostics = diagnostics,
    programs = program_table, assignments = assignments, reconstruction_error = reconstruction_error,
    all_partitions = partitions, bootstrap_seed = bootstrap_seed, bootstrap_units = bootstrap_units,
    functional_labels_used = TRUE, functional_labels_used_for_gene_selection = TRUE,
    rank_selection_uses_outcome_objective = FALSE)
}

# Validate and reuse stage 04's controlled expression and matched-scale plate
# null. The shared cache is read only and has an explicit primary-input signature.
cell_map <- fread(local_cell_mapping_path)
check(!anyDuplicated(cell_map$depmap_id) && nrow(cell_map) == 7, "Expected seven unique recipients")
metadata <- fread(tahoe_metadata_path)
metadata[, source_drug := drug]
metadata[, drug := panel_name(drug)]
metadata <- metadata[depmap_id %chin% cell_map$depmap_id]
metadata[, recipient := cell_map$depmap_name[match(depmap_id, cell_map$depmap_id)]]
metadata[, focal := drug %chin% local_drugs]
check(!anyDuplicated(metadata$observation_id) && all(metadata$dose_unit == "uM"), "Invalid Tahoe profile keys/units")
check(metadata[focal == TRUE, .N] == 147, "Expected 147 focal cell/drug/dose profiles")
check(all(metadata$treated_total > 0 & metadata$reference_total > 0 &
  metadata$treated_correction_factor > 0), "Invalid recorded normalization")
audit <- fread(gene_audit_path)
gene_table <- fread(file.path(tahoe_dir, "genes.tsv"))[order(gene_index)]
genes <- sort(audit[protein_coding & median_TPM >= 5 & fraction_ge_1_TPM >= 0.8 &
  gene %chin% gene_table$feature, gene])
check(!anyDuplicated(genes) && length(genes) > 1000, "Invalid expressed coding gene universe")
gene_index <- match(genes, gene_table$feature)
input_files <- c(gene_audit_path, tahoe_metadata_path, local_cell_mapping_path,
  file.path(tahoe_dir, "genes.tsv"), file.path(tahoe_dir, "pdex_plate_reference.npy"),
  file.path(tahoe_dir, "pdex_plate_reference_metadata.tsv"), file.path(cache_dir, unique(metadata$shard_path)))
check(all(file.exists(input_files)), "Missing primary molecular input")
prep_signature <- digest::digest(list(genes = genes, metadata = metadata,
  input_sizes = file.info(input_files)$size, input_mtime = file.info(input_files)$mtime,
  normalizer = "paired_trimmed_median_log_ratio_cp10k_v1"), algo = "sha256")
prep_path <- file.path(local_output_dir, "data/controlled_molecular_inputs.qs")
check(file.exists(prep_path), "Missing shared controlled-molecular preparation; run stage 04 first")
prepared <- qread(prep_path)
check(identical(prepared$signature, prep_signature), "Incompatible shared molecular preparation cache")
check(identical(rownames(prepared$response), metadata$observation_id) &&
  identical(colnames(prepared$response), genes) && all(is.finite(prepared$response)),
  "Shared molecular matrix does not align with current inputs")
message("Reusing validated stage-04 controlled-molecular preparation")
fwrite(prepared$checks, file.path(data_dir, "normalization_checks.tsv"), sep = "\t")
fwrite(metadata, file.path(data_dir, "molecular_metadata.tsv.gz"), sep = "\t")

# Match retrospective outcomes before the supervised gene filter.
# ALMANAC monotherapy uses assay signal / matched vehicle, without time-zero
# subtraction. Growth-inhibition percentages are different quantities and unused.
almanac_path <- file.path(work_dir, "almanac_conditionals/raw/ComboDrugGrowth_Nov2017.csv")
drug_mapping_path <- file.path(work_dir, "almanac_conditionals/repaired/results/drug_mapping.tsv")
drug_map <- fread(drug_mapping_path)
drug_map[, drug := panel_name(tahoe)]
drug_map <- drug_map[primary == TRUE & drug %chin% metadata$drug]
check(!anyDuplicated(drug_map$nsc), "Invalid primary ALMANAC identity mapping")
fwrite(drug_map, file.path(data_dir, "almanac_drug_identity_mapping.tsv"), sep = "\t")
almanac <- fread(almanac_path, select = c("NSC1", "NSC2", "CONC1", "CONC2", "CONCUNIT1", "VALID",
  "CELLNAME", "SCREENER", "STUDY", "TESTVALUE", "CONTROLVALUE"))
almanac <- almanac[VALID == "Y" & is.na(NSC2) & NSC1 %in% drug_map$nsc &
  CELLNAME %chin% cell_map$CELLNAME & CONC1 > 0 & is.finite(TESTVALUE) &
  is.finite(CONTROLVALUE) & CONTROLVALUE > 0]
check(all(almanac$CONCUNIT1 == "M"), "Unexpected ALMANAC concentration units")
almanac[, `:=`(drug = drug_map$drug[match(NSC1, drug_map$nsc)],
  depmap_id = cell_map$depmap_id[match(CELLNAME, cell_map$CELLNAME)], dose = CONC1 * 1e6,
  abundance_pct = TESTVALUE / CONTROLVALUE * 100)]
almanac_summary <- almanac[, .(abundance_pct = mean(abundance_pct),
  abundance_sd = sd(abundance_pct), n_records = .N, n_studies = uniqueN(STUDY),
  n_laboratories = uniqueN(SCREENER)), by = .(depmap_id, drug, dose)]
fwrite(almanac_summary, file.path(data_dir, "almanac_monotherapy_observed.tsv"), sep = "\t")

# PRISM supplies measured log2 abundance at a different dose grid. The release
# recommends MTS010 when available. Use that screen when it brackets a query;
# otherwise use an explicitly recorded bracketing screen. Never extrapolate.
prism_dir <- file.path(cache_dir, "prism")
prism_info_path <- file.path(prism_dir, "secondary-screen-replicate-collapsed-treatment-info.csv")
prism_matrix_path <- file.path(prism_dir, "secondary-screen-replicate-collapsed-logfold-change.csv")
prism_info <- fread(prism_info_path)
prism_info[, drug := panel_name(name)]
# The existing reviewed crosswalk supplies name aliases, without using outcome
# values. Direct name matches cover additional Tahoe drugs beyond that crosswalk.
function_mapping_path <- file.path(work_dir, "tahoe_prism_v2/artifacts/01_inputs/function_outcomes.tsv.gz")
function_mapping <- unique(fread(function_mapping_path, select = c("tahoe_drug", "exact_name_key")))
function_mapping[, `:=`(drug = panel_name(tahoe_drug), name_key = panel_name(exact_name_key))]
function_mapping <- unique(function_mapping[, .(drug, name_key)])
check(!anyDuplicated(function_mapping$name_key), "Ambiguous reviewed PRISM name crosswalk")
alias <- match(prism_info$drug, function_mapping$name_key)
prism_info[!is.na(alias), drug := function_mapping$drug[alias[!is.na(alias)]]]
prism_info <- prism_info[drug %chin% metadata$drug & dose > 0]
fwrite(prism_info[, .(name, drug, broad_id, screen_id, dose)],
  file.path(data_dir, "prism_drug_identity_mapping.tsv"), sep = "\t")
prism_header <- names(fread(prism_matrix_path, nrows = 0))
prism_info <- prism_info[column_name %chin% prism_header]
prism_matrix <- fread(prism_matrix_path, select = c(prism_header[1], prism_info$column_name))
setnames(prism_matrix, 1L, "depmap_id")
prism_matrix <- prism_matrix[depmap_id %chin% cell_map$depmap_id]
measurement_columns <- setdiff(names(prism_matrix), "depmap_id")
prism_matrix[, (measurement_columns) := lapply(.SD, as.numeric), .SDcols = measurement_columns]
prism_long <- melt(prism_matrix, id.vars = "depmap_id", variable.name = "column_name",
  value.name = "log2_fold_change")
prism_long <- merge(prism_long, prism_info[, .(column_name, drug, dose, screen_id, broad_id)],
  by = "column_name", all.x = TRUE, sort = FALSE)
prism_long <- prism_long[is.finite(log2_fold_change)]
prism_points <- prism_long[, .(log2_fold_change = mean(log2_fold_change), n_records = .N,
  broad_ids = paste(sort(unique(broad_id)), collapse = ",")), by = .(depmap_id, drug, dose, screen_id)]
fwrite(prism_points, file.path(data_dir, "prism_observed_dose_points.tsv"), sep = "\t")

sensitivity <- copy(metadata[, .(observation_id, depmap_id, recipient, drug, dose, focal)])
sensitivity[, `:=`(almanac_abundance_pct = NA_real_, almanac_status = "missing exact dose",
  almanac_n_records = NA_integer_, almanac_n_studies = NA_integer_,
  prism_abundance_pct = NA_real_, prism_status = "missing bracket",
  prism_lower_dose = NA_real_, prism_upper_dose = NA_real_, prism_lower_log2 = NA_real_,
  prism_upper_log2 = NA_real_, prism_screen = NA_character_)]
for (i in seq_len(nrow(sensitivity))) {
  if (i %% 1000 == 0) message("Aligned sensitivity for ", i, "/", nrow(sensitivity), " conditions")
  id <- sensitivity$depmap_id[i]
  drug_name <- sensitivity$drug[i]
  query_dose <- sensitivity$dose[i]
  a <- almanac_summary[depmap_id == id & drug == drug_name & abs(log(dose / query_dose)) < 1e-6]
  check(nrow(a) <= 1, "Ambiguous exact ALMANAC monotherapy query")
  if (nrow(a) == 1) {
    sensitivity[i, `:=`(almanac_abundance_pct = a$abundance_pct, almanac_status = "exact",
      almanac_n_records = a$n_records, almanac_n_studies = a$n_studies)]
  }
  p <- prism_points[depmap_id == id & drug == drug_name]
  screens <- sort(unique(p$screen_id))
  screens <- c(intersect("MTS010", screens), setdiff(screens, "MTS010"))
  for (screen in screens) {
    pp <- p[screen_id == screen][order(dose)]
    if (query_dose < min(pp$dose) || query_dose > max(pp$dose)) next
    exact <- which(abs(log(pp$dose / query_dose)) < 1e-6)
    lo <- if (length(exact)) exact[1] else max(which(pp$dose < query_dose))
    hi <- if (length(exact)) exact[1] else min(which(pp$dose > query_dose))
    log_abundance <- if (lo == hi) pp$log2_fold_change[lo] else approx(log(pp$dose[c(lo, hi)]),
      pp$log2_fold_change[c(lo, hi)], xout = log(query_dose))$y
    sensitivity[i, `:=`(prism_abundance_pct = 100 * 2^log_abundance,
      prism_status = if (lo == hi) "exact" else "interpolated",
      prism_lower_dose = pp$dose[lo], prism_upper_dose = pp$dose[hi],
      prism_lower_log2 = pp$log2_fold_change[lo], prism_upper_log2 = pp$log2_fold_change[hi],
      prism_screen = screen)]
    break
  }
}
check(!anyDuplicated(sensitivity$observation_id), "Duplicate sensitivity annotation")
fwrite(sensitivity, file.path(data_dir, "sensitivity_alignment.tsv"), sep = "\t")
fwrite(sensitivity[, .(conditions = .N), by = .(almanac_status, prism_status)],
  file.path(data_dir, "sensitivity_coverage.tsv"), sep = "\t")
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm = TRUE) / 25) * 25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm = TRUE) / 25) * 25)
message("Sensitivity bar maxima (% vehicle): ALMANAC ", almanac_bar_max, "; PRISM ", prism_bar_max)

# Select genes by two separate assay-specific outcome correlations. Each matched
# drug/dose observation has equal weight, with all doses pooled. Retain either-sign
# genes with nominal p <=0.05 in either assay, capped at 500 by
# the best qualifying assay p-value. BH adjusted values are diagnostics only. Molecular
# variation/control metrics remain diagnostics and never gate admission.
line_genes <- list()
gene_selection <- list()
outcome_test_tables <- list()
outcome_test_summary <- list()
gene_selection_signatures <- list()
condition_checks <- list()
for (id in cell_map$depmap_id) {
  rows <- which(metadata$depmap_id == id)
  x <- prepared$response[rows, , drop = FALSE]
  null <- prepared$null[[id]]
  baseline <- prepared$baseline[[id]]
  noise_variance <- matrixStats::colVars(null)
  abundance_bin <- pmin(10L, ceiling(rank(baseline, ties.method = "first") / length(baseline) * 10))
  noise_floor <- ave(noise_variance, abundance_bin, FUN = median)
  noise <- pmax(noise_variance, noise_floor)
  check(all(noise > 0), "Nonpositive abundance-stratified control variance")
  null_threshold <- as.numeric(matrixStats::colQuantiles(abs(null), probs = 0.95))
  exceeds <- sweep(abs(x), 2L, null_threshold, ">")
  supported_drugs <- colSums(rowsum(exceeds * 1L, metadata$drug[rows], reorder = FALSE) > 0)
  variance <- matrixStats::colVars(x)
  focus <- which(metadata$focal[rows])
  focal_variance <- matrixStats::colVars(x[focus, , drop = FALSE])
  table <- data.table(depmap_id = id, recipient = metadata$recipient[rows[1]], gene = genes,
    response_variance = variance, control_variance = noise,
    excess_variance = pmax(variance - noise, 0), focal_excess_variance = pmax(focal_variance - noise, 0),
    supported_drugs = supported_drugs, focal_supported_conditions = colSums(exceeds[focus, , drop = FALSE]),
    abundance_bin = abundance_bin)
  # Keep original-scale variability visible without imposing the old SD/noise gates.
  table[, `:=`(response_sd = sqrt(response_variance), variance_ratio = response_variance / control_variance)]
  matched_outcomes <- sensitivity[match(metadata$observation_id[rows], observation_id)]
  check(identical(matched_outcomes$observation_id, metadata$observation_id[rows]),
    "Outcome rows do not align with gene displacement")
  tests <- rbindlist(lapply(c("PRISM", "ALMANAC"), function(assay_name) {
    outcome_column <- if (assay_name == "PRISM") "prism_abundance_pct" else "almanac_abundance_pct"
    result <- gene_outcome_tests(x, matched_outcomes[[outcome_column]], metadata$drug[rows])
    result[, `:=`(depmap_id = id, recipient = metadata$recipient[rows[1]], assay = assay_name)]
    result
  }))
  tests[, q_value := p.adjust(p_value, method = "BH", n = sum(!is.na(p_value))), by = assay]
  tests[, `:=`(significant = !is.na(p_value) & p_value <= outcome_p_threshold &
      abs(correlation) >= outcome_minimum_abs_correlation,
    outcome_scale = "abundance (% vehicle)", dose_handling = "all matched doses pooled; equal condition weight",
    multiplicity_family = "diagnostic BH across estimable genes within this cell and assay",
    selection_threshold = sprintf("nominal two-sided p <=%g and |r| >=%g in the same assay; BH diagnostic only",
      outcome_p_threshold, outcome_minimum_abs_correlation))]
  outcome_test_tables[[id]] <- tests
  outcome_test_summary[[id]] <- tests[, .(candidate_genes = .N, tested_genes = sum(test_status == "tested"),
    significant_genes = sum(significant), matched_conditions = matched_conditions[1],
    matched_drugs = matched_drugs[1], non_estimable_genes = sum(test_status != "tested")),
    by = .(depmap_id, recipient, assay)]
  prism_tests <- tests[assay == "PRISM", .(gene, prism_correlation = correlation,
    prism_p_value = p_value, prism_q_value = q_value, prism_outcome_eligible = significant)]
  almanac_tests <- tests[assay == "ALMANAC", .(gene, almanac_correlation = correlation,
    almanac_p_value = p_value, almanac_q_value = q_value, almanac_outcome_eligible = significant)]
  table <- merge(table, prism_tests, by = "gene", all.x = TRUE, sort = FALSE)
  table <- merge(table, almanac_tests, by = "gene", all.x = TRUE, sort = FALSE)
  table[, best_outcome_p_value := pmin(prism_p_value, almanac_p_value, na.rm = TRUE)]
  table[is.na(prism_p_value) & is.na(almanac_p_value), best_outcome_p_value := NA_real_]
  table[, best_outcome_q_value := pmin(prism_q_value, almanac_q_value, na.rm = TRUE)]
  table[is.na(prism_q_value) & is.na(almanac_q_value), best_outcome_q_value := NA_real_]
  table[, strongest_outcome_correlation := pmax(abs(prism_correlation), abs(almanac_correlation), na.rm = TRUE)]
  table[is.na(prism_correlation) & is.na(almanac_correlation), strongest_outcome_correlation := NA_real_]
  table[, outcome_eligible := prism_outcome_eligible | almanac_outcome_eligible]
  # Cap ordering uses the qualifying assays, with absolute correlation breaking ties.
  table[, best_qualifying_p_value := pmin(fifelse(prism_outcome_eligible, prism_p_value, NA_real_),
    fifelse(almanac_outcome_eligible, almanac_p_value, NA_real_), na.rm = TRUE)]
  table[, strongest_qualifying_correlation := pmax(fifelse(prism_outcome_eligible, abs(prism_correlation), NA_real_),
    fifelse(almanac_outcome_eligible, abs(almanac_correlation), NA_real_), na.rm = TRUE)]
  table[outcome_eligible == FALSE, `:=`(best_qualifying_p_value = NA_real_, strongest_qualifying_correlation = NA_real_)]
  chosen <- sort(table[outcome_eligible == TRUE][order(best_qualifying_p_value, -strongest_qualifying_correlation, gene),
    head(gene, local_gene_limit)])
  check(length(chosen) >= 3L, paste("Fewer than three outcome-significant genes for", id,
    "at nominal p", outcome_p_threshold, "and |r| >=", outcome_minimum_abs_correlation,
    "-- no variance-selected substitute"))
  table[, selected := gene %chin% chosen]
  table[, selection_role := "outcome-significant in PRISM or ALMANAC; assays tested separately"]
  table[selected == FALSE, selection_role := "not selected"]
  table[, exclusion_reason := fcase(is.na(best_outcome_p_value), "no estimable outcome correlation",
    !outcome_eligible, "fails nominal p or absolute-correlation gate in both assays",
    selected, "selected", default = "significant outside gene cap")]
  line_genes[[id]] <- chosen
  gene_selection[[id]] <- table
  gene_selection_signatures[[id]] <- digest::digest(list(rule = "separate Pearson assays; either-sign nominal significance",
    p_threshold = outcome_p_threshold, minimum_abs_correlation = outcome_minimum_abs_correlation,
    cap = local_gene_limit, observations = metadata$observation_id[rows],
    outcomes = matched_outcomes[, .(observation_id, prism_abundance_pct, almanac_abundance_pct)],
    tests = tests[, .(gene, assay, correlation, p_value, q_value)], selected_genes = chosen), algo = "sha256")
  message(metadata$recipient[rows[1]], ": ", sum(table$outcome_eligible), " outcome-significant genes; retained ",
    length(chosen), "; PRISM/ALMANAC matched conditions ",
    paste(outcome_test_summary[[id]]$matched_conditions, collapse = "/"))
  null_energy <- mean(rowMeans(null^2))
  energy <- rowMeans(x^2)
  support <- energy / (energy + null_energy)
  metadata[rows, `:=`(condition_support = support, global_displacement = rowMeans(x),
    raw_response_rms = sqrt(energy), control_rms = sqrt(null_energy),
    response_control_rms_ratio = sqrt(energy / null_energy))]
  condition_checks[[id]] <- metadata[rows, .(observation_id, depmap_id, recipient, drug, dose, focal,
    n_plates, condition_support, raw_response_rms, control_rms, response_control_rms_ratio,
    global_displacement, contributing_plates)]
}
fwrite(rbindlist(gene_selection), file.path(data_dir, "gene_selection.tsv.gz"), sep = "\t")
selection_audit <- rbindlist(gene_selection)
fwrite(selection_audit[, .(genes = .N), by = .(recipient, exclusion_reason)],
  file.path(data_dir, "gene_filter_counts.tsv"), sep = "\t")
fwrite(rbindlist(outcome_test_tables), file.path(data_dir, "gene_outcome_correlations.tsv.gz"), sep = "\t")
fwrite(rbindlist(outcome_test_summary), file.path(data_dir, "gene_outcome_test_summary.tsv"), sep = "\t")
# Use a finite display floor for log10(0); every untransformed p-value is exported.
selection_audit[, displayed_minus_log10_p := -log10(pmax(best_outcome_p_value, .Machine$double.xmin))]
p <- ggplot(selection_audit, aes(displayed_minus_log10_p, fill = selected)) +
  geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.6, na.rm = TRUE) +
  geom_vline(xintercept = -log10(outcome_p_threshold), lty = 2) +
  facet_wrap(~recipient, ncol = 4) + scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#27818A")) +
  labs(x = "-log10(best outcome nominal p-value)", y = "Genes", fill = "Selected",
    title = "Outcome-informed gene selection: PRISM and ALMANAC tested separately",
    subtitle = sprintf("All doses pooled equally; nominal p <=%g AND |r| >=%g in the same assay; 500-gene cap",
      outcome_p_threshold, outcome_minimum_abs_correlation))
save_plot(p, "gene_outcome_filter", 11.5, 5.7)
fwrite(rbindlist(condition_checks), file.path(data_dir, "condition_signal_vs_control.tsv"), sep = "\t")

# Discover within cell lines across drugs and doses. These are candidate
# banks; biological reproducibility is unavailable from aggregated PDEX profiles.
discoveries <- list()
for (id in cell_map$depmap_id) {
  rows <- which(metadata$depmap_id == id)
  scope <- paste0("cell_", gsub("[^A-Za-z0-9]+", "_", metadata$recipient[rows[1]]))
  scope_model_dir <- file.path(output_dir, scope, "models")
  dir.create(scope_model_dir, recursive = TRUE, showWarnings = FALSE)
  bank_path <- file.path(scope_model_dir, "program_bank.qs")
  if (render_only) {
    check(file.exists(bank_path), paste("No saved molecular analysis to render for", scope))
    bank <- qread(bank_path)
    check(identical(bank$genes, line_genes[[id]]) &&
      identical(bank$metadata$observation_id, metadata$observation_id[rows]) &&
      identical(bank$ranks$rank, as.integer(local_rank_grid[local_rank_grid <=
        min(length(bank$genes) - 1L, length(rows) - 1L)])) &&
      nrow(bank$resampling) == nrow(bank$ranks) * local_resamples &&
      identical(bank$functional_labels_used, TRUE) &&
      identical(bank$functional_labels_used_for_gene_selection, TRUE) &&
      identical(bank$rank_selection_uses_outcome_objective, FALSE) &&
      identical(bank$gene_selection_signature, gene_selection_signatures[[id]]),
      paste("Saved bank settings, outcomes or gene-admission rule disagree for", scope))
    # Check the actual saved molecular structure against current primary data,
    # not just filenames. Render-only never changes weights or selected ranks.
    current_correlation <- gene_correlation(prepared$response[rows, bank$genes, drop = FALSE],
      group_weights(metadata$drug[rows]) * metadata$condition_support[rows])
    check(max(abs(current_correlation - bank$correlation)) < 1e-10,
      paste("Saved molecular correlation differs from current inputs for", scope))
    check(max(abs(cor(t(prepared$response[rows, bank$genes, drop = FALSE])) - bank$condition_correlation)) < 1e-10,
      paste("Saved sample correlation differs from current inputs for", scope))
    current_selection <- select_silhouette(bank$ranks)
    check(identical(bank$ranks$selected, current_selection$selected) &&
      identical(bank$ranks$selection_status, current_selection$selection_status) &&
      identical(bank$ranks$selection_rule, current_selection$selection_rule) &&
      identical(bank$ranks$minimum_size_qc_pass, current_selection$minimum_size_qc_pass) &&
      identical(bank$ranks$silhouette_candidate, current_selection$silhouette_candidate),
      paste("Saved cut does not follow the last-blue-dot mean gene silhouette rule for", scope))
    check(identical(bank$bootstrap_seed, seed + match(id, cell_map$depmap_id)),
      paste("Saved bootstrap seed differs for", scope))
    set.seed(bank$bootstrap_seed)
    expected_units <- replicate(local_resamples,
      sample(sort(unique(bank$metadata$drug)), uniqueN(bank$metadata$drug), replace = TRUE), simplify = FALSE)
    check(identical(bank$bootstrap_units, expected_units), paste("Saved resampling units differ for", scope))
    message("Rendering saved molecular analysis: ", scope)
  } else {
    message("Discovering ", scope, " over ", length(rows), " conditions")
    bank <- discover_local(prepared$response[rows, , drop = FALSE], prepared$null[[id]],
      metadata[rows], line_genes[[id]], "drug", scope, "cell")
    bank$gene_selection_signature <- gene_selection_signatures[[id]]
    bank$gene_selection_rule <- sprintf(
      "Separate Pearson assays; equal-weight pooled doses; either sign; p <=%g AND |r| >=%g in the same assay; cap %d",
      outcome_p_threshold, outcome_minimum_abs_correlation, local_gene_limit)
    qsave(bank, bank_path, preset = "fast")
  }
  discoveries[[scope]] <- bank
}
fwrite(rbindlist(lapply(discoveries, `[[`, "ranks")), file.path(data_dir, "rank_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(lapply(discoveries, `[[`, "resampling")), file.path(data_dir, "resampling_diagnostics.tsv"), sep = "\t")
fwrite(rbindlist(lapply(discoveries, `[[`, "diagnostics")), file.path(data_dir, "discovery_diagnostics.tsv"), sep = "\t")
program_summary <- rbindlist(lapply(discoveries, `[[`, "programs"))
if (!ncol(program_summary)) program_summary <- data.table(scope = character(), orientation = character(),
  program = character(), seed_genes = integer(), effective_weighted_genes = numeric(),
  mean_gene_correlation = numeric(), mean_control_gene_correlation = numeric(), response_rms = numeric(),
  control_rms = numeric(), response_control_rms_ratio = numeric(), marker_genes = character())
fwrite(program_summary, file.path(data_dir, "program_diagnostics.tsv"), sep = "\t")
freeze <- data.table(scope = names(discoveries), functional_labels_used = TRUE,
  functional_labels_used_for_gene_selection = TRUE, rank_selection_uses_outcome_objective = FALSE,
  bank_path = file.path(output_dir, names(discoveries), "models/program_bank.qs"))
freeze[, sha256 := vapply(bank_path, digest::digest, character(1), file = TRUE, algo = "sha256")]
fwrite(freeze, file.path(data_dir, "molecular_bank_freeze.tsv"), sep = "\t")
message("All seven outcome-informed banks saved; continuing descriptive figure rendering")
# Each scope owns a separate diagnostics PDF and PNG, with all cuts visible.
for (bank in discoveries) {
  dir.create(file.path(output_dir, bank$scope, "figures"), recursive = TRUE, showWarnings = FALSE)
  r <- bank$ranks
  p1 <- ggplot(r, aes(rank, mean_silhouette)) + geom_line() + geom_point(aes(color = minimum_size_qc_pass)) +
    scale_color_manual(values = c(`FALSE` = "#F8766D", `TRUE` = "#00BFC4")) +
    geom_vline(xintercept = local_minimum_rank - 0.5, lty = 3, color = "grey60") +
    geom_vline(data = r[selected == TRUE], aes(xintercept = rank), lty = 2, color = "#B24C35") +
    geom_point(data = r[selected == TRUE], shape = 1, size = 4, color = "#B24C35") +
    labs(x = "Candidate programs", y = "Mean gene silhouette",
      color = paste0("Size QC >=", local_minimum_genes),
      subtitle = "Select last blue dot with positive silhouette; ranks 2--5 barred")
  silhouette_plot <- p1 + labs(title = gsub("_", " ", bank$scope),
    caption = paste0("Red ring/dashed line: selected rank. Blue: every seed >=", local_minimum_genes,
      " genes.\nSelect the largest allowed size-qualified rank with mean gene silhouette >0."))
  ggsave(file.path(output_dir, bank$scope, "figures/silhouette.pdf"), silhouette_plot,
    width = 8, height = 4.5, device = cairo_pdf)
  ggsave(file.path(output_dir, bank$scope, "figures/silhouette.png"), silhouette_plot,
    width = 8, height = 4.5, dpi = 160, bg = "white")
  p2 <- ggplot(r, aes(rank, consensus_cophenetic)) + geom_line() + geom_point() +
    geom_rect(aes(xmin = 1.5, xmax = 5.5, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey80", alpha = 0.2, data = data.frame(unused = 1)) +
    geom_vline(data = r[selected == TRUE], aes(xintercept = rank), lty = 2, color = "#B24C35") +
    geom_point(data = r[selected == TRUE], size = 3, color = "#B24C35") +
    labs(x = "Candidate programs", y = "Consensus cophenetic correlation",
      subtitle = "Diagnostic only; rank selected on mean gene silhouette")
  cophenetic_plot <- p2 + labs(title = gsub("_", " ", bank$scope),
    caption = "24 whole-drug resamples; red: rank selected from the mean gene silhouette curve.")
  ggsave(file.path(output_dir, bank$scope, "figures/cophenetic.pdf"), cophenetic_plot,
    width = 8, height = 4.5, device = cairo_pdf)
  ggsave(file.path(output_dir, bank$scope, "figures/cophenetic.png"), cophenetic_plot,
    width = 8, height = 4.5, dpi = 160, bg = "white")
  p3 <- ggplot(r, aes(rank, median_resample_ari)) + geom_line() + geom_point() +
    geom_line(aes(y = q10_resample_ari), lty = 2, color = "grey50") +
    labs(x = "Candidate programs", y = "Resampling adjusted Rand", subtitle = "Solid: median; dashed: 10th percentile")
  p4 <- ggplot(r, aes(rank, pac)) + geom_line() + geom_point() +
    labs(x = "Candidate programs", y = "Ambiguous consensus fraction (PAC)", subtitle = "Gene-pair consensus in (0.1, 0.9)")
  composite <- (p1 + p2) / (p3 + p4) + plot_annotation(title = gsub("_", " ", bank$scope),
    subtitle = sprintf("%d conditions; %d drug/recipient resampling units; raw-tree cophenetic %.3f; selected k=%d",
      bank$diagnostics$conditions, bank$diagnostics$discovery_units,
      bank$diagnostics$raw_tree_cophenetic, bank$diagnostics$selected_rank),
    caption = paste("Gene admission uses retrospective outcomes; rank uses silhouette on admitted genes.",
      "\n24 whole-drug resamples with the gene filter fixed; collection stability, not experimental replication."))
  ggsave(file.path(output_dir, bank$scope, "figures/diagnostics.pdf"), composite,
    width = 9.8, height = 7, device = cairo_pdf)
  ggsave(file.path(output_dir, bank$scope, "figures/diagnostics.png"), composite,
    width = 9.8, height = 7, dpi = 160, bg = "white")
  scope_data_dir <- file.path(output_dir, bank$scope, "data")
  dir.create(scope_data_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(bank$ranks, file.path(scope_data_dir, "rank_diagnostics.tsv"), sep = "\t")
  fwrite(bank$resampling, file.path(scope_data_dir, "resampling_diagnostics.tsv"), sep = "\t")
  fwrite(rbindlist(lapply(seq_along(bank$bootstrap_units), function(j) {
    data.table(resample = j, draw = seq_along(bank$bootstrap_units[[j]]), drug = bank$bootstrap_units[[j]],
      bootstrap_seed = bank$bootstrap_seed)
  })), file.path(scope_data_dir, "resampling_units.tsv.gz"), sep = "\t")
  fwrite(bank$diagnostics, file.path(scope_data_dir, "discovery_diagnostics.tsv"), sep = "\t")
  cell_id <- unique(bank$metadata$depmap_id)
  check(length(cell_id) == 1L, "A cell-line bank contains multiple recipients")
  cell_selection <- gene_selection[[cell_id]]
  fwrite(cell_selection, file.path(scope_data_dir, "gene_selection.tsv.gz"), sep = "\t")
  fwrite(outcome_test_tables[[cell_id]], file.path(scope_data_dir, "gene_outcome_correlations.tsv.gz"), sep = "\t")
  fwrite(outcome_test_summary[[cell_id]], file.path(scope_data_dir, "gene_outcome_test_summary.tsv"), sep = "\t")
  fwrite(cell_selection[, .(genes = .N), by = exclusion_reason],
    file.path(scope_data_dir, "gene_filter_counts.tsv"), sep = "\t")
  if (is.null(bank$weights)) writeLines(c(bank$diagnostics$selection_status,
    "No dictionary was selected. See figures/diagnostics.pdf and data/rank_diagnostics.tsv.",
    "figures/displacement_unresolved.pdf shows full raw displacements; figures/focal_panel.pdf is the readable focus.",
    paste("Original responses remain in the shared preparation:", prep_path)),
    file.path(output_dir, bank$scope, "README.txt"))
  else writeLines(c(paste("Selected programs:", bank$diagnostics$selected_rank),
    sprintf("Genes: two-sided Pearson with PRISM or ALMANAC; p <=%g AND |r| >=%g in the same assay; both signs; cap %d.",
      outcome_p_threshold, outcome_minimum_abs_correlation, local_gene_limit),
    "Assays tested separately; matched doses pooled with equal condition weight; BH values diagnostic only.",
    "Rank selection: last blue dot with positive mean gene silhouette; ranks 2--5 barred; every seed >=12 genes.",
    "figures/silhouette.pdf shows the selection; cophenetic is a diagnostic only.",
    "figures/displacement_overview.pdf shows all raw displacements grouped by the selected programs.",
    "figures/focal_panel.pdf shows the seven-drug focus; displacement_pages.pdf covers all conditions.",
    "Within program: PRISM first, ALMANAC-only second, no data last; lower % vehicle first.",
    "Outcome-informed exploratory dictionary; the displayed outcomes were used for gene selection."),
    file.path(output_dir, bank$scope, "README.txt"))
}


# A measured zero is sensitive; NA is unavailable. The separate value strip
# retains that distinction even when the numerical bar is absent or tiny.
sensitivity_palette <- c(`<=50% vehicle` = "#762A83", `>50% vehicle` = "#1B9E77", `No value (NA)` = "#D9D9D9")
sensitivity_category <- function(values) {
  fifelse(is.na(values), "No value (NA)", fifelse(values <= 50, "<=50% vehicle", ">50% vehicle"))
}
sensitivity_track <- function(values) {
  ComplexHeatmap::anno_simple(sensitivity_category(values), col = sensitivity_palette,
    pch = fifelse(is.na(values), 4L, NA_integer_), pt_gp = grid::gpar(col = "#555555"),
    pt_size = grid::unit(min(1.5, 220 / length(values)), "mm"), height = grid::unit(4, "mm"), border = TRUE)
}
sensitivity_legend <- ComplexHeatmap::Legend(title = "Assay abundance", labels = names(sensitivity_palette),
  legend_gp = grid::gpar(fill = unname(sensitivity_palette)),
  title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))

# The added plot uses the original signed sample correlation on both axes,
# aligned to the existing displacement-map column order. Used for both statuses.
render_sample_correlation <- function(bank, identifiers, annotation, palette, destination, source_dir) {
  correlation <- bank$condition_correlation[identifiers, identifiers, drop = FALSE]
  check(identical(rownames(correlation), rownames(annotation)) &&
    identical(colnames(correlation), rownames(annotation)), "Sample correlation annotations are misaligned")
  programs <- if ("Program" %in% names(annotation)) annotation$Program else NULL
  hm <- draw_heatmap(correlation, cluster_rows = FALSE, cluster_columns = FALSE,
    row_anno_df = annotation, row_anno_cols = names(annotation),
    col_anno_df = annotation, col_anno_cols = names(annotation), palette_overrides = palette,
    ann_top_args = list(show_legend = FALSE), ann_row_args = list(show_annotation_name = FALSE), scale_rows = FALSE,
    row_split = programs, column_split = programs,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_row_names = FALSE, show_column_names = FALSE,
    row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
    value_range = c(-1, 0, 1), value_colors = c("#2166AC", "white", "#B2182B"), hm_name = "Pearson r",
    title = paste0(gsub("_", " ", bank$scope), " | sample/sample correlation\n",
      "Ward.D2 on Pearson distance across selected genes; same column order as displacement maps"))
  save_heatmap(hm, "condition_correlation", 8.5, 8, destination)
  fwrite(data.table(observation_id = identifiers, correlation),
    file.path(source_dir, "condition_correlation.tsv.gz"), sep = "\t")
}

# Render each cell-line bank. Every heatmap uses the shared
# ComplexHeatmap wrapper. Off-block values and signed original effects remain.
# Draw full-compendium overviews and readable pages from one frozen dictionary.
# Reused annotation/heatmap construction stays here beside its scientific owner.
render_local <- function(bank) {
  scope <- bank$scope
  scope_data_dir <- file.path(output_dir, scope, "data")
  scope_figure_dir <- file.path(output_dir, scope, "figures")
  for (path in c(scope_data_dir, scope_figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  program_ids <- colnames(bank$weights)
  k <- length(program_ids)
  md <- copy(bank$metadata)
  annotation <- sensitivity[match(md$observation_id, observation_id)]
  check(identical(md$observation_id, annotation$observation_id), "Sensitivity annotation ordering mismatch")
  # Within each frozen program block, show PRISM annotations first, then
  # ALMANAC-only annotations, then unavailable samples. Lower % vehicle means
  # greater sensitivity. ALMANAC breaks PRISM ties; molecular order breaks other ties.
  annotation[, sensitivity_sort_source := fcase(!is.na(prism_abundance_pct), "PRISM",
    !is.na(almanac_abundance_pct), "ALMANAC", default = "none")]
  annotation[, sensitivity_sort_abundance_pct := fcase(sensitivity_sort_source == "PRISM", prism_abundance_pct,
    sensitivity_sort_source == "ALMANAC", almanac_abundance_pct, default = NA_real_)]
  x <- prepared$response[match(md$observation_id, rownames(prepared$response)), bank$genes, drop = FALSE]
  scores <- x %*% bank$weights
  standardized <- sweep(scores, 2L, bank$score_sd, "/")
  dominant <- bank$condition_partition
  row_order <- order(bank$partition, match(seq_along(bank$genes), bank$tree$order))
  column_order <- unlist(lapply(seq_len(k), function(j) {
    rows <- which(dominant == j)
    rows <- rows[order(match(rows, bank$condition_tree$order))]
    ranks <- match(annotation$sensitivity_sort_source[rows], c("PRISM", "ALMANAC", "none"))
    rows <- rows[order(ranks, annotation$sensitivity_sort_abundance_pct[rows],
      annotation$almanac_abundance_pct[rows], na.last = TRUE)]
    rows
  }))
  check(setequal(column_order, seq_len(nrow(md))), "Incomplete condition ordering")
  md <- md[column_order]
  annotation <- annotation[column_order]
  scores <- scores[column_order, , drop = FALSE]
  dominant <- dominant[column_order]
  x <- x[column_order, , drop = FALSE]
  stopifnot(identical(rownames(x), md$observation_id))
  row_meta <- data.frame(Program = factor(program_ids[bank$partition[row_order]], levels = program_ids),
    row.names = bank$genes[row_order])
  col_meta <- data.frame(Program = factor(program_ids[dominant], levels = program_ids),
    Dose = sprintf("%g", md$dose), check.names = FALSE, row.names = md$observation_id)
  palette <- list(Program = setNames(hcl.colors(k, "Dark 3"), program_ids))
  render_sample_correlation(bank, md$observation_id, col_meta, palette, scope_figure_dir, scope_data_dir)
  selected_rows <- unlist(lapply(seq_len(k), function(j) {
    idx <- which(bank$partition == j)
    head(idx[order(bank$weights[idx, j], decreasing = TRUE)], local_display_genes_per_program)
  }))
  selected_rows <- row_order[row_order %in% selected_rows]
  representative <- t(x[, selected_rows, drop = FALSE])
  full_matrix <- t(x[, row_order, drop = FALSE])
  limit <- max(0.25, as.numeric(quantile(abs(full_matrix), 0.98)))
  # Transliterate display names for PDF text measurement; source names and all
  # identity joins retain the original Unicode drug strings.
  display_drugs <- gsub("\u03b1", "alpha", md$drug, fixed = TRUE)
  display_drugs <- gsub("\u03b2", "beta", display_drugs, fixed = TRUE)
  display_drugs <- gsub("\u03b3", "gamma", display_drugs, fixed = TRUE)
  labels <- paste(substr(display_drugs, 1, 38), sprintf("%g uM", md$dose))

  # Bar scale comes from actual observations and is shared across scopes within
  # each assay. Preserve abundance >100%; no clipping of large positive values.
  # NA graphics placeholders have an explicit match strip; exports retain NA.
  top_args <- function(indices) {
    av <- annotation$almanac_abundance_pct[indices]
    pv <- annotation$prism_abundance_pct[indices]
    av[is.na(av)] <- 0
    pv[is.na(pv)] <- 0
    list(`ALMANAC % vehicle` = ComplexHeatmap::anno_barplot(av,
      gp = grid::gpar(fill = "#B24C35", col = NA), ylim = c(0, almanac_bar_max),
      height = grid::unit(11, "mm"), axis_param = list(at = c(0, 100, almanac_bar_max))),
      `PRISM % vehicle` = ComplexHeatmap::anno_barplot(pv,
        gp = grid::gpar(fill = "#27818A", col = NA), ylim = c(0, prism_bar_max),
        height = grid::unit(11, "mm"), axis_param = list(at = c(0, 100, prism_bar_max))),
      `ALMANAC value` = sensitivity_track(annotation$almanac_abundance_pct[indices]),
      `PRISM value` = sensitivity_track(annotation$prism_abundance_pct[indices]),
      `RNA / control RMS` = ComplexHeatmap::anno_barplot(md$response_control_rms_ratio[indices],
        gp = grid::gpar(fill = "#777777", col = NA), height = grid::unit(8, "mm")),
      show_legend = c(Program = FALSE, Dose = TRUE))
  }
  # All maps share the frozen order, colors and original signed units. Neither
  # tiles nor off-block values are masked to improve the apparent diagonal.
  displacement_heatmap <- function(indices, gene_rows, show_labels, caption) {
    matrix <- t(x[indices, gene_rows, drop = FALSE])
    rm <- row_meta[rownames(matrix), , drop = FALSE]
    cm <- col_meta[indices, , drop = FALSE]
    stopifnot(identical(rownames(matrix), rownames(rm)), identical(colnames(matrix), rownames(cm)))
    draw_heatmap(matrix, col_anno_df = cm, col_anno_cols = names(cm), row_anno_df = rm,
      row_anno_cols = "Program", ann_top_args = top_args(indices), palette_overrides = palette,
      hm_name = "log2 delta", scale_rows = FALSE, value_range = c(-limit, 0, limit),
      value_colors = c("#2166AC", "white", "#B2182B"),
      heatmap_legend_param = list(at = c(-limit, 0, limit), labels = sprintf("%.2f", c(-limit, 0, limit))),
      show_row_names = show_labels, show_column_names = show_labels, column_labels = labels[indices],
      row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 5),
      row_split = rm$Program, column_split = cm$Program,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10),
      title = paste0(gsub("_", " ", scope), " | outcome-selected genes | ", caption, "\nSigned original displacement; ",
        "color saturates at +/-", sprintf("%.2f", limit), "; no row scaling"))
  }
  all_columns <- seq_len(nrow(md))
  overview <- displacement_heatmap(all_columns, row_order, FALSE,
    paste(ncol(full_matrix), "conditions;", nrow(full_matrix), "genes; silhouette-selected k=", k,
      "\nWithin program: PRISM sensitivity, ALMANAC-only sensitivity, then no data; lower % vehicle first"))
  save_heatmap(overview, "displacement_overview", if (nrow(md) > 60) 16 else 11.5, 9,
    scope_figure_dir, list(sensitivity_legend))
  # Detail pages show all conditions in the identical overview order, at most
  # 48 columns/page. Full unabridged drug names are in condition_order.tsv.
  pages <- split(all_columns, ceiling(all_columns / 48))
  page_height <- max(8.5, 5 + length(selected_rows) * 0.06)
  cairo_pdf(file.path(scope_figure_dir, "displacement_pages.pdf"), width = 13.5, height = page_height)
  for (page in seq_along(pages)) {
    hm <- displacement_heatmap(pages[[page]], selected_rows, TRUE,
      paste0("detail page ", page, "/", length(pages), "; top 12 seed genes/program; drug labels <=38 chars"))
    ComplexHeatmap::draw(hm, merge_legends = TRUE, annotation_legend_list = list(sensitivity_legend))
  }
  dev.off()
  preview <- displacement_heatmap(pages[[1]], selected_rows, TRUE,
    paste0("detail page 1/", length(pages), "; top 12 seed genes/program"))
  grDevices::png(file.path(scope_figure_dir, "displacement_page_01.png"), width = 13.5, height = page_height,
    units = "in", res = 160)
  ComplexHeatmap::draw(preview, merge_legends = TRUE, annotation_legend_list = list(sensitivity_legend))
  dev.off()
  if (bank$orientation == "cell") {
    focal_columns <- which(md$focal)
    hm <- displacement_heatmap(focal_columns, selected_rows, TRUE, "seven-drug focus; top 12 seed genes/program")
    save_heatmap(hm, "focal_panel", 11.5, page_height, scope_figure_dir, list(sensitivity_legend))
  }
  # Full gene correlation and resampling consensus share the gene block order.
  mark <- bank$genes[unlist(lapply(seq_len(k), function(j) {
    idx <- which(bank$partition == j)
    head(idx[order(bank$weights[idx, j], decreasing = TRUE)], 3L)
  }))]
  mark_table <- data.frame(gene = mark,
    Program = program_ids[bank$partition[match(mark, bank$genes)]])
  for (kind in c("gene_correlation", "resampling_consensus")) {
    matrix <- if (kind == "gene_correlation") bank$correlation else bank$consensus
    matrix <- matrix[row_order, row_order, drop = FALSE]
    stopifnot(identical(rownames(matrix), rownames(row_meta)) || is.null(rownames(matrix)))
    dimnames(matrix) <- list(bank$genes[row_order], bank$genes[row_order])
    hm <- draw_heatmap(matrix, col_anno_df = row_meta, col_anno_cols = "Program",
      row_anno_df = row_meta, row_anno_cols = "Program", palette_overrides = palette,
      ann_top_args = list(show_legend = FALSE), top_labels = mark_table,
      hm_name = if (kind == "gene_correlation") "Pearson r" else "co-clustering",
      scale_rows = FALSE, value_range = if (kind == "gene_correlation") c(-1, 0, 1) else c(0, 0.5, 1),
      value_colors = if (kind == "gene_correlation") c("#2166AC", "white", "#B2182B") else
        c("white", "#9ECAE1", "#08519C"),
      row_split = row_meta$Program, column_split = row_meta$Program,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 6),
      show_row_names = FALSE,
      title = paste(gsub("_", " ", scope), "|", gsub("_", " ", kind), "| all", length(bank$genes), "genes"))
    save_heatmap(hm, kind, 7.8, 7.2, scope_figure_dir)
  }
  # Program scores use the same all-drug column order and signed response units.
  score_map <- t(scores)
  score_limit <- max(0.1, as.numeric(quantile(abs(score_map), 0.98)))
  hm <- draw_heatmap(score_map, col_anno_df = col_meta, col_anno_cols = names(col_meta),
    ann_top_args = top_args(all_columns), palette_overrides = palette,
    hm_name = "weighted delta", scale_rows = FALSE, value_range = c(-score_limit, 0, score_limit),
    value_colors = c("#2166AC", "white", "#B2182B"),
    heatmap_legend_param = list(at = c(-score_limit, 0, score_limit), labels = sprintf("%.2f", c(-score_limit, 0, score_limit))),
    show_row_names = TRUE, show_column_names = nrow(md) <= 48, column_labels = labels,
    column_names_gp = grid::gpar(fontsize = 6), column_split = col_meta$Program,
    cluster_column_slices = FALSE,
    title = paste(gsub("_", " ", scope), "| original signed program scores |", nrow(md), "conditions"))
  save_heatmap(hm, "program_scores", if (nrow(md) > 60) 16 else 11.5, 6, scope_figure_dir, list(sensitivity_legend))
  # Export complete sources, rank partitions, and identical plotted identifiers.
  fwrite(cbind(data.table(gene = bank$genes), as.data.table(bank$weights)),
    file.path(scope_data_dir, "gene_weights.tsv.gz"), sep = "\t")
  corr <- bank$correlation[row_order, row_order, drop = FALSE]
  fwrite(cbind(data.table(gene = rownames(corr)), as.data.table(corr)),
    file.path(scope_data_dir, "gene_correlation.tsv.gz"), sep = "\t")
  fwrite(cbind(data.table(gene = rownames(full_matrix)), as.data.table(full_matrix)),
    file.path(scope_data_dir, "signed_displacements.tsv.gz"), sep = "\t")
  fwrite(cbind(md[, .(observation_id, recipient, drug, dose)], as.data.table(scores)),
    file.path(scope_data_dir, "program_scores.tsv"), sep = "\t")
  gene_meta <- bank$assignments[row_order]
  gene_meta[, `:=`(plotted_order = seq_len(.N), detail_label = gene %chin% rownames(representative))]
  fwrite(gene_meta, file.path(scope_data_dir, "gene_order.tsv"), sep = "\t")
  fwrite(cbind(md[, .(observation_id, recipient, drug, dose, focal)],
    data.table(display_label = labels, plotted_order = seq_len(nrow(md)),
      detail_page = ceiling(seq_len(nrow(md)) / 48),
      dominant_program = program_ids[dominant]),
    annotation[, .(almanac_abundance_pct, almanac_status, prism_abundance_pct, prism_status,
      sensitivity_sort_source, sensitivity_sort_abundance_pct)]),
    file.path(scope_data_dir, "condition_order.tsv"), sep = "\t")
  fwrite(bank$ranks, file.path(scope_data_dir, "rank_diagnostics.tsv"), sep = "\t")
  fwrite(bank$programs, file.path(scope_data_dir, "program_diagnostics.tsv"), sep = "\t")
  all_assignments <- rbindlist(lapply(names(bank$all_partitions), function(rank) {
    data.table(rank = as.integer(rank), gene = bank$genes, seed_cluster = bank$all_partitions[[rank]])
  }))
  fwrite(all_assignments, file.path(scope_data_dir, "all_rank_gene_assignments.tsv.gz"), sep = "\t")
  fwrite(rbindlist(lapply(names(bank$all_condition_partitions), function(rank) {
    data.table(rank = as.integer(rank), observation_id = bank$metadata$observation_id,
      sample_cluster = bank$all_condition_partitions[[rank]])
  })), file.path(scope_data_dir, "all_rank_condition_assignments.tsv.gz"), sep = "\t")
  invisible(NULL)
}

# If no allowed blue-point silhouette cut is supported, retain raw molecular inspection
# with sensitivity annotations. Dendrogram order is not a selected program cut.
render_unresolved <- function(bank) {
  md <- bank$metadata
  annotation <- sensitivity[match(md$observation_id, observation_id)]
  check(identical(annotation$observation_id, md$observation_id), "Unresolved-map sensitivity mismatch")
  genes <- bank$genes[bank$tree$order]
  response <- prepared$response[match(md$observation_id, rownames(prepared$response)), genes, drop = FALSE]
  condition_order <- bank$condition_tree$order
  md <- md[condition_order]
  annotation <- annotation[condition_order]
  x <- t(prepared$response[match(md$observation_id, rownames(prepared$response)), genes, drop = FALSE])
  limit <- max(0.25, as.numeric(quantile(abs(x), 0.98)))
  unresolved_map <- function(indices, gene_rows, show_labels, caption) {
    av <- annotation$almanac_abundance_pct[indices]
    pv <- annotation$prism_abundance_pct[indices]
    av[is.na(av)] <- 0
    pv[is.na(pv)] <- 0
    column_metadata <- data.frame(Dose = sprintf("%g", md$dose[indices]), row.names = md$observation_id[indices])
    annotation_args <- list(`ALMANAC % vehicle` = ComplexHeatmap::anno_barplot(av,
      gp = grid::gpar(fill = "#B24C35", col = NA), ylim = c(0, almanac_bar_max), height = grid::unit(11, "mm")),
      `PRISM % vehicle` = ComplexHeatmap::anno_barplot(pv,
        gp = grid::gpar(fill = "#27818A", col = NA), ylim = c(0, prism_bar_max), height = grid::unit(11, "mm")),
      `ALMANAC value` = sensitivity_track(annotation$almanac_abundance_pct[indices]),
      `PRISM value` = sensitivity_track(annotation$prism_abundance_pct[indices]))
    display_drugs <- gsub("\u03b1", "alpha", md$drug[indices], fixed = TRUE)
    display_drugs <- gsub("\u03b2", "beta", display_drugs, fixed = TRUE)
    display_drugs <- gsub("\u03b3", "gamma", display_drugs, fixed = TRUE)
    draw_heatmap(x[gene_rows, indices, drop = FALSE], cluster_rows = FALSE, cluster_columns = FALSE,
      scale_rows = FALSE, value_range = c(-limit, 0, limit), value_colors = c("#2166AC", "white", "#B2182B"),
      heatmap_legend_param = list(at = c(-limit, 0, limit), labels = sprintf("%.2f", c(-limit, 0, limit))),
      hm_name = "log2 delta", show_row_names = show_labels, show_column_names = show_labels,
      column_labels = paste(substr(display_drugs, 1, 38), sprintf("%g uM", md$dose[indices])),
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
      col_anno_df = column_metadata, col_anno_cols = names(column_metadata), ann_top_args = annotation_args,
      title = paste0(gsub("_", " ", bank$scope), " | No selected program rank\n", caption,
        "\nColor saturates at +/-", sprintf("%.2f", limit), " in original signed displacement units"))
  }
  heatmap <- unresolved_map(seq_len(ncol(x)), seq_len(nrow(x)), FALSE, "raw clustered overview; no row scaling")
  destination <- file.path(output_dir, bank$scope, "figures")
  scope_data_dir <- file.path(output_dir, bank$scope, "data")
  sample_annotation <- data.frame(Dose = sprintf("%g", md$dose), row.names = md$observation_id)
  render_sample_correlation(bank, md$observation_id, sample_annotation, list(), destination, scope_data_dir)
  save_heatmap(heatmap, "displacement_unresolved", 12.8, 8, destination, list(sensitivity_legend))
  # A readable 21-condition supplement uses 48 high-variance admitted genes,
  # retaining dendrogram order. This is a display subset, not a selected module.
  focal_columns <- which(md$focal)
  focal_columns <- focal_columns[order(md$drug[focal_columns], md$dose[focal_columns])]
  focal_genes <- head(genes[order(matrixStats::rowVars(x), decreasing = TRUE)], 48L)
  focal_gene_rows <- which(genes %chin% focal_genes)
  heatmap <- unresolved_map(focal_columns, focal_gene_rows, TRUE,
    "seven-drug focus; 48 high-variance admitted genes; no row scaling")
  save_heatmap(heatmap, "focal_panel", 10.8, 9, destination, list(sensitivity_legend))
  correlation <- bank$correlation[genes, genes, drop = FALSE]
  heatmap <- draw_heatmap(correlation, cluster_rows = FALSE, cluster_columns = FALSE,
    scale_rows = FALSE, value_range = c(-1, 0, 1), value_colors = c("#2166AC", "white", "#B2182B"),
    hm_name = "Pearson r", show_row_names = FALSE, show_column_names = FALSE,
    title = paste(gsub("_", " ", bank$scope), "| gene correlation; no selected cut"))
  save_heatmap(heatmap, "gene_correlation", 7.2, 7.2, destination)
  scope_data_dir <- file.path(output_dir, bank$scope, "data")
  fwrite(data.table(gene = genes, plotted_order = seq_along(genes), focal_display = genes %chin% focal_genes),
    file.path(scope_data_dir, "gene_order.tsv"), sep = "\t")
  fwrite(cbind(md[, .(observation_id, recipient, drug, dose, focal)],
    annotation[, .(almanac_abundance_pct, almanac_status, prism_abundance_pct, prism_status)]),
    file.path(scope_data_dir, "condition_order.tsv"), sep = "\t")
  fwrite(data.table(gene = rownames(x), x), file.path(scope_data_dir, "signed_displacements.tsv.gz"), sep = "\t")
  fwrite(data.table(gene = genes, correlation), file.path(scope_data_dir, "gene_correlation.tsv.gz"), sep = "\t")
  fwrite(rbindlist(lapply(names(bank$all_partitions), function(rank) {
    data.table(rank = as.integer(rank), gene = bank$genes, seed_cluster = bank$all_partitions[[rank]])
  })), file.path(scope_data_dir, "all_rank_gene_assignments.tsv.gz"), sep = "\t")
}
for (bank in discoveries) {
  if (is.null(bank$weights)) {
    message("Program map withheld for ", bank$scope, ": ", bank$diagnostics$selection_status)
    render_unresolved(bank)
    next
  }
  message("Rendering ", bank$scope)
  render_local(bank)
}

# Raw RNA signal and functional abundance are shown on separate descriptive
# scales. The panel is a complete 7 x 7 x 3 view, with no outcome-based sorting.
signal <- metadata[focal == TRUE, .(observation_id, raw_response_rms, response_control_rms_ratio)]
signal <- merge(sensitivity[focal == TRUE], signal, by = "observation_id", all.x = TRUE, sort = FALSE)
signal[, drug := factor(drug, levels = local_drugs)]
p <- ggplot(signal, aes(dose, response_control_rms_ratio, color = recipient, group = recipient)) +
  geom_hline(yintercept = 1, lty = 2, color = "grey60") + geom_line() + geom_point() +
  scale_x_log10(breaks = c(0.05, 0.5, 5)) + facet_wrap(~drug, ncol = 4) +
  labs(x = "Tahoe dose (uM)", y = "Whole-transcriptome RNA RMS / plate-control RMS",
    title = "Molecular support across the complete focal panel", color = "Recipient")
save_plot(p, "panel_molecular_support", 11.5, 5.7)
p <- ggplot(signal, aes(dose, prism_abundance_pct, color = recipient, group = recipient)) +
  geom_hline(yintercept = 100, lty = 2, color = "grey60") + geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  scale_x_log10(breaks = c(0.05, 0.5, 5)) + facet_wrap(~drug, ncol = 4) +
  labs(x = "Tahoe query dose (uM)", y = "PRISM abundance (% vehicle)",
    title = "Sensitivity annotations at Tahoe doses", subtitle = "Bracketed log-dose interpolation; no extrapolation",
    color = "Recipient")
save_plot(p, "panel_sensitivity", 11.5, 5.7)

# Program/null coherence exposes shared-control artifacts rather than claiming
# that a clean diagonal establishes biological validity.
programs <- rbindlist(lapply(discoveries, `[[`, "programs"))
if (nrow(programs) >= 3L) {
  association <- cor(programs$mean_control_gene_correlation, programs$mean_gene_correlation)
  p <- ggplot(programs, aes(mean_control_gene_correlation, mean_gene_correlation, color = scope)) +
    geom_point(alpha = 0.75) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    theme(aspect.ratio = 1) + labs(x = "Within-module correlation in vehicle pseudo-responses",
      y = "Within-module correlation in treated responses", color = "Cell-line bank",
      title = sprintf("Control structure versus program coherence | Pearson r=%.2f", association))
  save_plot(p, "program_control_structure", 5.8, 4.8)
} else {
  message("Program/control association figure unavailable: fewer than three selected programs")
}

# Record exploratory dose trends and plate confounding for original program
# scores. With three doses these correlations are descriptive, not tests.
dose_rows <- list()
plate_rows <- list()
for (bank in discoveries) {
  if (is.null(bank$weights)) next
  for (j in seq_len(ncol(bank$scores))) {
    md <- copy(bank$metadata)
    md[, `:=`(value = bank$scores[, j], scope = bank$scope, program = colnames(bank$scores)[j])]
    dose_rows[[length(dose_rows) + 1L]] <- md[, .(
      n_doses = uniqueN(dose), log_dose_correlation = if (uniqueN(dose) >= 3 && sd(value) > 0)
        cor(log(dose), value) else NA_real_, minimum_score = min(value), maximum_score = max(value)),
      by = .(scope, program, recipient, drug)]
    md[, primary_plate := sub(",.*", "", contributing_plates)]
    fit <- lm(value ~ factor(primary_plate), data = md)
    plate_rows[[length(plate_rows) + 1L]] <- data.table(scope = bank$scope,
      program = colnames(bank$scores)[j], primary_plate_r_squared = summary(fit)$r.squared,
      interpretation = "Drug/dose and plate are confounded; association is not a causal batch estimate")
  }
}
dose_summary <- if (length(dose_rows)) rbindlist(dose_rows) else data.table(scope = character(), program = character(),
  recipient = character(), drug = character(), n_doses = integer(), log_dose_correlation = numeric(),
  minimum_score = numeric(), maximum_score = numeric())
plate_summary <- if (length(plate_rows)) rbindlist(plate_rows) else data.table(scope = character(), program = character(),
  primary_plate_r_squared = numeric(), interpretation = character())
fwrite(dose_summary, file.path(data_dir, "program_dose_trends.tsv"), sep = "\t")
fwrite(plate_summary, file.path(data_dir, "program_plate_associations.tsv"), sep = "\t")

# The index points to each cell's own outputs; shared tables summarize the run.
artifact_index <- rbindlist(lapply(discoveries, function(bank) {
  artifacts <- c("silhouette.pdf", "cophenetic.pdf", "diagnostics.pdf", "gene_correlation.pdf", "condition_correlation.pdf",
    "resampling_consensus.pdf", "displacement_overview.pdf", "displacement_pages.pdf", "program_scores.pdf", "focal_panel.pdf")
  if (is.null(bank$weights)) artifacts <- c("silhouette.pdf", "cophenetic.pdf", "diagnostics.pdf",
    "gene_correlation.pdf", "condition_correlation.pdf", "displacement_unresolved.pdf", "focal_panel.pdf")
  data.table(recipient = bank$metadata$recipient[1], scope = bank$scope,
    selected_rank = bank$diagnostics$selected_rank, selection_status = bank$diagnostics$selection_status,
    genes = length(bank$genes),
    artifact = artifacts, path = file.path(output_dir, bank$scope, "figures", artifacts))
}))
fwrite(artifact_index, file.path(output_dir, "artifact_index.tsv"), sep = "\t")
writeLines(c("Cell-line-first Tahoe displacement pilot: 379 cached drugs, 7651 profiles, seven recipients.",
  "Creating script: /work/src/program_pilot/06_mockup_sample_programs.R",
  "Each cell_* folder owns its models, data and figures. artifact_index.tsv lists figure paths.",
  sprintf("Gene admission: separate PRISM/ALMANAC tests; p <=%g AND |r| >=%g in the same assay, either sign.",
    outcome_p_threshold, outcome_minimum_abs_correlation),
  "All matched drug/dose rows pooled with equal weight; no dose averaging, splitting or adjustment.",
  "BH adjusted values are diagnostic only; keep up to 500 eligible genes by best qualifying assay p-value.",
  "No displacement SD, variance-ratio, drug-recurrence or abundance-bin gate is used.",
  "Rank: last blue dot in mean gene silhouette; largest allowed rank with positive silhouette and all seeds >=12 genes.",
  "Ranks 2--5 barred; scan all feasible cuts from 2--80, including cuts failing size QC.",
  "Same Ward.D2 method on sqrt(2*(1-r)), now clustering sample/sample Pearson correlations across selected genes.",
  "Genes map to their strongest positive weighted Pearson association with sample-cluster indicators.",
  "Gene silhouette, minimum gene count and overlapping gene-seed score weights retain the stage-05 definitions.",
  "condition_correlation.pdf/png adds the signed sample/sample matrix in the same column order as existing maps.",
  "silhouette.pdf/png: selection curve; cophenetic.pdf/png: diagnostic only; diagnostics.pdf/png: four panels.",
  "displacement_pages.pdf: readable pages with sensitivity bars where observed/bracketed.",
  "Within each program: PRISM first, ALMANAC-only second, no data last; ascending % vehicle (most sensitive first).",
  "Sensitivity: purple <=50% vehicle; green >50%; grey cross means no value. Measured zero stays purple.",
  "gene_correlation.pdf/png: all retained genes in block order; off-block values are preserved.",
  "data/gene_selection.tsv.gz includes every candidate gene and its exclusion reason.",
  "Candidate programs only: control structure remains and resampling stability varies.",
  "Bootstrap draws use fixed per-cell seeds and are saved independently of the rank grid in resampling_units.tsv.gz.",
  "Outcomes are retrospective training/screening evidence used for gene admission; these maps are not held-out validation.",
  "Unsupported selections have raw annotated maps but no program dictionary.",
  "The previous experiments remain unchanged in ../04_mockup_local_programs and ../05_mockup_outcome_programs."),
  file.path(output_dir, "README.txt"))

# License provenance is an inventory, not a commercial-use clearance. No new
# resource was downloaded. R package copyleft obligations are recorded.
packages <- c("data.table", "qs", "ComplexHeatmap", "circlize", "ggplot2", "cluster", "clue", "patchwork", "matrixStats")
package_inventory <- rbindlist(lapply(packages, function(package) {
  description <- packageDescription(package)
  data.table(resource = package, version = description$Version, license = description$License,
    provenance = find.package(package), commercial_review = if (grepl("GPL", description$License))
      "Copyleft: review distribution obligations" else "Review license terms before distribution")
}))
fwrite(package_inventory, file.path(data_dir, "package_licenses.tsv"), sep = "\t")
writeLines(c("Tahoe PDEX, PRISM secondary screen, ALMANAC: existing cached research resources.",
  "No commercial-use clearance was established by this exploratory run.",
  "Resource origins and normalization are retained in the input manifests and source README files."),
  file.path(output_dir, "resource_provenance.txt"))
write_manifest("06_mockup_sample_programs", c(input_files, prep_path, almanac_path, drug_mapping_path,
  prism_info_path, prism_matrix_path, function_mapping_path), list(drugs = sort(unique(metadata$drug)),
    focal_supplement = local_drugs, recipients = cell_map$depmap_name,
    gene_limit_per_line = local_gene_limit,
    gene_selection = "separate PRISM and ALMANAC Pearson tests; either sign; nominal p threshold; BH diagnostic only",
    outcome_p_threshold = outcome_p_threshold, outcome_minimum_abs_correlation = outcome_minimum_abs_correlation,
    gene_cap_order = "best qualifying-assay p-value, largest qualifying absolute correlation, gene symbol",
    outcome_correlation_dose_handling = "all matched drug/dose rows pooled with equal condition weight; no dose adjustment",
    resamples = local_resamples, ranks = local_rank_grid,
    programs = "Ward.D2 sample/sample Pearson clustering; associated gene seeds; unchanged overlapping gene weights",
    rank_rule = "last blue dot in mean gene silhouette on outcome-selected genes: largest allowed positive size-qualified cut",
    minimum_selectable_rank = local_minimum_rank, minimum_seed_genes = local_minimum_genes,
    sensitivity = "exact ALMANAC; PRISM bracketing log-dose/log-abundance interpolation; no extrapolation",
    within_program_display_order = "PRISM ascending % vehicle, ALMANAC-only ascending % vehicle, then no data; ALMANAC breaks PRISM ties",
    render_only = render_only,
    frozen_banks = freeze[, .(scope, sha256)],
    scientific_status = "retrospective outcome-informed exploratory programs; nominal screening; no independent evaluation"))
message("Complete: seven outcome-informed banks, full-panel maps, silhouette diagnostics and gene-outcome source tables")
