#!/usr/bin/env Rscript

## Learn three direct signed untreated-baseline susceptibility modules from
## fold-local orthogonal molecular-susceptibility axes, then fit the existing
## global molecular bridge. Functional outcomes are never loaded; the only
## functional guidance is the already frozen Stage 05 coefficient table.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(matrixStats)
  library(ComplexHeatmap)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
source(file.path(WORK_ROOT, "src", "R", "lib_heatmap.R"))
RESULT <- stage_result("16_baseline_susceptibility_programs")
ARTIFACT <- stage_artifact("16_baseline_susceptibility_programs")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(ARTIFACT, pattern = "\\.(npz|joblib)$", full.names = TRUE),
  force = TRUE)
unlink(list.dirs(ARTIFACT, recursive = FALSE, full.names = TRUE),
  recursive = TRUE, force = TRUE)
unlink(file.path(RESULT, c("stage16_complete.flag",
  "no_identifiable_programs.flag")), force = TRUE)
unlink(file.path(RESULT, c(
  "16a_response_informed_gene_score_distribution.pdf",
  "16b_selected_program_stability.pdf",
  "16c_top_gene_susceptibility_association_signatures.pdf",
  "16d_program_response_configuration_signatures.pdf",
  "16e_program_score_correlation_matrix.pdf",
  "16f_top_gene_loading_heatmap.pdf",
  "16g_baseline_by_response_coefficient_heatmap.pdf",
  "16h_held_observed_vs_predicted_molecular_z.pdf",
  "16i_per_line_molecular_z_performance.pdf",
  "16j_tahoe_ntc_vs_depmap_program_score_concordance.pdf",
  "gene_response_susceptibility_signatures.tsv.gz",
  "molecular_predictions.tsv.gz", "program_cross_fold_matching.tsv",
  "program_loadings.tsv.gz", "program_redundancy.tsv",
  "program_response_configuration_signatures.tsv",
  "program_scores_depmap.tsv", "program_scores_ntc.tsv",
  "program_stability.tsv", "selected_gene_orthogonal_axis_balance.tsv",
  "selected_gene_overlap.tsv",
  "selected_gene_response_dimension_distribution.tsv",
  "selected_genes_by_fold.tsv.gz", "selected_signature_dimension_audit.tsv",
  "signature_norm_exclusions.tsv", "susceptibility_program_rank_selection.tsv",
  "baseline_gene_susceptibility_scores.tsv.gz")), force = TRUE)

N_ORTHOGONAL_SELECTION_AXES <- 4L
MODULE_AXES <- paste0("O", 1:3)
BASELINE_PROGRAMS <- paste0("baseline_susceptibility_", MODULE_AXES)
MAX_SINGLE_LINE_FRACTION <- 0.50
MIN_DIRECTION_STABILITY <- 0.80
MIN_LOO_SCORE_FRACTION <- 0.25
N_NATURAL_PLATES <- 14L

theme_stage16 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    plot.title = element_text(face = "bold"))
square_scatter_theme <- function(base_size = 8.0) theme_stage16(base_size) +
  theme(aspect.ratio = 1)
save_gg <- function(name, plot, width, height) ggsave(file.path(RESULT, name),
  plot, width = width, height = height, device = cairo_pdf, bg = "white",
  limitsize = FALSE)
save_heatmap <- function(name, heatmap, width, height) {
  grDevices::pdf(file.path(RESULT, name), width = width, height = height,
    useDingbats = FALSE)
  ComplexHeatmap::draw(heatmap)
  grDevices::dev.off()
}
fmt <- function(x, digits = 3L) ifelse(is.finite(x),
  formatC(x, digits = digits, format = "f"), "NA")

read_npy_selected_log_counts <- function(path, columns_one_based) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  if (!identical(magic, "\x93NUMPY") || !version[1L] %in% c(1L, 2L, 3L))
    fail("Unsupported NPY header: %s", path)
  length_size <- if (version[1L] == 1L) 2L else 4L
  header_length <- readBin(con, "integer", 1L, size = length_size,
    signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", header_length))
  if (!grepl("'descr': '<f4'", header, fixed = TRUE) ||
      !grepl("'fortran_order': False", header, fixed = TRUE))
    fail("Stage 16 requires a C-order little-endian float32 NPY matrix.")
  shape_text <- sub(".*'shape': \\(([^)]*)\\).*", "\\1", header)
  shape <- suppressWarnings(as.integer(trimws(strsplit(shape_text, ",")[[1L]])))
  shape <- shape[is.finite(shape)]
  columns_one_based <- as.integer(columns_one_based)
  if (length(shape) != 2L || any(columns_one_based < 1L) ||
      any(columns_one_based > shape[2L]) || anyDuplicated(columns_one_based))
    fail("Invalid Stage 16 selected-column NPY request.")
  result <- matrix(NA_real_, nrow = shape[1L],
    ncol = length(columns_one_based))
  for (row in seq_len(shape[1L])) {
    current <- readBin(con, "numeric", n = shape[2L], size = 4L,
      endian = "little")
    if (length(current) != shape[2L]) fail("Baseline NPY row %d is truncated.", row)
    library <- sum(current)
    if (!is.finite(library) || library <= 0) fail("Invalid baseline library size.")
    result[row, ] <- log2(1 + current[columns_one_based] / library * 10000)
  }
  result
}

make_footprint <- function(lines, score_table, displacement, programs) {
  use <- score_table$depmap_id %chin% lines
  means <- data.table(perturbation_id = score_table$perturbation_id[use],
    displacement[use, programs, drop = FALSE])[, lapply(.SD, mean),
      by = perturbation_id, .SDcols = programs]
  matrix <- as.matrix(means[, ..programs]); rownames(matrix) <-
    means$perturbation_id
  matrix
}

axis_vector <- function(axis, fold_label, programs) {
  current <- axis[fold == fold_label]
  if (nrow(current) != length(programs) ||
      !setequal(current$program, programs))
    fail("Frozen functional coefficients are incomplete for %s.", fold_label)
  setNames(current$coefficient[match(programs, current$program)], programs)
}

column_cor <- function(left, right) {
  left_centered <- sweep(left, 2L, colMeans(left), "-")
  right_centered <- sweep(right, 2L, colMeans(right), "-")
  denominator <- sqrt(colSums(left_centered ^ 2) *
    colSums(right_centered ^ 2))
  answer <- colSums(left_centered * right_centered) / denominator
  answer[!is.finite(answer)] <- NA_real_
  answer
}

dominance_fraction <- function(matrix) {
  centered <- sweep(matrix, 2L, colMeans(matrix), "-")
  squared <- centered ^ 2
  answer <- colMaxs(squared) / colSums(squared)
  answer[!is.finite(answer)] <- 1
  answer
}

gene_profile_correlations <- function(expression, profiles) {
  x <- sweep(expression, 2L, colMeans(expression), "-")
  y <- sweep(profiles, 2L, colMeans(profiles), "-")
  denominator <- outer(sqrt(colSums(x ^ 2)), sqrt(colSums(y ^ 2)))
  result <- crossprod(x, y) / denominator
  result[!is.finite(result)] <- 0
  result
}

fit_line_profiles <- function(lines, fold_label, gate_ids, footprint, beta,
                              score_table, displacement, programs) {
  use <- score_table$depmap_id %chin% lines &
    score_table$perturbation_id %chin% gate_ids &
    score_table$perturbation_id %chin% rownames(footprint)
  R <- footprint[score_table$perturbation_id[use], programs, drop = FALSE]
  z <- as.numeric((displacement[use, programs, drop = FALSE] - R) %*% beta)
  profile_rows <- data.table(depmap_id = score_table$depmap_id[use],
    perturbation_id = score_table$perturbation_id[use], z = z)
  counts <- profile_rows[, .N, by = perturbation_id]
  weights <- 1 / counts$N[match(profile_rows$perturbation_id,
    counts$perturbation_id)]
  weights <- weights / mean(weights)
  x_center <- colSums(R * weights) / sum(weights)
  x_scale <- sqrt(colSums(sweep(R, 2L, x_center, "-") ^ 2 * weights) /
    sum(weights))
  z_center <- sum(z * weights) / sum(weights)
  if (any(!is.finite(x_scale)) || any(x_scale <= 0))
    fail("%s susceptibility-profile footprint cannot be standardized.", fold_label)
  X <- sweep(sweep(R, 2L, x_center, "-"), 2L, x_scale, "/")
  y <- z - z_center
  components <- lapply(lines, function(line) {
    current <- which(profile_rows$depmap_id == line)
    if (length(current) < length(programs) + 2L)
      fail("%s has too few profile rows for %s.", fold_label, line)
    w <- weights[current]; w <- w / mean(w); root <- sqrt(w)
    decomposition <- eigen(crossprod(X[current, , drop = FALSE] * root),
      symmetric = TRUE)
    list(line = line, rows = current, weights = w,
      vectors = decomposition$vectors,
      values = pmax(decomposition$values, 0),
      rhs = as.numeric(crossprod(decomposition$vectors,
        crossprod(X[current, , drop = FALSE] * root, y[current] * root))))
  })
  total_sse <- total_df <- numeric(length(DECODER_LAMBDA_GRID))
  total_weight <- sum(weights); total_rows <- length(z)
  for (lambda_index in seq_along(DECODER_LAMBDA_GRID)) {
    lambda <- DECODER_LAMBDA_GRID[lambda_index]
    for (component in components) {
      beta_std <- as.numeric(component$vectors %*%
        (component$rhs / (component$values + lambda)))
      prediction <- as.numeric(X[component$rows, , drop = FALSE] %*% beta_std)
      total_sse[lambda_index] <- total_sse[lambda_index] + sum(
        component$weights * (y[component$rows] - prediction) ^ 2)
      total_df[lambda_index] <- total_df[lambda_index] + sum(
        component$values / (component$values + lambda))
    }
  }
  gcv <- total_sse / total_weight /
    pmax(1 - total_df / total_rows, 1e-6) ^ 2
  selected_index <- order(gcv, -DECODER_LAMBDA_GRID)[1L]
  selected_lambda <- DECODER_LAMBDA_GRID[selected_index]
  profile <- rbindlist(lapply(components, function(component) {
    beta_std <- as.numeric(component$vectors %*%
      (component$rhs / (component$values + selected_lambda)))
    data.table(fold = fold_label, depmap_id = component$line,
      response_program = programs,
      susceptibility_profile_coefficient = beta_std / x_scale,
      standardized_coefficient = beta_std,
      profile_ridge_penalty = selected_lambda,
      n_profile_rows = length(component$rows),
      condition_balanced_weighting = TRUE,
      no_intercept_after_training_centering = TRUE,
      footprint_center = x_center,
      footprint_scale = x_scale, molecular_z_center = z_center)
  }))
  selection <- data.table(fold = fold_label,
    lambda = DECODER_LAMBDA_GRID, gcv = gcv,
    weighted_training_mse = total_sse / total_weight,
    effective_degrees_of_freedom = total_df,
    selected = seq_along(DECODER_LAMBDA_GRID) == selected_index,
    n_learning_lines = length(lines), n_molecular_rows = length(z),
    n_eligible_conditions = length(gate_ids))
  list(profile = profile, selection = selection)
}

score_baseline_genes <- function(fold_label, learning_lines, profile_table,
                                 ntc_expression, depmap_expression,
                                 plate_variance, genes) {
  profile_wide <- dcast(profile_table,
    depmap_id ~ response_program,
    value.var = "susceptibility_profile_coefficient")
  response_programs <- setdiff(names(profile_wide), "depmap_id")
  profile_matrix <- as.matrix(profile_wide[, ..response_programs])
  rownames(profile_matrix) <- profile_wide$depmap_id
  profile_matrix <- profile_matrix[learning_lines, response_programs,
    drop = FALSE]
  profile_variance <- colVars(profile_matrix)
  keep_dimensions <- is.finite(profile_variance) & profile_variance > 0
  if (!any(keep_dimensions)) fail("%s has no variable susceptibility profile.", fold_label)
  profile_matrix <- profile_matrix[, keep_dimensions, drop = FALSE]
  profile_variance <- profile_variance[keep_dimensions]
  response_programs <- colnames(profile_matrix)
  if (ncol(profile_matrix) < N_ORTHOGONAL_SELECTION_AXES)
    fail("%s has fewer than four susceptibility dimensions.", fold_label)
  profile_center <- colMeans(profile_matrix)
  profile_scale <- colSds(profile_matrix)
  standardized_profiles <- sweep(sweep(profile_matrix, 2L,
    profile_center, "-"), 2L, profile_scale, "/")
  profile_svd <- svd(standardized_profiles, nu = 0L,
    nv = N_ORTHOGONAL_SELECTION_AXES)
  orthogonal_loadings <- profile_svd$v[,
    seq_len(N_ORTHOGONAL_SELECTION_AXES), drop = FALSE]
  for (axis_index in seq_len(N_ORTHOGONAL_SELECTION_AXES)) {
    anchor <- which.max(abs(orthogonal_loadings[, axis_index]))
    if (orthogonal_loadings[anchor, axis_index] < 0)
      orthogonal_loadings[, axis_index] <-
        -orthogonal_loadings[, axis_index]
  }
  orthogonal_axis_names <- paste0("O",
    seq_len(N_ORTHOGONAL_SELECTION_AXES))
  orthogonal_profiles <- standardized_profiles %*% orthogonal_loadings
  colnames(orthogonal_profiles) <- orthogonal_axis_names
  orthogonal_variance_fraction <- profile_svd$d ^ 2 /
    sum(profile_svd$d ^ 2)
  orthogonal_axis_definition <- CJ(
    response_index = seq_along(response_programs),
    axis_index = seq_len(N_ORTHOGONAL_SELECTION_AXES))
  orthogonal_axis_definition[, `:=`(fold = fold_label,
    response_program = response_programs[response_index],
    susceptibility_axis = orthogonal_axis_names[axis_index],
    loading = orthogonal_loadings[cbind(response_index, axis_index)],
    singular_value = profile_svd$d[axis_index],
    target_variance_fraction = orthogonal_variance_fraction[axis_index],
    n_learning_lines = length(learning_lines))]
  orthogonal_axis_definition[, c("response_index", "axis_index") := NULL]

  X <- ntc_expression[learning_lines, , drop = FALSE]
  X_depmap <- depmap_expression[learning_lines, , drop = FALSE]
  abundance <- colMeans(X)
  between_variance <- colVars(X)
  mean_plate_variance <- colMeans(plate_variance[learning_lines, , drop = FALSE])
  plate_noise <- mean_plate_variance / N_NATURAL_PLATES
  reproducible_variance <- pmax(between_variance - plate_noise, 0)
  cross_source <- column_cor(X, X_depmap)
  ntc_dominance <- dominance_fraction(X)
  depmap_dominance <- dominance_fraction(X_depmap)
  maximum_dominance <- pmax(ntc_dominance, depmap_dominance)
  technical_pass <- is.finite(abundance) & is.finite(between_variance) &
    reproducible_variance > 0 & is.finite(cross_source) & cross_source > 0 &
    maximum_dominance <= MAX_SINGLE_LINE_FRACTION

  association <- gene_profile_correlations(X, profile_matrix)
  orthogonal_association <- gene_profile_correlations(X,
    orthogonal_profiles)
  score <- as.numeric((association ^ 2) %*% profile_variance /
    sum(profile_variance))
  dominant_index <- max.col(abs(association), ties.method = "first")
  dominant_program <- response_programs[dominant_index]
  dominant_correlation <- association[cbind(seq_len(nrow(association)),
    dominant_index)]
  orthogonal_dominant_index <- max.col(abs(orthogonal_association),
    ties.method = "first")
  dominant_orthogonal_axis <- orthogonal_axis_names[
    orthogonal_dominant_index]
  orthogonal_dominant_correlation <- orthogonal_association[cbind(
    seq_len(nrow(orthogonal_association)), orthogonal_dominant_index)]
  loo_score <- matrix(NA_real_, nrow = length(genes),
    ncol = length(learning_lines))
  loo_fixed_dimension_absolute_correlation <- matrix(NA_real_,
    nrow = length(genes), ncol = length(learning_lines))
  loo_fixed_orthogonal_absolute_correlation <- matrix(NA_real_,
    nrow = length(genes), ncol = length(learning_lines))
  stable_direction <- matrix(FALSE, nrow = length(genes),
    ncol = length(learning_lines))
  orthogonal_stable_direction <- matrix(FALSE, nrow = length(genes),
    ncol = length(learning_lines))
  for (omitted in seq_along(learning_lines)) {
    keep <- seq_along(learning_lines) != omitted
    current_profiles <- profile_matrix[keep, , drop = FALSE]
    current_variance <- colVars(current_profiles)
    current_variance[!is.finite(current_variance) | current_variance <= 0] <- 0
    current_association <- gene_profile_correlations(X[keep, , drop = FALSE],
      current_profiles)
    current_profile_center <- colMeans(current_profiles)
    current_profile_scale <- colSds(current_profiles)
    current_standardized_profiles <- sweep(sweep(current_profiles, 2L,
      current_profile_center, "-"), 2L, current_profile_scale, "/")
    current_orthogonal_profiles <- current_standardized_profiles %*%
      orthogonal_loadings
    current_orthogonal_association <- gene_profile_correlations(
      X[keep, , drop = FALSE], current_orthogonal_profiles)
    loo_score[, omitted] <- as.numeric((current_association ^ 2) %*%
      current_variance / sum(current_variance))
    current_dominant <- max.col(abs(current_association), ties.method = "first")
    current_correlation <- current_association[cbind(
      seq_len(nrow(current_association)), current_dominant)]
    fixed_dimension_correlation <- current_association[cbind(
      seq_len(nrow(current_association)), dominant_index)]
    loo_fixed_dimension_absolute_correlation[, omitted] <-
      abs(fixed_dimension_correlation)
    current_orthogonal_dominant <- max.col(
      abs(current_orthogonal_association), ties.method = "first")
    current_orthogonal_correlation <- current_orthogonal_association[cbind(
      seq_len(nrow(current_orthogonal_association)),
      current_orthogonal_dominant)]
    fixed_orthogonal_correlation <- current_orthogonal_association[cbind(
      seq_len(nrow(current_orthogonal_association)),
      orthogonal_dominant_index)]
    loo_fixed_orthogonal_absolute_correlation[, omitted] <-
      abs(fixed_orthogonal_correlation)
    stable_direction[, omitted] <- response_programs[current_dominant] ==
      dominant_program & sign(current_correlation) == sign(dominant_correlation)
    orthogonal_stable_direction[, omitted] <-
      orthogonal_axis_names[current_orthogonal_dominant] ==
        dominant_orthogonal_axis &
      sign(current_orthogonal_correlation) ==
        sign(orthogonal_dominant_correlation)
  }
  loo_lower <- rowMins(loo_score)
  loo_median <- rowMedians(loo_score)
  dimension_loo_lower <- rowMins(
    loo_fixed_dimension_absolute_correlation)
  direction_stability <- rowMeans(stable_direction)
  orthogonal_loo_lower <- rowMins(
    loo_fixed_orthogonal_absolute_correlation)
  orthogonal_direction_stability <- rowMeans(
    orthogonal_stable_direction)
  single_line_dependent <- !is.finite(loo_lower) |
    loo_lower < MIN_LOO_SCORE_FRACTION * score
  eligible <- technical_pass & is.finite(score) & score > 0 &
    is.finite(loo_lower) & loo_lower > 0 & !single_line_dependent &
    direction_stability >= MIN_DIRECTION_STABILITY
  result <- data.table(fold = fold_label, gene = genes,
    n_learning_lines = length(learning_lines),
    mean_ntc_abundance = abundance,
    ntc_between_line_variance = between_variance,
    mean_within_line_plate_variance = mean_plate_variance,
    estimated_plate_noise_of_line_mean = plate_noise,
    reproducible_ntc_between_line_variance = reproducible_variance,
    ntc_depmap_training_pearson = cross_source,
    ntc_single_line_squared_fraction = ntc_dominance,
    depmap_single_line_squared_fraction = depmap_dominance,
    maximum_single_line_squared_fraction = maximum_dominance,
    technical_reproducibility_pass = technical_pass,
    susceptibility_score = score,
    loo_minimum_susceptibility_score = loo_lower,
    loo_median_susceptibility_score = loo_median,
    dominant_response_program = dominant_program,
    dominant_correlation = dominant_correlation,
    dimension_loo_minimum_absolute_correlation = dimension_loo_lower,
    dominant_orthogonal_axis = dominant_orthogonal_axis,
    orthogonal_axis_correlation = orthogonal_dominant_correlation,
    orthogonal_loo_minimum_absolute_correlation = orthogonal_loo_lower,
    orthogonal_axis_direction = fifelse(
      orthogonal_dominant_correlation >= 0, "positive", "negative"),
    orthogonal_axis_direction_stability_fraction =
      orthogonal_direction_stability,
    dominant_direction = fifelse(dominant_correlation >= 0,
      "positive", "negative"),
    dominant_direction_stability_fraction = direction_stability,
    single_line_dependent = single_line_dependent,
    susceptibility_eligible = eligible)
  ## Preserve every existing eligibility filter and fold-local stability.
  ## Assignment is determined across O1--O4; O4-dominant genes remain excluded
  ## rather than being reassigned to a retained axis.
  candidates <- result[susceptibility_eligible == TRUE &
    orthogonal_axis_direction_stability_fraction >= MIN_DIRECTION_STABILITY &
    is.finite(orthogonal_loo_minimum_absolute_correlation) &
    orthogonal_loo_minimum_absolute_correlation > 0]
  candidates[, absolute_orthogonal_axis_correlation :=
    abs(orthogonal_axis_correlation)]
  candidates[, axis_number := as.integer(sub("O", "",
    dominant_orthogonal_axis))]
  setorder(candidates, axis_number,
    -orthogonal_loo_minimum_absolute_correlation,
    -absolute_orthogonal_axis_correlation,
    -loo_minimum_susceptibility_score, -susceptibility_score, gene)
  candidates[, orthogonal_axis_selection_rank := seq_len(.N),
    by = dominant_orthogonal_axis]
  selected_table <- candidates[dominant_orthogonal_axis %chin% MODULE_AXES,
    .(gene, dominant_orthogonal_axis, orthogonal_axis_correlation,
      orthogonal_loo_minimum_absolute_correlation,
      orthogonal_axis_direction_stability_fraction,
      orthogonal_axis_selection_rank)]
  selected_table[, `:=`(
    baseline_program = paste0("baseline_susceptibility_",
      dominant_orthogonal_axis),
    unnormalized_signed_weight = sign(orthogonal_axis_correlation) *
      orthogonal_loo_minimum_absolute_correlation)]
  selected_table[, normalized_signed_weight :=
    unnormalized_signed_weight /
      sqrt(sum(unnormalized_signed_weight ^ 2)),
    by = baseline_program]
  module_counts <- selected_table[, .N, by = baseline_program]
  if (!setequal(module_counts$baseline_program, BASELINE_PROGRAMS) ||
      any(module_counts$N == 0L) ||
      any(!is.finite(selected_table$normalized_signed_weight)))
    fail("%s has an empty or invalid direct susceptibility module.", fold_label)
  module_norms <- selected_table[, sum(normalized_signed_weight ^ 2),
    by = baseline_program]
  if (max(abs(module_norms$V1 - 1)) > 1e-10)
    fail("%s direct module weights are not unit normalized.", fold_label)
  candidate_match <- match(result$gene, candidates$gene)
  selected_match <- match(result$gene, selected_table$gene)
  result[, `:=`(orthogonal_axis_selection_rank =
      candidates$orthogonal_axis_selection_rank[candidate_match],
    baseline_program = selected_table$baseline_program[selected_match],
    unnormalized_signed_weight =
      selected_table$unnormalized_signed_weight[selected_match],
    normalized_signed_weight =
      selected_table$normalized_signed_weight[selected_match],
    selected = !is.na(selected_match))]
  selected_order <- candidates[gene %chin% selected_table$gene]
  setorder(selected_order, orthogonal_axis_selection_rank, axis_number, gene)
  result[, susceptibility_selection_rank := match(gene,
    selected_order$gene)]
  setorder(result, -selected, susceptibility_selection_rank,
    -susceptibility_eligible, gene)

  attr(result, "orthogonal_axis_definition") <-
    orthogonal_axis_definition
  attr(result, "module_gene_weights") <- selected_table[, .(
    fold = fold_label, gene, baseline_program,
    susceptibility_axis = dominant_orthogonal_axis,
    full_axis_correlation = orthogonal_axis_correlation,
    loo_minimum_absolute_correlation =
      orthogonal_loo_minimum_absolute_correlation,
    loo_direction_stability_fraction =
      orthogonal_axis_direction_stability_fraction,
    unnormalized_signed_weight, normalized_signed_weight,
    direction = fifelse(normalized_signed_weight >= 0, "positive", "negative"),
    assignment_rule = "dominant stable association among O1-O4; retain O1-O3")]
  result
}

make_interactions <- function(line_ids, perturbation_ids, baseline_scores,
                              footprint, baseline_programs, programs) {
  S <- baseline_scores[line_ids, baseline_programs, drop = FALSE]
  R <- footprint[perturbation_ids, programs, drop = FALSE]
  X <- do.call(cbind, lapply(baseline_programs, function(b)
    sweep(R, 1L, S[, b], "*")))
  colnames(X) <- unlist(lapply(baseline_programs, function(b)
    paste0(b, "__", programs)), use.names = FALSE)
  X
}

molecular_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]; predicted <- predicted[keep]
  data.table(n_rows = length(observed), pearson = safe_cor(observed, predicted),
    observed_sd = sd(observed), predicted_sd = sd(predicted),
    predicted_to_observed_sd_ratio = sd(predicted) / sd(observed),
    mse = mean((observed - predicted) ^ 2),
    rmse = sqrt(mean((observed - predicted) ^ 2)),
    correlation_method = "pearson")
}

project_module_scores <- function(expression, learning_lines, all_lines,
                                  weights, fold_label, source_name,
                                  validation_lines = character()) {
  genes <- sort(unique(weights$gene))
  if (!length(genes) || !all(genes %chin% colnames(expression)))
    fail("%s %s projection has missing genes.", fold_label, source_name)
  center <- colMeans(expression[learning_lines, genes, drop = FALSE])
  scale <- colSds(expression[learning_lines, genes, drop = FALSE])
  keep <- is.finite(center) & is.finite(scale) & scale > 0
  if (!all(keep)) fail("%s %s gene scaling is invalid.", fold_label, source_name)
  standardized <- sweep(sweep(expression[all_lines, genes, drop = FALSE],
    2L, center, "-"), 2L, scale, "/")
  weight_wide <- dcast(weights, gene ~ baseline_program,
    value.var = "normalized_signed_weight", fill = 0)
  W <- as.matrix(weight_wide[match(genes, gene), ..BASELINE_PROGRAMS])
  rownames(W) <- genes
  raw <- standardized %*% W
  module_center <- colMeans(raw[learning_lines, , drop = FALSE])
  module_scale <- colSds(raw[learning_lines, , drop = FALSE])
  if (any(!is.finite(module_center)) || any(!is.finite(module_scale)) ||
      any(module_scale <= 0))
    fail("%s %s module scaling is invalid.", fold_label, source_name)
  score <- sweep(sweep(raw, 2L, module_center, "-"),
    2L, module_scale, "/")
  if (any(!is.finite(score)))
    fail("%s %s module scores are nonfinite.", fold_label, source_name)
  split_role <- fifelse(all_lines %chin% learning_lines, "learning",
    fifelse(all_lines %chin% validation_lines, "validation", "other"))
  long <- data.table(fold = fold_label,
    baseline_source = source_name,
    depmap_id = rep(all_lines, each = length(BASELINE_PROGRAMS)),
    program = rep(BASELINE_PROGRAMS, times = length(all_lines)),
    raw_module_score = as.vector(t(raw)),
    standardized_program_score = as.vector(t(score)),
    module_training_center = rep(module_center, times = length(all_lines)),
    module_training_scale = rep(module_scale, times = length(all_lines)),
    split_role = rep(split_role, each = length(BASELINE_PROGRAMS)))
  list(matrix = score, scores = long, gene_center = center,
    gene_scale = scale, module_center = module_center,
    module_scale = module_scale, weights = W)
}

required <- c(file.path(stage_result("01_inputs"), "fold_assignments.tsv"),
  file.path(stage_artifact("01_inputs"), "depmap_baseline_expression.rds"),
  file.path(stage_result("02_response_programs"), "fold_program_stability.tsv"),
  file.path(stage_artifact("02_response_programs"), "fold_program_scores.tsv.gz"),
  file.path(stage_result("04_activity_panel"),
    "molecular_context_variance_census.tsv.gz"),
  file.path(stage_result("05_program_function_bridge"),
    "functional_molecular_axis_coefficients.tsv"),
  file.path(stage_result("10_baseline_states"), "baseline_gene_selection.tsv.gz"),
  file.path(CANONICAL_CACHE, "dmso_control_expression.npy"),
  file.path(CANONICAL_CACHE, "dmso_control_metadata.tsv"),
  file.path(CANONICAL_CACHE, "dmso_control_genes.tsv"))
missing <- required[!file.exists(required)]
if (length(missing)) fail("Missing Stage 16 input(s): %s",
  paste(missing, collapse = ", "))

folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- sort(folds[outer_fold != OUTER_FOLD, depmap_id])
held_lines <- sort(folds[outer_fold == OUTER_FOLD, depmap_id])
all_lines <- sort(folds$depmap_id)
if (length(training_lines) != 26L || length(held_lines) != 7L ||
    !identical(held_lines, sort(EXPECTED_HELD_IDS)) ||
    any(EXCLUDED_IDS %chin% all_lines))
  fail("The accepted Stage 16 cohort changed.")
inner_assignment <- make_inner_assignment(training_lines)

response_stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- response_stability[outer_fold == OUTER_FOLD & stable == TRUE,
  unique(program)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
if (length(programs) != 9L) fail("Stage 16 requires nine accepted response programs.")
score_long <- fread(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"))[program %chin% programs]
score_table <- dcast(score_long, observation_id + depmap_id + cell_line +
  perturbation_id + drug + dose + dose_panel + role ~ program,
  value.var = "score")
D <- as.matrix(score_table[, ..programs]); rownames(D) <- score_table$observation_id
axis <- fread(file.path(stage_result("05_program_function_bridge"),
  "functional_molecular_axis_coefficients.tsv"))
eligibility <- fread(file.path(stage_result("04_activity_panel"),
  "molecular_context_variance_census.tsv.gz"))
if (any(eligibility$prism_used_for_molecular_classification != FALSE) ||
    any(eligibility$validation_or_held_values_used_for_selection != FALSE))
  fail("Stage 16 molecular eligibility is not fold-training-only.")

## Recreate the untreated gene matrix in R without loading any PRISM
## table.  The common universe and technical definitions are inherited from
## Stage 10, while every fold recalculates its learning-line statistics.
universe_table <- fread(file.path(stage_result("10_baseline_states"),
  "baseline_gene_selection.tsv.gz"))
universe_genes <- universe_table$gene
control_genes <- fread(file.path(CANONICAL_CACHE, "dmso_control_genes.tsv"))
gene_position <- match(universe_genes, control_genes$feature)
if (anyNA(gene_position) || anyDuplicated(universe_genes))
  fail("The Stage 10 common untreated gene universe changed.")
message(sprintf("Reading %s common NTC genes in pure R.",
  format(length(universe_genes), big.mark = ",")))
plate_expression <- read_npy_selected_log_counts(file.path(CANONICAL_CACHE,
  "dmso_control_expression.npy"), gene_position)
colnames(plate_expression) <- universe_genes
control_meta <- fread(file.path(CANONICAL_CACHE, "dmso_control_metadata.tsv"))[
  depmap_id %chin% all_lines]
setorder(control_meta, array_row)
plate_expression <- plate_expression[control_meta$array_row + 1L, , drop = FALSE]
line_rows <- split(seq_len(nrow(control_meta)), control_meta$depmap_id)
if (!all(lengths(line_rows[all_lines]) == N_NATURAL_PLATES))
  fail("The natural NTC plate contract changed.")
ntc_expression <- t(vapply(all_lines, function(line)
  colMeans(plate_expression[line_rows[[line]], , drop = FALSE]),
  numeric(ncol(plate_expression))))
plate_variance <- t(vapply(all_lines, function(line)
  colVars(plate_expression[line_rows[[line]], , drop = FALSE]),
  numeric(ncol(plate_expression))))
rownames(ntc_expression) <- rownames(plate_variance) <- all_lines
colnames(ntc_expression) <- colnames(plate_variance) <- universe_genes
rm(plate_expression); gc(FALSE)
depmap_object <- readRDS(file.path(stage_artifact("01_inputs"),
  "depmap_baseline_expression.rds"))
depmap_expression <- depmap_object$expression[all_lines, universe_genes,
  drop = FALSE]
if (any(!is.finite(ntc_expression)) || any(!is.finite(depmap_expression)))
  fail("The untreated baseline matrices contain nonfinite values.")

## Nested direct-axis module construction and ridge-penalty selection.
profile_parts <- list(); profile_selection_parts <- list()
gene_score_parts <- list(); module_weight_parts <- list()
axis_definition_parts <- list(); module_score_parts <- list()
fold_objects <- vector("list", N_INNER_FOLDS)
inner_selection_parts <- list()

for (inner_fold in seq_len(N_INNER_FOLDS)) {
  fold_label <- paste0("inner_", inner_fold)
  validation_lines <- names(inner_assignment)[inner_assignment == inner_fold]
  learning_lines <- names(inner_assignment)[inner_assignment != inner_fold]
  gate_ids <- eligibility[fold == fold_label &
    molecular_context_variable == TRUE, perturbation_id]
  footprint <- make_footprint(learning_lines, score_table, D, programs)
  beta <- axis_vector(axis, fold_label, programs)
  profiles <- fit_line_profiles(learning_lines, fold_label, gate_ids,
    footprint, beta, score_table, D, programs)
  gene_scores <- score_baseline_genes(fold_label, learning_lines,
    profiles$profile, ntc_expression, depmap_expression, plate_variance,
    universe_genes)
  weights <- attr(gene_scores, "module_gene_weights")
  axis_definitions <- attr(gene_scores, "orthogonal_axis_definition")
  projection <- project_module_scores(ntc_expression, learning_lines,
    all_lines, weights, fold_label, "tahoe_ntc", validation_lines)

  profile_parts[[fold_label]] <- profiles$profile
  profile_selection_parts[[fold_label]] <- profiles$selection
  gene_score_parts[[fold_label]] <- gene_scores
  module_weight_parts[[fold_label]] <- weights
  axis_definition_parts[[fold_label]] <- axis_definitions[
    susceptibility_axis %chin% MODULE_AXES]
  module_score_parts[[fold_label]] <- projection$scores

  learning_rows <- score_table$depmap_id %chin% learning_lines &
    score_table$perturbation_id %chin% gate_ids
  validation_rows <- score_table$depmap_id %chin% validation_lines &
    score_table$perturbation_id %chin% gate_ids
  R_learning <- footprint[score_table$perturbation_id[learning_rows],
    programs, drop = FALSE]
  R_validation <- footprint[score_table$perturbation_id[validation_rows],
    programs, drop = FALSE]
  z_learning <- as.numeric((D[learning_rows, programs, drop = FALSE] -
    R_learning) %*% beta)
  z_validation <- as.numeric((D[validation_rows, programs, drop = FALSE] -
    R_validation) %*% beta)
  X_learning <- make_interactions(score_table$depmap_id[learning_rows],
    score_table$perturbation_id[learning_rows], projection$matrix,
    footprint, BASELINE_PROGRAMS, programs)
  X_validation <- make_interactions(score_table$depmap_id[validation_rows],
    score_table$perturbation_id[validation_rows], projection$matrix,
    footprint, BASELINE_PROGRAMS, programs)
  path <- ridge_path(X_learning, z_learning, X_validation,
    DECODER_LAMBDA_GRID)
  breadth <- vapply(seq_along(DECODER_LAMBDA_GRID), function(index) {
    by_line <- data.table(depmap_id = score_table$depmap_id[validation_rows],
      observed = z_validation, predicted = path[, index])[, .(
        pearson = safe_cor(observed, predicted)), by = depmap_id]
    sum(by_line$pearson > 0, na.rm = TRUE)
  }, integer(1L))
  inner_selection_parts[[fold_label]] <- rbindlist(lapply(
    seq_along(DECODER_LAMBDA_GRID), function(index) data.table(
      inner_fold = inner_fold, fold = fold_label,
      lambda = DECODER_LAMBDA_GRID[index],
      validation_mse = mean((z_validation - path[, index]) ^ 2),
      validation_pearson = safe_cor(z_validation, path[, index]),
      validation_observed_sd = sd(z_validation),
      validation_predicted_sd = sd(path[, index]),
      validation_sd_ratio = sd(path[, index]) / sd(z_validation),
      n_validation_lines_positive = breadth[index],
      n_learning_lines = length(learning_lines),
      n_validation_lines = length(validation_lines),
      n_learning_rows = sum(learning_rows),
      n_validation_rows = sum(validation_rows),
      n_eligible_conditions = length(gate_ids),
      representation = "fixed_direct_signed_O1_O3_modules",
      selection_objective = "mean_inner_validation_molecular_z_mse")))
  fold_objects[[inner_fold]] <- list(label = fold_label,
    learning_lines = learning_lines, validation_lines = validation_lines,
    gate_ids = gate_ids, footprint = footprint, beta = beta,
    projection = projection)
}

inner_selection <- rbindlist(inner_selection_parts)
aggregate_selection <- inner_selection[, .(
  mean_validation_mse = mean(validation_mse),
  mean_validation_pearson = mean(validation_pearson),
  mean_validation_sd_ratio = mean(validation_sd_ratio),
  total_validation_lines_positive = sum(n_validation_lines_positive)),
  by = lambda]
chosen_lambda <- aggregate_selection[order(mean_validation_mse, -lambda),
  lambda[1L]]
inner_selection[, selected := lambda == chosen_lambda]
aggregate_selection[, selected := lambda == chosen_lambda]
inner_selection_output <- rbindlist(list(inner_selection,
  aggregate_selection[, .(inner_fold = NA_integer_, fold = "aggregate",
    lambda, validation_mse = mean_validation_mse,
    validation_pearson = mean_validation_pearson,
    validation_observed_sd = NA_real_, validation_predicted_sd = NA_real_,
    validation_sd_ratio = mean_validation_sd_ratio,
    n_validation_lines_positive = total_validation_lines_positive,
    n_learning_lines = 26L, n_validation_lines = 0L,
    n_learning_rows = NA_integer_, n_validation_rows = NA_integer_,
    n_eligible_conditions = NA_integer_,
    representation = "fixed_direct_signed_O1_O3_modules",
    selection_objective = "mean_inner_validation_molecular_z_mse",
    selected)]), use.names = TRUE, fill = TRUE)
fwrite(inner_selection_output, file.path(RESULT,
  "susceptibility_model_inner_selection.tsv"), sep = "\t")
fwrite(inner_selection_output, file.path(RESULT, "model_selection.tsv"),
  sep = "\t")

## Compact frozen record only; the prior factorization branch is not executable.
previous_onmf <- data.table(
  method = "orthogonal_gene_selection_then_susceptibility_signature_oNMF",
  candidate_rank = 2:4,
  mean_inner_validation_mse = c(0.159589807092948, 0.157738899482000,
    0.159593676311728),
  mean_inner_validation_pearson = c(0.350273893159559,
    0.343281564096995, 0.328276594542723),
  mean_inner_validation_sd_ratio = c(0.360989554047971,
    0.292966294865143, 0.281505094036461),
  total_inner_validation_lines_positive = c(22L, 21L, 20L),
  outer_selected_signature_effective_rank = 2.97564117234999,
  outer_first_signature_axis_variance_fraction = 0.621416824534232,
  final_model_available = FALSE,
  preservation = "compact frozen comparator; not executed by Stage 16")
fwrite(previous_onmf, file.path(RESULT,
  "previous_onmf_susceptibility_experiments.tsv"), sep = "\t")

## Repeat all response-informed definitions using only the outer-training lines.
outer_gate <- eligibility[fold == "outer_1" &
  molecular_context_variable == TRUE, perturbation_id]
outer_footprint <- make_footprint(training_lines, score_table, D, programs)
outer_beta <- axis_vector(axis, "outer_training", programs)
outer_profiles <- fit_line_profiles(training_lines, "outer_training",
  outer_gate, outer_footprint, outer_beta, score_table, D, programs)
outer_gene_scores <- score_baseline_genes("outer_training", training_lines,
  outer_profiles$profile, ntc_expression, depmap_expression, plate_variance,
  universe_genes)
outer_weights <- attr(outer_gene_scores, "module_gene_weights")
outer_axis_definitions <- attr(outer_gene_scores,
  "orthogonal_axis_definition")[susceptibility_axis %chin% MODULE_AXES]
ntc_projection <- project_module_scores(ntc_expression, training_lines,
  all_lines, outer_weights, "outer_training", "tahoe_ntc", held_lines)
depmap_projection <- project_module_scores(depmap_expression, training_lines,
  all_lines, outer_weights, "outer_training", "depmap_static", held_lines)

profile_parts$outer_training <- outer_profiles$profile
profile_selection_parts$outer_training <- outer_profiles$selection
gene_score_parts$outer_training <- outer_gene_scores
module_weight_parts$outer_training <- outer_weights
axis_definition_parts$outer_training <- outer_axis_definitions
module_score_parts$outer_training <- ntc_projection$scores

profiles_all <- rbindlist(profile_parts)
profile_selection_all <- rbindlist(profile_selection_parts)
gene_scores_all <- rbindlist(gene_score_parts, use.names = TRUE, fill = TRUE)
weights_all <- rbindlist(module_weight_parts)
axis_definitions_all <- rbindlist(axis_definition_parts)
axis_definitions_all[, baseline_program := paste0(
  "baseline_susceptibility_", susceptibility_axis)]
setorder(axis_definitions_all, fold, susceptibility_axis, response_program)
fwrite(profiles_all, file.path(RESULT,
  "susceptibility_profile_by_training_line.tsv"), sep = "\t")
fwrite(profile_selection_all, file.path(RESULT,
  "susceptibility_profile_ridge_selection.tsv"), sep = "\t")
fwrite(gene_scores_all, file.path(ARTIFACT,
  "baseline_gene_susceptibility_scores.tsv.gz"), sep = "\t",
  compress = "gzip")
fwrite(axis_definitions_all, file.path(RESULT,
  "orthogonal_susceptibility_axis_definitions.tsv"), sep = "\t")
fwrite(weights_all, file.path(RESULT,
  "baseline_susceptibility_module_gene_weights.tsv.gz"), sep = "\t",
  compress = "gzip")

add_score_labels <- function(table) {
  table[, `:=`(cell_line = folds$cell_line[match(depmap_id, folds$depmap_id)],
    split_role = fifelse(depmap_id %chin% training_lines,
      "training", "held_out"))]
  table
}
ntc_scores <- add_score_labels(copy(ntc_projection$scores))
depmap_scores <- add_score_labels(copy(depmap_projection$scores))
fwrite(ntc_scores, file.path(RESULT,
  "baseline_susceptibility_module_scores_ntc.tsv"), sep = "\t")
fwrite(depmap_scores, file.path(RESULT,
  "baseline_susceptibility_module_scores_depmap.tsv"), sep = "\t")
fwrite(rbindlist(list(ntc_scores, depmap_scores)), file.path(RESULT,
  "susceptibility_program_scores.tsv.gz"), sep = "\t", compress = "gzip")

top_genes <- outer_weights[order(baseline_program, direction,
  -abs(normalized_signed_weight), gene), head(.SD, 25L),
  by = .(baseline_program, direction)]
top_genes[, rank_within_direction := seq_len(.N),
  by = .(baseline_program, direction)]
fwrite(top_genes, file.path(RESULT,
  "baseline_susceptibility_module_top_genes.tsv"), sep = "\t")

## Cross-fold response-axis, score, and disjoint-gene-set stability.
fold_labels <- names(module_weight_parts)
stability_parts <- list(); stability_id <- 0L
for (left in seq_along(fold_labels)) for (right in left:length(fold_labels)) {
  fold_left <- fold_labels[left]; fold_right <- fold_labels[right]
  for (axis_name in MODULE_AXES) {
    current_program <- paste0("baseline_susceptibility_", axis_name)
    loading_left <- axis_definitions_all[fold == fold_left &
      susceptibility_axis == axis_name][match(programs, response_program), loading]
    loading_right <- axis_definitions_all[fold == fold_right &
      susceptibility_axis == axis_name][match(programs, response_program), loading]
    score_left <- module_score_parts[[fold_left]][program == current_program][
      match(all_lines, depmap_id), standardized_program_score]
    score_right <- module_score_parts[[fold_right]][program == current_program][
      match(all_lines, depmap_id), standardized_program_score]
    genes_left <- module_weight_parts[[fold_left]][baseline_program ==
      current_program, gene]
    genes_right <- module_weight_parts[[fold_right]][baseline_program ==
      current_program, gene]
    overlap <- intersect(genes_left, genes_right)
    stability_id <- stability_id + 1L
    stability_parts[[stability_id]] <- data.table(
      fold_1 = fold_left, fold_2 = fold_right,
      baseline_program = current_program,
      susceptibility_axis = axis_name,
      response_configuration_cosine = sum(loading_left * loading_right) /
        sqrt(sum(loading_left ^ 2) * sum(loading_right ^ 2)),
      module_score_pearson = safe_cor(score_left, score_right),
      n_genes_fold_1 = length(genes_left), n_genes_fold_2 = length(genes_right),
      n_gene_overlap = length(overlap),
      gene_set_jaccard = length(overlap) / length(union(genes_left, genes_right)),
      correlation_method = "pearson")
  }
}
module_stability <- rbindlist(stability_parts)
fwrite(module_stability, file.path(RESULT,
  "baseline_susceptibility_module_stability.tsv"), sep = "\t")

score_cor_parts <- list(); score_cor_id <- 0L
source_matrices <- list(tahoe_ntc = ntc_projection$matrix,
  depmap_static = depmap_projection$matrix)
role_lines <- list(training = training_lines, held_out = held_lines,
  all = all_lines)
for (source_name in names(source_matrices)) for (role_name in names(role_lines)) {
  matrix_now <- source_matrices[[source_name]][role_lines[[role_name]],
    BASELINE_PROGRAMS, drop = FALSE]
  for (left in seq_len(length(BASELINE_PROGRAMS) - 1L))
    for (right in (left + 1L):length(BASELINE_PROGRAMS)) {
      p1 <- BASELINE_PROGRAMS[left]; p2 <- BASELINE_PROGRAMS[right]
      genes_1 <- outer_weights[baseline_program == p1, gene]
      genes_2 <- outer_weights[baseline_program == p2, gene]
      score_cor_id <- score_cor_id + 1L
      score_cor_parts[[score_cor_id]] <- data.table(
        baseline_source = source_name, role = role_name,
        program_1 = p1, program_2 = p2,
        n_lines = nrow(matrix_now), pearson = safe_cor(matrix_now[, p1],
          matrix_now[, p2]), n_gene_overlap = length(intersect(genes_1, genes_2)),
        gene_set_jaccard = length(intersect(genes_1, genes_2)) /
          length(union(genes_1, genes_2)), correlation_method = "pearson")
    }
}
module_correlations <- rbindlist(score_cor_parts)
fwrite(module_correlations, file.path(RESULT,
  "baseline_susceptibility_module_correlations.tsv"), sep = "\t")

## Compact pathway interpretation of the fixed outer modules.
pathways <- fread(file.path(WORK_ROOT, "data", "01_curated_pathways",
  "canonical_pathways_compendium.txt"))
pathway_sets <- strsplit(pathways$genes_str, ";", fixed = TRUE)
module_universe <- unique(outer_weights$gene)
pathway_parts <- list()
for (program in BASELINE_PROGRAMS) for (direction in c("positive", "negative")) {
  members <- outer_weights[baseline_program == program &
    get("direction") == direction, gene]
  current <- rbindlist(lapply(seq_len(nrow(pathways)), function(index) {
    pathway_genes <- intersect(pathway_sets[[index]], module_universe)
    overlap <- intersect(members, pathway_genes)
    if (length(pathway_genes) < 3L || !length(overlap)) return(NULL)
    data.table(baseline_program = program, direction = direction,
      pathway = pathways$name[index], pathway_source = pathways$source[index],
      pathway_size_in_module_universe = length(pathway_genes),
      module_direction_gene_count = length(members), overlap_count = length(overlap),
      overlap_genes = paste(sort(overlap), collapse = ";"),
      enrichment_p = phyper(length(overlap) - 1L, length(members),
        length(module_universe) - length(members), length(pathway_genes),
        lower.tail = FALSE))
  }), fill = TRUE)
  if (nrow(current)) {
    current[, enrichment_q := p.adjust(enrichment_p, "BH")]
    setorder(current, enrichment_q, enrichment_p, -overlap_count, pathway)
    pathway_parts[[paste(program, direction)]] <- current[1:min(.N, 10L)]
  }
}
pathway_annotations <- rbindlist(pathway_parts, fill = TRUE)
fwrite(pathway_annotations, file.path(RESULT,
  "baseline_susceptibility_module_pathways.tsv"), sep = "\t")

## One outer global molecular conditioner; only lambda was selected above.
outer_train_rows <- score_table$depmap_id %chin% training_lines &
  score_table$perturbation_id %chin% outer_gate
outer_held_rows <- score_table$depmap_id %chin% held_lines
R_train <- outer_footprint[score_table$perturbation_id[outer_train_rows],
  programs, drop = FALSE]
z_train <- as.numeric((D[outer_train_rows, programs, drop = FALSE] -
  R_train) %*% outer_beta)
X_train <- make_interactions(score_table$depmap_id[outer_train_rows],
  score_table$perturbation_id[outer_train_rows], ntc_projection$matrix,
  outer_footprint, BASELINE_PROGRAMS, programs)
susceptibility_fit <- fit_ridge(X_train, z_train, chosen_lambda)
raw_coefficients <- raw_ridge_coefficients(susceptibility_fit)
feature_map <- tstrsplit(colnames(X_train), "__", fixed = TRUE)
coefficient_table <- data.table(baseline_program = feature_map[[1L]],
  response_program = feature_map[[2L]], interaction = colnames(X_train),
  coefficient = raw_coefficients$beta[colnames(X_train)],
  intercept = raw_coefficients$intercept, standardized_coefficient = 0,
  ridge_penalty = chosen_lambda, n_training_lines = length(training_lines),
  n_training_rows = sum(outer_train_rows), baseline_source_fit = "tahoe_ntc")
coefficient_table$standardized_coefficient[susceptibility_fit$keep] <-
  as.numeric(susceptibility_fit$beta_standardized)
fwrite(coefficient_table, file.path(RESULT,
  "susceptibility_model_coefficients.tsv"), sep = "\t")

make_held_predictions <- function(baseline, source_name) {
  held_index <- which(outer_held_rows)
  R_held <- outer_footprint[score_table$perturbation_id[held_index],
    programs, drop = FALSE]
  observed_context <- D[held_index, programs, drop = FALSE] - R_held
  observed_z <- as.numeric(observed_context %*% outer_beta)
  X_held <- make_interactions(score_table$depmap_id[held_index],
    score_table$perturbation_id[held_index], baseline, outer_footprint,
    BASELINE_PROGRAMS, programs)
  eligible <- score_table$perturbation_id[held_index] %chin% outer_gate
  ungated <- predict_ridge(susceptibility_fit, X_held)
  predicted_z <- fifelse(eligible, ungated, 0)
  predicted_context <- outer(predicted_z / sum(outer_beta ^ 2), outer_beta)
  Dhat <- R_held + predicted_context
  result <- data.table(observation_id = score_table$observation_id[held_index],
    depmap_id = score_table$depmap_id[held_index],
    cell_line = score_table$cell_line[held_index],
    perturbation_id = score_table$perturbation_id[held_index],
    drug = score_table$drug[held_index], dose = score_table$dose[held_index],
    dose_panel = score_table$dose_panel[held_index],
    baseline_source = source_name, molecular_eligible = eligible,
    observed_z = observed_z, predicted_z = predicted_z,
    ungated_predicted_z = ungated,
    reconstruction_projection_error = max(abs(as.numeric(
      (Dhat - R_held) %*% outer_beta) - predicted_z)))
  for (program in programs) {
    set(result, j = paste0("reference_", program), value = R_held[, program])
    set(result, j = paste0("observed_context_", program),
      value = observed_context[, program])
    set(result, j = paste0("reconstructed_displacement_", program),
      value = Dhat[, program])
  }
  result
}
held_predictions <- rbindlist(list(make_held_predictions(
  ntc_projection$matrix, "tahoe_ntc"), make_held_predictions(
  depmap_projection$matrix, "depmap_static")))
activity <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
held_predictions[, activity_supported := activity$activity_supported[match(
  perturbation_id, activity$perturbation_id)]]
fwrite(held_predictions, file.path(RESULT,
  "held_susceptibility_molecular_predictions.tsv.gz"), sep = "\t",
  compress = "gzip")

## Stage-owned direct-axis diagnostics.
configuration_matrix <- dcast(outer_axis_definitions,
  susceptibility_axis ~ response_program, value.var = "loading")
configuration_values <- as.matrix(configuration_matrix[, ..programs])
rownames(configuration_values) <- paste0("baseline_susceptibility_",
  configuration_matrix$susceptibility_axis)
configuration_values <- configuration_values[BASELINE_PROGRAMS, programs,
  drop = FALSE]
configuration_limit <- max(abs(configuration_values))
save_heatmap("16a_O1_O3_response_configuration_signatures.pdf",
  draw_heatmap(configuration_values, scale_rows = FALSE,
    value_range = c(-configuration_limit, 0, configuration_limit),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Direct susceptibility-axis response configurations"), 7.0, 3.8)

top_order <- copy(top_genes)
top_order[, direction_order := fifelse(direction == "positive", 1L, 2L)]
setorder(top_order, baseline_program, direction_order, rank_within_direction)
top_matrix <- dcast(outer_weights[gene %chin% top_order$gene],
  gene ~ baseline_program, value.var = "normalized_signed_weight", fill = 0)
top_values <- as.matrix(top_matrix[match(top_order$gene, gene),
  ..BASELINE_PROGRAMS]); rownames(top_values) <- top_order$gene
if (!identical(rownames(top_values), top_order$gene))
  fail("Top-gene heatmap annotation order is not synchronized.")
top_limit <- max(abs(top_values))
save_heatmap("16b_positive_negative_top_genes_by_module.pdf",
  draw_heatmap(top_values,
    row_anno_df = data.frame(Program = top_order$baseline_program,
      Direction = top_order$direction),
    row_anno_cols = c("Program", "Direction"), scale_rows = FALSE,
    value_range = c(-top_limit, 0, top_limit),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 4.8),
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Positive and negative genes in direct signed modules"),
  7.2, min(9.5, 3.5 + 0.07 * nrow(top_values)))

if (nrow(pathway_annotations)) {
  pathway_plot <- copy(pathway_annotations)
  pathway_plot[, module_direction := paste(
    sub("baseline_susceptibility_", "", baseline_program), direction,
    sep = "_")]
  pathway_plot[, pathway_label := make.unique(pathway)]
  pathway_wide <- dcast(pathway_plot, pathway_label ~ module_direction,
    value.var = "enrichment_q", fill = 1)
  pathway_columns <- setdiff(names(pathway_wide), "pathway_label")
  pathway_values <- -log10(pmax(as.matrix(pathway_wide[, ..pathway_columns]),
    .Machine$double.xmin)); rownames(pathway_values) <- pathway_wide$pathway_label
  save_heatmap("16c_pathway_enrichment_by_module.pdf",
    draw_heatmap(pathway_values, scale_rows = FALSE,
      value_range = c(0, max(pathway_values) / 2, max(pathway_values)),
      value_colors = c("white", "#67A9CF", "#B2182B"),
      show_row_names = TRUE, show_column_names = TRUE,
      row_names_gp = grid::gpar(fontsize = 4.5),
      cluster_rows = FALSE, cluster_columns = FALSE,
      title = "Direct-module pathway enrichment (-log10 FDR)"),
    8.0, min(9.5, 3.5 + 0.13 * nrow(pathway_values)))
} else {
  save_gg("16c_pathway_enrichment_by_module.pdf",
    ggplot() + annotate("text", 0, 0, label = "No pathway overlaps") +
      xlim(-1, 1) + ylim(-1, 1) + theme_void(), 5.0, 3.5)
}

score_correlation_matrix <- cor(ntc_projection$matrix[training_lines,
  BASELINE_PROGRAMS, drop = FALSE], method = "pearson")
save_heatmap("16d_baseline_module_score_correlations.pdf",
  draw_heatmap(score_correlation_matrix, scale_rows = FALSE,
    value_range = c(-1, 0, 1),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Outer-training baseline-module score correlations"), 5.2, 4.8)

inner_selected <- inner_selection[selected == TRUE]
p_inner <- ggplot(inner_selected,
    aes(factor(inner_fold), validation_pearson)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_col(fill = "#2166AC", width = 0.7) +
  labs(x = "Inner fold", y = "Molecular susceptibility Pearson r",
    title = "Inner-fold validation at selected ridge penalty") +
  theme_stage16()
save_gg("16e_inner_fold_molecular_validation.pdf", p_inner, 5.4, 4.0)

plot_held <- held_predictions[baseline_source == "tahoe_ntc" &
  molecular_eligible == TRUE]
p_z <- ggplot(plot_held, aes(observed_z, predicted_z)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.4, size = 0.9, color = "#2166AC") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.65, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("r=%s; SD ratio=%s; n=%d",
      fmt(safe_cor(plot_held$observed_z, plot_held$predicted_z)),
      fmt(sd(plot_held$predicted_z) / sd(plot_held$observed_z)),
      nrow(plot_held)), size = 2.35) +
  labs(x = "Observed molecular susceptibility",
    y = "Predicted molecular susceptibility",
    title = "Held-out direct-axis molecular prediction") +
  square_scatter_theme()
save_gg("16f_held_observed_vs_predicted_molecular_z.pdf", p_z, 5.0, 4.8)

per_line <- plot_held[, .(n_rows = .N,
  pearson = safe_cor(observed_z, predicted_z),
  predicted_to_observed_sd_ratio = sd(predicted_z) / sd(observed_z)),
  by = .(depmap_id, cell_line)]
p_line <- ggplot(per_line, aes(reorder(cell_line, pearson), pearson)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_segment(aes(xend = cell_line, y = 0, yend = pearson),
    color = "grey65", linewidth = 0.55) +
  geom_point(size = 2.1, color = "#2166AC") + coord_flip() +
  labs(x = NULL, y = "Molecular susceptibility Pearson r",
    title = "Held-out performance by cell line") + theme_stage16()
save_gg("16g_molecular_performance_by_held_line.pdf", p_line, 5.8, 4.5)

per_perturbation <- plot_held[, .(n_rows = .N,
  pearson = safe_cor(observed_z, predicted_z), drug = first(drug),
  dose = first(dose)), by = perturbation_id]
plot_perturbation <- per_perturbation[is.finite(pearson)][order(pearson)]
plot_perturbation[, order_id := seq_len(.N)]
p_perturbation <- ggplot(plot_perturbation, aes(order_id, pearson)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey55") +
  geom_segment(aes(xend = order_id, y = 0, yend = pearson),
    linewidth = 0.25, color = "grey65") +
  geom_point(size = 1.1, color = "#B2182B") +
  labs(x = "Perturbations sorted by performance",
    y = "Held-line Pearson r",
    title = "Molecular performance by perturbation") + theme_stage16()
save_gg("16h_molecular_performance_by_perturbation.pdf", p_perturbation,
  7.0, 4.2)

score_concordance <- merge(ntc_scores[, .(depmap_id, cell_line, split_role,
  program, tahoe_ntc = standardized_program_score)],
  depmap_scores[, .(depmap_id, program,
    depmap_static = standardized_program_score)],
  by = c("depmap_id", "program"))
p_concordance <- ggplot(score_concordance,
    aes(tahoe_ntc, depmap_static, color = program, shape = split_role)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(alpha = 0.8, size = 1.15) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55, color = "black") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
    label = sprintf("pooled r=%s", fmt(safe_cor(
      score_concordance$tahoe_ntc, score_concordance$depmap_static))),
    size = 2.35) +
  labs(x = "Tahoe NTC module score", y = "DepMap module score",
    color = "Module", shape = NULL,
    title = "Fixed NTC modules projected into DepMap") +
  square_scatter_theme()
save_gg("16i_ntc_vs_depmap_module_scores.pdf", p_concordance, 5.4, 5.0)

coefficient_wide <- dcast(coefficient_table,
  baseline_program ~ response_program, value.var = "coefficient")
coefficient_values <- as.matrix(coefficient_wide[, ..programs]); rownames(
  coefficient_values) <- coefficient_wide$baseline_program
coefficient_values <- coefficient_values[BASELINE_PROGRAMS, programs,
  drop = FALSE]
coefficient_limit <- max(abs(coefficient_values))
save_heatmap("16j_baseline_by_response_coefficient_heatmap.pdf",
  draw_heatmap(coefficient_values, scale_rows = FALSE,
    value_range = c(-coefficient_limit, 0, coefficient_limit),
    value_colors = c("#2166AC", "white", "#B2182B"),
    show_row_names = TRUE, show_column_names = TRUE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    title = "Baseline-module by reference-response coefficients"), 7.0, 3.8)

pdfs <- c("16a_O1_O3_response_configuration_signatures.pdf",
  "16b_positive_negative_top_genes_by_module.pdf",
  "16c_pathway_enrichment_by_module.pdf",
  "16d_baseline_module_score_correlations.pdf",
  "16e_inner_fold_molecular_validation.pdf",
  "16f_held_observed_vs_predicted_molecular_z.pdf",
  "16g_molecular_performance_by_held_line.pdf",
  "16h_molecular_performance_by_perturbation.pdf",
  "16i_ntc_vs_depmap_module_scores.pdf",
  "16j_baseline_by_response_coefficient_heatmap.pdf")
pdf_page_count <- function(path) {
  content <- readLines(path, warn = FALSE, encoding = "bytes")
  sum(grepl("/Type /Page([^s]|$)", content, perl = TRUE))
}
figures <- data.table(pdf = pdfs,
  producing_script = "scripts/16_learn_baseline_susceptibility_programs.R")
figures[, `:=`(exists = file.exists(file.path(RESULT, pdf)),
  pages = vapply(file.path(RESULT, pdf), pdf_page_count, integer(1L)))]
figures[, passed := exists & pages == 1L]
fwrite(figures, file.path(RESULT, "FIGURE_INDEX.tsv"), sep = "\t")

ntc_eligible <- held_predictions[baseline_source == "tahoe_ntc" &
  molecular_eligible == TRUE]
depmap_eligible <- held_predictions[baseline_source == "depmap_static" &
  molecular_eligible == TRUE]
outer_counts <- outer_weights[, .N, by = baseline_program]
stage16_summary <- data.table(
  representation = "direct_signed_O1_O3_susceptibility_modules",
  selected_susceptibility_ridge_penalty = chosen_lambda,
  n_baseline_programs = length(BASELINE_PROGRAMS),
  n_outer_module_genes = nrow(outer_weights),
  n_outer_O1_genes = outer_counts[baseline_program == BASELINE_PROGRAMS[1L], N],
  n_outer_O2_genes = outer_counts[baseline_program == BASELINE_PROGRAMS[2L], N],
  n_outer_O3_genes = outer_counts[baseline_program == BASELINE_PROGRAMS[3L], N],
  n_outer_eligible_genes = outer_gene_scores[susceptibility_eligible == TRUE, .N],
  n_outer_eligible_conditions = length(outer_gate),
  inner_selected_mean_mse = aggregate_selection[selected == TRUE,
    mean_validation_mse],
  inner_selected_mean_pearson = aggregate_selection[selected == TRUE,
    mean_validation_pearson],
  ntc_held_eligible_molecular_z_pearson = safe_cor(ntc_eligible$observed_z,
    ntc_eligible$predicted_z),
  ntc_held_eligible_molecular_z_sd_ratio = sd(ntc_eligible$predicted_z) /
    sd(ntc_eligible$observed_z),
  depmap_held_eligible_molecular_z_pearson = safe_cor(
    depmap_eligible$observed_z, depmap_eligible$predicted_z),
  depmap_held_eligible_molecular_z_sd_ratio = sd(depmap_eligible$predicted_z) /
    sd(depmap_eligible$observed_z),
  n_ntc_positive_held_lines = sum(per_line$pearson > 0, na.rm = TRUE),
  response_target = "function-informed molecular z from frozen Stage 05 coefficients",
  prism_outcome_table_loaded = FALSE, depmap_used_for_selection = FALSE)
fwrite(stage16_summary, file.path(RESULT, "stage16_summary.tsv"), sep = "\t")

failure_modes <- data.table(failure_mode = c(
  "empty_direct_axis_module", "invalid_module_variance_or_scaling",
  "weak_held_molecular_susceptibility", "insufficient_prediction_amplitude",
  "limited_positive_held_line_breadth"),
  triggered = c(FALSE, FALSE,
    stage16_summary$ntc_held_eligible_molecular_z_pearson <= 0.333,
    stage16_summary$ntc_held_eligible_molecular_z_sd_ratio < 0.25,
    stage16_summary$n_ntc_positive_held_lines < 5L),
  detail = c("all O1-O3 modules contain eligible genes",
    "training-only NTC and DepMap gene/module scaling succeeded",
    sprintf("NTC eligible r=%s; prespecified target >0.333",
      fmt(stage16_summary$ntc_held_eligible_molecular_z_pearson)),
    sprintf("NTC eligible SD ratio=%s; prespecified target >=0.25",
      fmt(stage16_summary$ntc_held_eligible_molecular_z_sd_ratio)),
    sprintf("%d of seven held lines have positive r; target >=5",
      stage16_summary$n_ntc_positive_held_lines)))
fwrite(failure_modes, file.path(RESULT, "failure_mode_summary.tsv"), sep = "\t")

saveRDS(list(label = "direct_signed_O1_O3_baseline_susceptibility_modules",
  selected_lambda = chosen_lambda, baseline_programs = BASELINE_PROGRAMS,
  response_programs = programs, module_gene_weights = outer_weights,
  response_axis_definitions = outer_axis_definitions,
  susceptibility_fit = susceptibility_fit, raw_coefficients = raw_coefficients,
  functional_molecular_axis = outer_beta, molecular_eligibility = outer_gate,
  ntc_gene_center = ntc_projection$gene_center,
  ntc_gene_scale = ntc_projection$gene_scale,
  ntc_module_center = ntc_projection$module_center,
  ntc_module_scale = ntc_projection$module_scale,
  depmap_gene_center = depmap_projection$gene_center,
  depmap_gene_scale = depmap_projection$gene_scale,
  depmap_module_center = depmap_projection$module_center,
  depmap_module_scale = depmap_projection$module_scale,
  held_out_used_for_selection_or_fitting = FALSE,
  prism_outcomes_loaded = FALSE), file.path(ARTIFACT,
    "baseline_susceptibility_program_model.rds"))

projection_error <- max(held_predictions$reconstruction_projection_error)
axis_norms <- axis_definitions_all[, sqrt(sum(loading ^ 2)),
  by = .(fold, susceptibility_axis)]
axis_orientations <- axis_definitions_all[order(response_program),
  loading[which.max(abs(loading))] > 0,
  by = .(fold, susceptibility_axis)]
module_norms <- weights_all[, sum(normalized_signed_weight ^ 2),
  by = .(fold, baseline_program)]
within_fold_overlap <- weights_all[, {
  sets <- split(gene, baseline_program)
  overlaps <- combn(names(sets), 2L, function(pair)
    length(intersect(sets[[pair[1L]]], sets[[pair[2L]]])))
  .(maximum_overlap = max(overlaps))
}, by = fold]
checks <- data.table(check = c("fixed_26_7_split_and_exclusions",
  "exactly_nine_response_programs", "only_O1_O3_baseline_programs",
  "no_onmf_executed", "no_prism_outcome_table_loaded",
  "functional_guidance_only_frozen_molecular_axis",
  "validation_and_held_excluded_from_definitions",
  "deterministic_unit_orthogonal_axis_orientation",
  "signed_unit_norm_gene_modules", "gene_sets_disjoint",
  "no_identity_features", "one_global_molecular_conditioner",
  "context_ineligible_exactly_zero", "reconstruction_exact_molecular_z",
  "depmap_fixed_ntc_definitions", "pearson_correlations_only",
  "all_diagnostic_pdfs_one_page"),
  value = c(length(training_lines) == 26L && length(held_lines) == 7L &&
      !any(EXCLUDED_IDS %chin% all_lines), length(programs) == 9L,
    setequal(unique(weights_all$baseline_program), BASELINE_PROGRAMS) &&
      setequal(unique(axis_definitions_all$susceptibility_axis), MODULE_AXES),
    TRUE, TRUE, TRUE,
    all(gene_scores_all[fold != "outer_training", n_learning_lines] <= 21L) &&
      outer_gene_scores$n_learning_lines[1L] == 26L,
    max(abs(axis_norms$V1 - 1)) <= 1e-10 && all(axis_orientations$V1),
    max(abs(module_norms$V1 - 1)) <= 1e-10 &&
      all(weights_all$direction %chin% c("positive", "negative")),
    all(within_fold_overlap$maximum_overlap == 0L),
    !any(grepl("drug|perturbation|cell_line", coefficient_table$interaction)),
    TRUE, all(held_predictions[molecular_eligible == FALSE, predicted_z] == 0),
    projection_error <= 1e-10,
    identical(dim(ntc_projection$weights), dim(depmap_projection$weights)) &&
      max(abs(ntc_projection$weights - depmap_projection$weights)) == 0,
    all(c(module_correlations$correlation_method,
      module_stability$correlation_method) == "pearson"),
    all(figures$passed)),
  detail = c("accepted fixed cohort", paste(programs, collapse = ","),
    paste(BASELINE_PROGRAMS, collapse = ","),
    "direct signed modules; no factorization call",
    "Stage 16 loads no PRISM outcome table",
    "Stage 05 coefficients enter only the molecular susceptibility target",
    "all gene, axis, weight, scaling, and ridge fitting are fold-learning-only",
    "unit response configurations; largest absolute response loading positive",
    "signed minimum absolute LOO correlation; unit L2 per module",
    "dominant-axis assignment produces disjoint sets; O4 is not reassigned",
    "only baseline-module by reference-response interactions",
    "one ridge fit shared across lines and perturbations",
    "reference footprint retained for ineligible conditions",
    sprintf("maximum reconstruction error %.3g", projection_error),
    "DepMap changes only training-line platform center and scale",
    "every reported correlation explicitly uses Pearson",
    "one page per Stage 16 PDF"))
fwrite(checks, file.path(RESULT, "implementation_checks.tsv"), sep = "\t")
fwrite(checks, file.path(RESULT, "stage16_required_checks.tsv"), sep = "\t")
if (!all(checks$value)) fail("A required Stage 16 check failed.")
write_artifact_manifest(ARTIFACT,
  "scripts/16_learn_baseline_susceptibility_programs.R")
writeLines("complete", file.path(RESULT, "stage16_complete.flag"))
message(sprintf(paste0("Stage 16 complete: direct O1-O3 modules; lambda %.6g; ",
  "NTC held eligible r=%s; SD ratio=%s."), chosen_lambda,
  fmt(stage16_summary$ntc_held_eligible_molecular_z_pearson),
  fmt(stage16_summary$ntc_held_eligible_molecular_z_sd_ratio)))
