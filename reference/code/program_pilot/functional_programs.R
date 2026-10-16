# Shared regression, clustering, and paired plotting for stages 11--15.
# The panel's aliases are explicit; no salt or synonym inference is performed.
panel_name <- function(x) {
  out <- canonical_drug(x)
  out[out == "rapamycin"] <- "sirolimus"
  out
}

# The mixed-correlation experiment changes the NMF representation only.
# Both correlation matrices use the existing per-view observation weights.
# The same signed two-channel encoding retains negative correlations.
functional_nmf_input <- function(support, displacement_correlation, include_channels = TRUE) {
  loading <- support$loadings
  correlation <- gene_correlation(t(loading), rep(1, ncol(loading)))
  representation <- support$settings$clustering_representation
  if (!is.null(representation)) {
    check(identical(representation, "mean_functional_displacement_correlation"),
      "Unsupported functional clustering representation")
    check(identical(dimnames(correlation), dimnames(displacement_correlation)),
      "Functional and displacement correlation genes differ")
    correlation <- (correlation + displacement_correlation) / 2
  }
  if (!include_channels) return(list(correlation = correlation, channels = NULL))
  if (is.null(representation)) {
    channels <- cbind(pmax(loading, 0), pmax(-loading, 0))
    labels <- colnames(loading)
  } else {
    channels <- cbind(pmax(correlation, 0), pmax(-correlation, 0))
    labels <- rownames(correlation)
  }
  colnames(channels) <- c(paste0("positive:", labels), paste0("negative:", labels))
  list(correlation = correlation, channels = channels)
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

# Reused by discovery and cached-bank validation. Select the maximum positive
# mean gene silhouette among cuts with every seed at least local_minimum_genes
# large. Exact ties select the smaller rank; single-fit NMF has no resampling.
select_silhouette <- function(rank_table, rank_override = NULL, selection_rule = "maximum_silhouette") {
  check(selection_rule %in% c("maximum_silhouette", "last_blue_dot", "silhouette_elbow"),
    "Unsupported rank selection rule")
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
  out[, silhouette_candidate := rank_allowed & minimum_size_qc_pass & is.finite(mean_silhouette) & mean_silhouette > 0]
  candidates <- out[silhouette_candidate == TRUE]
  candidates <- if (selection_rule == "last_blue_dot") candidates[order(-rank)] else
    candidates[order(-mean_silhouette, rank)]
  rule_label <- if (selection_rule == "last_blue_dot")
    "largest rank with positive mean gene silhouette and minimum seed size (last blue dot)" else
    "maximum positive mean gene silhouette with minimum seed size; smaller rank breaks exact ties"
  selected_rank <- if (nrow(candidates)) candidates$rank[1L] else NA_integer_
  status <- if (nrow(candidates)) "selected" else
    "unresolved: no allowed size-qualified cut with positive mean gene silhouette"
  if (selection_rule == "silhouette_elbow") {
    rule_label <- "silhouette elbow: best continuous two-segment fit after the eligible peak; descending then flatter; minimum seed size"
    selected_rank <- NA_integer_
    status <- "unresolved: no size-qualified descending-to-flatter silhouette elbow"
    out[, `:=`(elbow_sse = NA_real_, elbow_slope_before = NA_real_, elbow_slope_after = NA_real_,
      elbow_fitted_silhouette = NA_real_, elbow_peak_rank = NA_integer_)]
    if (nrow(candidates)) {
      # Fit the complete descending/tail curve, including size-QC-failing cuts,
      # but permit selection only at a positive, size-qualified observed cut.
      peak <- candidates[order(-mean_silhouette, rank), rank][1L]
      curve <- out[is.finite(mean_silhouette) & rank >= peak]
      out[, elbow_peak_rank := peak]
      for (break_rank in candidates$rank) {
        if (sum(curve$rank <= break_rank) < 3L || sum(curve$rank > break_rank) < 3L) next
        fit_data <- copy(curve)
        fit_data[, hinge := pmax(rank - break_rank, 0)]
        fit <- lm(mean_silhouette ~ rank + hinge, data = fit_data)
        before <- unname(coef(fit)["rank"])
        after <- before + unname(coef(fit)["hinge"])
        if (!is.finite(before) || !is.finite(after) || before >= 0 || after <= before) next
        out[rank == break_rank, `:=`(elbow_sse = sum(residuals(fit)^2),
          elbow_slope_before = before, elbow_slope_after = after)]
      }
      knees <- out[is.finite(elbow_sse)][order(elbow_sse, rank)]
      if (nrow(knees)) {
        selected_rank <- knees$rank[1L]
        status <- "selected"
        curve[, hinge := pmax(rank - selected_rank, 0)]
        fit <- lm(mean_silhouette ~ rank + hinge, data = curve)
        out[match(curve$rank, rank), elbow_fitted_silhouette := as.numeric(fitted(fit))]
      }
    }
  }
  out[, `:=`(selected = !is.na(selected_rank) & rank == selected_rank,
    selection_status = status,
    selection_rule = rule_label,
    cophenetic_change_to_next_rank = c(diff(consensus_cophenetic), NA_real_))]
  if (!is.null(rank_override)) {
    check(length(rank_override) == 1L && rank_override %in% out$rank, 'Requested rank is not a completed fit')
    out[, `:=`(selected = rank == rank_override, selection_status = 'selected: explicit user rank',
      selection_rule = paste('User-selected rank', rank_override))]
  }
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

# Reused by fitting and cache validation. Learn one sensitivity regression per
# gene and assay; score held-out drugs using training-only centering/scaling.
# Loading = (baseline squared error - gene squared error) / training outcome
# variance. The baseline is the training outcome mean. Unmeasured outcomes
# never become fitting columns. Zero training variance cannot estimate a slope
# and explicitly returns the baseline, with status recorded for that gene/fold.
functional_support <- function(x, metadata, outcomes, allow_zero = FALSE, sensitivity_sign = FALSE,
  unit_column = "drug") {
  check(identical(rownames(x), metadata$observation_id) &&
    identical(outcomes$observation_id, metadata$observation_id), "Functional-support row alignment mismatch")
  loadings <- list()
  observations <- list()
  coefficients <- list()
  folds <- list()
  assay_coverage <- list()
  for (assay_name in c("PRISM", "ALMANAC")) {
    outcome_column <- if (assay_name == "PRISM") "prism_abundance_pct" else "almanac_abundance_pct"
    matched <- which(!is.na(outcomes[[outcome_column]]))
    md <- metadata[matched]
    y <- outcomes[[outcome_column]][matched]
    drugs <- sort(unique(md[[unit_column]]))
    eligible <- length(matched) >= functional_folds && length(drugs) >= functional_folds
    assay_coverage[[assay_name]] <- data.table(assay = assay_name, conditions = length(matched),
      drugs = length(drugs), included = eligible,
      reason = if (eligible) "included" else "fewer than five measured drugs for whole-drug folds")
    if (!eligible) next
    check(length(drugs) >= functional_folds && all(is.finite(y)), "Insufficient drugs or invalid functional outcome")
    fold <- 1L + (match(md[[unit_column]], drugs) - 1L) %% functional_folds
    predicted <- matrix(NA_real_, nrow = length(matched), ncol = ncol(x))
    baseline <- normalization <- rep(NA_real_, length(matched))
    for (fold_id in seq_len(functional_folds)) {
      train <- which(fold != fold_id)
      test <- which(fold == fold_id)
      check(!any(md[[unit_column]][test] %in% md[[unit_column]][train]), "Group leakage in functional folds")
      y_mean <- mean(y[train])
      y_sd <- sqrt(mean((y[train] - y_mean)^2))
      check(is.finite(y_sd) && y_sd > 0, "Constant outcome in a functional training fold")
      baseline[test] <- y_mean
      normalization[test] <- y_sd^2
      x_train <- x[matched[train], , drop = FALSE]
      x_mean <- colMeans(x_train)
      x_sd <- matrixStats::colSds(x_train)
      y_standardized <- (y[train] - y_mean) / y_sd
      slopes <- numeric(ncol(x))
      status <- rep("ridge fit", ncol(x))
      for (j in seq_len(ncol(x))) {
        if (x_sd[j] == 0) {
          predicted[test, j] <- y_mean
          status[j] <- "constant training gene: baseline prediction"
          next
        }
        z <- (x_train[, j] - x_mean[j]) / x_sd[j]
        # glmnet requires two design columns. The constant column is exactly
        # zero and excluded; the only active predictor is this gene.
        fit <- glmnet::glmnet(cbind(gene = z, constant = 0), y_standardized,
          alpha = 0, lambda = functional_ridge_lambda, standardize = FALSE, intercept = FALSE,
          thresh = 1e-10, exclude = 2L)
        slopes[j] <- as.matrix(coef(fit))["gene", 1L]
        predicted[test, j] <- y_mean + y_sd * slopes[j] *
          (x[matched[test], j] - x_mean[j]) / x_sd[j]
      }
      coefficients[[paste(assay_name, fold_id)]] <- data.table(assay = assay_name, fold = fold_id,
        gene = colnames(x), training_conditions = length(train), training_drugs = uniqueN(md$drug[train]),
        training_gene_mean = x_mean, training_gene_sd = x_sd, training_outcome_mean = y_mean,
        training_outcome_sd = y_sd, standardized_beta = slopes, ridge_lambda = functional_ridge_lambda,
        fit_status = status)
    }
    check(all(is.finite(predicted)) && all(is.finite(normalization)) && all(normalization > 0),
      "Nonfinite functional predictions or scaling")
    residual <- sweep(predicted, 1L, y, "-")
    gain <- sweep(sweep(-residual^2, 1L, (y - baseline)^2, "+"), 1L, normalization, "/")
    loading <- t(gain)
    colnames(loading) <- paste0(assay_name, ":", md$observation_id)
    rownames(loading) <- colnames(x)
    loadings[[assay_name]] <- loading
    observations[[assay_name]] <- data.table(assay = assay_name, observation_id = md$observation_id,
      drug = md$drug, dose = md$dose, fold = fold, outcome = y, baseline_prediction = baseline,
      training_outcome_variance = normalization)
    folds[[assay_name]] <- data.table(assay = assay_name, observation_id = md$observation_id, drug = md$drug, fold = fold)
    if (unit_column != "drug") {
      observations[[assay_name]][, regression_unit := md[[unit_column]]]
      folds[[assay_name]][, regression_unit := md[[unit_column]]]
    }
  }
  loading <- do.call(cbind, loadings)
  check(length(loadings) > 0L, "Neither assay has enough drugs for functional regression")
  observation_table <- rbindlist(observations)
  raw_prediction_gain <- NULL
  if (sensitivity_sign) {
    # Multiply the complete existing gain by the measured response sign.
    # Negative gains are preserved; no value is clipped or discarded.
    check(identical(colnames(loading), paste0(observation_table$assay, ":", observation_table$observation_id)),
      "Response signs do not align with loading columns")
    raw_prediction_gain <- loading
    observation_table[, response_sign := sign(50 - outcome)]
    loading <- sweep(loading, 2L, observation_table$response_sign, "*")
  }
  check(all(is.finite(loading)) && (allow_zero || all(rowSums(loading^2) > 0)),
    "Invalid or entirely zero functional gene profile")
  channels <- cbind(pmax(loading, 0), pmax(-loading, 0))
  colnames(channels) <- c(paste0("positive:", colnames(loading)), paste0("negative:", colnames(loading)))
  result <- list(loadings = loading, channels = channels, observations = observation_table,
    coefficients = rbindlist(coefficients), folds = rbindlist(folds), assay_coverage = rbindlist(assay_coverage),
    settings = list(folds = functional_folds, ridge_lambda = functional_ridge_lambda,
      glmnet_version = as.character(packageVersion("glmnet")), normalization = "training outcome variance"))
  if (sensitivity_sign) {
    result$raw_prediction_gain <- raw_prediction_gain
    result$settings$sensitivity_sign <- "sign(50 - measured abundance percent vehicle); full gain retained"
  }
  if (unit_column != "drug") {
    result$settings$unit_column <- unit_column
    setnames(result$assay_coverage, "drugs", "regression_units")
    result$assay_coverage[included == FALSE, reason := "fewer than five measured cell lines for whole-cell folds"]
    result$coefficients[, training_units := {
      assay_value <- .BY$assay
      fold_value <- .BY$fold
      uniqueN(result$folds[assay == assay_value & fold != fold_value, regression_unit])
    }, by = .(assay, fold)]
  }
  result
}

# A top-level worker receives only this explicitly exported numerical context.
# It never serializes the caller's expression matrices or complete model bank.
# Commit each model before reporting completion or returning through the socket.
fit_functional_nmf_rank <- function(rank) {
  context <- nmf_rank_context
  path <- file.path(context$checkpoint_dir, sprintf('rank_%03d.qs', rank))
  if (file.exists(path)) {
    result <- qs::qread(path)
    stopifnot(identical(result$fit_signature, context$fit_signature), result$diagnostics$rank == rank)
    result$model <- NULL
    result$checkpoint_path <- path
    message(context$scope, '; reused saved NMF k=', rank)
    return(result)
  }
  fit_seed <- context$seed_base + rank
  model <- RcppML::nmf(context$channels, k=rank, seed=fit_seed, tol=context$tolerance,
    maxit=context$max_iterations, threads=1L, resource='cpu', verbose=FALSE)
  stopifnot(all(is.finite(model@w)), all(model@w >= 0), all(is.finite(model@h)),
    all(model@h >= 0), all(is.finite(model@d)), all(model@d >= 0))
  raw_partition <- max.col(model@w, ties.method='first')
  factor_order <- c(unique(raw_partition[context$tree_order]), setdiff(seq_len(rank), raw_partition))
  partition <- setNames(match(raw_partition, factor_order), context$genes)
  sizes <- tabulate(partition, nbins=rank)
  mean_silhouette <- if (data.table::uniqueN(partition) >= 2L)
    mean(cluster::silhouette(partition, context$distance)[, 'sil_width']) else NA_real_
  rank_row <- data.table::data.table(scope=context$scope, orientation=context$orientation, rank=rank,
    minimum_size=min(sizes), maximum_size=max(sizes), mean_silhouette=mean_silhouette,
    minimum_size_qc_pass=min(sizes) >= context$minimum_genes,
    nmf_fits=1L, nmf_seed=fit_seed, nmf_iterations=model@misc$iter,
    nmf_final_tolerance=model@misc$tol, nmf_converged=model@misc$tol <= context$tolerance,
    nmf_relative_squared_error=model@misc$loss / context$input_energy,
    consensus_cophenetic=NA_real_, pac=NA_real_, median_resample_ari=NA_real_, q10_resample_ari=NA_real_,
    stability_status='Unavailable: one NMF fit per rank; no repeats or resampling')
  result <- list(model=model[factor_order], partition=partition, diagnostics=rank_row,
    fit_signature=context$fit_signature)
  temporary <- tempfile(pattern=paste0(basename(path), '.'), tmpdir=dirname(path))
  on.exit(unlink(temporary), add=TRUE)
  qs::qsave(result, temporary, preset='fast')
  stopifnot(file.rename(temporary, path))
  data.table::fwrite(rank_row, sub('\\.qs$', '.tsv', path), sep='\t')
  message(context$scope, '; saved NMF k=', rank, '; minimum genes=', min(sizes),
    '; silhouette=', signif(mean_silhouette, 8), '; checkpoint=', path)
  result$model <- NULL
  result$checkpoint_path <- path
  result
}

# Shared Ward.D2 cut diagnostics: the hierarchy is supplied once by the caller.
# Each cut is persisted before its completion is announced.
fit_ward_cuts <- function(distance, tree, grid, scope, orientation, nmf_input_signature, checkpoint_dir) {
  genes <- tree$labels
  results <- setNames(vector("list", length(grid)), as.character(grid))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  for (rank in grid) {
    path <- file.path(checkpoint_dir, sprintf("rank_%03d.qs", rank))
    fit_signature <- digest::digest(list(algorithm = "ward.D2", nmf_input_signature,
      genes, scope, orientation, rank, local_minimum_genes, R.version.string), algo = "sha256")
    if (file.exists(path)) {
      result <- qread(path)
      check(identical(result$fit_signature, fit_signature), "Saved Ward.D2 cut inputs changed")
    } else {
      partition <- ordered_partition(tree, rank)
      sizes <- tabulate(partition, nbins = rank)
      rank_row <- data.table(scope = scope, orientation = orientation, rank = rank,
        minimum_size = min(sizes), maximum_size = max(sizes),
        mean_silhouette = mean(cluster::silhouette(partition, distance)[, "sil_width"]),
        minimum_size_qc_pass = min(sizes) >= local_minimum_genes,
        nmf_fits = 0L, nmf_seed = NA_integer_, nmf_iterations = NA_integer_,
        nmf_final_tolerance = NA_real_, nmf_converged = NA, nmf_relative_squared_error = NA_real_,
        consensus_cophenetic = NA_real_, pac = NA_real_, median_resample_ari = NA_real_,
        q10_resample_ari = NA_real_, stability_status = "Unavailable: one Ward.D2 tree; no resampling")
      result <- list(model = NULL, partition = partition, diagnostics = rank_row, fit_signature = fit_signature)
      temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
      qsave(result, temporary, preset = "fast")
      check(file.rename(temporary, path), "Could not commit Ward.D2 cut")
      fwrite(rank_row, sub("\\.qs$", ".tsv", path), sep = "\t")
    }
    results[[as.character(rank)]] <- result
    message(scope, "; saved Ward.D2 cut k=", rank, "; silhouette=",
      signif(result$diagnostics$mean_silhouette, 8))
  }
  results
}

# Freeze programs after outcome-informed gene admission. Use one NMF fit per
# rank, or cuts of one Ward.D2 hierarchy, without repeats or subsampling.
discover_local <- function(x, null, metadata, genes, unit_column, scope, orientation, support = NULL,
  sensitivity_sign = FALSE, previous = NULL, checkpoint_dir = file.path(output_dir, scope, 'models'),
  rank_override = NULL, shared_structure = NULL) {
  x <- x[, genes, drop = FALSE]
  null <- null[, genes, drop = FALSE]
  units <- metadata[[unit_column]]
  weights <- group_weights(units) * metadata$condition_support
  outcomes <- sensitivity[match(metadata$observation_id, observation_id)]
  if (is.null(support)) support <- functional_support(x, metadata, outcomes, sensitivity_sign = sensitivity_sign)
  check(identical(rownames(support$loadings), genes), "Functional loading gene order differs from NMF genes")
  mixed_correlation <- identical(support$settings$clustering_representation,
    "mean_functional_displacement_correlation")
  ward_clustering <- identical(support$settings$clustering_algorithm, "ward.D2")
  check(is.null(support$settings$clustering_algorithm) || ward_clustering,
    "Unsupported functional clustering algorithm")
  channels <- support$channels
  nmf_input_signature <- digest::digest(channels, algo = "sha256")
  nmf_settings <- if (ward_clustering) list(fits_per_rank = 0L, resamples = 0L,
    algorithm = "ward.D2", input_transform = functional_input_transform) else
    list(fits_per_rank = 1L, resamples = 0L, threads_per_fit = 1L, workers = nmf_workers,
    tolerance = nmf_tolerance, max_iterations = nmf_max_iterations,
    package_version = as.character(packageVersion("RcppML")), input_transform = functional_input_transform)
  linkage_label <- if (ward_clustering) "hclust ward.D2 on mean functional-loading/displacement gene correlations" else
    "single-fit NMF of gene/condition functional-support profiles"
  grid <- local_rank_grid[local_rank_grid < min(ncol(x), nrow(x))]
  check(length(grid) > 0, paste("No feasible ranks in", scope))
  reused_grid <- integer()
  if (!is.null(shared_structure)) {
    check(is.null(previous) && ward_clustering &&
      all(genes %in% shared_structure$genes) && identical(shared_structure$tree$labels, shared_structure$genes) &&
      identical(shared_structure$tree$method, "ward.D2") &&
      all(grid %in% shared_structure$ranks$rank) &&
      !is.na(shared_structure$diagnostics$selected_rank), "Invalid fixed global Ward.D2 structure")
    correlation <- gene_correlation(x, weights)
    clustering_correlation <- shared_structure$clustering_correlation[genes, genes, drop = FALSE]
    nmf_input_signature <- shared_structure$nmf_input_signature
    support$channels <- NULL # The shared correlation input is stored once globally.
    reused_grid <- grid
  } else if (!is.null(previous)) {
    check(identical(previous$genes, genes) && identical(previous$scope, scope) &&
      identical(previous$orientation, orientation) &&
      identical(as.data.frame(previous$metadata), as.data.frame(metadata)) &&
      identical(previous$nmf_input_signature, nmf_input_signature) &&
      identical(previous$nmf_settings, nmf_settings), "Cached NMF inputs or fit settings changed")
    check(all(previous$ranks$rank %in% grid) &&
      all(previous$ranks$nmf_fits == if (ward_clustering) 0L else 1L),
      "Cached sweep is not a compatible extension")
    if (!ward_clustering) check(all(previous$ranks$nmf_seed ==
      seed + 100000L * match(metadata$depmap_id[1], cell_map$depmap_id) + previous$ranks$rank),
      "Cached NMF seeds changed")
    check(!is.null(previous$correlation) && !is.null(previous$clustering_correlation),
      "Restore cached correlation matrices before extending NMF")
    reused_grid <- previous$ranks$rank
    correlation <- previous$correlation
    clustering_correlation <- previous$clustering_correlation
  } else {
    correlation <- gene_correlation(x, weights)
    learning_input <- functional_nmf_input(support, correlation)
    clustering_correlation <- learning_input$correlation
    if (mixed_correlation) {
      channels <- learning_input$channels
      support$channels <- channels
      nmf_input_signature <- digest::digest(channels, algo = "sha256")
    }
    rm(learning_input)
  }
  # The same Pearson-distance silhouette now evaluates the learning input.
  # Retain raw molecular correlations above for the original comparison maps.
  distance <- if (!is.null(shared_structure$distance)) shared_structure$distance else
    as.dist(sqrt(2 * (1 - clustering_correlation)))
  tree <- if (!is.null(shared_structure)) shared_structure$tree else
    if (is.null(previous)) hclust(distance, method = "ward.D2") else previous$tree
  source_structure <- if (!is.null(shared_structure)) shared_structure else previous
  new_grid <- setdiff(grid, reused_grid)
  rank_results <- setNames(vector("list", length(grid)), as.character(grid))
  for (rank in reused_grid) {
    saved_rank <- rank
    saved_model <- if (!is.null(source_structure$nmf_model) && rank == source_structure$diagnostics$selected_rank)
      source_structure$nmf_model else NULL
    saved_diagnostics <- copy(source_structure$ranks[rank == saved_rank])
    set(saved_diagnostics, j = "scope", value = scope)
    set(saved_diagnostics, j = "orientation", value = orientation)
    rank_results[[as.character(rank)]] <- list(model = saved_model,
      partition = source_structure$all_partitions[[as.character(rank)]][genes], diagnostics = saved_diagnostics)
  }
  if (length(reused_grid)) message(scope, "; reusing ", length(reused_grid),
    if (ward_clustering) " saved Ward.D2 cuts; calculating " else " completed single-fit NMF ranks; fitting ",
    length(new_grid), " additional ranks")
  if (length(new_grid) && ward_clustering) {
    rank_results[as.character(new_grid)] <- fit_ward_cuts(distance, tree, new_grid,
      scope, orientation, nmf_input_signature, checkpoint_dir)
  }
  if (length(new_grid) && !ward_clustering) {
    dir.create(checkpoint_dir, recursive=TRUE, showWarnings=FALSE)
    seed_base <- seed + 100000L * match(metadata$depmap_id[1], cell_map$depmap_id)
    nmf_rank_context <- list(channels=channels, genes=genes, tree_order=tree$order, distance=distance,
      seed_base=seed_base, tolerance=nmf_tolerance, max_iterations=nmf_max_iterations,
      minimum_genes=local_minimum_genes, input_energy=sum(channels^2), scope=scope,
      orientation=orientation, checkpoint_dir=checkpoint_dir,
      fit_signature=digest::digest(list(nmf_input_signature, seed_base, nmf_tolerance,
        nmf_max_iterations, nmf_settings$package_version, genes, tree$order, scope, orientation,
        local_minimum_genes), algo='sha256'))
    # Fresh R workers avoid forking native numerical-library thread pools.
    nmf_cluster <- parallel::makePSOCKcluster(min(nmf_workers, length(new_grid)), outfile = "")
    on.exit(parallel::stopCluster(nmf_cluster), add = TRUE)
    invisible(parallel::clusterEvalQ(nmf_cluster, {
      suppressPackageStartupMessages({library(data.table); library(cluster); library(RcppML)})
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
      data.table::setDTthreads(1L)
    }))
    parallel::clusterExport(nmf_cluster, c('nmf_rank_context', 'fit_functional_nmf_rank'),
      envir = environment())
    new_results <- parallel::parLapplyLB(nmf_cluster, new_grid, 'fit_functional_nmf_rank', chunk.size=1L)
    rank_results[as.character(new_grid)] <- new_results
  }
  check(!any(vapply(rank_results, inherits, logical(1), "try-error")), "A single-rank NMF worker failed")
  models <- lapply(rank_results, `[[`, "model")
  partitions <- lapply(rank_results, `[[`, "partition")
  names(models) <- names(partitions) <- as.character(grid)
  rank_table <- rbindlist(lapply(rank_results, `[[`, "diagnostics"), use.names = TRUE, fill = TRUE)
  unit_ids <- sort(unique(units))
  bootstrap_seed <- NULL
  bootstrap_units <- list()
  agreement <- data.table(scope = character(), orientation = character(), rank = integer(), resample = integer(),
    adjusted_rand = numeric(), sampled_units = integer(), omitted_units = integer())
  selection_rule <- if (is.null(support$settings$rank_selection_rule)) "maximum_silhouette" else
    support$settings$rank_selection_rule
  rank_table <- select_silhouette(rank_table, rank_override, selection_rule)
  selected_rank <- if (any(rank_table$selected)) rank_table[selected == TRUE, rank] else NA_integer_
  if (!is.null(shared_structure)) check(identical(selected_rank, shared_structure$diagnostics$selected_rank),
    "Fixed global rank changed in a cell-line view")
  selection_status <- rank_table$selection_status[1]
  # An unresolved scientific criterion is a reported abstention, not an input
  # error or permission to select the largest rank. Other scopes still complete.
  if (is.na(selected_rank)) {
    null_correlation <- cor(null)
    diagnostics <- data.table(scope = scope, orientation = orientation, conditions = nrow(x),
      discovery_units = length(unit_ids), genes = length(genes), selected_rank = NA_integer_,
      selection_status = selection_status, linkage = linkage_label,
      raw_tree_cophenetic = cor(as.numeric(distance), as.numeric(cophenetic(tree))),
      correlation_with_control_structure = cor(correlation[upper.tri(correlation)],
        null_correlation[upper.tri(null_correlation)]), centering_adjusted_rand = NA_real_,
      mean_reconstruction_relative_error = NA_real_, null_profiles = nrow(null),
      experimental_replicate_agreement_available = FALSE)
    return(list(scope = scope, orientation = orientation, genes = genes, correlation = correlation,
      tree = tree, clustering_correlation = clustering_correlation, functional_support = support,
      nmf_settings = nmf_settings, nmf_input_signature = nmf_input_signature,
      metadata = metadata, ranks = rank_table, resampling = agreement,
      diagnostics = diagnostics, programs = data.table(), weights = NULL,
      all_partitions = partitions, bootstrap_seed = bootstrap_seed, bootstrap_units = bootstrap_units,
      functional_labels_used = TRUE, functional_labels_used_for_gene_selection = TRUE,
      rank_selection_uses_outcome_objective = TRUE))
  }
  if (!ward_clustering) {
    selected_result <- rank_results[[as.character(selected_rank)]]
    if (is.null(selected_result$checkpoint_path)) {
      selected_path <- file.path(checkpoint_dir, sprintf('rank_%03d.qs', selected_rank))
      if (file.exists(selected_path)) {
        cached <- qs::qread(selected_path)
        expected_signature <- digest::digest(list(nmf_input_signature,
          seed + 100000L * match(metadata$depmap_id[1], cell_map$depmap_id), nmf_tolerance,
          nmf_max_iterations, nmf_settings$package_version, genes, tree$order, scope, orientation,
          local_minimum_genes), algo='sha256')
        check(identical(cached$fit_signature, expected_signature), 'Selected rank checkpoint settings changed')
        models[[as.character(selected_rank)]] <- cached$model
        rm(cached)
      }
    }
    if (is.null(models[[as.character(selected_rank)]]) && !is.null(selected_result$checkpoint_path))
      models[[as.character(selected_rank)]] <- qs::qread(selected_result$checkpoint_path)$model
    check(!is.null(models[[as.character(selected_rank)]]) ||
      (!is.null(rank_override) && !is.null(previous$all_partitions[[as.character(selected_rank)]])),
      "Extended sweep selected an unavailable cached model; rank-selection settings must stay unchanged")
  }
  partition <- partitions[[as.character(selected_rank)]]
  active_programs <- if (is.null(shared_structure)) seq_len(selected_rank) else sort(unique(partition))
  program_ids <- sprintf("P%02d", active_programs)
  if (!is.null(shared_structure)) partition <- setNames(match(partition, active_programs), genes)
  program_count <- length(program_ids)
  gene_tree_order <- if (is.null(shared_structure)) tree$order else
    match(shared_structure$genes[tree$order][shared_structure$genes[tree$order] %in% genes], genes)
  # Seed modules define profiles; every admitted gene can have positive weight
  # in several programs. Positive correlation preserves co-displacement direction.
  seeds <- sapply(seq_len(program_count), function(j) rowMeans(x[, partition == j, drop = FALSE]))
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
  # A centering comparison would require a second NMF fit and is not performed.
  q <- qr(membership)
  check(q$rank == program_count, "Dependent local program weight columns")
  projected <- x %*% qr.Q(q)
  response_energy <- rowSums(x^2)
  reconstruction_error <- pmax(response_energy - rowSums(projected^2), 0) / response_energy
  # Control comparisons are descriptive; vehicle contrasts share reference plates.
  diagnostics <- data.table(scope = scope, orientation = orientation,
    conditions = nrow(x), discovery_units = length(unit_ids), genes = length(genes),
    selected_rank = selected_rank, selection_status = selection_status,
    linkage = linkage_label,
    raw_tree_cophenetic = cor(as.numeric(distance), as.numeric(cophenetic(tree))),
    correlation_with_control_structure = cor(correlation[upper.tri(correlation)],
      null_correlation[upper.tri(null_correlation)]),
    centering_adjusted_rand = NA_real_,
    mean_reconstruction_relative_error = mean(reconstruction_error),
    null_profiles = nrow(null), experimental_replicate_agreement_available = FALSE)
  program_table <- rbindlist(lapply(seq_len(program_count), function(j) {
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
    tree_order = match(seq_along(genes), gene_tree_order))
  list(scope = scope, orientation = orientation, genes = genes, correlation = correlation,
    rank_override = rank_override,
    global_structure = if (is.null(shared_structure)) NULL else
      list(path = shared_structure$path, signature = shared_structure$input_signature,
        selected_rank = selected_rank, active_programs = active_programs,
        program_colors = setNames(hcl.colors(selected_rank, "Dark 3"), sprintf("P%02d", seq_len(selected_rank)))),
    gene_tree_order = gene_tree_order,
    consensus = NULL, tree = tree, partition = partition, nmf_model = models[[as.character(selected_rank)]],
    clustering_correlation = clustering_correlation, functional_support = support,
    nmf_settings = nmf_settings, nmf_input_signature = nmf_input_signature,
    weights = membership, scores = scores, null_scores = null_scores, score_sd = score_sd,
    metadata = metadata, ranks = rank_table, resampling = agreement, diagnostics = diagnostics,
    programs = program_table, assignments = assignments, reconstruction_error = reconstruction_error,
    all_partitions = partitions, bootstrap_seed = bootstrap_seed, bootstrap_units = bootstrap_units,
    functional_labels_used = TRUE, functional_labels_used_for_gene_selection = TRUE,
    rank_selection_uses_outcome_objective = TRUE)
}

# Include completed diagnostics recovered from an interrupted sweep, even when
# its model checkpoint was lost. Precision and checkpoint availability are explicit.
merge_completed_rank_diagnostics <- function(bank) {
  path <- file.path(output_dir, bank$scope, 'data/completed_rank_diagnostics.tsv')
  ranks <- copy(bank$ranks)
  if (!file.exists(path)) return(ranks)
  completed <- fread(path)
  extra <- completed[!rank %in% ranks$rank]
  if (!'diagnostic_source' %in% names(ranks)) ranks[, diagnostic_source := 'saved model; full precision']
  ranks <- rbindlist(list(ranks, extra), use.names=TRUE, fill=TRUE)[order(rank)]
  selection_rule <- if (is.null(bank$functional_support$settings$rank_selection_rule)) "maximum_silhouette" else
    bank$functional_support$settings$rank_selection_rule
  ranks <- select_silhouette(ranks, bank$rank_override, selection_rule)
  check(identical(ranks[selected == TRUE, rank], bank$diagnostics$selected_rank),
    'Completed rank diagnostics select a different model; load that checkpoint before plotting')
  ranks
}

# Reused for every cell as soon as its bank has been saved or validated.
# Each scope owns a separate diagnostics PDF and PNG, with all cuts visible.
render_diagnostics <- function(bank, rank_plots_only = FALSE) {
  bank$ranks <- merge_completed_rank_diagnostics(bank)
  display_scope <- if (bank$orientation == "pooled") "Combined cell-line panel" else gsub("_", " ", bank$scope)
  if (!is.null(bank$display_name)) display_scope <- bank$display_name
  if (!is.null(bank$global_structure)) display_scope <- paste(display_scope, "| fixed global programs")
  dir.create(file.path(output_dir, bank$scope, "figures"), recursive = TRUE, showWarnings = FALSE)
  r <- bank$ranks
  ward_clustering <- identical(bank$functional_support$settings$clustering_algorithm, "ward.D2")
  last_blue_dot <- identical(bank$functional_support$settings$rank_selection_rule, "last_blue_dot")
  silhouette_elbow <- identical(bank$functional_support$settings$rank_selection_rule, "silhouette_elbow")
  p1 <- ggplot(r, aes(rank, mean_silhouette)) + geom_line() + geom_point(aes(color = minimum_size_qc_pass)) +
    scale_color_manual(values = c(`FALSE` = "#F8766D", `TRUE` = "#00BFC4")) +
    geom_vline(xintercept = local_minimum_rank - 0.5, lty = 3, color = "grey60") +
    geom_vline(data = r[selected == TRUE], aes(xintercept = rank), lty = 2, color = "#B24C35") +
    geom_point(data = r[selected == TRUE], shape = 1, size = 4, color = "#B24C35") +
    labs(x = "Candidate programs", y = "Mean gene silhouette",
      color = paste0("Size QC >=", local_minimum_genes),
      subtitle = if (is.null(bank$rank_override) && silhouette_elbow)
        sprintf("Select the steep-to-flatter silhouette elbow among size-qualified cuts %d--%d", min(r$rank), max(r$rank)) else
        if (is.null(bank$rank_override) && last_blue_dot)
        sprintf("Select largest positive-silhouette, size-qualified cut among %d--%d", min(r$rank), max(r$rank)) else
        if (is.null(bank$rank_override))
        sprintf("Select maximum positive mean silhouette among size-qualified ranks %d--%d", min(r$rank), max(r$rank)) else
        sprintf("User-selected rank %d; completed ranks %d--%d", bank$rank_override, min(r$rank), max(r$rank)))
  silhouette_plot <- p1 + labs(title = display_scope,
    caption = paste0("Red ring/dashed line: selected rank. Blue: every seed >=", local_minimum_genes,
      " genes.\n", if (is.null(bank$rank_override) && silhouette_elbow)
        "Select the bend of a continuous two-segment fit after the eligible peak; dashed grey: fitted decline and tail." else
        if (is.null(bank$rank_override) && last_blue_dot)
        "Select the largest rank with positive mean silhouette and minimum seed size (last blue dot)." else
        if (is.null(bank$rank_override))
        "Select the maximum positive mean silhouette; exact ties choose the smaller rank." else
        paste('Explicit user selection:', bank$rank_override, 'programs.'),
      if ('diagnostic_source' %in% names(r) && any(grepl('console', r$diagnostic_source)))
        '\nExtended completed ranks recovered from console scores (3 significant figures); unfinished ranks omitted.' else ''))
  if (silhouette_elbow) {
    fitted_curve <- r[is.finite(elbow_fitted_silhouette)]
    p1 <- p1 + geom_line(data = fitted_curve, aes(y = elbow_fitted_silhouette), lty = 2, color = "grey45")
    silhouette_plot <- silhouette_plot + geom_line(data = fitted_curve,
      aes(y = elbow_fitted_silhouette), lty = 2, color = "grey45")
  }
  ggsave(file.path(output_dir, bank$scope, "figures/silhouette.pdf"), silhouette_plot,
    width = 8, height = 4.5, device = cairo_pdf)
  ggsave(file.path(output_dir, bank$scope, "figures/silhouette.png"), silhouette_plot,
    width = 8, height = 4.5, dpi = 160, bg = "white")
  unavailable <- ggplot() + annotate("text", x = 0, y = 0, size = 4,
    label = if (ward_clustering) "Unavailable for one Ward.D2 tree\nNo resampling" else
      "Unavailable for one NMF fit per rank\nNo repeat fits or resampling") + theme_void()
  p2 <- unavailable + labs(title = "Consensus cophenetic")
  cophenetic_plot <- p2 + labs(title = display_scope,
    caption = if (is.null(bank$rank_override) && silhouette_elbow)
      "Rank selection uses the descending-to-flatter silhouette elbow; no stability resampling is performed." else
      if (is.null(bank$rank_override) && last_blue_dot)
      "Rank selection uses the last blue dot; no stability resampling is performed." else
      if (is.null(bank$rank_override))
      "Rank selection maximizes mean gene silhouette; no stability resampling is performed." else
      paste('User-selected rank', bank$rank_override, '; no stability resampling is performed.'))
  ggsave(file.path(output_dir, bank$scope, "figures/cophenetic.pdf"), cophenetic_plot,
    width = 8, height = 4.5, device = cairo_pdf)
  ggsave(file.path(output_dir, bank$scope, "figures/cophenetic.png"), cophenetic_plot,
    width = 8, height = 4.5, dpi = 160, bg = "white")
  p3 <- unavailable + labs(title = "Resampling adjusted Rand")
  p4 <- unavailable + labs(title = "Consensus ambiguity (PAC)")
  composite <- (p1 + p2) / (p3 + p4) + plot_annotation(title = display_scope,
    subtitle = sprintf("%d conditions; %d drug/recipient resampling units; raw-tree cophenetic %.3f; selected k=%d",
      bank$diagnostics$conditions, bank$diagnostics$discovery_units,
      bank$diagnostics$raw_tree_cophenetic, bank$diagnostics$selected_rank),
    caption = if (ward_clustering)
      "Gene admission uses retrospective outcomes; silhouette evaluates the mean correlation distance.\nOne Ward.D2 tree; no repeats, subsampling or resampling stability diagnostics." else
      paste("Gene admission uses retrospective outcomes; silhouette evaluates functional-loading profiles.",
        "\nOne NMF fit per rank; no repeated fits, subsampling or resampling stability diagnostics."))
  ggsave(file.path(output_dir, bank$scope, "figures/diagnostics.pdf"), composite,
    width = 9.8, height = 7, device = cairo_pdf)
  ggsave(file.path(output_dir, bank$scope, "figures/diagnostics.png"), composite,
    width = 9.8, height = 7, dpi = 160, bg = "white")
  scope_data_dir <- file.path(output_dir, bank$scope, "data")
  dir.create(scope_data_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(bank$ranks, file.path(scope_data_dir, "rank_diagnostics.tsv"), sep = "\t")
  fwrite(bank$resampling, file.path(scope_data_dir, "resampling_diagnostics.tsv"), sep = "\t")
  fwrite(data.table(resample = integer(), draw = integer(), drug = character(), bootstrap_seed = integer()),
    file.path(scope_data_dir, "resampling_units.tsv.gz"), sep = "\t")
  fwrite(bank$diagnostics, file.path(scope_data_dir, "discovery_diagnostics.tsv"), sep = "\t")
  if (rank_plots_only) return(invisible(NULL))
  cell_id <- if (bank$orientation == "pooled") "pooled" else
    if (bank$orientation == "drug") unique(bank$metadata$drug) else unique(bank$metadata$depmap_id)
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
    sprintf("Genes: two-sided Pearson with PRISM or ALMANAC; p <=%g AND |r| >=%g in the same assay; both signs; no gene cap.",
      outcome_p_threshold, outcome_minimum_abs_correlation),
    "Assays tested separately; matched doses pooled with equal condition weight; BH values diagnostic only.",
    paste("Rank selection:", r$selection_rule[1L]),
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

# Render each cell-line bank. Every heatmap uses the shared
# ComplexHeatmap wrapper. Off-block values and signed original effects remain.
# Draw full-compendium overviews and readable pages from one frozen dictionary.
# Reused annotation/heatmap construction stays here beside its scientific owner.
render_local <- function(bank, detail_pages = TRUE, cluster_within_program = TRUE) {
  bank$ranks <- merge_completed_rank_diagnostics(bank)
  scope <- bank$scope
  display_scope <- if (bank$orientation == "pooled") "Combined cell-line panel" else gsub("_", " ", scope)
  scope_data_dir <- file.path(output_dir, scope, "data")
  scope_figure_dir <- file.path(output_dir, scope, "figures")
  for (path in c(scope_data_dir, scope_figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  program_ids <- colnames(bank$weights)
  k <- length(program_ids)
  many_programs <- k > 24L
  program_legend <- if (many_programs) list(annotation_legend_param = list(
    Program = list(ncol = ceiling(k / 20), labels_gp = grid::gpar(fontsize = 7)))) else list()
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
  dominant <- max.col(abs(standardized), ties.method = "first")
  gene_tree_order <- if (is.null(bank$gene_tree_order)) bank$tree$order else bank$gene_tree_order
  row_order <- order(bank$partition, match(seq_along(bank$genes), gene_tree_order))
  column_order <- unlist(lapply(seq_len(k), function(j) {
    rows <- which(dominant == j)
    if (cluster_within_program && length(rows) > 1)
      rows <- rows[hclust(dist(standardized[rows, , drop = FALSE]), method = "average")$order]
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
  palette <- list(Program = if (is.null(bank$global_structure)) setNames(hcl.colors(k, "Dark 3"), program_ids) else
    bank$global_structure$program_colors[program_ids])
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
  if (bank$orientation %in% c("pooled", "drug")) labels <- paste(md$recipient, labels, sep = " | ")

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
      row_anno_cols = "Program", ann_row_args = program_legend,
      ann_top_args = top_args(indices), palette_overrides = palette,
      hm_name = "log2 delta", scale_rows = FALSE, value_range = c(-limit, 0, limit),
      value_colors = c("#2166AC", "white", "#B2182B"),
      heatmap_legend_param = list(at = c(-limit, 0, limit), labels = sprintf("%.2f", c(-limit, 0, limit))),
      show_row_names = show_labels, show_column_names = show_labels, column_labels = labels[indices],
      row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 5),
      row_split = rm$Program, column_split = cm$Program,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      raster_resize_mat = bank$orientation == 'pooled' && !identical(bank$raster_resize_matrix, FALSE) && !show_labels,
      row_title_gp = grid::gpar(fontsize = if (many_programs) 5 else 8),
      row_title_rot = if (many_programs) 0 else 90,
      row_gap = grid::unit(if (many_programs) 2 else 1, 'mm'),
      column_title_gp = grid::gpar(fontsize = 10),
      title = paste0(display_scope, " | outcome-selected genes | ", caption, "\nSigned original displacement; ",
        "color saturates at +/-", sprintf("%.2f", limit), "; no row scaling",
        if (bank$orientation == 'pooled' && !show_labels && !identical(bank$raster_resize_matrix, FALSE)) '; overview pixels average adjacent values' else ''))
  }
  all_columns <- seq_len(nrow(md))
  overview <- displacement_heatmap(all_columns, row_order, FALSE,
    paste(ncol(full_matrix), "conditions;", nrow(full_matrix), "genes; selected k=", k,
      "\nWithin program: PRISM sensitivity, ALMANAC-only sensitivity, then no data; lower % vehicle first"))
  save_heatmap(overview, "displacement_overview", if (nrow(md) > 60) 16 else 11.5,
    if (many_programs) max(13, 4 + .16 * k) else 9,
    scope_figure_dir, list(sensitivity_legend))
  # Detail pages show all conditions in the identical overview order, at most
  # 48 columns/page. Full unabridged drug names are in condition_order.tsv.
  page_height <- max(8.5, 5 + length(selected_rows) * 0.06)
  # Return the same detail renderer when pooling so summary/loading plots can
  # be saved before the hundreds of detail pages. Cell defaults are unchanged.
  render_details <- function() {
    if (bank$orientation == 'pooled')
      fwrite(cbind(data.table(gene = rownames(full_matrix)), as.data.table(full_matrix)),
        file.path(scope_data_dir, 'signed_displacements.tsv.gz'), sep='\t')
    pages <- split(all_columns, ceiling(all_columns / 48))
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
    if (bank$orientation %in% c("cell", "drug")) {
      focal_columns <- if (bank$orientation == "drug") head(all_columns, 48L) else which(md$focal)
      hm <- displacement_heatmap(focal_columns, selected_rows, TRUE,
        if (bank$orientation == "drug") "cell-line focus; first 48 conditions; top 12 seed genes/program" else
          "seven-drug focus; top 12 seed genes/program")
      save_heatmap(hm, "focal_panel", 11.5, page_height, scope_figure_dir, list(sensitivity_legend))
    }
    invisible(NULL)
  }
  if (detail_pages) render_details()
  # Full gene correlation and resampling consensus share the gene block order.
  mark <- bank$genes[unlist(lapply(seq_len(k), function(j) {
    idx <- which(bank$partition == j)
    head(idx[order(bank$weights[idx, j], decreasing = TRUE)], 3L)
  }))]
  mark_table <- data.frame(gene = mark,
    Program = program_ids[bank$partition[match(mark, bank$genes)]])
  for (kind in "gene_correlation") {
    matrix <- if (kind == "gene_correlation") bank$correlation else bank$consensus
    matrix <- matrix[row_order, row_order, drop = FALSE]
    stopifnot(identical(rownames(matrix), rownames(row_meta)) || is.null(rownames(matrix)))
    dimnames(matrix) <- list(bank$genes[row_order], bank$genes[row_order])
    hm <- draw_heatmap(matrix, col_anno_df = row_meta, col_anno_cols = "Program",
      row_anno_df = row_meta, row_anno_cols = "Program", ann_row_args = program_legend, palette_overrides = palette,
      ann_top_args = list(show_legend = FALSE), top_labels = mark_table,
      hm_name = if (kind == "gene_correlation") "Pearson r" else "co-clustering",
      scale_rows = FALSE, value_range = if (kind == "gene_correlation") c(-1, 0, 1) else c(0, 0.5, 1),
      value_colors = if (kind == "gene_correlation") c("#2166AC", "white", "#B2182B") else
        c("white", "#9ECAE1", "#08519C"),
      row_split = row_meta$Program, column_split = row_meta$Program,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      raster_resize_mat = bank$orientation == 'pooled' && !identical(bank$raster_resize_matrix, FALSE),
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 6),
      show_row_names = FALSE,
      title = paste(display_scope, "|", gsub("_", " ", kind), "| all", length(bank$genes), "genes"))
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
    title = paste(display_scope, "| original signed program scores |", nrow(md), "conditions"))
  save_heatmap(hm, "program_scores", if (nrow(md) > 60) 16 else 11.5, 6, scope_figure_dir, list(sensitivity_legend))
  # Export complete sources, rank partitions, and identical plotted identifiers.
  fwrite(cbind(data.table(gene = bank$genes), as.data.table(bank$weights)),
    file.path(scope_data_dir, "gene_weights.tsv.gz"), sep = "\t")
  # The exact correlation source is already the exported signed displacement
  # matrix plus these weights. Avoid a second quadratic-size copy per cell.
  fwrite(data.table(observation_id = bank$metadata$observation_id,
    weight = group_weights(bank$metadata[[if (!is.null(bank$unit_column)) bank$unit_column else
      if (bank$orientation == "pooled") "cell_drug" else
      if (bank$orientation == "drug") "depmap_id" else "drug"]]) *
      bank$metadata$condition_support), file.path(scope_data_dir, "gene_correlation_weights.tsv"), sep = "\t")
  if (bank$orientation != 'pooled')
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
  invisible(if (detail_pages) NULL else render_details)
}

# If no allowed blue-point silhouette cut is supported, retain raw molecular inspection
# with sensitivity annotations. Dendrogram order is not a selected program cut.
render_unresolved <- function(bank, column_order = NULL) {
  md <- bank$metadata
  annotation <- sensitivity[match(md$observation_id, observation_id)]
  check(identical(annotation$observation_id, md$observation_id), "Unresolved-map sensitivity mismatch")
  genes <- bank$genes[bank$tree$order]
  response <- prepared$response[match(md$observation_id, rownames(prepared$response)), genes, drop = FALSE]
  condition_order <- if (is.null(column_order)) hclust(dist(response), method = "average")$order else column_order
  check(setequal(condition_order, seq_len(nrow(response))), "Incomplete unresolved column order")
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
      column_labels = paste(if (bank$orientation == "drug") md$recipient[indices] else substr(display_drugs, 1, 38),
        sprintf("%g uM", md$dose[indices])),
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
      col_anno_df = column_metadata, col_anno_cols = names(column_metadata), ann_top_args = annotation_args,
      title = paste0(gsub("_", " ", bank$scope), " | No selected program rank\n", caption,
        "\nColor saturates at +/-", sprintf("%.2f", limit), " in original signed displacement units"))
  }
  heatmap <- unresolved_map(seq_len(ncol(x)), seq_len(nrow(x)), FALSE, "raw clustered overview; no row scaling")
  destination <- file.path(output_dir, bank$scope, "figures")
  save_heatmap(heatmap, "displacement_unresolved", 12.8, 8, destination, list(sensitivity_legend))
  # A readable 21-condition supplement uses 48 high-variance admitted genes,
  # retaining dendrogram order. This is a display subset, not a selected module.
  focal_columns <- if (bank$orientation == "drug") head(seq_len(nrow(md)), 48L) else which(md$focal)
  focal_columns <- focal_columns[order(md$drug[focal_columns], md$dose[focal_columns])]
  focal_genes <- head(genes[order(matrixStats::rowVars(x), decreasing = TRUE)], 48L)
  focal_gene_rows <- which(genes %chin% focal_genes)
  heatmap <- unresolved_map(focal_columns, focal_gene_rows, TRUE,
    if (bank$orientation == "drug") "cell-line focus; 48 high-variance admitted genes; no row scaling" else
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

# Reused for each completed cell. Read the existing displacement source order,
# never recluster either axis. Missing functional observations stay grey.
render_functional_loadings <- function(bank) {
  display_scope <- if (bank$orientation == "pooled") "Combined cell-line panel" else gsub("_", " ", bank$scope)
  source_dir <- file.path(output_dir, bank$scope, "data")
  destination <- file.path(output_dir, bank$scope, "figures")
  gene_order <- fread(file.path(source_dir, "gene_order.tsv"))
  condition_order <- fread(file.path(source_dir, "condition_order.tsv"))
  check(setequal(gene_order$gene, bank$genes) && !anyDuplicated(gene_order$gene) &&
    setequal(condition_order$observation_id, bank$metadata$observation_id), "Incomplete loading/displacement order")
  program_ids <- if (is.null(bank$weights)) character() else colnames(bank$weights)
  program_legend <- if (length(program_ids) > 24L) list(annotation_legend_param = list(
    Program = list(ncol = ceiling(length(program_ids) / 20), labels_gp = grid::gpar(fontsize = 7)))) else list()
  row_annotation <- if (is.null(bank$weights)) NULL else data.frame(
    Program = factor(program_ids[bank$partition[match(gene_order$gene, bank$genes)]], levels = program_ids),
    row.names = gene_order$gene)
  column_annotation <- if (is.null(bank$weights)) NULL else data.frame(
    Program = factor(condition_order$dominant_program, levels = program_ids),
    row.names = condition_order$observation_id)
  palette <- if (is.null(bank$weights)) NULL else list(Program = if (is.null(bank$global_structure))
    setNames(hcl.colors(length(program_ids), "Dark 3"), program_ids) else bank$global_structure$program_colors[program_ids])
  matrices <- list()
  for (assay_name in c("PRISM", "ALMANAC")) {
    loading <- matrix(NA_real_, nrow = nrow(gene_order), ncol = nrow(condition_order),
      dimnames = list(gene_order$gene, condition_order$observation_id))
    observed <- bank$functional_support$observations[assay == assay_name, observation_id]
    if (length(observed)) {
      loading[, observed] <- bank$functional_support$loadings[gene_order$gene, paste0(assay_name, ":", observed), drop = FALSE]
    }
    matrices[[assay_name]] <- loading
  }
  # The measured loading values are exactly the finite entries in the maps.
  # Avoid materializing another full copy of the NA-padded pooled panels.
  limit <- max(0.1, as.numeric(quantile(abs(bank$functional_support$loadings), .98)))
  for (focal_only in c(FALSE, TRUE)) {
    rows <- if (!focal_only) seq_len(nrow(gene_order)) else
      if (is.null(bank$weights)) which(gene_order$focal_display) else which(gene_order$detail_label)
    columns <- if (focal_only && bank$orientation == "drug") head(seq_len(nrow(condition_order)), 48L) else
      if (focal_only) which(condition_order$focal) else seq_len(nrow(condition_order))
    # A pooled focus contains 39 cells: preserve its order and page labels
    # rather than squeezing hundreds of named columns into one small panel.
    column_pages <- if (focal_only && bank$orientation == "pooled" && length(columns) > 48L)
      split(columns, ceiling(seq_along(columns) / 48L)) else list(columns)
    make_maps <- function(columns) lapply(c("PRISM", "ALMANAC"), function(assay_name) {
      matrix <- matrices[[assay_name]][rows, columns, drop = FALSE]
      row_meta <- if (is.null(row_annotation)) NULL else row_annotation[rows, , drop = FALSE]
      col_meta <- if (is.null(column_annotation)) NULL else column_annotation[columns, , drop = FALSE]
      draw_heatmap(matrix, col_anno_df = col_meta, row_anno_df = row_meta,
        col_anno_cols = if (is.null(col_meta)) NULL else "Program",
        row_anno_cols = if (is.null(row_meta)) NULL else "Program", ann_row_args = program_legend,
        palette_overrides = palette,
        ann_top_args = list(show_legend = FALSE),
        row_split = if (is.null(row_meta)) NULL else row_meta$Program,
        column_split = if (is.null(col_meta)) NULL else col_meta$Program,
        cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        raster_resize_mat = bank$orientation == 'pooled' && !identical(bank$raster_resize_matrix, FALSE) && !focal_only,
        row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 6),
        scale_rows = FALSE, show_row_names = focal_only, show_column_names = focal_only,
        row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
        column_labels = if (bank$orientation %in% c("pooled", "drug"))
          paste(condition_order$recipient[columns], condition_order$drug[columns],
            sprintf("%g uM", condition_order$dose[columns]), sep = " | ") else
          paste(condition_order$drug[columns], sprintf("%g uM", condition_order$dose[columns])),
        hm_name = paste(assay_name, "gain"), na_col = "#D9D9D9", value_range = c(-limit, 0, limit),
        value_colors = c("#2166AC", "white", "#B2182B"),
        title = paste0(display_scope, " | ", assay_name,
          " functional loadings | same displacement order\n",
          if (is.null(bank$functional_support$raw_prediction_gain))
            "Positive: improves prediction; negative: worsens; grey: no outcome" else
            "Original gain x sign(50 - measured abundance); no clipping; grey: no outcome",
          if (bank$orientation == 'pooled' && !focal_only && !identical(bank$raster_resize_matrix, FALSE)) '\nOverview pixels average adjacent values; full values in source tables' else ''))
    })
    width <- if (focal_only) 11.5 else 16
    height <- if (focal_only) 2 * max(8.5, 5 + length(rows) * 0.06) else 15
    name <- if (focal_only) "functional_loadings_focal" else "functional_loadings"
    if (length(column_pages) == 1L) {
      maps <- make_maps(column_pages[[1L]])
      combined <- ComplexHeatmap::`%v%`(maps[[1L]], maps[[2L]])
      save_heatmap(combined, name, width, height, destination)
    } else {
      cairo_pdf(file.path(destination, paste0(name, ".pdf")), width = width, height = height)
      for (page in column_pages) {
        maps <- make_maps(page)
        ComplexHeatmap::draw(ComplexHeatmap::`%v%`(maps[[1L]], maps[[2L]]), merge_legends = TRUE)
      }
      dev.off()
      maps <- make_maps(column_pages[[1L]])
      grDevices::png(file.path(destination, paste0(name, ".png")), width = width, height = height,
        units = "in", res = 160)
      ComplexHeatmap::draw(ComplexHeatmap::`%v%`(maps[[1L]], maps[[2L]]), merge_legends = TRUE)
      dev.off()
    }
  }
  for (assay_name in names(matrices)) {
    loading <- matrices[[assay_name]]
    fwrite(cbind(data.table(gene=rownames(loading)), as.data.table(loading)),
      file.path(source_dir, paste0(tolower(assay_name), '_functional_loadings.tsv.gz')), sep='\t')
  }
  fwrite(bank$functional_support$observations, file.path(source_dir, "functional_prediction_context.tsv"), sep = "\t")
  fwrite(bank$functional_support$coefficients, file.path(source_dir, "functional_ridge_coefficients.tsv.gz"), sep = "\t")
  fwrite(bank$functional_support$folds, file.path(source_dir, "functional_drug_folds.tsv"), sep = "\t")
  if (!is.null(bank$functional_support$raw_prediction_gain)) {
    raw_gain <- bank$functional_support$raw_prediction_gain
    fwrite(cbind(data.table(gene = rownames(raw_gain)), as.data.table(raw_gain)),
      file.path(source_dir, "original_prediction_gain.tsv.gz"), sep = "\t")
  }
  if (is.null(bank$global_structure)) {
    fwrite(cbind(data.table(gene = bank$genes), as.data.table(bank$functional_support$channels)),
      file.path(source_dir, "nmf_functional_input.tsv.gz"), sep = "\t")
  } else writeLines(bank$global_structure$path, file.path(source_dir, "global_clustering_input.txt"))
}


# Same outcome matching used for cohort coverage and per-cell modeling.
align_panel_outcomes <- function(metadata, cell_map, data_dir) {
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
  sensitivity
}
# Aggregate counts are normalized with the original paired trimmed-median factor.
composition_factor <- function(treated, reference) {
  check(all(is.finite(treated)) && all(is.finite(reference)) &&
    all(treated >= 0) && all(reference >= 0), "Invalid count profile for normalization")
  tp <- treated / sum(treated)
  rp <- reference / sum(reference)
  valid <- tp > 0 & rp > 0 & tp * 10000 >= 0.1 & rp * 10000 >= 0.1
  check(sum(valid) >= 100, "Insufficient expressed genes for paired normalization")
  m <- log2(tp[valid]) - log2(rp[valid])
  a <- (log2(tp[valid]) + log2(rp[valid])) / 2
  mr <- quantile(m, c(0.3, 0.7))
  ar <- quantile(a, c(0.05, 0.95))
  keep <- m >= mr[1] & m <= mr[2] & a >= ar[1] & a <= ar[2]
  check(sum(keep) >= 100, "Insufficient genes after paired normalization trimming")
  2^median(m[keep])
}

# Derived correlations and signed channels have exact reconstruction sources.
# Save the fit, raw loadings, and hashes once rather than duplicating quadratic
# matrices across every cell. Rehydration verifies all reconstructed values.
save_functional_bank <- function(bank, path) {
  bank$derived_matrix_signatures <- list(correlation=digest::digest(bank$correlation,algo='sha256'),
    clustering_correlation=digest::digest(bank$clustering_correlation,algo='sha256'))
  bank$correlation <- NULL
  bank$clustering_correlation <- NULL
  bank$functional_support$channels <- NULL
  temporary <- tempfile(pattern=paste0(basename(path), '.'), tmpdir=dirname(path))
  on.exit(unlink(temporary), add=TRUE)
  qsave(bank,temporary,preset='high')
  check(file.rename(temporary, path), 'Could not commit saved program bank')
}
restore_functional_bank <- function(bank, response, shared_structure = NULL) {
  if(is.null(bank$derived_matrix_signatures)) return(bank)
  md <- bank$metadata
  x <- response[match(md$observation_id,rownames(response)),bank$genes,drop=FALSE]
  unit <- if (!is.null(bank$unit_column)) bank$unit_column else
    if(bank$orientation=='pooled') 'cell_drug' else if(bank$orientation=='drug') 'depmap_id' else 'drug'
  bank$correlation <- gene_correlation(x,group_weights(md[[unit]])*md$condition_support)
  if (is.null(bank$global_structure)) {
    learning_input <- functional_nmf_input(bank$functional_support, bank$correlation)
  } else {
    if (is.null(shared_structure)) shared_structure <- qread(bank$global_structure$path)
    check(identical(shared_structure$input_signature, bank$global_structure$signature) &&
      all(bank$genes %in% shared_structure$genes) &&
      identical(shared_structure$diagnostics$selected_rank, bank$global_structure$selected_rank),
      "Fixed global clustering reference changed")
    learning_input <- list(correlation = shared_structure$clustering_correlation[bank$genes, bank$genes, drop = FALSE], channels = NULL)
  }
  bank$clustering_correlation <- learning_input$correlation
  check(identical(digest::digest(bank$correlation,algo='sha256'),bank$derived_matrix_signatures$correlation) &&
    identical(digest::digest(bank$clustering_correlation,algo='sha256'),bank$derived_matrix_signatures$clustering_correlation),
    'Derived correlation reconstruction differs from saved fit')
  bank$functional_support$channels <- learning_input$channels
  bank
}

# Additional comparison maps: retain the learned gene programs and raw values,
# but cluster complete signed program-score vectors to order drug/dose columns.
# SD scaling gives each program equal units; signs and relative amplitudes remain.
render_program_score_clustered_maps <- function(bank, response) {
  check(!is.null(bank$weights), "Program-score ordering requires selected programs")
  md <- copy(bank$metadata)
  x <- response[match(md$observation_id, rownames(response)), bank$genes, drop = FALSE]
  check(identical(rownames(x), md$observation_id) && all(is.finite(x)),
    "Program-score map displacement rows do not align")
  scores <- x %*% bank$weights
  check(max(abs(scores - bank$scores)) < 1e-10 && all(bank$score_sd > 0),
    "Program-score map differs from saved program scores")
  standardized <- sweep(scores, 2L, bank$score_sd, "/")
  tree <- hclust(dist(standardized), method = "ward.D2")
  gene_tree_order <- if (is.null(bank$gene_tree_order)) bank$tree$order else bank$gene_tree_order
  rows <- order(bank$partition, match(seq_along(bank$genes), gene_tree_order))
  programs <- colnames(bank$weights)
  many_programs <- length(programs) > 24L
  map_height <- if (many_programs) max(13, 4 + .16 * length(programs)) else 9
  program_legend <- if (many_programs) list(annotation_legend_param = list(
    Program = list(ncol = ceiling(length(programs) / 20), labels_gp = grid::gpar(fontsize = 7)))) else list()
  row_annotation <- data.frame(Program = factor(programs[bank$partition[rows]], levels = programs),
    row.names = bank$genes[rows])
  column_annotation <- data.frame(Dose = sprintf("%g", md$dose), row.names = md$observation_id)
  palette <- list(Program = if (is.null(bank$global_structure)) setNames(hcl.colors(length(programs), "Dark 3"), programs) else
    bank$global_structure$program_colors[programs])
  annotation <- sensitivity[match(md$observation_id, observation_id)]
  check(identical(annotation$observation_id, md$observation_id), "Program-score map assay rows do not align")
  source_dir <- file.path(output_dir, bank$scope, "data")
  destination <- file.path(output_dir, bank$scope, "figures")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  order_table <- cbind(data.table(column_position = seq_len(nrow(md))), md[tree$order],
    annotation[tree$order, .(prism_abundance_pct, almanac_abundance_pct)])
  fwrite(order_table, file.path(source_dir, "condition_order_program_score_clustered.tsv"), sep = "\t")
  qsave(list(tree = tree, score_scaling = "divide each signed score by its across-condition SD",
    standardized_scores_sha256 = digest::digest(standardized, algo = "sha256"),
    gene_order = bank$genes[rows], column_order = md$observation_id[tree$order]),
    file.path(output_dir, bank$scope, "models/program_score_column_tree.qs"), preset = "high")
  display_scope <- gsub("_", " ", bank$scope)
  # Matrix columns and every annotation start in saved bank order. ComplexHeatmap
  # applies this exact tree without dendrogram reordering or program-based splits.
  make_map <- function(matrix, name, limit, caption, assay_bars = FALSE, dendrogram = TRUE) {
    check(identical(rownames(matrix), rownames(row_annotation)) &&
      identical(colnames(matrix), rownames(column_annotation)), "Program-score heatmap annotation mismatch")
    top <- list(show_legend = FALSE)
    if (assay_bars) {
      av <- annotation$almanac_abundance_pct
      pv <- annotation$prism_abundance_pct
      av[is.na(av)] <- 0
      pv[is.na(pv)] <- 0
      top <- c(top, list(`ALMANAC % vehicle` = ComplexHeatmap::anno_barplot(av,
        gp = grid::gpar(fill = "#B24C35", col = NA), ylim = c(0, almanac_bar_max), height = grid::unit(11, "mm")),
        `PRISM % vehicle` = ComplexHeatmap::anno_barplot(pv,
          gp = grid::gpar(fill = "#27818A", col = NA), ylim = c(0, prism_bar_max), height = grid::unit(11, "mm")),
        `ALMANAC value` = sensitivity_track(annotation$almanac_abundance_pct),
        `PRISM value` = sensitivity_track(annotation$prism_abundance_pct)))
    }
    draw_heatmap(matrix, row_anno_df = row_annotation, row_anno_cols = "Program", ann_row_args = program_legend,
      col_anno_df = column_annotation, col_anno_cols = "Dose", ann_top_args = top,
      palette_overrides = palette, hm_name = name, scale_rows = FALSE,
      value_range = c(-limit, 0, limit), value_colors = c("#2166AC", "white", "#B2182B"),
      na_col = "#D9D9D9", row_split = row_annotation$Program,
      cluster_rows = FALSE, cluster_row_slices = FALSE, cluster_columns = tree,
      column_dend_reorder = FALSE, show_column_dend = dendrogram,
      show_row_names = FALSE, show_column_names = FALSE,
      row_title_gp = grid::gpar(fontsize = if (many_programs) 5 else 8),
      row_title_rot = if (many_programs) 0 else 90, row_gap = grid::unit(if (many_programs) 2 else 1, 'mm'),
      title = paste0(display_scope, " | ", caption,
        "\nColumns: Ward.D2 on complete signed program-score vectors / program SD"))
  }
  displacement <- t(x[, rows, drop = FALSE])
  displacement_limit <- max(0.25, as.numeric(quantile(abs(displacement), .98)))
  map <- make_map(displacement, "log2 delta", displacement_limit,
    "signed original displacements; gene programs unchanged", assay_bars = TRUE)
  save_heatmap(map, "displacement_program_score_clustered", 16, map_height, destination, list(sensitivity_legend))
  # A fixed global gene union can include genes never fitted in this drug.
  # Such display entries remain NA; derive the color scale from fitted entries.
  fitted_loadings <- bank$functional_support$loadings[is.finite(bank$functional_support$loadings)]
  loading_limit <- if (length(fitted_loadings))
    max(0.1, as.numeric(quantile(abs(fitted_loadings), .98))) else 0.1
  loading_maps <- lapply(c("PRISM", "ALMANAC"), function(assay_name) {
    loading <- matrix(NA_real_, nrow = length(rows), ncol = nrow(md),
      dimnames = list(bank$genes[rows], md$observation_id))
    observed <- bank$functional_support$observations[assay == assay_name, observation_id]
    if (length(observed)) loading[, observed] <- bank$functional_support$loadings[
      bank$genes[rows], paste0(assay_name, ":", observed), drop = FALSE]
    make_map(loading, paste(assay_name, "gain"), loading_limit,
      paste(assay_name, "signed functional loadings; same displacement order; grey: no outcome"),
      dendrogram = assay_name == "PRISM")
  })
  combined <- ComplexHeatmap::`%v%`(loading_maps[[1L]], loading_maps[[2L]])
  save_heatmap(combined, "functional_loadings_program_score_clustered", 16,
    if (many_programs) 2 * map_height else 15, destination)
  message("Program-score clustered maps ready: ", bank$scope)
  invisible(list(gene_order = bank$genes[rows], column_order = md$observation_id[tree$order]))
}

# Exact script-15 admission test, reused for excluded-family reconstruction.
regression_admission <- function(support, id, recipient, p_threshold, minimum_loading_iqr,
  unit_column = "drug") {
  candidates <- rownames(support$loadings)
  regression_tests <- rbindlist(lapply(unique(support$observations$assay), function(assay_name) {
    columns <- which(support$observations$assay == assay_name)
    gain <- support$raw_prediction_gain[, columns, drop = FALSE]
    loading <- support$loadings[, columns, drop = FALSE]
    drug <- support$observations[[if (unit_column == "drug") "drug" else "regression_unit"]][columns]
    drug_sums <- rowsum(t(gain), drug, reorder = FALSE)
    drug_counts <- as.numeric(table(factor(drug, levels = rownames(drug_sums))))
    drug_gain <- t(sweep(drug_sums, 1L, drug_counts, "/"))
    mean_gain <- rowMeans(drug_gain)
    gain_se <- matrixStats::rowSds(drug_gain) / sqrt(ncol(drug_gain))
    estimable <- is.finite(gain_se) & gain_se > 0
    statistic <- p_value <- rep(NA_real_, length(candidates))
    statistic[estimable] <- mean_gain[estimable] / gain_se[estimable]
    p_value[estimable] <- stats::pt(statistic[estimable], df = ncol(drug_gain) - 1L, lower.tail = FALSE)
    loading_iqr <- matrixStats::rowIQRs(loading)
    result <- data.table(depmap_id = id, recipient = recipient,
      gene = candidates, assay = assay_name, matched_conditions = ncol(gain), matched_drugs = ncol(drug_gain),
      mean_gain_equal_condition = rowMeans(gain), mean_gain_equal_drug = mean_gain,
      drug_gain_standard_error = gain_se, statistic = statistic, p_value = p_value,
      test_status = fifelse(estimable, "tested", "constant drug-level gain: not estimable"),
      loading_iqr = loading_iqr, minimum_loading_iqr = minimum_loading_iqr)
    result[, q_value := p.adjust(p_value, method = "BH")]
    result[, `:=`(regression_significant = !is.na(p_value) & p_value <= p_threshold,
      loading_variation_eligible = loading_iqr >= minimum_loading_iqr)]
    result[, selected_in_assay := regression_significant & loading_variation_eligible]
    result
  }))
  if (unit_column != "drug") {
    setnames(regression_tests, c("matched_drugs", "mean_gain_equal_drug", "drug_gain_standard_error"),
      c("matched_cell_lines", "mean_gain_equal_cell", "cell_gain_standard_error"))
    regression_tests[test_status == "constant drug-level gain: not estimable",
      test_status := "constant cell-level gain: not estimable"]
  }
  regression_tests
} 
