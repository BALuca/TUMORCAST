# Nonlinear/calibration comparison on the existing frozen family-transfer bank.
# Invoked through canonical script 20 with --family-transfer=MEK|ERBB --nonlinear.
# Only decoder fitting changes; each completed fold and final model is committed.
fit_program_svr <- function(x, y, kernel, cost, gamma) {
  mu <- colMeans(x)
  sigma <- apply(x, 2L, sd)
  active <- which(sigma > 1e-12)
  stopifnot(length(active) >= 2L, all(is.finite(x)), all(is.finite(y)), sd(y) > 0)
  z <- sweep(sweep(x[, active, drop = FALSE], 2L, mu[active]), 2L, sigma[active], '/')
  ym <- mean(y)
  ys <- sd(y)
  model <- e1071::svm(z, (y - ym) / ys, type = 'eps-regression', kernel = kernel,
    cost = cost, gamma = gamma, epsilon = .1, scale = FALSE, cachesize = 256,
    tolerance = if (kernel == 'linear') .01 else .001)
  list(model = model, active = active, x_mean = mu, x_sd = sigma, y_mean = ym, y_sd = ys)
}

predict_program_svr <- function(f, x) {
  stopifnot(requireNamespace('e1071', quietly = TRUE))
  z <- sweep(sweep(x[, f$active, drop = FALSE], 2L, f$x_mean[f$active]), 2L, f$x_sd[f$active], '/')
  as.numeric(stats::predict(f$model, z)) * f$y_sd + f$y_mean
}

# Same coefficients predict absolute effects and differences between recipients
# or treatments. The contrast views are linear transformations, not new samples.
context_response_views <- function(x, y, metadata) {
  weights <- observation_weights(metadata)
  center_group <- function(mat, group) {
    totals <- rowsum(mat, group, reorder = FALSE)
    counts <- tabulate(match(group, rownames(totals)), nbins = nrow(totals))
    means <- totals / counts
    mat - means[match(group, rownames(means)), , drop = FALSE]
  }
  combined <- cbind(x, outcome = y)
  views <- list(combined,
    center_group(combined, paste(metadata$drug, metadata$dose)),
    center_group(combined, paste(metadata$depmap_id, metadata$dose)))
  stacked <- do.call(rbind, views)
  list(x = stacked[, seq_len(ncol(x)), drop = FALSE], y = stacked[, ncol(stacked)],
    weights = rep(weights / 3, 3L))
}

# Baseline is centered across distinct fitting recipients, while displacement
# is scaled without centering so all interaction features vanish at zero input.
context_response_design <- function(b, d, substrate, scaling) {
  bz <- sweep(sweep(b, 2L, scaling$b_mean, '-'), 2L, scaling$b_sd, '/')
  dz <- sweep(d, 2L, scaling$d_sd, '/')
  interactions <- do.call(cbind, lapply(seq_len(ncol(dz)), function(k) bz * dz[, k]))
  colnames(interactions) <- unlist(lapply(colnames(d), function(k)
    paste(colnames(b), k, sep = ':')), use.names = FALSE)
  cbind(SVR_gain = substrate, interactions)
}

# Self-contained decoder: reusable for any response scored in this frozen bank.
predict_context_response <- function(model, b, d) {
  substrate <- predict_program_svr(model$svr, d) - model$svr_zero
  x <- context_response_design(b, d, substrate, model$scaling)
  coefficients <- as.numeric(stats::coef(model$fit, s = model$lambda))[-1L]
  shared <- coefficients[1L] * substrate
  contextual <- as.numeric(x[, -1L, drop = FALSE] %*% coefficients[-1L])
  list(predicted = shared + contextual, shared = shared, contextual = contextual)
}

# Retrospective failure tracing. Replays frozen fits; never refits or selects
# features on the examples. Existing discordant cases retain their exact IDs.
audit_contextual_failures <- function(family) {
  root <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
  save_table <- function(x, name) data.table::fwrite(x,
    file.path(root, 'data', paste0('contextual_audit_', name, '.tsv')), sep = '\t')
  settings <- qs::qread(file.path(root, 'settings.qs'))
  prepared <- qs::qread(file.path(root, 'data/prepared_scores.qs'))
  bank <- qs::qread(file.path(root, 'models/program_bank.qs'))
  model <- qs::qread(file.path(root, 'models/contextual_actual.qs'))
  substrate <- qs::qread(file.path(root, 'models/nonlinear_RBF_SVR_D.qs'))
  panel <- qs::qread('/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
  labels <- data.table::as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
    depmap_id %in% rownames(prepared$baseline)]
  train <- labels[!depmap_id %in% settings$held_ids & !drug %in% settings$excluded_drugs]
  y <- log2(train$prism_abundance_pct / 100)
  d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  b <- prepared$baseline[train$depmap_id, , drop = FALSE]
  folds <- data.table::fread(file.path(root, 'data/nonlinear_training_folds.tsv'))
  fold <- folds$fold[match(train$depmap_id, folds$depmap_id)]
  stopifnot(!anyNA(fold), setequal(unique(train$depmap_id), model$training_ids))
  candidate_grid <- data.table::fread(file.path(root, 'data/nonlinear_validation.tsv'))
  candidate <- candidate_grid[model == 'RBF_SVR_D' & cost == substrate$selected$cost &
    gamma_multiplier == substrate$selected$gamma_multiplier, candidate]
  stopifnot(length(candidate) == 1L)
  penalties <- qs::qread(file.path(root, 'contextual_settings.qs'))$penalties
  stopifnot(max(abs(model$fit$lambda / penalties - 1)) < 1e-12)
  # Correlations explicitly separate recipients at one treatment from drugs
  # at one recipient/dose. Native spread is in percentage points of vehicle.
  describe_prediction <- function(md, observed, prediction, shared) {
    z <- data.table::copy(md[, .(depmap_id, drug, dose)])
    z[, `:=`(observed = observed, prediction = prediction, shared = shared,
      observed_pct = 100 * 2^observed, predicted_pct = 100 * 2^prediction)]
    z[, `:=`(yc = observed - mean(observed), pc = prediction - mean(prediction),
      rc = observed - shared - mean(observed - shared),
      cc = prediction - shared - mean(prediction - shared)), by = .(drug, dose)]
    within <- z[, .(r = cor_safe(observed, prediction)), by = .(depmap_id, dose)]
    spreads <- z[, .(observed_sd_pp = sd(observed_pct), predicted_sd_pp = sd(predicted_pct)), by = .(drug, dose)]
    data.table::data.table(across_cells_r = cor_safe(z$yc, z$pc),
      across_drugs_fixed_dose_r = median(within$r, na.rm = TRUE),
      context_residual_r = cor_safe(z$rc, z$cc),
      observed_sd_pp = mean(spreads$observed_sd_pp, na.rm = TRUE),
      predicted_sd_pp = mean(spreads$predicted_sd_pp, na.rm = TRUE),
      median_absolute_context_pp = median(abs(100 * 2^prediction - 100 * 2^shared)),
      sensitive = sum(observed <= -1), sensitive_called = sum(observed <= -1 & prediction <= -1))
  }
  metrics <- list()
  oof_predictions <- list()
  coefficient_tables <- list()
  for (k in 1:3) {
    message(family, ' | replaying context fold ', k)
    tr <- which(fold != k)
    saved <- qs::qread(file.path(root, 'models', sprintf('contextual_actual_fold%d.qs', k)))
    va <- saved$validation_rows
    inner_s <- rep(NA_real_, length(tr))
    for (j in 1:2) {
      inner <- qs::qread(file.path(root, 'models', sprintf('contextual_substrate_fold%d_inner%d.qs', k, j)))
      stopifnot(!length(intersect(inner$training_rows, va)), !length(intersect(inner$validation_rows, va)))
      iv <- inner$validation_rows
      zero <- predict_program_svr(inner$model, matrix(0, 1L, ncol(d)))
      inner_s[match(iv, tr)] <- predict_program_svr(inner$model, d[iv, , drop = FALSE]) - zero
    }
    outer <- qs::qread(file.path(root, 'models', sprintf('nonlinear_RBF_SVR_D_candidate%d_fold%d.qs', candidate, k)))
    validation_s <- predict_program_svr(outer$model, d[va, , drop = FALSE]) -
      predict_program_svr(outer$model, matrix(0, 1L, ncol(d)))
    stopifnot(all(is.finite(inner_s)))
    xx <- context_response_design(b[tr, ], d[tr, ], inner_s, saved$scaling)
    vx <- context_response_design(b[va, ], d[va, ], validation_s, saved$scaling)
    coefficients <- as.matrix(stats::coef(saved$fit, s = penalties))[-1L, , drop = FALSE]
    pp_train <- as.matrix(stats::predict(saved$fit, newx = xx, s = penalties))
    pp_validation <- as.matrix(stats::predict(saved$fit, newx = vx, s = penalties))
    replay <- qs::qread(file.path(root, 'models', sprintf('contextual_predictions_actual_fold%d.qs', k)))
    stopifnot(max(abs(pp_validation - replay$prediction)) < 1e-10)
    for (j in seq_along(penalties)) {
      for (cohort in c('Fitting cells', 'Excluded development cells')) {
        rows <- if (cohort == 'Fitting cells') tr else va
        pp <- if (cohort == 'Fitting cells') pp_train[, j] else pp_validation[, j]
        ss <- coefficients[1L, j] * if (cohort == 'Fitting cells') inner_s else validation_s
        metrics[[paste(k, j, cohort)]] <- cbind(data.table::data.table(fold = k,
          lambda = penalties[j], cohort = cohort), describe_prediction(train[rows], y[rows], pp, ss))
      }
      oof_predictions[[paste(k, j)]] <- cbind(train[va, .(observation_id, depmap_id, recipient, drug, dose)],
        data.table::data.table(lambda = penalties[j], predicted = pp_validation[, j],
          shared = coefficients[1L, j] * validation_s, prediction_source = 'Whole-cell development OOF'))
      coefficient_tables[[paste(k, j)]] <- data.table::data.table(fold = k, lambda = penalties[j],
        feature = rownames(coefficients)[-1L], coefficient = coefficients[-1L, j])
    }
  }
  metrics <- data.table::rbindlist(metrics)
  save_table(metrics, 'fold_metrics')
  coefficients <- data.table::rbindlist(coefficient_tables)
  stability <- coefficients[, .(positive_folds = sum(coefficient > 0),
    negative_folds = sum(coefficient < 0), coefficient_min = min(coefficient),
    coefficient_max = max(coefficient)), by = .(lambda, feature)]
  save_table(stability, 'coefficient_stability')
  oof <- data.table::rbindlist(oof_predictions)
  save_table(oof, 'development_predictions')
  # Baseline geometry: biological recipients, rather than drug rows, determine
  # independent context sample size. Report collinearity and expression-source agreement.
  ids <- model$training_ids
  bz <- sweep(sweep(prepared$baseline[ids, ], 2L, model$scaling$b_mean, '-'), 2L, model$scaling$b_sd, '/')
  eigenvalues <- svd(bz, nu = 0L, nv = 0L)$d^2
  geometry <- data.table::data.table(recipients = length(ids), outcome_rows = nrow(train),
    interaction_features = 1600L, baseline_rank = sum(eigenvalues > max(eigenvalues) * 1e-10),
    baseline_pc1_fraction = eigenvalues[1L] / sum(eigenvalues),
    baseline_effective_dimension = sum(eigenvalues)^2 / sum(eigenvalues^2),
    median_absolute_baseline_program_correlation = median(abs(cor(bz)[upper.tri(cor(bz))])))
  save_table(geometry, 'baseline_geometry')
  reference_root <- '/work/results/program_pilot/20_multidrug_baseline_response_function_MEK_holdout'
  cases <- data.table::fread(file.path(reference_root, 'data/diagnostic_geometry_example_resolution.tsv'))
  cases[, case := sprintf('C%02d', seq_len(.N))]
  case_ids <- unique(c(sub('__.*', '', cases$observation_i), sub('__.*', '', cases$observation_j)))
  endpoints <- unique(c(cases$observation_i, cases$observation_j))
  profiles <- list()
  vehicle_scores <- matrix(NA_real_, nrow(prepared$baseline), ncol(d), dimnames = dimnames(prepared$baseline))
  vehicle_gene_scores <- matrix(NA_real_, nrow(prepared$baseline), length(bank$genes),
    dimnames = list(rownames(prepared$baseline), bank$genes))
  for (id in rownames(vehicle_scores)) {
    molecular <- qs::qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
    stopifnot(length(molecular$baseline[[id]]) == length(molecular$genes),
      identical(colnames(molecular$response), molecular$genes), all(bank$genes %in% molecular$genes))
    vehicle_baseline <- stats::setNames(molecular$baseline[[id]], molecular$genes)
    vehicle_gene_scores[id, ] <- vehicle_baseline[bank$genes]
    vehicle_scores[id, ] <- as.numeric(vehicle_baseline[bank$genes] %*% bank$weights)
    if (id %in% case_ids) {
      obs <- intersect(endpoints, rownames(molecular$response))
      profiles[[id]] <- list(response = molecular$response[obs, , drop = FALSE],
        null = molecular$null[[id]], baseline = vehicle_baseline,
        metadata = molecular$metadata[match(obs, observation_id)])
      expected <- prepared$measured[match(obs, prepared$metadata$observation_id), , drop = FALSE]
      stopifnot(max(abs(as.matrix(profiles[[id]]$response[, bank$genes, drop = FALSE] %*% bank$weights) - expected)) < 1e-12)
    }
    rm(molecular)
  }
  agreement <- data.table::data.table(program = colnames(d),
    across_cell_ccle_vehicle_r = vapply(seq_len(ncol(d)), function(k)
      cor_safe(prepared$baseline[ids, k], vehicle_scores[ids, k]), numeric(1)))
  save_table(agreement, 'baseline_source_agreement')
  # Inspect named baseline genes as measurements, without claiming causality.
  context_genes <- c('ERBB2', 'EGFR', 'MET', 'ALK', 'ROS1', 'BRAF', 'KRAS', 'NRAS', 'STK11',
    'PIK3CA', 'MAP2K1', 'MAP2K2', 'MTOR', 'PSMB5',
    'BCL2', 'BCL2L1', 'MCL1', 'BAX', 'BAK1', 'BCL2L11', 'BID', 'BBC3', 'PMAIP1', 'CASP3', 'CASP8',
    'APAF1', 'XIAP', 'BIRC5', 'MKI67', 'PCNA', 'E2F1', 'MYC', 'CDKN1A', 'CDKN2A', 'RB1', 'TP53')
  expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
  header <- names(data.table::fread(expression_path, nrows = 0L))
  symbols <- sub(' \\(ENSG.*', '', header[-1L])
  selected <- which(symbols %in% unique(c(bank$genes, context_genes)))
  expression_table <- data.table::fread(expression_path, select = c(1L, selected + 1L))
  expression <- as.matrix(expression_table[, -1L])
  dimnames(expression) <- list(expression_table[[1L]], symbols[selected])
  context_genes <- intersect(context_genes, colnames(expression))
  source_agreement <- data.table::data.table(gene = bank$genes,
    across_cell_ccle_vehicle_r = vapply(bank$genes, function(gene)
      cor_safe(expression[ids, gene], vehicle_gene_scores[ids, gene]), numeric(1)))
  save_table(source_agreement, 'baseline_gene_source_agreement')
  # Per-treatment functional coverage reveals genuine extrapolation cases,
  # rather than interpreting absence of sensitive fitting examples as calibration failure.
  coverage <- train[, .(training_recipients = data.table::uniqueN(depmap_id),
    training_sensitive = sum(prism_abundance_pct <= 50), training_min_pct = min(prism_abundance_pct),
    training_max_pct = max(prism_abundance_pct)), by = .(drug, dose)]
  held_coverage <- merge(labels[depmap_id %in% settings$held_ids], coverage,
    by = c('drug', 'dose'), all.x = TRUE, sort = FALSE)
  held_coverage[, `:=`(outside_training_range = prism_abundance_pct < training_min_pct |
      prism_abundance_pct > training_max_pct,
    sensitive_without_training_sensitive = prism_abundance_pct <= 50 & training_sensitive == 0)]
  save_table(held_coverage, 'held_treatment_coverage')
  mutation_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv'
  mutation_header <- names(data.table::fread(mutation_path, nrows = 0L))
  mutation_symbols <- sub(' \\(.*', '', mutation_header[-1L])
  driver_genes <- c('BRAF', 'KRAS', 'NRAS', 'EGFR', 'ERBB2', 'PIK3CA', 'STK11', 'TP53', 'RB1', 'PTEN')
  mutation_selected <- which(mutation_symbols %in% driver_genes)
  mutation <- data.table::fread(mutation_path, select = c(1L, mutation_selected + 1L))
  data.table::setnames(mutation, c('depmap_id', mutation_symbols[mutation_selected]))
  stopifnot(all(driver_genes %in% names(mutation)), all(case_ids %in% mutation$depmap_id))
  mutation_long <- data.table::melt(mutation[depmap_id %in% rownames(prepared$baseline)],
    id.vars = 'depmap_id', variable.name = 'gene', value.name = 'hotspot_flag')
  mutation_long[, `:=`(cell = panel$cell_map$depmap_name[match(depmap_id, panel$cell_map$depmap_id)],
    cohort = ifelse(depmap_id %in% ids, 'Training cells',
      ifelse(depmap_id %in% settings$held_ids, 'Held cells', 'No functional training data')))]
  stopifnot(!anyNA(mutation_long$hotspot_flag), all(mutation_long$hotspot_flag %in% 0:1))
  save_table(mutation_long, 'driver_hotspots')
  save_table(mutation_long[, .(recipients = .N, positive = sum(hotspot_flag)), by = .(gene, cohort)],
    'driver_coverage')
  final <- data.table::fread(file.path(root, 'data/contextual_components.tsv'))[variant == 'actual']
  final[, prediction_source := 'Held-cell frozen prediction']
  prediction_lookup <- data.table::rbindlist(list(oof[lambda == model$lambda],
    final[, .(observation_id, depmap_id, recipient, drug, dose, predicted, shared, prediction_source)]), fill = TRUE)
  stopifnot(!anyDuplicated(prediction_lookup$observation_id))
  case_summary <- program_values <- gene_values <- baseline_gene_values <- neighbors <- list()
  for (k in seq_len(nrow(cases))) {
    pair <- cases[k]
    obs <- c(pair$observation_i, pair$observation_j)
    cell_ids <- sub('__.*', '', obs)
    left <- profiles[[cell_ids[1L]]]
    right <- profiles[[cell_ids[2L]]]
    common <- intersect(colnames(left$response), colnames(right$response))
    gi <- left$response[obs[1L], common]
    gj <- right$response[obs[2L], common]
    mi <- match(obs, prepared$metadata$observation_id)
    di <- prepared$measured[mi, , drop = FALSE]
    bi <- prepared$baseline[cell_ids, , drop = FALSE]
    vi <- vehicle_scores[cell_ids, , drop = FALSE]
    ni <- as.matrix(left$null[, bank$genes, drop = FALSE] %*% bank$weights)
    nj <- as.matrix(right$null[, bank$genes, drop = FALSE] %*% bank$weights)
    null_pairs <- data.table::CJ(a = seq_len(nrow(ni)), z = seq_len(nrow(nj)))
    if (cell_ids[1L] == cell_ids[2L]) null_pairs <- null_pairs[a != z]
    null_difference <- sqrt(rowMeans((ni[null_pairs$a, ] - nj[null_pairs$z, ])^2))
    response_difference <- sqrt(mean((di[1L, ] - di[2L, ])^2))
    gene_sd <- sqrt(apply(left$null[, common, drop = FALSE], 2L, var) +
      apply(right$null[, common, drop = FALSE], 2L, var))
    gene_difference <- gi - gj
    gene_z <- gene_difference / pmax(gene_sd, 1e-12)
    pp <- prediction_lookup[match(obs, observation_id)]
    outcomes <- labels[match(obs, observation_id)]
    stopifnot(!anyNA(outcomes$observation_id))
    case_summary[[k]] <- data.table::data.table(case = pair$case, comparison = pair$comparison,
      cohort = pair$cohort, example_type = pair$example_type, cell_i = pair$cell_i, cell_j = pair$cell_j,
      drug_i = pair$drug_i, drug_j = pair$drug_j, dose = pair$dose,
      observed_i = outcomes$prism_abundance_pct[1L], observed_j = outcomes$prism_abundance_pct[2L],
      predicted_i = 100 * 2^pp$predicted[1L], predicted_j = 100 * 2^pp$predicted[2L],
      shared_i = 100 * 2^pp$shared[1L], shared_j = 100 * 2^pp$shared[2L],
      context_i_pp = 100 * (2^pp$predicted[1L] - 2^pp$shared[1L]),
      context_j_pp = 100 * (2^pp$predicted[2L] - 2^pp$shared[2L]),
      prediction_source_i = pp$prediction_source[1L], prediction_source_j = pp$prediction_source[2L],
      program_difference = response_difference, vehicle_difference_median = median(null_difference),
      vehicle_difference_95 = quantile(null_difference, .95),
      program_difference_vehicle_percentile = mean(null_difference <= response_difference),
      raw_gene_pearson = cor(gi, gj), genes_above_two_vehicle_sd = sum(abs(gene_z) > 2),
      assayed_common_genes = length(common),
      ccle_baseline_rms_z = sqrt(mean(((bi[1L, ] - bi[2L, ]) / model$scaling$b_sd)^2)),
      vehicle_baseline_rms = sqrt(mean((vi[1L, ] - vi[2L, ])^2)),
      prism_status_i = outcomes$prism_status[1L], prism_status_j = outcomes$prism_status[2L],
      prism_screen_i = outcomes$prism_screen[1L], prism_screen_j = outcomes$prism_screen[2L],
      prism_lower_pct_i = 100 * 2^outcomes$prism_lower_log2[1L],
      prism_upper_pct_i = 100 * 2^outcomes$prism_upper_log2[1L],
      prism_lower_pct_j = 100 * 2^outcomes$prism_lower_log2[2L],
      prism_upper_pct_j = 100 * 2^outcomes$prism_upper_log2[2L],
      almanac_i = outcomes$almanac_abundance_pct[1L], almanac_j = outcomes$almanac_abundance_pct[2L])
    stopifnot(all(is.finite(c(pp$predicted, pp$shared))))
    program_values[[k]] <- data.table::data.table(case = pair$case, program = colnames(d),
      displacement_i = di[1L, ], displacement_j = di[2L, ],
      baseline_ccle_i = bi[1L, ], baseline_ccle_j = bi[2L, ],
      baseline_ccle_z_i = (bi[1L, ] - model$scaling$b_mean) / model$scaling$b_sd,
      baseline_ccle_z_j = (bi[2L, ] - model$scaling$b_mean) / model$scaling$b_sd,
      baseline_ccle_difference_z = (bi[1L, ] - bi[2L, ]) / model$scaling$b_sd,
      baseline_vehicle_i = vi[1L, ], baseline_vehicle_j = vi[2L, ])
    top <- order(abs(gene_z), decreasing = TRUE)[seq_len(min(50L, length(common)))]
    gene_values[[k]] <- data.table::data.table(case = pair$case, gene = common[top],
      displacement_i = gi[top], displacement_j = gj[top], difference = gene_difference[top],
      difference_vehicle_sd = gene_z[top], in_bank = common[top] %in% bank$genes,
      program = colnames(d)[bank$partition[match(common[top], bank$genes)]])
    context_genes_pair <- intersect(context_genes, common)
    baseline_gene_values[[k]] <- data.table::data.table(case = pair$case, gene = context_genes_pair,
      ccle_i = expression[cell_ids[1L], context_genes_pair], ccle_j = expression[cell_ids[2L], context_genes_pair],
      vehicle_i = left$baseline[context_genes_pair], vehicle_j = right$baseline[context_genes_pair],
      in_bank = context_genes_pair %in% bank$genes)
    # Same-treatment neighbors expose whether response- or baseline-nearest
    # training recipients actually have compatible abundance. Never fit here.
    for (side in 1:2) {
      md <- outcomes[side]
      ii <- which(train$drug == md$drug & train$dose == md$dose & train$depmap_id != md$depmap_id)
      if (!length(ii)) next
      dd <- sqrt(rowMeans(sweep(sweep(d[ii, , drop = FALSE], 2L, di[side, ], '-'),
        2L, model$scaling$d_sd, '/')^2))
      bb <- sqrt(rowMeans(sweep(sweep(b[ii, , drop = FALSE], 2L, bi[side, ], '-'),
        2L, model$scaling$b_sd, '/')^2))
      for (space in c('Displacement', 'Baseline', 'Both')) {
        distance <- switch(space, Displacement = dd, Baseline = bb, Both = sqrt(dd^2 + bb^2))
        nn <- order(distance)[seq_len(min(3L, length(ii)))]
        neighbors[[paste(k, side, space)]] <- data.table::data.table(case = pair$case, side = side,
          query_cell = md$recipient, query_abundance_pct = md$prism_abundance_pct, space = space,
          neighbor_cell = train$recipient[ii[nn]], neighbor_abundance_pct = train$prism_abundance_pct[ii[nn]],
          displacement_distance = dd[nn], baseline_distance = bb[nn], rank = seq_along(nn))
      }
    }
  }
  save_table(data.table::rbindlist(case_summary), 'cases')
  save_table(data.table::rbindlist(program_values), 'case_programs')
  save_table(data.table::rbindlist(gene_values), 'case_response_genes')
  save_table(data.table::rbindlist(baseline_gene_values), 'case_baseline_genes')
  save_table(data.table::rbindlist(neighbors), 'case_training_neighbors')
  save_table(data.table::data.table(program = colnames(d), genes = as.integer(colSums(bank$weights > 0))), 'program_sizes')
  focal_driver <- labels[drug == 'vemurafenib' & dose == .5,
    .(depmap_id, recipient, drug, dose, prism_abundance_pct, almanac_abundance_pct)]
  focal_driver <- merge(focal_driver, mutation_long[gene == 'BRAF', .(depmap_id, hotspot_flag, cohort)],
    by = 'depmap_id', all.x = TRUE, sort = FALSE)
  save_table(focal_driver, 'vemurafenib_braf_context')
  provenance <- c('/work/src/program_pilot/nonlinear_response.R',
    '/work/src/program_pilot/family_transfer_plots.R', file.path(root, 'models/contextual_actual.qs'),
    file.path(reference_root, 'data/diagnostic_geometry_example_resolution.tsv'), expression_path, mutation_path)
  save_table(data.table::data.table(path = provenance,
    sha256 = vapply(provenance, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'provenance')
  source('/work/src/program_pilot/family_transfer_plots.R')
  plot_contextual_audit(root)
  print(geometry)
  print(data.table::rbindlist(case_summary)[example_type == 'Similar programs, different outcomes',
    .(case, cell_i, cell_j, drug_i, drug_j, observed_i, observed_j, predicted_i, predicted_j,
      context_i_pp, context_j_pp, genes_above_two_vehicle_sd, program_difference_vehicle_percentile)])
}

run_contextual_response <- function(family) {
  # Small predeclared grid; fitting and selection balance three comparison axes.
  penalty_grid <- c(1000, 100, 10, 1, .1, .01)
  root <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
  data_dir <- file.path(root, 'data')
  model_dir <- file.path(root, 'models')
  write_table <- function(x, name) data.table::fwrite(x,
    file.path(data_dir, paste0('contextual_', name, '.tsv')), sep = '\t')
  status <- function(stage) {
    message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', family, ' | ', stage)
    data.table::fwrite(data.table::data.table(stage = stage, time = as.character(Sys.time()), pid = Sys.getpid()),
      file.path(root, 'contextual_status.tsv'), sep = '\t')
  }
  commit <- function(object, path) {
    temporary <- tempfile(tmpdir = dirname(path))
    qs::qsave(object, temporary, preset = 'high')
    stopifnot(file.rename(temporary, path))
  }
  settings <- qs::qread(file.path(root, 'settings.qs'))
  prepared <- qs::qread(file.path(data_dir, 'prepared_scores.qs'))
  base_path <- file.path(model_dir, 'nonlinear_RBF_SVR_D.qs')
  base <- qs::qread(base_path)
  panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
  panel <- qs::qread(panel_path)
  old_manifest <- data.table::fread(file.path(data_dir, 'nonlinear_input_hashes.tsv'))
  unchanged <- old_manifest[!grepl('/nonlinear_response.R$', path)]
  stopifnot(all(unchanged$sha256 == unname(vapply(unchanged$path, digest::digest,
    character(1), file = TRUE, algo = 'sha256'))))
  sources <- c('/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R',
    panel_path, base_path, file.path(root, 'settings.qs'), file.path(data_dir, 'prepared_scores.qs'))
  manifest <- data.table::data.table(path = sources,
    sha256 = unname(vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256')))
  signature <- digest::digest(list(manifest, penalty_grid, 'three_equal_views_v1'), algo = 'sha256')
  write_table(manifest, 'input_hashes')
  capture.output(sessionInfo(), file = file.path(root, 'contextual_sessionInfo.txt'))
  labels <- data.table::as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
    depmap_id %in% rownames(prepared$baseline)]
  train <- labels[!depmap_id %in% settings$held_ids & !drug %in% settings$excluded_drugs]
  query <- labels[depmap_id %in% settings$held_ids, .(observation_id, depmap_id, recipient, drug, dose)]
  y <- log2(train$prism_abundance_pct / 100)
  d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  dt <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  b <- prepared$baseline[train$depmap_id, , drop = FALSE]
  bt <- prepared$baseline[query$depmap_id, , drop = FALSE]
  programs <- colnames(d)
  stopifnot(ncol(d) == 40L, all(is.finite(d)), all(is.finite(dt)), all(is.finite(b)),
    setequal(unique(train$depmap_id), base$training_ids), !anyDuplicated(train$observation_id),
    !any(train$depmap_id %in% settings$held_ids), !any(train$drug %in% settings$excluded_drugs))
  folds <- data.table::fread(file.path(data_dir, 'nonlinear_training_folds.tsv'))
  fold <- folds$fold[match(train$depmap_id, folds$depmap_id)]
  stopifnot(!anyNA(fold), identical(sort(unique(fold)), 1:3))
  validation <- data.table::fread(file.path(data_dir, 'nonlinear_validation.tsv'))
  selected <- validation[model == 'RBF_SVR_D' & cost == base$selected$cost &
    gamma_multiplier == base$selected$gamma_multiplier]
  stopifnot(nrow(selected) == 1L)
  candidate <- selected$candidate
  gamma <- base$selected$gamma_multiplier / ncol(d)
  svr_zero <- predict_program_svr(base$model, matrix(0, 1L, ncol(d)))
  # Treatment means use molecular observations from fitting recipients only,
  # including unlabelled observations; no evaluation outcomes or recipients.
  mean_input <- function(ids, metadata) {
    ii <- which(prepared$metadata$depmap_id %in% ids)
    keys <- paste(prepared$metadata$drug[ii], prepared$metadata$dose[ii])
    totals <- rowsum(prepared$measured[ii, , drop = FALSE], keys, reorder = FALSE)
    counts <- tabulate(match(keys, rownames(totals)), nbins = nrow(totals))
    requested <- match(paste(metadata$drug, metadata$dose), rownames(totals))
    stopifnot(!length(intersect(ids, settings$held_ids)))
    # Rare treatments can lack any other recipient in a nested split. Missing
    # means are explicitly excluded from that comparator, never imputed.
    result <- (totals / counts)[requested, , drop = FALSE]
    attr(result, 'supported') <- !is.na(requested)
    result
  }
  scaling_for <- function(rows) {
    ids <- unique(train$depmap_id[rows])
    bb <- prepared$baseline[ids, , drop = FALSE]
    scaling <- list(b_mean = colMeans(bb), b_sd = apply(bb, 2L, sd), d_sd = apply(d[rows, ], 2L, sd))
    stopifnot(all(scaling$b_sd > 0), all(scaling$d_sd > 0))
    scaling
  }
  fit_context_path <- function(x, yy, md) {
    views <- context_response_views(x, yy, md)
    glmnet::glmnet(views$x, views$y, family = 'gaussian', alpha = 0,
      lambda = penalty_grid, weights = views$weights, intercept = FALSE, standardize = FALSE,
      penalty.factor = c(0, rep(1, ncol(x) - 1L)), thresh = 1e-6)
  }
  validation_score <- function(prediction, yy, md) {
    views <- context_response_views(as.matrix(prediction), yy, md)
    colSums(abs(views$x - views$y) * views$weights) / sum(views$weights)
  }
  oof <- list(actual = rep(NA_real_, nrow(train)), average = rep(NA_real_, nrow(train)))
  mean_oof <- matrix(NA_real_, nrow(train), ncol(d), dimnames = dimnames(d))
  losses <- list(actual = rep(0, length(penalty_grid)), average = rep(0, length(penalty_grid)))
  cv_tables <- list()
  commit(list(signature = signature, family = family, penalties = penalty_grid,
    frozen_svr_parameters = base$selected, training_ids = base$training_ids,
    held_ids = settings$held_ids, excluded_drugs = settings$excluded_drugs,
    conditioning = 'all 1600 cross-program products; no intercept; anchored SVR with unpenalized gain',
    fitting = 'equal absolute, fixed-drug/dose and fixed-cell/dose contrast views',
    validation = 'three whole-cell folds; nested two-fold whole-cell substrate cross-fitting',
    evaluation = 'retrospective; previously revealed outcomes'), file.path(root, 'contextual_settings.qs'))
  for (k in 1:3) {
    status(paste('context validation fold', k, '/3'))
    tr <- which(fold != k)
    va <- which(fold == k)
    fold_base <- qs::qread(file.path(model_dir,
      sprintf('nonlinear_RBF_SVR_D_candidate%d_fold%d.qs', candidate, k)))
    stopifnot(identical(fold_base$signature, base$signature), identical(fold_base$validation_rows, va))
    raw <- predict_program_svr(fold_base$model, d[va, , drop = FALSE])
    stopifnot(max(abs(raw - fold_base$prediction)) < 1e-10,
      max(abs(raw - base$oof[va])) < 1e-10)
    zero <- predict_program_svr(fold_base$model, matrix(0, 1L, ncol(d)))
    vd_mean <- mean_input(unique(train$depmap_id[tr]), train[va])
    oof$actual[va] <- raw - zero
    supported <- attr(vd_mean, 'supported')
    oof$average[va[supported]] <- predict_program_svr(fold_base$model, vd_mean[supported, , drop = FALSE]) - zero
    mean_oof[va, ] <- vd_mean
    # Nested substrates prevent validation recipients influencing contextual
    # training through the shared SVR. Each completed fit is saved immediately.
    ids <- sort(unique(train$depmap_id[tr]))
    inner_map <- setNames(rep(1:2, length.out = length(ids)), ids)
    inner_s <- list(actual = rep(NA_real_, length(tr)), average = rep(NA_real_, length(tr)))
    inner_mean <- matrix(NA_real_, length(tr), ncol(d), dimnames = list(NULL, programs))
    for (j in 1:2) {
      it <- tr[inner_map[train$depmap_id[tr]] != j]
      iv <- tr[inner_map[train$depmap_id[tr]] == j]
      path <- file.path(model_dir, sprintf('contextual_substrate_fold%d_inner%d.qs', k, j))
      if (file.exists(path)) {
        saved <- qs::qread(path)
        stopifnot(identical(saved$signature, signature), identical(saved$training_rows, it),
          identical(saved$validation_rows, iv))
        f <- saved$model
      } else {
        status(paste('nested substrate', k, j))
        f <- fit_program_svr(d[it, , drop = FALSE], y[it], 'radial', base$selected$cost, gamma)
        commit(list(signature = signature, training_rows = it, validation_rows = iv, model = f), path)
      }
      ii <- match(iv, tr)
      dm <- mean_input(unique(train$depmap_id[it]), train[iv])
      zero_inner <- predict_program_svr(f, matrix(0, 1L, ncol(d)))
      inner_s$actual[ii] <- predict_program_svr(f, d[iv, , drop = FALSE]) - zero_inner
      supported <- attr(dm, 'supported')
      inner_s$average[ii[supported]] <- predict_program_svr(f, dm[supported, , drop = FALSE]) - zero_inner
      inner_mean[ii, ] <- dm
    }
    scaling <- scaling_for(tr)
    for (variant in c('actual', 'average')) {
      path <- file.path(model_dir, sprintf('contextual_%s_fold%d.qs', variant, k))
      input_train <- if (variant == 'actual') d[tr, , drop = FALSE] else inner_mean
      input_validation <- if (variant == 'actual') d[va, , drop = FALSE] else vd_mean
      fit_rows <- which(complete.cases(input_train) & is.finite(inner_s[[variant]]))
      validation_rows <- which(complete.cases(input_validation) & is.finite(oof[[variant]][va]))
      write_table(data.table::data.table(variant = variant, fold = k,
        training_rows = length(fit_rows), excluded_training_rows = length(tr) - length(fit_rows),
        validation_rows = length(validation_rows), excluded_validation_rows = length(va) - length(validation_rows)),
        paste0('coverage_', variant, '_fold', k))
      xx <- context_response_design(b[tr[fit_rows], ], input_train[fit_rows, , drop = FALSE],
        inner_s[[variant]][fit_rows], scaling)
      vx <- context_response_design(b[va[validation_rows], ], input_validation[validation_rows, , drop = FALSE],
        oof[[variant]][va[validation_rows]], scaling)
      if (file.exists(path)) {
        saved <- qs::qread(path)
        stopifnot(identical(saved$signature, signature), identical(saved$validation_rows, va[validation_rows]))
        fit <- saved$fit
      } else {
        status(paste('ridge context', variant, 'fold', k))
        fit <- fit_context_path(xx, y[tr[fit_rows]], train[tr[fit_rows]])
        commit(list(signature = signature, fit = fit, scaling = scaling, validation_rows = va[validation_rows]), path)
      }
      pp <- as.matrix(stats::predict(fit, newx = vx, s = penalty_grid))
      score <- validation_score(pp, y[va[validation_rows]], train[va[validation_rows]])
      losses[[variant]] <- losses[[variant]] + score / 3
      cv_tables[[paste(variant, k)]] <- data.table::data.table(variant = variant, fold = k,
        lambda = penalty_grid, balanced_absolute_log2_error = score)
      write_table(data.table::rbindlist(cv_tables), 'validation')
      commit(list(signature = signature, validation_rows = va[validation_rows], prediction = pp),
        file.path(model_dir, sprintf('contextual_predictions_%s_fold%d.qs', variant, k)))
    }
  }
  stopifnot(all(is.finite(oof$actual)), all(is.finite(oof$average[complete.cases(mean_oof)])))
  scaling <- scaling_for(seq_len(nrow(train)))
  query_mean <- mean_input(base$training_ids, query)
  stopifnot(all(is.finite(query_mean)))
  mean_training_rows <- which(complete.cases(mean_oof) & is.finite(oof$average))
  write_table(train[-mean_training_rows, .(observation_id, depmap_id, drug, dose)], 'average_training_exclusions')
  predictions <- list()
  components <- list()
  for (variant in c('actual', 'average')) {
    status(paste('final context', variant))
    input <- if (variant == 'actual') d else mean_oof
    fit_rows <- if (variant == 'actual') seq_len(nrow(train)) else mean_training_rows
    xx <- context_response_design(b[fit_rows, ], input[fit_rows, , drop = FALSE], oof[[variant]][fit_rows], scaling)
    best <- penalty_grid[which.min(losses[[variant]])]
    path <- file.path(model_dir, paste0('contextual_', variant, '.qs'))
    if (file.exists(path)) {
      model <- qs::qread(path)
      stopifnot(identical(model$signature, signature), identical(model$lambda, best))
    } else {
      model <- list(signature = signature, fit = fit_context_path(xx, y[fit_rows], train[fit_rows]), lambda = best,
        scaling = scaling, svr = base$model, svr_zero = svr_zero, training_ids = base$training_ids)
      commit(model, path)
    }
    result <- predict_context_response(model, bt, if (variant == 'actual') dt else query_mean)
    predictions[[paste0('context_', variant)]] <- result$predicted
    components[[variant]] <- result
    zero_check <- predict_context_response(model, bt, dt * 0)$predicted
    stopifnot(max(abs(zero_check)) < 1e-12, all(is.finite(result$predicted)))
    coefficients <- as.numeric(stats::coef(model$fit, s = best))[-1L]
    coefficient_table <- data.table::data.table(baseline_program = rep(programs, times = length(programs)),
      displaced_program = rep(programs, each = length(programs)), coefficient = coefficients[-1L])
    write_table(coefficient_table, paste0('coefficients_', variant))
    write_table(data.table::data.table(variant = variant, lambda = best, svr_gain = coefficients[1L],
      zero_input_max_abs = max(abs(zero_check)), max_abs_context_coefficient = max(abs(coefficients[-1L]))),
      paste0('settings_', variant))
  }
  # Calibrate the anchored response-only comparator on the same three views.
  for (variant in c('actual', 'average')) {
    fit_rows <- if (variant == 'actual') seq_len(nrow(train)) else mean_training_rows
    views <- context_response_views(matrix(oof[[variant]][fit_rows], ncol = 1L), y[fit_rows], train[fit_rows])
    gain <- sum(views$weights * views$x[, 1L] * views$y) / sum(views$weights * views$x[, 1L]^2)
    vv <- predict_program_svr(base$model, if (variant == 'actual') dt else query_mean) - svr_zero
    predictions[[paste0('anchored_', variant)]] <- gain * vv
    commit(list(signature = signature, svr = base$model, svr_zero = svr_zero, gain = gain),
      file.path(model_dir, paste0('contextual_anchored_', variant, '.qs')))
  }
  prior <- data.table::fread(file.path(data_dir, 'nonlinear_predictions.tsv.gz'))[
    model == 'calibrated_RBF_SVR_D']
  stopifnot(setequal(prior$observation_id, query$observation_id))
  predictions$previous_SVR <- prior$predicted[match(query$observation_id, prior$observation_id)]
  # Strong direct-baseline reference for known treatments. This diagnostic has
  # treatment intercepts and drug-specific baseline slopes, unlike the decoder.
  drug_levels <- sort(unique(train$drug))
  condition_levels <- sort(unique(paste(train$drug, train$dose)))
  baseline_design <- function(bb, md, scaling) {
    bz <- sweep(sweep(bb, 2L, scaling$b_mean, '-'), 2L, scaling$b_sd, '/')
    n <- nrow(bz)
    p <- ncol(bz)
    drug_index <- match(md$drug, drug_levels)
    condition_index <- match(paste(md$drug, md$dose), condition_levels)
    stopifnot(!anyNA(drug_index), !anyNA(condition_index))
    slopes <- Matrix::sparseMatrix(i = rep(seq_len(n), times = p),
      j = (rep(drug_index, times = p) - 1L) * p + rep(seq_len(p), each = n),
      x = as.vector(bz), dims = c(n, p * length(drug_levels)))
    intercepts <- Matrix::sparseMatrix(i = seq_len(n), j = condition_index, x = 1,
      dims = c(n, length(condition_levels)))
    cbind(slopes, intercepts)
  }
  baseline_loss <- rep(0, length(penalty_grid))
  for (k in 1:3) {
    tr <- which(fold != k)
    va <- which(fold == k & paste(train$drug, train$dose) %in%
      paste(train$drug[tr], train$dose[tr]))
    sc <- scaling_for(tr)
    xx <- baseline_design(b[tr, ], train[tr], sc)
    vx <- baseline_design(b[va, ], train[va], sc)
    f <- glmnet::glmnet(xx, y[tr], alpha = 0, lambda = penalty_grid,
      weights = observation_weights(train[tr]), intercept = TRUE, standardize = TRUE)
    pp <- as.matrix(stats::predict(f, newx = vx, s = penalty_grid))
    baseline_loss <- baseline_loss + validation_score(pp, y[va], train[va]) / 3
    commit(list(signature = signature, fit = f, scaling = sc, validation_rows = va),
      file.path(model_dir, sprintf('contextual_baseline_fold%d.qs', k)))
    write_table(data.table::data.table(fold = k, validation_rows = length(va),
      unsupported_rows = sum(fold == k) - length(va)), paste0('baseline_coverage_fold', k))
  }
  baseline_fit <- glmnet::glmnet(baseline_design(b, train, scaling), y, alpha = 0,
    lambda = penalty_grid, weights = observation_weights(train), intercept = TRUE, standardize = TRUE)
  baseline_lambda <- penalty_grid[which.min(baseline_loss)]
  baseline_rows <- which(paste(query$drug, query$dose) %in% condition_levels)
  baseline_prediction <- as.numeric(stats::predict(baseline_fit,
    newx = baseline_design(bt[baseline_rows, ], query[baseline_rows], scaling), s = baseline_lambda))
  commit(list(signature = signature, fit = baseline_fit, lambda = baseline_lambda, scaling = scaling,
    drug_levels = drug_levels, condition_levels = condition_levels, training_ids = base$training_ids),
    file.path(model_dir, 'contextual_baseline_reference.qs'))
  pred <- data.table::rbindlist(lapply(names(predictions), function(name) {
    tab <- data.table::copy(query)
    tab[, `:=`(model = name, predicted = predictions[[name]])]
    tab
  }))
  baseline_table <- data.table::copy(query[baseline_rows])
  baseline_table[, `:=`(model = 'baseline_reference', predicted = baseline_prediction)]
  pred <- rbind(pred, baseline_table)
  write_table(query[-baseline_rows], 'baseline_unsupported_treatments')
  # Write and hash outcome-free predictions before attaching evaluation labels.
  prediction_path <- file.path(data_dir, 'contextual_predictions.tsv.gz')
  data.table::fwrite(pred, prediction_path, sep = '\t')
  write_table(data.table::data.table(path = prediction_path,
    sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')
  component_table <- data.table::rbindlist(lapply(names(components), function(name) {
    tab <- data.table::copy(query)
    tab[, `:=`(variant = name, shared = components[[name]]$shared,
      contextual = components[[name]]$contextual, predicted = components[[name]]$predicted)]
    tab
  }))
  write_table(component_table, 'components')
  status('predictions saved; reporting and plotting')
  source('/work/src/program_pilot/family_transfer_plots.R')
  report_contextual_response(root, panel)
  status('COMPLETE: contextual decoder')
}

run_nonlinear_response <- function(family) {
  # Fixed small grid; all selection is inside development recipients.
  cv_folds <- 3L
  svm_grid <- data.table::CJ(cost = c(1, 10), gamma_multiplier = c(.25, 1))
  root <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
  data_dir <- file.path(root, 'data')
  model_dir <- file.path(root, 'models')
  status_path <- file.path(root, 'nonlinear_status.tsv')
  status <- function(stage) {
    message(format(Sys.time()), ' | ', family, ' | ', stage)
    data.table::fwrite(data.table::data.table(stage = stage, time = as.character(Sys.time()), pid = Sys.getpid()),
      status_path, sep = '\t')
  }
  commit <- function(object, path) {
    temporary <- tempfile(tmpdir = dirname(path))
    qs::qsave(object, temporary, preset = 'high')
    stopifnot(file.rename(temporary, path))
  }
  write_table <- function(x, name) data.table::fwrite(x,
    file.path(data_dir, paste0('nonlinear_', name, '.tsv')), sep = '\t')
  settings <- qs::qread(file.path(root, 'settings.qs'))
  prepared <- qs::qread(file.path(data_dir, 'prepared_scores.qs'))
  original <- qs::qread(file.path(model_dir, 'functional_decoder.qs'))
  panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
  panel <- qs::qread(panel_path)
  sources <- c('/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R',
    panel_path, file.path(root, c('settings.qs', 'data/prepared_scores.qs', 'models/functional_decoder.qs')))
  manifest <- data.table::data.table(path = sources,
    sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256'))
  signature <- digest::digest(list(manifest, cv_folds, svm_grid, lambda_grid,
    epsilon = .1, seed = 20261010L), algo = 'sha256')
  write_table(manifest, 'input_hashes')
  capture.output(sessionInfo(), file = file.path(root, 'nonlinear_sessionInfo.txt'))
  labels <- data.table::as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
    depmap_id %in% rownames(prepared$baseline)]
  labels[, observed := log2(prism_abundance_pct / 100)]
  train <- labels[!depmap_id %in% settings$held_ids & !drug %in% settings$excluded_drugs]
  test <- labels[depmap_id %in% settings$held_ids]
  stopifnot(setequal(unique(train$depmap_id), original$training_ids),
    !any(train$drug %in% settings$excluded_drugs), !any(train$depmap_id %in% settings$held_ids),
    !anyDuplicated(train$observation_id), !anyDuplicated(test$observation_id))
  d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  dt <- prepared$measured[match(test$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  b <- prepared$baseline[train$depmap_id, , drop = FALSE]
  bt <- prepared$baseline[test$depmap_id, , drop = FALSE]
  colnames(b) <- paste0('B_', colnames(b))
  colnames(bt) <- colnames(b)
  features <- list(D = d, DB = cbind(d, b))
  queries <- list(D = dt, DB = cbind(dt, bt))
  b_old <- sweep(sweep(b, 2L, original$baseline_mean), 2L, original$baseline_sd, '/')
  bt_old <- sweep(sweep(bt, 2L, original$baseline_mean), 2L, original$baseline_sd, '/')
  old_design <- cbind(d, d * b_old)
  old_query <- cbind(dt, dt * bt_old)
  ids <- sort(unique(train$depmap_id))
  set.seed(20261010L)
  fold_map <- setNames(sample(rep(seq_len(cv_folds), length.out = length(ids))), ids)
  fold <- unname(fold_map[train$depmap_id])
  write_table(data.table::data.table(depmap_id = ids, fold = unname(fold_map)), 'training_folds')
  commit(list(signature = signature, family = family, excluded_cells = settings$held_ids,
    excluded_drugs = settings$excluded_drugs, training_ids = ids, grid = svm_grid,
    cv_folds = cv_folds, selection = 'weighted absolute log2 error in held-out training cells',
    svm_row_weights = 'uniform: e1071 regression has no case weights',
    solver_tolerance = c(linear = .01, radial = .001),
    seed = 20261010L), file.path(root, 'nonlinear_settings.qs'))
  predictions <- list(original_ridge = predict_ridge(original$models$measured_D_conditioned, old_query))
  # Training-cell out-of-fold predictions calibrate the original frozen penalty.
  ridge_oof <- numeric(nrow(train))
  ridge_db_oof <- matrix(NA_real_, nrow(train), length(lambda_grid))
  for (k in seq_len(cv_folds)) {
    tr <- which(fold != k)
    va <- which(fold == k)
    f <- fit_path(old_design[tr, , drop = FALSE], train$observed[tr], observation_weights(train[tr]), FALSE)
    ridge_oof[va] <- pred_path(f, old_design[va, , drop = FALSE], original$models$measured_D_conditioned$lambda)
    g <- fit_path(features$DB[tr, , drop = FALSE], train$observed[tr], rep(1, length(tr)), TRUE)
    ridge_db_oof[va, ] <- pred_path(g, features$DB[va, , drop = FALSE])
  }
  w <- observation_weights(train)
  calibration <- lm(train$observed ~ ridge_oof, weights = w)
  predictions$calibrated_ridge <- unname(coef(calibration)[1L] + coef(calibration)[2L] * predictions$original_ridge)
  commit(list(signature = signature, calibration = calibration, oof = ridge_oof),
    file.path(model_dir, 'nonlinear_ridge_calibration.qs'))
  ridge_loss <- colSums(abs(ridge_db_oof - train$observed) * w) / sum(w)
  selected_lambda <- lambda_grid[which.min(ridge_loss)]
  ridge_db <- fit_path(features$DB, train$observed, rep(1, nrow(train)), TRUE)
  ridge_db$lambda <- selected_lambda
  predictions$ridge_DB <- predict_ridge(ridge_db, queries$DB)
  commit(list(signature = signature, model = ridge_db), file.path(model_dir, 'nonlinear_ridge_DB.qs'))
  cv_results <- list()
  calibrators <- list(original_ridge = coef(calibration))
  for (method in c('linear_SVR_DB', 'RBF_SVR_D', 'RBF_SVR_DB')) {
    design <- if (method == 'RBF_SVR_D') 'D' else 'DB'
    kernel <- if (method == 'linear_SVR_DB') 'linear' else 'radial'
    grid <- if (kernel == 'linear') data.table::data.table(cost = c(1, 10), gamma_multiplier = 0) else svm_grid
    xx <- features[[design]]
    oof <- matrix(NA_real_, nrow(train), nrow(grid))
    for (candidate in seq_len(nrow(grid))) {
      gamma <- if (kernel == 'linear') 1 / ncol(xx) else grid$gamma_multiplier[candidate] / ncol(xx)
      for (k in seq_len(cv_folds)) {
        status(paste(method, 'candidate', candidate, '/', nrow(grid), 'fold', k, '/', cv_folds))
        path <- file.path(model_dir, sprintf('nonlinear_%s_candidate%d_fold%d.qs', method, candidate, k))
        tr <- which(fold != k)
        va <- which(fold == k)
        if (file.exists(path)) {
          saved <- qs::qread(path)
          stopifnot(identical(saved$signature, signature), identical(saved$validation_rows, va))
          oof[va, candidate] <- saved$prediction
        } else {
          f <- fit_program_svr(xx[tr, , drop = FALSE], train$observed[tr], kernel, grid$cost[candidate], gamma)
          pp <- predict_program_svr(f, xx[va, , drop = FALSE])
          commit(list(signature = signature, model = f, validation_rows = va, prediction = pp), path)
          oof[va, candidate] <- pp
        }
      }
      cv_results[[paste(method, candidate)]] <- data.table::data.table(model = method, candidate = candidate,
        cost = grid$cost[candidate], gamma_multiplier = grid$gamma_multiplier[candidate],
        validation_absolute_log2_error = sum(abs(oof[, candidate] - train$observed) * w) / sum(w))
      write_table(data.table::rbindlist(cv_results), 'validation')
    }
    error <- colSums(abs(oof - train$observed) * w) / sum(w)
    best <- which.min(error)
    status(paste('fitting selected', method, '| candidate', best))
    gamma <- if (kernel == 'linear') 1 / ncol(xx) else grid$gamma_multiplier[best] / ncol(xx)
    final_path <- file.path(model_dir, paste0('nonlinear_', method, '.qs'))
    if (file.exists(final_path)) {
      saved <- qs::qread(final_path)
      stopifnot(identical(saved$signature, signature))
      f <- saved$model
    } else {
      f <- fit_program_svr(xx, train$observed, kernel, grid$cost[best], gamma)
      commit(list(signature = signature, model = f, selected = grid[best], training_ids = ids,
        training_drugs = unique(train$drug), oof = oof[, best]), final_path)
    }
    predictions[[method]] <- predict_program_svr(f, queries[[design]])
    if (kernel == 'radial') {
      oof_prediction <- oof[, best]
      cal <- lm(train$observed ~ oof_prediction, weights = w)
      predictions[[paste0('calibrated_', method)]] <- unname(coef(cal)[1L] + coef(cal)[2L] * predictions[[method]])
      calibrators[[method]] <- coef(cal)
      commit(list(signature = signature, calibration = cal),
        file.path(model_dir, paste0('nonlinear_', method, '_calibration.qs')))
    }
    partial <- data.table::rbindlist(lapply(names(predictions), function(name) {
      z <- test[, .(observation_id, depmap_id, drug, dose, observed)]
      z[, `:=`(model = name, predicted = predictions[[name]])]
      z
    }))
    data.table::fwrite(partial, file.path(data_dir, 'nonlinear_predictions.tsv.gz'), sep = '\t')
  }
  # Persist and hash predictions before reading evaluation outcomes for reporting.
  pred <- data.table::rbindlist(lapply(names(predictions), function(name) {
    z <- test[, .(observation_id, depmap_id, drug, dose, observed)]
    z[, `:=`(model = name, predicted = predictions[[name]])]
    z
  }))
  stopifnot(all(is.finite(pred$predicted)))
  pred[, `:=`(observed_pct = 100 * 2^observed, predicted_pct = 100 * 2^predicted,
    cell = panel$cell_map$depmap_name[match(depmap_id, panel$cell_map$depmap_id)],
    drug_set = ifelse(drug %in% settings$excluded_drugs, 'Excluded family', 'Seen drugs'))]
  data.table::fwrite(pred, file.path(data_dir, 'nonlinear_predictions.tsv.gz'), sep = '\t')
  write_table(data.table::data.table(path = file.path(data_dir, 'nonlinear_predictions.tsv.gz'),
    sha256 = digest::digest(file = file.path(data_dir, 'nonlinear_predictions.tsv.gz'), algo = 'sha256')), 'prediction_hash')
  per_cell <- pred[, .(conditions = .N, pearson_log2 = cor_safe(observed, predicted),
    pearson_pct = cor_safe(observed_pct, predicted_pct), sensitive = sum(observed_pct <= 50),
    sensitive_correct = sum(observed_pct <= 50 & predicted_pct <= 50), resistant = sum(observed_pct > 50),
    resistant_correct = sum(observed_pct > 50 & predicted_pct > 50),
    median_absolute_error_percentage_points = median(abs(predicted_pct - observed_pct))), by = .(model, drug_set, cell)]
  write_table(per_cell, 'by_cell')
  overall <- per_cell[, .(median_within_cell_r_log2 = median(pearson_log2),
    sensitive = sum(sensitive), sensitive_correct = sum(sensitive_correct), resistant = sum(resistant),
    resistant_correct = sum(resistant_correct),
    median_cell_absolute_error_percentage_points = median(median_absolute_error_percentage_points)), by = .(model, drug_set)]
  write_table(overall, 'summary')
  pred[, `:=`(yy = observed - mean(observed), pp = predicted - mean(predicted)), by = .(model, drug, dose)]
  write_table(pred[, .(pearson = cor_safe(yy, pp), observed_sd = sd(yy), predicted_sd = sd(pp)),
    by = .(model, drug_set)], 'within_treatment')
  write_table(data.table::rbindlist(lapply(names(calibrators), function(name)
    data.table::data.table(model = name, intercept = calibrators[[name]][1L], slope = calibrators[[name]][2L]))), 'calibration')
  print(overall)
  source('/work/src/program_pilot/family_transfer_plots.R')
  plot_nonlinear_response(root)
  status('COMPLETE: nonlinear comparison')
}
