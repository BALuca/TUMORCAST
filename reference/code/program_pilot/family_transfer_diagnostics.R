# Retrospective failure localization using frozen family banks. No discovery or
# held-cell labels enter fitting. Drug-specific probes explicitly admit family
# outcomes from training cells and are NOT family-transfer validation.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet)})
source('/work/src/program_pilot/response_models.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(2L)
lambda_grid <- c(100, 10, 1, .1, .01)
within_cell_only <- identical(commandArgs(TRUE), '--within-cell')
geometry_only <- identical(commandArgs(TRUE), '--geometry')
stopifnot(length(commandArgs(TRUE)) == 0L || within_cell_only || geometry_only)
panel <- qread('/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')

# Audit many-to-one and one-to-many response relationships without fitting models.
# Neighbors use molecular coordinates only, with no outcome-dependent selection.
audit_response_geometry <- function(labels, d, b, training_rows, root, bank) {
  coordinate_sd <- apply(d[training_rows, , drop = FALSE], 2L, sd)
  stopifnot(all(is.finite(coordinate_sd)), all(coordinate_sd > 0))
  z <- sweep(d, 2L, coordinate_sd, '/')
  groups <- list(
    same_treatment_different_cells = split(seq_len(nrow(labels)),
      paste(labels$new_cell, labels$drug, labels$dose)),
    same_cell_different_drugs = split(seq_len(nrow(labels)),
      paste(labels$depmap_id, labels$dose)))
  tables <- list()
  for (comparison in names(groups)) {
    parts <- lapply(groups[[comparison]], function(ii) {
      if (length(ii) < 2L) return(NULL)
      distances <- as.matrix(dist(z[ii, , drop = FALSE])) / sqrt(ncol(z))
      diag(distances) <- Inf
      nearest <- max.col(-distances, ties.method = 'first')
      endpoints <- which(upper.tri(distances), arr.ind = TRUE)
      left <- endpoints[, 1L]
      right <- endpoints[, 2L]
      data.table(i = ii[left], j = ii[right], displacement_distance = distances[endpoints],
        molecular_neighbor = nearest[left] == right | nearest[right] == left)
    })
    pairs <- rbindlist(parts)
    stopifnot(nrow(pairs) > 0L)
    a <- z[pairs$i, , drop = FALSE]
    c <- z[pairs$j, , drop = FALSE]
    na <- sqrt(rowSums(a^2))
    nc <- sqrt(rowSums(c^2))
    pairs[, `:=`(comparison = comparison,
      cohort = ifelse(labels$new_cell[i], 'Held cells', 'Development cells'),
      direction_cosine = rowSums(a * c) / (na * nc),
      magnitude_ratio = pmin(na, nc) / pmax(na, nc),
      raw_program_rms = sqrt(rowMeans((d[i, , drop = FALSE] - d[j, , drop = FALSE])^2)),
      baseline_distance = sqrt(rowMeans((b[i, , drop = FALSE] - b[j, , drop = FALSE])^2)),
      abundance_i = labels$prism_abundance_pct[i], abundance_j = labels$prism_abundance_pct[j])]
    pairs[, `:=`(outcome_gap_pp = abs(abundance_i - abundance_j),
      sensitive_resistant_discordant = (abundance_i <= 50) != (abundance_j <= 50),
      strong_discordant = pmin(abundance_i, abundance_j) <= 30 & pmax(abundance_i, abundance_j) >= 80,
      outcome_group = ifelse(pmax(abundance_i, abundance_j) <= 50, 'Both sensitive',
        ifelse(pmin(abundance_i, abundance_j) > 50, 'Both resistant', 'Mixed')))]
    pairs[, `:=`(matched_direction_amplitude = direction_cosine >= .9 & magnitude_ratio >= .8,
      function_matched = outcome_gap_pp <= 10,
      observation_i = labels$observation_id[i], observation_j = labels$observation_id[j],
      cell_i = labels$recipient[i], cell_j = labels$recipient[j],
      drug_i = labels$drug[i], drug_j = labels$drug[j], dose = labels$dose[i])]
    stopifnot(all(is.finite(pairs$displacement_distance)), all(is.finite(pairs$outcome_gap_pp)),
      !anyDuplicated(pairs[, .(i, j)]))
    tables[[comparison]] <- pairs
  }
  pairs <- rbindlist(tables)
  save <- function(x, name) fwrite(x, file.path(root, 'data', paste0('diagnostic_geometry_', name)), sep = '\t')
  overlap <- copy(labels[is.finite(almanac_abundance_pct)])
  overlap[, `:=`(prism_centered = prism_abundance_pct - mean(prism_abundance_pct),
    almanac_centered = almanac_abundance_pct - mean(almanac_abundance_pct)), by = .(drug, dose)]
  save(overlap, 'assay_overlap.tsv')
  save(overlap[, .(conditions = .N, prism_sensitive = sum(prism_abundance_pct <= 50),
    almanac_sensitive = sum(almanac_abundance_pct <= 50),
    same_sensitive_class = sum((prism_abundance_pct <= 50) == (almanac_abundance_pct <= 50)),
    pearson_pct = cor(prism_abundance_pct, almanac_abundance_pct),
    within_treatment_pearson_pct = cor(prism_centered, almanac_centered),
    median_assay_gap_pp = median(abs(prism_abundance_pct - almanac_abundance_pct)))], 'assay_summary.tsv')
  # Use same-assay vehicle expression for an estimated post-treatment state.
  # The response uses plate-matched references; the cached baseline is the mean
  # across vehicle plates, so baseline + response is a state proxy, not exact
  # reconstructed treated expression. Never add CCLE TPM to these coordinates.
  ids <- unique(labels$depmap_id)
  vehicle <- matrix(NA_real_, length(ids), ncol(d), dimnames = list(ids, colnames(d)))
  for (id in ids) {
    molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
    stopifnot(identical(molecular$normalizer, 'paired_trimmed_median_log_ratio_cp10k_v1'))
    gi <- match(bank$genes, molecular$genes)
    stopifnot(!anyNA(gi))
    vehicle[id, ] <- as.numeric(molecular$baseline[[id]][gi] %*% bank$weights)
  }
  baseline_z <- sweep(vehicle[labels$depmap_id, , drop = FALSE], 2L, coordinate_sd, '/')
  post_z <- baseline_z + z
  stopifnot(all(is.finite(post_z)))
  pairs[, `:=`(baseline_expression_distance = NA_real_, post_state_distance = NA_real_)]
  for (start in seq.int(1L, nrow(pairs), by = 20000L)) {
    ii <- seq.int(start, min(start + 19999L, nrow(pairs)))
    a <- pairs$i[ii]
    c <- pairs$j[ii]
    pairs[ii, `:=`(
      baseline_expression_distance = sqrt(rowMeans((baseline_z[a, ] - baseline_z[c, ])^2)),
      post_state_distance = sqrt(rowMeans((post_z[a, ] - post_z[c, ])^2)))]
  }
  # Adding a shared baseline cannot change within-cell Euclidean distances.
  within <- pairs[comparison == 'same_cell_different_drugs']
  stopifnot(max(abs(within$post_state_distance - within$displacement_distance)) < 1e-12,
    max(abs(within$baseline_expression_distance)) == 0)
  pairs[, unit := ifelse(comparison == 'same_cell_different_drugs', cell_i, drug_i)]
  pairs[, stratum := paste(comparison, cohort, unit, dose, sep = '__')]
  save(pairs, 'pairs.tsv.gz')
  # Compare molecular-distance deciles within each fixed biological comparison.
  # Summarize each stratum first so large cell/drug groups do not dominate.
  long <- melt(pairs[, .(comparison, cohort, dose, unit, stratum, i, j,
    displacement_distance, post_state_distance, outcome_gap_pp, sensitive_resistant_discordant)],
    id.vars = c('comparison', 'cohort', 'dose', 'unit', 'stratum', 'i', 'j',
      'outcome_gap_pp', 'sensitive_resistant_discordant'),
    measure.vars = c('displacement_distance', 'post_state_distance'),
    variable.name = 'representation', value.name = 'distance')
  long[, decile := pmin(10L, ceiling(10 * frank(distance, ties.method = 'first') / .N)),
    by = .(stratum, representation)]
  strata <- long[, .(pairs = .N, median_gap_pp = median(outcome_gap_pp),
    discordant_fraction = mean(sensitive_resistant_discordant), median_distance = median(distance)),
    by = .(comparison, cohort, dose, unit, stratum, representation, decile)]
  save(strata, 'distance_strata.tsv')
  trend <- strata[, .(strata = .N, pairs = sum(pairs), median_gap_pp = median(median_gap_pp),
    q25_gap_pp = quantile(median_gap_pp, .25), q75_gap_pp = quantile(median_gap_pp, .75),
    median_discordant_fraction = median(discordant_fraction)),
    by = .(comparison, cohort, representation, decile)]
  save(trend, 'distance_trend.tsv')
  dose_trend <- strata[, .(strata = .N, median_gap_pp = median(median_gap_pp)),
    by = .(comparison, cohort, dose, representation, decile)]
  save(dose_trend, 'distance_by_dose.tsv')
  # A molecular nearest neighbour for each observation; no outcome is used.
  edges <- rbindlist(list(long[, .(comparison, cohort, unit, stratum, representation,
    observation = i, neighbour = j, distance, outcome_gap_pp, sensitive_resistant_discordant)],
    long[, .(comparison, cohort, unit, stratum, representation,
      observation = j, neighbour = i, distance, outcome_gap_pp, sensitive_resistant_discordant)]))
  setorder(edges, stratum, representation, observation, distance, neighbour)
  nearest <- edges[, .SD[1L], by = .(stratum, representation, observation)]
  save(nearest, 'nearest_observations.tsv.gz')
  neighbor_summary <- nearest[, .(observations = .N, median_gap_pp = median(outcome_gap_pp),
    discordant = sum(sensitive_resistant_discordant), discordant_fraction = mean(sensitive_resistant_discordant)),
    by = .(comparison, cohort, representation)]
  save(neighbor_summary, 'nearest_summary.tsv')
  # Show the matched-direction/amplitude rule in every cell and drug, not just
  # selected counterexamples. Also retain all-pair denominators for context.
  units <- rbindlist(lapply(c(FALSE, TRUE), function(matched) {
    pairs[!matched | matched_direction_amplitude,
      .(pairs = .N, discordant = sum(sensitive_resistant_discordant),
        discordant_fraction = mean(sensitive_resistant_discordant),
        strong_discordant = sum(strong_discordant), median_gap_pp = median(outcome_gap_pp)),
      by = .(comparison, cohort, unit)][, selection := if (matched) 'Matched profiles' else 'All pairs']
  }))
  save(units, 'by_cell_drug.tsv')
  print(neighbor_summary)
  summary <- rbindlist(lapply(c('all', 'nearest', 'direction_amplitude'), function(selection) {
    use <- switch(selection, all = rep(TRUE, nrow(pairs)), nearest = pairs$molecular_neighbor,
      direction_amplitude = pairs$matched_direction_amplitude)
    pairs[which(use), .(pairs = .N, median_outcome_gap_pp = median(outcome_gap_pp),
      discordant_pairs = sum(sensitive_resistant_discordant), strong_discordant_pairs = sum(strong_discordant),
      median_displacement_distance = median(displacement_distance),
      median_direction_cosine = median(direction_cosine), median_magnitude_ratio = median(magnitude_ratio)),
      by = .(comparison, cohort)][, selection := selection]
  }))
  save(summary, 'summary.tsv')
  function_matched <- pairs[function_matched == TRUE,
    .(pairs = .N, median_displacement_distance = median(displacement_distance),
      median_direction_cosine = median(direction_cosine), opposite_direction_pairs = sum(direction_cosine < 0),
      strong_pattern_matches = sum(matched_direction_amplitude)), by = .(comparison, cohort, outcome_group)]
  save(function_matched, 'function_matched.tsv')
  examples <- pairs[molecular_neighbor & strong_discordant & pmax(abundance_i, abundance_j) <= 150]
  setorder(examples, displacement_distance)
  examples <- examples[, head(.SD, 5L), by = .(comparison, cohort)]
  save(examples, 'discordant_examples.tsv')
  convergent <- pairs[function_matched == TRUE & outcome_group == 'Both sensitive' & direction_cosine < 0]
  setorder(convergent, outcome_gap_pp, -displacement_distance)
  convergent <- convergent[, head(.SD, 5L), by = .(comparison, cohort)]
  save(convergent, 'convergent_examples.tsv')
  # For named examples, test resolution against empirical vehicle contrasts.
  # These controls provide a background reference, not treated-replicate error.
  selected <- rbindlist(list(examples[, example_type := 'Similar programs, different outcomes'],
    convergent[, example_type := 'Different programs, similar outcomes']))
  endpoints <- unique(c(selected$i, selected$j))
  profiles <- list()
  for (id in unique(labels$depmap_id[endpoints])) {
    path <- panel$input_paths[match(id, panel$cell_map$depmap_id)]
    molecular <- qread(path)
    obs <- labels$observation_id[endpoints[labels$depmap_id[endpoints] == id]]
    stopifnot(all(obs %in% rownames(molecular$response)), all(bank$genes %in% colnames(molecular$response)))
    gene_response <- molecular$response[obs, bank$genes, drop = FALSE]
    expected <- d[match(obs, labels$observation_id), , drop = FALSE]
    stopifnot(max(abs(as.matrix(gene_response %*% bank$weights) - expected)) < 1e-12)
    null_genes <- molecular$null[[id]][, bank$genes, drop = FALSE]
    profiles[[id]] <- list(response = gene_response, null_genes = null_genes,
      null_programs = as.matrix(null_genes %*% bank$weights))
  }
  resolution <- rbindlist(lapply(seq_len(nrow(selected)), function(k) {
    pair <- selected[k]
    id_i <- labels$depmap_id[pair$i]
    id_j <- labels$depmap_id[pair$j]
    p_i <- profiles[[id_i]]
    p_j <- profiles[[id_j]]
    g_i <- p_i$response[pair$observation_i, ]
    g_j <- p_j$response[pair$observation_j, ]
    null_endpoints <- CJ(a = seq_len(nrow(p_i$null_programs)), c = seq_len(nrow(p_j$null_programs)))
    if (id_i == id_j) null_endpoints <- null_endpoints[a != c]
    stopifnot(nrow(null_endpoints) > 0L)
    program_background <- sqrt(rowMeans((p_i$null_programs[null_endpoints$a, , drop = FALSE] -
      p_j$null_programs[null_endpoints$c, , drop = FALSE])^2))
    gene_background <- sqrt(rowMeans((p_i$null_genes[null_endpoints$a, , drop = FALSE] -
      p_j$null_genes[null_endpoints$c, , drop = FALSE])^2))
    gene_difference <- sqrt(mean((g_i - g_j)^2))
    data.table(observation_i = pair$observation_i, observation_j = pair$observation_j,
      example_type = pair$example_type, raw_gene_pearson = cor(g_i, g_j),
      raw_gene_rms = gene_difference, gene_vehicle_difference_95 = quantile(gene_background, .95),
      program_vehicle_difference_95 = quantile(program_background, .95),
      program_difference_over_vehicle_95 = pair$raw_program_rms / quantile(program_background, .95),
      gene_difference_over_vehicle_95 = gene_difference / quantile(gene_background, .95))
  }))
  save(merge(selected, resolution, by = c('observation_i', 'observation_j', 'example_type')),
    'example_resolution.tsv')
  print(summary)
  print(function_matched)
  print(examples[, .(comparison, cohort, cell_i, cell_j, drug_i, drug_j, dose,
    abundance_i, abundance_j, displacement_distance, direction_cosine, magnitude_ratio)])
}

# Remove treatment means inside each fitting fold, including validation targets.
# Query predictions represent relative sensitivity, not absolute abundance.
fit_within_treatment <- function(x, y, metadata) {
  fold <- unname(inner_map[metadata$depmap_id])
  loss <- rep(0, length(lambda_grid))
  for (k in sort(unique(fold))) {
    tr <- which(fold != k)
    va <- which(fold == k)
    conditions <- unique(metadata$condition[tr])
    mx <- t(vapply(conditions, function(a) colMeans(x[tr[metadata$condition[tr] == a], , drop = FALSE]), numeric(ncol(x))))
    my <- vapply(conditions, function(a) mean(y[tr[metadata$condition[tr] == a]]), numeric(1))
    stopifnot(all(metadata$condition[va] %in% conditions))
    tx <- x[tr, , drop = FALSE] - mx[match(metadata$condition[tr], conditions), , drop = FALSE]
    vx <- x[va, , drop = FALSE] - mx[match(metadata$condition[va], conditions), , drop = FALSE]
    ty <- y[tr] - my[match(metadata$condition[tr], conditions)]
    vy <- y[va] - my[match(metadata$condition[va], conditions)]
    f <- fit_path(tx, ty, observation_weights(metadata[tr]), FALSE)
    loss <- loss + colSums((pred_path(f, vx) - vy)^2 * observation_weights(metadata[va]))
  }
  centered_x <- x
  centered_y <- y
  for (condition in unique(metadata$condition)) {
    ii <- which(metadata$condition == condition)
    centered_x[ii, ] <- sweep(x[ii, , drop = FALSE], 2, colMeans(x[ii, , drop = FALSE]))
    centered_y[ii] <- y[ii] - mean(y[ii])
  }
  f <- fit_path(centered_x, centered_y, observation_weights(metadata), FALSE)
  f$lambda <- lambda_grid[which.min(loss)]
  f
}

for (family in c('MEK', 'ERBB')) {
  message('Diagnosing ', family)
  root <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
  data_dir <- file.path(root, 'data')
  save_table <- function(x, name) fwrite(x, file.path(data_dir, paste0('diagnostic_', name, '.tsv')), sep = '\t')
  settings <- qread(file.path(root, 'settings.qs'))
  bank <- qread(file.path(root, 'models/program_bank.qs'))
  prepared <- qread(file.path(data_dir, 'prepared_scores.qs'))
  decoder <- qread(file.path(root, 'models/functional_decoder.qs'))
  programs <- colnames(prepared$measured)
  labels <- as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
    depmap_id %in% rownames(prepared$baseline)]
  labels[, `:=`(observed = log2(prism_abundance_pct / 100), condition = paste(drug, dose),
    new_cell = depmap_id %in% settings$held_ids, new_family = drug %in% settings$excluded_drugs)]
  d <- prepared$measured[match(labels$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  b <- sweep(sweep(prepared$baseline[labels$depmap_id, , drop = FALSE], 2, decoder$baseline_mean), 2, decoder$baseline_sd, '/')
  x <- cbind(d, d * b)
  tr <- which(!labels$new_cell & !labels$new_family)
  stopifnot(setequal(unique(labels$depmap_id[tr]), decoder$training_ids), all(is.finite(x)))
  if (geometry_only) {
    audit_response_geometry(labels, d, b, tr, root, bank)
    source('/work/src/program_pilot/family_transfer_plots.R')
    plot_systematic_response_geometry(root)
    # Show audited pairs in the MEK bank without refitting the representation.
    if (family == 'MEK') {
      plot_response_geometry(root)
      plot_response_distance_scatter(root)
    }
    next
  }
  ids <- sort(unique(labels$depmap_id[tr]))
  set.seed(20261010L)
  inner_map <- setNames(sample(rep(1:5, length.out = length(ids))), ids)
  labels[, predicted := predict_ridge(decoder$models$measured_D_conditioned, x)]
  if (within_cell_only) {
    # Compare shared and cell-specific mappings on identical held-out drugs.
    # Discovery is frozen and used these outcomes: conditional diagnostic only.
    development <- copy(labels[tr])
    development[, `:=`(shared_D = NA_real_, shared_conditioned = NA_real_, cell_specific_D = NA_real_)]
    drugs <- sort(unique(development$drug))
    drug_folds <- setNames(1L + (seq_along(drugs) - 1L) %% 5L, drugs)
    original_inner_map <- inner_map
    for (outer in 1:5) {
      message(family, ' | whole-drug diagnostic fold ', outer, '/5')
      training_rows <- which(drug_folds[development$drug] != outer)
      evaluation_rows <- which(drug_folds[development$drug] == outer)
      inner_map <- original_inner_map
      for (design in c('shared_D', 'shared_conditioned')) {
        features <- if (design == 'shared_D') d[tr, , drop = FALSE] else x[tr, , drop = FALSE]
        f <- fit_ridge(features[training_rows, , drop = FALSE], development$observed[training_rows],
          development$depmap_id[training_rows], observation_weights(development[training_rows]), FALSE)
        development[evaluation_rows, (design) := predict_ridge(f, features[evaluation_rows, , drop = FALSE])]
      }
      for (id in ids) {
        train_cell <- training_rows[development$depmap_id[training_rows] == id]
        eval_cell <- evaluation_rows[development$depmap_id[evaluation_rows] == id]
        if (!length(eval_cell)) next
        local_drugs <- sort(unique(development$drug[train_cell]))
        stopifnot(length(local_drugs) >= 5L,
          !any(development$drug[eval_cell] %in% development$drug[train_cell]))
        inner_map <- setNames(1L + (seq_along(local_drugs) - 1L) %% 5L, local_drugs)
        f <- fit_ridge(d[tr[train_cell], , drop = FALSE], development$observed[train_cell],
          development$drug[train_cell], observation_weights(development[train_cell]), FALSE)
        development[eval_cell, cell_specific_D := predict_ridge(f, d[tr[eval_cell], , drop = FALSE])]
      }
    }
    inner_map <- original_inner_map
    stopifnot(all(is.finite(as.matrix(development[, .(shared_D, shared_conditioned, cell_specific_D)]))))
    development[, cell := panel$cell_map$depmap_name[match(depmap_id, panel$cell_map$depmap_id)]]
    fwrite(development, file.path(data_dir, 'diagnostic_within_cell_predictions.tsv.gz'), sep = '\t')
    long_cell <- melt(development, id.vars = c('observation_id', 'depmap_id', 'cell', 'drug', 'dose', 'observed'),
      measure.vars = c('shared_D', 'shared_conditioned', 'cell_specific_D'), variable.name = 'model', value.name = 'p')
    cell_summary <- long_cell[, .(observations = .N, pearson = cor_safe(observed, p),
      median_absolute_error_percentage_points = median(abs(100 * 2^observed - 100 * 2^p))), by = .(cell, model)]
    save_table(cell_summary, 'within_cell_model_comparison')
    long_cell[, `:=`(yc = observed - mean(observed), pc = p - mean(p)), by = .(model, drug, dose)]
    save_table(long_cell[, .(within_drug_dose_r = cor_safe(yc, pc), observed_sd = sd(yc), predicted_sd = sd(pc)),
      by = model], 'within_cell_models_across_cells')
    print(cell_summary[, .(median_within_cell_r = median(pearson)), by = model])
    source('/work/src/program_pilot/family_transfer_plots.R')
    plot_within_cell_models(root)
    next
  }
  labels[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)),
    by = .(new_cell, new_family, condition)]
  quadrants <- labels[, .(rows = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug),
    overall_r = cor_safe(observed, predicted), within_treatment_r = cor_safe(yc, pc),
    observed_within_sd = sd(yc), predicted_within_sd = sd(pc)), by = .(new_cell, new_family)]
  save_table(quadrants, 'transfer_quadrants')
  print(quadrants)
  save_table(labels[, .(pearson = cor_safe(observed, predicted)), by = .(new_cell, new_family, depmap_id)], 'within_cell')
  training <- labels[tr]
  outcome_between <- 1 - sum(training$yc^2) / sum((training$observed - mean(training$observed))^2)
  feature_variance <- rbindlist(lapply(seq_along(programs), function(j) {
    v <- d[tr, j]
    residual <- v - ave(v, training$condition, FUN = mean)
    data.table(program = programs[j], treatment_mean_fraction = 1 - sum(residual^2) / sum((v - mean(v))^2))
  }))
  save_table(feature_variance, 'program_variance')
  save_table(data.table(outcome_treatment_mean_fraction = outcome_between,
    median_program_treatment_mean_fraction = median(feature_variance$treatment_mean_fraction),
    existing_lambda = decoder$models$measured_D_conditioned$lambda), 'training_objective')

  centered_model <- fit_within_treatment(x[tr, , drop = FALSE], training$observed, training)
  qsave(centered_model, file.path(root, 'models/diagnostic_within_treatment_decoder.qs'))
  evaluation <- copy(labels[new_cell & new_family])
  te <- which(labels$new_cell & labels$new_family)
  evaluation[, `:=`(existing = predicted, centered_training = predict_ridge(centered_model, x[te, , drop = FALSE]))]
  evaluation[, drug_specific := NA_real_]
  evaluation[, `:=`(drug_specific_displacement = NA_real_, drug_specific_baseline = NA_real_)]
  associations <- list()
  for (drug_name in unique(evaluation$drug)) {
    local_tr <- which(!labels$new_cell & labels$drug == drug_name)
    local_te <- which(labels$new_cell & labels$drug == drug_name)
    # Additional functional family data are allowed ONLY for this diagnostic.
    f <- fit_ridge(x[local_tr, , drop = FALSE], labels$observed[local_tr], labels$depmap_id[local_tr],
      observation_weights(labels[local_tr]), TRUE)
    evaluation[drug == drug_name, drug_specific := predict_ridge(f, x[local_te, , drop = FALSE])]
    for (feature_set in c('displacement', 'baseline')) {
      features <- if (feature_set == 'displacement') d else b
      ff <- fit_ridge(features[local_tr, , drop = FALSE], labels$observed[local_tr], labels$depmap_id[local_tr],
        observation_weights(labels[local_tr]), TRUE)
      evaluation[drug == drug_name, (paste0('drug_specific_', feature_set)) := predict_ridge(ff, features[local_te, , drop = FALSE])]
    }
    rtrain <- rtest <- numeric(length(programs))
    for (j in seq_along(programs)) {
      train_y <- labels$observed[local_tr] - ave(labels$observed[local_tr], labels$dose[local_tr], FUN = mean)
      test_y <- labels$observed[local_te] - ave(labels$observed[local_te], labels$dose[local_te], FUN = mean)
      train_d <- d[local_tr, j] - ave(d[local_tr, j], labels$dose[local_tr], FUN = mean)
      test_d <- d[local_te, j] - ave(d[local_te, j], labels$dose[local_te], FUN = mean)
      rtrain[j] <- cor_safe(train_y, train_d)
      rtest[j] <- cor_safe(test_y, test_d)
    }
    associations[[drug_name]] <- data.table(drug = drug_name, program = programs, training_r = rtrain,
      held_cell_r = rtest, strongest_training_program = seq_along(programs) == which.max(abs(rtrain)))
  }
  qsave(evaluation, file.path(data_dir, 'diagnostic_predictions.qs'))
  save_table(rbindlist(associations), 'individual_program_associations')
  long <- melt(evaluation, id.vars = c('drug', 'dose', 'depmap_id', 'observed'),
    measure.vars = c('existing', 'centered_training', 'drug_specific', 'drug_specific_displacement', 'drug_specific_baseline'),
    variable.name = 'diagnostic', value.name = 'p')
  long[, `:=`(yc = observed - mean(observed), pc = p - mean(p)), by = .(drug, dose, diagnostic)]
  comparison <- long[, .(within_dose_r = cor_safe(yc, pc), observed_sd = sd(yc), predicted_sd = sd(pc)), by = .(drug, diagnostic)]
  save_table(comparison, 'decoder_comparison')
  save_table(long[diagnostic != 'centered_training', .(
    observed_min_pct = min(100 * 2^observed), observed_max_pct = max(100 * 2^observed),
    predicted_min_pct = min(100 * 2^p), predicted_max_pct = max(100 * 2^p)),
    by = .(drug, dose, diagnostic)], 'assay_ranges')
  print(dcast(comparison, drug ~ diagnostic, value.var = 'within_dose_r'))

  # Trace the change from per-cell admitted genes to a fixed full-gene mean.
  score_audit <- list()
  for (id in bank$training_ids) {
    i <- match(id, panel$cell_map$depmap_id)
    molecular <- qread(panel$input_paths[i])
    local <- qread(file.path(root, 'models', paste0('source_', id, '.qs')))
    md <- as.data.table(molecular$metadata)[!drug %in% settings$excluded_drugs]
    xx <- molecular$response[match(md$observation_id, rownames(molecular$response)), bank$genes, drop = FALSE]
    yy <- labels$observed[match(md$observation_id, labels$observation_id)]
    for (j in seq_along(programs)) {
      members <- bank$genes[bank$partition == j]
      admitted <- intersect(members, local$genes)
      full <- rowMeans(xx[, members, drop = FALSE])
      subgroup <- if(length(admitted) >= 3L) rowMeans(xx[, admitted, drop = FALSE]) else rep(NA_real_, nrow(xx))
      ok <- is.finite(yy)
      centered_genes <- scale(xx[, members, drop = FALSE], center = TRUE, scale = FALSE)
      denom <- mean(colMeans(centered_genes^2))
      score_audit[[paste(id, j)]] <- data.table(depmap_id = id, program = programs[j], genes = length(members),
        admitted_genes = length(admitted), coherent_variance_fraction = mean((full - mean(full))^2) / denom,
        full_outcome_r = if(sum(ok) >= 3L) cor_safe(full[ok], yy[ok]) else NA_real_,
        admitted_outcome_r = if(sum(ok) >= 3L && length(admitted) >= 3L) cor_safe(subgroup[ok], yy[ok]) else NA_real_,
        full_admitted_r = if(length(admitted) >= 3L) cor_safe(full, subgroup) else NA_real_)
    }
  }
  save_table(rbindlist(score_audit), 'hard_mean_audit')
  save_table(data.table(path = c('/work/src/program_pilot/family_transfer_diagnostics.R',
    file.path(data_dir, 'prepared_scores.qs'), file.path(root, 'models/functional_decoder.qs')),
    sha256 = vapply(c('/work/src/program_pilot/family_transfer_diagnostics.R',
      file.path(data_dir, 'prepared_scores.qs'), file.path(root, 'models/functional_decoder.qs')),
      digest::digest, character(1), file = TRUE, algo = 'sha256')), 'provenance')
  message('Completed diagnostics: ', family)
}
