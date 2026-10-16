#!/usr/bin/env Rscript
# Test whether a small molecular interaction model transfers to HGF-like partners.
# Drug nomination and program selection use molecular data only. The exact B1
# monotherapy predictions are protected; no HGF functional outcome is fitted.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(digest)
  library(ggplot2)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

# Fixed, small experiment: six same-direction partners, three opposite controls,
# five HGF-supported programs and their 15 symmetric pair products.
workspace <- '/work'
root <- file.path(workspace, 'results/program_pilot/75_hgf_pattern_partner_transfer')
reference <- file.path(workspace, 'results/program_pilot/74_v20_lapatinib_hgf_frozen_query')
original <- file.path(reference, 'inputs/original')
anchor <- 'Lapatinib ditosylate'
programs <- sprintf('P%02d', 1:40)
nomination_dose <- 2.5
positive_partners <- 6L
negative_partners <- 3L
interaction_programs <- 5L
minimum_sign_agreement <- .75
minimum_partner_rows <- 60L
minimum_partner_cells <- 25L
lambda_grid <- 10^seq(1, -3, length.out = 20L)
elastic_alpha <- 1
query_doses <- c(1, 2.5)
seed <- 75012L
set.seed(seed)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
status('Reading the frozen original40 inputs and exact B1 models')
setwd(original)
source('modular_pair_pilot/model.R')
stopifnot(identical(colnames(PREP$W), programs), nrow(PREP$W) == 8842L,
  !AU %in% PAIR$depmap_id, !AU %in% SINGLE$depmap_id)
models <- lapply(1:3, function(k) readRDS(file.path(reference, 'inputs/results', paste0('fold_', k, '.rds')))$models$full)
stopifnot(all(vapply(models, function(x) !AU %in% x$training_single_ids, logical(1))))
molecular <- as.data.table(PREP$m)
all_pairs <- as.data.table(PAIR)
stopifnot(nrow(molecular) == 4410L, uniqueN(molecular$depmap_id) == 42L,
  uniqueN(molecular$drug) == 35L, nrow(all_pairs) == 50948L,
  max(abs(all_pairs$I - (all_pairs$y - all_pairs$yA - all_pairs$yB))) < 1e-10)

# This fixed molecular metric measures variation across receiving cells within
# drug/dose. It is not fitted using functional outcomes or Wilson measurements.
centered <- molecular[, lapply(.SD, function(x) x - mean(x)), by = .(drug, dose), .SDcols = programs]
response_scale <- vapply(centered[, programs, with = FALSE], sd, numeric(1))
stopifnot(all(is.finite(response_scale)), all(response_scale > 0))
hgf_path <- file.path(reference, 'data/hgf_transported_inputs.tsv')
hgf <- fread(hgf_path)
stopifnot(nrow(hgf) == 13L, sum(hgf$context == 'Equal-context mean') == 1L)
source_hgf <- as.matrix(hgf[context != 'Equal-context mean', programs, with = FALSE])
mean_hgf <- as.numeric(hgf[context == 'Equal-context mean', programs, with = FALSE])
scaled_hgf <- sweep(source_hgf, 2L, response_scale, '/')
scaled_mean_hgf <- mean_hgf / response_scale

# Nominate using individually measured single responses interpolated to 2.5uM.
# Characterized single interventions are legal inputs for excluded-pair transfer.
lo <- molecular[dose == .5, c('depmap_id', 'drug', programs), with = FALSE]
hi <- molecular[dose == 5, c('depmap_id', 'drug', programs), with = FALSE]
nomination <- merge(lo, hi, by = c('depmap_id', 'drug'), suffixes = c('_lo', '_hi'))
stopifnot(nrow(nomination) == 42L * 35L, !anyDuplicated(nomination[, .(depmap_id, drug)]))
fraction <- log(nomination_dose / .5) / log(5 / .5)
native <- (1 - fraction) * as.matrix(nomination[, paste0(programs, '_lo'), with = FALSE]) +
  fraction * as.matrix(nomination[, paste0(programs, '_hi'), with = FALSE])
scaled <- sweep(native, 2L, response_scale, '/')
cosine <- as.numeric(scaled %*% scaled_mean_hgf) /
  sqrt(rowSums(scaled^2) * sum(scaled_mean_hgf^2))
stopifnot(all(is.finite(cosine)))
nomination[, hgf_cosine := cosine]
ranking <- nomination[, .(median_cell_cosine = median(hgf_cosine),
  positive_cells = sum(hgf_cosine > 0), negative_cells = sum(hgf_cosine < 0), molecular_cells = .N), by = drug]
means <- rowsum(scaled, nomination$drug)
means <- means / as.numeric(table(nomination$drug)[rownames(means)])
context_cosine <- means %*% t(scaled_hgf) /
  outer(sqrt(rowSums(means^2)), sqrt(rowSums(scaled_hgf^2)))
ranking[, median_source_context_cosine := apply(context_cosine, 1L, median)[match(drug, rownames(context_cosine))]]
lap_pairs <- copy(all_pairs[drugA == anchor | drugB == anchor])
lap_pairs[, partner := fifelse(drugA == anchor, drugB, drugA)]
coverage <- lap_pairs[, .(combination_rows = .N, combination_cells = uniqueN(depmap_id)), by = .(drug = partner)]
ranking <- merge(ranking, coverage, by = 'drug', all.x = TRUE)
ranking[, eligible_coverage := !is.na(combination_rows) & combination_rows >= minimum_partner_rows &
  combination_cells >= minimum_partner_cells & drug != anchor]
same <- ranking[eligible_coverage & median_cell_cosine > .1 &
  positive_cells / molecular_cells >= minimum_sign_agreement & median_source_context_cosine > 0]
setorder(same, -median_source_context_cosine, drug)
opposite <- ranking[eligible_coverage & median_cell_cosine < -.1 &
  negative_cells / molecular_cells >= minimum_sign_agreement & median_source_context_cosine < 0]
setorder(opposite, median_source_context_cosine, drug)
stopifnot(nrow(same) >= positive_partners, nrow(opposite) >= negative_partners)
chosen <- rbind(same[seq_len(positive_partners)][, role := 'HGF-pattern neighbor'],
  opposite[seq_len(negative_partners)][, role := 'Opposite-pattern control'])
partners <- chosen$drug
stopifnot(length(partners) == 9L, !anyDuplicated(partners))
write_table(ranking[order(-median_source_context_cosine)], 'all_drug_molecular_rankings')
write_table(chosen, 'selected_partners')
write_table(cbind(nomination[, .(depmap_id, drug, hgf_cosine)], as.data.table(native)), 'nomination_single_responses')

# HGF-supported dimensions are selected before any combination outcome fitting.
agreement <- pmax(colMeans(scaled_hgf > 0), colMeans(scaled_hgf < 0))
program_support <- data.table(program = programs, sign_agreement = agreement,
  median_absolute_scaled_hgf = apply(abs(scaled_hgf), 2L, median),
  median_scaled_hgf = apply(scaled_hgf, 2L, median))
eligible <- program_support[sign_agreement >= minimum_sign_agreement]
setorder(eligible, -median_absolute_scaled_hgf, program)
stopifnot(nrow(eligible) >= interaction_programs)
selected_programs <- eligible$program[seq_len(interaction_programs)]
program_support[, selected := program %in% selected_programs]
write_table(program_support, 'hgf_program_support')
weights_path <- file.path(workspace,
  'almanac_conditionals/state_conditioned_bridge/program_first/results/frozen_program_weights.tsv')
membership <- fread(weights_path)
stopifnot(setequal(membership$gene, rownames(PREP$W)))
write_table(membership[program %in% selected_programs], 'selected_program_gene_membership')

# Use one outcome-free, frozen ensemble of the three molecular operators. Functional
# monotherapy heads are not involved in the interaction response coordinates.
d <- copy(all_pairs[drugA %in% c(anchor, partners) & drugB %in% c(anchor, partners)])
d[, row_id := .I]
stopifnot(nrow(d) > 1000L, uniqueN(d$depmap_id) >= 25L, !AU %in% d$depmap_id)
query <- unique(rbind(d[, .(depmap_id, drug = drugA, dose = doseA)],
  d[, .(depmap_id, drug = drugB, dose = doseB)]))
setorder(query, depmap_id, drug, dose)
query[, query_id := .I]
response <- matrix(0, nrow(query), 40L, dimnames = list(NULL, programs))
for (model in models) response <- response + base_predict(PREP, model$operators, query) / length(models)
query_key <- function(cell, drug, dose) paste(cell, drug, sprintf('%.10g', dose), sep = '|')
keys <- query_key(query$depmap_id, query$drug, query$dose)
stopifnot(!anyDuplicated(keys))
ia <- match(query_key(d$depmap_id, d$drugA, d$doseA), keys)
ib <- match(query_key(d$depmap_id, d$drugB, d$doseB), keys)
stopifnot(!anyNA(c(ia, ib)), all(is.finite(response)))
A <- sweep(response[ia, , drop = FALSE], 2L, response_scale, '/')
B <- sweep(response[ib, , drop = FALSE], 2L, response_scale, '/')
baseline <- getmatrix(PREP$raw, d$depmap_id, programs)
d[, partner := fifelse(drugA == anchor, drugB, fifelse(drugB == anchor, drugA, NA_character_))]

# Conservatively purge reused single-control settings. The prepared archive lacks
# study/date/sample IDs, so the screen/cell/drug/dose key over-purges rather than
# treating untraceable controls as independent observations.
reference_a <- paste(d$SCREENER, d$depmap_id, d$drugA, sprintf('%.10g', d$doseA), sep = '|')
reference_b <- paste(d$SCREENER, d$depmap_id, d$drugB, sprintf('%.10g', d$doseB), sep = '|')
training_rows <- function(excluded, test_rows, available = seq_len(nrow(d))) {
  tr <- available[!d$drugA[available] %in% excluded & !d$drugB[available] %in% excluded]
  test_refs <- unique(c(reference_a[test_rows], reference_b[test_rows]))
  tr <- tr[!reference_a[tr] %in% test_refs & !reference_b[tr] %in% test_refs]
  stopifnot(!length(intersect(tr, test_rows)),
    !any(excluded %in% c(d$drugA[tr], d$drugB[tr])),
    !length(intersect(c(reference_a[tr], reference_b[tr]), test_refs)))
  tr
}
pair_weights <- function(rows) {
  group <- paste(d$pair[rows], d$depmap_id[rows], sep = '|')
  cell_counts <- d[rows, .(cells = uniqueN(depmap_id)), by = pair]
  w <- 1 / as.numeric(table(group)[group]) / cell_counts$cells[match(d$pair[rows], cell_counts$pair)]
  w / sum(w)
}
ij <- which(upper.tri(matrix(0, interaction_programs, interaction_programs), diag = TRUE), arr.ind = TRUE)
program_index <- match(selected_programs, programs)
interaction_basis <- function(a, b) {
  stopifnot(nrow(a) == nrow(b), ncol(a) == 40L, ncol(b) == 40L)
  x <- vapply(seq_len(nrow(ij)), function(k) {
    i <- program_index[ij[k, 1]]
    j <- program_index[ij[k, 2]]
    (a[, i] * b[, j] + a[, j] * b[, i]) / 2
  }, numeric(nrow(a)))
  colnames(x) <- paste(selected_programs[ij[, 1]], selected_programs[ij[, 2]], sep = '__')
  x
}
pattern_basis <- interaction_basis(A, B)
strength_basis <- matrix(sqrt(rowSums(A^2)) * sqrt(rowSums(B^2)), ncol = 1L,
  dimnames = list(NULL, 'response_strength_product'))
stopifnot(max(abs(pattern_basis - interaction_basis(B, A))) < 1e-12,
  max(abs(interaction_basis(A * 0, B))) == 0)

# Baseline modulation is estimated from training receiver baselines only. Two PCs
# express a small context model; they do not replace the named response dimensions.
context_fit <- function(rows) {
  first <- rows[!duplicated(d$depmap_id[rows])]
  stopifnot(length(first) >= 10L)
  pc <- prcomp(baseline[first, , drop = FALSE], center = TRUE, scale. = TRUE, rank. = 2L)
  stopifnot(all(is.finite(pc$sdev[1:2])), all(pc$sdev[1:2] > 0))
  list(pc = pc, sd = pc$sdev[1:2], ids = d$depmap_id[first])
}
context_predict <- function(fit, b) {
  z <- predict(fit$pc, newdata = b)[, 1:2, drop = FALSE]
  tanh(sweep(z, 2L, fit$sd, '/'))
}
design_matrix <- function(base, gate) {
  x <- cbind(base, base * gate[, 1], base * gate[, 2])
  colnames(x) <- c(colnames(base), paste0(colnames(base), '__baseline_PC1'),
    paste0(colnames(base), '__baseline_PC2'))
  x
}

# These helpers are reused across every excluded partner and inner split.
fit_path <- function(train, test, family) {
  stopifnot(length(train) >= 100L, length(test) > 0L)
  base <- if (family == 'Program interactions') pattern_basis else strength_basis
  context <- context_fit(train)
  gate <- context_predict(context, baseline)
  x <- design_matrix(base, gate)
  w <- pair_weights(train)
  rms <- sqrt(colSums(x[train, , drop = FALSE]^2 * w))
  constant <- which(rms < 1e-10)
  rms[constant] <- 1
  xx <- sweep(x, 2L, rms, '/')
  penalty <- c(rep(1, ncol(base)), rep(4, 2L * ncol(base)))
  fit <- glmnet(xx[train, , drop = FALSE], d$I[train], weights = w,
    alpha = elastic_alpha, lambda = lambda_grid, intercept = FALSE,
    standardize = FALSE, penalty.factor = penalty, thresh = 1e-8)
  prediction <- as.matrix(predict(fit, xx[test, , drop = FALSE], s = lambda_grid))
  stopifnot(all(is.finite(prediction)))
  list(fit = fit, context = context, rms = rms, family = family,
    prediction = prediction, constant_columns = colnames(x)[constant],
    training_rows = train, feature_names = colnames(x))
}
select_fit <- function(available, excluded, family) {
  # Whole-cell inner folds retain viable lapatinib comparisons after outer
  # single-control purging. Outer evaluation still excludes an entire partner.
  # No held-cell combination outcome enters its inner fitting subset.
  cells <- sort(unique(d$depmap_id[available]))
  set.seed(seed)
  cell_fold <- setNames(sample(rep(1:3, length.out = length(cells))), cells)
  losses <- list()
  for (k in 1:3) {
    held_cells <- names(cell_fold)[cell_fold == k]
    test <- available[d$depmap_id[available] %in% held_cells & !is.na(d$partner[available])]
    stopifnot(length(test) >= 30L)
    train_available <- available[!d$depmap_id[available] %in% held_cells]
    train <- training_rows(excluded, test, train_available)
    stopifnot(!any(d$depmap_id[train] %in% held_cells))
    fitted <- fit_path(train, test, family)
    test_w <- pair_weights(test)
    error <- colSums((fitted$prediction - d$I[test])^2 * test_w)
    losses[[k]] <- data.table(inner_cell_fold = k, lambda = lambda_grid, objective = error,
      training_rows = length(train), test_rows = length(test))
    losses[[k]] <- rbind(losses[[k]], data.table(inner_cell_fold = k,
      lambda = Inf, objective = sum(d$I[test]^2 * test_w), training_rows = length(train), test_rows = length(test)))
  }
  table <- rbindlist(losses)
  mean_loss <- table[, .(objective = mean(objective)), by = lambda]
  best <- min(mean_loss$objective)
  eligible <- mean_loss[objective <= best * 1.01]
  chosen_lambda <- max(eligible$lambda)
  list(lambda = chosen_lambda, tuning = table)
}
predict_fitted <- function(model, a, b, state) {
  base <- if (model$family == 'Program interactions') interaction_basis(a, b) else
    matrix(sqrt(rowSums(a^2)) * sqrt(rowSums(b^2)), ncol = 1L,
      dimnames = list(NULL, 'response_strength_product'))
  x <- design_matrix(base, context_predict(model$context, state))
  x <- sweep(x, 2L, model$rms, '/')
  as.numeric(x %*% model$coefficients)
}
fit_selected <- function(train, test, family, lambda) {
  fitted <- fit_path(train, test, family)
  if (is.infinite(lambda)) {
    fitted$coefficients <- setNames(rep(0, length(fitted$feature_names)), fitted$feature_names)
  } else {
    cc <- as.matrix(coef(fitted$fit, s = lambda))[-1, 1]
    fitted$coefficients <- setNames(cc, fitted$feature_names)
  }
  fitted$lambda <- lambda
  fitted$fit <- NULL
  fitted$prediction <- NULL
  fitted
}

# Nomination, program basis and all query inputs are saved before held evaluation.
qsave(list(partners = chosen, programs = selected_programs, response_scale = response_scale,
  lambda_grid = lambda_grid, alpha = elastic_alpha, program_basis = '15 symmetric products; two baseline gates',
  inner_selection = 'Three complete-recipient folds, evaluated on lapatinib combinations',
  nomination_dose = nomination_dose, seed = seed,
  limitation = 'Conditional on frozen molecular nomination and previously inspected resources'),
  file.path(root, 'models/settings.qs'))
qsave(list(metadata = d, A = A, B = B, baseline = baseline, query = query, response = response),
  file.path(root, 'data/prepared_inputs.qs'))
write_table(d[, .(row_id, depmap_id, CELLNAME, SCREENER, pair, drugA, drugB, doseA, doseB, partner)],
  'combination_metadata')
status(paste('Nominated', length(partners), 'partners;', nrow(d), 'combination rows;',
  length(selected_programs), 'HGF-supported programs'))
print(chosen[, .(drug, role, median_source_context_cosine, combination_rows, combination_cells)])
print(program_support[selected == TRUE])

# Whole partner drug exclusions, with conservative single-reference purging.
# The primary held rows are lapatinib combinations; other selected pairs teach g.
predictions <- audits <- tuning <- list()
families <- c('Program interactions', 'Response strength control')
for (partner in partners) {
  test <- which(!is.na(d$partner) & d$partner == partner)
  train <- training_rows(partner, test)
  for (family in families) {
    status(paste('Fitting excluded partner', partner, '|', family))
    selection <- select_fit(train, partner, family)
    fitted <- fit_selected(train, test, family, selection$lambda)
    prediction <- predict_fitted(fitted, A[test, , drop = FALSE], B[test, , drop = FALSE], baseline[test, , drop = FALSE])
    stopifnot(max(abs(prediction - predict_fitted(fitted, B[test, , drop = FALSE],
      A[test, , drop = FALSE], baseline[test, , drop = FALSE]))) < 1e-10)
    tag <- paste(gsub('[^A-Za-z0-9]+', '_', partner), gsub(' ', '_', family), sep = '__')
    qsave(fitted, file.path(root, 'models', paste0('held_', tag, '.qs')))
    predictions[[tag]] <- data.table(row_id = d$row_id[test], excluded_partner = partner,
      family = family, predicted_g = prediction)
    tuning[[tag]] <- copy(selection$tuning)[, `:=`(excluded_partner = partner, family = family)]
    audits[[tag]] <- data.table(excluded_partner = partner, family = family,
      training_rows = length(train), training_pairs = uniqueN(d$pair[train]),
      test_rows = length(test), test_cells = uniqueN(d$depmap_id[test]),
      lambda = selection$lambda, retained_terms = sum(abs(fitted$coefficients) > 1e-10),
      shared_single_settings = length(intersect(c(reference_a[train], reference_b[train]),
        c(reference_a[test], reference_b[test]))))
    write_table(rbindlist(predictions), 'held_predictions_frozen')
    write_table(rbindlist(audits), 'fit_audit')
  }
}

# Final parameter selection uses complete-cell development folds only.
final <- list()
for (family in families) {
  status(paste('Final development fit |', family))
  selection <- select_fit(seq_len(nrow(d)), character(), family)
  fitted <- fit_selected(seq_len(nrow(d)), seq_len(nrow(d)), family, selection$lambda)
  final[[family]] <- fitted
  qsave(fitted, file.path(root, 'models', paste0('final_', gsub(' ', '_', family), '.qs')))
  tuning[[paste0('final_', family)]] <- copy(selection$tuning)[,
    `:=`(excluded_partner = 'Full development fit', family = family)]
  write_table(data.table(family = family, term = names(fitted$coefficients),
    coefficient_scaled_feature = fitted$coefficients, feature_rms = fitted$rms,
    coefficient_original_product = fitted$coefficients / fitted$rms),
    paste0('final_coefficients_', gsub(' ', '_', family)))
}
write_table(rbindlist(tuning), 'inner_parameter_selection')

# Freeze the HGF prediction before attaching combination evaluation labels.
au_baseline <- getmatrix(PREP$raw, AU, programs)
hgf_predictions <- list()
support <- list()
for (dose in query_doses) {
  l_ensemble <- Reduce('+', lapply(models, function(model)
    uppred(PREP, model$operators[[anchor]], AU, dose))) / length(models)
  aa <- matrix(rep(as.numeric(l_ensemble / response_scale), each = nrow(hgf)), nrow(hgf), 40L)
  bb <- sweep(as.matrix(hgf[, programs, with = FALSE]), 2L, response_scale, '/')
  state <- au_baseline[rep(1L, nrow(hgf)), , drop = FALSE]
  for (family in families) {
    g <- predict_fitted(final[[family]], aa, bb, state)
    stopifnot(max(abs(predict_fitted(final[[family]], aa * 0, bb, state))) == 0,
      max(abs(g - predict_fitted(final[[family]], bb, aa, state))) < 1e-10)
    for (k in seq_along(models)) {
      model <- models[[k]]
      l <- uppred(PREP, model$operators[[anchor]], AU, dose)
      y_l <- predict_single(model, l, au_baseline, 'secondary')
      y_h <- predict_single(model, as.matrix(hgf[, programs, with = FALSE]), state, 'secondary')
      hgf_predictions[[length(hgf_predictions) + 1L]] <- data.table(family = family,
        original_fold = k, dose = dose, hgf_context = hgf$context,
        lapatinib_alone_pct = 100 * 2^y_l, hgf_modifier_log2 = g,
        drug_effect_pct_of_hgf_environment = 100 * 2^(y_l + g),
        hgf_alone_extrapolated_pct = 100 * 2^y_h,
        combined_extrapolated_pct = 100 * 2^(y_l + y_h + g))
    }
    gate <- context_predict(final[[family]]$context, au_baseline)
    support[[length(support) + 1L]] <- data.table(family = family, dose = dose,
      baseline_PC1_gate = gate[1, 1], baseline_PC2_gate = gate[1, 2],
      query_response_metric = 'Native ensemble displacement / frozen within-drug-dose molecular SD',
      transfer = 'ALMANAC interaction applied to PRISM; no assay-specific interaction calibration')
  }
}
hgf_prediction <- rbindlist(hgf_predictions)
stopifnot(all(is.finite(hgf_prediction$drug_effect_pct_of_hgf_environment)))
write_table(hgf_prediction, 'hgf_predictions_frozen')
write_table(rbindlist(support), 'hgf_query_support')
frozen_paths <- file.path(root, 'data', c('held_predictions_frozen.tsv', 'hgf_predictions_frozen.tsv'))
write_table(data.table(path = frozen_paths,
  sha256 = vapply(frozen_paths, digest, character(1), file = TRUE, algo = 'sha256')), 'prediction_freeze_hashes')

# Evaluate departures from measured constituent effects, not pooled drug killing.
held <- merge(rbindlist(predictions), d, by = 'row_id', all.x = TRUE, sort = FALSE)
stopifnot(nrow(held) == sum(!is.na(d$partner)) * length(families), !anyNA(held$I))
held[, predicted_abundance_pct := 100 * 2^(yA + yB + predicted_g)]
held[, independent_abundance_pct := 100 * 2^(yA + yB)]
held[, observed_abundance_pct := 100 * 2^y]
held[, observed_departure_pp := observed_abundance_pct - independent_abundance_pct]
held[, predicted_departure_pp := predicted_abundance_pct - independent_abundance_pct]
safe_cor <- function(x, y) if (length(x) >= 4L && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
metrics <- held[, .(rows = .N, cells = uniqueN(depmap_id),
  pearson_interaction_log2 = safe_cor(I, predicted_g),
  pearson_departure_pp = safe_cor(observed_departure_pp, predicted_departure_pp),
  observed_positive_modification = sum(I > .25),
  observed_negative_modification = sum(I < -.25),
  correct_direction_large_effect = sum(sign(I) == sign(predicted_g) & abs(I) > .25),
  median_observed_departure_pp = median(observed_departure_pp),
  median_predicted_departure_pp = median(predicted_departure_pp)), by = .(family, excluded_partner)]
pooled <- held[, .(rows = .N, cells = uniqueN(depmap_id), partners = uniqueN(excluded_partner),
  pearson_interaction_log2 = safe_cor(I, predicted_g),
  pearson_departure_pp = safe_cor(observed_departure_pp, predicted_departure_pp)), by = family]
within_treatment <- held[, .(cells = uniqueN(depmap_id),
  pearson_interaction_log2 = safe_cor(I, predicted_g)),
  by = .(family, excluded_partner, doseA, doseB, SCREENER)]
write_table(held, 'held_predictions_evaluated')
write_table(metrics, 'held_partner_metrics')
write_table(pooled, 'pooled_held_metrics')
write_table(within_treatment, 'fixed_pair_dose_across_cell_metrics')
write_table(rbindlist(audits), 'fit_audit')

# Exact protected monotherapy replay: same441 recipients and original predictions.
mono <- fread(file.path(reference, 'inputs/results/predictions.tsv'))[model == 'full' & dose == 2.5]
stopifnot(nrow(mono) == 441L, !anyDuplicated(mono$depmap_id))
for (k in 1:3) {
  ii <- which(mono$fold == k)
  model <- models[[k]]
  dl <- uppred(PREP, model$operators[[anchor]], mono$depmap_id[ii], mono$dose[ii])
  ss <- getmatrix(PREP$raw, mono$depmap_id[ii], programs)
  stopifnot(max(abs(predict_single(model, dl, ss, 'secondary') - mono$predicted[ii])) < 1e-10)
}
write_table(mono, 'protected_441_monotherapy_predictions')
write_table(data.table(cells = 441L, pearson_log2 = cor(mono$observed, mono$predicted),
  observed_sensitive = sum(mono$observed <= -1), observed_above_50pct = sum(mono$observed > -1)),
  'protected_monotherapy_metric')

# Essential plots: held interaction predictions and the HGF treatment contrast.
labels <- copy(metrics)
labels[, panel := paste0(sub(' \\(.*', '', excluded_partner), '\nn=', rows,
  '; r=', ifelse(is.na(pearson_interaction_log2), 'undefined', sprintf('%.2f', pearson_interaction_log2)))]
plot_rows <- merge(held[family == 'Program interactions'], labels[family == 'Program interactions',
  .(excluded_partner, panel)], by = 'excluded_partner', all.x = TRUE)
p <- ggplot(plot_rows, aes(I, predicted_g)) + geom_point(size = .9, alpha = .5, color = '#2166AC') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
  facet_wrap(~panel, ncol = 3, scales = 'free') + theme_classic(base_size = 9) +
  theme(aspect.ratio = 1, strip.background = element_blank(), strip.text = element_text(size = 8),
    plot.caption = element_text(hjust = 0, size = 8)) +
  labs(x = 'Observed departure from separate effects (log2)', y = 'Predicted departure (log2)',
    title = 'Lapatinib combinations: each partner excluded from interaction training',
    caption = 'Whole partner exclusions; reused single-control settings purged.\nKnown single-response inputs; retrospective development test, conditional on the frozen molecular shortlist.')
ggsave(file.path(root, 'figures/held_partner_interactions.pdf'), p, width = 8, height = 8.8, device = cairo_pdf)
ggsave(file.path(root, 'figures/held_partner_interactions.png'), p, width = 8, height = 8.8, dpi = 150)
q <- hgf_prediction[hgf_context == 'Equal-context mean' & family == 'Program interactions']
q <- rbind(q[, .(original_fold, dose, condition = 'Lapatinib alone', remaining = lapatinib_alone_pct)],
  q[, .(original_fold, dose, condition = 'Lapatinib in HGF', remaining = drug_effect_pct_of_hgf_environment)])
q[, condition := factor(condition, levels = c('Lapatinib alone', 'Lapatinib in HGF'))]
q[, exposure := paste(dose, 'µM')]
hp <- ggplot(q, aes(condition, remaining, group = factor(original_fold), color = factor(original_fold))) +
  geom_line(linewidth = .7) + geom_point(size = 2) + facet_wrap(~exposure, nrow = 1) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom',
    plot.caption = element_text(hjust = 0, size = 8)) +
  scale_color_manual(values = c('#2166AC', '#B2182B', '#477A41'), name = 'Original B1 fit') +
  labs(x = NULL, y = 'Drug-treated abundance (% corresponding environment)',
    title = 'AU565: HGF modifier learned from other drug combinations',
    caption = 'HGF coordinates: equal mean of all12 transported source contexts.\nDrug-alone control is untreated; HGF control is HGF alone. No HGF functional label fitted; assay transfer uncalibrated.')
ggsave(file.path(root, 'figures/hgf_conditional_prediction.pdf'), hp, width = 7.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(root, 'figures/hgf_conditional_prediction.png'), hp, width = 7.5, height = 4.5, dpi = 150)

inputs <- c(file.path(original, 'modular_pair_pilot/inputs.rds'), hgf_path, weights_path,
  file.path(original, 'modular_pair_pilot/model.R'), file.path(original, 'modular_pair_pilot/bank_operator.R'),
  file.path(workspace, 'src/program_pilot/75_hgf_pattern_partner_transfer.R'))
write_table(data.table(path = inputs,
  sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')), 'source_hashes')
writeLines(capture.output(sessionInfo()), file.path(root, 'session_info.txt'))
writeLines(c('Original40 molecular-pattern partner-transfer experiment. Original B1 monotherapy protected.',
  'Nomination is molecular-only; closest drugs are pattern neighbors, not established HGF agonists.',
  'Five HGF-supported programs;15 symmetric products with two training-baseline gates; sparse regularization.',
  'Whole partner outcomes excluded; conservative screen/cell/drug/dose single-control purging.',
  'Interaction inputs use frozen ensemble baseline-predicted displacements. Known partner single RNA used for nomination.',
  'Primary held comparison: departures from matched separate-agent effects for lapatinib combinations.',
  'HGF remains outside functional fitting;12 transported source contexts and mean preserved.',
  'Applying ALMANAC g to the PRISM readout assumes interaction transfer across assays; no new calibration.',
  'Individual HGF-alone and absolute combined predictions are extrapolations. Conditional treatment contrasts are primary.',
  'Previously inspected outcomes and fixed molecular bank make this retrospective development evidence.'),
  file.path(root, 'REPORT.txt'))
print(pooled)
print(metrics[, .(excluded_partner, family, rows, pearson_interaction_log2)])
print(hgf_prediction[hgf_context == 'Equal-context mean',
  .(family, original_fold, dose, lapatinib_alone_pct, hgf_modifier_log2, drug_effect_pct_of_hgf_environment)])
status('Complete')
