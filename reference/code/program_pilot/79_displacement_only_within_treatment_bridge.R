#!/usr/bin/env Rscript
# Strict molecular bridge: baseline is used only by molecular displacement models.
# Functional sensitivity is fitted to within-treatment contrasts, separately from
# average-effect calibration. Molecular mean-only models are not candidates.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(glmnet)
  library(digest)
  library(ggplot2)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

root <- '/work/results/program_pilot/79_displacement_only_within_treatment_bridge'
previous <- '/work/results/program_pilot/78_native_dose_prism_functional_interpolation'
reference <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query'
original <- file.path(reference, 'inputs/original')
script_path <- '/work/src/program_pilot/79_displacement_only_within_treatment_bridge.R'
native_doses <- c(.05, .5, 5)
lambda_grid <- 10^seq(0, -3, length.out = 12L)
molecular_lambdas <- c(100, 10, 1, .1)
average_lambda <- .1
nonlinear_selection_margin <- .02
seed <- 77012L
anchor <- 'Lapatinib ditosylate'
candidates <- c('Linear displacement', 'Nonlinear displacement')
set.seed(seed)
for (folder in c('data', 'models', 'figures', 'models/molecular_cache')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
options(error = function() {
  status('FAILED; see run log')
  traceback(5)
  quit(save = 'no', status = 1L)
})
pearson <- function(x, y) {
  if (length(x) < 4L || sd(x) < 1e-10 || sd(y) < 1e-10) return(NA_real_)
  cor(x, y)
}

status('Loading frozen v78 dose alignment; no baseline variables in the functional design')
setwd(original)
source('modular_pair_pilot/bank_operator.R')
archive <- readRDS('modular_pair_pilot/inputs.rds')
PREP <- archive$prep
DRUGS <- archive$drugs
settings78 <- qread(file.path(previous, 'models/settings.qs'))
selected_drugs <- settings78$selected_drugs
DRUGS <- DRUGS[DRUGS %in% selected_drugs]
PREP$m <- PREP$m[PREP$m$drug %in% selected_drugs, ]
stopifnot(length(DRUGS) == 10L, !AU %in% PREP$m$depmap_id)
rm(archive)

# The upstream ridge models learn baseline -> measured program displacement.
# Selection assesses centered molecular responses within each native dose. A
# dose's common footprint cannot improve this criterion; no mean-only candidate
# is offered. Shrinkage is still allowed and is not evidence of genuine signal.
upchoose <- function(prep, m) {
  stopifnot(length(unique(m$depmap_id)) >= 6L)
  fold <- fmget(foldmap(m$depmap_id, 3L, 20261101), m$depmap_id)
  observed <- as.matrix(m[P])
  losses <- numeric(length(molecular_lambdas))
  for (j in seq_along(molecular_lambdas)) {
    for (k in 1:3) {
      train <- fold != k
      fit <- upfit(prep, m[train, ], 'programs', molecular_lambdas[j])
      predicted <- uppred(prep, fit, m$depmap_id[!train], m$dose[!train])
      train_residual <- observed[train, , drop = FALSE]
      for (dose in native_doses) {
        ii <- which(m$dose[train] == dose)
        train_residual[ii, ] <- sweep(train_residual[ii, , drop = FALSE], 2L,
          colMeans(train_residual[ii, , drop = FALSE]))
      }
      scale <- apply(train_residual, 2L, sd)
      scale <- pmax(scale, quantile(scale, .1), 1e-5)
      error <- observed[!train, , drop = FALSE] - predicted
      for (dose in native_doses) {
        ii <- which(m$dose[!train] == dose)
        error[ii, ] <- sweep(error[ii, , drop = FALSE], 2L, colMeans(error[ii, , drop = FALSE]))
      }
      losses[j] <- losses[j] + mean(sweep(error, 2L, scale, '/')^2) / 3
    }
  }
  pick <- which(losses <= min(losses) * 1.01)[1]
  fit <- upfit(prep, m, 'programs', molecular_lambdas[pick])
  fit$selection <- data.frame(kind = 'programs', lambda = molecular_lambdas,
    within_dose_internal_objective = losses)
  fit
}

# Reuse v78's functional measurements and interpolation, without reading the raw
# PRISM curves again. All original measured endpoints and all441 B1 cells remain.
SINGLE <- as.data.frame(fread(file.path(previous, 'data/individual_training_outcomes.tsv')))
d <- fread(file.path(previous, 'data/combination_training_outcomes.tsv'))
evaluation <- fread(file.path(previous, 'data/held_single_predictions_frozen.tsv'))[
  family == 'Joint singles and combinations', .(single_id, depmap_id, drug, dose, observed_log2)]
setorder(evaluation, single_id)
fold_table <- fread(file.path(previous, 'data/recipient_folds.tsv'))
recipient_fold <- setNames(fold_table$fold, fold_table$depmap_id)
stopifnot(nrow(d) == 5640L, uniqueN(d$pair) == 42L, !anyDuplicated(evaluation$single_id),
  sum(evaluation$drug == anchor & evaluation$dose == 2.5) == 441L,
  !AU %in% c(SINGLE$depmap_id, d$depmap_id))
PAIR <- as.data.frame(d)
entities <- unique(rbind(as.data.table(SINGLE)[, .(depmap_id, drug)],
  d[, .(depmap_id, drug = drugA)], d[, .(depmap_id, drug = drugB)], evaluation[, .(depmap_id, drug)]))
native_query <- entities[, .(dose = native_doses), by = .(depmap_id, drug)]
key <- function(cell, drug, dose) paste(cell, drug, sprintf('%.10g', dose), sep = '|')
native_key <- key(native_query$depmap_id, native_query$drug, native_query$dose)
indices <- function(cell, drug) {
  ix <- vapply(native_doses, function(dose) match(key(cell, drug, dose), native_key), integer(length(cell)))
  stopifnot(!anyNA(ix))
  ix
}
weights <- function(dose) {
  stopifnot(all(dose >= .05 - 1e-10), all(dose <= 5 + 1e-10))
  w <- matrix(0, length(dose), 3L)
  for (i in seq_along(dose)) {
    lo <- max(which(native_doses <= dose[i] + 1e-10))
    hi <- min(which(native_doses >= dose[i] - 1e-10))
    if (lo == hi) w[i, lo] <- 1 else {
      fraction <- log(dose[i] / native_doses[lo]) / log(native_doses[hi] / native_doses[lo])
      w[i, lo] <- 1 - fraction
      w[i, hi] <- fraction
    }
  }
  stopifnot(max(abs(rowSums(w) - 1)) < 1e-12)
  w
}
is <- indices(SINGLE$depmap_id, SINGLE$drug)
ia <- indices(d$depmap_id, d$drugA)
ib <- indices(d$depmap_id, d$drugB)
ie <- indices(evaluation$depmap_id, evaluation$drug)
ws <- weights(SINGLE$dose)
wa <- weights(d$doseA)
wb <- weights(d$doseB)
we <- weights(evaluation$dose)

# Nonlinear sensitivity uses every program pair, rather than restricting nonlinear
# effects to five HGF-nominated dimensions. Both candidates have no baseline gates.
ij <- which(upper.tri(matrix(0, 40L, 40L), diag = TRUE), arr.ind = TRUE)
basis_names <- c(P, paste(P[ij[, 1]], P[ij[, 2]], sep = '*'))
assay_deviations <- c('primary', 'nci_FF', 'nci_FG')
feature_names <- c(basis_names, unlist(lapply(assay_deviations, function(a) paste(a, P, sep = ':'))))
linear_columns <- c(seq_len(40L), length(basis_names) + seq_len(120L))
penalty <- c(rep(1, 40L), rep(4, nrow(ij)), rep(8, 120L))
response_basis <- function(response, context) {
  z <- sweep(response, 2L, context$response_scale, '/')
  quadratic <- z[, ij[, 1], drop = FALSE] * z[, ij[, 2], drop = FALSE]
  result <- cbind(z, quadratic)
  colnames(result) <- basis_names
  result
}
design <- function(basis, assay) {
  result <- cbind(basis, do.call(cbind, lapply(assay_deviations, function(a)
    basis[, seq_len(40L), drop = FALSE] * as.numeric(assay == a))))
  colnames(result) <- feature_names
  result
}
interpolate <- function(basis, ix, w) {
  result <- basis[ix[, 1], , drop = FALSE] * w[, 1]
  for (i in 2:3) result <- result + basis[ix[, i], , drop = FALSE] * w[, i]
  result
}
sum_basis <- function(a, b) {
  result <- a + b
  result[, 40L + seq_len(nrow(ij))] <- result[, 40L + seq_len(nrow(ij))] +
    a[, ij[, 1], drop = FALSE] * b[, ij[, 2], drop = FALSE] +
    a[, ij[, 2], drop = FALSE] * b[, ij[, 1], drop = FALSE]
  result
}

ns <- nrow(SINGLE)
np <- nrow(d)
pair_group <- paste(d$SCREENER, d$pair, d$doseA, d$doseB, sep = '|')
metadata <- rbindlist(list(
  as.data.table(SINGLE)[, .(depmap_id, source, group, control_key, kind = 'single', y)],
  d[, .(depmap_id, source = 'combination_total', group = pair_group,
    control_key = paste0('combination:', row_id), kind = 'combination', y)],
  d[, .(depmap_id, source = 'combination_departure', group = pair_group,
    control_key = paste0('departure:', row_id), kind = 'departure', y = I)]))
metadata[, center_group := paste(kind, group, sep = '|')]
source_totals <- settings78$source_totals
row_weights <- function(rows) {
  md <- metadata[rows]
  w <- numeric(nrow(md))
  for (src in unique(md$source)) {
    ii <- which(md$source == src)
    z <- md[ii]
    unit <- 1 / as.numeric(table(z$control_key)[z$control_key])
    totals <- tapply(unit, z$group, sum)
    w[ii] <- source_totals[src] * unit / as.numeric(totals[z$group]) / uniqueN(z$group)
  }
  stopifnot(all(w > 0), all(is.finite(w)))
  w / sum(w)
}
training_rows <- function(si, pi) c(si, ns + pi, ns + np + pi)
ref_a <- paste(d$SCREENER, d$depmap_id, d$drugA, sprintf('%.10g', d$doseA), sep = '|')
ref_b <- paste(d$SCREENER, d$depmap_id, d$drugB, sprintf('%.10g', d$doseB), sep = '|')
pair_treatments <- unique(d[, .(drugA, drugB, doseA, doseB)])
treatment_key <- function(a, b, da, db) paste(a, b, da, db, sep = '|')
pair_reference_index <- match(treatment_key(d$drugA, d$drugB, d$doseA, d$doseB),
  treatment_key(pair_treatments$drugA, pair_treatments$drugB, pair_treatments$doseA, pair_treatments$doseB))

prepare <- function(excluded) {
  molecular_ids <- sort(intersect(excluded, unique(PREP$m$depmap_id)))
  cache_hash <- digest(list(molecular_ids, native_key, molecular_lambdas), algo = 'sha256')
  cache_path <- file.path(root, 'models/molecular_cache', paste0(cache_hash, '.qs'))
  if (file.exists(cache_path)) {
    desc <- qread(cache_path)
    stopifnot(identical(desc$cache_hash, cache_hash))
  } else {
    desc <- build_descriptors(PREP, as.data.frame(native_query), as.data.frame(d[0]),
      excluded = excluded, crossfit = TRUE)
    desc$cache_hash <- cache_hash
    qsave(desc, cache_path)
  }
  stopifnot(all(vapply(desc$operators, function(op)
    op$kind == 'programs' && !any(op$training_ids %in% excluded), logical(1))))
  molecular <- as.data.table(PREP$m)[!depmap_id %in% excluded]
  centered <- molecular[, lapply(.SD, function(x) x - mean(x)), by = .(drug, dose), .SDcols = P]
  response_scale <- vapply(centered[, P, with = FALSE], sd, numeric(1))
  stopifnot(all(is.finite(response_scale)), all(response_scale > 0))
  ref_ids <- sort(unique(molecular$depmap_id))
  context <- list(response_scale = response_scale, reference_cells = ref_ids,
    excluded_molecular_recipients = excluded)
  basis <- response_basis(desc$X, context)
  bs <- interpolate(basis, is, ws)
  ba <- interpolate(basis, ia, wa)
  bb <- interpolate(basis, ib, wb)
  xsum <- design(sum_basis(ba, bb), paste0('nci_', d$SCREENER))
  xa <- design(ba, paste0('nci_', d$SCREENER))
  xb <- design(bb, paste0('nci_', d$SCREENER))
  xs <- design(bs, SINGLE$assay)
  xe <- design(interpolate(basis, ie, we), rep('secondary', nrow(evaluation)))

  # Reference inputs are molecular predictions in the training RNA-cell population,
  # not baseline features in the functional decoder. Their moments are constant
  # across receiving cells at a fixed treatment and dose, including new partners.
  ref_query <- CJ(depmap_id = ref_ids, drug = selected_drugs, dose = native_doses)
  ref_response <- base_predict(PREP, desc$operators, as.data.frame(ref_query))
  ref_basis <- response_basis(ref_response, context)
  ref_key <- key(ref_query$depmap_id, ref_query$drug, ref_query$dose)
  ref_by_drug <- lapply(selected_drugs, function(dr) lapply(native_doses, function(dose) {
    ix <- match(key(ref_ids, dr, dose), ref_key)
    stopifnot(!anyNA(ix))
    ref_basis[ix, , drop = FALSE]
  }))
  names(ref_by_drug) <- selected_drugs
  reference_at <- function(drug, dose) {
    w <- weights(dose)
    result <- matrix(0, length(ref_ids), length(basis_names))
    for (i in 1:3) result <- result + ref_by_drug[[drug]][[i]] * w[1, i]
    result
  }
  single_reference <- function(drug, dose) {
    unique_treatments <- unique(data.table(drug, dose))
    mean_matrix <- t(vapply(seq_len(nrow(unique_treatments)), function(i)
      colMeans(reference_at(unique_treatments$drug[i], unique_treatments$dose[i])), numeric(length(basis_names))))
    mean_matrix[match(paste(drug, dose), paste(unique_treatments$drug, unique_treatments$dose)), , drop = FALSE]
  }
  rxs <- design(single_reference(SINGLE$drug, SINGLE$dose), SINGLE$assay)
  rxa <- design(single_reference(d$drugA, d$doseA), paste0('nci_', d$SCREENER))
  rxb <- design(single_reference(d$drugB, d$doseB), paste0('nci_', d$SCREENER))
  ref_pair <- t(vapply(seq_len(nrow(pair_treatments)), function(i) {
    t <- pair_treatments[i]
    colMeans(sum_basis(reference_at(t$drugA, t$doseA), reference_at(t$drugB, t$doseB)))
  }, numeric(length(basis_names))))
  rxsum <- design(ref_pair[pair_reference_index, , drop = FALSE], paste0('nci_', d$SCREENER))
  rxe <- design(single_reference(evaluation$drug, evaluation$dose), rep('secondary', nrow(evaluation)))
  departure <- xsum - xa - xb
  reference_departure <- rxsum - rxa - rxb
  stopifnot(max(abs(departure[, linear_columns])) < 1e-8,
    max(abs(reference_departure[, linear_columns])) < 1e-8,
    max(abs(design(response_basis(matrix(0, 1L, 40L), context), 'secondary'))) == 0)
  list(x = rbind(xs, xsum, departure), rx = rbind(rxs, rxsum, reference_departure),
    xe = xe, rxe = rxe, xa = xa, xb = xb, xsum = xsum, rxa = rxa, rxb = rxb, rxsum = rxsum,
    context = context, operators = desc$operators, ref_by_drug = ref_by_drug)
}

# Sensitivity coefficients see only weighted within-context contrasts. Means are
# estimated using training rows exclusively. Average calibration is a second fit
# on group means and cannot change sensitivity coefficients or model selection.
fit_path <- function(prepared, rows, candidate) {
  columns <- if (candidate == candidates[1]) linear_columns else seq_along(feature_names)
  x <- prepared$x[rows, columns, drop = FALSE]
  rx <- prepared$rx[rows, columns, drop = FALSE]
  w <- row_weights(rows)
  group <- match(metadata$center_group[rows], unique(metadata$center_group[rows]))
  group_weight <- as.numeric(rowsum(w, group, reorder = FALSE))
  means <- rowsum(x * w, group, reorder = FALSE) / group_weight
  reference_means <- rowsum(rx * w, group, reorder = FALSE) / group_weight
  observed <- metadata$y[rows]
  ymean <- as.numeric(rowsum(observed * w, group, reorder = FALSE)) / group_weight
  centered_x <- x - means[group, , drop = FALSE]
  centered_y <- observed - ymean[group]
  stopifnot(max(abs(rowsum(centered_x * w, group, reorder = FALSE))) < 1e-8,
    max(abs(as.numeric(rowsum(centered_y * w, group, reorder = FALSE)))) < 1e-8)
  rms <- sqrt(as.numeric(crossprod(w, centered_x^2)))
  keep <- is.finite(rms) & rms > 1e-10
  yscale <- sqrt(sum(w * centered_y^2))
  stopifnot(any(keep), is.finite(yscale), yscale > 0)
  fit <- glmnet(centered_x[, keep, drop = FALSE] %*% Diagonal(x = 1 / rms[keep]),
    centered_y / yscale, alpha = 0, lambda = lambda_grid, weights = w,
    intercept = FALSE, standardize = FALSE, penalty.factor = penalty[columns][keep],
    thresh = 1e-7, maxit = 100000L)
  stopifnot(fit$jerr == 0L)
  beta <- matrix(0, length(feature_names), length(lambda_grid), dimnames = list(feature_names, NULL))
  beta[columns[keep], ] <- sweep(as.matrix(coef(fit, s = lambda_grid))[-1, , drop = FALSE] * yscale,
    1L, rms[keep], '/')
  adjusted_means <- matrix(ymean, length(ymean), length(lambda_grid)) -
    (means - reference_means) %*% beta[columns, , drop = FALSE]
  mean_rms <- sqrt(as.numeric(crossprod(group_weight, reference_means^2)))
  mean_keep <- is.finite(mean_rms) & mean_rms > 1e-10
  mean_scale <- sqrt(sum(group_weight * rowMeans(adjusted_means^2)))
  average_fit <- glmnet(reference_means[, mean_keep, drop = FALSE] %*%
      Diagonal(x = 1 / mean_rms[mean_keep]), adjusted_means / mean_scale,
    family = 'mgaussian', alpha = 0, lambda = average_lambda, weights = group_weight,
    intercept = FALSE, standardize = FALSE, penalty.factor = penalty[columns][mean_keep],
    thresh = 1e-7, maxit = 100000L)
  stopifnot(average_fit$jerr == 0L)
  gamma <- beta * 0
  gamma[columns[mean_keep], ] <- sweep(vapply(coef(average_fit), function(b)
    as.numeric(b[-1, 1]), numeric(sum(mean_keep))) * mean_scale, 1L, mean_rms[mean_keep], '/')
  list(beta = beta, gamma = gamma, training_rows = rows,
    zero_sensitivity_columns = feature_names[columns[!keep]])
}
predict_path <- function(x, rx, fit) (x - rx) %*% fit$beta + rx %*% fit$gamma
selection_score <- function(rows, predicted) {
  md <- metadata[rows]
  scores <- vapply(seq_len(ncol(predicted)), function(j) {
    z <- md[, .(source, center_group, y)]
    z[, predicted := predicted[, j]]
    r <- z[, .(r = pearson(y, predicted)), by = .(source, center_group)]
    # Undefined correlations contribute zero; no samples or low-variation groups
    # are removed from fitting/evaluation and no arbitrary variation cutoff is used.
    r[!is.finite(r), r := 0]
    by_source <- r[, .(r = mean(r)), by = source]
    sum(by_source$r * source_totals[by_source$source]) / sum(source_totals[by_source$source])
  }, numeric(1))
  scores
}

jobs <- c(lapply(1:3, function(k) list(tag = paste0('recipient_', k), regime = 'recipient', fold = k,
  held_cells = names(recipient_fold)[recipient_fold == k], held_partner = NA_character_)),
  lapply(setdiff(selected_drugs, anchor), function(dr) list(tag = paste0('partner_', gsub('[^A-Za-z0-9]+', '_', dr)),
    regime = 'partner', fold = NA_integer_, held_cells = character(), held_partner = dr)),
  list(list(tag = 'final', regime = 'final', fold = NA_integer_, held_cells = character(), held_partner = NA_character_)))
single_predictions <- pair_predictions <- selection_tables <- operator_tables <- audits <- list()
final_models <- list()
run_contract <- digest(list(digest(script_path, file = TRUE, algo = 'sha256'),
  molecular_lambdas, lambda_grid, average_lambda, nonlinear_selection_margin, feature_names), algo = 'sha256')
for (job in jobs) {
  result_path <- file.path(root, 'models', paste0(job$tag, '_completed.qs'))
  if (file.exists(result_path)) {
    result <- qread(result_path)
    stopifnot(identical(result$run_contract, run_contract))
    status(paste(job$tag, '| restored completed fit'))
  } else {
    status(paste(job$tag, '| fitting molecular-only bridge; functional baseline is excluded'))
    si <- which(!SINGLE$depmap_id %in% job$held_cells)
    pi <- which(!d$depmap_id %in% job$held_cells)
    test_pair <- if (job$regime == 'recipient') which(d$depmap_id %in% job$held_cells) else
      if (job$regime == 'partner') which(d$partner == job$held_partner & !is.na(d$partner)) else integer()
    if (job$regime == 'partner') {
      pi <- pi[!d$drugA[pi] %in% job$held_partner & !d$drugB[pi] %in% job$held_partner]
      held_refs <- unique(c(ref_a[test_pair], ref_b[test_pair]))
      pi <- pi[!ref_a[pi] %in% held_refs & !ref_b[pi] %in% held_refs]
      si <- si[!SINGLE$control_key[si] %in% held_refs]
      stopifnot(!any(job$held_partner %in% c(d$drugA[pi], d$drugB[pi])),
        !length(intersect(c(ref_a[pi], ref_b[pi], SINGLE$control_key[si]), held_refs)))
    }
    excluded <- union(AU, job$held_cells)
    ids <- sort(unique(c(SINGLE$depmap_id[si], d$depmap_id[pi])))
    inner_fold <- foldmap(ids, 3L, seed + 1L)
    scores <- matrix(0, length(lambda_grid), length(candidates), dimnames = list(NULL, candidates))
    for (k in 1:3) {
      status(paste(job$tag, '| training validation', k, 'of3'))
      held <- names(inner_fold)[inner_fold == k]
      inner_si <- si[!SINGLE$depmap_id[si] %in% held]
      inner_pi <- pi[!d$depmap_id[pi] %in% held]
      val_si <- si[SINGLE$depmap_id[si] %in% held]
      val_pi <- pi[d$depmap_id[pi] %in% held]
      prepared <- prepare(union(excluded, held))
      val <- training_rows(val_si, val_pi)
      for (candidate in candidates) {
        fitted <- fit_path(prepared, training_rows(inner_si, inner_pi), candidate)
        # Calibration is context-constant and cannot affect these within-context
        # correlations. Select using sensitivity predictions alone explicitly.
        predicted <- prepared$x[val, , drop = FALSE] %*% fitted$beta
        scores[, candidate] <- scores[, candidate] + selection_score(val, predicted) / 3
      }
      rm(prepared, fitted, predicted)
      gc(verbose = FALSE)
    }
    chosen_lambda <- vapply(candidates, function(candidate) {
      s <- scores[, candidate]
      which(s >= max(s) - .01)[1]
    }, integer(1))
    best <- vapply(candidates, function(candidate) scores[chosen_lambda[candidate], candidate], numeric(1))
    chosen_candidate <- if (best[2] > best[1] + nonlinear_selection_margin) candidates[2] else candidates[1]
    selection <- rbindlist(lapply(candidates, function(candidate) data.table(job = job$tag, candidate,
      lambda = lambda_grid, training_validation_within_treatment_pearson = scores[, candidate],
      selected_lambda = seq_along(lambda_grid) == chosen_lambda[candidate],
      selected_architecture = candidate == chosen_candidate)))
    status(paste(job$tag, '| final training fit; selected', chosen_candidate))
    prepared <- prepare(excluded)
    out_single <- out_pair <- models <- list()
    for (candidate in candidates) {
      fitted <- fit_path(prepared, training_rows(si, pi), candidate)
      j <- chosen_lambda[candidate]
      model <- list(candidate = candidate, beta = fitted$beta[, j], gamma = fitted$gamma[, j],
        context = prepared$context, operators = prepared$operators, lambda = lambda_grid[j],
        average_lambda = average_lambda, feature_names = feature_names, native_doses = native_doses,
        held_cells = job$held_cells, held_partner = job$held_partner,
        training_single_rows = si, training_pair_rows = pi, reference_basis = prepared$ref_by_drug,
        functional_inputs = 'Displacements and molecular-reference moments only; no baseline inputs',
        selected_from_training = candidate == chosen_candidate)
      model_path <- file.path(root, 'models', paste0(job$tag, '__', gsub(' ', '_', candidate), '.qs'))
      qsave(model, model_path)
      models[[candidate]] <- model
      if (job$regime == 'recipient') {
        test <- which(evaluation$depmap_id %in% job$held_cells)
        out <- copy(evaluation[test])
        out[, `:=`(fold = job$fold, family = candidate,
          predicted_log2 = as.numeric(predict_path(prepared$xe[test, , drop = FALSE],
            prepared$rxe[test, , drop = FALSE], fitted)[, j]),
          sensitivity_log2 = as.numeric((prepared$xe[test, , drop = FALSE] -
            prepared$rxe[test, , drop = FALSE]) %*% model$beta),
          average_log2 = as.numeric(prepared$rxe[test, , drop = FALSE] %*% model$gamma))]
        out_single[[candidate]] <- out
      }
      if (length(test_pair)) {
        out <- d[test_pair, .(row_id, depmap_id, CELLNAME, pair, drugA, drugB, doseA, doseB, SCREENER,
          partner, observed_log2 = y, observed_departure_log2 = I)]
        out[, `:=`(regime = job$regime, excluded_partner = job$held_partner, family = candidate,
          predicted_log2 = as.numeric(predict_path(prepared$xsum[test_pair, , drop = FALSE],
            prepared$rxsum[test_pair, , drop = FALSE], fitted)[, j]),
          predicted_singleA_log2 = as.numeric(predict_path(prepared$xa[test_pair, , drop = FALSE],
            prepared$rxa[test_pair, , drop = FALSE], fitted)[, j]),
          predicted_singleB_log2 = as.numeric(predict_path(prepared$xb[test_pair, , drop = FALSE],
            prepared$rxb[test_pair, , drop = FALSE], fitted)[, j]),
          average_log2 = as.numeric(prepared$rxsum[test_pair, , drop = FALSE] %*% model$gamma))]
        out[, predicted_departure_log2 := predicted_log2 - predicted_singleA_log2 - predicted_singleB_log2]
        if (candidate == candidates[1]) stopifnot(max(abs(out$predicted_departure_log2)) < 1e-8)
        out_pair[[candidate]] <- out
      }
    }
    if (length(out_single)) {
      selected <- copy(out_single[[chosen_candidate]])
      selected[, family := 'Selected from training']
      average <- copy(selected)
      average[, `:=`(family = 'Average effects only', predicted_log2 = average_log2, sensitivity_log2 = 0)]
      out_single <- c(out_single, list(selected = selected, average = average))
    }
    if (length(out_pair)) {
      selected <- copy(out_pair[[chosen_candidate]])
      selected[, family := 'Selected from training']
      out_pair <- c(out_pair, list(selected = selected))
    }
    operator_audit <- rbindlist(lapply(names(prepared$operators), function(dr) {
      op <- prepared$operators[[dr]]
      data.table(job = job$tag, drug = dr, kind = op$kind, lambda = op$lambda,
        training_molecular_cells = length(op$training_ids), shared_held_cells = sum(op$training_ids %in% excluded))
    }))
    result <- list(run_contract = run_contract, single = rbindlist(out_single), pair = rbindlist(out_pair),
      selection = selection, operators = operator_audit, selected = chosen_candidate,
      models = if (job$regime == 'final') models else NULL)
    qsave(result, result_path)
    rm(prepared, fitted)
    gc(verbose = FALSE)
  }
  single_predictions[[job$tag]] <- result$single
  pair_predictions[[job$tag]] <- result$pair
  selection_tables[[job$tag]] <- result$selection
  operator_tables[[job$tag]] <- result$operators
  if (job$regime == 'final') final_models <- result$models
  if (length(single_predictions)) write_table(rbindlist(single_predictions), 'held_single_predictions_frozen')
  if (length(pair_predictions)) write_table(rbindlist(pair_predictions), 'held_combination_predictions_frozen')
  write_table(rbindlist(selection_tables), 'training_parameter_selection')
  write_table(rbindlist(operator_tables), 'molecular_operator_audit')
}

status('Applying the displacement-only bridge to HGF; common untreated denominator')
hgf <- fread(file.path(reference, 'data/hgf_transported_inputs.tsv'))
h <- as.matrix(hgf[, P, with = FALSE])
hgf_predictions <- list()
for (candidate in candidates) {
  model <- final_models[[candidate]]
  for (dose in c(1, 2.5)) {
    w <- weights(dose)
    lap <- combined <- matrix(0, nrow(h), length(basis_names))
    ref_lap <- ref_combined <- matrix(0, length(model$context$reference_cells), length(basis_names))
    for (i in 1:3) {
      l <- uppred(PREP, model$operators[[anchor]], AU, native_doses[i])
      l <- l[rep(1L, nrow(h)), , drop = FALSE]
      lap <- lap + response_basis(l, model$context) * w[1, i]
      combined <- combined + response_basis(l + h, model$context) * w[1, i]
      ref_lap <- ref_lap + model$reference_basis[[anchor]][[i]] * w[1, i]
    }
    fl <- fh <- flh <- numeric(nrow(h))
    for (j in seq_len(nrow(h))) {
      rb <- matrix(0, nrow(ref_lap), length(basis_names))
      hb <- response_basis(h[j, , drop = FALSE], model$context)
      for (i in 1:3) {
        repeated_h <- hb[rep(1L, nrow(ref_lap)), , drop = FALSE]
        rb <- rb + sum_basis(model$reference_basis[[anchor]][[i]], repeated_h) * w[1, i]
      }
      x_l <- design(lap[j, , drop = FALSE], 'secondary')
      rx_l <- design(matrix(colMeans(ref_lap), 1L), 'secondary')
      x_h <- design(hb, 'secondary')
      x_lh <- design(combined[j, , drop = FALSE], 'secondary')
      rx_lh <- design(matrix(colMeans(rb), 1L), 'secondary')
      fl[j] <- (x_l - rx_l) %*% model$beta + rx_l %*% model$gamma
      fh[j] <- x_h %*% model$gamma
      flh[j] <- (x_lh - rx_lh) %*% model$beta + rx_lh %*% model$gamma
    }
    hgf_predictions[[paste(candidate, dose)]] <- data.table(family = candidate, dose,
      hgf_context = hgf$context, lapatinib_log2 = fl, hgf_log2 = fh, combination_log2 = flh,
      lapatinib_abundance_pct = 100 * 2^fl, hgf_abundance_pct = 100 * 2^fh,
      combination_abundance_pct = 100 * 2^flh, interaction_log2 = flh - fl - fh,
      absolute_hgf_effect_pp = 100 * (2^flh - 2^fl))
  }
}
write_table(rbindlist(hgf_predictions), 'au565_hgf_predictions_frozen')
prediction_paths <- file.path(root, 'data', c('held_single_predictions_frozen.tsv',
  'held_combination_predictions_frozen.tsv', 'au565_hgf_predictions_frozen.tsv'))
write_table(data.table(path = prediction_paths, sha256 = vapply(prediction_paths,
  digest, character(1), file = TRUE, algo = 'sha256')), 'prediction_freeze_hashes')

# Evaluation starts after model selection and serialized predictions are frozen.
status('Evaluating fixed-treatment, held-cell correlations')
single <- rbindlist(single_predictions)
pair <- rbindlist(pair_predictions)
stopifnot(!anyDuplicated(single[, .(family, single_id)]))
single_metrics <- single[, .(cells = uniqueN(depmap_id), pearson_log2 = pearson(observed_log2, predicted_log2),
  observed_below_50pct = sum(observed_log2 <= -1), predicted_below_50pct = sum(predicted_log2 <= -1)),
  by = .(family, drug, dose)]
within_fold <- single[, .(cells = uniqueN(depmap_id), pearson_log2 = pearson(observed_log2, predicted_log2),
  predicted_sd_log2 = sd(predicted_log2)), by = .(family, fold, drug, dose)]
write_table(single_metrics, 'fixed_drug_dose_single_metrics')
write_table(within_fold, 'fixed_drug_dose_fold_metrics')
combo_context <- pair[, .(cells = uniqueN(depmap_id), pearson_abundance = pearson(observed_log2, predicted_log2),
  pearson_departure = pearson(observed_departure_log2, predicted_departure_log2)),
  by = .(family, regime, SCREENER, pair, doseA, doseB)]
combo_summary <- combo_context[, .(contexts = .N,
  finite_abundance = sum(is.finite(pearson_abundance)),
  median_pearson_abundance = median(pearson_abundance, na.rm = TRUE),
  finite_departure = sum(is.finite(pearson_departure)),
  median_pearson_departure = median(pearson_departure, na.rm = TRUE)), by = .(family, regime)]
pooled_combo <- pair[, .(conditions = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
  pearson_abundance = pearson(observed_log2, predicted_log2),
  pearson_departure = pearson(observed_departure_log2, predicted_departure_log2)), by = .(family, regime)]
write_table(combo_context, 'fixed_pair_dose_across_cell_metrics')
write_table(combo_summary, 'fixed_pair_dose_summary')
write_table(pooled_combo, 'pooled_combination_metrics')
b1 <- single_metrics[drug == anchor & dose == 2.5]
stopifnot(all(b1$cells == 441L))
write_table(b1, 'b1_441_cell_comparison')
comparison78 <- fread(file.path(previous, 'data/fixed_drug_dose_single_metrics.tsv'))[
  family == 'Joint singles and combinations']
comparison78[, family := 'v78 baseline-conditioned']
write_table(rbindlist(list(single_metrics, comparison78), use.names = TRUE, fill = TRUE), 'comparison_with_v78')

p <- ggplot(single[family == 'Selected from training' & dose == 2.5],
  aes(x = observed_log2, y = predicted_log2)) +
  geom_point(size = .6, alpha = .45) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .45, color = '#286b83') +
  facet_wrap(~drug, ncol = 5, scales = 'free') +
  theme_bw(base_size = 9) + theme(aspect.ratio = 1, strip.text = element_text(size = 8)) +
  labs(x = 'Observed log2 abundance / untreated', y = 'Predicted log2 abundance / untreated',
    title = 'Displacement-only bridge: held-out cells at 2.5 µM')
annotations <- single_metrics[family == 'Selected from training' & dose == 2.5,
  .(drug, label = sprintf('r = %.3f; n = %d', pearson_log2, cells))]
p <- p + geom_text(data = annotations, aes(x = -Inf, y = Inf, label = label),
  inherit.aes = FALSE, hjust = -.05, vjust = 1.3, size = 2.7)
ggsave(file.path(root, 'figures/all_drugs_2_5uM_across_cells.pdf'), p, width = 12, height = 5.4)
ggsave(file.path(root, 'figures/all_drugs_2_5uM_across_cells.png'), p, width = 12, height = 5.4, dpi = 150)

selected_context <- combo_context[family == 'Selected from training']
figure_pair <- copy(pair[family == 'Selected from training'])
figure_pair[, comparison := ifelse(regime == 'recipient', 'Held-out cells', 'Held-out partners')]
cp <- ggplot(figure_pair, aes(observed_log2, predicted_log2)) +
  geom_point(size = .5, alpha = .25) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5) +
  facet_wrap(~comparison, scales = 'free') + theme_bw(base_size = 10) + theme(aspect.ratio = 1) +
  labs(x = 'Observed combination log2 abundance', y = 'Predicted combination log2 abundance',
    title = 'Mixed treatment conditions; assess fixed-treatment correlations separately')
cp_annotations <- pooled_combo[family == 'Selected from training',
  .(comparison = ifelse(regime == 'recipient', 'Held-out cells', 'Held-out partners'),
    label = sprintf('Pooled r = %.3f', pearson_abundance))]
cp <- cp + geom_text(data = cp_annotations, aes(x = -Inf, y = Inf, label = label),
  inherit.aes = FALSE, hjust = -.05, vjust = 1.3, size = 3)
ggsave(file.path(root, 'figures/combinations_vs_ground_truth.pdf'), cp, width = 7.5, height = 4)
ggsave(file.path(root, 'figures/combinations_vs_ground_truth.png'), cp, width = 7.5, height = 4, dpi = 150)

input_paths <- c(script_path, file.path(previous, 'models/settings.qs'),
  file.path(previous, 'data/individual_training_outcomes.tsv'),
  file.path(previous, 'data/combination_training_outcomes.tsv'),
  file.path(previous, 'data/recipient_folds.tsv'), file.path(original, 'modular_pair_pilot/inputs.rds'),
  file.path(reference, 'data/hgf_transported_inputs.tsv'))
write_table(data.table(path = input_paths, sha256 = vapply(input_paths, digest,
  character(1), file = TRUE, algo = 'sha256')), 'source_manifest')
qsave(list(run_contract = run_contract, selected_drugs = selected_drugs, native_doses = native_doses,
  molecular_lambdas = molecular_lambdas, molecular_selection = 'Within-native-dose molecular contrasts; no mean candidate',
  functional_inputs = 'Program displacement only; baseline occurs only upstream',
  sensitivity_objective = 'Weighted within-drug/dose or pair/dose contrasts',
  model_selection = 'Inner held-cell mean within-treatment Pearson; nonlinear improvement must exceed .02',
  average_lambda = average_lambda, source_totals = source_totals, feature_names = feature_names,
  retrospective = TRUE, no_functional_HGF_or_Wilson_training = TRUE), file.path(root, 'models/settings.qs'))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
capture.output(system('git -C /work rev-parse HEAD', intern = TRUE), file = file.path(root, 'git_revision.txt'))
capture.output(cat('V79: strict baseline -> molecular displacement -> function.\n',
  'No baseline covariates or gates in the functional regression.\n',
  'Mean-only molecular models excluded; shrinkage chosen on within-dose RNA differences.\n',
  'Sensitivity fitted using training-only within-treatment centering. Average calibration cannot change sensitivity coefficients.\n',
  'Functional model/regularization selection uses inner within-treatment Pearson, never average effects.\n',
  'Nonlinear candidate adds all820 displacement products; no baseline bypass.\n',
  'Whole receiving-cell and whole-partner exclusions; v78 dose alignment and all441 B1 points retained.\n',
  'Average-only pooled correlations can reflect fold offsets; inspect within-fold controls.\n',
  'Retrospective development; HGF/Wilson function not fitted.\n'),
  print(b1), print(single_metrics[dose == 2.5]), print(combo_summary), print(pooled_combo),
  print(rbindlist(hgf_predictions)[hgf_context == 'Equal-context mean']),
  file = file.path(root, 'RESULTS.txt'))
status('Complete')
