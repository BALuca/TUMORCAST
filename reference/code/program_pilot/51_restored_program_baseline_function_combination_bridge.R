#!/usr/bin/env Rscript
# Frozen34 hard programs: baseline -> response -> function, and two-response interactions.
# The43 bank used full-panel outcomes; all current validation is retrospective.
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/functional_programs.R')
source('/work/src/program_pilot/response_models.R')
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(20261011L)
options(warn = 1)

root <- '/work/results/program_pilot/51_restored_program_baseline_function_combination_bridge'
data_dir <- file.path(root, 'data')
model_dir <- file.path(root, 'models')
figure_dir <- file.path(root, 'figures')
for (path in c(root, data_dir, model_dir, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
script_path <- '/work/src/program_pilot/51_restored_program_baseline_function_combination_bridge.R'
bank_path <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
combo_path <- '/work/almanac_broad9_2026-10-07/results/prepared.rds'
control_path <- '/work/almanac_broad9_2026-10-07/structural_diagnosis/results/raw_eligible_records.tsv.gz'
decoder_root <- '/work/results/program_pilot/50_restore_calibrated_svr_current_program_bank'
families <- c('ERBB', 'MEK')
focal_drugs <- c('lapatinib', 'crizotinib', 'alpelisib', 'trametinib', 'gefitinib', 'everolimus', 'bortezomib')
target_genes <- list(lapatinib = c('ERBB2', 'EGFR'), crizotinib = c('MET', 'ALK', 'ROS1'),
  alpelisib = 'PIK3CA', trametinib = c('MAP2K1', 'MAP2K2'), gefitinib = 'EGFR',
  everolimus = 'MTOR', bortezomib = 'PSMB5')
lambda_grid <- c(100, 10, 1, .1, .01)
held_ids <- c('ACH-000138', 'ACH-000178', 'ACH-000493', 'ACH-000552', 'ACH-000669',
  'ACH-000875', 'ACH-000927', 'ACH-000750', 'ACH-000463')
write_table <- function(x, name) fwrite(as.data.table(x), file.path(data_dir, name), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(8); quit(save = 'no', status = 1L)})
status('preparing frozen programs, baseline and original combination records')
paths <- c(bank_path, panel_path, expression_path, combo_path, control_path,
  file.path(decoder_root, 'models', paste0(families, '_current_decoder.qs')))
check(all(file.exists(paths)), 'Missing bridge input')
manifest <- data.table(path = paths,
  sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
write_table(manifest, 'input_hashes.tsv')
signature <- digest::digest(list(manifest, held_ids, target_genes, lambda_grid,
  design = 'unpenalized dose intercepts + baseline + baseline by log10dose;34 hard means;all measured combo partners'),
  algo = 'sha256')
qsave(list(signature = signature, held_ids = held_ids, focal_drugs = focal_drugs,
  targets = target_genes, lambda_grid = lambda_grid, combination_primary_decoder = 'ERBB'),
  file.path(root, 'settings.qs'))
bank <- qread(bank_path)
panel <- qread(panel_path)
combo <- readRDS(combo_path)
md <- as.data.table(bank$metadata)
measured <- bank$scores
programs <- colnames(measured)
rank <- length(programs)
check(rank == 34L && identical(rownames(measured), md$observation_id) && all(is.finite(measured)),
  'Frozen43 scores are invalid')
check(identical(rownames(bank$weights), bank$genes) && identical(colnames(bank$weights), programs),
  'Gene/program identities changed')
check(max(abs(colSums(bank$weights) - 1)) < 1e-12 && all(rowSums(bank$weights > 0) == 1),
  'Expected normalized hard cluster means')
check(!anyDuplicated(md$observation_id) && !anyDuplicated(md[, .(depmap_id, drug, dose)]),
  'Duplicated molecular conditions')
# A frozen bank can score additional perturbations without relearning programs.
#43 saved its105 discovery drugs;20's combination test also uses other partners.
needed_drugs <- unique(c(focal_drugs, panel_name(c(combo$d$drugA, combo$d$drugB))))
extra_drugs <- intersect(setdiff(needed_drugs, md$drug), panel$metadata$drug)
extra_md <- as.data.table(panel$metadata)[drug %chin% extra_drugs]
if (nrow(extra_md)) {
  projection_path <- file.path(data_dir, 'additional_partner_projection.qs')
  if (file.exists(projection_path)) {
    projected <- qread(projection_path)
    check(identical(projected$signature, signature) && identical(projected$metadata$observation_id,
      extra_md$observation_id), 'Stale additional partner projection')
    extra_matrix <- projected$scores
  } else {
  extra_scores <- list()
  for (id in unique(extra_md$depmap_id)) {
    native <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
    ids <- extra_md[depmap_id == id, observation_id]
    check(all(ids %chin% rownames(native$response)), 'Additional partner signatures missing')
    extra_scores[[id]] <- as.matrix(native$response[ids, bank$genes, drop = FALSE] %*% bank$weights)
  }
  extra_matrix <- do.call(rbind, extra_scores)[extra_md$observation_id, , drop = FALSE]
  qsave(list(signature = signature, metadata = extra_md, scores = extra_matrix), projection_path)
  rm(extra_scores, native)
  }
  measured <- rbind(measured, extra_matrix)
  md <- rbind(md, extra_md, fill = TRUE)
  rm(extra_matrix)
  check(identical(rownames(measured), md$observation_id) && all(is.finite(measured)), 'Additional partner scoring failed')
}
write_table(data.table(drug = needed_drugs, original43_scores = needed_drugs %chin% bank$metadata$drug,
  projected = needed_drugs %chin% md$drug), 'combination_partner_coverage.tsv')
decoders <- setNames(lapply(families, function(f) qread(file.path(decoder_root, 'models',
  paste0(f, '_current_decoder.qs')))), families)
decode <- function(family, x) {
  f <- decoders[[family]]
  cf <- unname(coef(f$calibration))
  cf[1L] + cf[2L] * predict_program_svr(f$fitted, x)
}
for (f in families) check(setequal(decoders[[f]]$excluded_cells, held_ids), 'Decoder cell split changed')
write_table(data.table(family = families, zero_response_predicted_pct = vapply(families,
  function(f) 100 * 2^decode(f, matrix(0, 1L, rank)), numeric(1))), 'decoder_zero_input.tsv')

# Project pretreatment expression onto the identical complete gene memberships.
baseline_path <- file.path(data_dir, 'baseline_programs.qs')
if (file.exists(baseline_path)) {
  prepared <- qread(baseline_path)
  check(identical(prepared$signature, signature), 'Stale baseline cache')
  baseline <- prepared$baseline
  targets <- prepared$targets
} else {
  header <- names(fread(expression_path, nrows = 0L))
  symbols <- sub(' \\(ENSG.*', '', header[-1L])
  selected <- which(symbols %chin% unique(c(bank$genes, unlist(target_genes))))
  expression_table <- fread(expression_path, select = c(1L, selected + 1L))
  expression <- as.matrix(expression_table[, -1L])
  rownames(expression) <- expression_table[[1L]]
  colnames(expression) <- symbols[selected]
  check(!anyDuplicated(colnames(expression)) && all(bank$genes %chin% colnames(expression)) &&
    all(unlist(target_genes) %chin% colnames(expression)), 'Missing/duplicated baseline genes')
  ids <- intersect(unique(md$depmap_id), rownames(expression))
  baseline <- as.matrix(expression[ids, bank$genes, drop = FALSE] %*% bank$weights)
  targets <- expression[ids, unique(unlist(target_genes)), drop = FALSE]
  qsave(list(signature = signature, baseline = baseline, targets = targets), baseline_path)
  rm(expression, expression_table)
  invisible(gc())
}
check(all(is.finite(baseline)) && all(held_ids %chin% rownames(baseline)), 'Invalid baseline coverage')
# Check the scoring operation against saved43 native responses without rereading the panel.
first_cell <- md$depmap_id[1L]
native <- qread(panel$input_paths[match(first_cell, panel$cell_map$depmap_id)])
ii <- which(md$depmap_id == first_cell)[1:3]
replayed <- as.matrix(native$response[md$observation_id[ii], bank$genes, drop = FALSE] %*% bank$weights)
check(max(abs(replayed - measured[ii, ])) < 1e-10, 'Hard scores differ from43')
rm(native, replayed)
inner_ids <- sort(setdiff(rownames(baseline), held_ids))
inner_map <- setNames(sample(rep(1:5, length.out = length(inner_ids))), inner_ids)
write_table(data.table(depmap_id = names(inner_map), fold = unname(inner_map)), 'baseline_inner_folds.tsv')
write_table(panel$cell_map[, .(depmap_id, depmap_name, held = depmap_id %chin% held_ids,
  baseline_available = depmap_id %chin% rownames(baseline))], 'cell_split.tsv')
baseline_mu <- colMeans(baseline[inner_ids, , drop = FALSE])
baseline_sd <- apply(baseline[inner_ids, , drop = FALSE], 2L, sd)
check(all(baseline_sd > 0), 'Constant baseline program')
baseline_z <- sweep(sweep(baseline, 2L, baseline_mu, '-'), 2L, baseline_sd, '/')

# Dose rows remain independent. Shared baseline slopes and baseline-by-logdose
# slopes borrow information across doses; dose intercepts preserve average effects.
baseline_design <- function(drug_name, rows) {
  ids <- md$depmap_id[rows]
  dose_levels <- sort(unique(md[drug == drug_name, dose]))
  dd <- outer(md$dose[rows], dose_levels, '==') * 1
  colnames(dd) <- paste0('dose_', dose_levels)
  dd <- dd[, -1L, drop = FALSE] # The intercept represents the first dose.
  bb <- baseline[ids, , drop = FALSE]
  tt <- if (drug_name %in% names(target_genes)) targets[ids, target_genes[[drug_name]], drop = FALSE] else NULL
  bx <- cbind(bb, tt)
  slope <- bx * (log10(md$dose[rows]) - mean(log10(dose_levels)))
  colnames(slope) <- paste0(colnames(bx), '_logdose')
  cbind(dd, bx, slope)
}
# These two fits reuse20's target scaling/prediction contract and glmnet ridge.
# The only extension is unpenalized dose contrasts, shared by both baseline heads.
fit_dose_path <- function(x, y) {
  multi <- is.matrix(y)
  center <- if (multi) colMeans(y) else 0
  scale <- if (multi) apply(y, 2L, sd) else 1
  check(all(is.finite(x)) && all(is.finite(y)) && all(scale > 1e-12), 'Invalid baseline fit inputs')
  target <- if (multi) sweep(sweep(y, 2L, center, '-'), 2L, scale, '/') else y
  penalty <- as.numeric(!startsWith(colnames(x), 'dose_'))
  fit <- glmnet(x, target, family = if (multi) 'mgaussian' else 'gaussian', alpha = 0,
    lambda = lambda_grid, standardize = TRUE, intercept = TRUE, penalty.factor = penalty, thresh = 1e-7)
  list(fit = fit, center = center, scale = scale, multi = multi)
}
fit_dose_ridge <- function(x, y, ids) {
  fold <- unname(inner_map[ids])
  check(!anyNA(fold) && length(unique(fold)) >= 3L, 'Invalid baseline inner folds')
  loss <- numeric(length(lambda_grid))
  for (k in sort(unique(fold))) {
    tr <- which(fold != k)
    va <- which(fold == k)
    fit <- fit_dose_path(x[tr, , drop = FALSE], if (is.matrix(y)) y[tr, , drop = FALSE] else y[tr])
    p <- pred_path(fit, x[va, , drop = FALSE])
    if (is.matrix(y)) {
      for (j in seq_along(lambda_grid)) {
        loss[j] <- loss[j] + sum(rowMeans(sweep(p[, , j] - y[va, , drop = FALSE], 2L, fit$scale, '/')^2))
      }
    } else loss <- loss + colSums((p - y[va])^2)
  }
  fit <- fit_dose_path(x, y)
  fit$lambda <- lambda_grid[which.min(loss)]
  fit$selection_loss <- loss
  fit$training_ids <- unique(ids)
  fit
}
cd <- as.data.table(combo$d)
check(!anyDuplicated(cd$source_row) && max(abs(cd$interaction - (cd$y - cd$yA - cd$yB))) < 1e-10,
  'Combination endpoint identities failed')
combination_drugs <- unique(panel_name(c(cd$drugA, cd$drugB)))
fit_drugs <- unique(c(focal_drugs, intersect(combination_drugs, md$drug)))
check(all(focal_drugs %chin% md$drug), 'Missing focal drug')
predicted_response <- matrix(NA_real_, nrow(md), rank, dimnames = dimnames(measured))
mean_response <- predicted_response
direct_function <- rep(NA_real_, nrow(md))
labels <- as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
  observation_id %chin% md$observation_id & depmap_id %chin% rownames(baseline)]
labels[, `:=`(measured_pct = prism_abundance_pct, observed_log2 = log2(prism_abundance_pct / 100),
  split = fifelse(depmap_id %chin% held_ids, 'Held cells', 'Training cells'))]
support <- list()
for (drug_name in fit_drugs) {
  rows <- which(md$drug == drug_name & md$depmap_id %chin% rownames(baseline))
  ids <- md$depmap_id[rows]
  tr <- which(!ids %chin% held_ids)
  te <- which(ids %chin% held_ids)
  check(uniqueN(ids[tr]) >= 8L && length(te) > 0L, paste('Insufficient molecular support:', drug_name))
  x <- baseline_design(drug_name, rows)
  y <- measured[rows, , drop = FALSE]
  path <- file.path(model_dir, paste0('baseline_', gsub('[^a-z0-9]+', '_', drug_name), '.qs'))
  if (file.exists(path)) {
    saved <- qread(path)
    check(identical(saved$signature, signature), 'Stale baseline model')
    f <- saved$model
  } else {
    status(paste('baseline -> displacement:', drug_name))
    f <- fit_dose_ridge(x[tr, , drop = FALSE], y[tr, , drop = FALSE], ids[tr])
    qsave(list(signature = signature, model = f, drug = drug_name,
      training_observations = md$observation_id[rows[tr]], feature_names = colnames(x)), path)
  }
  check(!any(f$training_ids %chin% held_ids), 'Held cell in baseline fit')
  predicted_response[rows, ] <- predict_ridge(f, x)
  for (dose_value in sort(unique(md$dose[rows]))) {
    dtr <- tr[md$dose[rows[tr]] == dose_value]
    dest <- rows[md$dose[rows] == dose_value]
    mean_response[dest, ] <- matrix(colMeans(y[dtr, , drop = FALSE]), length(dest), rank, byrow = TRUE)
  }
  # A fair direct baseline->function control uses the identical pooled-dose design.
  direct_lambda <- NA_real_
  if (drug_name %chin% focal_drugs) {
  ld <- labels[drug == drug_name]
  ltr <- which(!ld$depmap_id %chin% held_ids)
  check(length(ltr) >= 8L, paste('Insufficient direct-function support:', drug_name))
  lx <- x[match(ld$observation_id, md$observation_id[rows]), , drop = FALSE]
  direct_path <- file.path(model_dir, sub('baseline_', 'direct_', basename(path), fixed = TRUE))
  if (file.exists(direct_path)) {
    ds <- qread(direct_path)
    check(identical(ds$signature, signature), 'Stale direct model')
    df <- ds$model
  } else {
    df <- fit_dose_ridge(lx[ltr, , drop = FALSE], ld$observed_log2[ltr], ld$depmap_id[ltr])
    qsave(list(signature = signature, model = df, drug = drug_name,
      training_observations = ld$observation_id[ltr]), direct_path)
  }
  direct_function[match(ld$observation_id, md$observation_id)] <- predict_ridge(df, lx)
  direct_lambda <- df$lambda
  }
  support[[drug_name]] <- data.table(drug = drug_name, focal = drug_name %chin% focal_drugs,
    molecular_train_cells = uniqueN(ids[tr]), molecular_held_cells = uniqueN(ids[te]),
    native_rows = length(rows), baseline_lambda = f$lambda, direct_lambda = direct_lambda)
  write_table(rbindlist(support), 'baseline_training_support.tsv')
  qsave(list(signature = signature, prediction = predicted_response[rows, ],
    mean_response = mean_response[rows, ], observation_id = md$observation_id[rows]),
    file.path(data_dir, paste0('response_', gsub('[^a-z0-9]+', '_', drug_name), '.qs')))
  status(paste('saved baseline response/function models:', drug_name))
}
qsave(list(signature = signature, metadata = md, actual = measured, predicted = predicted_response,
  reference = mean_response, baseline = baseline, targets = targets), file.path(data_dir, 'bridge_responses.qs'))

# Molecular diagnostics remove the TRAINING treatment mean from both vectors.
# This isolates recipient variation rather than crediting the average dose shape.
molecular_predictions <- rbindlist(lapply(focal_drugs, function(drug_name) {
  rows <- which(md$drug == drug_name & md$depmap_id %chin% rownames(baseline))
  rbindlist(lapply(seq_along(rows), function(k) {
    i <- rows[k]
    data.table(observation_id = md$observation_id[i], depmap_id = md$depmap_id[i], drug = drug_name,
      dose = md$dose[i], split = fifelse(md$depmap_id[i] %chin% held_ids, 'Held cells', 'Training cells'),
      program = programs, observed = measured[i, ], predicted = predicted_response[i, ], reference = mean_response[i, ])
  }))
}))
write_table(molecular_predictions, 'molecular_predictions.tsv.gz')
molecular_metrics <- molecular_predictions[, {
  context <- paste(program, dose)
  observed_contrast <- observed - ave(observed, context, FUN = mean)
  predicted_contrast <- predicted - ave(predicted, context, FUN = mean)
  .(cells = uniqueN(depmap_id), profile_r = cor_safe(observed, predicted),
    training_mean_residual_r = cor_safe(observed - reference, predicted - reference),
    across_cell_program_treatment_r = cor_safe(observed_contrast, predicted_contrast))
}, by = .(split, drug)]
write_table(molecular_metrics, 'molecular_by_drug.tsv')
write_table(molecular_predictions[, .(cells = uniqueN(depmap_id), r = cor_safe(observed, predicted)),
  by = .(split, drug, dose, program)], 'molecular_within_treatment_program.tsv')
variation <- molecular_predictions[, .(observed_sd = sd(observed), predicted_sd = sd(predicted)),
  by = .(split, drug, dose, program)]
variation[, predicted_to_observed_sd := predicted_sd / observed_sd]
write_table(variation, 'molecular_prediction_variation.tsv')
write_table(variation[, .(median_predicted_to_observed_sd = median(predicted_to_observed_sd)),
  by = .(split, drug)], 'molecular_variation_by_drug.tsv')

# Freeze functional predictions for BOTH pre-existing family decoders. Their
# single-drug family exclusions remain intact. No response->function refitting.
lp <- labels[drug %chin% focal_drugs]
li <- match(lp$observation_id, md$observation_id)
check(!anyNA(li) && all(is.finite(predicted_response[li, ])), 'Missing full-chain responses')
function_predictions <- list()
for (family in families) {
  for (mode in c('Measured response', 'Predicted response', 'Mean response')) {
    xx <- switch(mode, 'Measured response' = measured[li, , drop = FALSE],
      'Predicted response' = predicted_response[li, , drop = FALSE],
      'Mean response' = mean_response[li, , drop = FALSE])
    p <- copy(lp)
    p[, `:=`(decoder = family, model = mode, predicted_log2 = decode(family, xx),
      drug_scope = fifelse(drug %chin% decoders[[family]]$excluded_drugs, 'Excluded family', 'Familiar drug'))]
    function_predictions[[length(function_predictions) + 1L]] <- p
  }
  p <- copy(lp)
  p[, `:=`(decoder = family, model = 'Direct baseline', predicted_log2 = direct_function[li],
    drug_scope = fifelse(drug %chin% decoders[[family]]$excluded_drugs, 'Excluded family', 'Familiar drug'))]
  function_predictions[[length(function_predictions) + 1L]] <- p
  means <- lp[split == 'Training cells', .(reference_y = mean(observed_log2)), by = .(drug, dose)]
  p <- merge(copy(lp), means, by = c('drug', 'dose'), all.x = TRUE, sort = FALSE)
  p[, `:=`(decoder = family, model = 'Treatment mean', predicted_log2 = reference_y,
    drug_scope = fifelse(drug %chin% decoders[[family]]$excluded_drugs, 'Excluded family', 'Familiar drug'))]
  function_predictions[[length(function_predictions) + 1L]] <- p
}
fp <- rbindlist(function_predictions, fill = TRUE)
fp[, predicted_pct := 100 * 2^predicted_log2]
check(all(is.finite(fp$predicted_pct)), 'Nonfinite functional prediction')
write_table(fp[, .(decoder, model, observation_id, predicted_log2, predicted_pct)], 'function_predictions_frozen.tsv.gz')
write_table(data.table(sha256 = digest::digest(file = file.path(data_dir, 'function_predictions_frozen.tsv.gz'),
  algo = 'sha256')), 'function_prediction_hash.tsv')
write_table(fp, 'function_predictions.tsv.gz')

# The grouping axis is stated in every table: individual drugs, fixed treatments,
# and different drugs at one fixed dose within each recipient.
function_summary <- function(z) {
  sensitive <- z$measured_pct <= 50
  fixed <- paste(z$drug, z$dose)
  oy <- z$observed_log2 - ave(z$observed_log2, fixed, FUN = mean)
  py <- z$predicted_log2 - ave(z$predicted_log2, fixed, FUN = mean)
  cross <- z[, .(r = cor_safe(observed_log2, predicted_log2)), by = .(depmap_id, dose)]
  data.table(rows = nrow(z), cells = uniqueN(z$depmap_id), drugs = uniqueN(z$drug),
    pooled_r = cor_safe(z$observed_log2, z$predicted_log2),
    across_cells_fixed_treatment_r = if (sd(py) > 1e-10) cor_safe(oy, py) else NA_real_,
    median_across_drugs_same_cell_dose_r = median(cross$r, na.rm = TRUE),
    sensitive = sum(sensitive), recovered = sum(sensitive & z$predicted_pct <= 50),
    resistant = sum(!sensitive), false_sensitive = sum(!sensitive & z$predicted_pct <= 50),
    measured_sensitive_median_pct = median(z$measured_pct[sensitive]),
    predicted_sensitive_median_pct = median(z$predicted_pct[sensitive]))
}
fm <- fp[, function_summary(.SD), by = .(decoder, model, split)]
write_table(fm, 'function_summary.tsv')
fd <- fp[, function_summary(.SD), by = .(decoder, model, split, drug_scope)]
write_table(fd, 'function_by_training_scope.tsv')
drug_metrics <- fp[, .(rows = .N, cells = uniqueN(depmap_id), pooled_r = cor_safe(observed_log2, predicted_log2),
  sensitive = sum(measured_pct <= 50), recovered = sum(measured_pct <= 50 & predicted_pct <= 50),
  resistant = sum(measured_pct > 50), false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50)),
  by = .(decoder, model, split, drug, drug_scope)]
write_table(drug_metrics, 'function_by_drug.tsv')
write_table(fp[, .(rows = .N, r = cor_safe(observed_log2, predicted_log2)),
  by = .(decoder, model, split, drug, dose)], 'function_within_treatment.tsv')
write_table(fp[, .(drugs = uniqueN(drug), r = cor_safe(observed_log2, predicted_log2)),
  by = .(decoder, model, split, depmap_id, dose)], 'function_across_drugs_same_cell_dose.tsv')
write_table(fp[, .(measured_cell_range_pp = max(measured_pct) - min(measured_pct),
  predicted_cell_range_pp = max(predicted_pct) - min(predicted_pct)),
  by = .(decoder, model, split, drug, dose)], 'function_cell_sensitivity_range.tsv')
status('baseline -> response -> function finished; preparing combination transfer')

# Native scores and baseline predictions share the SAME logdose interpolation.
# Exact doses are used directly; outside-range records are explicitly excluded.
query <- unique(rbind(cd[, .(depmap_id, drug = drugA, dose = doseA)],
  cd[, .(depmap_id, drug = drugB, dose = doseB)]))
query[, `:=`(molecular_drug = panel_name(drug), matching = 'no measured cell/drug',
  lower_dose = NA_real_, upper_dose = NA_real_, alpha = NA_real_)]
actual_query <- predicted_query <- matrix(NA_real_, nrow(query), rank, dimnames = list(NULL, programs))
for (i in seq_len(nrow(query))) {
  rows <- which(md$depmap_id == query$depmap_id[i] & md$drug == query$molecular_drug[i])
  if (!length(rows)) next
  rows <- rows[order(md$dose[rows])]
  doses <- md$dose[rows]
  target <- query$dose[i]
  exact <- which(abs(log(target / doses)) < 1e-8)
  if (length(exact)) {
    lo <- hi <- exact[1L]
    alpha <- 0
    matching <- 'exact'
  } else {
    query[i, matching := 'outside measured dose range']
    if (target < min(doses) || target > max(doses)) next
    lo <- max(which(doses < target))
    hi <- min(which(doses > target))
    alpha <- log(target / doses[lo]) / log(doses[hi] / doses[lo])
    matching <- 'interpolated'
  }
  actual_query[i, ] <- (1 - alpha) * measured[rows[lo], ] + alpha * measured[rows[hi], ]
  predicted_query[i, ] <- (1 - alpha) * predicted_response[rows[lo], ] + alpha * predicted_response[rows[hi], ]
  set(query, i, c('lower_dose', 'upper_dose', 'alpha', 'matching'), list(doses[lo], doses[hi], alpha, matching))
}
write_table(query, 'combination_dose_matching.tsv')
key <- function(cell, drug, dose) paste(cell, drug, dose, sep = '__')
qkey <- key(query$depmap_id, query$drug, query$dose)
qa <- match(key(cd$depmap_id, cd$drugA, cd$doseA), qkey)
qb <- match(key(cd$depmap_id, cd$drugB, cd$doseB), qkey)
check(!anyNA(qa) && !anyNA(qb), 'Missing combination query')
keep <- is.finite(actual_query[qa, 1L]) & is.finite(actual_query[qb, 1L]) &
  cd$depmap_id %chin% rownames(baseline)
write_table(cd[!keep], 'combination_excluded_rows.tsv.gz')
cd <- cd[keep]
qa <- qa[keep]
qb <- qb[keep]
a <- actual_query[qa, , drop = FALSE]
bb <- actual_query[qb, , drop = FALSE]
ap <- predicted_query[qa, , drop = FALSE]
bp <- predicted_query[qb, , drop = FALSE]
check(all(is.finite(ap)) && all(is.finite(bp)), 'Missing baseline-only partner responses')
legacy <- new.env(parent = globalenv())
sys.source('/work/almanac_broad9_2026-10-07/model.R', envir = legacy)
legacy$P <- programs
products <- legacy$bilinear(a, bb)
predicted_products <- legacy$bilinear(ap, bp)
check(max(abs(products - legacy$bilinear(bb, a))) < 1e-10 &&
  max(abs(legacy$bilinear(a * 0, bb))) < 1e-12, 'Interaction symmetry/zero boundary failed')
# A treatment-average control is learned from development singles only. Dose
# interpolation has the same logdose rule; out-of-range contexts stay unsupported.
predict_treatment_average <- function(means, drug, dose) {
  answer <- rep(NA_real_, length(drug))
  for (i in seq_along(drug)) {
    d <- means[drug_name == drug[i]][order(dose_value)]
    if (!nrow(d)) next
    exact <- which(abs(log(dose[i] / d$dose_value)) < 1e-8)
    if (length(exact)) answer[i] <- d$observed[exact[1L]]
    else if (nrow(d) >= 2L && dose[i] >= min(d$dose_value) && dose[i] <= max(d$dose_value)) {
      answer[i] <- approx(log(d$dose_value), d$observed, xout = log(dose[i]), rule = 1L)$y
    }
  }
  answer
}
identity_levels <- sort(unique(c(cd$idA, cd$idB)))
identity <- matrix(0, nrow(cd), length(identity_levels), dimnames = list(NULL, identity_levels))
identity[cbind(seq_len(nrow(cd)), match(cd$idA, identity_levels))] <- 1
identity[cbind(seq_len(nrow(cd)), match(cd$idB, identity_levels))] <- 1
potency <- cbind(cd$yA + cd$yB, abs(cd$yA - cd$yB), cd$yA * cd$yB)
design <- list(molecular_products = products,
  conditioned_products = cbind(products, a * bb * baseline_z[cd$depmap_id, , drop = FALSE]),
  identity_potency = cbind(identity, potency))
predicted_design <- list(molecular_products = predicted_products,
  conditioned_products = cbind(predicted_products, ap * bp * baseline_z[cd$depmap_id, , drop = FALSE]))
raw <- fread(control_path, select = c('source_row', 'single_groupA', 'single_groupB'))
refs <- unique(rbind(raw[, .(source_row, ref = single_groupA)], raw[, .(source_row, ref = single_groupB)]))
refs <- refs[source_row %in% cd$source_row]
refs[, row := match(source_row, cd$source_row)]
check(all(cd$source_row %in% refs$source_row) && !anyNA(refs$ref), 'Missing ALMANAC control references')
purge <- function(train, test) setdiff(train, refs[ref %chin% refs[row %in% test, ref], row])
base_train <- which(!cd$depmap_id %chin% held_ids)
base_test <- which(cd$depmap_id %chin% held_ids)
splits <- list(held_cells = list(train = purge(base_train, base_test), test = base_test))
for (k in 1:5) {
  test <- base_test[cd$outer_fold[base_test] == k]
  train <- base_train[cd$outer_fold[base_train] != k]
  splits[[paste0('held_cells_pairs_', k)]] <- list(train = purge(train, test), test = test)
}
# New combination heads learn ALMANAC log2 interaction residuals. Their inputs
# are two responses in one frozen molecular coordinate system, never pair IDs.
dp_list <- audits <- list()
for (split_name in names(splits)) {
  tr <- splits[[split_name]]$train
  te <- splits[[split_name]]$test
  check(length(tr) > 100L && length(te) > 0L &&
    !length(intersect(cd$depmap_id[tr], cd$depmap_id[te])) &&
    !length(intersect(refs[row %in% tr, ref], refs[row %in% te, ref])), 'Combination split/control leakage')
  if (split_name != 'held_cells') check(!length(intersect(cd$pair[tr], cd$pair[te])), 'Held pair leakage')
  # Assay calibration uses ONLY development singles. Distinct screen outcomes
  # remain separate; duplicates of the same cell/drug/dose/outcome are removed.
  sa <- data.table(depmap_id = cd$depmap_id[tr], drug = cd$idA[tr], dose = cd$doseA[tr],
    observed = cd$yA[tr], query_index = qa[tr], decoded = decode('ERBB', a[tr, , drop = FALSE]))
  sb <- data.table(depmap_id = cd$depmap_id[tr], drug = cd$idB[tr], dose = cd$doseB[tr],
    observed = cd$yB[tr], query_index = qb[tr], decoded = decode('ERBB', bb[tr, , drop = FALSE]))
  singles <- unique(rbind(sa, sb), by = c('depmap_id', 'drug', 'dose', 'observed'))
  calibration <- lm(observed ~ decoded, data = singles, weights = observation_weights(singles))
  cf <- unname(coef(calibration))
  check(length(cf) == 2L && all(is.finite(cf)), 'Invalid single-assay calibration')
  molecular_additive <- 2 * cf[1] + cf[2] * (decode('ERBB', a[te, , drop = FALSE]) +
    decode('ERBB', bb[te, , drop = FALSE]))
  baseline_additive <- 2 * cf[1] + cf[2] * (decode('ERBB', ap[te, , drop = FALSE]) +
    decode('ERBB', bp[te, , drop = FALSE]))
  qsave(list(signature = signature, calibration = calibration, training_singles = singles,
    training_ids = unique(cd$depmap_id[tr]), training_pairs = unique(cd$pair[tr])),
    file.path(model_dir, paste0('single_assay_calibration_', split_name, '.qs')))
  # Learn the ALMANAC single endpoint on matched development controls, rather
  # than asking a PRISM-trained function decoder to transfer across assays.
  # Reuse50's exact C10, gamma1/p and epsilon0.1; no hyperparameter search.
  sx <- actual_query[singles$query_index, , drop = FALSE]
  single_signature <- digest::digest(list(signature, singles, cost = 10, gamma = 1 / rank,
    calibration = 'leave-one-development-cell-out equal-drug affine'), algo = 'sha256')
  single_path <- file.path(model_dir, paste0('almanac_single_decoder_', split_name, '.qs'))
  if (file.exists(single_path)) {
    single_decoder <- qread(single_path)
    check(identical(single_decoder$signature, single_signature), 'Stale ALMANAC single decoder')
  } else {
    status(paste('ALMANAC matched single-function head:', split_name))
    oof <- rep(NA_real_, nrow(singles))
    for (cell in sort(unique(singles$depmap_id))) {
      it <- which(singles$depmap_id != cell)
      va <- which(singles$depmap_id == cell)
      fold_path <- file.path(model_dir, paste0('almanac_single_', split_name, '_', cell, '.qs'))
      if (file.exists(fold_path)) {
        fs <- qread(fold_path)
        check(identical(fs$signature, single_signature), 'Stale ALMANAC single fold')
      } else {
        fitted <- fit_program_svr(sx[it, , drop = FALSE], singles$observed[it], 'radial', 10, 1 / rank)
        fs <- list(signature = single_signature, fitted = fitted, training_ids = unique(singles$depmap_id[it]),
          held_cell = cell, validation_prediction = predict_program_svr(fitted, sx[va, , drop = FALSE]))
        qsave(fs, fold_path)
      }
      check(!cell %in% fs$training_ids, 'ALMANAC inner calibration cell leakage')
      oof[va] <- fs$validation_prediction
    }
    check(all(is.finite(oof)), 'Incomplete ALMANAC calibration predictions')
    sc <- lm(observed ~ oof, data = singles, weights = observation_weights(singles))
    sf <- fit_program_svr(sx, singles$observed, 'radial', 10, 1 / rank)
    single_decoder <- list(signature = single_signature, fitted = sf, calibration = sc,
      training_ids = unique(singles$depmap_id), training_pairs = unique(cd$pair[tr]),
      training_singles = singles, oof = oof, cost = 10, gamma = 1 / rank)
    qsave(single_decoder, single_path)
  }
  check(!any(single_decoder$training_ids %chin% held_ids), 'Held cell used in ALMANAC single head')
  single_coefficients <- unname(coef(single_decoder$calibration))
  check(length(single_coefficients) == 2L && all(is.finite(single_coefficients)), 'Invalid ALMANAC head calibration')
  predict_single <- function(x) single_coefficients[1] + single_coefficients[2] *
    predict_program_svr(single_decoder$fitted, x)
  matched_additive <- predict_single(a[te, , drop = FALSE]) + predict_single(bb[te, , drop = FALSE])
  matched_baseline_additive <- predict_single(ap[te, , drop = FALSE]) + predict_single(bp[te, , drop = FALSE])
  single_test <- rbind(
    data.table(depmap_id = cd$depmap_id[te], drug = cd$idA[te], dose = cd$doseA[te], observed_log2 = cd$yA[te],
      prism_transferred = cf[1] + cf[2] * decode('ERBB', a[te, , drop = FALSE]),
      almanac_matched = predict_single(a[te, , drop = FALSE]),
      baseline_almanac = predict_single(ap[te, , drop = FALSE])),
    data.table(depmap_id = cd$depmap_id[te], drug = cd$idB[te], dose = cd$doseB[te], observed_log2 = cd$yB[te],
      prism_transferred = cf[1] + cf[2] * decode('ERBB', bb[te, , drop = FALSE]),
      almanac_matched = predict_single(bb[te, , drop = FALSE]),
      baseline_almanac = predict_single(bp[te, , drop = FALSE])))
  single_test <- unique(single_test, by = c('depmap_id', 'drug', 'dose', 'observed_log2'))
  write_table(single_test, paste0('almanac_single_predictions_', split_name, '.tsv.gz'))
  make_prediction <- function(name, interaction, additive, input) {
    p <- cd[te, .(source_row, depmap_id, pair, CELLNAME, drugA, drugB, doseA, doseB,
      observed_interaction = interaction, observed_log2 = y, single_a_log2 = yA, single_b_log2 = yB)]
    p[, `:=`(split = split_name, model = name, input = input,
      predicted_interaction = as.numeric(interaction), predicted_log2 = as.numeric(additive + interaction))]
    p
  }
  for (name in names(design)) {
    xx <- design[[name]]
    path <- file.path(model_dir, paste0('combination_', split_name, '_', name, '.qs'))
    if (file.exists(path)) {
      saved <- qread(path)
      check(identical(saved$signature, signature), 'Stale combination model')
      fit <- saved$model
    } else {
      status(paste('combination:', split_name, name))
      loss <- numeric(length(lambda_grid))
      w <- observation_weights(cd[tr], 'pair')
      for (cell in sort(unique(cd$depmap_id[tr]))) {
        va <- tr[cd$depmap_id[tr] == cell]
        it <- purge(setdiff(tr, va), va)
        inner <- fit_path(xx[it, , drop = FALSE], cd$interaction[it], observation_weights(cd[it], 'pair'),
          intercept = name == 'identity_potency')
        pp <- pred_path(inner, xx[va, , drop = FALSE])
        loss <- loss + colSums((pp - cd$interaction[va])^2 * w[match(va, tr)])
      }
      fit <- fit_path(xx[tr, , drop = FALSE], cd$interaction[tr], w, intercept = name == 'identity_potency')
      fit$lambda <- lambda_grid[which.min(loss)]
      fit$training_ids <- unique(cd$depmap_id[tr])
      fit$training_pairs <- unique(cd$pair[tr])
      fit$selection_loss <- loss
      qsave(list(signature = signature, model = fit, training_source_rows = cd$source_row[tr]), path)
    }
    check(!any(fit$training_ids %chin% held_ids), 'Held combination cell used in fitting')
    ip <- predict_ridge(fit, xx[te, , drop = FALSE])
    dp_list[[length(dp_list) + 1L]] <- make_prediction(name, ip, cd$yA[te] + cd$yB[te], 'Measured single outcomes')
    if (name != 'identity_potency') {
      check(max(abs(predict_ridge(fit, xx[te[1:2], , drop = FALSE] * 0))) < 1e-12,
        'Nonzero interaction at zero partner')
      ipp <- predict_ridge(fit, predicted_design[[name]][te, , drop = FALSE])
      dp_list[[length(dp_list) + 1L]] <- make_prediction(name, ip, molecular_additive, 'Measured molecular responses')
      dp_list[[length(dp_list) + 1L]] <- make_prediction(name, ipp, baseline_additive, 'Baseline only')
      dp_list[[length(dp_list) + 1L]] <- make_prediction(name, ip, matched_additive,
        'Measured responses | ALMANAC head')
      dp_list[[length(dp_list) + 1L]] <- make_prediction(name, ipp, matched_baseline_additive,
        'Baseline only | ALMANAC head')
    }
    write_table(rbindlist(dp_list), 'combination_partial_predictions.tsv.gz')
    status(paste('saved combination model:', split_name, name))
  }
  dp_list[[length(dp_list) + 1L]] <- make_prediction('independent_effects', rep(0, length(te)),
    cd$yA[te] + cd$yB[te], 'Measured single outcomes')
  dp_list[[length(dp_list) + 1L]] <- make_prediction('independent_effects', rep(0, length(te)),
    molecular_additive, 'Measured molecular responses')
  dp_list[[length(dp_list) + 1L]] <- make_prediction('independent_effects', rep(0, length(te)),
    baseline_additive, 'Baseline only')
  dp_list[[length(dp_list) + 1L]] <- make_prediction('independent_effects', rep(0, length(te)),
    matched_additive, 'Measured responses | ALMANAC head')
  dp_list[[length(dp_list) + 1L]] <- make_prediction('independent_effects', rep(0, length(te)),
    matched_baseline_additive, 'Baseline only | ALMANAC head')
  cell_means <- singles[, .(observed = mean(observed)), by = .(drug, dose, depmap_id)]
  treatment_means <- cell_means[, .(observed = mean(observed)), by = .(drug_name = drug, dose_value = dose)]
  treatment_additive <- predict_treatment_average(treatment_means, cd$idA[te], cd$doseA[te]) +
    predict_treatment_average(treatment_means, cd$idB[te], cd$doseB[te])
  p <- make_prediction('independent_effects', rep(0, length(te)), treatment_additive, 'Treatment-average singles')
  write_table(p[!is.finite(predicted_log2)], paste0('treatment_average_unsupported_', split_name, '.tsv.gz'))
  dp_list[[length(dp_list) + 1L]] <- p[is.finite(predicted_log2)]
  audits[[split_name]] <- data.table(split = split_name, train_rows = length(tr), test_rows = length(te),
    train_cells = uniqueN(cd$depmap_id[tr]), test_cells = uniqueN(cd$depmap_id[te]),
    train_pairs = uniqueN(cd$pair[tr]), test_pairs = uniqueN(cd$pair[te]),
    overlapping_pairs = length(intersect(cd$pair[tr], cd$pair[te])), shared_controls = 0L,
    assay_calibration_intercept = cf[1], assay_calibration_slope = cf[2])
  write_table(rbindlist(audits), 'combination_split_audit.tsv')
}
dp <- rbindlist(dp_list)
dp[, `:=`(scope = fifelse(split == 'held_cells', 'Held cells', 'Held cells + pairs'),
  measured_pct = 100 * 2^observed_log2, predicted_pct = 100 * 2^predicted_log2,
  observed_bliss_pp = 100 * (2^(single_a_log2 + single_b_log2) - 2^observed_log2),
  conditional_predicted_bliss_pp = 100 * 2^(single_a_log2 + single_b_log2) * (1 - 2^predicted_interaction),
  predicted_bliss_pp = 100 * (2^(predicted_log2 - predicted_interaction) - 2^predicted_log2))]
write_table(dp[, .(source_row, split, model, input, predicted_interaction, predicted_log2)],
  'combination_predictions_frozen.tsv.gz')
write_table(data.table(sha256 = digest::digest(file = file.path(data_dir, 'combination_predictions_frozen.tsv.gz'),
  algo = 'sha256')), 'combination_prediction_hash.tsv')
write_table(dp, 'combination_predictions.tsv.gz')
cm <- dp[, .(rows = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
  endpoint_r = cor_safe(observed_log2, predicted_log2), interaction_r = cor_safe(observed_interaction, predicted_interaction),
  bliss_pp_r = cor_safe(observed_bliss_pp, predicted_bliss_pp),
  sensitive = sum(measured_pct <= 50), recovered = sum(measured_pct <= 50 & predicted_pct <= 50),
  resistant = sum(measured_pct > 50), false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50),
  median_absolute_error_pp = median(abs(measured_pct - predicted_pct))), by = .(scope, input, model)]
write_table(cm, 'combination_summary.tsv')
common_comparison <- dp[scope == 'Held cells + pairs' & input == 'Treatment-average singles', source_row]
write_table(dp[scope == 'Held cells + pairs' & source_row %in% common_comparison,
  .(rows = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
    endpoint_r = cor_safe(observed_log2, predicted_log2), sensitive = sum(measured_pct <= 50),
    recovered = sum(measured_pct <= 50 & predicted_pct <= 50), false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50)),
  by = .(input, model)], 'combination_common_treatment_mean_comparison.tsv')
write_table(dp[, .(rows = .N, r = cor_safe(observed_log2, predicted_log2),
  interaction_r = cor_safe(observed_interaction, predicted_interaction)),
  by = .(scope, input, model, pair, depmap_id)], 'combination_by_pair_cell.tsv')
single_results <- rbindlist(lapply(names(splits), function(s) {
  p <- fread(file.path(data_dir, paste0('almanac_single_predictions_', s, '.tsv.gz')))
  p[, split := s]
  p
}))
single_long <- melt(single_results, id.vars = c('split', 'depmap_id', 'drug', 'dose', 'observed_log2'),
  variable.name = 'model', value.name = 'predicted_log2')
single_long[, `:=`(measured_pct = 100 * 2^observed_log2, predicted_pct = 100 * 2^predicted_log2,
  scope = fifelse(split == 'held_cells', 'Held cells', 'Held cells + pairs'))]
single_metrics <- single_long[, function_summary(.SD), by = .(split, model)]
write_table(single_metrics, 'almanac_single_function_summary.tsv')
write_table(single_long, 'almanac_single_function_predictions.tsv.gz')
single_training <- rbindlist(lapply(names(splits), function(s) {
  f <- qread(file.path(model_dir, paste0('almanac_single_decoder_', s, '.qs')))
  p <- copy(f$training_singles)
  p[, `:=`(split = s, actual_training_prediction = predict_program_svr(f$fitted,
    actual_query[query_index, , drop = FALSE]), raw_cell_oof = f$oof)]
  p
}))
write_table(single_training[, .(controls = .N, cells = uniqueN(depmap_id),
  own_fit_r = cor_safe(observed, actual_training_prediction), raw_whole_cell_oof_r = cor_safe(observed, raw_cell_oof)),
  by = split], 'almanac_single_training_validation.tsv')
# Define informative interaction settings from measured SINGLE effects only.
# These are descriptive evaluations, never used to tune/select any fitted model.
dp[, interaction_context := fifelse(100 * 2^(single_a_log2 + single_b_log2) < 20,
  'Independent effect below20%: limited interaction headroom', 'Independent effect at least20%')]
write_table(dp[, .(rows = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
  interaction_r = cor_safe(observed_interaction, predicted_interaction),
  bliss_pp_r = cor_safe(observed_bliss_pp, predicted_bliss_pp),
  meaningful_synergy = sum(observed_bliss_pp >= 10),
  recovered_synergy = sum(observed_bliss_pp >= 10 & predicted_bliss_pp >= 10),
  false_synergy = sum(observed_bliss_pp < 10 & predicted_bliss_pp >= 10)),
  by = .(scope, input, model, interaction_context)], 'combination_interaction_contexts.tsv')
status('all bridge models and predictions saved; generating consolidated figures')

# One plotting function handles genuinely repeated observed/predicted comparisons.
save_plot <- function(p, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, '.pdf')), p, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(figure_dir, paste0(name, '.png')), p, width = width, height = height, dpi = 130, limitsize = FALSE)
}
agreement <- function(z, x, y, facets, xlab, ylab, title) {
  stats <- z[, .(r = cor_safe(get(x), get(y)), n = .N), by = facets]
  pretty_model <- function(values) {
    values <- as.character(values)
    labels <- c(molecular_products = 'Molecular interaction', conditioned_products = 'Baseline-conditioned interaction',
      independent_effects = 'Independent effects', prism_transferred = 'Transferred PRISM decoder',
      almanac_matched = 'ALMANAC decoder', baseline_almanac = 'Predicted response -> ALMANAC decoder')
    answer <- unname(labels[values])
    answer[is.na(answer)] <- values[is.na(answer)]
    answer
  }
  p <- ggplot(z, aes(x = .data[[x]], y = .data[[y]])) + geom_point(size = .65, alpha = .4) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
    facet_wrap(as.formula(paste('~', paste(facets, collapse = '+'))), scales = 'free',
      ncol = if ('drug' %in% facets) 7L else 3L, labeller = labeller(model = pretty_model)) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = sprintf('r=%.2f; n=%d', r, n)),
      hjust = -.05, vjust = 1.3, inherit.aes = FALSE, size = 2.7) +
    labs(title = title, x = xlab, y = ylab) + theme_classic(base_size = 9) +
    theme(aspect.ratio = 1, strip.text = element_text(size = 8))
  if (x == 'observed_log2') p <- p + geom_vline(xintercept = -1, linetype = 2, color = '#757575', linewidth = .3) +
    geom_hline(yintercept = -1, linetype = 2, color = '#757575', linewidth = .3)
  p
}
mp <- copy(molecular_predictions[split == 'Held cells'])
mp[, `:=`(observed_cell_specific = observed - mean(observed),
  predicted_cell_specific = predicted - mean(predicted)), by = .(drug, dose, program)]
p <- agreement(mp, 'observed_cell_specific', 'predicted_cell_specific', 'drug',
  'Measured response: within-treatment/program cell contrast', 'Predicted response: same cell contrast',
  'Baseline -> recipient-specific molecular response | nine held cells')
save_plot(p, 'baseline_cell_specific_response', 18, 3.4)
baseline_plot <- p
primary <- fp[decoder == 'ERBB' & split == 'Held cells' & model %chin%
  c('Measured response', 'Predicted response', 'Mean response', 'Direct baseline')]
p <- agreement(primary, 'observed_log2', 'predicted_log2', c('model', 'drug'),
  'Measured log2(% vehicle / 100)', 'Predicted log2(% vehicle / 100)',
  'Full bridge | each drug shown | nine held cells | ERBB-excluded decoder')
save_plot(p, 'full_chain_all_focal_drugs', 18, 12.5)
function_plot <- p
training <- fp[decoder == 'ERBB' & split == 'Training cells' & model %chin%
  c('Measured response', 'Predicted response', 'Mean response', 'Direct baseline')]
p <- agreement(training, 'observed_log2', 'predicted_log2', c('model', 'drug'),
  'Measured log2(% vehicle / 100)', 'Predicted log2(% vehicle / 100)',
  'Training recipients | baseline/direct are in-sample | ERBB drugs excluded from decoder training')
save_plot(p, 'training_chain_all_focal_drugs', 18, 12.5)
cross <- fread(file.path(data_dir, 'function_across_drugs_same_cell_dose.tsv'))
p <- ggplot(cross[decoder == 'ERBB' & split == 'Held cells'], aes(model, r, color = factor(dose))) +
  geom_hline(yintercept = 0, color = 'grey75') + geom_boxplot(aes(group = model), outlier.shape = NA,
    color = 'grey40', fill = 'grey95') + geom_point(position = position_jitter(width = .12, seed = 20261011L)) +
  labs(title = 'Across drugs in the same recipient and dose | held cells', x = NULL,
    y = 'Pearson r across the seven focal drugs', color = 'Dose (uM)') + theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, 'across_drugs_same_cell_dose', 7.5, 4)

# Actual and baseline-predicted heatmaps share the exact sensitivity order and
# row tree. Above each, plot observed sensitivity as a bar, with drug/dose/cell.
held <- lp[split == 'Held cells'][order(drug, measured_pct, dose, depmap_id)]
rows <- match(held$observation_id, md$observation_id)
row_order <- hclust(dist(t(measured[rows, , drop = FALSE])), method = 'ward.D2')$order
annotations <- data.frame(Dose = factor(held$dose), Cell = factor(held$recipient), row.names = held$observation_id)
annotation_args <- list(`PRISM % vehicle` = ComplexHeatmap::anno_barplot(held$measured_pct, border = FALSE,
  gp = grid::gpar(fill = ifelse(held$measured_pct <= 50, '#B2182B', '#4D4D4D')),
  height = grid::unit(1.2, 'cm')), annotation_name_gp = grid::gpar(fontsize = 8))
heatmaps <- list()
for (mode in c('Actual', 'Predicted')) {
  x <- if (mode == 'Actual') measured[rows, ] else predicted_response[rows, ]
  x <- t(x)[row_order, , drop = FALSE]
  check(identical(colnames(x), held$observation_id), 'Heatmap annotation alignment failed')
  ht <- draw_heatmap(x, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    col_anno_df = annotations, col_anno_cols = c('Dose', 'Cell'), ann_top_args = annotation_args,
    show_row_names = TRUE, show_column_names = FALSE, hm_name = 'Displacement',
    value_range = c(-1, 0, 1) * max(abs(c(measured[rows, ], predicted_response[rows, ]))),
    value_colors = c('#2166AC', '#F7F7F7', '#B2182B'),
    title = focal_drugs,
    column_split = factor(held$drug, focal_drugs), cluster_column_slices = FALSE,
    row_names_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10))
  heatmaps[[mode]] <- ht
  pdf(file.path(figure_dir, paste0(tolower(mode), '_focal_response_heatmap.pdf')), width = 13, height = 7)
  ComplexHeatmap::draw(ht, column_title = paste(mode, 'program response | same gene bank/order | held cells'))
  dev.off()
  png(file.path(figure_dir, paste0(tolower(mode), '_focal_response_heatmap.png')), width = 1690, height = 910, res = 130)
  ComplexHeatmap::draw(ht, column_title = paste(mode, 'program response | same gene bank/order | held cells'))
  dev.off()
}
cp <- dp[scope == 'Held cells + pairs' &
  input %chin% c('Measured single outcomes', 'Measured responses | ALMANAC head', 'Baseline only | ALMANAC head') &
  model %chin% c('molecular_products', 'conditioned_products', 'independent_effects')]
p <- agreement(cp, 'observed_log2', 'predicted_log2', c('input', 'model'),
  'Measured combination log2(% vehicle / 100)', 'Predicted combination log2(% vehicle / 100)',
  'Combination transfer | held recipients AND held pairs | endpoints')
save_plot(p, 'combination_full_bridge', 10, 10)
combination_plot <- p
p <- agreement(cp[input == 'Measured single outcomes' & model != 'independent_effects'],
  'observed_bliss_pp', 'predicted_bliss_pp', 'model', 'Measured Bliss interaction (percentage points)',
  'Predicted Bliss interaction (percentage points)', 'Molecular interaction beyond measured single-agent effects')
save_plot(p, 'combination_interaction', 7.5, 3.8)
interaction_plot <- p
single_plot <- agreement(single_long[split == 'held_cells'], 'observed_log2', 'predicted_log2', 'model',
  'Measured single log2(% vehicle / 100)', 'Predicted single log2(% vehicle / 100)',
  'ALMANAC single-effect bridge | three held recipients | dashed lines:50% vehicle')
save_plot(single_plot, 'almanac_single_function_stage', 10, 3.8)
overview <- fm[decoder == 'ERBB' & split == 'Held cells']
ov <- melt(overview[, .(model, across_cells_fixed_treatment_r, median_across_drugs_same_cell_dose_r)],
  id.vars = 'model', variable.name = 'axis', value.name = 'r')
ov[, axis := fifelse(axis == 'across_cells_fixed_treatment_r', 'Across cells at identical drug/dose',
  'Across drugs at identical cell/dose')]
ov <- ov[is.finite(r)] # Mean controls have no cross-cell prediction, so r is undefined.
p <- ggplot(ov, aes(model, r, fill = axis)) + geom_col(position = 'dodge', width = .65) +
  geom_hline(yintercept = 0, color = 'grey60') + coord_flip() +
  labs(title = 'Where the baseline -> response -> function bridge succeeds or fails', x = NULL, y = 'Pearson r', fill = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', legend.text = element_text(size = 8))
save_plot(p, 'bridge_overview', 8, 4)
overview_plot <- p
# A short document connects the two stages; full7-drug details remain in their
# readable larger PDFs, while this document has fixed paper-sized pages.
pdf(file.path(figure_dir, 'bridge_report.pdf'), width = 11, height = 8.5)
print(overview_plot)
print(combination_plot)
print(interaction_plot)
print(single_plot)
ComplexHeatmap::draw(heatmaps[['Actual']], column_title = 'Actual molecular responses | held cells')
ComplexHeatmap::draw(heatmaps[['Predicted']], column_title = 'Baseline-predicted molecular responses | held cells')
dev.off()

# Meaningful replay, split and schema checks are recorded alongside the results.
reference50 <- fread(file.path(decoder_root, 'data', 'evaluated_predictions.tsv.gz'))
reference50 <- reference50[bank == '43 current / 34 programs' & calibration == 'Calibrated' &
  evaluation == 'Held cells' & drug %chin% focal_drugs]
replay <- fp[split == 'Held cells' & model == 'Measured response']
expected <- reference50$predicted_log2[match(paste(replay$decoder, replay$observation_id),
  paste(reference50$family, reference50$observation_id))]
check(!anyNA(expected) && max(abs(expected - replay$predicted_log2)) < 1e-10,
  'Frozen50 decoder predictions changed')
write_table(data.table(rows = nrow(replay), maximum_absolute_difference = max(abs(expected - replay$predicted_log2))),
  'frozen_decoder_replay.tsv')
original_support <- fread('/work/results/program_pilot/20_multidrug_baseline_response_function_union_programs/data/D_combination_predictions.tsv.gz')
original_rows <- unique(original_support$source_row)
check(setequal(cd$source_row[base_test], original_rows), 'Combination test coverage differs from20')
write_table(data.table(original20_test_rows = length(original_rows), current51_test_rows = length(base_test),
  exact_source_row_match = TRUE), 'original_combination_coverage_replay.tsv')
write_table(data.table(check = c('43 hard scoring replay', 'baseline held IDs excluded',
  'combination held cells excluded', 'held pair folds disjoint', 'shared controls purged',
  'symmetric bilinear representation', 'zero interaction for zero partner', 'finite full-chain predictions',
  'exact50 measured-response prediction replay', 'exact20 held combination source rows'),
  passed = TRUE), 'acceptance.tsv')
code_paths <- c(script_path, '/work/src/program_pilot/response_models.R',
  '/work/src/program_pilot/nonlinear_response.R', '/work/almanac_broad9_2026-10-07/model.R', '/work/src/R/lib_heatmap.R')
write_table(data.table(path = code_paths,
  sha256 = vapply(code_paths, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'code_hashes.tsv')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
report <- c('51: frozen34-program baseline -> displacement -> function -> combination bridge',
  'Retrospective:43 discovery used all-panel outcomes; the original nine held recipients are unchanged.',
  'Baseline response models use pooled native dose rows, shared baseline slopes and baseline-by-logdose slopes.',
  'Each response model and direct comparator selects ridge regularization only in grouped fitting-cell CV.',
  'Measured/predicted responses use exactly the unchanged calibrated50 ERBB/MEK functional heads.',
  'Training-cell baseline/direct predictions are in-sample, not independent validation.',
  'For the functional head, excluded-family drugs in training cells were not fitted; drug_scope identifies them.',
  'The ERBB-excluded decoder is the predeclared primary full-chain/combination head; both heads have single-drug tables.',
  'Across-drug correlation compares different drugs at fixed cell AND dose; fixed-treatment correlation compares cells.',
  'Full baseline combinations use neither measured molecular responses nor single outcomes of test recipients.',
  'ALMANAC single calibration is trained on development singles only; interaction fitting is separate.',
  'ALMANAC-matched SVR uses the same fixed C10/gamma1p settings; grouped-cell OOF affine calibration uses development singles.',
  'The original transferred-PRISM single head is retained as the matched-assay comparator.',
  'Treatment-average single controls use development-cell means and logdose interpolation, with no extrapolation.',
  'Their supported contexts and the matched-context comparison are explicitly recorded.',
  'Interaction-only diagnostics use measured single outcomes to isolate the molecular interaction head.',
  'Predicted Bliss points use each model own predicted additive endpoint; conditional conversion with measured singles is separate.',
  'Drug/drug combination evidence does not validate HGF+lapatinib rescue. Decoder zero-input output is recorded, not altered.',
  '', 'Held focal molecular results:', capture.output(print(molecular_metrics[split == 'Held cells'])),
  '', 'Held focal function:', capture.output(print(fm[split == 'Held cells'])),
  '', 'Combination endpoints and interactions:', capture.output(print(cm)),
  '', 'Combination split audit:', capture.output(print(rbindlist(audits))))
report <- c(report, '', 'Single-function stages in ALMANAC:', capture.output(print(single_metrics)),
  '', 'Baseline response variation retained:',
  capture.output(print(fread(file.path(data_dir, 'molecular_variation_by_drug.tsv'))[split == 'Held cells'])),
  '', 'Interpretation:',
  'The seven-drug cross-drug ordering is largely reproduced by treatment-average responses.',
  'Recipient-specific sensitivity remains weak; baseline-response fits shrink most recipient variation toward treatment means.',
  'ALMANAC-matched single-function learning improves measured-response endpoint prediction over transferred PRISM.',
  'Baseline-predicted responses expose a second transfer problem: their distribution differs from measured responses.',
  'Molecular interaction signal is modest and does not improve overall endpoints over independent effects.',
  'This implements the complete drug bridge and its controls; the current results do not establish the full grant wager.')
writeLines(report, file.path(root, 'REPORT.txt'))
status('COMPLETE')
