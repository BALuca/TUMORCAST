#!/usr/bin/env Rscript
# One joint objective learns individual effects and pair modifications.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(glmnet)
  library(digest)
})
root <- '/work/results/program_pilot/68_joint_tahoe_prism_combination_learning'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
single_path <- '/work/results/program_pilot/62_separate_response_interaction_bridge/data/single_predictions.tsv.gz'
script_path <- '/work/src/program_pilot/68_joint_tahoe_prism_combination_learning.R'
lambda_grid <- 10^seq(1, -3, length.out = 24L)
alpha <- .2
interaction_penalties <- c(1, 4)
set.seed(68011L)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn = 1)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
}
status <- function(x) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'STATUS.txt'))
}
options(error = function() {
  status('FAILED; see run log')
  traceback(6)
  quit(save = 'no', status = 1L)
})
status('Projecting both assays through the same measured genes and original40 weights')
panel <- qread(panel_path)
molecular <- qread(file.path(prior, 'data/matched_molecular.qs'))
bank <- qread(file.path(prior, 'data/representation_scores.qs'))
W <- qread(file.path(prior, 'models/original40_weights.qs'))
common <- intersect(intersect(rownames(W), rownames(molecular$displacement)), panel$genes)
stopifnot(length(common) > 500L, ncol(W) == 40L)
coverage <- data.table(program = colnames(W), measured_genes = colSums(W[common, ] != 0),
  all_genes = colSums(W != 0), absolute_weight_fraction = colSums(abs(W[common, ])) / colSums(abs(W)))
write_table(coverage, 'common_gene_coverage')
write_table(data.table(gene = common), 'common_genes')
bt_profiles <- bank$metadata[time_h == 3, profile_id]
bt_scores <- t(as.matrix(crossprod(W[common, ], molecular$displacement[common, bt_profiles, drop = FALSE])))
labels <- copy(panel$sensitivity[is.finite(prism_abundance_pct) & prism_abundance_pct > 0])
stopifnot(!anyDuplicated(labels$observation_id))
tp_parts <- list()
input_paths <- character()
for (i in seq_along(panel$input_paths)) {
  id <- panel$cell_map$depmap_id[i]
  ids <- labels[depmap_id == id, observation_id]
  if (!length(ids)) next
  input_paths <- c(input_paths, panel$input_paths[i])
  obj <- qread(panel$input_paths[i])
  idx <- match(ids, rownames(obj$response))
  stopifnot(!anyNA(idx), all(common %in% colnames(obj$response)),
    identical(obj$normalizer, 'paired_trimmed_median_log_ratio_cp10k_v1'))
  tp_parts[[id]] <- as.matrix(obj$response[idx, common, drop = FALSE] %*% W[common, ])
}
tp_scores <- do.call(rbind, tp_parts)
rownames(tp_scores) <- unlist(lapply(tp_parts, rownames), use.names = FALSE)
tp_scores <- tp_scores[labels$observation_id, , drop = FALSE]
stopifnot(identical(rownames(tp_scores), labels$observation_id), all(is.finite(tp_scores)))
rm(obj, tp_parts, molecular)

# Different assays have different units. RMS uses only individual RNA, preserves
# zero, and is frozen before any outcome fitting. BT20 deviations are shrunk.
bt_rms <- sqrt(colMeans(bt_scores^2))
tp_rms <- sqrt(colMeans(tp_scores^2))
stopifnot(all(bt_rms > 0), all(tp_rms > 0))
bt <- sweep(bt_scores, 2, bt_rms, '/')
tp <- sweep(tp_scores, 2, tp_rms, '/')
smooth_basis <- function(x) {
  sat <- tanh(x)
  colnames(sat) <- paste0(colnames(x), '_sat')
  cbind(x, sat)
}
cross_features <- function(a, b) {
  ij <- which(upper.tri(matrix(0, ncol(a), ncol(a)), diag = TRUE), arr.ind = TRUE)
  x <- vapply(seq_len(nrow(ij)), function(k) {
    i <- ij[k, 1]
    j <- ij[k, 2]
    (a[, i] * b[, j] + a[, j] * b[, i]) / 2
  }, numeric(nrow(a)))
  colnames(x) <- paste0('g:', colnames(a)[ij[, 1]], '__', colnames(a)[ij[, 2]])
  x
}
bt_basis <- smooth_basis(bt)
tp_basis <- smooth_basis(tp)
singles <- fread(single_path)[representation == 'Original40' & time_h == 3]
singles <- singles[, .(sample_id, source, drug, nominal_dose_uM, biological_replicate, profile_id, ratio)]
stopifnot(!anyDuplicated(paste(singles$source, singles$sample_id)), all(singles$ratio > 0))
sx <- bt_basis[singles$profile_id, , drop = FALSE]
obs <- fread(file.path(prior, 'data/matched_molecular_links.tsv.gz'))[time_h == 3 & molecular_supported == TRUE]
setorder(obs, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)
stopifnot(nrow(obs) == 496L, !anyDuplicated(obs$observation_id), all(obs$ratio_AB > 0))
a <- bt_basis[obs$profile_A, , drop = FALSE]
b <- bt_basis[obs$profile_B, , drop = FALSE]
gx <- cross_features(a, b)
stopifnot(ncol(gx) == 3240L, max(abs(gx - cross_features(b, a))) == 0,
  max(abs(cross_features(a * 0, b))) == 0, max(abs(smooth_basis(bt * 0))) == 0)
partners <- sort(unique(obs$partner))
models <- c('Joint_Tahoe_PRISM', 'Additive_Tahoe_PRISM', 'Joint_BT20_only')
inputs <- c(panel_path, single_path, file.path(prior, c('data/matched_molecular_links.tsv.gz',
  'data/matched_molecular.qs', 'data/representation_scores.qs', 'models/original40_weights.qs')), input_paths)
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')),
  'input_manifest')
write_table(labels, 'Tahoe_PRISM_training_labels')
write_table(singles, 'BT20_training_single_labels')
write_table(data.table(source = c('Tahoe_PRISM_singles', 'BT20_singles', 'BT20_combinations'),
  observations = c(nrow(labels), nrow(singles), nrow(obs)),
  cells = c(uniqueN(labels$depmap_id), 1, 1),
  interventions_or_pairs = c(uniqueN(labels$drug), uniqueN(singles$drug), uniqueN(paste(obs$anchor, obs$partner)))),
  'data_counts')
qsave(list(common_genes = common, bt_scores = bt_scores, tp_scores = tp_scores, bt_rms = bt_rms,
  tp_rms = tp_rms, lambda_grid = lambda_grid, alpha = alpha, interaction_penalties = interaction_penalties,
  bt20_deviation_penalty = 4, original_weight_sha256 = bank$original_weight_sha256,
  gene_coverage = coverage, single_basis = 'x and tanh(x)', primary_RNA_time_h = 3,
  objective = 'joint singles and combinations; equal source totals; no measured single offsets'),
  file.path(root, 'models/settings.qs'))

# Source weights prevent 9,000 external rows from overwhelming the BT20 assay.
unit_weights <- function(group) {
  w <- 1 / as.numeric(table(group)[group])
  w / sum(w)
}
tp_weights <- unit_weights(paste(labels$depmap_id, labels$drug))
single_weights <- unit_weights(paste(singles$source, singles$drug))
tp_y <- log2(labels$prism_abundance_pct / 100)
single_y <- log2(singles$ratio)
combo_y <- log2(obs$ratio_AB)
center_context <- function(x, md) {
  ans <- x
  groups <- interaction(md$anchor, md$nominal_A_uM, md$biological_replicate, drop = TRUE)
  for (ii in split(seq_len(nrow(md)), groups)) {
    if (is.matrix(x)) ans[ii, ] <- sweep(x[ii, , drop = FALSE], 2, colMeans(x[ii, , drop = FALSE]), '-')
    else ans[ii] <- x[ii] - mean(x[ii])
  }
  ans
}
objective <- function(y, pred, md) {
  w <- unit_weights(paste(md$anchor, md$partner, md$biological_replicate))
  .5 * colSums((pred - y)^2 * w) + .5 * colSums((center_context(pred, md) - center_context(y, md))^2 * w)
}

# Sparse rows for singles contain exactly zero interaction terms. The same f
# coefficients occur once for singles and as f(a)+f(b) in combination rows.
designs <- list()
for (model in models) {
  external <- model != 'Joint_BT20_only'
  joint <- model != 'Additive_Tahoe_PRISM'
  fs <- sx
  fc <- a + b
  ft <- tp_basis
  colnames(fs) <- colnames(fc) <- colnames(ft) <- paste0('f:', colnames(sx))
  if (external) {
    dev_s <- sx
    dev_c <- a + b
    dev_t <- tp_basis * 0
    colnames(dev_s) <- colnames(dev_c) <- colnames(dev_t) <- paste0('BT20:', colnames(sx))
    fs <- cbind(fs, dev_s)
    fc <- cbind(fc, dev_c)
    ft <- cbind(ft, dev_t)
  }
  fs <- Matrix(fs, sparse = TRUE)
  ft <- Matrix(ft, sparse = TRUE)
  if (joint) {
    z <- Matrix(0, nrow(sx), ncol(gx), sparse = TRUE, dimnames = list(NULL, colnames(gx)))
    fs <- cbind(fs, z)
    ft <- cbind(ft, Matrix(0, nrow(tp), ncol(gx), sparse = TRUE, dimnames = list(NULL, colnames(gx))))
    fc <- cbind(fc, gx)
  }
  designs[[model]] <- list(s = fs, t = ft, c = fc, external = external, joint = joint)
}

# Fit paths are cached for each unordered set of excluded partners, so the
# same inner exclusion is never fitted twice. Each final model is flushed.
cache <- new.env(parent = emptyenv())
fit_joint <- function(model, excluded, multiplier) {
  key <- paste(model, paste(sort(excluded), collapse = '+'), multiplier, sep = '|')
  if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
  d <- designs[[model]]
  tr <- which(!obs$partner %chin% excluded)
  w <- unit_weights(paste(obs$anchor[tr], obs$partner[tr], obs$biological_replicate[tr]))
  cx <- d$c[tr, , drop = FALSE]
  cy <- combo_y[tr]
  xx <- rbind(d$s, Matrix(cx, sparse = TRUE), Matrix(center_context(cx, obs[tr]), sparse = TRUE))
  yy <- c(single_y, cy, center_context(cy, obs[tr]))
  ww <- c(single_weights, .5 * w, .5 * w)
  if (d$external) {
    xx <- rbind(d$t, xx)
    yy <- c(tp_y, yy)
    ww <- c(tp_weights, ww)
  }
  ww <- ww / sum(ww)
  rms <- sqrt(as.numeric(crossprod(ww, xx^2)))
  keep <- rms > 1e-10 & is.finite(rms)
  xx <- xx[, keep, drop = FALSE] %*% Diagonal(x = 1 / rms[keep])
  block <- colnames(d$c)[keep]
  penalty <- ifelse(startsWith(block, 'g:'), multiplier, ifelse(startsWith(block, 'BT20:'), 4, 1))
  yscale <- sqrt(sum(ww * yy^2))
  fit <- glmnet(xx, yy / yscale, weights = ww, alpha = alpha, lambda = lambda_grid,
    penalty.factor = penalty, intercept = FALSE, standardize = FALSE, thresh = 1e-6, maxit = 30000L)
  stopifnot(fit$jerr == 0L)
  coefficients <- matrix(0, ncol(d$c), length(lambda_grid), dimnames = list(colnames(d$c), NULL))
  coefficients[keep, ] <- sweep(as.matrix(coef(fit, s = lambda_grid))[-1, , drop = FALSE] * yscale,
    1, rms[keep], '/')
  out <- list(coefficients = coefficients, prediction = d$c %*% coefficients,
    training_rows = tr, rms = rms, yscale = yscale)
  assign(key, out, cache)
  out
}

predictions <- single_predictions <- tuning <- validation <- coefficients <- list()
counter <- 0L
for (held in partners) {
  status(paste('Joint fitting: excluding all', held, 'combinations'))
  train_partners <- setdiff(partners, held)
  te <- which(obs$partner == held)
  tr <- which(obs$partner != held)
  source_ids <- function(x) unique(unlist(strsplit(x, ';', fixed = TRUE)))
  overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]), source_ids(obs$source_sample_ids_AB[te]))
  stopifnot(length(overlap) == 0L)
  for (model in models) {
    mult <- if (designs[[model]]$joint) interaction_penalties else 1
    losses <- matrix(0, length(lambda_grid), length(mult))
    for (j in seq_along(mult)) {
      for (inner in train_partners) {
        candidate <- fit_joint(model, c(held, inner), mult[j])
        va <- which(obs$partner == inner)
        losses[, j] <- losses[, j] + objective(combo_y[va], candidate$prediction[va, , drop = FALSE], obs[va]) /
          length(train_partners)
      }
    }
    choice <- which(losses == min(losses), arr.ind = TRUE)[1, ]
    fitted <- fit_joint(model, held, mult[choice[2]])
    beta <- fitted$coefficients[, choice[1]]
    d <- designs[[model]]
    pred <- as.numeric(d$c %*% beta)
    g_cols <- startsWith(names(beta), 'g:')
    gpred <- if (any(g_cols)) as.numeric(d$c[, g_cols, drop = FALSE] %*% beta[g_cols]) else rep(0, nrow(obs))
    fcols <- startsWith(names(beta), 'f:')
    devcols <- startsWith(names(beta), 'BT20:')
    single_beta <- beta[fcols] + if (any(devcols)) beta[devcols] else 0
    fa <- as.numeric(a %*% single_beta)
    fb <- as.numeric(b %*% single_beta)
    stopifnot(max(abs(pred - fa - fb - gpred)) < 1e-8, all(is.finite(2^pred)))
    saved <- list(model = model, held_partner = held, beta = beta, lambda = lambda_grid[choice[1]],
      interaction_penalty = mult[choice[2]], alpha = alpha,
      training_combination_ids = obs$observation_id[tr], training_single_ids = singles$sample_id,
      training_PRISM_ids = if (d$external) labels$observation_id else character(),
      original_weight_sha256 = bank$original_weight_sha256)
    path <- file.path(root, 'models', paste0(held, '_', model, '.qs'))
    qsave(saved, path)
    replay <- as.numeric(d$c %*% qread(path)$beta)
    stopifnot(max(abs(pred - replay)) < 1e-10)
    counter <- counter + 1L
    out <- obs[, .(observation_id, condition_id, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)]
    out[, `:=`(model = model, fold = held, split = fifelse(partner == held, 'held', 'training'),
      predicted_log2 = pred, predicted_A_log2 = fa, predicted_B_log2 = fb, predicted_g = gpred)]
    # Prediction-only artifacts are written and hashed before outcomes are joined.
    write_table(out, paste0('predictions_', held, '_', model))
    prediction_hash <- digest(out, algo = 'sha256')
    out[, `:=`(ratio_AB = obs$ratio_AB, ratio_A = obs$ratio_A, ratio_B = obs$ratio_B)]
    predictions[[counter]] <- out
    sp <- singles[, .(sample_id, source, drug, nominal_dose_uM, biological_replicate, ratio)]
    sp[, `:=`(model = model, fold = held, predicted_log2 = as.numeric(d$s %*% beta))]
    tp_pred <- if (d$external) as.numeric(d$t %*% beta) else rep(NA_real_, nrow(labels))
    if (d$external) {
      ext <- labels[, .(sample_id = observation_id, source = 'Tahoe_PRISM', drug, nominal_dose_uM = dose,
        biological_replicate = NA_integer_, ratio = prism_abundance_pct / 100)]
      ext[, `:=`(model = model, fold = held, predicted_log2 = tp_pred)]
      sp <- rbind(sp, ext)
    }
    single_predictions[[counter]] <- sp
    grid <- CJ(interaction_penalty = mult, lambda_index = seq_along(lambda_grid))
    grid[, `:=`(lambda = lambda_grid[lambda_index], inner_objective = as.vector(losses),
      selected = lambda_index == choice[1] & interaction_penalty == mult[choice[2]], model = model, fold = held)]
    tuning[[counter]] <- grid
    validation[[counter]] <- data.table(model = model, fold = held, source_overlap = length(overlap),
      replay_difference = max(abs(pred - replay)), prediction_sha256 = prediction_hash,
      training_combinations = length(tr), training_conditions = uniqueN(obs$condition_id[tr]),
      training_pairs = uniqueN(paste(obs$anchor[tr], obs$partner[tr])),
      active_single_terms = sum(beta[!g_cols] != 0), active_interactions = sum(beta[g_cols] != 0),
      lambda = saved$lambda, interaction_penalty = saved$interaction_penalty)
    coefficients[[counter]] <- data.table(model = model, fold = held, term = names(beta), coefficient = beta)
    write_table(rbindlist(predictions), 'predictions')
    write_table(rbindlist(validation), 'validation')
    write_table(rbindlist(tuning), 'tuning')
    message('Saved ', model, ' | ', held, ' | active interaction terms: ', sum(beta[g_cols] != 0))
  }
}
write_table(rbindlist(single_predictions), 'single_predictions')
write_table(rbindlist(coefficients), 'coefficients')
write_table(data.table(path = script_path, sha256 = digest(file = script_path, algo = 'sha256')), 'script_manifest')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
writeLines(system('git -C /work rev-parse HEAD', intern = TRUE), file.path(root, 'git_revision.txt'))
status('Joint models complete; predictions saved for all six excluded partners')
