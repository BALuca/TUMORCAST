#!/usr/bin/env Rscript
# Target constituent experiments remain available. Only their joint outcomes
# are excluded, across both anchors, all doses and biological experiments.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(e1071)
  library(digest)
})
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(61011L)
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
lambda_grid <- 10^seq(3, -4, length.out = 22L)
svr_grid <- CJ(cost = c(1, 10, 100), gamma_factor = c(.25, 1))
rna_times <- c(3, 24)
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', name), sep = '\t', na = 'NA')
status <- function(x) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', x)
  writeLines(x, file.path(root, 'MODEL_STATUS.txt'))
}
links_path <- file.path(root, 'data/matched_molecular_links.tsv.gz')
molecular_path <- file.path(root, 'data/matched_molecular.qs')
scores_path <- file.path(root, 'data/representation_scores.qs')
stopifnot(file.exists(links_path), file.exists(molecular_path), file.exists(scores_path))
all_links <- fread(links_path)
molecular <- qread(molecular_path)
bank <- qread(scores_path)
stopifnot(is.matrix(bank$original40), ncol(bank$original40) == 40L,
  all(is.finite(bank$original40)), is.matrix(bank$biology))
score_sets <- list(Original40 = t(bank$original40), Biology = t(bank$biology))
# Both constituent experiments, including the excluded partner, are available.
# Unsupervised reference PCA uses those individual measurements, not combination
# outcomes. This reference is twelve axes, not a complete-gene oracle.
gene_pca <- prcomp(t(molecular$displacement), center = TRUE, scale. = FALSE, rank. = 12L)
score_sets$GenePC12 <- t(gene_pca$x[, 1:12, drop = FALSE])
qsave(gene_pca, file.path(root, 'models/gene_reference_pca.qs'))
qsave(list(lambda_grid = lambda_grid, svr_grid = svr_grid, rna_times = rna_times,
  primary_time_h = 3, exclusion = 'entire partner joint outcomes; individuals available',
  representation = 'original40 signed weights; measured-gene contribution without renormalization'),
  file.path(root, 'models/model_settings.qs'))

# Same biological-replicate and measured-single controls for every model.
make_controls <- function(z) {
  a <- z$ratio_A
  b <- z$ratio_B
  cbind(single_A = a, single_B = b, single_product = a * b,
    single_A_sq = a^2, single_B_sq = b^2,
    logdose_A = log10(z$nominal_A_uM), logdose_B = log10(z$nominal_B_uM),
    model.matrix(~anchor + factor(biological_replicate), z)[, -1, drop = FALSE])
}
make_response <- function(scores, z, quadratic = FALSE) {
  stopifnot(all(z$profile_A %in% colnames(scores)), all(z$profile_B %in% colnames(scores)))
  a <- t(scores[, z$profile_A, drop = FALSE])
  b <- t(scores[, z$profile_B, drop = FALSE])
  sums <- a + b
  colnames(sums) <- paste0('sum_', colnames(a))
  if (quadratic) {
    ij <- which(upper.tri(matrix(0, ncol(a), ncol(a)), diag = TRUE), arr.ind = TRUE)
    product <- vapply(seq_len(nrow(ij)), function(k) {
      i <- ij[k, 1]
      j <- ij[k, 2]
      (a[, i] * b[, j] + a[, j] * b[, i]) / 2
    }, numeric(nrow(a)))
    colnames(product) <- paste0('cross_', colnames(a)[ij[, 1]], '__', colnames(a)[ij[, 2]])
    cbind(sums, product)
  } else {
    difference <- abs(a - b)
    product <- a * b
    colnames(difference) <- paste0('absdifference_', colnames(a))
    colnames(product) <- paste0('product_', colnames(a))
    cbind(sums, difference, product)
  }
}
safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 4L || sd(a[ok]) < 1e-10 || sd(b[ok]) < 1e-10) return(NA_real_)
  cor(a[ok], b[ok])
}
fit_ridge <- function(x, y, metadata) {
  mean_x <- colMeans(x)
  sd_x <- apply(x, 2, sd)
  keep <- is.finite(sd_x) & sd_x > 1e-10
  xx <- sweep(sweep(x[, keep, drop = FALSE], 2, mean_x[keep], '-'), 2, sd_x[keep], '/')
  # Equal total weight per partner/anchor/experiment; dose density cannot
  # increase that context's total contribution to the ridge fitting objective.
  groups <- interaction(metadata$partner, metadata$anchor, metadata$biological_replicate, drop = TRUE)
  w <- 1 / as.numeric(table(groups)[as.character(groups)])
  f <- glmnet(xx, y, alpha = 0, lambda = lambda_grid, weights = w,
    standardize = FALSE, intercept = TRUE, thresh = 1e-8, maxit = 100000L)
  list(fit = f, mean = mean_x[keep], sd = sd_x[keep], keep = keep)
}
predict_ridge <- function(f, x) {
  xx <- sweep(sweep(x[, f$keep, drop = FALSE], 2, f$mean, '-'), 2, f$sd, '/')
  as.matrix(predict(f$fit, newx = xx, s = lambda_grid))
}
balanced_loss <- function(y, pred, md) {
  # Internal fitting criterion only; the scientific report uses effect sizes,
  # dose curves, cell loss, and consistency across biological experiments.
  err <- (pred - y)^2
  groups <- interaction(md$anchor, md$biological_replicate, drop = TRUE)
  mean(vapply(split(err, groups), mean, numeric(1)))
}

predictions <- tuning <- capacity <- coefficients <- list()
counter <- 0L
for (tt in rna_times) {
  obs <- copy(all_links[time_h == tt & molecular_supported == TRUE])
  setorder(obs, partner, anchor, biological_replicate, nominal_A_uM, nominal_B_uM)
  stopifnot(nrow(obs) > 100L, !anyDuplicated(obs$observation_id), all(is.finite(obs$excess_inhibition)))
  partners <- sort(unique(obs$partner))
  stopifnot(length(partners) == 6L, uniqueN(obs$anchor) == 2L)
  controls <- make_controls(obs)
  class <- fcase(obs$partner %chin% c('lapatinib', 'neratinib'), 'ERBB',
    obs$partner %chin% c('selumetinib', 'trametinib'), 'MEK',
    obs$partner == 'saracatinib', 'SRC', obs$partner == 'nvp-tae684', 'TAE684_multikinase')
  moa <- model.matrix(~factor(class) - 1)
  designs <- list(Singles_ridge = controls, Singles_SVR = controls,
    Singles_MOA_SVR = cbind(controls, moa),
    Original40_ridge = cbind(controls, make_response(score_sets$Original40, obs, quadratic = TRUE)),
    Original40_SVR = cbind(controls, make_response(score_sets$Original40, obs)),
    Biology_SVR = cbind(controls, make_response(score_sets$Biology, obs)),
    GenePC12_SVR = cbind(controls, make_response(score_sets$GenePC12, obs)))
  stopifnot(all(vapply(designs, function(x) all(is.finite(x)), logical(1))))
  # Only training indices expose outcomes to a fitting operation below.
  outcomes <- obs$excess_inhibition
  for (excluded in partners) {
    tr <- which(obs$partner != excluded)
    te <- which(obs$partner == excluded)
    stopifnot(!excluded %in% obs$partner[tr],
      all(obs$profile_A[te] %in% colnames(score_sets$Original40)),
      all(obs$profile_B[te] %in% colnames(score_sets$Original40)))
    inner_partners <- partners[partners != excluded]
    for (model_name in names(designs)) {
      counter <- counter + 1L
      status(paste(tt, 'h | excluded', excluded, '|', model_name))
      x <- designs[[model_name]]
      ridge_model <- grepl('_ridge$', model_name)
      candidates <- if (ridge_model) length(lambda_grid) else nrow(svr_grid)
      oof <- matrix(NA_real_, length(tr), candidates)
      losses <- matrix(NA_real_, length(inner_partners), candidates)
      for (j in seq_along(inner_partners)) {
        ii <- which(obs$partner[tr] != inner_partners[j])
        vv <- which(obs$partner[tr] == inner_partners[j])
        if (ridge_model) {
          f <- fit_ridge(x[tr[ii], , drop = FALSE], outcomes[tr[ii]], obs[tr[ii]])
          pp <- predict_ridge(f, x[tr[vv], , drop = FALSE])
        } else {
          pp <- vapply(seq_len(candidates), function(k) {
            f <- fit_program_svr(x[tr[ii], , drop = FALSE], outcomes[tr[ii]],
              'radial', svr_grid$cost[k], svr_grid$gamma_factor[k] / ncol(x))
            predict_program_svr(f, x[tr[vv], , drop = FALSE])
          }, numeric(length(vv)))
        }
        oof[vv, ] <- pp
        losses[j, ] <- vapply(seq_len(candidates), function(k) {
          balanced_loss(outcomes[tr[vv]], pp[, k], obs[tr[vv]])
        }, numeric(1))
      }
      stopifnot(all(is.finite(oof)), all(is.finite(losses)))
      # Minimize held-training-partner loss. No one-SE preference for a nearly
      # constant fit; amplitude calibration remains visible in diagnostics.
      selected <- which.min(colMeans(losses))
      tuning[[counter]] <- data.table(time_h = tt, held_partner = excluded, model = model_name,
        candidate = seq_len(candidates), lambda = if (ridge_model) lambda_grid else NA_real_,
        cost = if (ridge_model) NA_real_ else svr_grid$cost,
        gamma_factor = if (ridge_model) NA_real_ else svr_grid$gamma_factor,
        internal_grouped_loss = colMeans(losses), selected = seq_len(candidates) == selected)
      if (ridge_model) {
        f <- fit_ridge(x[tr, , drop = FALSE], outcomes[tr], obs[tr])
        p_all <- predict_ridge(f, x)
        pp <- p_all[, selected]
        cf <- as.matrix(coef(f$fit, s = lambda_grid[selected]))
        coefficients[[counter]] <- data.table(time_h = tt, held_partner = excluded, model = model_name,
          term = rownames(cf), standardized_coefficient = cf[, 1])
        capacity[[counter]] <- rbindlist(lapply(seq_len(candidates), function(k) {
          data.table(time_h = tt, held_partner = excluded, model = model_name,
            candidate = k, lambda = lambda_grid[k], selected = k == selected,
            training_r = safe_cor(outcomes[tr], p_all[tr, k]),
            inner_partner_r = safe_cor(outcomes[tr], oof[, k]),
            observed_effect_p95 = quantile(outcomes[tr], .95),
            fitted_effect_p95 = quantile(p_all[tr, k], .95))
        }))
      } else {
        f <- fit_program_svr(x[tr, , drop = FALSE], outcomes[tr], 'radial',
          svr_grid$cost[selected], svr_grid$gamma_factor[selected] / ncol(x))
        pp <- predict_program_svr(f, x)
        capacity[[counter]] <- data.table(time_h = tt, held_partner = excluded, model = model_name,
          candidate = selected, lambda = NA_real_, selected = TRUE,
          training_r = safe_cor(outcomes[tr], pp[tr]),
          inner_partner_r = safe_cor(outcomes[tr], oof[, selected]),
          observed_effect_p95 = quantile(outcomes[tr], .95),
          fitted_effect_p95 = quantile(pp[tr], .95))
      }
      saved <- list(fit = f, selection = tuning[[counter]][selected == TRUE],
        columns = colnames(x), training_partners = inner_partners,
        target_individual_profiles_available = unique(unlist(obs[te, .(profile_A, profile_B)])))
      key <- paste0(tt, 'h_', excluded, '_', model_name)
      qsave(saved, file.path(root, 'models', paste0(key, '.qs')), preset = 'high')
      pred <- rbind(
        data.table(observation_id = obs$observation_id[te], split = 'held_partner', predicted_extra_pp = pp[te]),
        data.table(observation_id = obs$observation_id[tr], split = 'training_fit', predicted_extra_pp = pp[tr]),
        data.table(observation_id = obs$observation_id[tr], split = 'inner_partner_OOF', predicted_extra_pp = oof[, selected]))
      pred[, `:=`(time_h = tt, held_partner = excluded, model = model_name)]
      # Commit predictions before joining any target combination measurements.
      write_table(pred, paste0('prediction_', key, '.tsv.gz'))
      predictions[[counter]] <- pred
      write_table(rbindlist(tuning), 'model_training_selection.tsv')
      write_table(rbindlist(capacity, fill = TRUE), 'model_training_capacity.tsv')
    }
  }
  counter <- counter + 1L
  predictions[[counter]] <- obs[, .(observation_id, time_h = tt, held_partner = partner,
    model = 'Independent_singles', split = 'held_partner', predicted_extra_pp = 0)]
}
pred <- rbindlist(predictions, use.names = TRUE)
write_table(pred, 'predictions_before_outcome_join.tsv.gz')
writeLines(digest(pred, algo = 'sha256'), file.path(root, 'models/prediction_sha256.txt'))
write_table(rbindlist(coefficients), 'original40_ridge_coefficients.tsv.gz')

# Evaluation only begins after all primary and diagnostic predictions are frozen.
observed <- fread(links_path)[molecular_supported == TRUE]
observed[, `:=`(observed_extra_pp = excess_inhibition,
  relative_A = ratio_A, relative_B = ratio_B, relative_AB = ratio_AB,
  dose_anchor = dose_A_uM, dose_partner = dose_B_uM)]
pred <- merge(pred, observed, by = c('observation_id', 'time_h'), all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(pred$observed_extra_pp), all(is.finite(pred$predicted_extra_pp)))
pred[, predicted_relative_AB := relative_A * relative_B - predicted_extra_pp / 100]
pred[, impossible_negative_cell_count := predicted_relative_AB <= 0]
pred[, predicted_GR_AB := NA_real_]
pred[predicted_relative_AB > 0,
  predicted_GR_AB := 2^(log2(predicted_relative_AB * ncontrol72_AB / n0_AB) /
    log2(ncontrol72_AB / n0_AB)) - 1]
write_table(pred, 'predictions.tsv.gz')

# Average technical-normalized biological observations only for a separate
# condition summary. Preserve every biological experiment in the main table.
held <- pred[split == 'held_partner']
condition <- held[, .(observed_extra_pp = mean(observed_extra_pp),
  predicted_extra_pp = mean(predicted_extra_pp),
  relative_A = mean(relative_A), relative_B = mean(relative_B), relative_AB = mean(relative_AB),
  predicted_relative_AB = mean(predicted_relative_AB),
  biological_repeats = uniqueN(biological_replicate),
  observed_all_repeat_positive = all(observed_extra_pp > 0),
  combo_loss_all_repeats = all(GR_AB < 0),
  neither_single_loss_all_repeats = all(GR_A >= 0 & GR_B >= 0)),
  by = .(time_h, model, anchor, partner, nominal_A_uM, nominal_B_uM)]
write_table(condition, 'held_condition_predictions.tsv.gz')
summary <- condition[, .(dose_contexts = .N,
  observed_extra_median_pp = median(observed_extra_pp), predicted_extra_median_pp = median(predicted_extra_pp),
  extra_effect_r = safe_cor(observed_extra_pp, predicted_extra_pp),
  combination_count_r = safe_cor(relative_AB, predicted_relative_AB),
  observed_effect_at_least10pp = sum(observed_extra_pp >= 10),
  predicted_effect_at_least10pp = sum(predicted_extra_pp >= 10),
  both_effects_at_least10pp = sum(observed_extra_pp >= 10 & predicted_extra_pp >= 10),
  observed_combo_specific_loss_contexts = sum(combo_loss_all_repeats & neither_single_loss_all_repeats)),
  by = .(time_h, model, anchor, partner)]
write_table(summary, 'held_biological_summary.tsv')
inputs <- c(links_path, molecular_path, scores_path,
  '/work/src/program_pilot/61_matched_combination_models.R',
  '/work/src/program_pilot/nonlinear_response.R')
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')),
  'model_source_manifest.tsv')
writeLines(capture.output(sessionInfo()), file.path(root, 'model_sessionInfo.txt'))
status('COMPLETE: predictions frozen, outcomes joined, all partner/anchor results saved')
print(summary[time_h == 3 & model == 'Original40_SVR'])
