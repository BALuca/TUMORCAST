#!/usr/bin/env Rscript
# Post-primary input ablation: forecast combined cell counts without constituent
# functional measurements. All five designs and both RNA times are fixed below.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(e1071)
  library(digest)
})
source("/work/src/program_pilot/nonlinear_response.R")
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(61012L)
root <- "/work/results/program_pilot/61_individual_responses_combination_transfer"
data_dir <- file.path(root, "data")
svr_grid <- CJ(cost = c(1, 10, 100), gamma_factor = c(.25, 1))
rna_times <- c(3, 24)
model_names <- c("ContextOnly_SVR", "Original40_RNAonly_SVR", "Biology_RNAonly_SVR",
                 "GenePC12_RNAonly_SVR", "Singles_logcounts_SVR")
write_table <- function(x, name) fwrite(x, file.path(data_dir, name), sep = "\t", na = "NA")
status <- function(x) {
  message(format(Sys.time(), tz = "UTC", usetz = TRUE), " | ", x)
  writeLines(x, file.path(root, "ABLATION_STATUS.txt"))
}
# Freeze the post-primary design before loading any combination outcomes.
writeLines(c(
  "Post-primary input ablation, initiated because functional-only controls performed strongly in the main analysis.",
  "This is a retrospective diagnostic, not a preregistered or newly blinded confirmatory analysis.",
  "Fixed models: ContextOnly_SVR; Original40_RNAonly_SVR; Biology_RNAonly_SVR; GenePC12_RNAonly_SVR; Singles_logcounts_SVR.",
  "Target: log2(combined treated cell count / source vehicle count); all actual counts are positive.",
  "Context features: log10 doses, anchor identity, biological-repeat indicators. No partner identity or constituent function enters RNA-only designs.",
  "RNA features: symmetric sums, absolute differences and elementwise products of individual drug scores.",
  "Original40 uses frozen signed measured-gene contributions; Biology uses fixed biological scores; GenePC12 uses the already frozen unsupervised individual-RNA PCA.",
  "Singles_logcounts_SVR adds measured log2 individual-drug relative counts and their symmetric sum, difference and product; it is an oracle-input comparator.",
  "Outer exclusions: each of all six partners, both anchors, every dose and biological replicate. All individual responses remain available.",
  "Inner selection: five whole-partner folds, cost1/10/100 and gamma .25/1 divided by feature count, unchanged established SVR scaling and epsilon.",
  "Primary RNA time3h and secondary24h; all models at both times are reported. No selection among outer results.",
  "Predictions of positive relative combination counts are written before outcome joining. Extra inhibition is calculated only afterward using measured constituent function.",
  "Direct combined-count prediction is essential: posthoc observed/predicted extra effects share the measured constituent-product term and cannot alone establish molecular information.",
  paste0("Frozen before fitting at ", Sys.time())
), file.path(root, "ABLATION_DESIGN.txt"))
qsave(list(models = model_names, rna_times = rna_times, svr_grid = svr_grid, seed = 61012L,
           primary_time_h = 3, target = "log2_ratio_AB", retrospective_post_primary = TRUE),
      file.path(root, "models/ablation_settings.qs"))
links_path <- file.path(data_dir, "matched_molecular_links.tsv.gz")
scores_path <- file.path(data_dir, "representation_scores.qs")
pca_path <- file.path(root, "models/gene_reference_pca.qs")
links <- fread(links_path)
bank <- qread(scores_path)
pca <- qread(pca_path)
score_sets <- list(Original40 = t(bank$original40), Biology = t(bank$biology), GenePC12 = t(pca$x[, 1:12, drop = FALSE]))
stopifnot(nrow(score_sets$Original40) == 40L,
          all(vapply(score_sets, function(x) all(is.finite(x)), logical(1))))

# Reuse the primary model's symmetric representation without functional covariates.
make_response <- function(scores, z) {
  stopifnot(all(z$profile_A %in% colnames(scores)), all(z$profile_B %in% colnames(scores)))
  a <- t(scores[, z$profile_A, drop = FALSE])
  b <- t(scores[, z$profile_B, drop = FALSE])
  sums <- a + b
  difference <- abs(a - b)
  product <- a * b
  colnames(sums) <- paste0("sum_", colnames(a))
  colnames(difference) <- paste0("absdifference_", colnames(a))
  colnames(product) <- paste0("product_", colnames(a))
  cbind(sums, difference, product)
}
balanced_loss <- function(y, pred, md) {
  # Internal selection only; report independent dose curves and assay-scale outcomes.
  error <- (pred - y)^2
  groups <- interaction(md$anchor, md$biological_replicate, drop = TRUE)
  mean(vapply(split(error, groups), mean, numeric(1)))
}
safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 4L || sd(a[ok]) < 1e-10 || sd(b[ok]) < 1e-10) return(NA_real_)
  cor(a[ok], b[ok])
}
predictions <- tuning <- list()
counter <- 0L
for (tt in rna_times) {
  obs <- copy(links[time_h == tt & molecular_supported == TRUE])
  setorder(obs, partner, anchor, biological_replicate, nominal_A_uM, nominal_B_uM)
  stopifnot(all(obs$ratio_AB > 0), all(obs$ratio_A > 0), all(obs$ratio_B > 0), !anyDuplicated(obs$observation_id))
  partners <- sort(unique(obs$partner))
  stopifnot(length(partners) == 6L, uniqueN(obs$anchor) == 2L)
  context <- cbind(logdose_A = log10(obs$nominal_A_uM), logdose_B = log10(obs$nominal_B_uM),
                   model.matrix(~anchor + factor(biological_replicate), obs)[, -1, drop = FALSE])
  log_A <- log2(obs$ratio_A)
  log_B <- log2(obs$ratio_B)
  designs <- list(ContextOnly_SVR = context,
                  Original40_RNAonly_SVR = cbind(context, make_response(score_sets$Original40, obs)),
                  Biology_RNAonly_SVR = cbind(context, make_response(score_sets$Biology, obs)),
                  GenePC12_RNAonly_SVR = cbind(context, make_response(score_sets$GenePC12, obs)),
                  Singles_logcounts_SVR = cbind(context, log_A = log_A, log_B = log_B,
                                                log_sum = log_A + log_B, log_absdifference = abs(log_A - log_B),
                                                log_product = log_A * log_B))
  stopifnot(identical(names(designs), model_names), all(vapply(designs, function(x) all(is.finite(x)), logical(1))))
  # Only this target enters fitting, and only at training indices in each fold.
  y <- log2(obs$ratio_AB)
  for (held in partners) {
    tr <- which(obs$partner != held)
    te <- which(obs$partner == held)
    inner_partners <- partners[partners != held]
    for (model_name in model_names) {
      counter <- counter + 1L
      status(paste(tt, "h | held", held, "|", model_name))
      x <- designs[[model_name]]
      oof <- matrix(NA_real_, length(tr), nrow(svr_grid))
      losses <- matrix(NA_real_, length(inner_partners), nrow(svr_grid))
      for (j in seq_along(inner_partners)) {
        ii <- which(obs$partner[tr] != inner_partners[j])
        vv <- which(obs$partner[tr] == inner_partners[j])
        pp <- vapply(seq_len(nrow(svr_grid)), function(k) {
          fit <- fit_program_svr(x[tr[ii], , drop = FALSE], y[tr[ii]], "radial",
                                  svr_grid$cost[k], svr_grid$gamma_factor[k] / ncol(x))
          predict_program_svr(fit, x[tr[vv], , drop = FALSE])
        }, numeric(length(vv)))
        oof[vv, ] <- pp
        losses[j, ] <- vapply(seq_len(nrow(svr_grid)), function(k) balanced_loss(y[tr[vv]], pp[, k], obs[tr[vv]]), numeric(1))
      }
      stopifnot(all(is.finite(oof)), all(is.finite(losses)))
      chosen <- which.min(colMeans(losses))
      tuning[[counter]] <- data.table(time_h = tt, held_partner = held, model = model_name,
                                      candidate = seq_len(nrow(svr_grid)), cost = svr_grid$cost,
                                      gamma_factor = svr_grid$gamma_factor, internal_grouped_loss = colMeans(losses),
                                      selected = seq_len(nrow(svr_grid)) == chosen)
      fit <- fit_program_svr(x[tr, , drop = FALSE], y[tr], "radial", svr_grid$cost[chosen], svr_grid$gamma_factor[chosen] / ncol(x))
      pp <- predict_program_svr(fit, x)
      key <- paste0(tt, "h_", held, "_", model_name)
      qsave(list(fit = fit, selection = tuning[[counter]][selected == TRUE], columns = colnames(x),
                 training_partners = inner_partners, target = "log2_ratio_AB", input_design = model_name),
            file.path(root, "models", paste0("ablation_", key, ".qs")), preset = "high")
      pred <- rbind(data.table(observation_id = obs$observation_id[te], split = "held_partner", predicted_log2_relative_AB = pp[te]),
                    data.table(observation_id = obs$observation_id[tr], split = "training_fit", predicted_log2_relative_AB = pp[tr]),
                    data.table(observation_id = obs$observation_id[tr], split = "inner_partner_OOF", predicted_log2_relative_AB = oof[, chosen]))
      pred[, `:=`(time_h = tt, held_partner = held, model = model_name,
                    predicted_relative_AB = 2^predicted_log2_relative_AB)]
      stopifnot(all(is.finite(pred$predicted_relative_AB)), all(pred$predicted_relative_AB > 0))
      write_table(pred, paste0("ablation_prediction_", key, ".tsv.gz"))
      predictions[[counter]] <- pred
      write_table(rbindlist(tuning), "ablation_training_selection.tsv")
    }
  }
}
pred <- rbindlist(predictions)
write_table(pred, "ablation_predictions_before_outcome_join.tsv.gz")
writeLines(digest(pred, algo = "sha256"), file.path(root, "models/ablation_prediction_sha256.txt"))

# Evaluation starts only after every design/time/fold prediction is frozen.
observed <- fread(links_path)[molecular_supported == TRUE]
observed[, `:=`(observed_extra_pp = excess_inhibition, relative_A = ratio_A, relative_B = ratio_B,
                relative_AB = ratio_AB, observed_relative_AB = ratio_AB,
                dose_anchor = dose_A_uM, dose_partner = dose_B_uM)]
pred <- merge(pred, observed, by = c("observation_id", "time_h"), all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(pred$relative_AB), all(is.finite(pred$predicted_relative_AB)))
pred[, `:=`(predicted_extra_pp = 100 * (relative_A * relative_B - predicted_relative_AB),
             predicted_GR_AB = 2^(log2(predicted_relative_AB * ncontrol72_AB / n0_AB) / log2(ncontrol72_AB / n0_AB)) - 1,
             impossible_negative_cell_count = FALSE)]
write_table(pred, "ablation_predictions.tsv.gz")
condition <- pred[split == "held_partner", .(observed_relative_AB = mean(relative_AB),
                                             predicted_relative_AB = mean(predicted_relative_AB),
                                             observed_extra_pp = mean(observed_extra_pp), predicted_extra_pp = mean(predicted_extra_pp),
                                             observed_GR_AB = mean(GR_AB), predicted_GR_AB = mean(predicted_GR_AB),
                                             biological_repeats = uniqueN(biological_replicate),
                                             observed_net_loss_all_repeats = all(GR_AB < 0),
                                             predicted_net_loss_all_repeats = all(predicted_GR_AB < 0)),
                  by = .(time_h, model, anchor, partner, condition_id, nominal_A_uM, nominal_B_uM)]
write_table(condition, "ablation_held_condition_predictions.tsv.gz")
summary <- condition[, .(dose_contexts = .N, combination_count_r = safe_cor(observed_relative_AB, predicted_relative_AB),
                         extra_effect_r = safe_cor(observed_extra_pp, predicted_extra_pp),
                         observed_count_median_pct = 100 * median(observed_relative_AB),
                         predicted_count_median_pct = 100 * median(predicted_relative_AB),
                         observed_extra_median_pp = median(observed_extra_pp), predicted_extra_median_pp = median(predicted_extra_pp),
                         observed_effect_at_least10pp = sum(observed_extra_pp >= 10),
                         predicted_effect_at_least10pp = sum(predicted_extra_pp >= 10),
                         both_effects_at_least10pp = sum(observed_extra_pp >= 10 & predicted_extra_pp >= 10),
                         observed_net_loss_all_repeats = sum(observed_net_loss_all_repeats),
                         predicted_net_loss_all_repeats = sum(predicted_net_loss_all_repeats)),
                     by = .(time_h, model, anchor, partner)]
write_table(summary, "ablation_held_biological_summary.tsv")
overall <- condition[, .(dose_contexts = .N,
                         combination_count_r = safe_cor(observed_relative_AB, predicted_relative_AB),
                         extra_effect_r = safe_cor(observed_extra_pp, predicted_extra_pp),
                         observed_extra_median_pp = median(observed_extra_pp), predicted_extra_median_pp = median(predicted_extra_pp),
                         observed_count_median_pct = 100 * median(observed_relative_AB), predicted_count_median_pct = 100 * median(predicted_relative_AB)),
                     by = .(time_h, model)]
write_table(overall, "ablation_overall_summary.tsv")
writeLines(c(
  "Input-ablation interpretation: all fixed models and both RNA times are retained.",
  "Read direct held-combination cell-count prediction first. An RNA contribution requires evidence beyond the same context-only comparator.",
  "The context-only comparator contains doses, anchor identity and experiment indicators; it contains neither individual drug function nor RNA.",
  "Posthoc extra inhibition subtracts predicted or observed combination counts from the same measured constituent-product reference.",
  "Consequently a high extra-effect correlation can coexist with poor direct combined-count agreement. It must not be interpreted alone as a molecular signal.",
  "This ablation was defined after inspection of the primary controls. Its results are diagnostic and retrospective, not preregistered confirmation.",
  "All evaluations are within BT20 and use six held partners with two previously observed anchors; 198 dose contexts do not imply 198 independent interventions.",
  capture.output(print(overall))
), file.path(root, "ABLATION_INTERPRETATION.txt"))
inputs <- c(links_path, scores_path, pca_path, "/work/src/program_pilot/61_molecular_only_ablation.R",
            "/work/src/program_pilot/nonlinear_response.R")
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = "sha256")),
            "ablation_source_manifest.tsv")
writeLines(capture.output(sessionInfo()), file.path(root, "ablation_sessionInfo.txt"))
status("COMPLETE: all fixed input-ablation designs predicted and evaluated")
print(overall)
