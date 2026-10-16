#!/usr/bin/env Rscript
# Audit the existing explicit interaction fits. No model is fitted or selected.
suppressPackageStartupMessages({library(data.table); library(qs)})
root <- '/work/results/program_pilot/67_full_shrunken_program_interactions'
prior <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data'
RhpcBLASctl::blas_set_num_threads(1L)
setDTthreads(2L)
obs <- fread(file.path(prior, 'matched_molecular_links.tsv.gz'))[time_h == 3 & molecular_supported == TRUE]
setorder(obs, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)
bank <- qread(file.path(prior, 'representation_scores.qs'))
settings <- qread(file.path(root, 'models/settings.qs'))
response <- sweep(bank$original40, 2, settings$rms, '/')
a <- response[obs$profile_A, , drop = FALSE]
b <- response[obs$profile_B, , drop = FALSE]
ij <- which(upper.tri(matrix(0, 40, 40), diag = TRUE), arr.ind = TRUE)
x <- vapply(seq_len(nrow(ij)), function(k) {
  i <- ij[k, 1]
  j <- ij[k, 2]
  (a[, i] * b[, j] + a[, j] * b[, i]) / 2
}, numeric(nrow(a)))
colnames(x) <- paste(colnames(a)[ij[, 1]], colnames(a)[ij[, 2]], sep = '__')
checks <- list()
for (held in sort(unique(obs$partner))) {
  fit <- qread(file.path(root, 'models', paste0(held, '_Full.qs')))
  stopifnot(identical(colnames(x), fit$columns), !held %in% fit$training_partners)
  tr <- which(obs$partner %in% fit$training_partners)
  te <- which(obs$partner == held)
  scaled <- sweep(x[, fit$keep, drop = FALSE], 2, fit$sx, '/')
  # Training rows constrain only their row space. Held feature components
  # orthogonal to it cannot be identified from these training outcomes.
  decomp <- svd(scaled[tr, , drop = FALSE], nu = 0)
  keep <- decomp$d > max(decomp$d) * 1e-8
  v <- decomp$v[, keep, drop = FALSE]
  supported <- scaled[te, , drop = FALSE] %*% v %*% t(v)
  unsupported <- scaled[te, , drop = FALSE] - supported
  beta <- fit$raw_coefficients[fit$keep] * fit$sx
  total_prediction <- as.numeric(scaled[te, , drop = FALSE] %*% beta)
  unsupported_prediction <- as.numeric(unsupported %*% beta)
  out <- obs[te, .(condition_id, anchor, partner)]
  out[, `:=`(outside_training_span_fraction = sqrt(rowSums(unsupported^2) / rowSums(scaled[te, , drop = FALSE]^2)),
    predicted_g = total_prediction, unsupported_g = unsupported_prediction,
    training_rank = sum(keep), candidate_terms = ncol(scaled),
    effective_dimensions_95pct = which(cumsum(decomp$d^2) / sum(decomp$d^2) >= .95)[1])]
  checks[[held]] <- out[, .(outside_training_span_fraction = mean(outside_training_span_fraction),
    predicted_g = mean(predicted_g), unsupported_g = mean(unsupported_g),
    training_rank = unique(training_rank), candidate_terms = unique(candidate_terms),
    effective_dimensions_95pct = unique(effective_dimensions_95pct)), by = .(condition_id, anchor, partner)]
}
span <- rbindlist(checks)
stopifnot(nrow(span) == 198L, !anyDuplicated(span$condition_id))
fwrite(span, file.path(root, 'data/interaction_feature_support.tsv.gz'), sep = '\t')
summary <- span[, .(conditions = .N, training_rank = unique(training_rank),
  dimensions_95pct = unique(effective_dimensions_95pct),
  median_outside_span_fraction = median(outside_training_span_fraction),
  rms_unsupported_prediction = sqrt(mean(unsupported_g^2)),
  rms_prediction = sqrt(mean(predicted_g^2))), by = partner]
fwrite(summary, file.path(root, 'data/interaction_feature_support_summary.tsv.gz'), sep = '\t')
print(summary)
co <- fread(file.path(root, 'data/coefficients.tsv.gz'))[model == 'Full']
recurrence <- co[, .(selected_folds = sum(abs(coefficient) > 1e-10),
  positive_folds = sum(coefficient > 1e-10), negative_folds = sum(coefficient < -1e-10)), by = term]
fwrite(recurrence, file.path(root, 'data/interaction_term_recurrence.tsv.gz'), sep = '\t')
# Biological repeats assess reproducibility of the observed departure itself.
held_predictions <- fread(file.path(root, 'data/predictions.tsv.gz'))[model == 'Full' & split == 'held']
held_predictions[, observed_g := log2(ratio_AB / (ratio_A * ratio_B))]
replicates <- dcast(held_predictions, condition_id + anchor + partner ~ biological_replicate, value.var = 'observed_g')
replicate_columns <- setdiff(names(replicates), c('condition_id', 'anchor', 'partner'))
repeat_results <- lapply(combn(replicate_columns, 2, simplify = FALSE), function(pair) {
  available <- replicates[is.finite(get(pair[1])) & is.finite(get(pair[2]))]
  data.table(repeat_a = pair[1], repeat_b = pair[2], conditions = nrow(available),
    observed_departure_r = cor(available[[pair[1]]], available[[pair[2]]]))
})
repeat_results <- rbindlist(repeat_results)
fwrite(repeat_results, file.path(root, 'data/observed_interaction_repeat_agreement.tsv.gz'), sep = '\t')
print(repeat_results)
