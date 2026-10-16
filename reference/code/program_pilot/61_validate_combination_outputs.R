#!/usr/bin/env Rscript
# Validate the completed workflow against its scientific contracts. This script
# replays frozen models and scores; it neither refits nor selects a model.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(digest)
})
source('/work/src/program_pilot/nonlinear_response.R')
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
data_dir <- file.path(root, 'data')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
links <- fread(file.path(data_dir, 'matched_molecular_links.tsv.gz'))
scores <- qread(file.path(data_dir, 'representation_scores.qs'))
molecular <- qread(file.path(data_dir, 'matched_molecular.qs'))
weights <- qread(file.path(root, 'models/original40_weights.qs'))
primary <- fread(file.path(data_dir, 'predictions.tsv.gz'))
ablation <- fread(file.path(data_dir, 'ablation_predictions.tsv.gz'))
checks <- list()

# A failed numerical or exclusion check stops execution before a success report.
record <- function(name, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = passed, evidence = evidence)
  if (!isTRUE(passed)) stop(name, ': ', evidence)
}
supported <- links[molecular_supported == TRUE]
record('Matched observations are unique within RNA time',
  !anyDuplicated(supported[, .(observation_id, time_h)]),
  paste(nrow(supported), 'observation-time rows'))
record('All target individuals have measured RNA profiles',
  all(c(supported$profile_A, supported$profile_B) %in% rownames(scores$original40)),
  'Both individual interventions remain available in each combination exclusion')
record('Dose contexts have actual biological repeats',
  all(supported[, uniqueN(biological_replicate), by = .(time_h, condition_id)]$V1 >= 2L),
  'At least two source biological experiments for each supported nominal context')

# Recompute original signed scores from measured members without renormalization.
genes <- intersect(rownames(molecular$displacement), rownames(weights))
reprojected <- as.matrix(t(molecular$displacement[genes, , drop = FALSE]) %*% weights[genes, , drop = FALSE])
reference <- scores$original40[rownames(reprojected), colnames(reprojected), drop = FALSE]
score_difference <- max(abs(reprojected - reference))
record('Original40 scores replay from frozen signed weights', score_difference < 1e-10,
  paste('Maximum score discrepancy:', format(score_difference, digits = 5)))

# Check the endpoint definition against the saved source-normalized counts.
extra_difference <- max(abs(supported$excess_inhibition -
  100 * (supported$ratio_A * supported$ratio_B - supported$ratio_AB)))
record('Extra effect uses the recorded individual and combined counts', extra_difference < 1e-8,
  paste('Maximum discrepancy in percentage points:', format(extra_difference, digits = 5)))
wells <- fread(file.path(data_dir, 'matched_functional_wells.tsv.gz'))
recomputed_gr <- 2^(log2(wells$n72 / wells$n0) / log2(wells$ncontrol72 / wells$n0)) - 1
gr_difference <- max(abs(wells$recomputed_GR - recomputed_gr))
record('Each physical well GR uses its initial and vehicle counts', gr_difference < 1e-8,
  paste('Maximum GR discrepancy:', format(gr_difference, digits = 5)))
# GR is nonlinear: the reported biological-observation value is the mean of
# physical-well GRs, not GR calculated after averaging the normalization counts.
well_gr <- wells[arm == 'combination', .(source_mean_GR = mean(recomputed_GR)),
  by = .(condition_id = matched_condition_id, biological_replicate)]
gr_check <- merge(supported, well_gr, by = c('condition_id', 'biological_replicate'), all.x = TRUE)
gr_mean_difference <- max(abs(gr_check$GR_AB - gr_check$source_mean_GR))
record('Biological-observation GR preserves the mean of source well GRs',
  !anyNA(gr_check$source_mean_GR) && gr_mean_difference < 1e-8,
  paste('Maximum aggregated GR discrepancy:', format(gr_mean_difference, digits = 5)))

# All held predictions must exclude the entire partner, at both anchors and all doses.
for (kind in c('primary', 'ablation')) {
  p <- if (kind == 'primary') primary else ablation
  held <- p[split == 'held_partner']
  record(paste(kind, 'held outcomes are exactly the excluded partner'),
    all(held$partner == held$held_partner), paste(nrow(held), 'held predictions'))
  record(paste(kind, 'training rows exclude the target partner'),
    all(p[split %chin% c('training_fit', 'inner_partner_OOF'), partner != held_partner]),
    'Checked the stored partner of every training and inner-validation prediction')
  sizes <- held[, .(observations = .N, contexts = uniqueN(condition_id), partners = uniqueN(partner)),
    by = .(time_h, model)]
  record(paste(kind, 'each design reports all supported contexts'),
    all(sizes$observations == 496L & sizes$contexts == 198L & sizes$partners == 6L),
    '496 biological-repeat observations; 198 dose contexts; six excluded partners per design/time')
}

# Replay both direct-count designs for every outer fold. Checking the input
# column names also prevents constituent functional effects entering RNA-only fits.
replay <- list()
for (tt in c(3, 24)) {
  obs <- copy(supported[time_h == tt])
  setorder(obs, partner, anchor, biological_replicate, nominal_A_uM, nominal_B_uM)
  context <- cbind(logdose_A = log10(obs$nominal_A_uM), logdose_B = log10(obs$nominal_B_uM),
    model.matrix(~anchor + factor(biological_replicate), obs)[, -1, drop = FALSE])
  a <- scores$original40[obs$profile_A, , drop = FALSE]
  b <- scores$original40[obs$profile_B, , drop = FALSE]
  sums <- a + b
  difference <- abs(a - b)
  product <- a * b
  colnames(sums) <- paste0('sum_', colnames(a))
  colnames(difference) <- paste0('absdifference_', colnames(a))
  colnames(product) <- paste0('product_', colnames(a))
  designs <- list(ContextOnly_SVR = context,
    Original40_RNAonly_SVR = cbind(context, sums, difference, product))
  for (partner_name in sort(unique(obs$partner))) {
    for (model_name in names(designs)) {
      path <- file.path(root, 'models', paste0('ablation_', tt, 'h_', partner_name, '_', model_name, '.qs'))
      model <- qread(path)
      stopifnot(!partner_name %in% model$training_partners,
        setequal(model$training_partners, setdiff(obs$partner, partner_name)),
        identical(model$columns, colnames(designs[[model_name]])))
      te <- which(obs$partner == partner_name)
      pp <- 2^predict_program_svr(model$fit, designs[[model_name]][te, , drop = FALSE])
      saved <- ablation[time_h == tt & model == model_name & held_partner == partner_name & split == 'held_partner']
      saved <- saved[match(obs$observation_id[te], observation_id)]
      stopifnot(!anyNA(saved$observation_id))
      replay[[length(replay) + 1L]] <- data.table(time_h = tt, partner = partner_name,
        model = model_name, observations = length(te), max_count_fraction_discrepancy = max(abs(pp - saved$predicted_relative_AB)))
    }
  }
}
replay <- rbindlist(replay)
record('All24 original40/context direct-count model replays match',
  all(replay$max_count_fraction_discrepancy < 1e-8),
  paste('Largest discrepancy:', format(max(replay$max_count_fraction_discrepancy), digits = 5)))
fwrite(replay, file.path(data_dir, 'validation_model_replay.tsv'), sep = '\t')

# Invalid predictions in the primary extra-effect regression remain visible;
# direct log-count models are positive by construction and are never clipped.
physical <- rbindlist(list(
  primary[split == 'held_partner', .(workflow = 'Primary extra-effect model', observations = .N,
    nonpositive_predictions = sum(predicted_relative_AB <= 0)), by = .(time_h, model)],
  ablation[split == 'held_partner', .(workflow = 'Direct log-count input ablation', observations = .N,
    nonpositive_predictions = sum(predicted_relative_AB <= 0)), by = .(time_h, model)]))
fwrite(physical, file.path(data_dir, 'validation_physical_prediction_counts.tsv'), sep = '\t')
record('Direct-count models are finite and physically positive',
  all(is.finite(ablation$predicted_relative_AB) & ablation$predicted_relative_AB > 0),
  'Primary regression negative predictions are retained separately, never silently clipped')
fwrite(rbindlist(checks), file.path(data_dir, 'validation_checks.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'validation_sessionInfo.txt'))
writeLines(c('Completed source/count, original-score, exclusion and frozen-model replay checks.',
  'These checks establish implementation fidelity, not biological predictive validity.',
  paste('Completed UTC:', format(Sys.time(), tz = 'UTC'))), file.path(root, 'VALIDATION.txt'))
message('Completed validation of workflow61 outputs.')
