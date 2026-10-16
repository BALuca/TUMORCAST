#!/usr/bin/env Rscript
# Verify the export and replay the actual frozen models, without fitting.
entry <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
source(file.path(dirname(normalizePath(entry)), 'config.R'))
check_packages(c('data.table', 'digest'))
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
dir.create(regenerated_root, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(file.path(bundle_root, 'provenance/file_manifest.tsv'))
paths <- file.path(bundle_root, manifest$relative_path)
stopifnot(all(file.exists(paths)))
actual <- vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256')
stopifnot(identical(unname(actual), manifest$sha256))

archive <- readRDS(file.path(input_root, 'molecular_functional_inputs/inputs.rds'))
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread(file.path(input_root, 'training_outcomes/selected_doses.tsv'))
drugs <- bank_selection$drug
contexts <- fread(file.path(result_root, 'single_drug_bridge/linear/data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
source(file.path(bundle_root, 'code/single_drug_prediction_functions.R'))

baseline_saved <- fread(file.path(result_root, 'single_drug_bridge/linear/data/baseline_predictions_frozen.tsv'))
measured_saved <- fread(file.path(result_root, 'single_drug_bridge/linear/data/measured_predictions_frozen.tsv'))
pair_saved <- fread(file.path(result_root, 'single_drug_bridge/linear/data/combination_predictions_frozen.tsv'))
folds <- fread(file.path(result_root, 'single_drug_bridge/linear/data/recipient_folds.tsv'))
checks <- list()
held_displacements <- list()
row_key <- function(x) paste(x$depmap_id, x$drug, sprintf('%.10g', x$dose), sep = '|')
molecular <- as.data.table(archive$prep$m)
stopifnot(!anyDuplicated(row_key(molecular)))

for (job in sort(unique(baseline_saved$fold))) {
  label <- if (job > nrow(folds)) 'final' else sprintf('loo_%02d', job)
  model <- readRDS(file.path(result_root, 'single_drug_bridge/linear/models', paste0(label, '.rds')))
  z <- baseline_saved[fold == job]
  replay <- predict_single(model, z)
  discrepancy <- max(abs(-replay$sensitivity - z$predicted_log2))
  stopifnot(discrepancy < 1e-10)
  checks[[paste0(label, '_baseline')]] <- data.table(
    component = 'single-drug bridge baseline to displacement to function', fit = label,
    rows = nrow(z), maximum_absolute_difference = discrepancy)
  held <- which(z$evaluation == 'held_cell')
  held_displacements[[label]] <- cbind(z[held, .(depmap_id, drug, dose, fold)],
                                      as.data.table(replay$D[held, , drop = FALSE]))

  z <- measured_saved[fold == job]
  if (nrow(z)) {
    ix <- match(row_key(z), row_key(molecular))
    stopifnot(!anyNA(ix))
    replay <- predict_single(model, z, as.matrix(molecular[ix, ..programs]))
    discrepancy <- max(abs(-replay$sensitivity - z$predicted_log2))
    stopifnot(discrepancy < 1e-10)
    checks[[paste0(label, '_measured')]] <- data.table(
      component = 'single-drug bridge measured displacement to function', fit = label,
      rows = nrow(z), maximum_absolute_difference = discrepancy)
  }
  z <- pair_saved[fold == job]
  if (nrow(z)) {
    replay <- predict_pair(model, z)
    discrepancy <- max(abs(-replay$sensitivity - z$predicted_log2),
                       abs(-replay$departure - z$predicted_departure_log2))
    stopifnot(discrepancy < 1e-10)
    checks[[paste0(label, '_pairs')]] <- data.table(
      component = 'single-drug bridge original drug-pair model', fit = label,
      rows = nrow(z), maximum_absolute_difference = discrepancy)
  }
}
predicted_D <- rbindlist(held_displacements)
setnames(predicted_D, tail(names(predicted_D), 40), programs)
fwrite(predicted_D, file.path(regenerated_root, 'held_predicted_displacements.tsv'), sep = '\t')

# Replay the rescue with precisely the frozen molecular transformations used in ligand rescue.
model <- readRDS(file.path(input_root, 'ligand_response_inputs/models/final.rds'))
reference <- readRDS(file.path(result_root, 'single_drug_bridge/linear/models/final.rds'))
stopifnot(identical(model$par$level, reference$par$level),
          identical(model$par$log_amplitude, reference$par$log_amplitude),
          identical(model$context$W, reference$context$W),
          identical(model$rna$lap_receptors, reference$rna$lap_receptors))
source(file.path(bundle_root, 'code/ligand_query_prediction_functions.R'))
keep <- abs(model$template$correlations[, 'Lapatinib ditosylate']) >= .5
stopifnot(sum(keep) == 7)
source(file.path(bundle_root, 'code/bounded_rescue_functions.R'))
rescue <- readRDS(file.path(result_root, 'ligand_rescue/models/final.rds'))
profiles <- fread(file.path(input_root, 'ligand_response_inputs/data/time_matched_HGF_program_inputs.tsv'))
saved <- fread(file.path(result_root, 'ligand_rescue/data/AU565_HGF_predictions_frozen.tsv'))
Ha <- profiles[CellLine == 'AU565']
queries <- list()
for (dose in c(.5, 5)) {
  meta <- data.table(depmap_id = 'ACH-000248', drug = 'Lapatinib ditosylate', dose = dose)
  lap <- predict_single(model, meta)
  ff <- make_features(meta[rep(1, nrow(Ha))], lap$D[rep(1, nrow(Ha)), , drop = FALSE],
                      as.matrix(Ha[, ..programs]))
  rr <- predict_rho(rescue$g, ff$relative_opposition, rep(dose, nrow(Ha)),
                    ff$drug_activity > 0 & ff$ligand_norm > 0)
  cp <- compose(rep(2^(-lap$sensitivity), nrow(Ha)), rr$rho, rr$eta)
  queries[[as.character(dose)]] <- data.table(
    CellLine = 'AU565', dose = dose, source_profile = Ha$source_profile,
    lapatinib_pct_control = 100 * 2^(-lap$sensitivity), rho = rr$rho,
    predicted_combination_pct_control = 100 * cp$abundance_ratio,
    relative_opposition = ff$relative_opposition)
}
query <- rbindlist(queries)
ix <- match(paste(saved$dose, saved$source_profile), paste(query$dose, query$source_profile))
stopifnot(!anyNA(ix))
discrepancy <- max(abs(saved$predicted_combination_pct_control - query$predicted_combination_pct_control[ix]),
                   abs(saved$rho - query$rho[ix]))
stopifnot(discrepancy < 1e-10)
checks[['ligand rescue_HGF']] <- data.table(component = 'ligand rescue frozen HGF query', fit = 'final',
                                 rows = nrow(saved), maximum_absolute_difference = discrepancy)
fwrite(query, file.path(regenerated_root, 'HGF_queries_replayed.tsv'), sep = '\t')

# Reconstruct the held-cell interaction predictions and the end-to-end diagnostic.
outcomes <- fread(file.path(input_root, 'ligand_response_inputs/data/external_combination_outcomes.tsv'))[
  study == 'elife32271']
outcomes[, observation_id := .I]
ligands <- fread(file.path(input_root, 'ligand_response_inputs/data/external_ligand_programs.tsv'))
ix <- match(paste(outcomes$depmap_id, outcomes$gene), paste(ligands$depmap_id, ligands$gene))
stopifnot(!anyNA(ix))
single <- predict_single(model, outcomes[, .(depmap_id, drug, dose)])
ff <- make_features(outcomes, single$D, as.matrix(ligands[ix, ..programs]))
combo_saved <- fread(file.path(result_root, 'ligand_rescue/data/combination_predictions_frozen.tsv'))
combo_replay <- list()
for (i in seq_len(nrow(combo_saved))) {
  row <- combo_saved[i]
  k <- match(row$observation_id, outcomes$observation_id)
  label <- if (row$evaluation == 'held_cell') paste0('held_cell_', row$depmap_id) else 'final'
  rr <- predict_rho(rescue$folds[[label]], ff$relative_opposition[k], outcomes$dose[k],
                    ff$drug_activity[k] > 0 & ff$ligand_norm[k] > 0)
  conditional <- compose(outcomes$observed_drug_pct_vehicle[k] / 100, rr$rho, rr$eta)
  complete <- compose(2^(-single$sensitivity[k] * 3/5), rr$rho, rr$eta)
  vh <- 2^outcomes$observed_ligand_log2[k]
  combo_replay[[i]] <- data.table(
    predicted_departure_log2 = conditional$departure_log2,
    predicted_combination_pct_vehicle = 100 * vh * conditional$abundance_ratio,
    full_bridge_combination_pct_vehicle = 100 * vh * complete$abundance_ratio)
}
replay <- rbindlist(combo_replay)
discrepancy <- max(vapply(names(replay), function(n) max(abs(replay[[n]] - combo_saved[[n]])), numeric(1)))
stopifnot(discrepancy < 1e-10)
checks[['ligand rescue_combinations']] <- data.table(component = 'ligand rescue combination predictions', fit = 'all saved folds',
                                          rows = nrow(combo_saved), maximum_absolute_difference = discrepancy)
fwrite(rbindlist(checks), file.path(regenerated_root, 'prediction_replay_checks.tsv'), sep = '\t')
writeLines(capture.output(sessionInfo()), file.path(regenerated_root, 'replay_session_info.txt'))
cat('All packaged file hashes and frozen prediction replays agree. No models refitted.\n')
