#!/usr/bin/env Rscript
# Audit the exact B1 regression and probe HGF using saved fits only.
# The summed-response probe adds an explicit assumption; it is not a fitted pair model.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(digest)
  library(ggplot2)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

# Configuration: fixed source archive, recipient, exposures and all HGF contexts.
workspace <- '/work'
root <- file.path(workspace, 'results/program_pilot/74_v20_lapatinib_hgf_frozen_query')
archive <- file.path(workspace, 'almanac_combination_training_2026-10-03/TUMORCAST_lapatinib_ablation_2026-10-04.zip')
query_path <- file.path(workspace, 'grant_rebuild_2026-10-04/queries/AU565_component_queries.tsv')
bank_path <- file.path(workspace, 'almanac_conditionals/state_conditioned_bridge/program_first/bridge/prepared.rds')
recipient <- 'ACH-000248'
doses <- c(1, 2.5)
programs <- sprintf('P%02d', 1:40)
stopifnot(file.exists(archive), file.exists(query_path), file.exists(bank_path))
for (folder in c('inputs', 'data', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')

# Preserve the original version-specific implementation and model objects.
unzip(archive, exdir = file.path(root, 'inputs'))
original <- file.path(root, 'inputs/original')
saved_results <- file.path(root, 'inputs/results')
setwd(original)
source('modular_pair_pilot/model.R')
stopifnot(identical(colnames(PREP$W), programs), recipient %in% PREP$raw$depmap_id,
  !recipient %in% SINGLE$depmap_id)
baseline <- getmatrix(PREP$raw, recipient, programs)
reference_bank <- readRDS(bank_path)
weight_difference <- max(abs(PREP$W - reference_bank$W[rownames(PREP$W), programs]))
baseline_difference <- max(abs(baseline -
  as.matrix(reference_bank$raw[depmap_id == recipient, programs, with = FALSE])))
common <- match(PREP$m$sample_id, reference_bank$m$sample_id)
keep <- which(!is.na(common))
stopifnot(length(keep) == 252L)
molecular_difference <- max(abs(as.matrix(PREP$m[keep, programs]) -
  as.matrix(reference_bank$m[common[keep], programs, with = FALSE])))
stopifnot(weight_difference < 1e-12, baseline_difference < 1e-12, molecular_difference < 1e-8)

# Exact replay uses each held recipient's saved fold and original upstream operator.
saved <- as.data.table(read_tsv(file.path(saved_results, 'predictions.tsv'), show_col_types = FALSE))
models <- lapply(1:3, function(k) readRDS(file.path(saved_results, paste0('fold_', k, '.rds')))$models$full)
replay <- vector('list', 3L)
for (k in 1:3) {
  model <- models[[k]]
  rows <- saved[model == 'full' & fold == k]
  stopifnot(!recipient %in% model$training_single_ids,
    length(model$training_pair_ids) == 0L, all(model$par$Q == 0), model$cfg$family == 'null',
    !any(rows$depmap_id %in% model$training_single_ids))
  x <- uppred(PREP, model$operators[['Lapatinib ditosylate']], rows$depmap_id, rows$dose)
  s <- getmatrix(PREP$raw, rows$depmap_id, programs)
  replayed <- predict_single(model, x, s, 'secondary')
  difference <- max(abs(replayed - rows$predicted))
  stopifnot(difference < 1e-10)
  replay[[k]] <- data.table(fold = k, rows = nrow(rows), maximum_absolute_difference = difference)
}
write_table(rbindlist(replay), 'held_prediction_replay')
metrics <- saved[dose == 2.5 & model %in% c('full', 'refit_no_gates', 'refit_mean_displacement', 'direct_baseline'),
  .(n = .N, pearson_log2 = cor(observed, predicted),
    sensitive_observed = sum(observed <= -1), resistant_observed = sum(observed > -1),
    sensitive_called = sum(observed <= -1 & predicted <= -1),
    false_sensitive = sum(observed > -1 & predicted <= -1)), by = model]
stopifnot(all(metrics$n == 441L),
  abs(metrics[model == 'full', pearson_log2] - .5089411777536913) < 1e-12,
  abs(metrics[model == 'refit_no_gates', pearson_log2] - .2666291495842529) < 1e-12)
write_table(metrics, 'exact_b1_correlations')

# These HGF vectors already include the original AU565 baseline transport to
# corrected treated-minus-reference log2(1+CP10k), then the same frozen W.
# No direct projection of raw CytoSig logFC or HGF-based calibration is performed.
queries <- fread(query_path)
hgf <- queries[query == 'HGF' & type == 'ligand']
stopifnot(nrow(hgf) == 12L, !anyDuplicated(hgf$context), all(is.na(hgf$dose)),
  all(is.finite(as.matrix(hgf[, programs, with = FALSE]))))
mean_hgf <- data.table(context = 'Equal-context mean')
mean_hgf[, (programs) := as.list(colMeans(as.matrix(hgf[, programs, with = FALSE])))]
hgf <- rbindlist(list(hgf[, c('context', programs), with = FALSE], mean_hgf), use.names = TRUE)
write_table(hgf, 'hgf_transported_inputs')

# Query both possible calculations while keeping their meanings explicit.
# Archived rule: separate f(L)+f(H), because g=0. Molecular probe: f(L+H).
predictions <- vector('list', length(models) * length(doses) * nrow(hgf))
support <- vector('list', length(predictions))
vectors <- vector('list', length(models) * length(doses))
counter <- 0L
for (k in seq_along(models)) {
  model <- models[[k]]
  for (dose_value in doses) {
    l <- uppred(PREP, model$operators[['Lapatinib ditosylate']], recipient, dose_value)
    vectors[[(k - 1L) * length(doses) + match(dose_value, doses)]] <-
      cbind(data.table(fold = k, dose = dose_value), as.data.table(l))
    f0 <- predict_single(model, matrix(0, 1, 40), baseline, 'secondary')
    fl <- predict_single(model, l, baseline, 'secondary')
    stopifnot(f0 == 0)
    # Training-range flags describe functional input extrapolation, not calibrated uncertainty.
    training_rows <- SINGLE %>% filter(depmap_id %in% model$training_single_ids,
      assay != 'secondary' | drug %in% ANCHORS)
    tx <- base_predict(PREP, model$operators, training_rows)
    minima <- apply(tx, 2L, min)
    maxima <- apply(tx, 2L, max)
    for (i in seq_len(nrow(hgf))) {
      counter <- counter + 1L
      h <- as.matrix(hgf[i, programs, with = FALSE])
      fh <- predict_single(model, h, baseline, 'secondary')
      fs <- predict_single(model, l + h, baseline, 'secondary')
      # Exact regression nonlinearity supplies any interaction in the sum probe.
      interaction <- fs - fl - fh + f0
      base_rows <- data.table(fold = k, dose = dose_value, hgf_context = hgf$context[i],
        rule = rep(c('Archived separate effects', 'Summed molecular responses: extrapolation'), each = 4L),
        condition = rep(c('Untreated', 'HGF', 'Lapatinib', 'Lapatinib + HGF'), 2L),
        log2_abundance = c(f0, fh, fl, fl + fh, f0, fh, fl, fs))
      base_rows[, abundance_pct := 100 * 2^log2_abundance]
      base_rows[, drug_effect_in_hgf_log2 := rep(c(fl - f0, fs - fh), each = 4L)]
      base_rows[, normalized_rescue_log2 := rep(c(0, interaction), each = 4L)]
      base_rows[, drug_effect_pct_of_environment := 100 * 2^drug_effect_in_hgf_log2]
      stopifnot(all(is.finite(base_rows$abundance_pct)))
      predictions[[counter]] <- base_rows
      support[[counter]] <- data.table(fold = k, dose = dose_value, hgf_context = hgf$context[i],
        input = c('Lapatinib', 'HGF', 'Lapatinib + HGF'),
        reference = 'Deployed operators on functional-training rows; not cross-fitted descriptors',
        coordinates_outside_training_range = c(sum(l < minima | l > maxima),
          sum(h < minima | h > maxima), sum(l + h < minima | l + h > maxima)),
        normalized_rescue_log2 = interaction)
    }
  }
}
prediction <- rbindlist(predictions)
stopifnot(nrow(prediction) == 3L * 2L * 13L * 8L,
  !anyDuplicated(prediction[, .(fold, dose, hgf_context, rule, condition)]))
write_table(prediction, 'au565_all_queries')
write_table(rbindlist(support), 'functional_input_support')
write_table(rbindlist(vectors), 'au565_lapatinib_displacements')
write_table(prediction[hgf_context == 'Equal-context mean'], 'au565_mean_hgf_queries')
summary <- prediction[hgf_context != 'Equal-context mean' & condition == 'Lapatinib + HGF',
  .(fold_context_queries = .N, positive_normalized_rescue = sum(normalized_rescue_log2 > 1e-10),
    median_rescue_log2 = median(normalized_rescue_log2), minimum_rescue_log2 = min(normalized_rescue_log2),
    maximum_rescue_log2 = max(normalized_rescue_log2)), by = .(rule, dose)]
write_table(summary, 'hgf_context_summary')

# A compact diagnostic shows the exact saved folds, not a selected rescue curve.
plot_rows <- prediction[hgf_context == 'Equal-context mean' & dose == 1]
plot_rows[, condition := factor(condition, levels = c('Untreated', 'HGF', 'Lapatinib', 'Lapatinib + HGF'),
  labels = c('Untreated', 'HGF', 'Lapatinib', 'Lapatinib\n+ HGF'))]
plot_rows[, rule := factor(rule, levels = c('Archived separate effects', 'Summed molecular responses: extrapolation'),
  labels = c('Archived rule: separate effects', 'Probe: summed molecular responses'))]
p <- ggplot(plot_rows, aes(condition, abundance_pct, group = factor(fold), color = factor(fold))) +
  geom_hline(yintercept = 100, linetype = 2, color = '#999999', linewidth = .35) +
  geom_line(linewidth = .55) + geom_point(size = 2) + facet_wrap(~rule, nrow = 1) +
  scale_color_manual(values = c('#1B6CA8', '#A62738', '#477A41'), name = 'Saved fold') +
  labs(x = NULL, y = 'Predicted abundance (% untreated)',
    title = 'Exact v20 B1 fit queried for AU565',
    subtitle = 'Lapatinib 1 µM; equal mean of all 12 transported HGF source contexts',
    caption = 'Frozen model calculations. No matched combination validation.\nFold spread is refit sensitivity; summed responses are an added molecular assumption.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom',
    plot.caption = element_text(hjust = 0, size = 8), strip.background = element_blank(),
    strip.text = element_text(face = 'bold', size = 9))
ggsave(file.path(root, 'figures/au565_frozen_hgf_query.pdf'), p, width = 7.4, height = 4.6, device = cairo_pdf)
ggsave(file.path(root, 'figures/au565_frozen_hgf_query.png'), p, width = 7.4, height = 4.6, dpi = 160)

# The B1-matched dose is an equally explicit second view of the same frozen query.
plot_rows_25 <- prediction[hgf_context == 'Equal-context mean' & dose == 2.5]
plot_rows_25[, condition := factor(condition, levels = c('Untreated', 'HGF', 'Lapatinib', 'Lapatinib + HGF'),
  labels = c('Untreated', 'HGF', 'Lapatinib', 'Lapatinib\n+ HGF'))]
plot_rows_25[, rule := factor(rule, levels = c('Archived separate effects', 'Summed molecular responses: extrapolation'),
  labels = c('Archived rule: separate effects', 'Probe: summed molecular responses'))]
p25 <- p %+% plot_rows_25 +
  labs(subtitle = 'Lapatinib 2.5 µM; equal mean of all 12 transported HGF source contexts')
ggsave(file.path(root, 'figures/au565_b1_dose_hgf_query.pdf'), p25, width = 7.4, height = 4.6, device = cairo_pdf)
ggsave(file.path(root, 'figures/au565_b1_dose_hgf_query.png'), p25, width = 7.4, height = 4.6, dpi = 160)

sources <- c(archive, query_path, bank_path,
  file.path(original, 'modular_pair_pilot/model.R'),
  file.path(original, 'modular_pair_pilot/bank_operator.R'),
  file.path(workspace, 'src/program_pilot/74_v20_lapatinib_hgf_frozen_query.R'))
write_table(data.table(path = sources, sha256 = vapply(sources, digest, character(1), file = TRUE, algo = 'sha256')),
  'source_hashes')
output_paths <- file.path(root, 'data', c('au565_all_queries.tsv', 'hgf_transported_inputs.tsv',
  'au565_lapatinib_displacements.tsv'))
write_table(data.table(path = output_paths,
  sha256 = vapply(output_paths, digest, character(1), file = TRUE, algo = 'sha256')), 'query_hashes')
write_table(data.table(check = c('Weights match original40', 'AU565 baseline matches original40',
  'Common molecular scores match original40'), maximum_difference = c(weight_difference,
    baseline_difference, molecular_difference)), 'coordinate_checks')
writeLines(capture.output(sessionInfo()), file.path(root, 'session_info.txt'))
writeLines(c('Exact v20 B1 model; saved inference only. No refitting or outcome-based selection.',
  'Original rule: f(L,s)+f(H,s), with g=0; normalized HGF rescue is exactly zero.',
  'Alternative query f(L+H,s) assumes additive molecular responses and is unvalidated.',
  'Intrinsic baseline s is fixed across conditions. HGF transport is not measured AU565 RNA.',
  'No HGF-treated 441-cell functional cohort is present; no combination Pearson r can be computed.',
  'Source vectors: grant_rebuild_2026-10-04/queries/EXTERNAL_PROTOCOL.txt and',
  'almanac_combination_training_2026-10-03/SOURCES.txt document transfer and resource terms.'),
  file.path(root, 'REPORT.txt'))
print(metrics)
print(prediction[hgf_context == 'Equal-context mean' & condition == 'Lapatinib + HGF',
  .(fold, dose, rule, abundance_pct, drug_effect_pct_of_environment, normalized_rescue_log2)])
print(summary)
