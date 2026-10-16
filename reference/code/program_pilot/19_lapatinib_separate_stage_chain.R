#!/usr/bin/env Rscript
# Connect independently fitted baseline->displacement and measured-displacement->function models.
# Reuse script 18's completed hard-score fits without any refitting or model selection.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(ggplot2)})
data.table::setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(2L)
RhpcBLASctl::omp_set_num_threads(2L)
source_dir <- '/work/results/program_pilot/18_lapatinib_chain_diagnostics'
out <- '/work/results/program_pilot/19_lapatinib_separate_stage_chain'
for (sub in c('', 'data', 'models', 'figures')) dir.create(file.path(out, sub), recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(status = 'running', completed_folds = 0L), file.path(out, 'run_status.tsv'), sep = '\t')
model_paths <- file.path(source_dir, 'models', sprintf('diagnostic_fold_%02d.qs', 1:5))
models <- lapply(model_paths, qread)
stopifnot(all(vapply(models, function(m) identical(m$scoring, 'hard_cluster_mean'), logical(1))))
labels_path <- '/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge/joint40/prism_screen_audit/preferred_screen_functional_labels.tsv'
expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
labels <- fread(labels_path)[drug == 'Lapatinib ditosylate' & dose == 2.5,
  .(depmap_id, fold, observed = corrected_y)]
genes <- sort(unique(c(unlist(lapply(models, function(m) rownames(m$weights))), 'ERBB2', 'MET')))
headers <- names(fread(expression_path, nrows = 0))
symbols <- sub(' \\(ENSG.*', '', headers[-1L])
ix <- which(symbols %chin% genes)
tab <- fread(expression_path, select = c(1L, ix + 1L))
expression <- as.matrix(tab[, -1L])
rownames(expression) <- tab[[1L]]
colnames(expression) <- symbols[ix]
stopifnot(!anyDuplicated(colnames(expression)), all(genes %chin% colnames(expression)))
paths <- c(model_paths, labels_path, expression_path, '/work/src/program_pilot/19_lapatinib_separate_stage_chain.R')
hashes <- data.table(path = paths, sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(hashes, file.path(out, 'data/input_hashes.tsv'), sep = '\t')

# Apply the saved glmnet model and its stored response transformation exactly.
predict_saved <- function(fit, x) {
  raw <- predict(fit$fit, newx = x, s = fit$lambda)
  if (!fit$multi) return(as.numeric(raw))
  z <- matrix(raw, nrow(x), length(fit$center))
  sweep(sweep(z, 2L, fit$scale, '*'), 2L, fit$center, '+')
}
design <- function(d, b, kind) {
  if (kind == 'linear') return(d)
  if (kind == 'quadratic') return(cbind(d, d^2))
  stopifnot(kind == 'conditioned')
  cbind(d, d^2, d * b)
}
rows <- checks <- list()
for (k in 1:5) {
  m <- models[[k]]
  training_ids <- labels[fold != k, depmap_id]
  testing <- labels[fold == k]
  baseline <- expression[, rownames(m$weights), drop = FALSE] %*% m$weights
  bmu <- colMeans(baseline[training_ids, , drop = FALSE])
  bsd <- apply(baseline[training_ids, , drop = FALSE], 2L, sd)
  b <- sweep(sweep(baseline, 2L, bmu, '-'), 2L, bsd, '/')
  x <- list(receptors = expression[, c('ERBB2', 'MET'), drop = FALSE],
    baseline_programs = baseline,
    baseline_programs_receptors = cbind(baseline, expression[, c('ERBB2', 'MET'), drop = FALSE]))
  measured <- m$measured_displacements
  paired <- testing[depmap_id %chin% rownames(measured)]
  stage_links <- list()
  for (kind in c('linear', 'quadratic', 'conditioned')) {
    decoder <- m$models[[paste0('measured_', kind)]]
    stopifnot(!length(intersect(decoder$training_ids, m$excluded_ids)),
      all(decoder$training_ids %chin% rownames(measured)))
    # Measured-input prediction isolates the functional link on the identical paired test cells.
    observed_input <- predict_saved(decoder, design(measured[paired$depmap_id, , drop = FALSE],
      b[paired$depmap_id, , drop = FALSE], kind))
    rows[[length(rows) + 1L]] <- data.table(fold = k, cohort = 'paired', decoder = kind,
      response_input = 'measured', depmap_id = paired$depmap_id, observed = paired$observed, predicted = observed_input)
    for (upstream_name in names(x)) {
      upstream <- m$models[[paste0('upstream_', upstream_name)]]
      stopifnot(!length(intersect(upstream$training_ids, m$excluded_ids)))
      displacement <- predict_saved(upstream, x[[upstream_name]][testing$depmap_id, , drop = FALSE])
      predicted <- predict_saved(decoder, design(displacement, b[testing$depmap_id, , drop = FALSE], kind))
      zero <- predict_saved(decoder, design(displacement * 0, b[testing$depmap_id, , drop = FALSE], kind))
      stopifnot(all(is.finite(predicted)), max(abs(zero)) < 1e-10)
      z <- data.table(fold = k, cohort = 'all_475', decoder = kind, response_input = upstream_name,
        depmap_id = testing$depmap_id, observed = testing$observed, predicted = predicted)
      rows[[length(rows) + 1L]] <- z
      pz <- copy(z[depmap_id %chin% paired$depmap_id])
      pz[, cohort := 'paired']
      rows[[length(rows) + 1L]] <- pz
      stage_links[[paste(upstream_name, kind, sep = '__')]] <- list(upstream = upstream, decoder = decoder)
    }
  }
  qsave(list(fold = k, weights = m$weights, stages = stage_links,
    baseline_center = bmu, baseline_scale = bsd, excluded_ids = m$excluded_ids,
    source_hash = hashes$sha256[k], scoring = 'hard_cluster_mean'),
    file.path(out, sprintf('models/fold_%02d.qs', k)))
  checks[[k]] <- data.table(fold = k, upstream_train = length(m$models$upstream_baseline_programs_receptors$training_ids),
    decoder_train = length(m$models$measured_conditioned$training_ids), test = nrow(testing), paired_test = nrow(paired),
    held_recipient_excluded = TRUE, decoder_trained_on_measured_displacements = TRUE, zero_response_anchor = TRUE)
  fwrite(rbindlist(rows), file.path(out, 'data/predictions.tsv'), sep = '\t')
  fwrite(rbindlist(checks), file.path(out, 'data/stage_training_checks.tsv'), sep = '\t')
  fwrite(data.table(status = 'running', completed_folds = k), file.path(out, 'run_status.tsv'), sep = '\t')
  message('Separate-stage chain saved: fold ', k, '/5')
}
pred <- rbindlist(rows)
stopifnot(all(pred[, .N, by = .(cohort, decoder, response_input, depmap_id)]$N == 1L))
metrics <- pred[, .(n = .N, pearson = cor(observed, predicted), rmse = sqrt(mean((observed - predicted)^2))),
  by = .(cohort, decoder, response_input)]
fwrite(metrics, file.path(out, 'data/metrics.tsv'), sep = '\t')
# Independent replay verifies that composing the models did not change the measured-input decoder.
reference <- fread(file.path(source_dir, 'data/diagnostic_functional_predictions.tsv'))[cohort == 'paired_fit_and_test']
for (kind in c('linear', 'quadratic', 'conditioned')) {
  expected <- reference[model == paste0('measured_', kind)]
  actual <- pred[cohort == 'paired' & decoder == kind & response_input == 'measured']
  stopifnot(max(abs(actual$predicted - expected$predicted[match(actual$depmap_id, expected$depmap_id)])) < 1e-10)
}
# Main diagnostic: same decoder, paired cells, measured versus predicted displacement.
main <- pred[cohort == 'paired' & decoder == 'conditioned' &
  response_input %chin% c('measured', 'baseline_programs_receptors')]
main[, panel := ifelse(response_input == 'measured', 'Measured displacement -> function', 'Baseline -> displacement -> function')]
annotation <- main[, .(label = sprintf('n=%d; r=%.3f; RMSE=%.3f', .N, cor(observed, predicted),
  sqrt(mean((observed - predicted)^2)))), by = panel]
main <- merge(main, annotation, by = 'panel')
main[, facet := paste(panel, label, sep = '\n')]
p <- ggplot(main, aes(observed, predicted)) + geom_point(size = 1.3, alpha = .7) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
  facet_wrap(~facet, nrow = 1, scales = 'free') + theme_classic(base_size = 10) + theme(aspect.ratio = 1) +
  labs(x = 'Measured log2 abundance', y = 'Held-out predicted log2 abundance',
    title = 'Separately trained stages | identical measured-response decoder')
ggsave(file.path(out, 'figures/separate_stages_paired.pdf'), p, width = 8.3, height = 4.5)
ggsave(file.path(out, 'figures/separate_stages_paired.png'), p, width = 8.3, height = 4.5, dpi = 150)
report <- c('SEPARATELY TRAINED HARD-SCORE CHAIN', '',
  'Baseline -> displacement fitted on measured molecular responses only.',
  'Displacement -> function fitted on measured displacements plus functional outcomes only.',
  'Models are reused from the completed hard-score fits in script 18; no refitting or new selection.',
  'Conditioned decoder also receives untreated hard-program baseline scores.',
  'Same decoder evaluated with measured and baseline-predicted displacement on identical paired held-out cells.',
  'Functional decoder has only 26-30 paired training recipients per fold; one sensitive recipient across the paired cohort.',
  'All-475 prediction also uses this paired-trained decoder, not the previous 475-cell functional fit.',
  '5 uM molecular / 2.5 uM functional exposure proxy; retrospective outcomes.', '',
  capture.output(print(metrics)), '', 'TRAINING CHECKS', capture.output(print(rbindlist(checks))))
writeLines(report, file.path(out, 'REPORT.txt'))
fwrite(data.table(status = 'complete', completed_folds = 5L), file.path(out, 'run_status.tsv'), sep = '\t')
print(metrics)
