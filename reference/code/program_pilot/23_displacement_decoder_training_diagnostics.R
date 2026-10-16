#!/usr/bin/env Rscript
# Replay frozen decoders to locate failure before held-cell generalization.
# No program discovery, model fitting or parameter selection occurs here.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
})
source('/work/src/program_pilot/nonlinear_response.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/23_displacement_decoder_training_diagnostics'
dir.create(file.path(root, 'data'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, 'figures'), showWarnings = FALSE)
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
training_values <- summaries <- neighbors <- query_values <- hashes <- list()
for (family in c('MEK', 'ERBB')) {
  parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
  paths <- file.path(parent, c('settings.qs', 'data/prepared_scores.qs', 'models/nonlinear_RBF_SVR_D.qs'))
  hashes[[family]] <- data.table(path = paths,
    sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
  settings <- qread(paths[1L])
  prepared <- qread(paths[2L])
  saved <- qread(paths[3L])
  stopifnot(!length(intersect(settings$held_ids, saved$training_ids)))
  train <- as.data.table(panel$sensitivity)[depmap_id %in% saved$training_ids &
    !drug %in% settings$excluded_drugs & is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  stopifnot(setequal(unique(train$depmap_id), saved$training_ids), !anyDuplicated(train$observation_id))
  d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  stopifnot(all(is.finite(d)))
  zero <- predict_program_svr(saved$model, matrix(0, 1L, ncol(d)))
  raw <- predict_program_svr(saved$model, d)
  train[, `:=`(family = family, raw_prediction_pct = 100 * 2^raw,
    predicted_pct = 100 * 2^(raw - zero))]
  train[, `:=`(observed_pct = prism_abundance_pct,
    observed_class = ifelse(prism_abundance_pct <= 50, 'Sensitive', 'Resistant'))]
  stopifnot(all(is.finite(train$predicted_pct)))
  selected <- train[drug %in% c('anastrozole', 'ixazomib', 'vemurafenib') & dose == .5]
  training_values[[family]] <- selected
  summaries[[family]] <- selected[, .(training_cells = .N, sensitive = sum(observed_pct <= 50),
    sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
    resistant = sum(observed_pct > 50), false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50),
    observed_median_pct = median(observed_pct), predicted_median_pct = median(predicted_pct),
    raw_prediction_median_pct = median(raw_prediction_pct), pearson_r = cor(observed_pct, predicted_pct),
    zero_prediction_pct = 100 * 2^zero), by = .(family, drug, dose)]

  # Neighbors use exactly the standardization and RBF kernel of the frozen SVR.
  f <- saved$model
  z <- sweep(sweep(d[, f$active, drop = FALSE], 2L, f$x_mean[f$active], '-'), 2L, f$x_sd[f$active], '/')
  for (drug_name in c('anastrozole', 'ixazomib', 'vemurafenib')) {
    cell_name <- if (drug_name == 'vemurafenib') 'LOX IMVI' else 'CFPAC-1'
    query <- as.data.table(panel$sensitivity)[recipient == cell_name & drug == drug_name & dose == .5]
    stopifnot(nrow(query) == 1L, query$depmap_id %in% settings$held_ids)
    dq <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
    zq <- as.numeric((dq[, f$active] - f$x_mean[f$active]) / f$x_sd[f$active])
    squared_distance <- rowSums(sweep(z, 2L, zq, '-')^2)
    distances <- copy(train)
    distances[, `:=`(standardized_profile_distance = sqrt(squared_distance / ncol(z)),
      kernel_similarity = exp(-f$model$gamma * squared_distance), query_cell = cell_name,
      query_drug = drug_name, query_observed_pct = query$prism_abundance_pct)]
    nearest <- distances[order(standardized_profile_distance)][1:8]
    nearest[, `:=`(neighbor_scope = 'All training treatments', rank = seq_len(.N))]
    same <- distances[drug == drug_name & dose == .5][order(standardized_profile_distance)][1:5]
    same[, `:=`(neighbor_scope = 'Same drug and dose', rank = seq_len(.N))]
    neighbors[[paste(family, drug_name)]] <- rbindlist(list(nearest, same))
    query_values[[paste(family, drug_name)]] <- data.table(family = family, cell = cell_name,
      drug = drug_name, dose = .5, observed_pct = query$prism_abundance_pct,
      predicted_pct = 100 * 2^(predict_program_svr(f, dq) - zero))
  }
}
summary <- rbindlist(summaries)
values <- rbindlist(training_values)
fwrite(values, file.path(root, 'data/training_predictions.tsv'), sep = '\t')
fwrite(summary, file.path(root, 'data/training_summary.tsv'), sep = '\t')
fwrite(rbindlist(neighbors), file.path(root, 'data/query_training_neighbors.tsv'), sep = '\t')
fwrite(rbindlist(query_values), file.path(root, 'data/query_predictions.tsv'), sep = '\t')
manifest <- rbindlist(c(hashes, list(data.table(path = c(panel_path,
  '/work/src/program_pilot/23_displacement_decoder_training_diagnostics.R',
  '/work/src/program_pilot/nonlinear_response.R'), sha256 = vapply(c(panel_path,
  '/work/src/program_pilot/23_displacement_decoder_training_diagnostics.R',
  '/work/src/program_pilot/nonlinear_response.R'), digest::digest, character(1), file = TRUE, algo = 'sha256')))))
fwrite(manifest, file.path(root, 'data/input_hashes.tsv'), sep = '\t')

values[, panel := paste(family, drug, '0.5 uM', sep = ' | ')]
summary[, `:=`(panel = paste(family, drug, '0.5 uM', sep = ' | '),
  annotation = sprintf('Sensitive called: %d/%d\nMedian observed/predicted: %.0f%% / %.0f%%\nr=%.2f',
    sensitive_called, sensitive, observed_median_pct, predicted_median_pct, pearson_r))]
p <- ggplot(values, aes(observed_pct, predicted_pct, color = observed_class)) +
  geom_point(size = 2) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_vline(xintercept = 50, color = 'grey65', linetype = 2) +
  geom_hline(yintercept = 50, color = 'grey65', linetype = 2) +
  geom_text(data = summary, aes(x = Inf, y = Inf, label = annotation), inherit.aes = FALSE,
    hjust = 1.03, vjust = 1.1, size = 2.5) +
  facet_wrap(~panel, ncol = 3, scales = 'free') +
  scale_color_manual(values = c(Sensitive = '#2166AC', Resistant = '#B2182B')) +
  theme_classic(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom') +
  labs(x = 'Observed training-cell abundance (% vehicle)',
    y = 'Frozen SVR prediction (% vehicle)', color = 'Observed response',
    title = 'The displacement decoder misses sensitivity in its fitting cells',
    subtitle = 'Frozen models replayed; each point is a training cell at the named drug and dose',
    caption = 'Fitting observations, not validation. Predictions use the existing zero-response anchor; axes scale independently.')
ggsave(file.path(root, 'figures/training_decoder_failure.pdf'), p, width = 10.5, height = 7)
ggsave(file.path(root, 'figures/training_decoder_failure.png'), p, width = 10.5, height = 7, dpi = 140)
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
fwrite(data.table(stage = 'COMPLETE', models_refitted = 0L, time = as.character(Sys.time())),
  file.path(root, 'run_status.tsv'), sep = '\t')
print(summary)
