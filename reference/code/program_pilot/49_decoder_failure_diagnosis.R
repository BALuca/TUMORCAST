#!/usr/bin/env Rscript
# Diagnose existing decoder failures. No model fits or parameter selection.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(patchwork)})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/nonlinear_response.R')
source('/work/src/program_pilot/functional_programs.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/49_decoder_failure_diagnosis'
reference <- '/work/results/program_pilot/46_family_holdout_context_signal_pilot'
input_path <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr/data/frozen_comparison_inputs.qs'
for (p in file.path(root, c('data', 'figures'))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
status <- function(s) message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', s)
status('comparing cached predictions on identical family observations')
input <- qread(input_path)
md <- input$metadata
x <- input$scores[['43 averaged / 34 programs']]
projection <- qread(file.path(reference, 'data/family_projection_inputs.qs'))
fm <- projection$metadata
fx <- projection$scores
fam <- fread(file.path(reference, 'data/excluded_families.tsv'))
old <- fread('/work/results/program_pilot/45_whole_cell_regularized_program_svr/data/evaluated_predictions.tsv')
old <- old[bank == '43 averaged / 34 programs' & drug %chin% fam$drug & evaluation == 'Held cells']
old[, `:=`(family = fam$family[match(drug, fam$drug)], model = '45 family seen')]
unseen <- fread(file.path(reference, 'data/evaluated_family_predictions.tsv'))[recipient_set == 'Held cells']
unseen[, model := '46 family excluded']
bounded <- fread('/work/results/program_pilot/48_bounded_program_interaction_function_regression/data/evaluated_predictions.tsv')
bounded <- bounded[model == '48 Program interactions' & evaluation == 'Excluded family / held cells' &
  is.finite(measured_pct) & measured_pct > 0]
bounded[, model := '48 family excluded']
cols <- c('observation_id', 'depmap_id', 'recipient', 'family', 'drug', 'dose', 'measured_pct', 'predicted_log2', 'model')
comparisons <- rbindlist(list(old[, ..cols], unseen[is.finite(measured_pct) & measured_pct > 0, ..cols], bounded[, ..cols]))
check(comparisons[, all(uniqueN(model) == 3L), by = observation_id]$V1 |> all(), 'Matched predictions differ in coverage')
check(comparisons[, all(uniqueN(measured_pct) == 1L), by = observation_id]$V1 |> all(), 'Assay labels differ')
comparisons[, predicted_pct := 100 * 2^predicted_log2]
result <- comparisons[, .(observations = .N, sensitive = sum(measured_pct <= 50), resistant = sum(measured_pct > 50),
  recovered = sum(measured_pct <= 50 & predicted_pct <= 50), false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50),
  r_all_doses = cor(log2(measured_pct / 100), predicted_log2),
  measured_sensitive_median_pct = median(measured_pct[measured_pct <= 50]),
  predicted_sensitive_median_pct = median(predicted_pct[measured_pct <= 50])), by = .(family, drug, model)]
fwrite(comparisons, file.path(root, 'data/matched_held_cell_predictions.tsv'), sep = '\t')
fwrite(result, file.path(root, 'data/matched_all_drug_results.tsv'), sep = '\t')
status('replaying exact46 fits and tracing actual molecular neighbors')
audits <- neighbor_tables <- list()
for (family_name in c('ERBB', 'MEK')) {
  saved <- qread(file.path(reference, 'models', paste0(family_name, '_selected_decoder.qs')))
  fit <- saved$fitted
  tr <- match(saved$training_ids, md$observation_id)
  check(!anyNA(tr), 'Training identifiers missing')
  check(max(abs(colMeans(x[tr, , drop = FALSE]) - fit$x_mean)) < 1e-12, 'Training feature means differ')
  check(max(abs(apply(x[tr, , drop = FALSE], 2, sd) - fit$x_sd)) < 1e-12, 'Training feature scales differ')
  check(max(abs(100 * 2^md$observed_log2[tr] - md$measured_pct[tr])) < 1e-8, 'Response scale mismatch')
  pred <- predict_program_svr(fit, x[tr, , drop = FALSE])
  a <- copy(md[tr, .(observation_id, depmap_id, recipient, drug, dose, measured_pct, observed_log2)])
  alpha <- rep(0, length(tr))
  alpha[fit$model$index] <- drop(fit$model$coefs)
  a[, `:=`(family = family_name, predicted_pct = 100 * 2^pred,
    normalized_residual = (observed_log2 - pred) / fit$y_sd,
    dual_coefficient = alpha, at_negative_bound = alpha <= -fit$model$cost * (1 - 1e-5))]
  audits[[family_name]] <- a
  train_z <- sweep(sweep(x[tr, fit$active, drop = FALSE], 2, fit$x_mean[fit$active], '-'),
    2, fit$x_sd[fit$active], '/')
  query <- fm[family == family_name & is.finite(measured_pct) & measured_pct > 0 & recipient_set != 'Other cells']
  ii <- match(query$observation_id, fm$observation_id)
  query_z <- sweep(sweep(fx[ii, fit$active, drop = FALSE], 2, fit$x_mean[fit$active], '-'),
    2, fit$x_sd[fit$active], '/')
  distance2 <- pmax(outer(rowSums(query_z^2), rowSums(train_z^2), '+') - 2 * tcrossprod(query_z, train_z), 0)
  nearest <- t(apply(distance2, 1, function(d) order(d)[1:3]))
  n <- nrow(query)
  query[, `:=`(nearest_training_id = md$observation_id[tr[nearest[, 1]]],
    nearest_training_cell = md$recipient[tr[nearest[, 1]]],
    nearest_training_drug = md$drug[tr[nearest[, 1]]], nearest_training_dose = md$dose[tr[nearest[, 1]]],
    nearest_training_pct = md$measured_pct[tr[nearest[, 1]]], nearest_distance = sqrt(distance2[cbind(1:n, nearest[, 1])]),
    mean_three_training_pct = rowMeans(matrix(md$measured_pct[tr[nearest]], n, 3)),
    fraction_three_training_sensitive = rowMeans(matrix(md$measured_pct[tr[nearest]] <= 50, n, 3)))]
  for (j in 1:3) query[, (paste0('training_neighbor_', j)) := md$observation_id[tr[nearest[, j]]]]
  # These same-treatment neighbors were NOT part of the46 decoder training set.
  query[, (paste0('same_context_neighbor_', 1:3)) := NA_character_]
  query[, (c('mean_three_same_context_pct', 'fraction_three_same_context_sensitive',
    'same_context_sensitive_fraction', 'same_context_mean_pct', 'same_context_nearest_distance')) := NA_real_]
  for (k in seq_len(n)) {
    peers <- which(query$recipient_set == 'Existing cells' & query$drug == query$drug[k] &
      query$dose == query$dose[k] & query$depmap_id != query$depmap_id[k])
    check(length(peers) >= 3L, 'Insufficient comparison peers')
    ds <- rowSums(sweep(query_z[peers, , drop = FALSE], 2, query_z[k, ], '-')^2)
    closest <- peers[order(ds)[1:3]]
    for (j in 1:3) set(query, k, paste0('same_context_neighbor_', j), query$observation_id[closest[j]])
    set(query, k, 'mean_three_same_context_pct', mean(query$measured_pct[closest]))
    set(query, k, 'fraction_three_same_context_sensitive', mean(query$measured_pct[closest] <= 50))
    set(query, k, 'same_context_sensitive_fraction', mean(query$measured_pct[peers] <= 50))
    set(query, k, 'same_context_mean_pct', mean(query$measured_pct[peers]))
    set(query, k, 'same_context_nearest_distance', sqrt(min(ds)))
  }
  query[, `:=`(predicted_pct = 100 * 2^predict_program_svr(fit, fx[ii, , drop = FALSE]),
    neighbor_sensitive_enrichment = fraction_three_same_context_sensitive - same_context_sensitive_fraction)]
  neighbor_tables[[family_name]] <- query
  status(paste(family_name, 'NEIGHBOR AUDIT SAVED'))
  fwrite(rbindlist(neighbor_tables, fill = TRUE), file.path(root, 'data/all_family_neighbor_audit.tsv'), sep = '\t')
}
audit <- rbindlist(audits)
fwrite(audit, file.path(root, 'data/exact_training_loss_audit.tsv'), sep = '\t')
training_summary <- audit[, .(observations = .N, sensitive = sum(measured_pct <= 50),
  missed_sensitive = sum(measured_pct <= 50 & predicted_pct > 50),
  fraction_missed_sensitive_at_negative_bound = mean(at_negative_bound[measured_pct <= 50 & predicted_pct > 50]),
  median_normalized_residual_missed_sensitive = median(normalized_residual[measured_pct <= 50 & predicted_pct > 50])), by = family]
fwrite(training_summary, file.path(root, 'data/training_loss_summary.tsv'), sep = '\t')
print(training_summary)
nn <- rbindlist(neighbor_tables, fill = TRUE)
nn[, severe_miss := measured_pct <= 40 & predicted_pct >= 60]
nn_summary <- nn[recipient_set == 'Held cells' & severe_miss,
  .(cases = .N, median_measured_pct = median(measured_pct), median_predicted_pct = median(predicted_pct),
    median_nearest_training_pct = median(nearest_training_pct),
    median_three_same_context_pct = median(mean_three_same_context_pct),
    mean_same_context_neighbor_enrichment = mean(neighbor_sensitive_enrichment)), by = family]
fwrite(nn_summary, file.path(root, 'data/severe_miss_neighbor_summary.tsv'), sep = '\t')
print(nn_summary)

# Inspect existing candidates rather than fitting more algorithms. The case replay
# is diagnostic only; evaluated-family outcomes do not select a candidate here.
cv_audit <- case_replay <- cv_neighbor_tables <- list()
grid <- fread(file.path(reference, 'data/whole_cell_validation.tsv'))
for (family_name in c('ERBB', 'MEK')) {
  saved <- qread(file.path(reference, 'models', paste0(family_name, '_selected_decoder.qs')))
  tr <- match(saved$training_ids, md$observation_id)
  severe <- nn[family == family_name & recipient_set == 'Held cells' & severe_miss]
  sx <- fx[match(severe$observation_id, fm$observation_id), , drop = FALSE]
  for (j in 1:6) {
    oof <- rep(NA_real_, length(tr))
    replay <- matrix(NA_real_, nrow(severe), 3)
    for (k in 1:3) {
      fold <- qread(file.path(reference, 'models', sprintf('%s_candidate_%02d_fold_%d.qs', family_name, j, k)))
      oof[match(fold$validation_ids, md$observation_id[tr])] <- fold$prediction
      replay[, k] <- predict_program_svr(fold$fitted, sx)
      if (j == saved$chosen$candidate) {
        f <- fold$fitted
        ti <- match(fold$fitting_ids, md$observation_id)
        vi <- match(fold$validation_ids, md$observation_id)
        check(!any(md$depmap_id[ti] %chin% md$depmap_id[vi]), 'Neighbor audit overlaps validation cells')
        tz <- sweep(sweep(x[ti, f$active, drop = FALSE], 2, f$x_mean[f$active], '-'),
          2, f$x_sd[f$active], '/')
        vz <- sweep(sweep(x[vi, f$active, drop = FALSE], 2, f$x_mean[f$active], '-'),
          2, f$x_sd[f$active], '/')
        ds <- pmax(outer(rowSums(vz^2), rowSums(tz^2), '+') - 2 * tcrossprod(vz, tz), 0)
        closest <- t(apply(ds, 1, function(d) order(d)[1:3]))
        neighbors <- matrix(md$measured_pct[ti[closest]], length(vi), 3)
        cv_neighbor_tables[[paste(family_name, k)]] <- cbind(md[vi,
          .(observation_id, recipient, drug, dose, measured_pct)],
          data.table(family = family_name, fold = k, predicted_pct = 100 * 2^fold$prediction,
            sensitive_neighbor_fraction = rowMeans(neighbors <= 50),
            mean_three_neighbor_pct = rowMeans(neighbors),
            training_sensitive_fraction = mean(md$measured_pct[ti] <= 50)))
      }
    }
    setting <- grid[family == family_name & candidate == j]
    cv_audit[[paste(family_name, j)]] <- cbind(setting, data.table(selected = j == saved$chosen$candidate,
      sensitive = sum(md$measured_pct[tr] <= 50), resistant = sum(md$measured_pct[tr] > 50),
      recovered_sensitive = sum(md$measured_pct[tr] <= 50 & oof <= -1),
      false_sensitive = sum(md$measured_pct[tr] > 50 & oof <= -1),
      median_sensitive_prediction_pct = median(100 * 2^oof[md$measured_pct[tr] <= 50])))
    case_replay[[paste(family_name, j)]] <- cbind(severe[, .(family, observation_id, recipient, drug, dose, measured_pct)],
      data.table(candidate = j, cost = setting$cost, gamma_multiplier = setting$gamma_multiplier,
        selected = j == saved$chosen$candidate, prediction_from_saved_fold_models_pct = 100 * 2^rowMeans(replay)))
  }
}
fwrite(rbindlist(cv_audit), file.path(root, 'data/existing_candidate_selection_audit.tsv'), sep = '\t')
fwrite(rbindlist(case_replay), file.path(root, 'data/failure_replay_existing_fold_models.tsv'), sep = '\t')
print(rbindlist(cv_audit))
cv_neighbors <- rbindlist(cv_neighbor_tables)
fwrite(cv_neighbors, file.path(root, 'data/whole_cell_validation_neighbor_coherence.tsv'), sep = '\t')
coherence <- cv_neighbors[, .(sensitive_queries = sum(measured_pct <= 50),
  missed_sensitive_queries = sum(measured_pct <= 50 & predicted_pct > 50),
  missed_sensitive_with_three_sensitive_neighbors = sum(measured_pct <= 50 & predicted_pct > 50 & sensitive_neighbor_fraction == 1),
  resistant_queries_with_three_sensitive_neighbors = sum(measured_pct > 50 & sensitive_neighbor_fraction == 1),
  mean_sensitive_neighbor_fraction_for_sensitive_queries = mean(sensitive_neighbor_fraction[measured_pct <= 50]),
  mean_sensitive_neighbor_fraction_for_resistant_queries = mean(sensitive_neighbor_fraction[measured_pct > 50])), by = family]
fwrite(coherence, file.path(root, 'data/whole_cell_validation_neighbor_summary.tsv'), sep = '\t')
print(coherence)

# A direct local decoder check isolates information recovered from these inputs.
# Fixed three-neighbor averaging is descriptive, not a tuned new workflow.
neighbor_check <- copy(cv_neighbors)
neighbor_check[, `:=`(observed = log2(measured_pct / 100),
  neighbor = log2(mean_three_neighbor_pct / 100), svr = log2(predicted_pct / 100))]
neighbor_check[, `:=`(yc = observed - mean(observed), nc = neighbor - mean(neighbor),
  sc = svr - mean(svr)), by = .(family, drug, dose, fold)]
direct_signal <- neighbor_check[, .(rows = .N, sensitive = sum(measured_pct <= 50),
  resistant = sum(measured_pct > 50),
  SVR_sensitive_recovered = sum(measured_pct <= 50 & predicted_pct <= 50),
  neighbor_sensitive_recovered = sum(measured_pct <= 50 & mean_three_neighbor_pct <= 50),
  SVR_false_sensitive = sum(measured_pct > 50 & predicted_pct <= 50),
  neighbor_false_sensitive = sum(measured_pct > 50 & mean_three_neighbor_pct <= 50),
  SVR_all_dose_r = cor(observed, svr), neighbor_all_dose_r = cor(observed, neighbor),
  SVR_same_treatment_fold_centered_r = cor(yc, sc),
  neighbor_same_treatment_fold_centered_r = cor(yc, nc),
  SVR_sensitive_median_pct = median(predicted_pct[measured_pct <= 50]),
  neighbor_sensitive_median_pct = median(mean_three_neighbor_pct[measured_pct <= 50])), by = family]
fwrite(direct_signal, file.path(root, 'data/direct_neighbor_signal_check.tsv'), sep = '\t')
print(direct_signal)

# Compare the same held-cell outcomes: no new function fit is hidden in these views.
result[, label := paste(family, drug, sep = ' | ')]
result[, label := factor(label, levels = rev(unique(result[order(family, drug), label])))]
colours <- c('45 family seen' = '#4A4A4A', '46 family excluded' = '#B24C35', '48 family excluded' = '#27818A')
p_order <- ggplot(result, aes(r_all_doses, label, colour = model)) +
  geom_vline(xintercept = 0, colour = 'grey65', linetype = 2) +
  geom_point(position = position_dodge(width = .5), size = 2.5) +
  scale_colour_manual(values = colours) + scale_x_continuous(limits = c(-1, 1)) +
  labs(title = 'Same drugs, same nine held cells, same34 molecular programs',
    subtitle = '45 learned from familiar cells receiving these families;46/48 excluded the entire family.\nERBB45/46: same SVR C10 and gamma1/34. MEK also changes regularization; not an isolated family effect.',
    x = 'Pearson r across all doses', y = NULL, colour = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', axis.text.y = element_text(size = 9))
ggsave(file.path(root, 'figures/family_seen_vs_excluded.pdf'), p_order, width = 10, height = 5.5)
ggsave(file.path(root, 'figures/family_seen_vs_excluded.png'), p_order, width = 10, height = 5.5, dpi = 160)
rates <- melt(result[, .(label, model, sensitivity_recovered = recovered / sensitive,
  resistant_called_sensitive = false_sensitive / resistant)], id.vars = c('label', 'model'),
  variable.name = 'measurement', value.name = 'fraction')
p_recovery <- ggplot(rates, aes(fraction, label, colour = model)) +
  geom_point(position = position_dodge(width = .5), size = 2.4, na.rm = TRUE) +
  facet_wrap(~measurement, nrow = 1) + scale_colour_manual(values = colours) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  labs(title = 'Recovered sensitivity and false sensitivity calls',
    subtitle = 'All doses of each drug in the same nine held cells. Threshold:50% vehicle.\nBlank sensitivity row: no measured sensitive observations for that drug.',
    x = 'Fraction of dose-cell observations', y = NULL, colour = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', strip.text = element_text(size = 9))
ggsave(file.path(root, 'figures/sensitivity_recovery_vs_false_calls.pdf'), p_recovery, width = 11, height = 5.5)
ggsave(file.path(root, 'figures/sensitivity_recovery_vs_false_calls.png'), p_recovery, width = 11, height = 5.5, dpi = 160)

# Show the worst clear misses, with actual program profiles and outcome bars.
cases <- nn[recipient_set == 'Held cells' & severe_miss][order(-(predicted_pct - measured_pct))]
cases <- cases[!duplicated(drug)][1:min(6L, .N)]
check(nrow(cases) > 0, 'No severe miss cases')
cases[, case := sprintf('C%02d', .I)]
fwrite(cases, file.path(root, 'data/illustrated_failure_cases.tsv'), sep = '\t')
local_tables <- list()
for (k in seq_len(nrow(cases))) {
  saved <- qread(file.path(reference, 'models', paste0(cases$family[k], '_selected_decoder.qs')))
  tr <- match(saved$training_ids, md$observation_id)
  fit <- saved$fitted
  ztr <- sweep(sweep(x[tr, fit$active, drop = FALSE], 2, fit$x_mean[fit$active], '-'),
    2, fit$x_sd[fit$active], '/')
  qx <- fx[match(cases$observation_id[k], fm$observation_id), fit$active, drop = FALSE]
  zq <- sweep(sweep(qx, 2, fit$x_mean[fit$active], '-'), 2, fit$x_sd[fit$active], '/')
  ds <- rowSums(sweep(ztr, 2, drop(zq), '-')^2)
  ii <- order(ds)[1:30]
  local_tables[[k]] <- cbind(data.table(case = cases$case[k], neighbor_rank = 1:30,
    distance = sqrt(ds[ii]), kernel_similarity = exp(-fit$model$gamma * ds[ii])),
    audit[family == cases$family[k]][match(md$observation_id[tr[ii]], observation_id)])
  check(nrow(local_tables[[k]]) == 30L &&
    identical(local_tables[[k]]$observation_id, md$observation_id[tr[ii]]), 'Local trace rows differ')
}
local <- rbindlist(local_tables, fill = TRUE)
fwrite(local, file.path(root, 'data/nearest_thirty_training_profiles.tsv'), sep = '\t')

# Named baseline context for the actual failures, from an already produced cache.
# Hotspot flags are measurements to investigate, not proof of the cause.
mutations <- fread('/work/results/program_pilot/22_compressed_baseline_mutation_conditioning_MEK_holdout/data/mutation_features.tsv')
context_ids <- unique(c(cases$depmap_id,
  sub('__.*', '', unlist(cases[, paste0('same_context_neighbor_', 1:3), with = FALSE]))))
fwrite(mutations[depmap_id %chin% context_ids & gene %chin% c('BRAF', 'KRAS', 'NRAS', 'NF1')],
  file.path(root, 'data/failure_case_baseline_hotspot_context.tsv'), sep = '\t')
profiles <- annotations <- list()
for (k in seq_len(nrow(cases))) {
  family_name <- cases$family[k]
  fit <- qread(file.path(reference, 'models', paste0(family_name, '_selected_decoder.qs')))$fitted
  ids <- c(cases$observation_id[k], unlist(cases[k, paste0('training_neighbor_', 1:3), with = FALSE]),
    unlist(cases[k, paste0('same_context_neighbor_', 1:3), with = FALSE]))
  role <- c('Missed sensitive', rep('Decoder training neighbor', 3), rep('Excluded family neighbor', 3))
  raw <- matrix(NA_real_, 7, ncol(x), dimnames = list(paste(cases$case[k], 1:7, sep = '_'), colnames(x)))
  m <- match(ids, md$observation_id)
  raw[!is.na(m), ] <- x[m[!is.na(m)], , drop = FALSE]
  f <- match(ids, fm$observation_id)
  raw[!is.na(f), ] <- fx[f[!is.na(f)], , drop = FALSE]
  meta <- rbindlist(list(md[, .(observation_id, recipient, drug, dose, measured_pct)],
    fm[, .(observation_id, recipient, drug, dose, measured_pct)]))
  meta <- unique(meta, by = 'observation_id')[match(ids, observation_id)]
  check(all(is.finite(raw)) && identical(meta$observation_id, unname(ids)), 'Failure profiles do not align')
  pred <- predict_program_svr(fit, raw)
  scaled <- sweep(sweep(raw, 2, fit$x_mean, '-'), 2, fit$x_sd, '/')
  profiles[[k]] <- t(scaled)
  annotations[[k]] <- data.table(column_id = rownames(raw), case = cases$case[k], role = role,
    label = paste(meta$recipient, meta$drug, meta$dose, sep = ' | '), observation_id = ids,
    observed_pct = meta$measured_pct, predicted_pct = 100 * 2^pred)
}
plotted <- do.call(cbind, profiles)
ann <- rbindlist(annotations)
check(identical(colnames(plotted), ann$column_id), 'Failure plot annotations differ')
fwrite(ann, file.path(root, 'data/failure_heatmap_annotations.tsv'), sep = '\t')
fwrite(data.table(program = rownames(plotted), as.data.table(plotted)),
  file.path(root, 'data/failure_heatmap_svr_scaled_programs.tsv'), sep = '\t')
cm <- data.frame(Role = ann$role, row.names = ann$column_id)
maximum <- ceiling(max(c(ann$observed_pct, ann$predicted_pct)) / 50) * 50
top <- list(show_legend = FALSE,
  'Measured % vehicle' = ComplexHeatmap::anno_barplot(ann$observed_pct, ylim = c(0, maximum),
    gp = grid::gpar(fill = '#27818A', col = NA), height = grid::unit(15, 'mm')),
  'Predicted % vehicle' = ComplexHeatmap::anno_barplot(ann$predicted_pct, ylim = c(0, maximum),
    gp = grid::gpar(fill = '#B24C35', col = NA), height = grid::unit(15, 'mm')),
  'Assay class' = sensitivity_track(ann$observed_pct))
hm <- draw_heatmap(plotted, scale_rows = FALSE, cluster_rows = hclust(dist(plotted), method = 'ward.D2'),
  cluster_columns = FALSE, cluster_column_slices = FALSE,
  column_split = factor(ann$case, levels = cases$case), col_anno_df = cm, col_anno_cols = 'Role', ann_top_args = top,
  show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7), show_column_names = TRUE,
  column_labels = ann$label, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 5),
  value_range = c(-3, 0, 3), value_colors = c('#2166AC', 'white', '#B2182B'),
  hm_name = 'SVR-scaled program', title = 'Actual missed responses and their nearest molecular neighbors\nEach case: query;3 decoder-training neighbors;3 familiar same-drug/dose neighbors excluded from fit')
save_heatmap(hm, 'failure_cases_and_actual_neighbors', 15, 10, file.path(root, 'figures'), list(sensitivity_legend))
for (case_id in c('C01', cases[recipient == 'LOX IMVI' & drug == 'trametinib', case])) {
  idx <- which(ann$case == case_id)
  m <- plotted[, idx, drop = FALSE]
  a <- ann[idx]
  metadata <- cm[idx, , drop = FALSE]
  top_case <- list(show_legend = FALSE,
    'Measured % vehicle' = ComplexHeatmap::anno_barplot(a$observed_pct, ylim = c(0, maximum),
      gp = grid::gpar(fill = '#27818A', col = NA), height = grid::unit(16, 'mm')),
    'Predicted % vehicle' = ComplexHeatmap::anno_barplot(a$predicted_pct, ylim = c(0, maximum),
      gp = grid::gpar(fill = '#B24C35', col = NA), height = grid::unit(16, 'mm')),
    'Assay class' = sensitivity_track(a$observed_pct))
  colour_limit <- max(abs(m))
  hm_case <- draw_heatmap(m, scale_rows = FALSE, cluster_rows = hclust(dist(m), method = 'ward.D2'),
    cluster_columns = FALSE, col_anno_df = metadata, col_anno_cols = 'Role', ann_top_args = top_case,
    show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 8), show_column_names = TRUE,
    column_labels = a$label, column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 7),
    value_range = c(-colour_limit, 0, colour_limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'SVR-scaled program', column_title_gp = grid::gpar(fontsize = 11),
    title = paste0(case_id, ' | ', a$label[1], '\n',
      'Query | 3 training neighbors | 3 same-treatment peers'))
  save_heatmap(hm_case, paste0(case_id, '_actual_failure_trace'), 8, 9, file.path(root, 'figures'), list(sensitivity_legend))
}

# Keep the scientific diagnosis beside its evidence; do not propose another fit
# as a substitute for demonstrating which part of the current decoder fails.
report <- c(
  '49: DIAGNOSIS OF THE EXISTING MOLECULAR-DISPLACEMENT TO VIABILITY DECODER',
  '',
  'No new model was fitted or selected. This replays saved 45/46 models and',
  'predictions, compares 48, and traces actual molecular neighbors. All doses',
  'remain in the data. Sensitive means measured viability <=50% of vehicle.',
  '',
  '1. The selected SVR already misses sensitive observations in its own fit.',
  sprintf('   %s: %d/%d sensitive fitting observations predicted resistant.',
    training_summary$family, training_summary$missed_sensitive, training_summary$sensitive),
  '   Exact response transformation and training-only feature means/scales',
  '   agree with the original inputs. These checks did not identify a target',
  '   sign, scaling or observation-join error.',
  '',
  '2. Local program geometry contains phenotype information in whole-cell CV.',
  '   All observations from the query cell are excluded from its fold fit.',
  '   Distances use that fit\'s standardized program inputs, without labels.',
  sprintf('   %s: sensitive queries have %.1f%% sensitive neighbors on average;',
    coherence$family, 100 * coherence$mean_sensitive_neighbor_fraction_for_sensitive_queries),
  sprintf('   resistant queries have %.1f%%. Among missed sensitive queries,',
    100 * coherence$mean_sensitive_neighbor_fraction_for_resistant_queries),
  sprintf('   %s has %d with all three nearest training profiles sensitive.',
    coherence$family, coherence$missed_sensitive_with_three_sensitive_neighbors),
  sprintf('   %s also has %d resistant queries with all three neighbors sensitive.',
    coherence$family, coherence$resistant_queries_with_three_sensitive_neighbors),
  '   Thus local information exists, but neighborhoods are not perfectly',
  '   phenotype-consistent. This is a diagnostic, not a validated new predictor.',
  '',
  '3. Hs 766T / TAK-733 / 0.5 uM: a concrete fitting failure.',
  '   Measured viability 5.48%; selected SVR prediction 92.78%. Its three',
  '   nearest actual fitting profiles have measured viability 19.50%,17.21%,',
  '   19.62%, but the model predicts 95.06%,88.66%,89.66% for those training',
  '   observations themselves. Their kernel similarities to the query are',
  '   0.9973-0.9978; the surrounding 30 profiles include 12 sensitive and',
  '   18 resistant observations. Local sensitive structure is smoothed into',
  '   a mixed neighborhood. All six existing candidate fold settings still',
  '   predict roughly 90-93% for this query. Changing saved settings alone',
  '   does not rescue it.',
  '',
  '4. LOX IMVI / trametinib / 0.5 uM: a different ambiguity.',
  '   Measured viability 24.38%; predicted 94.69%. Its three nearest fitting',
  '   profiles average 99.73% measured viability; three familiar cells at',
  '   the same drug/dose average 92.39%. Thus these local profiles do not',
  '   provide the same sensitive-neighbor evidence as the TAK-733 case.',
  '   The cached baseline annotation marks LOX IMVI BRAF-hotspot positive,',
  '   while those three same-treatment neighbors are negative. This is an',
  '   observed context difference, not evidence that BRAF causes the mismatch.',
  '',
  '5. Current candidate-selection objective can reward the wrong compromise.',
  '   MEK selected C=1/gamma-multiplier=0.1 with same-treatment-centered CV',
  '   correlation 0.04787, versus 0.04604 for C=10/multiplier=1. Across all',
  '   doses the corresponding correlations are 0.436 and 0.544. Sensitive',
  '   recovery is 327/1635 versus 555/1635; false sensitive calls increase',
  '   from 106/4937 to 338/4937. The tiny centered-correlation difference',
  '   hides a material sensitivity/specificity tradeoff. It does not justify',
  '   choosing a model after inspecting held-family outcomes.',
  '',
  'Interpretation: a simple calibration shift is insufficient. The global',
  'decoder suppresses some local phenotype information even in training,',
  'and other failures involve similar displacements with different outcomes.',
  'These must be addressed separately. The 47/48 model changes did not',
  'establish a consistent rescue and are not accepted as an improved bridge.',
  '',
  'Figure guide: C01_actual_failure_trace.pdf shows the first failure;',
  'C04_actual_failure_trace.pdf shows the second. Columns are query, three',
  'actual fit neighbors, three same-drug/dose familiar-cell peers excluded',
  'from the function fit. Top bars show measured and predicted % vehicle.',
  'Colors show the actual standardized 34-program displacement inputs.',
  'The heatmap is not rescaled gene-by-gene or converted to sensitivity.',
  'failure_cases_and_actual_neighbors.pdf contains five drug-specific cases.',
  '',
  'Limits: these are retrospective diagnostics. Program bank 43 was discovered',
  'using outcomes across the full panel, including now-evaluated cells/drugs.',
  'Decoder whole-cell exclusion therefore does not make discovery independent.',
  '45 versus 46 changes family exposure; for MEK it also changes C/gamma,',
  'so that MEK comparison does not isolate family exclusion. No baseline-to-',
  'predicted-displacement model is tested here, and no mutation rescue is fitted.'
)
report <- c(report, '', 'Direct local-decoder check (fixed three-neighbor average)',
  'Whole-cell validation with the same saved fitting/query splits; no labels',
  'enter neighbor selection. Averaging measured training-neighbor viability',
  'improves pooled ordering and sensitivity recovery, but does not restore',
  'across-cell ordering within drug/dose. Centering includes fold, matching46.',
  sprintf('   %s: pooled r SVR %.3f; neighbors %.3f. Same-treatment r %.3f; %.3f.',
    direct_signal$family, direct_signal$SVR_all_dose_r, direct_signal$neighbor_all_dose_r,
    direct_signal$SVR_same_treatment_fold_centered_r,
    direct_signal$neighbor_same_treatment_fold_centered_r),
  sprintf('   %s: sensitive recovered SVR %d/%d; neighbors %d/%d.',
    direct_signal$family, direct_signal$SVR_sensitive_recovered, direct_signal$sensitive,
    direct_signal$neighbor_sensitive_recovered, direct_signal$sensitive),
  sprintf('   %s: false sensitive SVR %d/%d; neighbors %d/%d.',
    direct_signal$family, direct_signal$SVR_false_sensitive, direct_signal$resistant,
    direct_signal$neighbor_false_sensitive, direct_signal$resistant),
  'This limits the earlier inference: local phenotype coherence by itself',
  'does not show that cell-specific sensitivity is encoded in this geometry.',
  'Making the current kernel more local cannot be assumed to fix that gap.',
  'Before a new decoder, compare matched same-treatment molecular contrasts',
  'in41 per-drug programs against43 global hard means. Preserve all doses.',
  'This separates loss during bank aggregation from failure of the decoder',
  'or an outcome/context mismatch. Discovery labels remain retrospective.')
writeLines(report, file.path(root, 'REPORT.txt'))
provenance_paths <- c(input_path, file.path(reference, 'data/family_projection_inputs.qs'),
  file.path(reference, 'models/ERBB_selected_decoder.qs'),
  file.path(reference, 'models/MEK_selected_decoder.qs'),
  '/work/src/program_pilot/49_decoder_failure_diagnosis.R')
fwrite(data.table(path = provenance_paths, md5 = unname(tools::md5sum(provenance_paths))),
  file.path(root, 'data/diagnosis_provenance.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('diagnosis tables and figures complete; no model fitted')
