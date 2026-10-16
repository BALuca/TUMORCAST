#!/usr/bin/env Rscript
# Validate the controlled HGF-input replacement and compare its frozen outputs.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(digest); library(glmnet)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/135_v128_with_corrected_HGF_inputs/linear'
parent <- '/work/results/program_pilot/107_within_pair_departure_interaction/linear'
programs <- sprintf('P%02d', 1:40)
source_weight <- .2
interaction_penalty <- .05
assay_penalty <- .01
external_within_weight <- .8
primary_endpoint <- 'normalized_cell_count'
challenge_genes <- c('HGF', 'NRG1')
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
drugs <- bank_selection$drug
contexts <- fread(file.path(root, 'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]
context_key <- function(drug, dose) paste(drug, sprintf('%.10g', dose), sep = '|')
context_index <- function(drug, dose) match(context_key(drug, dose), context_key(contexts$drug, contexts$dose))
ligands <- fread(file.path(root, 'data/time_matched_mema_molecular_inputs.tsv'))
external_ligands <- fread(file.path(root, 'data/external_ligand_programs.tsv'))
external_molecular <- fread(file.path(root, 'data/external_exact_drug_molecular_observations.tsv'))
external <- fread(file.path(root, 'data/external_combination_outcomes.tsv'))
mutation <- fread('/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv',
                  select = c('V1', 'BRAF (673)'))
setnames(mutation, c('depmap_id', 'BRAF_hotspot'))
training_labels <- fread(file.path(root, 'data/mema_nonchallenge_training_outcomes.tsv'))
pairs <- fread(file.path(root, 'data/combination_training_outcomes.tsv'))
original_final <- qread(file.path(root, 'models/final_normalized_cell_count.qs'))
profiles <- fread(file.path(root, 'data/time_matched_HGF_program_inputs.tsv'))
source_views <- profiles[!grepl('^Leave_out_|^Primary_', source_profile)]

# Read only the exact function definitions; never execute the complete fit again.
required <- c('predict_rna', 'signature_node', 'pair_features', 'predict_single', 'predict_pair',
              'mema_features', 'center_geometry', 'external_features', 'fit_interaction',
              'select_penalty', 'fit_selected', 'predict_external_departure')
expressions <- parse('/work/src/program_pilot/135_v128_with_corrected_HGF_inputs.R')
found <- character()
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name('<-')) &&
      is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% required) {
    eval(expression, envir = .GlobalEnv)
    found <- c(found, as.character(expression[[2]]))
  }
}
stopifnot(setequal(required, found))
penalty_audit <- list()
query_hgf <- function(model, Dh, profile) {
  ctx <- model$context
  lap <- predict_single(model, data.table(depmap_id = 'ACH-000248', drug = 'Lapatinib ditosylate', dose = 5))
  Qh <- cbind(1, sweep(sweep(sweep(Dh, 2, ctx$fm), 2, ctx$fs, '/') %*% ctx$fv, 2, ctx$fscale, '/'))
  scale_map <- solve(crossprod(ctx$Q[1:10, ]) / 10 + diag(c(0, 1, 1, 1)),
                     crossprod(ctx$Q[1:10, ], log(ctx$sigma[1:10, ])) / 10)
  minimum <- matrix(apply(ctx$sigma[1:10, ], 2, min), nrow(Dh), 40, byrow = TRUE)
  maximum <- matrix(apply(ctx$sigma[1:10, ], 2, max), nrow(Dh), 40, byrow = TRUE)
  hscale <- pmin(pmax(exp(Qh %*% scale_map), minimum), maximum)
  H40 <- Dh / hscale
  H <- cbind(H40, H40 %*% ctx$W)
  L <- lap$program_activity[rep(1L, nrow(Dh)), , drop = FALSE]
  product <- pair_features(L, H, ctx$W, model$interaction_basis)
  interaction <- -as.numeric(sweep(product, 2, ctx$interaction_scale, '/') %*% model$par$Gamma)
  standalone <- -((H40 %*% ctx$W) * as.numeric(exp(Qh %*% model$par$log_amplitude)))[, 1] *
    model$ligand_single_conversion
  data.table(source_profile = profile, lapatinib_pct = 100 * 2^(-lap$sensitivity),
    HGF_alone_pct = 100 * 2^standalone, interaction_log2 = interaction,
    combination_pct_vehicle = 100 * 2^(-lap$sensitivity + standalone + interaction))
}
mf <- mema_features(original_final)
phi <- pair_features(mf$A, mf$B, original_final$context$W, original_final$interaction_basis)
reverse <- pair_features(mf$B, mf$A, original_final$context$W, original_final$interaction_basis)
zero <- pair_features(mf$A, mf$B * 0, original_final$context$W, original_final$interaction_basis)
old <- qread(file.path(parent, 'models/final.qs'))
stopifnot(max(abs(phi - reverse)) < 1e-12, max(abs(zero)) < 1e-12,
          identical(old$rna, original_final$rna), identical(old$par$level, original_final$par$level),
          identical(old$par$log_amplitude, original_final$par$log_amplitude),
          identical(old$template, original_final$template), !any(training_labels$gene %in% challenge_genes))
write_table(data.table(check = c('symmetry', 'zero_input', 'single_drug_unchanged', 'no_MEMA_HGF_NRG1_training'),
  maximum_difference = c(max(abs(phi - reverse)), max(abs(zero)), 0, 0)), 'structural_validation')
au <- profiles[CellLine == 'AU565']
source_queries <- query_hgf(original_final, as.matrix(au[, ..programs]), au$source_profile)
write_table(source_queries, 'HGF_source_study_queries_frozen_model')
primary_saved <- fread(file.path(root, 'data/au565_hgf_5uM_predictions_frozen.tsv'))[
  hgf_context == 'Study-balanced HGF (same recipe as training)' & assumed_HGF_route == 1]
stopifnot(abs(primary_saved$combination_pct_control -
  source_queries[source_profile == 'Primary_study_mean', combination_pct_vehicle]) < 1e-8)

# Source dependence is assessed in both the query and all training recipients.
original_ligands <- copy(ligands)
original_external_ligands <- copy(external_ligands)
original_source_views <- copy(source_views)
exclusions <- list()
for (profile in unique(profiles$source_profile)[grepl('^Leave_out_', unique(profiles$source_profile))]) {
  ligands <- copy(original_ligands)
  external_ligands <- copy(original_external_ligands)
  source_views <- original_source_views[source_profile != sub('^Leave_out_', '', profile)]
  replacement <- profiles[source_profile == profile]
  for (cell in unique(replacement$CellLine)) {
    values <- replacement[CellLine == cell, ..programs]
    if (cell %in% ligands$CellLine) ligands[CellLine == cell & gene == 'HGF', (programs) := values]
    if (cell %in% external_ligands$CellLine) external_ligands[CellLine == cell & gene == 'HGF', (programs) := values]
  }
  fit <- fit_selected(original_final, mema_features(original_final), primary_endpoint)
  exclusions[[profile]] <- query_hgf(fit$model, as.matrix(replacement[CellLine == 'AU565', ..programs]), profile)
  qsave(fit$model, file.path(root, 'models', paste0('RNA_source_', profile, '.qs')))
  cat('Source exclusion complete:', profile, '\n')
}
write_table(rbindlist(exclusions), 'HGF_source_study_exclusion_refits')

versions <- c(v128 = '128_signature_modulation_bypass_interaction',
              v131 = '131_precision_weighted_HGF_bridge', v135 = '135_v128_with_corrected_HGF_inputs')
metrics <- scatter <- bars <- list()
for (version in names(versions)) {
  folder <- file.path('/work/results/program_pilot', versions[[version]], 'linear/data')
  d <- fread(file.path(folder, 'HER2_external_held_cell_correlations.tsv'))[gene == 'NRG1']
  d[, version := version]
  metrics[[version]] <- d
  d <- fread(file.path(folder, 'external_held_cell_predictions_frozen.tsv'))[
    study == 'elife32271' & gene == 'NRG1']
  d[, version := version]
  scatter[[version]] <- d
  d <- fread(file.path(folder, 'au565_hgf_5uM_predictions_frozen.tsv'))[
    hgf_context == 'Study-balanced HGF (same recipe as training)' & assumed_HGF_route == 1]
  bars[[version]] <- data.table(version, treatment = c('Lapatinib', 'Lapatinib + HGF'),
    abundance_pct = c(d$lapatinib_pct_control, d$combination_pct_control))
}
metrics <- rbindlist(metrics)
write_table(metrics, 'comparison_NRG1_held_cell_correlations')
write_table(rbindlist(bars), 'comparison_AU565_HGF_predictions')
scatter <- rbindlist(scatter)
scatter[, panel := factor(paste0(dose, ' µM'), levels = c('0.5 µM', '5 µM'))]
annotations <- scatter[, .(label = sprintf('r = %.2f; n = %d',
  cor(observed_departure_log2, predicted_departure_log2), .N)), by = .(version, panel)]
p <- ggplot(scatter, aes(observed_departure_log2, predicted_departure_log2)) +
  geom_point(colour = '#0072B2', size = 1.6) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .4, colour = 'grey35') +
  geom_text(data = annotations, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -.03, vjust = 1.1, size = 2.7) +
  facet_wrap(vars(panel, version), ncol = 3, scales = 'free') + theme_bw(base_size = 9) +
  theme(aspect.ratio = 1) + labs(x = 'Measured interaction (log2)', y = 'Predicted interaction (log2)',
    title = 'NRG1 + lapatinib: held-out cell predictions',
    caption = 'Six cell lines per dose. Previously examined outcomes; retrospective comparison.')
ggsave(file.path(root, 'figures/NRG1_held_cell_version_comparison.pdf'), p, width = 7.1, height = 5.3, device = cairo_pdf)
ggsave(file.path(root, 'figures/NRG1_held_cell_version_comparison.png'), p, width = 7.1, height = 5.3, dpi = 150)
print(metrics)
print(rbindlist(bars))
print(rbindlist(exclusions))
