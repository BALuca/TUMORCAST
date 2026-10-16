#!/usr/bin/env Rscript
# Explain the existing partner holdout and high count correlations. Audit source
# overlaps and replay a lapatinib fit after changing only its excluded outcomes.
# This diagnostic does not select models or change the completed predictions.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(digest)
})
source('/work/src/program_pilot/nonlinear_response.R')
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
data_dir <- file.path(root, 'data')
figure_dir <- file.path(root, 'figures')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(61013L)
obs <- fread(file.path(data_dir, 'matched_molecular_links.tsv.gz'))[
  time_h == 3 & molecular_supported == TRUE]
setorder(obs, partner, anchor, biological_replicate, nominal_A_uM, nominal_B_uM)
bank <- qread(file.path(data_dir, 'representation_scores.qs'))
primary <- fread(file.path(data_dir, 'predictions.tsv.gz'))[
  time_h == 3 & split == 'held_partner']
stopifnot(nrow(obs) == 496L, !anyDuplicated(obs$observation_id))

# Physical-well provenance is checked independently of the fitted metadata.
sample_ids <- function(x) unique(unlist(strsplit(x, ';', fixed = TRUE)))
audit <- rbindlist(lapply(sort(unique(obs$partner)), function(held) {
  tr <- obs[partner != held]
  te <- obs[partner == held]
  train_combo <- sample_ids(tr$source_sample_ids_AB)
  held_combo <- sample_ids(te$source_sample_ids_AB)
  train_singles <- sample_ids(c(tr$source_sample_ids_A, tr$source_sample_ids_B))
  held_singles <- sample_ids(c(te$source_sample_ids_A, te$source_sample_ids_B))
  data.table(held_partner = held, training_pairs = uniqueN(tr$pair),
    training_contexts = uniqueN(tr$condition_id), training_observations = nrow(tr),
    held_pairs = uniqueN(te$pair), held_contexts = uniqueN(te$condition_id), held_observations = nrow(te),
    shared_train_held_combination_wells = length(intersect(train_combo, held_combo)),
    held_combination_wells_in_any_single_input = length(intersect(held_combo, c(train_singles, held_singles))),
    shared_individual_function_wells = length(intersect(train_singles, held_singles)),
    shared_individual_RNA_profiles = length(intersect(c(tr$profile_A, tr$profile_B), c(te$profile_A, te$profile_B))),
    shared_biological_experiments = length(intersect(tr$biological_replicate, te$biological_replicate)))
}))
stopifnot(all(audit$shared_train_held_combination_wells == 0L),
  all(audit$held_combination_wells_in_any_single_input == 0L))
fwrite(audit, file.path(data_dir, 'leakage_source_overlap_audit.tsv'), sep = '\t')

# Rebuild the exact primary input design: combined outcomes are not input fields.
make_design <- function(z) {
  controls <- cbind(single_A = z$ratio_A, single_B = z$ratio_B,
    single_product = z$ratio_A * z$ratio_B, single_A_sq = z$ratio_A^2, single_B_sq = z$ratio_B^2,
    logdose_A = log10(z$nominal_A_uM), logdose_B = log10(z$nominal_B_uM),
    model.matrix(~anchor + factor(biological_replicate), z)[, -1, drop = FALSE])
  a <- bank$original40[z$profile_A, , drop = FALSE]
  b <- bank$original40[z$profile_B, , drop = FALSE]
  sum_response <- a + b
  difference <- abs(a - b)
  product <- a * b
  colnames(sum_response) <- paste0('sum_', colnames(a))
  colnames(difference) <- paste0('absdifference_', colnames(a))
  colnames(product) <- paste0('product_', colnames(a))
  cbind(controls, sum_response, difference, product)
}
held <- 'lapatinib'
tr <- which(obs$partner != held)
te <- which(obs$partner == held)
x <- make_design(obs)
saved <- qread(file.path(root, 'models/3h_lapatinib_Original40_SVR.qs'))
stopifnot(identical(colnames(x), saved$columns), !held %in% saved$training_partners)
y <- obs$excess_inhibition
fit_original <- fit_program_svr(x[tr, , drop = FALSE], y[tr], 'radial',
  saved$selection$cost, saved$selection$gamma_factor / ncol(x))
p_original <- predict_program_svr(fit_original, x[te, , drop = FALSE])
p_saved <- predict_program_svr(saved$fit, x[te, , drop = FALSE])
stopifnot(max(abs(p_original - p_saved)) < 1e-8)

# Intervention: replace every excluded combination count with a different number
# and recompute its target. Preserve individual measurements and all training rows.
altered <- copy(obs)
altered[te, ratio_AB := seq(0.05, 2, length.out = length(te))]
altered[te, excess_inhibition := 100 * (ratio_A * ratio_B - ratio_AB)]
x_altered <- make_design(altered)
stopifnot(identical(x, x_altered), identical(altered$excess_inhibition[tr], y[tr]))
fit_altered <- fit_program_svr(x_altered[tr, , drop = FALSE], altered$excess_inhibition[tr], 'radial',
  saved$selection$cost, saved$selection$gamma_factor / ncol(x))
p_altered <- predict_program_svr(fit_altered, x_altered[te, , drop = FALSE])
change <- max(abs(p_original - p_altered))
stopifnot(change < 1e-8)
fwrite(data.table(held_partner = held, held_observations_changed = length(te),
  maximum_saved_model_replay_difference_pp = max(abs(p_original - p_saved)),
  maximum_prediction_change_after_held_outcome_replacement_pp = change,
  interpretation = 'Frozen-hyperparameter refit unchanged when only held combination outcomes change'),
  file.path(data_dir, 'leakage_held_outcome_intervention.tsv'), sep = '\t')

# Biological means retain the mean of per-experiment products, not a product of means.
all_counts <- primary[model %chin% c('Independent_singles', 'Singles_SVR', 'Original40_SVR'),
  .(observed_pct = 100 * mean(relative_AB), predicted_pct = 100 * mean(predicted_relative_AB),
    repeats = uniqueN(biological_replicate)),
  by = .(model, anchor, partner, nominal_A_uM, nominal_B_uM)]
count_labels <- c(Independent_singles = 'Multiply measured singles\nNO MODEL TRAINED',
  Singles_SVR = 'Learned correction\nSingle function; no RNA',
  Original40_SVR = 'Learned correction\nSingle function + original40 RNA')
all_counts[, design := factor(model, levels = names(count_labels), labels = unname(count_labels))]
count_stats <- all_counts[, .(r = cor(observed_pct, predicted_pct), contexts = .N), by = design]
count_stats[, label := sprintf('Pearson r = %.3f\nn = %d dose conditions', r, contexts)]
fwrite(all_counts, file.path(data_dir, 'leakage_walkthrough_count_plot.tsv.gz'), sep = '\t')
fwrite(count_stats, file.path(data_dir, 'leakage_walkthrough_count_summary.tsv'), sep = '\t')
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  aspect.ratio = 1, legend.position = 'bottom', plot.title.position = 'plot'))
anchor_colors <- c(alpelisib = '#0072B2', dactolisib = '#D55E00')
count_plot <- ggplot(all_counts, aes(observed_pct, predicted_pct)) +
  geom_point(aes(color = anchor), alpha = .72, size = 1.8) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .65) +
  geom_text(data = count_stats, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.03, vjust = 1.2, size = 3.8) +
  facet_wrap(~design, nrow = 1, scales = 'free') + scale_color_manual(values = anchor_colors) +
  labs(title = 'A high combination-count correlation already exists before any learning',
    subtitle = 'All six partner exclusions pooled; 3 h primary model; observed individual function is supplied',
    x = 'Measured combination cell count (% of vehicle)', y = 'Predicted combination cell count (% of vehicle)',
    color = 'First drug', caption = paste('198 conditions, averaged across 2–3 biological experiments.',
      'The left panel uses only measured individual-drug fractions A × B.',
      'Lines are fitted trends; axes have independent scales.'))
ggsave(file.path(figure_dir, 'leakage_walkthrough_total_counts.pdf'), count_plot,
  width = 11.5, height = 5.2, device = cairo_pdf)
ggsave(file.path(figure_dir, 'leakage_walkthrough_total_counts.png'), count_plot, width = 11.5, height = 5.2, dpi = 160)

# The fitted primary target is extra inhibition, displayed separately from counts.
effect <- primary[model %chin% c('Singles_SVR', 'Original40_SVR'),
  .(observed_pp = mean(observed_extra_pp), predicted_pp = mean(predicted_extra_pp)),
  by = .(model, anchor, partner, nominal_A_uM, nominal_B_uM)]
effect[, design := factor(model, levels = c('Singles_SVR', 'Original40_SVR'),
  labels = c('Single function; no RNA', 'Single function + original40 RNA'))]
effect_stats <- effect[, .(r = cor(observed_pp, predicted_pp), contexts = .N), by = design]
effect_stats[, label := sprintf('Pearson r = %.3f\nn = %d dose conditions', r, contexts)]
fwrite(effect, file.path(data_dir, 'leakage_walkthrough_learned_target_plot.tsv.gz'), sep = '\t')
effect_plot <- ggplot(effect, aes(observed_pp, predicted_pp)) +
  geom_point(aes(color = anchor), alpha = .72, size = 1.8) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .65) +
  geom_hline(yintercept = 0, color = 'grey75', linewidth = .3) +
  geom_vline(xintercept = 0, color = 'grey75', linewidth = .3) +
  geom_text(data = effect_stats, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.03, vjust = 1.2, size = 3.8) +
  facet_wrap(~design, nrow = 1, scales = 'free') + scale_color_manual(values = anchor_colors) +
  labs(title = 'What the primary regression actually learned: departure from the single-drug expectation',
    subtitle = 'Target = 100 × (measured single A × measured single B − measured combination)',
    x = 'Measured additional inhibition (percentage points)',
    y = 'Predicted additional inhibition (percentage points)', color = 'First drug',
    caption = 'The RNA model supplies a predicted correction to the measured-single product; adding RNA does not improve this test.')
ggsave(file.path(figure_dir, 'leakage_walkthrough_learned_target.pdf'), effect_plot,
  width = 9.0, height = 5.6, device = cairo_pdf)
ggsave(file.path(figure_dir, 'leakage_walkthrough_learned_target.png'), effect_plot, width = 9.0, height = 5.6, dpi = 160)

writeLines(c('Leakage walkthrough: source-level overlap, held-outcome intervention, and endpoint decomposition.',
  'No held combination well is used in training combinations or any individual-function input.',
  'Individual first-drug measurements, RNA profiles and biological experiments are deliberately shared.',
  'Only combination outcomes for the excluded partner are withheld. This is neither new-cell nor new-study validation.',
  paste('Lapatinib saved-model replay discrepancy (pp):', max(abs(p_original - p_saved))),
  paste('Prediction change after replacing all64 lapatinib combination outcomes (pp):', change),
  'The intervention reuses frozen hyperparameters and checks the exact final fit; it does not independently rerun tuning.',
  'Primary predicted count = measured singleA × measured singleB − predicted additional inhibition/100.',
  'Thus observed individual functional outcomes contribute most of the high total-count correlation.',
  'This is permitted under an individuals-observed/pair-unobserved task, but not evidence of a molecular-only or baseline-only bridge.'),
  file.path(root, 'LEAKAGE_WALKTHROUGH.txt'))
capture.output(sessionInfo(), file = file.path(root, 'leakage_walkthrough_sessionInfo.txt'))
message('Completed leakage provenance, frozen-model intervention, and walkthrough plots.')
