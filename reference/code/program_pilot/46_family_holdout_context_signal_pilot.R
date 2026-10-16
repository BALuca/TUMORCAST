#!/usr/bin/env Rscript
# Retrospective decoder-family transfer using the fixed43 hard-mean program bank.
# Cover every original ERBB/MEK family member; missing assays remain explicit.
suppressPackageStartupMessages({library(data.table); library(qs); library(Matrix); library(ggplot2); library(patchwork)})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/nonlinear_response.R')
source('/work/src/program_pilot/functional_programs.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
set.seed(20261010L)
options(warn = 1)

root <- '/work/results/program_pilot/46_family_holdout_context_signal_pilot'
input_path <- '/work/results/program_pilot/44_compare_program_banks_molecular_function_svr/data/frozen_comparison_inputs.qs'
bank_path <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
fold_path <- '/work/results/program_pilot/45_whole_cell_regularized_program_svr/data/whole_cell_folds.tsv'
grid_path <- '/work/results/program_pilot/45_whole_cell_regularized_program_svr/data/regularization_grid.tsv'
families <- c('ERBB', 'MEK')
workers <- 3L
# Descriptive contrast flag uses existing-cell assays only; no context is hidden.
minimum_cells <- 18L
minimum_per_class <- 4L
minimum_iqr_pct <- 25
for (p in file.path(root, c('data', 'models', 'figures'))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading43 bank and original family definitions')
input <- qread(input_path)
bank <- qread(bank_path)
panel <- qread(panel_path)
md <- copy(input$metadata)
x <- input$scores[['43 averaged / 34 programs']]
folds <- fread(fold_path)
grid <- fread(grid_path)
training_cells <- sort(unique(md[split == 'Training cells', depmap_id]))
held_cells <- input$held_ids
specifications <- lapply(families, function(family) qread(paste0(
  '/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout/settings.qs')))
names(specifications) <- families
family_map <- rbindlist(lapply(families, function(family) {
  s <- specifications[[family]]
  check(setequal(s$held_ids, held_cells), 'Original family split differs from44')
  data.table(family = family, drug = s$excluded_drugs)
}))
check(!anyDuplicated(family_map$drug) && length(unique(bank$partition)) == 34L, 'Invalid family map or43 bank')
fwrite(family_map, file.path(root, 'data/excluded_families.tsv'), sep = '\t')
md[, fold := folds$fold[match(depmap_id, folds$depmap_id)]]
source_files <- c(input_path, bank_path, panel_path, fold_path, grid_path,
  paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', families, '_holdout/settings.qs'))
hashes <- data.table(path = source_files,
  sha256 = vapply(source_files, digest::digest, character(1), file = TRUE, algo = 'sha256'))
signature <- digest::digest(list(hashes, input$input_hashes, family_map), algo = 'sha256')
fwrite(hashes, file.path(root, 'data/input_hashes.tsv'), sep = '\t')

# Score all17 compounds, including those outside43's original105 plotted drugs.
projection_path <- file.path(root, 'data/family_projection_inputs.qs')
if (file.exists(projection_path)) {
  projection <- qread(projection_path)
  check(identical(projection$signature, signature), 'Family projection inputs changed')
  fm <- projection$metadata
  fx <- projection$scores
} else {
  fm <- copy(panel$metadata[drug %chin% family_map$drug])
  ann <- panel$sensitivity[match(fm$observation_id, observation_id)]
  check(identical(ann$observation_id, fm$observation_id), 'Family molecular/assay rows do not align')
  fm[, `:=`(family = family_map$family[match(drug, family_map$drug)],
    measured_pct = ann$prism_abundance_pct, almanac_pct = ann$almanac_abundance_pct,
    recipient_set = fifelse(depmap_id %chin% training_cells, 'Existing cells',
      fifelse(depmap_id %chin% held_cells, 'Held cells', 'Other cells')))]
  fx <- matrix(NA_real_, nrow(fm), ncol(x), dimnames = list(fm$observation_id, colnames(x)))
  previous <- match(fm$observation_id, bank$metadata$observation_id)
  fx[!is.na(previous), ] <- bank$scores[previous[!is.na(previous)], , drop = FALSE]
  for (id in unique(fm$depmap_id[is.na(previous)])) {
    ii <- which(fm$depmap_id == id & is.na(previous))
    molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
    fx[ii, ] <- as.matrix(molecular$response[fm$observation_id[ii], bank$genes, drop = FALSE] %*% bank$weights)
  }
  check(all(is.finite(fx)), 'Missing/nonfinite family program scores')
  qsave(list(signature = signature, metadata = fm, scores = fx), projection_path, preset = 'high')
}
coverage <- fm[, .(molecular_columns = .N, molecular_cells = uniqueN(depmap_id),
  prism_columns = sum(is.finite(measured_pct) & measured_pct > 0),
  almanac_columns = sum(is.finite(almanac_pct))), by = .(family, drug)]
coverage[, assay_status := fifelse(prism_columns > 0, 'Matched PRISM available', 'No matched PRISM: molecular projection only')]
fwrite(coverage, file.path(root, 'data/all_family_drug_coverage.tsv'), sep = '\t')
status(paste('all', nrow(coverage), 'family compounds projected;', sum(coverage$prism_columns > 0), 'have matched PRISM'))

safe_cor <- function(a, b) if (length(a) >= 3L && sd(a) > 1e-12 && sd(b) > 1e-12) cor(a, b) else NA_real_
# Identical same-treatment validation criterion to45, now excluding each family.
cv_measure <- function(z, pred) {
  z <- copy(z)
  z[, predicted := pred]
  z[, `:=`(yc = observed_log2 - mean(observed_log2), pc = predicted - mean(predicted)), by = .(condition, fold)]
  data.table(pooled_r = safe_cor(z$observed_log2, z$predicted), within_treatment_r = safe_cor(z$yc, z$pc))
}
cv_scores <- selected <- predictions <- list()
for (family in families) {
  excluded <- specifications[[family]]$excluded_drugs
  train <- which(md$depmap_id %chin% training_cells & !md$drug %chin% excluded)
  check(!any(md$depmap_id[train] %chin% held_cells) && !any(md$drug[train] %chin% excluded), 'Decoder fitting leakage')
  check(uniqueN(md$depmap_id[train]) == 26L && all(!is.na(md$fold[train])), 'Invalid training recipients/folds')
  metrics <- list()
  for (j in seq_len(nrow(grid))) {
    setting <- grid[j]
    status(paste(family, 'whole-cell CV candidate', j, '/', nrow(grid)))
    results <- parallel::mclapply(1:3, function(fold_id) {
      RhpcBLASctl::blas_set_num_threads(1L)
      fit_rows <- train[md$fold[train] != fold_id]
      val_rows <- train[md$fold[train] == fold_id]
      path <- file.path(root, 'models', sprintf('%s_candidate_%02d_fold_%d.qs', family, j, fold_id))
      if (file.exists(path)) {
        s <- qread(path)
        check(identical(s$signature, signature), 'Stale CV checkpoint')
        return(s)
      }
      check(!any(md$depmap_id[fit_rows] %chin% md$depmap_id[val_rows]), 'Validation cells entered fitting')
      f <- fit_program_svr(x[fit_rows, , drop = FALSE], md$observed_log2[fit_rows],
        'radial', setting$cost, setting$gamma_multiplier / ncol(x))
      s <- list(signature = signature, fitted = f, setting = setting,
        fitting_ids = md$observation_id[fit_rows], fitting_cells = unique(md$depmap_id[fit_rows]),
        validation_ids = md$observation_id[val_rows], validation_cells = unique(md$depmap_id[val_rows]),
        prediction = predict_program_svr(f, x[val_rows, , drop = FALSE]), excluded_drugs = excluded)
      qsave(s, path, preset = 'high')
      s
    }, mc.cores = workers, mc.preschedule = FALSE)
    check(all(vapply(results, function(s) is.list(s) && !inherits(s, 'try-error'), logical(1))), 'CV worker failed')
    oof <- rep(NA_real_, length(train))
    for (s in results) oof[match(s$validation_ids, md$observation_id[train])] <- s$prediction
    check(all(is.finite(oof)), 'Incomplete CV predictions')
    score <- cbind(data.table(family = family, candidate = j, cost = setting$cost,
      gamma_multiplier = setting$gamma_multiplier), cv_measure(md[train], oof))
    metrics[[j]] <- score
    cv_scores[[paste(family, j)]] <- score
    fwrite(rbindlist(cv_scores), file.path(root, 'data/whole_cell_validation.tsv'), sep = '\t')
  }
  chosen <- rbindlist(metrics)[order(-within_treatment_r, cost, gamma_multiplier)][1]
  selected[[family]] <- chosen
  fwrite(rbindlist(selected), file.path(root, 'data/selected_settings.tsv'), sep = '\t')
  status(paste(family, 'SELECTED C', chosen$cost, 'gamma', chosen$gamma_multiplier, '/p'))
  path <- file.path(root, 'models', paste0(family, '_selected_decoder.qs'))
  if (file.exists(path)) {
    s <- qread(path)
    check(identical(s$signature, signature) && identical(s$chosen$candidate, chosen$candidate), 'Stale final decoder')
    f <- s$fitted
  } else {
    f <- fit_program_svr(x[train, , drop = FALSE], md$observed_log2[train],
      'radial', chosen$cost, chosen$gamma_multiplier / ncol(x))
    qsave(list(signature = signature, fitted = f, chosen = chosen, training_ids = md$observation_id[train],
      excluded_drugs = excluded, excluded_cells = held_cells), path, preset = 'high')
  }
  ii <- which(fm$family == family)
  pred <- predict_program_svr(f, fx[ii, , drop = FALSE])
  p <- fm[ii, .(observation_id, family, drug, dose, depmap_id, recipient, recipient_set)]
  p[, `:=`(predicted_log2 = pred, predicted_pct = 100 * 2^pred)]
  check(all(is.finite(p$predicted_pct)), 'Nonfinite family prediction')
  fwrite(p, file.path(root, 'data', paste0(family, '_outcome_free_predictions.tsv')), sep = '\t')
  fwrite(data.table(sha256 = digest::digest(p, algo = 'sha256')),
    file.path(root, 'data', paste0(family, '_prediction_hash.tsv')), sep = '\t')
  predictions[[family]] <- p
  status(paste(family, 'MODEL AND ALL FAMILY PREDICTIONS SAVED'))
}
z <- rbindlist(predictions)
z[, measured_pct := fm$measured_pct[match(observation_id, fm$observation_id)]]
fwrite(z, file.path(root, 'data/evaluated_family_predictions.tsv'), sep = '\t')

# Report every setting, with an explicit contrast flag and sensitive/resistant counts.
context <- z[recipient_set != 'Other cells', {
  v <- .SD[is.finite(measured_pct) & measured_pct > 0]
  n <- nrow(v)
  r <- s <- integer()
  if (n) {
    r <- which(v$measured_pct >= 60)
    s <- which(v$measured_pct <= 40)
  }
  pair_count <- 0L
  ordering <- NA_real_
  if (length(r) && length(s)) {
    valid <- outer(v$measured_pct[r], v$measured_pct[s], '-') >= 25
    gap <- outer(v$predicted_log2[r], v$predicted_log2[s], '-')
    pair_count <- sum(valid)
    if (pair_count) ordering <- mean((gap[valid] > 0) + .5 * (gap[valid] == 0))
  }
  .(measured_cells = n, sensitive = sum(v$measured_pct <= 50), resistant = sum(v$measured_pct > 50),
    clearly_sensitive = length(s), clearly_resistant = length(r),
    abundance_iqr_pct = if (n) IQR(v$measured_pct) else NA_real_,
    pearson_log2 = if (n >= 3) safe_cor(log2(v$measured_pct / 100), v$predicted_log2) else NA_real_,
    pearson_pct = if (n >= 3) safe_cor(v$measured_pct, v$predicted_pct) else NA_real_,
    separated_pairs = pair_count, correct_pair_order_fraction = ordering)
}, by = .(family, drug, dose, recipient_set)]
flags <- context[recipient_set == 'Existing cells', .(family, drug, dose,
  existing_cell_contrast = measured_cells >= minimum_cells & sensitive >= minimum_per_class &
    resistant >= minimum_per_class & abundance_iqr_pct >= minimum_iqr_pct)]
context <- merge(context, flags, by = c('family', 'drug', 'dose'), all.x = TRUE)
fwrite(context, file.path(root, 'data/all_drug_dose_results.tsv'), sep = '\t')
v <- z[is.finite(measured_pct) & measured_pct > 0 & recipient_set != 'Other cells']
v[, observed_log2 := log2(measured_pct / 100)]
v[, `:=`(yc = observed_log2 - mean(observed_log2), pc = predicted_log2 - mean(predicted_log2)),
  by = .(family, drug, dose, recipient_set)]
drug_summary <- v[, .(observations = .N, cells = uniqueN(depmap_id),
  pooled_pearson_log2 = safe_cor(observed_log2, predicted_log2),
  within_dose_pearson_log2 = safe_cor(yc, pc)), by = .(family, drug, recipient_set)]
fwrite(drug_summary, file.path(root, 'data/all_drug_results.tsv'), sep = '\t')
family_summary <- v[, .(observations = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug),
  pooled_pearson_log2 = safe_cor(observed_log2, predicted_log2),
  within_drug_dose_pearson_log2 = safe_cor(yc, pc)), by = .(family, recipient_set)]
fwrite(family_summary, file.path(root, 'data/family_summary.tsv'), sep = '\t')
print(family_summary)
status('all17 compounds evaluated; rendering all-drug figures')

# Matched rows across two heatmaps: correlations and sensitive/resistant ordering.
order_table <- unique(context[, .(family, drug, dose, existing_cell_contrast)])[order(family, drug, dose)]
order_table[, label := paste0(drug, ' | ', dose, ' uM')]
row_annotation <- data.frame(Family = order_table$family,
  Contrast = ifelse(order_table$existing_cell_contrast, 'Present in existing cells', 'Limited / no assay'),
  row.names = order_table$label)
matrices <- lapply(c('pearson_log2', 'correct_pair_order_fraction'), function(metric) {
  m <- matrix(NA_real_, nrow(order_table), 2, dimnames = list(order_table$label, c('Existing cells', 'Held cells')))
  for (j in seq_len(2)) {
    q <- context[recipient_set == colnames(m)[j]]
    ii <- match(paste(q$drug, q$dose), paste(order_table$drug, order_table$dose))
    m[ii, j] <- q[[metric]]
  }
  m
})
palette <- list(Family = c(ERBB = '#B24C35', MEK = '#27818A'),
  Contrast = c('Present in existing cells' = '#5B8A3C', 'Limited / no assay' = '#D9D9D9'))
hm1 <- draw_heatmap(matrices[[1]], scale_rows = FALSE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7),
  row_anno_df = row_annotation, row_anno_cols = c('Family', 'Contrast'), palette_overrides = palette,
  ann_row_args = list(show_annotation_name = FALSE), show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 8), column_names_rot = 45,
  value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'), na_col = '#D9D9D9',
  hm_name = 'Prediction r', title = 'Fixed drug/dose Pearson r\nlog2 viability; all family compounds')
hm2 <- draw_heatmap(matrices[[2]], scale_rows = FALSE, show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 7), row_names_side = 'right', row_names_max_width = grid::unit(7, 'cm'),
  show_column_names = TRUE, column_names_gp = grid::gpar(fontsize = 8), column_names_rot = 45,
  value_range = c(0, .5, 1), value_colors = c('#2166AC', 'white', '#B2182B'), na_col = '#D9D9D9',
  hm_name = 'Correct order', title = 'Clear sensitive/resistant pairs\nGrey: no eligible pairs / no assay')
save_heatmap(hm1 + hm2, 'all_family_context_results', 10, max(6, .15 * nrow(order_table) + 1.7),
  file.path(root, 'figures'))

# A labeled point plot exposes every measured context and its class balance.
wide <- dcast(context, family + drug + dose ~ recipient_set,
  value.var = c('pearson_log2', 'sensitive', 'resistant', 'measured_cells'))
wide <- wide[get('measured_cells_Existing cells') > 0 | get('measured_cells_Held cells') > 0]
wide[, label := sprintf('%s | %g uM   [S/R %d/%d ; %d/%d]', drug, dose,
  get('sensitive_Existing cells'), get('resistant_Existing cells'),
  get('sensitive_Held cells'), get('resistant_Held cells'))]
long <- merge(context, wide[, .(family, drug, dose, label)], by = c('family', 'drug', 'dose'))
long[, label := factor(label, levels = rev(wide[order(family, drug, dose), label]))]
wide[, label := factor(label, levels = levels(long$label))]
p <- ggplot(long, aes(pearson_log2, label, colour = recipient_set)) +
  geom_segment(data = wide, aes(x = get('pearson_log2_Existing cells'), xend = get('pearson_log2_Held cells'),
    y = label, yend = label), inherit.aes = FALSE, colour = 'grey75', linewidth = .5) +
  geom_vline(xintercept = 0, colour = 'grey55', linetype = 2) + geom_point(size = 2.2) +
  facet_grid(family ~ ., scales = 'free_y', space = 'free_y') +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
  scale_colour_manual(values = c('Existing cells' = '#27818A', 'Held cells' = '#B24C35')) +
  labs(title = 'Does the response-to-function map transfer to an excluded family?',
    subtitle = 'Every drug/dose with PRISM. S/R counts: existing cells ; held cells (50% vehicle cutoff).\nEntire family excluded from decoder fitting and selection; fixed43 program bank.',
    x = 'Prediction Pearson r across cells at the same drug and dose (log2 viability)', y = NULL, colour = NULL,
    caption = 'Seven compounds have no matched assay; see coverage table. Program discovery used all outcomes: retrospective comparison.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9))
ggsave(file.path(root, 'figures/all_assayed_contexts_readable.pdf'), p, width = 11, height = 8.5)
ggsave(file.path(root, 'figures/all_assayed_contexts_readable.png'), p, width = 11, height = 8.5, dpi = 160)

# Lead with the full, dose-inclusive prediction task; show the narrower question beside it.
drug_plot_data <- melt(drug_summary, id.vars = c('family', 'drug', 'recipient_set', 'observations', 'cells'),
  measure.vars = c('pooled_pearson_log2', 'within_dose_pearson_log2'),
  variable.name = 'comparison', value.name = 'pearson_r')
drug_plot_data[, comparison := factor(comparison,
  levels = c('pooled_pearson_log2', 'within_dose_pearson_log2'),
  labels = c('All doses together', 'Cell differences after subtracting dose averages'))]
drug_plot_data[, drug_label := paste(family, drug, sep = ' | ')]
drug_order <- unique(drug_plot_data[order(family, drug), drug_label])
drug_plot_data[, drug_label := factor(drug_label, levels = rev(drug_order))]
overview <- ggplot(drug_plot_data, aes(pearson_r, drug_label, colour = recipient_set)) +
  geom_vline(xintercept = 0, colour = 'grey65', linetype = 2) +
  geom_point(position = position_dodge(width = .4), size = 2.5) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
  scale_colour_manual(values = c('Existing cells' = '#27818A', 'Held cells' = '#B24C35'),
    labels = c('Existing cells' = '26 familiar cells; excluded drugs', 'Held cells' = '9 held-out cells; excluded drugs')) +
  labs(title = 'Measured program displacements predict some excluded-drug responses',
    subtitle = 'Left: all three doses contribute. Right: the additional distinction among cells at a given dose.\nAll assayed ERBB/MEK drugs; each family excluded from its decoder fitting and selection.',
    x = 'Pearson r: predicted versus measured log2 viability', y = NULL, colour = NULL,
    caption = 'Input: 34 measured program displacements. Seven other compounds lack matched PRISM.\nFixed bank43 was discovered using all outcomes: this is retrospective decoder transfer, not an independent full-pipeline test.') +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom', axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9),
    strip.text = element_text(size = 9), plot.caption = element_text(hjust = 0, size = 8))
ggsave(file.path(root, 'figures/all_drugs_dose_inclusive_overview.pdf'), overview, width = 11.5, height = 5.3)
ggsave(file.path(root, 'figures/all_drugs_dose_inclusive_overview.png'), overview, width = 11.5, height = 5.3, dpi = 160)

# A single multipage PDF connects the overview to actual predictions, dose coloured.
# Plot percentages on log2 axes so strong killing and near-vehicle responses both remain visible.
family_pages <- list()
dose_palette <- c('0.05' = '#009E73', '0.5' = '#E69F00', '5' = '#0072B2')
for (family_name in families) {
  d <- copy(v[family == family_name])
  d[, dose_label := factor(as.character(dose), levels = names(dose_palette))]
  drugs <- sort(unique(d$drug))
  panel_keys <- as.vector(rbind(paste(drugs, 'Existing cells'), paste(drugs, 'Held cells')))
  labels <- drug_summary[family == family_name]
  labels[, panel_key := paste(drug, recipient_set)]
  labels[, panel_label := sprintf('%s | %s\nAll-dose r = %.2f; %d cells, %d dose observations',
    drug, fifelse(recipient_set == 'Existing cells', 'familiar cells', 'held-out cells'),
    pooled_pearson_log2, cells, observations)]
  d[, panel := factor(labels$panel_label[match(paste(drug, recipient_set), labels$panel_key)],
    levels = labels$panel_label[match(panel_keys, labels$panel_key)])]
  p <- ggplot(d, aes(measured_pct, predicted_pct)) +
    geom_hline(yintercept = 50, colour = 'grey80', linetype = 2) +
    geom_vline(xintercept = 50, colour = 'grey80', linetype = 2) +
    geom_point(aes(colour = dose_label), size = 1.5, alpha = .8) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = 'grey35', linewidth = .55) +
    scale_x_continuous(trans = 'log2', breaks = c(1, 5, 10, 25, 50, 100, 200)) +
    scale_y_continuous(trans = 'log2', breaks = c(1, 5, 10, 25, 50, 100, 200)) +
    scale_colour_manual(values = dose_palette, drop = FALSE) +
    facet_wrap(~ panel, ncol = if (family_name == 'ERBB') 4 else 2, scales = 'free') +
    labs(title = paste(family_name, '| all doses together, every assayed drug'),
      subtitle = 'Measured 34-program displacements to viability. Entire family excluded from decoder fitting/selection.\nDoses remain independent observations; dashed lines mark 50% vehicle.',
      x = 'Measured PRISM % vehicle (log2 scale)', y = 'Predicted % vehicle (log2 scale)', colour = 'Dose (uM)',
      caption = 'Panels use independent x/y ranges. Correlations are descriptive; repeated doses are not independent cells.\nBank43 used all outcomes during discovery. Baseline-to-displacement is not tested here.') +
    theme_classic(base_size = 8) + theme(aspect.ratio = 1, legend.position = 'bottom',
      strip.text = element_text(size = 7), plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 9), plot.caption = element_text(hjust = 0, size = 8))
  family_pages[[family_name]] <- p
  page_width <- if (family_name == 'ERBB') 12 else 8
  page_height <- if (family_name == 'ERBB') 12 else 10
  ggsave(file.path(root, 'figures', paste0(family_name, '_all_drugs_dose_inclusive.pdf')),
    p, width = page_width, height = page_height)
  ggsave(file.path(root, 'figures', paste0(family_name, '_all_drugs_dose_inclusive.png')),
    p, width = page_width, height = page_height, dpi = 160)
}
grDevices::cairo_pdf(file.path(root, 'figures/combined_family_transfer.pdf'), width = 11.5, height = 11.5)
print(overview)
for (p in family_pages) print(p)
grDevices::dev.off()

# One six-panel prediction page and one34-program map for every compound.
program_order <- unique(bank$partition[bank$tree$labels[bank$tree$order]])
programs <- sprintf('P%02d', program_order)
row_meta <- data.frame(Program = factor(programs, levels = programs), row.names = programs)
program_colors <- setNames(hcl.colors(34, 'Dark 3'), sprintf('P%02d', 1:34))
for (k in seq_len(nrow(coverage))) {
  family <- coverage$family[k]
  drug_name <- coverage$drug[k]
  slug <- paste0(family, '_', gsub('[^A-Za-z0-9]+', '_', drug_name))
  q <- z[drug == drug_name]
  plots <- list()
  for (recipient_group in c('Existing cells', 'Held cells')) for (dose_value in sort(unique(q$dose))) {
    d <- q[recipient_set == recipient_group & dose == dose_value & is.finite(measured_pct) & measured_pct > 0]
    m <- context[drug == drug_name & recipient_set == recipient_group & dose == dose_value]
    if (nrow(d)) {
      plots[[paste(recipient_group, dose_value)]] <- ggplot(d, aes(log2(measured_pct / 100), predicted_log2)) +
        geom_hline(yintercept = -1, colour = 'grey80', linetype = 2) +
        geom_vline(xintercept = -1, colour = 'grey80', linetype = 2) +
        geom_point(size = 1.7, colour = '#27818A') +
        geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .55, colour = '#B24C35') +
        labs(title = paste(recipient_group, '|', dose_value, 'uM'),
          subtitle = sprintf('r = %.2f | n=%d; S/R=%d/%d', m$pearson_log2, m$measured_cells, m$sensitive, m$resistant),
          x = 'Measured log2(% vehicle /100)', y = 'Predicted log2(% vehicle /100)') +
        theme_classic(base_size = 8) + theme(aspect.ratio = 1, plot.title = element_text(size = 9))
    } else {
      plots[[paste(recipient_group, dose_value)]] <- ggplot() +
        annotate('text', x = 0, y = 0, label = 'No matched PRISM assay', size = 3) +
        labs(title = paste(recipient_group, '|', dose_value, 'uM')) +
        theme_void(base_size = 8) + theme(aspect.ratio = 1)
    }
  }
  p <- wrap_plots(plots, ncol = 3) + plot_annotation(title = paste(family, '|', drug_name),
    subtitle = paste0('Entire ', family, ' family excluded from decoder fitting/selection; bank43 fixed\n',
      'Existing cells: new-family test. Held cells: new-family + new-cell test. Retrospective bank discovery.'))
  ggsave(file.path(root, 'figures', paste0(slug, '_predictions.pdf')), p, width = 10, height = 6.5)
  ggsave(file.path(root, 'figures', paste0(slug, '_predictions.png')), p, width = 10, height = 6.5, dpi = 160)
  ii <- match(q$observation_id, fm$observation_id)
  group <- factor(q$recipient_set, levels = c('Existing cells', 'Held cells', 'Other cells'))
  co <- order(group, q$dose, q$measured_pct, q$recipient, na.last = TRUE)
  plotted <- t(fx[ii[co], programs, drop = FALSE])
  cm <- data.frame(Set = group[co], Dose = factor(q$dose[co]), row.names = q$observation_id[co])
  check(identical(colnames(plotted), rownames(cm)), 'Program map annotations are misaligned')
  observed <- q$measured_pct[co]
  predicted <- q$predicted_pct[co]
  bar_max <- max(125, ceiling(max(c(observed, predicted), na.rm = TRUE) / 25) * 25)
  observed_bar <- observed
  observed_bar[is.na(observed_bar)] <- 0
  top <- list(show_legend = FALSE,
    'PRISM % vehicle' = ComplexHeatmap::anno_barplot(observed_bar, ylim = c(0, bar_max),
      gp = grid::gpar(fill = '#27818A', col = NA), height = grid::unit(13, 'mm')),
    'Predicted % vehicle' = ComplexHeatmap::anno_barplot(predicted, ylim = c(0, bar_max),
      gp = grid::gpar(fill = '#B24C35', col = NA), height = grid::unit(13, 'mm')),
    'PRISM value' = sensitivity_track(observed))
  hm <- draw_heatmap(plotted, scale_rows = FALSE, row_anno_df = row_meta, row_anno_cols = 'Program',
    col_anno_df = cm, col_anno_cols = c('Set', 'Dose'), ann_top_args = top,
    palette_overrides = list(Program = program_colors), ann_row_args = list(show_legend = FALSE),
    column_split = interaction(cm$Set, cm$Dose, drop = TRUE, lex.order = TRUE),
    cluster_column_slices = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7), show_column_names = TRUE,
    column_labels = q$recipient[co], column_names_gp = grid::gpar(fontsize = 5), column_names_rot = 90,
    value_range = c(-.3, 0, .3), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'Mean log2 delta', title = paste0(drug_name, ' | fixed43 programs: all assigned-gene means\n',
      'Columns: recipient set, dose, observed sensitivity. Grey PRISM track = missing assay.'))
  save_heatmap(hm, paste0(slug, '_program_displacements'), 14, 7.5, file.path(root, 'figures'), list(sensitivity_legend))
  status(paste('FIGURES READY:', drug_name, '|', k, '/', nrow(coverage)))
}
writeLines(c('Fixed43bank:34programs,8429genes,hard means of all assigned genes.',
  'Separate ERBB/MEK decoders exclude the entire original family and nine held cells.',
  'SVR setting selected by3 whole-cell folds using only26 existing cells and non-family drugs.',
  'All17 family compounds and51dosecontexts reported;7 compounds have no matched PRISM.',
  'Existing cells are excluded-family evaluation endpoints, not decoder training predictions.',
  'Contrast flag: >=18existingcells,>=4sensitive/4resistant at50%vehicle, IQR>=25percentagepoints.',
  'Flag is descriptive only; no setting is dropped or used to tune its family decoder.',
  'Clear pair ordering:<=40vs>=60%vehicle, separation>=25percentagepoints.',
  'Displacements are measured. Baseline-to-displacement and direct baseline comparison not tested.',
  '43discovery used allcell/drug outcomes; prior test results seen. Retrospective decoder transfer only.'),
  file.path(root, 'README.txt'))
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('complete')
