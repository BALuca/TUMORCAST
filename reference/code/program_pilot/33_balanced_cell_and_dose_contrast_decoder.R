#!/usr/bin/env Rscript
# Learn recipient and dose sensitivity contrasts in the existing frozen programs.
suppressPackageStartupMessages({library(data.table); library(qs); library(glmnet); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

# Fixed signal gates use fitting recipients only; no tuning on outer outcomes.
root <- '/work/results/program_pilot/33_balanced_cell_and_dose_contrast_decoder'
families <- c('ERBB', 'MEK')
minimum_cells <- 18L
minimum_each_class <- 4L
sensitive_limit <- 40
resistant_limit <- 60
minimum_iqr_pp <- 25
minimum_pair_gap_pp <- 25
lambda_grid <- 10^seq(2, -6, by = -1)
objective_mass <- c(Cell = .45, Dose = .45, Anchor = .10)
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(6); quit(save = 'no', status = 1L)})

# Construct real within-context pairs; drug identities never enter the features.
# Row r is clearly resistant; row s is clearly sensitive. Dose pairs stay within
# a recipient and require the sensitive observation to have the higher dose.
contrast_pairs <- function(md, cell_keys, dose_drugs) {
  parts <- list()
  for (key in cell_keys) {
    ii <- which(md$condition == key)
    rr <- ii[md$observed_pct[ii] >= resistant_limit]
    ss <- ii[md$observed_pct[ii] <= sensitive_limit]
    if (!length(rr) || !length(ss)) next
    z <- CJ(r = rr, s = ss)
    z <- z[md$observed_pct[r] - md$observed_pct[s] >= minimum_pair_gap_pp]
    z[, `:=`(axis = 'Cell', unit = key, recipient = 'across recipients')]
    parts[[paste('Cell', key)]] <- z
  }
  keys <- unique(md[drug %in% dose_drugs, .(drug, depmap_id)])
  for (j in seq_len(nrow(keys))) {
    ii <- which(md$drug == keys$drug[j] & md$depmap_id == keys$depmap_id[j])
    rr <- ii[md$observed_pct[ii] >= resistant_limit]
    ss <- ii[md$observed_pct[ii] <= sensitive_limit]
    if (!length(rr) || !length(ss)) next
    z <- CJ(r = rr, s = ss)
    z <- z[md$dose[s] > md$dose[r] & md$observed_pct[r] - md$observed_pct[s] >= minimum_pair_gap_pp]
    z[, `:=`(axis = 'Dose', unit = keys$drug[j], recipient = keys$depmap_id[j])]
    parts[[paste('Dose', keys$drug[j], keys$depmap_id[j])]] <- z
  }
  if (!length(parts)) return(data.table(r = integer(), s = integer(), axis = character(),
    unit = character(), recipient = character()))
  pairs <- rbindlist(parts)
  pairs[, n_pairs := .N, by = .(axis, unit, recipient)]
  pairs[, n_recipients := uniqueN(recipient), by = .(axis, unit)]
  pairs[, base_weight := 1 / (n_pairs * n_recipients)]
  pairs[, base_weight := base_weight / sum(base_weight), by = axis]
  pairs
}

# Censor unsupported cell contrasts and unsupported dose trajectories independently.
screen_signal <- function(md) {
  cells <- md[, .(cells = uniqueN(depmap_id), clear_sensitive = sum(observed_pct <= sensitive_limit),
    clear_resistant = sum(observed_pct >= resistant_limit), iqr_pp = IQR(observed_pct)), by = .(drug, dose, condition)]
  candidates <- cells[cells >= minimum_cells & clear_sensitive >= minimum_each_class &
    clear_resistant >= minimum_each_class & iqr_pp >= minimum_iqr_pp, condition]
  pairs <- contrast_pairs(md, candidates, unique(md$drug))
  participants <- pairs[axis == 'Cell', .(paired_sensitive = uniqueN(s), paired_resistant = uniqueN(r)), by = unit]
  cells <- merge(cells, participants, by.x = 'condition', by.y = 'unit', all.x = TRUE)
  cells[is.na(paired_sensitive), `:=`(paired_sensitive = 0L, paired_resistant = 0L)]
  cells[, eligible := condition %in% candidates & paired_sensitive >= minimum_each_class & paired_resistant >= minimum_each_class]
  cells[, reason := fcase(cells < minimum_cells, 'fewer than 18 recipients',
    clear_sensitive < minimum_each_class, 'fewer than four clear sensitive recipients',
    clear_resistant < minimum_each_class, 'fewer than four clear resistant recipients',
    iqr_pp < minimum_iqr_pp, 'abundance IQR below 25 percentage points',
    paired_sensitive < minimum_each_class | paired_resistant < minimum_each_class,
    'too few recipients in pairs separated by 25 percentage points', default = 'eligible')]
  curves <- md[, .(doses = uniqueN(dose)), by = .(drug, depmap_id)]
  doses <- curves[, .(cells_with_two_doses = sum(doses >= 2L)), by = drug]
  transitions <- pairs[axis == 'Dose', .(transition_cells = uniqueN(recipient)), by = .(drug = unit)]
  doses <- merge(doses, transitions, by = 'drug', all.x = TRUE)
  doses[is.na(transition_cells), transition_cells := 0L]
  doses[, eligible := cells_with_two_doses >= minimum_cells & transition_cells >= minimum_each_class]
  doses[, reason := fcase(cells_with_two_doses < minimum_cells, 'fewer than 18 paired dose trajectories',
    transition_cells < minimum_each_class, 'fewer than four clear sensitive/resistant dose transitions', default = 'eligible')]
  list(cells = cells, doses = doses, cell_keys = cells[eligible == TRUE, condition],
    dose_drugs = doses[eligible == TRUE, drug])
}

# Represent the two contrast losses and small absolute anchor as ordinary glmnet
# rows. Pair differences preserve the exact shared 80-feature, zero-anchored model.
objective_blocks <- function(md, x, screen, training_energy = NULL) {
  pairs <- contrast_pairs(md, screen$cell_keys, screen$dose_drugs)
  stopifnot(all(c('Cell', 'Dose') %in% pairs$axis))
  dose_pairs <- pairs[axis == 'Dose']
  clear <- md$observed_pct <= sensitive_limit | md$observed_pct >= resistant_limit
  anchor_rows <- sort(unique(c(which(md$condition %in% screen$cell_keys & clear), dose_pairs$r, dose_pairs$s)))
  stopifnot(length(anchor_rows) > 0L)
  anchor <- data.table(row = anchor_rows, condition = md$condition[anchor_rows])
  anchor[, base_weight := 1 / .N, by = condition]
  anchor[, base_weight := base_weight / sum(base_weight)]
  xx <- rbind(x[pairs$r, , drop = FALSE] - x[pairs$s, , drop = FALSE], x[anchor$row, , drop = FALSE])
  target <- c(md$observed[pairs$r] - md$observed[pairs$s], md$observed[anchor$row])
  axis <- c(pairs$axis, rep('Anchor', nrow(anchor)))
  base_weight <- c(pairs$base_weight, anchor$base_weight)
  energy <- if (is.null(training_energy)) setNames(vapply(names(objective_mass), function(a)
    sum(base_weight[axis == a] * target[axis == a]^2), numeric(1)), names(objective_mass)) else training_energy
  stopifnot(all(is.finite(xx)), all(is.finite(target)), all(energy > 0), all(is.finite(energy)))
  weights <- objective_mass[axis] * base_weight / energy[axis]
  list(x = xx, target = target, axis = axis, base_weight = base_weight, weights = weights,
    energy = energy, pairs = pairs, anchor_rows = anchor_rows)
}

status('loading frozen scores; screening training signal')
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
stopifnot(!anyDuplicated(labels$observation_id))
all_predictions <- all_audits <- selection_tables <- source_paths <- capacity_tables <- list()
model_names <- c('Original absolute fit', 'Filtered absolute fit', 'Balanced cell + dose fit')

for (family_name in families) {
  parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family_name, '_holdout')
  paths <- file.path(parent, c('settings.qs', 'data/prepared_scores.qs', 'models/program_bank.qs', 'models/functional_decoder.qs'))
  settings <- qread(paths[1L])
  prepared <- qread(paths[2L])
  bank <- qread(paths[3L])
  old <- qread(paths[4L])
  stopifnot(identical(prepared$signature, settings$signature), identical(bank$signature, settings$signature),
    identical(old$signature, settings$signature), length(colnames(bank$weights)) == 40L,
    !any(old$training_ids %in% settings$held_ids), !any(bank$training_ids %in% settings$held_ids),
    !any(old$training_drugs %in% settings$excluded_drugs),
    identical(prepared$metadata$observation_id, rownames(prepared$measured)))
  source_paths[[family_name]] <- paths
  programs <- colnames(bank$weights)
  train <- copy(labels[depmap_id %in% old$training_ids & drug %in% old$training_drugs &
    is.finite(prism_abundance_pct) & prism_abundance_pct > 0])
  train[, `:=`(observed_pct = prism_abundance_pct, observed = log2(prism_abundance_pct / 100),
    condition = paste(drug, dose, sep = '__'))]
  stopifnot(uniqueN(train$depmap_id) == 26L, !anyDuplicated(train[, .(depmap_id, drug, dose)]))
  d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  b <- sweep(sweep(prepared$baseline[train$depmap_id, , drop = FALSE], 2L, old$baseline_mean, '-'),
    2L, old$baseline_sd, '/')
  x <- cbind(d, d * b)
  colnames(x) <- c(paste0('D_', programs), paste0('DB_', programs))
  stopifnot(all(is.finite(x)), ncol(x) == 80L, nrow(x) == nrow(train))
  screen <- screen_signal(train)
  blocks <- objective_blocks(train, x, screen)
  screen$cells[, family := family_name]
  screen$doses[, family := family_name]
  save_table(screen$cells, paste0(family_name, '_training_cell_context_screen'))
  save_table(screen$doses, paste0(family_name, '_training_dose_drug_screen'))
  save_table(data.table(family = family_name, axis = names(objective_mass),
    target_energy = blocks$energy, objective_mass = objective_mass), paste0(family_name, '_objective_scaling'))
  message(family_name, ': ', length(screen$cell_keys), ' cell contexts; ', length(screen$dose_drugs),
    ' dose-response drugs; ', sum(blocks$axis == 'Cell'), ' cell pairs; ', sum(blocks$axis == 'Dose'), ' dose pairs')

  # Whole-cell folds screen their own fitting outcomes. The frozen bank itself
  # was discovered on the full development cohort, so these are selection folds.
  ids <- sort(unique(train$depmap_id))
  set.seed(20261010L)
  inner_map <- setNames(sample(rep(1:5, length.out = length(ids))), ids)
  fold <- unname(inner_map[train$depmap_id])
  save_table(data.table(family = family_name, depmap_id = ids, fold = inner_map[ids]), paste0(family_name, '_cell_folds'))
  loss <- matrix(0, 2L, length(lambda_grid), dimnames = list(c('absolute', 'balanced'), lambda_grid))
  fold_paths <- fold_predictions <- fold_screens <- list()
  for (k in 1:5) {
    status(paste(family_name, 'fitting selection fold', k, '/ 5'))
    tr <- which(fold != k)
    va <- which(fold == k)
    stopifnot(!length(intersect(train$depmap_id[tr], train$depmap_id[va])))
    sc <- screen_signal(train[tr])
    tb <- objective_blocks(train[tr], x[tr, , drop = FALSE], sc)
    vb <- objective_blocks(train[va], x[va, , drop = FALSE], sc, tb$energy)
    ai <- which(tb$axis == 'Anchor')
    models <- list(absolute = fit_path(tb$x[ai, , drop = FALSE], tb$target[ai], tb$base_weight[ai], FALSE),
      balanced = fit_path(tb$x, tb$target, tb$weights, FALSE))
    for (name in names(models)) {
      pp <- pred_path(models[[name]], vb$x)
      if (name == 'absolute') {
        vv <- which(vb$axis == 'Anchor')
        loss[name, ] <- loss[name, ] + colSums((pp[vv, , drop = FALSE] - vb$target[vv])^2 * vb$base_weight[vv]) / 5
      } else loss[name, ] <- loss[name, ] + colSums((pp - vb$target)^2 * vb$weights) / 5
    }
    # Original objective, refitted without validation recipients at its frozen
    # original penalty, supplies a comparable development-fold reference.
    reference <- fit_path(x[tr, , drop = FALSE], train$observed[tr], observation_weights(train[tr]), FALSE)
    reference$lambda <- old$models$measured_D_conditioned$lambda
    fold_paths[[as.character(k)]] <- models
    fold_predictions[[as.character(k)]] <- list(raw = lapply(models, function(f) pred_path(f, x[va, , drop = FALSE])),
      reference = predict_ridge(reference, x[va, , drop = FALSE]), rows = va, screen = sc)
    fold_screens[[as.character(k)]] <- data.table(family = family_name, fold = k,
      fitting_cells = uniqueN(train$depmap_id[tr]), validation_cells = uniqueN(train$depmap_id[va]),
      cell_contexts = length(sc$cell_keys), dose_drugs = length(sc$dose_drugs),
      fitting_cell_pairs = sum(tb$axis == 'Cell'), validation_cell_pairs = sum(vb$axis == 'Cell'),
      fitting_dose_pairs = sum(tb$axis == 'Dose'), validation_dose_pairs = sum(vb$axis == 'Dose'))
    qsave(list(models = models, reference = reference, screen = sc, training_ids = unique(train$depmap_id[tr]),
      validation_ids = unique(train$depmap_id[va]), energy = tb$energy, source_signature = settings$signature),
      file.path(root, 'models', paste0(family_name, '_fold_', k, '.qs')))
  }
  save_table(rbindlist(fold_screens), paste0(family_name, '_fold_signal_audit'))
  selection <- rbindlist(lapply(rownames(loss), function(name) data.table(family = family_name,
    objective = name, lambda = lambda_grid, internal_objective_loss = loss[name, ])))
  selected <- apply(loss, 1L, which.min)
  selection[, selected := lambda == lambda_grid[selected[objective]]]
  save_table(selection, paste0(family_name, '_penalty_selection'))
  selection_tables[[family_name]] <- selection
  ai <- which(blocks$axis == 'Anchor')
  fits <- list(absolute = fit_path(blocks$x[ai, , drop = FALSE], blocks$target[ai], blocks$base_weight[ai], FALSE),
    balanced = fit_path(blocks$x, blocks$target, blocks$weights, FALSE))
  for (name in names(fits)) {
    fits[[name]]$lambda <- lambda_grid[selected[name]]
    stopifnot(max(abs(predict_ridge(fits[[name]], x * 0))) < 1e-12)
    qsave(list(model = fits[[name]], screen = screen, objective_mass = objective_mass, energy = blocks$energy,
      baseline_mean = old$baseline_mean, baseline_sd = old$baseline_sd, programs = programs,
      training_ids = ids, excluded_ids = settings$held_ids, excluded_drugs = settings$excluded_drugs,
      source_signature = settings$signature), file.path(root, 'models', paste0(family_name, '_', name, '_decoder.qs')))
    message(family_name, ' ', name, ': selected lambda ', fits[[name]]$lambda)
  }
  coefficients <- rbindlist(lapply(names(fits), function(name) data.table(family = family_name, objective = name,
    feature = colnames(x), coefficient = as.numeric(coef(fits[[name]]$fit, s = fits[[name]]$lambda))[-1L])))
  save_table(coefficients, paste0(family_name, '_coefficients'))

  # The already-fitted paths distinguish shrinkage from model capacity without
  # fitting extra architectures or choosing a penalty using external outcomes.
  path_predictions <- pred_path(fits$balanced, x)
  cell_rows <- which(train$condition %in% screen$cell_keys)
  capacity <- rbindlist(lapply(seq_along(lambda_grid), function(j) {
    z <- train[cell_rows, .(condition, observed)]
    z[, predicted := path_predictions[cell_rows, j]]
    z[, `:=`(yc = observed - mean(observed), pc = predicted - mean(predicted)), by = condition]
    pp <- blocks$pairs
    pp[, predicted_gap := path_predictions[r, j] - path_predictions[s, j]]
    data.table(family = family_name, lambda = lambda_grid[j],
      selected = j == selected['balanced'], training_cell_r = cor_safe(z$yc, z$pc),
      cell_pair_correct_pct = 100 * weighted.mean(as.numeric(pp[axis == 'Cell', predicted_gap] > 0),
        pp[axis == 'Cell', base_weight]),
      dose_pair_correct_pct = 100 * weighted.mean(as.numeric(pp[axis == 'Dose', predicted_gap] > 0),
        pp[axis == 'Dose', base_weight]))
  }))
  save_table(capacity, paste0(family_name, '_training_penalty_path'))
  capacity_tables[[family_name]] <- capacity

  # Produce predictions before reading outer-recipient outcomes for evaluation.
  held <- copy(as.data.table(prepared$metadata)[depmap_id %in% settings$held_ids & drug %in% old$training_drugs])
  held[, condition := paste(drug, dose, sep = '__')]
  held <- held[condition %in% screen$cell_keys | drug %in% screen$dose_drugs]
  hd <- prepared$measured[match(held$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  hb <- sweep(sweep(prepared$baseline[held$depmap_id, , drop = FALSE], 2L, old$baseline_mean, '-'),
    2L, old$baseline_sd, '/')
  hx <- cbind(hd, hd * hb)
  parts <- list()
  for (cohort in c('Fitting cells', 'Held-out cells')) {
    md <- if (cohort == 'Fitting cells') train[, .(observation_id, depmap_id, drug, dose, condition)] else held
    xx <- if (cohort == 'Fitting cells') x else hx
    for (j in seq_along(model_names)) {
      f <- if (j == 1L) old$models$measured_D_conditioned else fits[[c('absolute', 'balanced')[j - 1L]]]
      tab <- copy(md)
      tab[, `:=`(family = family_name, cohort = cohort, model = model_names[j], fold = 0L,
        cell_context = condition %in% screen$cell_keys, dose_drug = drug %in% screen$dose_drugs,
        predicted = predict_ridge(f, xx))]
      parts[[paste(cohort, j)]] <- tab[cell_context | dose_drug]
    }
  }
  for (k in 1:5) {
    saved <- fold_predictions[[as.character(k)]]
    md <- train[saved$rows, .(observation_id, depmap_id, drug, dose, condition)]
    for (j in seq_along(model_names)) {
      pp <- if (j == 1L) saved$reference else saved$raw[[c('absolute', 'balanced')[j - 1L]]][, selected[j - 1L]]
      tab <- copy(md)
      tab[, `:=`(family = family_name, cohort = 'Penalty-selection cell folds', model = model_names[j], fold = k,
        cell_context = condition %in% saved$screen$cell_keys, dose_drug = drug %in% saved$screen$dose_drugs,
        predicted = pp)]
      parts[[paste('fold', k, j)]] <- tab[cell_context | dose_drug]
    }
  }
  prediction <- rbindlist(parts, use.names = TRUE)
  save_table(prediction, paste0(family_name, '_frozen_predictions'))
  save_table(data.table(path = file.path(root, 'data', paste0(family_name, '_frozen_predictions.tsv')),
    sha256 = digest::digest(file = file.path(root, 'data', paste0(family_name, '_frozen_predictions.tsv')), algo = 'sha256')),
    paste0(family_name, '_prediction_hash'))
  all_predictions[[family_name]] <- prediction
  all_audits[[family_name]] <- data.table(family = family_name, cells = length(ids), training_rows = nrow(train),
    cell_contexts = length(screen$cell_keys), cell_drugs = uniqueN(screen$cells[eligible == TRUE, drug]),
    dose_drugs = length(screen$dose_drugs), clear_anchor_rows = length(blocks$anchor_rows),
    cell_pairs = sum(blocks$axis == 'Cell'), dose_pairs = sum(blocks$axis == 'Dose'),
    excluded_cell_overlap = length(intersect(ids, settings$held_ids)),
    excluded_drug_overlap = length(intersect(train$drug, settings$excluded_drugs)), zero_response_effect = 0)
}

status('evaluating both axes on training and excluded recipients')
prediction <- rbindlist(all_predictions)
prediction <- merge(prediction, labels[, .(observation_id, cell = recipient, observed_pct = prism_abundance_pct,
  prism_status, prism_lower_dose, prism_upper_dose)], by = 'observation_id', all.x = TRUE)
excluded <- prediction[!is.finite(observed_pct) | observed_pct <= 0]
save_table(excluded, 'missing_or_nonpositive_endpoints')
prediction <- prediction[is.finite(observed_pct) & observed_pct > 0]
prediction[, `:=`(observed = log2(observed_pct / 100), predicted_pct = 100 * 2^predicted)]
stopifnot(all(is.finite(prediction$predicted_pct)), !anyDuplicated(prediction[, .(family, cohort, model, observation_id)]))
prediction[, `:=`(observed_centered = observed - mean(observed), predicted_centered = predicted - mean(predicted)),
  by = .(family, cohort, model, fold, condition)]
save_table(prediction, 'evaluated_predictions')
save_table(rbindlist(all_audits), 'training_signal_audit')
condition_metrics <- prediction[cell_context == TRUE, .(cells = uniqueN(depmap_id),
  pearson = cor_safe(observed, predicted), sensitive = sum(observed_pct <= 50),
  resistant = sum(observed_pct > 50), clear_sensitive = sum(observed_pct <= sensitive_limit),
  clear_resistant = sum(observed_pct >= resistant_limit),
  sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)),
  by = .(family, cohort, model, fold, drug, dose, condition)]
condition_metrics[, status := fifelse(is.finite(pearson), 'Evaluable', 'Fewer than three recipients or zero variation')]
save_table(condition_metrics, 'fixed_drug_dose_metrics')

# Pairs are evaluation contrasts, not independent biological samples. Repeated
# recipients and shared endpoints do not increase the independent cell count.
pair_values <- list()
groups <- unique(prediction[, .(family, cohort, model, fold)])
for (j in seq_len(nrow(groups))) {
  g <- groups[j]
  z <- prediction[family == g$family & cohort == g$cohort & model == g$model & fold == g$fold]
  pairs <- contrast_pairs(z, unique(z[cell_context == TRUE, condition]), unique(z[dose_drug == TRUE, drug]))
  pairs[, `:=`(family = g$family, cohort = g$cohort, model = g$model, fold = g$fold,
    drug = z$drug[r], dose_resistant = z$dose[r], dose_sensitive = z$dose[s],
    cell_resistant = z$cell[r], cell_sensitive = z$cell[s],
    observed_gap_log2 = z$observed[r] - z$observed[s],
    predicted_gap_log2 = z$predicted[r] - z$predicted[s])]
  pair_values[[as.character(j)]] <- pairs
}
pairs <- rbindlist(pair_values)
save_table(pairs, 'cell_and_dose_contrast_values')
pair_metrics <- pairs[, .(contrasts = .N, contexts = uniqueN(unit),
  correctly_ordered_pct = 100 * weighted.mean(as.numeric(predicted_gap_log2 > 0), base_weight),
  observed_gap_log2 = weighted.mean(observed_gap_log2, base_weight),
  predicted_gap_log2 = weighted.mean(predicted_gap_log2, base_weight)), by = .(family, cohort, model, axis)]
pair_metrics[, `:=`(observed_resistant_sensitive_ratio = 2^observed_gap_log2,
  predicted_resistant_sensitive_ratio = 2^predicted_gap_log2)]
save_table(pair_metrics, 'contrast_order_and_magnitude_metrics')
cell_metrics <- prediction[cell_context == TRUE, .(cells = uniqueN(depmap_id), observations = .N,
  pearson = cor_safe(observed_centered, predicted_centered)), by = .(family, cohort, model)]
cell_metrics <- merge(cell_metrics, condition_metrics[, .(evaluable_contexts = sum(is.finite(pearson)),
  median_context_pearson = median(pearson[is.finite(pearson)]), sensitive = sum(sensitive),
  sensitive_called = sum(sensitive_called), resistant = sum(resistant), false_sensitive = sum(false_sensitive)),
  by = .(family, cohort, model)], by = c('family', 'cohort', 'model'))
save_table(cell_metrics, 'across_cell_summary')

status('exporting compact objective diagnostics')
colors <- setNames(c('#777777', '#D55E00', '#0072B2'), model_names)
prediction[, model := factor(model, levels = model_names)]
cell_metrics[, model := factor(model, levels = model_names)]
pair_metrics[, model := factor(model, levels = model_names)]
pair_metrics[, cohort := factor(cohort, levels = c('Fitting cells', 'Penalty-selection cell folds', 'Held-out cells'))]
plots <- list()
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height, device = grDevices::cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 150, bg = 'white')
  plots[[name]] <<- p
}
for (cohort_name in c('Fitting cells', 'Held-out cells')) {
  z <- prediction[cohort == cohort_name & cell_context == TRUE]
  m <- cell_metrics[cohort == cohort_name]
  m[, label := sprintf('Across-cell r = %.2f\n%d cells; %d informative settings\nSensitive recovered: %d/%d',
    pearson, cells, evaluable_contexts, sensitive_called, sensitive)]
  p <- ggplot(z, aes(observed_centered, predicted_centered)) +
    geom_point(aes(color = model), size = 1, alpha = .5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#222222', linewidth = .45) +
    geom_text(data = m, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.1, size = 2.9) +
    facet_wrap(vars(family, model), ncol = 3, scales = 'free') + scale_color_manual(values = colors, guide = 'none') +
    theme_classic(base_size = 10) + theme(aspect.ratio = 1) +
    labs(x = 'Measured log2 abundance: drug/dose mean removed', y = 'Predicted log2 abundance: drug/dose mean removed',
      title = paste(cohort_name, '| fixed-treatment recipient sensitivity'),
      subtitle = 'Identical observations across models; signal gates use training recipients only. Axes scale independently.',
      caption = if (cohort_name == 'Fitting cells') 'Apparent training fit. This is the first check, not evidence of transfer.' else
        'Nine original excluded recipients on known drugs. Previously inspected endpoints: retrospective evaluation.')
  save_plot(p, if (cohort_name == 'Fitting cells') '01_training_fixed_treatment_fit' else '02_held_cells_fixed_treatment_fit', 11, 8)
}
p <- ggplot(pair_metrics, aes(model, correctly_ordered_pct, color = axis, group = axis)) +
  geom_hline(yintercept = 50, linetype = 2, color = 'grey70') + geom_line(linewidth = .5) + geom_point(size = 2) +
  facet_wrap(vars(family, cohort), ncol = 3) + theme_classic(base_size = 10) +
  scale_color_manual(values = c(Cell = '#0072B2', Dose = '#D55E00')) +
  scale_y_continuous(limits = c(0, 100)) + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = 'bottom') +
  labs(x = NULL, y = 'Clear sensitivity contrasts correctly ordered (%)', color = 'Contrast axis',
    title = 'Did the objective learn both recipient differences and dose changes?',
    subtitle = 'Cell: same drug/dose, different recipients. Dose: same drug/recipient, higher versus lower dose.',
    caption = 'Contrasts share samples and are not independent replicates. Selection folds condition on the frozen development bank.')
save_plot(p, '03_cell_and_dose_contrast_ordering', 12, 8)
capacity <- rbindlist(capacity_tables)
p <- ggplot(capacity, aes(lambda, training_cell_r, color = family, group = family)) +
  geom_line(linewidth = .6) + geom_point(size = 1.8) +
  geom_point(data = capacity[selected == TRUE], shape = 8, size = 3.5) +
  scale_x_log10(breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) + theme_classic(base_size = 10) +
  scale_color_manual(values = c(ERBB = '#0072B2', MEK = '#D55E00')) + theme(legend.position = 'bottom') +
  labs(x = 'Ridge penalty (smaller = less shrinkage)', y = 'Training across-cell r at fixed drug/dose', color = 'Frozen bank',
    title = 'Can this same model fit training cell differences with weaker shrinkage?',
    subtitle = 'Stars mark penalties chosen in whole-cell development folds. All points reuse the existing fitted paths.',
    caption = 'Apparent training capacity only; the most flexible fit is not selected using held-recipient results.')
save_plot(p, '05_training_penalty_path', 7, 4.8)

# Keep every informative context visible; unknown correlations remain grey/NA.
for (family_name in families) {
  cm <- condition_metrics[family == family_name & cohort %in% c('Fitting cells', 'Held-out cells')]
  row_info <- unique(cm[, .(drug, dose, condition)])[order(drug, dose)]
  rows <- row_info$condition
  cols <- unlist(lapply(c('Fitting cells', 'Held-out cells'), function(cohort) paste(cohort, model_names, sep = ' | ')))
  mat <- matrix(NA_real_, length(rows), length(cols), dimnames = list(rows, cols))
  text_mat <- matrix('NA', length(rows), length(cols), dimnames = dimnames(mat))
  ij <- cbind(match(cm$condition, rows), match(paste(cm$cohort, cm$model, sep = ' | '), cols))
  stopifnot(!anyDuplicated(cm[, .(condition, cohort, model)]), !anyNA(ij))
  mat[ij] <- cm$pearson
  text_mat[ij] <- sprintf('%.2f\nn=%d', cm$pearson, cm$cells)
  ca <- data.frame(Cohort = rep(c('Fitting cells', 'Held-out cells'), each = 3), row.names = cols)
  stopifnot(identical(rownames(ca), colnames(mat)), identical(dimnames(mat), dimnames(text_mat)))
  h <- draw_heatmap(mat, col_anno_df = ca, col_anno_cols = 'Cohort', scale_rows = FALSE,
    value_range = c(-1, 0, 1), value_colors = c('#B2182B', 'white', '#2166AC'), hm_name = 'Pearson r',
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = TRUE,
    row_labels = paste(row_info$drug, row_info$dose, 'µM'),
    column_labels = rep(c('Original absolute', 'Filtered absolute', 'Cell + dose objective'), 2),
    column_names_rot = 35, row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 9),
    title = paste(family_name, '| every retained fixed-dose cell context'),
    cell_fun = function(j, i, x, y, width, height, fill) grid::grid.text(text_mat[i, j], x, y,
      gp = grid::gpar(fontsize = 7, col = if (is.finite(mat[i, j]) && abs(mat[i, j]) > .65) 'white' else '#222222')))
  grDevices::cairo_pdf(file.path(root, 'figures', paste0('04_', family_name, '_all_cell_contexts.pdf')),
    width = 11, height = 2.5 + nrow(mat) * .24)
  ComplexHeatmap::draw(h)
  dev.off()
  png(file.path(root, 'figures', paste0('04_', family_name, '_all_cell_contexts.png')),
    width = 11, height = 2.5 + nrow(mat) * .24, units = 'in', res = 150)
  ComplexHeatmap::draw(h)
  dev.off()
}
grDevices::cairo_pdf(file.path(root, 'objective_diagnostics.pdf'), width = 12, height = 9, onefile = TRUE)
for (p in plots) print(p)
dev.off()
sources <- unique(c(panel_path, unlist(source_paths), '/work/src/program_pilot/33_balanced_cell_and_dose_contrast_decoder.R',
  '/work/src/program_pilot/response_models.R', '/work/src/R/lib_heatmap.R'))
save_table(data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'input_hashes')
qsave(list(minimum_cells = minimum_cells, minimum_each_class = minimum_each_class,
  sensitive_limit = sensitive_limit, resistant_limit = resistant_limit, minimum_iqr_pp = minimum_iqr_pp,
  minimum_pair_gap_pp = minimum_pair_gap_pp, objective_mass = objective_mass, lambda_grid = lambda_grid,
  seed = 20261010L, families = families, representation_changed = FALSE, retrospective = TRUE), file.path(root, 'settings.qs'))
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(system('git -C /work rev-parse HEAD', intern = TRUE), file.path(root, 'git_revision.txt'))
writeLines(c('BALANCED CELL AND DOSE CONTRAST OBJECTIVE', '',
  'Start with objective_diagnostics.pdf: training fit, excluded cells, both contrast axes.',
  '04_*_all_cell_contexts.pdf shows every training-defined fixed-dose context, including failures.', '',
  'Inputs stay unchanged: 40 hard-mean measured displacements + 40 matched baseline interactions.',
  'Output stays log2(PRISM abundance / 100), with zero predicted effect at zero displacement.',
  '45% cell-contrast loss + 45% dose-contrast loss + 10% absolute anchor; ridge regularization.',
  'Each contrast axis is normalized by fitting target energy and balanced across contexts.',
  'Cell contexts: >=18 recipients, >=4 clear sensitive <=40%, >=4 clear resistant >=60%, IQR >=25 pp.',
  'Pairs must differ by >=25 pp. Dose drugs require >=18 multi-dose recipients and >=4 clear transitions.',
  'Dose contrasts compare higher-dose sensitive versus lower-dose resistant in the same recipient.',
  'No-signal fixed-dose contexts have zero cell-contrast weight but may inform qualified dose contrasts.',
  'All gates are repeated inside five whole-recipient fitting folds; penalty uses training-development validation only.',
  'Original absolute fit and filtered absolute fit distinguish screening from changing the objective.', '',
  'LIMITS',
  'Clear phenotype contrasts are signal gates, not replicate-based proof that assay noise is absent.',
  'PRISM endpoints can be bracket-interpolated. All provenance is retained in evaluated_predictions.tsv.',
  'Selection folds reuse the development-discovered bank and select penalties; they are not independent validation.',
  'The nine outer recipients were excluded from discovery and functional fitting; previously seen endpoints are retrospective.',
  'Relative ordering, contrast size and absolute sensitivity calibration are separate readouts.',
  'Measured D is used throughout. Baseline-to-D and combinations have not been tested.', '',
  'TRAINING SIGNAL', capture.output(print(rbindlist(all_audits))), '',
  'ACROSS-CELL RESULTS', capture.output(print(cell_metrics)), '',
  'CONTRAST RESULTS', capture.output(print(pair_metrics))), file.path(root, 'READING_GUIDE.txt'))
status('COMPLETE')
print(cell_metrics)
print(pair_metrics)
