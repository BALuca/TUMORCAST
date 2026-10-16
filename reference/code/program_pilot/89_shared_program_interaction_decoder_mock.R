#!/usr/bin/env Rscript
# A measured-RNA-only mock: one shared elastic-net decoder with explicit program
# products, fitted to cell-to-cell sensitivity contrasts within each drug/dose.
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(qs)
  library(ggplot2)
  library(digest)
})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

# Fixed configuration: dose choices are the user's upfront retrospective choices.
root <- '/work/results/program_pilot/89_shared_program_interaction_decoder_mock'
v84 <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
v87 <- '/work/results/program_pilot/87_selected_dose_joint_program_bridge'
script_path <- '/work/src/program_pilot/89_shared_program_interaction_decoder_mock.R'
programs <- sprintf('P%02d', 1:40)
alpha <- .2
lambda_grid <- 10^seq(log10(1), log10(.003), length.out = 30L)
seed <- 89012L
set.seed(seed)
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
pearson <- function(x, y) if (length(x) >= 4L && sd(x) > 1e-12 && sd(y) > 1e-12) cor(x, y) else NA_real_

points_path <- file.path(v84, 'data/observed_program_viability_points.tsv.gz')
selection_path <- file.path(v87, 'data/selected_doses.tsv')
fold_path <- file.path(v87, 'data/recipient_folds.tsv')
selected <- fread(selection_path)
points <- merge(fread(points_path), selected[, .(drug, dose)], by = c('drug', 'dose'))
wide <- dcast(points, depmap_id + drug + dose + observed_log2_abundance ~ program,
  value.var = 'observed_displacement')
wide <- merge(wide, fread(fold_path)[, .(depmap_id, fold)], by = 'depmap_id', all.x = TRUE)
setorder(wide, drug, depmap_id)
stopifnot(!anyNA(wide$fold), !anyDuplicated(wide[, .(depmap_id, drug, dose)]),
  nrow(wide) == sum(selected$cells), uniqueN(wide$drug) == 9L,
  !'Vemurafenib' %in% wide$drug, all(is.finite(wide$observed_log2_abundance)),
  all(is.finite(as.matrix(wide[, ..programs]))))
write_table(selected, 'selected_doses')
write_table(wide[, .(depmap_id, drug, dose, fold)], 'measured_profile_folds')
cat('Input: ', nrow(wide), ' measured profiles, ', uniqueN(wide$depmap_id), ' recipients, nine drugs.\n', sep = '')

# Build products BEFORE group centering. This retains context from the complete
# displacement pattern, including interactions between mean movements and cell
# variation. No drug indicator, baseline value, or predicted RNA enters X.
displacements <- as.matrix(wide[, ..programs])
pair_index <- t(combn(seq_along(programs), 2L))
products <- displacements[, pair_index[, 1], drop = FALSE] *
  displacements[, pair_index[, 2], drop = FALSE]
colnames(products) <- paste(programs[pair_index[, 1]], programs[pair_index[, 2]], sep = ':')
design <- list(linear = displacements, interactions = cbind(displacements, products))
feature_key <- rbind(data.table(feature = programs, type = 'main', programA = programs, programB = NA_character_),
  data.table(feature = colnames(products), type = 'interaction',
    programA = programs[pair_index[, 1]], programB = programs[pair_index[, 2]]))
stopifnot(ncol(design$interactions) == 820L, !anyDuplicated(feature_key$feature))
write_table(feature_key, 'feature_definitions')

# Reused for inner selection, outer exclusion and the final all-data mock.
# Group centering uses ONLY the recipients being fitted, after feature creation.
# Each drug/dose receives equal total weight. glmnet standardizes the centered
# columns and returns coefficients in the original physical feature units.
fit_decoder <- function(rows, x) {
  training <- wide[rows]
  groups <- unique(training$drug)
  g <- match(training$drug, groups)
  counts <- tabulate(g, nbins = length(groups))
  stopifnot(length(groups) == 9L, all(counts >= 2L))
  xt <- x[rows, , drop = FALSE]
  feature_mean <- rowsum(xt, g, reorder = FALSE) / counts
  outcome_mean <- as.numeric(rowsum(training$observed_log2_abundance, g, reorder = FALSE)) / counts
  xt <- xt - feature_mean[g, , drop = FALSE]
  yt <- training$observed_log2_abundance - outcome_mean[g]
  fit <- glmnet(xt, yt, alpha = alpha, lambda = lambda_grid, weights = 1 / counts[g],
    intercept = FALSE, standardize = TRUE, thresh = 1e-8, maxit = 100000L)
  stopifnot(fit$jerr == 0L, all(is.finite(as.matrix(coef(fit)))))
  list(fit = fit, groups = groups, feature_mean = feature_mean,
    feature_sd = sqrt(colMeans(xt^2)), outcome_mean = outcome_mean,
    training_cells = unique(training$depmap_id))
}

# Select regularization by mean SIGNED Pearson across drugs in inner held-cell
# predictions. For this sensitivity objective, fold offsets are removed only
# while scoring. Constant predictions are ineligible, so an average-only fit
# cannot be selected. No outer held-out outcome selects lambda or alpha.
select_lambda <- function(rows, x, job) {
  cells <- sort(unique(wide$depmap_id[rows]))
  # Both feature sets use identical inner recipient splits for a fair comparison.
  outer_job <- sub('_(linear|interactions)$', '', job)
  set.seed(seed + match(outer_job, c(paste0('fold_', 1:3), 'all_data')))
  inner <- setNames(sample(rep(1:3, length.out = length(cells))), cells)
  held_predictions <- matrix(NA_real_, length(rows), length(lambda_grid))
  for (j in 1:3) {
    held <- which(inner[wide$depmap_id[rows]] == j)
    training_rows <- rows[-held]
    model <- fit_decoder(training_rows, x)
    stopifnot(!any(wide$depmap_id[rows[held]] %in% model$training_cells))
    held_predictions[held, ] <- predict(model$fit, newx = x[rows[held], , drop = FALSE], s = lambda_grid)
  }
  stopifnot(all(is.finite(held_predictions)))
  group <- paste(wide$drug[rows], inner[wide$depmap_id[rows]], sep = '|')
  residual_y <- wide$observed_log2_abundance[rows]
  for (g in unique(group)) {
    ii <- which(group == g)
    residual_y[ii] <- residual_y[ii] - mean(residual_y[ii])
    held_predictions[ii, ] <- sweep(held_predictions[ii, , drop = FALSE],
      2L, colMeans(held_predictions[ii, , drop = FALSE]))
  }
  drug_correlations <- vapply(unique(wide$drug[rows]), function(drug) {
    ii <- which(wide$drug[rows] == drug)
    vapply(seq_along(lambda_grid), function(l) pearson(residual_y[ii], held_predictions[ii, l]), numeric(1))
  }, numeric(length(lambda_grid)))
  eligible <- rowSums(is.finite(drug_correlations)) == 9L
  score <- rep(-Inf, length(lambda_grid))
  score[eligible] <- rowMeans(drug_correlations[eligible, , drop = FALSE])
  stopifnot(any(eligible))
  selected_index <- which.max(score)
  write_table(data.table(job, lambda = lambda_grid, eligible, inner_mean_signed_pearson = score,
    selected = seq_along(lambda_grid) == selected_index), paste0(job, '_regularization_selection'))
  list(lambda = lambda_grid[selected_index], inner_cells = data.table(depmap_id = cells, inner_fold = inner[cells]))
}

predictions <- coefficients <- fitting_summary <- list()
jobs <- c(paste0('fold_', 1:3), 'all_data')
for (job in jobs) {
  k <- if (job == 'all_data') NA_integer_ else as.integer(sub('fold_', '', job))
  train <- if (is.na(k)) seq_len(nrow(wide)) else which(wide$fold != k)
  test <- if (is.na(k)) seq_len(nrow(wide)) else which(wide$fold == k)
  if (!is.na(k)) stopifnot(!length(intersect(wide$depmap_id[train], wide$depmap_id[test])))
  for (family in names(design)) {
    tag <- paste(job, family, sep = '_')
    chosen <- select_lambda(train, design[[family]], tag)
    model <- fit_decoder(train, design[[family]])
    beta <- as.numeric(coef(model$fit, s = chosen$lambda))[-1L]
    names(beta) <- colnames(design[[family]])
    stopifnot(sum(abs(beta) > 1e-10) > 0L)
    out <- wide[test, .(depmap_id, drug, dose, fold, observed_log2 = observed_log2_abundance)]
    out[, `:=`(decoder = family, evaluation = if (is.na(k)) 'fitted' else 'held_cell',
      predicted_function_score = as.numeric(design[[family]][test, , drop = FALSE] %*% beta))]
    predictions[[tag]] <- out
    coefficients[[tag]] <- data.table(job, decoder = family, feature = names(beta),
      coefficient = beta, effect_per_feature_SD = beta * model$feature_sd,
      active = abs(beta) > 1e-10)
    fitting_summary[[tag]] <- data.table(job, decoder = family, lambda = chosen$lambda,
      fitting_cells = length(model$training_cells), main_terms = sum(abs(beta[programs]) > 1e-10),
      interaction_terms = sum(abs(beta[setdiff(names(beta), programs)]) > 1e-10))
    qsave(list(decoder = family, coefficient = beta, lambda = chosen$lambda, alpha = alpha,
      programs = programs, pair_index = if (family == 'interactions') pair_index else NULL,
      training_cells = model$training_cells, group_feature_mean = model$feature_mean,
      group_outcome_mean = model$outcome_mean, group_names = model$groups,
      mean_terms_fit = FALSE, baseline_inputs = FALSE, drug_identity_inputs = FALSE,
      functional_score_only = TRUE, retrospective_dose_selection = TRUE),
      file.path(root, 'models', paste0(tag, '.qs')))
    write_table(chosen$inner_cells, paste0(tag, '_inner_cell_folds'))
    cat('Finished ', tag, ': lambda=', signif(chosen$lambda, 3),
      ', main terms=', sum(abs(beta[programs]) > 1e-10),
      ', products=', sum(abs(beta[setdiff(names(beta), programs)]) > 1e-10), '\n', sep = '')
  }
}
pred <- rbindlist(predictions)
write_table(pred, 'predictions_frozen')
write_table(rbindlist(fitting_summary), 'fitting_summary')
coef_table <- merge(rbindlist(coefficients), feature_key, by = 'feature', all.x = TRUE)
write_table(coef_table, 'coefficients')

# Raw-score correlations expose fold-level calibration differences. Centered
# correlations isolate sensitivity patterns, and are explicitly labelled as such.
# Neither uses fitted or held outcomes to modify saved model coefficients.
pred[, observed_variation := observed_log2 - mean(observed_log2), by = .(drug, decoder, evaluation, fold)]
pred[, predicted_variation := predicted_function_score - mean(predicted_function_score),
  by = .(drug, decoder, evaluation, fold)]
metrics <- pred[, .(cells = uniqueN(depmap_id), raw_score_r = pearson(observed_log2, predicted_function_score),
  within_fold_r = pearson(observed_variation, predicted_variation)), by = .(drug, dose, decoder, evaluation)]
write_table(metrics, 'per_drug_correlations')
comparison <- dcast(metrics[evaluation == 'held_cell'], drug + dose + cells ~ decoder, value.var = 'within_fold_r')
write_table(comparison, 'held_cell_decoder_comparison')
print(comparison)

# Export strongest all-data associations and their coefficient consistency.
# Three exclusion fits describe stability; they do not establish significance.
stability <- coef_table[job != 'all_data', .(active_folds = sum(active),
  positive_folds = sum(coefficient > 1e-10), negative_folds = sum(coefficient < -1e-10)), by = .(decoder, feature)]
signatures <- merge(coef_table[job == 'all_data' & decoder == 'interactions' & active],
  stability, by = c('decoder', 'feature'))
signatures <- signatures[order(-abs(effect_per_feature_SD))]
write_table(signatures, 'selected_program_signature_terms')

# Same square layout for linear and interaction fits, so the comparison is direct.
held <- pred[evaluation == 'held_cell']
held[, panel := paste0(drug, ' | ', dose, ' µM')]
labels <- metrics[evaluation == 'held_cell', .(decoder, panel = paste0(drug, ' | ', dose, ' µM'),
  label = sprintf('r = %.3f; n = %d', within_fold_r, cells))]
plots <- list()
for (family in names(design)) {
  z <- held[decoder == family]
  z[, fold_label := factor(fold)]
  p <- ggplot(z, aes(observed_variation, predicted_variation)) +
    geom_point(aes(color = fold_label), size = 1.3, alpha = .85) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
    geom_text(data = labels[decoder == family], aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -.05, vjust = 1.3, size = 2.8) +
    facet_wrap(~panel, ncol = 3, scales = 'free') +
    scale_color_manual(values = c('#2378a0', '#bd5a28', '#72419d'), name = 'Held-cell fold') +
    theme_bw(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom', strip.text = element_text(size = 8)) +
    labs(title = paste0('v89 measured-RNA shared decoder: ', family),
      x = 'Observed variation in PRISM log2 abundance', y = 'Predicted functional-score variation',
      caption = paste0('One selected dose per drug; each cell line excluded from coefficient fitting and regularization selection.\n',
        'Both axes have drug/dose/fold means removed to isolate sensitivity variation. Scores are not absolute viability predictions.\n',
        'Dose selection used all observed outcomes upfront; these are retrospective associations.'))
  plots[[family]] <- p
  ggsave(file.path(root, 'figures', paste0(family, '_held_cell_sensitivity.png')), p, width = 9, height = 10, dpi = 150)
  ggsave(file.path(root, 'figures', paste0(family, '_held_cell_sensitivity.pdf')), p, width = 9, height = 10, device = cairo_pdf)
}
cairo_pdf(file.path(root, 'figures/shared_decoder_comparison.pdf'), width = 9, height = 10, onefile = TRUE)
for (p in plots) print(p)
dev.off()
writeLines(c('Measured Tahoe RNA only; nine selected drug/dose contexts; no baseline or imputed displacement inputs.',
  '40 shared main effects and 780 explicit shared products; no drug-specific decoder coefficients.',
  'Features and targets are centered per drug/dose inside fitting cells only, after product construction.',
  'Elastic net alpha=0.2 fixed; lambda selected by inner whole-cell mean signed within-drug Pearson.',
  'All-zero/constant-score candidates cannot be selected. Absolute drug averages do not enter the fitting objective.',
  'Raw-score and within-fold Pearson are both exported. Figures use within-fold contrasts and do not demonstrate abundance calibration.',
  'Dose selection used all observed outcomes upfront; outer cell exclusions are coefficient-fitting exclusions, not independent confirmation.'),
  file.path(root, 'SCOPE.txt'))
paths <- c(points_path, selection_path, fold_path, script_path)
write_table(data.table(path = paths, sha256 = vapply(paths, digest, character(1), file = TRUE, algo = 'sha256')), 'input_hashes')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
cat('Complete.\n')
