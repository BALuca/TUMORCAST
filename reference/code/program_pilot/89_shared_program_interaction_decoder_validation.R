#!/usr/bin/env Rscript
# Replay saved predictions and show the single all-data fit separately from the
# held-cell figures. This does not refit or select any coefficients.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
})
root <- '/work/results/program_pilot/89_shared_program_interaction_decoder_mock'
source_root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
programs <- sprintf('P%02d', 1:40)
selected <- fread(file.path(root, 'data/selected_doses.tsv'))
points <- merge(fread(file.path(source_root, 'data/observed_program_viability_points.tsv.gz')),
  selected[, .(drug, dose)], by = c('drug', 'dose'))
wide <- dcast(points, depmap_id + drug + dose ~ program, value.var = 'observed_displacement')
pred <- fread(file.path(root, 'data/predictions_frozen.tsv'))
stopifnot(!anyDuplicated(pred[, .(depmap_id, drug, dose, decoder, evaluation)]),
  all(is.finite(pred$predicted_function_score)), uniqueN(wide$depmap_id) == 30L)
replay <- list()
for (job in c(paste0('fold_', 1:3), 'all_data')) {
  for (family in c('linear', 'interactions')) {
    tag <- paste(job, family, sep = '_')
    model <- qread(file.path(root, 'models', paste0(tag, '.qs')))
    if (job == 'all_data') {
      rows <- pred[evaluation == 'fitted' & decoder == family]
    } else {
      k <- as.integer(sub('fold_', '', job))
      rows <- pred[evaluation == 'held_cell' & decoder == family & fold == k]
      stopifnot(!any(rows$depmap_id %in% model$training_cells))
    }
    matched <- merge(rows, wide, by = c('depmap_id', 'drug', 'dose'), all.x = TRUE, sort = FALSE)
    x <- as.matrix(matched[, ..programs])
    if (family == 'interactions') {
      pairs <- model$pair_index
      products <- x[, pairs[, 1], drop = FALSE] * x[, pairs[, 2], drop = FALSE]
      colnames(products) <- paste(programs[pairs[, 1]], programs[pairs[, 2]], sep = ':')
      x <- cbind(x, products)
    }
    stopifnot(identical(colnames(x), names(model$coefficient)),
      !model$baseline_inputs, !model$drug_identity_inputs, !model$mean_terms_fit)
    difference <- max(abs(as.numeric(x %*% model$coefficient) - matched$predicted_function_score))
    stopifnot(difference < 1e-10)
    a <- fread(file.path(root, 'data', paste0(job, '_linear_inner_cell_folds.tsv')))
    b <- fread(file.path(root, 'data', paste0(job, '_interactions_inner_cell_folds.tsv')))
    stopifnot(identical(a, b), !anyNA(a$inner_fold))
    if (job != 'all_data') stopifnot(!any(rows$depmap_id %in% a$depmap_id))
    replay[[tag]] <- data.table(job, decoder = family, rows = nrow(rows),
      maximum_replay_difference = difference, identical_inner_recipient_splits = TRUE)
  }
}
fwrite(rbindlist(replay), file.path(root, 'data/prediction_replay_validation.tsv'), sep = '\t')

# Show the all-data fit with ordinary within-drug Pearson, explicitly marked as
# fitted observations. Global treatment means are removed for readable axes.
fitted <- pred[decoder == 'interactions' & evaluation == 'fitted']
fitted[, `:=`(observed_variation = observed_log2 - mean(observed_log2),
  predicted_variation = predicted_function_score - mean(predicted_function_score)), by = drug]
fitted[, panel := paste0(drug, ' | ', dose, ' µM')]
labels <- fitted[, .(label = sprintf('Fitted r = %.3f; n = %d',
  cor(observed_log2, predicted_function_score), .N)), by = panel]
p <- ggplot(fitted, aes(observed_variation, predicted_variation)) +
  geom_point(size = 1.5, color = '#2378a0', alpha = .85) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.3, size = 2.8) +
  facet_wrap(~panel, ncol = 3, scales = 'free') +
  theme_bw(base_size = 10) + theme(aspect.ratio = 1, strip.text = element_text(size = 8)) +
  labs(title = 'v89: one shared interaction decoder fitted to all measured observations',
    x = 'Observed variation in PRISM log2 abundance', y = 'Fitted functional-score variation',
    caption = paste0('Fitted observations: all 30 recipients contributed to these coefficients. This figure does not evaluate held-out predictions.\n',
      'One common function across all nine drugs; five selected program products; each treatment mean removed for plotting.\n',
      'Scores describe sensitivity variation and are not calibrated absolute viability predictions.'))
ggsave(file.path(root, 'figures/interactions_all_data_fit.png'), p, width = 9, height = 9.5, dpi = 150)
ggsave(file.path(root, 'figures/interactions_all_data_fit.pdf'), p, width = 9, height = 9.5, device = cairo_pdf)
cat('Verified: saved prediction replay, recipient exclusions, identical inner folds, shared molecular-only inputs.\n')
