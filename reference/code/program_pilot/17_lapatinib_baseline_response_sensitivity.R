#!/usr/bin/env Rscript
# Recipient-held-out lapatinib chain using training-only outcome-informed programs.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/functional_programs.R')
data.table::setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(2L)
RhpcBLASctl::omp_set_num_threads(2L)
options(warn = 1)
set.seed(20261010L)

# Fixed, pre-evaluation choices. Molecular/function exposures are explicitly a proxy.
rank <- 40L
upstream_lambda <- 0.01
functional_lambda <- 0.1
bootstrap_repeats <- 1000L
output_dir <- '/work/results/program_pilot/17_lapatinib_baseline_response_sensitivity'
data_dir <- file.path(output_dir, 'data')
model_dir <- file.path(output_dir, 'models')
figure_dir <- file.path(output_dir, 'figures')
for (path in c(output_dir, data_dir, model_dir, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
parent_dir <- '/work/results/program_pilot/15_mean_functional_displacement_correlation_programs_ward'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
labels_path <- '/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge/joint40/prism_screen_audit/preferred_screen_functional_labels.tsv'
expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
script_path <- '/work/src/program_pilot/17_lapatinib_baseline_response_sensitivity.R'
panel <- qread(panel_path)
bank_paths <- file.path(parent_dir, paste0('cell_', gsub('[^A-Za-z0-9]+', '_', panel$cell_map$depmap_name)), 'models/program_bank.qs')
check(all(file.exists(c(bank_paths, panel$input_paths, labels_path, expression_path))), 'Required source absent')
paths <- c(panel_path, labels_path, expression_path, script_path, bank_paths, panel$input_paths)
manifest <- data.table(path = paths, sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(manifest, file.path(data_dir, 'input_hashes.tsv'), sep = '\t')
signature <- digest::digest(manifest, algo = 'sha256')
status <- function(phase, fold = 0L) {
  fwrite(data.table(status = phase, completed_folds = fold, total_folds = 5L, time = as.character(Sys.time())),
    file.path(output_dir, 'run_status.tsv'), sep = '\t')
  message(format(Sys.time()), ' | ', phase, ' | completed folds: ', fold)
}
status('preparing corrected labels and baseline expression')

# Labels and existing recipient folds are read as-is; no response-based filtering.
labels <- fread(labels_path)[drug == 'Lapatinib ditosylate' & dose == 2.5,
  .(depmap_id, fold, observed = corrected_y, screen_id, dose)]
check(!anyDuplicated(labels$depmap_id) && all(is.finite(labels$observed)), 'Invalid functional cohort')
header <- names(fread(expression_path, nrows = 0))
symbols <- sub(' \\(ENSG.*', '', header[-1L])
candidate_genes <- sort(unique(c(panel$genes, 'ERBB2', 'MET')))
gene_indices <- which(symbols %chin% candidate_genes)
check(!anyDuplicated(symbols[gene_indices]), 'Ambiguous CCLE gene symbols')
expression_table <- fread(expression_path, select = c(1L, gene_indices + 1L))
expression <- as.matrix(expression_table[, -1L])
rownames(expression) <- expression_table[[1L]]
colnames(expression) <- symbols[gene_indices]
check(!anyDuplicated(rownames(expression)) && all(is.finite(expression)), 'Invalid CCLE baseline')
coverage <- labels[, .(depmap_id, fold, baseline_available = depmap_id %chin% rownames(expression),
  molecular_panel = depmap_id %chin% panel$cell_map$depmap_id)]
fwrite(coverage, file.path(data_dir, 'cohort_coverage.tsv'), sep = '\t')
check(all(coverage$baseline_available), 'A corrected-label recipient lacks CCLE baseline')
check(all(c('ERBB2', 'MET') %chin% colnames(expression)), 'Receptor genes missing')
source_coverage <- panel$cell_map[, .(depmap_id, depmap_name,
  baseline_available = depmap_id %chin% rownames(expression))]
fwrite(source_coverage, file.path(data_dir, 'discovery_baseline_coverage.tsv'), sep = '\t')
all_source_ids <- sort(unique(c(labels$depmap_id, panel$cell_map$depmap_id)))
all_ids <- intersect(all_source_ids, rownames(expression))
fold_map <- setNames(labels$fold, labels$depmap_id)
extra_ids <- setdiff(all_source_ids, names(fold_map))
if (length(extra_ids)) fold_map <- c(fold_map, setNames(rep(1:5, length.out = length(extra_ids)), extra_ids))
fwrite(data.table(depmap_id = all_ids, fold = unname(fold_map[all_ids]),
  functional_available = all_ids %chin% labels$depmap_id,
  discovery_available = all_ids %chin% panel$cell_map$depmap_id), file.path(data_dir, 'recipient_folds.tsv'), sep = '\t')
receptors <- expression[all_ids, c('ERBB2', 'MET'), drop = FALSE]

# Reconstruct each exact local relationship once; preserve the source hashes.
# Only these small admitted-gene matrices are held in memory, not whole-panel RNA.
local <- molecular_rows <- list()
for (i in seq_len(nrow(panel$cell_map))) {
  id <- panel$cell_map$depmap_id[i]
  bank <- qread(bank_paths[i])
  molecular <- qread(panel$input_paths[i])
  md <- bank$metadata
  x <- molecular$response[match(md$observation_id, rownames(molecular$response)), bank$genes, drop = FALSE]
  w <- group_weights(md$drug) * md$condition_support
  displacement_correlation <- gene_correlation(x, w)
  relationship <- functional_nmf_input(bank$functional_support, displacement_correlation, include_channels = FALSE)$correlation
  check(identical(digest::digest(relationship, algo = 'sha256'), bank$derived_matrix_signatures$clustering_correlation),
    paste('Local matrix reconstruction differs:', id))
  local[[id]] <- list(genes = bank$genes, relationship = relationship,
    observation_id = md$observation_id, weight = w / sum(w))
  mi <- which(molecular$metadata$drug == 'lapatinib' & molecular$metadata$dose == 5)
  check(length(mi) == 1L, paste('Expected one 5 uM lapatinib profile:', id))
  molecular_rows[[id]] <- molecular$response[mi, , drop = FALSE]
  message('Source relationship ready: ', id, ' (', i, '/', nrow(panel$cell_map), ')')
}
molecular_gene <- do.call(rbind, molecular_rows)
rownames(molecular_gene) <- names(molecular_rows)
rm(bank, molecular, x, relationship, displacement_correlation, molecular_rows, expression_table)
invisible(gc())

# Common, fixed-penalty ridge operators. No held-out label selects a parameter.
fit_upstream <- function(ids, scores) {
  glmnet(receptors[ids, , drop = FALSE], scores[ids, , drop = FALSE], family = 'mgaussian',
    alpha = 0, lambda = upstream_lambda, standardize = TRUE, standardize.response = TRUE,
    intercept = TRUE, thresh = 1e-8)
}
predict_upstream <- function(fit, ids) {
  z <- predict(fit, newx = receptors[ids, , drop = FALSE], s = upstream_lambda)[, , 1L, drop = FALSE]
  matrix(z, nrow = length(ids), ncol = rank, dimnames = list(ids, sprintf('P%02d', 1:rank)))
}
ridge_function <- function(x, y, intercept = FALSE) {
  rms <- sqrt(colMeans(x^2))
  keep <- is.finite(rms) & rms > 1e-10
  check(sum(keep) >= 2L, 'Insufficient varying functional predictors')
  z <- sweep(x[, keep, drop = FALSE], 2L, rms[keep], '/')
  fit <- glmnet(z, y, alpha = 0, lambda = functional_lambda, standardize = FALSE,
    intercept = intercept, thresh = 1e-8)
  list(fit = fit, scale = rms, keep = keep)
}
predict_function <- function(fit, x) {
  as.numeric(predict(fit$fit, newx = sweep(x[, fit$keep, drop = FALSE], 2L,
    fit$scale[fit$keep], '/'), s = functional_lambda))
}
functional_design <- function(displacement, baseline, conditioned) {
  answer <- cbind(displacement, displacement^2)
  if (conditioned) answer <- cbind(answer, displacement * baseline)
  answer
}
safe_cor <- function(x, y) {
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  cor(x, y)
}
save_scatter <- function(d, filename, title, xlab, ylab) {
  r <- safe_cor(d$observed, d$predicted)
  p <- ggplot(d, aes(observed, predicted)) + geom_point(size = 1.5, alpha = .7) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .55) +
    labs(title = paste(strwrap(title, width = 43), collapse = '\n'),
      subtitle = sprintf('n = %d; Pearson r = %.3f', nrow(d), r), x = xlab, y = ylab) +
    theme_classic(base_size = 10) + theme(aspect.ratio = 1, plot.title = element_text(size = 10))
  ggsave(file.path(figure_dir, paste0(filename, '.pdf')), p, width = 4.3, height = 4.3)
  ggsave(file.path(figure_dir, paste0(filename, '.png')), p, width = 4.3, height = 4.3, dpi = 160)
}

all_function <- all_molecular <- audits <- list()
for (fold in 1:5) {
  status(paste('building training-only rank-40 programs; fold', fold), fold - 1L)
  held_ids <- names(fold_map)[fold_map == fold]
  discovery_ids <- setdiff(names(local), held_ids)
  genes <- sort(unique(unlist(lapply(local[discovery_ids], `[[`, 'genes'))))
  check(all(genes %chin% colnames(expression)), 'A selected gene lacks baseline expression; no silent dropping')
  total <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  count <- matrix(0L, length(genes), length(genes))
  for (id in discovery_ids) {
    ix <- match(local[[id]]$genes, genes)
    total[ix, ix] <- total[ix, ix] + local[[id]]$relationship
    count[ix, ix] <- count[ix, ix] + 1L
  }
  correlation <- total / pmax(count, 1L)
  check(all(is.finite(correlation)) && max(abs(correlation - t(correlation))) < 1e-12, 'Invalid training aggregate')
  distance <- as.dist(sqrt(2 * (1 - pmin(pmax(correlation, -1), 1))))
  tree <- hclust(distance, method = 'ward.D2')
  partition <- ordered_partition(tree, rank)
  programs <- sprintf('P%02d', 1:rank)
  no_edge_fraction <- mean(count[upper.tri(count)] == 0L)
  model_path <- file.path(model_dir, sprintf('fold_%02d.qs', fold))
  representation <- list(fold = fold, genes = genes, partition = partition, tree = tree,
    discovery_ids = discovery_ids, excluded_ids = held_ids, signature = signature,
    no_edge_fraction = no_edge_fraction, rank = rank)
  qsave(representation, model_path)
  fwrite(data.table(gene = genes, program = programs[partition]),
    file.path(data_dir, sprintf('fold_%02d_gene_assignments.tsv', fold)), sep = '\t')
  rm(total, count, correlation, distance)
  invisible(gc())

  # Fit the existing positive squared correlation weights once on training cells.
  # Streaming sufficient statistics avoid concatenating every cell's full matrix.
  seed_weights <- sparseMatrix(i = seq_along(genes), j = partition,
    x = 1 / tabulate(partition, nbins = rank)[partition], dims = c(length(genes), rank))
  sx <- sx2 <- numeric(length(genes))
  ss <- ss2 <- numeric(rank)
  cross <- matrix(0, length(genes), rank)
  mass <- 0
  for (id in discovery_ids) {
    molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
    x <- molecular$response[match(local[[id]]$observation_id, rownames(molecular$response)), genes, drop = FALSE]
    w <- local[[id]]$weight / length(discovery_ids)
    seeds <- as.matrix(x %*% seed_weights)
    sx <- sx + colSums(x * w)
    sx2 <- sx2 + colSums(x^2 * w)
    ss <- ss + colSums(seeds * w)
    ss2 <- ss2 + colSums(seeds^2 * w)
    cross <- cross + crossprod(x, seeds * w)
    mass <- mass + sum(w)
  }
  check(abs(mass - 1) < 1e-10, 'Training projection weights do not sum to one')
  covariance <- cross - outer(sx, ss)
  scales <- sqrt(outer(sx2 - sx^2, ss2 - ss^2))
  weights <- pmax(covariance / scales, 0)^2
  weights <- sweep(weights, 2L, colSums(weights), '/')
  dimnames(weights) <- list(genes, programs)
  check(all(is.finite(weights)) && max(abs(colSums(weights) - 1)) < 1e-10, 'Invalid frozen expression projection')
  measured <- molecular_gene[, genes, drop = FALSE] %*% weights
  baseline <- expression[all_ids, genes, drop = FALSE] %*% weights
  representation$weights <- weights
  qsave(representation, model_path)
  fwrite(cbind(data.table(gene = genes), as.data.table(weights)),
    file.path(data_dir, sprintf('fold_%02d_gene_weights.tsv.gz', fold)), sep = '\t')
  # Representation inspection precedes functional fitting. Every program is shown.
  heat <- draw_heatmap(t(measured[discovery_ids, , drop = FALSE]), scale_rows = FALSE,
    cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE,
    hm_name = 'Displacement', title = sprintf('Fold %d training cells | 5 uM lapatinib | frozen scores', fold))
  grDevices::pdf(file.path(figure_dir, sprintf('fold_%02d_training_displacements.pdf', fold)), width = 8, height = 7)
  ComplexHeatmap::draw(heat)
  grDevices::dev.off()
  fwrite(data.table(program = programs, as.data.table(cor(measured[discovery_ids, , drop = FALSE]))),
    file.path(data_dir, sprintf('fold_%02d_program_correlations.tsv', fold)), sep = '\t')
  rm(molecular, x, seeds, sx, sx2, ss, ss2, covariance, scales, cross)
  invisible(gc())
  status(paste('fitting molecular and functional stages; fold', fold), fold - 1L)

  # Molecular cross-validation: test cells have no influence on the representation.
  molecular_training_ids <- intersect(discovery_ids, all_ids)
  upstream <- fit_upstream(molecular_training_ids, measured)
  predicted <- predict_upstream(upstream, all_ids)
  held_molecular <- intersect(intersect(held_ids, rownames(measured)), all_ids)
  response_mean <- colMeans(measured[molecular_training_ids, , drop = FALSE])
  response_sd <- apply(measured[molecular_training_ids, , drop = FALSE], 2L, sd)
  check(all(response_sd > 0), 'Constant training program response')
  mp <- rbindlist(lapply(held_molecular, function(id) data.table(fold = fold, depmap_id = id,
    program = programs, observed = measured[id, ], predicted = predicted[id, ],
    training_mean = response_mean, training_sd = response_sd)))
  all_molecular[[fold]] <- mp
  fwrite(mp, file.path(data_dir, sprintf('fold_%02d_molecular_predictions.tsv', fold)), sep = '\t')
  # Upstream out-of-fold inputs for paired training recipients; all doses stay grouped.
  crossfit <- predicted
  for (inner in setdiff(1:5, fold)) {
    target <- molecular_training_ids[fold_map[molecular_training_ids] == inner]
    if (!length(target)) next
    inner_train <- setdiff(molecular_training_ids, target)
    check(length(inner_train) >= 5L, 'Insufficient upstream training recipients')
    inner_fit <- fit_upstream(inner_train, measured)
    crossfit[target, ] <- predict_upstream(inner_fit, target)
  }
  training <- labels[!depmap_id %chin% held_ids]
  testing <- labels[depmap_id %chin% held_ids]
  check(!length(intersect(training$depmap_id, held_ids)), 'Functional recipient leakage')
  baseline_mean <- colMeans(baseline[training$depmap_id, , drop = FALSE])
  baseline_sd <- apply(baseline[training$depmap_id, , drop = FALSE], 2L, sd)
  check(all(baseline_sd > 0), 'Constant baseline coordinate')
  b <- sweep(sweep(baseline, 2L, baseline_mean, '-'), 2L, baseline_sd, '/')
  dt <- crossfit[training$depmap_id, , drop = FALSE]
  de <- predicted[testing$depmap_id, , drop = FALSE]
  bt <- b[training$depmap_id, , drop = FALSE]
  be <- b[testing$depmap_id, , drop = FALSE]
  fits <- predictions <- list()
  predictions$training_mean <- rep(mean(training$observed), nrow(testing))
  for (conditioned in c(FALSE, TRUE)) {
    name <- if (conditioned) 'baseline_conditioned_chain' else 'no_functional_baseline'
    fits[[name]] <- ridge_function(functional_design(dt, bt, conditioned), training$observed)
    predictions[[name]] <- predict_function(fits[[name]], functional_design(de, be, conditioned))
    check(max(abs(predict_function(fits[[name]], functional_design(de * 0, be, conditioned)))) < 1e-10,
      'Functional zero-response anchor failed')
  }
  direct_x <- cbind(b, receptors)
  fits$direct_baseline <- ridge_function(direct_x[training$depmap_id, , drop = FALSE], training$observed, intercept = TRUE)
  predictions$direct_baseline <- predict_function(fits$direct_baseline, direct_x[testing$depmap_id, , drop = FALSE])
  fp <- rbindlist(lapply(names(predictions), function(name) data.table(fold = fold,
    depmap_id = testing$depmap_id, model = name, observed = testing$observed, predicted = predictions[[name]],
    measured_response_available = testing$depmap_id %chin% rownames(measured))))
  paired_test <- testing[depmap_id %chin% rownames(measured)]
  if (nrow(paired_test)) {
    for (conditioned in c(FALSE, TRUE)) {
      name <- if (conditioned) 'baseline_conditioned_chain' else 'no_functional_baseline'
      pred <- predict_function(fits[[name]], functional_design(measured[paired_test$depmap_id, , drop = FALSE],
        b[paired_test$depmap_id, , drop = FALSE], conditioned))
      fp <- rbind(fp, data.table(fold = fold, depmap_id = paired_test$depmap_id,
        model = paste0(name, '_measured_substitution'), observed = paired_test$observed,
        predicted = pred, measured_response_available = TRUE))
    }
  }
  check(all(is.finite(fp$predicted)), 'Nonfinite functional prediction')
  all_function[[fold]] <- fp
  representation$upstream <- upstream
  representation$functional_models <- fits
  representation$baseline_center <- baseline_mean
  representation$baseline_scale <- baseline_sd
  representation$predicted_displacements <- predicted
  representation$crossfit_displacements <- crossfit
  representation$measured_displacements <- measured
  representation$functional_training_ids <- training$depmap_id
  qsave(representation, model_path)
  fwrite(fp, file.path(data_dir, sprintf('fold_%02d_functional_predictions.tsv', fold)), sep = '\t')
  fwrite(rbindlist(all_function), file.path(data_dir, 'functional_predictions.tsv'), sep = '\t')
  fwrite(rbindlist(all_molecular), file.path(data_dir, 'molecular_predictions.tsv'), sep = '\t')
  audits[[fold]] <- data.table(fold = fold, discovery_cells = length(discovery_ids),
    excluded_discovery_cells = length(held_molecular), training_genes = length(genes), programs = rank,
    unobserved_gene_pair_fraction = no_edge_fraction, functional_train = nrow(training), functional_test = nrow(testing),
    measured_test = nrow(paired_test), held_cell_in_discovery = FALSE, zero_response_anchor = TRUE)
  fwrite(rbindlist(audits), file.path(data_dir, 'fold_audit.tsv'), sep = '\t')
  save_scatter(fp[model == 'baseline_conditioned_chain'], sprintf('fold_%02d_sensitivity', fold),
    paste('Lapatinib sensitivity | held-out fold', fold), 'Measured log2 abundance', 'Predicted log2 abundance')
  pm <- copy(mp)
  pm[, `:=`(observed = (observed - training_mean) / training_sd,
    predicted = (predicted - training_mean) / training_sd)]
  save_scatter(pm, sprintf('fold_%02d_molecular_response', fold), paste('Molecular response | held-out fold', fold),
    'Measured displacement (training SD)', 'Predicted displacement (training SD)')
  status('fold saved with plots', fold)
  rm(representation, weights, measured, baseline, predicted, crossfit, fits, upstream)
  invisible(gc())
}

# All recipients, paired comparisons and uncertainty; no selection of a favorable program.
fp <- rbindlist(all_function)
mp <- rbindlist(all_molecular)
metrics <- fp[, .(n = .N, pearson = safe_cor(observed, predicted),
  rmse = sqrt(mean((observed - predicted)^2)), prediction_sd = sd(predicted)), by = model]
fwrite(metrics, file.path(data_dir, 'functional_metrics.tsv'), sep = '\t')
fwrite(fp[, .(n = .N, pearson = safe_cor(observed, predicted), rmse = sqrt(mean((observed - predicted)^2))),
  by = .(fold, model)], file.path(data_dir, 'functional_metrics_by_fold.tsv'), sep = '\t')
paired <- fp[measured_response_available == TRUE]
fwrite(paired[, .(n = .N, pearson = safe_cor(observed, predicted), rmse = sqrt(mean((observed - predicted)^2))),
  by = model], file.path(data_dir, 'paired_stage_metrics.tsv'), sep = '\t')
mm <- mp[, .(n = .N, pearson = safe_cor(observed, predicted),
  model_mse = mean((observed - predicted)^2), mean_response_mse = mean((observed - training_mean)^2)),
  by = .(fold, program)]
fwrite(mm, file.path(data_dir, 'molecular_metrics_by_fold_program.tsv'), sep = '\t')
ms <- mp[, .(recipients = uniqueN(depmap_id), standardized_model_mse = mean(((observed - predicted) / training_sd)^2),
  standardized_mean_response_mse = mean(((observed - training_mean) / training_sd)^2)), by = fold]
fwrite(ms, file.path(data_dir, 'molecular_metrics_by_fold.tsv'), sep = '\t')
association <- merge(mp, labels[, .(depmap_id, abundance = observed)], by = 'depmap_id')
association <- association[, .(n = .N, pearson = safe_cor(observed, abundance)), by = .(fold, program)]
fwrite(association, file.path(data_dir, 'measured_response_function_associations.tsv'), sep = '\t')
association_plot <- ggplot(association, aes(program, pearson)) + geom_point(size = 1) +
  geom_hline(yintercept = 0, color = 'grey70', linewidth = .3) + coord_flip() + facet_wrap(~fold, nrow = 1) +
  labs(x = 'Fold-local program', y = 'Pearson r: measured displacement vs abundance',
    title = 'All measured response–function associations',
    subtitle = 'Excluded recipients only; different programs in each fold; small per-fold sample sizes') +
  theme_classic(base_size = 9)
ggsave(file.path(figure_dir, 'measured_response_function_associations.pdf'), association_plot, width = 9, height = 7)
wide <- dcast(fp[!grepl('substitution', model)], depmap_id + fold + observed ~ model, value.var = 'predicted')
set.seed(20261010L)
boot <- rbindlist(lapply(seq_len(bootstrap_repeats), function(i) {
  z <- wide[sample.int(nrow(wide), replace = TRUE)]
  data.table(draw = i, baseline_mse_gain = mean((z$observed - z$no_functional_baseline)^2 -
    (z$observed - z$baseline_conditioned_chain)^2),
    direct_baseline_mse_gain = mean((z$observed - z$direct_baseline)^2 - (z$observed - z$baseline_conditioned_chain)^2))
}))
fwrite(boot, file.path(data_dir, 'paired_recipient_bootstrap.tsv'), sep = '\t')
for (name in unique(fp$model)) {
  save_scatter(fp[model == name], paste0('sensitivity_', name), gsub('_', ' ', name),
    'Measured log2 abundance | 2.5 uM', 'Predicted log2 abundance | 5 uM RNA proxy')
}
p <- ggplot(metrics, aes(reorder(model, pearson), pearson)) + geom_point(size = 2) + coord_flip() +
  labs(x = NULL, y = 'Held-out Pearson correlation', title = 'Lapatinib | new training-fold programs') +
  theme_classic(base_size = 10)
ggsave(file.path(figure_dir, 'functional_model_comparison.pdf'), p, width = 7, height = 3.2)
capture.output(sessionInfo(), file = file.path(data_dir, 'sessionInfo.txt'))
check(all(fp[!grepl('substitution', model), .N, by = .(depmap_id, model)]$N == 1L), 'Repeated evaluation recipient')
check(all(fp[!grepl('substitution', model), uniqueN(depmap_id), by = model]$V1 == nrow(labels)), 'Incomplete evaluation')
writeLines(c('Lapatinib baseline -> response -> sensitivity: retrospective held-recipient evaluation.',
  'Five training-only Ward.D2 program sets at rank 40; program IDs are fold-local.',
  'Each fold excludes held recipients from gene selection, relationship aggregation, weights and model fitting.',
  '5 uM molecular profiles predict corrected provider-preferred PRISM outcomes at 2.5 uM: dose proxy.',
  'Ridge lambda .01 upstream; .1 downstream. Baseline interactions are program-matched B*D.',
  'Measured substitution uses the SAME held-out decoder; it is a stage diagnostic, not a separately optimized model.',
  'Bootstrap intervals condition on fitted predictions; previously explored outcomes remain retrospective.'),
  file.path(output_dir, 'README.txt'))
status('complete', 5L)
print(metrics)
