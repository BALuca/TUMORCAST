#!/usr/bin/env Rscript
# Shared multidrug function mapping, focal baseline conditioning, and pair transfer.
# All outer evaluation cells are excluded from outcome-informed program discovery.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(glmnet)
  library(ggplot2)
})
source('/work/src/program_pilot/config.R')
source('/work/src/program_pilot/functional_programs.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(2L)
RhpcBLASctl::omp_set_num_threads(2L)
options(warn = 1)
set.seed(20261010L)

# Fixed choices approved before this run; no rank or outcome-based drug selection.
program_bank_type <- 'step15_union'
rank <- 40L # Used only by the previous combined Ward bank; union size is learned below.
lambda_grid <- c(100, 10, 1, .1, .01)
held_ids <- c('ACH-000138', 'ACH-000178', 'ACH-000493', 'ACH-000552', 'ACH-000669',
  'ACH-000875', 'ACH-000927', 'ACH-000750', 'ACH-000463')
focal_drugs <- c('lapatinib', 'crizotinib', 'alpelisib', 'trametinib', 'gefitinib', 'everolimus', 'bortezomib')
target_genes <- list(lapatinib = c('ERBB2', 'EGFR'), crizotinib = c('MET', 'ALK', 'ROS1'),
  alpelisib = 'PIK3CA', trametinib = c('MAP2K1', 'MAP2K2'), gefitinib = 'EGFR',
  everolimus = 'MTOR', bortezomib = 'PSMB5')
# A separate evaluation mode reuses the exact discovery and model helpers.
# It exits before touching the historical union-bank outputs.
family_argument <- commandArgs(TRUE)
if (length(family_argument) == 2L && identical(family_argument[2L], '--contextual-audit') &&
    grepl('^--family-transfer=(MEK|ERBB)$', family_argument[1L])) {
  family <- sub('--family-transfer=', '', family_argument[1L], fixed = TRUE)
  source('/work/src/program_pilot/response_models.R')
  source('/work/src/program_pilot/nonlinear_response.R')
  audit_contextual_failures(family)
  quit(save = 'no', status = 0L)
}
if (length(family_argument) == 2L && identical(family_argument[2L], '--contextual') &&
    grepl('^--family-transfer=(MEK|ERBB)$', family_argument[1L])) {
  family <- sub('--family-transfer=', '', family_argument[1L], fixed = TRUE)
  source('/work/src/program_pilot/response_models.R')
  source('/work/src/program_pilot/nonlinear_response.R')
  run_contextual_response(family)
  quit(save = 'no', status = 0L)
}
if (length(family_argument) == 2L && identical(family_argument[2L], '--nonlinear') &&
    grepl('^--family-transfer=(MEK|ERBB)$', family_argument[1L])) {
  family <- sub('--family-transfer=', '', family_argument[1L], fixed = TRUE)
  source('/work/src/program_pilot/response_models.R')
  source('/work/src/program_pilot/nonlinear_response.R')
  run_nonlinear_response(family)
  quit(save = 'no', status = 0L)
}
if (length(family_argument) == 1L && grepl('^--family-transfer=(MEK|ERBB)$', family_argument)) {
  family <- sub('--family-transfer=', '', family_argument, fixed = TRUE)
  source('/work/src/program_pilot/response_models.R')
  source('/work/src/program_pilot/family_transfer.R')
  quit(save = 'no', status = 0L)
}
script_path <- '/work/src/program_pilot/20_multidrug_baseline_response_function.R'
out <- '/work/results/program_pilot/20_multidrug_baseline_response_function_union_programs'
data_dir <- file.path(out, 'data')
model_dir <- file.path(out, 'models')
figure_dir <- file.path(out, 'figures')
for (p in c(out, data_dir, model_dir, figure_dir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(data_dir, name), sep = '\t')
status <- function(stage) {
  message(format(Sys.time()), ' | ', stage)
  fwrite(data.table(stage = stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(out, 'run_status.tsv'), sep = '\t')
}
options(error = function() {
  status('FAILED: see run.log')
  traceback(10)
  quit(save = 'no', status = 1)
})
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
parent_dir <- '/work/results/program_pilot/15_mean_functional_displacement_correlation_programs_ward'
expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
combo_path <- '/work/almanac_broad9_2026-10-07/results/prepared.rds'
control_path <- '/work/almanac_broad9_2026-10-07/structural_diagnosis/results/raw_eligible_records.tsv.gz'
panel <- qread(panel_path)
cell_map <- as.data.table(panel$cell_map)
stopifnot(all(held_ids %chin% cell_map$depmap_id), !'ACH-000248' %chin% cell_map$depmap_id)
bank_paths <- file.path(parent_dir, paste0('cell_', gsub('[^A-Za-z0-9]+', '_', cell_map$depmap_name)),
  'models/program_bank.qs')
discovery_ids <- setdiff(cell_map$depmap_id, held_ids)
paths <- c(panel_path, expression_path, combo_path, control_path, bank_paths, panel$input_paths,
  '/work/src/program_pilot/functional_programs.R', '/work/src/program_pilot/utilities.R')
stopifnot(all(file.exists(paths)))
manifest <- data.table(path = paths, sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
write_table(manifest, 'input_hashes.tsv')
prep_signature <- digest::digest(list(manifest, held_ids, program_bank_type, rank, target_genes,
  scoring = 'hard_mean'), algo = 'sha256')
write_table(data.table(path = script_path, sha256 = digest::digest(file = script_path, algo = 'sha256')), 'code_hash.tsv')
analysis_sources <- c(script_path, '/work/almanac_broad9_2026-10-07/model.R', '/work/src/R/lib_heatmap.R')
write_table(data.table(path = analysis_sources,
  sha256 = vapply(analysis_sources, digest::digest, character(1), file = TRUE, algo = 'sha256')), 'analysis_code_hashes.tsv')
qsave(list(held_ids = held_ids, discovery_ids = discovery_ids, focal_drugs = focal_drugs,
  targets = target_genes, program_bank_type = program_bank_type, rank = rank,
  lambda_grid = lambda_grid, signature = prep_signature), file.path(out, 'settings.qs'))
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))

# Prepare the one frozen training-only bank and all native-dose hard scores.
# Reuse only a data/settings-hashed cache; model fitting below always runs afresh.
cache_path <- file.path(data_dir, 'prepared_scores.qs')
if (file.exists(cache_path)) {
  status('loading verified hard-score cache')
  prepared <- qread(cache_path)
  stopifnot(identical(prepared$signature, prep_signature))
  bank <- prepared$bank
  measured <- prepared$measured
  md <- prepared$metadata
  baseline <- prepared$baseline
  targets <- prepared$targets
} else {
  if (program_bank_type == 'step15_union') {
    status('assembling the union of training-cell script-15 programs')
    membership_rows <- source_rows <- list()
    for (id in discovery_ids) {
      i <- match(id, cell_map$depmap_id)
      local_bank <- qread(bank_paths[i])
      local_programs <- sort(unique(local_bank$partition))
      source_rows[[id]] <- data.table(depmap_id = id, cell_name = cell_map$depmap_name[i],
        programs = length(local_programs), genes = length(local_bank$genes),
        status = if (length(local_programs)) 'included' else 'no resolved script-15 programs')
      if (!length(local_programs)) next
      stopifnot(length(local_bank$partition) == length(local_bank$genes), !anyNA(local_bank$partition))
      prefix <- gsub('[^A-Za-z0-9]+', '_', cell_map$depmap_name[i])
      membership_rows[[id]] <- data.table(gene = local_bank$genes,
        program = paste0(prefix, '__P', sprintf('%02d', local_bank$partition)),
        source_depmap_id = id, source_cell = cell_map$depmap_name[i], source_program = local_bank$partition)
    }
    assignments <- rbindlist(membership_rows)
    source_audit <- rbindlist(source_rows)
    stopifnot(!any(assignments$source_depmap_id %chin% held_ids),
      !anyDuplicated(assignments[, .(gene, program)]))
    genes <- sort(unique(assignments$gene))
    programs <- unique(assignments$program)
    rank <- length(programs)
    sizes <- assignments[, .(genes = .N), by = .(program, source_depmap_id, source_cell, source_program)]
    sizes <- sizes[match(programs, program)]
    weights <- sparseMatrix(i = match(assignments$gene, genes), j = match(assignments$program, programs),
      x = 1 / sizes$genes[match(assignments$program, programs)], dims = c(length(genes), rank),
      dimnames = list(genes, programs))
    # Source programs may overlap across cells; never merge or recluster them.
    stopifnot(all(rowSums(weights > 0) >= 1L), max(abs(colSums(weights) - 1)) < 1e-12)
    bank <- list(genes = genes, weights = weights, program_membership = assignments,
      program_sources = sizes, discovery_ids = discovery_ids,
      contributing_ids = source_audit[programs > 0, depmap_id], excluded_ids = held_ids,
      signature = prep_signature, type = program_bank_type)
    qsave(bank, file.path(model_dir, 'program_bank.qs'))
    write_table(assignments, 'gene_assignments.tsv')
    write_table(sizes, 'program_sizes.tsv')
    write_table(source_audit, 'program_source_audit.tsv')
    message('Union bank: ', rank, ' programs from ', length(bank$contributing_ids),
      ' training cells; ', length(genes), ' unique genes; no reclustering')
  } else {
    stopifnot(program_bank_type == 'combined_ward')
  status('reconstructing training-cell gene relationships')
  local <- list()
  for (id in discovery_ids) {
    i <- match(id, cell_map$depmap_id)
    source_bank <- qread(bank_paths[i])
    molecular <- qread(panel$input_paths[i])
    x <- molecular$response[match(source_bank$metadata$observation_id, rownames(molecular$response)),
      source_bank$genes, drop = FALSE]
    w <- group_weights(source_bank$metadata$drug) * source_bank$metadata$condition_support
    relation <- functional_nmf_input(source_bank$functional_support, gene_correlation(x, w), FALSE)$correlation
    stopifnot(identical(digest::digest(relation, algo = 'sha256'),
      source_bank$derived_matrix_signatures$clustering_correlation))
    local[[id]] <- list(genes = source_bank$genes, relation = relation)
    message('Reconstructed ', id, ': ', length(source_bank$genes), ' genes')
  }
  genes <- sort(unique(unlist(lapply(local, `[[`, 'genes'))))
  total <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  count <- matrix(0L, length(genes), length(genes))
  for (z in local) {
    ix <- match(z$genes, genes)
    total[ix, ix] <- total[ix, ix] + z$relation
    count[ix, ix] <- count[ix, ix] + 1L
  }
  relationship <- total / pmax(count, 1L)
  stopifnot(all(is.finite(relationship)), max(abs(relationship - t(relationship))) < 1e-12)
  tree <- hclust(as.dist(sqrt(2 * (1 - pmin(pmax(relationship, -1), 1)))), method = 'ward.D2')
  partition <- ordered_partition(tree, rank)
  programs <- sprintf('P%02d', seq_len(rank))
  sizes <- tabulate(partition, nbins = rank)
  weights <- sparseMatrix(i = seq_along(genes), j = partition, x = 1 / sizes[partition],
    dims = c(length(genes), rank), dimnames = list(genes, programs))
  stopifnot(all(rowSums(weights > 0) == 1L), max(abs(colSums(weights) - 1)) < 1e-12)
  bank <- list(genes = genes, partition = partition, tree = tree, weights = weights,
    discovery_ids = discovery_ids, excluded_ids = held_ids, signature = prep_signature,
    missing_pair_fraction = mean(count[upper.tri(count)] == 0L))
  qsave(bank, file.path(model_dir, 'program_bank.qs'))
  qsave(list(correlation = relationship, count = count, genes = genes), file.path(data_dir, 'gene_relationships.qs'))
  write_table(data.table(gene = genes, program = programs[partition]), 'gene_assignments.tsv')
  write_table(data.table(program = programs, genes = sizes), 'program_sizes.tsv')
  # The exact complete relationship matrix is plotted in Ward order without another clustering.
  gene_order <- tree$order
  h <- draw_heatmap(relationship[gene_order, gene_order], scale_rows = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
    hm_name = 'Mean relation', title = 'Training-cell gene relationships | Ward order | missing pairs = 0')
  pdf(file.path(figure_dir, 'training_gene_relationships.pdf'), width = 7, height = 6.5)
  ComplexHeatmap::draw(h)
  dev.off()
  rm(local, total, count, relationship, h, source_bank, molecular, x)
  invisible(gc())
  }
  header <- names(fread(expression_path, nrows = 0))
  symbols <- sub(' \\(ENSG.*', '', header[-1L])
  selected <- which(symbols %chin% unique(c(genes, unlist(target_genes))))
  e <- fread(expression_path, select = c(1L, selected + 1L))
  expression <- as.matrix(e[, -1L])
  rownames(expression) <- e[[1L]]
  colnames(expression) <- symbols[selected]
  stopifnot(!anyDuplicated(colnames(expression)), all(genes %chin% colnames(expression)),
    all(unlist(target_genes) %chin% colnames(expression)))
  available_ids <- intersect(cell_map$depmap_id, rownames(expression))
  baseline <- as.matrix(expression[available_ids, genes, drop = FALSE] %*% weights)
  targets <- expression[available_ids, unique(unlist(target_genes)), drop = FALSE]
  score_list <- metadata_list <- list()
  for (i in seq_len(nrow(cell_map))) {
    molecular <- qread(panel$input_paths[i])
    metadata_list[[i]] <- as.data.table(molecular$metadata)
    stopifnot(identical(rownames(molecular$response), metadata_list[[i]]$observation_id))
    scores <- as.matrix(molecular$response[, genes, drop = FALSE] %*% weights)
    direct <- sapply(seq_len(rank), function(k)
      rowMeans(molecular$response[1:3, genes[which(weights[, k] > 0)], drop = FALSE]))
    stopifnot(max(abs(scores[1:3, ] - direct)) < 1e-10)
    score_list[[i]] <- scores
  }
  measured <- do.call(rbind, score_list)
  md <- rbindlist(metadata_list, fill = TRUE)
  stopifnot(!anyDuplicated(md$observation_id), identical(rownames(measured), md$observation_id), all(is.finite(measured)))
  prepared <- list(signature = prep_signature, bank = bank, measured = measured, metadata = md,
    baseline = baseline, targets = targets)
  qsave(prepared, cache_path)
  rm(expression, e, score_list, metadata_list, molecular, prepared)
  invisible(gc())
}
programs <- colnames(measured)
rank <- length(programs)
settings <- qread(file.path(out, 'settings.qs'))
settings$rank <- rank
qsave(settings, file.path(out, 'settings.qs'))
write_table(cell_map[, .(depmap_id, depmap_name, held_out = depmap_id %chin% held_ids,
  discovery = depmap_id %chin% discovery_ids, baseline_available = depmap_id %chin% rownames(baseline))], 'cell_split.tsv')
stopifnot(!length(intersect(bank$discovery_ids, held_ids)), all(is.finite(baseline)))
# A fixed deterministic inner assignment groups all observations of each cell.
inner_ids <- sort(setdiff(rownames(baseline), held_ids))
set.seed(20261010L)
inner_map <- setNames(sample(rep(1:5, length.out = length(inner_ids))), inner_ids)
write_table(data.table(depmap_id = names(inner_map), inner_fold = unname(inner_map)), 'inner_cell_folds.tsv')

# Reused ridge fitting for scalar and multivariate targets. Outcome normalization
# is fit independently in each inner training partition. Zero-anchor models do not
# center predictors or outcomes and have no intercept.
source('/work/src/program_pilot/response_models.R')
save_plot <- function(p, name, width = 7, height = 4) {
  ggsave(file.path(figure_dir, paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(figure_dir, paste0(name, '.png')), p, width = width, height = height, dpi = 150)
}
save_molecular_heatmap <- function(bm) {
  heat_data <- bm[model == 'baseline_programs_targets',
    .(improvement = 1 - mean(standardized_mse) / mean(reference_mse)), by = .(program, drug)]
  heat_wide <- dcast(heat_data, program ~ drug, value.var = 'improvement')
  heat_matrix <- as.matrix(heat_wide[, -1L])
  rownames(heat_matrix) <- heat_wide$program
  h <- draw_heatmap(heat_matrix, scale_rows = FALSE, cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE, hm_name = 'MSE improvement',
    title = 'Baseline prediction versus training mean')
  figure_height <- max(7, .12 * nrow(heat_matrix))
  pdf(file.path(figure_dir, 'B_program_prediction_improvement.pdf'), width = 8, height = figure_height)
  ComplexHeatmap::draw(h)
  dev.off()
  png(file.path(figure_dir, 'B_program_prediction_improvement.png'), width = 1200,
    height = ceiling(figure_height * 150), res = 150)
  ComplexHeatmap::draw(h)
  dev.off()
}
# Every prediction is saved before metrics are computed. Errors are summarized
# at drug-dose and cell levels, with separate cell- and drug-block bootstrap CIs.
evaluate <- function(predictions, prefix) {
  write_table(predictions, paste0(prefix, '_predictions.tsv.gz'))
  checks <- predictions[, .(n = .N, cells = uniqueN(depmap_id), pearson = cor_safe(observed, predicted),
    mse = mean((observed - predicted)^2), mean_mse = mean((observed - reference)^2)), by = .(model, drug, dose)]
  checks[, improvement_over_mean := 1 - mse / mean_mse]
  write_table(checks, paste0(prefix, '_per_drug_dose.tsv'))
  metrics <- predictions[, .(rows = .N, pearson = cor_safe(observed, predicted),
    rmse = sqrt(mean((observed - predicted)^2)), mse = mean((observed - predicted)^2),
    reference_mse = mean((observed - reference)^2)), by = model]
  within <- copy(predictions)
  within[, `:=`(observed_within = observed - mean(observed), predicted_within = predicted - mean(predicted)),
    by = .(model, drug, dose)]
  wm <- within[, .(within_drug_dose_pearson = cor_safe(observed_within, predicted_within)), by = model]
  metrics <- merge(metrics, wm, by = 'model')
  metrics <- merge(metrics, checks[, .(equal_condition_mse = mean(mse), equal_condition_reference_mse = mean(mean_mse)),
    by = model], by = 'model')
  write_table(metrics, paste0(prefix, '_metrics.tsv'))
  intervals <- list()
  for (unit in c('depmap_id', 'drug')) {
    errors <- predictions[, .(gain = mean((observed - reference)^2 - (observed - predicted)^2)), by = c('model', unit)]
    for (name in unique(errors$model)) {
      v <- errors[model == name, gain]
      set.seed(20261010L)
      boot <- replicate(1000, mean(sample(v, replace = TRUE)))
      intervals[[length(intervals) + 1L]] <- data.table(model = name, resampling_unit = unit,
        units = length(v), mse_gain = mean(v), lower = quantile(boot, .025), upper = quantile(boot, .975))
    }
  }
  write_table(rbindlist(intervals), paste0(prefix, '_uncertainty.tsv'))
  display_names <- c(measured_D = 'Measured response',
    baseline_only = 'Baseline-only residual model',
    measured_D_interaction_only = 'Centered response + baseline interaction',
    measured_D_conditioned = 'Measured response\nBaseline-conditioned decoder',
    direct_baseline_drug_dose = 'Direct baseline by drug + dose',
    drug_dose_mean = 'Training drug-dose mean',
    direct_focal_baseline_targets = 'Direct focal baseline + targets',
    baseline_programs__measured_D = 'Predicted response (baseline)\nResponse-only decoder',
    baseline_programs__measured_D_conditioned = 'Predicted response (baseline)\nBaseline-conditioned decoder',
    baseline_programs_targets__measured_D = 'Predicted response (baseline + targets)\nResponse-only decoder',
    baseline_programs_targets__measured_D_conditioned = 'Predicted response (baseline + targets)\nBaseline-conditioned decoder')
  lab <- metrics[, .(model, display = unname(display_names[model]),
    label = sprintf('r = %.2f; within = %.2f', pearson, within_drug_dose_pearson))]
  stopifnot(!anyNA(lab$display))
  pp <- merge(predictions, lab, by = 'model')
  pp[, panel := paste(display, label, sep = '\n')]
  p <- ggplot(pp, aes(observed, predicted)) + geom_point(size = .5, alpha = .3) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
    facet_wrap(~panel, scales = 'free', ncol = 3) +
    labs(x = 'Observed log2 abundance / vehicle', y = 'Predicted log2 abundance / vehicle') +
    theme_classic(base_size = 9) + theme(aspect.ratio = 1, strip.text = element_text(size = 8))
  save_plot(p, paste0(prefix, '_predictions'), 10, 3.3 * ceiling(uniqueN(pp$panel) / 3))
  print(metrics)
  invisible(metrics)
}

# E: descriptive diagnostics requested after evaluation. Reuse existing fits;
# do not select a new model, drug, rank or gene set using these results.
run_structure_diagnostics <- function() {
  status('E: within-drug displacement and script-16 scoring audit')
  sensitivity <- as.data.table(panel$sensitivity)
  sensitivity[, split := ifelse(depmap_id %chin% held_ids, 'Held out', 'Training')]
  focal <- sensitivity[drug %chin% focal_drugs & is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  focal[, log_abundance := log2(prism_abundance_pct / 100)]
  spread <- focal[, .(cells = uniqueN(depmap_id), sensitive_cells = sum(prism_abundance_pct <= 50),
    min_pct = min(prism_abundance_pct), median_pct = median(prism_abundance_pct),
    max_pct = max(prism_abundance_pct), sd_log2 = sd(log_abundance)), by = .(drug, dose, split)]
  write_table(spread[order(drug, dose, split)], 'E_focal_sensitivity_spread.tsv')
  focal[, cell_name := cell_map$depmap_name[match(depmap_id, cell_map$depmap_id)]]
  write_table(focal, 'E_focal_sensitivity_values.tsv')
  p <- ggplot(focal, aes(factor(dose), prism_abundance_pct, colour = split)) +
    geom_hline(yintercept = 50, linetype = 2, colour = '#777777') +
    geom_point(position = position_jitter(width = .16, height = 0, seed = 20), size = 1.6, alpha = .7) +
    facet_wrap(~drug, ncol = 3) + scale_y_log10() +
    scale_colour_manual(values = c('Training' = '#999999', 'Held out' = '#B2182B')) +
    labs(x = 'Dose (uM)', y = 'PRISM abundance (% vehicle; log scale)', colour = NULL,
      subtitle = 'Each point is one cell line; dashed line = 50% vehicle') + theme_classic(base_size = 10)
  save_plot(p, 'E_focal_sensitivity', 9, 7)

  score_rows <- list()
  for (drug_name in focal_drugs) {
    all_rows <- which(md$drug == drug_name & md$depmap_id %chin% rownames(baseline))
    training_rows <- all_rows[md$depmap_id[all_rows] %chin% discovery_ids]
    # Learn only display row order from training molecular data, identical in both panels.
    row_order <- hclust(dist(t(measured[training_rows, , drop = FALSE])), method = 'ward.D2')$order
    pdf(file.path(figure_dir, paste0('E_focal_', drug_name, '_displacements.pdf')),
      width = 13, height = max(11, .13 * rank))
    for (dose_value in sort(unique(md$dose[all_rows]))) {
      rows <- all_rows[md$dose[all_rows] == dose_value]
      ids <- md$depmap_id[rows]
      fits <- qread(file.path(model_dir, paste0('B_', drug_name, '_', dose_value, '.qs')))
      xx <- cbind(baseline[ids, , drop = FALSE], targets[ids, target_genes[[drug_name]], drop = FALSE])
      predicted <- predict_ridge(fits$models$baseline_programs_targets, xx)
      actual <- measured[rows, , drop = FALSE]
      reference <- colMeans(actual[!ids %chin% held_ids, , drop = FALSE])
      an <- sensitivity[match(md$observation_id[rows], observation_id)]
      stopifnot(identical(an$observation_id, md$observation_id[rows]))
      order_columns <- order(ids %chin% held_ids, an$prism_abundance_pct, ids, na.last = TRUE)
      names_columns <- paste0(cell_map$depmap_name[match(ids, cell_map$depmap_id)], ' [',
        ifelse(is.finite(an$prism_abundance_pct), sprintf('%.0f%%', an$prism_abundance_pct), 'NA'), ']')
      annotations <- data.frame(Split = ifelse(ids %chin% held_ids, 'Held out', 'Training'),
        row.names = names_columns)
      annotations <- annotations[order_columns, , drop = FALSE]
      for (kind in c('Raw displacement', 'Deviation from training drug-dose mean')) {
        observed_matrix <- if (kind == 'Raw displacement') actual else sweep(actual, 2L, reference, '-')
        predicted_matrix <- if (kind == 'Raw displacement') predicted else sweep(predicted, 2L, reference, '-')
        observed_matrix <- t(observed_matrix)[row_order, order_columns, drop = FALSE]
        predicted_matrix <- t(predicted_matrix)[row_order, order_columns, drop = FALSE]
        colnames(observed_matrix) <- colnames(predicted_matrix) <- names_columns[order_columns]
        rownames(observed_matrix) <- rownames(predicted_matrix) <- programs[row_order]
        stopifnot(identical(colnames(observed_matrix), rownames(annotations)))
        limit <- max(.01, quantile(abs(c(observed_matrix, predicted_matrix)), .99))
        display <- list()
        for (view in c('Measured', 'Predicted')) {
          matrix <- if (view == 'Measured') observed_matrix else predicted_matrix
          display[[view]] <- draw_heatmap(matrix, scale_rows = FALSE,
            col_anno_df = annotations, col_anno_cols = 'Split',
            palette_overrides = list(Split = c('Held out' = '#B2182B', 'Training' = '#999999')),
            cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
            show_column_names = view == 'Predicted', column_names_rot = 90,
            column_names_gp = grid::gpar(fontsize = 7), row_names_gp = grid::gpar(fontsize = 7),
            hm_name = paste(view, 'delta'), value_range = c(-limit, 0, limit),
            value_colors = c('#2166AC', 'white', '#B2182B'),
            title = paste(drug_name, dose_value, 'uM |', view, '|', kind,
              if (view == 'Predicted') '\nTraining = in-sample fit; held out = prediction; labels show % vehicle' else ''))
        }
        ComplexHeatmap::draw(ComplexHeatmap::`%v%`(display$Measured, display$Predicted), merge_legends = TRUE)
      }
      score_rows[[length(score_rows) + 1L]] <- rbindlist(lapply(seq_along(ids), function(j)
        data.table(drug = drug_name, dose = dose_value, depmap_id = ids[j], program = programs,
          split = ifelse(ids[j] %chin% held_ids, 'Held out', 'Training'),
          observed = actual[j, ], predicted = predicted[j, ], training_mean = reference,
          prism_abundance_pct = an$prism_abundance_pct[j])))
    }
    dev.off()
    message('Within-drug plots saved: ', drug_name)
  }
  write_table(rbindlist(score_rows), 'E_focal_displacement_values.tsv.gz')

  # Measure whether global scores already contain sensitivity-associated variation
  # within each cell. These are descriptive correlations, not new prediction tests.
  matched <- sensitivity[is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  matched[, outcome := log2(prism_abundance_pct / 100)]
  correlations <- rbindlist(lapply(unique(matched$depmap_id), function(id) {
    rows <- matched[depmap_id == id]
    x <- measured[match(rows$observation_id, md$observation_id), , drop = FALSE]
    data.table(depmap_id = id, cell_name = cell_map$depmap_name[match(id, cell_map$depmap_id)],
      split = ifelse(id %chin% held_ids, 'Held out', 'Training'), program = programs,
      n = nrow(rows), pearson = vapply(seq_len(rank), function(k) cor_safe(x[, k], rows$outcome), numeric(1)))
  }))
  write_table(correlations, 'E_global_program_sensitivity_by_cell.tsv')
  wide <- dcast(correlations, program ~ cell_name, value.var = 'pearson')
  cm <- as.matrix(wide[, -1L])
  rownames(cm) <- wide$program
  stopifnot(all(is.finite(cm)))
  hm <- draw_heatmap(cm, scale_rows = FALSE, cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE, hm_name = 'Pearson r',
    value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'),
    title = 'Shared-bank hard scores versus sensitivity within each cell\nDescriptive across drug/dose; no model fitting')
  pdf(file.path(figure_dir, 'E_global_program_sensitivity_by_cell.pdf'), width = 12, height = max(8, .12 * rank))
  ComplexHeatmap::draw(hm)
  dev.off()

  # Isolate AN3 CA score changes using the exact cached molecular measurements.
  id <- 'ACH-000940'
  local_bank <- qread('/work/results/program_pilot/16_global_mean_relationships_fixed_ward_programs/cell_AN3_CA/models/program_bank.qs')
  molecular <- qread(panel$input_paths[match(id, cell_map$depmap_id)])
  genes <- local_bank$genes
  selected_index <- match(genes, bank$genes)
  stopifnot(!anyNA(selected_index))
  selected_weights <- bank$weights[selected_index, , drop = FALSE]
  selected_count <- colSums(selected_weights > 0)
  total_count <- colSums(bank$weights > 0)
  write_table(data.table(program = programs, AN3_selected_genes = selected_count,
    full_union_genes = total_count, selected_fraction_of_hard_mean = selected_count / total_count),
    'E_AN3_CA_selected_gene_mass.tsv')
  yy <- matched[depmap_id == id]
  x <- molecular$response[match(yy$observation_id, rownames(molecular$response)), genes, drop = FALSE]
  restricted <- sapply(which(selected_count > 0), function(k) rowMeans(x[, selected_weights[, k] > 0, drop = FALSE]))
  local_hard <- sapply(sort(unique(local_bank$partition)), function(k) rowMeans(x[, local_bank$partition == k, drop = FALSE]))
  score_sets <- list(script16_cell_specific_weighted = x %*% local_bank$weights,
    script16_selected_genes_hard = local_hard,
    script20_selected_genes_only_hard = restricted,
    script20_current_bank_hard = measured[match(yy$observation_id, md$observation_id), , drop = FALSE])
  association <- rbindlist(lapply(names(score_sets), function(name) {
    xx <- score_sets[[name]]
    rr <- vapply(seq_len(ncol(xx)), function(k) cor_safe(xx[, k], yy$outcome), numeric(1))
    data.table(scoring = name, program_index = seq_along(rr), pearson = rr, n = nrow(xx))
  }))
  write_table(association, 'E_AN3_CA_scoring_associations.tsv')
  audit <- association[, .(programs = .N, max_abs_r = max(abs(pearson)), median_abs_r = median(abs(pearson))), by = scoring]
  write_table(audit, 'E_AN3_CA_scoring_summary.tsv')

  # Quantify the separation of questions: same cell across drugs versus same
  # drug/dose across cells, using existing decoder predictions only.
  decoder <- qread(file.path(model_dir, 'A_shared_decoder.qs'))
  eligible <- matched[depmap_id %chin% rownames(baseline)]
  xx <- measured[match(eligible$observation_id, md$observation_id), , drop = FALSE]
  eligible[, predicted := predict_ridge(decoder$models$measured_D, xx)]
  by_cell <- eligible[, .(n = .N, sensitivity_min = min(prism_abundance_pct),
    sensitivity_max = max(prism_abundance_pct), r_across_drugs_doses = cor_safe(outcome, predicted)),
    by = .(depmap_id, split)]
  by_cell[, cell_name := cell_map$depmap_name[match(depmap_id, cell_map$depmap_id)]]
  write_table(by_cell, 'E_existing_decoder_within_cell.tsv')
  train <- eligible[!depmap_id %chin% held_ids]
  train[, condition := factor(paste(drug, dose, sep = '__'))]
  variance <- data.table(model = c('drug_dose_only', 'cell_only', 'drug_dose_plus_cell'),
    training_descriptive_R2 = c(summary(lm(outcome ~ condition, data = train))$r.squared,
      summary(lm(outcome ~ factor(depmap_id), data = train))$r.squared,
      summary(lm(outcome ~ condition + factor(depmap_id), data = train))$r.squared))
  write_table(variance, 'E_training_outcome_variance.tsv')
  writeLines(c('Within-drug and representation audit: descriptive; no new model fitting',
    'Script 20 main decoder: D or D + B*D, shared slopes, no intercept, no drug-dose or cell fixed effects.',
    'Balancing observations changes weights; it does not remove drug or cell effects.',
    'The direct-baseline comparator has drug-dose intercepts and drug-specific baseline slopes.',
    'Baseline-to-D models are fitted separately per drug/dose; held cells are excluded.',
    sprintf('Script 16 AN3 CA uses %d locally outcome-selected genes and %d active globally assigned programs.',
      length(genes), ncol(local_bank$weights)),
    sprintf('Script 20 bank %s: %d programs over %d genes; median AN3-selected weight mass %.1f%%.',
      program_bank_type, rank, length(bank$genes), 100 * median(selected_count / total_count)),
    'Script 16 heatmap displays individual genes; its column scores use cell-specific soft weights.',
    'Script 20 uses common hard means. Global assignments alone do not make these scores identical.',
    'Within-cell sensitivity structure across drugs is not the same test as within-drug transfer across cells.',
    'AN3 CA gene-screen associations below are descriptive and outcome-selected, not validation.',
    '', capture.output(print(audit)), '', 'Existing decoder within AN3 CA:',
    capture.output(print(by_cell[depmap_id == id])), '', 'Training outcome variance:',
    capture.output(print(variance)), '', 'Focal sensitivity support:', capture.output(print(spread[order(drug, dose, split)]))),
    file.path(out, 'WITHIN_DRUG_DIAGNOSTICS.txt'))
  status('COMPLETE: within-drug plots and script-16 audit saved; no models refitted')
}
# F: predict variation across cells at fixed drug/dose. All means, baseline
# scaling and ridge penalties are estimated using training cells only.
run_drug_dose_adjusted <- function() {
  status('F: drug-dose-adjusted sensitivity and its baseline/molecular determinants')
  labels <- as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
    depmap_id %chin% rownames(baseline)]
  labels[, `:=`(observed = log2(prism_abundance_pct / 100), condition = paste(drug, dose, sep = '__'))]
  train <- which(!labels$depmap_id %chin% held_ids)
  test <- which(labels$depmap_id %chin% held_ids)
  response <- measured[match(labels$observation_id, md$observation_id), , drop = FALSE]
  b0 <- baseline[labels$depmap_id, , drop = FALSE]
  stopifnot(all(is.finite(response)), !anyDuplicated(labels$observation_id))
  # This transform is reused at inner validation, final training and chain inference.
  fit_adjustment <- function(rows) {
    groups <- split(rows, labels$condition[rows])
    dm <- do.call(rbind, lapply(groups, function(i) colMeans(response[i, , drop = FALSE])))
    ym <- vapply(groups, function(i) mean(labels$observed[i]), numeric(1))
    cells <- unique(labels$depmap_id[rows])
    mu <- colMeans(baseline[cells, , drop = FALSE])
    sc <- apply(baseline[cells, , drop = FALSE], 2L, sd)
    stopifnot(all(sc > 0), identical(rownames(dm), names(ym)))
    list(displacement_mean = dm, sensitivity_mean = ym, baseline_mean = mu,
      baseline_sd = sc, training_ids = cells, training_observations = labels$observation_id[rows])
  }
  transform_inputs <- function(adjustment, rows, input_response = response[rows, , drop = FALSE]) {
    k <- match(labels$condition[rows], rownames(adjustment$displacement_mean))
    stopifnot(!anyNA(k))
    list(d = input_response - adjustment$displacement_mean[k, , drop = FALSE],
      b = sweep(sweep(b0[rows, , drop = FALSE], 2L, adjustment$baseline_mean, '-'),
        2L, adjustment$baseline_sd, '/'), reference = unname(adjustment$sensitivity_mean[k]))
  }
  make_design <- function(z, name) {
    if (name == 'baseline_only') return(z$b)
    if (name == 'measured_D') return(z$d)
    if (name == 'measured_D_interaction_only') return(cbind(z$d, z$d * z$b))
    stopifnot(name == 'measured_D_conditioned')
    cbind(z$d, z$b, z$d * z$b)
  }
  kinds <- c('baseline_only', 'measured_D', 'measured_D_interaction_only', 'measured_D_conditioned')
  fold <- unname(inner_map[labels$depmap_id[train]])
  stopifnot(!anyNA(fold))
  w <- observation_weights(labels[train])
  loss <- matrix(0, length(kinds), length(lambda_grid), dimnames = list(kinds, lambda_grid))
  cv_predictions <- cv_checks <- list()
  for (k in sort(unique(fold))) {
    it <- train[fold != k]
    va <- train[fold == k]
    adjustment <- fit_adjustment(it)
    xt <- transform_inputs(adjustment, it)
    xv <- transform_inputs(adjustment, va)
    yt <- labels$observed[it] - xt$reference
    yv <- labels$observed[va] - xv$reference
    stopifnot(!any(labels$depmap_id[va] %chin% adjustment$training_ids))
    cv_checks[[as.character(k)]] <- data.table(inner_fold = k, train_cells = length(adjustment$training_ids),
      validation_cells = uniqueN(labels$depmap_id[va]), overlapping_cells = 0L)
    for (name in kinds) {
      f <- fit_path(make_design(xt, name), yt, observation_weights(labels[it]), intercept = TRUE)
      pp <- pred_path(f, make_design(xv, name))
      loss[name, ] <- loss[name, ] + colSums((pp - yv)^2 * w[fold == k])
      cv_predictions[[paste(k, name)]] <- rbindlist(lapply(seq_along(lambda_grid), function(j)
        data.table(observation_id = labels$observation_id[va], depmap_id = labels$depmap_id[va],
          inner_fold = k, model = name, lambda = lambda_grid[j], observed_residual = yv,
          predicted_residual = pp[, j], reference = xv$reference)))
    }
    message('F: completed inner cell fold ', k)
  }
  write_table(rbindlist(cv_checks), 'F_inner_split_checks.tsv')
  write_table(rbindlist(cv_predictions), 'F_inner_predictions.tsv.gz')
  adjustment <- fit_adjustment(train)
  stopifnot(!any(adjustment$training_ids %chin% held_ids))
  xt <- transform_inputs(adjustment, train)
  xv <- transform_inputs(adjustment, test)
  labels[, reference := unname(adjustment$sensitivity_mean[match(condition, names(adjustment$sensitivity_mean))])]
  check_means <- labels[train, .(mean_residual = mean(observed - reference)), by = condition]
  stopifnot(max(abs(check_means$mean_residual)) < 1e-12)
  write_table(check_means, 'F_training_residual_mean_checks.tsv')
  write_table(data.table(condition = names(adjustment$sensitivity_mean),
    sensitivity_mean_log2 = unname(adjustment$sensitivity_mean)), 'F_training_drug_dose_means.tsv')
  write_table(cbind(data.table(condition = rownames(adjustment$displacement_mean)),
    as.data.table(adjustment$displacement_mean)), 'F_training_displacement_means.tsv.gz')
  fits <- predictions <- coefficients <- list()
  for (name in kinds) {
    f <- fit_path(make_design(xt, name), labels$observed[train] - xt$reference, w, intercept = TRUE)
    f$lambda <- lambda_grid[which.min(loss[name, ])]
    f$inner_mse <- loss[name, ] / sum(w)
    f$training_ids <- adjustment$training_ids
    fits[[name]] <- f
    z <- labels[test, .(observation_id, depmap_id, drug, dose, observed, reference)]
    z[, `:=`(model = name, predicted = reference + predict_ridge(f, make_design(xv, name)))]
    predictions[[name]] <- z
    qsave(list(model = f, adjustment = adjustment), file.path(model_dir, paste0('F_', name, '.qs')))
    coef <- as.numeric(coef(f$fit, s = f$lambda))
    terms <- if (name == 'baseline_only') paste0('baseline:', programs) else if (name == 'measured_D')
      paste0('displacement:', programs) else if (name == 'measured_D_interaction_only')
        c(paste0('displacement:', programs), paste0('interaction:', programs)) else
          c(paste0('displacement:', programs), paste0('baseline:', programs), paste0('interaction:', programs))
    stopifnot(length(coef) == length(terms) + 1L)
    coefficients[[name]] <- data.table(model = name, term = c('(Intercept)', terms), coefficient = coef)
    write_table(rbindlist(predictions), 'F_partial_predictions.tsv.gz')
  }
  predictions$drug_dose_mean <- labels[test, .(observation_id, depmap_id, drug, dose, observed, reference,
    model = 'drug_dose_mean', predicted = reference)]
  # The previous direct baseline comparator is retained unchanged for reference.
  old <- fread(file.path(data_dir, 'A_shared_function_predictions.tsv.gz'))
  predictions$direct_baseline_drug_dose <- old[model == 'direct_baseline_drug_dose']
  fp <- rbindlist(predictions)
  metrics <- evaluate(fp, 'F_adjusted_function')
  write_table(rbindlist(coefficients), 'F_function_coefficients.tsv')
  qsave(list(models = fits, adjustment = adjustment), file.path(model_dir, 'F_adjusted_decoder.qs'))
  write_table(rbindlist(lapply(kinds, function(name) data.table(model = name, lambda = lambda_grid,
    inner_mse = loss[name, ] / sum(w), selected = lambda_grid == fits[[name]]$lambda))), 'F_selected_penalties.tsv')

  # Freeze the new measured-response decoder and substitute the existing baseline
  # response predictions, with the SAME training-only drug-dose centering.
  previous_chain <- fread(file.path(data_dir, 'C_full_chain_predictions.tsv.gz'))
  matched_ids <- unique(previous_chain$observation_id)
  rows <- which(labels$observation_id %chin% matched_ids)
  stopifnot(all(rows %in% test))
  bp <- fread(file.path(data_dir, 'B_molecular_predictions.tsv.gz'))[model == 'baseline_programs_targets']
  bp <- dcast(bp, depmap_id + drug + dose ~ program, value.var = 'predicted')
  key <- function(x) paste(x$depmap_id, x$drug, x$dose, sep = '__')
  ix <- match(key(labels[rows]), key(bp))
  stopifnot(!anyNA(ix))
  input_response <- as.matrix(bp[ix, ..programs])
  transformed <- transform_inputs(adjustment, rows, input_response)
  chain <- list()
  for (name in c('measured_D', 'measured_D_conditioned')) {
    z <- labels[rows, .(observation_id, depmap_id, drug, dose, observed, reference)]
    z[, `:=`(model = paste0('baseline_programs_targets__', name),
      predicted = reference + predict_ridge(fits[[name]], make_design(transformed, name)))]
    chain[[name]] <- z
  }
  cp <- rbind(rbindlist(chain), fp[observation_id %chin% matched_ids],
    previous_chain[model == 'direct_focal_baseline_targets'])
  stopifnot(!anyDuplicated(cp[, .(model, observation_id)]), length(unique(cp[, .N, by = model]$N)) == 1L)
  chain_metrics <- evaluate(cp, 'F_adjusted_chain')
  # Explicit residual plots prevent between-drug potency from inflating the visual association.
  residual <- copy(fp)
  residual[, `:=`(observed_residual = observed - reference, predicted_residual = predicted - reference)]
  residual_summary <- residual[, .(pearson = cor_safe(observed_residual, predicted_residual),
    mse = mean((observed_residual - predicted_residual)^2)), by = model]
  write_table(residual_summary, 'F_residual_metrics.tsv')
  p <- ggplot(residual, aes(observed_residual, predicted_residual)) + geom_point(size = .5, alpha = .3) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = '#B2182B', linewidth = .5) +
    facet_wrap(~model, scales = 'free', ncol = 3, labeller = label_wrap_gen(25)) +
    labs(x = 'Observed deviation from training drug-dose mean (log2)',
      y = 'Predicted deviation from training drug-dose mean (log2)') +
    theme_classic(base_size = 9) + theme(aspect.ratio = 1)
  save_plot(p, 'F_sensitivity_residual_predictions', 10, 7)
  comparison <- merge(fread(file.path(data_dir, 'A_shared_function_metrics.tsv')),
    metrics, by = 'model', suffixes = c('_original', '_adjusted'))
  write_table(comparison, 'F_adjustment_comparison.tsv')
  writeLines(c('Drug-dose-adjusted functional learning on the script-15 union bank',
    'Training target: log2(PRISM abundance / vehicle) minus its training drug-dose mean.',
    'Displacement features: native program scores minus their training drug-dose means.',
    'Means and baseline scaling are refit inside every inner cell fold; held cells never estimate them.',
    'Final prediction adds the training sensitivity mean back to the predicted residual.',
    'Baseline-only tests broad cell susceptibility; D tests molecular deviations.',
    'Interaction-only preserves D + B*D; conditioned additionally includes baseline main effects.',
    'No cell-ID coefficients are fitted: unseen-cell differences must be inferred from baseline.',
    'Same bank, held cells, upstream models, observations and ridge grid as the preceding union-bank run.',
    'Coefficients are regularized associations, not causal gene/program effects.',
    '', 'Measured response:', capture.output(print(metrics)), '', 'Complete chain:',
    capture.output(print(chain_metrics))), file.path(out, 'DRUG_DOSE_ADJUSTED_REPORT.txt'))
  status('COMPLETE: drug-dose-adjusted function and full chain saved')
}
if ('--adjusted-only' %in% commandArgs(TRUE)) {
  stopifnot(identical(commandArgs(TRUE), '--adjusted-only'))
  run_drug_dose_adjusted()
  quit(save = 'no', status = 0)
}

if ('--diagnostics-only' %in% commandArgs(TRUE)) {
  stopifnot(identical(commandArgs(TRUE), '--diagnostics-only'))
  run_structure_diagnostics()
  quit(save = 'no', status = 0)
}

# A: the same shared measured-response decoder is used subsequently in C.
report_only <- '--report-only' %in% commandArgs(TRUE)
combination_only <- '--combination-only' %in% commandArgs(TRUE) || report_only
stopifnot(all(commandArgs(TRUE) %in% c('--combination-only', '--report-only')))
if (!combination_only) {
status('A: shared measured-response functional models')
labels <- as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
labels[, observed := log2(prism_abundance_pct / 100)]
labels[, baseline_available := depmap_id %chin% rownames(baseline)]
write_table(labels[, .(rows = .N), by = .(depmap_id, baseline_available, prism_status)], 'prism_coverage.tsv')
labels <- labels[baseline_available == TRUE]
stopifnot(!anyDuplicated(labels$observation_id), all(labels$observation_id %chin% md$observation_id))
labels[, condition := paste(drug, dose, sep = '__')]
tr <- which(!labels$depmap_id %chin% held_ids)
te <- which(labels$depmap_id %chin% held_ids)
stopifnot(length(tr) > 0, length(te) > 0)
write_table(labels[, .(rows = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug)),
  by = .(held_out = depmap_id %chin% held_ids)], 'prism_split_counts.tsv')
d <- measured[match(labels$observation_id, md$observation_id), , drop = FALSE]
bmu <- colMeans(baseline[unique(labels$depmap_id[tr]), , drop = FALSE])
bsd <- apply(baseline[unique(labels$depmap_id[tr]), , drop = FALSE], 2L, sd)
stopifnot(all(bsd > 0))
b <- sweep(sweep(baseline, 2L, bmu, '-'), 2L, bsd, '/')
bl <- b[labels$depmap_id, , drop = FALSE]
condition_levels <- sort(unique(labels$condition[tr]))
stopifnot(all(labels$condition[te] %chin% condition_levels))
condition_design <- sparseMatrix(i = seq_len(nrow(labels)), j = match(labels$condition, condition_levels), x = 1,
  dims = c(nrow(labels), length(condition_levels)))
drug_levels <- sort(unique(labels$drug[tr]))
# Strong direct comparator: drug-specific baseline effects plus drug-dose intercepts.
interaction_design <- sparseMatrix(i = rep(seq_len(nrow(labels)), rank),
  j = rep((match(labels$drug, drug_levels) - 1L) * rank, rank) + rep(seq_len(rank), each = nrow(labels)),
  x = as.vector(bl), dims = c(nrow(labels), length(drug_levels) * rank))
a_design <- list(measured_D = d, measured_D_conditioned = cbind(d, d * bl),
  direct_baseline_drug_dose = cbind(condition_design, interaction_design))
means <- labels[tr, .(reference = mean(observed)), by = condition]
labels[, reference := means$reference[match(condition, means$condition)]]
w <- observation_weights(labels[tr])
a_models <- list()
a_predictions <- list()
for (name in names(a_design)) {
  status(paste('A: fitting', name))
  xx <- a_design[[name]]
  f <- fit_ridge(xx[tr, , drop = FALSE], labels$observed[tr], labels$depmap_id[tr], w,
    intercept = name == 'direct_baseline_drug_dose')
  if (name != 'direct_baseline_drug_dose') stopifnot(max(abs(predict_ridge(f, xx[te[1:2], , drop = FALSE] * 0))) < 1e-12)
  a_models[[name]] <- f
  qsave(f, file.path(model_dir, paste0('A_', name, '.qs')))
  a_predictions[[name]] <- labels[te, .(observation_id, depmap_id, drug, dose, observed, reference)]
  a_predictions[[name]][, `:=`(model = name, predicted = predict_ridge(f, xx[te, , drop = FALSE]))]
  write_table(rbindlist(a_predictions), 'A_partial_predictions.tsv.gz')
}
a_predictions[['drug_dose_mean']] <- labels[te, .(observation_id, depmap_id, drug, dose, observed, reference,
  model = 'drug_dose_mean', predicted = reference)]
a_predictions <- rbindlist(a_predictions)
a_metrics <- evaluate(a_predictions, 'A_shared_function')
qsave(list(models = a_models, baseline_center = bmu, baseline_scale = bsd), file.path(model_dir, 'A_shared_decoder.qs'))
# Inspect all programs in the training molecular signatures before interpreting fits.
focus_train <- which(md$drug %chin% focal_drugs & md$depmap_id %chin% discovery_ids)
h <- draw_heatmap(t(measured[focus_train, , drop = FALSE]), scale_rows = FALSE,
  cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = FALSE,
  hm_name = 'Displacement', title = paste(rank, 'hard programs | seven focal drugs | training cells'))
pdf(file.path(figure_dir, 'training_focal_displacements.pdf'), width = 11, height = max(6, .12 * rank))
ComplexHeatmap::draw(h)
dev.off()

# B/C: fit molecular conditioning only for focal drugs, using all training molecular
# cells with CCLE baseline. Direct function controls use that drug-dose's PRISM rows.
status('B/C: focal baseline models and frozen decoder composition')
b_predictions <- c_predictions <- support <- list()
for (drug_name in focal_drugs) {
  doses <- sort(unique(md[drug == drug_name, dose]))
  stopifnot(length(doses) > 0)
  for (dose_value in doses) {
    rows <- which(md$drug == drug_name & md$dose == dose_value & md$depmap_id %chin% rownames(baseline))
    ids <- md$depmap_id[rows]
    stopifnot(!anyDuplicated(ids))
    itr <- which(!ids %chin% held_ids)
    ite <- which(ids %chin% held_ids)
    stopifnot(length(itr) >= 8L, length(ite) > 0)
    y <- measured[rows, , drop = FALSE]
    train_mean <- colMeans(y[itr, , drop = FALSE])
    train_sd <- apply(y[itr, , drop = FALSE], 2L, sd)
    feature <- list(baseline_programs = baseline[ids, , drop = FALSE],
      baseline_programs_targets = cbind(baseline[ids, , drop = FALSE], targets[ids, target_genes[[drug_name]], drop = FALSE]))
    models <- list()
    for (name in names(feature)) {
      f <- fit_ridge(feature[[name]][itr, , drop = FALSE], y[itr, , drop = FALSE], ids[itr])
      prediction <- predict_ridge(f, feature[[name]][ite, , drop = FALSE])
      models[[name]] <- f
      bp <- rbindlist(lapply(seq_along(ite), function(j) data.table(drug = drug_name, dose = dose_value,
        depmap_id = ids[ite[j]], program = programs, model = name, observed = y[ite[j], ],
        predicted = prediction[j, ], reference = train_mean, training_sd = train_sd)))
      b_predictions[[length(b_predictions) + 1L]] <- bp
      test_labels <- labels[drug == drug_name & dose == dose_value & depmap_id %chin% ids[ite]]
      if (nrow(test_labels)) {
        query <- match(test_labels$depmap_id, ids[ite])
        predicted_d <- prediction[query, , drop = FALSE]
        base <- b[test_labels$depmap_id, , drop = FALSE]
        for (decoder in c('measured_D', 'measured_D_conditioned')) {
          xx <- if (decoder == 'measured_D') predicted_d else cbind(predicted_d, predicted_d * base)
          cp <- test_labels[, .(observation_id, depmap_id, drug, dose, observed, reference)]
          cp[, `:=`(model = paste(name, decoder, sep = '__'), predicted = predict_ridge(a_models[[decoder]], xx))]
          c_predictions[[length(c_predictions) + 1L]] <- cp
        }
      }
    }
    # Direct per-drug/dose baseline prediction, with exactly matched test rows.
    lt <- labels[drug == drug_name & dose == dose_value & !depmap_id %chin% held_ids]
    le <- labels[drug == drug_name & dose == dose_value & depmap_id %chin% held_ids]
    if (nrow(lt) >= 5L && nrow(le)) {
      bx <- cbind(baseline, targets[, target_genes[[drug_name]], drop = FALSE])
      f <- fit_ridge(bx[lt$depmap_id, , drop = FALSE], lt$observed, lt$depmap_id)
      cp <- le[, .(observation_id, depmap_id, drug, dose, observed, reference)]
      cp[, `:=`(model = 'direct_focal_baseline_targets', predicted = predict_ridge(f, bx[le$depmap_id, , drop = FALSE]))]
      c_predictions[[length(c_predictions) + 1L]] <- cp
      models[['direct_function']] <- f
    }
    key <- paste(drug_name, dose_value, sep = '_')
    qsave(list(models = models, drug = drug_name, dose = dose_value, training_ids = ids[itr], test_ids = ids[ite]),
      file.path(model_dir, paste0('B_', key, '.qs')))
    support[[key]] <- data.table(drug = drug_name, dose = dose_value, molecular_train = length(itr),
      molecular_test = length(ite), function_train = nrow(lt), function_test = nrow(le))
    write_table(rbindlist(b_predictions), 'B_molecular_predictions.tsv.gz')
    write_table(rbindlist(c_predictions), 'C_partial_predictions.tsv.gz')
    status(paste('B/C: finished', key))
  }
}
bp <- rbindlist(b_predictions)
write_table(rbindlist(support), 'B_training_support.tsv')
bm <- bp[, .(cells = uniqueN(depmap_id), standardized_mse = mean(((observed - predicted) / training_sd)^2),
  reference_mse = mean(((observed - reference) / training_sd)^2), pearson = cor_safe(observed, predicted)),
  by = .(model, drug, dose, program)]
bm[, improvement := 1 - standardized_mse / reference_mse]
write_table(bm, 'B_molecular_metrics.tsv')
write_table(bp[, .(standardized_mse = mean(((observed - predicted) / training_sd)^2),
  reference_mse = mean(((observed - reference) / training_sd)^2)), by = .(model, drug)], 'B_by_drug.tsv')
save_molecular_heatmap(bm)
cp <- rbindlist(c_predictions)
common <- Reduce(intersect, split(cp$observation_id, cp$model))
write_table(cp[!observation_id %chin% common], 'C_unmatched_rows.tsv.gz')
cp <- cp[observation_id %chin% common]
cp <- rbind(cp, a_predictions[observation_id %chin% common], fill = TRUE)
c_metrics <- evaluate(cp, 'C_full_chain')
status('A/B/C complete; starting combination transfer')
} else {
  status('resuming D using the saved A/B/C results and frozen bank')
  decoder <- qread(file.path(model_dir, 'A_shared_decoder.qs'))
  stopifnot(all(vapply(decoder$models, function(f) !any(f$training_ids %chin% held_ids), logical(1))))
  b <- sweep(sweep(baseline, 2L, decoder$baseline_center, '-'), 2L, decoder$baseline_scale, '/')
  a_metrics <- fread(file.path(data_dir, 'A_shared_function_metrics.tsv'))
  c_metrics <- fread(file.path(data_dir, 'C_full_chain_metrics.tsv'))
  if (report_only) {
    status('regenerating plots from frozen predictions; no model fitting')
    a_metrics <- evaluate(fread(file.path(data_dir, 'A_shared_function_predictions.tsv.gz')), 'A_shared_function')
    c_metrics <- evaluate(fread(file.path(data_dir, 'C_full_chain_predictions.tsv.gz')), 'C_full_chain')
    save_molecular_heatmap(fread(file.path(data_dir, 'B_molecular_metrics.tsv')))
  }
}

# D: reuse pairing, endpoints, pair folds and symmetric-product machinery; replace
# the historical nine-program features with the new forty hard cluster means.
if (!report_only) {
combo <- readRDS(combo_path)
cd <- as.data.table(combo$d)
stopifnot(max(abs(cd$interaction - (cd$y - cd$yA - cd$yB))) < 1e-10)
query <- unique(rbind(cd[, .(depmap_id, drug = drugA, dose = doseA)],
  cd[, .(depmap_id, drug = drugB, dose = doseB)]))
query[, molecular_drug := panel_name(drug)]
query[, `:=`(query_id = .I, lower_dose = NA_real_, upper_dose = NA_real_, alpha = NA_real_,
  matching = 'no measured cell/drug')]
zs <- matrix(NA_real_, nrow(query), rank, dimnames = list(NULL, programs))
for (i in seq_len(nrow(query))) {
  rows <- which(md$depmap_id == query$depmap_id[i] & md$drug == query$molecular_drug[i])
  if (!length(rows)) next
  rows <- rows[order(md$dose[rows])]
  doses <- md$dose[rows]
  target <- query$dose[i]
  exact <- which(abs(log(target / doses)) < 1e-8)
  if (length(exact)) {
    lo <- hi <- exact[1L]
    alpha <- 0
    matching <- 'exact'
  } else {
    query[i, matching := 'outside measured cell/drug dose range']
    if (target < min(doses) || target > max(doses)) next
    lo <- max(which(doses < target))
    hi <- min(which(doses > target))
    alpha <- log(target / doses[lo]) / log(doses[hi] / doses[lo])
    matching <- 'interpolated'
  }
  zs[i, ] <- (1 - alpha) * measured[rows[lo], ] + alpha * measured[rows[hi], ]
  set(query, i, c('lower_dose', 'upper_dose', 'alpha', 'matching'), list(doses[lo], doses[hi], alpha, matching))
}
write_table(query, 'D_molecular_dose_matching.tsv')
key <- function(cell, drug, dose) paste(cell, drug, dose, sep = '__')
qkey <- key(query$depmap_id, query$drug, query$dose)
qa <- match(key(cd$depmap_id, cd$drugA, cd$doseA), qkey)
qb <- match(key(cd$depmap_id, cd$drugB, cd$doseB), qkey)
stopifnot(!anyNA(qa), !anyNA(qb))
keep <- is.finite(zs[qa, 1]) & is.finite(zs[qb, 1]) & cd$depmap_id %chin% rownames(baseline)
write_table(cd[!keep], 'D_excluded_rows.tsv.gz')
a <- zs[qa[keep], , drop = FALSE]
z <- zs[qb[keep], , drop = FALSE]
cd <- cd[keep]
legacy <- new.env(parent = globalenv())
sys.source('/work/almanac_broad9_2026-10-07/model.R', envir = legacy)
legacy$P <- programs
products <- legacy$bilinear(a, z)
stopifnot(max(abs(products - legacy$bilinear(z, a))) < 1e-10,
  max(abs(legacy$bilinear(a * 0, z))) < 1e-12)
identity_levels <- sort(unique(c(cd$idA, cd$idB)))
identity <- matrix(0, nrow(cd), length(identity_levels), dimnames = list(NULL, identity_levels))
identity[cbind(seq_len(nrow(cd)), match(cd$idA, identity_levels))] <- 1
identity[cbind(seq_len(nrow(cd)), match(cd$idB, identity_levels))] <- 1
# Symmetric single-agent potency controls use actual ALMANAC single outcomes,
# giving a stringent comparator on the same information available for the residual.
potency <- cbind(cd$yA + cd$yB, abs(cd$yA - cd$yB), cd$yA * cd$yB)
d_design <- list(molecular_products = products,
  molecular_products_conditioned = cbind(products, a * z * b[cd$depmap_id, , drop = FALSE]),
  identity_potency = cbind(identity, potency))
raw <- fread(control_path, select = c('source_row', 'single_groupA', 'single_groupB'))
refs <- unique(rbind(raw[, .(source_row, ref = single_groupA)], raw[, .(source_row, ref = single_groupB)]))
refs <- refs[source_row %in% cd$source_row]
refs[, row := match(source_row, cd$source_row)]
stopifnot(all(cd$source_row %in% refs$source_row), !anyNA(refs$ref))
purge <- function(train, test) setdiff(train, refs[ref %chin% refs[row %in% test, ref], row])
base_train <- which(!cd$depmap_id %chin% held_ids)
base_test <- which(cd$depmap_id %chin% held_ids)
splits <- list(held_cells = list(train = purge(base_train, base_test), test = base_test))
for (k in 1:5) {
  test <- base_test[cd$outer_fold[base_test] == k]
  train <- base_train[cd$outer_fold[base_train] != k]
  splits[[paste0('held_cells_pairs_', k)]] <- list(train = purge(train, test), test = test)
}
d_predictions <- d_audits <- list()
for (split_name in names(splits)) {
  train <- splits[[split_name]]$train
  test <- splits[[split_name]]$test
  stopifnot(!length(intersect(cd$depmap_id[train], cd$depmap_id[test])), length(train) > 100L,
    !length(intersect(refs[row %in% train, ref], refs[row %in% test, ref])))
  if (split_name != 'held_cells') stopifnot(!length(intersect(cd$pair[train], cd$pair[test])))
  w <- observation_weights(cd[train], 'pair')
  split_models <- list()
  for (name in names(d_design)) {
    status(paste('D:', split_name, name))
    xx <- d_design[[name]]
    intercept <- name == 'identity_potency'
    loss <- numeric(length(lambda_grid))
    for (cell in sort(unique(cd$depmap_id[train]))) {
      va <- train[cd$depmap_id[train] == cell]
      it <- purge(setdiff(train, va), va)
      stopifnot(length(it) > 100L)
      inner <- fit_path(xx[it, , drop = FALSE], cd$interaction[it], observation_weights(cd[it], 'pair'), intercept)
      pp <- pred_path(inner, xx[va, , drop = FALSE])
      loss <- loss + colSums((pp - cd$interaction[va])^2 * w[match(va, train)])
    }
    f <- fit_path(xx[train, , drop = FALSE], cd$interaction[train], w, intercept)
    f$lambda <- lambda_grid[which.min(loss)]
    f$inner_mse <- loss / sum(w)
    f$training_ids <- unique(cd$depmap_id[train])
    f$training_pairs <- unique(cd$pair[train])
    if (!intercept) stopifnot(max(abs(predict_ridge(f, xx[test[1:2], , drop = FALSE] * 0))) < 1e-12)
    split_models[[name]] <- f
    dp <- cd[test, .(source_row, depmap_id, pair, CELLNAME, drugA, drugB, doseA, doseB,
      observed = interaction, combination_observed = y, additive = yA + yB)]
    dp[, `:=`(split = split_name, model = name, predicted = predict_ridge(f, xx[test, , drop = FALSE]))]
    d_predictions[[length(d_predictions) + 1L]] <- dp
  }
  dp <- cd[test, .(source_row, depmap_id, pair, CELLNAME, drugA, drugB, doseA, doseB,
    observed = interaction, combination_observed = y, additive = yA + yB)]
  dp[, `:=`(split = split_name, model = 'independent_effects', predicted = 0)]
  d_predictions[[length(d_predictions) + 1L]] <- dp
  qsave(split_models, file.path(model_dir, paste0('D_', split_name, '.qs')))
  write_table(rbindlist(d_predictions), 'D_partial_predictions.tsv.gz')
  d_audits[[split_name]] <- data.table(split = split_name, training_rows = length(train), test_rows = length(test),
    training_cells = uniqueN(cd$depmap_id[train]), test_cells = uniqueN(cd$depmap_id[test]),
    training_pairs = uniqueN(cd$pair[train]), test_pairs = uniqueN(cd$pair[test]),
    overlapping_pairs = length(intersect(unique(cd$pair[train]), unique(cd$pair[test]))), shared_controls = 0L)
  write_table(rbindlist(d_audits), 'D_split_audit.tsv')
}
dp <- rbindlist(d_predictions)
dp[, scope := ifelse(split == 'held_cells', 'held_cells', 'held_cells_and_pairs')]
dp[, combination_predicted := additive + predicted]
write_table(dp, 'D_combination_predictions.tsv.gz')
dmetrics <- dp[, .(rows = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
  interaction_pearson = cor_safe(observed, predicted), interaction_mse = mean((observed - predicted)^2),
  zero_interaction_mse = mean(observed^2), combination_pearson = cor_safe(combination_observed, combination_predicted)),
  by = .(scope, model)]
write_table(dmetrics, 'D_combination_metrics.tsv')
} else {
  dp <- fread(file.path(data_dir, 'D_combination_predictions.tsv.gz'))
  dmetrics <- fread(file.path(data_dir, 'D_combination_metrics.tsv'))
}
pe <- dp[, .(mse = mean((observed - predicted)^2), zero_mse = mean(observed^2)), by = .(scope, model, pair, depmap_id)]
write_table(pe, 'D_pair_cell_errors.tsv')
write_table(pe[, .(equal_pair_cell_mse = mean(mse), zero_mse = mean(zero_mse)), by = .(scope, model)],
  'D_balanced_metrics.tsv')
di <- list()
for (unit in c('pair', 'depmap_id')) {
  err <- pe[, .(gain = mean(zero_mse - mse)), by = c('scope', 'model', unit)]
  for (sc in unique(err$scope)) for (name in unique(err$model)) {
    v <- err[scope == sc & model == name, gain]
    set.seed(20261010L)
    boot <- replicate(1000, mean(sample(v, replace = TRUE)))
    di[[length(di) + 1L]] <- data.table(scope = sc, model = name, resampling_unit = unit, units = length(v),
      mse_gain = mean(v), lower = quantile(boot, .025), upper = quantile(boot, .975))
  }
}
write_table(rbindlist(di), 'D_uncertainty.tsv')
plot_d <- merge(dp, dmetrics[, .(scope, model, r = interaction_pearson)], by = c('scope', 'model'))
p <- ggplot(plot_d, aes(observed, predicted)) + geom_point(size = .4, alpha = .2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B2182B', linewidth = .5) +
  facet_grid(scope ~ model, scales = 'free', labeller = label_wrap_gen(width = 20)) +
  geom_text(data = unique(plot_d[, .(scope, model, r)]),
    aes(x = -Inf, y = Inf, label = sprintf('r = %.2f', r)), hjust = -.1, vjust = 1.2, inherit.aes = FALSE, size = 3) +
  labs(x = 'Observed interaction (log2 residual)', y = 'Predicted interaction') +
  theme_classic(base_size = 9) + theme(aspect.ratio = 1, strip.text = element_text(size = 8))
save_plot(p, 'D_interaction_transfer', 12, 6.5)
print(dmetrics)

# Compact report: observed results, explicit endpoint/split limitations, and files.
bm_summary <- fread(file.path(data_dir, 'B_by_drug.tsv'))[model == 'baseline_programs_targets']
bm_summary[, improvement := 1 - standardized_mse / reference_mse]
report <- c(paste('Multidrug baseline → displacement → function:', rank, 'hard programs |', program_bank_type),
  sprintf('Discovery: %d training cells; %d genes; %d held cells excluded.',
    length(discovery_ids), length(bank$genes), length(held_ids)),
  'No outer held-cell outcome informed discovery, penalty selection or fitting.',
  'All validation is retrospective; historical test outcomes were previously inspected.',
  'Inner penalty validation conditions on the frozen outer program bank.',
  'PRISM: log2 abundance/vehicle, cached bracket interpolation, doses retained separately.',
  'Baseline: CCLE expression hard means; molecular response: original corrected displacement hard means.',
  'A: shared function mapping fitted on measured displacement; no drug labels in molecular decoders.',
  'B: each focal drug/dose learns its own baseline response map.',
  'C: identical A decoder receives predicted rather than measured displacement; no refitting.',
  'D: ALMANAC interaction residual; total predictions add OBSERVED matched single-agent effects.',
  'D is a measured-single-response interaction test, not a baseline-only combination prediction.',
  'Cell and drug/pair uncertainty are reported separately; three ALMANAC test cells limit inference.',
  '', 'Interpretation',
  sprintf('Measured D: within-drug/dose r = %.3f; RMSE = %.3f versus drug-dose mean %.3f.',
    a_metrics[model == 'measured_D', within_drug_dose_pearson],
    a_metrics[model == 'measured_D', rmse], a_metrics[model == 'drug_dose_mean', rmse]),
  sprintf('Baseline + targets improves molecular MSE by %.1f%% to %.1f%% across focal drugs.',
    100 * min(bm_summary$improvement), 100 * max(bm_summary$improvement)),
  'Large pooled functional correlations mostly reflect between-drug/dose differences; consult within-condition results.',
  'Compare measured and predicted displacement on the identical held-out rows to isolate each link.',
  'Compare molecular combination models against both independence and identity/potency controls.',
  '', 'A: shared function', capture.output(print(a_metrics)), '', 'C: complete chain',
  capture.output(print(c_metrics)), '', 'D: interaction transfer', capture.output(print(dmetrics)))
writeLines(report, file.path(out, 'REPORT.txt'))
status('COMPLETE: A/B/C/D results and plots saved')
run_structure_diagnostics()
if (!report_only) run_drug_dose_adjusted()
