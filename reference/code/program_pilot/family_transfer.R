# Script-20 evaluation mode: script-15 discovery, global Ward40, excluded family
# and excluded recipients. Invoked by --family-transfer=MEK or =ERBB.
# The workflow is sequential and checkpoints each source and fitted model.
functional_folds <- 5L
functional_ridge_lambda <- .1
outcome_p_threshold <- .05
regression_gain_p_threshold <- .05
minimum_loading_iqr <- .05
rank <- 40L
out <- file.path(work_dir, 'results/program_pilot',
  paste0('20_multidrug_baseline_response_function_', family, '_holdout'))
data_dir <- file.path(out, 'data')
model_dir <- file.path(out, 'models')
figure_dir <- file.path(out, 'figures')
for (path in c(out, data_dir, model_dir, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(data_dir, name), sep = '\t')
status <- function(stage) {
  message(format(Sys.time()), ' | ', family, ' | ', stage)
  fwrite(data.table(stage = stage, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(out, 'run_status.tsv'), sep = '\t')
}
commit <- function(object, path) {
  temporary <- tempfile(tmpdir = dirname(path))
  qsave(object, temporary, preset = 'high')
  check(file.rename(temporary, path), paste('Cannot commit', path))
}
options(error = function() { status('FAILED: see run.log'); traceback(8); quit(save = 'no', status = 1L) })
status('preparing family exclusion')
panel_path <- file.path(work_dir, 'results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
annotation_path <- file.path(cache_dir, 'prism/secondary-screen-replicate-collapsed-treatment-info.csv')
expression_path <- file.path(work_dir, 'outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv')
panel <- qread(panel_path)
cell_map <- as.data.table(panel$cell_map)
all_drugs <- sort(unique(panel$metadata$drug))
annotation <- unique(fread(annotation_path)[, .(name, moa, target)])
annotation[, drug := panel_name(name)]
if (family == 'MEK') {
  explicit_family <- c('trametinib', 'trametinib (dmso_tf solvate)', 'cobimetinib', 'binimetinib')
  annotated <- annotation[grepl('(^|, )MEK inhibitor', moa) |
    grepl('(^|, )(MAP2K1|MAP2K2)(,|$)', target), drug]
} else {
  explicit_family <- c('afatinib', 'almonertinib (hydrochloride)', 'almonertinib (mesylate)',
    'erlotinib', 'gefitinib', 'lapatinib', 'neratinib', 'neratinib (maleate)',
    'osimertinib (mesylate)', 'simotinib', 'tucatinib', 'vandetanib')
  annotated <- annotation[grepl('(^|, )EGFR inhibitor', moa) |
    grepl('(^|, )(EGFR|ERBB2)(,|$)', target), drug]
}
excluded_drugs <- sort(intersect(all_drugs, unique(c(explicit_family, annotated))))
check(all(explicit_family %chin% excluded_drugs), 'Explicit family compound absent from panel')
write_table(merge(data.table(drug = excluded_drugs), annotation, by = 'drug', all.x = TRUE), 'excluded_family.tsv')
training_ids <- setdiff(cell_map$depmap_id, held_ids)
source_files <- c(panel_path, panel$input_paths, annotation_path, expression_path,
  file.path(code_dir, c('functional_programs.R', 'response_models.R', 'family_transfer.R', 'utilities.R')))
manifest <- data.table(path = source_files,
  sha256 = vapply(source_files, digest::digest, character(1), file = TRUE, algo = 'sha256'))
signature <- digest::digest(list(manifest, family, excluded_drugs, held_ids,
  p = outcome_p_threshold, gain_p = regression_gain_p_threshold, iqr = minimum_loading_iqr,
  ridge = functional_ridge_lambda, rank = rank), algo = 'sha256')
write_table(manifest, 'input_hashes.tsv')
commit(list(family = family, excluded_drugs = excluded_drugs, held_ids = held_ids,
  rank = rank, signature = signature, lambda_grid = lambda_grid), file.path(out, 'settings.qs'))
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))

# Recompute discovery from permitted rows; never use all-drug regression caches.
local_banks <- list()
audit <- list()
for (id in training_ids) {
  i <- match(id, cell_map$depmap_id)
  recipient <- cell_map$depmap_name[i]
  path <- file.path(model_dir, paste0('source_', id, '.qs'))
  if (file.exists(path)) {
    local <- qread(path)
    check(identical(local$signature, signature), paste('Stale source cache:', id))
  } else {
    status(paste('discovering', recipient, match(id, training_ids), '/', length(training_ids)))
    molecular <- qread(panel$input_paths[i])
    metadata <- as.data.table(molecular$metadata)[!drug %chin% excluded_drugs]
    x <- molecular$response[match(metadata$observation_id, rownames(molecular$response)), , drop = FALSE]
    check(identical(rownames(x), metadata$observation_id) && all(is.finite(x)), 'Invalid discovery rows')
    outcomes <- as.data.table(panel$sensitivity)[match(metadata$observation_id, observation_id)]
    check(!any(metadata$depmap_id %chin% held_ids) && !any(outcomes$drug %chin% excluded_drugs), 'Discovery leakage')
    tests <- rbindlist(lapply(c('PRISM', 'ALMANAC'), function(assay) {
      column <- if (assay == 'PRISM') 'prism_abundance_pct' else 'almanac_abundance_pct'
      test <- gene_outcome_tests(x, outcomes[[column]], metadata$drug)
      test[, assay := assay]
      test
    }))
    candidates <- sort(unique(tests[is.finite(p_value) & p_value <= outcome_p_threshold, gene]))
    check(length(candidates) >= 3L, paste('Fewer than three Pearson candidates:', recipient))
    support <- functional_support(x[, candidates, drop = FALSE], metadata, outcomes, sensitivity_sign = TRUE)
    admission <- regression_admission(support, id, recipient, regression_gain_p_threshold, minimum_loading_iqr)
    genes <- sort(unique(admission[selected_in_assay == TRUE, gene]))
    # Insufficient sources are explicit; no replacement genes or altered cutoff.
    local <- list(signature = signature, genes = genes, correlation = NULL,
      screen = tests, admission = admission, metadata = metadata,
      coefficients = support$coefficients[gene %chin% genes], observations = support$observations)
    if (length(genes) >= 3L) {
      support$loadings <- support$loadings[genes, , drop = FALSE]
      support$settings$clustering_representation <- 'mean_functional_displacement_correlation'
      null_energy <- mean(rowMeans(molecular$null[[id]]^2))
      energy <- rowMeans(x^2)
      condition_support <- energy / (energy + null_energy)
      displacement_correlation <- gene_correlation(x[, genes, drop = FALSE],
        group_weights(metadata$drug) * condition_support)
      local$correlation <- functional_nmf_input(support, displacement_correlation, FALSE)$correlation
      local$loadings <- support$loadings
    }
    commit(local, path)
    rm(molecular, x, support)
  }
  eligible <- length(local$genes) >= 3L
  audit[[id]] <- data.table(depmap_id = id, cell = recipient, genes = length(local$genes),
    conditions = nrow(local$metadata), drugs = uniqueN(local$metadata$drug),
    status = if (eligible) 'included' else 'fewer than three admitted genes: excluded')
  if (eligible) local_banks[[id]] <- local[c('genes', 'correlation')]
  write_table(rbindlist(audit), 'source_discovery_audit.tsv')
  message('Committed ', recipient, ': ', length(local$genes), ' admitted genes')
}

bank_path <- file.path(model_dir, 'program_bank.qs')
if (file.exists(bank_path)) {
  bank <- qread(bank_path)
  check(identical(bank$signature, signature), 'Stale global bank')
} else {
  status('aggregating observed gene pairs and cutting Ward tree at 40')
  genes <- sort(unique(unlist(lapply(local_banks, `[[`, 'genes'))))
  total <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  count <- matrix(0L, length(genes), length(genes))
  for (local in local_banks) {
    ix <- match(local$genes, genes)
    total[ix, ix] <- total[ix, ix] + local$correlation
    count[ix, ix] <- count[ix, ix] + 1L
  }
  relation <- total / pmax(count, 1L)
  check(max(abs(relation - t(relation))) < 1e-12 && all(is.finite(relation)), 'Invalid mean relationships')
  tree <- hclust(as.dist(sqrt(2 * (1 - pmin(pmax(relation, -1), 1)))), method = 'ward.D2')
  partition <- ordered_partition(tree, rank)
  programs <- sprintf('P%02d', seq_len(rank))
  size <- tabulate(partition, rank)
  weights <- sparseMatrix(i = seq_along(genes), j = partition, x = 1 / size[partition],
    dims = c(length(genes), rank), dimnames = list(genes, programs))
  bank <- list(genes = genes, partition = partition, tree = tree, weights = weights,
    signature = signature, training_ids = training_ids, excluded_ids = held_ids,
    excluded_drugs = excluded_drugs, missing_pair_fraction = mean(count[upper.tri(count)] == 0L))
  commit(bank, bank_path)
  # Raw unobserved pairs remain missing; zero is solely the Ward input convention.
  raw_relation <- relation
  raw_relation[count == 0L] <- NA_real_
  commit(list(genes = genes, mean_relationship = raw_relation, contributing_cells = count),
    file.path(data_dir, 'averaged_relationships.qs'))
  h <- draw_heatmap(relation[tree$order, tree$order], scale_rows = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
    hm_name = 'Mean relation', title = paste(family, 'excluded | Ward 40'))
  pdf(file.path(figure_dir, 'training_gene_relationships.pdf'), width = 7, height = 6.5)
  ComplexHeatmap::draw(h)
  dev.off()
  rm(total, count, relation, raw_relation, h)
}
rm(local_banks)
invisible(gc())
programs <- colnames(bank$weights)
write_table(data.table(gene = bank$genes, program = programs[bank$partition]), 'gene_assignments.tsv')
write_table(data.table(program = programs, genes = as.integer(colSums(bank$weights > 0))), 'program_sizes.tsv')
# Membership overlap reports how excluding data changes the exploratory bank.
reference <- fread(file.path(work_dir,
  'results/program_pilot/16_global_mean_relationships_fixed_ward_programs/data/global_gene_assignments.tsv'))
overlap <- merge(data.table(gene = bank$genes, program = programs[bank$partition]), reference, by = 'gene')
write_table(overlap[, .(shared_genes = .N), by = .(program, global_program)], 'reference40_overlap.tsv')

score_path <- file.path(data_dir, 'prepared_scores.qs')
if (file.exists(score_path)) {
  prepared <- qread(score_path)
  check(identical(prepared$signature, signature), 'Stale score cache')
} else {
  status('scoring the frozen forty programs in every molecular sample')
  scores <- metadata_parts <- list()
  for (i in seq_len(nrow(cell_map))) {
    molecular <- qread(panel$input_paths[i])
    metadata_parts[[i]] <- as.data.table(molecular$metadata)[, .(observation_id, depmap_id, drug, dose)]
    scores[[i]] <- as.matrix(molecular$response[, bank$genes, drop = FALSE] %*% bank$weights)
    check(identical(rownames(scores[[i]]), metadata_parts[[i]]$observation_id), 'Score order mismatch')
  }
  header <- names(fread(expression_path, nrows = 0L))
  symbols <- sub(' \\(ENSG.*', '', header[-1L])
  selected <- which(symbols %chin% unique(c(bank$genes, unlist(target_genes))))
  expression_table <- fread(expression_path, select = c(1L, selected + 1L))
  expression <- as.matrix(expression_table[, -1L])
  dimnames(expression) <- list(expression_table[[1L]], symbols[selected])
  ids <- intersect(cell_map$depmap_id, rownames(expression))
  prepared <- list(signature = signature, measured = do.call(rbind, scores),
    metadata = rbindlist(metadata_parts),
    baseline = as.matrix(expression[ids, bank$genes, drop = FALSE] %*% bank$weights),
    targets = expression[ids, unique(unlist(target_genes)), drop = FALSE])
  commit(prepared, score_path)
  rm(scores, metadata_parts, expression, expression_table, molecular)
}
measured <- prepared$measured
md <- prepared$metadata
baseline <- prepared$baseline
targets <- prepared$targets
check(all(is.finite(measured)) && all(is.finite(baseline)), 'Nonfinite program scores')
rm(prepared)
# Only development outcomes are read into a fitting table. Evaluation outcomes
# are kept in a separate table and do not enter any centering or model selection.
labels <- as.data.table(panel$sensitivity)[is.finite(prism_abundance_pct) & prism_abundance_pct > 0 &
  depmap_id %chin% rownames(baseline)]
labels[, observed := log2(prism_abundance_pct / 100)]
labels[, condition := paste(drug, dose, sep = '__')]
train <- labels[!depmap_id %chin% held_ids & !drug %chin% excluded_drugs]
test <- labels[depmap_id %chin% held_ids & drug %chin% excluded_drugs]
check(nrow(train) > 0L && nrow(test) > 0L, 'Empty functional split')
check(!any(train$drug %chin% excluded_drugs) && !any(train$depmap_id %chin% held_ids), 'Functional leakage')
write_table(rbindlist(list(train[, .(split = 'training', rows = .N, cells = uniqueN(depmap_id)), by = drug],
  test[, .(split = 'test', rows = .N, cells = uniqueN(depmap_id)), by = drug])), 'functional_split.tsv')
write_table(test[, .(observation_id, depmap_id, drug, dose, prism_status, prism_lower_dose,
  prism_upper_dose, prism_screen)], 'test_dose_provenance.tsv')
inner_ids <- sort(unique(train$depmap_id))
set.seed(20261010L)
inner_map <- setNames(sample(rep(1:5, length.out = length(inner_ids))), inner_ids)
# Same grouped ridge code as script 20. Validation is conditional on the frozen
# outer-training bank; outer test cells and family never enter parameter selection.
bmu <- colMeans(baseline[inner_ids, , drop = FALSE])
bsd <- apply(baseline[inner_ids, , drop = FALSE], 2L, sd)
check(all(bsd > 0), 'Constant baseline program')
b <- sweep(sweep(baseline, 2L, bmu, '-'), 2L, bsd, '/')
dtrain <- measured[match(train$observation_id, md$observation_id), , drop = FALSE]
dtest <- measured[match(test$observation_id, md$observation_id), , drop = FALSE]
btrain <- b[train$depmap_id, , drop = FALSE]
btest <- b[test$depmap_id, , drop = FALSE]
# Treatment-average molecular responses use training recipients only. Their
# functional outcomes remain excluded, including for these query compounds.
means <- vapply(seq_len(nrow(test)), function(i) {
  ii <- which(md$drug == test$drug[i] & md$dose == test$dose[i] & !md$depmap_id %chin% held_ids)
  check(length(ii) > 0L, 'No training-cell molecular reference')
  colMeans(measured[ii, , drop = FALSE])
}, numeric(rank))
dmean <- t(means)
colnames(dmean) <- programs
swapped <- seq_len(nrow(test))
set.seed(20261010L)
for (condition in unique(test$condition)) {
  ii <- which(test$condition == condition)
  swapped[ii] <- sample(ii)
}
designs <- list(measured_D = dtrain, measured_D_conditioned = cbind(dtrain, dtrain * btrain),
  response_magnitude = cbind(magnitude = sqrt(rowMeans(dtrain^2)),
    btrain * sqrt(rowMeans(dtrain^2))))
decoder_path <- file.path(model_dir, 'functional_decoder.qs')
models <- list()
prediction_parts <- list()
for (name in names(designs)) {
  status(paste('fitting shared functional model:', name))
  path <- file.path(model_dir, paste0('decoder_', name, '.qs'))
  if (file.exists(path)) {
    saved <- qread(path)
    check(identical(saved$signature, signature), 'Stale decoder')
    f <- saved$model
  } else {
    f <- fit_ridge(designs[[name]], train$observed, train$depmap_id,
      observation_weights(train), intercept = FALSE)
    commit(list(signature = signature, model = f), path)
  }
  models[[name]] <- f
  inputs <- if (name == 'response_magnitude') list(measured = dtest) else
    list(measured = dtest, training_mean = dmean, swapped_recipient = dtest[swapped, , drop = FALSE])
  for (input in names(inputs)) {
    dd <- inputs[[input]]
    xx <- switch(name, measured_D = dd, measured_D_conditioned = cbind(dd, dd * btest),
      response_magnitude = cbind(magnitude = sqrt(rowMeans(dd^2)), btest * sqrt(rowMeans(dd^2))))
    pred <- test[, .(observation_id, depmap_id, drug, dose, observed)]
    pred[, `:=`(model = name, response_input = input, predicted = predict_ridge(f, xx))]
    prediction_parts[[paste(name, input)]] <- pred
  }
  write_table(rbindlist(prediction_parts), 'functional_predictions.tsv.gz')
}
commit(list(signature = signature, models = models, baseline_mean = bmu, baseline_sd = bsd,
  training_ids = unique(train$depmap_id), training_drugs = unique(train$drug)), decoder_path)

# Minimal reporting: native dose curves and within-dose associations, no local
# discovery plots or additional cluster-count search. Preserve response magnitude.
pred <- rbindlist(prediction_parts)
pred[, `:=`(observed_pct = 100 * 2^observed, predicted_pct = 100 * 2^predicted)]
pred[, cell := cell_map$depmap_name[match(depmap_id, cell_map$depmap_id)]]
write_table(pred, 'functional_predictions.tsv.gz')
per_condition <- pred[, .(cells = .N, pearson = cor_safe(observed, predicted),
  observed_min_pct = min(observed_pct), observed_max_pct = max(observed_pct),
  predicted_min_pct = min(predicted_pct), predicted_max_pct = max(predicted_pct),
  median_absolute_error_percentage_points = median(abs(predicted_pct - observed_pct))),
  by = .(model, response_input, drug, dose)]
write_table(per_condition, 'functional_within_dose.tsv')
centered <- copy(pred)
centered[, `:=`(observed_within = observed - mean(observed), predicted_within = predicted - mean(predicted)),
  by = .(model, response_input, drug, dose)]
write_table(centered[, .(cells = uniqueN(depmap_id), observations = .N,
  within_dose_pearson = cor_safe(observed_within, predicted_within),
  observed_cell_spread = sd(observed_within), predicted_cell_spread = sd(predicted_within)),
  by = .(model, response_input, drug)], 'functional_by_drug.tsv')
curves <- pred[model == 'measured_D_conditioned' & response_input %chin% c('measured', 'training_mean')]
observed <- unique(curves[, .(drug, dose, cell, observed_pct)])
plot_rows <- rbind(curves[, .(drug, dose, cell, abundance = predicted_pct,
  series = ifelse(response_input == 'measured', 'Individual molecular response', 'Average molecular response'))],
  observed[, .(drug, dose, cell, abundance = observed_pct, series = 'Observed PRISM')])
for (drug_name in unique(plot_rows$drug)) {
  p <- ggplot(plot_rows[drug == drug_name], aes(dose, abundance, colour = series, group = series)) +
    geom_line(linewidth = .45) + geom_point(size = 1.3) + scale_x_log10(breaks = c(.05, .5, 5)) +
    facet_wrap(~cell, ncol = 3) + theme_classic(base_size = 9) +
    theme(legend.position = 'bottom') + labs(x = 'Dose (uM)', y = 'Abundance (% vehicle)', colour = NULL,
      title = paste(drug_name, '| family and recipients excluded'),
      subtitle = 'PRISM bracket-interpolated to molecular doses; same fixed forty programs')
  ggsave(file.path(figure_dir, paste0(gsub('[^a-z0-9]+', '_', drug_name), '_dose_curves.pdf')),
    p, width = 9, height = 7)
}
source(file.path(code_dir, 'family_transfer_plots.R'))
plot_family_transfer(pred, figure_dir, data_dir)
status('COMPLETE: measured-response family transfer')
message('Saved results: ', out)
