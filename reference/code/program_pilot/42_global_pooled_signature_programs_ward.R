#!/usr/bin/env Rscript
# Outcome-informed descriptive pilot. Gene admission uses retrospective PRISM
# and ALMANAC correlations plus held-out regression gain and loading variation.
# Ward.D2 clusters the 50/50 mean of functional-loading and displacement correlations;
# full-data gene screening makes this retrospective, not independent validation.
# Run from /work. The molecular-only 04 script/results remain a separate experiment.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(cluster)
  library(patchwork)
  library(RcppML)
  library(glmnet)
})
source("/work/src/program_pilot/config.R")
arguments <- commandArgs(trailingOnly = TRUE)
check(length(arguments) == 0L || identical(arguments, "--render-only"),
  "Supported invocation: Rscript 42_global_pooled_signature_programs_ward.R [--render-only]")
render_only <- identical(arguments, "--render-only")
# One Ward.D2 hierarchy supplies every candidate cluster count.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
RhpcBLASctl::blas_set_num_threads(8L)
RhpcBLASctl::omp_set_num_threads(1L)
local_resamples <- 0L # One Ward.D2 tree, without bootstraps or repeated fits.
local_minimum_rank <- 2L
local_rank_grid <- 2:80 # Expanded search requested for this Ward.D2 version.
global_rank_override <- 60L # User-selected global cut; reuse the saved hierarchy.
local_rank_selection <- "silhouette_elbow" # Stop at the steep-to-flatter bend; every seed still requires >=12 genes.
functional_folds <- 5L # Hold whole cell lines together; every dose remains a separate row.
functional_ridge_lambda <- 0.1 # Fixed mild ridge penalty on training-standardized predictor and outcome.
functional_input_transform <- "50/50 mean gene-gene Pearson correlation of signed functional loadings and original displacements; Ward.D2 on sqrt(2*(1-C))"
regression_gain_p_threshold <- 0.05 # One-sided nominal test: positive held-out gain across whole cell lines.
minimum_loading_iqr <- 0.05 # Signed gain spread, in units of training outcome variance; robust to isolated outliers.
# This separately requested experiment keeps its outputs and gene rule local.
# Test each drug/dose row equally; pool doses without dose adjustment or averaging.
outcome_p_threshold <- 0.05 # Two-sided nominal threshold; BH values remain diagnostic only.
outcome_minimum_abs_correlation <- 0 # Earlier lenient rule has no effect-size floor.
local_gene_limit <- Inf # Retain every significant gene in the existing expressed coding universe.
output_dir <- file.path(work_dir, "results/program_pilot/42_global_pooled_signature_programs_ward")
data_dir <- file.path(output_dir, "data")
figure_dir <- file.path(output_dir, "figures")
model_dir <- file.path(output_dir, "models")
for (path in c(output_dir, data_dir, figure_dir, model_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
set.seed(seed)
options(warn = 1)
theme_set(theme_classic(base_size = 10))
message("Outcome-selected displacement pilot; seed ", seed, "; results: ", output_dir)

source(file.path(code_dir, "functional_programs.R"))
# All eligible cells come from the explicit panel preparation, including assay-only lines.
prep_path <- file.path(work_dir, "results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs")
check(file.exists(prep_path), "Run prepare_functional_panel.R before the full panel")
panel <- qread(prep_path)
cell_map <- panel$cell_map
metadata <- copy(panel$metadata)
sensitivity <- copy(panel$sensitivity)
genes <- panel$genes
input_files <- c(gene_audit_path, panel$input_paths)
responses <- null_list <- baseline_list <- checks <- list()
for (id in cell_map$depmap_id) {
  molecular <- qread(panel$input_paths[match(id, cell_map$depmap_id)])
  rows <- which(metadata$depmap_id == id)
  response <- molecular$response[match(metadata$observation_id[rows], rownames(molecular$response)), , drop = FALSE]
  check(identical(colnames(response), genes) && all(is.finite(response)), "Invalid prepared panel response")
  responses[[id]] <- response
  null_list[[id]] <- molecular$null[[id]]
  baseline_list[[id]] <- molecular$baseline[[id]]
  checks[[id]] <- data.table(depmap_id=id, profiles=nrow(response), null_profiles=nrow(null_list[[id]]), finite=TRUE)
}
prepared <- list(response=do.call(rbind,responses), null=null_list, baseline=baseline_list, checks=rbindlist(checks))
metadata <- metadata[match(rownames(prepared$response), observation_id)]
rm(responses, molecular, response, null_list, baseline_list)
invisible(gc())
fwrite(prepared$checks, file.path(data_dir, "normalization_checks.tsv"), sep="\t")
fwrite(cell_map, file.path(data_dir, "panel_cell_mapping.tsv"), sep="\t")
fwrite(metadata, file.path(data_dir, "molecular_metadata.tsv.gz"), sep="\t")
fwrite(sensitivity, file.path(data_dir, "sensitivity_alignment.tsv"), sep="\t")
almanac_path <- file.path(work_dir, "almanac_conditionals/raw/ComboDrugGrowth_Nov2017.csv")
drug_mapping_path <- file.path(work_dir, "almanac_conditionals/repaired/results/drug_mapping.tsv")
prism_info_path <- file.path(cache_dir, "prism/secondary-screen-replicate-collapsed-treatment-info.csv")
prism_matrix_path <- file.path(cache_dir, "prism/secondary-screen-replicate-collapsed-logfold-change.csv")
function_mapping_path <- file.path(work_dir, "tahoe_prism_v2/artifacts/01_inputs/function_outcomes.tsv.gz")
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm = TRUE) / 25) * 25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm = TRUE) / 25) * 25)
message("Sensitivity bar maxima (% vehicle): ALMANAC ", almanac_bar_max, "; PRISM ", prism_bar_max)

setDTthreads(1L)

# Preserve script-15 molecular/control support separately for each source cell.
for (id in cell_map$depmap_id) {
  rows <- which(metadata$depmap_id == id)
  energy <- rowMeans(prepared$response[rows, , drop = FALSE]^2)
  null_energy <- mean(rowMeans(prepared$null[[id]]^2))
  metadata[rows, `:=`(condition_support = energy / (energy + null_energy),
    global_displacement = rowMeans(prepared$response[rows, , drop = FALSE]),
    raw_response_rms = sqrt(energy), control_rms = sqrt(null_energy),
    response_control_rms_ratio = sqrt(energy / null_energy))]
}
check(!anyDuplicated(metadata$observation_id) && all(is.finite(metadata$condition_support)),
  "Invalid molecular input identities or control support")
fwrite(metadata, file.path(data_dir, "molecular_metadata.tsv.gz"), sep = "\t")

# One bank pools every original cell-line/drug/dose signature. Drug-level
# coverage controls only which per-drug views are exported, never admission to
# the pooled input. Missing assay outcomes remain missing, with no imputation.
coverage <- sensitivity[, .(molecular_conditions = .N, cell_lines = uniqueN(depmap_id),
  prism_cells = uniqueN(depmap_id[is.finite(prism_abundance_pct)]),
  almanac_cells = uniqueN(depmap_id[is.finite(almanac_abundance_pct)])), by = drug]
coverage[, eligible := prism_cells >= functional_folds | almanac_cells >= functional_folds]
fwrite(coverage, file.path(data_dir, "drug_assay_coverage.tsv"), sep = "\t")
plot_drugs <- coverage[eligible == TRUE, drug]
drugs <- sort(unique(metadata$drug))
# After the global fit, show focal drugs first and the remaining measured drugs.
plot_drugs <- c(intersect(local_drugs, plot_drugs), sort(setdiff(plot_drugs, local_drugs)))
message('GLOBAL INPUT: ', nrow(metadata), ' signatures; ', length(drugs), ' drugs; ', nrow(cell_map), ' cell lines')
gene_selection <- outcome_test_tables <- outcome_test_summary <- discoveries <- list()
run_rows <- list()
qsave(list(parent_script = '15_mean_functional_displacement_correlation_programs_ward.R',
  grouping = 'all cell-line/drug/dose signatures pooled', regression_units = 'whole cell lines; doses remain separate',
  ranks = local_rank_grid, rank_selection = local_rank_selection, rank_override = global_rank_override,
  minimum_genes = local_minimum_genes,
  outcome_p = outcome_p_threshold, gain_p = regression_gain_p_threshold, loading_iqr = minimum_loading_iqr,
  ridge_lambda = functional_ridge_lambda, eligible_drugs = drugs,
  scientific_status = 'retrospective discovery across all measured cell lines; no independent validation'),
  file.path(output_dir, 'settings.qs'))

for (drug_name in 'pooled') {
  started <- Sys.time()
  rows <- seq_len(nrow(metadata))
  md <- metadata[rows]
  scope <- '.'
  scope_data_dir <- file.path(output_dir, scope, 'data')
  scope_model_dir <- file.path(output_dir, scope, 'models')
  scope_figure_dir <- file.path(output_dir, scope, 'figures')
  for (path in c(scope_data_dir, scope_model_dir, scope_figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  message('START ', scope, ' | ', uniqueN(md$depmap_id), ' cells | ', nrow(md), ' cell/dose columns')
  x <- prepared$response[rows, , drop = FALSE]
  outcomes <- sensitivity[match(md$observation_id, observation_id)]
  check(identical(rownames(x), md$observation_id) && identical(outcomes$observation_id, md$observation_id),
    'Drug input/outcome alignment failed')
  # Each cell contributes its existing vehicle pseudo-responses. These pooled
  # controls supply the same coherence diagnostics and score construction.
  null <- do.call(rbind, prepared$null[unique(md$depmap_id)])
  tests <- rbindlist(lapply(c('PRISM', 'ALMANAC'), function(assay_name) {
    column <- if (assay_name == 'PRISM') 'prism_abundance_pct' else 'almanac_abundance_pct'
    tested <- gene_outcome_tests(x, outcomes[[column]], md$depmap_id)
    setnames(tested, 'matched_drugs', 'matched_cell_lines')
    tested[, `:=`(drug = drug_name, assay = assay_name)]
    tested
  }))
  tests[, q_value := p.adjust(p_value, method = 'BH'), by = assay]
  tests[, significant := !is.na(p_value) & p_value <= outcome_p_threshold &
    abs(correlation) >= outcome_minimum_abs_correlation]
  outcome_test_tables[[drug_name]] <- tests
  outcome_test_summary[[drug_name]] <- tests[, .(candidate_genes = .N,
    tested_genes = sum(test_status == 'tested'), significant_genes = sum(significant),
    matched_conditions = matched_conditions[1L], matched_cell_lines = matched_cell_lines[1L]), by = .(drug, assay)]
  candidates <- sort(unique(tests[significant == TRUE, gene]))
  selection <- data.table(drug = drug_name, gene = genes, outcome_eligible = genes %chin% candidates,
    response_variance = matrixStats::colVars(x), control_variance = matrixStats::colVars(null))
  selection[, `:=`(selected = FALSE, exclusion_reason = fifelse(outcome_eligible,
    'awaiting regression admission', 'fails nominal outcome correlation in both assays'))]
  fwrite(tests, file.path(scope_data_dir, 'gene_outcome_correlations.tsv.gz'), sep = '\t')
  fwrite(outcome_test_summary[[drug_name]], file.path(scope_data_dir, 'gene_outcome_test_summary.tsv'), sep = '\t')
  if (length(candidates) < 3L) {
    selection[outcome_eligible == TRUE, exclusion_reason := 'fewer than three candidate genes; no program fit']
    gene_selection[[drug_name]] <- selection
    fwrite(selection, file.path(scope_data_dir, 'gene_selection.tsv.gz'), sep = '\t')
    run_rows[[drug_name]] <- data.table(drug = drug_name, status = 'fewer than three outcome candidates',
      candidates = length(candidates), selected_genes = 0L, selected_rank = NA_integer_, seconds = as.numeric(difftime(Sys.time(), started, units='secs')))
    fwrite(rbindlist(run_rows), file.path(output_dir, 'global_status.tsv'), sep = '\t')
    next
  }
  # One regression per gene and assay across pooled signatures. Keep all drugs
  # and doses from each held-out cell together; earlier local fits are not reused.
  signature <- digest::digest(list(x = x[, candidates, drop = FALSE], metadata = md,
    outcomes = outcomes, candidates = candidates, folds = functional_folds,
    lambda = functional_ridge_lambda, unit = 'depmap_id', sign = TRUE), algo = 'sha256')
  regression_path <- file.path(scope_model_dir, 'functional_regression.qs')
  if (file.exists(regression_path)) {
    regression <- qread(regression_path)
    check(identical(regression$input_signature, signature), paste('Regression inputs changed for', scope))
    support <- regression$support
  } else {
    message(scope, ' | fitting ', length(candidates), ' gene candidates; whole-cell folds')
    support <- functional_support(x[, candidates, drop = FALSE], md, outcomes,
      allow_zero = TRUE, sensitivity_sign = TRUE, unit_column = 'depmap_id')
    support$channels <- NULL
    qsave(list(input_signature = signature, support = support), regression_path, preset = 'high')
  }
  fold_cells <- merge(support$folds, md[, .(observation_id, depmap_id)], by = 'observation_id')
  check(all(fold_cells[, .(n = uniqueN(fold)), by = .(assay, depmap_id)]$n == 1L),
    'Doses from one cell were split across regression folds')
  regression_tests <- regression_admission(support, drug_name, drug_name,
    regression_gain_p_threshold, minimum_loading_iqr, unit_column = 'depmap_id')
  regression_tests[, c('depmap_id', 'recipient') := NULL]
  regression_tests[, drug := drug_name]
  admission <- regression_tests[, .(regression_significant = any(regression_significant),
    loading_variation_eligible = any(loading_variation_eligible),
    functional_eligible = any(selected_in_assay)), by = gene]
  chosen <- candidates[candidates %chin% admission[functional_eligible == TRUE, gene]]
  selection <- merge(selection, admission, by = 'gene', all.x = TRUE, sort = FALSE)
  selection[, selected := gene %chin% chosen]
  selection[, exclusion_reason := fcase(!outcome_eligible, 'fails nominal outcome correlation in both assays',
    !regression_significant, 'no significant positive held-out cell-mean gain',
    !loading_variation_eligible, 'signed loading IQR below threshold',
    !functional_eligible, 'gain and variation fail to qualify in the same assay',
    selected, 'selected', default = 'invalid selection state')]
  gene_selection[[drug_name]] <- selection
  fwrite(selection, file.path(scope_data_dir, 'gene_selection.tsv.gz'), sep = '\t')
  fwrite(regression_tests, file.path(scope_data_dir, 'gene_regression_admission.tsv.gz'), sep = '\t')
  gain_plot <- ggplot(regression_tests, aes(-log10(pmax(p_value, .Machine$double.xmin)), fill = selected_in_assay)) +
    geom_histogram(bins = 35L, position = 'identity', alpha = .6, na.rm = TRUE) +
    geom_vline(xintercept = -log10(regression_gain_p_threshold), lty = 2) + facet_wrap(~assay, ncol = 1L) +
    scale_fill_manual(values = c(`FALSE` = 'grey70', `TRUE` = '#27818A')) +
    labs(x = '-log10(held-out cell-mean gain nominal p)', y = 'Candidate genes', fill = 'Pass both gates')
  variation_plot <- ggplot(regression_tests, aes(loading_iqr, fill = selected_in_assay)) +
    geom_histogram(bins = 35L, position = 'identity', alpha = .6) +
    geom_vline(xintercept = minimum_loading_iqr, lty = 2) + facet_wrap(~assay, ncol = 1L) +
    scale_fill_manual(values = c(`FALSE` = 'grey70', `TRUE` = '#27818A')) +
    labs(x = 'Signed-loading IQR (normalized gain units)', y = 'Candidate genes', fill = 'Pass both gates')
  admission_plot <- gain_plot + variation_plot + plot_layout(guides = 'collect') +
    plot_annotation(title = paste(drug_name, '| regression gene admission'),
      subtitle = sprintf('Retain %d / %d genes; same-assay gain p <=%g and signed-loading IQR >=%g',
        length(chosen), length(candidates), regression_gain_p_threshold, minimum_loading_iqr))
  for (extension in c('pdf', 'png')) ggsave(file.path(scope_figure_dir, paste0('regression_gene_filter.', extension)),
    admission_plot, width = 10, height = if (uniqueN(regression_tests$assay) == 2L) 6.2 else 4.2,
    device = if (extension == 'pdf') cairo_pdf else 'png', dpi = 160, bg = 'white')
  if (length(chosen) < 3L) {
    run_rows[[drug_name]] <- data.table(drug = drug_name, status = 'fewer than three regression-admitted genes',
      candidates = length(candidates), selected_genes = length(chosen), selected_rank = NA_integer_,
      seconds = as.numeric(difftime(Sys.time(), started, units='secs')))
    fwrite(rbindlist(run_rows), file.path(output_dir, 'global_status.tsv'), sep = '\t')
    message(scope, ' | no supported program fit; admission diagnostics saved')
    next
  }
  support$loadings <- support$loadings[chosen, , drop = FALSE]
  support$raw_prediction_gain <- support$raw_prediction_gain[chosen, , drop = FALSE]
  support$coefficients <- support$coefficients[gene %chin% chosen]
  support$channels <- cbind(pmax(support$loadings, 0), pmax(-support$loadings, 0))
  colnames(support$channels) <- c(paste0('positive:', colnames(support$loadings)),
    paste0('negative:', colnames(support$loadings)))
  support$settings$clustering_representation <- 'mean_functional_displacement_correlation'
  support$settings$clustering_algorithm <- 'ward.D2'
  support$settings$rank_selection_rule <- local_rank_selection
  support$settings$gene_admission <- list(p_threshold = regression_gain_p_threshold,
    minimum_loading_iqr = minimum_loading_iqr, test = 'one-sided held-out cell-mean gain t test; same-assay IQR')
  bank_signature <- digest::digest(list(signature, chosen, regression_tests, local_rank_grid,
    local_rank_selection, local_minimum_genes, functional_input_transform), algo = 'sha256')
  bank_path <- file.path(scope_model_dir, 'program_bank.qs')
  if (file.exists(bank_path)) {
    bank <- qread(bank_path)
    check(identical(bank$gene_selection_signature, bank_signature), paste('Saved bank inputs changed for', scope))
    bank <- restore_functional_bank(bank, prepared$response)
  } else {
    check(!render_only, paste('No saved bank for', scope))
    message(scope, ' | Ward.D2 clustering ', length(chosen), ' genes')
    bank <- discover_local(x, null, md, chosen, 'depmap_id', scope, 'pooled',
      support = support, sensitivity_sign = TRUE)
    bank$unit_column <- 'depmap_id'
    bank$gene_selection_signature <- bank_signature
    bank$gene_selection_rule <- 'Script-15 thresholds; whole cell lines replace whole drugs as regression units'
    save_functional_bank(bank, bank_path)
  }
  if (!identical(bank$diagnostics$selected_rank, global_rank_override)) {
    check(!render_only, 'Run without --render-only to apply the requested global rank')
    previous_path <- file.path(model_dir, sprintf('program_bank_rank%02d.qs', bank$diagnostics$selected_rank))
    if (!file.exists(previous_path)) check(file.copy(bank_path, previous_path), 'Could not preserve previous selected bank')
    message('Selecting global rank ', global_rank_override, ' from saved Ward cuts; no regressions or hierarchy refitted')
    saved_partitions <- bank$all_partitions
    saved_signature <- bank$gene_selection_signature
    bank <- discover_local(x, null, md, chosen, 'depmap_id', scope, 'pooled',
      support = bank$functional_support, sensitivity_sign = TRUE, previous = bank,
      rank_override = global_rank_override)
    check(identical(bank$all_partitions, saved_partitions), 'Rank selection altered the saved Ward partitions')
    bank$unit_column <- 'depmap_id'
    bank$gene_selection_signature <- saved_signature
    bank$gene_selection_rule <- 'Global script-41 screening; user-selected rank 60 of the saved Ward hierarchy'
    save_functional_bank(bank, bank_path)
  }
  check(uniqueN(bank$metadata$drug) == length(drugs) && identical(bank$metadata$observation_id, md$observation_id),
    'Global bank lost drugs or input columns')
  # The installed ComplexHeatmap matrix resizer has an out-of-bounds bug for
  # this wide input. Standard device rasterization preserves plotted values.
  bank$raster_resize_matrix <- FALSE
  render_diagnostics(bank)
  fwrite(data.table(stage = 'global clustering saved; rendering', genes = length(chosen),
    rank = bank$diagnostics$selected_rank, signatures = nrow(md)),
    file.path(output_dir, 'run_status.tsv'), sep = '\t')
  # Each local object is only a view of the frozen global bank, never a refit.
  view_status <- list()
  for (view_drug in plot_drugs) {
    idx <- which(bank$metadata$drug == view_drug)
    view <- bank
    view$scope <- paste0('drug_', gsub('[^A-Za-z0-9]+', '_', view_drug))
    view$orientation <- 'drug'
    view$metadata <- bank$metadata[idx]
    view$scores <- if (is.null(bank$scores)) NULL else bank$scores[idx, , drop = FALSE]
    cols <- which(bank$functional_support$observations$observation_id %chin% view$metadata$observation_id)
    view$functional_support$observations <- bank$functional_support$observations[cols]
    view$functional_support$loadings <- bank$functional_support$loadings[, cols, drop = FALSE]
    view$functional_support$raw_prediction_gain <- bank$functional_support$raw_prediction_gain[, cols, drop = FALSE]
    view$functional_support$folds <- bank$functional_support$folds[observation_id %chin% view$metadata$observation_id]
    # This small reference suppresses repeated export of the global quadratic
    # correlation input while preserving globally fixed program colors.
    if (!is.null(bank$weights)) view$global_structure <- list(path = bank_path,
      program_colors = setNames(hcl.colors(ncol(bank$weights), 'Dark 3'), colnames(bank$weights)))
    view_data <- file.path(output_dir, view$scope, 'data')
    dir.create(view_data, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, view$scope, 'models'), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, view$scope, 'figures'), recursive = TRUE, showWarnings = FALSE)
    if (is.null(view$weights)) render_unresolved(view) else {
      render_program_score_clustered_maps(view, prepared$response)
      render_local(view)
    }
    render_functional_loadings(view)
    check(identical(view$genes, bank$genes) && identical(view$partition, bank$partition) &&
      identical(view$weights, bank$weights), 'Per-drug view changed the global programs')
    writeLines(c(paste('Fixed global bank:', bank_path),
      'No gene screen, regression, rank selection or clustering of genes was repeated for this drug.'),
      file.path(output_dir, view$scope, 'README.txt'))
    view_status[[view_drug]] <- data.table(drug = view_drug, status = 'plots complete',
      genes = length(bank$genes), rank = bank$diagnostics$selected_rank, signatures = length(idx))
    fwrite(rbindlist(view_status), file.path(output_dir, 'per_drug_plot_status.tsv'), sep = '\t')
    message('FIXED GLOBAL PROGRAM FIGURES READY: ', file.path(output_dir, view$scope, 'figures'))
    rm(view)
  }
  if (is.null(bank$weights)) {
    overview_order <- order(outcomes$prism_abundance_pct, outcomes$almanac_abundance_pct, na.last = TRUE)
    render_unresolved(bank, column_order = overview_order)
  } else {
    # Sorting the full overview by sensitivity within program avoids a quadratic
    # 42,679-column distance. The per-drug views retain all original clustering.
    render_local(bank, detail_pages = FALSE, cluster_within_program = FALSE)
  }
  render_functional_loadings(bank)
  writeLines(c(readLines(file.path(output_dir, 'README.txt')),
    'ONE GLOBAL FIT: gene screen, gene regression and Ward.D2 all use the pooled signatures.',
    'All 42,679 original molecular signatures retained; no drug/cell centering or dose averaging.',
    'Five whole-cell folds keep all drugs and doses from a held-out cell together.',
    'Global overview sorts sensitivity within program; per-drug maps retain original column clustering.',
    'Per-drug figures use the same global gene partition, weights and program-score scaling.',
    'Outcome-informed retrospective discovery across all cells; not independent validation.'),
    file.path(output_dir, 'README.txt'))
  message('GLOBAL FIGURES READY: ', figure_dir)
  discoveries[[drug_name]] <- bank[c('scope', 'ranks', 'diagnostics', 'programs')]
  run_rows[[drug_name]] <- data.table(drug = drug_name, status = 'plots complete', candidates = length(candidates),
    selected_genes = length(chosen), selected_rank = bank$diagnostics$selected_rank,
    seconds = as.numeric(difftime(Sys.time(), started, units='secs')))
  fwrite(rbindlist(run_rows), file.path(output_dir, 'global_status.tsv'), sep = '\t')
  message('FIGURES READY: ', scope_figure_dir)
  rm(support, bank, x, null)
  invisible(gc())
}
fwrite(rbindlist(gene_selection), file.path(data_dir, 'gene_selection.tsv.gz'), sep = '\t')
fwrite(rbindlist(outcome_test_summary), file.path(data_dir, 'gene_outcome_test_summary.tsv'), sep = '\t')
if (length(discoveries)) {
  fwrite(rbindlist(lapply(discoveries, `[[`, 'ranks')), file.path(data_dir, 'rank_diagnostics.tsv'), sep = '\t')
  fwrite(rbindlist(lapply(discoveries, `[[`, 'diagnostics')), file.path(data_dir, 'discovery_diagnostics.tsv'), sep = '\t')
  fwrite(rbindlist(lapply(discoveries, `[[`, 'programs'), fill = TRUE), file.path(data_dir, 'program_diagnostics.tsv'), sep = '\t')
}
figures <- list.files(output_dir, pattern = '[.]pdf$', recursive = TRUE, full.names = TRUE)
fwrite(data.table(path = figures), file.path(output_dir, 'artifact_index.tsv'), sep = '\t')
write_manifest('42_global_pooled_signature_programs_ward', c(input_files, prep_path),
  list(parent = 'script 15', grouping = 'all cell-line/drug/dose signatures pooled', regression_units = 'cell lines', drugs = drugs,
    functional_folds = functional_folds, ridge_lambda = functional_ridge_lambda,
    rank_grid = local_rank_grid, rank_selection = local_rank_selection, rank_override = global_rank_override,
    outcome_p = outcome_p_threshold, gain_p = regression_gain_p_threshold, loading_iqr = minimum_loading_iqr))
writeLines(capture.output(sessionInfo()), file.path(output_dir, 'sessionInfo.txt'))
message('COMPLETE: one global bank across ', length(drugs), ' drugs; results ', output_dir)
fwrite(data.table(stage = 'complete', drugs = length(drugs), plotted_drugs = length(plot_drugs)),
  file.path(output_dir, 'run_status.tsv'), sep = '\t')
