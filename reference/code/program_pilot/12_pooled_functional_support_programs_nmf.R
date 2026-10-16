#!/usr/bin/env Rscript
# Shared genes across cell lines; columns are measured cell/drug/dose/assay
# signed functional loadings and gene admission from completed stage 13.
# Retrospective analysis: admission and loadings both use functional outcomes.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(cluster)
  library(patchwork)
  library(RcppML)
  library(glmnet)
})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
arguments <- commandArgs(trailingOnly = TRUE)
check(length(arguments) == 0L || identical(arguments, '--render-only'),
  'Usage: Rscript 12_pooled_functional_support_programs_nmf.R [--render-only]')
render_only <- identical(arguments, '--render-only')
parent_dir <- file.path(work_dir, 'results/program_pilot/13_response_signed_functional_support_programs_nmf')
molecular_index_path <- file.path(work_dir,
  'results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
output_dir <- file.path(work_dir, 'results/program_pilot/12_pooled_functional_support_programs_nmf')
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
model_dir <- file.path(output_dir, 'models')
for (path in c(data_dir, figure_dir, model_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
nmf_tolerance <- 1e-4
nmf_max_iterations <- 100L
nmf_workers <- 4L # Memory bound for the pooled matrix; one start per rank, unchanged.
local_resamples <- 0L
local_minimum_rank <- 2L
local_rank_grid <- 2:80 # Extended pooled sweep; reuse completed ranks without refitting.
selected_rank_override <- 20L # Explicit user selection from the saved NMF partitions.
if (render_only) {
  saved_bank <- qread(file.path(model_dir, 'program_bank.qs'))
  local_rank_grid <- saved_bank$ranks$rank
  rm(saved_bank)
  message('Rendering completed ranks only: ', min(local_rank_grid), '--', max(local_rank_grid))
}
functional_folds <- 5L
functional_ridge_lambda <- 0.1
functional_input_transform <- 'full cross-fitted gain times sign(50 - assay abundance); separate positive/negative channels; no clipping'
outcome_p_threshold <- 0.05
outcome_minimum_abs_correlation <- 0
local_gene_limit <- Inf
regression_gain_p_threshold <- 0.05
minimum_loading_iqr <- 0.05
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
data.table::setDTthreads(2L)
options(warn = 1)
theme_set(theme_classic(base_size = 10))

# Require the completed panel, not an opportunistic subset of finished cells.
completion_path <- file.path(parent_dir, 'data/panel_completed.tsv')
check(file.exists(completion_path), 'Stage 13 panel is not complete; finish stage 13 before pooling')
completed <- fread(completion_path)
cell_map <- fread(file.path(parent_dir, 'data/panel_cell_mapping.tsv'))
check(setequal(completed$depmap_id, cell_map$depmap_id) && all(completed$status == 'complete'),
  'Incomplete per-cell cohort')
metadata <- fread(file.path(parent_dir, 'data/molecular_metadata.tsv.gz'))
sensitivity <- fread(file.path(parent_dir, 'data/sensitivity_alignment.tsv'))
selection <- fread(file.path(parent_dir, 'data/gene_selection.tsv.gz'))
check(all(c('functional_eligible', 'regression_significant', 'loading_variation_eligible') %in% names(selection)) &&
  all(selection[selected == TRUE, outcome_eligible & functional_eligible]),
  'Parent selection lacks completed regression/loading admission')
genes <- sort(unique(selection[selected == TRUE, gene]))
check(length(genes) >= 3L && !anyDuplicated(genes), 'Invalid shared gene union')
fwrite(selection[, .(selected_in_cells = sum(selected), tested_in_cells = .N), by = gene][gene %chin% genes],
  file.path(data_dir, 'pooled_gene_admission.tsv'), sep='\t')
message('Pooling ', length(genes), ' stage-13 retained genes across ', nrow(cell_map),
  ' cell lines; ranks ', min(local_rank_grid), '--', max(local_rank_grid))
panel <- qread(molecular_index_path)
check(setequal(panel$cell_map$depmap_id, cell_map$depmap_id) && all(genes %in% panel$genes),
  'Prepared molecular panel does not cover the pooled union')
# Preserve the exact regression outcomes from qs, rather than round-tripping
# interpolated values through a display TSV when filling missing gene fits.
check(setequal(sensitivity$observation_id, panel$sensitivity$observation_id), 'Different parent outcome rows')
sensitivity <- copy(panel$sensitivity)
reuse_audit <- list()
fwrite(data.table(status='running', phase='union-gene regressions', cells_ready=0L, cells=nrow(cell_map)),
  file.path(output_dir, 'run_status.tsv'), sep='\t')

# Compute union-gene loadings even where that gene did not pass the local screen.
# Missing assay outcomes remain absent columns; they are never filled with zero.
loading_parts <- raw_gain_parts <- response_parts <- null_parts <- observation_parts <- coefficient_parts <- fold_parts <- list()
for (id in cell_map$depmap_id) {
  input_path <- panel$input_paths[match(id, panel$cell_map$depmap_id)]
  check(file.exists(input_path), paste('Missing controlled expression for', id))
  molecular <- qread(input_path)
  md <- metadata[depmap_id == id]
  x <- molecular$response[match(md$observation_id, rownames(molecular$response)), genes, drop=FALSE]
  check(identical(rownames(x), md$observation_id) && all(is.finite(x)), 'Pooled expression alignment mismatch')
  outcomes <- sensitivity[match(md$observation_id, observation_id)]
  loading_path <- file.path(data_dir, paste0('loadings_', id, '.qs'))
  scope_name <- paste0('cell_', gsub('[^A-Za-z0-9]+', '_', md$recipient[1L]))
  regression_path <- file.path(parent_dir, scope_name, 'models/functional_regression.qs')
  check(file.exists(regression_path), paste('Missing stage-13 regression cache for', scope_name))
  signature <- digest::digest(list(genes, x, outcomes, functional_folds, functional_ridge_lambda,
    sensitivity_sign=TRUE, glmnet_version=as.character(packageVersion('glmnet')),
    parent_regression=digest::digest(regression_path, file=TRUE, algo='sha256')), algo='sha256')
  if (file.exists(loading_path)) {
    support <- qread(loading_path)
    check(identical(support$input_signature, signature), paste('Incompatible pooled loadings for', id))
  } else {
    check(!render_only, paste('No cached union-gene loadings for', id))
    cached <- qread(regression_path)$support
    check(identical(cached$settings$folds, functional_folds) &&
      identical(cached$settings$ridge_lambda, functional_ridge_lambda) &&
      identical(cached$settings$glmnet_version, as.character(packageVersion('glmnet'))) &&
      identical(cached$settings$sensitivity_sign, 'sign(50 - measured abundance percent vehicle); full gain retained'),
      paste('Parent regression settings differ for', scope_name))
    reused <- genes[genes %in% rownames(cached$loadings)]
    missing <- setdiff(genes, reused)
    message(scope_name, ': reusing ', length(reused), ' gene regressions; fitting ', length(missing), ' missing union genes')
    additional <- if (length(missing)) functional_support(x[, missing, drop=FALSE], md, outcomes,
      allow_zero=TRUE, sensitivity_sign=TRUE) else NULL
    if (!is.null(additional)) {
      check(identical(as.data.frame(additional$observations), as.data.frame(cached$observations)) &&
        identical(as.data.frame(additional$folds), as.data.frame(cached$folds)) &&
        identical(additional$settings, cached$settings),
        paste('New and cached regression contexts differ for', scope_name))
    }
    loading <- raw_gain <- matrix(NA_real_, nrow=length(genes), ncol=ncol(cached$loadings),
      dimnames=list(genes, colnames(cached$loadings)))
    loading[reused, ] <- cached$loadings[reused, , drop=FALSE]
    raw_gain[reused, ] <- cached$raw_prediction_gain[reused, , drop=FALSE]
    coefficients <- cached$coefficients[gene %chin% reused]
    if (length(missing)) {
      loading[missing, ] <- additional$loadings[missing, , drop=FALSE]
      raw_gain[missing, ] <- additional$raw_prediction_gain[missing, , drop=FALSE]
      coefficients <- rbindlist(list(coefficients, additional$coefficients))
    }
    coefficients <- coefficients[order(match(assay, c('PRISM', 'ALMANAC')), fold, match(gene, genes))]
    check(all(is.finite(loading)) && all(is.finite(raw_gain)), 'Uncomputed union-gene loading')
    support <- list(loadings=loading, raw_prediction_gain=raw_gain, observations=copy(cached$observations),
      coefficients=coefficients, folds=copy(cached$folds), assay_coverage=copy(cached$assay_coverage),
      settings=cached$settings, reused_genes=reused, newly_fitted_genes=missing)
    support$input_signature <- signature
    support$channels <- NULL # Deterministically reconstructed once after concatenation.
    qsave(support, loading_path, preset='high')
    rm(cached, additional, loading, raw_gain, coefficients)
  }
  check(identical(rownames(support$loadings), genes) &&
    identical(support$loadings, sweep(support$raw_prediction_gain, 2L, support$observations$response_sign, '*')),
    'Different gene rows or sensitivity sign across cells')
  support$observations[, depmap_id := id]
  support$coefficients[, depmap_id := id]
  support$folds[, depmap_id := id]
  loading_parts[[id]] <- support$loadings
  raw_gain_parts[[id]] <- support$raw_prediction_gain
  observation_parts[[id]] <- support$observations
  coefficient_parts[[id]] <- support$coefficients
  fold_parts[[id]] <- support$folds
  response_parts[[id]] <- x
  null_parts[[id]] <- molecular$null[[id]][, genes, drop=FALSE]
  reuse_audit[[id]] <- data.table(depmap_id=id, recipient=md$recipient[1L],
    union_genes=length(genes), reused_regressions=length(support$reused_genes),
    newly_fitted_regressions=length(support$newly_fitted_genes), loading_columns=ncol(support$loadings))
  fwrite(rbindlist(reuse_audit), file.path(data_dir, 'pooled_regression_reuse.tsv'), sep='\t')
  fwrite(data.table(status='running', phase='union-gene regressions', cells_ready=length(reuse_audit), cells=nrow(cell_map)),
    file.path(output_dir, 'run_status.tsv'), sep='\t')
  message('Union-gene loadings ready: ', cell_map$depmap_name[match(id, cell_map$depmap_id)])
  rm(molecular, support, x)
  invisible(gc())
}
loadings <- do.call(cbind, loading_parts)
check(!anyDuplicated(colnames(loadings)) && all(is.finite(loadings)), 'Invalid pooled functional columns')
channels <- cbind(pmax(loadings, 0), pmax(-loadings, 0))
colnames(channels) <- c(paste0('positive:', colnames(loadings)), paste0('negative:', colnames(loadings)))
support <- list(loadings=loadings, raw_prediction_gain=do.call(cbind, raw_gain_parts),
  channels=channels, observations=rbindlist(observation_parts),
  coefficients=rbindlist(coefficient_parts), folds=rbindlist(fold_parts),
  settings=list(folds=functional_folds, ridge_lambda=functional_ridge_lambda,
    glmnet_version=as.character(packageVersion('glmnet')), normalization='training outcome variance',
    sensitivity_sign='sign(50 - measured abundance percent vehicle); full gain retained',
    gene_admission=list(parent='13', p_threshold=regression_gain_p_threshold, minimum_loading_iqr=minimum_loading_iqr,
      rule='same-assay prediction gain and loading variation in at least one cell; union across cells')))
prepared <- list(response=do.call(rbind, response_parts))
metadata <- metadata[match(rownames(prepared$response), observation_id)]
metadata[, cell_drug := paste(depmap_id, drug, sep=':')]
null <- do.call(rbind, null_parts)
rm(loading_parts, raw_gain_parts, response_parts, null_parts, observation_parts, coefficient_parts, fold_parts, loadings, channels)
invisible(gc())
fwrite(metadata, file.path(data_dir, 'molecular_metadata.tsv.gz'), sep='\t')
fwrite(sensitivity, file.path(data_dir, 'sensitivity_alignment.tsv'), sep='\t')
fwrite(support$observations, file.path(data_dir, 'loading_columns.tsv'), sep='\t')
almanac_bar_max <- max(125, ceiling(max(sensitivity$almanac_abundance_pct, na.rm=TRUE)/25)*25)
prism_bar_max <- max(125, ceiling(max(sensitivity$prism_abundance_pct, na.rm=TRUE)/25)*25)
gene_selection <- list(pooled=selection)
outcome_test_tables <- list(pooled=fread(file.path(parent_dir, 'data/gene_outcome_correlations.tsv.gz')))
outcome_test_summary <- list(pooled=fread(file.path(parent_dir, 'data/gene_outcome_test_summary.tsv')))
prep_path <- molecular_index_path
input_signature <- digest::digest(list(genes, metadata, support$loadings,
  local_rank_grid, local_minimum_rank, local_minimum_genes, nmf_tolerance, nmf_max_iterations,
  functional_input_transform), algo='sha256')
scope <- '.' # One pooled analysis: figures/data live directly in this numbered folder.
bank_path <- file.path(model_dir, 'program_bank.qs')
if (file.exists(bank_path)) {
  bank <- qread(bank_path)
  cached_signature <- digest::digest(list(genes, metadata, support$loadings,
    bank$ranks$rank, local_minimum_rank, local_minimum_genes, nmf_tolerance, nmf_max_iterations,
    functional_input_transform), algo='sha256')
  check(identical(bank$pooled_input_signature, cached_signature), 'Pooled inputs or fit settings changed')
  check(all(bank$ranks$rank %in% local_rank_grid), 'Cached ranks are outside the requested sweep')
  bank <- restore_functional_bank(bank, prepared$response)
  if (render_only) {
    message('Recalculating mean gene silhouettes for every saved partition')
    loading_distance <- as.dist(sqrt(2 * (1 - bank$clustering_correlation)))
    for (saved_rank in bank$ranks$rank) {
      saved_partition <- bank$all_partitions[[as.character(saved_rank)]]
      value <- if (uniqueN(saved_partition) >= 2L)
        mean(cluster::silhouette(saved_partition, loading_distance)[, 'sil_width']) else NA_real_
      bank$ranks[rank == saved_rank, mean_silhouette := value]
    }
    bank$ranks <- select_silhouette(bank$ranks, selected_rank_override)
    message('Silhouettes recalculated; requested rank ', selected_rank_override)
    rm(loading_distance)
  }
  if (!identical(bank$diagnostics$selected_rank, selected_rank_override) ||
      !identical(bank$rank_override, selected_rank_override)) {
    check(selected_rank_override %in% as.integer(names(bank$all_partitions)), 'Requested partition is not saved')
    saved_partition <- bank$all_partitions[[as.character(selected_rank_override)]]
    message('Rebuilding plot profiles from the saved rank-', selected_rank_override, ' partition; no NMF fit')
    bank <- discover_local(prepared$response, null, metadata, genes, 'cell_drug', scope, 'pooled',
      support=support, sensitivity_sign=TRUE, previous=bank, rank_override=selected_rank_override)
    check(identical(bank$partition, saved_partition), 'Selected saved partition changed')
    bank$pooled_input_signature <- input_signature
    save_functional_bank(bank, bank_path)
    rm(saved_partition)
  }
  if (!identical(bank$ranks$rank, local_rank_grid)) {
    check(!render_only, 'Extend the saved rank sweep before rendering only')
    fwrite(data.table(status='running', phase='extending combined NMF to rank 80',
      cells_ready=nrow(cell_map), cells=nrow(cell_map)), file.path(output_dir, 'run_status.tsv'), sep='\t')
    reused_ranks <- copy(bank$ranks)
    previous_partitions <- bank$all_partitions
    bank <- discover_local(prepared$response, null, metadata, genes, 'cell_drug', scope, 'pooled',
      support=support, sensitivity_sign=TRUE, previous=bank, rank_override=selected_rank_override)
    # Verify actual saved partitions and fit diagnostics, rather than fitting old ranks again.
    unchanged_columns <- c('rank', 'minimum_size', 'maximum_size', 'mean_silhouette', 'nmf_fits', 'nmf_seed',
      'nmf_iterations', 'nmf_final_tolerance', 'nmf_converged', 'nmf_relative_squared_error')
    check(identical(as.data.frame(reused_ranks[, ..unchanged_columns]),
      as.data.frame(bank$ranks[rank %in% reused_ranks$rank, ..unchanged_columns])) &&
      identical(previous_partitions, bank$all_partitions[names(previous_partitions)]),
      'Previously completed NMF ranks changed during extension')
    fwrite(data.table(rank=local_rank_grid, fit_source=ifelse(local_rank_grid %in% reused_ranks$rank,
      'reused completed fit', 'new single fit'), nmf_fits=1L), file.path(data_dir, 'rank_fit_reuse.tsv'), sep='\t')
    bank$pooled_input_signature <- input_signature
    save_functional_bank(bank, bank_path)
    rm(reused_ranks, previous_partitions)
  }
} else {
  check(!render_only, 'No saved pooled fit to render')
  fwrite(data.table(status='running', phase='combined NMF ranks 2--80', cells_ready=nrow(cell_map), cells=nrow(cell_map)),
    file.path(output_dir, 'run_status.tsv'), sep='\t')
  message('Fitting combined NMF: ', length(genes), ' genes x ', ncol(support$loadings), ' signed loading columns')
  bank <- discover_local(prepared$response, null, metadata, genes, 'cell_drug', scope, 'pooled',
    support=support, sensitivity_sign=TRUE, rank_override=selected_rank_override)
  bank$pooled_input_signature <- input_signature
  save_functional_bank(bank, bank_path)
}
rm(support, null)
invisible(gc())
bank$ranks <- merge_completed_rank_diagnostics(bank)
fwrite(bank$ranks, file.path(data_dir, 'rank_diagnostics.tsv'), sep='\t')
fwrite(data.table(status='running', phase='rendering combined overview figures',
  cells_ready=nrow(cell_map), cells=nrow(cell_map)), file.path(output_dir, 'run_status.tsv'), sep='\t')
render_diagnostics(bank)
render_details <- NULL
if (is.null(bank$weights)) render_unresolved(bank) else render_details <- render_local(bank, detail_pages=FALSE)
render_functional_loadings(bank)
message('Combined overview, diagnostics, scores and loading figures ready: ', figure_dir)
fwrite(data.table(status='running', phase='full displacement detail pages; overview figures ready',
  cells_ready=nrow(cell_map), cells=nrow(cell_map)), file.path(output_dir, 'run_status.tsv'), sep='\t')
if (!is.null(render_details)) render_details()
fwrite(bank$ranks, file.path(data_dir, 'rank_diagnostics.tsv'), sep='\t')
fwrite(bank$diagnostics, file.path(data_dir, 'discovery_diagnostics.tsv'), sep='\t')
writeLines(c('Shared gene rows: union of every retained gene in completed stage 13, without a cap.',
  'Columns: separate measured cell line x drug x dose x assay functional loadings.',
  'Regressions remain per gene, cell line and assay; five whole-drug folds; ridge lambda 0.1.',
  'Loadings retain full original gain times sign(50 - measured abundance); no clipping.',
  'Gene admission: same-assay held-out drug-mean gain p <=0.05 and signed-loading IQR >=0.05 in any cell.',
  sprintf('NMF fits positive/negative loading channels directly; one start per rank; completed ranks %d--%d.',
    min(bank$ranks$rank), max(bank$ranks$rank)),
  sprintf('Rank: explicit user selection, %d programs; all completed silhouette diagnostics retained.', selected_rank_override),
  'Displacement and loading heatmaps have the same exported gene and condition order.',
  'Retrospective outcome-informed analysis; no independent evaluation is claimed.'), file.path(output_dir,'README.txt'))
writeLines(capture.output(sessionInfo()), file.path(data_dir,'sessionInfo.txt'))
fwrite(data.table(script='12_pooled_functional_support_programs_nmf', cells=nrow(cell_map), genes=length(genes),
  loading_columns=ncol(bank$functional_support$loadings), selected_rank=bank$diagnostics$selected_rank,
  code_sha256=digest::digest(file.path(code_dir,'12_pooled_functional_support_programs_nmf.R'),file=TRUE,algo='sha256'),
  input_signature=input_signature), file.path(data_dir,'run_manifest.tsv'),sep='\t')
figures <- c('silhouette.pdf', 'cophenetic.pdf', 'diagnostics.pdf', 'gene_correlation.pdf',
  'functional_loadings.pdf', 'functional_loadings_focal.pdf')
figures <- c(figures, if (is.null(bank$weights)) c('displacement_unresolved.pdf', 'focal_panel.pdf') else
  c('displacement_overview.pdf', 'displacement_pages.pdf', 'program_scores.pdf'))
fwrite(data.table(artifact=figures, path=file.path(figure_dir, figures)), file.path(output_dir, 'artifact_index.tsv'), sep='\t')
check(all(file.exists(file.path(figure_dir, figures))), 'A required combined figure is missing')
fwrite(data.table(status='complete', phase='all combined figures ready', cells_ready=nrow(cell_map), cells=nrow(cell_map)),
  file.path(output_dir, 'run_status.tsv'), sep='\t')
message('Pooled figures ready: ', file.path(output_dir,scope,'figures'))
