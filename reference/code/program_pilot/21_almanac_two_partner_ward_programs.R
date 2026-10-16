#!/usr/bin/env Rscript
# Retrospective ALMANAC interaction study from two measured same-cell signatures.
# Run from /work with Rscript (the workspace .Rprofile activates its R library).
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
  library(cluster)
})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))

# Use the same gene admission, regression penalty and Ward settings as stage 15.
output_dir <- file.path(work_dir, 'results/program_pilot/21_almanac_two_partner_ward_programs')
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
model_dir <- file.path(output_dir, 'models')
functional_folds <- 5L
functional_ridge_lambda <- 0.1
outcome_p_threshold <- 0.05
regression_gain_p_threshold <- 0.05
minimum_loading_iqr <- 0.05
local_minimum_rank <- 2L
local_rank_grid <- 2:80
local_minimum_genes <- 12L
local_resamples <- 0L
local_rank_selection <- "silhouette_elbow"
functional_input_transform <- "50/50 mean gene-gene Pearson correlation of signed functional loadings and original displacements; Ward.D2 on sqrt(2*(1-C))"
settings <- list(seed = seed, functional_folds = functional_folds, ridge_lambda = functional_ridge_lambda,
  screen_p = outcome_p_threshold, admission_p = regression_gain_p_threshold, loading_iqr = minimum_loading_iqr,
  ranks = local_rank_grid, minimum_seed_genes = local_minimum_genes, rank_rule = local_rank_selection,
  endpoint = "unclipped Bliss excess inhibition, percentage points",
  regression_inputs = "raw displacement of the same gene for partner A and partner B")
arguments <- commandArgs(trailingOnly = TRUE)
check(all(arguments %in% c('--prepare-only', '--cell') | grepl('^ACH-[0-9]+$', arguments)),
  'Supported arguments: --prepare-only or --cell ACH-000681')
selected_cell <- NULL
if ('--cell' %in% arguments) {
  check(length(arguments) == 2L && arguments[1] == '--cell', 'Use --cell followed by one DepMap ID')
  selected_cell <- arguments[2]
}
check(!('--prepare-only' %in% arguments) || identical(arguments, '--prepare-only'),
  '--prepare-only must be used alone')
for (path in c(output_dir, data_dir, figure_dir, model_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
data.table::setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
Sys.setenv(OMP_NUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1', MKL_NUM_THREADS = '1')
set.seed(seed)
theme_set(theme_classic(base_size = 10))
options(warn = 1)

panel_path <- file.path(work_dir,
  'results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
raw_path <- file.path(work_dir, 'almanac_conditionals/raw/ComboDrugGrowth_Nov2017.csv')
mapping_path <- file.path(work_dir, 'almanac_conditionals/repaired/results/drug_mapping.tsv')
script_path <- file.path(code_dir, '21_almanac_two_partner_ward_programs.R')
check(all(file.exists(c(panel_path, raw_path, mapping_path))), 'Missing ALMANAC or prepared molecular inputs')
panel <- qread(panel_path)
cell_map <- copy(panel$cell_map[almanac_present == TRUE])
metadata <- copy(panel$metadata[depmap_id %chin% cell_map$depmap_id])
check(!anyDuplicated(cell_map$depmap_id) && !anyDuplicated(cell_map$CELLNAME) && !anyNA(cell_map$CELLNAME),
  'Invalid ALMANAC cell mapping')
check(!anyDuplicated(metadata$observation_id) && all(is.finite(metadata$dose)) && all(metadata$dose > 0),
  'Invalid molecular metadata identifiers or doses')
mapping <- fread(mapping_path)[primary == TRUE, .(nsc, identity, drug = panel_name(tahoe), nci_name, mapping)]
check(!anyDuplicated(mapping$nsc) && !anyNA(mapping), 'Invalid primary NSC mapping')
input_paths <- c(panel_path, raw_path, mapping_path, panel$input_paths[match(cell_map$depmap_id, panel$cell_map$depmap_id)])
check(all(file.exists(input_paths)), 'Missing prepared cell matrix')
input_manifest <- data.table(path = input_paths,
  sha256 = vapply(input_paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
input_signature <- digest::digest(list(inputs = input_manifest, settings = settings), algo = 'sha256')
fwrite(input_manifest, file.path(data_dir, 'input_manifest.tsv'), sep = '\t')
fwrite(mapping, file.path(data_dir, 'drug_identity_mapping.tsv'), sep = '\t')
fwrite(cell_map, file.path(data_dir, 'cell_mapping.tsv'), sep = '\t')
qsave(settings, file.path(model_dir, 'settings.qs'))
prepared_path <- file.path(model_dir, 'matched_combinations.qs')

# Reconstruct singles in their assay context before aggregating combination rows.
# Nonnegative raw signal includes zero survival; negative signals and invalid
# controls are excluded explicitly. Survival above vehicle is retained unclipped.
if (file.exists(prepared_path)) {
  prepared <- qread(prepared_path)
  check(identical(prepared$input_signature, input_signature), 'Cached combination input/settings mismatch')
  combinations <- prepared$combinations
  queries <- prepared$queries
  message('Reused matched ALMANAC endpoints and recipient-specific dose audit')
} else {
  raw_columns <- c('COMBODRUGSEQ', 'SCREENER', 'STUDY', 'TESTDATE', 'CELLNAME', 'NSC1', 'NSC2',
    'SAMPLE1', 'SAMPLE2', 'CONC1', 'CONC2', 'CONCUNIT1', 'CONCUNIT2', 'TESTVALUE', 'CONTROLVALUE', 'VALID')
  raw <- fread(raw_path, select = raw_columns)
  check(!anyDuplicated(raw$COMBODRUGSEQ), 'Duplicate raw ALMANAC record ID')
  raw <- raw[CELLNAME %chin% cell_map$CELLNAME & NSC1 %in% mapping$nsc]
  raw[, valid_endpoint := VALID == 'Y' & is.finite(TESTVALUE) & TESTVALUE >= 0 &
    is.finite(CONTROLVALUE) & CONTROLVALUE > 0 & CONCUNIT1 == 'M']
  fwrite(raw[, .(records = .N, invalid_endpoints = sum(!valid_endpoint), zero_signal = sum(TESTVALUE == 0)),
    by = .(CELLNAME, SCREENER)], file.path(data_dir, 'raw_endpoint_audit.tsv'), sep = '\t')
  raw[, survival := TESTVALUE / CONTROLVALUE]
  context <- c('SCREENER', 'STUDY', 'TESTDATE', 'CELLNAME')
  singles <- raw[is.na(NSC2) & valid_endpoint == TRUE, .(survival = mean(survival), single_records = .N),
    by = c(context, 'NSC1', 'SAMPLE1', 'CONC1')]
  pairs <- raw[!is.na(NSC2) & NSC2 %in% mapping$nsc & NSC1 != NSC2]
  pairs[, valid_endpoint := valid_endpoint & CONCUNIT2 == 'M']
  single_a <- copy(singles)
  setnames(single_a, c('survival', 'single_records'), c('survival_a', 'single_records_a'))
  pairs <- merge(pairs, single_a, by = c(context, 'NSC1', 'SAMPLE1', 'CONC1'), all.x = TRUE, sort = FALSE)
  single_b <- copy(singles)
  setnames(single_b, c('NSC1', 'SAMPLE1', 'CONC1', 'survival', 'single_records'),
    c('NSC2', 'SAMPLE2', 'CONC2', 'survival_b', 'single_records_b'))
  pairs <- merge(pairs, single_b, by = c(context, 'NSC2', 'SAMPLE2', 'CONC2'), all.x = TRUE, sort = FALSE)
  check(!anyDuplicated(pairs$COMBODRUGSEQ), 'Matched-single joins multiplied combination records')
  pairs[, reason := fcase(!valid_endpoint, 'invalid combination endpoint or concentration units',
    !is.finite(survival_a) | !is.finite(survival_b), 'missing valid same-context single', default = 'eligible endpoint')]
  fwrite(pairs[, .(records = .N), by = .(CELLNAME, reason)],
    file.path(data_dir, 'matched_single_exclusions.tsv'), sep = '\t')
  pairs <- pairs[reason == 'eligible endpoint']
  check(all(is.finite(pairs$CONC1) & pairs$CONC1 > 0 & is.finite(pairs$CONC2) & pairs$CONC2 > 0),
    'Invalid positive combination concentrations')
  # Canonical NSC order is for storage and display only; models are symmetric.
  reverse <- which(pairs$NSC1 > pairs$NSC2)
  for (fields in list(c('NSC1', 'NSC2'), c('SAMPLE1', 'SAMPLE2'), c('CONC1', 'CONC2'),
    c('survival_a', 'survival_b'), c('single_records_a', 'single_records_b'))) {
    saved <- pairs[[fields[1]]][reverse]
    set(pairs, reverse, fields[1], pairs[[fields[2]]][reverse])
    set(pairs, reverse, fields[2], saved)
  }
  pairs[, `:=`(depmap_id = cell_map$depmap_id[match(CELLNAME, cell_map$CELLNAME)],
    drug_a = mapping$drug[match(NSC1, mapping$nsc)], drug_b = mapping$drug[match(NSC2, mapping$nsc)],
    dose_a = CONC1 * 1e6, dose_b = CONC2 * 1e6,
    pair = paste(NSC1, NSC2, sep = '::'), bliss_pp = 100 * (survival_a * survival_b - survival))]
  fwrite(pairs[drug_a == drug_b, .(records = .N), by = .(drug_a, NSC1, NSC2)],
    file.path(data_dir, 'same_molecular_compound_exclusions.tsv'), sep = '\t')
  pairs <- pairs[drug_a != drug_b]
  fwrite(pairs[, .(COMBODRUGSEQ, depmap_id, pair, SCREENER, STUDY, TESTDATE, NSC1, NSC2,
    dose_a, dose_b, survival_a, survival_b, survival, bliss_pp)],
    file.path(data_dir, 'raw_bliss_reconstruction.tsv.gz'), sep = '\t')
  combinations <- pairs[, .(bliss_pp = mean(bliss_pp), single_a_pct = 100 * mean(survival_a),
    single_b_pct = 100 * mean(survival_b), expected_pct = 100 * mean(survival_a * survival_b),
    combination_pct = 100 * mean(survival), raw_records = .N, studies = uniqueN(STUDY)),
    by = .(depmap_id, CELLNAME, pair, NSC1, NSC2, drug_a, drug_b, dose_a, dose_b, SCREENER)]
  check(max(abs(combinations$bliss_pp - (combinations$expected_pct - combinations$combination_pct))) < 1e-10,
    'Bliss endpoint reconstruction mismatch')
  queries <- unique(rbind(combinations[, .(depmap_id, drug = drug_a, dose = dose_a)],
    combinations[, .(depmap_id, drug = drug_b, dose = dose_b)]))
  setorder(queries, depmap_id, drug, dose)
  queries[, `:=`(query_id = .I, status = 'no same-cell molecular drug', observation_lo = NA_character_,
    observation_hi = NA_character_, dose_lo = NA_real_, dose_hi = NA_real_, fraction = NA_real_)]
  for (i in seq_len(nrow(queries))) {
    query <- queries[i]
    md <- metadata[depmap_id == query$depmap_id & drug == query$drug][order(dose)]
    if (!nrow(md)) next
    check(!anyDuplicated(md$dose), 'Ambiguous same-cell molecular dose')
    exact <- which(abs(log(md$dose / query$dose)) < 1e-6)
    if (length(exact)) {
      check(length(exact) == 1L, 'Ambiguous exact molecular query')
      lo <- hi <- exact
      interpolation_fraction <- 0
      query_status <- 'exact'
    } else {
      if (query$dose < min(md$dose) || query$dose > max(md$dose)) {
        queries[i, status := 'outside same-cell measured dose range']
        next
      }
      lo <- max(which(md$dose < query$dose))
      hi <- min(which(md$dose > query$dose))
      interpolation_fraction <- log(query$dose / md$dose[lo]) / log(md$dose[hi] / md$dose[lo])
      query_status <- 'log-dose interpolated'
    }
    queries[i, `:=`(status = query_status, observation_lo = md$observation_id[lo], observation_hi = md$observation_id[hi],
      dose_lo = md$dose[lo], dose_hi = md$dose[hi], fraction = interpolation_fraction)]
  }
  a_keys <- paste(combinations$depmap_id, combinations$drug_a, combinations$dose_a)
  b_keys <- paste(combinations$depmap_id, combinations$drug_b, combinations$dose_b)
  query_keys <- paste(queries$depmap_id, queries$drug, queries$dose)
  check(!anyDuplicated(query_keys), 'Duplicate molecular query key')
  combinations[, `:=`(query_a = queries$query_id[match(a_keys, query_keys)],
    query_b = queries$query_id[match(b_keys, query_keys)])]
  check(!anyNA(combinations[, .(query_a, query_b)]), 'Missing query join')
  combinations[, `:=`(status_a = queries$status[query_a], status_b = queries$status[query_b])]
  fwrite(combinations[, .(conditions = .N, pairs = uniqueN(pair)), by = .(CELLNAME, status_a, status_b)],
    file.path(data_dir, 'molecular_coverage_exclusions.tsv'), sep = '\t')
  combinations <- combinations[status_a %chin% c('exact', 'log-dose interpolated') &
    status_b %chin% c('exact', 'log-dose interpolated')]
  setorder(combinations, depmap_id, pair, dose_a, dose_b, SCREENER)
  combinations[, `:=`(combination_id = sprintf('combo_%06d', .I),
    recipient = cell_map$depmap_name[match(depmap_id, cell_map$depmap_id)])]
  check(nrow(combinations) > 0 && all(is.finite(combinations$bliss_pp)), 'No finite eligible combination outcomes')
  check(!anyDuplicated(combinations[, .(depmap_id, pair, dose_a, dose_b, SCREENER)]),
    'Duplicate aggregated combination condition')
  qsave(list(input_signature = input_signature, combinations = combinations, queries = queries), prepared_path)
  rm(raw, singles, single_a, single_b, pairs)
  invisible(gc())
}
fwrite(queries, file.path(data_dir, 'dose_query_audit.tsv'), sep = '\t')
fwrite(combinations, file.path(data_dir, 'matched_combinations.tsv.gz'), sep = '\t')
coverage <- combinations[, .(conditions = .N, pairs = uniqueN(pair),
  drugs = uniqueN(c(NSC1, NSC2)), raw_records = sum(raw_records),
  positive_interaction = sum(bliss_pp > 0), negative_interaction = sum(bliss_pp < 0),
  zero_interaction = sum(bliss_pp == 0)), by = .(depmap_id, recipient)]
fwrite(coverage, file.path(data_dir, 'combination_coverage.tsv'), sep = '\t')
print(coverage)
if ('--prepare-only' %in% arguments) quit(save = 'no')

# The stage-15 per-gene regression with two raw partner columns. Outcome scaling,
# ridge penalty, complete-intervention folds and prediction-gain definition are
# unchanged. A combination pair replaces a single drug as the held-out unit.
functional_support_pair <- function(a, b, md) {
  check(identical(dimnames(a), dimnames(b)) && identical(rownames(a), md$combination_id),
    'Partner signatures and combination outcomes do not align')
  pairs <- sort(unique(md$pair))
  check(length(pairs) >= functional_folds && all(is.finite(md$bliss_pp)),
    'Insufficient pairs or invalid interaction endpoint')
  fold <- 1L + (match(md$pair, pairs) - 1L) %% functional_folds
  predicted <- matrix(NA_real_, nrow(a), ncol(a), dimnames = dimnames(a))
  baseline <- normalization <- rep(NA_real_, nrow(a))
  coefficients <- list()
  for (k in seq_len(functional_folds)) {
    train <- which(fold != k)
    test <- which(fold == k)
    check(!any(md$pair[test] %in% md$pair[train]), 'Pair leakage in functional folds')
    y_mean <- mean(md$bliss_pp[train])
    y_sd <- sqrt(mean((md$bliss_pp[train] - y_mean)^2))
    check(is.finite(y_sd) && y_sd > 0, 'Constant training interaction outcome')
    baseline[test] <- y_mean
    normalization[test] <- y_sd^2
    mean_a <- colMeans(a[train, , drop = FALSE])
    mean_b <- colMeans(b[train, , drop = FALSE])
    sd_a <- matrixStats::colSds(a[train, , drop = FALSE])
    sd_b <- matrixStats::colSds(b[train, , drop = FALSE])
    beta_a <- beta_b <- numeric(ncol(a))
    status <- rep('ridge fit', ncol(a))
    for (j in seq_len(ncol(a))) {
      scale_a <- if (sd_a[j] > 0) sd_a[j] else 1
      scale_b <- if (sd_b[j] > 0) sd_b[j] else 1
      z <- cbind(A = (a[train, j] - mean_a[j]) / scale_a,
        B = (b[train, j] - mean_b[j]) / scale_b)
      inactive <- which(c(sd_a[j], sd_b[j]) == 0)
      if (length(inactive) == 2L) {
        predicted[test, j] <- y_mean
        status[j] <- 'constant training gene in both partners: baseline prediction'
        next
      }
      fit <- glmnet::glmnet(z, (md$bliss_pp[train] - y_mean) / y_sd, alpha = 0,
        lambda = functional_ridge_lambda, standardize = FALSE, intercept = FALSE,
        thresh = 1e-10, exclude = inactive)
      beta <- as.matrix(coef(fit))[-1L, 1L]
      beta_a[j] <- beta[1L]
      beta_b[j] <- beta[2L]
      z_test <- cbind(A = (a[test, j] - mean_a[j]) / scale_a,
        B = (b[test, j] - mean_b[j]) / scale_b)
      if (length(inactive)) z_test[, inactive] <- 0
      predicted[test, j] <- y_mean + y_sd * as.numeric(z_test %*% beta)
      if (length(inactive)) status[j] <- 'one constant training partner: that coefficient excluded'
    }
    coefficients[[k]] <- data.table(assay = 'ALMANAC', fold = k, gene = colnames(a),
      training_conditions = length(train), training_drugs = uniqueN(md$pair[train]),
      training_gene_mean_a = mean_a, training_gene_mean_b = mean_b,
      training_gene_sd_a = sd_a, training_gene_sd_b = sd_b,
      training_outcome_mean = y_mean, training_outcome_sd = y_sd,
      standardized_beta_a = beta_a, standardized_beta_b = beta_b,
      ridge_lambda = functional_ridge_lambda, fit_status = status)
    message(md$recipient[1L], ': joint gene regression fold ', k, '/', functional_folds, '; ', ncol(a), ' genes')
  }
  check(all(is.finite(predicted)) && all(normalization > 0), 'Nonfinite functional predictions or scaling')
  residual <- sweep(predicted, 1L, md$bliss_pp, '-')
  gain <- sweep(sweep(-residual^2, 1L, (md$bliss_pp - baseline)^2, '+'), 1L, normalization, '/')
  raw_gain <- t(gain)
  colnames(raw_gain) <- paste0('ALMANAC:', md$combination_id)
  loading <- sweep(raw_gain, 2L, sign(md$bliss_pp), '*')
  channels <- cbind(pmax(loading, 0), pmax(-loading, 0))
  colnames(channels) <- c(paste0('positive:', colnames(loading)), paste0('negative:', colnames(loading)))
  observations <- md[, .(assay = 'ALMANAC', observation_id = combination_id, drug = pair,
    dose = dose_a, dose_a, dose_b, outcome = bliss_pp)]
  observations[, `:=`(fold = fold, baseline_prediction = baseline,
    training_outcome_variance = normalization, response_sign = sign(outcome))]
  list(loadings = loading, raw_prediction_gain = raw_gain, channels = channels,
    predictions = t(predicted), observations = observations, coefficients = rbindlist(coefficients),
    folds = observations[, .(assay, observation_id, drug, fold)],
    assay_coverage = data.table(assay = 'ALMANAC', conditions = nrow(md), drugs = length(pairs), included = TRUE),
    settings = list(folds = functional_folds, ridge_lambda = functional_ridge_lambda,
      glmnet_version = as.character(packageVersion('glmnet')), normalization = 'training outcome variance',
      sensitivity_sign = 'sign(Bliss excess inhibition); full gain retained',
      regression_inputs = 'raw displacement A and raw displacement B'))
}

# Combination plots retain stage-15 program blocks, signed colors, control bars,
# top 12 seed genes in detail pages and its Ward program-score column ordering.
# Only the partner panel stacking and the functional endpoint change.
render_combination <- function(bank, a, b, md, destination) {
  data <- file.path(destination, 'data')
  figures <- file.path(destination, 'figures')
  programs <- if (is.null(bank$weights)) character() else colnames(bank$weights)
  if (is.null(bank$weights)) {
    row_order <- bank$tree$order
    representative_rows <- head(row_order, 72L)
    column_order <- hclust(dist(cbind(a[, bank$genes, drop = FALSE], b[, bank$genes, drop = FALSE])),
      method = 'ward.D2')$order
    row_meta <- NULL
    palette <- NULL
    col_tree <- NULL
  } else {
    row_order <- order(bank$partition, match(seq_along(bank$genes), bank$gene_tree_order))
    score_a <- a[, bank$genes, drop = FALSE] %*% bank$weights
    score_b <- b[, bank$genes, drop = FALSE] %*% bank$weights
    standardized <- cbind(sweep(score_a, 2L, bank$score_sd, '/'), sweep(score_b, 2L, bank$score_sd, '/'))
    col_tree <- hclust(dist(standardized), method = 'ward.D2')
    column_order <- col_tree$order
    row_meta <- data.frame(Program = factor(programs[bank$partition[row_order]], levels = programs),
      row.names = bank$genes[row_order])
    palette <- list(Program = setNames(hcl.colors(length(programs), 'Dark 3'), programs))
    representatives <- unlist(lapply(seq_along(programs), function(j) {
      rows <- which(bank$partition == j)
      head(rows[order(bank$weights[rows, j], decreasing = TRUE)], local_display_genes_per_program)
    }))
    representative_rows <- row_order[row_order %in% representatives]
    qsave(list(tree = col_tree, gene_order = bank$genes[row_order],
      column_order = md$combination_id[column_order], score_scaling = 'separate A/B signed scores divided by program SD'),
      file.path(destination, 'models/program_score_column_tree.qs'))
    fwrite(data.table(combination_id = md$combination_id,
      setNames(as.data.table(score_a), paste0('A_', programs)),
      setNames(as.data.table(score_b), paste0('B_', programs))), file.path(data, 'partner_program_scores.tsv'), sep = '\t')
    fwrite(data.table(gene = bank$genes, as.data.table(bank$weights)), file.path(data, 'gene_weights.tsv.gz'), sep = '\t')
    fwrite(bank$programs, file.path(data, 'program_diagnostics.tsv'), sep = '\t')
  }
  a <- a[, bank$genes, drop = FALSE]
  b <- b[, bank$genes, drop = FALSE]
  limit <- max(0.25, as.numeric(quantile(abs(c(a, b)), 0.98)))
  load_limit <- max(0.1, as.numeric(quantile(abs(bank$functional_support$loadings), 0.98)))
  md <- md[column_order]
  a <- a[column_order, , drop = FALSE]
  b <- b[column_order, , drop = FALSE]
  col_meta <- data.frame(`A dose` = sprintf('%g', md$dose_a), `B dose` = sprintf('%g', md$dose_b),
    check.names = FALSE, row.names = md$combination_id)
  rownames(a) <- rownames(b) <- md$combination_id
  check(identical(rownames(a), rownames(col_meta)) && identical(rownames(a), rownames(b)),
    'Stacked partner panel column alignment failed')
  bliss_range <- range(c(0, md$bliss_pp))
  loading <- bank$functional_support$loadings[, paste0('ALMANAC:', md$combination_id), drop = FALSE]
  colnames(loading) <- md$combination_id
  # The two outcome annotations from the original plot become the combination
  # Bliss bar and its continuous value strip. Separate molecular-support bars
  # are associated with the appropriate partner panel.
  top_args <- function(columns, partner) {
    support <- if (partner == 'A') md$rna_control_ratio_a[columns] else md$rna_control_ratio_b[columns]
    list(`ALMANAC Bliss pp` = ComplexHeatmap::anno_barplot(md$bliss_pp[columns],
      gp = grid::gpar(fill = '#B24C35', col = NA), ylim = bliss_range,
      height = grid::unit(11, 'mm')),
      `Bliss value` = ComplexHeatmap::anno_simple(md$bliss_pp[columns],
        col = circlize::colorRamp2(c(bliss_range[1L], 0, bliss_range[2L]), c('#2166AC', 'white', '#B2182B'))),
      `RNA / control RMS` = ComplexHeatmap::anno_barplot(support,
        gp = grid::gpar(fill = '#777777', col = NA), height = grid::unit(8, 'mm')),
      show_legend = FALSE)
  }
  make_panel <- function(matrix, rows, columns, partner, show_labels, caption, loading_map = FALSE) {
    gene_names <- bank$genes[rows]
    matrix <- if (loading_map) matrix[gene_names, columns, drop = FALSE] else t(matrix[columns, rows, drop = FALSE])
    rm <- if (is.null(row_meta)) NULL else row_meta[gene_names, , drop = FALSE]
    cm <- col_meta[columns, , drop = FALSE]
    check(identical(colnames(matrix), rownames(cm)), 'Map annotations do not match combinations')
    labels <- paste(md$drug_a[columns], sprintf('%g uM', md$dose_a[columns]), '+',
      md$drug_b[columns], sprintf('%g uM', md$dose_b[columns]))
    draw_heatmap(matrix, col_anno_df = cm, col_anno_cols = names(cm), row_anno_df = rm,
      row_anno_cols = if (is.null(rm)) NULL else 'Program',
      ann_top_args = if (loading_map) list(show_legend = FALSE) else top_args(columns, partner),
      palette_overrides = palette, hm_name = if (loading_map) 'signed loading' else paste(partner, 'log2 delta'),
      scale_rows = FALSE, value_range = if (loading_map) c(-load_limit, 0, load_limit) else c(-limit, 0, limit),
      value_colors = c('#2166AC', 'white', '#B2182B'),
      show_row_names = show_labels, show_column_names = show_labels && partner == 'B', column_labels = labels,
      row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 5),
      row_split = if (is.null(rm)) NULL else rm$Program,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 10),
      row_title = if (is.null(rm)) paste('Partner', partner) else
        paste(if (loading_map) 'Joint' else partner, levels(droplevels(rm$Program))),
      row_title_rot = 0, title = paste(md$recipient[1L], '|', caption),
      height = grid::unit(if (show_labels) max(40, length(rows) * 2.1) else 72, 'mm'))
  }
  columns <- seq_len(nrow(md))
  overview <- make_panel(a, row_order, columns, 'A', FALSE, 'original signed displacement | partner A') %v%
    make_panel(b, row_order, columns, 'B', FALSE, 'original signed displacement | partner B')
  save_heatmap(overview, 'displacement_overview', 16, 10.8, figures)
  save_heatmap(overview, 'displacement_program_score_clustered', 16, 10.8, figures)
  joint_map <- make_panel(loading, row_order, columns, 'joint', FALSE,
    'signed functional loadings | same combination order', loading_map = TRUE)
  save_heatmap(joint_map, 'functional_loadings_program_score_clustered', 16, 6, figures)
  if (!is.null(bank$weights)) {
    correlation <- bank$correlation[row_order, row_order, drop = FALSE]
    hm <- draw_heatmap(correlation, col_anno_df = row_meta, col_anno_cols = 'Program',
      row_anno_df = row_meta, row_anno_cols = 'Program', palette_overrides = palette,
      ann_top_args = list(show_legend = FALSE), hm_name = 'Pearson r', scale_rows = FALSE,
      value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'),
      row_split = row_meta$Program, column_split = row_meta$Program,
      cluster_row_slices = FALSE, cluster_column_slices = FALSE,
      row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 6),
      show_row_names = FALSE, title = paste(md$recipient[1L], '| gene correlation | all admitted genes'))
    save_heatmap(hm, 'gene_correlation', 7.8, 7.2, figures)
    score_meta <- data.frame(Program = factor(programs, levels = programs), row.names = programs)
    score_limit <- max(0.1, as.numeric(quantile(abs(c(score_a, score_b)), 0.98)))
    score_maps <- lapply(c('A', 'B'), function(partner) {
      matrix <- t(if (partner == 'A') score_a else score_b)[, column_order, drop = FALSE]
      colnames(matrix) <- md$combination_id
      draw_heatmap(matrix, col_anno_df = col_meta, col_anno_cols = names(col_meta),
        row_anno_df = score_meta, row_anno_cols = 'Program', palette_overrides = palette,
        ann_top_args = top_args(columns, partner), hm_name = paste(partner, 'weighted delta'),
        scale_rows = FALSE, value_range = c(-score_limit, 0, score_limit),
        value_colors = c('#2166AC', 'white', '#B2182B'), show_row_names = TRUE,
        row_title = paste('Partner', partner), row_title_rot = 90,
        title = paste(md$recipient[1L], '| original signed program scores'), height = grid::unit(35, 'mm'))
    })
    save_heatmap(score_maps[[1L]] %v% score_maps[[2L]], 'program_scores', 16, 8, figures)
  }
  # Detail pages use the original 48-column limit and representative-gene rule.
  pages <- split(columns, ceiling(columns / 48L))
  # Two sets of original assay/control tracks, row-block gaps and combination
  # labels need additional height above the two signed signature panels.
  page_height <- 6.8 + 2 * max(40, length(representative_rows) * 2.1) / 25.4
  cairo_pdf(file.path(figures, 'displacement_pages.pdf'), width = 13.5, height = page_height)
  for (page in seq_along(pages)) {
    maps <- make_panel(a, representative_rows, pages[[page]], 'A', TRUE,
      paste('partner A | detail page', page, '/', length(pages))) %v%
      make_panel(b, representative_rows, pages[[page]], 'B', TRUE, 'partner B | same combinations')
    ComplexHeatmap::draw(maps, merge_legends = TRUE)
  }
  dev.off()
  preview <- make_panel(a, representative_rows, pages[[1L]], 'A', TRUE, 'partner A | detail page 1') %v%
    make_panel(b, representative_rows, pages[[1L]], 'B', TRUE, 'partner B | same combinations')
  grDevices::png(file.path(figures, 'displacement_page_01.png'), width = 13.5, height = page_height,
    units = 'in', res = 160)
  ComplexHeatmap::draw(preview, merge_legends = TRUE)
  dev.off()
  fwrite(md[, .(combination_id, pair, drug_a, drug_b, dose_a, dose_b, SCREENER, bliss_pp)],
    file.path(data, 'condition_order_program_score_clustered.tsv'), sep = '\t')
  fwrite(data.table(gene = bank$genes[row_order], plotted_order = seq_along(row_order),
    detail_label = row_order %in% representative_rows), file.path(data, 'gene_order.tsv'), sep = '\t')
  qsave(list(partner_a = a, partner_b = b, column_order = md$combination_id),
    file.path(destination, 'models/plotted_partner_signatures.qs'))
  fwrite(bank$ranks, file.path(data, 'rank_diagnostics.tsv'), sep = '\t')
  writeLines(c('Same stage-15 Ward method; two raw partner predictors; ALMANAC Bliss endpoint.',
    'Partner A is above B with identical combination columns and gene order.',
    'Original signed displacement units, blue-white-red colors, program blocks and control bars.',
    'Columns use the original Ward.D2 program-score ordering extended to both separate score vectors.',
    'Full-data screening and gene admission remain retrospective, as in stage 15.',
    bank$diagnostics$selection_status), file.path(destination, 'README.txt'))
}

# Cache signatures cover scientific calculations; changes to labels/layout alone
# do not rerun gene regressions or Ward cuts. The original shared code is hashed.
calculation_signature <- digest::digest(list(settings = settings, regression = body(functional_support_pair),
  shared_code = digest::digest(file.path(code_dir, 'functional_programs.R'), file = TRUE, algo = 'sha256')),
  algo = 'sha256')
original_paths <- c(file.path(code_dir, '15_mean_functional_displacement_correlation_programs_ward.R'),
  file.path(code_dir, 'functional_programs.R'), file.path(code_dir, 'config.R'), file.path(work_dir, 'src/R/lib_heatmap.R'))
original_hashes <- data.table(path = original_paths,
  sha256 = vapply(original_paths, digest::digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(original_hashes, file.path(data_dir, 'original_pipeline_hashes.tsv'), sep = '\t')
cells_to_run <- cell_map$depmap_id
if (!is.null(selected_cell)) {
  check(selected_cell %chin% cells_to_run, 'Requested cell has no ALMANAC molecular overlap')
  cells_to_run <- selected_cell
}
inventories <- correlation_tables <- list()
for (id in cells_to_run) {
  md <- copy(combinations[depmap_id == id])
  scope <- paste0('cell_', gsub('[^A-Za-z0-9]+', '_', md$recipient[1L]))
  destination <- file.path(output_dir, scope)
  for (path in c(destination, file.path(destination, 'data'), file.path(destination, 'figures'),
    file.path(destination, 'models'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
  check(identical(colnames(molecular$response), panel$genes) && all(is.finite(molecular$response)),
    'Prepared gene universe or molecular matrix mismatch')
  q <- queries[depmap_id == id & query_id %in% c(md$query_a, md$query_b)]
  lo <- match(q$observation_lo, rownames(molecular$response))
  hi <- match(q$observation_hi, rownames(molecular$response))
  check(!anyNA(c(lo, hi)) && all(q$fraction >= 0 & q$fraction <= 1), 'Invalid molecular interpolation bounds')
  query_x <- molecular$response[lo, , drop = FALSE] * (1 - q$fraction) +
    molecular$response[hi, , drop = FALSE] * q$fraction
  a <- query_x[match(md$query_a, q$query_id), , drop = FALSE]
  b <- query_x[match(md$query_b, q$query_id), , drop = FALSE]
  rownames(a) <- rownames(b) <- md$combination_id
  check(all(is.finite(a)) && all(is.finite(b)), 'Nonfinite partner displacement')
  qsave(list(queries = q, displacement = query_x), file.path(destination, 'models/partner_queries.qs'))
  null <- molecular$null[[id]]
  null_energy <- mean(rowMeans(null^2))
  energy_a <- rowMeans(a^2)
  energy_b <- rowMeans(b^2)
  md[, `:=`(rna_control_ratio_a = sqrt(energy_a / null_energy),
    rna_control_ratio_b = sqrt(energy_b / null_energy))]
  # Apply the original per-panel outcome screen to the two constituent signatures.
  tests <- rbindlist(lapply(c('A', 'B'), function(partner) {
    result <- gene_outcome_tests(if (partner == 'A') a else b, md$bliss_pp, md$pair)
    result[, `:=`(partner = partner, recipient = md$recipient[1L], outcome_scale = 'Bliss excess inhibition pp')]
    result
  }))
  tests[, q_value := p.adjust(p_value, method = 'BH'), by = partner]
  tests[, significant := !is.na(p_value) & p_value <= outcome_p_threshold]
  candidates <- sort(unique(tests[significant == TRUE, gene]))
  check(length(candidates) >= 3L, paste('Fewer than three outcome-screened genes for', scope))
  fwrite(tests, file.path(destination, 'data/gene_outcome_correlations.tsv.gz'), sep = '\t')
  regression_signature <- digest::digest(list(calculation_signature, input_signature, id,
    metadata = md, candidates = candidates, a = a[, candidates, drop = FALSE], b = b[, candidates, drop = FALSE]),
    algo = 'sha256')
  regression_path <- file.path(destination, 'models/functional_regression.qs')
  if (file.exists(regression_path)) {
    regression <- qread(regression_path)
    check(identical(regression$input_signature, regression_signature), 'Saved combination regression inputs changed')
    support <- regression$support
    message('Reused joint raw-partner gene regressions: ', scope)
  } else {
    support <- functional_support_pair(a[, candidates, drop = FALSE], b[, candidates, drop = FALSE], md)
    qsave(list(input_signature = regression_signature, support = support), regression_path, preset = 'high')
  }
  # Predictions are committed before the exploratory representation and report.
  fwrite(data.table(path = regression_path,
    sha256 = digest::digest(regression_path, file = TRUE, algo = 'sha256'), full_data_screening = TRUE,
    interpretation = 'retrospective cross-fitted per-gene predictions; same stage-15 screening'),
    file.path(destination, 'data/regression_prediction_freeze.tsv'), sep = '\t')
  admission <- regression_admission(support, id, md$recipient[1L], regression_gain_p_threshold, minimum_loading_iqr)
  chosen <- candidates[candidates %chin% admission[selected_in_assay == TRUE, gene]]
  fwrite(admission, file.path(destination, 'data/gene_regression_admission.tsv.gz'), sep = '\t')
  check(length(chosen) >= 3L, paste('Fewer than three genes pass original regression admission for', scope))
  check(identical(rownames(support$loadings), candidates) && all(is.finite(support$loadings)),
    'Functional loading candidate matrix is invalid')
  # Scientific correlations are reported in interaction units, without reporting
  # optimization losses as evidence. They retain the stage-15 retrospective status.
  predictions <- t(support$predictions)
  prediction_r <- as.numeric(cor(md$bliss_pp, predictions))
  gene_correlations <- data.table(depmap_id = id, recipient = md$recipient[1L], gene = candidates,
    selected = candidates %chin% chosen, conditions = nrow(md), independent_pairs = uniqueN(md$pair),
    pearson_r = prediction_r, observed_minimum_pp = min(md$bliss_pp), observed_maximum_pp = max(md$bliss_pp),
    predicted_minimum_pp = apply(predictions, 2L, min), predicted_maximum_pp = apply(predictions, 2L, max),
    interpretation = 'retrospective per-gene joint prediction across combinations within one cell')
  correlation_tables[[id]] <- gene_correlations
  fwrite(gene_correlations, file.path(destination, 'data/interaction_prediction_correlations.tsv.gz'), sep = '\t')
  support$loadings <- support$loadings[chosen, , drop = FALSE]
  support$raw_prediction_gain <- support$raw_prediction_gain[chosen, , drop = FALSE]
  support$predictions <- support$predictions[chosen, , drop = FALSE]
  support$coefficients <- support$coefficients[gene %chin% chosen]
  support$channels <- cbind(pmax(support$loadings, 0), pmax(-support$loadings, 0))
  support$settings$clustering_representation <- 'mean_functional_displacement_correlation'
  support$settings$clustering_algorithm <- 'ward.D2'
  support$settings$rank_selection_rule <- local_rank_selection
  # The original displacement-correlation calculation receives both panels as
  # measured observations. No averaging, products or row transformations are used.
  a_md <- copy(md)
  b_md <- copy(md)
  a_md[, `:=`(observation_id = paste0(combination_id, ':A'), drug = pair, dose = dose_a,
    condition_support = energy_a / (energy_a + null_energy), partner = 'A')]
  b_md[, `:=`(observation_id = paste0(combination_id, ':B'), drug = pair, dose = dose_b,
    condition_support = energy_b / (energy_b + null_energy), partner = 'B')]
  discovery_metadata <- rbind(a_md, b_md)
  x <- rbind(a, b)
  rownames(x) <- discovery_metadata$observation_id
  sensitivity <- discovery_metadata[, .(observation_id, bliss_pp)]
  bank_signature <- digest::digest(list(regression_signature, chosen, support$settings,
    x = x[, chosen, drop = FALSE], null = null[, chosen, drop = FALSE], discovery_metadata), algo = 'sha256')
  bank_path <- file.path(destination, 'models/program_bank.qs')
  if (file.exists(bank_path)) {
    bank <- qread(bank_path)
    check(identical(bank$combination_input_signature, bank_signature), 'Saved combination Ward bank inputs changed')
    bank <- restore_functional_bank(bank, x)
    message('Reused original Ward bank: ', scope)
  } else {
    bank <- discover_local(x, null, discovery_metadata, chosen, 'drug', scope, 'cell', support = support,
      sensitivity_sign = TRUE)
    bank$combination_input_signature <- bank_signature
    save_functional_bank(bank, bank_path)
  }
  render_diagnostics(bank, rank_plots_only = TRUE)
  render_combination(bank, a, b, md, destination)
  inventories[[id]] <- data.table(depmap_id = id, recipient = md$recipient[1L],
    conditions = nrow(md), pairs = uniqueN(md$pair), candidate_genes = length(candidates),
    admitted_genes = length(chosen), selected_programs = bank$diagnostics$selected_rank,
    selection_status = bank$diagnostics$selection_status,
    median_admitted_gene_prediction_r = median(gene_correlations[selected == TRUE, pearson_r]))
  fwrite(rbindlist(inventories), file.path(data_dir, 'combination_program_inventory.tsv'), sep = '\t')
  rm(a, b, x, query_x, molecular, support, bank, predictions)
  invisible(gc())
  message('Combination analysis ready: ', scope)
}
fwrite(rbindlist(correlation_tables), file.path(data_dir, 'interaction_prediction_correlations.tsv.gz'), sep = '\t')
current_hashes <- vapply(original_paths, digest::digest, character(1), file = TRUE, algo = 'sha256')
check(identical(unname(current_hashes), original_hashes$sha256), 'Original stage-15 files changed during adaptation')
write_manifest('21_almanac_two_partner_ward_programs', input_paths, list(settings = settings,
  calculation_signature = calculation_signature, original_files_unchanged = TRUE,
  interpretation = 'same retrospective stage-15 method, joint raw-partner gene regression and Bliss outcome'))
writeLines(c('ALMANAC combination adaptation of the stage-15 Ward method.',
  'Run from /work: Rscript src/program_pilot/21_almanac_two_partner_ward_programs.R',
  'The original script and shared calculations are unchanged.',
  'Two original signed partner displacements are the per-gene ridge predictors.',
  'No sum/product features, program decoder or nested program-bank validation was added.',
  'Partner A/B heatmaps are stacked with identical combination and gene orders.',
  'Functional outcome: unclipped Bliss excess inhibition in percentage points.',
  'Full-data outcome screening and gene admission are retrospective, matching stage 15.',
  'Positive/negative counts describe observed interaction direction, not significance or clinical response.',
  'Prediction correlations compare combinations within each cell; they do not measure cross-cell transport.',
  'Baseline-to-displacement transport is untested; both signatures are measured within the recipient.'),
  file.path(output_dir, 'README.txt'))
capture.output(sessionInfo(), file = file.path(output_dir, 'sessionInfo.txt'))
message('Completed the original-method combination panel: ', length(inventories), ' cells')

