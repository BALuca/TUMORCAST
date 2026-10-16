#!/usr/bin/env Rscript
# Restore the historical joint bridge and symmetric nonlinear combination model.
# The current frozen bank is the only scientific substitution; selected historical
# penalties, recipients, outcomes and outer folds are held fixed.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})
source('/work/src/program_pilot/utilities.R')
source('/work/src/program_pilot/functional_programs.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn = 1)
set.seed(20260929L)

root <- '/work/results/program_pilot/53_historical_joint_and_nonlinear_bridge_current_bank'
bank_path <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs'
historical_bridge <- '/work/almanac_conditionals/state_conditioned_bridge/program_first/bridge'
historical_combo <- '/work/almanac_broad9_2026-10-07'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
stage <- commandArgs(trailingOnly = TRUE)
if (!length(stage)) stage <- 'all'
check(length(stage) == 1L && stage %in% c('all', 'combination', 'joint', 'prepare', 'report'), 'Unknown stage')
for (path in file.path(root, c('data', 'models', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', name), sep = '\t')
status <- function(message_text) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', message_text)
  fwrite(data.table(stage = message_text, time = as.character(Sys.time()), pid = Sys.getpid()),
    file.path(root, paste0('status_', stage, '.tsv')), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(6); quit(save = 'no', status = 1L)})
bank <- qread(bank_path)
genes <- bank$genes
W <- bank$weights[genes, , drop = FALSE]
P <- colnames(W)
K <- ncol(W)
check(K == 34L && length(genes) == 8429L && all(abs(colSums(W) - 1) < 1e-12), 'Current bank changed')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 150)
}
cor_safe <- function(x, y) {
  if (length(x) < 3L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y)
}
# Pair-balanced historical scientific summaries give each pair equal weight.
pair_weights <- function(z) {
  counts <- copy(z[, .(pair, depmap_id)])
  counts[, row_count := .N, by = .(pair, depmap_id)]
  counts[, cell_count := uniqueN(depmap_id), by = pair]
  w <- 1 / (counts$row_count * counts$cell_count)
  w / sum(w)
}
weighted_cor <- function(x, y, w) {
  x <- x - sum(w * x)
  y <- y - sum(w * y)
  denominator <- sqrt(sum(w * x^2) * sum(w * y^2))
  if (denominator < 1e-15) return(NA_real_)
  sum(w * x * y) / denominator
}

if (stage %in% c('all', 'combination')) {
  status('restoring557 historical combination rows, original42 splits and selected penalties')
  historical <- readRDS(file.path(historical_combo, 'cross_context_enrichment/results/models_and_splits.rds'))
  d <- copy(historical$d)
  queries <- fread(file.path(historical_combo, 'results/observed_input_displacements.tsv.gz'))
  settings <- fread(file.path(historical_combo, 'bank40_comparison/results/penalties.tsv'))
  check(nrow(d) == 557L && !anyDuplicated(d$source_row), 'Historical combination rows changed')
  combo_signature <- digest::digest(list(W, d$source_row, historical$splits, settings), algo = 'sha256')
  query_path <- file.path(root, 'data/combination_query_scores.qs')
  if (file.exists(query_path)) {
    cached <- qread(query_path)
    check(identical(cached$signature, combo_signature), 'Stale combination inputs')
    query_scores <- cached$scores
  } else {
    panel <- qread(panel_path)
    query_scores <- matrix(NA_real_, nrow(queries), K, dimnames = list(as.character(queries$query_id), P))
    for (id in unique(queries$depmap_id)) {
      cell_index <- match(id, panel$cell_map$depmap_id)
      check(!is.na(cell_index), 'Historical combination recipient missing from molecular cache')
      molecular <- qread(panel$input_paths[cell_index])
      check(all(genes %in% colnames(molecular$response)), 'Incomplete current-bank gene coverage')
      native <- as.matrix(molecular$response[, genes, drop = FALSE] %*% W)
      md <- copy(molecular$metadata)
      md[, drug := panel_name(drug)]
      for (i in which(queries$depmap_id == id)) {
        lo <- which(md$drug == panel_name(queries$drug[i]) & abs(md$dose - queries$lo[i]) < 1e-8)
        hi <- which(md$drug == panel_name(queries$drug[i]) & abs(md$dose - queries$hi[i]) < 1e-8)
        check(length(lo) == 1L && length(hi) == 1L,
          paste('Historical query cannot be projected:', id, queries$drug[i], queries$lo[i], queries$hi[i]))
        query_scores[i, ] <- (1 - queries$alpha[i]) * native[lo, ] + queries$alpha[i] * native[hi, ]
      }
      status(paste('combination inputs projected:', id))
    }
    check(all(is.finite(query_scores)), 'Incomplete projected combination inputs')
    qsave(list(signature = combo_signature, scores = query_scores, queries = queries), query_path)
  }
  A <- query_scores[match(d$iA, queries$query_id), , drop = FALSE]
  B <- query_scores[match(d$iB, queries$query_id), , drop = FALSE]
  check(all(is.finite(A)) && all(is.finite(B)), 'Combination query join failed')
  raw <- fread(file.path(historical_combo, 'structural_diagnosis/results/raw_eligible_records.tsv.gz'))
  refs <- unique(rbind(raw[, .(source_row, ref = single_groupA)], raw[, .(source_row, ref = single_groupB)]))
  refs <- refs[source_row %in% d$source_row]
  refs[, row := match(source_row, d$source_row)]
  purge <- function(tr, te) setdiff(tr, refs[ref %in% refs[row %in% te, ref], row])
  product_indices <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  products <- vapply(seq_len(nrow(product_indices)), function(j) {
    u <- product_indices[j, 1]
    v <- product_indices[j, 2]
    if (u == v) A[, u] * B[, v] else A[, u] * B[, v] + A[, v] * B[, u]
  }, numeric(nrow(A)))
  distance2 <- function(x, y) pmax(outer(rowSums(x^2), rowSums(y^2), '+') - 2 * tcrossprod(x, y), 0)
  training_weights <- function(z) {
    w <- numeric(nrow(z))
    for (stratum_name in unique(z$stratum)) {
      ix <- which(z$stratum == stratum_name)
      w[ix] <- pair_weights(z[ix]) / uniqueN(z$stratum)
    }
    w
  }
  # Exact historical symmetric kernels and zero-either-input magnitude gate.
  kernel_matrix <- function(a, b, train_a, train_b, bandwidth, kind) {
    amplitude <- sqrt(rowSums(a^2) * rowSums(b^2))
    amplitude <- amplitude / (1 + amplitude)
    training_amplitude <- sqrt(rowSums(train_a^2) * rowSums(train_b^2))
    training_amplitude <- training_amplitude / (1 + training_amplitude)
    if (kind == 'individual') {
      similarity <- (exp(-distance2(a, train_a) / bandwidth) + exp(-distance2(b, train_b) / bandwidth) +
        exp(-distance2(a, train_b) / bandwidth) + exp(-distance2(b, train_a) / bandwidth)) / 4
    } else {
      direct <- distance2(a, train_a) + distance2(b, train_b)
      swapped <- distance2(a, train_b) + distance2(b, train_a)
      similarity <- (exp(-direct / bandwidth) + exp(-swapped / bandwidth)) / 2
    }
    similarity * outer(amplitude, training_amplitude)
  }
  # Replay the earlier nonlinear solver on its old native inputs before using
  # the substituted bank. This checks the model implementation independently.
  old_scores <- fread(file.path(historical_combo, 'bank40_comparison/results/displacements40.tsv.gz'))
  old_programs <- grep('^P[0-9]+$', names(old_scores), value = TRUE)
  old_a <- as.matrix(old_scores[match(d$iA, query_id), ..old_programs])
  old_b <- as.matrix(old_scores[match(d$iB, query_id), ..old_programs])
  old_predictions <- fread(file.path(historical_combo, 'bank40_comparison/results/predictions.tsv.gz'))
  replay_checks <- list()
  for (split_name in names(historical$splits)[seq_len(2L)]) {
    split <- historical$splits[[split_name]]
    tr <- purge(split$tr, split$te)
    te <- split$te
    fold_match <- if (is.na(split$pair_fold)) is.na(settings$pair_fold) else settings$pair_fold == split$pair_fold
    lambda <- settings[task == split$task & held_cell == split$cell & fold_match & model == 'flexible', lambda]
    w <- training_weights(d[tr])
    scale <- pmax(sqrt(colSums((old_a[tr, , drop = FALSE]^2 + old_b[tr, , drop = FALSE]^2) * w) / 2), 1e-12)
    aa <- sweep(old_a, 2L, scale, '/')
    bb <- sweep(old_b, 2L, scale, '/')
    direct <- distance2(aa[tr, ], aa[tr, ]) + distance2(bb[tr, ], bb[tr, ])
    swapped <- distance2(aa[tr, ], bb[tr, ]) + distance2(bb[tr, ], aa[tr, ])
    bandwidth <- max(median(pmin(direct, swapped)[upper.tri(direct)]), 1e-8)
    kernel <- kernel_matrix(aa[tr, ], bb[tr, ], aa[tr, ], bb[tr, ], bandwidth, 'flexible')
    sw <- sqrt(w)
    eig <- eigen(kernel * outer(sw, sw), symmetric = TRUE)
    coefficients <- sw * as.numeric(eig$vectors %*%
      (crossprod(eig$vectors, d$interaction[tr] * sw) / (pmax(eig$values, 0) + lambda)))
    replay <- as.numeric(kernel_matrix(aa[te, ], bb[te, ], aa[tr, ], bb[tr, ], bandwidth, 'flexible') %*% coefficients)
    expected <- old_predictions[task == split$task & held_cell == split$cell & model == 'flexible']
    if (!is.na(split$pair_fold)) expected <- expected[pair_fold == split$pair_fold]
    error <- max(abs(replay - expected$prediction[match(d$source_row[te], expected$source_row)]))
    check(error < 1e-9, 'Historical nonlinear decoder does not replay')
    replay_checks[[length(replay_checks) + 1L]] <- data.table(check = paste0('historical_kernel_', split_name), max_error = error)
  }
  write_table(rbindlist(replay_checks), 'combination_solver_acceptance.tsv')
  for (split_name in names(historical$splits)) {
    split <- historical$splits[[split_name]]
    te <- split$te
    tr <- purge(split$tr, te)
    check(!split$cell %in% d$depmap_id[tr], 'Held combination cell enters fitting')
    if (split$task == 'recipient_pair') check(!length(intersect(d$pair[tr], d$pair[te])), 'Held pairs enter fitting')
    for (kind in c('raw', 'flexible', 'individual')) {
      path <- file.path(root, 'models', paste0('combination_', split_name, '_', kind, '.qs'))
      if (file.exists(path)) {
        fit <- qread(path)
        check(identical(fit$signature, combo_signature), 'Stale combination model')
      } else {
        fold_match <- if (is.na(split$pair_fold)) is.na(settings$pair_fold) else settings$pair_fold == split$pair_fold
        selected <- settings[task == split$task & held_cell == split$cell & fold_match & model == kind]
        check(nrow(selected) == 1L, 'Historical penalty join is ambiguous')
        lambda <- selected$lambda
        w <- training_weights(d[tr])
        sw <- sqrt(w)
        if (kind == 'raw') {
          scale <- pmax(sqrt(colSums(products[tr, , drop = FALSE]^2 * w)), 1e-12)
          z <- sweep(products[tr, , drop = FALSE], 2L, scale, '/')
          t <- z * sw
          eig <- eigen(tcrossprod(t), symmetric = TRUE)
          beta <- as.numeric(crossprod(t, eig$vectors %*%
            (crossprod(eig$vectors, d$interaction[tr] * sw) / (pmax(eig$values, 0) + lambda))))
          prediction <- as.numeric(sweep(products[te, , drop = FALSE], 2L, scale, '/') %*% beta)
          model <- list(scale = scale, beta = beta)
          zero_error <- 0
          symmetry_error <- 0
        } else {
          scale <- pmax(sqrt(colSums((A[tr, , drop = FALSE]^2 + B[tr, , drop = FALSE]^2) * w) / 2), 1e-12)
          a <- sweep(A, 2L, scale, '/')
          b <- sweep(B, 2L, scale, '/')
          direct <- distance2(a[tr, ], a[tr, ]) + distance2(b[tr, ], b[tr, ])
          swapped <- distance2(a[tr, ], b[tr, ]) + distance2(b[tr, ], a[tr, ])
          bandwidth <- max(median(pmin(direct, swapped)[upper.tri(direct)]), 1e-8)
          kernel <- kernel_matrix(a[tr, ], b[tr, ], a[tr, ], b[tr, ], bandwidth, kind)
          eig <- eigen(kernel * outer(sw, sw), symmetric = TRUE)
          coefficients <- sw * as.numeric(eig$vectors %*%
            (crossprod(eig$vectors, d$interaction[tr] * sw) / (pmax(eig$values, 0) + lambda)))
          kt <- kernel_matrix(a[te, ], b[te, ], a[tr, ], b[tr, ], bandwidth, kind)
          prediction <- as.numeric(kt %*% coefficients)
          reverse <- kernel_matrix(b[te, ], a[te, ], a[tr, ], b[tr, ], bandwidth, kind)
          symmetry_error <- max(abs(kt - reverse))
          zero_error <- max(abs(kernel_matrix(a[te, ] * 0, b[te, ], a[tr, ], b[tr, ], bandwidth, kind)))
          check(symmetry_error < 1e-12 && zero_error == 0, 'Historical kernel composition boundary failed')
          model <- list(scale = scale, bandwidth = bandwidth, train_a = a[tr, ], train_b = b[tr, ], coefficients = coefficients)
        }
        fit <- list(signature = combo_signature, model = model, kind = kind, lambda = lambda,
          training_source_rows = d$source_row[tr], held_source_rows = d$source_row[te],
          prediction = prediction, split = split, symmetry_error = symmetry_error, zero_error = zero_error)
        qsave(fit, path)
      }
      z <- copy(d[te])
      z[, `:=`(task = split$task, held_cell = split$cell, pair_fold = split$pair_fold,
        model = kind, bank = 'Current43 / 34 programs', prediction = fit$prediction)]
      fwrite(z[, .(source_row, task, held_cell, pair_fold, model, prediction)],
        file.path(root, 'data', paste0('combination_', split_name, '_', kind, '_frozen.tsv')), sep = '\t')
      write_table(z, paste0('combination_', split_name, '_', kind, '.tsv'))
      status(paste('MODEL AND PREDICTIONS SAVED:', split_name, kind))
    }
  }
  paths <- list.files(file.path(root, 'data'), '^combination_.*_(raw|flexible|individual)\\.tsv$', full.names = TRUE)
  current <- rbindlist(lapply(paths, fread))
  previous <- fread(file.path(historical_combo, 'bank40_comparison/results/predictions.tsv.gz'))
  previous[, bank := 'Historical Tahoe/CytoSig / 40 programs']
  combined <- rbindlist(list(previous, current), fill = TRUE)
  write_table(combined, 'combination_predictions.tsv.gz')
  metrics <- combined[stratum == 'enriched', {
    w <- pair_weights(.SD)
    .(measurements = .N, cells = uniqueN(depmap_id), pairs = uniqueN(pair),
      interaction_r = weighted_cor(prediction, interaction, w))
  }, by = .(bank, task, model)]
  by_cell <- combined[stratum == 'enriched', .(measurements = .N, pairs = uniqueN(pair),
    interaction_r = weighted_cor(prediction, interaction, pair_weights(.SD))), by = .(bank, task, held_cell, model)]
  write_table(metrics, 'combination_summary.tsv')
  write_table(by_cell, 'combination_by_cell.tsv')
  z <- combined[stratum == 'enriched' & task == 'recipient_pair' & held_cell == 'ACH-000552']
  z[, bank := factor(bank, levels = c('Historical Tahoe/CytoSig / 40 programs', 'Current43 / 34 programs'))]
  annotations <- by_cell[task == 'recipient_pair' & held_cell == 'ACH-000552']
  annotations[, bank := factor(bank, levels = levels(z$bank))]
  annotations[, label := sprintf('Pair-balanced r = %.2f', interaction_r)]
  plot <- ggplot(z, aes(interaction, prediction)) + geom_point(size = 1, alpha = .65, color = '#24587b') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = '#ac4035') +
    geom_text(data = annotations, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.03, vjust = 1.3, size = 3) + facet_grid(bank ~ model, scales = 'free') +
    labs(x = 'Measured combination interaction (log2)', y = 'Predicted interaction (log2)',
      title = 'Same 119 HT29 measurements; cell and evaluated pairs held out',
      subtitle = 'Original folds, reference purging and selected penalties; frozen bank substituted') +
    theme_classic(base_size = 10) + theme(aspect.ratio = 1, strip.background = element_blank())
  save_plot(plot, 'historical_current_combination_comparison', 10, 6.7)
  names_by_id <- setNames(d$CELLNAME, d$depmap_id)
  all_cells <- combined[stratum == 'enriched' & task == 'recipient_pair' & model == 'flexible']
  all_cells[, `:=`(cell = unname(names_by_id[held_cell]),
    bank = factor(bank, levels = levels(z$bank)))]
  annotations <- copy(by_cell[task == 'recipient_pair' & model == 'flexible'])
  annotations[, `:=`(cell = unname(names_by_id[held_cell]), bank = factor(bank, levels = levels(z$bank)),
    label = sprintf('r = %.2f; n = %d', interaction_r, measurements))]
  plot <- ggplot(all_cells, aes(interaction, prediction)) +
    geom_point(size = .9, alpha = .6, color = '#24587b') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = '#ac4035') +
    geom_text(data = annotations, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.02, vjust = 1.3, size = 2.7) + facet_grid(bank ~ cell, scales = 'free') +
    labs(x = 'Measured combination interaction (log2)', y = 'Predicted interaction (log2)',
      title = 'All seven recipients: cell and evaluated pairs held out',
      subtitle = 'Original signal-enriched cohort; symmetric nonlinear joint-profile decoder') +
    theme_classic(base_size = 9) + theme(aspect.ratio = 1, strip.background = element_blank())
  save_plot(plot, 'historical_current_combinations_all_cells', 15, 5.1)
  status('combination restoration complete')
  print(metrics)
  if (stage == 'combination') quit(save = 'no', status = 0L)
}

if (stage %in% c('all', 'prepare', 'joint', 'report')) {
  joint_source <- file.path(historical_bridge, 'joint40/linear_decoder/artifacts/inputs.rds')
  original <- readRDS(joint_source)
  joint_signature <- digest::digest(list(W, original$M, original$F, original$folds,
    original$receptors, 'historical-selected-feedback-and-ridge'), algo = 'sha256')
  input_path <- file.path(root, 'data/joint_inputs.qs')
  if (file.exists(input_path)) {
    prepared <- qread(input_path)
    check(identical(prepared$signature, joint_signature), 'Stale joint inputs')
    a <- prepared$inputs
  } else {
    status('projecting all8429 genes onto the unchanged270 drug and37 ligand rows')
    suppressPackageStartupMessages({library(jsonlite); library(arrow)})
    g <- readRDS('/work/almanac_conditionals/state_conditioned_bridge/sensitivity_selected/artifacts/eligible_gene_measurements.rds')
    refs <- readRDS(file.path(historical_bridge, 'artifacts/reference_genes.rds'))
    common <- intersect(genes, colnames(g$response))
    additional <- setdiff(genes, colnames(g$response))
    write_table(data.table(gene = genes, archived_gene_input = genes %in% common), 'gene_input_coverage.tsv')
    recovery_path <- file.path(root, 'data/recovered_gene_inputs.qs')
    if (file.exists(recovery_path)) {
      recovery <- qread(recovery_path)
      check(identical(recovery$signature, joint_signature), 'Stale recovered gene inputs')
      response <- recovery$response
      C <- recovery$C
      PC <- recovery$PC
      md <- recovery$metadata
    } else {
    # The archived gene cache lacks99 current-bank genes. Recover them from the
    # primary R-readable count files; never drop genes or alter program weights.
    expressions <- parse(file.path(historical_bridge, 'scripts/10_prepare_bridge.R'))
    helper_index <- which(vapply(expressions, function(x)
      is.call(x) && as.character(x[[1]]) %in% c('=', '<-') && identical(as.character(x[[2]]), 'read_npy_rows'), logical(1)))
    check(length(helper_index) == 1L, 'Historical R matrix-row reader missing')
    eval(expressions[[helper_index]], envir = environment())
    primary_dir <- '/work/cache/tahoe_prism_full_chain/tahoe_pdex_selected'
    gene_table <- fread(file.path(primary_dir, 'genes.tsv'))
    extra_index <- match(additional, gene_table$feature)
    check(!anyNA(extra_index), 'Current bank genes missing in primary molecular source')
    reference_counts <- read_npy_matrix(file.path(primary_dir, 'pdex_plate_reference.npy'))
    reference_metadata <- fread(file.path(primary_dir, 'pdex_plate_reference_metadata.tsv'))
    original_rows <- match(paste(refs$cm$depmap_id, refs$cm$plate),
      paste(reference_metadata$depmap_id, reference_metadata$plate))
    check(!anyNA(original_rows), 'Original reference plates changed')
    C <- matrix(NA_real_, nrow(refs$C), length(genes), dimnames = list(NULL, genes))
    C[, common] <- refs$C[, common, drop = FALSE]
    C[, additional] <- reference_counts[original_rows, extra_index, drop = FALSE]
    PC <- matrix(NA_real_, nrow(refs$PC), length(genes), dimnames = list(NULL, genes))
    PC[, common] <- refs$PC[, common, drop = FALSE]
    response <- matrix(NA_real_, nrow(g$response), length(genes), dimnames = list(g$metadata$sample_id, genes))
    response[, common] <- g$response[, common, drop = FALSE]
    md <- copy(g$metadata)
    metadata <- fread('/work/almanac_conditionals/state_conditioned_bridge/program_first/results/full_tahoe_profile_metadata.tsv.gz')
    primary <- md[source == 'original cache']
    primary[, original_row := match(sample_id, md$sample_id)]
    key <- function(x) paste(x$depmap_id, canonical_drug(x$drug), sprintf('%.12g', x$dose))
    source_index <- match(key(primary), key(metadata))
    check(!anyNA(source_index), 'Original drug rows absent from primary metadata')
    primary_metadata <- metadata[source_index]
    for (shard in unique(primary_metadata$shard_path)) {
      ix <- which(primary_metadata$shard_path == shard)
      delta <- read_npy_rows(file.path('/work/cache/tahoe_prism_full_chain', shard), primary_metadata$shard_row[ix])
      for (j in seq_along(ix)) {
        r <- primary_metadata[ix[j]]
        plates <- trimws(strsplit(r$contributing_plates, ',', fixed = TRUE)[[1]])
        ri <- match(paste(r$depmap_id, plates), paste(reference_metadata$depmap_id, reference_metadata$plate))
        check(!anyNA(ri), 'Original treated reference plates missing')
        ref <- colMeans(reference_counts[ri, extra_index, drop = FALSE])
        treated <- pmax(delta[j, extra_index] + ref, 0)
        response[primary$original_row[ix[j]], additional] <-
          log2(1 + treated / r$treated_total * 1e4 / r$treated_correction_factor) -
          log2(1 + ref / r$reference_total * 1e4)
      }
    }
    rm(delta, reference_counts)
    status('original count genes restored; recovering public count genes and reference plates')
    mapping <- unique(fread('/work/almanac_conditionals/diagnosis/public_tahoe_cell_metadata.tsv')[,
      .(cell_line = Cell_ID_Cellosaur, depmap_id = Cell_ID_DepMap)])
    public_keys <- paste(refs$pcm$depmap_id, refs$pcm$plate)
    remaining <- unique(public_keys)
    response_parts <- list()
    response_metadata <- list()
    public_wanted <- key(md[source == 'public recovered'])
    focal_files <- readLines('/work/almanac_conditionals/state_conditioned_bridge/mechanism_diagnosis/focal_parquet_chunks.txt')
    files <- unique(c(focal_files, list.files('/work/almanac_conditionals/diagnosis/pdex_projection_chunks',
      pattern = 'json$', full.names = TRUE)))
    reader <- arrow::ParquetFileReader$create('/work/almanac_conditionals/diagnosis/pdex_selected_sparse.parquet')
    for (file in files) {
      if (!length(remaining) && length(response_parts)) {
        recovered <- unique(key(rbindlist(response_metadata)))
        if (all(public_wanted %in% recovered)) break
      }
      description <- jsonlite::fromJSON(file)
      if (!length(description$profiles)) next
      profiles <- as.data.table(description$profiles)
      profiles[, depmap_id := mapping$depmap_id[match(cell_line, mapping$cell_line)]]
      profiles[, plate_key := paste(depmap_id, plate)]
      profiles <- profiles[plate_key %in% remaining | key(profiles) %in% public_wanted]
      if (!nrow(profiles)) next
      tab <- as.data.table(reader$ReadRowGroup(as.integer(description$rowgroup), c(0L, 1L, 2L, 3L, 11L, 12L)))
      for (j in seq_len(nrow(profiles))) {
        r <- profiles[j]
        profile_key <- paste(r$depmap_id, r$plate, r$drug, r$dose)
        if (!r$plate_key %in% remaining && profile_key %in% names(response_parts)) next
        v <- tab[cell_line == r$cell_line & plate == r$plate & grepl(paste0("('", r$drug, "',"), target, fixed = TRUE)]
        targets <- unique(v$target)
        if (length(targets) > 1L) {
          doses <- as.numeric(sub(".*', ([0-9.eE+-]+),.*", '\\1', targets))
          check(any(abs(doses - r$dose) < 1e-8), 'Public target dose mismatch')
          v <- v[target == targets[which.min(abs(doses - r$dose))]]
        }
        check(nrow(v) == nrow(gene_table) && !anyDuplicated(v$feature), 'Incomplete public count profile')
        gi <- match(additional, v$feature)
        check(!anyNA(gi), 'Missing public genes')
        if (r$plate_key %in% remaining) {
          pi <- which(public_keys == r$plate_key)
          PC[pi, additional] <- rep(v$ref_mean[gi], each = length(pi))
          check(abs(sum(v$ref_mean) - refs$pcm$total[pi[1]]) < 1e-3, 'Public reference totals changed')
          remaining <- setdiff(remaining, r$plate_key)
        }
        if (key(r) %in% public_wanted && !profile_key %in% names(response_parts)) {
          dx <- log2(1 + v$target_mean[gi] / r$target_total * 1e4 / r$factor) -
            log2(1 + v$ref_mean[gi] / r$ref_total * 1e4)
          response_parts[[profile_key]] <- dx
          response_metadata[[profile_key]] <- r[, .(depmap_id, drug, dose, plate)]
        }
      }
      if (length(response_parts) %% 20L == 0L || !length(remaining)) {
        status(paste('public genes:', length(response_parts), 'drug plates;', length(remaining), 'reference plates remain'))
      }
    }
    check(!length(remaining), 'Public baseline reference coverage incomplete')
    public_md <- rbindlist(response_metadata)
    public_response <- do.call(rbind, response_parts)
    public_md[, response_row := .I]
    public_md[, condition := key(public_md)]
    for (i in which(md$source == 'public recovered')) {
      ix <- public_md[condition == key(md[i]), response_row]
      check(length(ix) > 0, 'Public drug condition absent')
      plates <- sort(public_md$plate[ix])
      required_plates <- sort(as.integer(strsplit(md$contributing_plates[i], ',', fixed = TRUE)[[1]]))
      check(identical(plates, required_plates), 'Public drug plate aggregation changed')
      response[i, additional] <- colMeans(public_response[ix, , drop = FALSE])
    }
    check(all(is.finite(response)) && all(is.finite(C)) && all(is.finite(PC)), 'Gene reconstruction incomplete')
    qsave(list(response = response, C = C, PC = PC, metadata = md, signature = joint_signature),
      file.path(root, 'data/recovered_gene_inputs.qs'))
    }
    # Apply the original baseline rank convention to the new native programs.
    feature <- function(v) {
      ranks <- (rank(v, na.last = 'keep') - .5) / sum(is.finite(v)) - .5
      ranks[!is.finite(ranks)] <- 0
      as.numeric(ranks %*% W) / colSums(abs(W))
    }
    expression_path <- '/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv'
    expression <- fread(expression_path)
    symbols <- sub(' \\(ENSG.*', '', names(expression)[-1L])
    check(!anyDuplicated(symbols[symbols %in% genes]) && all(genes %in% symbols), 'CCLE gene coverage changed')
    ids <- expression[[1L]]
    x <- as.matrix(expression[, match(genes, symbols) + 1L, with = FALSE])
    rownames(x) <- ids
    totals <- rowSums(pmax(2^as.matrix(expression[, -1L]) - 1, 0))
    names(totals) <- ids
    CP <- sweep(pmax(2^x - 1, 0), 1L, totals[ids], '/') * 1e4
    colnames(CP) <- genes
    baseline <- t(apply(x, 1L, feature))
    rownames(baseline) <- ids
    colnames(baseline) <- P
    m <- original$M[source_type == 'drug']
    check(identical(m$sample_id, md$sample_id), 'Historical molecular row order changed')
    xm <- matrix(NA_real_, nrow(m), K, dimnames = list(m$sample_id, P))
    for (i in seq_len(nrow(m))) {
      r <- md[i]
      control_md <- if (r$source == 'original cache') refs$cm else refs$pcm
      control <- if (r$source == 'original cache') C else PC
      exclude <- strsplit(r$contributing_plates, ',', fixed = TRUE)[[1]]
      ix <- which(control_md$depmap_id == r$depmap_id & !as.character(control_md$plate) %in% exclude)
      check(length(ix) > 0L, 'Independent baseline has no reference plates')
      v <- log2(1 + colSums(control[ix, , drop = FALSE]) / sum(control_md$total[ix]) * 1e4)
      xm[i, ] <- feature(v)
    }
    for (id in unique(m$depmap_id)) baseline[id, ] <- colMeans(xm[m$depmap_id == id, , drop = FALSE])
    # Reuse the historical R expression readers for the same fixed37 ligand pairs.
    # Only their gene projection changes; no new study or contrast is admitted.
    ligand_dir <- file.path(historical_bridge, 'ligand_conditioning_trial')
    ligand_env <- new.env(parent = globalenv())
    ligand_env$G <- file.path(ligand_dir, 'inputs/geo')
    ligand_env$sm <- fread(file.path(ligand_dir, 'inputs/sample_metadata.tsv'))
    ligand_env$genes <- genes
    ligand_env$mapcache <- list()
    ligand_env$exprcache <- list()
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    ligand_env$len <- fread(file.path(ligand_env$G, 'GSE119606_HSCs_vehicle_PDGF_TGF_raw.tsv.gz'))[,
      .(gene = GeneName, len = CodingLength)][is.finite(len) & len > 0, .(len = median(len)), by = gene]
    ligand_env$ensmap <- AnnotationDbi::mapIds(org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db, 'ENSEMBL'), keytype = 'ENSEMBL', column = 'SYMBOL', multiVals = 'asNA')
    ligand_expressions <- parse(file.path(ligand_dir, 'scripts/04_prepare.R'))
    for (name in c('readtab', 'getmap', 'getexpr')) {
      ix <- which(vapply(ligand_expressions, function(e)
        is.call(e) && as.character(e[[1]]) %in% c('=', '<-') && identical(as.character(e[[2]]), name), logical(1)))
      check(length(ix) == 1L, 'Historical ligand reader missing')
      eval(ligand_expressions[[ix]], envir = ligand_env)
    }
    pairs <- fread(file.path(ligand_dir, 'inputs/curated_pairs.tsv'))
    sample_metadata <- ligand_env$sm
    cyto_profiles <- unique(unlist(strsplit(pairs[pair_id %in% original$ligand_meta$pair_id, profiles], ';', fixed = TRUE)))
    cyto <- fread('/work/data/01_CytoSig_data/diff.merge.gz', select = c('V1', cyto_profiles))
    setnames(cyto, 1L, 'gene')
    records <- list()
    xl <- matrix(NA_real_, nrow(original$ligand_meta), K, dimnames = list(original$ligand_meta$pair_id, P))
    for (i in seq_len(nrow(original$ligand_meta))) {
      pair_id_value <- original$ligand_meta$pair_id[i]
      r <- pairs[pair_id == pair_id_value]
      check(nrow(r) == 1L, 'Historical ligand pair is missing')
      ex <- ligand_env$getexpr(r)
      covered <- intersect(genes, rownames(ex$X))
      controls <- strsplit(r$controls, ';', fixed = TRUE)[[1]]
      check(all(controls %in% colnames(ex$X)), 'Historical ligand controls missing')
      baseline_gene <- setNames(rep(NA_real_, length(genes)), genes)
      baseline_gene[covered] <- rowMeans(ex$X[covered, controls, drop = FALSE])
      profiles <- strsplit(r$profiles, ';', fixed = TRUE)[[1]]
      fc <- rowMeans(as.matrix(cyto[, ..profiles]), na.rm = TRUE)[match(genes, cyto$gene)]
      names(fc) <- genes
      records[[pair_id_value]] <- list(baseline = baseline_gene, fc = fc,
        agreement = original$ligand_records[[pair_id_value]]$agreement,
        coverage = length(covered) / length(genes), mode = ex$mode)
      xl[i, ] <- feature(baseline_gene)
    }
    # Keep receptor covariates byte-identical; these are a separate upstream block.
    a <- original
    a$W <- W
    a$genes <- genes
    a$CP <- CP
    a$XM <- cbind(rbind(xm, xl), original$XM[, original$receptors, drop = FALSE])
    a$XF <- cbind(baseline[original$F$depmap_id, , drop = FALSE], original$XF[, original$receptors, drop = FALSE])
    a$baseline_features <- cbind(baseline[rownames(original$baseline_features), , drop = FALSE],
      original$baseline_features[, original$receptors, drop = FALSE])
    a$features <- c(P, original$receptors)
    colnames(a$XM) <- colnames(a$XF) <- colnames(a$baseline_features) <- a$features
    a$Ydrug <- as.matrix(response %*% W)
    a$ligand_records <- records
    a$MX <- matrix(NA_real_, nrow(a$F), K, dimnames = list(NULL, P))
    for (drug_name in unique(a$F$drug)) {
      for (id in intersect(unique(a$F$depmap_id), m$depmap_id)) {
        ix <- which(m$intervention == drug_name & m$depmap_id == id)
        if (!length(ix)) next
        ix <- ix[order(m$dose[ix])]
        fi <- which(a$F$drug == drug_name & a$F$depmap_id == id)
        for (j in seq_len(K)) a$MX[fi, j] <- approx(log(m$dose[ix]), a$Ydrug[ix, j], log(a$F$dose[fi]), rule = 1)$y
      }
    }
    a$matched <- complete.cases(a$MX)
    check(identical(a$M, original$M) && identical(a$F, original$F) && identical(a$folds, original$folds),
      'Historical recipients, outcomes or folds changed')
    check(identical(a$matched, original$matched), 'Paired measured-function anchor membership changed')
    check(all(is.finite(a$XM)) && all(is.finite(a$XF)) && all(is.finite(a$Ydrug)), 'Nonfinite joint inputs')
    qsave(list(signature = joint_signature, inputs = a), input_path)
    write_table(a$M, 'joint_molecular_rows.tsv')
    write_table(a$F, 'joint_functional_rows.tsv')
    write_table(data.table(ligand = a$ligand_meta$ligand, pair_id = names(records),
      coverage = vapply(records, function(r) r$coverage, numeric(1)),
      missing_fc = vapply(records, function(r) sum(!is.finite(r$fc)), integer(1))), 'joint_ligand_coverage.tsv')
    status('complete historical joint inputs saved')
  }
  if (stage == 'prepare') quit(save = 'no', status = 0L)
}

if (stage %in% c('all', 'joint', 'report')) {
  # Reuse the original numerical objective and its linear-head correction.
  # All functions use this explicit environment; historical files stay unchanged.
  legacy <- new.env(parent = globalenv())
  sys.source(file.path(historical_bridge, 'joint40/scripts/02_lib.R'), envir = legacy)
  linear_expressions <- parse(file.path(historical_bridge, 'joint40/linear_decoder/scripts/02_lib.R'))
  for (name in c('unpack_head', 'head', 'head_grad', 'head_features')) {
    ix <- which(vapply(linear_expressions, function(e)
      is.call(e) && as.character(e[[1]]) %in% c('=', '<-') && identical(as.character(e[[2]]), name), logical(1)))
    check(length(ix) == 1L, 'Historical linear head definition missing')
    eval(linear_expressions[[ix]], envir = legacy)
  }
  legacy$fit_quadratic_implementation <- legacy$fit_model
  legacy$IJ <- matrix(integer(), 0L, 2L)
  legacy$diagq <- logical()
  checks <- list()
  original_rows <- list()
  for (fold in seq_len(5L)) {
    saved <- readRDS(file.path(historical_bridge, 'joint40/linear_decoder/artifacts', paste0('outer', fold, '.rds')))
    fi <- which(original$F$block %in% saved$model$excluded)
    rows <- copy(original$F[fi])
    rows[, intervention := drug]
    delta <- legacy$predict_displacement(saved$model, original$XF[fi, , drop = FALSE], rows)
    replay <- legacy$predict_viability(saved$model, delta)
    expected <- saved$fun[model == 'joint', predicted]
    error <- max(abs(replay - expected))
    check(error < 1e-9, 'Historical joint prediction replay failed')
    checks[[length(checks) + 1L]] <- data.table(check = paste0('historical_outer', fold, '_replay'), max_error = error)
    prior <- copy(saved$fun)
    prior[, bank := 'Historical Tahoe/CytoSig / 40 programs']
    original_rows[[fold]] <- prior
  }
  legacy$a <- a
  legacy$K <- K
  legacy$IJ <- matrix(integer(), 0L, 2L)
  legacy$diagq <- logical()
  legacy$INTERVENTIONS <- sort(unique(a$M$intervention))
  legacy$knots <- setNames(lapply(legacy$INTERVENTIONS,
    function(i) sort(unique(a$M[intervention == i, dose]))), legacy$INTERVENTIONS)
  legacy$.ctx <- new.env(parent = emptyenv())
  legacy$.warm <- new.env(parent = emptyenv())
  # The saved linear wrapper hard-codes40 in an assertion. Call its underlying
  # unchanged optimizer with the same linear functions, now with34 native inputs.
  fit_joint <- function(context, gamma) {
    f <- legacy$fit_quadratic_implementation(context, gamma)
    f$decoder$kind <- 'full_program_linear'
    check(length(f$decoder$coefficients) == K, 'Functional head is not linear in native programs')
    check(f$convergence == 0L, 'Joint optimization did not converge')
    check(all(is.finite(f$b)) && all(is.finite(f$decoder$coefficients)), 'Nonfinite joint model')
    check(abs(legacy$predict_viability(f, matrix(0, 1L, K))) < 1e-12, 'Functional vehicle anchor lost')
    f
  }
  score_joint <- function(f, held, label, fold) {
    fi <- which(a$F$block %in% held)
    rows <- copy(a$F[fi])
    rows[, intervention := drug]
    delta <- legacy$predict_displacement(f, a$XF[fi, , drop = FALSE], rows)
    prediction <- legacy$predict_viability(f, delta)
    frozen <- rows[, .(depmap_id, drug, dose, block)]
    frozen[, `:=`(predicted = prediction, model = label, fold = fold)]
    write_table(frozen, paste0('joint_outer', fold, '_', label, '_frozen.tsv'))
    rows[, `:=`(predicted = prediction, model = label, fold = fold,
      bank = 'Current43 / 34 programs', measured_prediction = NA_real_)]
    paired <- a$matched[fi]
    rows$measured_prediction[paired] <- legacy$predict_viability(f, a$MX[fi[paired], , drop = FALSE])
    mu <- a$F[!block %in% held, .(mean_y = mean(y)), by = .(drug, dose)]
    rows[, mean_prediction := mu$mean_y[match(paste(drug, dose), paste(mu$drug, mu$dose))]]
    mi <- which(a$M$block %in% held)
    predicted_molecular <- legacy$predict_displacement(f, a$XM[mi, , drop = FALSE], a$M[mi])
    write_table(rows, paste0('joint_outer', fold, '_', label, '_functional.tsv'))
    qsave(list(metadata = a$M[mi], prediction = predicted_molecular,
      observed = legacy$context(held)$Y[mi, , drop = FALSE], signature = joint_signature),
      file.path(root, 'data', paste0('joint_outer', fold, '_', label, '_molecular.qs')))
    rows
  }
  if (stage != 'report') {
    status('checking joint objective gradients and starting five historical outer folds')
    first <- readRDS(file.path(historical_bridge, 'joint40/linear_decoder/artifacts/outer1.rds'))
    c <- legacy$context(first$model$excluded)
    warm <- legacy$warm_head(c)
    theta <- c(as.numeric(c$b0), warm$v)
    gamma <- first$model$gamma
    analytic <- legacy$objective(theta, c, gamma)$gradient
    set.seed(20261011L)
    ii <- unique(c(sample(seq_along(theta), 12L), tail(seq_along(theta), K)[seq_len(8L)]))
    finite <- vapply(ii, function(i) {
      plus <- minus <- theta
      plus[i] <- plus[i] + 1e-5
      minus[i] <- minus[i] - 1e-5
      (legacy$objective(plus, c, gamma)$value - legacy$objective(minus, c, gamma)$value) / 2e-5
    }, numeric(1))
    gradient_error <- max(abs(analytic[ii] - finite))
    check(gradient_error < 1e-5, 'Reused joint objective gradient failed')
    checks[[length(checks) + 1L]] <- data.table(check = 'current_joint_gradient', max_error = gradient_error)
    write_table(rbindlist(checks), 'acceptance.tsv')
    fold_paths <- parallel::mclapply(seq_len(5L), function(fold) {
      historical_fit <- readRDS(file.path(historical_bridge, 'joint40/linear_decoder/artifacts', paste0('outer', fold, '.rds')))
      held <- historical_fit$model$excluded
      context <- legacy$context(held)
      check(!any(a$F$depmap_id[context$fi] %in% a$F$depmap_id[a$F$block %in% held]), 'Functional held cells enter fitting')
      for (label in c('joint', 'staged')) {
        gamma <- if (label == 'joint') historical_fit$model$gamma else 0
        path <- file.path(root, 'models', paste0('joint_outer', fold, '_', label, '.qs'))
        status(paste('joint fold', fold, label, 'gamma', gamma, 'started'))
        if (file.exists(path)) {
          fitted <- qread(path)
          check(identical(fitted$signature, joint_signature), 'Stale joint model checkpoint')
          f <- fitted$model
        } else {
          f <- fit_joint(context, gamma)
          qsave(list(signature = joint_signature, model = f, held_blocks = held,
            fitting_functional_rows = context$fi, fitting_molecular_rows = context$mi), path)
        }
        rows <- score_joint(f, held, label, fold)
        status(paste('MODEL AND PREDICTIONS SAVED: joint fold', fold, label))
        ann <- rows[, .(r = cor_safe(100 * 2^y, 100 * 2^predicted)), by = drug]
        ann[, text := sprintf('r = %.2f', r)]
        plot <- ggplot(rows, aes(100 * 2^y, 100 * 2^predicted)) +
          geom_point(size = .8, alpha = .55, color = '#24587b') +
          geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#ac4035', linewidth = .5) +
          geom_text(data = ann, aes(x = -Inf, y = Inf, label = text), inherit.aes = FALSE,
            hjust = -.03, vjust = 1.3, size = 3) + facet_wrap(~drug, scales = 'free') +
          labs(x = 'Measured PRISM remaining abundance (%)', y = 'Predicted remaining abundance (%)',
            title = paste('Current bank;', label, 'model; held-out fold', fold)) +
          theme_classic(base_size = 10) + theme(aspect.ratio = 1)
        save_plot(plot, paste0('joint_outer', fold, '_', label), 7.1, 3.8)
      }
      file.path(root, 'models', paste0('joint_outer', fold, '_joint.qs'))
    }, mc.cores = 3L, mc.preschedule = FALSE)
    check(all(vapply(fold_paths, function(path) is.character(path) && length(path) == 1L && file.exists(path), logical(1))),
      'A historical outer fold failed')
    historical_full <- readRDS(file.path(historical_bridge, 'joint40/linear_decoder/artifacts/frozen_model.rds'))
    full_context <- legacy$context(character())
    for (label in c('joint', 'staged')) {
      path <- file.path(root, 'models', paste0('joint_full_', label, '.qs'))
      if (file.exists(path)) {
        fitted <- qread(path)
        check(identical(fitted$signature, joint_signature), 'Stale full joint model')
        f <- fitted$model
      } else {
        gamma <- if (label == 'joint') historical_full$gamma else 0
        f <- fit_joint(full_context, gamma)
        qsave(list(signature = joint_signature, model = f), path)
      }
      rows <- copy(a$F)
      rows[, intervention := drug]
      rows[, predicted := legacy$predict_viability(f, legacy$predict_displacement(f, a$XF, rows))]
      rows[, `:=`(model = label, bank = 'Current43 / 34 programs', evaluation = 'Own fitting data')]
      write_table(rows, paste0('joint_full_', label, '_functional.tsv'))
      status(paste('full', label, 'model saved'))
    }
  }
  # Summarize repeated dose rows separately from unique recipient outcomes.
  paths <- list.files(file.path(root, 'data'), '^joint_outer[1-5]_(joint|staged)_functional\\.tsv$', full.names = TRUE)
  current <- rbindlist(lapply(paths, fread))
  previous <- rbindlist(original_rows)
  previous[, model := fifelse(model == 'staged40', 'staged', model)]
  combined <- rbindlist(list(previous, current), fill = TRUE)
  write_table(combined, 'joint_held_functional_predictions.tsv.gz')
  metrics <- combined[, {
    sensitive <- y <= -1
    .(rows = .N, cells = uniqueN(depmap_id), pearson_log2 = cor_safe(y, predicted),
      sensitive_rows = sum(sensitive), recovered_sensitive_rows = sum(sensitive & predicted <= -1),
      resistant_rows = sum(!sensitive), false_sensitive_rows = sum(!sensitive & predicted <= -1),
      measured_sensitive_median_pct = median(100 * 2^y[sensitive]),
      predicted_sensitive_median_pct = median(100 * 2^predicted[sensitive]))
  }, by = .(bank, model, drug)]
  write_table(metrics, 'joint_held_functional_summary.tsv')
  by_dose <- combined[, .(cells = uniqueN(depmap_id), pearson_log2 = cor_safe(y, predicted),
    sensitive = sum(y <= -1), recovered = sum(y <= -1 & predicted <= -1),
    resistant = sum(y > -1), false_sensitive = sum(y > -1 & predicted <= -1)), by = .(bank, model, drug, dose)]
  write_table(by_dose, 'joint_held_functional_by_dose.tsv')
  # Dose averages are a secondary overview, not a replacement for independent
  # dose outcomes; the primary per-dose figure and source table retain all doses.
  consensus <- combined[, .(observed = mean(y), predicted = mean(predicted)), by = .(bank, model, drug, depmap_id)]
  write_table(consensus, 'joint_unique_recipient_predictions.tsv.gz')
  summary <- consensus[, .(cells = .N, r = cor_safe(observed, predicted), sensitive = sum(observed <= -1),
    recovered = sum(observed <= -1 & predicted <= -1), resistant = sum(observed > -1),
    false_sensitive = sum(observed > -1 & predicted <= -1)), by = .(bank, model, drug)]
  write_table(summary, 'joint_unique_recipient_summary.tsv')
  summary[, label := sprintf('r = %.2f; sensitive recovered %d/%d', r, recovered, sensitive)]
  plot <- ggplot(consensus, aes(observed, predicted)) + geom_point(size = .85, alpha = .6, color = '#24587b') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = '#ac4035') +
    geom_text(data = summary, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.02, vjust = 1.3, size = 2.8) + facet_grid(bank + model ~ drug, scales = 'free') +
    labs(x = 'Measured PRISM effect (log2 remaining fraction)', y = 'Predicted effect (log2 remaining fraction)',
      title = 'Historical joint architecture: original bank versus current bank',
      subtitle = 'Same recipients and whole-recipient/study folds; dose-mean predictions give one point per cell') +
    theme_classic(base_size = 10) + theme(aspect.ratio = 1, strip.background = element_blank())
  save_plot(plot, 'historical_current_joint_bridge', 8.3, 12)
  by_dose[, panel := paste(drug, dose, 'uM')]
  combined[, panel := paste(drug, dose, 'uM')]
  by_dose[, label := sprintf('r = %.2f; recovered %d/%d', pearson_log2, recovered, sensitive)]
  plot <- ggplot(combined, aes(y, predicted)) + geom_point(size = .75, alpha = .55, color = '#24587b') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = '#ac4035') +
    geom_text(data = by_dose, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.02, vjust = 1.3, size = 2.6) + facet_grid(bank + model ~ panel, scales = 'free') +
    labs(x = 'Measured PRISM effect (log2 remaining fraction)', y = 'Predicted effect (log2 remaining fraction)',
      title = 'Across cell lines receiving the same drug and dose',
      subtitle = '475 held-out recipients per panel; original five recipient/study folds and functional endpoints') +
    theme_classic(base_size = 9) + theme(aspect.ratio = 1, strip.background = element_blank(),
      strip.text = element_text(size = 8))
  save_plot(plot, 'historical_current_joint_bridge_by_dose', 16, 10.2)

  # Molecular variation and measured-input decoding localize the two links.
  measured_input <- current[is.finite(measured_prediction)]
  measured_summary <- measured_input[, .(cells = uniqueN(depmap_id), rows = .N,
    predicted_response_r = cor_safe(y, predicted), measured_response_r = cor_safe(y, measured_prediction)),
    by = .(model, drug, dose)]
  write_table(measured_summary, 'joint_paired_function_stage.tsv')
  molecular_rows <- list()
  for (fold in seq_len(5L)) {
    for (label in c('joint', 'staged')) {
      x <- qread(file.path(root, 'data', paste0('joint_outer', fold, '_', label, '_molecular.qs')))
      keep <- which(x$metadata$source_type == 'drug')
      z <- copy(x$metadata[keep])
      z[, signature_row := .I]
      z <- z[rep(seq_len(nrow(z)), each = K)]
      z[, `:=`(program = rep(P, length(keep)), observed = as.vector(t(x$observed[keep, , drop = FALSE])),
        predicted = as.vector(t(x$prediction[keep, , drop = FALSE])), model = label, fold = fold)]
      molecular_rows[[length(molecular_rows) + 1L]] <- z
    }
  }
  molecular <- rbindlist(molecular_rows)
  write_table(molecular, 'joint_held_molecular_predictions.tsv.gz')
  molecular_summary <- molecular[, .(cells = uniqueN(depmap_id), r = cor_safe(observed, predicted),
    predicted_to_observed_sd = sd(predicted) / sd(observed)), by = .(model, intervention, dose, program)]
  write_table(molecular_summary, 'joint_molecular_variation_by_program.tsv')
  write_table(molecular_summary[, .(median_program_r = median(r, na.rm = TRUE),
    median_predicted_to_observed_sd = median(predicted_to_observed_sd, na.rm = TRUE)),
    by = .(model, intervention)], 'joint_molecular_variation_summary.tsv')
  full_paths <- file.path(root, 'data', paste0('joint_full_', c('joint', 'staged'), '_functional.tsv'))
  own <- rbindlist(lapply(full_paths, fread))
  write_table(own[, .(cells = uniqueN(depmap_id), r = cor_safe(y, predicted), sensitive = sum(y <= -1),
    recovered = sum(y <= -1 & predicted <= -1), resistant = sum(y > -1),
    false_sensitive = sum(y > -1 & predicted <= -1)), by = .(model, drug, dose)], 'joint_own_fitting_summary.tsv')
  own_annotations <- own[, .(r = cor_safe(y, predicted)), by = .(model, drug, dose)]
  own_annotations[, label := sprintf('r = %.2f', r)]
  train_plot <- ggplot(own, aes(y, predicted)) + geom_point(size = .7, alpha = .55, color = '#24587b') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .5, color = '#ac4035') +
    geom_text(data = own_annotations, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.02, vjust = 1.3, size = 2.8) +
    facet_grid(model ~ drug + dose, scales = 'free') + labs(x = 'Measured PRISM effect (log2)',
      y = 'Predicted effect (log2)', title = 'Own fitting data: full developmental models; in-sample') +
    theme_classic(base_size = 9) + theme(aspect.ratio = 1)
  save_plot(train_plot, 'joint_own_fitting_data', 15, 5.8)
  # All paired cell-line profiles are shown, ordered by the common2.5uM endpoint.
  # Actual and predicted native molecular doses retain an identical display order.
  for (drug_name in unique(a$F$drug)) {
    z <- molecular[model == 'joint' & intervention == drug_name]
    assay <- a$F[drug == drug_name & dose == 2.5, .(depmap_id, y)]
    z <- z[depmap_id %in% assay$depmap_id]
    columns <- unique(z[, .(sample_id, depmap_id, dose)])
    columns[, measured_pct := 100 * 2^assay$y[match(depmap_id, assay$depmap_id)]]
    setorder(columns, dose, measured_pct, depmap_id)
    actual <- dcast(z, sample_id ~ program, value.var = 'observed')[match(columns$sample_id, sample_id)]
    predicted <- dcast(z, sample_id ~ program, value.var = 'predicted')[match(columns$sample_id, sample_id)]
    actual <- t(as.matrix(actual[, ..P]))
    predicted <- t(as.matrix(predicted[, ..P]))
    colnames(actual) <- colnames(predicted) <- columns$sample_id
    row_order <- hclust(dist(actual), method = 'ward.D2')$order
    annotations <- data.frame(Dose = factor(columns$dose), row.names = columns$sample_id)
    args <- list(`PRISM2.5uM % vehicle` = ComplexHeatmap::anno_barplot(columns$measured_pct,
      border = FALSE, gp = grid::gpar(fill = ifelse(columns$measured_pct <= 50, '#B2182B', '#4D4D4D')),
      height = grid::unit(1, 'cm')), annotation_name_gp = grid::gpar(fontsize = 8))
    range <- max(abs(c(actual, predicted)))
    filename <- paste0('joint_', gsub('[^a-z0-9]+', '_', tolower(drug_name)), '_molecular_heatmaps.pdf')
    grDevices::pdf(file.path(root, 'figures', filename), width = 11, height = 6)
    for (mode in c('Measured', 'Held-out predicted')) {
      matrix <- if (mode == 'Measured') actual else predicted
      heatmap <- draw_heatmap(matrix[row_order, , drop = FALSE], scale_rows = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE,
        col_anno_df = annotations, col_anno_cols = 'Dose', ann_top_args = args,
        hm_name = 'Displacement', value_range = c(-range, 0, range),
        value_colors = c('#2166AC', '#FFFFFF', '#B2182B'), column_split = factor(columns$dose),
        title = paste(drug_name, mode, '| same order and scale | full34 native program means'))
      ComplexHeatmap::draw(heatmap)
      preview_path <- file.path(root, 'figures', sub('\\.pdf$', paste0('_', gsub('[^a-z]+', '_', tolower(mode)), '.png'), filename))
      grDevices::png(preview_path, width = 1500, height = 820, res = 135)
      ComplexHeatmap::draw(heatmap)
      grDevices::dev.off()
    }
    grDevices::dev.off()
  }
  status('joint restoration complete')
  print(metrics)
}

if (stage %in% c('all', 'joint', 'report')) {
  # Preserve the exact configuration and source provenance with the report.
  inputs <- c(bank_path, joint_source, panel_path,
    file.path(historical_bridge, 'joint40/scripts/02_lib.R'),
    file.path(historical_bridge, 'joint40/linear_decoder/scripts/02_lib.R'),
    file.path(historical_combo, 'cross_context_enrichment/results/models_and_splits.rds'),
    file.path(historical_combo, 'bank40_comparison/results/penalties.tsv'))
  write_table(data.table(path = inputs, sha256 = vapply(inputs, digest::digest, character(1), file = TRUE, algo = 'sha256')),
    'input_hashes.tsv')
  frozen <- list.files(file.path(root, 'data'), '_frozen\\.tsv$', full.names = TRUE)
  write_table(data.table(path = frozen, sha256 = vapply(frozen, digest::digest, character(1), file = TRUE, algo = 'sha256')),
    'prediction_hashes.tsv')
  qsave(list(bank = bank_path, programs = K, genes = length(genes), joint_signature = joint_signature,
    feedback = .25, upstream_ridge = .03, decoder_ridge = .1, held_folds = a$folds,
    settings_source = 'Original fold-specific selected settings; no new tuning'), file.path(root, 'settings.qs'))
  writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
  lines <- c('53: historical joint and nonlinear combination architectures on current43 native34 programs',
    'Same original labels, recipients, study/recipient folds, receptor inputs, dose interpolation and selected settings.',
    'All8429 bank genes retained. The99 absent from old gene caches were recovered from primary counts and GEO inputs in R.',
    'Joint functional input is predicted program displacement alone; sensitivity gradients update the upstream molecular mapping.',
    'Functional zero input returns100% vehicle. Gamma0 staged control uses the same bank and rows.',
    'Combination model uses the original symmetric nonlinear kernels and zero-partner gate; no baseline prediction is used in that test.',
    'Combination target is interaction beyond the measured single-agent effects, not the complete baseline-only combination endpoint.',
    'Current43 discovery used full-panel functional labels. These results are retrospective and conditional on that frozen bank.',
    'No HGF rescue validation or active-reference replacement is claimed.', '',
    'Joint held-recipient/study summaries (pooled doses):', capture.output(print(metrics)), '',
    'Same-drug/same-dose across-cell summaries:', capture.output(print(by_dose)), '',
    'Molecular stage:', capture.output(print(fread(file.path(root, 'data/joint_molecular_variation_summary.tsv')))), '',
    'Combination interaction summaries:', capture.output(print(fread(file.path(root, 'data/combination_summary.tsv')))), '',
    'Interpretation: historical correlation is largely retained; sensitivity magnitude remains underpredicted.',
    'At2.5uM neither bank identifies the43 lapatinib-sensitive or107 crizotinib-sensitive cells at the50% threshold.',
    'Per-dose results are essential: crizotinib pooled correlation is much stronger than across-cell prediction at fixed dose.')
  writeLines(lines, file.path(root, 'REPORT.txt'))
  figures <- file.path(root, 'figures', c('historical_current_joint_bridge_by_dose.pdf',
    'joint_own_fitting_data.pdf', 'historical_current_combinations_all_cells.pdf',
    'historical_current_combination_comparison.pdf', 'joint_lapatinib_ditosylate_molecular_heatmaps.pdf',
    'joint_crizotinib_molecular_heatmaps.pdf'))
  check(all(file.exists(figures)), 'Restoration figure set incomplete')
  write_table(data.table(file = figures, contents = c('Held-out same-drug/dose old/current comparison',
    'Current models own fitting data', 'All seven held recipients: nonlinear interaction',
    'HT29: nonlinear and product model controls', 'Lapatinib actual/predicted displacement in identical order',
    'Crizotinib actual/predicted displacement in identical order')), 'figure_index.tsv')
  status('COMPLETE: models, source tables, comparison figures and report saved')
}
