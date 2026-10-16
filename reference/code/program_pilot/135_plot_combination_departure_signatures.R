#!/usr/bin/env Rscript
# Descriptive v106-style correlations for the existing v135 training data.
# No functional model is fitted, and no new sensitivity signatures are claimed.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs)})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
root <- '/work/results/program_pilot/135_v128_with_corrected_HGF_inputs/linear'
out <- file.path(root, 'figures/combination_departure_signatures')
data_out <- file.path(root, 'data/combination_departure_signatures')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(data_out, recursive = TRUE, showWarnings = FALSE)
programs <- sprintf('P%02d', 1:40)
model <- qread(file.path(root, 'models/final_normalized_cell_count.qs'))
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
bank_selection <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
drugs <- bank_selection$drug
contexts <- fread(file.path(root, 'data/molecular_contexts.tsv'))
extra_contexts <- contexts[11:nrow(contexts)]

# Load frozen prediction functions only; these do not refit RNA or function.
required <- c('context_key', 'context_index', 'predict_rna', 'signature_node', 'predict_single')
expressions <- parse('/work/src/program_pilot/135_v128_with_corrected_HGF_inputs.R')
found <- character()
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name('<-')) &&
      is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% required) {
    eval(expression, envir = .GlobalEnv)
    found <- c(found, as.character(expression[[2]]))
  }
}
stopifnot(setequal(required, found))

# The same zero-referenced molecular calibration used by g; no outcome enters it.
scale_ligands <- function(D) {
  ctx <- model$context
  Q <- cbind(1, sweep(sweep(sweep(D, 2, ctx$fm), 2, ctx$fs, '/') %*% ctx$fv, 2, ctx$fscale, '/'))
  map <- solve(crossprod(ctx$Q[1:10, ]) / 10 + diag(c(0, 1, 1, 1)),
               crossprod(ctx$Q[1:10, ], log(ctx$sigma[1:10, ])) / 10)
  low <- matrix(apply(ctx$sigma[1:10, ], 2, min), nrow(D), 40, byrow = TRUE)
  high <- matrix(apply(ctx$sigma[1:10, ], 2, max), nrow(D), 40, byrow = TRUE)
  D / pmin(pmax(exp(Q %*% map), low), high)
}
short_drug <- function(x) {
  x <- sub(' ditosylate| \\(hydrochloride\\)', '', x)
  sub('5-Azacytidine', 'Azacitidine', x, fixed = TRUE)
}

pairs <- fread(file.path(root, 'data/combination_training_outcomes.tsv'))
stopifnot(max(abs(pairs$I - (pairs$y - pairs$yA - pairs$yB))) < 1e-10)
A <- predict_single(model, pairs[, .(depmap_id, drug = drugA, dose = doseA)])$program_activity[, 1:40]
B <- predict_single(model, pairs[, .(depmap_id, drug = drugB, dose = doseB)])$program_activity[, 1:40]
pairs[, context := paste0(short_drug(drugA), ' ', doseA, ' + ', short_drug(drugB), ' ', doseB)]
pairs[, `:=`(outcome = I, sample = CELLNAME)]

external <- fread(file.path(root, 'data/external_combination_outcomes.tsv'))
predicted <- fread(file.path(root, 'data/external_predicted_drug_programs.tsv'))
key <- function(d) paste(d$depmap_id, d$gene, d$drug, d$dose, sep = '|')
stopifnot(!anyDuplicated(key(predicted)))
index <- match(key(external), key(predicted))
stopifnot(!anyNA(index))
EA <- as.matrix(predicted[index, ..programs])
for (dr in unique(external$drug)) {
  rows <- which(external$drug == dr)
  if (dr == 'Lapatinib ditosylate') {
    scale <- model$context$sigma[context_index(external$drug[rows], external$dose[rows]), , drop = FALSE]
    EA[rows, ] <- EA[rows, , drop = FALSE] / scale
  } else EA[rows, ] <- sweep(EA[rows, , drop = FALSE], 2, model$external_rna[[dr]]$sigma, '/')
}
ligand <- fread(file.path(root, 'data/external_ligand_programs.tsv'))
index <- match(paste(external$depmap_id, external$gene), paste(ligand$depmap_id, ligand$gene))
stopifnot(!anyNA(index))
EB <- scale_ligands(as.matrix(ligand[index, ..programs]))
external[, `:=`(context = paste0(gene, ' + ', short_drug(drug), ' ', dose),
                outcome = observed_departure_log2, sample = CellLine)]
stopifnot(max(abs(external$outcome - (external$observed_drug_arm_increment_log2 - external$observed_ligand_log2))) < 1e-10)

# Each physical observation occurs once. HGF source alternatives are not counted
# as additional cells: the study-balanced transported profile is shown here.
blocks <- list(ALMANAC = list(meta = pairs, A = A, B = B,
  titles = c('Drug A programs', 'Drug B programs'),
  title = 'ALMANAC | across cells at each fixed drug pair and dose',
  note = 'Outcome: log2[V(A+B)/(V(A)V(B))]. Positive = protection; negative = excess killing.\nA/B follow the displayed drug order. Displacements predicted in g units; doses in µM.'))
for (dr in c('PD-184352', 'PLX-4720', 'Lapatinib ditosylate')) {
  rows <- which(external$drug == dr)
  id <- if (dr == 'Lapatinib ditosylate') 'HER2' else dr
  note <- if (id == 'HER2')
    'Outcome: log2[V(drug+ligand)/(V(drug)V(ligand))], final abundance at 72 h.\nNRG1: six cells; EGF: two cells (grey, insufficient for correlation). Drug predicted; ligand transported.' else
    'Outcome: reported log2 growth interaction; net-growth readout, not final cell abundance.\nSix cells per ligand at a fixed drug dose. Drug predicted; ligand transported; no combination RNA measured.'
  blocks[[id]] <- list(meta = external[rows], A = EA[rows, , drop = FALSE], B = EB[rows, , drop = FALSE],
    titles = c('Drug displacement', 'Ligand displacement'),
    title = paste0(id, ' | across cells within each ligand + drug/dose'), note = note)
}
rows <- which(external$study == 'Straussman2012')
across <- copy(external[rows])
across[, context := paste0(CellLine, ' | ', drug, ' ', dose)]
across[, sample := gene]
blocks[['Growth_across_ligands']] <- list(meta = across, B = EB[rows, , drop = FALSE],
  titles = 'Ligand displacement', title = 'RAF/MEK experiments | across ligands within each cell and drug',
  note = 'Outcome: reported net-growth interaction. Each column compares 13 ligands within one cell/drug.\nTransported ligand programs in g molecular units; fixed drug vector has no across-ligand variation.')

mema <- fread(file.path(root, 'data/mema_nonchallenge_training_outcomes.tsv'))[endpoint == 'normalized_cell_count']
ml <- fread(file.path(root, 'data/mema_molecular_inputs.tsv'))
index <- match(paste(mema$CellLine, mema$gene), paste(ml$CellLine, ml$gene))
stopifnot(!anyNA(index))
mema[, `:=`(context = CellLine, sample = gene, outcome = y - mean(y)), by = CellLine]
blocks[['MEMA_relative_only']] <- list(meta = mema, B = scale_ligands(as.matrix(ml[index, ..programs])),
  titles = 'Ligand displacement', title = 'MEMA | relative ligand response, not an absolute interaction',
  note = 'Outcome: cell-centred log2[drug/vehicle under the same ligand]; lapatinib 0.75 µM, 72 h.\n26 training ligands per cell. No ligand-free control: absolute departure cannot be calculated. HGF/NRG1 excluded.')

correlations <- observations <- list()
for (id in names(blocks)) {
  block <- blocks[[id]]
  meta <- block$meta
  stopifnot(!anyDuplicated(meta[, .(context, sample)]))
  observations[[id]] <- meta[, .(dataset = id, context, sample, outcome)]
  parts <- intersect(c('A', 'B'), names(block))
  for (part in parts) {
    X <- block[[part]]
    stopifnot(nrow(X) == nrow(meta), ncol(X) == 40, all(is.finite(X)))
    for (context in unique(meta$context)) {
      ii <- which(meta$context == context)
      y <- meta$outcome[ii]
      r <- vapply(seq_len(40), function(j) {
        x <- X[ii, j]
        if (length(y) < 4 || sd(y) < 1e-10 || sd(x) < 1e-10) return(NA_real_)
        cor(x, y)
      }, numeric(1))
      correlations[[paste(id, part, context)]] <- data.table(dataset = id, partner = part, context,
        program = programs, n = length(ii), outcome_SD = sd(y), r = r)
    }
  }
}
corr <- rbindlist(correlations)
fwrite(corr, file.path(data_out, 'program_departure_correlations.tsv'), sep = '\t')
fwrite(rbindlist(observations), file.path(data_out, 'observed_outcomes.tsv'), sep = '\t')
fwrite(corr[, .(n = first(n), outcome_SD = first(outcome_SD), finite_programs = sum(is.finite(r))),
  by = .(dataset, partner, context)], file.path(data_out, 'correlation_coverage.tsv'), sep = '\t')

# One row order across every page, learned descriptively from all finite columns.
wide <- dcast(corr, program ~ dataset + partner + context, value.var = 'r')
all_cor <- as.matrix(wide[match(programs, program), -1])
keep <- colSums(is.finite(all_cor)) == 40L
stopifnot(any(keep))
row_order <- hclust(dist(all_cor[, keep, drop = FALSE]), method = 'ward.D2')$order
fwrite(data.table(program = programs[row_order], order = seq_len(40)), file.path(data_out, 'shared_program_order.tsv'), sep = '\t')
collapsed_tables <- list()
render <- function(id) {
  block <- blocks[[id]]
  parts <- intersect(c('A', 'B'), names(block))
  columns <- unique(block$meta$context)
  matrices <- lapply(parts, function(part) {
    x <- dcast(corr[dataset == id & partner == part], program ~ context, value.var = 'r')
    M <- as.matrix(x[match(programs, program), ..columns])
    dimnames(M) <- list(programs, columns)
    M
  })
  joint <- do.call(rbind, matrices)
  valid <- which(colSums(is.finite(joint)) == nrow(joint))
  invalid <- setdiff(seq_along(columns), valid)
  group <- rep(NA_integer_, length(columns))
  if (length(valid) >= 2) {
    tree <- hclust(dist(t(joint[, valid, drop = FALSE])), method = 'ward.D2')
    order <- c(valid[tree$order], invalid)
    if (length(valid) >= 4) group[valid] <- cutree(tree, k = 2)
  } else order <- seq_along(columns)
  fwrite(data.table(context = columns[order], group = group[order], order = seq_along(order)),
    file.path(data_out, paste0(id, '_column_groups.tsv')), sep = '\t')
  ht <- NULL
  for (j in seq_along(parts)) {
    M <- matrices[[j]][row_order, order, drop = FALSE]
    coverage <- corr[dataset == id & partner == parts[j], .(n = first(n)), by = context]
    annotation <- data.frame(N = as.character(coverage$n[match(columns[order], coverage$context)]),
                             row.names = columns[order])
    n_levels <- sort(unique(annotation$N))
    palette <- list(N = setNames(gray.colors(length(n_levels), start = .78, end = .32), n_levels))
    if (any(!is.na(group))) {
      annotation$Group <- paste0('Group ', group[order])
      palette$Group <- c('Group 1' = '#D58B36', 'Group 2' = '#775DA6')
    }
    stopifnot(identical(colnames(M), rownames(annotation)))
    hm <- draw_heatmap(M, scale_rows = FALSE, value_range = c(-1, 0, 1),
      value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = paste0('Pearson r ', parts[j]),
      col_anno_df = annotation, col_anno_cols = names(annotation), palette_overrides = palette,
      na_col = '#D0D0D0', show_heatmap_legend = j == 1,
      heatmap_legend_param = list(title = 'Pearson r'),
      show_row_names = j == 1, row_names_side = 'left', row_names_gp = grid::gpar(fontsize = 8),
      show_column_names = TRUE, column_names_rot = 55, column_names_gp = grid::gpar(fontsize = 7),
      column_names_max_height = grid::unit(65, 'mm'), width = grid::unit(max(28, 6 * ncol(M)), 'mm'),
      height = grid::unit(145, 'mm'), title = if (ncol(M) <= 4) gsub(' ', '\n', block$titles[j]) else block$titles[j],
      column_title_gp = grid::gpar(fontsize = 10))
    ht <- if (is.null(ht)) hm else ht + hm
    if (any(!is.na(group))) {
      C <- vapply(1:2, function(k) rowMeans(matrices[[j]][, which(group == k), drop = FALSE]), numeric(40))
      dimnames(C) <- list(programs, c('Group 1', 'Group 2'))
      collapsed_tables[[paste(id, parts[j])]] <<- data.table(dataset = id, partner = parts[j], program = programs, C)
      summary <- draw_heatmap(C[row_order, , drop = FALSE], scale_rows = FALSE,
        value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'),
        hm_name = paste0('Mean r ', parts[j]), show_heatmap_legend = FALSE,
        show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 55,
        column_names_gp = grid::gpar(fontsize = 7), width = grid::unit(16, 'mm'),
        height = grid::unit(145, 'mm'), title = 'Means', column_title_gp = grid::gpar(fontsize = 9))
      ht <- ht + summary
    }
  }
  ComplexHeatmap::draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'right',
    column_title = block$title, column_title_gp = grid::gpar(fontsize = 12),
    padding = grid::unit(c(19, 5, 4, 4), 'mm'), ht_gap = grid::unit(5, 'mm'))
  grid::grid.text(paste0(block$note, '\nDescriptive correlations, not held-out prediction. Two-group means summarize clustering; they are not fitted routes.'),
    x = .5, y = grid::unit(3, 'mm'), just = c('centre', 'bottom'), gp = grid::gpar(fontsize = 7.5))
}
book <- file.path(out, 'all_combination_departure_heatmaps.pdf')
cairo_pdf(book, width = 11.7, height = 10.5, onefile = TRUE)
for (id in names(blocks)) render(id)
dev.off()
for (id in names(blocks)) {
  columns <- uniqueN(blocks[[id]]$meta$context)
  parts <- length(intersect(c('A', 'B'), names(blocks[[id]])))
  width <- max(6.7, min(11.7, (parts * (max(28, columns * 6) + ifelse(columns >= 4, 21, 0)) + 60) / 25.4))
  cairo_pdf(file.path(out, paste0(id, '.pdf')), width = width, height = 10.5)
  render(id)
  dev.off()
  png(file.path(out, paste0(id, '.png')), width = width, height = 10.5, units = 'in', res = 140, type = 'cairo')
  render(id)
  dev.off()
}
fwrite(rbindlist(collapsed_tables), file.path(data_out, 'descriptive_collapsed_correlations.tsv'), sep = '\t')
writeLines(c('Existing v135 inputs and measured functional outcomes; no model fitting.',
  '40 programme rows; Pearson r, no row z-scoring. Identical row order on all pages.',
  'Drug vectors are inferred from baseline; ligand vectors are transported CytoSig responses.',
  'Two groups are descriptive hierarchical-clustering summaries, not validated signatures.',
  'HGF alternate source views are not counted as independent observations.',
  'ALMANAC and HER2 count departures, Straussman net-growth departures, and MEMA relative outcomes are separated.',
  'Grey entries are undefined (fewer than four independent observations or no variation).',
  'Wilson outcomes are absent. Previously inspected data make this exploratory, not validation.'),
  file.path(data_out, 'READ_ME.txt'))
cat('Created:', book, '\n')
