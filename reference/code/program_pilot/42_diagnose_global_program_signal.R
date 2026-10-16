#!/usr/bin/env Rscript
# Descriptive audit of the frozen global bank. No fitting, selection or refitting.
suppressPackageStartupMessages({library(data.table); library(qs)})
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(4L)
root <- '/work/results/program_pilot/42_global_pooled_signature_programs_ward'
bank <- qread(file.path(root, 'models/program_bank.qs'))
panel <- qread('/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs')
md <- copy(bank$metadata)
a <- panel$sensitivity[match(md$observation_id, observation_id)]
stopifnot(identical(md$observation_id, a$observation_id), nrow(bank$scores) == nrow(md))
md[, abundance := a$prism_abundance_pct]
md[, context := paste(drug, dose, sep = '@')]
md[, score_row := .I]
programs <- colnames(bank$scores)
keep <- which(is.finite(md$abundance))
x <- bank$scores[keep, , drop = FALSE]
d <- md[keep]
y <- d$abundance
stopifnot(!anyDuplicated(d[, paste(depmap_id, context)]))
# Drug-dose and cell-line centering here is descriptive, not a fitted predictor.
within_x <- x
within_y <- y
for (idx in split(seq_len(nrow(d)), d$context)) {
  within_x[idx, ] <- sweep(x[idx, , drop = FALSE], 2, colMeans(x[idx, , drop = FALSE]), '-')
  within_y[idx] <- y[idx] - mean(y[idx])
}
cell_x <- x
cell_y <- y
for (idx in split(seq_len(nrow(d)), d$depmap_id)) {
  cell_x[idx, ] <- sweep(x[idx, , drop = FALSE], 2, colMeans(x[idx, , drop = FALSE]), '-')
  cell_y[idx] <- y[idx] - mean(y[idx])
}
cx <- rowsum(x, d$context) / as.numeric(table(d$context)[rownames(rowsum(x, d$context))])
cy <- tapply(y, d$context, mean)[rownames(cx)]
correlations <- data.table(program = programs,
  pooled_r = as.numeric(cor(x, y)),
  between_drug_dose_means_r = as.numeric(cor(cx, cy)),
  across_cells_within_drug_dose_r = as.numeric(cor(within_x, within_y)),
  across_drugs_within_cell_r = as.numeric(cor(cell_x, cell_y)),
  fraction_score_variation_within_drug_dose = colSums(within_x^2) / colSums(scale(x, scale = FALSE)^2))
fwrite(correlations, file.path(root, 'data/signal_axis_correlations.tsv'), sep = '\t')
# Assay contrasts and program associations within each exact drug and dose.
context_tables <- list()
for (context in unique(d$context)) {
  idx <- which(d$context == context)
  if (length(idx) < 6L || sd(y[idx]) == 0) next
  r <- as.numeric(cor(x[idx, , drop = FALSE], y[idx]))
  context_tables[[context]] <- data.table(drug = d$drug[idx[1]], dose = d$dose[idx[1]],
    cells = length(idx), sensitive = sum(y[idx] <= 50), resistant = sum(y[idx] > 50),
    abundance_q10 = quantile(y[idx], .1), abundance_q90 = quantile(y[idx], .9),
    program = programs, r = r)
}
ctx <- rbindlist(context_tables)
fwrite(ctx, file.path(root, 'data/signal_within_drug_dose_correlations.tsv'), sep = '\t')
# Program redundancy, within-cluster agreement and score mixing from other genes.
relationships <- qread(file.path(root, 'data/global_gene_relationships_rank60.qs'))
w <- bank$weights
coherence <- rbindlist(lapply(seq_along(programs), function(j) {
  idx <- which(bank$partition == j)
  outside <- which(bank$partition != j)
  vals <- lapply(relationships$matrices, function(m) m[idx, idx, drop = FALSE][upper.tri(m[idx, idx, drop = FALSE])])
  data.table(program = programs[j], genes = length(idx),
    weight_inside_seed_fraction = sum(w[idx, j]) / sum(w[, j]),
    displacement_within_median_r = median(vals$gene_correlation),
    loading_within_median_r = median(vals$gene_functional_loading_correlation),
    displacement_outside_median_r = median(relationships$matrices$gene_correlation[idx, outside]))
}))
fwrite(coherence, file.path(root, 'data/signal_program_coherence.tsv'), sep = '\t')
# For the requested odd-cell audit: identify low displacement/control ratios
# and disagreement between the two assays without relabeling or dropping data.
qc <- md[, .(signatures = .N, median_response_control_ratio = median(response_control_rms_ratio),
  fraction_response_below_control = mean(response_control_rms_ratio < 1),
  median_displacement_rms = median(raw_response_rms)), by = .(depmap_id, recipient)]
fwrite(qc, file.path(root, 'data/signal_cell_molecular_quality.tsv'), sep = '\t')
matched <- which(is.finite(a$prism_abundance_pct) & is.finite(a$almanac_abundance_pct))
score_cor <- cor(x)
report <- capture.output({
  cat('PRISM: ', nrow(d), ' signatures; ', uniqueN(d$depmap_id), ' cells; ', uniqueN(d$drug), ' drugs\n', sep = '')
  cat('Median / max absolute correlation by comparison axis:\n')
  print(correlations[, lapply(.SD, function(v) c(median = median(abs(v)), max = max(abs(v)))),
    .SDcols = c('pooled_r', 'between_drug_dose_means_r', 'across_cells_within_drug_dose_r', 'across_drugs_within_cell_r')])
  cat('Programs with highest pooled association:\n')
  print(correlations[order(-abs(pooled_r))][1:8])
  cat('Program coherence and weight fractions (quantiles):\n')
  print(coherence[, lapply(.SD, quantile, probs = c(0, .25, .5, .75, 1)), .SDcols = 3:6])
  cat('Median pairwise program-score correlation: ', median(score_cor[upper.tri(score_cor)]), '\n')
  cat('Top score-PC variance fraction: ', max(eigen(score_cor, symmetric = TRUE, only.values = TRUE)$values) / length(programs), '\n')
  cat('Matched PRISM+ALMANAC: ', length(matched), '; conflicting sensitive/resistant labels: ',
    sum((a$prism_abundance_pct[matched] <= 50) != (a$almanac_abundance_pct[matched] <= 50)), '\n')
  cat('Lowest molecular/control signal cells:\n')
  print(qc[order(median_response_control_ratio)][1:8])
  cat('Focal drug-dose contrasts and maximum observed association (descriptive, 60 programs):\n')
  print(ctx[drug %in% c('everolimus', 'sirolimus', 'erlotinib', 'gefitinib', 'lapatinib', 'neratinib', 'trametinib'),
    .(cells = first(cells), sensitive = first(sensitive), resistant = first(resistant),
      q10 = first(abundance_q10), q90 = first(abundance_q90), max_abs_r = max(abs(r), na.rm = TRUE)), by = .(drug, dose)])
})
writeLines(report, file.path(root, 'signal_diagnosis.txt'))
cat(paste(report, collapse = '\n'), '\n')

# Isolate scoring from discovery: compare the same fixed clusters using only
# their own gene means. This diagnostic does not replace the saved bank/scores.
hard_weights <- Matrix::sparseMatrix(i = seq_along(bank$genes), j = bank$partition,
  x = 1 / tabulate(bank$partition, length(programs))[bank$partition],
  dims = c(length(bank$genes), length(programs)))
hard <- matrix(NA_real_, nrow(md), length(programs), dimnames = list(md$observation_id, programs))
for (id in unique(md$depmap_id)) {
  molecular <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
  idx <- which(md$depmap_id == id)
  values <- molecular$response[md$observation_id[idx], bank$genes, drop = FALSE]
  hard[idx, ] <- as.matrix(values %*% hard_weights)
}
stopifnot(all(is.finite(hard)))
hx <- hard[keep, , drop = FALSE]
hwithin <- hx
for (idx in split(seq_len(nrow(d)), d$context)) {
  hwithin[idx, ] <- sweep(hx[idx, , drop = FALSE], 2, colMeans(hx[idx, , drop = FALSE]), '-')
}
hard_check <- data.table(program = programs, pooled_r = as.numeric(cor(hx, y)),
  across_cells_within_drug_dose_r = as.numeric(cor(hwithin, within_y)))
fwrite(hard_check, file.path(root, 'data/signal_hard_mean_diagnostic.tsv'), sep = '\t')
hard_context <- rbindlist(lapply(unique(ctx[, paste(drug, dose, sep = '@')]), function(context) {
  idx <- which(d$context == context)
  data.table(drug = d$drug[idx[1]], dose = d$dose[idx[1]], program = programs,
    hard_mean_r = as.numeric(cor(hx[idx, , drop = FALSE], y[idx])))
}))
paired <- merge(ctx, hard_context, by = c('drug', 'dose', 'program'))
fwrite(paired, file.path(root, 'data/signal_soft_vs_hard_by_context.tsv'), sep = '\t')
extra <- capture.output({
  cat('\nSCORING ISOLATION: unchanged clusters; own-gene means only\n')
  hc <- cor(hx)
  cat('Median score-score r: ', median(hc[upper.tri(hc)]), '\n')
  cat('Maximum absolute pooled r: ', max(abs(hard_check$pooled_r)), '\n')
  cat('Maximum absolute across-cell within-drug-dose r: ', max(abs(hard_check$across_cells_within_drug_dose_r)), '\n')
  cat('Within-context focal comparisons, same data and same clusters:\n')
  print(paired[drug %in% c('sirolimus', 'trametinib', 'lapatinib'),
    .(sensitive = first(sensitive), resistant = first(resistant),
      best_soft_abs_r = max(abs(r)), best_hard_abs_r = max(abs(hard_mean_r))), by = .(drug, dose)])
})
cat(paste(extra, collapse = '\n'), '\n')
write(extra, file.path(root, 'signal_diagnosis.txt'), append = TRUE)
