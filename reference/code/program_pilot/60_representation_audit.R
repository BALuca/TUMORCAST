#!/usr/bin/env Rscript
# Molecular-only audit: frozen program banks versus genes and prior biological sets.
# No functional outcome is used to choose genes, programs, contrasts or thresholds.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(limma)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/60_defensible_response_composition'
prior <- '/work/results/program_pilot/59_mechanistic_evidence_datasets/data'
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
out <- function(x, name) fwrite(x, file.path(root, 'data', paste0('representation_', name, '.tsv.gz')), sep = '\t')
plot_save <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0('representation_', name, '.pdf')), p,
    device = cairo_pdf, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0('representation_', name, '.png')), p, width = width, height = height, dpi = 130)
}
cor_safe <- function(x, y) {
  good <- is.finite(x) & is.finite(y)
  if (sum(good) < 4L || sd(x[good]) < 1e-10 || sd(y[good]) < 1e-10) return(NA_real_)
  cor(x[good], y[good])
}

# Existing discovery results remain frozen. The same hard assigned-gene means
# are used for the three empirical banks; all prior biological sets are retained.
paths <- c(bank43 = '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs',
  bank16 = '/work/results/program_pilot/16_global_mean_relationships_fixed_ward_programs/data/global_gene_assignments.tsv',
  bank42 = '/work/results/program_pilot/42_global_pooled_signature_programs_ward/models/program_bank.qs',
  biology = '/work/results/program_pilot/57_erbb_met_hgf_program_biology/data/biological_readout_membership.tsv.gz',
  landmark = file.path(prior, 'matched_gene_annotation.tsv'))
stopifnot(all(file.exists(paths)))
membership <- list()
for (bn in c('bank43', 'bank42')) {
  obj <- qread(paths[[bn]])
  stopifnot(!anyDuplicated(obj$genes), all(obj$genes %in% names(obj$partition)))
  membership[[bn]] <- data.table(bank = bn, gene = obj$genes,
    feature = sprintf('P%02d', obj$partition[obj$genes]))
}
a16 <- fread(paths[['bank16']])
stopifnot(!anyDuplicated(a16$gene), uniqueN(a16$global_program) == 40L)
membership[['bank16']] <- a16[, .(bank = 'bank16', gene, feature = global_program)]
membership[['biology']] <- unique(fread(paths[['biology']])[, .(bank = 'biology', gene, feature = readout)])
members <- rbindlist(membership)
members[, weight := 1 / .N, by = .(bank, feature)]
members[, feature_id := paste(bank, feature, sep = '__')]
stopifnot(!anyDuplicated(members[, .(bank, feature, gene)]))
out(members, 'membership')
fwrite(members[, .(representation = bank, feature, gene, weight)],
  file.path(root, 'data/representation_membership.tsv'), sep = '\t')
weights <- lapply(split(members, by = 'bank'), function(z) {
  genes <- sort(unique(z$gene))
  features <- sort(unique(z$feature))
  sparseMatrix(i = match(z$gene, genes), j = match(z$feature, features), x = z$weight,
    dims = c(length(genes), length(features)), dimnames = list(genes, features))
})
qsave(weights, file.path(root, 'models/representation_frozen_weights.qs'))
out(data.table(bank = names(paths), path = unname(paths),
  sha256 = vapply(paths, digest, character(1), file = TRUE, algo = 'sha256')), 'source_manifest')
landmark <- unique(fread(paths[['landmark']])$gene_symbol)
stopifnot(length(landmark) == 978L)
message('Frozen weights written: ', paste(names(weights), lengths(lapply(weights, colnames)), collapse = '; '))

# Score each set by its measured genes, explicitly reporting missing coverage.
# The sparse-L1000 experiment keeps only the 978 directly measured landmarks;
# inferred L1000 genes are not substituted for measurements.
all_scores <- all_contrasts <- gene_summaries <- gene_coherence <- list()
replicate_agreement <- sparse_agreement <- coverage_parts <- list()
contrast_gene_list <- list()
for (cell in c('MCF7', 'LNCaP')) {
  obj <- qread(file.path(prior, paste0('diaz_', cell, '_molecular.qs')))
  x <- obj$log2_one_plus_TMM_CP10K
  delta <- obj$displacement
  md <- as.data.table(obj$samples)
  stopifnot(identical(colnames(x), md$sample_id), identical(dimnames(x), dimnames(delta)), all(is.finite(x)))
  # Detection only, before contrasts: CP10K >= 0.1 (CPM >= 1) in >=3 samples.
  detected <- rowSums(x >= log2(1.1)) >= 3L
  out(data.table(cell_line = cell, gene = rownames(x), detected = detected), paste0(cell, '_gene_universe'))
  score_parts <- list()
  for (bn in names(weights)) {
    ww <- weights[[bn]]
    full_genes <- intersect(rownames(ww), rownames(x))
    sparse_genes <- intersect(full_genes, landmark)
    full_w <- ww[full_genes, , drop = FALSE]
    sparse_w <- ww[sparse_genes, , drop = FALSE]
    full_mass <- colSums(full_w)
    sparse_mass <- colSums(sparse_w)
    n_full <- colSums(full_w != 0)
    n_sparse <- colSums(sparse_w != 0)
    full_scores <- as.matrix(crossprod(full_w, delta[full_genes, , drop = FALSE])) / full_mass
    sparse_scores <- as.matrix(crossprod(sparse_w, delta[sparse_genes, , drop = FALSE])) / sparse_mass
    sparse_scores[sparse_mass == 0, ] <- NA_real_
    full_expression <- as.matrix(crossprod(full_w, x[full_genes, , drop = FALSE])) / full_mass
    rownames(full_scores) <- rownames(sparse_scores) <- rownames(full_expression) <- paste(bn, colnames(ww), sep = '__')
    coverage_parts[[paste(cell, bn)]] <- data.table(cell_line = cell, bank = bn, feature = colnames(ww),
      feature_id = rownames(full_scores), total_genes = colSums(ww != 0), full_genes = n_full,
      landmark_genes = n_sparse, full_fraction = full_mass, landmark_fraction = sparse_mass)
    score_parts[[bn]] <- full_expression
    for (mode in c('full', 'landmark')) {
      mat <- if (mode == 'full') full_scores else sparse_scores
      lt <- as.data.table(as.table(mat))
      setnames(lt, c('feature_id', 'sample_id', 'value'))
      lt[, `:=`(bank = bn, measurement = mode)]
      all_scores[[paste(cell, bn, mode)]] <- merge(lt, md, by = 'sample_id')
    }
    # Measurement validation uses molecular conditions only, averaging the
    # deposited triplicates. No cell-count value enters these decisions.
    eval_md <- md[arm != 'DMSO' & time_h > 0 & time_h <= 24]
    condition_levels <- unique(eval_md$condition_id)
    full_means <- vapply(condition_levels, function(cc) rowMeans(full_scores[, md$condition_id == cc, drop = FALSE]),
      numeric(nrow(full_scores)))
    sparse_means <- vapply(condition_levels, function(cc) rowMeans(sparse_scores[, md$condition_id == cc, drop = FALSE]),
      numeric(nrow(sparse_scores)))
    sparse_agreement[[paste(cell, bn)]] <- rbindlist(lapply(seq_len(nrow(full_scores)), function(j) {
      aa <- full_means[j, ]
      bb <- sparse_means[j, ]
      informative <- abs(aa) >= .05 & is.finite(bb)
      slope <- if (all(is.finite(bb)) && sd(aa) > 1e-10) coef(lm(bb ~ aa))[2] else NA_real_
      data.table(cell_line = cell, bank = bn, feature = colnames(ww)[j], feature_id = rownames(full_scores)[j],
        n_conditions = length(aa), landmark_genes = n_sparse[j], r = cor_safe(aa, bb), slope = slope,
        informative_conditions = sum(informative),
        direction_agreement = if (sum(informative) > 0) mean(sign(aa[informative]) == sign(bb[informative])) else NA_real_)
    }))
  }
  feature_expression <- do.call(rbind, score_parts)
  # Replicate agreement reflects the deposited wells, not independent studies.
  reps <- unique(md$replicate)
  stopifnot(length(reps) == 3L)
  for (pair in combn(reps, 2L, simplify = FALSE)) {
    aa <- md[replicate == pair[1] & arm != 'DMSO' & time_h <= 24]
    bb <- md[replicate == pair[2] & arm != 'DMSO' & time_h <= 24]
    bb <- bb[match(aa$condition_id, bb$condition_id)]
    stopifnot(identical(aa$condition_id, bb$condition_id))
    for (bn in c('genes', names(weights))) {
      mat <- if (bn == 'genes') delta[detected, , drop = FALSE] else
        do.call(rbind, lapply(all_scores[paste(cell, bn, 'full')], function(z) {
          ss <- dcast(z, feature_id ~ sample_id, value.var = 'value')
          mm <- as.matrix(ss[, -1]); rownames(mm) <- ss$feature_id; mm
        }))
      rr <- vapply(seq_len(nrow(aa)), function(j) cor_safe(mat[, aa$sample_id[j]], mat[, bb$sample_id[j]]), numeric(1))
      replicate_agreement[[paste(cell, bn, paste(pair, collapse = '-'))]] <- data.table(cell_line = cell,
        bank = bn, condition_id = aa$condition_id, arm = aa$arm, time_h = aa$time_h,
        replicate_pair = paste(pair, collapse = '-'), r = rr)
    }
  }

  # Saturated condition means with within-condition residuals. The interaction
  # is AB-A-B+vehicle on log2(1+CP10K), not Bliss and not a causal decomposition.
  group <- factor(paste(md$arm, md$time_h, sep = '_'))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  cc <- list()
  cm <- list()
  for (tt in c(0, 3, 6, 9, 12, 24)) {
    for (pair in c('TM', 'TW', 'MW')) {
      aa <- strsplit(pair, '')[[1]]
      terms <- paste(c(pair, aa, 'DMSO'), tt, sep = '_')
      stopifnot(all(terms %in% colnames(design)))
      v <- setNames(numeric(ncol(design)), colnames(design))
      v[terms] <- c(1, -1, -1, 1)
      id <- paste(pair, tt, sep = '_')
      cc[[id]] <- v
      cm[[id]] <- data.table(contrast_id = id, arm = pair, time_h = tt, contrast_type = 'pair_interaction')
    }
    if (cell == 'LNCaP') {
      terms <- paste(c('TMW', 'TM', 'TW', 'MW', 'T', 'M', 'W', 'DMSO'), tt, sep = '_')
      v <- setNames(numeric(ncol(design)), colnames(design))
      v[terms] <- c(1, -1, -1, -1, 1, 1, 1, -1)
      id <- paste('TMW', tt, sep = '_')
      cc[[id]] <- v
      cm[[id]] <- data.table(contrast_id = id, arm = 'TMW', time_h = tt, contrast_type = 'three_way_interaction')
    }
  }
  contrast_matrix <- do.call(cbind, cc)
  contrast_meta <- rbindlist(cm)
  for (type in c('gene', 'feature')) {
    yy <- if (type == 'gene') x[detected, , drop = FALSE] else feature_expression
    if (type == 'feature') {
      supported <- rowSums(is.finite(yy)) == ncol(yy)
      out(data.table(cell_line = cell, feature_id = rownames(yy), analyzed = supported,
        reason = ifelse(supported, 'Measured gene mean available', 'No measured member genes in full RNA')),
        paste0(cell, '_feature_availability'))
      yy <- yy[supported, , drop = FALSE]
    }
    fit <- eBayes(contrasts.fit(lmFit(yy, design), contrast_matrix), trend = TRUE, robust = TRUE)
    lt <- rbindlist(lapply(seq_len(ncol(contrast_matrix)), function(j) {
      data.table(id = rownames(yy), contrast_id = colnames(contrast_matrix)[j],
        effect = fit$coefficients[, j], standard_error = fit$stdev.unscaled[, j] * sqrt(fit$s2.post),
        p_value = fit$p.value[, j], q_value = p.adjust(fit$p.value[, j], 'BH'))
    }))
    lt <- merge(lt, contrast_meta, by = 'contrast_id')
    lt[, `:=`(cell_line = cell, uncertainty = 'limma moderated within-deposited-well variation; biological independence unresolved')]
    if (type == 'gene') {
      setnames(lt, 'id', 'gene')
      lt[, clear_interaction := q_value <= .05 & abs(effect) >= .25]
      out(lt, paste0(cell, '_gene_interactions'))
      contrast_gene_list[[cell]] <- lt
      gene_summaries[[cell]] <- lt[, .(genes_tested = .N, clear_genes = sum(clear_interaction),
        positive = sum(clear_interaction & effect > 0), negative = sum(clear_interaction & effect < 0)),
        by = .(cell_line, contrast_id, arm, time_h, contrast_type)]
    } else {
      setnames(lt, 'id', 'feature_id')
      lt <- merge(lt, unique(members[, .(bank, feature, feature_id)]), by = 'feature_id')
      all_contrasts[[cell]] <- lt
    }
  }
  gc()
  message('Completed molecular contrasts: ', cell)
}
coverage <- rbindlist(coverage_parts)
agreement <- rbindlist(sparse_agreement)
agreement[, pass_in_cell := landmark_genes >= 5L & informative_conditions >= 5L &
  is.finite(r) & r >= .7 & is.finite(slope) & slope >= .5 & slope <= 2 &
  is.finite(direction_agreement) & direction_agreement >= .75]
eligibility <- agreement[, .(cells_checked = .N, eligible_landmark = all(pass_in_cell),
  minimum_r = min(r), minimum_sign_agreement = min(direction_agreement), minimum_landmarks = min(landmark_genes)),
  by = .(bank, feature, feature_id)]
eligibility[, rule := 'Both cells: >=5 landmarks; >=5 full-response conditions |delta|>=0.05; r>=0.7; slope0.5-2; sign>=75%']
fwrite(eligibility[, .(representation = bank, feature, eligible_LINCS = eligible_landmark,
  minimum_r, minimum_sign_agreement, minimum_landmarks, rule)],
  file.path(root, 'data/representation_eligibility.tsv'), sep = '\t')
out(coverage, 'coverage')
out(agreement, 'landmark_agreement')
out(agreement[!is.finite(r)], 'unavailable_landmark_correlations')
out(eligibility, 'landmark_eligibility')
out(rbindlist(replicate_agreement), 'replicate_agreement')
scores <- rbindlist(all_scores)
out(scores, 'scores')
out(scores[measurement == 'full', .(mean_displacement = mean(value),
  well_SE = sd(value) / sqrt(.N), RNA_wells = .N),
  by = .(cell_line, bank, feature_id, arm, time_h)], 'condition_scores')
fc <- rbindlist(all_contrasts)
out(fc, 'feature_interactions')
gs <- rbindlist(gene_summaries)
out(gs, 'gene_interaction_counts')

# Directly locate cancellation within empirical programs. Every detected member
# contributes; sign heterogeneity is not hidden by retaining only favorable genes.
coherence <- list()
for (cell in names(contrast_gene_list)) {
  genes <- contrast_gene_list[[cell]]
  joined <- merge(genes, members, by = 'gene', allow.cartesian = TRUE)
  zz <- joined[, .(detected_members = .N, clear_members = sum(clear_interaction),
    positive_clear = sum(clear_interaction & effect > 0), negative_clear = sum(clear_interaction & effect < 0),
    mean_absolute_gene_interaction = mean(abs(effect)), mean_gene_interaction = mean(effect),
    cancellation_fraction = 1 - abs(mean(effect)) / mean(abs(effect))),
    by = .(cell_line, arm, time_h, contrast_id, bank, feature, feature_id)]
  zz[, dominant_direction_fraction := pmax(positive_clear, negative_clear) / pmax(1L, clear_members)]
  coherence[[cell]] <- merge(zz, fc[cell_line == cell, .(feature_id, contrast_id,
    program_interaction = effect, program_q = q_value, program_SE = standard_error)], by = c('feature_id', 'contrast_id'))
}
coherence <- rbindlist(coherence)
out(coherence, 'within_program_cancellation')

# Independently measured nuclei counts retain original percent-vehicle units.
fn <- fread(file.path(prior, 'diaz_functional_condition_summaries.tsv.gz'))[assay == 'InCell']
out(fn, 'independent_nuclei_counts')
count_plot <- ggplot(fn[time_h <= 24], aes(time_h, value, color = arm, group = arm)) +
  geom_hline(yintercept = 100, color = 'grey75') + geom_line() + geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .35, na.rm = TRUE) +
  facet_wrap(~cell_line, nrow = 1) + theme_classic(base_size = 10) +
  labs(x = 'Hours after treatment', y = 'Counted nuclei (% vehicle)',
    title = 'Independent functional measurements', color = 'Treatment',
    caption = 'Source means ± SD; separate wells from RNA. No functional value selected a gene or program.')
gene_plot <- ggplot(gs[arm != 'TMW'], aes(time_h, clear_genes, color = arm, group = arm)) +
  geom_line() + geom_point(size = 1.5) + facet_wrap(~cell_line, nrow = 1) + theme_classic(base_size = 10) +
  labs(x = 'Hours after treatment', y = 'Genes with combination interaction', color = 'Combination',
    title = 'Joint RNA responses depart from the sum of single responses',
    caption = 'Interaction = AB − A − B + vehicle, log2(1+CP10K). |effect|≥0.25 and within-well BH q≤0.05.\nTriplicates quantify the deposited wells; independent biological replication is unresolved.')
plot_save(count_plot / gene_plot, 'function_and_gene_interactions', 10, 8)

p <- ggplot(agreement, aes(bank, r, color = cell_line)) +
  geom_hline(yintercept = .7, linetype = 2, color = 'grey50') +
  geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .5, seed = 1), size = 1.8, na.rm = TRUE) +
  theme_classic(base_size = 11) + labs(x = 'Frozen representation', y = 'Full RNA versus 978-gene score: Pearson r',
    title = 'Do directly measured LINCS genes preserve each program?', color = 'Cell',
    caption = 'One point per program/readout; correlations across treatment/time condition means.\nEligibility additionally requires ≥5 landmarks, slope 0.5–2 and ≥75% direction agreement in both cells.')
plot_save(p, 'landmark_preservation', 8, 5)

cancel_plot <- coherence[time_h == 24 & arm != 'TMW' & clear_members >= 10]
cancel_cor <- cancel_plot[, .(label = sprintf('Pearson r = %.2f',
  cor_safe(dominant_direction_fraction, cancellation_fraction))), by = .(cell_line, arm)]
p <- ggplot(cancel_plot,
  aes(dominant_direction_fraction, cancellation_fraction, color = bank)) +
  geom_point(size = 1.6, alpha = .75) + geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x,
    se = FALSE, color = 'grey35') +
  geom_text(data = cancel_cor, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.25, size = 2.7) +
  facet_grid(cell_line ~ arm) + theme_classic(base_size = 10) + theme(aspect.ratio = 1) +
  labs(x = 'Fraction of clear member genes sharing the dominant sign', y = 'Fraction cancelled by averaging',
    title = 'Gene directions and cancellation within each frozen program',
    caption = '24 h; programs with ≥10 clear member genes. Cancellation = 1 − |mean gene effect| / mean |gene effect|.\nThis diagnoses averaging; it does not establish that those genes cause the functional response.')
plot_save(p, 'gene_cancellation', 11, 7)

# The higher-dose experiment is a separate batch. Check the repeated, identical
# single-drug conditions before interpreting stronger-dose controls as causal.
dose <- qread(file.path(prior, 'diaz_dose_expression.qs'))
dose_md <- as.data.table(dose$samples)
dose_x <- dose$displacement_CP10K
stopifnot(identical(colnames(dose_x), dose_md$sample_id))
old <- qread(file.path(prior, 'diaz_MCF7_molecular.qs'))
old_md <- as.data.table(old$samples)
dose_fn <- fread(file.path(prior, 'diaz_dose_matched_function.tsv'))
cross_batch <- list()
dose_readouts <- list()
for (bn in names(weights)) {
  ww <- weights[[bn]]
  gg <- intersect(rownames(ww), rownames(dose_x))
  ww <- ww[gg, , drop = FALSE]
  ds <- as.matrix(crossprod(ww, dose_x[gg, , drop = FALSE])) / colSums(ww)
  lt <- as.data.table(as.table(ds))
  setnames(lt, c('feature', 'sample_id', 'value'))
  lt[, bank := bn]
  dose_readouts[[bn]] <- merge(lt, dose_md[, .(sample_id, drug, dose_uM, condition_id)], by = 'sample_id')
}
dose_readouts <- rbindlist(dose_readouts)
out(dose_readouts, 'higher_dose_readouts')
common <- intersect(rownames(old$displacement), rownames(dose_x))
detect <- rowSums(old$log2_one_plus_TMM_CP10K[common, ] >= log2(1.1)) >= 3L
common <- common[detect]
for (drug_name in c('Tamoxifen', 'Mefloquine')) {
  arm_code <- if (drug_name == 'Tamoxifen') 'T' else 'M'
  dose_value <- if (drug_name == 'Tamoxifen') 20 else 10
  aa <- rowMeans(old$displacement[common, old_md$arm == arm_code & old_md$time_h == 24, drop = FALSE])
  bb <- rowMeans(dose_x[common, dose_md$drug == drug_name & dose_md$dose_uM == dose_value, drop = FALSE])
  cross_batch[[drug_name]] <- data.table(drug = drug_name, dose_uM = dose_value, time_h = 24,
    genes = length(common), gene_displacement_r = cor_safe(aa, bb),
    original_nuclei_pct = fn[cell_line == 'MCF7' & arm == arm_code & time_h == 24, value],
    later_nuclei_pct = dose_fn[drug == drug_name & dose_uM == dose_value, mean],
    interpretation = 'Separate batch about two years later; same-dose functional response differs; no pooled significance test')
}
cross_batch <- rbindlist(cross_batch)
out(cross_batch, 'higher_dose_batch_check')
dr <- dose_readouts[bank == 'biology' & drug != 'DMSO', .(value = mean(value), se = sd(value) / sqrt(.N)),
  by = .(feature, drug, dose_uM)]
reference <- scores[cell_line == 'MCF7' & bank == 'biology' & measurement == 'full' &
  arm == 'TM' & time_h == 24, .(value = mean(value)), by = feature_id]
reference[, feature := sub('^biology__', '', feature_id)]
p <- ggplot(dr, aes(dose_uM, value, color = drug, group = drug)) +
  geom_hline(data = reference, aes(yintercept = value), linetype = 2, color = 'grey40') +
  geom_line() + geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = .3) +
  facet_wrap(~feature, scales = 'free_y', ncol = 3) + theme_classic(base_size = 9) +
  labs(x = 'Single-drug dose (uM)', y = 'Mean gene displacement: log2(1+CP10K)', color = 'Single drug',
    title = 'Higher doses can be compared descriptively; the batches are not equivalent',
    caption = 'Dashed: original tamoxifen20 + mefloquine10 combination at24h. Points: later dose batch mean +/- well SE.\nSame-dose nuclei counts differ across batches: T20 108% vs68%; M10 118% vs86%. These are not matched functional controls.')
plot_save(p, 'higher_dose_batch_control', 11, 8)

# Compact fixed-readout heatmap; all nine previously defined sets are shown.
hh <- fc[bank == 'biology' & arm != 'TMW' & time_h %in% c(3, 6, 12, 24)]
hh[, column := paste(cell_line, arm, time_h, sep = '__')]
wide <- dcast(hh, feature ~ column, value.var = 'effect')
hm <- as.matrix(wide[, -1])
rownames(hm) <- wide$feature
ca <- unique(hh[, .(column, cell_line, arm, time_h)])
setorder(ca, cell_line, arm, time_h)
hm <- hm[, ca$column, drop = FALSE]
ro <- hclust(dist(hm), method = 'ward.D2')$order
hm <- hm[ro, , drop = FALSE]
stopifnot(identical(colnames(hm), ca$column))
lim <- max(abs(hm))
readout_labels <- c(Immediate_early = 'Immediate early', ERK_output = 'ERK output',
  MTORC1_SIGNALING = 'mTORC1 response', APOPTOSIS = 'Apoptosis genes',
  EPITHELIAL_MESENCHYMAL_TRANSITION = 'EMT genes', GLYCOLYSIS = 'Glycolysis',
  UNFOLDED_PROTEIN_RESPONSE = 'Unfolded protein response', E2F_TARGETS = 'E2F targets',
  MYC_TARGETS_V1 = 'MYC targets')
ht <- draw_heatmap(hm, col_anno_df = as.data.frame(ca), col_anno_cols = c('cell_line', 'arm', 'time_h'),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = FALSE,
  scale_rows = FALSE, value_range = c(-lim, 0, lim), value_colors = c('#286B9E', 'white', '#B43B32'),
  column_split = factor(paste(ca$cell_line, ca$arm), levels = unique(paste(ca$cell_line, ca$arm))),
  row_labels = readout_labels[rownames(hm)], row_names_gp = grid::gpar(fontsize = 9),
  title = 'Fixed biological readouts: combination interaction (AB − A − B + vehicle)',
  hm_name = 'Mean gene\ninteraction')
cairo_pdf(file.path(root, 'figures/representation_biological_interactions.pdf'), width = 10, height = 4.5)
ComplexHeatmap::draw(ht)
dev.off()
png(file.path(root, 'figures/representation_biological_interactions.png'), width = 1400, height = 630, res = 140)
ComplexHeatmap::draw(ht)
dev.off()

summary_rows <- eligibility[, .(features = .N, eligible_landmark = sum(eligible_landmark)), by = bank]
out(summary_rows, 'eligibility_summary')
out(coherence[time_h == 24 & arm != 'TMW', .(programs = .N,
  programs_with_10_clear_genes = sum(clear_members >= 10),
  mixed_sign_programs = sum(clear_members >= 10 & dominant_direction_fraction < .8),
  cancelled_half = sum(clear_members >= 10 & cancellation_fraction > .5)), by = .(cell_line, arm, bank)],
  'cancellation_summary')
writeLines(c('Molecular representation audit; all functional outcomes excluded from molecular selection.',
  'Banks16,42,43 are frozen historical outcome-informed banks; biology is the nine prior independent readouts.',
  'Full RNA scores are unscaled mean log2(1+TMM CP10K) differences against matched-time vehicle.',
  'Pair interaction is AB-A-B+vehicle on this molecular scale. It is not a functional Bliss score.',
  'All pair/time contrasts are reported; the LNCaP triple remains a separate three-way boundary test.',
  'RNA triplicates quantify deposited-well agreement; biological independence is unresolved.',
  'Landmark eligibility is frozen from molecular-only agreement in both Diaz cells; functional labels never enter.',
  capture.output(print(summary_rows)),
  'Molecular interaction effect size is AB-A-B+vehicle, not a log fold change and not functional synergy.',
  'At24h TM extra UPR mean response: MCF7 +0.363 (51 positive /10 negative clear genes), LNCaP +0.152 (28/9).',
  'E2F interaction is positive in TM, but its actual displacement is negative: suppression is LESS THAN ADDITIVE.',
  'This distinction prevents interpreting a positive E2F interaction as increased proliferation.',
  'Large changes within a program can have opposite signs and cancel. Cancellation does not itself prove lost causal information.',
  'Time-zero pair controls: MCF7 all0 clear genes; LNCaP TM/TW0, MW131. Full counts at all times are exported.',
  'Higher-dose controls were collected in a different batch about two years later and are not pooled for significance.',
  capture.output(print(cross_batch)),
  'Same-dose gene response correlations are high across batches, while nuclei counts differ materially.',
  'Therefore this later dose batch cannot decisively distinguish composition biology from stronger single exposure.'),
  file.path(root, 'REPRESENTATION_INTERPRETATION.txt'))
out(data.table(path = '/work/src/program_pilot/60_representation_audit.R',
  sha256 = digest(file = '/work/src/program_pilot/60_representation_audit.R', algo = 'sha256')), 'script_hash')
capture.output(sessionInfo(), file = file.path(root, 'representation_sessionInfo.txt'))
message('Completed representation audit.')
