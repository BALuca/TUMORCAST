#!/usr/bin/env Rscript
# Reuse the original grant bank exactly; no gene/program selection on outcomes.
# Frozen signed weights are retained, including on partially measured assays.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(digest)
  library(jsonlite)
  library(xml2)
  library(ggplot2)
  library(patchwork)
})
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
original <- '/work/almanac_conditionals/state_conditioned_bridge/program_first'
prior59 <- '/work/results/program_pilot/59_mechanistic_evidence_datasets/data'
prior60 <- '/work/results/program_pilot/60_defensible_response_composition/data'
for (folder in c('data', 'models', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
sources <- c(weights = file.path(original, 'results/frozen_program_weights.tsv'),
  freeze = file.path(original, 'results/program_bank_freeze.json'),
  terms = file.path(original, 'results/program_theme_top_terms.tsv'),
  themes = file.path(original, 'bridge/inspection/results/selected_program_themes.tsv'),
  gene_summary = file.path(original, 'results/program_gene_summary.tsv'),
  historical_prepared = file.path(original, 'all_programs_prepared.rds'),
  historical_genes = '/work/almanac_conditionals/state_conditioned_bridge/sensitivity_selected/artifacts/eligible_gene_measurements.rds',
  biological_sets = '/work/results/program_pilot/57_erbb_met_hgf_program_biology/data/biological_readout_membership.tsv.gz',
  landmark = file.path(prior59, 'matched_gene_annotation.tsv'),
  h1792 = file.path(prior60, 'h1792_molecular_replication.qs'),
  grant = '/work/ERC_StG_2027_B1_v20.docx')
stopifnot(all(file.exists(sources)))
freeze <- fromJSON(sources[['freeze']])
hash <- digest(file = sources[['weights']], algo = 'sha256')
stopifnot(identical(hash, freeze[['frozen_program_weights.tsv']]), freeze$rank == 40L,
  identical(freeze$sensitivity_used, FALSE), identical(freeze$AU565_used, FALSE))
membership <- fread(sources[['weights']])
programs <- sprintf('P%02d', 1:40)
genes <- sort(unique(membership$gene))
stopifnot(length(genes) == 8842L, setequal(membership$program, programs), !anyDuplicated(membership$gene),
  all(is.finite(membership$weight)))
W <- sparseMatrix(i = match(membership$gene, genes), j = match(membership$program, programs),
  x = membership$weight, dims = c(length(genes), length(programs)), dimnames = list(genes, programs))
qsave(W, file.path(root, 'models/original40_weights.qs'))
write_table(membership, 'original40_membership')

# The program identifiers and annotations come from the original frozen bank.
# Previously curated labels take precedence; remaining labels use its original
# all-eligible-background enrichment, without rerunning or selecting on outcomes.
themes <- fread(sources[['themes']])[, .(program, curated_theme = theme_label, curated_q = theme_q)]
terms <- fread(sources[['terms']])[background == 'all_eligible']
top <- terms[order(q, p), .SD[1L], by = program]
annotations <- merge(data.table(program = programs), themes, by = 'program', all.x = TRUE)
annotations <- merge(annotations, top[, .(program, source_term = term, source_collection = collection,
  source_q = q, source_overlap = overlap)], by = 'program', all.x = TRUE)
annotations <- merge(annotations, fread(sources[['gene_summary']]), by = 'program', all.x = TRUE)
annotations[, theme := fifelse(!is.na(curated_theme), curated_theme,
  fifelse(is.finite(source_q) & source_q <= .05, gsub('_', ' ', sub('^(REACTOME_|HALLMARK_)', '', source_term)),
    'No clear enrichment'))]
annotations[, display_name := paste(program, theme, sep = ': ')]
annotations[, annotation_source := fifelse(!is.na(curated_theme), sources[['themes']], sources[['terms']])]
write_table(annotations, 'original40_program_annotations')

# Reproduce the actual archived 40-program scores, not a newly reconstructed
# interpretation of their weights. This checks signs, magnitude and orientation.
old <- readRDS(sources[['historical_prepared']])
old_gene <- readRDS(sources[['historical_genes']])
stopifnot(all(genes %in% colnames(old_gene$response)), all(genes %in% colnames(old_gene$baseline)))
replayed <- as.matrix(old_gene$response[, genes, drop = FALSE] %*% W)
row_index <- match(old$m$sample_id, old_gene$metadata$sample_id)
stopifnot(!anyNA(row_index))
stored <- as.matrix(old$m[, programs, with = FALSE])
maximum_error <- max(abs(replayed[row_index, , drop = FALSE] - stored))
weight_error <- max(abs(as.matrix(W) - old$W[genes, programs]))
stopifnot(maximum_error < 1e-8, weight_error < 1e-12)
write_table(data.table(check = c('Exact archived weights', 'Archived single-drug projections'),
  maximum_absolute_difference = c(weight_error, maximum_error),
  rows = c(length(genes), nrow(stored)), columns = 40L), 'original40_historical_replay')
rm(old, old_gene, replayed, stored)
gc()

# Preserve visible grant text and exact source provenance, excluding Word's
# embedded Zotero field instructions from the extraction.
doc <- read_xml(unz(sources[['grant']], 'word/document.xml'))
paragraphs <- xml_find_all(doc, './/w:p', xml_ns(doc))
visible <- vapply(paragraphs, function(p) paste(xml_text(xml_find_all(p, './/w:t', xml_ns(doc))), collapse = ''), character(1))
selected <- grep('40 response programs|441 cell lines|baseline-dependent prediction|Removing untreated', visible,
  ignore.case = TRUE)
write_table(data.table(paragraph = selected, text = visible[selected]), 'original40_grant_claim')
write_table(data.table(source = names(sources), path = unname(sources),
  sha256 = vapply(sources, digest, character(1), file = TRUE, algo = 'sha256')), 'original40_source_manifest')

# Full and subset projections use exactly the same W. Missing genes contribute
# no observed term; they are NOT assumed biologically unchanged. Do not rescale
# the surviving weights to make a partial measurement look like a complete one.
landmarks <- unique(fread(sources[['landmark']])$gene_symbol)
stopifnot(length(landmarks) == 978L)
project <- function(delta, dataset, supplied_genes = rownames(delta)) {
  stopifnot(is.matrix(delta), !anyDuplicated(rownames(delta)), !anyDuplicated(colnames(delta)),
    all(is.finite(delta)))
  measured <- intersect(intersect(genes, supplied_genes), rownames(delta))
  stopifnot(all(measured %in% rownames(delta)), length(measured) > 0L)
  ww <- W[measured, , drop = FALSE]
  scores <- t(as.matrix(crossprod(ww, delta[measured, , drop = FALSE])))
  coverage <- data.table(dataset = dataset, program = programs,
    total_genes = colSums(W != 0), measured_genes = colSums(ww != 0),
    absolute_weight_fraction = colSums(abs(ww)) / colSums(abs(W)),
    squared_weight_fraction = colSums(ww ^ 2) / colSums(W ^ 2),
    positive_genes_measured = colSums(ww > 0), negative_genes_measured = colSums(ww < 0))
  coverage[, gene_fraction := measured_genes / total_genes]
  coverage[, measurement := fifelse(measured_genes == total_genes, 'Complete frozen projection',
    'Partial measured contribution; unmeasured contribution unknown')]
  list(scores = scores, coverage = coverage)
}
cor_safe <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 4L || sd(a[ok]) < 1e-10 || sd(b[ok]) < 1e-10) return(NA_real_)
  cor(a[ok], b[ok])
}
coverage_list <- score_list <- fidelity_list <- list()
for (cell in c('MCF7', 'LNCaP')) {
  obj <- qread(file.path(prior59, paste0('diaz_', cell, '_molecular.qs')))
  delta <- obj$displacement
  md <- as.data.table(obj$samples)
  full <- project(delta, paste0('Diaz_', cell, '_full'))
  limited <- project(delta, paste0('Diaz_', cell, '_978'), landmarks)
  coverage_list[[paste0(cell, '_full')]] <- full$coverage
  coverage_list[[paste0(cell, '_978')]] <- limited$coverage
  for (mode in c('full', '978')) {
    ss <- if (mode == 'full') full$scores else limited$scores
    lt <- as.data.table(as.table(ss))
    setnames(lt, c('sample_id', 'program', 'score'))
    lt[, `:=`(dataset = paste0('Diaz_', cell), measurement = mode)]
    score_list[[paste(cell, mode)]] <- merge(lt, md[, .(sample_id, condition_id, cell_line, arm, time_h)], by = 'sample_id')
  }
  cc <- unique(md[arm != 'DMSO' & time_h > 0 & time_h <= 24, condition_id])
  aa <- t(vapply(cc, function(id) colMeans(full$scores[md$condition_id == id, , drop = FALSE]), numeric(40)))
  bb <- t(vapply(cc, function(id) colMeans(limited$scores[md$condition_id == id, , drop = FALSE]), numeric(40)))
  fidelity_list[[cell]] <- rbindlist(lapply(seq_along(programs), function(j) {
    active <- abs(aa[, j]) >= .1
    data.table(dataset = paste0('Diaz_', cell), program = programs[j], independent_conditions = length(cc),
      full_subset_r = cor_safe(aa[, j], bb[, j]), slope_subset_on_full = coef(lm(bb[, j] ~ aa[, j]))[2],
      active_conditions = sum(active), sign_agreement = if (sum(active) > 0L)
        mean(sign(aa[active, j]) == sign(bb[active, j])) else NA_real_,
      median_absolute_full = median(abs(aa[, j])), median_absolute_subset = median(abs(bb[, j])))
  }))
}

# H1792 has two drug conditions with replicate wells: insufficient independent
# conditions for a per-program transport correlation. Export both actual arms.
h1792 <- qread(sources[['h1792']])
for (source in c('RNA', 'Tahoe')) {
  delta <- h1792[[paste0(source, '_delta')]]
  full <- project(delta, paste0('H1792_', source, '_full'))
  limited <- project(delta, paste0('H1792_', source, '_978'), landmarks)
  coverage_list[[paste(source, 'full')]] <- full$coverage
  coverage_list[[paste(source, '978')]] <- limited$coverage
  md <- as.data.table(h1792[[paste0(source, '_metadata')]])
  identifier <- if (source == 'RNA') 'sample_id' else 'observation_id'
  stopifnot(identical(rownames(full$scores), md[[identifier]]))
  for (mode in c('full', '978')) {
    ss <- if (mode == 'full') full$scores else limited$scores
    lt <- as.data.table(as.table(ss))
    setnames(lt, c('profile_id', 'program', 'score'))
    lt[, measurement := mode]
    md2 <- copy(md)
    setnames(md2, identifier, 'profile_id')
    write_table(merge(lt, md2, by = 'profile_id'), paste0('original40_H1792_', source, '_', mode, '_scores'))
  }
}
write_table(rbindlist(coverage_list), 'original40_external_coverage')
fwrite(rbindlist(score_list), file.path(root, 'data/original40_Diaz_scores.tsv.gz'), sep = '\t')
fidelity <- rbindlist(fidelity_list)
fidelity[, interpretation := 'Molecular subset audit only; all40 programs remain in the primary model']
write_table(fidelity, 'original40_full_vs_landmark_fidelity')
message('Original bank identity, replay and external molecular projection complete.')

# Score the actual new combination-source measurements. These are L1000
# fluorescence responses, not RNA-seq CP10K; weights transfer, absolute units do
# not. No held combination outcome is loaded anywhere in this script.
input_path <- file.path(root, 'data/matched_molecular.qs')
stopifnot(file.exists(input_path))
molecular <- qread(input_path)
delta <- molecular$displacement
md <- as.data.table(molecular$metadata)
stopifnot(identical(colnames(delta), md$profile_id), nrow(delta) == 978L,
  all(rownames(delta) %in% landmarks), all(is.finite(delta)))
projection <- project(delta, 'BT20_measured_L1000')
biological_members <- unique(fread(sources[['biological_sets']])[, .(feature = readout, gene)])
biological_features <- sort(unique(biological_members$feature))
biological <- matrix(NA_real_, nrow(md), length(biological_features),
  dimnames = list(md$profile_id, biological_features))
biology_coverage <- list()
for (j in seq_along(biological_features)) {
  feature_name <- biological_features[j]
  gg <- biological_members[feature == feature_name, gene]
  available <- intersect(gg, rownames(delta))
  stopifnot(length(available) > 0L)
  biological[, j] <- colMeans(delta[available, , drop = FALSE])
  biology_coverage[[feature_name]] <- data.table(feature = feature_name, total_genes = length(gg),
    measured_genes = length(available), measured_fraction = length(available) / length(gg),
    measured_members = paste(available, collapse = ';'),
    interpretation = 'Mean of directly measured gene members; transcriptional readout, not kinase or protein activity')
}
write_table(projection$coverage, 'original40_BT20_coverage')
write_table(rbindlist(biology_coverage), 'biological_BT20_coverage')
write_table(biological_members, 'biological_membership')
representation <- list(original40 = projection$scores, biology = biological, metadata = md,
  coverage = projection$coverage, biology_coverage = rbindlist(biology_coverage),
  program_annotations = annotations, original_weight_sha256 = hash,
  units = 'Weighted log2(1+quantile-normalized measured L1000 fluorescence) treated-minus-own-plate-vehicle; partial measured-gene projection, not original CP10K units',
  gene_displacement_source = input_path,
  primary_rule = 'All40 original signed programs retained; no outcome-based selection or eligibility filter')
stopifnot(identical(rownames(representation$original40), md$profile_id),
  identical(rownames(representation$biology), md$profile_id), ncol(representation$original40) == 40L,
  all(is.finite(representation$original40)), all(is.finite(representation$biology)))
qsave(representation, file.path(root, 'data/representation_scores.qs'))
for (bn in c('original40', 'biology')) {
  lt <- as.data.table(as.table(representation[[bn]]))
  setnames(lt, c('profile_id', 'feature', 'score'))
  lt[, representation := bn]
  fwrite(merge(lt, md, by = 'profile_id'), file.path(root, 'data', paste0('representation_', bn, '_scores.tsv.gz')), sep = '\t')
}
write_table(data.table(path = c(input_path, '/work/src/program_pilot/61_original_bank_projection.R'),
  sha256 = vapply(c(input_path, '/work/src/program_pilot/61_original_bank_projection.R'),
    digest, character(1), file = TRUE, algo = 'sha256')), 'original40_run_manifest')
writeLines(c('Original grant40: frozen Tahoe/CytoSig signed orthogonalNMF bank,8842genes. No sensitivity used in discovery.',
  paste('OriginalweightSHA256:', hash),
  sprintf('Archived projection reproduced across270profilesx40programs; maximum discrepancy %.3g.', maximum_error),
  'Original scoring uses corrected log2(1+CP10K) treated-minus-reference, multiplied by unchanged signed weights.',
  'All40 programs retained in the primary experiment. Missing members are unknown, not imputed unchanged.',
  'BT20 L1000 measurements have different fluorescence units; partial weighted scores are calibrated only within that assay.',
  'FullRNA versus978 molecular audit reports correlation, sign agreement and amplitude attenuation separately; no outcome-guided feature filtering.',
  'Biological fallback uses all9 prior fixed readouts, recording their measured members; no inferred phosphorylation.',
  'Program labels are original identifiers plus original enrichment/curated annotations; no new outcome-selected biological story.'),
  file.path(root, 'ORIGINAL40_PROJECTION.txt'))
capture.output(sessionInfo(), file = file.path(root, 'original40_sessionInfo.txt'))
message('Completed all40 original-bank and independent biological projections for ', nrow(md), ' molecular profiles.')

# Molecular diagnostics only: the figures below never load combination outcomes.
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0('original40_', name, '.pdf')), p,
    device = cairo_pdf, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0('original40_', name, '.png')), p,
    width = width, height = height, dpi = 130)
}
plot_names <- setNames(paste(annotations$program, strtrim(annotations$theme, 33), sep = ': '), annotations$program)
fidelity_plot <- merge(fidelity, projection$coverage[, .(program, measured_genes, absolute_weight_fraction)], by = 'program')
fidelity_plot[, label := factor(plot_names[program], levels = rev(plot_names[programs]))]
fidelity_plot[, dataset := sub('Diaz_', '', dataset)]
fidelity_long <- melt(fidelity_plot, id.vars = c('program', 'label', 'dataset'),
  measure.vars = c('full_subset_r', 'slope_subset_on_full'), variable.name = 'metric', value.name = 'value')
fidelity_long[, metric := factor(metric, levels = c('full_subset_r', 'slope_subset_on_full'),
  labels = c('Shape: full-to-978 Pearson correlation', 'Amplitude: 978-on-full regression slope'))]
p <- ggplot(fidelity_long, aes(value, label, color = dataset)) +
  geom_vline(xintercept = 0, color = 'grey80') +
  geom_point(position = position_dodge(width = .45), size = 1.5) +
  facet_wrap(~metric, scales = 'free_x', nrow = 1) + theme_classic(base_size = 8.5) +
  theme(legend.position = 'bottom') +
  labs(x = 'Molecular measurement agreement', y = NULL, color = 'Independent RNA experiment',
    title = 'Original grant bank: all 40 programs, full RNA versus directly measured LINCS genes',
    caption = 'Same immutable signed weights in both projections; missing members are not renormalized.\nCondition means across doses/times; no function used. Good shape agreement can coexist with strong amplitude loss. All40 remain in the primary model.')
save_plot(p, 'full_to_landmark_fidelity', 11, 10)

# Compare the two studies on the SAME measured genes so unequal gene coverage
# cannot itself create apparent disagreement between their program responses.
h <- qread(sources[['h1792']])
rna_md <- as.data.table(h$RNA_metadata)
tahoe_md <- as.data.table(h$Tahoe_metadata)
common <- Reduce(intersect, list(genes, rownames(h$RNA_delta), rownames(h$Tahoe_delta)))
rna_scores <- project(h$RNA_delta, 'H1792_RNA_common', common)$scores
tahoe_scores <- project(h$Tahoe_delta, 'H1792_Tahoe_common', common)$scores
treated <- rna_md$treatment == 'trametinib'
vehicle <- rna_md$treatment == 'DMSO'
stopifnot(sum(treated) == 3L, sum(vehicle) == 3L)
rna_mean <- colMeans(rna_scores[treated, , drop = FALSE]) - colMeans(rna_scores[vehicle, , drop = FALSE])
rna_se <- sqrt(apply(rna_scores[treated, , drop = FALSE], 2, var) / sum(treated) +
  apply(rna_scores[vehicle, , drop = FALSE], 2, var) / sum(vehicle))
erk_genes <- c('PHLDA1', 'SPRY2', 'DUSP6', 'ETV4', 'ETV5')
erk_membership <- membership[gene %in% erk_genes]
stopifnot(setequal(erk_membership$gene, erk_genes), all(erk_genes %in% common))
erk_programs <- sort(unique(erk_membership$program))
erk_membership[, directly_measured_LINCS := gene %in% landmarks]
write_table(erk_membership, 'original40_ERK_gene_membership')
comparison <- rbindlist(lapply(which(tahoe_md$dose == .5), function(i) {
  data.table(program = programs, independent_RNA = rna_mean, independent_RNA_SE = rna_se,
    Tahoe = as.numeric(tahoe_scores[i, programs]), formulation = tahoe_md$formulation[i],
    RNA_time_h = 16, Tahoe_time_h = 24, dose_uM = .5,
    source_profile = tahoe_md$observation_id[i], common_genes = length(common),
    ERK_gene_member = programs %in% erk_programs)
}))
write_table(comparison, 'original40_H1792_common_gene_replication')
agreement <- comparison[, .(programs = .N, common_genes = common_genes[1],
  program_response_r = cor_safe(independent_RNA, Tahoe),
  direction_agreement = mean(sign(independent_RNA) == sign(Tahoe))), by = formulation]
write_table(agreement, 'original40_H1792_replication_summary')
agreement[, text := sprintf('All40: r = %.2f\nSame sign: %d/40', program_response_r, round(direction_agreement * 40))]
p <- ggplot(comparison, aes(independent_RNA, Tahoe)) +
  geom_hline(yintercept = 0, color = 'grey80') + geom_vline(xintercept = 0, color = 'grey80') +
  geom_point(aes(color = ERK_gene_member), size = 2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey35') +
  geom_text(data = comparison[ERK_gene_member == TRUE], aes(label = program),
    hjust = -.15, vjust = -.6, size = 2.6, check_overlap = TRUE) +
  geom_text(data = agreement, aes(x = -Inf, y = Inf, label = text), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.2, size = 3) +
  facet_wrap(~formulation, scales = 'free', nrow = 1) +
  scale_color_manual(values = c('FALSE' = '#777777', 'TRUE' = '#B34D24'),
    labels = c('Other original programs', 'Contains a prespecified ERK-output gene')) +
  theme_classic(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom') +
  labs(x = 'Independent H1792 RNA response, 16h', y = 'Tahoe H1792 RNA response, 24h', color = NULL,
    title = 'Trametinib 0.5 uM: independent replication across the complete original bank',
    caption = paste0('Both projections use the same ', length(common), ' measured genes and unchanged signed weights.\nRNA-seq triplicates versus Tahoe pooled response; collection time, culture and formulation differ. This is molecular replication, not functional prediction.'))
save_plot(p, 'H1792_replication_all40', 9, 5.5)

# Trace every prespecified ERK gene into its assigned program. Positive program
# values are not synonymous with pathway activation because weights are signed.
rna_gene_mean <- rowMeans(h$RNA_delta[erk_genes, treated, drop = FALSE]) -
  rowMeans(h$RNA_delta[erk_genes, vehicle, drop = FALSE])
rna_gene_se <- sqrt(apply(h$RNA_delta[erk_genes, treated, drop = FALSE], 1, var) / sum(treated) +
  apply(h$RNA_delta[erk_genes, vehicle, drop = FALSE], 1, var) / sum(vehicle))
gene_trace <- data.table(gene = erk_genes, response = rna_gene_mean, response_SE = rna_gene_se,
  source = 'Independent RNA16h')
for (i in which(tahoe_md$dose == .5)) {
  gene_trace <- rbind(gene_trace, data.table(gene = erk_genes,
    response = h$Tahoe_delta[erk_genes, i], response_SE = NA_real_,
    source = paste('Tahoe24h', tahoe_md$formulation[i])))
}
gene_trace <- merge(gene_trace, erk_membership[, .(gene, program, weight, directly_measured_LINCS)], by = 'gene')
gene_trace[, signed_contribution := response * weight]
write_table(gene_trace, 'original40_H1792_ERK_gene_trace')
source_names <- unique(gene_trace$source)
component_rows <- list()
for (src in source_names) {
  for (pr in erk_programs) {
    gene_part <- gene_trace[source == src & program == pr, sum(signed_contribution)]
    full_part <- if (src == 'Independent RNA16h') rna_mean[pr] else {
      fm <- sub('^Tahoe24h ', '', src)
      comparison[formulation == fm & program == pr, Tahoe]
    }
    component_rows[[paste(src, pr)]] <- data.table(source = src, program = pr,
      contribution = c('Prespecified ERK genes', 'All other measured members', 'Full program'),
      value = c(gene_part, full_part - gene_part, full_part))
  }
}
components <- rbindlist(component_rows)
write_table(components, 'original40_H1792_ERK_program_decomposition')
components[, contribution := factor(contribution,
  levels = c('Prespecified ERK genes', 'All other measured members', 'Full program'))]
program_info <- annotations[program %in% erk_programs]
program_labels <- setNames(paste(program_info$program, strtrim(program_info$theme, 24), sep = ': '), program_info$program)
components[, program_label := factor(program_labels[program], levels = program_labels[erk_programs])]
gene_trace[, gene_label := factor(paste0(gene, ' -> ', program, ' [', ifelse(weight > 0, '+', '-'), ' weight]'),
  levels = unique(paste0(gene, ' -> ', program, ' [', ifelse(weight > 0, '+', '-'), ' weight]')))]
p_gene <- ggplot(gene_trace, aes(response, gene_label, color = source)) +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_errorbar(aes(xmin = response - response_SE, xmax = response + response_SE),
    orientation = 'y', position = position_dodge(width = .55), width = .2, na.rm = TRUE) +
  geom_point(position = position_dodge(width = .55), size = 2) + theme_classic(base_size = 10) +
  labs(x = 'Measured gene response: log2(1+CP10K)', y = NULL, color = NULL,
    title = 'Trametinib 0.5 uM: five prespecified ERK-output genes and original assignments') +
  theme(legend.position = 'bottom')
p_component <- ggplot(components, aes(contribution, value, color = source, group = source)) +
  geom_hline(yintercept = 0, color = 'grey75') +
  geom_point(position = position_dodge(width = .5), size = 1.8) +
  facet_wrap(~program_label, scales = 'free_y', nrow = 1) + theme_classic(base_size = 9) +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(color = 'none') +
  labs(x = NULL, y = 'Signed contribution to original program score', color = NULL,
    title = 'The four programs containing those genes: ERK contribution versus other members')
p <- (p_gene / p_component) + plot_layout(heights = c(1, 1.2), guides = 'collect') +
  plot_annotation(caption = 'All five genes were fixed before this experiment; programs selected solely by membership, not outcomes or agreement.\nIndependent RNA error bars are well SE including vehicle variation; Tahoe aggregated profiles have no replicate SE here.\nETV4 has a negative original weight: lower ETV4 gives a positive contribution. Program theme names are enrichment labels, not direct pathway activity.')
p <- p & theme(legend.position = 'bottom')
save_plot(p, 'H1792_ERK_gene_to_program_trace', 12, 8)
message('Completed original-bank measurement and H1792 mechanism diagnostic figures.')
