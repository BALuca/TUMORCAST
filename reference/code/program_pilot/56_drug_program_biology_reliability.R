#!/usr/bin/env Rscript
# Read bank43 as drug biology: repeatability, changing programs, genes and decoder capture.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/56_drug_program_biology_reliability'
prior <- '/work/results/program_pilot/55_current_bank_within_treatment_signal_audit'
bank_root <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs'
panel_root <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data'
prism_root <- '/work/cache/tahoe_prism_full_chain/prism'
focus_drugs <- c('afatinib', 'erlotinib', 'gefitinib', 'lapatinib', 'trametinib',
  'selumetinib', 'cobimetinib', 'binimetinib', 'alpelisib', 'ribociclib',
  'palbociclib', 'bortezomib', 'ixazomib')
for (directory in file.path(root, c('data', 'figures'))) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
cor_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) < 1e-12 || sd(y[ok]) < 1e-12) return(NA_real_)
  cor(x[ok], y[ok])
}
paths <- c(file.path(prior, 'data/diagnostic_predictions.tsv.gz'),
  file.path(bank_root, 'models/program_bank.qs'),
  file.path(panel_root, 'molecular_inputs/panel_index.qs'),
  file.path(prism_root, 'secondary-screen-logfold-change.csv'),
  file.path(prism_root, 'secondary-screen-replicate-treatment-info.csv'),
  '/work/grant_rebuild_2026-10-04/solution_evidence/hallmark_membership.tsv',
  '/work/data/01_curated_pathways/ReactomePathways.gmt.zip',
  '/work/outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv')
hashes <- vapply(paths, function(p) digest(file = p, algo = 'sha256'), character(1))
save_table(data.table(path = paths, sha256 = hashes), 'input_manifest')
writeLines(capture.output(sessionInfo()), file.path(root, 'session_info.txt'))
bank <- qread(paths[2])
panel <- qread(paths[3])
all_rows <- fread(paths[1])
# Use one saved model per drug: MEK exclusion where the drug is seen, otherwise ERBB exclusion.
choice <- unique(all_rows[, .(drug, family)])[order(drug, family != 'MEK')][, .SD[1L], by = drug]
rows <- merge(all_rows, choice, by = c('drug', 'family'))
stopifnot(!anyDuplicated(rows[, .(evaluation, observation_id)]))
md <- unique(rows[, .(observation_id, drug, depmap_id, recipient, dose,
  measured_pct, observed_log2, prism_screen, prism_lower_dose, prism_upper_dose)])
stopifnot(!anyDuplicated(md$observation_id), all(md$observation_id %in% rownames(bank$scores)))
relations <- merge(fread(file.path(prior, 'data/drug_program_relationships.tsv')), choice, by = c('drug', 'family'))
selected <- rbind(relations[, .SD[which.max(abs(pooled_r))], by = drug][, view := 'Full dose curve'],
  relations[, .SD[which.max(abs(within_dose_r))], by = drug][, view := 'Differences between cells at fixed dose'])
save_table(selected, 'development_selected_programs')

# Actual membership enrichment, with the program-bank universe as the background.
hallmark <- fread(paths[6])
connection <- unz(paths[7], 'ReactomePathways.gmt')
gmt <- strsplit(readLines(connection), '\t', fixed = TRUE)
close(connection)
reactome <- rbindlist(lapply(gmt, function(z) data.table(collection = 'Reactome', term = z[1], gene = z[-c(1,2)])))
membership <- unique(rbind(hallmark, reactome)[gene %in% bank$genes])
membership <- membership[, if (.N >= 8L & .N <= 1500L) .SD, by = .(collection, term)]
gene_sets <- split(membership$gene, paste(membership$collection, membership$term, sep = '|'))
enrichment <- rbindlist(lapply(colnames(bank$weights), function(program) {
  genes <- rownames(bank$weights)[as.numeric(bank$weights[, program]) > 0]
  rbindlist(lapply(names(gene_sets), function(key) {
    set <- gene_sets[[key]]
    overlap <- intersect(genes, set)
    data.table(program, pathway = key, program_genes = length(genes), pathway_genes = length(set),
      overlap = length(overlap), expected = length(genes) * length(set) / length(bank$genes),
      p_value = phyper(length(overlap) - 1L, length(set), length(bank$genes) - length(set), length(genes), lower.tail = FALSE),
      genes = paste(sort(overlap), collapse = ', '))
  }))
}))
enrichment[, q_value := p.adjust(p_value, 'BH'), by = program]
enrichment[, fold_enrichment := overlap / expected]
save_table(enrichment, 'program_pathway_enrichment')
save_table(enrichment[overlap >= 3L][order(q_value, -fold_enrichment)][, head(.SD, 5L), by = program], 'program_top_themes')
message('Program biology annotations saved for ', ncol(bank$weights), ' programs')

# Match released replicate columns through the original reviewed broad-ID mapping.
drug_mapping <- unique(fread(file.path(panel_root, 'prism_drug_identity_mapping.tsv'))[, .(broad_id, screen_id, drug)])
stopifnot(!anyDuplicated(drug_mapping[, .(broad_id, screen_id)]))
rep_info <- fread(paths[5])
rep_info <- merge(rep_info, drug_mapping, by = c('broad_id', 'screen_id'))
header <- names(fread(paths[4], nrows = 0L))
rep_info <- rep_info[drug %in% md$drug & column_name %in% header]
stopifnot(!anyDuplicated(rep_info$column_name))
raw <- fread(paths[4], select = c(header[1], rep_info$column_name))
setnames(raw, 1L, 'depmap_id')
raw <- raw[depmap_id %in% md$depmap_id]
stopifnot(!anyDuplicated(raw$depmap_id), all(md$depmap_id %in% raw$depmap_id))
measure_columns <- setdiff(names(raw), 'depmap_id')
raw[, (measure_columns) := lapply(.SD, as.numeric), .SDcols = measure_columns]
long <- melt(raw, id.vars = 'depmap_id', variable.name = 'column_name', value.name = 'rep_log2', variable.factor = FALSE)
long <- merge(long[is.finite(rep_log2)], rep_info[, .(column_name, drug, dose, screen_id, detection_plate, broad_id)], by = 'column_name')
long[, dose_key := sprintf('%.8g', dose)]
low <- merge(md, long, by.x = c('drug', 'depmap_id', 'prism_screen'), by.y = c('drug', 'depmap_id', 'screen_id'),
  allow.cartesian = TRUE, suffixes = c('', '_source'))
low <- low[abs(log(dose_source / prism_lower_dose)) < 1e-7]
high <- merge(md[, .(observation_id, drug, depmap_id, prism_screen, prism_upper_dose)], long,
  by.x = c('drug', 'depmap_id', 'prism_screen'), by.y = c('drug', 'depmap_id', 'screen_id'), allow.cartesian = TRUE)
high <- high[abs(log(dose / prism_upper_dose)) < 1e-7,
  .(observation_id, detection_plate, broad_id, high_log2 = rep_log2)]
replicates <- merge(low, high, by = c('observation_id', 'detection_plate', 'broad_id'))
stopifnot(!anyDuplicated(replicates[, .(observation_id, detection_plate, broad_id)]))
replicates[, fraction := fifelse(prism_lower_dose == prism_upper_dose, 0,
  log(dose / prism_lower_dose) / log(prism_upper_dose / prism_lower_dose))]
replicates[, value := rep_log2 + fraction * (high_log2 - rep_log2)]
replicates[, replicate := sub('.*_X([0-9]+).*', '\\1', detection_plate)]
save_table(replicates[, .(observation_id, drug, depmap_id, recipient, dose, detection_plate,
  broad_id, replicate, value, abundance_pct = 100 * 2^value)], 'matched_prism_replicates')
rep_summary <- replicates[, .(replicates = .N, replicate_mean_log2 = mean(value), replicate_sd_log2 = sd(value),
  replicate_min_pct = min(100 * 2^value), replicate_max_pct = max(100 * 2^value),
  all_replicates_below_50 = all(100 * 2^value <= 50), all_replicates_above_50 = all(100 * 2^value > 50)), by = observation_id]
md <- merge(md, rep_summary, by = 'observation_id', all.x = TRUE)
save_table(md, 'observation_reliability')
save_table(md[, .(observations = .N, with_two_or_more_replicates = sum(replicates >= 2, na.rm = TRUE),
  median_replicate_sd_log2 = median(replicate_sd_log2, na.rm = TRUE)), by = drug], 'drug_replicate_coverage')
save_table(md[is.na(replicates)], 'missing_matched_replicates')
message('PRISM repeat support saved: ', sum(md$replicates >= 2L, na.rm = TRUE), '/', nrow(md), ' observations')

# Across-cell replicate agreement is evaluated within each drug/dose, preserving dose curves separately.
rep_cell <- replicates[, .(value = mean(value)), by = .(drug, dose, depmap_id, replicate)]
rep_pairs <- merge(rep_cell, rep_cell, by = c('drug', 'dose', 'depmap_id'), allow.cartesian = TRUE)
rep_pairs <- rep_pairs[replicate.x < replicate.y]
rep_agreement <- rep_pairs[, .(cells = .N, r = cor_safe(value.x, value.y)), by = .(drug, dose, replicate.x, replicate.y)]
save_table(rep_agreement, 'within_dose_replicate_agreement')

# Independent vehicle pseudo-responses bound control variation of each actual program score.
null_tables <- focus_gene_rows <- list()
for (id in unique(md$depmap_id)) {
  path <- panel$input_paths[match(id, panel$cell_map$depmap_id)]
  molecular <- qread(path)
  null <- as.matrix(molecular$null[[id]][, bank$genes, drop = FALSE] %*% bank$weights)
  null_tables[[id]] <- data.table(depmap_id = id, program = colnames(null),
    null_sd = apply(null, 2L, sd), null_abs_95 = apply(abs(null), 2L, quantile, probs = .95))
  ids <- md[depmap_id == id & drug %in% focus_drugs, observation_id]
  if (length(ids)) {
    x <- molecular$response[ids, bank$genes, drop = FALSE]
    stopifnot(max(abs(as.matrix(x %*% bank$weights) - bank$scores[ids, , drop = FALSE])) < 1e-10)
    focus_gene_rows[[id]] <- x
  }
}
null_summary <- rbindlist(null_tables)
save_table(null_summary, 'program_vehicle_variation')
focus_genes <- do.call(rbind, focus_gene_rows)
qsave(focus_genes, file.path(root, 'data/focus_drug_gene_displacements.qs'), preset = 'high')

# Retrospective low/high response quartiles at each dose, with replicate range shown.
# Groups are descriptive and use existing collapsed outcomes, so their replicate
# consistency is technical support, not independent confirmation of selection.
groups <- copy(md)
groups[, response_group := fcase(observed_log2 <= quantile(observed_log2, .25), 'Lower abundance',
  observed_log2 >= quantile(observed_log2, .75), 'Higher abundance', default = 'Middle'), by = .(drug, dose)]
score_long <- melt(as.data.table(bank$scores[md$observation_id, , drop = FALSE], keep.rownames = 'observation_id'),
  id.vars = 'observation_id', variable.name = 'program', value.name = 'displacement', variable.factor = FALSE)
score_long <- merge(score_long, groups, by = 'observation_id')
score_long <- merge(score_long, null_summary, by = c('depmap_id', 'program'))
score_long[, above_control := abs(displacement) > null_abs_95]
save_table(score_long, 'all_drug_program_observations')
contrasts <- score_long[response_group != 'Middle', .(
  lower_cells = sum(response_group == 'Lower abundance'), higher_cells = sum(response_group == 'Higher abundance'),
  lower_abundance_pct = median(measured_pct[response_group == 'Lower abundance']),
  higher_abundance_pct = median(measured_pct[response_group == 'Higher abundance']),
  lower_program = mean(displacement[response_group == 'Lower abundance']),
  higher_program = mean(displacement[response_group == 'Higher abundance']),
  control_sd = sqrt(mean(null_sd^2)),
  lower_above_control = mean(above_control[response_group == 'Lower abundance']),
  higher_above_control = mean(above_control[response_group == 'Higher abundance'])), by = .(drug, dose, program)]
contrasts[, `:=`(abundance_gap_pp = higher_abundance_pct - lower_abundance_pct,
  program_difference = lower_program - higher_program)]
save_table(contrasts, 'drug_dose_program_contrasts')

# All-gene signed contrasts identify the genes that actually move inside a named program.
gene_contrasts <- list()
for (drug_name in intersect(focus_drugs, md$drug)) {
  for (d in sort(unique(groups[drug == drug_name, dose]))) {
    g <- groups[drug == drug_name & dose == d]
    lo <- g[response_group == 'Lower abundance', observation_id]
    hi <- g[response_group == 'Higher abundance', observation_id]
    if (!length(lo) || !length(hi)) next
    delta <- colMeans(focus_genes[lo, , drop = FALSE]) - colMeans(focus_genes[hi, , drop = FALSE])
    gene_contrasts[[paste(drug_name, d)]] <- data.table(drug = drug_name, dose = d, gene = names(delta),
      program = colnames(bank$weights)[bank$partition], lower_minus_higher = as.numeric(delta))
  }
}
gene_contrasts <- rbindlist(gene_contrasts)
save_table(gene_contrasts, 'focus_drug_gene_contrasts')
gene_coherence <- gene_contrasts[, .(genes = .N, mean_gene_difference = mean(lower_minus_higher),
  median_gene_difference = median(lower_minus_higher),
  fraction_same_direction = mean(sign(lower_minus_higher) == sign(mean(lower_minus_higher))),
  changing_genes = paste(gene[order(-abs(lower_minus_higher))][seq_len(min(12L, .N))], collapse = ', ')),
  by = .(drug, dose, program)]
save_table(gene_coherence, 'focus_program_gene_coherence')
selected_contrasts <- merge(contrasts, unique(selected[, .(drug, program, view)]), by = c('drug', 'program'))
selected_contrasts <- merge(selected_contrasts, gene_coherence, by = c('drug', 'dose', 'program'), all.x = TRUE)
save_table(selected_contrasts, 'selected_program_biological_contrasts')

# Published ERK output genes and outcome-independent Hallmark sets provide
# mechanistic reference readouts. They are diagnostics, not fitted replacements.
# Compare the same genes with their reconstruction from the hard program means:
# each gene is replaced by the mean of its assigned43 program. No gene re-selection.
anchors <- list(ERK_output = c('PHLDA1', 'SPRY2', 'DUSP6', 'DUSP4', 'ETV4', 'ETV5'))
for (term_name in c('HALLMARK_E2F_TARGETS', 'HALLMARK_MYC_TARGETS_V1',
  'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_APOPTOSIS')) {
  anchors[[term_name]] <- hallmark[term == term_name, gene]
}
anchor_rows <- anchor_members <- list()
for (name in names(anchors)) {
  original <- anchors[[name]]
  genes <- intersect(original, bank$genes)
  stopifnot(length(genes) >= 3L)
  program <- colnames(bank$weights)[bank$partition[match(genes, bank$genes)]]
  anchor_members[[name]] <- data.table(readout = name, gene = original,
    available = original %in% genes, program = program[match(original, genes)])
  measured <- rowMeans(focus_genes[, genes, drop = FALSE])
  reconstructed <- rowMeans(bank$scores[rownames(focus_genes), program, drop = FALSE])
  anchor_rows[[name]] <- data.table(observation_id = names(measured), readout = name,
    measured_gene_score = as.numeric(measured), program_reconstruction = as.numeric(reconstructed))
}
anchor_data <- merge(rbindlist(anchor_rows), groups, by = 'observation_id')
save_table(rbindlist(anchor_members), 'mechanism_reference_membership')
save_table(anchor_data, 'mechanism_reference_observations')
anchor_associations <- anchor_data[, .(cells = .N,
  gene_readout_r = cor_safe(measured_gene_score, observed_log2),
  program_reconstruction_r = cor_safe(program_reconstruction, observed_log2),
  gene_program_agreement = cor_safe(measured_gene_score, program_reconstruction),
  lower_gene_readout = mean(measured_gene_score[response_group == 'Lower abundance']),
  higher_gene_readout = mean(measured_gene_score[response_group == 'Higher abundance'])),
  by = .(drug, dose, readout)]
save_table(anchor_associations, 'mechanism_reference_associations')

# Biological hypothesis prompted by ERK suppression in both response groups:
# a coexisting PI3K driver may uncouple ERK suppression from growth response.
# This is an exploratory genotype comparison, not an independent validation.
# Hotspot absence does not establish that the full pathway is wild type.
mutation_header <- names(fread(paths[8], nrows = 0L))
mutation_genes <- sub(' \\(.*', '', mutation_header)
mutation_columns <- which(mutation_genes %in% c('BRAF', 'KRAS', 'NRAS', 'PIK3CA'))
mutation <- fread(paths[8], select = c(1L, mutation_columns))
setnames(mutation, c('depmap_id', mutation_genes[mutation_columns]))
stopifnot(!anyDuplicated(mutation$depmap_id), all(md$depmap_id %in% mutation$depmap_id))
genotype <- merge(anchor_data[drug %in% c('trametinib', 'cobimetinib') & readout == 'ERK_output'],
  mutation, by = 'depmap_id')
genotype[, RAS_RAF_hotspot := BRAF == 1L | KRAS == 1L | NRAS == 1L]
genotype_predictions <- rows[drug %in% c('trametinib', 'cobimetinib') &
  evaluation %in% c('Whole-cell validation', 'Held cells'),
  .(observation_id, evaluation, predicted_pct)]
stopifnot(!anyDuplicated(genotype_predictions$observation_id))
genotype <- merge(genotype, genotype_predictions, by = 'observation_id')
save_table(genotype, 'mek_genotype_biology')
genotype_summary <- genotype[RAS_RAF_hotspot == TRUE, .(cells = .N,
  measured_abundance = median(measured_pct), predicted_abundance = median(predicted_pct),
  ERK_gene_readout = median(measured_gene_score),
  minimum_replicates = min(replicates, na.rm = TRUE),
  recipients = paste(recipient, collapse = ', ')), by = .(drug, dose, PIK3CA)]
save_table(genotype_summary, 'mek_genotype_summary')
g <- copy(genotype[RAS_RAF_hotspot == TRUE & dose == .5])
g[, group := factor(PIK3CA, levels = 0:1, labels = c('No PIK3CA hotspot\n(n=20)', 'PIK3CA hotspot\n(n=5)'))]
genotype_plots <- list()
for (view in c('measured_pct', 'measured_gene_score', 'predicted_pct')) {
  z <- copy(g)
  z[, value := get(view)]
  label <- switch(view, measured_pct = 'Measured abundance (% vehicle)',
    measured_gene_score = 'ERK gene displacement', predicted_pct = 'Shared-model abundance (% vehicle)')
  genotype_plots[[view]] <- ggplot(z, aes(group, value, color = group)) +
    geom_boxplot(outlier.shape = NA, width = .55, color = 'grey55') +
    geom_point(aes(shape = evaluation), position = position_jitter(width = .09, height = 0, seed = 3), size = 1.8) +
    facet_wrap(~drug, ncol = 1, scales = 'free_y') +
    scale_color_manual(values = c('#247F8A', '#B95720')) +
    labs(x = NULL, y = label) + theme_classic(base_size = 10) +
    theme(legend.position = 'none', axis.text.x = element_text(size = 8))
}
genotype_figure <- wrap_plots(genotype_plots, nrow = 1) +
  plot_annotation(title = 'MEK inhibition: ERK response, genotype and survival at0.5uM',
    subtitle = 'RAS/RAF-hotspot cells only. Same cells in all panels; each point is one cell.',
    caption = 'Exploratory association: five co-mutant cells, including a sensitive exception (LS180).\nPredictions exclude the cell from decoder fitting; bank43 discovery used full-panel outcomes. Hotspot absence is not full wild-type status.')
ggsave(file.path(root, 'figures/MEK_genotype_response.pdf'), genotype_figure, width = 11, height = 7)
ggsave(file.path(root, 'figures/MEK_genotype_response.png'), genotype_figure, width = 11, height = 7, dpi = 170)

# Influence check: remove each cell's whole dose curve, without choosing a new
# program, and re-evaluate the development association. A single-cell-driven
# correlation is not evidence of a general biological response relationship.
influence <- list()
for (i in seq_len(nrow(selected))) {
  s <- selected[i]
  d <- rows[drug == s$drug & evaluation == 'Own fitting data']
  d$score <- bank$scores[d$observation_id, s$program]
  center <- s$view == 'Differences between cells at fixed dose'
  deleted <- rbindlist(lapply(unique(d$depmap_id), function(id) {
    z <- copy(d[depmap_id != id])
    if (center) z[, `:=`(score = score - mean(score), observed_log2 = observed_log2 - mean(observed_log2)), by = dose]
    data.table(omitted_cell = id, r = cor_safe(z$score, z$observed_log2))
  }))
  deleted[, `:=`(drug = s$drug, program = s$program, view = s$view)]
  influence[[i]] <- deleted
}
influence <- rbindlist(influence)
save_table(influence, 'whole_cell_deletion_influence')
save_table(influence[, .(minimum_r = min(r, na.rm = TRUE), maximum_r = max(r, na.rm = TRUE),
  weakest_abs_r = min(abs(r), na.rm = TRUE), influential_cell = omitted_cell[which.min(abs(r))]),
  by = .(drug, program, view)], 'program_association_stability')

# Saved decoder: evaluate the SAME phenotype quartiles, including technical support.
model_rows <- merge(rows, groups[, .(observation_id, response_group, replicates,
  all_replicates_below_50, all_replicates_above_50)], by = 'observation_id')
capture <- model_rows[response_group != 'Middle', .(lower_cells = sum(response_group == 'Lower abundance'),
  higher_cells = sum(response_group == 'Higher abundance'),
  measured_low_pct = median(measured_pct[response_group == 'Lower abundance']),
  measured_high_pct = median(measured_pct[response_group == 'Higher abundance']),
  predicted_low_pct = median(predicted_pct[response_group == 'Lower abundance']),
  predicted_high_pct = median(predicted_pct[response_group == 'Higher abundance']),
  raw_curve_r = cor_safe(observed_log2, predicted_log2)), by = .(family, drug, dose, evaluation)]
save_table(capture, 'shared_model_same_biological_contrasts')
cv_capture <- copy(rows[evaluation == 'Whole-cell validation'])
cv_capture[, `:=`(yc = observed_log2 - mean(observed_log2),
  local_c = local_prediction - mean(local_prediction), shared_c = predicted_log2 - mean(predicted_log2)),
  by = .(drug, dose, fold)]
save_table(cv_capture[, .(cells = .N, local_program_r = cor_safe(yc, local_c),
  shared_model_r = cor_safe(yc, shared_c)), by = .(drug, dose)], 'same_dose_excluded_cell_capture')
message('Drug contrasts, genes, control variation and saved-model capture saved')

# Named biological readouts, with the actual genes and phenotype visible together.
# Marker lists are fixed by mechanism; they are not selected by outcome correlation.
source('/work/src/R/lib_heatmap.R')
marker_sets <- list(
  ERK = c('PHLDA1', 'SPRY2', 'DUSP6', 'DUSP4', 'ETV4', 'ETV5'),
  Cell_cycle = c('MCM2', 'MCM5', 'PCNA', 'TOP2A', 'MKI67', 'CCNA2', 'CCNB1', 'CCNE2', 'RRM2', 'CDC6', 'TYMS', 'CDKN1A'),
  ER_stress = c('ATF4', 'DDIT3', 'TRIB3', 'ASNS', 'CHAC1', 'PMAIP1', 'HSPA5', 'XBP1', 'HERPUD1', 'DNAJB9'))
marker_drugs <- data.table(drug = c('trametinib', 'cobimetinib', 'afatinib', 'palbociclib', 'ribociclib', 'bortezomib', 'ixazomib'),
  marker = c('ERK', 'ERK', 'ERK', 'Cell_cycle', 'Cell_cycle', 'ER_stress', 'ER_stress'))
pdf(file.path(root, 'figures/named_biological_readouts.pdf'), width = 16, height = 7.5)
for (i in seq_len(nrow(marker_drugs))) {
  drug_name <- marker_drugs$drug[i]
  if (!drug_name %in% md$drug) next
  z <- md[drug == drug_name][order(dose, measured_pct)]
  genes <- intersect(marker_sets[[marker_drugs$marker[i]]], bank$genes)
  missing <- setdiff(marker_sets[[marker_drugs$marker[i]]], genes)
  matrix <- t(focus_genes[z$observation_id, genes, drop = FALSE])
  cm <- data.frame(Dose = factor(z$dose), row.names = z$observation_id)
  rm <- data.frame(Program = factor(colnames(bank$weights)[bank$partition[match(genes, bank$genes)]]), row.names = genes)
  stopifnot(identical(colnames(matrix), rownames(cm)), identical(rownames(matrix), rownames(rm)))
  limit <- max(.1, as.numeric(quantile(abs(matrix), .98)))
  hm <- draw_heatmap(matrix, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    col_anno_df = cm, col_anno_cols = 'Dose', row_anno_df = rm, row_anno_cols = 'Program',
    ann_top_args = list(`PRISM % vehicle` = ComplexHeatmap::anno_barplot(z$measured_pct,
      ylim = c(0, max(125, z$measured_pct)), height = grid::unit(14, 'mm'), gp = grid::gpar(fill = '#247F8A', col = NA))),
    column_split = cm$Dose, cluster_column_slices = FALSE,
    show_row_names = TRUE, show_column_names = TRUE, column_labels = z$recipient,
    column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 5), row_names_gp = grid::gpar(fontsize = 9),
    value_range = c(-limit, 0, limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'Gene log2 delta', raster_resize_mat = FALSE,
    title = paste0(drug_name, ' | named ', marker_drugs$marker[i], ' readouts\n',
      'All cells; doses retained; columns ordered by PRISM abundance within dose. Original gene displacements, no row scaling.',
      if (length(missing)) paste0('\nAbsent from bank universe: ', paste(missing, collapse = ', ')) else ''))
  ComplexHeatmap::draw(hm)
  png(file.path(root, 'figures', paste0(drug_name, '_named_readouts.png')), width = 1920, height = 900, res = 120)
  ComplexHeatmap::draw(hm)
  dev.off()
}
dev.off()

# Compact drug dossiers use the actual selected program and show every matched cell.
pdf(file.path(root, 'figures/drug_biology_dossiers.pdf'), width = 11, height = 8)
for (drug_name in intersect(focus_drugs, md$drug)) {
  programs <- unique(selected[drug == drug_name, program])
  for (program_name in programs) {
    z <- score_long[drug == drug_name & program == program_name]
    z[, cohort := fifelse(depmap_id %in% unique(rows[evaluation == 'Held cells', depmap_id]), 'Held cells', 'Development cells')]
    ann <- z[, .(r = cor_safe(displacement, observed_log2)), by = dose]
    p <- ggplot(z, aes(displacement, measured_pct, color = cohort, shape = cohort)) +
      geom_linerange(aes(ymin = replicate_min_pct, ymax = replicate_max_pct), alpha = .45, linewidth = .4, na.rm = TRUE) +
      geom_point(size = 1.7) + geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, linewidth = .4) +
      facet_wrap(~dose, scales = 'free_x', nrow = 1) +
      geom_text(data = ann, aes(x = -Inf, y = Inf, label = sprintf('r(log2 abundance) = %.2f', r)),
        inherit.aes = FALSE, hjust = -.05, vjust = 1.2, size = 3) +
      scale_color_manual(values = c('Development cells' = '#4C626D', 'Held cells' = '#B95720')) +
      labs(x = paste(program_name, 'mean gene displacement'), y = 'PRISM abundance (% vehicle)',
        title = paste(drug_name, '|', program_name),
        subtitle = 'Columns: dose (uM). Vertical ranges: matched PRISM technical replicates. All matched cells shown.') +
      theme_classic(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom')
    cs <- selected_contrasts[drug == drug_name & program == program_name]
    cs <- unique(cs, by = c('drug', 'dose', 'program'))
    theme_rows <- enrichment[program == program_name & overlap >= 3 & q_value <= .05][order(q_value, -fold_enrichment)][1:3]
    theme_rows <- theme_rows[!is.na(pathway)]
    theme_text <- if (nrow(theme_rows)) paste0(sub('^[^|]+\\|', '', theme_rows$pathway),
      ' (', theme_rows$overlap, ' genes; q=', format(theme_rows$q_value, digits = 2), ')') else
      'No pathway enrichment supported at BH q <=0.05; do not assign a mechanism from nominal top terms.'
    biology <- c(paste0('Program: ', cs$genes[1], ' genes. Descriptive pathway enrichments:'),
      theme_text, '',
      paste0(cs$dose, ' uM: lower/higher abundance ', round(cs$lower_abundance_pct), '/', round(cs$higher_abundance_pct),
        '%; program difference ', round(cs$program_difference, 3),
        '; genes with same direction ', round(100 * cs$fraction_same_direction), '%'),
      '', 'Largest gene changes in these contrasts:',
      paste0(cs$dose, ' uM: ', cs$changing_genes), '',
      'Quartiles and bank discovery use these outcomes: descriptive biology, not independent validation.',
      'Vehicle variation is a control reference; PRISM repeats are technical, not biological replication.')
    annotation <- ggplot() + annotate('text', x = 0, y = 1,
      label = paste(biology, collapse = '\n'), hjust = 0, vjust = 1, size = 2.7) +
      xlim(0, 1) + ylim(0, 1) + theme_void()
    figure <- p / annotation + plot_layout(heights = c(1.2, 1))
    print(figure)
    ggsave(file.path(root, 'figures', paste0(drug_name, '_', program_name, '.png')), figure, width = 11, height = 8, dpi = 160)
  }
}
dev.off()
stopifnot(identical(hashes, vapply(paths, function(p) digest(file = p, algo = 'sha256'), character(1))))
save_table(data.table(script = '/work/src/program_pilot/56_drug_program_biology_reliability.R',
  sha256 = digest(file = '/work/src/program_pilot/56_drug_program_biology_reliability.R', algo = 'sha256')),
  'script_provenance')
writeLines(c('56: Drug biology and technical measurement support in bank43.',
  paste('Drugs:', uniqueN(md$drug), '| Cells:', uniqueN(md$depmap_id), '| Observations:', nrow(md)),
  paste('Observations with >=2 paired PRISM replicates:', sum(md$replicates >= 2, na.rm = TRUE)),
  'See data tables for all drugs and biological dossier figures for the declared mechanism groups.',
  'No new decoder fitted. Program selection uses development cells; biology/replicate checks remain retrospective.',
  'Source files hash-identical after execution.'), file.path(root, 'REPORT.txt'))
message('Completed ', root)
