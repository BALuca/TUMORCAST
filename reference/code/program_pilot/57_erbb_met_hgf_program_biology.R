#!/usr/bin/env Rscript
# Biological diagnosis of the current bank along ERBB, HGF-MET and downstream responses.
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
root <- '/work/results/program_pilot/57_erbb_met_hgf_program_biology'
prior <- '/work/results/program_pilot/56_drug_program_biology_reliability/data'
panel_root <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data'
prism_root <- '/work/cache/tahoe_prism_full_chain/prism'
bank_path <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs'
cyto_path <- '/work/data/01_CytoSig_data/diff.merge.gz'
drug_roles <- data.table(
  drug = c('lapatinib', 'afatinib', 'neratinib', 'erlotinib', 'crizotinib (hydrochloride)',
    'crizotinib', 'capmatinib', 'alpelisib', 'trametinib', 'cobimetinib'),
  role = c(rep('ERBB inhibition', 4), rep('MET inhibition', 3), 'PI3K inhibition', rep('MEK inhibition', 2)))
for (d in file.path(root, c('data', 'figures'))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
cor_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L || sd(x[ok]) < 1e-12 || sd(y[ok]) < 1e-12) return(NA_real_)
  cor(x[ok], y[ok])
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 160)
}
paths <- c(bank_path, file.path(panel_root, 'molecular_inputs/panel_index.qs'), cyto_path,
  '/work/data/01_CytoSig_data/meta_info.txt', file.path(prior, 'program_pathway_enrichment.tsv.gz'),
  file.path(prior, 'observation_reliability.tsv.gz'), file.path(prior, 'matched_prism_replicates.tsv.gz'),
  '/work/grant_rebuild_2026-10-04/solution_evidence/hallmark_membership.tsv')
stopifnot(all(file.exists(paths)))
save_table(data.table(path = paths, sha256 = vapply(paths, function(p) digest(file = p, algo = 'sha256'), character(1))),
  'input_manifest')
bank <- qread(bank_path)
panel <- qread(paths[2])
md <- merge(panel$metadata[drug %in% drug_roles$drug], drug_roles, by = 'drug')
stopifnot(setequal(md$drug, drug_roles$drug), !anyDuplicated(md$observation_id))
sens <- copy(panel$sensitivity[observation_id %in% md$observation_id])
names_to_keep <- c('observation_id', 'prism_abundance_pct', 'prism_lower_dose', 'prism_upper_dose', 'prism_screen')
md <- merge(md, sens[, ..names_to_keep], by = 'observation_id')
md[, measured_pct := prism_abundance_pct]
md[, observed_log2 := log2(measured_pct / 100)]

# Extend the existing PRISM salt identity explicitly; no other source identity changes.
pi <- fread(file.path(prism_root, 'secondary-screen-replicate-collapsed-treatment-info.csv'))[tolower(name) == 'crizotinib']
ph <- names(fread(file.path(prism_root, 'secondary-screen-replicate-collapsed-logfold-change.csv'), nrows = 0L))
px <- fread(file.path(prism_root, 'secondary-screen-replicate-collapsed-logfold-change.csv'), select = c(ph[1], pi$column_name))
setnames(px, 1L, 'depmap_id')
px <- px[depmap_id %in% md$depmap_id]
points <- melt(px, id.vars = 'depmap_id', variable.name = 'column_name', value.name = 'value', variable.factor = FALSE)
points <- merge(points, pi[, .(column_name, dose, screen_id, broad_id)], by = 'column_name')
points <- points[is.finite(value), .(value = mean(value)), by = .(depmap_id, dose, screen_id)]
for (i in which(md$drug == 'crizotinib (hydrochloride)')) {
  id <- md$depmap_id[i]
  query_dose <- md$dose[i]
  p <- points[depmap_id == id]
  screens <- c(intersect('MTS010', unique(p$screen_id)), setdiff(sort(unique(p$screen_id)), 'MTS010'))
  for (screen in screens) {
    z <- p[screen_id == screen][order(dose)]
    if (query_dose < min(z$dose) || query_dose > max(z$dose)) next
    exact <- which(abs(log(z$dose / query_dose)) < 1e-6)
    lo <- if (length(exact)) exact[1] else max(which(z$dose < query_dose))
    hi <- if (length(exact)) exact[1] else min(which(z$dose > query_dose))
    y <- if (lo == hi) z$value[lo] else approx(log(z$dose[c(lo, hi)]), z$value[c(lo, hi)], xout = log(query_dose))$y
    md[i, `:=`(measured_pct = 100 * 2^y, observed_log2 = y, prism_lower_dose = z$dose[lo],
      prism_upper_dose = z$dose[hi], prism_screen = screen)]
    break
  }
}
save_table(data.table(tahoe_drug = 'crizotinib (hydrochloride)', prism_drug = 'crizotinib',
  rationale = 'Same active crizotinib moiety; distinct molecular preparation retained; (S)-crizotinib excluded'), 'explicit_salt_alias')
md[, cohort := 'Development cells']
old_pred <- fread('/work/results/program_pilot/55_current_bank_within_treatment_signal_audit/data/diagnostic_predictions.tsv.gz')
held_ids <- unique(old_pred[evaluation == 'Held cells', depmap_id])
md[depmap_id %in% held_ids, cohort := 'Originally held cells']
md[!is.finite(measured_pct), cohort := 'No PRISM outcome']

# Reuse56 technical repeats and extend exactly the same interpolation to the salt panel.
replicates <- fread(paths[7])[observation_id %in% md$observation_id]
ri <- fread(file.path(prism_root, 'secondary-screen-replicate-treatment-info.csv'))[tolower(name) == 'crizotinib']
rh <- names(fread(file.path(prism_root, 'secondary-screen-logfold-change.csv'), nrows = 0L))
ri <- ri[column_name %in% rh]
rx <- fread(file.path(prism_root, 'secondary-screen-logfold-change.csv'), select = c(rh[1], ri$column_name))
setnames(rx, 1L, 'depmap_id')
rx <- rx[depmap_id %in% md$depmap_id]
rcols <- setdiff(names(rx), 'depmap_id')
rx[, (rcols) := lapply(.SD, as.numeric), .SDcols = rcols]
rl <- melt(rx, id.vars = 'depmap_id', variable.name = 'column_name', value.name = 'value', variable.factor = FALSE)
rl <- merge(rl[is.finite(value)], ri[, .(column_name, dose, screen_id, broad_id, detection_plate)], by = 'column_name')
extra_reps <- list()
for (i in which(md$drug == 'crizotinib (hydrochloride)' & is.finite(md$observed_log2))) {
  z <- md[i]
  a <- rl[depmap_id == z$depmap_id & screen_id == z$prism_screen & abs(log(dose / z$prism_lower_dose)) < 1e-7]
  b <- rl[depmap_id == z$depmap_id & screen_id == z$prism_screen & abs(log(dose / z$prism_upper_dose)) < 1e-7]
  ab <- merge(a, b, by = c('detection_plate', 'broad_id'))
  if (!nrow(ab)) next
  fraction <- if (z$prism_lower_dose == z$prism_upper_dose) 0 else
    log(z$dose / z$prism_lower_dose) / log(z$prism_upper_dose / z$prism_lower_dose)
  extra_reps[[length(extra_reps) + 1L]] <- ab[, .(observation_id = z$observation_id, drug = z$drug,
    depmap_id = z$depmap_id, recipient = z$recipient, dose = z$dose, detection_plate, broad_id,
    replicate = sub('.*_X([0-9]+).*', '\\1', detection_plate), value = value.x + fraction * (value.y - value.x))]
}
replicates <- rbind(replicates, rbindlist(extra_reps), fill = TRUE)
replicates[, abundance_pct := 100 * 2^value]
save_table(replicates, 'matched_prism_replicates')
rep_summary <- replicates[, .(repeats = .N, repeat_low = min(abundance_pct), repeat_high = max(abundance_pct),
  all_repeats_sensitive = all(abundance_pct <= 50)), by = observation_id]
md <- merge(md, rep_summary, by = 'observation_id', all.x = TRUE)
rc <- replicates[, .(value = mean(value)), by = .(drug, dose, depmap_id, replicate)]
rp <- merge(rc, rc, by = c('drug', 'dose', 'depmap_id'), allow.cartesian = TRUE)[replicate.x < replicate.y]
save_table(rp[, .(cells = .N, r = cor_safe(value.x, value.y)), by = .(drug, dose, replicate.x, replicate.y)],
  'within_dose_replicate_agreement')

# Project the unchanged hard program means, including previously unscored salts.
gene_rows <- baseline_rows <- list()
for (id in unique(md$depmap_id)) {
  obj <- qread(panel$input_paths[match(id, panel$cell_map$depmap_id)])
  ids <- md[depmap_id == id, observation_id]
  stopifnot(all(ids %in% rownames(obj$response)), identical(colnames(obj$response), obj$genes))
  gene_rows[[id]] <- obj$response[ids, , drop = FALSE]
  baseline_rows[[id]] <- data.table(depmap_id = id, gene = obj$genes, expression = obj$baseline[[id]])
}
gene_delta <- do.call(rbind, gene_rows)
scores <- as.matrix(gene_delta[, bank$genes, drop = FALSE] %*% bank$weights)
common <- intersect(rownames(scores), rownames(bank$scores))
agreement <- max(abs(scores[common, , drop = FALSE] - bank$scores[common, , drop = FALSE]))
stopifnot(agreement < 1e-10)
baseline <- rbindlist(baseline_rows)
receptors <- dcast(baseline[gene %in% c('ERBB2', 'EGFR', 'MET', 'ALK', 'ROS1')], depmap_id ~ gene, value.var = 'expression')
stopifnot(all(c('ERBB2', 'EGFR', 'MET') %in% names(receptors)))
md <- merge(md, receptors, by = 'depmap_id')
setorder(md, drug, dose, measured_pct, depmap_id)
save_table(md, 'drug_observations')
save_table(md[, .(profiles = .N, cells = uniqueN(depmap_id), with_PRISM = sum(is.finite(measured_pct)),
  sensitive = sum(measured_pct <= 50, na.rm = TRUE), repeat_supported_sensitive = sum(all_repeats_sensitive & repeats >= 2, na.rm = TRUE)),
  by = .(drug, dose)], 'drug_coverage_and_sensitive_counts')
score_long <- melt(as.data.table(scores, keep.rownames = 'observation_id'), id.vars = 'observation_id',
  variable.name = 'program', value.name = 'displacement', variable.factor = FALSE)
score_long <- merge(score_long, md, by = 'observation_id')
save_table(score_long, 'drug_program_observations')
message('Drug inputs ready: ', nrow(md), ' profiles; exact frozen score agreement = ', agreement)

# Crizotinib free base and salt are paired molecular preparations, not biological replicates.
salt_pairs <- merge(md[drug == 'crizotinib', .(depmap_id, dose, free_id = observation_id)],
  md[drug == 'crizotinib (hydrochloride)', .(depmap_id, dose, salt_id = observation_id)], by = c('depmap_id', 'dose'))
salt_agreement <- rbindlist(lapply(seq_len(nrow(salt_pairs)), function(i) {
  z <- salt_pairs[i]
  data.table(z, gene_r = cor_safe(gene_delta[z$free_id, bank$genes], gene_delta[z$salt_id, bank$genes]),
    program_r = cor_safe(scores[z$free_id, ], scores[z$salt_id, ]))
}))
save_table(salt_agreement, 'crizotinib_preparation_agreement')

# Existing CytoSig contrasts; retain every selected experiment and missing gene explicitly.
cm <- fread(paths[4], skip = 1L, header = FALSE, col.names = c('profile_id', 'treatment', 'condition', 'duration', 'dose'))
cm <- cm[toupper(treatment) %in% c('HGF', 'EGF', 'NRG1')]
cm[, study := sub('.*@((GSE[0-9]+|E-MEXP-[0-9]+|E-MTAB-[0-9]+)).*', '\\1', profile_id)]
stopifnot(all(grepl('^(GSE|E-MEXP|E-MTAB)', cm$study)), !anyDuplicated(cm$profile_id))
cm[, study_family := fifelse(study %in% c('GSE26732', 'GSE26735'), 'GSE26736', study)]
cm[, duration_original := duration]
cm[, duration_note := 'As recorded by CytoSig']
# The GEO/CytoSig48h metadata conflicts with the original RNA experiment methods.
# All three papers using these same RNA samples specify16h; keep the original too.
cm[study == 'GSE135552', `:=`(duration = '16h',
  duration_note = 'Corrected from48h to16h using RNA methods in PMID31940827 and PMC7699574')]
cm[, hours := suppressWarnings(as.numeric(sub('[^0-9.].*$', '', duration)))]
cm[grepl('m$', duration), hours := hours / 60]
cm[grepl('d$', duration), hours := hours * 24]
cm[, time_group := fcase(hours <= 6, 'Early <=6h', hours > 6, 'Later >6h', default = 'Time unspecified')]
cyhead <- names(fread(cyto_path, nrows = 0L))
stopifnot(all(cm$profile_id %in% cyhead))
cy <- fread(cyto_path, select = c(cyhead[1], cm$profile_id), showProgress = FALSE)
setnames(cy, 1L, 'gene')
stopifnot(!anyDuplicated(cy$gene))
cy_matrix <- t(as.matrix(cy[, -1L]))
colnames(cy_matrix) <- cy$gene
cy_scores <- matrix(NA_real_, nrow(cm), ncol(scores), dimnames = list(cm$profile_id, colnames(scores)))
coverage <- list()
for (p in colnames(scores)) {
  genes <- bank$genes[as.numeric(bank$weights[, p]) > 0]
  present <- intersect(genes, colnames(cy_matrix))
  z <- cy_matrix[cm$profile_id, present, drop = FALSE]
  n <- rowSums(is.finite(z))
  cy_scores[, p] <- rowMeans(z, na.rm = TRUE)
  cy_scores[n < .8 * length(genes), p] <- NA_real_
  coverage[[p]] <- data.table(profile_id = cm$profile_id, program = p, observed_genes = n,
    program_genes = length(genes), fraction = n / length(genes))
}
save_table(rbindlist(coverage), 'ligand_program_gene_coverage')
cy_long <- melt(as.data.table(cy_scores, keep.rownames = 'profile_id'), id.vars = 'profile_id',
  variable.name = 'program', value.name = 'displacement', variable.factor = FALSE)
cy_long <- merge(cy_long, cm, by = 'profile_id')
save_table(cm, 'ligand_experiments')
save_table(cy_long, 'ligand_program_observations')
cy_bg_genes <- intersect(bank$genes, colnames(cy_matrix))
cy_background <- data.table(profile_id = cm$profile_id,
  gene_mean = rowMeans(cy_matrix[cm$profile_id, cy_bg_genes, drop = FALSE], na.rm = TRUE),
  gene_median = apply(cy_matrix[cm$profile_id, cy_bg_genes, drop = FALSE], 1L, median, na.rm = TRUE))
save_table(merge(cy_background, cm, by = 'profile_id'), 'ligand_genome_background')
# Average experimental variants within each study first; studies then have equal weight.
cy_study <- cy_long[, .(displacement = mean(displacement, na.rm = TRUE)), by = .(treatment, study_family, time_group, program)]
cy_summary <- cy_study[, .(studies = sum(is.finite(displacement)), mean = mean(displacement, na.rm = TRUE),
  positive_studies = sum(displacement > 0, na.rm = TRUE), negative_studies = sum(displacement < 0, na.rm = TRUE)),
  by = .(treatment, time_group, program)]
save_table(cy_summary, 'ligand_study_balanced_program_summary')
message('Ligand inputs ready: ', nrow(cm), ' contrasts, including ', cm[treatment == 'HGF', .N], ' HGF')

# Fixed biological references locate loss at admission, program averaging or decoder.
hm <- fread(paths[8])
readouts <- list(ERK_output = c('PHLDA1', 'SPRY2', 'DUSP6', 'DUSP4', 'ETV4', 'ETV5'),
  Immediate_early = c('FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'EGR2', 'IER2', 'IER3'))
for (term_name in c('HALLMARK_E2F_TARGETS', 'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MTORC1_SIGNALING',
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_APOPTOSIS',
  'HALLMARK_GLYCOLYSIS')) {
  readouts[[sub('HALLMARK_', '', term_name)]] <- hm[term == term_name, gene]
}
readout_members <- drug_readouts <- ligand_readouts <- list()
for (name in names(readouts)) {
  original <- readouts[[name]]
  tg <- intersect(original, colnames(gene_delta))
  bg <- intersect(original, bank$genes)
  programs <- colnames(scores)[bank$partition[match(bg, bank$genes)]]
  readout_members[[name]] <- data.table(readout = name, gene = original, in_Tahoe_input = original %in% tg,
    in_bank = original %in% bg, program = programs[match(original, bg)])
  drug_readouts[[name]] <- data.table(observation_id = rownames(gene_delta), readout = name,
    measured_genes = rowMeans(gene_delta[, tg, drop = FALSE]),
    bank_genes = rowMeans(gene_delta[, bg, drop = FALSE]),
    program_reconstruction = rowMeans(scores[, programs, drop = FALSE]))
  lg <- intersect(original, colnames(cy_matrix))
  bgl <- intersect(bg, colnames(cy_matrix))
  ligand_readouts[[name]] <- data.table(profile_id = cm$profile_id, readout = name,
    measured_genes = rowMeans(cy_matrix[cm$profile_id, lg, drop = FALSE], na.rm = TRUE),
    bank_genes = rowMeans(cy_matrix[cm$profile_id, bgl, drop = FALSE], na.rm = TRUE),
    program_reconstruction = rowMeans(cy_scores[, programs, drop = FALSE], na.rm = TRUE))
}
readout_members <- rbindlist(readout_members)
dr <- merge(rbindlist(drug_readouts), md, by = 'observation_id')
lr <- merge(rbindlist(ligand_readouts), cm, by = 'profile_id')
save_table(readout_members, 'biological_readout_membership')
save_table(dr, 'drug_biological_readouts')
save_table(lr, 'ligand_biological_readouts')
save_table(dr[, .(cells = .N, gene_program_r = cor_safe(bank_genes, program_reconstruction),
  gene_fate_r = cor_safe(measured_genes, observed_log2), program_fate_r = cor_safe(program_reconstruction, observed_log2),
  gene_mean = mean(measured_genes), program_mean = mean(program_reconstruction)), by = .(drug, dose, readout)],
  'readout_recovery_and_function')
receptor_relations <- rbindlist(lapply(c('ERBB2', 'EGFR', 'MET'), function(g) {
  dr[, .(r_baseline_response = cor_safe(get(g), measured_genes), r_baseline_fate = cor_safe(get(g), observed_log2)),
    by = .(drug, dose, readout)][, receptor := g]
}))
save_table(receptor_relations, 'baseline_receptor_associations')
save_table(dr[recipient != 'BT-474' & drug == 'lapatinib' & readout == 'ERK_output',
  .(cells = .N, baseline_response_r = cor_safe(ERBB2, measured_genes),
    baseline_fate_r = cor_safe(ERBB2, observed_log2)), by = dose], 'lapatinib_receptor_without_BT474')

# Drug-level evidence retains dose transitions AND differences among cells at each dose.
associations <- score_long[is.finite(observed_log2), .(cells = .N, r = cor_safe(displacement, observed_log2),
  mean_delta = mean(displacement), mean_sensitive = mean(displacement[measured_pct <= 50]),
  mean_resistant = mean(displacement[measured_pct > 50])), by = .(drug, dose, program)]
save_table(associations, 'program_function_by_dose')
centered <- copy(score_long[is.finite(observed_log2)])
centered[, `:=`(dx = displacement - mean(displacement), dy = observed_log2 - mean(observed_log2)), by = .(drug, dose, program)]
summary <- centered[, .(profiles = .N, cells = uniqueN(depmap_id), full_curve_r = cor_safe(displacement, observed_log2),
  within_dose_r = cor_safe(dx, dy)), by = .(drug, program)]
save_table(summary, 'program_function_full_curve_and_within_dose')
cs <- score_long[, .(mean_delta = mean(displacement)), by = .(drug, dose, program)]
cs <- merge(cs, cy_summary[treatment == 'HGF'], by = 'program', allow.cartesian = TRUE)
cs[, opposing_direction := mean_delta * mean < 0]
save_table(cs, 'hgf_drug_program_direction_comparison')

# Saved shared decoder on exactly these observations; salt extrapolations were not fitted previously.
choice <- unique(old_pred[, .(drug, family)])[order(drug, family != 'MEK')][, .SD[1L], by = drug]
pred <- merge(old_pred, choice, by = c('drug', 'family'))
pred <- pred[observation_id %in% md$observation_id]
save_table(pred, 'saved_decoder_predictions')
save_table(pred[, .(profiles = .N, r = cor_safe(observed_log2, predicted_log2),
  measured_sensitive = sum(measured_pct <= 50), predicted_sensitive = sum(predicted_pct <= 50)),
  by = .(drug, dose, evaluation)], 'saved_decoder_capture')

# Paired biological pictures. Original signed values, no automatic row standardization.
source('/work/src/R/lib_heatmap.R')
heat_pdf_png <- function(h, name, width, height) {
  pdf(file.path(root, 'figures', paste0(name, '.pdf')), width = width, height = height)
  ComplexHeatmap::draw(h, heatmap_legend_side = 'right', annotation_legend_side = 'right',
    padding = grid::unit(c(3, 10, 3, 3), 'mm'))
  dev.off()
  png(file.path(root, 'figures', paste0(name, '.png')), width = width * 150, height = height * 150, res = 150)
  ComplexHeatmap::draw(h, heatmap_legend_side = 'right', annotation_legend_side = 'right',
    padding = grid::unit(c(3, 10, 3, 3), 'mm'))
  dev.off()
}
program_order <- hclust(dist(t(scores)), method = 'ward.D2')$order
enrichment <- fread(paths[5])
theme_top <- enrichment[q_value <= .05 & overlap >= 3][order(q_value)][, .SD[1L], by = program]
program_labels <- setNames(colnames(scores), colnames(scores))
for (p in theme_top$program) {
  label <- sub('^[^|]+\\|', '', theme_top[program == p, pathway])
  label <- gsub('HALLMARK_|_', ' ', label)
  program_labels[p] <- paste(p, substr(label, 1, 37), sep = ' | ')
}
# Short, gene-supported theme labels; complete enrichment statistics remain in56.
short_themes <- c(P03 = 'RNA metabolism', P04 = 'Integrated stress response', P05 = 'Translation',
  P08 = 'Calmodulin-related', P09 = 'MYC targets', P12 = 'Rho / cytoskeleton', P13 = 'MYC / ribosome biogenesis',
  P17 = 'Immediate early / NF-kB', P18 = 'Cell cycle', P22 = 'Cilium trafficking',
  P30 = 'Immediate early / NF-kB', P32 = 'Respiration')
program_labels[names(short_themes)] <- paste(names(short_themes), short_themes, sep = ' | ')
markers <- unique(c(readouts$ERK_output, readouts$Immediate_early,
  'MCM2', 'MCM5', 'PCNA', 'TOP2A', 'MKI67', 'CCNA2', 'CCNB1', 'RRM2', 'CDC6', 'TYMS', 'CDKN1A',
  'ATF4', 'DDIT3', 'ASNS', 'CHAC1', 'TRIB3', 'SLC7A5', 'ITGA5', 'VIM', 'FN1', 'CDH1', 'AXL', 'MET', 'ERBB2',
  'HK2', 'LDHA', 'SLC2A1', 'ADH1B'))
markers <- intersect(markers, colnames(gene_delta))
marker_program <- colnames(scores)[bank$partition[match(markers, bank$genes)]]
marker_program[is.na(marker_program)] <- 'Not in bank43'
for (drug_name in drug_roles$drug) {
  z <- md[drug == drug_name][order(dose, measured_pct, depmap_id)]
  ca <- data.frame(Dose = factor(z$dose), Cohort = z$cohort, row.names = z$observation_id)
  a <- t(scores[z$observation_id, program_order, drop = FALSE])
  stopifnot(identical(colnames(a), rownames(ca)))
  bar <- if (any(is.finite(z$measured_pct))) list(`PRISM % vehicle` = ComplexHeatmap::anno_barplot(
    z$measured_pct, ylim = c(0, max(130, z$measured_pct, na.rm = TRUE)), height = grid::unit(14, 'mm'),
    gp = grid::gpar(fill = '#317E8C', col = NA))) else list()
  h <- draw_heatmap(a, scale_rows = FALSE, col_anno_df = ca, col_anno_cols = c('Dose', 'Cohort'), ann_top_args = bar,
    cluster_rows = FALSE, cluster_columns = FALSE, column_split = ca$Dose, cluster_column_slices = FALSE,
    show_row_names = TRUE, row_labels = program_labels[rownames(a)], row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = TRUE, column_labels = z$recipient, column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = 5.5), value_range = c(-.3, 0, .3),
    value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'Program delta',
    title = paste0(drug_name, ' | frozen bank43\nWithin each dose, cells ordered by PRISM abundance; missing outcomes last. Colors saturate at +/-0.3 log2 delta.'))
  slug <- gsub('[^a-z0-9]+', '_', drug_name)
  heat_pdf_png(h, paste0(slug, '_programs'), 16, 9)
  a <- t(gene_delta[z$observation_id, markers, drop = FALSE])
  ra <- data.frame(Program = marker_program, row.names = markers)
  h <- draw_heatmap(a, scale_rows = FALSE, col_anno_df = ca, col_anno_cols = c('Dose', 'Cohort'), ann_top_args = bar,
    row_anno_df = ra, row_anno_cols = 'Program', cluster_rows = FALSE, cluster_columns = FALSE,
    column_split = ca$Dose, cluster_column_slices = FALSE, show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8), show_column_names = TRUE, column_labels = z$recipient,
    column_names_rot = 90, column_names_gp = grid::gpar(fontsize = 5.5), value_range = c(-1, 0, 1),
    value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'Gene delta',
    title = paste0(drug_name, ' | named downstream response genes\nSame cells and ordering as program heatmap. Original log2 deltas; colors saturate at +/-1.'))
  heat_pdf_png(h, paste0(slug, '_genes'), 16, 10)
  message('Drug figures saved: ', drug_name)
}

# All HGF experiments on a shared order, without hiding cell/time heterogeneity.
hz <- cm[treatment == 'HGF'][order(study_family, condition, hours, study)]
hca <- data.frame(Study = hz$study_family, Time = hz$time_group, row.names = hz$profile_id)
hl <- paste(hz$condition, ifelse(nzchar(hz$duration), hz$duration, '?time'), sep = ' | ')
short_cells <- c('Lymphatic Endothelial' = 'Lymphatic EC', 'Primary Renal Proximal Tubule Epithelial' = 'Renal tubule',
  'Umbilical Vein Endothelial' = 'Umbilical EC', 'Bxpc3 Pancreatic Cancer' = 'BxPC3',
  'Hnscc Cell Line Detroit562' = 'Detroit562', 'Hnscc Cell Line Fadu' = 'FaDu', 'Hnscc Cell Line Scc9' = 'SCC9',
  'Prostate Adenocarcinoma Cell Line Pc3M' = 'PC3M', 'Mda-Mb435' = 'MDA-MB435',
  'Expanded Islets' = 'Islets', 'H1975' = 'H1975')
stopifnot(all(hz$condition %in% names(short_cells)))
hl <- paste(short_cells[hz$condition], ifelse(nzchar(hz$duration), hz$duration, '?time'), sep = ' | ')
hl[hz$study %in% c('GSE26732', 'GSE26735')] <- paste(hl[hz$study %in% c('GSE26732', 'GSE26735')],
  ifelse(hz$study[hz$study %in% c('GSE26732', 'GSE26735')] == 'GSE26732', 'Affy', 'Illumina'))
h <- draw_heatmap(t(cy_scores[hz$profile_id, program_order, drop = FALSE]), scale_rows = FALSE,
  col_anno_df = hca, col_anno_cols = c('Study', 'Time'), column_split = hca$Study,
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE, show_row_names = TRUE,
  row_labels = program_labels[colnames(scores)[program_order]], row_names_gp = grid::gpar(fontsize = 8),
  show_column_names = TRUE, column_labels = hl, column_names_rot = 55, column_names_gp = grid::gpar(fontsize = 7),
  value_range = c(-.3, 0, .3), value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'CytoSig delta',
  title = 'HGF responses through bank43 | every available experiment\nStudy and time retained; grey means <80% gene coverage. GSE26732/35 share one experiment.')
heat_pdf_png(h, 'HGF_all_programs', 13, 10)
mg <- intersect(markers, colnames(cy_matrix))
h <- draw_heatmap(t(cy_matrix[hz$profile_id, mg, drop = FALSE]), scale_rows = FALSE,
  col_anno_df = hca, col_anno_cols = c('Study', 'Time'), column_split = hca$Study,
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE, show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8), show_column_names = TRUE, column_labels = hl,
  column_names_rot = 55, column_names_gp = grid::gpar(fontsize = 7), value_range = c(-1, 0, 1),
  value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'CytoSig gene delta',
  title = 'HGF | named growth-response genes\nOriginal source log2 changes; not quantitatively calibrated to Tahoe displacements.')
heat_pdf_png(h, 'HGF_named_genes', 13, 10)

# Mechanistic summary: doses and cell-specific responses remain points, not only averages.
theme_pub <- theme_classic(base_size = 10) + theme(legend.position = 'bottom')
focal_drugs <- c('lapatinib', 'crizotinib (hydrochloride)', 'capmatinib', 'alpelisib')
readout_view <- c('ERK_output', 'MYC_TARGETS_V1', 'E2F_TARGETS', 'MTORC1_SIGNALING')
dp <- dr[drug %in% focal_drugs & readout %in% readout_view]
dp[, sensitivity := fifelse(is.na(measured_pct), 'Unmeasured', fifelse(measured_pct <= 50, '<=50% vehicle', '>50% vehicle'))]
p <- ggplot(dp, aes(factor(dose), measured_genes, color = sensitivity)) +
  geom_hline(yintercept = 0, color = 'grey65', linewidth = .3) +
  geom_point(position = position_jitter(width = .13, height = 0, seed = 57), size = 1.2, alpha = .75) +
  facet_grid(readout ~ drug, scales = 'free_y') +
  scale_color_manual(values = c('<=50% vehicle' = '#B2462D', '>50% vehicle' = '#237D92', 'Unmeasured' = 'grey55')) +
  labs(x = 'Dose (uM)', y = 'Mean original gene displacement', color = 'PRISM abundance',
    title = 'Do branch inhibitors change the expected downstream readouts?',
    subtitle = 'Each point is a cell. Expression readouts do not directly measure kinase phosphorylation.') + theme_pub
save_plot(p, 'Y_downstream_drug_responses', 12, 9)
lview <- lr[readout %in% readout_view]
p <- ggplot(lview, aes(treatment, measured_genes, color = time_group)) +
  geom_hline(yintercept = 0, color = 'grey65', linewidth = .3) +
  geom_point(position = position_jitter(width = .18, height = 0, seed = 57), size = 1.8, alpha = .8) +
  facet_wrap(~readout, scales = 'free_y', nrow = 1) +
  labs(x = NULL, y = 'Mean source gene log2 change', color = 'Exposure',
    title = 'HGF, EGF and NRG1: shared output or context-dependent response?',
    subtitle = 'Every point is a source contrast, not an independent study; cells, platforms and times differ.') + theme_pub
save_plot(p, 'Y_ligand_downstream_responses', 12, 4)

# Baseline receptor context and gene-readout retention are separate links.
bp <- dr[drug %in% c('lapatinib', 'crizotinib (hydrochloride)') & readout == 'ERK_output']
bp[, receptor_expression := fifelse(drug == 'lapatinib', ERBB2, MET)]
bp[, receptor_label := fifelse(drug == 'lapatinib', 'ERBB2', 'MET')]
bp[, label := fifelse(recipient == 'BT-474' | receptor_expression >= quantile(receptor_expression, .9), recipient, '')]
ann <- bp[, .(r = cor_safe(receptor_expression, measured_genes)), by = .(drug, dose)]
p <- ggplot(bp, aes(receptor_expression, measured_genes)) + geom_point(aes(color = measured_pct), size = 1.8) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#444444', linewidth = .5) +
  geom_text(aes(label = label), check_overlap = TRUE, size = 2.5, vjust = -1) +
  facet_grid(drug ~ dose, scales = 'free') +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = sprintf('r = %.2f', r)),
    inherit.aes = FALSE, hjust = -.05, vjust = 1.3, size = 3) +
  scale_color_gradient(low = '#B2462D', high = '#237D92', na.value = 'grey60') +
  labs(x = 'Baseline ERBB2 (lapatinib) / MET (crizotinib), log2 CP10K', y = 'ERK target-gene displacement',
    color = 'PRISM % vehicle', title = 'Does receptor context explain who responds molecularly?') +
  theme_pub + theme(aspect.ratio = 1)
save_plot(p, 'baseline_receptor_to_response', 10, 7)
ret <- dr[drug %in% c('lapatinib', 'crizotinib (hydrochloride)', 'alpelisib') & readout %in% c('ERK_output', 'E2F_TARGETS')]
ret_ann <- ret[, .(r = cor_safe(bank_genes, program_reconstruction)), by = .(drug, readout)]
p <- ggplot(ret, aes(bank_genes, program_reconstruction, color = factor(dose))) + geom_point(size = 1.4, alpha = .7) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .4) +
  facet_grid(readout ~ drug, scales = 'free') +
  geom_text(data = ret_ann, aes(x = -Inf, y = Inf, label = sprintf('r = %.2f', r)), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.3, size = 3) +
  labs(x = 'Original reference genes retained in bank43', y = 'Same genes replaced by their program means', color = 'Dose (uM)',
    title = 'Which named biological responses survive program averaging?') + theme_pub + theme(aspect.ratio = 1)
save_plot(p, 'biological_signal_retention', 10, 7)

# Direct Y-oriented comparison, all34 rows, with a distinct legend for each source.
# BT474 is shown separately because a panel average can erase a context-specific ERBB response.
overview <- score_long[drug != 'crizotinib', .(displacement = mean(displacement), cells = .N),
  by = .(drug, dose, program)]
overview[, context := paste(drug, dose, 'uM')]
bt <- score_long[recipient == 'BT-474' & drug == 'lapatinib',
  .(drug = 'lapatinib | BT474', dose, program, displacement, cells = 1L,
    context = paste('lapatinib | BT474', dose, 'uM'))]
overview <- rbind(overview, bt)
save_table(overview, 'program_drug_overview')
om <- dcast(overview, program ~ context, value.var = 'displacement')
oa <- as.matrix(om[, -1L])
rownames(oa) <- om$program
oa <- oa[colnames(scores)[program_order], , drop = FALSE]
oc <- unique(overview[, .(context, drug, dose)])
setorder(oc, drug, dose)
oa <- oa[, oc$context, drop = FALSE]
oc <- data.frame(Drug = oc$drug, row.names = oc$context)
left <- draw_heatmap(oa, scale_rows = FALSE, col_anno_df = oc, col_anno_cols = 'Drug',
  column_split = oc$Drug, cluster_column_slices = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, row_labels = program_labels[rownames(oa)], row_names_gp = grid::gpar(fontsize = 7),
  row_names_side = 'left',
  show_column_names = TRUE, column_labels = sub('.* ', '', sub(' uM$', '', colnames(oa))),
  column_names_rot = 0, column_names_gp = grid::gpar(fontsize = 7),
  value_range = c(-.15, 0, .15), value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'Tahoe mean delta',
  title = 'Drugs: mean across39 cells, plus separate BT474 columns\nDoses (uM) retained; colors saturate at +/-0.15')
right_matrix <- t(cy_scores[hz$profile_id, colnames(scores)[program_order], drop = FALSE])
right <- draw_heatmap(right_matrix, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = TRUE, column_labels = hl, column_names_rot = 55,
  column_names_gp = grid::gpar(fontsize = 6.5), value_range = c(-.3, 0, .3),
  value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'HGF source delta',
  title = 'HGF: individual contexts\nSeparate source scale, +/-0.3')
heat_pdf_png(left + right, 'Y_programs_across_drugs_and_HGF', 20, 9)

# A readable focal contrast for the raw ERK genes and their reconstruction from program means.
focal_gene <- rbind(
  dr[recipient == 'BT-474' & drug == 'lapatinib' & readout == 'ERK_output',
    .(context = paste('BT474 lapatinib', dose, 'uM'), source = 'Tahoe', bank_genes, program_reconstruction)],
  lr[treatment == 'HGF' & study == 'GSE129075' & readout == 'ERK_output',
    .(context = paste('BxPC3 HGF', duration), source = 'CytoSig', bank_genes, program_reconstruction)],
  lr[treatment == 'HGF' & study == 'GSE135552' & readout == 'ERK_output',
    .(context = paste(sub('Hnscc Cell Line ', '', condition), 'HGF16h'), source = 'CytoSig', bank_genes, program_reconstruction)])
save_table(focal_gene, 'focal_ERK_signal_retention')
fg <- melt(focal_gene, id.vars = c('context', 'source'), variable.name = 'representation', value.name = 'value')
fg[, context := factor(context, levels = unique(focal_gene$context))]
fg[, representation := factor(representation, levels = c('bank_genes', 'program_reconstruction'),
  labels = c('Original ERK genes in bank43', 'Same genes replaced by program means'))]
p <- ggplot(fg, aes(context, value, fill = representation)) +
  geom_hline(yintercept = 0, color = 'grey50', linewidth = .4) +
  geom_col(position = position_dodge(.8), width = .75) + facet_grid(. ~ source, scales = 'free_x', space = 'free_x') +
  scale_fill_manual(values = c('#1F7391', '#D89539')) +
  labs(x = NULL, y = 'Mean signed response (each source scale)', fill = NULL,
    title = 'Expected ERK-response directions can disappear in the program means',
    subtitle = 'Identical five genes on both sides: PHLDA1, SPRY2, DUSP6, ETV4, ETV5. DUSP4 is outside the bank.') +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p, 'focal_ERK_signal_retention', 11, 4.5)

# Gene-level evidence for the response themes, separate from a pathway label.
theme_programs <- c('P05', 'P09', 'P13', 'P17', 'P30')
gene_evidence <- list()
for (p in theme_programs) {
  genes <- bank$genes[as.numeric(bank$weights[, p]) > 0]
  z <- md[drug %in% c('lapatinib', 'afatinib', 'crizotinib (hydrochloride)', 'capmatinib', 'alpelisib')]
  g <- as.data.table(gene_delta[z$observation_id, genes, drop = FALSE], keep.rownames = 'observation_id')
  g <- melt(g, id.vars = 'observation_id', variable.name = 'gene', value.name = 'delta', variable.factor = FALSE)
  g <- merge(g, z[, .(observation_id, drug, recipient, dose, measured_pct)], by = 'observation_id')
  g[, program := p]
  gene_evidence[[p]] <- g
}
ge <- rbindlist(gene_evidence)
gene_summary <- ge[, .(mean_delta = mean(delta), gene_function_r = cor_safe(delta, log2(measured_pct / 100))),
  by = .(drug, dose, program, gene)]
save_table(gene_summary, 'themed_program_gene_evidence')
save_table(gene_summary[, .(genes = .N, mean_delta = mean(mean_delta), median_delta = median(mean_delta),
  same_direction_fraction = mean(sign(mean_delta) == sign(mean(mean_delta))),
  most_changed_genes = paste(head(gene[order(-abs(mean_delta))], 10), collapse = ', ')),
  by = .(drug, dose, program)], 'themed_program_gene_coherence')
ligand_gene_summary <- rbindlist(lapply(theme_programs, function(p) {
  genes <- intersect(bank$genes[as.numeric(bank$weights[, p]) > 0], colnames(cy_matrix))
  a <- as.data.table(cy_matrix[hz$profile_id, genes, drop = FALSE], keep.rownames = 'profile_id')
  a <- melt(a, id.vars = 'profile_id', variable.name = 'gene', value.name = 'delta', variable.factor = FALSE)
  merge(a[, program := p], hz, by = 'profile_id')
}))
save_table(ligand_gene_summary, 'HGF_themed_program_genes')
save_table(ligand_gene_summary[, .(observed_genes = sum(is.finite(delta)), mean_delta = mean(delta, na.rm = TRUE),
  median_delta = median(delta, na.rm = TRUE), same_direction_fraction = mean(sign(delta) == sign(mean(delta, na.rm = TRUE)), na.rm = TRUE),
  most_changed_genes = paste(head(gene[order(-abs(delta))], 10), collapse = ', ')),
  by = .(profile_id, condition, duration, study, program)], 'HGF_themed_program_gene_coherence')
source_notes <- data.table(
  topic = c('Y mechanism', 'CytoSig units', 'HGF pancreatic expression', 'HGF receptor context', 'HGF RNA exposure correction',
    'Shared HGF experiment', 'ERK reference genes'),
  url = c('https://pmc.ncbi.nlm.nih.gov/articles/PMC3724525/',
    'https://www.nature.com/articles/s41592-021-01274-5',
    'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3692977',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC7013520/',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC7699574/',
    'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26736',
    'https://aacrjournals.org/clincancerres/article/23/6/1471/80485/Clinically-Viable-Gene-Expression-Assays-with'))
save_table(source_notes, 'biological_source_notes')
extra_paths <- c(panel$input_paths,
  file.path(prism_root, c('secondary-screen-replicate-collapsed-treatment-info.csv',
    'secondary-screen-replicate-collapsed-logfold-change.csv', 'secondary-screen-replicate-treatment-info.csv',
    'secondary-screen-logfold-change.csv')),
  '/work/results/program_pilot/55_current_bank_within_treatment_signal_audit/data/diagnostic_predictions.tsv.gz',
  '/work/src/R/lib_heatmap.R')
save_table(data.table(path = extra_paths, sha256 = vapply(extra_paths, function(p) digest(file = p, algo = 'sha256'), character(1))),
  'additional_input_manifest')
figure_index <- data.table(
  file = c('Y_programs_across_drugs_and_HGF.pdf', 'focal_ERK_signal_retention.pdf',
    'baseline_receptor_to_response.pdf', 'HGF_all_programs.pdf', 'HGF_named_genes.pdf',
    'Y_downstream_drug_responses.pdf', 'Y_ligand_downstream_responses.pdf', 'biological_signal_retention.pdf'),
  purpose = c('All34 programs: drug dose means plus individual HGF experiments; same program order',
    'Exact five-gene ERK response versus their program-mean reconstruction',
    'Baseline receptor to molecular response with PRISM abundance', 'Every HGF experiment and program',
    'Named HGF response genes, same experiment order', 'Four drug branches, dose, downstream biology and sensitivity',
    'HGF, EGF and NRG1 readouts with early/late exposure', 'Readout preservation across all drug cells and doses'))
save_table(figure_index, 'figure_index')

writeLines(capture.output(sessionInfo()), file.path(root, 'session_info.txt'))
save_table(data.table(script = '/work/src/program_pilot/57_erbb_met_hgf_program_biology.R',
  sha256 = digest(file = '/work/src/program_pilot/57_erbb_met_hgf_program_biology.R', algo = 'sha256')), 'script_provenance')
writeLines(c('57: ERBB / MET / HGF biology using frozen bank43, no refitting.',
  paste('Drug profiles:', nrow(md), '; ligand contrasts:', nrow(cm), '; exact saved-score error:', agreement),
  'Drug and ligand heatmaps use original signed response values, with explicit saturation limits.',
  'No matched drug-plus-HGF expression measurements; cross-source pattern agreement cannot establish the Y topology.',
  'Bank43 is outcome-informed; all phenotype analyses here are retrospective diagnostics.',
  'See INTERPRETATION.txt for the inspected biological conclusions.'), file.path(root, 'REPORT.txt'))
message('Completed ', root)
