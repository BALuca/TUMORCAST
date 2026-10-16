#!/usr/bin/env Rscript
# Measured external bypass evidence for analysis 58. No fitting or gene selection.
# Reuse the frozen historical 54-gene definition and its verified GEO reader cache.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source('/work/src/R/lib_heatmap.R')
out <- '/work/results/program_pilot/58_molecular_function_information_audit'
old <- '/work/tahoe_prism_v2/results/35_frozen_program_bypass_validation'
cache <- '/work/tahoe_prism_v2/artifacts/33_braf_erbb_programs/verification_inputs.rds'
dir.create(file.path(out, 'data'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, 'figures'), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(out, 'data', paste0('external_', name, '.tsv')), sep = '\t')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(out, 'figures', paste0('external_', name, '.pdf')), p,
         width = width, height = height, device = cairo_pdf, bg = 'white')
  ggsave(file.path(out, 'figures', paste0('external_', name, '.png')), p,
         width = width, height = height, dpi = 180, bg = 'white')
}
v <- readRDS(cache)
membership <- fread(file.path(old, 'frozen_program_genes_exact_copy.tsv'))
stopifnot(nrow(membership) == 54L, !anyDuplicated(membership$gene))
genes <- membership$gene
program_labels <- c(P1 = 'P1: proliferation suppression', P2 = 'P2: adaptive response', P3 = 'P3: shared output')
write_table(membership, 'frozen_membership')

# Gene-level time courses preserve both the drug-induced suppression relative
# to vehicle and the ligand increment relative to each pretreated background.
ma <- as.data.table(v$full_time)
stopifnot(identical(colnames(v$absolute127), ma$gsm))
da <- as.data.table(v$absolute127[genes, , drop = FALSE], keep.rownames = 'gene')
ga <- melt(da, id.vars = 'gene', variable.name = 'gsm', value.name = 'delta_vehicle')
ga <- merge(ga, ma[, .(gsm, arm, egf_h, biological_set)], by = 'gsm')
ga <- merge(ga, membership[, .(gene, program, orientation)], by = 'gene')
ga[, signed_delta := delta_vehicle * orientation]
gene_day <- ga[, .(delta_vehicle = mean(delta_vehicle), signed_delta = mean(signed_delta),
                   n_technical = .N), by = .(gene, program, arm, egf_h, biological_set)]
write_table(gene_day, 'egf_gene_timecourses')
program_day <- gene_day[, .(score = mean(signed_delta, na.rm = TRUE),
                            measured_genes = sum(is.finite(signed_delta))),
                       by = .(program, arm, egf_h, biological_set)]
program_day[, program_label := factor(program, levels = names(program_labels), labels = program_labels)]
write_table(program_day, 'egf_program_timecourses_vehicle_reference')
p <- ggplot(program_day[arm != 'DMSO'], aes(egf_h, score, color = arm,
             group = interaction(arm, biological_set), linetype = factor(biological_set))) +
  geom_hline(yintercept = 0, color = 'grey75') + geom_line(linewidth = .7) + geom_point(size = 1.7) +
  facet_wrap(~program_label, nrow = 1, scales = 'free_y') +
  scale_color_manual(values = c('Vem' = '#BB6530', 'Vem + Cobi' = '#63539B')) +
  theme_classic(base_size = 11) + theme(legend.position = 'bottom') +
  labs(x = 'Hours after EGF addition', y = 'Signed RNA change from vehicle', color = 'Drug background',
       linetype = 'Biological day', title = 'Measured ligand response during BRAF inhibition',
       subtitle = 'A375: 24 h drug pretreatment; two biological days, two technical wells per condition')
save_plot(p, 'egf_bypass_timecourses', 9, 4)

# Quantify ligand induction and its blockade for every measured program member.
mi <- as.data.table(v$s127)
di <- as.data.table(v$d127[genes, , drop = FALSE], keep.rownames = 'gene')
gi <- melt(di, id.vars = 'gene', variable.name = 'gsm', value.name = 'increment')
gi <- merge(gi, mi[, .(gsm, arm, egf_h, biological_set)], by = 'gsm')
gi <- merge(gi, membership[, .(gene, program, orientation)], by = 'gene')
gi[, signed_increment := increment * orientation]
gd <- gi[, .(score = mean(signed_increment)), by = .(gene, program, arm, egf_h, biological_set)]
setorder(gd, gene, program, arm, biological_set, egf_h)
auc <- gd[, .(auc = sum(diff(egf_h) * (head(score, -1) + tail(score, -1)) / 2)),
          by = .(gene, program, arm, biological_set)]
aw <- dcast(auc, gene + program + biological_set ~ arm, value.var = 'auc')
setnames(aw, c('Vem', 'Vem + Cobi'), c('vem_auc', 'vem_cobi_auc'))
aw[, blockade := vem_auc - vem_cobi_auc]
write_table(aw, 'egf_gene_auc')
coherence <- aw[is.finite(vem_auc) & is.finite(vem_cobi_auc),
  .(genes = .N, induced_genes = sum(vem_auc > 0), blocked_genes = sum(blockade > 0),
    induced_and_blocked = sum(vem_auc > 0 & blockade > 0),
    mean_vem_auc = mean(vem_auc), mean_vem_cobi_auc = mean(vem_cobi_auc)), by = .(program, biological_set)]
write_table(coherence, 'egf_gene_coherence')
repeat_auc <- dcast(aw, gene + program ~ biological_set, value.var = c('vem_auc', 'blockade'))
repeat_stats <- repeat_auc[, .(n = sum(complete.cases(vem_auc_3, vem_auc_4)),
  induction_day_r = cor(vem_auc_3, vem_auc_4, use = 'complete.obs'),
  blockade_day_r = cor(blockade_3, blockade_4, use = 'complete.obs')), by = program]
write_table(repeat_stats, 'egf_day_reproducibility')

# The full frozen membership is displayed, including the one unavailable gene.
gd[, profile := paste(arm, biological_set, egf_h, sep = '__')]
col_meta <- unique(gd[, .(profile, arm, biological_set, egf_h)])
col_meta[, arm := factor(arm, levels = c('Vem', 'Vem + Cobi'))]
setorder(col_meta, arm, biological_set, egf_h)
hm_dt <- dcast(gd, gene ~ profile, value.var = 'score')
hm <- as.matrix(hm_dt[match(genes, gene), col_meta$profile, with = FALSE])
rownames(hm) <- genes
ann <- data.frame(Background = col_meta$arm, Day = factor(col_meta$biological_set),
                  Hour = factor(col_meta$egf_h), row.names = col_meta$profile)
stopifnot(identical(colnames(hm), rownames(ann)))
limit <- unname(quantile(abs(hm), .98, na.rm = TRUE))
h <- draw_heatmap(hm, col_anno_df = ann, col_anno_cols = names(ann),
  scale_rows = FALSE, value_range = c(-limit, 0, limit),
  value_colors = c('#2166AC', 'white', '#B2182B'), na_col = '#AAAAAA',
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
  cluster_column_slices = FALSE, row_split = factor(membership$program, levels = c('P1','P2','P3')),
  column_split = factor(paste(col_meta$arm, col_meta$biological_set),
                        levels = unique(paste(col_meta$arm, col_meta$biological_set))),
  show_column_names = FALSE, show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 7),
  hm_name = 'Signed RNA', heatmap_legend_param = list(title = 'EGF increment\nlog2 scale'))
for (ext in c('pdf', 'png')) {
  path <- file.path(out, 'figures', paste0('external_egf_all_frozen_genes.', ext))
  if (ext == 'pdf') cairo_pdf(path, 10, 9) else png(path, 10, 9, units = 'in', res = 180, type = 'cairo')
  ComplexHeatmap::draw(h, heatmap_legend_side = 'right', annotation_legend_side = 'right',
    column_title = 'Measured EGF response during BRAF inhibition; all frozen genes')
  dev.off()
}

# Verify historical five-line functional association from the measured arrays.
# Effects are signed changes from each cell's vehicle controls, not predictions.
m684 <- as.data.table(v$s684)
stopifnot(identical(colnames(v$d684), m684$gsm))
gene684 <- melt(as.data.table(v$d684[genes, , drop = FALSE], keep.rownames = 'gene'),
                id.vars = 'gene', variable.name = 'gsm', value.name = 'delta')
gene684 <- merge(gene684, m684[, .(gsm, cell_line, arm, replicate)], by = 'gsm')
gene684 <- merge(gene684, membership[, .(gene, program, orientation)], by = 'gene')
gene684[, signed_delta := delta * orientation]
scores <- gene684[, .(score = mean(signed_delta)), by = .(gsm, cell_line, arm, replicate, program)]
means <- scores[, .(score = mean(score), se = sd(score) / sqrt(.N), n = .N), by = .(cell_line, arm, program)]
fact <- dcast(means, cell_line + program ~ arm, value.var = 'score')
fact <- fact[complete.cases(fact)]
fact[, `:=`(interaction_suppression = -(plx_lap - plx - lapatinib + control),
             combo_suppression = -(plx_lap - control),
             plx_suppression = -(plx - control),
             lap_suppression = -(lapatinib - control),
             added_to_plx_suppression = -(plx_lap - plx))]
historical <- fread(file.path(old, 'gse68453_frozen_program_function.tsv'))
bliss <- unique(historical[, .(cell_line, bliss)])
fact <- merge(fact, bliss, by = 'cell_line')
check <- merge(fact[, .(cell_line, program, interaction_suppression)],
               historical[, .(cell_line, program, suppression)], by = c('cell_line', 'program'))
stopifnot(max(abs(check$interaction_suppression - check$suppression)) < 1e-10)
write_table(fact, 'combo_factorial')
write_table(means, 'combo_program_arm_means')
write_table(gene684, 'combo_individual_gene_measurements')
long <- melt(fact, id.vars = c('cell_line', 'program', 'bliss'),
  measure.vars = c('interaction_suppression', 'combo_suppression', 'plx_suppression',
                  'lap_suppression', 'added_to_plx_suppression'), variable.name = 'readout', value.name = 'value')
# All 120 permutations of five cell labels: no asymptotic significance claim.
perm <- as.matrix(expand.grid(rep(list(seq_len(5)), 5)))
perm <- perm[apply(perm, 1, function(z) length(unique(z)) == 5), , drop = FALSE]
stats <- long[, {
  stopifnot(.N == 5)
  r <- cor(value, bliss)
  null_r <- apply(perm, 1, function(ii) cor(value, bliss[ii]))
  loo_r <- vapply(seq_len(.N), function(i) cor(value[-i], bliss[-i]), numeric(1))
  list(n_cells = .N, r = r, exact_two_sided_p = mean(abs(null_r) >= abs(r) - 1e-12),
       loo_r_min = min(loo_r), loo_r_max = max(loo_r),
       r_without_A375 = cor(value[cell_line != 'A375'], bliss[cell_line != 'A375']))
}, by = .(program, readout)]
write_table(stats, 'combo_association_robustness')
loo <- long[, rbindlist(lapply(seq_len(.N), function(i)
  data.table(omitted_cell = cell_line[i], r = cor(value[-i], bliss[-i])))), by = .(program, readout)]
write_table(loo, 'combo_leave_one_cell_out')
arm_uncertainty <- dcast(means, cell_line + program ~ arm, value.var = 'se')
arm_uncertainty <- arm_uncertainty[complete.cases(arm_uncertainty)]
arm_uncertainty[, interaction_se := sqrt(control^2 + plx^2 + lapatinib^2 + plx_lap^2)]
fact_se <- merge(fact, arm_uncertainty[, .(cell_line, program, interaction_se)], by = c('cell_line', 'program'))
fact_se[, molecular_interaction_over_two_se := abs(interaction_suppression) > 2 * interaction_se]
write_table(fact_se, 'combo_molecular_interaction_uncertainty')
gene_means <- gene684[, .(score = mean(signed_delta)), by = .(gene, program, cell_line, arm)]
gene_fact <- dcast(gene_means, gene + program + cell_line ~ arm, value.var = 'score')
gene_fact <- merge(gene_fact[complete.cases(gene_fact)], bliss, by = 'cell_line')
gene_fact[, interaction_suppression := -(plx_lap - plx - lapatinib + control)]
gene_associations <- gene_fact[, .(n = .N, r = cor(interaction_suppression, bliss)), by = .(gene, program)]
write_table(gene_associations, 'combo_gene_associations')
write_table(gene_associations[, .(genes = .N, positive_genes = sum(r > 0), median_gene_r = median(r)), by = program],
            'combo_gene_coherence')

# All twelve lines have PLX-only RNA, extending the same frozen scores to a
# larger comparison without choosing genes or signs using functional outcomes.
all_benefit <- unique(fread('/work/tahoe_prism_v2/results/33_braf_erbb_programs/gse68453_expression_bliss_join.tsv')
                     [, .(cell_line, bliss)])
plx_all <- merge(means[arm == 'plx'], all_benefit, by = 'cell_line')
plx_all[, published_resistant := cell_line %chin% c('SKMEL28','SKMEL24','DM13','VMM5A','VMM18','DM331')]
plx_stats <- plx_all[, {
  ct <- cor.test(score, bliss)
  lr <- vapply(seq_len(.N), function(i) cor(score[-i], bliss[-i]), numeric(1))
  list(n = .N, r = unname(ct$estimate), p = ct$p.value, low = ct$conf.int[1], high = ct$conf.int[2],
       loo_min = min(lr), loo_max = max(lr), resistant_mean = mean(score[published_resistant]),
       sensitive_mean = mean(score[!published_resistant]))
}, by = program]
write_table(plx_all, 'twelve_line_plx_response_and_benefit')
write_table(plx_stats, 'twelve_line_plx_association')
plx_all[, program_label := factor(program, levels = names(program_labels), labels = program_labels)]
plx_stats[, program_label := factor(program, levels = names(program_labels), labels = program_labels)]
plx_stats[, label := sprintf('Pearson r = %.2f; n = %d', r, n)]
p <- ggplot(plx_all, aes(score, bliss, color = published_resistant)) + geom_point(size = 2) +
  geom_text(aes(label = cell_line), size = 2.5, vjust = -1, check_overlap = TRUE) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = '#406080', linewidth = .6) +
  geom_label(data = plx_stats, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
             hjust = -.02, vjust = 1.05, size = 2.7, label.size = 0) +
  facet_wrap(~program_label, nrow = 1, scales = 'free_x') + theme_classic(base_size = 11) +
  theme(aspect.ratio = 1, legend.position = 'bottom') +
  scale_y_continuous(expand = expansion(mult = c(.08, .3))) +
  labs(x = 'Measured 8 h PLX-only program displacement', y = 'Later added combination benefit (Bliss)',
       color = 'Published PLX resistance', title = 'A single-agent response compared across twelve melanoma lines',
       subtitle = 'Frozen historical gene sets; all twelve lines shown; retrospective association, not a fitted predictor')
save_plot(p, 'twelve_line_plx_response', 11, 4.8)
gp <- long[readout == 'interaction_suppression']
gp <- merge(gp, fact_se[, .(cell_line, program, interaction_se)], by = c('cell_line', 'program'))
gp[, program_label := factor(program, levels = names(program_labels), labels = program_labels)]
labels <- stats[readout == 'interaction_suppression']
labels[, program_label := factor(program, levels = names(program_labels), labels = program_labels)]
labels[, label := sprintf('r = %.2f; exact p = %.2f\nleave-one-cell-out r: %.2f to %.2f', r, exact_two_sided_p,
                          loo_r_min, loo_r_max)]
p <- ggplot(gp, aes(value, bliss)) + geom_point(size = 2) +
  geom_segment(aes(x = value - interaction_se, xend = value + interaction_se, yend = bliss),
               linewidth = .45, color = '#888888') +
  geom_text(aes(label = cell_line), size = 2.7, vjust = -1, check_overlap = TRUE) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#406080', linewidth = .6) +
  geom_label(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
             hjust = -.02, vjust = 1.05, size = 2.7, label.size = 0) +
  facet_wrap(~program_label, nrow = 1, scales = 'free_x') +
  scale_y_continuous(expand = expansion(mult = c(.08, .4))) +
  scale_x_continuous(expand = expansion(mult = c(.15, .2))) +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1) +
  labs(x = '8 h molecular suppression beyond additive drug responses', y = '3-day combination benefit (Bliss)',
       title = 'Does measured molecular interaction track functional interaction?',
       subtitle = 'Five melanoma factorials; horizontal bars = 1 RNA SE; benefit digitized historically; retrospective')
save_plot(p, 'combo_function_robustness', 11, 4.5)

# Preserve actual measured phosphopeptide evidence, clearly separate from HGF.
phospho <- fread(file.path(old, 'gerosa_phosphoproteomic_summary.tsv'))
write_table(phospho, 'egf_measured_phosphoproteins')
pg <- phospho[condition %chin% c('DMSO', '1 µM vem', '1 µM vem + EGF')]
pg[, condition := factor(condition, levels = c('DMSO', '1 µM vem', '1 µM vem + EGF'),
                         labels = c('Vehicle', 'Vem', 'Vem + EGF'))]
p <- ggplot(pg, aes(condition, mean_relative, group = analyte)) +
  geom_hline(yintercept = 1, linetype = 2, color = 'grey70') +
  geom_line(color = '#406080') + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_relative - se_relative, ymax = mean_relative + se_relative),
                width = .1, na.rm = TRUE) + facet_wrap(~analyte, nrow = 1) + theme_classic(base_size = 10) +
  labs(x = NULL, y = 'Measured phosphopeptide / vehicle', title = 'A measured protein trace of EGF bypass',
       subtitle = 'A375 targeted SRM; vehicle and vem n=4; EGF condition n=1; this is EGF, not HGF')
save_plot(p, 'egf_measured_phosphoproteins', 9, 3.6)

# Matched BxPC3 ligand comparison: both HGF and MSP activate motility, but the
# source study reports proliferation only for HGF. MST1 is the CytoSig label
# for secreted MSP here, not the unrelated STK4/MST1 kinase. These four columns
# are processed mean contrasts; no replicate-level uncertainty is available.
cyto_file <- '/work/data/01_CytoSig_data/diff.merge.gz'
cyto_meta_file <- '/work/data/01_CytoSig_data/meta_info.txt'
cm <- fread(cyto_meta_file, skip = 1L, header = FALSE,
            col.names = c('profile_id', 'treatment', 'condition', 'duration', 'dose'))
cm <- cm[grepl('GSE129075', profile_id) & treatment %chin% c('HGF', 'MST1')]
stopifnot(nrow(cm) == 4, !anyDuplicated(cm$profile_id))
cm[, ligand := fifelse(treatment == 'MST1', 'MSP', 'HGF')]
cm[, hours := as.numeric(sub('h$', '', duration))]
setorder(cm, hours, ligand)
cy_names <- names(fread(cyto_file, nrows = 0L))
cy <- fread(cyto_file, select = c(cy_names[1], cm$profile_id), showProgress = FALSE)
setnames(cy, 1L, 'gene')
stopifnot(!anyDuplicated(cy$gene))
cy_long <- melt(cy, id.vars = 'gene', variable.name = 'profile_id', value.name = 'log2_change')
cy_long <- merge(cy_long, cm, by = 'profile_id')
write_table(cm, 'hgf_msp_profiles')
write_table(cy_long, 'hgf_msp_all_genes')
readout_file <- '/work/results/program_pilot/57_erbb_met_hgf_program_biology/data/biological_readout_membership.tsv.gz'
readout_members <- fread(readout_file)[readout %chin% c('ERK_output', 'Immediate_early', 'E2F_TARGETS',
                                                       'MYC_TARGETS_V1', 'MTORC1_SIGNALING')]
cy_readouts <- merge(cy_long, readout_members[, .(gene, readout)], by = 'gene', allow.cartesian = TRUE)
readout_scores <- cy_readouts[, .(score = mean(log2_change, na.rm = TRUE),
  measured_genes = sum(is.finite(log2_change)), positive_genes = sum(log2_change > 0, na.rm = TRUE)),
  by = .(ligand, hours, readout)]
write_table(readout_scores, 'hgf_msp_fixed_readouts')
write_table(cy_readouts, 'hgf_msp_readout_genes')
readout_pair <- dcast(cy_readouts, gene + readout + hours ~ ligand, value.var = 'log2_change')
readout_coherence <- readout_pair[is.finite(HGF) & is.finite(MSP),
  .(jointly_measured_genes = .N, HGF_above_MSP = sum(HGF > MSP),
    mean_difference = mean(HGF - MSP)), by = .(readout, hours)]
write_table(readout_coherence, 'hgf_msp_readout_coherence')
paired <- dcast(cy_long, gene + hours ~ ligand, value.var = 'log2_change')
paired <- paired[is.finite(HGF) & is.finite(MSP)]
paired[, difference := HGF - MSP]
similarity <- paired[, .(genes = .N, pearson = cor(HGF, MSP),
  cosine = sum(HGF * MSP) / sqrt(sum(HGF^2) * sum(MSP^2)),
  hgf_median = median(HGF), msp_median = median(MSP),
  hgf_mean = mean(HGF), msp_mean = mean(MSP),
  hgf_abs_over_half_log2 = sum(abs(HGF) >= .5), msp_abs_over_half_log2 = sum(abs(MSP) >= .5),
  hgf_abs_mean = mean(abs(HGF)), msp_abs_mean = mean(abs(MSP))), by = hours]
write_table(similarity, 'hgf_msp_genome_similarity')
setorder(paired, hours, -difference)
write_table(paired, 'hgf_msp_gene_differences')
p <- ggplot(readout_scores, aes(factor(hours), score, color = ligand, group = ligand)) +
  geom_hline(yintercept = 0, color = 'grey75') + geom_line(linewidth = .7) + geom_point(size = 2) +
  facet_wrap(~readout, scales = 'free_y', nrow = 1) + theme_classic(base_size = 10) +
  scale_color_manual(values = c(HGF = '#B24A30', MSP = '#367E95')) + theme(legend.position = 'bottom') +
  labs(x = 'Hours after ligand', y = 'Mean gene log2 change from vehicle',
       title = 'Same cells, different ligands: can early RNA distinguish their consequences?',
       subtitle = 'BxPC3; HGF80ng/mL vsMSP100ng/mL. Paper: both increase migration; onlyHGF proliferation.')
save_plot(p, 'hgf_msp_fixed_readouts', 11, 3.8)
named_genes <- unique(c(readout_members[readout %chin% c('ERK_output', 'Immediate_early'), gene],
  c('MYC','CCND1','CCNE1','CDK2','CDK4','MKI67','PCNA','MCM2','MCM5','CDC6','E2F1','E2F2','RRM2','TYMS',
    'EIF4EBP1','SLC2A1','HK2','LDHA','SLC7A5','ASNS','DDIT4','HBEGF','EREG','AREG','MET','MST1R')))
ng <- cy[gene %chin% named_genes]
ng <- ng[match(intersect(named_genes, ng$gene), gene)]
matrix_ng <- as.matrix(ng[, cm$profile_id, with = FALSE])
rownames(matrix_ng) <- ng$gene
colnames(matrix_ng) <- paste(cm$ligand, cm$duration)
h <- draw_heatmap(matrix_ng, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 0,
  row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 9),
  value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'),
  hm_name = 'RNA log2FC', heatmap_legend_param = list(title = 'Log2 change\nfrom vehicle'))
for (ext in c('pdf', 'png')) {
  path <- file.path(out, 'figures', paste0('external_hgf_msp_named_genes.', ext))
  if (ext == 'pdf') cairo_pdf(path, 5, 8) else png(path, 5, 8, units = 'in', res = 180, type = 'cairo')
  ComplexHeatmap::draw(h, column_title = 'BxPC3: HGF vsMSP\nFixed named genes; original source changes')
  dev.off()
}
similarity[, label := sprintf('Pearson r=%.2f; n=%s genes', pearson, format(genes, big.mark = ','))]
p <- ggplot(paired, aes(MSP, HGF)) + geom_point(size = .4, alpha = .2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#406080', linewidth = .6) +
  geom_label(data = similarity, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
             hjust = -.02, vjust = 1.05, size = 3, label.size = 0) +
  facet_wrap(~hours, scales = 'free', labeller = label_both) + theme_classic(base_size = 11) +
  theme(aspect.ratio = 1) + labs(x = 'MSP response, log2 change', y = 'HGF response, log2 change',
    title = 'Genome-wide similarity does not establish equal function',
    subtitle = 'All jointly measured genes; no selection by response strength; processed mean contrasts')
save_plot(p, 'hgf_msp_genome_similarity', 8, 4.5)
writeLines(c('Matched BxPC3 molecular comparison: GSE129075.',
  'Primary paper DOI:10.1016/j.cyto.2019.154762.',
  'https://www.sciencedirect.com/science/article/abs/pii/S1043466619301917',
  'Source design: HGF80ng/mL, MSP100ng/mL, vehicle;1h and4hRNA, biologicaltriplicate.',
  'CytoSig MST1 maps to macrophage-stimulating protein in this experiment.',
  'The source publication reports both ligands increase migration, but onlyHGF increases proliferation.',
  'The exact proliferation assay and exposure duration could not be verified from accessible primary text;',
  'do not equate the1h/4hRNA measurements with the later reported functional endpoint.',
  'No numeric functional outcome is digitized or imputed here. This is not a drug-rescue experiment.',
  'The processed CytoSig contrasts lack individualreplicate uncertainty; differential thresholds describe',
  'effect magnitude only, not statistical significance. Readouts and named genes were fixed beforehand.',
  'Magnitude differences can reflect ligand potency/receptor competence; these concentrations are not',
  'demonstrated to be signaling-equivalent. This does not establish a pathway-specific causal decoder.',
  'Observed4hresponses distinguishHGF fromMSP: sustainedDUSP4,DUSP6,PHLDA1,SPRY2,MYC induction.',
  'HGF andMSP share onlymoderate genomecorrelation(~.57); HGF has farlarger andbroader changes.',
  'The4hHGFgenomemean is−.096(median−.070), versusMSP+.002(median0): positiveMYC/mTORreadouts',
  'are not a uniform positivegenomeoffset. Their gene-directioncoherence is nevertheless incomplete.',
  'E2F-targets fall afterHGF(−.328), whileMSPisflat. ETV4/5 fall despiteinduction ofotherERKfootprintgenes.',
  'ThusRNA differs meaningfully, but a genericgrowth/cell-cycleactivity score does not explain the',
  'reportedproliferation distinction. The result supportsresponse heterogeneity, not an identified',
  'functionalcode or drugrescue mechanism.'
), file.path(out, 'data', 'external_hgf_msp_provenance.txt'))
print(readout_scores)
print(similarity)

inputs <- c(cache, cyto_file, cyto_meta_file, readout_file, file.path(old, 'frozen_program_genes_exact_copy.tsv'),
            file.path(old, 'gse68453_frozen_program_function.tsv'),
            file.path(old, 'gerosa_phosphoproteomic_summary.tsv'))
write_table(data.table(path = inputs, md5 = unname(tools::md5sum(inputs))), 'input_manifest')
writeLines(capture.output(sessionInfo()), file.path(out, 'data', 'external_session_info.txt'))
writeLines(c(
  'External molecular bypass evidence: executed audit, 2026-10-11',
  '',
  'GSE127988 provides a real EGF bypass trace during BRAF inhibition in A375.',
  'Two independent biological days, two technical wells per condition. Fixed historical 54-gene programs.',
  'Adaptive P2: 16/19 and15/19 measured genes induced by EGF and attenuated by cobimetinib.',
  'Gene-level induction across-day r=.836; blockade r=.751. SharedP3 induction r=.883, blockade r=.805.',
  'RNA induction is transient and P3 remains suppressed relative to untreated vehicle: partial restoration,',
  'not normalization of the entire molecular state. P1 is signed toward proliferation suppression.',
  'Actual targeted phosphopeptide measurements independently show ERK1 .340→1.085 and ERK2 .156→1.318',
  'relative to vehicle when EGF is added to1uM vemurafenib. EGF protein condition has only one measurement.',
  '',
  'GSE68453 is a separate melanoma dataset with8h RNA and published3-day combination-benefit values.',
  'Historical benefit values were digitized; they are not newly quantified raw functional measurements.',
  'Five complete PLX4720/lapatinib factorials reproduce frozenP3 interaction-benefit r=.729.',
  'Leave-one-cell-out r=.726–.912; exact two-sided permutation p=.217. All14 member genes have positive r.',
  'However four interaction scores are .007–.041 in magnitude versus SE .208–.248. OnlyA375 exceeds2SE,',
  'and its interaction is less-than-additive suppression. There is no resolved super-additive RNA',
  'suppression across the remaining cells. A high correlation of uncertain point estimates was oversold.',
  'P2 r=.452 flips to−.450 whenA375 is omitted. ProliferationP1 does not correlate positively.',
  'Across all12 lines, PLX-only frozenP3 displacement versus later combination benefit is weak: r=.212,p=.508.',
  '',
  'HGF/MET: these EGF/BRAF results are not HGF rescue evidence. Analysis57 shows HGF-alone induction in',
  'unmatched recipients and separate lapatinib suppression. No matchedHGF+drug transcriptional factorial',
  'has been verified among these inputs. SavedWilson ±1 phosphorylation labels are excluded: hand-coded,',
  'not quantitative molecular measurements. Published phosphoprotein/functional HGF rescue is distinct',
  'from having a measured RNA bridge in our present dataset.',
  '',
  'Interpretation: earlyRNA can measure a reproducible intervention-dependent pathway restoration and',
  'blockade. These data do not yet establish a robust transferable RNA→functional-combination map.',
  'No new genes, signs, programs, or predictive models were selected using these outcomes.'
), file.path(out, 'data', 'external_INTERPRETATION.txt'))
print(coherence)
print(repeat_stats)
print(stats[readout == 'interaction_suppression'])
message('External measured-bypass audit complete: ', out)
