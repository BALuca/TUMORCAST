#!/usr/bin/env Rscript
# Complete BT20 combination and constituent-response views, followed by reserved
# prediction diagnostics. Expectations were frozen before models were evaluated.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(2L)
options(warn = 1)
root <- '/work/results/program_pilot/61_individual_responses_combination_transfer'
data_dir <- file.path(root, 'data')
figure_dir <- file.path(root, 'figures')
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- file.path(data_dir, c('matched_combination_observations.tsv.gz',
  'matched_functional_controls.tsv', 'matched_molecular.qs', 'representation_scores.qs',
  'biological_expectations_frozen.tsv'))
stopifnot(all(file.exists(inputs)))
anchors <- c('alpelisib', 'dactolisib')
partners <- c('lapatinib', 'neratinib', 'selumetinib', 'trametinib', 'saracatinib', 'nvp-tae684')
model_levels <- c('Independent_singles', 'Singles_ridge', 'Singles_SVR', 'Singles_MOA_SVR',
                  'Original40_ridge', 'Original40_SVR', 'Biology_SVR', 'GenePC12_SVR')
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))
save_table <- function(x, name) fwrite(x, file.path(data_dir, paste0('figure_', name, '.tsv.gz')), sep = '\t')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(figure_dir, paste0(name, '.png')), p, width = width, height = height, dpi = 150)
}
combo <- fread(inputs[1])
stopifnot(setequal(unique(combo$anchor), anchors), setequal(unique(combo$partner), partners),
          uniqueN(paste(combo$anchor, combo$partner)) == 12L,
          !anyDuplicated(combo$observation_id), all(is.finite(combo$GR_AB)))
combo[, `:=`(anchor = factor(anchor, levels = anchors), partner = factor(partner, levels = partners))]
combo[, dose_color := factor(as.numeric(sub('^[^:]+:', '', dose_key_A)))]
functional <- combo[, .(GR_AB = mean(GR_AB), GR_A = mean(GR_A), GR_B = mean(GR_B),
  SD_GR_AB = sd(GR_AB), excess_inhibition = mean(excess_inhibition),
  SD_extra = sd(excess_inhibition), ratio_AB = mean(ratio_AB), ratio_A = mean(ratio_A), ratio_B = mean(ratio_B),
  n_bio = .N, dose_A_uM = mean(dose_A_uM), dose_B_uM = mean(dose_B_uM),
  repeated_positive = all(excess_inhibition > 0),
  net_loss_repeats = sum(GR_AB < 0),
  single_to_combination_kill_repeats = sum(GR_AB < 0 & GR_A >= 0 & GR_B >= 0)),
  by = .(condition_id, anchor, partner, dose_key_A, dose_key_B)]
functional[, `:=`(dose_color = factor(as.numeric(sub('^[^:]+:', '', dose_key_A))),
  clear_extra_effect = excess_inhibition >= 10 & repeated_positive & n_bio >= 2)]
save_table(functional, 'complete_functional_condition_means')
summary <- functional[, .(dose_contexts = .N, three_bio_contexts = sum(n_bio == 3),
  clear_extra_contexts = sum(clear_extra_effect),
  clear_extra_two_repeat_contexts = sum(clear_extra_effect & n_bio == 2),
  clear_extra_three_repeat_contexts = sum(clear_extra_effect & n_bio == 3),
  all_bio_combo_killing = sum(net_loss_repeats == n_bio & n_bio >= 2),
  all_bio_new_killing_from_nonkilling_singles = sum(single_to_combination_kill_repeats == n_bio & n_bio >= 2),
  mean_extra_inhibition_pp = mean(excess_inhibition),
  median_extra_inhibition_pp = median(excess_inhibition)), by = .(anchor, partner)]
save_table(summary, 'all12_functional_summary')

# Lines show complete measured dose series at each anchor concentration. The
# partner-alone measurements are plotted separately and are never predicted RNA.
functional_curve <- function(tab, title, subtitle) {
  single_b <- tab[, .(GR_B = mean(GR_B)), by = .(anchor, partner, dose_B_uM)]
  ggplot(tab, aes(dose_B_uM, GR_AB, color = dose_color, group = dose_color)) +
    geom_hline(yintercept = 0, color = 'grey60', linewidth = .3) +
    geom_hline(yintercept = 1, color = 'grey60', linetype = 2, linewidth = .3) +
    geom_line(linewidth = .45) + geom_point(size = 1) +
    geom_line(data = single_b, aes(dose_B_uM, GR_B, group = 1), color = 'black',
              linetype = 2, linewidth = .6, inherit.aes = FALSE) +
    scale_x_log10() + facet_grid(partner ~ anchor) +
    labs(title = title, subtitle = subtitle, x = 'Partner concentration (µM)', y = 'GR: 1 unchanged growth; 0 arrest; <0 net cell loss',
      color = 'PI3K-anchor dose (µM)',
      caption = 'Black dashed curve: measured partner alone. Colored curves: measured combination at every available anchor dose.\nAll12 pairs and all measured doses are retained; panel counts are dose contexts, not independent biological experiments.') +
    theme(legend.position = 'bottom', legend.key.width = grid::unit(10, 'mm'),
      strip.text = element_text(size = 9), axis.text = element_text(size = 8))
}
overview <- functional_curve(functional,
  'BT20: the complete observed combination responses', 'Mean of available biological experiments; individual experiments follow in the PDF')
pdf(file.path(figure_dir, 'all12_observed_growth_curves.pdf'), width = 9.5, height = 12)
print(overview)
for (bio in sort(unique(combo$biological_replicate))) {
  print(functional_curve(combo[biological_replicate == bio],
    paste('BT20: biological experiment', bio), 'Technical wells averaged within experiment; same matched single-agent controls'))
}
dev.off()
ggsave(file.path(figure_dir, 'all12_observed_growth_curves.png'), overview, width = 9.5, height = 12, dpi = 150)

# Counts distinguish small growth-rate shifts from the qualitative appearance of
# net loss. Different pairs have different dose coverage, explicitly shown.
summary_long <- melt(summary, id.vars = c('anchor', 'partner', 'dose_contexts'),
  measure.vars = c('clear_extra_contexts', 'all_bio_combo_killing', 'all_bio_new_killing_from_nonkilling_singles'),
  variable.name = 'criterion', value.name = 'contexts')
summary_long[, criterion := factor(criterion,
  levels = c('clear_extra_contexts', 'all_bio_combo_killing', 'all_bio_new_killing_from_nonkilling_singles'),
  labels = c('Extra inhibition >=10 pp\npositive in every available repeat', 'Combination cell loss\nin every available repeat',
    'Combination-only net cell loss\nin every available repeat'))]
p <- ggplot(summary_long, aes(partner, contexts, fill = anchor)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_text(aes(label = paste0(contexts, '/', dose_contexts)), position = position_dodge(width = .8),
            vjust = -.3, size = 2.8) + facet_wrap(~criterion, ncol = 1, scales = 'free_y') +
  scale_y_continuous(expand = expansion(mult = c(0, .17))) + theme(legend.position = 'bottom') +
  labs(title = 'Which effects reproduce across biological experiments?',
    subtitle = 'Qualifying / all measured contexts; at least2 repeats required; available2- and3-repeat counts saved separately',
    x = NULL, y = 'Dose contexts', fill = 'PI3K anchor')
save_plot(p, 'all12_functional_consistency', 10.3, 8.2)
controls <- fread(inputs[2])
save_table(controls, 'control_counts')
p <- ggplot(controls, aes(factor(biological_replicate), provided_vehicle_count, color = factor(technical_replicate))) +
  geom_point(position = position_dodge(width = .25), size = 2.3) +
  geom_point(aes(y = provided_initial_count), position = position_dodge(width = .25), shape = 1, size = 2.3) +
  labs(title = 'The biological experiments started at different cell densities',
    subtitle = 'Filled: supplied72h vehicle count; open: supplied initial count; normalization is preserved',
    x = 'Biological experiment', y = 'Counted cells', color = 'Technical repeat') +
  theme(legend.position = 'bottom')
save_plot(p, 'functional_control_counts', 6.6, 4.3)
individual <- fread(file.path(data_dir, 'matched_same_screen_individual_function.tsv.gz'))
individual <- individual[drug_a %in% c(anchors, partners)]
individual[, drug_a := factor(drug_a, levels = c(anchors, partners))]
p <- ggplot(individual, aes(dose_A_uM, GR, color = factor(biological_replicate), group = biological_replicate)) +
  geom_hline(yintercept = 0, color = 'grey55', linewidth = .3) +
  geom_hline(yintercept = 1, color = 'grey55', linetype = 2, linewidth = .3) +
  geom_line(linewidth = .6) + geom_point(size = 1.4) + scale_x_log10() + facet_wrap(~drug_a, ncol = 4) +
  labs(title = 'The measured individual-drug responses supplied to the models',
    subtitle = 'Same combination screen; biological repeats remain separate; GR above1 indicates faster growth than the supplied vehicle reference',
    x = 'Drug concentration (µM)', y = 'GR: 0 growth arrest; negative net cell loss', color = 'Biological repeat') +
  theme(legend.position = 'bottom')
save_plot(p, 'individual_drug_functional_controls', 11, 6.4)

# Original signed program projection, not a replacement bank or hard-cluster mean.
mol <- qread(inputs[3])
scores <- qread(inputs[4])
stopifnot(all(c('original40', 'metadata') %in% names(scores)))
smd <- as.data.table(scores$metadata)
stopifnot(identical(rownames(scores$original40), smd$profile_id))
selected <- smd[drug %in% c(anchors, partners)]
stopifnot(setequal(selected$drug, c(anchors, partners)), all(c(3, 24) %in% selected$time_h))
selected[, drug := factor(drug, levels = c(anchors, partners))]
setorder(selected, time_h, drug, dose_uM)
program_matrix <- t(scores$original40[selected$profile_id, , drop = FALSE])
stopifnot(nrow(program_matrix) == 40L, all(is.finite(program_matrix)))
row_order <- hclust(dist(program_matrix), method = 'ward.D2')$order
program_matrix <- program_matrix[row_order, , drop = FALSE]
program_annotation <- fread(file.path(data_dir, 'original40_program_annotations.tsv'))
program_labels <- program_annotation$display_name[match(rownames(program_matrix), program_annotation$program)]
stopifnot(!anyNA(program_labels))
program_labels <- ifelse(nchar(program_labels) > 54, paste0(substr(program_labels, 1, 51), '...'), program_labels)
save_table(data.table(feature = rownames(program_matrix), display_order = seq_len(nrow(program_matrix))), 'program_display_order')
program_long <- as.data.table(as.table(program_matrix))
setnames(program_long, c('feature', 'profile_id', 'score'))
program_long <- merge(program_long, selected[, .(profile_id, drug, dose_uM, time_h)], by = 'profile_id')
save_table(program_long, 'original40_single_responses')
hm_limit <- quantile(abs(program_matrix), .98)
stopifnot(hm_limit > 0)
pdf(file.path(figure_dir, 'original40_individual_response_heatmaps.pdf'), width = 13, height = 9)
for (t in c(3, 24)) {
  tab <- selected[time_h == t]
  mm <- program_matrix[, tab$profile_id, drop = FALSE]
  colnames(mm) <- paste(as.character(tab$drug), tab$dose_uM, 'µM')
  ann <- data.frame(Drug = tab$drug, row.names = colnames(mm))
  stopifnot(identical(colnames(mm), rownames(ann)))
  hm <- draw_heatmap(mm, col_anno_df = ann, col_anno_cols = 'Drug', scale_rows = FALSE,
    value_range = c(-hm_limit, 0, hm_limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'Weighted\nresponse', show_row_names = TRUE, show_column_names = TRUE,
    row_labels = program_labels,
    cluster_rows = FALSE, cluster_columns = FALSE, column_split = tab$drug, cluster_column_slices = FALSE,
    column_labels = paste(tab$dose_uM, 'µM'),
    row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 6), column_names_rot = 70,
    column_title_gp = grid::gpar(fontsize = 10),
    title = paste('Original40 signed programs: BT20 individual treatments at', t, 'h\nObserved978-gene contributions only; identical row order and color scale at both times'))
  ComplexHeatmap::draw(hm)
  png(file.path(figure_dir, paste0('original40_individual_response_heatmap_', t, 'h.png')),
      width = 1950, height = 1350, res = 150)
  ComplexHeatmap::draw(hm)
  dev.off()
}
dev.off()

# These named transcript readouts were frozen from the primary study and general
# pathway roles before predictions. They are not phosphoprotein or nuclear FOXO measurements.
gene_groups <- list('Receptors / feedback' = c('EGFR', 'ERBB2', 'ERBB3', 'IGF1R', 'IRS1', 'IRS2', 'MET'),
  'ERK-response transcripts' = c('DUSP4', 'DUSP6', 'PHLDA1', 'SPRY2', 'ETV4', 'ETV5', 'FOS', 'JUN'),
  'Cell cycle / stress transcripts' = c('CDKN1A', 'CDKN1B', 'CCND1', 'MYC', 'BCL2L11', 'XBP1', 'ATF3', 'DDIT3', 'HSPA5'))
gene_annotation <- rbindlist(lapply(names(gene_groups), function(term) data.table(gene = gene_groups[[term]], group = term)))
gene_annotation[, measured := gene %in% rownames(mol$displacement)]
save_table(gene_annotation, 'biological_gene_coverage')
genes <- gene_annotation[measured == TRUE, gene]
gm <- mol$displacement[genes, selected$profile_id, drop = FALSE]
gl <- as.data.table(as.table(gm))
setnames(gl, c('gene', 'profile_id', 'displacement'))
gl <- merge(gl, selected[, .(profile_id, drug, dose_uM, time_h)], by = 'profile_id')
save_table(gl, 'named_single_gene_responses')
gene_limit <- quantile(abs(gm), .98)
pdf(file.path(figure_dir, 'named_gene_individual_response_heatmaps.pdf'), width = 13, height = 6.3)
for (t in c(3, 24)) {
  tab <- selected[time_h == t]
  mm <- gm[, tab$profile_id, drop = FALSE]
  colnames(mm) <- paste(as.character(tab$drug), tab$dose_uM, 'µM')
  gene_ann <- gene_annotation[match(rownames(mm), gene)]
  stopifnot(identical(rownames(mm), gene_ann$gene))
  gene_col_ann <- data.frame(Drug = tab$drug, row.names = colnames(mm))
  hm <- draw_heatmap(mm, col_anno_df = gene_col_ann, col_anno_cols = 'Drug',
    scale_rows = FALSE, value_range = c(-gene_limit, 0, gene_limit),
    value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'RNA\ndisplacement',
    show_row_names = TRUE, show_column_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = factor(gene_ann$group, levels = names(gene_groups)), cluster_row_slices = FALSE,
    column_split = tab$drug, cluster_column_slices = FALSE,
    column_labels = paste(tab$dose_uM, 'µM'),
    row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 6), column_names_rot = 70,
    column_title_gp = grid::gpar(fontsize = 10), row_title_gp = grid::gpar(fontsize = 8),
    title = paste('BT20 measured transcript readouts at', t, 'h\nLog2(1+fluorescence), quantile normalized, minus own-plate vehicle; no row scaling'))
  ComplexHeatmap::draw(hm)
  png(file.path(figure_dir, paste0('named_gene_individual_response_heatmap_', t, 'h.png')),
      width = 1950, height = 945, res = 150)
  ComplexHeatmap::draw(hm)
  dev.off()
}
dev.off()
anchor_genes <- gl[drug %in% anchors & gene %in% c('EGFR', 'ERBB3', 'DUSP4', 'DUSP6', 'CCND1', 'MYC')]
anchor_genes[, technical_SE := vapply(seq_len(.N), function(i)
  mol$displacement_assay_replicate_se[gene[i], profile_id[i]], numeric(1))]
save_table(anchor_genes, 'anchor_gene_time_courses')
p <- ggplot(anchor_genes, aes(dose_uM, displacement, color = factor(time_h), group = time_h)) +
  geom_hline(yintercept = 0, color = 'grey60', linewidth = .3) + geom_line() + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = displacement - 2 * technical_SE, ymax = displacement + 2 * technical_SE), width = .04) +
  scale_x_log10() + facet_grid(gene ~ drug, scales = 'free_y') + theme(legend.position = 'bottom') +
  labs(title = 'The two PI3K anchors have distinguishable measured molecular responses',
    subtitle = 'Prespecified receptor, ERK-response and cell-cycle transcripts; all measured concentrations',
    x = 'Anchor concentration (µM)', y = 'RNA displacement from vehicle', color = 'RNA collection (h)',
    caption = 'Bars are +/-2SE across technical RNA measurements, not biological-replicate confidence intervals.\nERBB3 induction and DUSP responses are transcript observations; no receptor phosphorylation or causal rescue is inferred.')
save_plot(p, 'anchor_gene_time_courses', 8.7, 10.7)

# This script is rerunnable as outputs arrive. Molecular/function figures are
# useful before model fitting ends; explicit status prevents partial work being
# mistaken for a finished prediction test.
prediction_path <- file.path(data_dir, 'predictions.tsv.gz')
if (!file.exists(prediction_path)) {
  writeLines('Observed function and molecular figures complete. Reserved predictions are not yet available; prediction diagnostics pending.',
             file.path(root, 'FIGURE_STATUS.txt'))
  message('Observed data figures complete; prediction diagnostics pending.')
  quit(save = 'no', status = 0)
}
pred <- fread(prediction_path)
required <- c('anchor', 'partner', 'time_h', 'model', 'split', 'condition_id', 'biological_replicate',
              'observed_extra_pp', 'predicted_extra_pp', 'relative_AB', 'predicted_relative_AB', 'GR_AB', 'predicted_GR_AB')
stopifnot(all(required %in% names(pred)), all(is.finite(pred$observed_extra_pp)),
          all(is.finite(pred$predicted_extra_pp)), all(pred$model %in% model_levels))
held <- pred[split == 'held_partner']
stopifnot(nrow(held) > 0L, !anyDuplicated(held[, .(condition_id, biological_replicate, time_h, model)]))
pm <- held[, .(observed_extra = mean(observed_extra_pp), predicted_extra = mean(predicted_extra_pp),
  observed_relative = mean(relative_AB), predicted_relative = mean(predicted_relative_AB),
  observed_GR = mean(GR_AB), predicted_GR = mean(predicted_GR_AB), n_bio = .N,
  observed_min_extra = min(observed_extra_pp), observed_max_extra = max(observed_extra_pp)),
  by = .(condition_id, anchor, partner, time_h, model)]
pm <- merge(pm, functional[, .(condition_id, dose_A_uM, dose_B_uM, clear_extra_effect)], by = 'condition_id')
pm[, `:=`(model = factor(model, levels = model_levels), partner = factor(partner, levels = partners),
           anchor = factor(anchor, levels = anchors))]
save_table(pm, 'held_partner_condition_predictions')
performance <- pm[, .(dose_contexts = .N, Pearson_r = if (sd(predicted_extra) > 0) cor(observed_extra, predicted_extra) else NA_real_,
  clear_extra_contexts = sum(clear_extra_effect),
  clear_extra_predicted_positive = sum(clear_extra_effect & predicted_extra > 0),
  clear_extra_predicted_ge10pp = sum(clear_extra_effect & predicted_extra >= 10),
  mean_observed_extra_pp = mean(observed_extra), mean_predicted_extra_pp = mean(predicted_extra),
  contexts_predicted_within_bio_range = sum(predicted_extra >= observed_min_extra & predicted_extra <= observed_max_extra)),
  by = .(time_h, model, anchor, partner)]
save_table(performance, 'held_partner_interpretable_results')

# Every reserved partner appears. Scatterplot points are dose-context means,
# avoiding presentation of technical wells as independent validation samples.
pdf(file.path(figure_dir, 'held_partner_observed_vs_predicted.pdf'), width = 12, height = 7.8)
for (t in c(3, 24)) {
  z <- pm[time_h == t]
  lab <- z[, .(label = sprintf('r = %s; n = %d contexts',
    if (sd(predicted_extra) > 0) format(round(cor(observed_extra, predicted_extra), 2), nsmall = 2) else 'NA', .N)), by = model]
  p <- ggplot(z, aes(observed_extra, predicted_extra, color = partner, shape = anchor)) +
    geom_point(size = 1, alpha = .55) +
    geom_smooth(aes(x = observed_extra, y = predicted_extra, group = 1), inherit.aes = FALSE,
      method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .45) +
    geom_hline(yintercept = 0, color = 'grey60', linewidth = .25) + geom_vline(xintercept = 0, color = 'grey60', linewidth = .25) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -.03, vjust = 1.2, size = 2.5) +
    facet_wrap(~model, ncol = 4, scales = 'free') + theme(aspect.ratio = 1, legend.position = 'bottom') +
    labs(title = paste('Excluded-partner combination prediction from', t, 'h individual responses'),
      subtitle = 'All combinations containing the partner excluded under both anchors; molecular and single-agent data remain available',
      x = 'Observed extra inhibition beyond measured singles (percentage points)',
      y = 'Predicted extra inhibition (percentage points)', color = 'Held partner', shape = 'Anchor',
      caption = 'Points average biological repeats within dose context. 3 h is primary; 24 h is prespecified secondary. Published outcomes make this retrospective.')
  print(p)
  ggsave(file.path(figure_dir, paste0('held_partner_observed_vs_predicted_', t, 'h.png')), p, width = 12, height = 7.8, dpi = 150)
}
dev.off()

# Consolidated pair-level comparisons retain each full measured effect range.
# The primary original program model is compared with an identical singles-only
# nonlinear model; no held-out result chooses the representation displayed here.
comparison <- pm[model %in% c('Singles_SVR', 'Original40_SVR')]
comparison_long <- rbindlist(list(
  unique(comparison[, .(condition_id, anchor, partner, time_h, dose_A_uM, dose_B_uM,
                        series = 'Observed', value = observed_extra)]),
  comparison[, .(condition_id, anchor, partner, time_h, dose_A_uM, dose_B_uM,
                  series = as.character(model), value = predicted_extra)]))
comparison_long[, dose_color := factor(signif(dose_A_uM, 3))]
save_table(comparison_long, 'complete_reserved_effect_curves')
pdf(file.path(figure_dir, 'all12_reserved_effect_curves.pdf'), width = 12.5, height = 12)
for (t in c(3, 24)) {
  for (anchor_name in anchors) {
    z <- comparison_long[time_h == t & anchor == anchor_name]
    p <- ggplot(z, aes(dose_B_uM, value, color = dose_color, group = dose_color)) +
      geom_hline(yintercept = 0, color = 'grey60', linewidth = .25) + geom_line(linewidth = .45) + geom_point(size = .9) +
      scale_x_log10() + facet_grid(partner ~ factor(series, levels = c('Observed', 'Singles_SVR', 'Original40_SVR'))) +
      theme(legend.position = 'bottom', strip.text = element_text(size = 8), axis.text = element_text(size = 8)) +
      labs(title = paste(anchor_name, '+ every held-out partner:', t, 'h inputs'),
        subtitle = 'Same dose contexts and axis units for observations, singles-only model and original40 response model',
        x = 'Partner concentration (µM)', y = 'Extra inhibition beyond measured singles (percentage points)',
        color = 'Anchor concentration (µM)',
        caption = 'Full dose curves; biological-repeat means. Positive values indicate greater inhibition than measured-single Bliss expectation.')
    print(p)
    ggsave(file.path(figure_dir, paste0('reserved_effect_curves_', anchor_name, '_', t, 'h.png')),
      p, width = 12.5, height = 12, dpi = 130)
  }
}
dev.off()

# Trace largest prediction failures to the measured individual responses without
# refitting. The rule is fixed: ten largest absolute errors, for each RNA time.
failure <- pm[model == 'Original40_SVR']
failure[, absolute_error_pp := abs(predicted_extra - observed_extra)]
setorder(failure, time_h, -absolute_error_pp)
failure <- failure[, head(.SD, 10), by = time_h]
links <- fread(file.path(data_dir, 'matched_molecular_links.tsv.gz'))
failure_links <- unique(links[, .(condition_id, time_h, profile_A, profile_B)])
failure <- merge(failure, failure_links, by = c('condition_id', 'time_h'), all.x = TRUE)
stopifnot(!anyNA(failure$profile_A), !anyNA(failure$profile_B))
failure[, `:=`(single_A_median_absolute_RNA_change = vapply(profile_A, function(id) median(abs(mol$displacement[, id])), numeric(1)),
               single_B_median_absolute_RNA_change = vapply(profile_B, function(id) median(abs(mol$displacement[, id])), numeric(1)))]
# Magnitudes describe the measured RNA inputs, not prediction-error performance.
save_table(failure, 'largest_reserved_failures_with_inputs')

# The most serious errors are traced to actual measured cell counts and singles.
# Net growth/loss remains visible even when excess-over-Bliss is numerically large.
failure_cases <- merge(failure[, .(condition_id, time_h, absolute_error_pp)],
  held[model == 'Original40_SVR'], by = c('condition_id', 'time_h'))
save_table(failure_cases, 'largest_failures_all_biological_repeats')
lap <- held[time_h == 3 & partner == 'lapatinib' & model %in% c('Singles_SVR', 'Original40_SVR')]
stopifnot(all(c('dose_A_uM', 'dose_B_uM', 'relative_A') %in% names(lap)))
lap_counts <- rbindlist(list(
  unique(lap[, .(condition_id, biological_replicate, anchor, dose_A_uM, dose_B_uM,
    series = 'Observed combination', percent_vehicle = 100 * relative_AB)]),
  unique(lap[, .(condition_id, biological_replicate, anchor, dose_A_uM, dose_B_uM,
    series = 'Anchor alone', percent_vehicle = 100 * relative_A)]),
  lap[, .(condition_id, biological_replicate, anchor, dose_A_uM, dose_B_uM,
    series = as.character(model), percent_vehicle = 100 * predicted_relative_AB)]))
lap_counts[, dose_color := factor(signif(dose_A_uM, 3))]
lap_counts <- merge(lap_counts, combo[, .(condition_id, biological_replicate,
  ratio_technical_sd_A, ratio_technical_sd_AB)], by = c('condition_id', 'biological_replicate'))
lap_counts[, technical_SD_percent := fifelse(series == 'Observed combination', 100 * ratio_technical_sd_AB,
  fifelse(series == 'Anchor alone', 100 * ratio_technical_sd_A, NA_real_))]
save_table(lap_counts, 'lapatinib_actual_counts_and_predictions')
pdf(file.path(figure_dir, 'lapatinib_functional_failure_trace.pdf'), width = 11.5, height = 9)
for (anchor_name in anchors) {
  z <- lap_counts[anchor == anchor_name]
  p <- ggplot(z, aes(dose_B_uM, percent_vehicle, color = dose_color, group = dose_color)) +
    geom_hline(yintercept = 100, color = 'grey65', linetype = 2) + geom_line(linewidth = .6) + geom_point(size = 1.2) +
    geom_errorbar(data = z[is.finite(technical_SD_percent)],
      aes(ymin = percent_vehicle - technical_SD_percent, ymax = percent_vehicle + technical_SD_percent), width = .035) +
    scale_x_log10() + facet_grid(biological_replicate ~ factor(series,
      levels = c('Anchor alone', 'Observed combination', 'Singles_SVR', 'Original40_SVR'))) +
    theme(legend.position = 'bottom', strip.text = element_text(size = 8)) +
    labs(title = paste('Measured counts:', anchor_name, '+ lapatinib'),
      subtitle = 'Every reserved dose and biological experiment; primary3h models. Anchor-alone curve is the matched single control.',
      x = 'Lapatinib concentration (µM)', y = 'Measured / predicted cell count (% of vehicle)', color = 'Anchor (µM)',
      caption = 'An extra-over-Bliss value does not necessarily imply growth arrest or killing. Bars: +/-1 technical-well SD, descriptive, not biological confidence intervals.')
  print(p)
  ggsave(file.path(figure_dir, paste0('lapatinib_failure_trace_', anchor_name, '.png')),
    p, width = 11.5, height = 9, dpi = 130)
}
dev.off()

# Post-primary input ablation is displayed separately. It predicts actual
# combination counts, avoiding algebraic reuse of single effects in an excess score.
ablation_path <- file.path(data_dir, 'ablation_predictions.tsv.gz')
if (file.exists(ablation_path)) {
  ab <- fread(ablation_path)[split == 'held_partner']
  stopifnot(all(is.finite(ab$predicted_relative_AB)), all(ab$predicted_relative_AB > 0))
  am <- ab[, .(observed_percent = 100 * mean(relative_AB), predicted_percent = 100 * mean(predicted_relative_AB),
    n_bio = .N, observed_min_percent = 100 * min(relative_AB), observed_max_percent = 100 * max(relative_AB),
    observed_GR = mean(GR_AB), predicted_GR = mean(predicted_GR_AB)),
    by = .(condition_id, anchor, partner, time_h, model)]
  save_table(am, 'input_ablation_actual_count_predictions')
  ap <- am[, .(dose_contexts = .N, Pearson_r = cor(observed_percent, predicted_percent),
    mean_observed_percent = mean(observed_percent), mean_predicted_percent = mean(predicted_percent),
    predicted_within_bio_range = sum(predicted_percent >= observed_min_percent & predicted_percent <= observed_max_percent),
    observed_net_loss = sum(observed_GR < 0), predicted_net_loss = sum(predicted_GR < 0),
    true_net_loss_predicted = sum(observed_GR < 0 & predicted_GR < 0),
    false_predicted_net_loss = sum(observed_GR >= 0 & predicted_GR < 0)), by = .(time_h, model, anchor, partner)]
  save_table(ap, 'input_ablation_interpretable_results')
  pdf(file.path(figure_dir, 'input_ablation_observed_vs_predicted_counts.pdf'), width = 13.5, height = 8.5)
  for (t in c(3, 24)) {
    z <- am[time_h == t]
    lab <- z[, .(label = sprintf('r = %.2f; n = %d contexts', cor(observed_percent, predicted_percent), .N)), by = model]
    p <- ggplot(z, aes(observed_percent, predicted_percent, color = partner, shape = anchor)) +
      geom_point(size = 1, alpha = .6) +
      geom_smooth(aes(x = observed_percent, y = predicted_percent, group = 1), inherit.aes = FALSE,
        method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
      geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
                hjust = -.03, vjust = 1.2, size = 2.7) +
      facet_wrap(~model, ncol = 3, scales = 'free') + theme(aspect.ratio = 1, legend.position = 'bottom') +
      labs(title = paste('Can individual RNA predict actual combined cell counts?', t, 'h inputs'),
        subtitle = 'Post-primary input ablation; identical excluded partners. RNA-only models receive no individual functional measurements.',
        x = 'Observed combined cell count (% of vehicle)', y = 'Predicted combined cell count (% of vehicle)',
        color = 'Held partner', shape = 'PI3K anchor',
        caption = 'ContextOnly uses context/dose; Singles_logcounts uses measured individual function. No candidate is selected from these results.')
    print(p)
    ggsave(file.path(figure_dir, paste0('input_ablation_counts_', t, 'h.png')), p, width = 13.5, height = 8.5, dpi = 140)
  }
  dev.off()

  # Decision figure: the controls, full study coverage and an observed failure
  # are shown together. It is deliberately not a positive-only grant panel.
  primary_counts <- pm[time_h == 3 & model %in% c('Independent_singles', 'Singles_SVR', 'Original40_SVR'),
    .(condition_id, anchor, partner, model = as.character(model),
      observed_percent = 100 * observed_relative, predicted_percent = 100 * predicted_relative)]
  ablation_counts <- am[time_h == 3 & model %in% c('ContextOnly_SVR', 'Original40_RNAonly_SVR', 'Singles_logcounts_SVR')]
  count_scatter <- function(tab, levels, labels, title, subtitle) {
    tab <- copy(tab)
    tab[, partner := factor(partner, levels = partners)]
    tab[, display := factor(model, levels = levels, labels = labels)]
    labels_dt <- tab[, .(label = sprintf('r = %.2f; %d contexts', cor(observed_percent, predicted_percent), .N)), by = display]
    ggplot(tab, aes(observed_percent, predicted_percent, color = partner, shape = anchor)) +
      geom_point(size = 1.2, alpha = .6) +
      geom_smooth(aes(x = observed_percent, y = predicted_percent, group = 1), inherit.aes = FALSE,
        method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
      geom_text(data = labels_dt, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -.02, vjust = 1.2, size = 2.8) +
      facet_wrap(~display, nrow = 1, scales = 'free') + theme(aspect.ratio = 1, legend.position = 'bottom') +
      labs(title = title, subtitle = subtitle, x = 'Observed combined cells (% of vehicle)',
        y = 'Predicted combined cells (% of vehicle)', color = 'Excluded partner', shape = 'Anchor')
  }
  p_primary <- count_scatter(primary_counts,
    c('Independent_singles', 'Singles_SVR', 'Original40_SVR'),
    c('Independent single effects', 'Learned from single function', 'Single function + original40 RNA'),
    'A  Primary test: does the original program bank improve combination prediction?',
    'Whole partner excluded under both PI3K anchors; all individual measurements available; 3 h RNA primary')
  p_ablation <- count_scatter(ablation_counts,
    c('ContextOnly_SVR', 'Original40_RNAonly_SVR', 'Singles_logcounts_SVR'),
    c('Context and dose only', 'Context + original40 RNA', 'Context + single function'),
    'B  Direct-count input ablation: RNA does not outperform context and dose',
    'Same exclusions; direct combined-count target; RNA model receives no individual functional measurements')
  decision_page1 <- (p_primary / p_ablation) + plot_annotation(
    title = 'Matched individual responses provide a test of composition, not yet evidence of added molecular prediction',
    subtitle = 'BT20: 2 PI3K anchors x 6 partners; 198 RNA-supported dose contexts; 496 biological observations; 1,488 combination wells',
    caption = 'Every point averages available biological experiments (100 contexts have 3; 98 have 2). Lines are fitted trends, not identity lines.\n3 h RNA-only r = 0.815 versus context-only 0.880; 24 h RNA-only r = 0.828 versus the same 0.880 control. All models/times appear in the accompanying full figures.\nThe primary target was extra inhibition; panels show actual combined cells to avoid attributing a shared single-drug reference to RNA.')
  supported <- functional[condition_id %in% pm$condition_id]
  stopifnot(nrow(supported) == 198L, sum(supported$n_bio) == 496L,
    sum(supported$n_bio == 2L) == 98L, sum(supported$n_bio == 3L) == 100L)
  supported_summary <- supported[, .(contexts = .N, clear_extra = sum(clear_extra_effect),
    new_killing = sum(single_to_combination_kill_repeats == n_bio)), by = .(anchor, partner)]
  save_table(supported_summary, 'supported_pair_biological_summary')
  bs <- melt(supported_summary, id.vars = c('anchor', 'partner', 'contexts'),
    measure.vars = c('clear_extra', 'new_killing'), variable.name = 'criterion', value.name = 'n')
  bs[, criterion := factor(criterion, levels = c('clear_extra', 'new_killing'),
    labels = c('Extra effect >=10 pp; same sign in every repeat', 'Combination-only net cell loss in every repeat'))]
  p_coverage <- ggplot(bs, aes(partner, n, fill = anchor)) +
    geom_col(position = position_dodge(width = .8), width = .7) +
    geom_text(aes(label = paste0(n, '/', contexts)), position = position_dodge(width = .8), vjust = -.3, size = 2.7) +
    facet_wrap(~criterion, ncol = 1, scales = 'free_y') +
    scale_y_continuous(expand = expansion(mult = c(0, .23))) + theme(legend.position = 'bottom') +
    labs(title = 'C  Repeated interactions are present, but extra inhibition and killing are different outcomes',
      x = NULL, y = 'Qualifying / measured contexts', fill = 'PI3K anchor')
  trace_low <- lap_counts[anchor == 'dactolisib' & dose_A_uM < .05]
  trace_low[, series := factor(series, levels = c('Anchor alone', 'Observed combination', 'Singles_SVR', 'Original40_SVR'),
    labels = c('Dactolisib alone', 'Measured combination', 'Single-function prediction', 'Original40 + function prediction'))]
  p_failure <- ggplot(trace_low, aes(dose_B_uM, percent_vehicle, color = series, group = series)) +
    geom_line(linewidth = .7) + geom_point(size = 1.6) + scale_x_log10() +
    geom_errorbar(data = trace_low[is.finite(technical_SD_percent)],
      aes(ymin = percent_vehicle - technical_SD_percent, ymax = percent_vehicle + technical_SD_percent), width = .035) +
    facet_wrap(~biological_replicate, nrow = 1, labeller = label_both) +
    theme(legend.position = 'bottom') +
    labs(title = 'D  Missed anchor-specific increase in cell counts after adding lapatinib',
      subtitle = 'Dactolisib ~0.04 µM; all lapatinib doses and repeats; bars: +/-1 technical SD; primary 3 h predictions',
      x = 'Lapatinib concentration (µM)', y = 'Cells (% of vehicle)', color = NULL)
  p_genes <- ggplot(anchor_genes[gene %in% c('ERBB3', 'DUSP4') & time_h == 24],
    aes(dose_uM, displacement, color = drug, group = drug)) + geom_hline(yintercept = 0, color = 'grey60') +
    geom_line(linewidth = .7) + geom_point(size = 1.6) + scale_x_log10() +
    facet_wrap(~gene, nrow = 1, scales = 'free_y') + theme(legend.position = 'bottom') +
    labs(title = 'E  The anchors do have distinguishable measured single-drug RNA responses',
      subtitle = '24 h transcript displacements; full genes, both times and technical uncertainty in anchor_gene_time_courses.pdf',
      x = 'Anchor concentration (µM)', y = 'Log2-fluorescence displacement', color = 'Anchor')
  decision_page2 <- (p_coverage / p_failure / p_genes) + plot_layout(heights = c(1.45, 1, 1)) + plot_annotation(
    title = 'Biological diagnosis: repeated joint effects, missed anchor-specific increase in counts',
    caption = 'Only the 198 exactly RNA-matched contexts are counted above. Every available biological repeat must satisfy the stated criterion.\nAt 10 µM lapatinib, combination-minus-dactolisib counts are +20.8, +1.7 and +27.6 percentage points; the second repeat has a small difference.\nThis error-selected trace establishes measured behavior, not causal protection. RNA readouts are transcripts, not phosphorylation. Initial densities differ by experiment.')
  pdf(file.path(figure_dir, 'grant_evidence_summary.pdf'), width = 11.5, height = 12)
  print(decision_page1); print(decision_page2); dev.off()
  ggsave(file.path(figure_dir, 'grant_evidence_summary.png'), decision_page1, width = 11.5, height = 12, dpi = 140)
  ggsave(file.path(figure_dir, 'grant_evidence_summary_biology.png'), decision_page2, width = 11.5, height = 12, dpi = 140)

  writeLines(c(
    'Biological interpretation: matched individual responses and combination transfer (workflow 61)',
    '',
    'The scientific test',
    'We used BT20 from Niepel et al. 2017 (doi:10.1038/s41467-017-01383-w): alpelisib or dactolisib combined with lapatinib, neratinib, selumetinib, trametinib, saracatinib or NVP-TAE684. All 12 pairs and full measured dose ranges are retained. A partner is withheld under both anchors at once; its individual molecular and functional measurements remain available. This is retrospective partner exclusion within one cell line, not new-cell transfer, HGF rescue prediction or blind replication of unknown combinations.',
    'Primary molecular input is 3 h, with 24 h reported separately. Exactly matched individual RNA supports 198 dose contexts and 496 biological observations (100 contexts have 3 repeats; 98 have 2), from 1,488 technical combination wells. Each biological experiment is summarized before analysis. Author-provided initial and vehicle counts are preserved.',
    '',
    'What biology is actually measured',
    'Individual-drug RNA contains the expected broad receptor-feedback and ERK-response distinctions. ERBB3 rises with both PI3K anchors; the dactolisib response appears at lower tested doses and remains high at 24 h. DUSP4 is suppressed by alpelisib at 3 h and by several ErbB/MEK partners, whereas dactolisib produces positive DUSP4 transcript displacement. These are measured RNA differences, not inferred phosphoprotein activity. The low-dose anchors are therefore not molecularly indistinguishable.',
    'Named readouts were specified from the primary biological study before model outcomes: receptor feedback, ERK-response genes, cell-cycle and stress genes. Only 13 are measured by the landmark panel; missing genes remain missing. FOXO nuclear localization, receptor phosphorylation and causal bypass were measured in separate experiments in the published study, not in our individual RNA table. We do not infer those measurements from FOXO, ERBB3 or DUSP transcript abundance.',
    'The exact historical original40 bank is used, with signed learned weights and all 40 programs. It is not the later hard-cluster bank. L1000 provides a partial weighted projection with approximately 10% of each program absolute weight observed at the median. External full-versus-landmark agreement and the H1792 ERK trace are separate representation audits, not evidence that these combination outcomes are predicted.',
    '',
    'What the functional data establish',
    sprintf('In the 198 RNA-supported contexts, %d have mean extra inhibition >=10 percentage points and positive extra inhibition in every available biological repeat. %d show combination-only net cell loss in every repeat: neither single causes net cell loss.',
      sum(supported$clear_extra_effect), sum(supported$single_to_combination_kill_repeats == supported$n_bio)),
    'The interaction label is not equivalent to killing. At alpelisib approximately 0.04 µM plus lapatinib 3.33 µM, large extra inhibition coexists with mean GR above 1; individual normalized counts above vehicle inflate the Bliss reference. All such observations are retained and shown as direct cell counts and GR.',
    'An anchor-specific counterexample is dactolisib approximately 0.04 µM plus lapatinib 10 µM: combination cell counts exceed dactolisib alone in all 3 biological experiments. Dactolisib-alone counts are 42.3%, 42.4%, 58.8% of vehicle; the combinations are 63.1%, 44.1%, 86.4%. Differences are +20.8, +1.7 and +27.6 percentage points: marked in two repeats, small in the second. Figures show technical-well SD, not biological confidence intervals. Both tested primary models miss the increase; the original40 model predicts net loss rather than the observed continued growth. This establishes an observed functional difference, not a causal protection mechanism.',
    'Biological experiment 3 starts at much lower cell density and has larger typical excess effects. We preserve it, display its individual measurements and do not treat technical repeats or dose contexts as independent cell lines. Different density is a documented assay difference, not a demonstrated cause of every discrepancy.',
    '',
    'What the model test establishes',
    'The primary original40 model does not improve on the model supplied with the same measured single-agent functional responses. The 3 h extra-effect correlations across 198 context means are 0.714 for original40 plus function and 0.810 for single function alone; these are not molecular-only results.',
    'The direct-count ablation is the cleaner attribution test: RNA-only original40 predicts combined counts with r=0.815 at 3 h and 0.828 at 24 h, but context and dose alone give r=0.880. Measured single function gives r=0.936. Biology and gene-PC RNA views also do not beat that context control. Hence an apparently strong aggregate molecular prediction does not establish added molecular information in this test.',
    'The matched-partner audit asks a sharper biological question: with anchor, both nominal doses and biological experiment fixed, does changing the partner alter the outcome? There are 369 matched contrasts across all 15 partner pairs. Of 154 measured differences at least 10 percentage points, 149 have the same direction in every available biological repeat. Partner-specific functional variation is therefore real. Original40 RNA-only recovers the direction in 104/154 and both the direction and at least 10-point magnitude in 60/154; its contrast correlation is 0.364 at 3 h and 0.217 at 24 h. Measured single function gives 126/154 correct directions and r=0.598. Original40 thus retains partial partner-specific information, substantially weaker than its pooled count correlation suggests. These contrasts share observations and use separate partner-exclusion fits; no independence-based significance claim is made. The context-only contrast is not a meaningful biological partner predictor because its two predictions differ only through separate fitted backgrounds.',
    'Do not use RNA-only extra-effect correlation as evidence of molecular transfer: observed and predicted extra effects share the same measured single-drug product. For dactolisib+lapatinib this gives approximately 0.943 extra-effect correlation while direct combined-count correlation is only approximately 0.114. Actual combined counts and GR expose the failure.',
    'The unconstrained primary extra-effect model can imply impossible negative cell counts: original40_SVR does so in 11 of 496 held biological observations at 3 h (minimum -7.85% of vehicle). These are displayed without silent clipping. The direct-log-count ablation predicts positive counts by construction.',
    '',
    'Implication for the grant',
    'This matched resource establishes that measurable individual responses and repeatable combination phenotypes can be joined in a clean partner-exclusion test. It does not yet establish the central claim that the original40 molecular representation adds transferable information about composition beyond dose/context or measured individual function. Visible RNA biology, positive excess labels and a high pooled correlation cannot substitute for that comparison.',
    'The concrete unresolved step is the decoding of anchor-specific interactions, illustrated by an increase in combination counts despite informative, distinguishable individual RNA. The current results do not prove that the RNA lacks all needed information, that the bank is sufficient, or that a particular bypass causes the failure. Those claims require independent evidence rather than retrospective pathway stories.',
    '',
    'Figure entry points',
    'figures/grant_evidence_summary.pdf: two-page decision figure; controls, coverage and biological failure.',
    'figures/all12_observed_growth_curves.pdf: every pair, mean and each biological experiment.',
    'figures/input_ablation_observed_vs_predicted_counts.pdf: all five input comparisons, both RNA times.',
    'figures/partner_contrasts_3h.pdf and partner_contrasts_24h.pdf: partner-specific cell-count differences after matching anchor, doses and biological experiment.',
    'figures/lapatinib_functional_failure_trace.pdf: every lapatinib dose, both anchors and all repeats.',
    'figures/anchor_gene_time_courses.pdf and named_gene_individual_response_heatmaps.pdf: measured biological readouts.',
    'figures/original40_individual_response_heatmaps.pdf: all original40 programs, both times, identical display order.',
    'All plotted values are provided in data/figure_*.tsv.gz; original source and model outputs remain separate.'
  ), file.path(root, 'BIOLOGICAL_INTERPRETATION.txt'))
}
figure_files <- list.files(figure_dir, pattern = '\\.pdf$', full.names = FALSE)
figure_index <- data.table(file = file.path('figures', figure_files),
  description = c(
    all12_functional_consistency = 'Repeat-consistent extra inhibition, net cell loss and combination-only net cell loss across all measured dose grids.',
    all12_observed_growth_curves = 'All 12 pairs: mean and each biological experiment; GR and measured partner-alone curves.',
    all12_reserved_effect_curves = 'Observed and reserved primary predictions for every pair and both times.',
    anchor_gene_time_courses = 'Measured ERBB3/EGFR, ERK-response and cell-cycle transcripts; both anchors, doses, times; technical uncertainty.',
    functional_control_counts = 'Supplied initial and vehicle cell counts by biological/technical experiment.',
    grant_evidence_summary = 'Two-page grant decision: full coverage, actual cell counts, primary and ablation controls, missed count increase and measured RNA.',
    held_partner_observed_vs_predicted = 'All eight primary models, both times; extra-effect scatter plots.',
    individual_drug_functional_controls = 'Actual single-drug GR by dose and biological experiment.',
    input_ablation_observed_vs_predicted_counts = 'Five input ablations: actual combined counts, both times, identical partner exclusions.',
    lapatinib_functional_failure_trace = 'Direct counts at every lapatinib dose under both anchors; all biological experiments and primary models.',
    named_gene_individual_response_heatmaps = 'All measured prespecified biological transcript readouts across eight drugs, doses and both times.',
    original40_H1792_ERK_gene_to_program_trace = 'Independent H1792 mechanistic gene-to-original40-program decomposition.',
    original40_H1792_replication_all40 = 'Independent H1792 RNA versus Tahoe original40 responses; formulations separate.',
    original40_full_to_landmark_fidelity = 'External full-gene versus observed landmark projections of the exact original40 bank.',
    original40_individual_response_heatmaps = 'All exact historical original40 signed program responses to single drugs; both times, fixed order.'
    ,training_inner_held_diagnostics = 'Direct-count input-ablation capacity at training, inner-partner and held-partner stages: original40 RNA-only, context and single-function models.'
    ,training_inner_held_fold_summary = 'Fold-level training and exclusion summaries for the direct-count input-ablation models.'
    ,partner_contrasts_3h = 'Primary RNA: measured versus predicted changes when the partner varies but anchor, doses and biological experiment are matched.'
    ,partner_contrasts_24h = 'Secondary RNA: measured versus predicted changes when the partner varies but anchor, doses and biological experiment are matched.'
  )[sub('\\.pdf$', '', figure_files)])
figure_index[is.na(description), description := 'Companion audit figure; see the matching source table and workflow report.']
fwrite(figure_index, file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
writeLines(c('Complete observed function, original40 and named-gene molecular figures, and reserved prediction figures generated.',
  'No representation or timepoint was selected using held-partner outcomes. All12 pairs remain visible.'), file.path(root, 'FIGURE_STATUS.txt'))
capture.output(sessionInfo(), file = file.path(data_dir, 'diagnostic_figures_session_info.txt'))
message('Complete combination diagnostic figures generated.')
