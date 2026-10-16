#!/usr/bin/env Rscript
# Descriptive audit of the exact matrix shown in scripts 35/36; no model fitting.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2)})
source('/work/src/program_pilot/response_models.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- '/work/results/program_pilot/40_constrained_program_sensitivity_classifier'
input_path <- '/work/results/program_pilot/35_all_training_program_displacement_heatmaps/data/ERBB_all_training_inputs.qs'
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout'
input <- qread(input_path)
md <- copy(input$annotations)
d <- input$displacement
programs <- colnames(d)
stopifnot(nrow(d) == 6260L, ncol(d) == 40L, identical(rownames(d), md$observation_id),
  uniqueN(md$depmap_id) == 26L, all(is.finite(d)), all(md$measured_pct > 0))
z <- as.data.table(d)
z[, `:=`(condition = md$condition, depmap_id = md$depmap_id, recipient = md$recipient,
  log2_abundance = log2(md$measured_pct / 100), abundance_pct = md$measured_pct)]
tables <- list()
for (program in programs) {
  x <- z[[program]]
  y <- z$log2_abundance
  within_context_x <- x - ave(x, z$condition)
  within_context_y <- y - ave(y, z$condition)
  within_cell_x <- x - ave(x, z$depmap_id)
  within_cell_y <- y - ave(y, z$depmap_id)
  context <- z[, .(displacement = mean(get(program)), log2_abundance = mean(log2_abundance)), by = condition]
  by_cell <- z[, .(program = program, observations = .N,
    log2_abundance_r = cor_safe(get(program), log2_abundance),
    abundance_pct_r = cor_safe(get(program), abundance_pct)), by = recipient]
  tables[[program]] <- list(cell = by_cell, axes = data.table(program = program,
    all_observations_r = cor_safe(x, y), across_treatment_mean_r = cor_safe(context$displacement, context$log2_abundance),
    within_treatment_across_cells_r = cor_safe(within_context_x, within_context_y),
    within_cell_across_treatments_r = cor_safe(within_cell_x, within_cell_y),
    median_cell_r = median(by_cell$log2_abundance_r),
    positive_cells = sum(by_cell$log2_abundance_r > 0), negative_cells = sum(by_cell$log2_abundance_r < 0)))
}
axes <- rbindlist(lapply(tables, `[[`, 'axes'))
by_cell <- rbindlist(lapply(tables, `[[`, 'cell'))
fwrite(axes, file.path(root, 'data', 'program_relationship_axes.tsv'), sep = '\t')
fwrite(by_cell, file.path(root, 'data', 'program_relationships_by_cell.tsv'), sep = '\t')

# The total-variance partition is descriptive and uses sample-weighted context means.
# It is not held-out prediction, replicate reliability or a causal decomposition.
partition <- data.table(outcome_scale = c('Abundance (% vehicle)', 'Log2 abundance'),
  variation_between_drug_dose_means_pct = c(
    100 * var(ave(z$abundance_pct, z$condition)) / var(z$abundance_pct),
    100 * var(ave(z$log2_abundance, z$condition)) / var(z$log2_abundance)))
fwrite(partition, file.path(root, 'data', 'training_outcome_variation_partition.tsv'), sep = '\t')

# Show EVERY program and EVERY training cell. The clustering organizes this
# descriptive association matrix, without changing any frozen gene program.
cells <- sort(unique(by_cell$recipient))
mat <- matrix(NA_real_, length(programs), length(cells), dimnames = list(programs, cells))
mat[cbind(match(by_cell$program, programs), match(by_cell$recipient, cells))] <- by_cell$log2_abundance_r
stopifnot(all(is.finite(mat)))
ro <- hclust(dist(mat), method = 'ward.D2')$order
co <- hclust(dist(t(mat)), method = 'ward.D2')$order
qsave(list(correlations = mat, row_order = ro, column_order = co),
  file.path(root, 'data', 'program_relationships_by_cell.qs'))
for (format in c('pdf', 'png')) {
  path <- file.path(root, 'figures', paste0('program_viability_relationships_by_cell.', format))
  if (format == 'pdf') grDevices::cairo_pdf(path, width = 8.5, height = 7.5) else {
    png(path, width = 8.5, height = 7.5, units = 'in', res = 180, type = 'cairo')
  }
  h <- draw_heatmap(mat[ro, co, drop = FALSE], scale_rows = FALSE,
    value_range = c(-1, 0, 1), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = 'Pearson r\n(log2 abundance)', cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8), column_names_rot = 60,
    title = 'Each program versus viability ACROSS TREATMENTS within each training cell',
    column_title_gp = grid::gpar(fontsize = 10))
  ComplexHeatmap::draw(h, heatmap_legend_side = 'right')
  dev.off()
}

# P14 is shown explicitly because the preceding loss audit already identified it
# as the strongest S/R mean difference; it is not selected using held outcomes.
plot_data <- data.table(displacement = d[, 'P14'], abundance_pct = md$measured_pct,
  condition = md$condition, dose = factor(md$dose))
centered <- copy(plot_data)
centered[, `:=`(displacement = displacement - mean(displacement),
  abundance_pct = abundance_pct - mean(abundance_pct)), by = condition]
plot_data[, axis := 'All treatments combined']
centered[, axis := 'Cells at the same drug and dose: centered']
plot_data <- rbindlist(list(plot_data, centered))
stats <- plot_data[, .(label = sprintf('Pearson r = %.3f', cor_safe(displacement, abundance_pct))), by = axis]
fwrite(stats, file.path(root, 'data', 'P14_relationship_plot_statistics.tsv'), sep = '\t')
p <- ggplot(plot_data, aes(displacement, abundance_pct, color = dose)) +
  geom_point(alpha = .15, size = .5) +
  geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = stats, aes(x = -Inf, y = Inf, label = label), color = 'black',
    inherit.aes = FALSE, hjust = -.03, vjust = 1.2, size = 3) +
  facet_wrap(~axis, scales = 'free', nrow = 1) + theme_classic(base_size = 10) +
  theme(aspect.ratio = 1, legend.position = 'bottom') +
  labs(x = 'P14 expression displacement (centered in right panel)',
    y = 'Abundance (% vehicle; centered in right panel)', color = 'Dose (uM)',
    title = 'The combined sensitivity pattern and recipient differences are different signals',
    caption = 'Exact 6,260 training observations. Centering removes each drug/dose mean from both axes; all observations retained.')
ggsave(file.path(root, 'figures', 'P14_relationship_comparison.pdf'), p,
  device = grDevices::cairo_pdf, width = 8, height = 4.5)
ggsave(file.path(root, 'figures', 'P14_relationship_comparison.png'), p, width = 8, height = 4.5, dpi = 180)

# Preserve the previously performed baseline-source audit as evidence; do not
# invent a new baseline-conditioning fix or repeat models already tried.
baseline_paths <- file.path(parent, 'data', c('contextual_audit_baseline_geometry.tsv',
  'contextual_audit_baseline_source_agreement.tsv', 'contextual_audit_baseline_gene_source_agreement.tsv'))
sources <- c(input_path, baseline_paths, '/work/src/program_pilot/40_program_relationship_audit.R',
  '/work/src/program_pilot/response_models.R', '/work/src/R/lib_heatmap.R')
fwrite(data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data', 'relationship_audit_input_hashes.tsv'), sep = '\t')
writeLines(capture.output(sessionInfo()), file.path(root, 'relationship_audit_sessionInfo.txt'))
print(axes[program %in% c('P14', 'P26', 'P31', 'P40')])
print(partition)
message('Relationship audit COMPLETE: no models fitted; outputs kept in existing analysis folder.')
