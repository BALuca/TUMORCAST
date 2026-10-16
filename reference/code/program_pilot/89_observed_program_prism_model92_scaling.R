#!/usr/bin/env Rscript
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
# Shared v92 per-program scaling, plotted across the original nine v89 drug/dose
# contexts used in v89. This is a plot variant, with no fitting changes.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
setDTthreads(1L)
root <- '/work/results/program_pilot/89_shared_program_interaction_decoder_mock'
figure_root <- file.path(root, 'figures/observed_programs_across_drugs')
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
points_path <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge/data/observed_program_viability_points.tsv.gz'
selection_path <- file.path(root, 'data/selected_doses.tsv')
script_path <- '/work/src/program_pilot/89_observed_program_prism_model92_scaling.R'
programs <- sprintf('P%02d', 1:40)
inputs <- c(points_path, selection_path)
hashes <- vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')
drug_names <- c('5-Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan',
  'Cytarabine (hydrochloride)', 'Lapatinib ditosylate', 'Topotecan (hydrochloride)',
  'Vandetanib', 'Vismodegib')
drug_labels <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
  'Lapatinib', 'Topotecan', 'Vandetanib', 'Vismodegib')
drug_colors <- setNames(c('#4E79A7', '#F28E2B', '#59A14F', '#E15759', '#76B7B2',
  '#B07AA1', '#9C755F', '#B89800', '#79706E'), drug_names)
pearson <- function(x, y) if (length(x) >= 3L && sd(x) > 1e-12 && sd(y) > 1e-12) cor(x, y) else NA_real_

selected <- fread(selection_path)
points <- merge(fread(points_path), selected[, .(drug, dose)], by = c('drug', 'dose'))
keys <- c('depmap_id', 'drug', 'dose', 'program')
stopifnot(setequal(points$drug, drug_names), setequal(points$program, programs),
  !anyDuplicated(points[, ..keys]), nrow(points) == sum(selected$cells) * 40L,
  all(points[, .N, by = program]$N == 256L),
  all(is.finite(points$observed_displacement)), all(is.finite(points$observed_log2_abundance)),
  !'Vemurafenib' %in% points$drug)
# Reconstruct the final v92 shared per-program scale from its training RNA.
# All ten model contexts (including crizotinib 5 uM) contribute to this scale;
# the plotted profiles remain exactly the nine contexts in the requested PDF.
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
model_selection <- fread('/work/results/program_pilot/92_shared_molecular_pattern_sensitivity_decoder/data/selected_doses.tsv')
m <- merge(as.data.table(archive$prep$m), model_selection, by = c('drug', 'dose'))[depmap_id != 'ACH-000248']
D <- as.matrix(m[, ..programs])
residual <- D
for (dr in unique(m$drug)) {
  ii <- which(m$drug == dr)
  residual[ii, ] <- sweep(D[ii, , drop = FALSE], 2L, colMeans(D[ii, , drop = FALSE]))
}
shared_scale <- sqrt(colMeans(residual^2))
stopifnot(all(is.finite(shared_scale)), all(shared_scale > 0))
scale_table <- data.table(program = programs, shared_within_treatment_RMS = shared_scale)
fwrite(scale_table, file.path(root, 'data/model92_shared_program_scales.tsv'), sep = '	')
points[, model_scale := shared_scale[match(program, programs)]]
points[, scaled_displacement := observed_displacement / model_scale]
# Scaling is before the six-axis projection; no drug/dose mean is subtracted
# from the plotted input, and PRISM labels retain their original assay scale.
points[, drug := factor(drug, levels = drug_names)]
legend_rows <- selected[match(drug_names, drug)]
legend_labels <- setNames(paste0(drug_labels, ' | ', legend_rows$dose, ' µM; n=', legend_rows$cells), drug_names)

# These residuals are used only for the descriptive correlation annotation.
# All plotted coordinates remain the original uncentered physical values.
points[, x_within_drug := observed_displacement - mean(observed_displacement), by = .(drug, program)]
points[, y_within_drug := observed_log2_abundance - mean(observed_log2_abundance), by = .(drug, program)]
metrics <- points[, .(profiles = .N,
  pooled_r = pearson(observed_displacement, observed_log2_abundance),
  within_drug_r = pearson(x_within_drug, y_within_drug)), by = program]
by_drug <- points[, .(cells = .N, r = pearson(observed_displacement, observed_log2_abundance)),
  by = .(drug, dose, program)]
fwrite(metrics, file.path(root, 'data/observed_program_across_drugs_model92_scaling_correlations.tsv'), sep = '\t')
fwrite(by_drug, file.path(root, 'data/observed_program_across_drugs_model92_scaling_each_drug_correlations.tsv'), sep = '\t')
fwrite(points[, .(depmap_id, drug, dose, program, contributing_plates,
  observed_displacement, model_scale, scaled_displacement, observed_log2_abundance)],
  file.path(root, 'data/observed_program_across_drugs_model92_scaling_plot_points.tsv.gz'), sep = '\t')

make_page <- function(selected_programs) {
  rows <- copy(points[program %in% selected_programs])
  rows[, program := factor(program, levels = selected_programs)]
  labels <- copy(metrics[program %in% selected_programs])
  labels[, program := factor(program, levels = selected_programs)]
  labels[, label := sprintf('Pooled r = %.2f\nWithin-drug r = %.2f', pooled_r, within_drug_r)]
  ggplot(rows, aes(scaled_displacement, observed_log2_abundance, color = drug)) +
    geom_point(size = 1, alpha = .75) +
    geom_smooth(aes(group = drug), method = 'lm', formula = y ~ x, se = FALSE, linewidth = .55) +
    geom_label(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.02, vjust = 1.03, size = 2.5, label.size = 0, fill = '#FFFFFFD9',
      label.padding = grid::unit(.06, 'lines')) +
    facet_wrap(~program, ncol = min(5L, length(selected_programs)), scales = 'free_x') +
    scale_color_manual(values = drug_colors, labels = legend_labels, drop = FALSE,
      name = 'Drug | included molecular dose; cell lines') +
    scale_y_continuous(expand = expansion(mult = c(.04, .14))) +
    guides(color = guide_legend(ncol = 3, byrow = TRUE,
      override.aes = list(size = 2.5, alpha = 1))) +
    theme_bw(base_size = 10) +
    theme(aspect.ratio = 1, legend.position = 'bottom', panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = '#F3F3F3'),
      plot.margin = margin(6, 8, 5, 5)) +
    labs(title = 'Model-scaled program displacement versus PRISM abundance: same nine contexts',
      subtitle = 'All nine contexts shown together. Each colored line compares cells receiving the same drug and dose.',
      x = 'Program displacement / shared within-treatment RMS (v92 scale)',
      y = 'Observed PRISM log2 abundance / vehicle',
      caption = paste0('Raw displacement divided by one shared v92 scale per program; drug means retained; PRISM unchanged. 256 profiles / 30 cells.\n',
        'Pooled r mixes drug differences and cell variation; within-drug r removes each drug/dose mean from both variables. ',
        'Abundance: 0 = 100%, −1 = 50%, −2 = 25% of vehicle.'))
}

cairo_pdf(file.path(figure_root, 'all_40_programs_across_included_drugs_model92_scaling.pdf'),
  width = 12.5, height = 7.8, onefile = TRUE)
for (page in 1:4) {
  these_programs <- programs[(page - 1L) * 10L + 1:10]
  p <- make_page(these_programs)
  print(p)
  ggsave(file.path(figure_root, paste0('programs_page_', page, '_across_included_drugs_model92_scaling.png')),
    p, width = 12.5, height = 7.8, dpi = 160)
}
dev.off()
stopifnot(identical(hashes, vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')))
fwrite(data.table(path = c(inputs, script_path),
  sha256 = vapply(c(inputs, script_path), digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/observed_program_across_drugs_model92_scaling_input_hashes.tsv'), sep = '\t')
cat('Complete: 40 program scatters, nine selected drug/dose contexts, 256 measured profiles; shared v92 scaling; original inputs unchanged.\n')
