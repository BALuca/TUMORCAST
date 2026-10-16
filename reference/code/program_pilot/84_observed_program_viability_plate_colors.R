#!/usr/bin/env Rscript
# A plate-colored view of the existing v84 observed-data scatters. Coordinates,
# Pearson annotations, matched profiles, and the measured RNA reference are frozen.
# Multi-plate source profiles remain one point with their pooled plate group named.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
setDTthreads(1L)

# This is an alternate figure view, so its outputs stay in the existing v84 root.
root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
figure_root <- file.path(root, 'figures/observed_program_viability')
points_path <- file.path(root, 'data/observed_program_viability_points.tsv.gz')
metrics_path <- file.path(root, 'data/observed_program_viability_pooled_correlations.tsv')
script_path <- '/work/src/program_pilot/84_observed_program_viability_plate_colors.R'
programs <- sprintf('P%02d', 1:40)
dose_shapes <- c('0.05' = 16, '0.5' = 17, '5' = 15)
# Blue, orange, and reddish purple from the colorblind-compatible Okabe-Ito palette.
# These three colors identify the three plate groups within each drug's figure.
plate_palette <- c('#0072B2', '#D55E00', '#CC79A7')
inputs <- c(points_path, metrics_path)
input_hashes <- vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')

points <- fread(points_path)
metrics <- fread(metrics_path)
keys <- c('depmap_id', 'drug', 'dose', 'program')
required <- c(keys, 'contributing_plates', 'observed_displacement', 'observed_log2_abundance')
stopifnot(all(required %in% names(points)), !anyDuplicated(points[, ..keys]),
  setequal(points$program, programs), uniqueN(points$drug) == 10L,
  setequal(points$dose, c(.05, .5, 5)),
  all(is.finite(points$observed_displacement)), all(is.finite(points$observed_log2_abundance)),
  !anyNA(points$contributing_plates), !any(points$contributing_plates == ''),
  !anyDuplicated(metrics[, .(drug, program)]),
  nrow(unique(points[, .(depmap_id, drug, dose)])) == 867L)
points[, dose_label := factor(as.character(dose), levels = names(dose_shapes))]
points[, plate_label := ifelse(grepl(',', contributing_plates, fixed = TRUE),
  paste0(gsub(',', '+', contributing_plates, fixed = TRUE), ' (pooled)'), contributing_plates)]

# Verify and export the exact plate/dose mapping used for colors. Sorting doses
# establishes a deterministic within-drug color order, without changing points.
mapping <- unique(points[, .(drug, dose, contributing_plates, plate_label)])
setorder(mapping, drug, dose)
stopifnot(all(mapping[, .N, by = .(drug, dose)]$N == 1L),
  all(mapping[, uniqueN(plate_label), by = drug]$V1 == 3L))
mapping[, color := plate_palette, by = drug]
fwrite(mapping, file.path(root, 'data/observed_program_viability_plate_color_key.tsv'), sep = '\t')

# Reused for each ten-program page. Each colored line is fitted across cells at
# one dose, with the color naming that dose's actual contributing plate group.
make_page <- function(drug_value, page) {
  these_programs <- programs[(page - 1L) * 10L + 1:10]
  rows <- copy(points[drug == drug_value & program %in% these_programs])
  rows[, program := factor(program, levels = these_programs)]
  labels <- copy(metrics[drug == drug_value & program %in% these_programs])
  labels[, program := factor(program, levels = these_programs)]
  labels[, label := sprintf('Pooled r = %.2f\nWithin-dose r = %.2f', pooled_r, within_dose_r)]
  drug_mapping <- mapping[drug == drug_value]
  plate_colors <- setNames(drug_mapping$color, drug_mapping$plate_label)
  rows[, plate_label := factor(plate_label, levels = names(plate_colors))]
  note <- if (drug_value == 'Lapatinib ditosylate') {
    'Lapatinib: dose and plate are completely confounded (0.05 / 0.5 / 5 µM = plates 1 / 2 / 3).'
  } else if (any(grepl(',', drug_mapping$contributing_plates, fixed = TRUE))) {
    'A label such as 6+14 (pooled) names an already-pooled profile; it is not separate plate measurements.'
  } else {
    'Within this drug, each dose comes from a different plate; the two labels cannot isolate their effects.'
  }
  ggplot(rows, aes(observed_displacement, observed_log2_abundance)) +
    geom_point(aes(color = plate_label, shape = dose_label), size = 1.35, alpha = .82) +
    geom_smooth(aes(color = plate_label, group = dose_label), method = 'lm', formula = y ~ x,
      se = FALSE, linewidth = .55) +
    geom_label(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
      hjust = -.02, vjust = 1.03, size = 2.4, label.size = 0, fill = '#FFFFFFCC',
      label.padding = grid::unit(.06, 'lines')) +
    facet_wrap(~program, ncol = 5, scales = 'free_x') +
    scale_color_manual(values = plate_colors, drop = FALSE, name = 'Contributing plate(s)') +
    scale_shape_manual(values = dose_shapes, drop = FALSE, name = 'RNA dose (µM)') +
    scale_y_continuous(expand = expansion(mult = c(.04, .16))) +
    guides(color = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 1, size = 2.4)),
      shape = guide_legend(order = 2, override.aes = list(color = '#333333', alpha = 1, size = 2.4))) +
    theme_bw(base_size = 10) +
    theme(aspect.ratio = 1, legend.position = 'bottom', legend.box = 'horizontal',
      panel.grid.minor = element_blank(), strip.background = element_rect(fill = '#F3F3F3'),
      plot.margin = margin(6, 8, 5, 5)) +
    labs(x = 'Measured program displacement (signed RNA score)',
      y = 'Observed PRISM log2 abundance / vehicle',
      title = paste0(drug_value, ': original v84 scatters with plate colors and dose shapes'),
      subtitle = paste0('All original measured points retained; programs ', these_programs[1], '–',
        these_programs[10], '. Lines compare cells at the same dose.'),
      caption = paste0(note, '\nPRISM is interpolated to the RNA dose; 0 = 100%, −1 = 50%, −2 = 25% of vehicle. ',
        'Within-dose r removes each dose mean from both axes.'))
}

# Save every drug and the combined forty-page book. Focus previews show the two
# programs named by the user in their original program-order pages.
drug_order <- c('Lapatinib ditosylate', 'Busulfan', setdiff(unique(points$drug),
  c('Lapatinib ditosylate', 'Busulfan')))
all_pages <- list()
index <- list()
for (drug_value in drug_order) {
  slug <- gsub('^_|_$', '', gsub('[^a-z0-9]+', '_', tolower(drug_value)))
  output_path <- file.path(figure_root, paste0(slug, '_all_40_programs_plate_colors.pdf'))
  drug_pages <- lapply(1:4, function(page) make_page(drug_value, page))
  # Cairo retains Unicode dose and minus symbols without PDF font substitutions.
  cairo_pdf(output_path, width = 12.5, height = 6.5, onefile = TRUE)
  for (p in drug_pages) print(p)
  dev.off()
  all_pages <- c(all_pages, drug_pages)
  if (drug_value == 'Lapatinib ditosylate') {
    ggsave(file.path(figure_root, 'lapatinib_ditosylate_page_3_plate_colors.png'),
      drug_pages[[3]], width = 12.5, height = 6.5, dpi = 160)
  }
  if (drug_value == 'Busulfan') {
    ggsave(file.path(figure_root, 'busulfan_page_2_plate_colors.png'),
      drug_pages[[2]], width = 12.5, height = 6.5, dpi = 160)
  }
  index[[drug_value]] <- data.table(drug = drug_value, pdf = output_path,
    pages = 4L, profiles = uniqueN(points[drug == drug_value, .(depmap_id, dose)]))
  cat('Saved ', drug_value, ': ', output_path, '\n', sep = '')
}
cairo_pdf(file.path(figure_root, 'all_ten_drugs_all_40_programs_plate_colors.pdf'),
  width = 12.5, height = 6.5, onefile = TRUE)
for (p in all_pages) print(p)
dev.off()
fwrite(rbindlist(index), file.path(root, 'data/observed_program_viability_plate_figure_index.tsv'), sep = '\t')

# Ensure this plot-only operation preserved both original data inputs byte-for-byte.
stopifnot(identical(input_hashes, vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')))
fwrite(data.table(path = c(inputs, script_path),
  sha256 = vapply(c(inputs, script_path), digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/observed_program_viability_plate_input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'data/observed_program_viability_plate_session.txt'))
cat('Complete: 400 scatters; 867 original cell–dose profiles; no filtering, refitting, or coordinate changes.\n')
