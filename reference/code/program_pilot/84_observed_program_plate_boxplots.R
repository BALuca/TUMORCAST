#!/usr/bin/env Rscript
# Boxplot view of the original v84 measured displacement table. Each observation
# is one cell-line profile; archived replicate pools stay labelled as pools.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
setDTthreads(1L)

# Plot variants belong to the existing analysis folder, with no model changes.
root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
figure_root <- file.path(root, 'figures/observed_program_viability')
points_path <- file.path(root, 'data/observed_program_viability_points.tsv.gz')
mapping_path <- file.path(root, 'data/observed_program_viability_plate_color_key.tsv')
script_path <- '/work/src/program_pilot/84_observed_program_plate_boxplots.R'
programs <- sprintf('P%02d', 1:40)
inputs <- c(points_path, mapping_path)
input_hashes <- vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')

points <- fread(points_path)
mapping <- fread(mapping_path)
keys <- c('depmap_id', 'drug', 'dose', 'program')
required <- c(keys, 'contributing_plates', 'observed_displacement')
stopifnot(all(required %in% names(points)), !anyDuplicated(points[, ..keys]),
  setequal(points$program, programs), uniqueN(points$drug) == 10L,
  setequal(points$dose, c(.05, .5, 5)),
  all(is.finite(points$observed_displacement)),
  !anyNA(points$contributing_plates), !any(points$contributing_plates == ''),
  nrow(unique(points[, .(depmap_id, drug, dose)])) == 867L,
  !anyDuplicated(mapping[, .(drug, dose, contributing_plates)]),
  all(c('plate_label', 'color') %in% names(mapping)))
points <- merge(points, mapping[, .(drug, dose, contributing_plates, plate_label, color)],
  by = c('drug', 'dose', 'contributing_plates'), all.x = TRUE, sort = FALSE)
stopifnot(nrow(points) == 34680L, !anyNA(points$plate_label), !anyNA(points$color))

# Export unscaled assay-score distribution summaries for all 1,200 boxes.
summaries <- points[, .(n_cells = .N, mean = mean(observed_displacement),
  median = median(observed_displacement), sd = sd(observed_displacement),
  q25 = unname(quantile(observed_displacement, .25)),
  q75 = unname(quantile(observed_displacement, .75)),
  minimum = min(observed_displacement), maximum = max(observed_displacement)),
  by = .(drug, dose, contributing_plates, plate_label, program)]
setorder(summaries, drug, program, dose)
stopifnot(nrow(summaries) == 1200L, all(summaries$n_cells >= 2L))
fwrite(summaries, file.path(root, 'data/observed_program_plate_boxplot_summary.tsv'), sep = '\t')

# Reused across every drug and program page. A common y scale within each
# program allows the three plate/dose groups to be compared directly.
make_page <- function(drug_value, selected_programs) {
  rows <- copy(points[drug == drug_value & program %in% selected_programs])
  counts <- copy(summaries[drug == drug_value & program %in% selected_programs])
  drug_mapping <- copy(mapping[drug == drug_value])
  setorder(drug_mapping, dose)
  drug_mapping[, x_label := paste0(as.character(dose), '\nplate ', plate_label)]
  plate_colors <- setNames(drug_mapping$color, drug_mapping$plate_label)
  x_labels <- setNames(drug_mapping$x_label, as.character(drug_mapping$dose))
  rows[, dose_group := factor(as.character(dose), levels = c('0.05', '0.5', '5'))]
  counts[, dose_group := factor(as.character(dose), levels = c('0.05', '0.5', '5'))]
  rows[, program := factor(program, levels = selected_programs)]
  counts[, program := factor(program, levels = selected_programs)]
  pooled_note <- if (any(grepl(',', drug_mapping$contributing_plates, fixed = TRUE))) {
    'A+B (pooled) is an archived combined replicate profile, with one observation per cell line.'
  } else {
    'Each dose has one contributing plate; dose and plate effects are confounded.'
  }
  ggplot(rows, aes(dose_group, observed_displacement)) +
    geom_hline(yintercept = 0, color = '#999999', linetype = 'dashed', linewidth = .3) +
    geom_boxplot(aes(fill = plate_label), width = .55, alpha = .5,
      outlier.shape = NA, linewidth = .45) +
    geom_point(aes(color = plate_label),
      position = position_jitter(width = .12, height = 0, seed = 84), size = .9, alpha = .7) +
    geom_point(data = counts, aes(y = mean), shape = 18, size = 2, color = '#111111') +
    geom_text(data = counts, aes(y = Inf, label = paste0('n=', n_cells)),
      vjust = 1.3, size = 2.5, color = '#333333') +
    facet_wrap(~program, ncol = min(5L, length(selected_programs)), scales = 'free_y') +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(.05, .13))) +
    scale_fill_manual(values = plate_colors, guide = 'none') +
    scale_color_manual(values = plate_colors, guide = 'none') +
    theme_bw(base_size = 10) +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(), strip.background = element_rect(fill = '#F3F3F3'),
      axis.text.x = element_text(size = 8, lineheight = .95),
      plot.margin = margin(6, 8, 5, 5)) +
    labs(title = paste0(drug_value, ': ', paste(range(selected_programs), collapse = '–')),
      subtitle = 'Original v84 observed molecular displacements across cell lines; each box is one plate/dose group.',
      x = 'Drug concentration (µM) and contributing plate(s)',
      y = 'Measured program displacement (signed RNA score)',
      caption = paste0('Boxes: median and interquartile range; whiskers: 1.5×IQR; points: cell lines; diamonds: means.\n',
        'Original matched-plate control, no rescaling or filtering. ', pooled_note))
}

# Produce the same four ten-program pages per drug as the earlier scatter view.
drug_order <- c('Lapatinib ditosylate', 'Busulfan', setdiff(unique(points$drug),
  c('Lapatinib ditosylate', 'Busulfan')))
all_pages <- list()
index <- list()
for (drug_value in drug_order) {
  slug <- gsub('^_|_$', '', gsub('[^a-z0-9]+', '_', tolower(drug_value)))
  output_path <- file.path(figure_root, paste0(slug, '_all_40_programs_plate_boxplots.pdf'))
  drug_pages <- lapply(1:4, function(page) make_page(drug_value, programs[(page - 1L) * 10L + 1:10]))
  cairo_pdf(output_path, width = 12.5, height = 7, onefile = TRUE)
  for (p in drug_pages) print(p)
  dev.off()
  all_pages <- c(all_pages, drug_pages)
  index[[drug_value]] <- data.table(drug = drug_value, pdf = output_path, pages = 4L,
    profiles = uniqueN(points[drug == drug_value, .(depmap_id, dose)]))
  cat('Saved ', drug_value, '\n', sep = '')
}
cairo_pdf(file.path(figure_root, 'all_ten_drugs_all_40_programs_plate_boxplots.pdf'),
  width = 12.5, height = 7, onefile = TRUE)
for (p in all_pages) print(p)
dev.off()

# Show the two programs highlighted by the user immediately in a compact view.
focus <- wrap_plots(make_page('Lapatinib ditosylate', 'P21'), make_page('Busulfan', 'P12'), ncol = 1)
ggsave(file.path(figure_root, 'lapatinib_P21_busulfan_P12_plate_boxplots.png'),
  focus, width = 8, height = 11, dpi = 160)
ggsave(file.path(figure_root, 'lapatinib_P21_busulfan_P12_plate_boxplots.pdf'),
  focus, width = 8, height = 11, device = cairo_pdf)
fwrite(rbindlist(index), file.path(root, 'data/observed_program_plate_boxplot_figure_index.tsv'), sep = '\t')

# Verify plot generation preserved the frozen input tables.
stopifnot(identical(input_hashes, vapply(inputs, digest, character(1), file = TRUE, algo = 'sha256')))
fwrite(data.table(path = c(inputs, script_path),
  sha256 = vapply(c(inputs, script_path), digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/observed_program_plate_boxplot_input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'data/observed_program_plate_boxplot_session.txt'))
cat('Complete: 400 program panels, 1,200 plate/dose boxes, 867 original cell–dose profiles.\n')
