#!/usr/bin/env Rscript
# Plate is the organizing unit: drugs sharing a plate are adjacent boxes.
# Original v84 replicate pools have separate x positions and are never unpooled.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
setDTthreads(1L)

root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
figure_root <- file.path(root, 'figures/observed_program_viability')
points_path <- file.path(root, 'data/observed_program_viability_points.tsv.gz')
script_path <- '/work/src/program_pilot/84_observed_program_boxplots_by_plate.R'
programs <- sprintf('P%02d', 1:40)
plate_order <- c('1', '2', '3', '4', '5', '6,14', '7', '7,13', '8', '8,13',
  '9', '9,13', '10', '11', '12')
drug_names <- c('5-Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan',
  'Cytarabine (hydrochloride)', 'Lapatinib ditosylate', 'Topotecan (hydrochloride)',
  'Vandetanib', 'Vemurafenib', 'Vismodegib')
drug_labels <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
  'Lapatinib', 'Topotecan', 'Vandetanib', 'Vemurafenib', 'Vismodegib')
drug_colors <- setNames(c('#4E79A7', '#F28E2B', '#59A14F', '#E15759', '#76B7B2',
  '#B07AA1', '#9C755F', '#EDC949', '#FF9DA7', '#79706E'), drug_names)
input_hash <- digest(points_path, file = TRUE, algo = 'sha256')

points <- fread(points_path)
keys <- c('depmap_id', 'drug', 'dose', 'program')
stopifnot(all(c(keys, 'contributing_plates', 'observed_displacement') %in% names(points)),
  !anyDuplicated(points[, ..keys]), nrow(points) == 34680L,
  setequal(points$program, programs), setequal(points$drug, drug_names),
  setequal(points$contributing_plates, plate_order),
  all(is.finite(points$observed_displacement)))
mapping <- unique(points[, .(contributing_plates, dose)])
stopifnot(!anyDuplicated(mapping$contributing_plates))
mapping[, plate_group := factor(contributing_plates, levels = plate_order)]
setorder(mapping, plate_group)
mapping[, axis_label := paste0(gsub(',', '+', contributing_plates, fixed = TRUE), '\n', as.character(dose))]
axis_labels <- setNames(mapping$axis_label, mapping$contributing_plates)
points[, plate_group := factor(contributing_plates, levels = plate_order)]
points[, drug := factor(drug, levels = drug_names)]
summaries <- points[, .(n_cells = .N, mean = mean(observed_displacement),
  median = median(observed_displacement), sd = sd(observed_displacement),
  q25 = unname(quantile(observed_displacement, .25)),
  q75 = unname(quantile(observed_displacement, .75))),
  by = .(contributing_plates, dose, drug, program)]
stopifnot(nrow(summaries) == 1200L)
fwrite(summaries, file.path(root, 'data/observed_program_boxplots_by_plate_summary.tsv'), sep = '\t')
fwrite(mapping, file.path(root, 'data/observed_program_boxplots_by_plate_key.tsv'), sep = '\t')

# Two wide program panels per page keep all plate groups and drug boxes readable.
# Each program has one shared y axis across every plate and drug on the page.
make_page <- function(selected_programs) {
  rows <- copy(points[program %in% selected_programs])
  rows[, program := factor(program, levels = selected_programs)]
  ggplot(rows, aes(plate_group, observed_displacement, fill = drug)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = '#999999', linewidth = .3) +
    geom_boxplot(position = position_dodge2(width = .82, preserve = 'single', padding = .12),
      width = .76, linewidth = .35, alpha = .85, outlier.size = .65, outlier.alpha = .6) +
    facet_wrap(~program, ncol = 1, scales = 'free_y') +
    scale_x_discrete(labels = axis_labels, drop = FALSE) +
    scale_fill_manual(values = drug_colors, labels = drug_labels, drop = FALSE, name = 'Drug') +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = '#F3F3F3'),
      axis.text.x = element_text(size = 9, lineheight = .95),
      legend.position = 'bottom', legend.key.height = grid::unit(.35, 'cm'),
      plot.margin = margin(6, 8, 5, 5)) +
    labs(title = 'v84 measured program displacements: drugs side by side within each plate',
      subtitle = 'Each box contains cell-line measurements for one drug on the indicated plate/dose group.',
      x = 'Contributing plate(s) (top label) / concentration in µM (bottom label)',
      y = 'Measured program displacement (signed RNA score)',
      caption = paste0('Boxes: median and interquartile range; whiskers: 1.5×IQR; points: outliers. Original matched-plate control; no rescaling or filtering.\n',
        'Groups marked A+B are already-pooled replicate profiles, shown separately from single plates. Dose and plate remain confounded.'))
}

book_path <- file.path(figure_root, 'all_drugs_40_programs_boxplots_by_plate.pdf')
cairo_pdf(book_path, width = 14, height = 7.5, onefile = TRUE)
for (page in 1:20) print(make_page(programs[(page - 1L) * 2L + 1:2]))
dev.off()
focus <- make_page(c('P21', 'P12'))
ggsave(file.path(figure_root, 'P21_P12_all_drugs_boxplots_by_plate.png'),
  focus, width = 14, height = 7.5, dpi = 160)
ggsave(file.path(figure_root, 'P21_P12_all_drugs_boxplots_by_plate.pdf'),
  focus, width = 14, height = 7.5, device = cairo_pdf)
stopifnot(identical(input_hash, digest(points_path, file = TRUE, algo = 'sha256')))
fwrite(data.table(path = c(points_path, script_path),
  sha256 = vapply(c(points_path, script_path), digest, character(1), file = TRUE, algo = 'sha256')),
  file.path(root, 'data/observed_program_boxplots_by_plate_input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'data/observed_program_boxplots_by_plate_session.txt'))
cat('Complete: 40 programs, all ten drugs, 15 explicit plate groups, 1,200 boxes; original inputs unchanged.\n')
