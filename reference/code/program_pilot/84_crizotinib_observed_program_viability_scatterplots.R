#!/usr/bin/env Rscript
# Crizotinib diagnostic using exactly the v84 measured-program scatter definition.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(qs);library(digest)})
setDTthreads(1L)
root <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
archive_path <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
annotations_path <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data/original40_program_annotations.tsv'
figure_root <- file.path(root,'figures/observed_program_viability/crizotinib')
dir.create(figure_root,recursive=TRUE,showWarnings=FALSE)
programs <- sprintf('P%02d',1:40)
native_doses <- c(.05,.5,5)
colors <- c('0.05'='#1879a4','0.5'='#df9b32','5'='#ae3769')
pearson <- function(x,y) if(length(x)<3L || sd(x)<1e-12 || sd(y)<1e-12) NA_real_ else cor(x,y)
archive <- readRDS(archive_path)
selected_drugs <- 'crizotinib'
molecular <- as.data.table(archive$prep$m)[drug==selected_drugs]
annotations <- fread(annotations_path)
keys <- c('depmap_id','drug','dose')
prism <- '/work/cache/tahoe_prism_full_chain/prism'
treatment_path <- file.path(prism,'secondary-screen-replicate-treatment-info.csv')
raw_path <- file.path(prism,'secondary-screen-logfold-change.csv')
treatment <- fread(treatment_path)[tolower(name)=='crizotinib' & perturbation_type=='experimental_treatment']
stopifnot(nrow(treatment)>0L,!anyDuplicated(treatment$column_name))
header <- names(fread(raw_path,nrows=0))
raw <- fread(raw_path,select=c(header[1],treatment$column_name))
setnames(raw,1L,'depmap_id')
raw <- raw[depmap_id %in% molecular$depmap_id]
long_raw <- melt(raw,id.vars='depmap_id',variable.name='column_name',value.name='y',variable.factor=FALSE)
long_raw <- merge(long_raw,treatment[,.(column_name,dose,screen_id,broad_id)],by='column_name')
curve <- long_raw[is.finite(y),.(y=mean(y),replicates=.N),by=.(depmap_id,screen_id,broad_id,dose)]
native <- curve[,{o=order(dose); .(dose=native_doses,
  y=approx(log(dose[o]),y[o],xout=log(native_doses),rule=1)$y,
  weight=approx(log(dose[o]),replicates[o],xout=log(native_doses),rule=1)$y)},by=.(depmap_id,screen_id,broad_id)]
outcomes <- native[is.finite(y)&is.finite(weight),.(y=weighted.mean(y,weight)),by=.(depmap_id,dose)]
outcomes[,drug:='crizotinib']
stopifnot(!anyDuplicated(molecular[,..keys]),!anyDuplicated(outcomes[,..keys]),all(is.finite(as.matrix(molecular[,..programs]))))
fwrite(curve,file.path(figure_root,'prism_measured_curve_points.tsv'),sep='\t')
fwrite(native,file.path(figure_root,'prism_native_dose_screen_estimates.tsv'),sep='\t')
# Retain every matched RNA profile. Unmatched profiles are reported explicitly;
# missing functional outcomes are never filled with predictions.
coverage_rows <- merge(molecular[, c(keys, 'contributing_plates'), with = FALSE],
  outcomes[, c(keys, 'y'), with = FALSE], by = keys, all.x = TRUE)
coverage <- coverage_rows[, .(rna_cells = uniqueN(depmap_id),
  matched_cells = sum(is.finite(y)), unmatched_cells = sum(!is.finite(y)),
  rna_plates = paste(sort(unique(contributing_plates)), collapse = ',')), by = .(drug, dose)]
matched <- merge(molecular, outcomes[, c(keys, 'y'), with = FALSE], by = keys)
stopifnot(nrow(matched) == sum(coverage$matched_cells), setequal(matched$drug, selected_drugs))
metadata <- as.data.table(archive$meta)[, .(depmap_id = ModelID, cell_line = CellLineName)]
stopifnot(!anyDuplicated(metadata$depmap_id), all(matched$depmap_id %in% metadata$depmap_id))
matched <- merge(matched, metadata, by = 'depmap_id')
setnames(matched, 'y', 'observed_log2_abundance')
matched[, observed_percent_vehicle := 100 * 2^observed_log2_abundance]
long <- melt(matched, id.vars = c(keys, 'cell_line', 'contributing_plates',
  'observed_log2_abundance', 'observed_percent_vehicle'), measure.vars = programs,
  variable.name = 'program', value.name = 'observed_displacement', variable.factor = FALSE)
long[, dose_label := factor(as.character(dose), levels = names(colors))]
long[, `:=`(centered_displacement = observed_displacement - mean(observed_displacement),
  centered_response = observed_log2_abundance - mean(observed_log2_abundance)),
  by = .(drug, dose, program)]

# Report both axes of association: pooling doses, and comparing cells after each
# dose mean is removed. Also retain all three separate-dose correlations.
pooled <- long[, .(cells = uniqueN(depmap_id), cell_dose_points = .N,
  pooled_r = pearson(observed_displacement, observed_log2_abundance),
  within_dose_r = pearson(centered_displacement, centered_response)), by = .(drug, program)]
by_dose <- long[, .(cells = .N,
  r = pearson(observed_displacement, observed_log2_abundance),
  mean_displacement = mean(observed_displacement), sd_displacement = sd(observed_displacement),
  mean_log2_abundance = mean(observed_log2_abundance),
  sensitive_below_50pct = sum(observed_log2_abundance < -1)), by = .(drug, dose, program)]

# An explicit dose contrast cancels a cell's constant response component. This
# table tests whether the observed change in a program accompanies the observed
# change in functional response; it is not used to train or select anything.
contrasts <- list()
for (i in 1:2) {
  lower <- long[dose == native_doses[i], .(depmap_id, drug, program,
    displacement_low = observed_displacement, response_low = observed_log2_abundance)]
  upper <- long[dose == native_doses[i + 1L], .(depmap_id, drug, program,
    displacement_high = observed_displacement, response_high = observed_log2_abundance)]
  pairs <- merge(lower, upper, by = c('depmap_id', 'drug', 'program'))
  pairs[, `:=`(dose_low = native_doses[i], dose_high = native_doses[i + 1L],
    displacement_change = displacement_high - displacement_low,
    response_change = response_high - response_low)]
  contrasts[[i]] <- pairs
}
contrast_points <- rbindlist(contrasts)
contrast_metrics <- contrast_points[, .(cells = .N,
  dose_change_r = pearson(displacement_change, response_change),
  mean_displacement_change = mean(displacement_change),
  mean_response_log2_change = mean(response_change)), by = .(drug, program, dose_low, dose_high)]

fwrite(coverage, file.path(figure_root, 'observed_program_viability_coverage.tsv'), sep = '\t')
fwrite(long, file.path(figure_root, 'observed_program_viability_points.tsv.gz'), sep = '\t')
fwrite(pooled, file.path(figure_root, 'observed_program_viability_pooled_correlations.tsv'), sep = '\t')
fwrite(by_dose, file.path(figure_root, 'observed_program_viability_each_dose_correlations.tsv'), sep = '\t')
fwrite(contrast_points, file.path(figure_root, 'observed_program_viability_dose_contrast_points.tsv.gz'), sep = '\t')
fwrite(contrast_metrics, file.path(figure_root, 'observed_program_viability_dose_contrast_correlations.tsv'), sep = '\t')
fwrite(annotations, file.path(figure_root, 'observed_program_viability_program_annotations.tsv'), sep = '\t')

# Ten square panels per page keep all forty programs readable, in their frozen
# numerical order. The three colored trend lines each compare cells at one dose.
# Axis ranges are independent; no unit matching or displacement scaling is used.
pages <- list()
drug_index <- list()
for (drug_value in selected_drugs) {
  slug <- gsub('^_|_$', '', gsub('[^a-z0-9]+', '_', tolower(drug_value)))
  drug_path <- file.path(figure_root, paste0(slug, '_all_40_programs.pdf'))
  drug_pages <- list()
  for (page in 1:4) {
    these_programs <- programs[(page - 1L) * 10L + 1:10]
    rows <- long[drug == drug_value & program %in% these_programs]
    rows[, program := factor(program, levels = these_programs)]
    labels <- pooled[drug == drug_value & program %in% these_programs]
    labels[, program := factor(program, levels = these_programs)]
    labels[, label := sprintf('Pooled r = %.2f\nWithin-dose r = %.2f', pooled_r, within_dose_r)]
    p <- ggplot(rows, aes(observed_displacement, observed_log2_abundance)) +
      geom_point(aes(color = dose_label), size = .8, alpha = .65) +
      geom_smooth(aes(color = dose_label), method = 'lm', formula = y ~ x,
        se = FALSE, linewidth = .45) +
      geom_label(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -.02, vjust = 1.03, size = 2.4, linewidth = 0, fill = '#FFFFFFCC',
        label.padding = grid::unit(.06, 'lines')) +
      facet_wrap(~program, ncol = 5, scales = 'free_x') +
      scale_color_manual(values = colors, drop = FALSE, name = 'Measured RNA dose (µM)') +
      scale_y_continuous(expand = expansion(mult = c(.04, .16))) +
      theme_bw(base_size = 10) +
      theme(aspect.ratio = 1, legend.position = 'bottom',
        panel.grid.minor = element_blank(), strip.background = element_rect(fill = '#F3F3F3'),
        plot.margin = margin(6, 8, 5, 5)) +
      labs(x = 'Measured program displacement (signed RNA score)',
        y = 'Observed PRISM log2 abundance / vehicle',
        title = paste0(drug_value, ': observed program changes versus functional response'),
        subtitle = paste0('Each point is one measured cell–dose profile; programs ',
          these_programs[1], '–', these_programs[10], '. Colored lines compare cells at the same dose.'),
        caption = paste0('PRISM curves interpolated to Tahoe doses; these are separate experiments, not paired wells. ',
          '0 = 100%, −1 = 50%, −2 = 25% of vehicle abundance.\n',
          'Pooled r mixes dose and cell effects; within-dose r subtracts the mean at each dose from both axes. ',
          'All matched profiles are included.'))
    drug_pages[[page]] <- p
    pages[[paste0(slug, '_', page)]] <- p
    if (drug_value == 'crizotinib') {
      ggsave(file.path(figure_root, paste0(slug, '_page_', page, '.png')),
        p, width = 12.5, height = 6.5, dpi = 160)
    }
  }
  grDevices::cairo_pdf(drug_path, width = 12.5, height = 6.5)
  for (p in drug_pages) print(p)
  dev.off()
  drug_index[[drug_value]] <- data.table(drug = drug_value, pdf = drug_path, pages = 4L,
    measured_cells = uniqueN(long[drug == drug_value]$depmap_id),
    cell_dose_profiles = nrow(matched[drug == drug_value]))
}
grDevices::cairo_pdf(file.path(figure_root, 'crizotinib_all_doses_all_40_programs.pdf'), width = 12.5, height = 6.5,
  family = "sans")
for (p in pages) print(p)
dev.off()
fwrite(rbindlist(drug_index), file.path(figure_root, 'observed_program_viability_figure_index.tsv'), sep = '\t')

print(coverage)
print(by_dose[order(-abs(r))][1:12])
writeLines(c('Crizotinib uses the identical measured displacement / interpolated PRISM definition as v84.',
  'All matched profiles retained. Dose/plate confounding remains; no model fitted.',capture.output(sessionInfo())),
  file.path(figure_root,'run_details.txt'))
