#!/usr/bin/env Rscript
# Show measured RNA progression against identically normalized inter-plate DMSO
# background. These plots do not estimate same-plate sampling error or fit models.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setDTthreads(1L)
root <- '/work/results/program_pilot/86_molecular_dose_and_plate_diagnosis'
prior <- '/work/results/program_pilot/84_full_panel_lincs_tahoe_joint_bridge'
observed <- fread(file.path(prior, 'data/observed_program_viability_points.tsv.gz'))
background <- fread(file.path(root, 'data/dmso_background_identical_normalizer.tsv'))
keys <- data.table(drug = c(rep('Lapatinib ditosylate', 3), 'Busulfan'),
  program = c('P11', 'P21', 'P40', 'P12'))
rows <- merge(observed, keys, by = c('drug', 'program'))
rows[, panel := paste(sub(' ditosylate', '', drug), program)]
rows[, dose_label := factor(as.character(dose), levels = c('0.05', '0.5', '5'))]
cohorts <- unique(rows[, .(depmap_id, drug, program, panel)])
bg <- merge(background, cohorts, by = c('depmap_id', 'program'), allow.cartesian = TRUE)
bands <- bg[, .(lower = quantile(displacement, .025), upper = quantile(displacement, .975)), by = panel]
means <- rows[, .(mean = mean(observed_displacement), sd = sd(observed_displacement), cells = .N),
  by = .(panel, dose_label)]
p <- ggplot(rows, aes(dose_label, observed_displacement)) +
  geom_rect(data = bands, aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
    inherit.aes = FALSE, fill = '#D8D8D8', alpha = .6) +
  geom_hline(yintercept = 0, linewidth = .3, color = '#666666') +
  geom_point(aes(color = dose_label), position = position_jitter(width = .08, height = 0, seed = 86),
    size = 1.3, alpha = .75) +
  geom_line(data = means, aes(y = mean, group = 1), color = 'black', linewidth = .5) +
  geom_errorbar(data = means, aes(y = mean, ymin = mean - sd, ymax = mean + sd),
    color = 'black', linewidth = .5, width = .13) +
  geom_point(data = means, aes(y = mean), color = 'black', shape = 18, size = 2.5) +
  facet_wrap(~panel, ncol = 2, scales = 'free_y') +
  scale_color_manual(values = c('0.05' = '#1879a4', '0.5' = '#df9b32', '5' = '#ae3769'),
    name = 'RNA dose (µM)') +
  theme_bw(base_size = 11) + theme(aspect.ratio = 1, legend.position = 'bottom',
    panel.grid.minor = element_blank()) +
  labs(x = 'Measured molecular dose (µM)', y = 'Measured program displacement',
    title = 'RNA dose progression and DMSO plate background',
    subtitle = 'Each colored point is one measured cell line; black markers and bars show mean ± SD across cells.',
    caption = paste0('Grey band: central 95% of DMSO plate-versus-other-13-plates contrasts in the same cell lines,\n',
      'using the identical full-gene normalizer and frozen program weights. This includes inter-plate control-state drift;\n',
      'it is not an estimate of error remaining after same-plate drug/control subtraction.'))
ggsave(file.path(root, 'figures/selected_program_dose_progression_vs_background.pdf'), p,
  width = 8.4, height = 9.1, device = cairo_pdf)
ggsave(file.path(root, 'figures/selected_program_dose_progression_vs_background.png'), p,
  width = 8.4, height = 9.1, dpi = 150)

# Compare variation across doses on identical recipient sets within each drug.
wide <- dcast(observed, depmap_id + drug + program ~ dose, value.var = 'observed_displacement')
setnames(wide, c('0.05', '0.5', '5'), c('low', 'middle', 'high'))
complete <- wide[is.finite(low) & is.finite(middle) & is.finite(high)]
spread <- complete[, .(cells = .N, sd_low = sd(low), sd_middle = sd(middle), sd_high = sd(high),
  mean_low = mean(low), mean_middle = mean(middle), mean_high = mean(high)), by = .(drug, program)]
spread[, `:=`(higher_high_dose_variation = sd_high > sd_low,
  increasing_variation = sd_low <= sd_middle & sd_middle <= sd_high,
  middle_mean_outside_endpoints = mean_middle < pmin(mean_low, mean_high) |
    mean_middle > pmax(mean_low, mean_high))]
fwrite(spread, file.path(root, 'data/same_recipient_dose_variation.tsv'), sep = '\t')
summary <- spread[, .(cells = cells[1], programs = .N,
  programs_with_more_high_dose_variation = sum(higher_high_dose_variation),
  programs_with_increasing_variation = sum(increasing_variation),
  programs_with_middle_excursion = sum(middle_mean_outside_endpoints)), by = drug]
fwrite(summary, file.path(root, 'data/same_recipient_dose_variation_summary.tsv'), sep = '\t')
print(summary)
