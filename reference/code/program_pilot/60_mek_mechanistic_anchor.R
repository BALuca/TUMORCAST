#!/usr/bin/env Rscript
# Orthogonal measured signaling and inhibitor-growth controls for the MEK pilot.
# Source studies and assay protocols remain separate; no RNA-based mechanism is
# invented and no combined-treatment RNA is inferred from these measurements.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
root <- '/work/results/program_pilot/60_defensible_response_composition'
prior <- '/work/results/program_pilot/59_mechanistic_evidence_datasets/data'
input_paths <- c(file.path(prior, 'mek_measurements.tsv.gz'),
                 file.path(prior, 'mek_growth_biological_replicates.tsv'))
stopifnot(all(file.exists(input_paths)))
for (folder in c('data', 'figures')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0('mek_anchor_', name, '.tsv.gz')), sep = '\t')
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0('mek_anchor_', name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0('mek_anchor_', name, '.png')), p, width = width, height = height, dpi = 150)
}
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))
m <- fread(input_paths[1])[study_id == 'vanBentum2025']
stopifnot(all(is.finite(m$value)))

# The displayed Luminex controls are copied into multiple ligand panels. We use
# only each source ligand panel's supplied controls, never pool those copies as
# independent replicates. Source lacks biological experiment IDs for these rows;
# points are source observations and uncertainty is descriptive mean +/- SD.
akt <- m[assay == 'Luminex_phosphoprotein']
stopifnot(all(akt$dose_A_uM %in% c(0, 1)), all(akt$ligand_duration_min == 20))
akt[, MEK_inhibitor := ifelse(dose_A_uM == 1, 'Selumetinib 1 µM', 'Vehicle')]
akt[, ligand_arm := ifelse(ligand == 'BSA', 'BSA control', 'Growth factor')]
akt[, ligand_arm := factor(ligand_arm, levels = c('BSA control', 'Growth factor'))]
akt_summary <- akt[, .(mean_log2 = mean(value), sd_log2 = sd(value), n_source_observations = .N),
                   by = .(cell_line, ligand_panel, ligand_arm, MEK_inhibitor)]
save_table(akt, 'pAKT_source')
save_table(akt_summary, 'pAKT_summary')
akt_wide <- dcast(akt_summary, cell_line + ligand_panel ~ ligand_arm + MEK_inhibitor, value.var = 'mean_log2')
stopifnot(!anyNA(akt_wide))
akt_effect <- akt_wide[, .(cell_line, ligand_panel,
  MEKi_effect_without_ligand_log2 = `BSA control_Selumetinib 1 µM` - `BSA control_Vehicle`,
  MEKi_effect_with_ligand_log2 = `Growth factor_Selumetinib 1 µM` - `Growth factor_Vehicle`,
  ligand_effect_without_MEKi_log2 = `Growth factor_Vehicle` - `BSA control_Vehicle`,
  ligand_effect_with_MEKi_log2 = `Growth factor_Selumetinib 1 µM` - `BSA control_Selumetinib 1 µM`)]
akt_effect[, interaction_log2 := MEKi_effect_with_ligand_log2 - MEKi_effect_without_ligand_log2]
akt_effect[, interpretation := 'Descriptive factorial contrast; biological replicate IDs are unavailable in display source']
save_table(akt_effect, 'pAKT_factorial_contrasts')
p <- ggplot(akt_summary, aes(ligand_arm, mean_log2, color = MEK_inhibitor, group = MEK_inhibitor)) +
  geom_line(linewidth = .6, position = position_dodge(width = .16)) +
  geom_errorbar(aes(ymin = mean_log2 - sd_log2, ymax = mean_log2 + sd_log2),
                position = position_dodge(width = .16), width = .08) +
  geom_point(size = 1.8, position = position_dodge(width = .16)) +
  facet_grid(cell_line ~ ligand_panel) + scale_color_manual(values = c('#D55E00', '#0072B2')) +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = 'Direct pAKT measurement: ligand responses depend on cell and MEK inhibition',
    subtitle = 'Selumetinib 1 µM: 3 h 40 min pretreatment, then 20 min growth factor; four full factorial comparisons',
    x = NULL, y = 'pAKT abundance change (log2)', color = NULL,
    caption = 'Means +/- SD of displayed source observations; repeated BSA controls are not independent evidence across panels.\nHGF is one of four ligands. A higher combined signal alone does not prove a positive interaction.')
save_plot(p, 'pAKT_factorials', 11, 8.2)

# Direct phosphosite assay: source intensities at four hours, independently of
# Luminex and of growth. Retain the ambiguous JNK1/JNK3 site assignment.
ph <- m[assay == 'SPIED_DIA_phosphoproteomics']
stopifnot(all(ph$dose_A_uM %in% c(0, 10)), all(ph$ligand_duration_min == 30))
ph[, treatment := factor(treatment, levels = c('DMSO+BSA', 'MEKi+BSA', 'DMSO+GFmix', 'MEKi+GFmix'))]
ph_summary <- ph[, .(mean_log10 = mean(value), sd_log10 = sd(value), n_biological_replicates = .N),
                 by = .(cell_line, feature, treatment)]
save_table(ph, 'phosphosite_source')
save_table(ph_summary, 'phosphosite_summary')
ph_wide <- dcast(ph_summary, cell_line + feature ~ treatment, value.var = 'mean_log10')
ph_effect <- ph_wide[, .(cell_line, feature,
  MEKi_effect_BSA_log10 = `MEKi+BSA` - `DMSO+BSA`,
  MEKi_effect_GFmix_log10 = `MEKi+GFmix` - `DMSO+GFmix`,
  interaction_log10 = `MEKi+GFmix` - `MEKi+BSA` - `DMSO+GFmix` + `DMSO+BSA`)]
save_table(ph_effect, 'phosphosite_factorial_contrasts')
p <- ggplot(ph, aes(treatment, value)) +
  geom_point(aes(color = treatment), position = position_jitter(width = .06, seed = 6011), size = 1.6) +
  geom_errorbar(data = ph_summary, aes(y = mean_log10, ymin = mean_log10 - sd_log10, ymax = mean_log10 + sd_log10),
                width = .12, color = 'black', inherit.aes = TRUE) +
  geom_point(data = ph_summary, aes(y = mean_log10), shape = 95, size = 5) +
  facet_grid(feature ~ cell_line, scales = 'free_y') + theme(legend.position = 'none',
    axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = 'ERK inhibition and alternative JNK signaling are measured separately',
    subtitle = 'Selumetinib 10 µM for 3.5 h, then 30 min EGF/HGF/FGF2/VEGF-C mix; independent phosphoproteomic assay',
    x = NULL, y = 'Phosphopeptide intensity (log10)',
    caption = 'Each point is a biological replicate (normally 3; DLD-1 GFmix-only has2). Lines show mean +/- SD.\nThe JNK phosphopeptide maps to JNK1 Y185 / JNK3 Y223; no isoform assignment is fabricated.')
save_plot(p, 'ERK_JNK_phosphosites', 10.3, 6.8)

# Functional comparison uses source biological-replicate means, not technical
# wells. Within-replicate vehicle scaling avoids conflating separate experiments.
growth <- fread(input_paths[2])[feature == 'doubling_time']
stopifnot(all(is.finite(growth$value)), all(growth$value > 0),
          !anyDuplicated(growth[, .(cell_line, treatment, dose_A_uM, dose_B_uM, replicate)]))
growth[, partner := sub('selumetinib\\+', '', treatment)]
baseline <- growth[dose_A_uM == 0 & dose_B_uM == 0,
  .(cell_line, treatment, replicate, vehicle_doubling_h = value)]
growth <- merge(growth, baseline, by = c('cell_line', 'treatment', 'replicate'))
growth[, relative_growth_rate_percent := 100 * vehicle_doubling_h / value]
growth[, partner_dose := factor(dose_B_uM)]
gs <- growth[, .(mean_doubling_h = mean(value), sd_doubling_h = sd(value),
  mean_relative_rate = mean(relative_growth_rate_percent), sd_relative_rate = sd(relative_growth_rate_percent),
  n_biological_replicates = uniqueN(replicate)), by = .(cell_line, partner, dose_A_uM, dose_B_uM, partner_dose)]
stopifnot(all(gs$n_biological_replicates == 3L))
save_table(growth, 'growth_biological_replicates')
save_table(gs, 'growth_summary')
p <- ggplot(gs, aes(dose_A_uM, mean_relative_rate, color = partner_dose, group = partner_dose)) +
  geom_hline(yintercept = 100, linetype = 2, color = 'grey60') +
  geom_line(linewidth = .65) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_relative_rate - sd_relative_rate, ymax = mean_relative_rate + sd_relative_rate), width = .04) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = .1), breaks = c(0, .2, 1, 5)) +
  facet_grid(cell_line ~ partner) + theme(legend.position = 'bottom') +
  labs(title = 'Does blocking the alternative route change growth?',
    subtitle = 'Complete single-agent and combination dose grids; three biological experiments per condition',
    x = 'Selumetinib (µM)', y = 'Growth rate (% of matched vehicle)', color = 'Partner dose (µM)',
    caption = 'Rate = vehicle doubling time / treated doubling time; mean +/- SD across biological experiments.\nGrowth protocols are separate from ligand-stimulated phosphoproteomics. Slower growth is not equivalent to cell death.')
save_plot(p, 'complete_growth_dose_grids', 10, 7)

# At1µM of each drug, show each biological experiment and the measured single
# controls. This dose is shared by both blocker grids; full curves remain above.
anchor <- growth[dose_A_uM %in% c(0, 1) & dose_B_uM %in% c(0, 1)]
anchor[, arm := factor(paste(dose_A_uM, dose_B_uM), levels = c('0 0', '1 0', '0 1', '1 1'),
                       labels = c('Vehicle', 'MEKi', 'Partner', 'MEKi + partner'))]
anchor_summary <- anchor[, .(mean_h = mean(value), sd_h = sd(value), n_biological_replicates = .N),
                        by = .(cell_line, partner, arm)]
save_table(anchor, 'one_uM_functional_factorial')
save_table(anchor_summary, 'one_uM_functional_summary')
aw <- dcast(anchor, cell_line + partner + replicate ~ arm, value.var = 'value')
aw[, `:=`(extra_doubling_hours_vs_MEKi = `MEKi + partner` - MEKi,
           extra_doubling_hours_vs_partner = `MEKi + partner` - Partner)]
ae <- aw[, .(mean_extra_h_vs_MEKi = mean(extra_doubling_hours_vs_MEKi),
  experiments_slower_than_MEKi = sum(extra_doubling_hours_vs_MEKi > 0),
  experiments_slower_than_partner = sum(extra_doubling_hours_vs_partner > 0), n = .N),
  by = .(cell_line, partner)]
save_table(ae, 'one_uM_functional_effects')
p <- ggplot(anchor, aes(arm, value, group = replicate)) +
  geom_line(color = 'grey65', linewidth = .4) + geom_point(aes(color = replicate), size = 2) +
  facet_grid(cell_line ~ partner, scales = 'free_y') + theme(legend.position = 'none',
    axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = 'Functional controls: how long cells take to double',
    subtitle = 'All single and combined treatments at 1 µM; each line connects one biological experiment',
    x = NULL, y = 'Doubling time (hours)',
    caption = 'Each drug-pair experiment has its own matched vehicle. Different control doubling times are preserved.\nThese comparisons show the consequence of blockade; they do not establish an HGF-rescue mechanism.')
save_plot(p, 'functional_controls_hours', 9.2, 6.7)

save_table(data.table(path = input_paths,
  sha256 = vapply(input_paths, digest, character(1), file = TRUE, algo = 'sha256')), 'input_manifest')
report <- c('Measured MEK biological anchor: signaling and growth remain separate', '',
  'Data: van Bentum2025, Nature Communications10.1038/s41467-025-59404-y.',
  'Direct pAKT: selumetinib1uM,3h40min then20min ligand. Four ligands and BSA controls are shown; HGF is not assumed to be synergistic.',
  'ERK/JNK phosphosites: selumetinib10uM,3.5h then30min GFmix; three biological replicates, except two in oneDLD-1arm.',
  'Growth: separate48h doubling-time analyses with three independent experiments and full single/combined dose grids.',
  'The drug differs from the planned Jaaks trametinib test; this is a same-pathway biological reference, not a matched validation outcome.',
  '', 'At1uM+1uM, added doubling time versus MEKi alone:',
  vapply(seq_len(nrow(ae)), function(i) sprintf('%s + %s: %.2f additional hours; slower than MEKi in %d/%d experiments; slower than partner in%d/%d.',
    ae$cell_line[i], ae$partner[i], ae$mean_extra_h_vs_MEKi[i], ae$experiments_slower_than_MEKi[i], ae$n[i],
    ae$experiments_slower_than_partner[i], ae$n[i]), character(1)),
  '', 'HGF pAKT factorial effects (log2 scale; descriptive, no invented biological-replicate pairing):',
  capture.output(print(akt_effect[ligand_panel == 'HGF', .(cell_line,
    MEKi_effect_without_ligand_log2, MEKi_effect_with_ligand_log2, interaction_log2)])),
  '', 'Use: these independent measurements identify which signaling direction a molecular representation should retain, and which cell contexts display a functional consequence of inhibiting an alternative route.',
  'Limits: no same-protocol RNA/functional chain is measured here. Positive pAKT does not itself establish survival rescue. Slower doubling does not establish killing.',
  'BSA controls occur repeatedly across display panels; no inferential p-values or inflated replicate counts are generated from those copies.')
writeLines(report, file.path(root, 'MEK_MECHANISTIC_ANCHOR.txt'))
capture.output(sessionInfo(), file = file.path(root, 'data/mek_anchor_session_info.txt'))
print(ae)
print(ph_effect)
message('MEK anchor complete.')
