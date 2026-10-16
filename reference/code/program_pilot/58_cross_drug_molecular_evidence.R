# Fixed biological readouts versus PRISM abundance across drugs.
# This is part of analysis58: molecular information before program compression.
# No readout, gene, drug or dose is selected on these functional outcomes.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(digest)
})
source('/work/src/R/lib_heatmap.R')

root <- '/work/results/program_pilot/58_molecular_function_information_audit'
prior56 <- '/work/results/program_pilot/56_drug_program_biology_reliability/data'
prior57 <- '/work/results/program_pilot/57_erbb_met_hgf_program_biology/data'
bootstrap_draws <- 199L
set.seed(5802)
dir.create(file.path(root, 'data'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, 'figures'), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0('cross_drug_', name, '.tsv.gz')), sep = '\t')
}
save_plot <- function(x, name, width, height) {
  path <- file.path(root, 'figures', paste0('cross_drug_', name))
  ggsave(paste0(path, '.pdf'), x, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(path, '.png'), x, width = width, height = height, dpi = 150, limitsize = FALSE)
}
paths <- c(file.path(prior56, c('mechanism_reference_observations.tsv.gz',
  'within_dose_replicate_agreement.tsv.gz', 'mechanism_reference_membership.tsv.gz')),
  file.path(prior57, c('drug_biological_readouts.tsv.gz',
  'within_dose_replicate_agreement.tsv.gz', 'biological_readout_membership.tsv.gz')))
stopifnot(all(file.exists(paths)))
save_table(data.table(path = paths, sha256 = vapply(paths, function(p) digest(file = p,
  algo = 'sha256'), character(1))), 'input_manifest')

# Retain every saved observation, including flat doses and contradictory cases.
# Earlier56 measured named genes present in bank43. Earlier57 additionally retained
# all measured named genes; use the common bank-gene subset for exact reconstruction comparisons.
a <- fread(paths[1])
b <- fread(paths[4])
stopifnot(!anyDuplicated(a[, .(observation_id, readout)]),
  !anyDuplicated(b[, .(observation_id, readout)]))
a[, readout := sub('HALLMARK_', '', readout)]
a <- a[, .(observation_id, drug, depmap_id, recipient, dose, readout,
  gene_score = measured_gene_score, all_gene_score = NA_real_,
  program_score = program_reconstruction, observed_log2, measured_pct,
  repeats = replicates, repeat_low = replicate_min_pct, repeat_high = replicate_max_pct,
  source = '56: named genes retained in bank43')]
b <- b[, .(observation_id, drug, depmap_id, recipient, dose, readout,
  gene_score = bank_genes, all_gene_score = measured_genes,
  program_score = program_reconstruction, observed_log2, measured_pct,
  repeats, repeat_low, repeat_high, source = '57: named genes before/after bank admission')]
overlap <- merge(a, b, by = c('observation_id', 'readout'))
stopifnot(max(abs(overlap$gene_score.x - overlap$gene_score.y)) < 1e-10,
  max(abs(overlap$program_score.x - overlap$program_score.y)) < 1e-10)
measurements <- rbind(b, a[!b, on = .(observation_id, readout)])
stopifnot(!anyDuplicated(measurements[, .(observation_id, readout)]))
save_table(measurements[!is.finite(observed_log2), .N, by = .(drug, readout)], 'missing_function')
measurements <- measurements[is.finite(observed_log2)]
stopifnot(all(is.finite(measurements$gene_score)), all(is.finite(measurements$program_score)),
  all(measurements$measured_pct > 0))
save_table(measurements, 'observations')

# Pearson r is summarized two ways: whole dose curves, and between-cell differences
# after removal of each drug-dose mean. All doses enter both estimates.
# Whole cells are resampled, carrying every dose together. Recompute dose means
# in every bootstrap sample. These are descriptive CIs, not independent validation.
cor_safe <- function(x, y) {
  if (length(x) < 5L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y)
}
cor_views <- function(d, column) {
  x <- d[[column]]
  y <- d$observed_log2
  xc <- x - ave(x, d$dose, FUN = mean)
  yc <- y - ave(y, d$dose, FUN = mean)
  c(full_dose_curve = cor_safe(x, y), dose_centered = cor_safe(xc, yc))
}
summary_rows <- list()
counter <- 0L
for (drug_name in sort(unique(measurements$drug))) {
  for (readout_name in sort(unique(measurements[drug == drug_name, readout]))) {
    d <- measurements[drug == drug_name & readout == readout_name]
    ids <- unique(d$depmap_id)
    indexes <- split(seq_len(nrow(d)), factor(d$depmap_id, levels = ids))
    boot_indexes <- replicate(bootstrap_draws,
      unlist(indexes[sample(seq_along(ids), length(ids), replace = TRUE)], use.names = FALSE), simplify = FALSE)
    for (column in c('gene_score', 'program_score', 'all_gene_score')) {
      if (all(is.na(d[[column]]))) next
      stopifnot(all(is.finite(d[[column]])))
      r <- cor_views(d, column)
      draws <- vapply(boot_indexes, function(i) cor_views(d[i], column), numeric(2))
      for (view in names(r)) {
        counter <- counter + 1L
        interval <- quantile(draws[view, ], c(.025, .975), na.rm = TRUE)
        summary_rows[[counter]] <- data.table(drug = drug_name, readout = readout_name,
          representation = column, view, n = nrow(d), cells = length(ids), doses = uniqueN(d$dose),
          r = r[[view]], lower = interval[1], upper = interval[2], bootstrap_draws,
          sensitive_observations = sum(d$measured_pct <= 50),
          repeated_observations = sum(d$repeats >= 2, na.rm = TRUE),
          gene_program_agreement = cor_safe(d$gene_score, d$program_score))
      }
    }
  }
  message('Readout audit finished: ', drug_name)
}
associations <- rbindlist(summary_rows)
save_table(associations, 'associations')

# Endpoint reproducibility is technical: it does not validate the separate Tahoe
# molecular assay, and low correlations in nearly flat contexts do not establish noise.
ra <- fread(paths[2])
rb <- fread(paths[5])
stopifnot(identical(names(ra), names(rb)))
reliability <- unique(rbind(ra, rb))
reliability <- reliability[drug %in% measurements$drug]
save_table(reliability, 'replicate_pairs')
repeat_summary <- reliability[, .(repeat_pairs = .N, min_cells = min(cells),
  median_repeat_r = median(r, na.rm = TRUE), min_repeat_r = min(r, na.rm = TRUE),
  max_repeat_r = max(r, na.rm = TRUE)), by = .(drug, dose)]
observations <- unique(measurements[, .(observation_id, drug, dose, depmap_id,
  measured_pct, repeats, repeat_low, repeat_high)])
counts <- observations[, .(observations = .N, sensitive = sum(measured_pct <= 50),
  repeated = sum(repeats >= 2, na.rm = TRUE),
  repeat_consistent_sensitive = sum(repeats >= 2 & repeat_high <= 50, na.rm = TRUE),
  repeat_consistent_above50 = sum(repeats >= 2 & repeat_low > 50, na.rm = TRUE),
  abundance_q10 = quantile(measured_pct, .1), abundance_q90 = quantile(measured_pct, .9)), by = .(drug, dose)]
save_table(merge(counts, repeat_summary, by = c('drug', 'dose'), all.x = TRUE), 'endpoint_reliability')

# Predeclared mechanism/drug pairs include all requested siblings, not best hits.
mechanisms <- data.table(drug = c('afatinib', 'lapatinib', 'neratinib', 'trametinib',
  'cobimetinib', 'palbociclib', 'ribociclib', 'bortezomib', 'ixazomib'),
  readout = c(rep('ERK_output', 5), rep('E2F_TARGETS', 2), rep('UNFOLDED_PROTEIN_RESPONSE', 2)),
  family = c(rep('ERBB', 3), rep('MEK', 2), rep('CDK4/6', 2), rep('Proteasome', 2)),
  expected_sign = c(rep(1, 7), -1, -1))
save_table(merge(mechanisms, associations, by = c('drug', 'readout')), 'prespecified_mechanisms')

# Compact forest plots expose all five commonly measured readouts for all14 drugs.
core <- c('ERK_output', 'E2F_TARGETS', 'MYC_TARGETS_V1', 'UNFOLDED_PROTEIN_RESPONSE', 'APOPTOSIS')
labels <- c(ERK_output = 'ERK output', E2F_TARGETS = 'E2F targets', MYC_TARGETS_V1 = 'MYC targets',
  UNFOLDED_PROTEIN_RESPONSE = 'Unfolded protein response', APOPTOSIS = 'Apoptosis-associated RNA')
forest <- associations[readout %in% core & representation != 'all_gene_score']
forest[, drug := factor(drug, levels = rev(sort(unique(drug))))]
forest[, readout := factor(readout, levels = core)]
forest[, representation := factor(representation, c('gene_score', 'program_score'),
  c('Named genes', 'Reconstructed from program means'))]
pd <- position_dodge(width = .55)
for (view_name in c('full_dose_curve', 'dose_centered')) {
  p <- ggplot(forest[view == view_name], aes(r, drug, color = representation)) +
    geom_vline(xintercept = 0, color = 'grey60', linewidth = .3) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = .25, position = pd, linewidth = .35) +
    geom_point(position = pd, size = 1.5) +
    facet_wrap(~readout, nrow = 1, labeller = as_labeller(labels)) +
    scale_color_manual(values = c('#1F78B4', '#D95F02')) +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-.5, 0, .5)) +
    labs(x = 'Pearson r with log2 PRISM relative abundance; 95% cell-block bootstrap interval', y = NULL,
      color = NULL, title = if (view_name == 'full_dose_curve')
        'All doses together: molecular response versus relative abundance' else
        'Across-cell differences after subtracting drug-dose means; all doses retained',
      subtitle = 'All14 drugs with saved readouts shown. Positive r: suppression accompanies lower abundance; stress readouts may run oppositely.',
      caption = 'Named genes are restricted to those present in bank43 for this matched compression comparison. Descriptive, retrospective associations; no predictor fitted.') +
    theme_bw(base_size = 10) + theme(legend.position = 'bottom', strip.text = element_text(size = 9),
      plot.margin = margin(8, 12, 8, 8), panel.grid.minor = element_blank())
  save_plot(p, paste0('all_drugs_', view_name), 15, 6.5)
}

# Full trajectories retain the dose transitions the user identified visually.
plot_data <- merge(measurements, mechanisms, by = c('drug', 'readout'))
plot_data[, panel := paste0(drug, '\n', labels[readout])]
plot_data[, panel := factor(panel, levels = unique(paste0(mechanisms$drug, '\n', labels[mechanisms$readout])))]
plot_data[, log10_dose := log10(dose)]
setorder(plot_data, drug, depmap_id, dose)
caption <- merge(mechanisms, associations[representation == 'gene_score' & view == 'full_dose_curve'],
  by = c('drug', 'readout'))
caption[, panel := paste0(drug, '\n', labels[readout])]
caption[, label := sprintf('r = %.2f [%.2f, %.2f]\n%d cells; %d cell-dose observations', r, lower, upper, cells, n)]
p <- ggplot(plot_data, aes(gene_score, observed_log2)) +
  geom_path(aes(group = depmap_id), color = 'grey75', linewidth = .25, alpha = .65) +
  geom_point(aes(color = log10_dose), size = 1.4, alpha = .8) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_hline(yintercept = -1, linetype = 'dashed', color = 'grey50', linewidth = .3) +
  geom_label(data = caption, aes(label = label), x = -Inf, y = -Inf, hjust = -.02, vjust = -.05,
    inherit.aes = FALSE, size = 2.6, linewidth = 0, fill = 'white') +
  facet_wrap(~panel, ncol = 3, scales = 'free') +
  scale_color_viridis_c(name = 'Dose (log10 µM)') +
  scale_y_continuous(breaks = c(-6, -4, -2, -1, 0, 1), labels = c('1.6', '6.3', '25', '50', '100', '200')) +
  labs(x = 'Named-gene molecular displacement (log2 expression change)',
    y = 'PRISM relative abundance (% vehicle; log2 spacing)',
    title = 'Do established molecular responses accompany functional effects across sibling drugs?',
    subtitle = 'Each line connects one cell across doses. Dashed line: 50% abundance. All cells and available doses retained.',
    caption = 'Fixed ERK, E2F and unfolded-protein-response readouts; genes retained in bank43. Associations are descriptive, not held-out prediction.\nGrey trajectories preserve dose information; they do not establish matched exposure timing between Tahoe and PRISM.') +
  theme_bw(base_size = 10) + theme(aspect.ratio = 1, legend.position = 'bottom',
    panel.grid.minor = element_blank(), plot.margin = margin(8, 10, 8, 8))
save_plot(p, 'prespecified_response_trajectories', 11, 12)

# Save the difference between removing bank admission and removing compression.
membership56 <- fread(paths[3])
membership57 <- fread(paths[6])
save_table(membership56, 'readout_membership_56')
save_table(membership57, 'readout_membership_57')
writeLines(c(
  'CROSS-DRUG MOLECULAR INFORMATION AUDIT',
  sprintf('%d drugs; %d unique cell-dose observations; %d cell lines.',
    uniqueN(measurements$drug), uniqueN(measurements$observation_id), uniqueN(measurements$depmap_id)),
  'No model was fitted and no readout, drug or dose was selected on these outcomes.',
  'All available five core fixed readouts across14 drugs are shown; other57 readouts remain in the association table.',
  'Full-dose correlation includes dose effects and differences between cells. Dose-centered correlation removes drug-dose means but retains every dose.',
  '199 bootstrap resamples carry whole cells and all their doses; intervals are retrospective and unadjusted across comparisons.',
  'Readouts in the common comparison use named genes retained in outcome-informed bank43;57 all_gene_score uses all measured named genes.',
  'Program reconstruction replaces each named gene by its assigned hard program mean. Poor reconstruction is not proof that all34 coordinates cannot encode a footprint.',
  'Technical PRISM repeats support endpoint precision; they do not establish independent biological replication or matched molecular/functional exposure timing.',
  'Relative abundance is not a direct cell-death measurement.',
  'No capmatinib PRISM endpoint is available in these saved inputs; missing_function.tsv.gz lists missing contexts.'),
  file.path(root, 'CROSS_DRUG_METHODS.txt'))
writeLines(capture.output(sessionInfo()), file.path(root, 'cross_drug_session_info.txt'))
save_table(data.table(script = '/work/src/program_pilot/58_cross_drug_molecular_evidence.R',
  sha256 = digest(file = '/work/src/program_pilot/58_cross_drug_molecular_evidence.R', algo = 'sha256'),
  seed = 5802L, bootstrap_draws), 'script_provenance')
message('Cross-drug molecular information audit completed.')
