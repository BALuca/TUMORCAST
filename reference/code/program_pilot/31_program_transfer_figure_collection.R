#!/usr/bin/env Rscript
# Explain the frozen ERBB-family transfer experiment; no model fitting or discovery.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})
source('/work/src/program_pilot/response_models.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)

# Keep original experiment, endpoint, and all evaluation observations unchanged.
drug_name <- 'neratinib'
external_program_count <- 8L
parent <- '/work/results/program_pilot/20_multidrug_baseline_response_function_ERBB_holdout'
root <- '/work/results/program_pilot/31_program_transfer_figure_collection'
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
for (path in file.path(root, c('data', 'figures'))) dir.create(path, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time()), ' | ', stage)
  fwrite(data.table(stage, time = as.character(Sys.time()), pid = Sys.getpid()), file.path(root, 'run_status.tsv'), sep = '\t')
}
options(error = function() {status('FAILED'); traceback(5); quit(save = 'no', status = 1L)})
status('loading original family-excluded results')
settings <- qread(file.path(parent, 'settings.qs'))
prepared <- qread(file.path(parent, 'data/prepared_scores.qs'))
decoder <- qread(file.path(parent, 'models/functional_decoder.qs'))
bank <- qread(file.path(parent, 'models/program_bank.qs'))
predictions <- fread(file.path(parent, 'data/functional_predictions.tsv.gz'))
uncertainty <- fread(file.path(parent, 'data/functional_by_drug_uncertainty.tsv'))
panel <- qread(panel_path)
programs <- colnames(bank$weights)
lambda_grid <- settings$lambda_grid
stopifnot(length(programs) == 40L, !any(bank$training_ids %in% settings$held_ids),
  !any(decoder$training_ids %in% settings$held_ids),
  !any(decoder$training_drugs %in% settings$excluded_drugs),
  identical(bank$excluded_drugs, settings$excluded_drugs),
  identical(prepared$signature, settings$signature))
sources <- c(panel_path, file.path(parent, c('settings.qs', 'models/program_bank.qs',
  'models/functional_decoder.qs', 'data/prepared_scores.qs', 'data/functional_predictions.tsv.gz',
  'data/functional_by_drug_uncertainty.tsv', 'data/test_dose_provenance.tsv')),
  '/work/src/program_pilot/31_program_transfer_figure_collection.R',
  '/work/src/program_pilot/response_models.R', '/work/src/R/lib_heatmap.R')
hashes <- data.table(path = sources, sha256 = vapply(sources, digest::digest, character(1), file = TRUE, algo = 'sha256'))
save_table(hashes, 'input_hashes')
writeLines(system('git -C /work rev-parse HEAD', intern = TRUE), file.path(root, 'git_revision.txt'))

# All comparison names describe the original, saved prediction variants.
comparisons <- data.table(
  comparison = c('Full programs + baseline', 'Programs without baseline', 'Response magnitude + baseline',
    'Mean response + baseline', 'Swapped response + baseline'),
  model = c('measured_D_conditioned', 'measured_D', 'response_magnitude',
    'measured_D_conditioned', 'measured_D_conditioned'),
  response_input = c('measured', 'measured', 'measured', 'training_mean', 'swapped_recipient'),
  expected_r = c(.66093286926805, .531784907634374, .194615179948016, .00411201422258724, -.0743647298012749))
display <- merge(predictions, comparisons, by = c('model', 'response_input'))
display[, comparison := factor(comparison, levels = comparisons$comparison)]
display[, `:=`(observed_centered = observed - mean(observed),
  predicted_centered = predicted - mean(predicted)), by = .(comparison, drug, dose)]
metrics <- display[, .(cells = uniqueN(depmap_id), observations = .N,
  pearson = cor_safe(observed_centered, predicted_centered),
  predicted_measured_spread_ratio = sd(predicted_centered) / sd(observed_centered)),
  by = .(comparison, model, response_input, drug)]
metrics <- merge(metrics, uncertainty[, .(model, response_input, drug,
  ci_low = cell_bootstrap_lower, ci_high = cell_bootstrap_upper)], by = c('model', 'response_input', 'drug'))
headline <- metrics[drug == drug_name]
stopifnot(all(abs(headline$pearson - comparisons$expected_r[match(headline$comparison, comparisons$comparison)]) < 1e-12),
  all(headline$cells == 9L), all(headline$observations == 27L))
save_table(metrics, 'all_family_metrics')
save_table(headline, 'headline_numbers')
save_table(display, 'all_comparison_plot_values')
save_table(fread(file.path(parent, 'data/test_dose_provenance.tsv')), 'dose_provenance')

# Replay and decompose the zero-anchored original ridge: D_p * (beta_p + gamma_p * B_p).
query <- predictions[drug == drug_name & model == 'measured_D_conditioned' & response_input == 'measured']
setorder(query, dose, depmap_id)
d <- prepared$measured[match(query$observation_id, prepared$metadata$observation_id), , drop = FALSE]
b <- sweep(sweep(prepared$baseline[query$depmap_id, , drop = FALSE], 2L,
  decoder$baseline_mean, '-'), 2L, decoder$baseline_sd, '/')
fit <- decoder$models$measured_D_conditioned
coefficients <- as.numeric(coef(fit$fit, s = fit$lambda))
beta <- coefficients[2:41]
gamma <- coefficients[42:81]
shared <- sweep(d, 2L, beta, '*')
context <- sweep(d * b, 2L, gamma, '*')
contributions <- shared + context
replay_difference <- max(abs(rowSums(contributions) - query$predicted))
stopifnot(length(coefficients) == 81L, abs(coefficients[1L]) < 1e-12, replay_difference < 1e-10,
  max(abs(predict_ridge(fit, cbind(d, d * b)) - query$predicted)) < 1e-10)

# Rank display programs using development-cell contribution variation only.
labels <- as.data.table(panel$sensitivity)
training <- labels[depmap_id %in% decoder$training_ids & drug %in% decoder$training_drugs &
  is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
td <- prepared$measured[match(training$observation_id, prepared$metadata$observation_id), , drop = FALSE]
tb <- sweep(sweep(prepared$baseline[training$depmap_id, , drop = FALSE], 2L,
  decoder$baseline_mean, '-'), 2L, decoder$baseline_sd, '/')
tc <- sweep(td, 2L, beta, '*') + sweep(td * tb, 2L, gamma, '*')
training_contribution <- data.table(drug = training$drug, dose = training$dose, tc)
training_contribution[, (programs) := lapply(.SD, function(x) x - mean(x)), by = .(drug, dose), .SDcols = programs]
importance <- data.table(program = programs, training_contribution_sd = vapply(training_contribution[, ..programs], sd, numeric(1)),
  genes = as.integer(colSums(bank$weights > 0)), beta = beta, baseline_modifier = gamma)
setorder(importance, -training_contribution_sd, program)
importance[, display_rank := seq_len(.N)]
top_programs <- importance$program[seq_len(external_program_count)]

# Annotate this frozen bank itself; names from the exploratory bank are not inherited.
status('annotating frozen programs and exporting molecular tables')
hallmark <- as.data.table(msigdbr::msigdbr(species = 'Homo sapiens', category = 'H'))
pathways <- lapply(split(hallmark$gene_symbol, hallmark$gs_name), function(x) intersect(unique(x), bank$genes))
pathways <- pathways[lengths(pathways) >= 5L]
enrichment <- rbindlist(lapply(programs, function(p) {
  genes <- bank$genes[bank$partition == match(p, programs)]
  z <- as.data.table(fgsea::fora(pathways, genes, bank$genes, minSize = 5L))
  z[, `:=`(program = p, overlap_genes = vapply(overlapGenes, paste, character(1), collapse = ';'))]
  z[, overlapGenes := NULL]
  z
}))
enrichment[, q_across_programs := p.adjust(pval, 'BH')]
save_table(enrichment, 'frozen_program_hallmark_enrichment')
annotation <- rbindlist(lapply(programs, function(p) {
  hits <- enrichment[program == p & overlap >= 3L & q_across_programs <= .05][order(q_across_programs, pval)]
  term <- if (nrow(hits)) sub('^HALLMARK_', '', hits$pathway[1L]) else 'No supported Hallmark label'
  genes <- bank$genes[bank$partition == match(p, programs)]
  data.table(program = p, hallmark = gsub('_', ' ', term), genes = length(genes),
    members = paste(genes, collapse = ';'), supported_terms = nrow(hits))
}))
save_table(annotation, 'frozen_program_annotations')
importance <- merge(importance, annotation[, .(program, hallmark)], by = 'program', sort = FALSE)
setorder(importance, display_rank)
save_table(importance, 'training_program_display_selection')
save_table(data.table(resource = 'Hallmark', package = 'msigdbr', version = as.character(packageVersion('msigdbr')),
  sha256 = digest::digest(pathways, algo = 'sha256'), provenance = 'existing package gene sets; no download',
  license_note = 'Existing msigdbr/fgsea dependencies; MSigDB source terms retained'), 'annotation_resource')
save_table(data.table(gene = bank$genes, program = programs[bank$partition]), 'gene_memberships')

# Cell and program ordering use molecular structure, keeping all annotation identities aligned.
td_sd <- apply(td, 2L, sd)
stopifnot(all(td_sd > 0))
program_order <- programs[hclust(as.dist(sqrt(pmax(2 * (1 - cor(td)), 0))), method = 'ward.D2')$order]
cell_ids <- sort(unique(query$depmap_id))
cell_profiles <- do.call(rbind, lapply(cell_ids, function(id) colMeans(d[query$depmap_id == id, , drop = FALSE])))
cell_order <- cell_ids[hclust(dist(sweep(cell_profiles, 2L, td_sd, '/')), method = 'ward.D2')$order]
column_order <- unlist(lapply(sort(unique(query$dose)), function(dose_value)
  match(cell_order, query[dose == dose_value, depmap_id]) + which(query$dose == dose_value)[1L] - 1L))
stopifnot(setequal(column_order, seq_len(nrow(query))))
query_order <- query[column_order]
ordered_ids <- query_order$observation_id
row_labels <- setNames(paste0(programs, ' (', colSums(bank$weights > 0), ' genes)'), programs)
top_annotation <- annotation[match(top_programs, program)]
top_labels <- setNames(ifelse(top_annotation$supported_terms > 0L,
  paste0(top_programs, ' | ', tools::toTitleCase(tolower(top_annotation$hallmark))),
  paste0(top_programs, ' (', top_annotation$genes, ' genes)')), top_programs)
mean_d <- t(vapply(seq_len(nrow(query)), function(i) {
  ix <- prepared$metadata$drug == drug_name & prepared$metadata$dose == query$dose[i] &
    prepared$metadata$depmap_id %in% bank$training_ids
  stopifnot(any(ix))
  colMeans(prepared$measured[ix, , drop = FALSE])
}, numeric(40)))
deviation <- sweep(d - mean_d, 2L, td_sd, '/')
program_table <- rbindlist(lapply(seq_along(programs), function(i) data.table(
  observation_id = query$observation_id, cell = query$cell, dose = query$dose, program = programs[i],
  displacement = d[, i], training_mean_displacement = mean_d[, i],
  displacement_deviation_training_sd = deviation[, i], baseline_training_z = b[, i],
  shared_contribution_log2 = shared[, i], baseline_contribution_log2 = context[, i],
  total_contribution_log2 = contributions[, i], observed_pct = query$observed_pct)))
save_table(program_table, 'program_response_and_contribution_values')
save_table(query_order[, .(observation_id, cell, dose)], 'heatmap_column_order')
save_table(data.table(program = program_order), 'heatmap_row_order')

verification <- data.table(check = c('quoted five correlations', 'ridge and contribution replay',
  'held recipients excluded from discovery', 'held recipients excluded from functional fitting',
  'ERBB family excluded from functional fitting', 'external programs selected on training only'),
  detail = c('all discrepancies <1e-12', sprintf('maximum discrepancy %.3g', replay_difference),
    'zero overlaps', 'zero overlaps', 'zero overlaps', 'SD of training program contributions within drug/dose'))
save_table(verification, 'verification')

# Reused export helpers checkpoint each completed PDF, PNG and its caption.
figures <- list()
register_figure <- function(stem, kind, object, width, height, caption) {
  figures[[stem]] <<- list(kind = kind, object = object, width = width, height = height, caption = caption)
  index <- rbindlist(lapply(names(figures), function(name) data.table(figure = name,
    audience = ifelse(grepl('^external', name), 'external', 'internal'),
    pdf = file.path(root, 'figures', paste0(name, '.pdf')), png = file.path(root, 'figures', paste0(name, '.png')),
    caption = figures[[name]]$caption)))
  fwrite(index, file.path(root, 'FIGURE_INDEX.tsv'), sep = '\t')
}
save_plot <- function(p, stem, width, height, caption) {
  ggsave(file.path(root, 'figures', paste0(stem, '.pdf')), p, width = width, height = height,
    device = grDevices::cairo_pdf, bg = 'white')
  ggsave(file.path(root, 'figures', paste0(stem, '.png')), p, width = width, height = height, dpi = 180, bg = 'white')
  register_figure(stem, 'plot', p, width, height, caption)
}
save_map <- function(h, stem, width, height, caption) {
  grDevices::cairo_pdf(file.path(root, 'figures', paste0(stem, '.pdf')), width = width, height = height)
  ComplexHeatmap::draw(h, heatmap_legend_side = 'right', annotation_legend_side = 'right')
  dev.off()
  png(file.path(root, 'figures', paste0(stem, '.png')), width = width, height = height, units = 'in', res = 180)
  ComplexHeatmap::draw(h, heatmap_legend_side = 'right', annotation_legend_side = 'right')
  dev.off()
  register_figure(stem, 'heatmap', h, width, height, caption)
}
plot_theme <- theme_classic(base_size = 10) + theme(aspect.ratio = 1,
  plot.title = element_text(size = 11, face = 'bold'), legend.position = 'bottom')
cell_colors <- setNames(c('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00',
  '#56B4E9', '#6A3D9A', '#666666', '#B15928'), sort(unique(query$cell)))
method_colors <- setNames(c('#0072B2', '#009E73', '#D55E00', '#777777', '#CC79A7'), comparisons$comparison)
caption_scope <- paste('Nine excluded cell lines; 27 cell-dose observations, not 27 independent recipients.',
  'ERBB-family data excluded from program discovery and functional fitting.',
  'PRISM values bracket-interpolated to 0.05, 0.5 and 5 µM.',
  'Retrospective illustrative selection; other family results are shown separately.')

status('exporting external comparison figures')
headline[, comparison := factor(comparison, levels = rev(comparisons$comparison))]
control_plot <- ggplot(headline, aes(pearson, comparison, color = comparison)) +
  geom_vline(xintercept = 0, color = 'grey75') +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = comparison), linewidth = .55) +
  geom_point(size = 2.6) + geom_text(aes(label = sprintf(' %.2f', pearson)), hjust = -.2, vjust = -1, size = 3) +
  scale_color_manual(values = method_colors, guide = 'none') +
  scale_x_continuous(limits = c(-.6, .95), breaks = c(-.5, 0, .5)) +
  theme_classic(base_size = 10) + labs(x = 'Across-cell Pearson r at fixed drug/dose', y = NULL,
    title = 'Specific response patterns support transfer', subtitle = 'Neratinib | recipient-bootstrap 95% intervals')
save_plot(control_plot, 'external_01_program_pattern_controls', 7.5, 3.4,
  paste('Exact original log2-scale correlations after centering observed and predicted abundance within each dose.',
    'Mean and swapped responses enter the identical baseline-conditioned decoder.',
    'Magnitude is a separately fitted comparator with baseline modulation.',
    'Recipient swap is one saved shuffle, not a permutation test.', caption_scope))

full_rows <- display[drug == drug_name & comparison == comparisons$comparison[1L]]
full_metric <- headline[comparison == comparisons$comparison[1L]]
scatter <- ggplot(full_rows, aes(observed_centered, predicted_centered)) +
  geom_hline(yintercept = 0, color = 'grey90') + geom_vline(xintercept = 0, color = 'grey90') +
  geom_point(aes(color = cell, shape = factor(dose)), size = 2.3) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#333333', linewidth = .55) +
  scale_color_manual(values = cell_colors) + plot_theme +
  labs(x = 'Measured log2 abundance\n(deviation from dose mean)',
    y = 'Predicted log2 abundance\n(deviation from dose mean)', color = NULL, shape = 'Dose (µM)',
    title = 'Transfer to excluded neratinib and recipients',
    subtitle = sprintf('r = %.2f; 95%% interval %.2f–%.2f | response scale is compressed',
      full_metric$pearson, full_metric$ci_low, full_metric$ci_high)) +
  guides(color = guide_legend(nrow = 3), shape = guide_legend(nrow = 1))
save_plot(scatter, 'external_02_held_recipient_transfer', 5.8, 6,
  paste('Both axes are centered separately at each dose; this measures across-cell response ordering.',
    'Independent axis scales expose the smaller predicted spread.', caption_scope))

# Native units show what the correlation cannot establish about response calibration.
high <- query[dose == 5]
high[, cell := factor(cell, levels = rev(high[order(observed_pct), cell]))]
native_points <- melt(high, id.vars = 'cell', measure.vars = c('observed_pct', 'predicted_pct'),
  variable.name = 'series', value.name = 'abundance')
native_points[, series := factor(series, levels = c('observed_pct', 'predicted_pct'),
  labels = c('Measured PRISM', 'Predicted from programs'))]
native_plot <- ggplot(high, aes(y = cell)) +
  geom_vline(xintercept = 50, linetype = 2, color = 'grey70') +
  geom_segment(aes(x = observed_pct, xend = predicted_pct, yend = cell), color = 'grey65') +
  geom_point(data = native_points, aes(x = abundance, color = series, shape = series), size = 2.5) +
  scale_color_manual(values = c('Measured PRISM' = '#0072B2', 'Predicted from programs' = '#D55E00')) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom') +
  labs(x = 'Abundance (% vehicle)', y = NULL, color = NULL, shape = NULL,
    title = 'Neratinib at 5 µM: ordering and absolute scale',
    subtitle = 'Observed range: 30–182%; predicted: 70–103%')
save_plot(native_plot, 'external_03_native_abundance_calibration', 6, 4.2,
  paste('Each row is one held recipient at the same dose. Dashed line: 50% vehicle.',
    'Two observed sensitive recipients are both missed; high correlation does not establish sensitivity calls.', caption_scope))

# Continuous abundance strips and program matrices share exactly the same column order.
make_map <- function(mat, selected_programs, title, hm_name, labels, limit) {
  ordered <- t(mat[column_order, match(selected_programs, programs), drop = FALSE])
  colnames(ordered) <- ordered_ids
  rownames(ordered) <- selected_programs
  annotation_df <- data.frame(Cell = query_order$cell, row.names = ordered_ids)
  stopifnot(identical(colnames(ordered), rownames(annotation_df)))
  bottom <- HeatmapAnnotation(
    `Measured %vehicle` = query_order$observed_pct,
    `Predicted %vehicle` = query_order$predicted_pct,
    col = list(`Measured %vehicle` = circlize::colorRamp2(c(30, 100, 182), c('#2166AC', 'white', '#B2182B')),
      `Predicted %vehicle` = circlize::colorRamp2(c(30, 100, 182), c('#2166AC', 'white', '#B2182B'))),
    annotation_name_gp = grid::gpar(fontsize = 8), simple_anno_size = grid::unit(3, 'mm'))
  draw_heatmap(ordered, col_anno_df = annotation_df, col_anno_cols = 'Cell',
    palette_overrides = list(Cell = cell_colors), scale_rows = FALSE,
    ann_top_args = list(Dose = ComplexHeatmap::anno_block(
      labels = c('0.05 µM', '0.5 µM', '5 µM'), gp = grid::gpar(fill = 'white', col = NA),
      labels_gp = grid::gpar(fontsize = 9), height = grid::unit(5, 'mm'))),
    value_range = c(-limit, 0, limit), value_colors = c('#2166AC', 'white', '#B2182B'),
    hm_name = hm_name, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_labels = labels[selected_programs], column_labels = query_order$cell,
    column_split = factor(query_order$dose, levels = c(.05, .5, 5), labels = c('0.05 µM', '0.5 µM', '5 µM')),
    cluster_column_slices = FALSE, bottom_annotation = bottom,
    row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 7),
    column_names_rot = 90, title = title,
    column_title_gp = grid::gpar(fontsize = 10), heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8)))
}
top_order <- program_order[program_order %in% top_programs]
h_top <- make_map(deviation, top_order, 'Neratinib: program pattern differs between recipients',
  'D deviation\n(training SD)', top_labels, max(abs(deviation[, match(top_programs, programs)])))
save_map(h_top, 'external_04_program_patterns', 10.5, 4.9,
  paste('Eight programs chosen by within-treatment functional contribution variation in training data only.',
    'Tiles: measured D minus the training-recipient mean for this drug/dose, divided by training D SD.',
    'Rows and common cell order use molecular structure; no sensitivity sorting.',
    'Hallmark annotations describe this frozen bank, not causal drivers or pathway activity.', caption_scope))

family_plot <- ggplot(metrics[comparison %in% comparisons$comparison[1:3]],
  aes(pearson, drug, color = comparison)) +
  geom_vline(xintercept = 0, color = 'grey75') +
  geom_point(position = position_dodge(width = .5), size = 2.1) +
  scale_color_manual(values = method_colors) + scale_x_continuous(limits = c(-.6, .9)) +
  theme_classic(base_size = 10) + theme(legend.position = 'bottom') +
  guides(color = guide_legend(ncol = 1)) +
  labs(x = 'Across-cell Pearson r at fixed drug/dose', y = NULL, color = NULL,
    title = 'Full ERBB-family evaluation: transfer is uneven',
    subtitle = 'All seven supported drugs; same nine held recipients and three doses')
save_plot(family_plot, 'external_05_full_family_context', 7.2, 4.8,
  paste('All supported excluded-family drugs are retained; neratinib was selected retrospectively as an illustration.',
    'Broad family transfer is not established. Each point uses centered log2 abundance.', caption_scope))

status('exporting detailed internal figures')
design_boxes <- data.table(
  x = c(1, 4, 7, 1, 4, 7), y = c(3, 3, 3, 1.5, 1.5, 1.5),
  label = c('DISCOVER PROGRAMS\n30 development cells\nERBB family excluded',
    'FREEZE REPRESENTATION\n40 hard gene clusters\n5,109 genes',
    'FIT SHARED FUNCTION\n26 cells × 90 drugs\n6,260 observations',
    'HELD INPUTS\n9 excluded recipients\nMeasured D + untreated B',
    'FROZEN FUNCTION\nSum over 40 programs:\nD × (beta + gamma × B)',
    'HELD EVALUATION\nNeratinib, 3 doses\nIndependent PRISM assay'))
design_plot <- ggplot(design_boxes) +
  geom_rect(aes(xmin = x - 1.3, xmax = x + 1.3, ymin = y - .45, ymax = y + .45),
    fill = '#F0F5F8', color = '#0072B2', linewidth = .5) +
  geom_text(aes(x, y, label = label), size = 3.4, lineheight = 1.15) +
  annotate('segment', x = c(2.3, 5.3, 2.3, 5.3), xend = c(2.7, 5.7, 2.7, 5.7),
    y = c(3, 3, 1.5, 1.5), yend = c(3, 3, 1.5, 1.5),
    arrow = grid::arrow(length = grid::unit(2, 'mm'))) +
  annotate('text', x = 4, y = .45, size = 3.2,
    label = paste('Readout: correlation across cells after removing each dose mean from log2 abundance.',
      'Controls: omit baseline, use response magnitude, use mean D, or swap D between recipients.',
      'Baseline modifies measured D here; baseline → D and combination prediction are not evaluated.', sep = '\n')) +
  scale_x_continuous(limits = c(-.4, 8.4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3.65), expand = c(0, 0)) + theme_void() +
  labs(title = 'What was learned, what was withheld, and what r = 0.66 measures',
    subtitle = 'The original ERBB-family excluded experiment; this collection changes reporting only') +
  theme(plot.title = element_text(size = 12, face = 'bold'), plot.subtitle = element_text(size = 10))
save_plot(design_plot, 'internal_00_experiment_and_readout', 10, 4.5,
  paste('The outcome-informed program bank excludes held recipients and the entire ERBB family.',
    'Baseline z-scores use functional training recipients. Functional coefficients are shared across training drugs.', caption_scope))
internal_scatter <- ggplot(display[drug == drug_name], aes(observed_centered, predicted_centered)) +
  geom_point(aes(color = cell, shape = factor(dose)), size = 1.7) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#333333', linewidth = .5) +
  geom_text(data = headline, aes(x = -Inf, y = Inf, label = sprintf('r = %.3f', pearson)),
    inherit.aes = FALSE, hjust = -.08, vjust = 1.2, size = 3) +
  facet_wrap(~comparison, scales = 'free', ncol = 3) + scale_color_manual(values = cell_colors) +
  plot_theme + labs(x = 'Measured log2 abundance minus dose mean', y = 'Predicted log2 abundance minus dose mean',
    title = 'Every quoted comparison, on exactly the same 27 observations', color = NULL, shape = 'Dose (µM)')
save_plot(internal_scatter, 'internal_01_all_control_scatterplots', 11, 8,
  paste('Independent x/y scales in every square panel; fit lines show association, not identity.', caption_scope))

h_all <- make_map(deviation, program_order, 'All 40 program responses | training-reference deviations',
  'D deviation\n(training SD)', row_labels, max(abs(deviation)))
save_map(h_all, 'internal_02_all_program_responses', 10.5, 9.5,
  paste('Complete frozen bank, including unsupported pathway labels. Same cell order at all doses.',
    'See source table for raw D, reference mean and scaling; no test-label feature selection.', caption_scope))
h_raw <- make_map(d, program_order, 'All 40 program responses | original signed displacement',
  'Mean gene D', row_labels, max(abs(d)))
save_map(h_raw, 'internal_03_raw_program_displacements', 10.5, 9.5,
  paste('Hard means of signed corrected expression displacement, without heatmap row scaling.', caption_scope))
h_contribution <- make_map(contributions, program_order, 'Functional contributions | sum equals the frozen prediction',
  'Log2 effect', row_labels, max(abs(contributions)))
save_map(h_contribution, 'internal_04_program_functional_contributions', 10.5, 9.5,
  paste('Per-program contribution D_p * (beta_p + gamma_p * baseline_z_p). Positive increases predicted abundance.',
    'These are fitted prediction contributions, not causal gene effects.', caption_scope))
h_baseline <- make_map(context, program_order, 'Contribution added by recipient baseline conditioning',
  'Log2 effect', row_labels, max(abs(context)))
save_map(h_baseline, 'internal_05_baseline_modulation_contributions', 10.5, 9.5,
  paste('Tiles are D_p * gamma_p * baseline_z_p. Background modifies measured responses;',
    'this does not predict those responses from baseline.', caption_scope))

curve <- display[drug == drug_name & comparison %in% comparisons$comparison[c(1, 2, 3, 4)],
  .(cell, dose, series = as.character(comparison), abundance = predicted_pct)]
curve <- rbind(curve, query[, .(cell, dose, series = 'Measured PRISM', abundance = observed_pct)])
curve_colors <- c(method_colors, 'Measured PRISM' = '#111111')
curve_plot <- ggplot(curve, aes(dose, abundance, color = series, group = series)) +
  geom_hline(yintercept = 50, linetype = 2, color = 'grey80') +
  geom_line(linewidth = .55) + geom_point(size = 1.4) +
  facet_wrap(~cell, ncol = 3) + scale_x_log10(breaks = c(.05, .5, 5)) +
  scale_color_manual(values = curve_colors) + plot_theme + guides(color = guide_legend(ncol = 2)) +
  labs(x = 'Dose (µM)', y = 'Abundance (% vehicle)', color = NULL,
    title = 'Neratinib dose curves: all nine recipients and controls')
save_plot(curve_plot, 'internal_06_all_recipient_dose_curves', 10, 10,
  paste('All facets share native abundance scale; dashed line: 50% vehicle.', caption_scope))

dose_metrics <- display[drug == drug_name, .(pearson_log2 = cor_safe(observed, predicted),
  pearson_pct = cor_safe(observed_pct, predicted_pct),
  sensitive = sum(observed_pct <= 50), sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50),
  measured_min_pct = min(observed_pct), measured_max_pct = max(observed_pct),
  predicted_min_pct = min(predicted_pct), predicted_max_pct = max(predicted_pct)), by = .(comparison, dose)]
save_table(dose_metrics, 'per_dose_native_metrics')
dose_plot <- ggplot(dose_metrics, aes(factor(dose), pearson_log2, color = comparison, group = comparison)) +
  geom_hline(yintercept = 0, color = 'grey80') + geom_line(linewidth = .5) + geom_point(size = 2) +
  scale_color_manual(values = method_colors) + plot_theme + guides(color = guide_legend(ncol = 1)) +
  labs(x = 'Dose (µM)', y = 'Across-cell Pearson r (log2 abundance)', color = NULL,
    title = 'Ordering at each individual dose | nine cells per point')
save_plot(dose_plot, 'internal_07_per_dose_comparisons', 7, 6,
  paste('Correlations within single doses; these differ from the headline pooled centered correlation.', caption_scope))

leave_one_out <- rbindlist(lapply(unique(full_rows$cell), function(removed) {
  z <- copy(full_rows[cell != removed])
  z[, `:=`(oc = observed - mean(observed), pc = predicted - mean(predicted)), by = dose]
  data.table(removed_cell = removed, pearson = cor_safe(z$oc, z$pc))
}))
save_table(leave_one_out, 'recipient_influence')
influence_plot <- ggplot(leave_one_out, aes(pearson, reorder(removed_cell, pearson))) +
  geom_vline(xintercept = full_metric$pearson, color = '#0072B2', linetype = 2) +
  geom_point(size = 2) + theme_classic(base_size = 10) +
  labs(x = 'Correlation after removing this recipient', y = 'Removed recipient',
    title = 'Does one cell explain the neratinib result?', subtitle = 'Dose means recomputed after each removal; frozen predictions')
save_plot(influence_plot, 'internal_08_recipient_influence', 6.2, 4.2,
  paste('Leave-one-recipient-out descriptive sensitivity analysis; no models are refitted.', caption_scope))

single_program <- program_table[program %in% top_programs[1:4]]
single_program[, dose_label := paste(dose, 'µM')]
single_metrics <- single_program[, .(label = sprintf('r = %.2f', cor_safe(displacement, observed_pct))), by = .(program, dose_label)]
save_table(single_metrics, 'display_program_associations')
program_plot <- ggplot(single_program, aes(displacement, observed_pct)) +
  geom_point(aes(color = cell), size = 1.6) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#333333', linewidth = .45) +
  geom_text(data = single_metrics, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.05, vjust = 1.2, size = 2.8) +
  facet_wrap(vars(program, dose_label), ncol = 3, scales = 'free') +
  scale_color_manual(values = cell_colors) + plot_theme +
  labs(x = 'Measured program displacement', y = 'Measured PRISM abundance (% vehicle)', color = NULL,
    title = 'Individual program associations | display selection uses training data')
save_plot(program_plot, 'internal_09_selected_program_associations', 10, 12,
  paste('Four programs selected by training contribution rank, not held association.',
    'Univariate exploratory associations do not reproduce or explain the full multivariate prediction alone.', caption_scope))

# An external overview joins evidence and controls while retaining native-unit calibration.
overview <- (control_plot + native_plot) / scatter + plot_layout(heights = c(1, 1.4)) +
  plot_annotation(title = 'Program patterns carry functional information across excluded interventions',
    subtitle = 'Neratinib illustrates recipient-specific transfer; absolute response remains under-calibrated.',
    caption = 'Retrospective example. ERBB family and nine recipients excluded from discovery/fitting.\nMeasured molecular inputs; baseline→response and combination prediction are not tested here.',
    theme = theme(plot.title = element_text(size = 14, face = 'bold'), plot.caption = element_text(size = 9)))
save_plot(overview, 'external_00_transfer_overview', 10, 11,
  paste('Overview of the specific-pattern controls, native-scale limitation and held transfer association.', caption_scope))

# Produce convenient collection PDFs in addition to compact individual exports.
status('assembling internal and external collections')
external_names <- c('external_00_transfer_overview', 'external_04_program_patterns',
  'external_05_full_family_context')
internal_names <- c('internal_00_experiment_and_readout', 'internal_01_all_control_scatterplots', 'internal_02_all_program_responses',
  'internal_03_raw_program_displacements', 'internal_04_program_functional_contributions',
  'internal_05_baseline_modulation_contributions', 'internal_06_all_recipient_dose_curves',
  'internal_07_per_dose_comparisons', 'internal_08_recipient_influence',
  'internal_09_selected_program_associations', 'external_05_full_family_context')
for (audience in c('external', 'internal')) {
  selected_figures <- if (audience == 'external') external_names else internal_names
  grDevices::cairo_pdf(file.path(root, paste0(audience, '_figure_collection.pdf')), width = 11, height = 12, onefile = TRUE)
  for (name in selected_figures) {
    figure <- figures[[name]]
    if (figure$kind == 'plot') print(figure$object) else
      ComplexHeatmap::draw(figure$object, heatmap_legend_side = 'right', annotation_legend_side = 'right')
  }
  dev.off()
}
qsave(list(family = 'ERBB', example_drug = drug_name, external_programs = top_programs,
  excluded_drugs = settings$excluded_drugs, excluded_recipients = settings$held_ids,
  external_collection_order = external_names, internal_collection_order = internal_names,
  hashes = hashes, retrospective = TRUE), file.path(root, 'settings.qs'))
for (audience in c('external', 'internal')) {
  names_for_caption <- names(figures)[grepl(paste0('^', audience), names(figures))]
  caption_text <- unlist(lapply(names_for_caption, function(name)
    c(name, figures[[name]]$caption, '')), use.names = FALSE)
  writeLines(caption_text, file.path(root, paste0(audience, '_figure_captions.txt')))
}
writeLines(capture.output(sessionInfo()), file.path(root, 'sessionInfo.txt'))
writeLines(c('PROGRAM TRANSFER FIGURE COLLECTION', '',
  'Start with external_figure_collection.pdf (3 pages). Internal collection contains 11 pages.',
  'Each individual figure has a PDF, PNG and caption in FIGURE_INDEX.tsv.', '',
  'WHAT WAS TRAINED: original 40 hard-mean programs; shared zero-anchored ridge effect D*(beta+gamma*baseline_z).',
  'Development function data: 6,260 observations, 26 cells, 90 drugs; ERBB family excluded.',
  'WHAT WAS EVALUATED: measured neratinib responses in nine excluded cells at three doses.',
  'HEADLINE METRIC: log2 abundance centered separately within each dose, across recipients.',
  'PROGRAM MAPS: original signed D or training-reference deviation; prediction contributions are labeled separately.',
  'The eight external display programs are selected on training contribution variation only.',
  'Rows and cells are ordered by molecular structure, not sensitivity.', '',
  'SUPPORTED: a retrospective example of program-specific response-to-function transfer.',
  'LIMITS: family-wide results are mixed; amplitude is compressed; sensitive calls are not recovered at 5 µM.',
  'Baseline conditions measured D here. Baseline->D, combinations, and HGF rescue are not established by these figures.',
  'The one recipient swap is descriptive; bootstrap intervals resample recipients with doses kept together.',
  'Program discovery uses development functional outcomes; it is not an outcome-free molecular discovery model.',
  'Hallmark enrichment is an annotation of gene composition, not causal mechanism or pathway activation.', '',
  'EXACT HEADLINE NUMBERS', capture.output(print(headline[, .(comparison, pearson, ci_low, ci_high)])), '',
  'VERIFICATION', capture.output(print(verification))), file.path(root, 'READING_GUIDE.txt'))
status('COMPLETE')
message('Figure collection: ', root)
