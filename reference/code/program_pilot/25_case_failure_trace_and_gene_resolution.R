#!/usr/bin/env Rscript
# Trace frozen decoder failures through labels, neighborhoods and gene resolution.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
})
source('/work/src/program_pilot/nonlinear_response.R')
source('/work/src/program_pilot/response_models.R')
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
set.seed(20261010L)
root <- '/work/results/program_pilot/25_case_failure_trace_and_gene_resolution'
for (p in file.path(root, c('data', 'models', 'figures'))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
status <- function(stage) {
  message(format(Sys.time(), tz = 'Europe/Zurich', usetz = TRUE), ' | ', stage)
  fwrite(data.table(stage = stage, time = as.character(Sys.time())), file.path(root, 'run_status.tsv'), sep = '\t')
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 140)
}
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
lineage_path <- '/work/outputs/dataset_inventory/cache/depmap/Model.csv'
lineage <- fread(lineage_path, select = c('ModelID', 'CellLineName', 'OncotreeLineage', 'OncotreePrimaryDisease'))
cases <- data.table(case = c('CFPAC low', 'SW900 low', 'H2347 low', 'BT474 low', 'CFPAC middle', 'Hs766T docetaxel'),
  recipient = c('CFPAC-1', 'SW 900', 'NCI-H2347', 'BT-474', 'CFPAC-1', 'Hs 766T'),
  drug = c(rep('ixazomib', 5L), 'docetaxel'), dose = c(rep(.05, 4L), .5, .5))
cases <- merge(cases, labels, by = c('recipient', 'drug', 'dose'), all.x = TRUE, sort = FALSE)
stopifnot(nrow(cases) == 6L, !anyNA(cases$observation_id))
cases[, `:=`(prism_lower_abundance = 100 * 2^prism_lower_log2,
  prism_upper_abundance = 100 * 2^prism_upper_log2)]
save_table(cases, 'case_labels')
ix <- labels[drug == 'ixazomib' & dose == .05 & is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
stopifnot(nrow(ix) == 35L, !anyDuplicated(ix$depmap_id))

# Read each molecular file once; retain only the targeted signatures and controls.
status('loading targeted gene profiles and vehicle contrasts')
profiles <- list()
input_paths <- panel$input_paths[match(ix$depmap_id, panel$cell_map$depmap_id)]
for (i in seq_len(nrow(ix))) {
  id <- ix$depmap_id[i]
  m <- qread(input_paths[i])
  obs <- unique(c(ix$observation_id[i], cases[depmap_id == id, observation_id]))
  stopifnot(all(obs %in% rownames(m$response)), identical(colnames(m$response), m$genes))
  profiles[[id]] <- list(response = m$response[obs, , drop = FALSE], null = m$null[[id]],
    baseline = setNames(m$baseline[[id]], m$genes), metadata = m$metadata[match(obs, observation_id)])
}
genes <- colnames(profiles[[1L]]$response)
stopifnot(all(vapply(profiles, function(p) identical(colnames(p$response), genes), logical(1))))
gene_response <- do.call(rbind, lapply(seq_len(nrow(ix)), function(i)
  profiles[[ix$depmap_id[i]]]$response[ix$observation_id[i], , drop = FALSE]))
gene_baseline <- do.call(rbind, lapply(ix$depmap_id, function(id) profiles[[id]]$baseline))
dimnames(gene_baseline) <- list(ix$depmap_id, genes)
save_table(rbindlist(lapply(profiles, function(p) p$metadata)), 'molecular_metadata')

neighbors <- kernel_summary <- noise_predictions <- noise_summary <- representation_predictions <- lineage_tables <- list()
representation_selection <- gene_associations <- program_values <- pair_summaries <- gene_differences <- list()
source_paths <- c(panel_path, lineage_path, input_paths, '/work/src/program_pilot/25_case_failure_trace_and_gene_resolution.R',
  '/work/src/program_pilot/nonlinear_response.R', '/work/src/program_pilot/response_models.R', '/work/src/R/lib_heatmap.R')
for (family in c('MEK', 'ERBB')) {
  status(paste(family, 'exact frozen-model trace'))
  parent <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_', family, '_holdout')
  capacity_root <- paste0('/work/results/program_pilot/24_ixazomib_decoder_capacity_and_molecular_signal_', family, '_holdout')
  settings <- qread(file.path(parent, 'settings.qs'))
  prepared <- qread(file.path(parent, 'data/prepared_scores.qs'))
  bank <- qread(file.path(parent, 'models/program_bank.qs'))
  original <- qread(file.path(parent, 'models/nonlinear_RBF_SVR_D.qs'))
  cell_lineages <- copy(lineage[ModelID %in% c(original$training_ids, settings$held_ids)])
  cell_lineages[, `:=`(family = family, cohort = ifelse(ModelID %in% original$training_ids, 'Fitting cells', 'Held cells'))]
  lineage_tables[[family]] <- cell_lineages
  fits <- list(Original = original$model,
    Local = qread(file.path(capacity_root, 'models/capacity_2.qs'))$model,
    Very_local = qread(file.path(capacity_root, 'models/capacity_3.qs'))$model)
  source_paths <- c(source_paths, file.path(parent, c('settings.qs', 'data/prepared_scores.qs',
    'models/program_bank.qs', 'models/nonlinear_RBF_SVR_D.qs')),
    file.path(capacity_root, c('models/capacity_2.qs', 'models/capacity_3.qs')))
  train <- labels[depmap_id %in% original$training_ids & !drug %in% settings$excluded_drugs &
    is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
  d <- prepared$measured[match(train$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  query_d <- prepared$measured[match(cases$observation_id, prepared$metadata$observation_id), , drop = FALSE]
  stopifnot(all(cases$depmap_id %in% settings$held_ids), !length(intersect(settings$held_ids, original$training_ids)))
  for (model_name in names(fits)) {
    f <- fits[[model_name]]
    z <- sweep(sweep(d[, f$active, drop = FALSE], 2L, f$x_mean[f$active], '-'), 2L, f$x_sd[f$active], '/')
    stopifnot(max(abs(z[f$model$index, , drop = FALSE] - f$model$SV)) < 1e-10)
    for (k in seq_len(nrow(cases))) {
      cc <- cases[k]
      q <- query_d[k, , drop = FALSE]
      qz <- as.numeric((q[, f$active] - f$x_mean[f$active]) / f$x_sd[f$active])
      squared <- rowSums(sweep(z, 2L, qz, '-')^2)
      kernels <- exp(-f$model$gamma * squared)
      contribution <- kernels[f$model$index] * as.numeric(f$model$coefs) * f$y_sd
      intercept <- f$y_mean - f$model$rho * f$y_sd
      direct <- intercept + sum(contribution)
      predicted <- predict_program_svr(f, q)
      stopifnot(abs(direct - predicted) < 1e-8)
      kernel_summary[[paste(family, model_name, k)]] <- data.table(family = family, model = model_name, case = cc$case,
        observed_pct = cc$prism_abundance_pct, predicted_pct = 100 * 2^predicted,
        no_neighbor_prediction_pct = 100 * 2^intercept, maximum_kernel_similarity = max(kernels),
        neighbors_above_05 = sum(kernels > .5), positive_log2_sum = sum(contribution[contribution > 0]),
        negative_log2_sum = sum(contribution[contribution < 0]), replay_difference = direct - predicted)
      if (model_name == 'Original') {
        for (scope in c('Same drug and dose', 'All training treatments')) {
          eligible <- if (scope == 'Same drug and dose') which(train$drug == cc$drug & train$dose == cc$dose) else seq_len(nrow(train))
          nearest <- eligible[order(squared[eligible])][1:5]
          nn <- copy(train[nearest, .(observation_id, recipient, drug, dose, prism_abundance_pct)])
          nn[, `:=`(family = family, case = cc$case, scope = scope, rank = seq_len(.N),
            distance = sqrt(squared[nearest] / ncol(z)), query_abundance_pct = cc$prism_abundance_pct)]
          neighbors[[paste(family, k, scope)]] <- nn
        }
      }
      null <- as.matrix(profiles[[cc$depmap_id]]$null[, bank$genes, drop = FALSE] %*% bank$weights)
      null <- sweep(null, 2L, colMeans(null), '-')
      perturbed <- sweep(null, 2L, as.numeric(q), '+')
      pp <- 100 * 2^predict_program_svr(f, perturbed)
      noise_predictions[[paste(family, model_name, k)]] <- data.table(family = family, model = model_name,
        case = cc$case, vehicle_contrast = rownames(null), predicted_pct = pp)
      noise_summary[[paste(family, model_name, k)]] <- data.table(family = family, model = model_name, case = cc$case,
        original_predicted_pct = 100 * 2^predicted, vehicle_perturbed_median_pct = median(pp),
        vehicle_perturbed_min_pct = min(pp), vehicle_perturbed_max_pct = max(pp),
        fraction_changed_call = mean((pp <= 50) != (100 * 2^predicted <= 50)))
    }
  }
  save_table(rbindlist(kernel_summary), 'kernel_trace')
  save_table(rbindlist(neighbors), 'training_neighbors')
  save_table(rbindlist(noise_summary), 'vehicle_perturbation_summary')

  # A controlled feature-resolution comparison at the same dose and same split.
  status(paste(family, 'gene resolution and baseline-source diagnostic'))
  rows <- match(ix$observation_id, prepared$metadata$observation_id)
  scores <- prepared$measured[rows, , drop = FALSE]
  reconstructed <- as.matrix(gene_response[, bank$genes, drop = FALSE] %*% bank$weights)
  stopifnot(max(abs(reconstructed - scores)) < 1e-12)
  baseline_tahoe <- as.matrix(gene_baseline[, bank$genes, drop = FALSE] %*% bank$weights)
  fitting <- which(ix$depmap_id %in% original$training_ids)
  held <- which(ix$depmap_id %in% settings$held_ids)
  stopifnot(length(fitting) == 26L, length(held) == 9L)
  fold_map <- fread(file.path(parent, 'data/nonlinear_training_folds.tsv'))
  fold <- fold_map$fold[match(ix$depmap_id[fitting], fold_map$depmap_id)]
  y <- log2(ix$prism_abundance_pct / 100)
  matrices <- list('40 displacement means' = scores, 'Bank genes before averaging' = gene_response[, bank$genes, drop = FALSE],
    'All measured genes' = gene_response, 'CCLE baseline means' = prepared$baseline[ix$depmap_id, , drop = FALSE],
    'Tahoe baseline means' = baseline_tahoe)
  for (kind in names(matrices)) {
    x <- matrices[[kind]]
    active <- apply(x[fitting, , drop = FALSE], 2L, sd) > 1e-10
    x <- x[, active, drop = FALSE]
    cv <- cv.glmnet(x[fitting, , drop = FALSE], y[fitting], alpha = 0,
      lambda = 10^seq(3, -8, length.out = 45L), foldid = fold, standardize = TRUE,
      type.measure = 'mae', keep = TRUE, thresh = 1e-8)
    qsave(list(fit = cv, features = colnames(x), training_ids = ix$depmap_id[fitting]),
      file.path(root, 'models', paste0(family, '_', match(kind, names(matrices)), '.qs')))
    pp <- as.numeric(predict(cv, newx = x, s = 'lambda.min'))
    pred <- ix[, .(observation_id, depmap_id, recipient, observed_pct = prism_abundance_pct)]
    pred[, `:=`(family = family, representation = kind, cohort = ifelse(seq_len(.N) %in% fitting, 'Fitting cells', 'Held cells'),
      predicted_pct = 100 * 2^pp)]
    representation_predictions[[paste(family, kind)]] <- pred
    representation_selection[[paste(family, kind)]] <- data.table(family = family, representation = kind,
      features = ncol(x), selected_lambda = cv$lambda.min)
  }
  # Display raw differences for the named discordant same-treatment pair.
  ids <- c('ACH-000138', 'ACH-000927')
  gi <- match(ids, ix$depmap_id)
  delta <- gene_response[gi[1L], ] - gene_response[gi[2L], ]
  noise_sd <- sqrt(apply(profiles[[ids[1L]]]$null, 2L, var) + apply(profiles[[ids[2L]]]$null, 2L, var))
  gene_z <- delta / pmax(noise_sd, 1e-12)
  top <- order(abs(gene_z), decreasing = TRUE)[1:50]
  gene_differences[[family]] <- data.table(family = family, gene = genes[top],
    cfpac_displacement = gene_response[gi[1L], top], bt474_displacement = gene_response[gi[2L], top],
    difference = delta[top], difference_vehicle_sd = gene_z[top], in_bank = genes[top] %in% bank$genes,
    program = colnames(scores)[bank$partition[match(genes[top], bank$genes)]])
  null_a <- as.matrix(profiles[[ids[1L]]]$null[, bank$genes, drop = FALSE] %*% bank$weights)
  null_b <- as.matrix(profiles[[ids[2L]]]$null[, bank$genes, drop = FALSE] %*% bank$weights)
  pairs <- CJ(a = seq_len(nrow(null_a)), b = seq_len(nrow(null_b)))
  background <- sqrt(rowMeans((null_a[pairs$a, ] - null_b[pairs$b, ])^2))
  gap <- sqrt(mean((scores[gi[1L], ] - scores[gi[2L], ])^2))
  pair_summaries[[family]] <- data.table(family = family, gene_profile_correlation = cor(gene_response[gi[1L], ], gene_response[gi[2L], ]),
    program_gap = gap, median_vehicle_gap = median(background), fraction_vehicle_gaps_below_observed = mean(background <= gap),
    genes_above_two_vehicle_sd = sum(abs(gene_z) > 2), largest50_outside_bank = sum(!genes[top] %in% bank$genes))
  program_values[[family]] <- data.table(family = family, program = colnames(scores),
    cfpac_displacement = scores[gi[1L], ], bt474_displacement = scores[gi[2L], ],
    cfpac_ccle = prepared$baseline[ids[1L], ], bt474_ccle = prepared$baseline[ids[2L], ],
    cfpac_tahoe = baseline_tahoe[gi[1L], ], bt474_tahoe = baseline_tahoe[gi[2L], ])
}
predictions <- rbindlist(representation_predictions)
prediction_path <- file.path(root, 'data/outcome_free_representation_predictions.tsv')
fwrite(predictions[, .(observation_id, depmap_id, recipient, family, representation, cohort, predicted_pct)],
  prediction_path, sep = '\t')
save_table(data.table(path = prediction_path, sha256 = digest::digest(file = prediction_path, algo = 'sha256')), 'prediction_hash')
save_table(predictions, 'representation_predictions')
save_table(rbindlist(representation_selection), 'representation_selection')
summary <- predictions[, .(n = .N, pearson_r = cor_safe(observed_pct, predicted_pct),
  sensitive = sum(observed_pct <= 50), sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
  resistant = sum(observed_pct > 50), false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)),
  by = .(family, representation, cohort)]
save_table(summary, 'representation_summary')
save_table(rbindlist(pair_summaries), 'cfpac_bt474_resolution')
save_table(rbindlist(program_values), 'cfpac_bt474_programs')
save_table(rbindlist(gene_differences), 'cfpac_bt474_genes')
save_table(rbindlist(noise_predictions), 'vehicle_perturbation_predictions')
save_table(rbindlist(lineage_tables), 'lineage_membership')
save_table(rbindlist(lineage_tables)[, .(cells = .N, names = paste(CellLineName, collapse = ', ')),
  by = .(family, cohort, OncotreeLineage)], 'lineage_coverage')
source_paths <- unique(source_paths)
save_table(data.table(path = source_paths, sha256 = vapply(source_paths, digest::digest,
  character(1), file = TRUE, algo = 'sha256')), 'input_hashes')

status('plotting case evidence')
trace <- rbindlist(kernel_summary)
case_observed <- unique(trace[, .(family, case, observed_pct)])
case_observed[, `:=`(model = 'Observed', predicted_pct = observed_pct)]
bars <- rbindlist(list(trace[, .(family, case, model, predicted_pct)],
  case_observed[, .(family, case, model, predicted_pct)]))
bars[, model := factor(model, levels = c('Observed', 'Original', 'Local', 'Very_local'))]
p <- ggplot(bars, aes(model, predicted_pct, fill = model)) + geom_col(width = .7) +
  geom_hline(yintercept = 50, color = 'grey55', linetype = 2) + facet_grid(family ~ case, scales = 'free_y') +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 50, hjust = 1), legend.position = 'none') +
  labs(x = NULL, y = 'Abundance (% vehicle)', title = 'Named failure cases: measured response and frozen predictions',
    subtitle = 'Low: ixazomib 0.05 uM; middle: 0.5 uM; docetaxel: 0.5 uM')
save_plot(p, 'case_predictions', 12, 6)
nn <- rbindlist(neighbors)[scope == 'Same drug and dose' & case %in% c('CFPAC low', 'BT474 low', 'SW900 low')]
p <- ggplot(nn, aes(rank, prism_abundance_pct)) + geom_point(size = 2, color = '#2166AC') +
  geom_text(aes(label = recipient), vjust = -1, size = 2.3) +
  geom_hline(aes(yintercept = query_abundance_pct), color = '#B2182B', linewidth = .6) +
  geom_hline(yintercept = 50, color = 'grey60', linetype = 2) + facet_grid(family ~ case, scales = 'free_y') +
  scale_x_continuous(breaks = 1:5, expand = expansion(add = .65)) +
  scale_y_continuous(expand = expansion(mult = c(.06, .13))) + theme_classic(base_size = 9) +
  labs(x = 'Nearest training profile (same drug and dose)', y = 'Measured abundance (% vehicle)',
    title = 'Do molecular neighbors have the right phenotype?', subtitle = 'Blue: nearest training cells; red line: the held cell')
save_plot(p, 'same_treatment_neighbors', 10, 6)
pdf(file.path(root, 'figures/gene_resolution_models.pdf'), width = 12, height = 5.5)
for (family_name in c('MEK', 'ERBB')) {
  values <- predictions[family == family_name]
  a <- summary[family == family_name]
  a[, label := sprintf('r=%.2f; sensitive %d/%d', pearson_r, sensitive_called, sensitive)]
  p <- ggplot(values, aes(observed_pct, predicted_pct)) + geom_point(aes(color = observed_pct <= 50), size = 1.4) +
    geom_smooth(aes(group = 1), method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .4) +
    geom_hline(yintercept = 50, color = 'grey70', linetype = 2) + geom_vline(xintercept = 50, color = 'grey70', linetype = 2) +
    geom_text(data = a, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.03, vjust = 1.2, size = 2.1) +
    facet_wrap(vars(cohort, representation), ncol = 5, scales = 'free') + theme_classic(base_size = 8) +
    theme(aspect.ratio = 1, legend.position = 'none') +
    labs(x = 'Observed abundance (% vehicle)', y = 'Predicted abundance (% vehicle)',
      title = paste(family_name, '| ixazomib 0.05 uM: genes versus program means'),
      subtitle = 'Identical training/held cells; training-only regularization selection; held results retrospective')
  print(p)
  save_plot(p, paste0('gene_resolution_', family_name), 12, 5.5)
}
dev.off()
noise <- rbindlist(noise_predictions)
p <- ggplot(noise[model != 'Very_local'], aes(model, predicted_pct)) + geom_boxplot(outlier.size = .6) +
  geom_point(data = trace[model != 'Very_local'], color = '#B2182B', size = 1.8) +
  geom_hline(yintercept = 50, color = 'grey60', linetype = 2) + facet_grid(family ~ case, scales = 'free_y') +
  theme_classic(base_size = 9) + labs(x = NULL, y = 'Predicted abundance (% vehicle)',
    title = 'Sensitivity to empirical vehicle contrasts',
    subtitle = 'Boxes: predictions after adding centered vehicle contrasts; red: original prediction. These are not confidence intervals.')
save_plot(p, 'vehicle_perturbations', 12, 6)

# Raw gene changes and program means for the same discordant pair. Genes are
# displayed because they differ in this pair; this is not predictive selection.
pdf(file.path(root, 'figures/cfpac_bt474_molecular_profiles.pdf'), width = 5.5, height = 8)
for (family_name in c('MEK', 'ERBB')) {
  tab <- copy(gene_differences[[family_name]])
  tab[, group := ifelse(is.na(program), 'Outside bank', program)]
  setorder(tab, group, gene)
  mm <- as.matrix(tab[, .(cfpac_displacement, bt474_displacement)])
  dimnames(mm) <- list(tab$gene, c('CFPAC-1 (25%)', 'BT-474 (156%)'))
  annotation <- data.frame(Program = tab$group, row.names = tab$gene)
  stopifnot(identical(rownames(mm), rownames(annotation)))
  limit <- max(abs(mm))
  hm <- draw_heatmap(mm, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    row_anno_df = annotation, row_anno_cols = 'Program', show_row_names = TRUE, show_column_names = TRUE,
    row_split = tab$group, row_title = NULL, row_names_gp = grid::gpar(fontsize = 8),
    value_range = c(-limit, 0, limit), value_colors = c('#2166AC', '#FFFFFF', '#B2182B'),
    hm_name = 'Displacement', title = paste(family_name, '| largest 50 vehicle-scaled gene differences'))
  ComplexHeatmap::draw(hm)
  tab <- program_values[[family_name]]
  mm <- as.matrix(tab[, .(cfpac_displacement, bt474_displacement)])
  dimnames(mm) <- list(tab$program, c('CFPAC-1 (25%)', 'BT-474 (156%)'))
  order <- hclust(dist(mm), method = 'ward.D2')$order
  limit <- max(abs(mm))
  hm <- draw_heatmap(mm[order, , drop = FALSE], scale_rows = FALSE, cluster_rows = FALSE,
    cluster_columns = FALSE, show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8), value_range = c(-limit, 0, limit),
    value_colors = c('#2166AC', '#FFFFFF', '#B2182B'), hm_name = 'Displacement',
    title = paste(family_name, '| the same signatures after program averaging'))
  ComplexHeatmap::draw(hm)
}
dev.off()
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
status('COMPLETE')
print(summary)
print(rbindlist(pair_summaries))
