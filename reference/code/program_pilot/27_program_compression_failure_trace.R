#!/usr/bin/env Rscript
# Locate response distinctions lost between genes, local programs and frozen banks.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(Matrix)
  library(ggplot2)
})
source('/work/src/R/lib_heatmap.R')
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
root <- '/work/results/program_pilot/27_program_compression_failure_trace'
for (directory in file.path(root, c('data', 'figures'))) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
panel_path <- '/work/results/program_pilot/11_panel_functional_support_programs_nmf/data/molecular_inputs/panel_index.qs'
local_root <- '/work/results/program_pilot/15_mean_functional_displacement_correlation_programs_ward'
global_path <- '/work/results/program_pilot/16_global_mean_relationships_fixed_ward_programs/data/global_gene_assignments.tsv'
audit_path <- '/work/results/program_pilot/20_multidrug_baseline_response_function_MEK_holdout/data/contextual_audit_cases.tsv'
panel <- qread(panel_path)
labels <- as.data.table(panel$sensitivity)
cases <- fread(audit_path)[case %in% c('C01', 'C19', 'C25')]
cases <- cases[match(c('C01', 'C19', 'C25'), case)]
stopifnot(nrow(cases) == 3L, !anyNA(cases$case))
manifest_paths <- c(panel_path, global_path, audit_path,
  '/work/src/program_pilot/27_program_compression_failure_trace.R', '/work/src/R/lib_heatmap.R')

# Hard projections reproduce the direct cluster means used by the predictor.
hard_weights <- function(genes, partition) {
  levels <- sort(unique(partition))
  program_names <- if (is.numeric(levels)) sprintf('P%02d', levels) else as.character(levels)
  partition <- as.integer(factor(partition, levels = levels))
  size <- tabulate(partition)
  sparseMatrix(i = seq_along(genes), j = partition, x = 1 / size[partition],
    dims = c(length(genes), length(size)),
    dimnames = list(genes, program_names))
}
global <- fread(global_path)
global_weights <- hard_weights(global$gene, global$global_program)
frozen <- lapply(c('MEK', 'ERBB'), function(family) {
  path <- paste0('/work/results/program_pilot/20_multidrug_baseline_response_function_',
    family, '_holdout/models/program_bank.qs')
  manifest_paths <<- c(manifest_paths, path)
  qread(path)
})
names(frozen) <- c('MEK', 'ERBB')
capacity <- lapply(c('MEK', 'ERBB'), function(family) {
  path <- paste0('/work/results/program_pilot/24_ixazomib_decoder_capacity_and_molecular_signal_',
    family, '_holdout/data/shared_capacity_predictions.tsv.gz')
  manifest_paths <<- c(manifest_paths, path)
  x <- fread(path)
  x[, family := family]
  x
})
capacity <- rbindlist(capacity)
profiles <- list()
summary_tables <- coverage_tables <- program_tables <- gene_tables <- prediction_tables <- list()
unavailable <- list()

# Vehicle contrasts are an empirical between-plate reference, not treated replicates.
gap_summary <- function(response, null_i, null_j, same_cell) {
  gap <- sqrt(mean((response[1L, ] - response[2L, ])^2))
  if (same_cell) {
    background <- as.vector(dist(null_i)) / sqrt(ncol(null_i))
  } else {
    background <- as.vector(sapply(seq_len(nrow(null_j)), function(j) {
      sqrt(rowMeans(sweep(null_i, 2L, null_j[j, ], '-')^2))
    }))
  }
  data.table(coordinates = ncol(response), gap = gap,
    median_vehicle_gap = median(background), gap_to_vehicle_median = gap / median(background),
    vehicle_percentile = mean(background < gap))
}
save_plot <- function(plot, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), plot, width = width, height = height)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), plot, width = width, height = height, dpi = 140)
}

for (case_index in seq_len(nrow(cases))) {
  cc <- cases[case_index]
  message(format(Sys.time()), ' | ', cc$case, ': loading exact case profiles')
  lab <- rbindlist(list(labels[recipient == cc$cell_i & drug == cc$drug_i & dose == cc$dose],
    labels[recipient == cc$cell_j & drug == cc$drug_j & dose == cc$dose]))
  stopifnot(nrow(lab) == 2L, all(is.finite(lab$prism_abundance_pct)))
  for (id in unique(lab$depmap_id)) {
    if (is.null(profiles[[id]])) {
      path <- panel$input_paths[match(id, panel$cell_map$depmap_id)]
      stopifnot(length(path) == 1L, !is.na(path))
      manifest_paths <- c(manifest_paths, path)
      profiles[[id]] <- qread(path)
    }
  }
  mi <- profiles[[lab$depmap_id[1L]]]
  mj <- profiles[[lab$depmap_id[2L]]]
  stopifnot(identical(mi$genes, mj$genes))
  response <- rbind(mi$response[lab$observation_id[1L], , drop = FALSE],
    mj$response[lab$observation_id[2L], , drop = FALSE])
  null_i <- mi$null[[lab$depmap_id[1L]]]
  null_j <- mj$null[[lab$depmap_id[2L]]]
  same_cell <- lab$depmap_id[1L] == lab$depmap_id[2L]
  delta <- response[2L, ] - response[1L, ]
  vehicle_sd <- if (same_cell) matrixStats::colSds(null_i) else
    sqrt((matrixStats::colVars(null_i) + matrixStats::colVars(null_j)) / 2)
  stopifnot(all(vehicle_sd >= 0), all(is.finite(vehicle_sd)))
  # Constant vehicle genes have no estimable SD ratio; retain raw contrasts but
  # exclude them from the vehicle-scaled display ranking rather than invent a floor.
  vehicle_scaled <- rep(NA_real_, length(delta))
  variable_vehicle <- vehicle_sd > 0
  vehicle_scaled[variable_vehicle] <- abs(delta[variable_vehicle]) / vehicle_sd[variable_vehicle]
  message(sum(!variable_vehicle), ' genes have constant vehicle contrasts; scaled contrast is NA')
  gene_tab <- data.table(case = cc$case, gene = mi$genes, displacement_i = response[1L, ],
    displacement_j = response[2L, ], difference = delta, vehicle_sd = vehicle_sd,
    difference_vehicle_sd = vehicle_scaled)

  representations <- list()
  local_banks <- list()
  for (recipient in unique(lab$recipient)) {
    scope <- paste0('cell_', gsub('[^A-Za-z0-9]+', '_', recipient))
    path <- file.path(local_root, scope, 'models/program_bank.qs')
    selection_path <- file.path(local_root, scope, 'data/gene_selection.tsv.gz')
    manifest_paths <- c(manifest_paths, path, selection_path)
    local <- qread(path)
    selection <- fread(selection_path)
    local_banks[[recipient]] <- local
    gene_tab[, (paste0('filter_', scope)) := selection$exclusion_reason[match(gene, selection$gene)]]
    if (is.null(local$weights)) {
      unavailable[[paste(cc$case, recipient)]] <- data.table(case = cc$case, recipient = recipient,
        reason = local$diagnostics$selection_status, admitted_genes = length(local$genes))
      message(recipient, ': no saved local scores; ', local$diagnostics$selection_status)
      next
    }
    own <- which(lab$recipient == recipient)
    replay <- response[own, local$genes, drop = FALSE] %*% local$weights
    stopifnot(max(abs(replay - local$scores[lab$observation_id[own], , drop = FALSE])) < 1e-10)
    representations[[paste0('15 soft | ', recipient)]] <- local$weights
    representations[[paste0('15 hard | ', recipient)]] <- hard_weights(local$genes, local$partition)
    gene_tab[, (paste0('local_', scope)) := as.character(local$partition[match(gene, local$genes)])]
  }
  representations[['16 global hard | descriptive']] <- global_weights
  representations[['20 global hard | MEK']] <- frozen$MEK$weights
  representations[['20 global hard | ERBB']] <- frozen$ERBB$weights

  metrics <- list(All_genes = cbind(data.table(case = cc$case, representation = 'All measured genes'),
    gap_summary(response, null_i, null_j, same_cell)))
  scores_list <- list()
  for (representation in names(representations)) {
    weights <- representations[[representation]]
    genes <- rownames(weights)
    stopifnot(all(genes %in% mi$genes), max(abs(colSums(weights) - 1)) < 1e-10)
    scores <- as.matrix(response[, genes, drop = FALSE] %*% weights)
    control_i <- as.matrix(null_i[, genes, drop = FALSE] %*% weights)
    control_j <- as.matrix(null_j[, genes, drop = FALSE] %*% weights)
    metrics[[representation]] <- cbind(data.table(case = cc$case, representation = representation),
      gap_summary(scores, control_i, control_j, same_cell))
    metrics[[paste0(representation, '_genes')]] <- cbind(data.table(case = cc$case,
      representation = paste0(representation, ' | genes before averaging')),
      gap_summary(response[, genes, drop = FALSE], null_i[, genes, drop = FALSE],
        null_j[, genes, drop = FALSE], same_cell))
    signed_contributions <- sweep(as.matrix(weights), 1L, delta[genes], '*')
    mean_difference <- colSums(signed_contributions)
    absolute_difference <- colSums(abs(signed_contributions))
    program_tables[[paste(cc$case, representation)]] <- data.table(case = cc$case,
      representation = representation, program = colnames(weights), score_i = scores[1L, ],
      score_j = scores[2L, ], difference = mean_difference,
      weighted_absolute_gene_difference = absolute_difference,
      signed_fraction_retained = abs(mean_difference) / absolute_difference,
      vehicle_sd = if (same_cell) matrixStats::colSds(control_i) else
        sqrt((matrixStats::colVars(control_i) + matrixStats::colVars(control_j)) / 2))
    ranked <- gene_tab[is.finite(difference_vehicle_sd)][order(-difference_vehicle_sd), gene]
    coverage_tables[[paste(cc$case, representation)]] <- data.table(case = cc$case,
      representation = representation, selected_genes = length(genes),
      top50_genes_retained = sum(head(ranked, 50L) %in% genes),
      fraction_squared_gene_contrast_retained = sum(delta[genes]^2) / sum(delta^2),
      median_signed_fraction_retained = median(abs(mean_difference) / absolute_difference))
    scores_list[[representation]] <- scores
  }
  summary_tables[[cc$case]] <- rbindlist(metrics)
  gene_tables[[cc$case]] <- gene_tab
  predictions <- capacity[observation_id %in% lab$observation_id,
    .(family, model, observation_id, recipient, drug, dose, observed_pct, predicted_pct)]
  predictions[, case := cc$case]
  prediction_tables[[cc$case]] <- predictions

  # Commit the first case before reading or analyzing the next one.
  fwrite(rbindlist(summary_tables), file.path(root, 'data/representation_gaps.tsv'), sep = '\t')
  fwrite(rbindlist(coverage_tables), file.path(root, 'data/gene_coverage_and_cancellation.tsv'), sep = '\t')
  fwrite(rbindlist(program_tables), file.path(root, 'data/program_contrasts.tsv'), sep = '\t')
  fwrite(rbindlist(gene_tables, fill = TRUE), file.path(root, 'data/gene_contrasts.tsv.gz'), sep = '\t')
  fwrite(rbindlist(prediction_tables), file.path(root, 'data/decoder_predictions.tsv'), sep = '\t')
  if (length(unavailable)) fwrite(rbindlist(unavailable), file.path(root, 'data/unavailable_local_scores.tsv'), sep = '\t')
  fwrite(lab, file.path(root, 'data', paste0(cc$case, '_labels.tsv')), sep = '\t')
  fwrite(rbindlist(list(mi$metadata[match(lab$observation_id[1L], observation_id)],
    mj$metadata[match(lab$observation_id[2L], observation_id)])),
    file.path(root, 'data', paste0(cc$case, '_molecular_metadata.tsv')), sep = '\t')

  display <- rbindlist(program_tables)[case == cc$case &
    (grepl('15 soft', representation) | grepl('20 global', representation))]
  display <- melt(display, id.vars = c('representation', 'program'),
    measure.vars = c('score_i', 'score_j'), variable.name = 'signature', value.name = 'displacement')
  display[, signature := ifelse(signature == 'score_i',
    sprintf('%s / %s (%.1f%%)', cc$cell_i, cc$drug_i, lab$prism_abundance_pct[1L]),
    sprintf('%s / %s (%.1f%%)', cc$cell_j, cc$drug_j, lab$prism_abundance_pct[2L]))]
  plot <- ggplot(display, aes(program, displacement, color = signature, group = signature)) +
    geom_line(linewidth = .5) + geom_point(size = 1.3) +
    facet_wrap(vars(representation), ncol = 1, scales = 'free_x') + theme_classic(base_size = 9) +
    theme(legend.position = 'top', legend.text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    labs(x = 'Program', y = 'Measured program displacement', color = NULL,
      title = paste(cc$case, '| original local scores and frozen global means'),
      subtitle = 'Every panel applies one fixed gene projection to both signatures')
  save_plot(plot, paste0(cc$case, '_program_profiles'), 8, 2 + length(unique(display$representation)) * 1.3)

  if (same_cell) {
    local <- local_banks[[1L]]
    matched <- labels[depmap_id == lab$depmap_id[1L] & is.finite(prism_abundance_pct) & prism_abundance_pct > 0]
    matched <- matched[observation_id %in% rownames(local$scores)]
    points <- cbind(matched[, .(observation_id, drug, dose, prism_abundance_pct)],
      as.data.table(local$scores[matched$observation_id, , drop = FALSE]))
    points <- melt(points, id.vars = c('observation_id', 'drug', 'dose', 'prism_abundance_pct'),
      variable.name = 'program', value.name = 'displacement')
    points[, highlighted := ifelse(observation_id %in% lab$observation_id, drug, 'Other drugs/doses')]
    associations <- points[, .(r = cor(displacement, prism_abundance_pct)), by = program]
    associations[, label := sprintf('Pearson r = %.2f', r)]
    associations[, case := cc$case]
    fwrite(associations, file.path(root, 'data', paste0(cc$case, '_local_program_associations.tsv')), sep = '\t')
    dose_associations <- points[, .(conditions = .N, r = cor(displacement, prism_abundance_pct)), by = .(program, dose)]
    fwrite(dose_associations, file.path(root, 'data', paste0(cc$case, '_local_program_associations_by_dose.tsv')), sep = '\t')
    score_correlation <- cor(local$scores)
    eigenvalues <- eigen(score_correlation, symmetric = TRUE, only.values = TRUE)$values
    fwrite(data.table(case = cc$case, programs = ncol(local$scores),
      first_component_variance_fraction = eigenvalues[1L] / sum(eigenvalues),
      effective_dimensions = sum(eigenvalues)^2 / sum(eigenvalues^2)),
      file.path(root, 'data', paste0(cc$case, '_local_program_redundancy.tsv')), sep = '\t')
    fwrite(cbind(data.table(program = colnames(score_correlation)), as.data.table(score_correlation)),
      file.path(root, 'data', paste0(cc$case, '_local_program_correlation.tsv')), sep = '\t')
    plot <- ggplot(points, aes(displacement, prism_abundance_pct)) +
      geom_point(aes(color = highlighted), size = 1.3, alpha = .6) +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'grey40', linewidth = .5) +
      geom_point(data = points[highlighted != 'Other drugs/doses'], aes(color = highlighted), size = 3) +
      geom_text(data = associations, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
        hjust = -.1, vjust = 1.2, size = 3) +
      facet_wrap(vars(program), ncol = 3, scales = 'free') + theme_classic(base_size = 9) +
      theme(aspect.ratio = 1, legend.position = 'bottom') +
      labs(x = 'Script-15 program displacement', y = 'PRISM abundance (% vehicle)', color = NULL,
        title = paste(cc$case, '| local association and the discordant pair'),
        subtitle = 'All matched drugs/doses in CFPAC-1; original outcome-informed local programs')
    save_plot(plot, paste0(cc$case, '_local_function_associations'), 9, 6.5)
  }

  # Show raw signed gene values; row selection is descriptive, never predictive.
  top <- gene_tab[is.finite(difference_vehicle_sd)][order(-difference_vehicle_sd)][1:40]
  reference <- local_banks[[1L]]
  group <- ifelse(top$gene %in% reference$genes,
    paste0('Local P', sprintf('%02d', reference$partition[match(top$gene, reference$genes)])), 'Outside local bank')
  order <- order(group, -top$difference_vehicle_sd)
  matrix <- t(response[, top$gene[order], drop = FALSE])
  colnames(matrix) <- c(sprintf('%s / %s (%.1f%%)', cc$cell_i, cc$drug_i, lab$prism_abundance_pct[1L]),
    sprintf('%s / %s (%.1f%%)', cc$cell_j, cc$drug_j, lab$prism_abundance_pct[2L]))
  annotation <- data.frame(Program = group[order], row.names = top$gene[order])
  stopifnot(identical(rownames(matrix), rownames(annotation)))
  limit <- max(abs(matrix))
  heatmap <- draw_heatmap(matrix, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    row_anno_df = annotation, row_anno_cols = 'Program', row_split = annotation$Program, row_title = NULL,
    show_row_names = TRUE, show_column_names = TRUE, row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8), value_range = c(-limit, 0, limit),
    value_colors = c('#2166AC', 'white', '#B2182B'), hm_name = 'Displacement',
    title = paste(cc$case, '| top 40 vehicle-scaled gene contrasts'))
  pdf(file.path(root, 'figures', paste0(cc$case, '_gene_contrasts.pdf')), width = 6.4, height = 7.5)
  ComplexHeatmap::draw(heatmap)
  dev.off()
  png(file.path(root, 'figures', paste0(cc$case, '_gene_contrasts.png')), width = 896, height = 1050, res = 140)
  ComplexHeatmap::draw(heatmap)
  dev.off()
  fwrite(data.table(stage = paste(cc$case, 'complete'), time = as.character(Sys.time())),
    file.path(root, 'run_status.tsv'), sep = '\t')
  print(summary_tables[[cc$case]])
  print(rbindlist(coverage_tables)[case == cc$case])
  message(cc$case, ' results and figures saved')
}
manifest_paths <- unique(manifest_paths)
fwrite(data.table(path = manifest_paths, sha256 = vapply(manifest_paths, digest::digest,
  character(1), file = TRUE, algo = 'sha256')), file.path(root, 'data/input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
fwrite(data.table(stage = 'COMPLETE', time = as.character(Sys.time())), file.path(root, 'run_status.tsv'), sep = '\t')
