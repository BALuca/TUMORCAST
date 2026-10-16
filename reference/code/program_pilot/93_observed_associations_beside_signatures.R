#!/usr/bin/env Rscript
# Compare observed correlations with the frozen decoder's two program signatures.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(qs); library(data.table)})
source('/work/src/R/lib_heatmap.R')

root <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes'
out <- file.path(root, 'figures/program_and_route_diagnostics/clustered')
clusters <- qread(file.path(root, 'data/correlation_cluster_inputs.qs'))
model <- qread(file.path(root, 'models/final.qs'))
selection <- fread(file.path(root, 'data/selected_doses.tsv'))
drug_names <- c('Azacytidine', 'Bendamustine', 'Bortezomib', 'Busulfan', 'Cytarabine',
               'Lapatinib', 'Topotecan', 'Vandetanib', 'Vismodegib', 'Crizotinib')
programs <- sprintf('P%02d', 1:40)
ro <- clusters$row_order
co <- clusters$column_order
observed <- clusters$matrices$observed_r[ro, co, drop=FALSE]

# Orient each route toward sensitivity using its PRISM readout sign. Flipping W
# and its output coefficient together leaves tanh predictions exactly unchanged.
orientation <- sign(-model$par$C[1, ])
stopifnot(all(orientation != 0), identical(dim(model$par$W), c(40L, 2L)))
signatures <- sweep(model$par$W, 2, orientation, '*')
dimnames(signatures) <- list(programs, c('Signature 1', 'Signature 2'))
signatures <- signatures[ro, , drop=FALSE]
stopifnot(identical(rownames(observed), rownames(signatures)))
stopifnot(max(abs(tanh(model$par$W) %*% model$par$C[1, ] -
                  tanh(sweep(model$par$W, 2, orientation, '*')) %*%
                  (model$par$C[1, ] * orientation))) < 1e-12)
row_groups <- paste0('Pattern ', clusters$row_groups[ro])
row_split <- factor(row_groups, levels=unique(row_groups))
annotation <- data.frame(Dose=paste0(selection$dose[co], ' µM'), row.names=colnames(observed))

# Same program rows across both panels; separate legends retain the actual units.
render <- function() {
  left <- draw_heatmap(observed, scale_rows=FALSE, value_range=c(-1, 0, 1),
    value_colors=c('#2166AC', 'white', '#B2182B'), hm_name='Observed Pearson r',
    col_anno_df=annotation, col_anno_cols='Dose',
    palette_overrides=list(Dose=c('0.05 µM'='#DCEAF2', '0.5 µM'='#6BAED6', '5 µM'='#08519C')),
    show_row_names=TRUE, row_names_side='left', show_column_names=TRUE,
    column_labels=drug_names[co], column_names_rot=45,
    row_names_gp=grid::gpar(fontsize=8), column_names_gp=grid::gpar(fontsize=8),
    cluster_rows=FALSE, cluster_columns=FALSE, row_split=row_split,
    cluster_row_slices=FALSE, row_gap=grid::unit(2, 'mm'),
    row_title_gp=grid::gpar(fontsize=8), width=grid::unit(105, 'mm'),
    column_title_gp=grid::gpar(fontsize=11), title='Measured program–sensitivity associations')
  limit <- max(abs(signatures))
  right <- draw_heatmap(signatures, scale_rows=FALSE, value_range=c(-limit, 0, limit),
    value_colors=c('#2166AC', 'white', '#B2182B'), hm_name='Signature weight',
    show_row_names=FALSE, show_column_names=TRUE, column_names_rot=45,
    column_names_gp=grid::gpar(fontsize=8), cluster_rows=FALSE, cluster_columns=FALSE,
    row_split=row_split, cluster_row_slices=FALSE, row_gap=grid::unit(2, 'mm'),
    row_title=NULL, width=grid::unit(24, 'mm'),
    column_title_gp=grid::gpar(fontsize=11), title='Learned signatures')
  drawn <- ComplexHeatmap::draw(left + right, heatmap_legend_side='bottom',
    annotation_legend_side='right', ht_gap=grid::unit(8, 'mm'),
    column_title='v93: observed associations beside two learned sensitivity signatures',
    column_title_gp=grid::gpar(fontsize=13),
    padding=grid::unit(c(16, 5, 5, 5), 'mm'))
  stopifnot(identical(as.integer(unlist(ComplexHeatmap::row_order(drawn, name='Observed Pearson r'))), 1:40))
  grid::grid.text('Same program order. Red = positive association / sensitivity-oriented route weight.\nLeft: correlations across cells within drug/dose. Right: fitted tanh input weights, not correlations or conditional betas.',
    x=.5, y=grid::unit(4, 'mm'), just=c('centre', 'bottom'), gp=grid::gpar(fontsize=8))
}
stem <- file.path(out, 'observed_correlations_beside_learned_signatures')
cairo_pdf(paste0(stem, '.pdf'), width=9, height=9)
render()
dev.off()
png(paste0(stem, '.png'), width=9, height=9, units='in', res=160, type='cairo')
render()
dev.off()
fwrite(data.table(program=programs[ro], signature_1=signatures[, 1], signature_2=signatures[, 2]),
       file.path(root, 'data/sensitivity_oriented_signature_loadings.tsv'), sep='\t')
writeLines(c('Frozen v93 fit; no refitting or selection.',
  'Left panel: Pearson correlation of measured displacement with measured -log2 PRISM abundance, within drug/dose.',
  'Right panel: W route loadings, each multiplied by sign(-C_PRISM) to orient toward sensitivity.',
  'Both signs were reversed for display; paired reversal of W and C leaves model predictions unchanged.',
  'Route weights feed tanh; calibration, drug mixture, saturation and assay readout determine their final effect.',
  'Rows and drug columns retain the observed-correlation Ward.D2 clustering order. No row z-scoring.'),
  paste0(stem, '_definitions.txt'))
cat(paste0(stem, '.pdf\n'))
