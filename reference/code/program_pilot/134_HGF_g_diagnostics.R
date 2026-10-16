#!/usr/bin/env Rscript
# Read-only model interrogation: no coefficient fitting and no outcome-based selection.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs)})
root <- '/work/results/program_pilot/134_HGF_g_diagnostics/linear'
input <- '/work/results/program_pilot/133_counts_and_growth_order_bridge/linear'
dir.create(file.path(root, 'data'), recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0(name, '.tsv')), sep = '\t')
model <- qread(file.path(input, 'models/final_normalized_cell_count.qs'))
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
receptors <- as.data.table(archive$prep$receptors)
P <- sprintf('P%02d', 1:40)
ctx <- model$context
lap <- model$rna$lap_receptors
X <- as.matrix(receptors[depmap_id == 'ACH-000248', .(ERBB2, MET)])
stopifnot(nrow(X) == 1L)
Z <- sweep(sweep(X, 2, lap$xmean), 2, lap$xsd, '/') %*% lap$B
L <- as.numeric(lap$mu + as.numeric(Z) * lap$sigma)
g <- 6L
stopifnot(rownames(ctx$mu)[g] == 'Lapatinib ditosylate')
profiles <- fread(file.path(input, 'data/time_matched_HGF_program_inputs.tsv'))[CellLine == 'AU565']
H <- as.matrix(profiles[, ..P])
Qh <- cbind(1, sweep(sweep(sweep(H, 2, ctx$fm), 2, ctx$fs, '/') %*% ctx$fv, 2, ctx$fscale, '/'))
scale_map <- solve(crossprod(ctx$Q[1:10, ]) / 10 + diag(c(0, 1, 1, 1)),
                   crossprod(ctx$Q[1:10, ], log(ctx$sigma[1:10, ])) / 10)
minimum <- matrix(apply(ctx$sigma[1:10, ], 2, min), nrow(H), 40, byrow = TRUE)
maximum <- matrix(apply(ctx$sigma[1:10, ], 2, max), nrow(H), 40, byrow = TRUE)
hscale <- pmin(pmax(exp(Qh %*% scale_map), minimum), maximum)
Lscaled <- L / ctx$sigma[g, ]
Hscaled <- H / hscale
Hcommon <- sweep(H, 2, ctx$sigma[g, ], '/')
pair_diagnostics <- rbindlist(lapply(seq_len(nrow(H)), function(i) {
  data.table(source_profile = profiles$source_profile[i],
    r_across_40_raw_programs = cor(L, H[i, ]),
    cosine_raw = sum(L * H[i, ]) / sqrt(sum(L^2) * sum(H[i, ]^2)),
    r_across_40_model_scaled_programs = cor(Lscaled, Hscaled[i, ]),
    opposite_sign_programs = sum(L * H[i, ] < 0),
    total_programs = 40L,
    lap_route1_absolute = sum(Lscaled * ctx$W[, 1]),
    HGF_route1_own_scale = sum(Hscaled[i, ] * ctx$W[, 1]),
    HGF_route1_lapatinib_scale = sum(Hcommon[i, ] * ctx$W[, 1]),
    lap_route2_absolute = sum(Lscaled * ctx$W[, 2]),
    HGF_route2_own_scale = sum(Hscaled[i, ] * ctx$W[, 2]),
    HGF_route2_lapatinib_scale = sum(Hcommon[i, ] * ctx$W[, 2]))
}))
write_table(pair_diagnostics, 'AU565_lapatinib_HGF_program_opposition')
programs <- rbindlist(lapply(seq_len(nrow(H)), function(i) {
  data.table(source_profile = profiles$source_profile[i], program = P,
    lapatinib_displacement = L, HGF_displacement = H[i, ],
    opposite_sign = L * H[i, ] < 0,
    lapatinib_route1_contribution = Lscaled * ctx$W[, 1],
    HGF_route1_own_scale_contribution = Hscaled[i, ] * ctx$W[, 1],
    HGF_route1_lapatinib_scale_contribution = Hcommon[i, ] * ctx$W[, 1])
}))
write_table(programs, 'AU565_program_level_opposition')
sources <- c('GSE135552', 'GSE26736', 'GSE71404', 'GSE98979')
source_rows <- match(sources, profiles$source_profile)
source_cor <- cor(t(H[source_rows, , drop = FALSE]))
dimnames(source_cor) <- list(sources, sources)
write_table(as.data.table(as.table(source_cor)), 'HGF_source_agreement_across_40_AU565_programs')
epithelial <- H[match(c('GSE135552', 'GSE71404', 'GSE98979'), profiles$source_profile), , drop = FALSE]
write_table(data.table(program = P, all_three_same_sign = apply(sign(epithelial), 2, function(z) length(unique(z)) == 1L)),
            'epithelial_study_program_sign_agreement')

# Summarize every overnight fit using its already-frozen result tables. All rows
# are retained; no best-source or best-model selection is performed here.
folders <- list.dirs('/work/results/program_pilot', recursive = FALSE, full.names = TRUE)
numbers <- suppressWarnings(as.integer(sub('_.*', '', basename(folders))))
folders <- folders[is.finite(numbers) & numbers >= 108 & numbers <= 133]
all_queries <- metrics <- list()
for (folder in folders) {
  version <- basename(folder)
  path <- file.path(folder, 'linear/data/au565_hgf_5uM_predictions_frozen.tsv')
  if (file.exists(path)) {
    d <- fread(path)
    d[, version := version]
    all_queries[[version]] <- d
  }
  for (name in c('mema_correlations', 'combination_correlations', 'HER2_external_held_cell_correlations')) {
    path <- file.path(folder, 'linear/data', paste0(name, '.tsv'))
    if (file.exists(path)) {
      d <- fread(path)
      d[, `:=`(version = version, result_table = name)]
      metrics[[paste(version, name)]] <- d
    }
  }
}
write_table(rbindlist(all_queries, fill = TRUE), 'all_iterations_frozen_HGF_queries')
write_table(rbindlist(metrics, fill = TRUE), 'all_iterations_validation_correlations')
write_table(data.table(coefficients = length(model$par$Gamma),
  molecular_prototypes = nrow(model$interaction_basis$centres),
  ridge_penalty = model$interaction_training$interaction_penalty,
  MEMA_scale = model$mema_assay_scale,
  true_count_external_observations = model$interaction_training$exact_count_rows), 'current_g_dimensions')
writeLines(c('No models refitted. No molecular source selected using outcomes.',
  'Raw programme correlations compare 40 coordinates within AU565, not predictive correlation across cell lines.',
  'HGF is transplanted from donor gene responses; AU565 HGF RNA is not measured.',
  'Program opposition cannot establish restoration of a specific phosphoprotein pathway.',
  'Current g learns centered MEMA and ALMANAC contrasts; only external count observations also calibrate means.',
  'Existing source-exclusion fits retain the original additional HGF coordinate loading; they are not a full lineage-only rediscovery.'),
  file.path(root, 'DIAGNOSTIC_SCOPE.txt'))
print(pair_diagnostics)
print(data.table(source = sources, source_cor))
cat('Same-sign programmes across three epithelial studies:', sum(apply(sign(epithelial), 2, function(z) length(unique(z)) == 1L)), '/40\n')
cat('Saved diagnostic tables:', root, '\n')
