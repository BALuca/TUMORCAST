#!/usr/bin/env Rscript
# Query existing frozen models; no refitting or choice of HGF route from outcomes.
# Report intervention contrasts against zero displacement under the SAME transform.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs)})
folders <- c('104_erbb2_met_lapatinib_shared_decoder',
             '105_erbb2_met_without_crizotinib_vismodegib')
hgf <- fread('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv')
archive <- readRDS('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds')
receptors <- as.data.table(archive$prep$receptors)
programs <- sprintf('P%02d', 1:40)
hmean <- as.numeric(hgf[context == 'Equal-context mean', ..programs])
all_results <- list()
for (i in seq_along(folders)) {
  root <- file.path('/work/results/program_pilot', folders[i], 'linear')
  model <- qread(file.path(root, 'models/final.qs'))
  ctx <- model$context
  selected <- fread('/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data/selected_doses.tsv')
  g <- match('Lapatinib ditosylate', selected$drug)
  lapmap <- model$rna$lap_receptors
  receptor_values <- as.matrix(receptors[depmap_id == 'ACH-000248', .(ERBB2, MET)])
  lap_z <- sweep(sweep(receptor_values, 2, lapmap$xmean), 2, lapmap$xsd, '/') %*% lapmap$B
  lap_D <- as.numeric(sweep(sweep(lap_z, 2, lapmap$sigma, '*'), 2, lapmap$mu, '+'))
  lap_route <- ctx$route[g]
  lap_amplitude <- exp(as.numeric(ctx$Q[g, , drop = FALSE] %*% model$par$log_amplitude))
  lap_projection <- as.numeric((lap_D / ctx$sigma[g, ]) %*% ctx$W)
  lap_effect <- lap_amplitude * lap_projection[lap_route]
  LA <- matrix(0, 1, 2)
  LA[1, lap_route] <- lap_effect
  Qh <- cbind(1, sweep(matrix((hmean - ctx$fm) / ctx$fs, nrow = 1) %*% ctx$fv,
                        2, ctx$fscale, '/'))
  scale_map <- solve(crossprod(ctx$Q) / 10 + diag(c(0, 1, 1, 1)),
                     crossprod(ctx$Q, log(ctx$sigma)) / 10)
  hscale <- pmin(pmax(as.numeric(exp(Qh %*% scale_map)), apply(ctx$sigma, 2, min)),
                apply(ctx$sigma, 2, max))
  h_amplitude <- exp(as.numeric(Qh %*% model$par$log_amplitude))
  for (route in 1:2) {
    D <- as.matrix(hgf[, ..programs])
    h_effect <- h_amplitude * as.numeric(sweep(D, 2, hscale, '/') %*% ctx$W[, route])
    # Verify centered-head contrast equals the zero-referenced projection.
    centered <- as.numeric(sweep(sweep(D, 2, hmean), 2, hscale, '/') %*% ctx$W[, route])
    zero_score <- sum((-hmean / hscale) * ctx$W[, route])
    stopifnot(max(abs(h_effect - h_amplitude * (centered - zero_score))) < 1e-12)
    HA <- matrix(0, nrow(D), 2)
    HA[, route] <- h_effect
    A <- LA[rep(1L, nrow(D)), , drop = FALSE]
    products <- cbind(A[, 1] * HA[, 1], A[, 1] * HA[, 2] + A[, 2] * HA[, 1], A[, 2] * HA[, 2])
    J <- as.numeric(sweep(products, 2, ctx$interaction_scale, '/') %*% model$par$Gamma)
    out <- data.table(version = 103L + i, cell = 'AU565', lapatinib_uM = 5,
      hgf_context = hgf$context, assumed_HGF_route = route,
      lapatinib_pct_zero_control = 100 * 2^(-lap_effect),
      hgf_alone_pct_zero_control = 100 * 2^(-h_effect),
      combination_pct_zero_control = 100 * 2^(-lap_effect - h_effect - J),
      hgf_addition_relative_pct = 100 * (2^(-h_effect - J) - 1),
      interaction_relative_pct = 100 * (2^(-J) - 1),
      interaction_log2 = -J,
      reference = 'Each intervention contrasted against zero displacement using its identical frozen transformation; not raw PRISM-head predictions',
      qualification = 'HGF route unlearned; both assumptions retained. No HGF or Wilson outcomes fitted.')
    all_results[[paste(i, route)]] <- out
  }
  data <- rbindlist(all_results)[version == 103L + i]
  fwrite(data, file.path(root, 'data/au565_hgf_5uM_zero_reference_query.tsv'), sep = '\t')
}
print(rbindlist(all_results)[hgf_context == 'Equal-context mean',
  .(version, assumed_HGF_route, lapatinib_pct_zero_control, hgf_alone_pct_zero_control,
    combination_pct_zero_control, hgf_addition_relative_pct, interaction_relative_pct)])
