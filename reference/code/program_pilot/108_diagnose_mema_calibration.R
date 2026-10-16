#!/usr/bin/env Rscript
# Diagnose calibration jitter without changing or refitting v108 predictions.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- '/work/results/program_pilot/108_mema_calibrated_shared_interaction/linear/data'
pred <- fread(file.path(root, 'mema_predictions_scored.tsv'))
audit <- fread(file.path(root, 'fit_audit.tsv'))
final_scale <- audit[fit == 'final', .(endpoint, final_scale = mema_assay_scale)]
held <- merge(pred[evaluation == 'held_gene'], final_scale, by = 'endpoint')
# Constant full-data scale is an attribution diagnostic, NOT a validated predictor.
held[, diagnostic_constant_scale := final_scale * (uncalibrated_drug_log2 + interaction_log2)]
held[, diagnostic_scale_only := mema_assay_scale * uncalibrated_drug_log2]
summary <- held[, .(
  ligands = .N,
  held_r = cor(observed_log2, predicted_log2),
  diagnostic_constant_scale_r = cor(observed_log2, diagnostic_constant_scale),
  diagnostic_scale_only_r = cor(observed_log2, diagnostic_scale_only),
  held_scale_min = min(mema_assay_scale),
  held_scale_max = max(mema_assay_scale),
  molecular_interaction_SD = sd(final_scale * interaction_log2),
  scale_only_SD = sd(diagnostic_scale_only),
  observed_SD = sd(observed_log2)
), by = .(endpoint, CellLine)]
fwrite(summary, file.path(root, 'calibration_jitter_diagnostic.tsv'), sep = '\t')
fwrite(held, file.path(root, 'calibration_jitter_attribution.tsv'), sep = '\t')
print(summary)
