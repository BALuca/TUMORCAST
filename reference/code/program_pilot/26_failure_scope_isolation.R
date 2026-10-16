#!/usr/bin/env Rscript
# Retrospective scope checks on frozen predictions; no fitting or model selection.
suppressPackageStartupMessages(library(data.table))
root <- '/work/results/program_pilot/26_failure_scope_isolation'
dir.create(root, recursive = TRUE, showWarnings = FALSE)
paths <- paste0('/work/results/program_pilot/24_ixazomib_decoder_capacity_and_molecular_signal_',
  c('MEK', 'ERBB'), '_holdout/data/shared_capacity_predictions.tsv.gz')
values <- rbindlist(lapply(seq_along(paths), function(i) {
  x <- fread(paths[i])
  x[, family := c('MEK', 'ERBB')[i]]
  x
}))
stopifnot(!anyDuplicated(values[, .(family, model, cohort, observation_id)]),
  all(is.finite(values$observed_pct)), all(values$observed_pct > 0),
  all(is.finite(values$predicted_pct)), all(values$predicted_pct > 0))
values[, `:=`(observed = log2(observed_pct / 100), predicted = log2(predicted_pct / 100))]
association <- function(a, b) {
  if (length(a) < 3L || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  cor(a, b)
}
summarize_scope <- function(x) {
  x <- copy(x)
  x[, `:=`(observed_centered = observed - mean(observed),
    predicted_centered = predicted - mean(predicted),
    observed_pct_centered = observed_pct - mean(observed_pct),
    predicted_pct_centered = predicted_pct - mean(predicted_pct)), by = .(family, model, cohort, drug, dose)]
  x[, .(observations = .N, cells = uniqueN(depmap_id), drugs = uniqueN(drug),
    overall_r = association(observed, predicted),
    within_treatment_r = association(observed_centered, predicted_centered),
    within_treatment_pct_r = association(observed_pct_centered, predicted_pct_centered),
    sensitive = sum(observed_pct <= 50), sensitive_called = sum(observed_pct <= 50 & predicted_pct <= 50),
    resistant = sum(observed_pct > 50), false_sensitive = sum(observed_pct > 50 & predicted_pct <= 50)),
    by = .(family, model, cohort)]
}
scopes <- list(All = values, Without_ixazomib = values[drug != 'ixazomib'],
  Without_six_traced_observations = values[!((drug == 'ixazomib' & dose == .05 &
    recipient %in% c('CFPAC-1', 'SW 900', 'NCI-H2347', 'BT-474')) |
    (drug == 'ixazomib' & dose == .5 & recipient == 'CFPAC-1') |
    (drug == 'docetaxel' & dose == .5 & recipient == 'Hs 766T'))])
scope_summary <- rbindlist(lapply(scopes, summarize_scope), idcol = 'scope')
fwrite(scope_summary, file.path(root, 'scope_summary.tsv'), sep = '\t')
held <- values[cohort == 'Held cells']
leave_cell <- rbindlist(lapply(unique(held$recipient), function(cell) {
  z <- summarize_scope(held[recipient != cell])
  z[, omitted_cell := cell]
  z
}))
fwrite(leave_cell, file.path(root, 'leave_one_cell_out.tsv'), sep = '\t')
by_treatment <- held[, .(cells = .N, sensitive = sum(observed_pct <= 50),
  resistant = sum(observed_pct > 50), observed_range_pct = diff(range(observed_pct)),
  r = association(observed, predicted)), by = .(family, model, drug, dose)]
fwrite(by_treatment, file.path(root, 'per_treatment.tsv'), sep = '\t')
by_cell <- held[, .(observations = .N, across_treatment_r = association(observed, predicted)),
  by = .(family, model, recipient)]
fwrite(by_cell, file.path(root, 'per_cell.tsv'), sep = '\t')
treatment_summary <- by_treatment[sensitive > 0 & resistant > 0,
  .(mixed_treatments = .N, positive_r = sum(r > 0, na.rm = TRUE),
    negative_r = sum(r < 0, na.rm = TRUE), median_r = median(r, na.rm = TRUE)), by = .(family, model)]
fwrite(treatment_summary, file.path(root, 'mixed_treatment_summary.tsv'), sep = '\t')
source_paths <- c(paths, '/work/src/program_pilot/26_failure_scope_isolation.R')
fwrite(data.table(path = source_paths, sha256 = vapply(source_paths, digest::digest,
  character(1), file = TRUE, algo = 'sha256')), file.path(root, 'input_hashes.tsv'), sep = '\t')
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
print(scope_summary[cohort == 'Held cells' & model != 'Very local SVR'])
print(leave_cell[, .(min_r = min(within_treatment_r), max_r = max(within_treatment_r)), by = .(family, model)])
print(treatment_summary)
print(by_cell[, .(min_r = min(across_treatment_r), max_r = max(across_treatment_r)), by = .(family, model)])
