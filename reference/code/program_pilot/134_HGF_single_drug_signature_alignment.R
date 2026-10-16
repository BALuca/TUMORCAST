#!/usr/bin/env Rscript
# Descriptive alignment against the already-fixed v106 single-drug signatures.
# The |r| >= 0.5 display subset is not used to fit or select any model.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs)})
root <- '/work/results/program_pilot/134_HGF_g_diagnostics/linear/data'
model <- qread('/work/results/program_pilot/106_lapatinib_supported_dose_interactions/linear/models/final.qs')
values <- fread(file.path(root, 'AU565_program_level_opposition.tsv'))
R <- model$template$correlations
primary <- values[source_profile == 'Primary_study_mean']
primary <- primary[match(rownames(R), program)]
stopifnot(identical(primary$program, rownames(R)))
primary[, r_with_lapatinib_sensitivity := R[, 'Lapatinib ditosylate']]
primary[, HGF_opposes_sensitivity := HGF_displacement * r_with_lapatinib_sensitivity < 0]
primary[, strong_association_display := abs(r_with_lapatinib_sensitivity) >= .5]
fwrite(primary, file.path(root, 'HGF_vs_lapatinib_sensitivity_programs.tsv'), sep = '\t')
alignment <- data.table(drug = colnames(R),
  r_across_programs_HGF_vs_sensitivity_pattern = as.numeric(cor(primary$HGF_displacement, R)),
  qualification = 'Across 40 program coordinates, not a held-out viability prediction correlation')
fwrite(alignment, file.path(root, 'HGF_alignment_with_individual_drug_signatures.tsv'), sep = '\t')
sources <- values[!grepl('Primary|Leave_out', source_profile)]
sources[, r_with_lapatinib_sensitivity := R[program, 'Lapatinib ditosylate']]
agreement <- sources[, .(sources_opposing_sensitivity = sum(HGF_displacement * r_with_lapatinib_sensitivity < 0),
                         total_sources = .N), by = program]
fwrite(agreement, file.path(root, 'HGF_signature_direction_source_agreement.tsv'), sep = '\t')
collapsed <- data.table(signature = colnames(model$template$collapsed),
  r_across_programs_HGF_vs_sensitivity_pattern = as.numeric(cor(primary$HGF_displacement, model$template$collapsed)))
fwrite(collapsed, file.path(root, 'HGF_alignment_with_collapsed_signatures.tsv'), sep = '\t')
print(primary[strong_association_display == TRUE,
  .(program, r_with_lapatinib_sensitivity, lapatinib_displacement, HGF_displacement, HGF_opposes_sensitivity)])
print(alignment)
print(collapsed)
