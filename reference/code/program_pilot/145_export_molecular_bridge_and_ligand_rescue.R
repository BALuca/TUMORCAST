#!/usr/bin/env Rscript
# Export the frozen v106 single-drug bridge and v144 rescue without changing fits.
# The only edits to runnable source copies are local paths and qs -> base R RDS I/O.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(qs); library(digest) })
setDTthreads(1L)
programs <- sprintf('P%02d', 1:40)

bundle <- '/work/results/program_pilot/145_molecular_bridge_and_ligand_rescue_portable'
support <- '/work/src/program_pilot/145_portable_bundle'
v106 <- '/work/results/program_pilot/106_lapatinib_supported_dose_interactions'
v144 <- '/work/results/program_pilot/144_bounded_program_rescue'
v95 <- '/work/results/program_pilot/95_leave_one_cell_out_correlation_signatures'
v93 <- '/work/results/program_pilot/93_joint_calibrated_sensitivity_routes/data'
v136 <- '/work/results/program_pilot/136_programwise_reversal_reinforcement/linear'
v131 <- '/work/results/program_pilot/131_precision_weighted_HGF_bridge/linear/data'
archive <- '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/inputs/original/modular_pair_pilot/inputs.rds'
src106 <- '/work/src/program_pilot/106_lapatinib_supported_dose_interactions.R'
src144 <- '/work/src/program_pilot/144_bounded_program_rescue.R'
src138 <- '/work/src/program_pilot/138_pooled_direction_constrained_rescue.R'
wilson <- '/work/almanac_combination_training_repair_2026-10-03/results/wilson_observed_source.tsv'
stopifnot(!dir.exists(bundle))
dir.create(bundle, recursive = TRUE)
provenance <- list()

# Copy files, retaining original source hashes for every copied artifact.
copy_file <- function(source, relative) {
  stopifnot(file.exists(source), !dir.exists(source))
  target <- file.path(bundle, relative)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(source, target, overwrite = FALSE))
  original_hash <- digest(source, file = TRUE, algo = 'sha256')
  stopifnot(identical(original_hash, digest(target, file = TRUE, algo = 'sha256')))
  provenance[[relative]] <<- data.table(relative_path = relative, original_path = source,
                                       original_sha256 = original_hash)
}
copy_tree <- function(source, relative) {
  paths <- list.files(source, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  paths <- paths[!file.info(paths)$isdir]
  stopifnot(length(paths) > 0)
  for (path in paths) copy_file(path, file.path(relative, substring(path, nchar(source) + 2)))
}
copy_tree(v106, 'results/single_drug_bridge')
copy_tree(v144, 'results/ligand_rescue')
copy_tree(file.path(v95, 'models'), 'inputs/frozen_response_templates/models')
copy_file(file.path(v95, 'data/recipient_folds.tsv'), 'inputs/frozen_response_templates/data/recipient_folds.tsv')
for (name in c('selected_doses.tsv', 'individual_training_outcomes.tsv', 'combination_training_outcomes.tsv')) {
  copy_file(file.path(v93, name), file.path('inputs/training_outcomes', name))
}
copy_file(archive, 'inputs/molecular_functional_inputs/inputs.rds')
copy_file('/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv',
          'inputs/molecular_functional_inputs/hgf_transported_inputs.tsv')
copy_file(file.path(v136, 'models/final.qs'), 'inputs/ligand_response_inputs/models/final.qs')
for (name in c('molecular_contexts.tsv', 'external_combination_outcomes.tsv', 'external_ligand_programs.tsv',
               'time_matched_HGF_program_inputs.tsv', 'HGF_time_context_eligibility.tsv',
               'time_matched_HGF_source_hashes.tsv', 'external_source_hashes.tsv')) {
  copy_file(file.path(v136, 'data', name), file.path('inputs/ligand_response_inputs/data', name))
}
copy_file(file.path(v136, 'data/single_signatures_beside_combination_benefit/AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv'),
          'inputs/ligand_response_inputs/data/AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv')
for (name in c('HGF_input_reconstruction.qs', 'HGF_source_platform_precision_weights.tsv',
               'corrected_HGF_source_gene_response.tsv', 'HGF_source_reconstruction_hashes.tsv')) {
  copy_file(file.path(v131, name), file.path('inputs/HGF_gene_sources', name))
}
copy_file(wilson, 'inputs/wilson/wilson_observed_source.tsv')
for (name in c('Wilson_original_Fig1.jpg', 'Wilson_original_Fig3.jpg', 'bioc.xml')) {
  copy_file(file.path('/work/tmp/wilson_primary_check', name), file.path('inputs/wilson', name))
}
copy_file('/work/TUMORCAST_grant_evidence_current/data/wilson_SD_source.md',
          'provenance/wilson_original_SD_source.txt')
copy_file('/work/almanac_conditionals/state_conditioned_bridge/program_first/README.md',
          'provenance/program_bank_original_notes.txt')
for (name in c('program_bank_freeze.json', 'program_quality.tsv', 'program_gene_summary.tsv',
               'final_program_stability.tsv', 'program_theme_top_terms.tsv')) {
  copy_file(file.path('/work/almanac_conditionals/state_conditioned_bridge/program_first/results', name),
            file.path('provenance/program_bank', name))
}
copy_file('/work/results/program_pilot/141_audit_v106_to_v140/data/lapatinib_matched_population_correlations.tsv',
          'provenance/lapatinib_matched_population_correlations.tsv')
for (part in c('B1', 'B2')) {
  copy_file(paste0('/work/ERC_StG_2027_', part, '_v22.docx'), paste0('grant_context/ERC_StG_2027_', part, '_v22.docx'))
}

sources <- c(src106, src144, src138, '/work/src/program_pilot/144_plot_grant_support.R',
             '/work/src/program_pilot/95_leave_one_cell_out_correlation_signatures.R',
             '/work/src/program_pilot/131_precision_weighted_HGF_bridge.R',
             '/work/src/program_pilot/131_prepare_precision_weighted_hgf.R')
for (path in sources) copy_file(path, file.path('provenance/original_code', basename(path)))
copy_file('/work/src/R/lib_heatmap.R', 'code/lib_heatmap.R')
for (path in list.files(support, full.names = TRUE)) copy_file(path, file.path('code', basename(path)))

# RDS copies preserve numerical objects and remove qs as a runtime requirement.
objects <- list.files(bundle, pattern = '\\.qs$', full.names = TRUE, recursive = TRUE)
for (path in objects) saveRDS(qread(path), sub('\\.qs$', '.rds', path), compress = 'gzip')

# Extract existing function bodies with their original comments; do not redefine their mathematics.
extract_functions <- function(source, wanted, relative) {
  lines <- readLines(source, warn = FALSE)
  expr <- parse(source, keep.source = TRUE)
  refs <- attr(expr, 'srcref')
  found <- character()
  output <- c('# Frozen functions copied from the original numbered script; see provenance/original_code.', '')
  for (i in seq_along(expr)) {
    e <- expr[[i]]
    if (is.call(e) && identical(e[[1]], as.name('<-')) && is.symbol(e[[2]]) &&
        as.character(e[[2]]) %in% wanted) {
      found <- c(found, as.character(e[[2]]))
      r <- as.integer(refs[[i]])
      output <- c(output, lines[seq.int(r[1], r[3])], '')
    }
  }
  stopifnot(setequal(found, wanted))
  output <- gsub('v106', 'Single-drug bridge', output, fixed = TRUE)
  output <- gsub('v144', 'Ligand rescue', output, fixed = TRUE)
  writeLines(output, file.path(bundle, relative))
}
extract_functions(src106, c('pearson', 'context_key', 'context_index', 'predict_rna', 'signature_node',
                            'pair_features', 'predict_single', 'predict_pair', 'render_heatmap'),
                  'code/single_drug_prediction_functions.R')
extract_functions(src138, c('context_key', 'context_index', 'predict_rna', 'signature_node', 'predict_single'),
                  'code/ligand_query_prediction_functions.R')
extract_functions(src144, c('make_features', 'compose', 'predict_rho'), 'code/bounded_rescue_functions.R')

bootstrap <- c("entry <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])",
               "source(file.path(dirname(normalizePath(entry)), 'config.R'))")
replace_path <- function(lines, old, relative) {
  replacement <- paste0("file.path(bundle_root, '", relative, "')")
  lines <- gsub(paste0("'", old, "'"), replacement, lines, fixed = TRUE)
  gsub(paste0('"', old, '"'), replacement, lines, fixed = TRUE)
}
portable_script <- function(source, replacements, relative) {
  lines <- readLines(source, warn = FALSE)
  lib <- grep('^\\.libPaths', lines)
  stopifnot(length(lib) == 1)
  lines <- append(lines[-lib], bootstrap, after = lib - 1)
  for (old in names(replacements)) lines <- replace_path(lines, old, replacements[[old]])
  lines <- gsub('library(qs);', '', lines, fixed = TRUE)
  lines <- gsub('qread(', 'readRDS(', lines, fixed = TRUE)
  lines <- gsub('qsave(', 'saveRDS(', lines, fixed = TRUE)
  lines <- gsub('.qs', '.rds', lines, fixed = TRUE)
  # Public script/figure names describe the model rather than an iteration number.
  lines <- gsub('v106_v144_grant_support', 'molecular_bridge_and_ligand_rescue_support', lines, fixed = TRUE)
  lines <- gsub('v106', 'Single-drug bridge', lines, fixed = TRUE)
  lines <- gsub('v144', 'Ligand rescue', lines, fixed = TRUE)
  stopifnot(!any(grepl('/work/|/root/', lines)))
  writeLines(lines, file.path(bundle, relative))
  invisible(parse(file.path(bundle, relative)))
}
portable_script(src106, setNames(c('regenerated/single_drug_bridge', 'inputs/frozen_response_templates', 'inputs/training_outcomes',
                                  'inputs/molecular_functional_inputs/inputs.rds', 'code/lib_heatmap.R',
                                  'inputs/molecular_functional_inputs/hgf_transported_inputs.tsv', 'code/refit_single_drug_bridge.R'),
                                c(v106, v95, v93, archive, '/work/src/R/lib_heatmap.R',
                                  '/work/results/program_pilot/74_v20_lapatinib_hgf_frozen_query/data/hgf_transported_inputs.tsv',
                                  src106)), 'code/refit_single_drug_bridge.R')
portable_script(src144, setNames(c('regenerated/ligand_rescue', 'inputs/ligand_response_inputs',
                                  'results/single_drug_bridge/linear/models/final.qs', 'inputs/molecular_functional_inputs/inputs.rds',
                                  'inputs/training_outcomes/selected_doses.tsv', 'code/ligand_query_prediction_functions.R',
                                  'inputs/wilson/wilson_observed_source.tsv', 'code/refit_ligand_rescue.R'),
                                c(v144, v136, file.path(v106, 'linear/models/final.qs'), archive,
                                  file.path(v93, 'selected_doses.tsv'), src138, wilson, src144)),
                'code/refit_ligand_rescue.R')
plot_source <- '/work/src/program_pilot/144_plot_grant_support.R'
portable_script(plot_source, setNames(c('regenerated/ligand_rescue',
  'results/single_drug_bridge/linear/data/baseline_predictions_frozen.tsv',
  'results/single_drug_bridge/linear/data/measured_predictions_frozen.tsv',
  'inputs/ligand_response_inputs/data/AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv',
  'code/recreate_existing_grant_figure.R'), c(v144,
  file.path(v106, 'linear/data/baseline_predictions_frozen.tsv'),
  file.path(v106, 'linear/data/measured_predictions_frozen.tsv'),
  file.path(v136, 'data/single_signatures_beside_combination_benefit/AU565_lapatinib_5uM_vs_HGF_raw_displacements.tsv'),
  plot_source)), 'code/recreate_existing_grant_figure.R')

lines <- readLines(src106)
first <- grep('^plot_single <- function', lines)
last <- grep("^write_table\\(data.table\\(assay=", lines) - 1L
stopifnot(length(first) == 1, length(last) == 1)
figure_lines <- gsub('v106', 'Single-drug bridge', lines[first:last], fixed = TRUE)
writeLines(figure_lines, file.path(bundle, 'code/single_drug_figures_from_tables.R'))
lines <- readLines(src144)
first <- grep('^annotation <- metrics', lines)
last <- grep('^# Provenance and modeling assumptions', lines) - 1L
lines <- lines[first:last]
lines <- replace_path(lines, wilson, 'inputs/wilson/wilson_observed_source.tsv')
writeLines(lines, file.path(bundle, 'code/ligand_rescue_figures_from_tables.R'))

# Inspectable TSV inputs supplement R objects for other agents and languages.
a <- readRDS(file.path(bundle, 'inputs/molecular_functional_inputs/inputs.rds'))
tabular <- file.path(bundle, 'inputs/molecular_functional_inputs/tabular')
dir.create(tabular)
for (name in c('m', 'raw', 'receptors')) fwrite(as.data.table(a$prep[[name]]), file.path(tabular, paste0(name, '.tsv')), sep = '\t')
fwrite(data.table(gene_symbol = rownames(a$prep$W), a$prep$W), file.path(tabular, 'gene_program_weights.tsv.gz'), sep = '\t')
m <- readRDS(file.path(bundle, 'results/single_drug_bridge/linear/models/final.rds'))
fwrite(data.table(program = programs, m$context$W), file.path(tabular, 'functional_signature_weights.tsv'), sep = '\t')
# Only the recipient rows required to reproduce frozen HGF projection are exported.
hgf <- fread(file.path(bundle, 'inputs/ligand_response_inputs/data/time_matched_HGF_program_inputs.tsv'))
ids_wanted <- unique(hgf$depmap_id)
connection <- file('/work/outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv', 'r')
header <- readLines(connection, 1)
expression_rows <- list()
repeat {
  z <- readLines(connection, 100, warn = FALSE)
  if (!length(z)) break
  ids <- gsub('"', '', sub(',.*', '', z))
  for (i in which(ids %in% ids_wanted)) expression_rows[[ids[i]]] <- z[i]
  if (length(expression_rows) == length(ids_wanted)) break
}
close(connection)
stopifnot(setequal(names(expression_rows), ids_wanted))
expression <- fread(text = paste(c(header, unlist(expression_rows, use.names = FALSE)), collapse = '\n'))
setnames(expression, 1, 'depmap_id')
fwrite(expression, file.path(bundle, 'inputs/HGF_gene_sources/recipient_expression.tsv.gz'), sep = '\t')
source_urls <- data.table(resource = c('Wilson Figure 1: 72h AU565 curve', 'Wilson Figure 3: long-term staining',
                                     'Wilson main paper, local BioC XML', 'NRG1 calibration study identifier'),
                         url = c('https://www.nature.com/articles/nature11249/figures/1',
                                 'https://www.nature.com/articles/nature11249/figures/3',
                                 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3724525/',
                                 'https://elifesciences.org/articles/32271'))
fwrite(source_urls, file.path(bundle, 'provenance/source_urls.tsv'), sep = '\t')
fwrite(data.table(CellLine = 'AU565', dose_uM = 5, duration_h = 72,
                  condition = c('Lapatinib', 'Lapatinib + HGF'), approximate_viability_pct = c(13, 33),
                  measurement_status = 'Approximate visual reading of Figure 1b; not extracted raw observations; no SD',
                  use = 'Context only; not used for fitting, selection, prediction replay, or historical plots'),
       file.path(bundle, 'inputs/wilson/72h_curve_visual_estimates.tsv'), sep = '\t')

packages <- c('data.table', 'ggplot2', 'glmnet', 'digest', 'nloptr', 'RhpcBLASctl', 'ComplexHeatmap',
              'circlize', 'dplyr', 'matrixStats', 'RColorBrewer', 'viridisLite', 'patchwork', 'qs')
package_info <- rbindlist(lapply(packages, function(p) {
  d <- packageDescription(p)
  data.table(package = p, version = as.character(d$Version), license = as.character(d$License))
}))
fwrite(package_info, file.path(bundle, 'provenance/package_versions.tsv'), sep = '\t')
writeLines(capture.output(sessionInfo()), file.path(bundle, 'provenance/export_session_info.txt'))
fwrite(rbindlist(provenance), file.path(bundle, 'provenance/original_sources.tsv'), sep = '\t')
writeLines(c(
  'TUMORCAST: molecular single-drug bridge and bounded ligand rescue',
  '',
  'START HERE',
  'Use R. Unpack the archive; commands work from any directory when the script path is absolute.',
  'Rscript --vanilla code/01_verify_and_replay.R',
  'Rscript --vanilla code/02_regenerate_figures.R',
  'Rscript --vanilla code/03_rebuild_HGF_program_inputs.R',
  'Paths above assume the current directory is the unpacked folder.',
  'Existing R packages must be available; package versions/licenses are recorded in provenance/.',
  'Scripts do not install packages, download data, or refer to the original server.',
  '',
  'CONTENTS',
  'results/single_drug_bridge/: original tables, figures, 62 recipient-excluded fits and final model.',
  'results/ligand_rescue/: original rescue tables, figures, fitted calibration and recipient folds.',
  'inputs/molecular_functional_inputs/: complete processed input archive; readable baseline, receptor, RNA and gene-weight TSVs.',
  'inputs/training_outcomes/ and inputs/frozen_response_templates/: exact dose/outcome inputs and frozen fold-specific upstream templates.',
  'inputs/ligand_response_inputs/: exact Ligand rescue molecular inputs, ligand outcomes and frozen v136 single-drug object.',
  'inputs/HGF_gene_sources/: source gene responses, platform precision and the 14 recipient expression profiles.',
  'inputs/wilson/: original figures, main-paper XML, historical printed long-term values and approximate 72h readings.',
  'code/: runnable portable copies, exact prediction functions, table-only plotting and optional refits.',
  'provenance/: original code snapshots, source paths/hashes, package versions and program-bank notes.',
  'grant_context/: unedited B1/B2 v22 documents for context; their older numerical claims are not updated.',
  'regenerated/: new outputs. Original results are never overwritten.',
  '',
  'OPTIONAL REFITS',
  'Rscript --vanilla code/refit_single_drug_bridge.R',
  'Rscript --vanilla code/refit_ligand_rescue.R',
  'These can take longer. Ligand rescue deliberately keeps the packaged frozen single-drug model.',
  'Single-drug bridge starts from the packaged fold-specific v95 signatures/RNA templates, as the original Single-drug bridge did.',
  'The molecular gene-program bank is frozen. Rebuilding discovery from raw Tahoe/CytoSig is outside this export.',
  'The export includes all processed inputs needed for these fits; full original public screens are not needed.',
  'Portable refit edits: paths and qs/RDS serialization only. Original implementation snapshots are preserved.',
  '',
  'IMPORTANT SCIENTIFIC CONTEXT FOR THE NEXT AGENT',
  'Single-drug bridge: baseline -> predicted 40-program displacement -> shared calibrated functional decoder.',
  'Lapatinib molecular conditioning uses untreated ERBB2/HER2 and MET; no direct baseline input to function.',
  'Eight functional drugs/doses; Crizotinib and Vismodegib excluded. The RNA bank retains ten drugs.',
  'Held-recipient lapatinib Pearson r: .5864 on 30 RNA-overlap cells, .3421 across 441 PRISM cells, fixed 5 uM.',
  'Those correlations use log2 control-relative abundance. Raw-percent correlations differ.',
  'Recipient exclusions are retrospective: dose selection and repeated inspection used revealed outcomes.',
  'Training across drugs is not proof of transfer to an unseen drug; other end-to-end drug results are weak.',
  'Single-drug bridge original ALMANAC interaction remains archived, but Ligand rescue replaces that rescue component.',
  'Ligand rescue: seven frozen lapatinib-associated programs; normalized ligand opposition -> bounded recovery fraction.',
  'Calibration uses 12 NRG1 observations, six cells, 0.5/5 uM, 72h; four EGF observations are external controls.',
  'HGF/Wilson outcomes do not enter calibration fitting. Repeated HGF-directed development makes this retrospective.',
  'The query assumes HGF-alone equals 100% and transports HGF profiles from other cancer cells into AU565.',
  'At 5 uM: 36.30997% lapatinib-alone -> 60.24105% with the mean HGF profile; 37.5743% depletion restored.',
  'At 0.5 uM: 90.13471% -> 99.91352%. This is NOT a measured five-day HGF result.',
  'A 72h-derived recovery fraction is applied to the five-day PRISM prediction without a learned time course.',
  'Conditional interaction evaluation uses measured drug/ligand controls; end-to-end diagnostics are separate.',
  'The constant-recovery comparator matches/exceeds conditional correlations; molecular added value is unresolved.',
  'HGF source sensitivity is substantial: individual source predictions can go in the opposite direction.',
  'Wilson 72h Figure 1b: partial AU565 rescue, approximately 13% -> 33% at 5 uM; readings are approximate.',
  'Wilson long-term Figure 3c: 1 uM; treated wells day 12, untreated image day 6. Control normalization unresolved.',
  'Existing Wilson comparison figures are historical illustrations, NOT dose/time-matched validation.',
  'Do not infer functional rescue merely from opposite molecular signatures.',
  'No new model architectures, fit selection, or grant panel arrangements were created by this export.',
  '',
  'VERIFICATION',
  '01 verifies SHA256 for all exported artifacts and replays baseline, measured-RNA, ALMANAC and rescue predictions.',
  'All replayed predictions must agree with the original tables to absolute tolerance 1e-10.',
  'The replay also exports all held-cell predicted displacements for further analysis.',
  '02 redraws the existing scatterplots, fold/signature heatmaps and rescue figures from frozen artifacts.',
  '03 rebuilds HGF program inputs from source gene responses and recipient RNA and checks the frozen inputs.'
), file.path(bundle, 'README.txt'))

# A column dictionary allows another agent to inspect the data without loading models.
tables <- list.files(bundle, pattern = '\\.tsv(\\.gz)?$', full.names = TRUE, recursive = TRUE)
dictionary <- rbindlist(lapply(tables, function(path) {
  x <- fread(path)
  data.table(relative_path = substring(path, nchar(bundle) + 2), rows = nrow(x),
             column = names(x), type = vapply(x, function(v) class(v)[1], character(1)))
}))
fwrite(dictionary, file.path(bundle, 'provenance/table_dictionary.tsv'), sep = '\t')
code_files <- list.files(file.path(bundle, 'code'), pattern = '\\.R$', full.names = TRUE)
for (path in code_files) {
  invisible(parse(path))
  stopifnot(!any(grepl('/work/|/root/', readLines(path))))
}
files <- list.files(bundle, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
files <- files[!file.info(files)$isdir]
manifest <- data.table(relative_path = substring(files, nchar(bundle) + 2), bytes = file.info(files)$size,
                       sha256 = vapply(files, digest, character(1), file = TRUE, algo = 'sha256'))
fwrite(manifest, file.path(bundle, 'provenance/file_manifest.tsv'), sep = '\t')
cat('Exported', nrow(manifest), 'artifacts;', round(sum(manifest$bytes) / 1024^2, 1), 'MiB.\n')
cat(bundle, '\n')
