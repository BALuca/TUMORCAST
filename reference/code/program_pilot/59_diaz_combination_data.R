#!/usr/bin/env Rscript
# Matched single/combo RNA and independent functional time courses, Diaz2020.
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(jsonlite)
  library(curl)
  library(qs)
  library(Matrix)
  library(digest)
})
root <- '/work/results/program_pilot/59_mechanistic_evidence_datasets'
for (folder in c('source', 'data')) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
setDTthreads(2L)
options(timeout = 180)
# Use exact public primary-source supplementary files. A failed download stops
# this adapter; HTML challenge pages must never be parsed as measurements.
download_manifest <- list()
download_source <- function(url, name) {
  path <- file.path(root, 'source', paste0('diaz_', name))
  if (!file.exists(path)) curl_download(url, path, quiet = TRUE)
  stopifnot(file.exists(path), file.info(path)$size > 0)
  if (grepl('\\.(zip|xlsx)$', path)) {
    con <- file(path, 'rb')
    magic <- readBin(con, 'raw', n = 2L)
    close(con)
    stopifnot(identical(magic, charToRaw('PK')))
  }
  download_manifest[[name]] <<- data.table(study_id = 'DIAZ2020', url = url, path = path,
    bytes = file.info(path)$size, sha256 = digest(file = path, algo = 'sha256'),
    downloaded_utc = format(file.info(path)$mtime, tz = 'UTC', usetz = TRUE))
  fwrite(rbindlist(download_manifest), file.path(root, 'data/diaz_download_manifest.tsv'), sep = '\t')
  path
}
for (number in c(1L, 2L, 3L, 4L, 5L)) {
  name <- paste0('elife-52707-supp', number, '-v2.xlsx')
  download_source(paste0('https://cdn.elifesciences.org/articles/52707/', name), name)
}
for (number in c(6L, 7L, 8L)) {
  name <- if (number == 6L) 'elife-52707-supp6-v2.csv.zip' else paste0('elife-52707-supp', number, '-v2.zip')
  download_source(paste0('https://cdn.elifesciences.org/articles/52707/', name), name)
}
fwrite(rbindlist(download_manifest), file.path(root, 'data/diaz_download_manifest.tsv'), sep = '\t')
save_table <- function(x, name) fwrite(x, file.path(root, 'data', paste0('diaz_', name, '.tsv.gz')), sep = '\t')
archive <- file.path(root, 'source/diaz_elife-52707-supp8-v2.zip')
entries <- unzip(archive, list = TRUE)$Name
needed <- entries[grepl('(dynsyn_count.csv|lncap_anno.csv|lncap_ucsc_raw_counts|pheno_ss.txt|dose_ss.txt)$', entries) &
  !grepl('__MACOSX', entries)]
extracted <- unzip(archive, files = needed, exdir = file.path(root, 'source/diaz_rna'))
source_file <- function(name) {
  hit <- extracted[basename(extracted) == name]
  stopifnot(length(hit) == 1L)
  hit
}

# Source matrices have gene symbols. LNCaP sample aliases map to sequenced run IDs
# through the deposited annotation; do not infer treatment from column order.
mcf <- fread(source_file('dynsyn_count.csv'))
setnames(mcf, 1L, 'gene')
stopifnot(!anyDuplicated(mcf$gene))
mcf_matrix <- as.matrix(mcf[, -'gene'])
rownames(mcf_matrix) <- mcf$gene
ln <- fread(source_file('lncap_ucsc_raw_counts'), skip = 1L)
ln_annotation <- fread(source_file('lncap_anno.csv'))
setnames(ln_annotation, 1L, 'source_run_alias')
rna_cols <- names(ln)[-(1:6)]
run_aliases <- sub('.*star_align/([^_]+)_.*', '\\1', rna_cols)
stopifnot(!anyDuplicated(run_aliases), setequal(run_aliases, ln_annotation$source_run_alias))
ln_matrix <- as.matrix(ln[, ..rna_cols])
rownames(ln_matrix) <- ln$Geneid
colnames(ln_matrix) <- ln_annotation$Sample[match(run_aliases, ln_annotation$source_run_alias)]
stopifnot(!anyDuplicated(rownames(ln_matrix)), !anyDuplicated(colnames(ln_matrix)))
counts <- list(MCF7 = mcf_matrix, LNCaP = ln_matrix)
sample_tables <- list()
expression <- list()
bank <- qread('/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs/models/program_bank.qs')
coverage_tables <- score_tables <- list()
for (cell in names(counts)) {
  x <- counts[[cell]]
  stopifnot(all(is.finite(x)), all(x >= 0), all(x == round(x)))
  source_ids <- colnames(x)
  stopifnot(all(grepl('^(DMSO|T|M|W|TM|TW|MW|TMW)_[0-9]+[ABC]$', source_ids)))
  md <- data.table(study_id = 'DIAZ2020', cell_line = cell, source_sample_id = source_ids,
    sample_id = paste('DIAZ2020', cell, source_ids, sep = '__'),
    arm = sub('_.*', '', source_ids), time_h = as.numeric(sub('.*_([0-9]+)[ABC]$', '\\1', source_ids)),
    replicate = sub('.*([ABC])$', '\\1', source_ids), assay = 'RNAseq',
    replicate_type = 'Reported triplicate; biological independence not resolved from sample alias',
    tamoxifen_uM = fifelse(grepl('T', sub('_.*', '', source_ids)), 20, 0),
    mefloquine_uM = fifelse(grepl('M', sub('_.*', '', source_ids)) & !grepl('DMSO', source_ids), 10, 0),
    withaferin_A_uM = fifelse(grepl('W', sub('_.*', '', source_ids)), 5, 0),
    vehicle = '0.3% DMSO', source_file = if (cell == 'MCF7') source_file('dynsyn_count.csv') else source_file('lncap_ucsc_raw_counts'))
  md[, condition_id := paste(study_id, cell_line, arm, time_h, sep = '__')]
  stopifnot(md[arm == 'DMSO', all(tamoxifen_uM == 0 & mefloquine_uM == 0 & withaferin_A_uM == 0)])
  colnames(x) <- md$sample_id
  # Save counts independently of a transparent, outcome-blind score projection.
  qsave(x, file.path(root, 'data', paste0('diaz_', cell, '_counts.qs')), preset = 'high')
  dge <- edgeR::DGEList(counts = x)
  dge <- edgeR::calcNormFactors(dge, method = 'TMM')
  log_expression <- log2(1 + edgeR::cpm(dge) / 100)
  delta <- log_expression
  for (t in unique(md$time_h)) {
    controls <- which(md$arm == 'DMSO' & md$time_h == t)
    indices <- which(md$time_h == t)
    stopifnot(length(controls) == 3L)
    delta[, indices] <- sweep(log_expression[, indices, drop = FALSE], 1L,
      rowMeans(log_expression[, controls, drop = FALSE]), '-')
  }
  expression[[cell]] <- list(log2_one_plus_TMM_CP10K = log_expression, displacement = delta, samples = md)
  qsave(expression[[cell]], file.path(root, 'data', paste0('diaz_', cell, '_molecular.qs')), preset = 'high')
  sample_tables[[cell]] <- md
  for (program in colnames(bank$weights)) {
    genes <- bank$genes[as.numeric(bank$weights[, program]) > 0]
    present <- intersect(genes, rownames(delta))
    coverage_tables[[paste(cell, program)]] <- data.table(study_id = 'DIAZ2020', cell_line = cell,
      program = program, program_genes = length(genes), measured_genes = length(present),
      fraction_measured = length(present)/length(genes))
    score_tables[[paste(cell, program)]] <- data.table(sample_id = md$sample_id, program = program,
      displacement = if (length(present)/length(genes) >= .8) colMeans(delta[present, , drop = FALSE]) else NA_real_,
      baseline_reference = 'same-cell same-time DMSO mean',
      score_status = 'Fixed bank43 projection; partial membership renormalized if coverage >=80%')
  }
}
samples <- rbindlist(sample_tables)
stopifnot(!anyDuplicated(samples$sample_id))
save_table(samples, 'samples')
save_table(rbindlist(coverage_tables), 'bank43_coverage')
save_table(rbindlist(score_tables), 'bank43_scores')

# Independent functional source tables contain means and SDs, not individual
# wells. Preserve the published assays separately: ATP, counted nuclei, metabolic
# real-time signal. Do not equate them, manufacture replicate rows or clamp >100%.
drug_code <- c('Tamoxifen'='T', 'Mefloquine'='M', 'Withaferin A'='W', 'Tamoxifen - Mefloquine'='TM')
functional_rows <- list()
pair_rows <- list()
for (cell in c('MCF7', 'LNCaP')) {
  number <- if (cell == 'MCF7') 3L else 4L
  path <- file.path(root, 'source', paste0('diaz_elife-52707-supp', number, '-v2.xlsx'))
  tb <- as.data.table(read_excel(path, sheet = 'Combinations'))
  setnames(tb, names(tb), gsub('[^A-Za-z0-9]+', '_', names(tb)))
  if ('Viability_Assay' %in% names(tb)) setnames(tb, 'Viability_Assay', 'Assay')
  if ('Drug1_StdDev' %in% names(tb)) setnames(tb, 'Drug1_StdDev', 'Drug1_SD')
  combo_sd <- if ('Observed_Inhibition_StdDev' %in% names(tb)) tb$Observed_Inhibition_StdDev else tb$Observed_Viability_StdDev
  arm1 <- unname(drug_code[trimws(tb$Drug1)])
  arm2 <- unname(drug_code[trimws(tb$Drug2)])
  stopifnot(!anyNA(arm1), !anyNA(arm2))
  combo <- paste0(arm1, arm2)
  pair_rows[[cell]] <- data.table(study_id = 'DIAZ2020', cell_line = cell,
    assay = tb$Assay, time_h = tb$Hours, arm_A = arm1, arm_B = arm2, arm_AB = combo,
    A_pct = tb$Drug_1_Viability, B_pct = tb$Drug_2_Viability, AB_pct = tb$Observed_Viability,
    A_SD_pct = tb$Drug1_SD, B_SD_pct = tb$Drug2_StdDev, AB_SD_pct = combo_sd,
    reported_EOB_pct = as.numeric(tb$Excess_Over_BLISS), reported_EOB_error = tb$Excess_Over_BLISS_Error,
    source_file = path, source_row = seq_len(nrow(tb)) + 1L)
  for (role in c('A', 'B', 'AB')) {
    pp <- pair_rows[[cell]]
    functional_rows[[paste(cell, role)]] <- pp[, .(study_id, cell_line, assay, time_h,
      arm = get(paste0('arm_', role)), value = get(paste0(role, '_pct')),
      sd = get(paste0(role, '_SD_pct')), unit = 'percent vehicle signal',
      measurement_type = 'Published mean and SD; individual wells not available in this table',
      source_role = role, source_file, source_row)]
  }
}
func <- rbindlist(functional_rows)
func[, feature := fcase(assay == 'InCell', 'counted_nuclei_relative_to_vehicle',
  assay == 'CTG', 'ATP_luminescence_relative_to_vehicle', assay == 'RTG', 'real_time_metabolic_signal',
  default = 'unrecognized')]
stopifnot(!any(func$feature == 'unrecognized'), all(is.finite(func$value)),
  all(is.na(func$sd) | (is.finite(func$sd) & func$sd >= 0)))
func[, uncertainty_status := fifelse(is.finite(sd), 'Source SD reported',
  'Source SD missing; retained as NA, not imputed')]
# Shared single-agent means recur in pair tables; retain source rows and give
# an explicit link ID instead of counting duplicate table entries as replicates.
func[, condition_id := paste(study_id, cell_line, arm, time_h, sep = '__')]
func[, measurement_id := paste(condition_id, assay, source_row, sep = '__')]
save_table(func, 'measurements')
pairs <- rbindlist(pair_rows)
pairs[, bliss_unclipped_pct := A_pct * B_pct / 100 - AB_pct]
pairs[, source_report_discrepancy_pct := reported_EOB_pct - bliss_unclipped_pct]
save_table(pairs, 'functional_factorials')

# Build joinable RNA->independent-function links without pretending that RNA
# and function were measured in the identical wells or at an identical endpoint.
rna_conditions <- samples[, .(rna_replicates = .N, rna_sample_ids = paste(sample_id, collapse = ';')),
  by = .(study_id, cell_line, arm, time_h, condition_id)]
# TM occurs both as a pair outcome and as a constituent of TMW. These are the
# same mean, sometimes with rounded SD or a missing SD in one display row.
function_unique <- func[, {
  stopifnot(uniqueN(value) == 1L)
  reported <- .SD[is.finite(sd)]
  reported <- reported[order(source_role != 'AB', source_row)]
  sd_values <- unique(reported$sd)
  round_agreement <- length(sd_values) > 0L && diff(range(sd_values)) <= .05
  list(value = value[1], sd = if (round_agreement) reported$sd[1] else NA_real_,
    source_SD_values = paste(sd_values, collapse = ';'),
    uncertainty_status = if (!length(sd_values)) 'SD absent' else if (!round_agreement)
      'Conflicting SD; not combined' else if (length(sd_values) > 1L)
      'Repeated SD differs only by rounding; direct pair-outcome SD preferred' else 'Source SD retained',
    functional_source_rows = paste(source_row, collapse = ';'))
}, by = .(condition_id, study_id, cell_line, arm, time_h, assay, feature)]
stopifnot(!anyDuplicated(function_unique[, .(condition_id, assay)]))
save_table(function_unique, 'functional_condition_summaries')
links <- merge(rna_conditions, function_unique,
  by = c('condition_id','study_id','cell_line','arm','time_h'), all.x = TRUE, allow.cartesian = TRUE)
links[, linkage := fifelse(is.finite(value), 'same study/cell/arm/dose/time; separate assay wells',
  'no reported functional mean; DMSO is reference rather than fabricated numeric data')]
save_table(links, 'same_time_measurement_links')
later <- merge(rna_conditions[time_h <= 12],
  function_unique[time_h >= 24, .(study_id, cell_line, arm, function_time_h = time_h, assay, feature, value, sd)],
  by = c('study_id','cell_line','arm'), allow.cartesian = TRUE)
later[, linkage := 'early RNA to later function, same study/cell/arm/dose; separate assay wells']
save_table(later, 'early_rna_later_function_links')

# Every tested pair/time is retained. LNCaP also contains the complete T/M/W
# triple factorial: reserve TMW as a composition test, never program discovery.
contrast_rows <- list()
for (cell in unique(samples$cell_line)) {
  for (t in sort(unique(samples[cell_line == cell, time_h]))) {
    available <- samples[cell_line == cell & time_h == t, unique(arm)]
    for (pair in list(c('T','M'), c('T','W'), c('M','W'))) {
      arms <- c('DMSO', pair, paste0(pair, collapse = ''))
      contrast_rows[[paste(cell, t, arms[4])]] <- data.table(study_id = 'DIAZ2020', cell_line = cell,
        time_h = t, contrast_id = paste('DIAZ2020',cell,t,arms[4],sep='__'),
        question = 'Does combined RNA contain a response beyond both singles, and does independent function change?',
        required_arms = paste(arms, collapse = ';'), available_arms = paste(intersect(arms, available), collapse = ';'),
        molecular_complete = all(arms %in% available),
        molecular_contrast = paste0(arms[4], ' - ', pair[1], ' - ', pair[2], ' + DMSO'),
        functional_contrast = 'Observed combination versus both measured singles; counted nuclei primary, ATP comparator',
        role = 'Pair-response analysis; no outcome-selected gene definition')
    }
    if (cell == 'LNCaP') {
      arms <- c('DMSO','T','M','W','TM','TW','MW','TMW')
      contrast_rows[[paste(cell,t,'TMW')]] <- data.table(study_id = 'DIAZ2020', cell_line = cell,
        time_h = t, contrast_id = paste('DIAZ2020',cell,t,'TMW',sep='__'),
        question = 'Does a third agent change the pair response beyond its measured single and pair effects?',
        required_arms = paste(arms,collapse=';'),available_arms=paste(intersect(arms,available),collapse=';'),
        molecular_complete = all(arms %in% available),
        molecular_contrast = 'TMW - TM - TW - MW + T + M + W - DMSO',
        functional_contrast = 'Same factorial applied to independent counted nuclei; retain ATP as separate assay',
        role = 'Reserved triple-composition test; no fitting/selection on TMW outcomes')
    }
  }
}
contrasts <- rbindlist(contrast_rows)
contrasts[, functional_mean_coverage := vapply(seq_len(.N), function(i) {
  needed <- setdiff(strsplit(required_arms[i], ';', fixed = TRUE)[[1]], 'DMSO')
  observed <- function_unique[cell_line == contrasts$cell_line[i] & time_h == contrasts$time_h[i] & assay == 'InCell', arm]
  all(needed %in% observed)
}, logical(1))]
contrasts[, `:=`(functional_vehicle = 'Normalized reference100; raw vehicle wells not supplied',
  functional_uncertainty = 'Published means only; some joint SD missing; no replicated factorial significance claim')]
save_table(contrasts, 'contrasts')
study <- data.table(study_id = 'DIAZ2020', accession = 'GSE149428', doi = '10.7554/eLife.52707',
  source_url = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7546737/',
  cells = paste(names(counts),collapse=';'), rna_samples = nrow(samples),
  independent_function = 'nuclei counts, ATP and metabolic signal; published means/SD',
  measured_joint_rna = TRUE, measured_triple_rna = TRUE, phosphoprotein = FALSE,
  limitation = 'Two cell models; high drug doses with broad/off-target effects; no proof of ERBB/MET mechanism; destructive assays in separate wells; incomplete joint RNA at48h',
  license = 'Article CC-BY4.0; source supplementary data attributed to original authors')
fwrite(study, file.path(root,'data/diaz_studies.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(root,'data/diaz_session_info.txt'))
save_table(data.table(script='/work/src/program_pilot/59_diaz_combination_data.R',
  sha256=digest(file='/work/src/program_pilot/59_diaz_combination_data.R',algo='sha256')), 'script_provenance')
message('Prepared Diaz RNA samples: ', nrow(samples), '; independent functional source rows: ', nrow(func),
  '; complete RNA factorials: ', sum(rbindlist(contrast_rows)$molecular_complete))
