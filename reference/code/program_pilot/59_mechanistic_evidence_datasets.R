#!/usr/bin/env Rscript
# Consolidate the source-specific 59 adapters into a scientific experiment index.
# This does not train a predictor, select favorable drugs, or merge study scales.
suppressPackageStartupMessages({library(data.table); library(qs); library(digest)})
root <- '/work/results/program_pilot/59_mechanistic_evidence_datasets'
data_dir <- file.path(root, 'data')
setDTthreads(2L)

manifests <- list()
for (prefix in c('diaz', 'mek', 'matched', 'hgf')) {
  x <- fread(file.path(data_dir, paste0(prefix, '_download_manifest.tsv')))
  if ('local_path' %in% names(x)) setnames(x, 'local_path', 'path')
  if ('source_url' %in% names(x)) setnames(x, 'source_url', 'url')
  if ('size_bytes' %in% names(x)) setnames(x, 'size_bytes', 'bytes')
  x[, adapter := prefix]
  stopifnot(all(file.exists(x$path)))
  manifests[[prefix]] <- x
}
dose_manifest <- fread(file.path(data_dir, 'diaz_dose_source_manifest.tsv'))
setnames(dose_manifest, c('file','size_bytes','primary_url'), c('path','bytes','url'))
dose_manifest[, adapter := 'diaz_dose']
stopifnot(all(file.exists(dose_manifest$path)))
manifests[['diaz_dose']] <- dose_manifest
manifest <- rbindlist(manifests, fill = TRUE)
fwrite(manifest, file.path(data_dir, 'source_manifest.tsv'), sep = '\t')

diaz_samples <- fread(file.path(data_dir, 'diaz_samples.tsv.gz'))
diaz_factorials <- fread(file.path(data_dir, 'diaz_functional_factorials.tsv.gz'))
diaz_contrasts <- fread(file.path(data_dir, 'diaz_contrasts.tsv.gz'))
diaz_coverage <- fread(file.path(data_dir, 'diaz_bank43_coverage.tsv.gz'))
diaz_scores <- fread(file.path(data_dir, 'diaz_bank43_scores.tsv.gz'))
dose_coverage <- fread(file.path(data_dir, 'diaz_dose_bank43_coverage.tsv'))
dose_scores <- fread(file.path(data_dir, 'diaz_dose_bank43_scores.tsv'))
mek <- fread(file.path(data_dir, 'mek_measurements.tsv.gz'))
matched <- fread(file.path(data_dir, 'matched_studies.tsv'))
matched_links <- fread(file.path(data_dir, 'matched_rna_growth_conditions.tsv.gz'))
stopifnot(!anyDuplicated(diaz_samples$sample_id),
  setequal(diaz_samples$sample_id, diaz_scores$sample_id),
  all(diaz_coverage$fraction_measured >= .8), all(is.finite(diaz_scores$displacement)),
  nrow(dose_scores) == 612L, all(is.finite(dose_scores$displacement)),
  nrow(dose_coverage) == 34L, all(dose_coverage$fraction_measured >= .8),
  all(file.exists(mek$source_file)))

# Full genetic-control dose curves, including null doses, are retained. Four
# source experiment IDs are matched between MYC and its LacZ control, not between
# unrelated assays. A published percentage is not a raw cell count.
myc <- mek[source_panel %chin% c('Fig8F', 'Fig8G', 'Fig8H')]
myc_keys <- c('cell_line','treatment','dose_A_uM','replicate','biological_state')
stopifnot(myc[, max(.N), by = myc_keys]$V1 |> max() == 1L)
myc_pairs <- dcast(myc, cell_line + treatment + dose_A_uM + replicate ~ biological_state,
  value.var = 'value')
stopifnot(all(is.finite(myc_pairs$LacZ_control)), all(is.finite(myc_pairs$MYC_overexpression)))
myc_pairs[, MYC_minus_LacZ_percentage_points := MYC_overexpression - LacZ_control]
myc_curves <- myc_pairs[, .(independent_experiments = .N,
  LacZ_mean_pct = mean(LacZ_control), MYC_mean_pct = mean(MYC_overexpression),
  effect_percentage_points = mean(MYC_minus_LacZ_percentage_points),
  effect_SD = sd(MYC_minus_LacZ_percentage_points),
  experiments_with_higher_MYC_viability = sum(MYC_minus_LacZ_percentage_points > 0)),
  by = .(cell_line, treatment, dose_uM = dose_A_uM)]
myc_curves[, `:=`(unit = 'percent genotype-matched vehicle; difference in percentage points',
  interpretation = 'Direct genetic intervention; whole curves retained, no selected-dose mechanism inference')]
fwrite(myc_pairs, file.path(data_dir, 'MYC_genetic_control_paired_measurements.tsv'), sep = '\t')
fwrite(myc_curves, file.path(data_dir, 'MYC_genetic_control_complete_curves.tsv'), sep = '\t')

# Show every pair/triple functional condition, with ATP and nuclei counts kept
# separate. The table is assay evidence, not a fitted composition prediction.
functional_check <- copy(diaz_factorials)
functional_check[, benefit_over_better_constituent_pct := pmin(A_pct, B_pct) - AB_pct]
functional_check[, measurement_level := 'Published means/SD; source-table repeated singles are not additional replicates']
fwrite(functional_check, file.path(data_dir, 'Diaz_complete_functional_controls.tsv'), sep = '\t')

inventory <- data.table(
  priority = c(1L, 2L, 3L, 4L, 5L),
  study = c('Diaz2020', 'Macaya2023', 'Niepel2017', 'vanBentum2025', 'Shaban/Kamashev2024'),
  purpose = c('Direct response-composition test with observed joint RNA and function',
    'Causal MYC resistance control and RNA/protein discrepancy',
    'Matched response-function test across cells and complete dose curves',
    'Independent measurement of MEK feedback and alternate-route blockade',
    'ERBB ligand rescue reference with RNA/protein/function'),
  prepared_measurements = c(
    sprintf('%s combination RNA samples plus18dose-control RNA samples; %s complete RNA factorials; independent nuclei/ATP means and SD', nrow(diaz_samples), sum(diaz_contrasts$molecular_complete)),
    '14,904 gene symbols x 9 RNA samples; TMT protein; MYC qPCR; 3 genetic-control drug/cell dose curves, 4 experiments each',
    sprintf('%s measured RNA samples x 978 genes; 109 drugs, 6 cells, 6 doses; %s matched RNA-time/function conditions', matched$RNA_samples, matched$matched_RNA_function_conditions),
    '3 cells: ligand x MEKi pAKT and ERK/JNK phosphosites; 2 cells: PI3Ki/JNKi x MEKi growth grids',
    '10,485 numeric RNA coefficients; 1,018 named protein coefficients; 3 printed growth summaries; 12 RNA run records'),
  crucial_controls = c('Vehicle, T/M/W, TM/TW/MW, LNCaP TMW; time course; additional MCF7 dose experiment',
    'MYC versus LacZ at every dose, each with own vehicle; drug singles and combinations',
    'Plate vehicles; initial and endpoint cell counts; growth biological replicates; both constituent singles',
    'Ligand alone, MEKi alone, both; alternate-route inhibitor alone and with MEKi',
    'Vehicle, lapatinib, EGF, lapatinib+EGF, serum, lapatinib+serum'),
  main_limit = c('Two cells, three high-dose drugs; no ERBB/MET mechanism; functional replicate wells unavailable; separate assay seeding densities; dose experiment later batch',
    'No RNA-seq in joint-treatment or MYC-genetic arms; RNA16h/protein48h/function72h and genetic function10days',
    '978 measured landmarks; cannot claim complete bank43 coverage or joint-combination RNA; RNA technical replication',
    'No matched RNA; different MEKi dose/media in molecular and growth experiments; do not join as exact condition pairs',
    'EGF not HGF; censored RNA; uneven RNA replication; full LFQ lacks verified arm mapping; no growth replicate table'),
  source_url = c('https://pmc.ncbi.nlm.nih.gov/articles/PMC7546737/',
    'https://www.nature.com/articles/s41467-023-41828-z',
    'https://www.nature.com/articles/s41467-017-01383-w',
    'https://www.nature.com/articles/s41467-025-59404-y',
    'https://pmc.ncbi.nlm.nih.gov/articles/PMC11472020/'))
fwrite(inventory, file.path(data_dir, 'study_inventory.tsv'), sep = '\t')

tests <- data.table(
  test = c('Observed composition', 'Compression failure', 'Dose or mechanism',
    'Causal regulator control', 'Measured bypass', 'Across-cell functional information'),
  resource = c('Diaz2020', 'Diaz2020', 'Diaz2020 MCF7 dose experiment',
    'Macaya2023', 'vanBentum2025', 'Niepel2017'),
  distinguish = c('RNA response predictable from singles versus new response appearing under joint exposure',
    'Functionally aligned gene information discarded by program averaging versus retained by programs',
    'Combination behaves like greater single-drug exposure versus a different molecular state at similar cell-count effect',
    'MYC is merely correlated with survival versus changing MYC changes survival',
    'Upstream ERK target engagement versus persistent alternative AKT/JNK signaling',
    'Expression mostly identifies drug/dose versus contains reproducible information about different cellular outcomes'),
  measured_truth = c('Independent nuclei counts for all singles, pairs and triple; ATP kept separate',
    'Same independent functional measurements, never labels reconstructed from RNA',
    'Actual dose-specific nuclei counts; use pheno_ss experiment matching RNA, not other dose curve',
    'Viability curves for forced MYC and LacZ in four experiments; construct and protein controls',
    'Direct pERK/pAKT/pJNK plus separate growth under PI3K/JNK blockade',
    '72h GR and initial/endpoint counts, with biological replicate IDs'),
  molecular_input = c('Observed RNA all arms; fixed bank43 hard means',
    'Identical RNA matrix as genes and fixed programs, matched genes and samples',
    'RNA dose curves and combination time course; each normalized to own batch vehicle',
    'MYC qPCR and independent protein; single-treatment RNA/proteome as auxiliary evidence',
    'Quantitative phosphoprotein source rows; no RNA proxy for phospho measurement',
    '3h/24h directly measured RNA, plate vehicles, independent baseline protein'),
  comparison = c('Joint minus single-A minus single-B plus vehicle, before later function; full third-order contrast in LNCaP',
    'Compare repeatable gene-level response contrasts to program-level contrasts across the same complete factorial',
    'Trace full dose trajectories; compare conditions with overlapping cell-count effects; report batch confounding',
    'MYC-minus-LacZ viability across every measured dose, paired by experiment',
    '(ligand+MEKi - MEKi) - (ligand - vehicle); then blocker incremental growth effect after its single effect',
    'Within-drug full dose trajectories across cells; retain quiet contexts and measure uncertainty from repeats'),
  evaluation_unit = c('Whole pair or triple, all its times together; LNCaP triple reserved from fitting',
    'Whole treatment arm; no selecting programs on the evaluated functional outcome',
    'Whole dose trajectory and experiment; never random dose-point split',
    'Independent experiment; full dose curve',
    'Biological replicate; average technical wells within replicate; no pseudo-pairing between protocols',
    'Entire cell or drug family; separate within-drug cell contrasts from across-drug effects'),
  failure_means = c('New joint RNA: missing response-composition term; accurate joint RNA but wrong function: downstream mapping/measurement gap',
    'Repeatable relevant genes but absent program contrast: aggregation loses measured information; absence at both levels alone does not prove no biology',
    'Dose trajectory reproduces pattern: no evidence for special combination program; batch prevents a causal between-experiment conclusion',
    'No reproducible genetic curve shift would weaken MYC mediation; MYC RNA-protein discrepancy alone cannot indict all transcription',
    'Phospho route without selective blocker functional increment is insufficient to name causal bypass; unmatched assay conditions constrain inference',
    'Reproducible function without matching molecular separation localizes insufficiency to measured RNA/time/features; weak function replicate agreement identifies endpoint uncertainty'),
  status = c('Matrices, controls, functional links and fixed bank scores prepared; composition prediction not fitted',
    'Same-sample gene matrices and all 34 program scores prepared',
    '18dose-control RNA samples and9matched function summaries prepared; separate batch retained',
    'Actual complete genetic curves extracted and summarized',
    'Actual molecular factorial and growth grids extracted; no invented exact condition join',
    '7,827 matched RNA-time/function conditions prepared; existing retrospective resource'))
fwrite(tests, file.path(data_dir, 'experimental_tests.tsv'), sep = '\t')

validation <- data.table(check = c('Diaz unique RNA sample IDs', 'Diaz RNA/program sample IDs agree',
  'Diaz all 34 programs >=80 percent genes in both cells', 'Diaz all program scores finite',
  'Source manifest paths exist', 'Genetic-control arms have complete experiment pairing',
  'Niepel matched RNA/function count'),
  observed = c(nrow(diaz_samples), uniqueN(diaz_scores$sample_id), min(diaz_coverage$fraction_measured),
    sum(is.finite(diaz_scores$displacement)), nrow(manifest), nrow(myc_pairs), sum(matched_links$valid_function_match)),
  status = 'PASS')
stopifnot(sum(matched_links$valid_function_match) == 7827L)
fwrite(validation, file.path(data_dir, 'validation.tsv'), sep = '\t')

writeLines(c(
  '59: MECHANISTIC EVIDENCE DATASETS', '',
  'Purpose: establish an independently measured biological contrast before asking whether our RNA programs capture it.',
  'The old logical error was using an RNA pattern both to name a mechanism and as evidence that the mechanism was captured.',
  'Here molecular states and function remain separate measurements, joined only where cell, compound and exposure permit.', '',
  'START WITH THE COMPLETE COMBINATION EXPERIMENT',
  'Diaz2020 gives vehicle, three singles, three pairs and (in LNCaP) the triple, with gene-level RNA and independently measured nuclei counts.',
  'The RNA matrix contains 282 samples across MCF7 and LNCaP. All three singles/pairs have RNA at 0,3,6,9,12,24h; LNCaP has all eight factorial arms.',
  'There are 42 complete RNA factorial contrasts. MCF7 has additional48h singles/control, but no48h pair RNA.',
  'Zero-hour RNA was collected within approximately30minutes after exposure: it is not pretreatment baseline.',
  'The fixed43 bank has34programs. Every program has >=80% gene coverage in both cells; its hard-mean displacements have been saved without outcome fitting.',
  'Counted nuclei is the primary functional readout. ATP/metabolic assays remain separate, allowing disagreement between cell abundance and metabolism to be seen.',
  'RNA and function are separate wells, seeded at8000 and10000cells/well respectively. Functional tables provide means/SD, with30missingSDs retained asNA.',
  'Repeated single-treatment source means are consolidated before joining to RNA; a rounded repeated SD is not an extra functional replicate.',
  'The full RNA triple factorial is observed, but raw vehicle wells and some joint SDs are absent. Normalized functional means support descriptive contrasts, not a replicated third-order significance test.',
  'High single-drug dose controls are a separate later experiment. They test the magnitude alternative, but cannot erase its batch confounding.',
  'The18dose-control RNA samples and9independent condition-level functional means are prepared. Fixed-bank scores use the same CP10K transform as the combination data, with each experiment ownDMSO controls.',
  'These are high-dose tamoxifen/mefloquine/withaferin-A experiments. They test response composition; they do not establish HGF/MET rescue.', '',
  'HOW THIS DISTINGUISHES FAILURE POINTS',
  '1. Compare measured joint RNA against its measured single responses. This directly reveals information that exists only after joint treatment.',
  '2. Ask whether repeatable gene-level contrasts survive the fixed program averages, using the identical samples.',
  '3. Relate earlier molecular states to later independent cell counts, preserving treatment arms as the units of evaluation.',
  'Inspect cell counts at the RNA collection time too. An earlier-than-final RNA sample is not automatically earlier than functional damage.',
  '4. Any composition model must be evaluated on a whole excluded pair/triple. A successful fit to its own joint observations is not transfer evidence.',
  '5. Similar functional effects at different doses are retained; quiet/noisy trajectories are identified from replicate uncertainty, not discarded because a model misses them.', '',
  'MECHANISM CONTROLS',
  'Macaya2023 supplies forced-MYC versus LacZ dose curves in four experiments, directly testing whether altering the candidate regulator changes drug sensitivity.',
  'Its source also supplies MYC RNA, protein and drug-combination viability. No RNA-seq exists in the genetic arms or combination, so these tables do not create a complete transcriptomic rescue experiment.',
  'The complete genetic curves are in data/MYC_genetic_control_complete_curves.tsv; every dose, including null doses, is retained.',
  'At12.5nM trametinib, H2009 mean viability/outgrowth is40.9% forLacZ and74.8% forMYC; all4experiments increase. This is a10daycrystal-violet assay, not direct cell counting.',
  'vanBentum2025 supplies ligand-by-MEK-inhibitor phosphoprotein factorials and growth with PI3K/JNK blockers. It measures the proposed alternate routes directly.',
  'The molecular and growth protocols differ in dose/media. They are mechanistic controls in the same cell models, not exact-condition RNA-to-function training pairs.', '',
  'BROADER MATCHED RNA/FUNCTION RESOURCE',
  'Niepel2017 supplies109drugs x6cells x6doses; 7827condition-by-RNA-time records match independent72hfunction.',
  'It contains978directly measured RNA genes, not a measured full transcriptome. RNA replicates are technical; functional biological repeats are separately labeled.',
  'BT20 combination wells and constituent RNA support partner comparisons, but combination RNA was not measured.',
  'Published GR, normalized net growth, nuclei percentages and ATP percentages are different measurements; the import keeps their names/units distinct.', '',
  'ERBB/HGF SCOPE',
  'The newly acquired lapatinib-plus-EGF/serum study contains RNA, protein and growth, but the RNA table censors nonsignificant entries and the full LFQ specimen key is unresolved.',
  'It is an auxiliary ERBB rescue reference, not the primary quantitative RNA benchmark. EGF is not HGF.',
  'Historical Wilson/Watson HGF experiments are indexed as references. Hand-coded +/-1 signaling is excluded from numeric mechanistic evidence.', '',
  'FILES',
  'data/study_inventory.tsv: which measurements each study actually supplies, source URLs and limitations.',
  'data/experimental_tests.tsv: the alternatives, independent truth, exact controls, exclusions and failure interpretation.',
  'data/source_manifest.tsv: all downloaded/reused originals, local paths and hashes.',
  'data/diaz_*: full RNA counts, normalized molecular matrices, program scores, condition joins, functional assays, factorial controls.',
  'data/diaz_dose_*: separately prepared dose-control experiment with its own metadata and function.',
  'data/mek_*: phosphoprotein, protein, qPCR, genetic/chemical intervention function and RNA sources.',
  'data/matched_*: same-study broad RNA/function and baseline protein sources.',
  'data/hgf_*: ERBB rescue sources, explicitly flagged censored/unmapped entries.',
  'data/validation.tsv: completed identity, coverage, pairing and source existence checks.', '',
  'figures/: complete genetic-control dose curves, independent cell-count/ATP time courses, and fixed-program RNA heatmaps.', '',
  'STATUS: acquisition/preparation and fixed-bank projection completed. No new molecular-to-function predictor has been fitted here.',
  'Program projection is a readable transformation of measured RNA, not a claim that the program bank passes these tests.',
  'All earlier scripts/results were preserved. The59adapters share one result folder and retain source-specific measurement units.'
), file.path(root, 'SCIENTIFIC_DESIGN.txt'))

files <- list.files(data_dir, full.names = TRUE)
files <- files[!basename(files) %in% c('file_registry.tsv')]
file_registry <- data.table(file = basename(files), path = files, bytes = file.info(files)$size)
file_registry[, format := sub('^.*\\.', '', file)]
fwrite(file_registry, file.path(data_dir, 'file_registry.tsv'), sep = '\t')
writeLines(capture.output(sessionInfo()), file.path(data_dir, 'consolidation_session_info.txt'))
message('Consolidated ', nrow(inventory), ' study groups, ', nrow(tests), ' scientific tests and ', nrow(manifest), ' source records.')
