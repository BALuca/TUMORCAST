#!/usr/bin/env Rscript
# Read-only audit of completed models. No fitting, tuning, or new HGF selection.
# Tables retain original endpoints and validation axes; unlike assays are not pooled.
.libPaths(c('/work/renv/library/R-4.3/x86_64-pc-linux-gnu', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(qs); library(digest)})
result_root <- '/work/results/program_pilot'
out <- file.path(result_root, '141_audit_v106_to_v140', 'data')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
roots <- list.dirs(result_root, recursive=FALSE)
roots <- roots[grepl('^1[0-4][0-9]_', basename(roots))]
versions <- as.integer(substr(basename(roots), 1, 3))
roots <- roots[versions %in% 106:140]
versions <- as.integer(substr(basename(roots), 1, 3))
stopifnot(identical(sort(versions), 106:140))
names(roots) <- versions
write_table <- function(x, name) fwrite(x, file.path(out, paste0(name, '.tsv')), sep='\t')
read_data <- function(v, name) fread(file.path(roots[as.character(v)], 'linear/data', paste0(name, '.tsv')))
same <- function(a,b) isTRUE(all.equal(a,b))
reference <- qread(file.path(roots['106'], 'linear/models/final.qs'))
drugs <- reference$selected$drug

# Inspect the original single-drug components, separately from changing interaction objects.
invariance <- rbindlist(lapply(106:140, function(v) {
  if (v %in% c(130,134)) return(data.table(version=v, status='diagnostic only'))
  parent_version <- if (v %in% c(139,140)) 136L else v
  model <- qread(file.path(roots[as.character(parent_version)], 'linear/models/final.qs'))
  stopifnot(all(drugs %in% rownames(model$context$mu)))
  data.table(version=v, status=if(v %in% c(139,140)) 'uses frozen v136 single-drug parent' else 'embedded single-drug model',
    same_level=same(reference$par$level,model$par$level),
    same_amplitude=same(reference$par$log_amplitude,model$par$log_amplitude),
    same_signatures=same(reference$context$W,model$context$W),
    same_lapatinib_RNA=same(reference$rna$lap_receptors,model$rna$lap_receptors),
    same_original_drug_means=same(reference$context$mu[drugs,],model$context$mu[drugs,]),
    same_original_drug_SD=same(reference$context$sigma[drugs,],model$context$sigma[drugs,]))
}), fill=TRUE)
write_table(invariance, 'single_drug_component_invariance')

# Compare exactly the same lapatinib recipients before comparing the larger deployment set.
m <- read_data(106,'measured_predictions_frozen')[drug=='Lapatinib ditosylate' & evaluation=='held_cell']
b <- read_data(106,'baseline_predictions_frozen')[drug=='Lapatinib ditosylate' & evaluation=='held_cell']
keys <- c('depmap_id','drug','dose')
stopifnot(!anyDuplicated(m[,..keys]), !anyDuplicated(b[,..keys]))
matched <- merge(m,b,by=keys,suffixes=c('_measured','_predicted'))
stopifnot(nrow(matched)==30, nrow(b)==441,
          max(abs(matched$observed_log2_measured-matched$observed_log2_predicted))<1e-9)
single <- data.table(
  molecular_input=c('measured displacement','baseline-predicted displacement','baseline-predicted displacement'),
  recipient_set=c('same 30 measured-RNA recipients','same 30 measured-RNA recipients','all 441 PRISM recipients'),
  cells=c(30,30,441),
  r_log2_abundance=c(cor(matched$observed_log2_measured,matched$predicted_log2_measured),
    cor(matched$observed_log2_predicted,matched$predicted_log2_predicted),cor(b$observed_log2,b$predicted_log2)),
  r_abundance=c(cor(2^matched$observed_log2_measured,2^matched$predicted_log2_measured),
    cor(2^matched$observed_log2_predicted,2^matched$predicted_log2_predicted),cor(2^b$observed_log2,2^b$predicted_log2)))
single[, comparison_axis := 'across held-out cells at fixed lapatinib dose 5 uM; retrospective protocol']
write_table(single,'lapatinib_matched_population_correlations')

# Archive every existing query and validation summary with its original schema and provenance.
query_tables <- c('au565_hgf_5uM_zero_reference_query','au565_hgf_5uM_predictions_frozen',
                  'AU565_HGF_predictions_frozen','added_HGF_displacement_query')
metric_tables <- c('combination_correlations','mema_correlations','HER2_external_held_cell_correlations',
  'external_held_cell_correlations','external_whole_compound_correlations',
  'across_ligand_held_correlations','lapatinib_held_calibration_correlations')
collected_queries <- list()
collected_metrics <- list()
inputs <- list()
for (v in 106:140) {
  for (name in c(query_tables,metric_tables)) {
    path <- file.path(roots[as.character(v)],'linear/data',paste0(name,'.tsv'))
    # Different historical versions intentionally produced different named summaries.
    if (!file.exists(path)) next
    tab <- fread(path)
    tab[, `:=`(audit_version=v, audit_source=path, audit_table=name)]
    tag <- paste(v,name,sep='_')
    if (name %in% query_tables) collected_queries[[tag]] <- tab else collected_metrics[[tag]] <- tab
    inputs[[tag]] <- data.table(path=path,sha256=digest(path,file=TRUE,algo='sha256'))
  }
}
write_table(rbindlist(collected_queries,fill=TRUE),'all_saved_HGF_queries_original_units')
write_table(rbindlist(collected_metrics,fill=TRUE),'all_saved_validation_summaries_original_axes')
write_table(rbindlist(inputs),'input_hashes')

# Decompose direct addition under the unchanged decoder; do not infer a causal calibration.
contrib <- read_data(139,'added_HGF_displacement_program_contributions')
query <- read_data(139,'added_HGF_displacement_query')
programs <- sprintf('P%02d',1:40)
rr <- reference$template$correlations[,'Lapatinib ditosylate']
cc <- as.numeric(contrib[dose==5 & source_profile=='Primary_study_mean',..programs])
stopifnot(length(cc)==40, identical(names(rr),programs))
per_program <- data.table(program=programs, observed_lapatinib_sensitivity_r=rr,
  HGF_added_log2_abundance_contribution=cc, strong_lapatinib_association=abs(rr)>=.5)
write_table(per_program,'direct_addition_program_contributions')
write_table(query,'direct_addition_predictions')
stopifnot(abs(sum(cc)-query[dose==5 & source_profile=='Primary_study_mean',HGF_increment_log2_abundance])<1e-10)

# Test the simplest fixed sensitivity-opposition score against the existing HER2 ligand outcomes.
# This is descriptive: full-data frozen signatures, no fit and no new held-out claim.
ligands <- read_data(136,'external_ligand_programs')
outcomes <- read_data(136,'external_combination_outcomes')[study=='elife32271']
ix <- match(paste(outcomes$depmap_id,outcomes$gene),paste(ligands$depmap_id,ligands$gene))
stopifnot(!anyNA(ix),nrow(outcomes)==16)
L <- as.matrix(ligands[ix,..programs])
route <- reference$context$route['Lapatinib ditosylate']
w <- reference$context$W[,route]
sigma <- reference$context$sigma['Lapatinib ditosylate',]
outcomes[, opposition_raw_signature := as.numeric(-L %*% w)]
outcomes[, opposition_decoder_units := as.numeric(-L %*% (w/sigma))]
write_table(outcomes,'HER2_fixed_opposition_scores_descriptive')
opposition <- outcomes[,.(cells=.N,
  r_raw_signature=if(.N>=4) cor(opposition_raw_signature,observed_departure_log2) else NA_real_,
  r_decoder_units=if(.N>=4) cor(opposition_decoder_units,observed_departure_log2) else NA_real_),by=.(gene,dose)]
write_table(opposition,'HER2_fixed_opposition_correlations_descriptive')

# Human-reviewed implementation history; diagnostics are explicitly distinguished from models.
changes <- c(
 'Fixed two sensitivity routes; absolute single-drug decoder; three ALMANAC interaction products',
 'Freeze single drugs; fit interaction to within-pair departures',
 'Add MEMA ligand contrasts and assay calibration',
 'Joint conditional MEMA means and ALMANAC departures',
 'Signature-by-program interactions: 80 coefficients',
 'Add RAF/MEK melanoma ligand combinations',
 'External combinations plus signature-by-program interactions',
 'Nonlinear molecular modifiers',
 'BRAF-conditioned external drug RNA; gene-level HGF averaging before transport',
 'Separate MEMA conditional contrasts from externally anchored departures',
 'Training-only interaction penalty selection',
 'Signed signature opposition and alignment channels',
 'Four time-matched cancer HGF study families',
 'Remove functional gain from interaction molecular inputs',
 'Source-robust views; fit ligand-alone conversion',
 'Add an HGF consensus molecular coordinate',
 'RNA-only unanchored bilinear interaction basis',
 'Add 16 lapatinib NRG1/EGF count combinations',
 'Seven RNA axes plus HGF coordinate; 36 bilinear terms',
 'Equal external-study weighting',
 'Prioritize within-treatment contrasts over means',
 'Select penalties using within-treatment contrasts',
 'Add centered sensitivity-score modulation channels',
 'Add direct HGF-coordinate interaction terms',
 'DIAGNOSTIC: original-platform HGF measurement reliability',
 'Correct noisy HGF platform preprocessing and precision weights',
 'Add HCC1954 lapatinib-HGF conditional outcome',
 'Separate count magnitude from net-growth ordering',
 'DIAGNOSTIC: cumulative interaction and HGF audit',
 'Return to v128 architecture with corrected HGF; omit extra coordinate',
 'Programwise reversal and reinforcement products',
 'Train those products on RAF/MEK only',
 'Pool two routes and constrain rescue/reinforcement signs',
 'Raw-ligand shared coefficient vector; drug multipliers',
 'Fixed-route signed program contributions; sparse nonnegative rescue weights')
stopifnot(length(changes)==35)
write_table(data.table(version=106:140,folder=unname(roots[as.character(106:140)]),change=changes),
            'version_change_manifest')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
print(single)
print(opposition)
print(per_program[,.(log2_increment=sum(HGF_added_log2_abundance_contribution)),by=strong_lapatinib_association])
cat('Audit complete. No models were fitted or changed.\n')
