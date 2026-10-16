#!/usr/bin/env Rscript
# Same-study BT20 constituent RNA and independent imaged-cell combination outcomes.
# Whole-partner exclusions apply to combinations only; individual drug data stay available.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(limma)
  library(digest)
})
setDTthreads(4)
root <- "/work/results/program_pilot/61_individual_responses_combination_transfer"
out <- file.path(root, "data")
source_root <- "/work/results/program_pilot/59_mechanistic_evidence_datasets/data"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
partners <- c("lapatinib", "neratinib", "selumetinib", "trametinib", "saracatinib", "nvp-tae684")
anchors <- c("alpelisib", "dactolisib")
write_table <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
nominal_dose <- function(x) {
  grid <- c(0.04, 0.12, 0.37, 1.11, 3.33, 10)
  vapply(x, function(v) {
    if (!is.finite(v)) return(NA_real_)
    if (v == 0) return(0)
    i <- which.min(abs(log(grid / v)))
    if (abs(log(grid[i] / v)) <= 0.035) grid[i] else NA_real_
  }, numeric(1))
}
format_dose <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 10)

# Normalize measured RNA first so representation scoring can proceed while function is prepared.
# This reads the existing compact qs data, never the full GCTX or inferred genes.
meta <- fread(file.path(source_root, "matched_samples.tsv.gz"))[cell_line == "BT20"]
raw_rna <- qread(file.path(source_root, "matched_measured_RNA_Level2.qs"))
stopifnot(all(meta$sample_id %in% colnames(raw_rna)), nrow(raw_rna) == 978L,
          all(meta$time_h %in% c(3, 24)), all(meta$n_vehicle >= 5))
raw_rna <- raw_rna[, meta$sample_id, drop = FALSE]
stopifnot(all(is.finite(raw_rna)), all(raw_rna >= 0))
expression <- limma::normalizeQuantiles(log2(1 + raw_rna))
dimnames(expression) <- dimnames(raw_rna)
baseline_plate <- vapply(unique(meta$det_plate), function(p) {
  rowMeans(expression[, meta$det_plate == p & meta$pert_type == "ctl_vehicle", drop = FALSE])
}, numeric(nrow(expression)))
rownames(baseline_plate) <- rownames(expression)
displacement <- expression - baseline_plate[, match(meta$det_plate, colnames(baseline_plate)), drop = FALSE]
dimnames(displacement) <- dimnames(expression)
baseline <- vapply(c(3, 24), function(t) {
  rowMeans(baseline_plate[, unique(meta[time_h == t, det_plate]), drop = FALSE])
}, numeric(nrow(expression)))
dimnames(baseline) <- list(rownames(expression), c("BT20|3h", "BT20|24h"))
meta[, `:=`(cell = cell_line, profile_id = paste(cell_line, drug, dose_uM, paste0(time_h, "h"), sep = "|"))]
profiles <- meta[pert_type == "trt_cp", .(n_samples = .N, n_plates = uniqueN(det_plate),
                                          sample_ids = paste(sample_id, collapse = ";"),
                                          replicate_type = "L1000_assay_replicates; source_methods_call_technical"),
                 by = .(profile_id, cell, drug, dose_uM, time_h)]
profile_expression <- vapply(profiles$profile_id, function(p) {
  rowMeans(expression[, meta$profile_id == p, drop = FALSE])
}, numeric(nrow(expression)))
profile_displacement <- vapply(profiles$profile_id, function(p) {
  rowMeans(displacement[, meta$profile_id == p, drop = FALSE])
}, numeric(nrow(expression)))
profile_se <- vapply(profiles$profile_id, function(p) {
  z <- displacement[, meta$profile_id == p, drop = FALSE]
  if (ncol(z) > 1L) apply(z, 1L, sd) / sqrt(ncol(z)) else rep(NA_real_, nrow(z))
}, numeric(nrow(expression)))
rownames(profile_expression) <- rownames(profile_displacement) <- rownames(profile_se) <- rownames(expression)
stopifnot(identical(colnames(profile_displacement), profiles$profile_id), all(is.finite(profile_displacement)))
qsave(list(expression = profile_expression, displacement = profile_displacement,
           baseline = baseline, metadata = profiles, displacement_assay_replicate_se = profile_se,
           normalization = "log2(1+measured_L1000_fluorescence); limma_quantile_normalization; subtract_own_plate_vehicle_mean"),
      file.path(out, "matched_molecular.qs"), preset = "high")
write_table(meta, "matched_RNA_samples.tsv.gz")
write_table(profiles, "matched_RNA_profiles.tsv.gz")
message("RNA ready: ", nrow(profile_displacement), " genes x ", ncol(profile_displacement), " profiles")

# Preserve every physical functional well and source normalization. Controls are
# recorded independently; no outcome-based removal of a biological repeat.
raw <- fread(file.path(source_root, "matched_combination_counts.tsv.gz"))
stopifnot(all(raw$cell_line == "BT20"), all(raw$biological_replicate %in% 1:3),
          all(raw$technical_replicate %in% 1:3), all(raw$n0 > 0), all(raw$ncontrol72 > raw$n0),
          all(raw$n72 > 0), !anyDuplicated(raw$sample_id))
stopifnot(max(abs(raw$recomputed_GR - raw$GR)) < 1e-4)
raw[, `:=`(cell = cell_line, nominal_A_uM = nominal_dose(dose_a_uM), nominal_B_uM = nominal_dose(dose_b_uM))]
# Unsupported intermediate doses retain their exact values and are never merged
# merely because both have nominal NA. Round-equivalent grid doses share a key.
raw[, `:=`(dose_key_A = fifelse(is.finite(nominal_A_uM), paste0("nominal:", format_dose(nominal_A_uM)),
                               paste0("actual:", format_dose(dose_a_uM))),
            dose_key_B = fifelse(is.finite(nominal_B_uM), paste0("nominal:", format_dose(nominal_B_uM)),
                               paste0("actual:", format_dose(dose_b_uM))))]
raw[, matched_condition_id := paste(cell, drug_a, dose_key_A, drug_b, dose_key_B, sep = "|")]
write_table(raw, "matched_functional_wells.tsv.gz")
controls <- raw[arm == "vehicle", .(physical_wells = .N, measured_vehicle_mean = mean(n72),
                                     measured_vehicle_sd = sd(n72), provided_vehicle_count = unique(ncontrol72),
                                     provided_initial_count = unique(n0)),
                by = .(biological_replicate, technical_replicate)]
controls[, `:=`(vehicle_cv = measured_vehicle_sd / measured_vehicle_mean,
                control_growth_fold = provided_vehicle_count / provided_initial_count,
                normalization_note = "source_provided_vehicle_denominator_preserved; not_replaced_by_mean_of_listed_control_wells")]
write_table(controls, "matched_functional_controls.tsv")
avg <- raw[arm != "vehicle", .(ratio = mean(relative_cell_count),
                                 ratio_technical_sd = sd(relative_cell_count), GR = mean(recomputed_GR),
                                 GR_technical_sd = sd(recomputed_GR), fold_vs_day0 = mean(fold_vs_day0),
                                 n0 = mean(n0), n72 = mean(n72), ncontrol72 = mean(ncontrol72),
                                 technical_n = .N, source_sample_ids = paste(sample_id, collapse = ";"),
                                 source_actual_doses_A = paste(sort(unique(dose_a_uM)), collapse = ";"),
                                 source_actual_doses_B = paste(sort(unique(dose_b_uM)), collapse = ";"),
                                 dose_A_uM = mean(dose_a_uM), dose_B_uM = mean(dose_b_uM),
                                 nominal_A_uM = unique(nominal_A_uM), nominal_B_uM = unique(nominal_B_uM)),
             by = .(matched_condition_id, cell, arm, drug_a, drug_b, dose_key_A, dose_key_B, biological_replicate)]
stopifnot(!anyDuplicated(avg[, .(matched_condition_id, biological_replicate)]), all(avg$technical_n %in% c(2L, 3L)))
single <- avg[arm == "single"]
write_table(single, "matched_same_screen_individual_function.tsv.gz")
combo <- avg[arm == "combination"]
stopifnot(all(combo$drug_a %in% anchors), all(combo$drug_b %in% partners))
setnames(combo, c("matched_condition_id", "drug_a", "drug_b", "ratio", "ratio_technical_sd", "GR", "GR_technical_sd",
                 "fold_vs_day0", "n0", "n72", "ncontrol72", "technical_n", "source_sample_ids"),
         c("condition_id", "anchor", "partner", "ratio_AB", "ratio_technical_sd_AB", "GR_AB", "GR_technical_sd_AB",
           "fold_vs_day0_AB", "n0_AB", "n72_AB", "ncontrol72_AB", "technical_n_AB", "source_sample_ids_AB"))
for (a in c("A", "B")) {
  z <- single[, .(cell, drug = drug_a, dose_key = dose_key_A, biological_replicate,
                   ratio, ratio_technical_sd, GR, GR_technical_sd, fold_vs_day0, n0, n72, ncontrol72,
                   technical_n, source_sample_ids)]
  setnames(z, c("drug", "dose_key"), c(if (a == "A") "anchor" else "partner", paste0("dose_key_", a)))
  value_cols <- c("ratio", "ratio_technical_sd", "GR", "GR_technical_sd", "fold_vs_day0", "n0", "n72", "ncontrol72",
                  "technical_n", "source_sample_ids")
  setnames(z, value_cols, paste0(value_cols, "_", a))
  combo <- merge(combo, z, by = c("cell", if (a == "A") "anchor" else "partner", paste0("dose_key_", a), "biological_replicate"),
                 all.x = TRUE, sort = FALSE)
}
combo[, `:=`(observation_id = paste(condition_id, paste0("bio", biological_replicate), sep = "|"),
              pair = paste(anchor, partner, sep = " + "),
              ratio_expected = ratio_A * ratio_B,
              excess_inhibition = 100 * (ratio_A * ratio_B - ratio_AB),
              viability_A = 100 * ratio_A, viability_B = 100 * ratio_B, viability_AB = 100 * ratio_AB,
              net_loss_A = fold_vs_day0_A < 1, net_loss_B = fold_vs_day0_B < 1, net_loss_AB = fold_vs_day0_AB < 1,
              complete_measured_constituents = is.finite(ratio_A) & is.finite(ratio_B),
              function_time_h = 72, function_assay = "imaged_cell_counts", exclusion_unit = partner)]
stopifnot(all(combo$complete_measured_constituents), !anyDuplicated(combo$observation_id))
write_table(combo, "matched_combination_observations.tsv.gz")
message("Function ready: ", nrow(combo), " condition x biological-repeat observations; ", uniqueN(combo$condition_id), " union contexts")

# Join each constituent to a directly measured nominal RNA dose. No interpolation
# or extrapolation, including when the source provides an intermediate dose.
links <- rbindlist(list(copy(combo)[, time_h := 3], copy(combo)[, time_h := 24]))
for (a in c("A", "B")) {
  z <- profiles[, .(cell, drug, dose_uM, time_h, profile_id, n_samples, n_plates)]
  setnames(z, c("drug", "dose_uM", "profile_id", "n_samples", "n_plates"),
           c(if (a == "A") "anchor" else "partner", paste0("nominal_", a, "_uM"),
             paste0("profile_", a), paste0("RNA_samples_", a), paste0("RNA_plates_", a)))
  links <- merge(links, z, by = c("cell", if (a == "A") "anchor" else "partner", paste0("nominal_", a, "_uM"), "time_h"),
                 all.x = TRUE, sort = FALSE)
}
links[, `:=`(molecular_supported = !is.na(profile_A) & !is.na(profile_B),
              molecular_matching_rule = "same_cell_same_drug_nominal_dose_within_3.5pct_rounding; no_interpolation",
              primary_time = time_h == 3)]
write_table(links, "matched_molecular_links.tsv.gz")
coverage <- links[, .(biological_repeat_observations = .N, union_dose_contexts = uniqueN(condition_id),
                       physical_combination_wells = sum(technical_n_AB), biological_replicates = uniqueN(biological_replicate)),
                  by = .(anchor, partner, time_h, molecular_supported)]
write_table(coverage, "matched_coverage.tsv")
exclusions <- CJ(held_partner = partners, observed_partner = partners)
exclusions[, `:=`(combination_role = fifelse(held_partner == observed_partner, "held_out", "training"),
                   individual_RNA_available = TRUE, individual_function_available = TRUE,
                   rule = "exclude_all_combinations_containing_held_partner_in_both_anchors_all_doses_all_replicates")]
write_table(exclusions, "matched_partner_exclusion_manifest.tsv")

# Retain the wider independent single-agent panel for single-response calibration.
# These assays are not substituted into the matched-combination Bliss comparator.
individual <- fread(file.path(source_root, "matched_single_agent_growth.tsv.gz"))[cell_line == "BT20"]
individual[, `:=`(cell = cell_line, nominal_dose_uM = nominal_dose(dose_uM),
                   availability_under_partner_exclusion = "all_single_drug_data_available_in_every_fold")]
write_table(individual, "matched_individual_function_panel.tsv.gz")
individual_links <- rbindlist(list(copy(individual)[, RNA_time_h := 3], copy(individual)[, RNA_time_h := 24]))
individual_links <- merge(individual_links, profiles[, .(cell, drug, nominal_dose_uM = dose_uM,
                                                         RNA_time_h = time_h, profile_id)],
                          by = c("cell", "drug", "nominal_dose_uM", "RNA_time_h"), all.x = TRUE, sort = FALSE)
write_table(individual_links, "matched_individual_function_RNA_links.tsv.gz")

# Genuine biological-repeat uncertainty is computed only within shared nominal contexts.
condition <- combo[, .(anchor = unique(anchor), partner = unique(partner), biological_replicates = .N,
                        bio_ids = paste(sort(biological_replicate), collapse = ";"),
                        extra_mean = mean(excess_inhibition), extra_sd = sd(excess_inhibition),
                        extra_se = sd(excess_inhibition) / sqrt(.N),
                        extra_min = min(excess_inhibition), extra_max = max(excess_inhibition),
                        all_repeats_positive = all(excess_inhibition > 0), all_repeats_ge10pp = all(excess_inhibition >= 10),
                        all_repeats_ge20pp = all(excess_inhibition >= 20),
                        combination_viability_mean = mean(viability_AB), GR_AB_mean = mean(GR_AB),
                        net_loss_repeats = sum(net_loss_AB),
                        fold_vs_day0_AB_mean = mean(fold_vs_day0_AB),
                        single_floor_repeats = sum(pmin(viability_A, viability_B) < 10),
                        nominal_A_uM = unique(nominal_A_uM), nominal_B_uM = unique(nominal_B_uM)), by = condition_id]
condition[, molecular_supported := is.finite(nominal_A_uM) & is.finite(nominal_B_uM)]
write_table(condition, "matched_biological_repeatability.tsv")
repeat_pairs <- merge(combo[, .(condition_id, biological_replicate, extra = excess_inhibition)],
                      combo[, .(condition_id, biological_replicate, extra = excess_inhibition)],
                      by = "condition_id", allow.cartesian = TRUE)[biological_replicate.x < biological_replicate.y]
repeat_agreement <- repeat_pairs[, .(shared_dose_contexts = .N, Pearson_r = cor(extra.x, extra.y),
                                     median_absolute_difference_pp = median(abs(extra.x - extra.y)),
                                     sign_agreement_fraction = mean((extra.x > 0) == (extra.y > 0))),
                                 by = .(biological_replicate.x, biological_replicate.y)]
write_table(repeat_agreement, "matched_repeat_agreement.tsv")
supported_repeat_pairs <- repeat_pairs[condition_id %chin% condition[molecular_supported == TRUE, condition_id]]
supported_repeat_agreement <- supported_repeat_pairs[, .(shared_dose_contexts = .N, Pearson_r = cor(extra.x, extra.y),
                                                         median_absolute_difference_pp = median(abs(extra.x - extra.y)),
                                                         sign_agreement_fraction = mean((extra.x > 0) == (extra.y > 0))),
                                                     by = .(biological_replicate.x, biological_replicate.y)]
write_table(supported_repeat_agreement, "matched_supported_repeat_agreement.tsv")
input_paths <- file.path(source_root, c("matched_samples.tsv.gz", "matched_measured_RNA_Level2.qs", "matched_combination_counts.tsv.gz",
                                       "matched_single_agent_growth.tsv.gz"))
manifest <- data.table(path = input_paths, bytes = file.info(input_paths)$size,
                       sha256 = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256"))
write_table(manifest, "matched_input_manifest.tsv")
writeLines(c(
  "Same-study single-drug RNA to combination function: BT20, Niepel et al. 2017, DOI 10.1038/s41467-017-01383-w.",
  "Primary molecular time3h; 24h is a prespecified secondary time. Function72h imaged-cell counts, not ATP.",
  "All 12 anchor+partner pairs retained. Whole-partner combination exclusion across both anchors/all doses/all repeats; all individual drug RNA and function remain available.",
  "Measured 978 L1000 genes only. Log2(1+fluorescence), limma quantile normalization, own-plate vehicle gene means subtracted.",
  "These RNA displacement units are not RNAseq log2CP10K. Program projections measure observed-gene contributions rather than complete original weighted scores.",
  "Source-provided n0 and vehicle counts are preserved. Mean of listed vehicle wells does not exactly equal the supplied normalization denominator.",
  "Technical wells are averaged first within biological repeat. Constituent single-drug ratios are matched from the same combination screen by biological repeat and dose.",
  "Primary extra inhibition in percentage points =100*(mean singleA ratio * mean singleB ratio - mean combination ratio).",
  "GR and treated/initial cell-count ratios are retained separately to distinguish reduced proliferation from net cell loss.",
  "Exact nominal RNA matching permits only3.5% source dose rounding; intermediate doses remain unsupported. No RNA dose interpolation or extrapolation.",
  "Biological repeat3 has much lower initial/control cell counts; all repeats remain, with control differences and genuine shared-context repeatability reported.",
  "Physical wells and dose contexts are not independent biological replicates. RNA assay replication is technical per source methods.",
  "Previously available outcomes make this a retrospective excluded-partner test, not new blind experimental validation."
), file.path(root, "MATCHED_DATA_DESIGN.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "matched_preparation_session_info.txt"))
print(coverage)
print(repeat_agreement)
