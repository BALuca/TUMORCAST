#!/usr/bin/env Rscript
# Assess prespecified focal-pair functional measurements, independent of model predictions.
# No fitting, filtering for favorable effects, or reinterpretation of wells as biological repeats.
suppressPackageStartupMessages(library(data.table))
root <- "/work/results/program_pilot/60_defensible_response_composition"
out <- file.path(root, "data")
x <- fread(file.path(out, "transfer_molecular_links.tsv.gz"))[pair == "pictilisib + trametinib"]
stopifnot(nrow(x) > 0L, all(x$complete_measured_constituents))
x[, `:=`(pictilisib_uM = fifelse(drug_A == "pictilisib", dose_A_uM, dose_B_uM),
           trametinib_uM = fifelse(drug_A == "trametinib", dose_A_uM, dose_B_uM),
           pictilisib_viability = fifelse(drug_A == "pictilisib", viability_A, viability_B),
           trametinib_viability = fifelse(drug_A == "trametinib", viability_A, viability_B))]
x[, dose_context := paste(cell, pictilisib_uM, trametinib_uM, sep = "|")]
# A repeated physical dose context can occur on several plates from one assay date.
# First average within date; then each date receives equal weight for condition summaries.
day <- x[, .(pictilisib_viability = mean(pictilisib_viability), trametinib_viability = mean(trametinib_viability),
              viability_AB = mean(viability_AB), bliss_expected_viability = mean(bliss_expected_viability),
              excess_inhibition = mean(excess_inhibition), physical_wells = .N, plates = uniqueN(plate_id),
              plate_ids = paste(unique(plate_id), collapse = ";"),
              molecular_supported = all(molecular_supported),
              RNA_A_matching_status = paste(unique(A_matching_status), collapse = ";"),
              RNA_B_matching_status = paste(unique(B_matching_status), collapse = ";")),
          by = .(dose_context, cell, pictilisib_uM, trametinib_uM, source_release, assay_day)]
condition <- day[, .(pictilisib_viability = mean(pictilisib_viability), trametinib_viability = mean(trametinib_viability),
                    viability_AB = mean(viability_AB), bliss_expected_viability = mean(bliss_expected_viability),
                    excess_inhibition = mean(excess_inhibition), physical_wells = sum(physical_wells),
                    plates = sum(plates), assay_dates = uniqueN(assay_day), releases = uniqueN(source_release),
                    source_releases = paste(unique(source_release), collapse = ";"),
                    molecular_supported = all(molecular_supported),
                    excess_min_across_assay_days = min(excess_inhibition),
                    excess_max_across_assay_days = max(excess_inhibition),
                    excess_day_range = diff(range(excess_inhibition)),
                    excess_day_sd = if (.N > 1L) sd(excess_inhibition) else NA_real_,
                    same_sign_across_assay_days = all(excess_inhibition > 0) | all(excess_inhibition <= 0),
                    all_days_extra_inhibition_at_least10pp = all(excess_inhibition >= 10)),
                by = .(dose_context, cell, pictilisib_uM, trametinib_uM)]
stopifnot(!anyDuplicated(condition$dose_context))
condition[, `:=`(single_agent_floor_below10pct = pmin(pictilisib_viability, trametinib_viability) < 10,
                  bliss_floor_below10pct = bliss_expected_viability < 10,
                  function_assay = "ATP_CellTiterGlo_percent_of_plate_DMSO", RNA_time_h = 24, function_time_h = 72,
                  repeat_unit = "assay_creation_date_proxy; source_does_not_identify_biological_replicates")]
fwrite(day, file.path(out, "transfer_focal_assay_day_measurements.tsv"), sep = "\t", na = "NA")
fwrite(condition, file.path(out, "transfer_focal_condition_reliability.tsv"), sep = "\t", na = "NA")
summary <- rbindlist(lapply(c("all_functional_doses", "RNA_supported_doses", "RNA_unsupported_doses"), function(s) {
  z <- if (s == "all_functional_doses") condition else if (s == "RNA_supported_doses") condition[molecular_supported == TRUE] else condition[molecular_supported == FALSE]
  z[, .(stratum = s, dose_contexts = .N, physical_wells = sum(physical_wells),
        repeated_across_assay_dates = sum(assay_dates > 1),
        repeated_across_screen_releases = sum(releases > 1),
        pictilisib_viability_min = min(pictilisib_viability), pictilisib_viability_max = max(pictilisib_viability),
        trametinib_viability_min = min(trametinib_viability), trametinib_viability_max = max(trametinib_viability),
        combination_viability_min = min(viability_AB), combination_viability_max = max(viability_AB),
        excess_inhibition_median = median(excess_inhibition),
        excess_inhibition_min = min(excess_inhibition), excess_inhibition_max = max(excess_inhibition),
        extra_inhibition_at_least10pp = sum(excess_inhibition >= 10),
        antagonism_at_least10pp = sum(excess_inhibition <= -10),
        single_agent_floor_below10pct = sum(single_agent_floor_below10pct),
        bliss_floor_below10pct = sum(bliss_floor_below10pct),
        repeated_same_sign = sum(same_sign_across_assay_days & assay_dates > 1),
        repeated_all_days_extra_at_least10pp = sum(all_days_extra_inhibition_at_least10pp & assay_dates > 1),
        repeated_median_day_difference_pp = if (any(assay_dates > 1)) median(excess_day_range[assay_dates > 1]) else NA_real_),
    by = cell]
}))
fwrite(summary, file.path(out, "transfer_focal_reliability.tsv"), sep = "\t", na = "NA")
# The matching proportions describe molecular support only. They do not assess
# whether interpolation reproduces the unmeasured molecular response at that dose.
support <- unique(x[, .(dose_context, cell, molecular_supported, A_matching_status, B_matching_status)])
support[, both_measured_dose := A_matching_status == "measured_dose_within_3.5pct_rounding" &
          B_matching_status == "measured_dose_within_3.5pct_rounding"]
support[, any_interpolated := grepl("interpolation", A_matching_status) | grepl("interpolation", B_matching_status)]
support_summary <- support[, .(dose_contexts = .N, supported = sum(molecular_supported),
                               both_measured = sum(both_measured_dose),
                               supported_with_interpolation = sum(molecular_supported & any_interpolated)), by = cell]
fwrite(support_summary, file.path(out, "transfer_focal_molecular_support.tsv"), sep = "\t")
controls <- fread(file.path(out, "transfer_plate_controls.tsv"))[plate_id %chin% x$plate_id]
fwrite(controls, file.path(out, "transfer_focal_plate_controls.tsv"), sep = "\t")
source_scripts <- c("/work/src/program_pilot/60_prepare_combination_transfer.R",
                    "/work/src/program_pilot/60_transfer_functional_reliability.R")
code_manifest <- data.table(path = source_scripts,
                            sha256 = vapply(source_scripts, digest::digest, character(1), file = TRUE, algo = "sha256"))
fwrite(code_manifest, file.path(out, "transfer_preparation_code_manifest.tsv"), sep = "\t")
writeLines(c(
  "Focal trametinib+pictilisib functional reliability audit; measured outcomes only, no model-based selection.",
  "Context means average physical wells within assay date, then weight assay dates equally. Dose contexts are not independent experiments.",
  "Extra inhibition = measured single-agent Bliss expected viability minus measured combination viability, in percentage points.",
  "A 10-percentage-point effect is a descriptive assay-scale threshold, not a formal significance test.",
  "These focal combinations occur only in the Original release. No original-versus-validation comparison is possible.",
  "Only BT20 has exact dose contexts repeated across assay dates. MCF7 repeats plates within one date per layout; HS578T and MDAMB231 also lack exact-context cross-date repeats. Reciprocal layouts use different dose grids.",
  "No focal dose context has an individual-drug viability or Bliss expected viability below 10%. Strong-single-drug floor saturation does not explain the focal failure.",
  "All RNA-supported focal contexts require at least one constituent's log-dose interpolation between measured RNA doses. None measures both constituent doses exactly.",
  "RNA is measured at 24h in LINCS; function is 72h ATP in Jaaks. Cells, drugs and doses are joined, but these are different studies, media and experimental days.",
  "Repeated-day ranges quantify observed assay reproducibility, not confidence intervals based on independent biological replicates.",
  "A molecular-supported dose subset omits the low-trametinib and/or low-pictilisib arms; all functional doses remain available and plotted."
), file.path(root, "TRANSFER_FUNCTIONAL_RELIABILITY.txt"))
print(summary)
print(support_summary)
