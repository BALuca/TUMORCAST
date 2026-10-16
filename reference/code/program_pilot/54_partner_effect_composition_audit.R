#!/usr/bin/env Rscript
# Audit saved molecular composition predictions through matched partner changes.
# All fitting is historical. This script constructs contrasts and plots them in R.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(digest)
})
setDTthreads(2)

# Fixed source cohort, nomination thresholds and reporting choices.
root <- "/work"
base <- file.path(root, "almanac_broad9_2026-10-07")
out <- file.path(root, "results/program_pilot/54_partner_effect_composition_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (folder in c("data", "figures")) dir.create(file.path(out, folder), showWarnings = FALSE)
strong_cutoff <- 0.5
weak_cutoff <- 0.25
primary_cell <- "ACH-000552"
source(file.path(base, "model.R")) # Existing weighted Pearson calculation; no fitting invoked.
files <- c(
  battery = file.path(base, "redesigned_battery/results/held_out_predictions.tsv.gz"),
  splits = file.path(base, "redesigned_battery/results/splits.rds"),
  split_coverage = file.path(base, "redesigned_battery/results/split_coverage.tsv"),
  prepared = file.path(base, "results/prepared.rds"),
  selection = file.path(base, "cross_context_enrichment/results/selection_manifest.tsv.gz"),
  discovery_cells = file.path(base, "cross_context_enrichment/results/discovery_recipients.tsv"),
  quality = file.path(base, "redesigned_battery/results/functional_measurement_quality.tsv"),
  flexible = file.path(base, "failure_isolation/results/predictions.tsv.gz"),
  individual = file.path(base, "failure_isolation/results/additive_control_predictions.tsv.gz"),
  norm = file.path(base, "failure_isolation/results/flexible_norm_predictions.tsv.gz")
)
stopifnot(all(file.exists(files)))
input_hashes <- vapply(files, function(f) digest(file = f, algo = "sha256"), character(1))
write_table <- function(x, name) fwrite(x, file.path(out, "data", name), sep = "\t", na = "NA")
write_table(data.table(input = names(files), path = files, sha256 = input_hashes), "input_manifest.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
git_revision <- system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE)
writeLines(git_revision, file.path(out, "git_revision.txt"))
message("Reading saved predictions and freezing outcome-independent partner selection.")

# Define eligible profiles using only metadata and the already-frozen discovery
# nomination. All weak candidates are retained, not only the old same-pair controls.
selection <- fread(files[["selection"]])
stopifnot(!anyDuplicated(selection$source_row),
          !any(c("y", "yA", "yB", "interaction") %in% names(selection)))
discovery_cells <- fread(files[["discovery_cells"]])
stopifnot(!any(discovery_cells$depmap_id %in% selection$depmap_id))
selection[, support := fcase(stratum == "enriched", "Strong in discovery",
                             stratum %in% c("weak_candidate", "weak_control"), "Weak in discovery",
                             default = "Not nominated")]

# Both orientations use their own wholly held anchor drug. A profile holds its
# cell, laboratory and anchor concentration fixed; partner doses remain distinct.
oriented <- rbindlist(lapply(c("A", "B"), function(side) {
  z <- copy(selection)
  other <- if (side == "A") "B" else "A"
  z[, `:=`(anchor = get(paste0("id", side)), anchor_name = get(paste0("drug", side)),
           anchor_dose = get(paste0("dose", side)), partner = get(paste0("id", other)),
           partner_name = get(paste0("drug", other)), partner_dose = get(paste0("dose", other)),
           side = side)]
  z
}))
oriented[, profile := paste(anchor, sprintf("%.8g", anchor_dose), depmap_id, SCREENER, sep = "|")]
stopifnot(!anyDuplicated(oriented[, .(source_row, anchor)]))
metadata <- oriented[, .(source_row, anchor, anchor_name, anchor_dose, partner, partner_name,
                        partner_dose, profile, depmap_id, CELLNAME, SCREENER, side, support,
                        discovery_median, discovery_n, discovery_scope, selected, dose_stratum)]

# Reusable contrast construction for the two distinct saved validation regimes.
# Selection depends only on independent discovery and measurement coverage.
nominate_contrasts <- function(meta) {
  enriched <- meta[support == "Strong in discovery"]
  weak <- meta[support == "Weak in discovery"]
  a <- enriched[, .(profile, anchor, source_enriched = source_row, partner_enriched = partner,
                     dose_enriched = partner_dose, discovery_enriched = discovery_median)]
  b <- weak[, .(profile, source_weak = source_row, partner_weak = partner,
               dose_weak = partner_dose, discovery_weak = discovery_median)]
  pairs <- merge(a, b, by = "profile", allow.cartesian = TRUE)
  pairs <- pairs[partner_enriched != partner_weak &
                   abs(discovery_enriched - discovery_weak) >= strong_cutoff]
  pairs[, enriched_first := discovery_enriched < discovery_weak]
  pairs[, `:=`(source_more = fifelse(enriched_first, source_enriched, source_weak),
                source_less = fifelse(enriched_first, source_weak, source_enriched),
                partner_more = fifelse(enriched_first, partner_enriched, partner_weak),
                partner_less = fifelse(enriched_first, partner_weak, partner_enriched),
                discovery_delta = -abs(discovery_enriched - discovery_weak))]
  pairs[, partner_pair := paste(pmin(partner_more, partner_less), pmax(partner_more, partner_less), sep = "::")]
  pairs[, contrast_id := paste(profile, source_more, source_less, sep = "|")]
  stopifnot(!anyDuplicated(pairs$contrast_id), all(pairs$partner_more != pairs$partner_less),
            all(pairs$discovery_delta <= -strong_cutoff))
  pairs
}
primary_selection <- nominate_contrasts(metadata)
stopifnot(nrow(primary_selection) > 0L)
coverage <- metadata[, .(rows = .N, partners = uniqueN(partner),
                         enriched = sum(support == "Strong in discovery"),
                         weak = sum(support == "Weak in discovery"),
                         nominated_partners = uniqueN(partner[support != "Not nominated"])),
                     by = .(profile, anchor, anchor_name, anchor_dose, depmap_id, CELLNAME, SCREENER)]
coverage[, eligible := profile %in% primary_selection$profile]
coverage[, priority_cell := depmap_id == primary_cell]
setorder(coverage, -priority_cell, -nominated_partners, -partners, anchor, anchor_dose, depmap_id, SCREENER)
# Six distinct anchors for a compact overview, chosen without test outcomes.
display <- coverage[eligible == TRUE][!duplicated(anchor)][seq_len(min(6L, uniqueN(coverage[eligible == TRUE, anchor])))]
write_table(coverage, "profile_coverage.tsv")
write_table(primary_selection, "primary_contrast_selection.tsv.gz")
write_table(display, "coverage_selected_display_profiles.tsv")
selection_files <- file.path(out, "data", c("primary_contrast_selection.tsv.gz", "coverage_selected_display_profiles.tsv"))
write_table(data.table(path = selection_files,
                      sha256 = vapply(selection_files, function(f) digest(file = f, algo = "sha256"), character(1))),
            "display_selection_freeze.tsv")

# Verify historical drug exclusions from saved split indices. No model is loaded
# to refit; prepared rows are used only to reconstruct split membership.
prepared <- readRDS(files[["prepared"]])
source_rows <- prepared$d
splits <- readRDS(files[["splits"]])
drug_splits <- splits[vapply(splits, function(x) identical(x$test, "drug"), logical(1))]
split_audit <- rbindlist(lapply(drug_splits, function(s) {
  train <- source_rows[s$tr]
  test <- source_rows[s$te]
  data.table(anchor = s$fold, train_rows = nrow(train), test_rows = nrow(test),
             anchor_in_training = sum(train$idA == s$fold | train$idB == s$fold),
             all_test_rows_contain_anchor = all(test$idA == s$fold | test$idB == s$fold))
}))
stopifnot(all(split_audit$anchor_in_training == 0L), all(split_audit$all_test_rows_contain_anchor))
reference_split_audit <- fread(files[["split_coverage"]])[test == "drug"]
stopifnot(all(reference_split_audit$shared_reference_groups == 0L))
write_table(split_audit, "saved_drug_exclusion_audit.tsv")
rm(prepared)

# Join evaluated outcomes only after freezing eligibility and displayed profiles.
battery <- fread(files[["battery"]])[test == "drug"]
stopifnot(nrow(battery) == 2L * nrow(selection),
          all(battery$fold == battery$idA | battery$fold == battery$idB),
          !anyDuplicated(battery[, .(source_row, fold)]),
          max(abs(battery$interaction - (battery$y - battery$yA - battery$yB))) < 1e-10)
for (s in drug_splits) {
  stopifnot(setequal(battery[fold == s$fold, source_row], source_rows$source_row[s$te]))
}
setnames(battery, "fold", "anchor")
random_models <- sprintf("random_%02d", 1:20)
primary_models <- c("molecular", "magnitude", "identity", "zero", random_models)
keep <- c("source_row", "anchor", "interaction", "y", "yA", "yB", primary_models)
primary <- merge(metadata, battery[, ..keep], by = c("source_row", "anchor"), all.x = TRUE)
stopifnot(nrow(primary) == nrow(metadata), !anyNA(primary[, ..primary_models]))
primary[, `:=`(anchor_single = fifelse(side == "A", yA, yB),
                partner_single = fifelse(side == "A", yB, yA), regime = "Whole anchor drug excluded")]
primary[, `:=`(observed_abundance_pct = 100 * 2^y,
                independent_abundance_pct = 100 * 2^(yA + yB),
                observed_pct_of_independence = 100 * 2^interaction)]
quality <- fread(files[["quality"]])
primary <- merge(primary, quality, by = "source_row", all.x = TRUE)
stopifnot(!anyNA(primary$distinct_studies))
write_table(primary, "whole_drug_partner_predictions.tsv.gz")

# Matched secondary comparison: existing flexible and common-gate individual
# readouts, with recipient and unordered pair excluded (not whole drugs).
flexible <- fread(files[["flexible"]])[task == "recipient_pair" & model %in% c("raw", "flexible")]
individual <- fread(files[["individual"]])[task == "recipient_pair"]
norm <- fread(files[["norm"]])[task == "recipient_pair"]
stopifnot(uniqueN(norm$model) == 1L)
norm[, model := "flexible_norm"]
secondary_long <- rbindlist(list(flexible[, .(source_row, model, prediction)],
                                 individual[, .(source_row, model, prediction)],
                                 norm[, .(source_row, model, prediction)]))
stopifnot(!anyDuplicated(secondary_long[, .(source_row, model)]))
secondary_wide <- dcast(secondary_long, source_row ~ model, value.var = "prediction")
secondary_models <- c("raw", "flexible", "additive_control", "flexible_norm")
stopifnot(!anyNA(secondary_wide[, ..secondary_models]))
secondary_meta <- metadata[source_row %in% secondary_wide$source_row]
secondary_selection <- nominate_contrasts(secondary_meta)
secondary <- merge(secondary_meta, secondary_wide, by = "source_row")
secondary_outcome <- unique(flexible[, .(source_row, interaction, y, yA, yB, pair_fold)])
stopifnot(!anyDuplicated(secondary_outcome$source_row))
secondary <- merge(secondary, secondary_outcome, by = "source_row")
secondary[, `:=`(anchor_single = fifelse(side == "A", yA, yB),
                  partner_single = fifelse(side == "A", yB, yA),
                  regime = "Cell and pair excluded; constituent drugs may be seen")]
write_table(secondary, "cell_pair_partner_predictions.tsv.gz")
write_table(secondary_selection, "secondary_contrast_selection.tsv.gz")

# Calculate partner differences without refitting. Different source rows can
# carry different measured singleton references; preserve those differences.
evaluate_contrasts <- function(predictions, selection, models) {
  keys <- c("profile", "source_row")
  first <- match(paste(selection$profile, selection$source_more), paste(predictions$profile, predictions$source_row))
  second <- match(paste(selection$profile, selection$source_less), paste(predictions$profile, predictions$source_row))
  stopifnot(!anyNA(first), !anyNA(second))
  a <- predictions[first]
  b <- predictions[second]
  stopifnot(all(a$anchor == b$anchor), all(a$anchor_dose == b$anchor_dose),
            all(a$depmap_id == b$depmap_id), all(a$SCREENER == b$SCREENER), all(a$partner != b$partner))
  z <- copy(selection)
  z[, `:=`(regime = a$regime, depmap_id = a$depmap_id, CELLNAME = a$CELLNAME,
            anchor_name = a$anchor_name, anchor_dose = a$anchor_dose, SCREENER = a$SCREENER,
            partner_more_name = a$partner_name, partner_more_dose = a$partner_dose,
            partner_less_name = b$partner_name, partner_less_dose = b$partner_dose,
            observed_more = a$interaction, observed_less = b$interaction,
            observed_delta = a$interaction - b$interaction,
            anchor_reference_difference = a$anchor_single - b$anchor_single,
            partner_single_difference = a$partner_single - b$partner_single,
            raw_combination_difference = a$y - b$y,
            exact_both_profiles = a$dose_stratum == "both exact" & b$dose_stratum == "both exact")]
  if ("pair_fold" %in% names(predictions)) z[, same_pair_fit := a$pair_fold == b$pair_fold]
  else z[, same_pair_fit := TRUE]
  for (model in models) z[, (model) := a[[model]] - b[[model]]]
  # Algebra independently removes the singleton differences from the raw change.
  stopifnot(max(abs(z$observed_delta - (z$raw_combination_difference -
               z$anchor_reference_difference - z$partner_single_difference))) < 1e-10)
  z[, dose_contrasts := .N, by = .(profile, partner_pair)]
  z[, partner_pairs := uniqueN(partner_pair), by = profile]
  z[, profiles_per_anchor := uniqueN(profile), by = anchor]
  z[, weight := 1 / (dose_contrasts * partner_pairs * profiles_per_anchor)]
  z[, weight := weight / sum(weight)]
  z
}
primary_contrasts <- evaluate_contrasts(primary, primary_selection, primary_models)
secondary_contrasts <- evaluate_contrasts(secondary, secondary_selection, secondary_models)
write_table(primary_contrasts, "whole_drug_partner_contrasts.tsv.gz")
write_table(secondary_contrasts, "cell_pair_partner_contrasts.tsv.gz")
message("Partner contrasts saved: ", nrow(primary_contrasts), " whole-drug; ",
        nrow(secondary_contrasts), " cell/pair.")

# Summaries distinguish actual signal, predicted direction and effect magnitude.
# Counts are dependent contrasts, not independent experimental replications.
summarize_contrasts <- function(z, models) {
  rbindlist(lapply(models, function(model) {
    observed <- z$observed_delta
    predicted <- z[[model]]
    substantial <- observed <= -strong_cutoff
    opposite <- observed >= strong_cutoff
    w <- z$weight / sum(z$weight)
    data.table(model = model, contrasts = nrow(z), profiles = uniqueN(z$profile),
               anchors = uniqueN(z$anchor), recipients = uniqueN(z$depmap_id),
               weighted_pearson = wcorr(observed, predicted, w),
               discovery_direction_agreement = sum(w * (observed < 0)),
               substantial_observed = sum(substantial),
               substantial_direction_correct = sum(substantial & predicted < 0),
               substantial_magnitude_recovered = sum(substantial & predicted <= -strong_cutoff),
               strong_opposite_observed = sum(opposite),
               false_substantial_prediction = sum(predicted <= -strong_cutoff & observed > -weak_cutoff),
               median_observed_strong_delta = median(observed[substantial]),
               median_predicted_on_strong_delta = median(predicted[substantial]),
               median_predicted_on_strong_ratio_pct = median(100 * 2^predicted[substantial]))
  }))
}
primary_summary <- summarize_contrasts(primary_contrasts, primary_models)
secondary_summary <- summarize_contrasts(secondary_contrasts, secondary_models)
write_table(primary_summary, "whole_drug_contrast_summary.tsv")
write_table(secondary_summary, "cell_pair_contrast_summary.tsv")

# A correct contrast can arise by assigning an effect to the wrong partner.
# Require the independently strong endpoint AND the weak endpoint to be captured
# separately, with no test-outcome filtering of the original evaluated cohort.
endpoint_recovery <- function(predictions, contrasts, models) {
  i <- match(paste(contrasts$profile, contrasts$source_enriched),
             paste(predictions$profile, predictions$source_row))
  j <- match(paste(contrasts$profile, contrasts$source_weak),
             paste(predictions$profile, predictions$source_row))
  strong_observed <- predictions$interaction[i]
  weak_observed <- predictions$interaction[j]
  confirmed <- abs(strong_observed) >= strong_cutoff &
    sign(strong_observed) == sign(contrasts$discovery_enriched) & abs(weak_observed) <= weak_cutoff
  rbindlist(lapply(models, function(model) {
    strong_predicted <- predictions[[model]][i]
    weak_predicted <- predictions[[model]][j]
    recovered <- strong_predicted * sign(strong_observed) >= strong_cutoff &
      abs(weak_predicted) <= weak_cutoff
    data.table(contrast_id = contrasts$contrast_id, profile = contrasts$profile, anchor = contrasts$anchor,
               depmap_id = contrasts$depmap_id, model = model, confirmed_strong_and_weak = confirmed,
               both_effects_recovered = confirmed & recovered,
               observed_strong = strong_observed, observed_weak = weak_observed,
               predicted_strong = strong_predicted, predicted_weak = weak_predicted)
  }))
}
primary_endpoint_calls <- endpoint_recovery(primary, primary_contrasts, c("molecular", "magnitude", "identity"))
secondary_endpoint_calls <- endpoint_recovery(secondary, secondary_contrasts, secondary_models)
write_table(primary_endpoint_calls, "whole_drug_two_endpoint_recovery.tsv.gz")
write_table(secondary_endpoint_calls, "cell_pair_two_endpoint_recovery.tsv.gz")
endpoint_summary <- rbindlist(list(
  primary_endpoint_calls[, .(regime = "Whole drug", confirmed = sum(confirmed_strong_and_weak),
                            recovered = sum(both_effects_recovered),
                            recovered_profiles = uniqueN(profile[both_effects_recovered]),
                            recovered_anchors = uniqueN(anchor[both_effects_recovered])), by = model],
  secondary_endpoint_calls[, .(regime = "Cell and pair", confirmed = sum(confirmed_strong_and_weak),
                              recovered = sum(both_effects_recovered),
                              recovered_profiles = uniqueN(profile[both_effects_recovered]),
                              recovered_anchors = uniqueN(anchor[both_effects_recovered])), by = model]))
write_table(endpoint_summary, "two_endpoint_recovery_summary.tsv")
by_anchor <- primary_contrasts[, {
  z <- copy(.SD)
  z[, anchor := .BY$anchor]
  summarize_contrasts(z, c("molecular", "magnitude", "identity"))
}, by = anchor]
write_table(by_anchor, "whole_drug_summary_by_anchor.tsv")
by_profile <- primary_contrasts[, {
  z <- copy(.SD)
  z[, profile := .BY$profile]
  summarize_contrasts(z, c("molecular", "magnitude", "identity"))
}, by = profile]
write_table(by_profile, "whole_drug_summary_by_profile.tsv")

# Profile-centered endpoint correlations quantify partner ordering while removing
# each fixed anchor/context's average interaction. Both models use identical rows.
profile_signal <- function(z, models) {
  eligible <- uniqueN(z$profile)
  rbindlist(lapply(models, function(model) {
    q <- z[, .(anchor, profile, partner, observed = interaction, predicted = get(model))]
    q[, `:=`(observed_centered = observed - mean(observed),
              predicted_centered = predicted - mean(predicted)), by = profile]
    q[, exposures_per_partner := .N, by = .(profile, partner)]
    q[, partners_per_profile := uniqueN(partner), by = profile]
    q[, profiles_per_anchor := uniqueN(profile), by = anchor]
    q[, w := 1 / (exposures_per_partner * partners_per_profile * profiles_per_anchor)]
    q[, w := w / sum(w)]
    data.table(model = model, rows = nrow(q), profiles = eligible,
               raw_pearson = wcorr(q$observed, q$predicted, q$w),
               profile_centered_pearson = wcorr(q$observed_centered, q$predicted_centered, q$w))
  }))
}
write_table(profile_signal(primary[profile %in% primary_contrasts$profile],
                           c("molecular", "magnitude", "identity")), "whole_drug_profile_signal.tsv")
write_table(profile_signal(secondary[profile %in% secondary_contrasts$profile], secondary_models),
            "cell_pair_profile_signal.tsv")

# Reference sensitivity is a diagnostic restriction, never a replacement cohort.
# The 0.25-log2 cutoff is the existing weak-effect scale, not a tuned threshold.
reference_sensitivity <- primary_contrasts[abs(anchor_reference_difference) <= weak_cutoff]
write_table(summarize_contrasts(reference_sensitivity, c("molecular", "magnitude", "identity")),
            "similar_anchor_reference_sensitivity.tsv")
write_table(primary_contrasts[, .(contrasts = .N, profiles = uniqueN(profile),
                                  median_anchor_reference_gap = median(abs(anchor_reference_difference)),
                                  gap_over_quarter_log2 = sum(abs(anchor_reference_difference) > weak_cutoff),
                                  same_saved_fit = all(same_pair_fit))], "contrast_reference_audit.tsv")

# Flush a numerical report before the complete plot collection is generated.
report <- c(
  "54: PARTNER-EFFECT COMPOSITION AUDIT — saved predictions; no model refitted",
  "Nine-program historical ALMANAC bank. All results remain retrospective.",
  "Question: does changing the partner change outcome beyond its individual effect,",
  "and do the saved molecular models capture that change in a fixed cell/anchor/dose/lab?",
  "",
  sprintf("Whole-drug test: %d independently nominated contrasts in %d profiles, %d anchor drugs and %d recipients.",
          nrow(primary_contrasts), uniqueN(primary_contrasts$profile), uniqueN(primary_contrasts$anchor),
          uniqueN(primary_contrasts$depmap_id)),
  "Each contrast uses two predictions from the SAME model excluding its entire anchor drug.",
  "Enriched versus weak nomination uses forty other cells, not evaluated-cell outcomes.",
  "All existing weak candidates enter; the prior analysis chose only38 same-pair weak controls.",
  "This is recovery of independently nominated known interactions, not phenotype-unseen discovery.",
  "I=log2(vAB)-log2(vA)-log2(vB). Delta=I(lower discovery effect)-I(higher discovery effect).",
  "Delta=-1 means half as much abundance relative to the corresponding independent effects.",
  "Strong contrast: delta<=-0.5. Predicted magnitude recovery uses the same threshold.",
  "Counts are dependent contrasts; doses/cells/partners and single-agent references are reused.",
  "",
  "WHOLE-DRUG RESULTS:", capture.output(print(primary_summary[model %in% c("molecular", "magnitude", "identity", "zero")])),
  "", "CELL-AND-PAIR RESULTS:", capture.output(print(secondary_summary)),
  "", "STRONG AND WEAK ENDPOINTS MUST BOTH BE RIGHT:", capture.output(print(endpoint_summary)),
  "The nonlinear joint/common-gate comparison has a DIFFERENT holdout and smaller cohort.",
  "It excludes the cell and each pair; constituent drugs can occur in fitting.",
  "Different pair folds may supply the two predictions; same_pair_fit is exported.",
  "No saved common-gate individual-signature control exists in the broad whole-drug test.",
  "Do not attribute differences between these two experiments solely to decoder architecture.",
  "",
  "FIGURES: composition_overview.pdf; coverage_selected_partner_profiles.pdf;",
  "all_whole_drug_partner_profiles.pdf; cell_pair_control_profiles.pdf.",
  "retrospective_partner_examples.pdf shows recurring successes and wrong-allocation examples,",
  "explicitly selected after inspecting outcomes; these are explanations, not fresh validation.",
  "Rows in profile figures retain every saved partner exposure, including unselected rows.",
  "Main examples are coverage-selected before evaluation outcomes are joined.",
  "All drugs, doses, cells, source IDs and numerical contrasts are in data/.",
  "",
  "Limits: measured molecular single responses, not baseline-predicted responses.",
  "Interaction outcomes use measured functional singles; this isolates composition.",
  "Most functional contexts have one source study; molecular dose interpolation remains.",
  "HT29 alone was excluded from original unsupervised program discovery; other cells are",
  "functional holdouts where indicated. No causal program-edge or rescue claim follows automatically."
)
writeLines(report, file.path(out, "REPORT.txt"))
message("Numerical report saved. Rendering compact overview and partner profiles.")

# Named shared plot styles, reused across each complete partner profile.
labels <- c(molecular = "Program products", magnitude = "Response magnitude", identity = "Drug identity",
            flexible = "Joint program patterns", additive_control = "Individual patterns + gate",
            flexible_norm = "Magnitude + gate", raw = "Program products", observed = "Observed",
            discovery = "Other-cell discovery")
colors <- c(molecular = "#0072B2", magnitude = "#999999", identity = "#E69F00", flexible = "#0072B2",
            additive_control = "#D55E00", flexible_norm = "#999999", raw = "#009E73",
            observed = "#202020", discovery = "#CC79A7")
plot_theme <- theme_classic(base_size = 10) + theme(
  plot.title = element_text(size = 12, face = "bold"), plot.subtitle = element_text(size = 9),
  legend.position = "bottom", legend.title = element_blank(),
  strip.background = element_blank(), strip.text = element_text(size = 10),
  plot.caption = element_text(size = 8, hjust = 0), plot.margin = margin(8, 12, 8, 8))

contrast_scatter <- function(z, models, title, subtitle) {
  d <- melt(z, id.vars = c("contrast_id", "profile", "anchor", "observed_delta", "weight"),
            measure.vars = models, variable.name = "model", value.name = "prediction")
  d[, panel := factor(labels[as.character(model)], levels = labels[models])]
  annotation <- d[, .(r = wcorr(observed_delta, prediction, weight / sum(weight))), by = panel]
  annotation[, label := sprintf("Pearson r = %.2f", r)]
  ggplot(d, aes(observed_delta, prediction)) +
    geom_hline(yintercept = 0, linewidth = .3, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = .3, color = "grey75") +
    geom_point(size = .85, alpha = .22, color = "#0072B2") +
    geom_smooth(aes(weight = weight), method = "lm", formula = y ~ x, se = FALSE,
                linewidth = .65, color = "#D55E00") +
    geom_text(data = annotation, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -.04, vjust = 1.2, size = 3) +
    facet_wrap(~panel, nrow = 1, scales = "free") +
    labs(title = title, subtitle = subtitle, x = "Observed change in interaction (log2)",
         y = "Predicted change in interaction (log2)",
         caption = "Same anchor drug, dose, cell and laboratory; different partner drugs. Negative = greater suppression.\nPanels use independent scales; inspect profile plots for prediction amplitude.") +
    plot_theme + theme(aspect.ratio = 1)
}
primary_scatter <- contrast_scatter(primary_contrasts, c("molecular", "magnitude", "identity"),
                                    "Can the model distinguish partners for an excluded drug?",
                                    "Entire anchor drug excluded from functional fitting; contrasts nominated in other cells")
secondary_scatter <- contrast_scatter(secondary_contrasts, c("flexible", "additive_control", "flexible_norm"),
                                      "Does joint reading improve on individual response patterns?",
                                      "Separate experiment: cell and each pair excluded; constituent drugs may be seen")
ggsave(file.path(out, "figures", "whole_drug_partner_contrasts.png"), primary_scatter,
       width = 10.5, height = 4.2, dpi = 150)
ggsave(file.path(out, "figures", "cell_pair_partner_contrasts.png"), secondary_scatter,
       width = 10.5, height = 4.2, dpi = 150)
pdf(file.path(out, "figures", "composition_overview.pdf"), width = 10.5, height = 4.2)
print(primary_scatter)
print(secondary_scatter)
dev.off()

# A profile shows actual interaction magnitude and every available partner dose.
# Discovery sets row order; outcomes never select which rows appear.
partner_plot <- function(z, models, title, suffix = "") {
  z <- copy(z)
  setorder(z, discovery_median, partner_name, partner_dose, na.last = TRUE)
  z[, row_id := seq_len(.N)]
  z[, partner_label := sprintf("%s | %g uM  [%s]", partner_name, partner_dose,
                              fcase(support == "Strong in discovery" & discovery_median < 0, "extra suppression",
                                    support == "Strong in discovery" & discovery_median > 0, "less suppression",
                                    support == "Weak in discovery", "weak", default = "unselected"))]
  d <- melt(z, id.vars = c("row_id", "partner_label"), measure.vars = c("interaction", models),
            variable.name = "model", value.name = "value")
  d[model == "interaction", model := "observed"]
  d[, model := factor(as.character(model), levels = c("observed", models))]
  position <- position_dodge(width = .65)
  ggplot(d, aes(value, factor(row_id, levels = rev(z$row_id)), color = model, shape = model)) +
    geom_vline(xintercept = 0, linewidth = .6, color = "grey60") +
    geom_vline(xintercept = c(-1, 1), linewidth = .25, linetype = 3, color = "grey75") +
    geom_point(position = position, size = 2, alpha = .9) +
    scale_y_discrete(labels = setNames(z$partner_label, z$row_id)) +
    scale_color_manual(values = colors, labels = labels, drop = FALSE) +
    scale_shape_manual(values = setNames(c(16, 17, 15, 18), c("observed", models)), labels = labels) +
    labs(title = title, subtitle = paste0(z$regime[1], suffix), y = NULL,
         x = "Effect beyond independent single-agent abundance (log2 ratio)\n-1 = half as many cells; 0 = expected; +1 = twice as many",
         caption = "Bracketed labels refer to other-cell discovery, not the outcome shown. Every saved partner dose is retained.\nMeasured single-agent effects define the expectation; less suppression than independence need not mean rescue.") +
    plot_theme + theme(axis.text.y = element_text(size = 8), legend.text = element_text(size = 9))
}

# Write one complete collection, avoiding a proliferation of folders per profile.
render_profiles <- function(data, profile_ids, models, filename) {
  pdf(file.path(out, "figures", filename), width = 10, height = 8.5, onefile = TRUE)
  on.exit(dev.off())
  for (i in seq_along(profile_ids)) {
    z <- copy(data[profile == profile_ids[i]])
    setorder(z, discovery_median, partner_name, partner_dose, na.last = TRUE)
    pages <- split(seq_len(nrow(z)), ceiling(seq_len(nrow(z)) / 21L))
    for (j in seq_along(pages)) {
      title <- sprintf("%s %g uM + partners | %s | lab %s", z$anchor_name[1], z$anchor_dose[1],
                       z$CELLNAME[1], z$SCREENER[1])
      suffix <- sprintf(" | profile %d/%d, page %d/%d", i, length(profile_ids), j, length(pages))
      print(partner_plot(z[pages[[j]]], models, title, suffix))
    }
    if (i %% 40L == 0L) message(filename, ": ", i, "/", length(profile_ids), " profiles rendered")
  }
}
render_profiles(primary, display$profile, c("molecular", "magnitude", "identity"),
                "coverage_selected_partner_profiles.pdf")
preview <- primary[profile == display$profile[1]]
setorder(preview, discovery_median, partner_name, partner_dose, na.last = TRUE)
g <- partner_plot(preview[seq_len(min(.N, 21L))], c("molecular", "magnitude", "identity"),
                  sprintf("%s %g uM + partners | %s", preview$anchor_name[1], preview$anchor_dose[1], preview$CELLNAME[1]))
ggsave(file.path(out, "figures", "partner_profile_preview.png"), g, width = 10, height = 8.5, dpi = 150)
render_profiles(primary, coverage[eligible == TRUE, profile], c("molecular", "magnitude", "identity"),
                "all_whole_drug_partner_profiles.pdf")
secondary_profiles <- coverage[profile %in% secondary_contrasts$profile, profile]
render_profiles(secondary, secondary_profiles, c("flexible", "additive_control", "flexible_norm"),
                "cell_pair_control_profiles.pdf")

# Inspectable retrospective examples include every covered cell and anchor dose
# for the displayed partner concentrations. No favorable recipients are removed.
example_specs <- list(
  list(anchor = "idarubicin", partners = c("bortezomib", "vismodegib"), doses = c(.1, .3),
       title = "Idarubicin: partial recovery of a partner-dependent effect",
       note = "Bortezomib yields less suppression than independence; vismodegib is mostly near independence."),
  list(anchor = "mitoxantrone", partners = c("plicamycin", "oxaliplatin"), doses = c(.1, .1),
       title = "Mitoxantrone: right ordering can come from the wrong allocation",
       note = "Observed extra suppression is concentrated in plicamycin; the model often invents an oxaliplatin effect.")
)
example_data <- list()
pdf(file.path(out, "figures", "retrospective_partner_examples.pdf"), width = 11.5, height = 7)
for (k in seq_along(example_specs)) {
  spec <- example_specs[[k]]
  z <- primary[anchor == spec$anchor &
                 ((partner == spec$partners[1] & partner_dose == spec$doses[1]) |
                    (partner == spec$partners[2] & partner_dose == spec$doses[2]))]
  z[, example := spec$title]
  z[, partner_label := factor(paste0(partner_name, "\n", partner_dose, " uM"),
                               levels = unique(paste0(partner_name, "\n", partner_dose, " uM")))]
  z[, anchor_label := paste0("Anchor ", anchor_dose, " uM")]
  example_data[[k]] <- z
  d <- melt(z, id.vars = c("CELLNAME", "anchor_label", "partner_label"),
            measure.vars = c("interaction", "molecular", "magnitude", "identity"),
            variable.name = "model", value.name = "value")
  d[model == "interaction", model := "observed"]
  d[, model := factor(model, levels = c("observed", "molecular", "magnitude", "identity"))]
  p <- ggplot(d, aes(partner_label, value, color = model, shape = model, group = model)) +
    geom_hline(yintercept = 0, linewidth = .4, color = "grey55") +
    geom_line(linewidth = .4, alpha = .8) + geom_point(size = 1.6) +
    facet_grid(anchor_label ~ CELLNAME, drop = FALSE) +
    scale_color_manual(values = colors, labels = labels) +
    scale_shape_manual(values = c(observed = 16, molecular = 17, magnitude = 15, identity = 18), labels = labels) +
    scale_y_continuous(breaks = -4:4, labels = c("1/16x", "1/8x", "1/4x", "1/2x", "1x", "2x", "4x", "8x", "16x")) +
    labs(title = spec$title, subtitle = paste0(spec$note, "\nEntire anchor drug excluded from fitting; all covered cells and anchor doses shown."),
         x = "Partner and concentration", y = "Abundance relative to independent effects\n(log2-spaced axis)",
         caption = "Retrospective example selected after inspecting outcomes. Missing points indicate no saved matched prediction.\n1x = independent expectation; above1x = less suppression; below1x = more suppression. This is not evidence of rescue versus either drug alone.") +
    plot_theme + theme(axis.text.x = element_text(size = 7), strip.text = element_text(size = 8),
                       panel.spacing.x = grid::unit(.8, "lines"))
  print(p)
  ggsave(file.path(out, "figures", paste0("retrospective_", spec$anchor, ".png")), p,
         width = 11.5, height = 7, dpi = 160)
}
dev.off()
write_table(rbindlist(example_data), "retrospective_example_profiles.tsv.gz")

# Completed-run provenance and unchanged-input checks.
after_hashes <- vapply(files, function(f) digest(file = f, algo = "sha256"), character(1))
stopifnot(identical(input_hashes, after_hashes))
write_table(data.table(check = c("historical inputs unchanged", "whole anchor excluded", "no same-partner contrasts",
                                "interaction subtraction reconstructed", "display selection excluded evaluation labels"),
                      passed = TRUE), "checks.tsv")
artifacts <- list.files(file.path(out, "figures"), full.names = TRUE)
write_table(data.table(path = artifacts, bytes = file.info(artifacts)$size), "figure_index.tsv")
interpretation <- c(
  "INTERPRETATION",
  "The biological partner differences recur: 92.7% agreement with the independently",
  "nominated direction, weighting anchors/profiles/partner pairs/doses as specified.",
  "The whole-drug molecular model recovers a >=0.5 log2 contrast in only 38 of",
  "3700 substantial observed contrasts. These are dependent comparisons, not trials.",
  "Its median predicted contrast corresponds to 99% of the other partner's normalized",
  "abundance, versus 41% observed: it largely flattens the partner effects.",
  "This is not complete absence of ordering signal: within-profile centered Pearson",
  "is 0.147 across all eligible endpoints. The selected contrast correlation is ~0.",
  "Requiring the strong endpoint AND the near-independent endpoint to be correct",
  "leaves 21 recoveries among 3257 confirmed comparisons, spanning 7 profiles and",
  "5 anchors; the magnitude and identity controls recover none under this criterion.",
  "",
  "A recurring partial success is idarubicin with bortezomib versus vismodegib.",
  "At 1 uM idarubicin in HT29, bortezomib gives 4.92x the independent expectation",
  "and the model predicts 1.76x; vismodegib is near independence. This is less",
  "suppression than independence, not demonstrated rescue versus monotherapy.",
  "The figure retains both anchor doses and all covered cells. It was selected",
  "retrospectively and supports a limited example, not a general success rate.",
  "Mitoxantrone illustrates the converse: observed extra suppression is concentrated",
  "in plicamycin, whereas predicted contrasts can arise from an invented oxaliplatin",
  "effect. A correct ranking alone does not establish correct composition.",
  "",
  "The nonlinear model has modest ordering signal in the separate cell/pair test:",
  "within-profile centered Pearson 0.295 versus 0.245 for the individual-signature",
  "control. Both endpoints are recovered in 11/77 versus 7/77 confirmed comparisons.",
  "Those results do not establish whole-drug transfer or resolve the earlier",
  "single-drug identifiability failure. Combination outcomes provide an additional",
  "constraint, but these saved fits do not demonstrate a uniquely identified map.",
  "",
  "Grant implication: a limited retrospective example of molecular partner-effect",
  "transfer exists; broadly reliable composition and the baseline-to-combination",
  "bridge remain unproven. This audit uses the historical nine-program bank and",
  "measured single-drug molecular responses, not the newer bank or predicted responses."
)
writeLines(c(report, "", interpretation, "", "Complete: tables and all partner-profile figures generated; input hashes unchanged."),
           file.path(out, "REPORT.txt"))
message("Completed ", out)
