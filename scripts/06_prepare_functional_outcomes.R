#!/usr/bin/env Rscript
# Reconstruct PRISM curves and matched ALMANAC constituent/combination outcomes.
# The nine-drug shortlist and previously exposed dose choices are retrospective.
# Crizotinib's 5-uM context was added separately, with broader recipient coverage.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(qs)
})

output_dir <- "results/06_prepare_functional_outcomes"
data_dir <- file.path(output_dir, "data")
prism_dir <- "data/raw/prism"
native_doses <- c(0.05, 0.5, 5)
selection <- data.frame(
  drug = c(
    "5-Azacytidine", "Bendamustine", "Bortezomib", "Busulfan",
    "Cytarabine (hydrochloride)", "Lapatinib ditosylate", "Topotecan (hydrochloride)",
    "Vandetanib", "Vismodegib", "crizotinib"
  ),
  dose = c(5, 0.5, 0.5, 5, 0.5, 5, 5, 0.05, 0.5, 5),
  prism_name = c(
    "azacitidine", "bendamustine", "bortezomib", "busulfan", "cytarabine",
    "lapatinib", "topotecan", "vandetanib", "vismodegib", "crizotinib"
  )
)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(output_dir, "run.log"), split = TRUE)
write_result <- function(x, name) {
  fwrite(x, file.path(data_dir, paste0(name, ".tsv")), sep = "\t")
}
baseline <- qread("data/processed/expression.qs")$result$baseline
cohort <- read.delim("data/raw/annotations/historical_functional_cell_cohort.tsv")
stopifnot(!anyDuplicated(cohort$ModelID), nrow(cohort) == 1059L)

# Work from replicate wells. Average replicate log2 responses within each
# screen/compound/dose, interpolate in log dose, then weight across screens by
# the interpolated number of observed replicates. No extrapolation is allowed.
treatment <- fread(file.path(prism_dir, "secondary-screen-replicate-treatment-info.csv")) %>%
  filter(perturbation_type == "experimental_treatment", tolower(name) %in% selection$prism_name) %>%
  mutate(drug = selection$drug[match(tolower(name), selection$prism_name)])
header <- names(fread(file.path(prism_dir, "secondary-screen-logfold-change.csv"), nrows = 0))
raw <- fread(file.path(prism_dir, "secondary-screen-logfold-change.csv"),
  select = c(header[1], treatment$column_name)
)
setnames(raw, 1L, "depmap_id")
long <- melt(raw,
  id.vars = "depmap_id", variable.name = "column_name", value.name = "y",
  variable.factor = FALSE
) %>%
  inner_join(select(treatment, column_name, drug, dose, screen_id, broad_id), by = "column_name") %>%
  filter(is.finite(y), depmap_id %in% rownames(baseline), depmap_id != "ACH-000248")
rm(raw)
points <- long %>%
  group_by(depmap_id, drug, screen_id, broad_id, dose) %>%
  summarise(y = mean(y), replicates = n(), .groups = "drop")
write_result(points, "prism_replicate_curve_points")

# The earlier nine contexts required a finite measured endpoint at 0.15625,
# 0.625 or 2.5 uM in the
# original solid-tumor cohort. The later crizotinib preparation did not impose it.
eligible <- points %>%
  filter(vapply(
    dose, function(value) min(abs(log(value / c(0.15625, 0.625, 2.5)))) < 1e-4,
    logical(1)
  ), depmap_id %in% cohort$ModelID) %>%
  distinct(depmap_id, drug)
points <- points %>%
  mutate(key = paste(depmap_id, drug, sep = "|")) %>%
  filter(key %in% paste(eligible$depmap_id, eligible$drug, sep = "|") | drug == "crizotinib")
curves <- split(points, paste(points$depmap_id, points$drug, points$screen_id, points$broad_id, sep = "|"))
native <- bind_rows(lapply(curves, function(curve) {
  curve <- curve[order(curve$dose), ]
  if (nrow(curve) < 2L) {
    return(NULL)
  }
  data.frame(
    depmap_id = curve$depmap_id[1], drug = curve$drug[1],
    screen_id = curve$screen_id[1], broad_id = curve$broad_id[1], dose = native_doses,
    y = approx(log(curve$dose), curve$y, xout = log(native_doses), rule = 1)$y,
    weight = approx(log(curve$dose), curve$replicates, xout = log(native_doses), rule = 1)$y
  )
}))
write_result(native, "prism_native_dose_screen_estimates")
native <- native %>%
  filter(is.finite(y), is.finite(weight)) %>%
  group_by(depmap_id, drug, dose) %>%
  summarise(y = weighted.mean(y, weight), .groups = "drop")
write_result(native, "prism_native_dose_outcomes")

# Reproduce the original upfront dose selection for the nine retained contexts:
# maximize the equal-program mean absolute cross-cell Pearson association.
# All recipients' labels enter this historical choice, so subsequent performance
# remains retrospective. Crizotinib's later fixed 5-uM context is separate.
prepared <- qread("data/processed/molecular_inputs.qs")$result
programs <- sprintf("P%02d", 1:40)
dose_profiles <- prepared$m %>%
  filter(drug %in% selection$drug[1:9], depmap_id != "ACH-000248") %>%
  inner_join(native, by = c("depmap_id", "drug", "dose"), relationship = "one-to-one")
profile_groups <- split(dose_profiles, paste(dose_profiles$drug, dose_profiles$dose, sep = "|"))
dose_correlations <- bind_rows(lapply(profile_groups, function(profile) {
  data.frame(
    drug = profile$drug[1], dose = profile$dose[1], program = programs,
    cells = nrow(profile), r = as.numeric(cor(as.matrix(profile[, programs]), -profile$y))
  )
}))
stopifnot(all(is.finite(dose_correlations$r)))
dose_scores <- dose_correlations %>%
  group_by(drug, dose) %>%
  summarise(
    cells = first(cells), mean_abs_r = mean(abs(r)), median_abs_r = median(abs(r)),
    strongest_abs_r = max(abs(r)), .groups = "drop"
  ) %>%
  arrange(drug, desc(mean_abs_r), dose)
chosen_doses <- dose_scores %>%
  group_by(drug) %>%
  slice_head(n = 1) %>%
  ungroup()
index <- match(selection$drug[1:9], chosen_doses$drug)
stopifnot(!anyNA(index), all(selection$dose[1:9] == chosen_doses$dose[index]))
selection$dose[1:9] <- chosen_doses$dose[index]
write_result(dose_correlations, "observed_program_viability_each_dose_correlations")
write_result(dose_scores, "upfront_dose_selection_scores")

prism_selected <- inner_join(native, selection[, c("drug", "dose")], by = c("drug", "dose")) %>%
  mutate(
    assay = "secondary", source = "prism_secondary",
    control_key = paste(assay, depmap_id, drug, sprintf("%.10g", dose), sep = "|")
  )

# Map molecular formulations to published ALMANAC NSCs using the same audited
# synonym/salt rules. Prodrugs and stereoisomers are not merged by these rules.
canonical <- function(x) gsub("[^a-z0-9]", "", tolower(x))
parent_name <- function(x) {
  x <- tolower(x)
  x <- gsub("\\(|\\)", "", x)
  x <- gsub(paste(c(
    " hydrochloride", " dihydrochloride", " ditosylate", " sulfate hydrate",
    " sulfate", " disodium", " hemihydrate", " trihydrate", " tosylate", " free base"
  ), collapse = "|"), "", x)
  canonical(x)
}
conditions <- fread("data/raw/tahoe/pdex_selected/perturbations.tsv")
availability <- fread("data/raw/tahoe/pdex_selected/response_availability.tsv")
discovery_cells <- qread("data/processed/expression.qs")$result$eligibility_cells
coverage <- availability %>%
  filter(available == 1, depmap_id %in% discovery_cells) %>%
  left_join(select(conditions, condition_index, drug), by = "condition_index") %>%
  count(drug, name = "profiles")
compound_names <- fread("data/raw/almanac/ComboCompoundNames_small.txt",
  header = FALSE,
  col.names = c("nsc", "nci_name")
) %>%
  mutate(parent = parent_name(nci_name), nci_key = canonical(nci_name))
aliases <- c(
  `5fluorouracil` = "fluorouracil", `5azacytidine` = "azacitidine",
  rapamycin = "sirolimus", epirubicin = "4epiadriamycin", retinoicacid = "tretinoin"
)
formulations <- data.frame(drug = unique(conditions$drug)) %>%
  mutate(parent = parent_name(drug), molecular_key = canonical(drug))
index <- match(formulations$parent, names(aliases))
formulations$parent[!is.na(index)] <- unname(aliases[index[!is.na(index)]])
drug_mapping <- inner_join(formulations, compound_names, by = "parent", relationship = "many-to-many") %>%
  left_join(coverage, by = "drug") %>%
  filter(!is.na(profiles)) %>%
  mutate(strict = molecular_key == nci_key, identity = parent) %>%
  arrange(nsc, desc(strict), desc(profiles), drug) %>%
  group_by(nsc) %>%
  slice_head(n = 1) %>%
  ungroup()
write_result(drug_mapping, "almanac_drug_mapping")

# Match cell names against the original DepMap source, with explicit audited aliases.
raw_almanac <- fread("data/raw/almanac/ComboDrugGrowth_Nov2017.csv")
stopifnot(!anyDuplicated(raw_almanac$COMBODRUGSEQ), all(raw_almanac$VALID == "Y"))
models <- fread("data/raw/depmap/Model.csv") %>%
  transmute(cell_key = canonical(CellLineName), depmap_id = ModelID)
cell_mapping <- raw_almanac %>%
  distinct(CELLNAME) %>%
  mutate(cell_key = canonical(gsub("/ATCC|\\(TB\\)", "", CELLNAME)))
cell_aliases <- c(
  ovcar3 = "nihovcar3", nciadrres = "nciadrres", u251 = "u251mg",
  `7860` = "786o", mdamb435 = "mdamb435s"
)
index <- match(cell_mapping$cell_key, names(cell_aliases))
cell_mapping$cell_key[!is.na(index)] <- unname(cell_aliases[index[!is.na(index)]])
cell_mapping <- left_join(cell_mapping, models, by = "cell_key")
stopifnot(!anyDuplicated(cell_mapping$CELLNAME))
write_result(cell_mapping, "almanac_cell_mapping")

# Signal/control gives abundance relative to matched control, not percent growth.
# Match constituent singles within site, study, date, cell, compound, sample and dose.
raw_almanac <- raw_almanac %>% mutate(v = TESTVALUE / CONTROLVALUE)
assay_keys <- c("SCREENER", "STUDY", "TESTDATE", "CELLNAME")
# Preserve original row sorting and grouped arithmetic before deduplicating controls.
a <- as.data.table(raw_almanac)
keys <- assay_keys
mp <- as.data.table(drug_mapping)
setnames(mp, "drug", "tahoe")
mp[, primary := TRUE]
lm <- as.data.table(cell_mapping)
lm[, observed_tahoe := depmap_id %in% discovery_cells]
AU <- "ACH-000248"
s <- a[is.na(NSC2) & CONCUNIT1 == "M" & is.finite(v) & v > 0,
  .(v = mean(v), pg = mean(PERCENTGROWTH), single_n = .N),
  by = c(keys, "NSC1", "SAMPLE1", "CONC1")]
p <- a[!is.na(NSC2) & NSC1 != NSC2 & NSC1 %in% mp$nsc & NSC2 %in% mp$nsc & CONCUNIT1 == "M" & CONCUNIT2 == "M"]
s1 <- copy(s)
setnames(s1, c("v", "pg", "single_n"), c("vA", "pgA", "nA"))
p <- merge(p, s1, by = c(keys, "NSC1", "SAMPLE1", "CONC1"), all.x = TRUE)
s2 <- copy(s)
setnames(s2, c("NSC1", "SAMPLE1", "CONC1", "v", "pg", "single_n"), c("NSC2", "SAMPLE2", "CONC2", "vB", "pgB", "nB"))
p <- merge(p, s2, by = c(keys, "NSC2", "SAMPLE2", "CONC2"), all.x = TRUE)
# Never borrow singles across studies/dates; study-level matching is audited against NCI ExpectedGrowth.
p <- p[is.finite(vA) & is.finite(vB) & is.finite(v) & v > 0]
p[, `:=`(y = log2(v), yA = log2(vA), yB = log2(vB), doseA = CONC1 * 1e6, doseB = CONC2 * 1e6)]
p[, `:=`(I = y - yA - yB, expected_pg = ifelse(pmin(pgA, pgB) < 0, pmin(pgA, pgB), pmin(pgA, 100) * pmin(pgB, 100) / 100))]
p <- p[abs(expected_pg - EXPECTEDGROWTH) < .02]
mp <- mp[primary == TRUE]
p <- merge(p, mp[, .(NSC1 = nsc, drugA = tahoe, idA = identity, strictA = strict)], by = "NSC1")
p <- merge(p, mp[, .(NSC2 = nsc, drugB = tahoe, idB = identity, strictB = strict)], by = "NSC2")
p <- merge(p, lm, by = "CELLNAME")
p <- p[!is.na(depmap_id) & depmap_id != AU & idA != idB]
p[, pair := paste(pmin(idA, idB), pmax(idA, idB), sep = "::")]
# Canonicalize identities AND doses before aggregating independent assay repeats.
reverse <- p$idA > p$idB
for (nm in list(c("drugA", "drugB"), c("idA", "idB"), c("doseA", "doseB"),
  c("yA", "yB"), c("vA", "vB"), c("pgA", "pgB"), c("strictA", "strictB"))) {
  tmp <- p[[nm[1]]][reverse]
  set(p, which(reverse), nm[1], p[[nm[2]]][reverse])
  set(p, which(reverse), nm[2], tmp)
}
# Preserve measured dose points. Collapse only true replicate line/pair/dose/laboratory records.
d <- p[, .(y = mean(y), yA = mean(yA), yB = mean(yB), I = mean(I),
  pg = mean(PERCENTGROWTH), expected_pg = mean(expected_pg), v = mean(v),
  vA = mean(vA), vB = mean(vB), replicates = .N),
  by = .(pair, depmap_id, CELLNAME, PANEL, SCREENER, drugA, drugB, idA, idB,
    doseA, doseB, strictA, strictB, observed_tahoe)]
pairs <- as.data.frame(d) %>%
  filter(
    depmap_id %in% cohort$ModelID, doseA >= 0.05 - 1e-10, doseA <= 5 + 1e-10,
    doseB >= 0.05 - 1e-10, doseB <= 5 + 1e-10
  )
write_result(pairs, "almanac_matched_pairs")
# Preserve the historical tabular precision boundary before interpolation.
# In particular, molar-to-micromolar conversion otherwise leaves join keys such
# as 0.3 differing from their nominal decimal value by one floating-point unit.
pairs <- as.data.frame(fread(file.path(data_dir, "almanac_matched_pairs.tsv")))
saveRDS(pairs, "data/processed/almanac_pairs.rds")

# Interpolate each site's complete 2D log-dose grid to the original selected doses.
bracket <- function(grid, target) {
  grid <- sort(unique(grid))
  exact <- which(abs(grid - target) < 1e-10)
  if (length(exact)) {
    return(data.frame(dose = grid[exact[1]], weight = 1))
  }
  if (target < min(grid) || target > max(grid)) {
    return(NULL)
  }
  lower <- max(grid[grid < target])
  upper <- min(grid[grid > target])
  fraction <- log(target / lower) / log(upper / lower)
  data.frame(dose = c(lower, upper), weight = c(1 - fraction, fraction))
}
source_pairs <- pairs %>% filter(drugA %in% selection$drug[1:9], drugB %in% selection$drug[1:9])
groups <- split(source_pairs, paste(source_pairs$depmap_id, source_pairs$SCREENER, source_pairs$pair, sep = "|"))
selected_pairs <- bind_rows(lapply(groups, function(z) {
  target_a <- selection$dose[match(z$drugA[1], selection$drug)]
  target_b <- selection$dose[match(z$drugB[1], selection$drug)]
  a <- bracket(z$doseA, target_a)
  b <- bracket(z$doseB, target_b)
  if (is.null(a) || is.null(b)) {
    return(NULL)
  }
  corners <- merge(a, b, by = NULL, suffixes = c("A", "B")) %>%
    mutate(weight = weightA * weightB) %>%
    left_join(select(z, doseA, doseB, y, yA, yB), by = c("doseA", "doseB"))
  if (anyNA(corners$y) || anyNA(corners$yA) || anyNA(corners$yB)) {
    return(NULL)
  }
  stopifnot(abs(sum(corners$weight) - 1) < 1e-12)
  z[1, c("depmap_id", "drugA", "drugB", "pair", "CELLNAME", "PANEL", "SCREENER")] %>%
    mutate(
      doseA = target_a, doseB = target_b, y = sum(corners$weight * corners$y),
      yA = sum(corners$weight * corners$yA), yB = sum(corners$weight * corners$yB), I = y - yA - yB
    )
}))
write_result(selected_pairs, "combination_training_outcomes")
arms <- bind_rows(
  selected_pairs %>% transmute(depmap_id, pair, SCREENER, arm = "A", drug = drugA, dose = doseA, y = yA),
  selected_pairs %>% transmute(depmap_id, pair, SCREENER, arm = "B", drug = drugB, dose = doseB, y = yB)
)
# The historical unique() rule compared outcomes at exact floating-point
# precision. Preserve its admitted assay arms explicitly, rather than letting
# rounding decide whether a reused control appears once or twice. This input
# contains identifiers only; every outcome below is reconstructed from raw data.
retained_arms <- read.delim("data/raw/annotations/retained_almanac_control_arms.tsv")
constituents <- retained_arms %>%
  left_join(arms,
    by = c("depmap_id", "pair", "SCREENER", "arm", "drug", "dose"),
    relationship = "many-to-one"
  ) %>%
  transmute(depmap_id, drug, dose, y, assay = paste0("nci_", SCREENER)) %>%
  mutate(source = "nci_single", control_key = paste(assay, depmap_id, drug, sprintf("%.10g", dose), sep = "|"))
stopifnot(!anyNA(constituents), nrow(constituents) == 352L)
single <- bind_rows(prism_selected, constituents) %>%
  mutate(group = paste(assay, drug, dose, sep = "|"), row_id = row_number())
write_result(single, "individual_training_outcomes")
write_result(selection[, c("drug", "dose")], "selected_doses")
fold_cells <- sort(unique(c(single$depmap_id, selected_pairs$depmap_id)))
write_result(data.frame(depmap_id = fold_cells, fold = seq_along(fold_cells)), "recipient_folds")

# Compare rebuilt outcome measurements only after raw-data reconstruction.
original_dir <- "reference/results/93_joint_calibrated_sensitivity_routes"
comparison <- list()
for (name in c("individual_training_outcomes", "combination_training_outcomes")) {
  original <- read.delim(file.path(original_dir, paste0(name, ".tsv")))
  rebuilt <- if (name == "individual_training_outcomes") single else selected_pairs
  keys <- if (name == "individual_training_outcomes") {
    c("depmap_id", "drug", "dose", "assay")
  } else {
    c("depmap_id", "pair", "SCREENER", "doseA", "doseB")
  }
  rebuilt <- rebuilt %>%
    group_by(across(all_of(keys))) %>%
    arrange(y, .by_group = TRUE) %>%
    mutate(repeat_index = row_number()) %>%
    ungroup()
  original <- original %>%
    group_by(across(all_of(keys))) %>%
    arrange(y, .by_group = TRUE) %>%
    mutate(repeat_index = row_number()) %>%
    ungroup()
  joined <- inner_join(rebuilt, original,
    by = c(keys, "repeat_index"),
    suffix = c("_new", "_old"), relationship = "one-to-one"
  )
  comparison[[name]] <- data.frame(
    table = name, rebuilt_rows = nrow(rebuilt), reference_rows = nrow(original),
    matched_rows = nrow(joined), maximum_y_difference = max(abs(joined$y_new - joined$y_old))
  )
}
comparison <- bind_rows(comparison)
write_result(comparison, "reference_comparison")
print(comparison)
if (any(comparison$rebuilt_rows != comparison$reference_rows) || any(comparison$maximum_y_difference > 1e-8)) {
  stop("Functional reconstruction differs from the reference; inspect before fitting")
}

# Package ordinary data frames for the retained model equations. This object is
# generated here, never copied from the portable bundle.
archive <- list(
  prep = prepared, single = as.data.frame(single), pair = as.data.frame(pairs),
  drugs = qread("data/processed/molecular_gene_effects.qs")$result$drugs, meta = cohort
)
saveRDS(archive, "data/processed/molecular_functional_inputs.rds")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
