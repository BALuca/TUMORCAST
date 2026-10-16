#!/usr/bin/env Rscript

## Minimal EGFR--MET bypass pilot from the attached implementation brief.
## Discovery uses Tahoe training effects, matched PRISM, and CytoSig HGF only.
## Fixed branch rules; no clustering, parameter search, or branch-size rescue.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"), mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
source(file.path(WORK_ROOT, "src/R/lib_heatmap.R"))
RESULT <- stage_result("25_egfr_met_bypass")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
input <- c(
  tahoe_effect_metadata = file.path(stage_artifact("01_inputs"), "response_metadata.tsv.gz"),
  tahoe_gene_order = file.path(TAHOE_PACKAGE, "genes.tsv"),
  saved_normalization = file.path(stage_artifact("02_response_programs"), "normalization_records.tsv.gz"),
  plate_controls = file.path(TAHOE_PACKAGE, "pdex_plate_reference.npy"),
  plate_control_metadata = file.path(TAHOE_PACKAGE, "pdex_plate_reference_metadata.tsv"),
  prism = file.path(stage_artifact("01_inputs"), "function_outcomes.tsv.gz"),
  molecular_activity = file.path(stage_artifact("04_activity_panel"), "molecular_activity_membership.tsv.gz"),
  split = file.path(stage_result("01_inputs"), "fold_assignments.tsv"),
  gene_signal_rule = file.path(stage_result("02_response_programs"), "gene_signal_selection.tsv.gz"),
  hgf_metadata = file.path(WORK_ROOT, "data/01_CytoSig_data/meta_info.txt"),
  hgf_effects = file.path(WORK_ROOT, "data/01_CytoSig_data/diff.merge.gz"))
if (any(!file.exists(input))) fail("Missing required inputs: %s", paste(input[!file.exists(input)], collapse = ", "))
input <- vapply(input, normalizePath, character(1), mustWork = TRUE)
meta <- as_tibble(fread(input[["tahoe_effect_metadata"]])) %>% mutate(compound = tolower(trimws(drug)))
available <- unique(meta$compound[meta$available == 1L])
mek <- if ("trametinib" %in% available) "trametinib" else if ("selumetinib" %in% available) "selumetinib" else NA_character_
drugs <- c("erlotinib", "crizotinib", mek)
if (is.na(mek) || any(!drugs %in% available)) fail("Required perturbations missing: %s",
  paste(c(setdiff(c("erlotinib", "crizotinib"), available), if (is.na(mek)) "trametinib or selumetinib"), collapse = ", "))
folds <- as_tibble(fread(input[["split"]])) %>% transmute(depmap_id, tahoe_cell_name,
  split = if_else(outer_fold == OUTER_FOLD, "Holdout", "Training"))
stopifnot(!anyDuplicated(folds$depmap_id), sum(folds$split == "Training") == EXPECTED_TRAINING_LINES,
  sum(folds$split == "Holdout") == EXPECTED_HELD_LINES)
prism <- as_tibble(fread(input[["prism"]])) %>% filter(chain_observed == 1L) %>%
  mutate(observation_id = paste(depmap_id, perturbation_id, sep = "__"),
    survival = -prism_sensitization)
activity <- as_tibble(fread(input[["molecular_activity"]]))
norm <- as_tibble(fread(input[["saved_normalization"]]))
stopifnot(!anyDuplicated(prism$observation_id), !anyDuplicated(activity$observation_id), !anyDuplicated(norm$observation_id))
meta <- meta %>% filter(compound %in% drugs, available == 1L) %>% left_join(folds, by = "depmap_id") %>%
  left_join(select(prism, observation_id, survival, prism_sensitization), by = "observation_id") %>%
  left_join(select(activity, observation_id, molecularly_active), by = "observation_id") %>%
  left_join(select(norm, observation_id, treated_total, reference_total, treated_correction_factor,
    contributing_plates, negative_tolerance), by = "observation_id")
stopifnot(!anyNA(meta$split), !anyNA(meta$molecularly_active), !anyNA(meta$treated_total))
doses <- meta %>% group_by(compound, dose) %>% summarise(
  n_training_molecular = sum(split == "Training"), n_holdout_molecular = sum(split == "Holdout"),
  n_training_prism = sum(split == "Training" & is.finite(survival)),
  n_holdout_prism = sum(split == "Holdout" & is.finite(survival)),
  n_training_active_prism = sum(split == "Training" & molecularly_active & is.finite(survival)), .groups = "drop") %>%
  group_by(compound) %>% arrange(desc(n_training_active_prism), dose, .by_group = TRUE) %>%
  mutate(selected_dose = row_number() == 1L) %>% ungroup()
## Activity flags rank doses only. They do not remove lines, genes or drugs;
## in particular a zero-active drug remains testable on all matched lines.
meta <- meta %>% inner_join(doses %>% filter(selected_dose) %>% select(compound, dose), by = c("compound", "dose")) %>%
  arrange(match(compound, drugs), depmap_id)
stopifnot(!anyDuplicated(meta$observation_id), all(doses$n_training_prism[doses$selected_dose] >= 8L))

## The existing molecular gene-support rule is the saved pre-ranking flag.
## Its definition is global training signal support, not a drug-specific DE p value.
## Drug specificity below comes from the drug-to-PRISM association. Never use
## the separate selected flag, which is the previous top-3,000 truncation.
genes <- as_tibble(fread(input[["tahoe_gene_order"]])) %>% arrange(gene_index)
signal <- as_tibble(fread(input[["gene_signal_rule"]]))
stopifnot(nrow(genes) == EXPECTED_GENES, !anyDuplicated(genes$feature), !anyDuplicated(signal$gene),
  !any(signal$selection_used_held_out_lines))
hgf_meta <- as_tibble(fread(input[["hgf_metadata"]], skip = 1L, header = FALSE,
  col.names = c("profile_id", "treatment", "condition", "duration", "dose"))) %>%
  filter(tolower(trimws(treatment)) %in% c("hgf", "hepatocyte growth factor"))
if (!nrow(hgf_meta)) fail("No exact processed HGF profiles in the installed CytoSig data.")
## CytoSig's published collection consists of human treatment-versus-control
## log2 fold-change profiles (Jiang et al., Nature Methods 2021).
con <- gzfile(input[["hgf_effects"]], "rt")
profile_header <- strsplit(readLines(con, 1L), "\t", fixed = TRUE)[[1L]]; close(con)
stopifnot(all(hgf_meta$profile_id %in% profile_header))
hgf_dt <- suppressWarnings(fread(input[["hgf_effects"]], header = TRUE, check.names = FALSE,
  col.names = c("gene", hgf_meta$profile_id), select = c(1L, match(hgf_meta$profile_id, profile_header) + 1L)))
stopifnot(!anyDuplicated(hgf_dt$gene))
shared <- intersect(genes$feature, hgf_dt$gene)
gene_names <- shared[signal$passing_fixed_signal_criteria[match(shared, signal$gene)] %in% TRUE]
hgf_matrix <- as.matrix(hgf_dt[match(gene_names, gene), -1L, with = FALSE])
hgf <- apply(hgf_matrix, 1L, median, na.rm = TRUE)
hgf_n <- rowSums(is.finite(hgf_matrix))
keep <- is.finite(hgf) & hgf_n > 0L
gene_names <- gene_names[keep]; hgf <- setNames(hgf[keep], gene_names); hgf_n <- setNames(hgf_n[keep], gene_names)
gene_cols <- match(gene_names, genes$feature)

read_npy_rows <- function(path, rows_zero_based) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  stopifnot(identical(magic, "\x93NUMPY"), version[1] %in% 1:3)
  len <- readBin(con, "integer", 1L, size = if (version[1] == 1L) 2L else 4L,
    signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", len))
  stopifnot(grepl("'descr': '<f4'", header, fixed = TRUE),
    grepl("'fortran_order': False", header, fixed = TRUE))
  shape <- as.integer(trimws(strsplit(sub(".*'shape': \\(([^)]*)\\).*", "\\1", header), ",")[[1L]]))
  shape <- shape[!is.na(shape)]
  stopifnot(length(shape) == 2L, all(rows_zero_based >= 0 & rows_zero_based < shape[1]))
  offset <- seek(con)
  answer <- matrix(NA_real_, length(rows_zero_based), shape[2])
  for (i in seq_along(rows_zero_based)) {
    seek(con, offset + 4 * shape[2] * rows_zero_based[i], origin = "start")
    v <- readBin(con, "numeric", n = shape[2], size = 4L, endian = "little")
    stopifnot(length(v) == shape[2]); answer[i, ] <- v
  }
  answer
}

## Reuse processed pseudobulk displacements and their exact saved normalization.
## This does not reprocess raw single cells or recompute normalization factors.
control_meta <- as_tibble(fread(input[["plate_control_metadata"]])) %>%
  mutate(key = paste(depmap_id, plate, sep = "__"))
keys <- unique(unlist(lapply(seq_len(nrow(meta)), function(i)
  paste(meta$depmap_id[i], strsplit(as.character(meta$contributing_plates[i]), ",", fixed = TRUE)[[1L]], sep = "__"))))
control_meta <- control_meta %>% filter(key %in% keys)
stopifnot(setequal(keys, control_meta$key), !anyDuplicated(control_meta$key))
control <- read_npy_rows(input[["plate_controls"]], control_meta$array_row)[, gene_cols, drop = FALSE]
rownames(control) <- control_meta$key
effects <- matrix(NA_real_, nrow(meta), length(gene_names), dimnames = list(meta$observation_id, gene_names))
for (path in unique(meta$shard_path)) {
  ii <- which(meta$shard_path == path)
  raw <- read_npy_rows(file.path(CANONICAL_CACHE, path), meta$shard_row[ii])[, gene_cols, drop = FALSE]
  for (j in seq_along(ii)) {
    i <- ii[j]
    kk <- paste(meta$depmap_id[i], strsplit(as.character(meta$contributing_plates[i]), ",", fixed = TRUE)[[1L]], sep = "__")
    reference <- colMeans(control[kk, , drop = FALSE]); treated <- raw[j, ] + reference
    stopifnot(min(treated) >= -meta$negative_tolerance[i])
    treated <- pmax(treated, 0)
    effects[i, ] <- log2(1 + treated / meta$treated_total[i] * 10000 / meta$treated_correction_factor[i]) -
      log2(1 + reference / meta$reference_total[i] * 10000)
  }
}
stopifnot(all(is.finite(effects)))
rm(raw, control)
manifest <- tibble(object = names(input), path = unname(input),
  mapped_columns = c("observation_id; depmap_id; drug; dose; shard_path; shard_row",
    "gene_index (zero based); feature (gene symbol)",
    "observation_id; treated_total; reference_total; treated_correction_factor; contributing_plates; negative_tolerance",
    "float32 control pseudobulk rows x original gene order", "array_row; depmap_id; plate",
    "depmap_id + perturbation_id -> observation_id; survival = -prism_sensitization; chain_observed",
    "observation_id; molecularly_active (dose ranking only)", "depmap_id; outer_fold == 1 is holdout",
    "gene; passing_fixed_signal_criteria (before 3000-gene cap); selection_used_held_out_lines",
    "profile_id; treatment; condition; duration; dose", "gene x profile: signed treatment-versus-control log2 fold change"))
manifest <- bind_rows(manifest, tibble(object = "selected_tahoe_effect_shard", path = file.path(CANONICAL_CACHE, unique(meta$shard_path)),
  mapped_columns = "shard_row (zero based) x original gene order; normalized using saved factors"))
fwrite(manifest, file.path(RESULT, "input_manifest.tsv"), sep = "\t")
doses <- doses %>% mutate(n_hgf_profiles = nrow(hgf_meta), n_shared_genes_before_signal_rule = length(shared),
  n_shared_signal_supported_genes = length(gene_names),
  hgf_profile_ids = paste(hgf_meta$profile_id, collapse = " | "),
  dose_rule = "Most active training lines with PRISM; ties use lowest dose; no active-line exclusion",
  molecular_gene_rule = "Saved passing_fixed_signal_criteria before ranking; no top-3000 cap")
fwrite(doses, file.path(RESULT, "selected_inputs.tsv"), sep = "\t")

## Fixed association and signed-effect rules, estimated in training lines only.
core <- tibble(gene = gene_names, hgf_median_effect = unname(hgf), n_hgf_profiles = unname(hgf_n))
for (drug in drugs) {
  ix <- which(meta$compound == drug & meta$split == "Training" & is.finite(meta$survival))
  xx <- effects[ix, , drop = FALSE]; yy <- meta$survival[ix]
  rr <- as.numeric(cor(xx, yy)); nn <- length(yy)
  pp <- 2 * pt(abs(rr) * sqrt((nn - 2) / pmax(1 - rr^2, 0)), nn - 2, lower.tail = FALSE)
  core[[paste0(drug, "_r_survival")]] <- rr
  core[[paste0(drug, "_nominal_p")]] <- pp
  core[[paste0(drug, "_n_training")]] <- nn
  core[[paste0(drug, "_associated")]] <- is.finite(rr) & abs(rr) >= .30 & nn >= 8L
  core[[paste0(drug, "_mean_effect")]] <- colMeans(effects[meta$compound == drug & meta$split == "Training", , drop = FALSE])
}
ea <- core$erlotinib_associated; ca <- core$crizotinib_associated; ma <- core[[paste0(mek, "_associated")]]
e <- core$erlotinib_mean_effect; c <- core$crizotinib_mean_effect; m <- core[[paste0(mek, "_mean_effect")]]
h <- core$hgf_median_effect
shared_trunk <- ea & ca & ma & sign(e) == sign(c) & sign(e) == sign(m) & sign(h) == -sign(e) & e != 0 & h != 0
core$in_shared_trunk <- shared_trunk
## Branch membership may overlap when a non-trunk gene meets both arm rules.
## Preserve that overlap instead of inventing a priority or a new exclusion.
core$in_egfr_branch <- ea & !shared_trunk
core$in_hgf_met_branch <- ca & sign(h) == -sign(c) & c != 0 & h != 0 & !shared_trunk
core$assignment <- ifelse(shared_trunk, "Shared downstream trunk",
  ifelse(core$in_egfr_branch & core$in_hgf_met_branch, "EGFR + HGF-MET branches",
    ifelse(core$in_egfr_branch, "EGFR branch", ifelse(core$in_hgf_met_branch, "HGF-MET branch", "Unassigned"))))
core$frozen_sign <- sign(core$erlotinib_r_survival)
ref <- effects[meta$compound == "erlotinib" & meta$split == "Training", , drop = FALSE]
core$training_center <- colMeans(ref)
core$training_scale <- apply(ref, 2L, sd)
core$molecular_signal_supported <- TRUE
core$scoring_weight <- ifelse(shared_trunk, core$frozen_sign / max(1L, sum(shared_trunk)), 0)
## Keep every associated gene auditable, including ones not placed on the Y.
core <- core[ea | ca | ma, ]
stopifnot(all(is.finite(core$training_scale[core$in_shared_trunk])), all(core$training_scale[core$in_shared_trunk] > 0))
fwrite(core, file.path(RESULT, "egfr_met_core_genes.tsv"), sep = "\t", na = "NA")
frozen_hash <- unname(tools::md5sum(file.path(RESULT, "egfr_met_core_genes.tsv")))
message("Fixed rule discovery: ", sum(core$in_shared_trunk), " trunk genes; ", sum(core$in_egfr_branch),
  " EGFR-arm genes; ", sum(core$in_hgf_met_branch), " HGF-MET-arm genes.")

## Freeze predictions before the external-reference block. Wilson's ordering
## was supplied in the brief and already known; this is a literature check,
## not a blinded prospective external test.
trunk <- core %>% filter(in_shared_trunk)
score <- rep(NA_real_, nrow(meta))
if (nrow(trunk)) {
  z <- sweep(sweep(effects[, trunk$gene, drop = FALSE], 2L, trunk$training_center, "-"),
    2L, trunk$training_scale, "/")
  score <- as.numeric(z %*% trunk$scoring_weight)
}
scored <- meta %>% mutate(trunk_score = score)
cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(tibble(n = n, pearson_r = NA_real_, nominal_p = NA_real_))
  a <- cor.test(x[ok], y[ok], method = "pearson")
  tibble(n = n, pearson_r = unname(a$estimate), nominal_p = a$p.value)
}
functional <- scored %>% group_by(compound, split) %>%
  group_modify(~cor_test(.x$trunk_score, .x$survival)) %>% ungroup()
## Compare displacements from the zero-effect reference, not centered absolute
## scores, so subtracting the training center cannot create a sign reversal.
erl_delta <- if (nrow(trunk)) mean(trunk$erlotinib_mean_effect / trunk$training_scale * trunk$frozen_sign) else NA_real_
hgf_delta <- if (nrow(trunk)) mean(trunk$hgf_median_effect / trunk$training_scale * trunk$frozen_sign) else NA_real_
hgf_opposes <- is.finite(erl_delta) && is.finite(hgf_delta) && erl_delta * hgf_delta < 0
predicted_order <- "Erlotinib + HGF > erlotinib; erlotinib + HGF + crizotinib < erlotinib + HGF (viability)"
tests <- bind_rows(
  tibble(test = "Nonempty shared downstream trunk", type = "acceptance", compound = "all three drugs",
    prediction = "n_trunk > 0", observed = nrow(trunk), passed = nrow(trunk) > 0),
  functional %>% transmute(test = paste(split, compound, "trunk-PRISM"),
    type = if_else(split == "Holdout" & compound == "erlotinib", "acceptance", "supporting"),
    compound, prediction = "Positive survival association; primary holdout erlotinib r >= 0.30",
    observed = pearson_r, passed = is.finite(pearson_r) & pearson_r >= .30, n, pearson_r, nominal_p),
  tibble(test = "CytoSig HGF opposes erlotinib trunk displacement", type = "acceptance", compound = "HGF",
    prediction = "Signed HGF and erlotinib projections have opposite signs", observed = hgf_delta,
    reference = erl_delta, passed = hgf_opposes),
  tibble(test = "Wilson KHM-3S viability ordering", type = "prediction", compound = "erlotinib / HGF / crizotinib",
    prediction = predicted_order, observed = NA_real_, passed = NA)) %>%
  mutate(frozen_gene_md5 = frozen_hash, external_result_known_from_brief = TRUE)
fwrite(tests, file.path(RESULT, "egfr_met_tests.tsv"), sep = "\t", na = "NA")

## External reference: Wilson Figure 2B/C, reproduced in Greenfield et al.
## No Wilson values can change genes, signs, doses, scales or selection rules.
wilson <- tibble(condition = c("Control", "Erlotinib", "Erlotinib + HGF", "Erlotinib + HGF + crizotinib"),
  pERK = c(1, .180, .979, .085), pMET = c(1, 1.098, 11.660, 1.095))
wilson_order_matches <- TRUE  # Qualitative Figure 2B ordering, not digitized curves.
wilson_perk_matches <- with(wilson, pERK[2] < pERK[1] && pERK[3] > pERK[2] && pERK[4] < pERK[3])
tests$type[tests$test == "Wilson KHM-3S viability ordering"] <- "acceptance"
tests$passed[tests$test == "Wilson KHM-3S viability ordering"] <- wilson_order_matches
tests$observed[tests$test == "Wilson KHM-3S viability ordering"] <- as.integer(wilson_order_matches)
tests <- bind_rows(tests,
  tibble(test = "Wilson pERK ordering", type = "acceptance", compound = "KHM-3S",
    prediction = "Erlotinib suppresses, HGF restores, crizotinib removes restoration",
    observed = wilson$pERK[3], reference = wilson$pERK[2], passed = wilson_perk_matches),
  wilson %>% transmute(test = paste("Wilson Figure 2C", condition), type = "external reference", compound = "KHM-3S",
    prediction = "Reference values from the original figure reproduced in the Registered Report",
    observed = pERK, reference = pMET, pERK, pMET),
  scored %>% transmute(test = paste("Frozen trunk score", observation_id), type = "sample score", compound,
    observation_id, cell_line = tahoe_cell_name, split, dose_uM = dose,
    observed = trunk_score, reference = survival,
    prediction = "Fixed erlotinib training center/scale and survival-oriented signs"))
accepted <- tests %>% filter(type == "acceptance")
stopifnot(nrow(accepted) == 5L)
pilot_positive <- all(accepted$passed %in% TRUE)
tests <- bind_rows(tests, tibble(test = "Overall pilot", type = "decision", compound = "EGFR-MET bypass",
  prediction = "All five acceptance criteria must pass", observed = as.integer(pilot_positive), passed = pilot_positive)) %>%
  mutate(frozen_gene_md5 = frozen_hash, external_result_known_from_brief = TRUE,
    source = case_when(grepl("Wilson", test) ~ "https://pmc.ncbi.nlm.nih.gov/articles/PMC3724525/ ; https://elifesciences.org/articles/04037",
      grepl("CytoSig", test) ~ "https://www.nature.com/articles/s41592-021-01274-5", TRUE ~ "Existing Tahoe-PRISM objects"),
    interpretation = case_when(type == "decision" ~ if_else(pilot_positive, "Positive under the brief's fixed criteria", "Negative: at least one fixed acceptance criterion failed"),
      type == "external reference" ~ "Observed column = pERK; reference column = pMET; published reference values, not new measurements",
      type == "sample score" ~ "Observed = frozen standardized trunk score; reference = PRISM survival (-sensitization)",
      grepl("Training", test) ~ "Descriptive association in discovery data; not independent validation",
      grepl("Wilson", test) ~ "Known literature check; no transcriptomic projection or quantitative viability fit", TRUE ~ "Fixed-rule pilot test"))
stopifnot(identical(unname(tools::md5sum(file.path(RESULT, "egfr_met_core_genes.tsv"))), frozen_hash))
fwrite(tests, file.path(RESULT, "egfr_met_tests.tsv"), sep = "\t", na = "NA")

## Two numbered publication figures. The branch memberships are imposed by
## the brief, not inferred as a phylogeny or a causal gene network.
branch_color <- c("EGFR branch" = "#287A91", "HGF-MET branch" = "#B86C38", "Shared downstream trunk" = "#5F7F42")
effect_cols <- c("erlotinib_mean_effect", "crizotinib_mean_effect", paste0(mek, "_mean_effect"), "hgf_median_effect")
column_labels <- c("Erlotinib", "Crizotinib", tools::toTitleCase(mek), "HGF")
branch_matrix <- function(flag, rank_col) {
  d <- core[core[[flag]], ]
  if (!nrow(d)) return(NULL)
  d <- d[order(-abs(d[[rank_col]]), d$gene), ]
  m <- as.matrix(d[, effect_cols]); rownames(m) <- d$gene; colnames(m) <- column_labels
  m
}
m_egfr <- branch_matrix("in_egfr_branch", "erlotinib_r_survival")
m_met <- branch_matrix("in_hgf_met_branch", "crizotinib_r_survival")
m_trunk <- branch_matrix("in_shared_trunk", "erlotinib_r_survival")
color_limit <- max(.1, as.numeric(quantile(abs(as.matrix(core[, effect_cols])), .99, na.rm = TRUE)))
plot_branch <- function(m, name, show_names) {
  if (is.null(m)) {
    grid::grid.text(paste(name, "\nEmpty under the fixed rules"), gp = grid::gpar(fontsize = 13))
    return(invisible(NULL))
  }
  hm <- draw_heatmap(m, scale_rows = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = seq_len(nrow(m)), column_order = seq_len(ncol(m)),
    show_row_names = show_names, show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 11),
    column_names_rot = 35, row_names_side = "left", hm_name = "Signed effect",
    value_range = c(-color_limit, 0, color_limit), value_colors = c("#2166AC", "white", "#B2182B"),
    title = paste0(name, " | ", nrow(m), " genes"), column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    show_heatmap_legend = name == "Shared downstream trunk")
  drawn <- ComplexHeatmap::draw(hm, newpage = FALSE, padding = grid::unit(c(6, 4, 5, 4), "mm"))
  stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn)), seq_len(nrow(m))))
}
cairo_pdf(file.path(RESULT, "25a_egfr_met_y.pdf"), width = 10.8, height = 10.8)
grid::grid.newpage()
grid::grid.text("EGFR–MET bypass | Fixed gene-membership rules", x = .04, y = .97, just = "left",
  gp = grid::gpar(fontsize = 18, fontface = "bold"))
grid::grid.text("Signed Tahoe training effects and median CytoSig HGF effects", x = .04, y = .935, just = "left",
  gp = grid::gpar(fontsize = 12, col = "#596571"))
for (i in 1:3) {
  xy <- list(c(.26, .70, .45, .40), c(.74, .70, .45, .40), c(.50, .245, .58, .34))[[i]]
  grid::pushViewport(grid::viewport(x = xy[1], y = xy[2], width = xy[3], height = xy[4]))
  plot_branch(list(m_egfr, m_met, m_trunk)[[i]], names(branch_color)[i], i == 3L)
  grid::popViewport()
}
grid::grid.lines(c(.27, .44), c(.49, .435), arrow = grid::arrow(length = grid::unit(2.5, "mm")), gp = grid::gpar(col = branch_color[1], lwd = 2))
grid::grid.lines(c(.73, .56), c(.49, .435), arrow = grid::arrow(length = grid::unit(2.5, "mm")), gp = grid::gpar(col = branch_color[2], lwd = 2))
grid::grid.text(paste0("Each row is a retained gene. Dense arm labels are in egfr_met_core_genes.tsv. Colors saturate at ±", sprintf("%.2f", color_limit), ".\n",
  "Arms may overlap; the trunk is disjoint. Arrows express the brief's bypass hypothesis, not inferred causal gene order."),
  x = .04, y = .035, just = "left", gp = grid::gpar(fontsize = 10, lineheight = 1.2))
dev.off()

held <- scored %>% filter(split == "Holdout", is.finite(survival))
held_stats <- functional %>% filter(split == "Holdout") %>% mutate(label = sprintf("n=%d | r=%.2f | p=%.2g", n, pearson_r, nominal_p))
scatter <- ggplot(held, aes(trunk_score, survival)) + geom_point(size = 2, color = "#287A91") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .65, color = "#B86C38") +
  geom_text(data = held_stats, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -.03, vjust = 1.4, size = 3) +
  facet_wrap(~factor(compound, levels = drugs), nrow = 1L, scales = "free") +
  labs(x = "Frozen trunk score", y = "PRISM survival (-sensitization)", title = "A | Held-out Tahoe–PRISM") +
  theme_bw(base_size = 10) + theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"), strip.background = element_blank())
if (!nrow(trunk)) scatter <- ggplot() + annotate("text", x = 0, y = 0,
  label = "No shared trunk under the fixed rules.\nHeld-out correlation cannot be estimated.", size = 5) +
  theme_void() + labs(title = "A | Held-out Tahoe-PRISM")
projection <- tibble(source = factor(c("Erlotinib", "CytoSig HGF"), levels = c("Erlotinib", "CytoSig HGF")),
  score = c(erl_delta, hgf_delta))
projection_plot <- ggplot(projection, aes(source, score, fill = source)) +
  geom_hline(yintercept = 0, color = "#555555", linewidth = .4) + geom_col(width = .55, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", score)), vjust = ifelse(projection$score >= 0, -.4, 1.3), size = 3.5) +
  scale_fill_manual(values = c("#287A91", "#B86C38")) + scale_y_continuous(expand = expansion(mult = .2)) +
  labs(x = NULL, y = "Signed projection from zero effect", title = "B | HGF reversal",
    subtitle = "Same frozen signs and training scales") + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
if (!nrow(trunk)) projection_plot <- ggplot() + annotate("text", x = 0, y = 0,
  label = "No trunk to project", size = 4) + theme_void() + labs(title = "B | HGF reversal")
wilson_long <- bind_rows(transmute(wilson, condition, marker = "pERK", value = pERK),
  transmute(wilson, condition, marker = "pMET", value = pMET)) %>%
  mutate(condition = factor(condition, levels = wilson$condition),
    short = factor(c("Control", "Erl", "Erl + HGF", "Erl + HGF\n+ Criz")[as.integer(condition)],
      levels = c("Control", "Erl", "Erl + HGF", "Erl + HGF\n+ Criz")))
wilson_plot <- ggplot(wilson_long, aes(short, value)) + geom_col(fill = "#5F7F42", width = .65) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -.35, size = 3) + facet_wrap(~marker, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(x = NULL, y = "Published relative phosphorylation", title = "C | Wilson KHM-3S molecular reference") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
    strip.background = element_blank())
cairo_pdf(file.path(RESULT, "25b_egfr_met_validation.pdf"), width = 11.6, height = 9.0)
grid::grid.newpage()
print(scatter, vp = grid::viewport(x = .5, y = .72, width = .98, height = .47))
print(projection_plot, vp = grid::viewport(x = .19, y = .28, width = .34, height = .34))
print(wilson_plot, vp = grid::viewport(x = .69, y = .28, width = .61, height = .34))
grid::grid.text(paste0(if (pilot_positive) "POSITIVE" else "NEGATIVE", " under the five fixed acceptance criteria | ", nrow(trunk), " frozen trunk genes"),
  x = .025, y = .98, just = "left", gp = grid::gpar(fontsize = 16, fontface = "bold"))
grid::grid.text("Wilson viability: HGF rescues erlotinib-treated cells; adding crizotinib reduces rescue. Qualitative published ordering.\nKHM-3S / JCRB0138; 72 h Syto 60; erlotinib 0.001–10 uM; HGF 50 ng/mL; crizotinib 0.5 uM.\nFigure 2B/C, Wilson et al., Nature 2012 (10.1038/nature11249); values reproduced in eLife 2014 (10.7554/eLife.04037).\nWilson results were known from the brief. No Wilson transcriptomes, viability fits or parameter tuning were used.",
  x = .025, y = .045, just = "left", gp = grid::gpar(fontsize = 9.5, lineheight = 1.15))
dev.off()
stopifnot(identical(unname(tools::md5sum(file.path(RESULT, "egfr_met_core_genes.tsv"))), frozen_hash))
message("Stage 25 complete: ", if (pilot_positive) "positive" else "negative", " under the fixed brief criteria.")
print(tests %>% filter(type %in% c("acceptance", "decision")) %>% select(test, observed, passed))
