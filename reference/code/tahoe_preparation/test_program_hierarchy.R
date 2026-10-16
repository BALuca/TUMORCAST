#!/usr/bin/env Rscript

## Exploratory hierarchy sensitivity analysis authorized after input validation.
## Thresholds below are diagnostic assumptions, not validated activity calls.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
if (length(commandArgs(trailingOnly = TRUE)))
  fail("This input-validation stage accepts no command-line arguments; use TAHOE_V2_ROOT.")

read_required <- function(path, columns, key) {
  if (!file.exists(path)) fail("Missing required input: %s", path)
  x <- as_tibble(fread(path))
  missing <- setdiff(union(columns, key), names(x))
  if (length(missing)) fail("%s lacks columns: %s", path,
    paste(missing, collapse = ", "))
  if (!nrow(x)) fail("Empty required input: %s", path)
  if (anyNA(x[, key, drop = FALSE]) || anyDuplicated(x[, key, drop = FALSE]))
    fail("Missing or duplicate key (%s): %s", paste(key, collapse = ", "), path)
  x
}
stability <- read_required(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"), c("outer_fold", "program", "stable"),
  c("outer_fold", "program", "repeat"))
if (any((stability %>% group_by(outer_fold, program) %>%
    summarise(n_status = n_distinct(stable), .groups = "drop"))$n_status != 1L))
  fail("Repeated stability records disagree on retained program membership.")
programs <- stability %>% filter(outer_fold == OUTER_FOLD, stable) %>%
  distinct(program) %>% pull(program)
if (length(programs) != EXPECTED_RESPONSE_PROGRAMS)
  fail("Expected nine retained programs; found %d.", length(programs))
scores <- read_required(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"), c("observation_id", "depmap_id", "cell_line",
    "perturbation_id", "drug", "dose", "dose_panel", "role", "outer_fold",
    "program", "score"), c("outer_fold", "observation_id", "program")) %>%
  filter(outer_fold == OUTER_FOLD, program %in% programs)
if (!nrow(scores) || any(!is.finite(scores$score)))
  fail("Frozen retained program scores are empty or nonfinite.")
if (any((scores %>% count(observation_id))$n != length(programs)))
  fail("An observation lacks one or more retained program scores.")
metadata <- scores %>% distinct(observation_id, depmap_id, cell_line,
  perturbation_id, drug, dose, dose_panel, role)
if (anyDuplicated(metadata$observation_id))
  fail("Program rows disagree on observation metadata.")
if (any(!is.finite(metadata$dose)) || any(metadata$dose <= 0) ||
    anyNA(metadata$dose_panel) || any(!metadata$dose_panel %in% 1:3))
  fail("Missing or invalid within-drug dose order.")
doses <- metadata %>% distinct(drug, dose, dose_panel) %>% arrange(drug, dose_panel)
dose_check <- doses %>% group_by(drug) %>% summarise(
  valid = n() == 3L && n_distinct(dose_panel) == 3L && all(diff(dose) > 0),
  .groups = "drop")
if (!all(dose_check$valid)) fail("Within-drug dose order is not three increasing doses.")
if (anyDuplicated(metadata[, c("depmap_id", "drug", "dose_panel")]))
  fail("Duplicate drug--cell-line--dose observation.")
folds <- read_required(file.path(stage_result("01_inputs"),
  "fold_assignments.tsv"), c("depmap_id", "outer_fold"), "depmap_id")
fold_position <- match(metadata$depmap_id, folds$depmap_id)
if (anyNA(fold_position)) fail("A scored cell line lacks an existing split assignment.")
expected_role <- ifelse(folds$outer_fold[fold_position] == OUTER_FOLD,
  "held_out", "training")
if (anyNA(metadata$role) || any(metadata$role != expected_role))
  fail("Score roles disagree with the existing development/holdout split.")
panel <- read_required(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"), c("perturbation_id", "activity_supported"),
  "perturbation_id")
panel_position <- match(metadata$perturbation_id, panel$perturbation_id)
if (anyNA(panel_position) || anyNA(panel$activity_supported[panel_position]))
  fail("Invalid join to the existing molecular-activity filter.")
metadata$activity_supported <- panel$activity_supported[panel_position]
outcomes <- read_required(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"), c("depmap_id", "perturbation_id", "tahoe_drug",
    "tahoe_dose_uM", "prism_sensitization", "outer_fold", "chain_observed"),
  c("depmap_id", "perturbation_id")) %>% filter(chain_observed == 1L) %>%
  mutate(observation_id = paste(depmap_id, perturbation_id, sep = "__"))
position <- match(outcomes$observation_id, metadata$observation_id)
if (!nrow(outcomes) || anyNA(position))
  fail("Existing PRISM matches are empty or lack frozen program scores.")
if (any(!is.finite(outcomes$prism_sensitization)) ||
    any(outcomes$tahoe_drug != metadata$drug[position]) ||
    any(outcomes$tahoe_dose_uM != metadata$dose[position]) ||
    any(outcomes$outer_fold != folds$outer_fold[fold_position[position]]))
  fail("Existing PRISM matches disagree with score identifiers, split, or outcome.")
outcomes$activity_supported <- metadata$activity_supported[position]
outcomes$role <- metadata$role[position]

## Print exact input census before exploratory recruitment calculations.
for (cohort in c("all_frozen", "activity_supported")) {
  m <- metadata
  f <- outcomes
  if (cohort == "activity_supported") {
    m <- filter(m, activity_supported)
    f <- filter(f, activity_supported)
  }
  for (split in c("all", "training", "held_out")) {
    current <- if (split == "all") m else filter(m, role == split)
    functional <- if (split == "all") f else filter(f, role == split)
    complete <- current %>% count(drug, depmap_id) %>% filter(n == 3L)
    cat(sprintf(paste0("%s / %s: programs=%d; drugs=%d; cell_lines=%d; ",
      "complete_three_dose_series=%d; matched_PRISM_observations=%d\n"),
      cohort, split, length(programs), n_distinct(current$drug),
      n_distinct(current$depmap_id), nrow(complete), nrow(functional)))
  }
}

suppressPackageStartupMessages(library(ggplot2))
RESULT <- stage_result("18_program_hierarchy")
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
wide <- dcast(as.data.table(scores), observation_id ~ program, value.var = "score")
D <- as.matrix(wide[match(metadata$observation_id, wide$observation_id), ..programs])
if (any(!is.finite(D))) fail("Invalid wide score matrix.")
train <- outcomes$role == "training"
X <- D[position, , drop = FALSE]
frozen <- readRDS(file.path(stage_artifact("05_program_function_bridge"),
  "function_decoder_model.rds"))
if (!identical(frozen$feature_names, programs)) fail("Frozen decoder features differ.")
base_fit <- fit_ridge(X[train, ], outcomes$prism_sensitization[train], frozen$lambda)
base_prediction <- predict_ridge(base_fit, X)
if (max(abs(base_prediction - predict_ridge(frozen, X))) > 1e-8)
  fail("Refitted current model does not reproduce the frozen decoder.")
saved <- fread(file.path(stage_artifact("05_program_function_bridge"),
  "measured_function_predictions.tsv.gz"))[evaluation_included == TRUE]
if (anyDuplicated(saved$observation_id)) fail("Duplicate frozen evaluation key.")
held_position <- match(outcomes$observation_id[!train], saved$observation_id)
if (anyNA(held_position) || max(abs(base_prediction[!train] -
    saved$measured_prediction[held_position])) > 1e-8)
  fail("Current held predictions do not reproduce Stage 05.")

## Signed thresholds in units of each program's training SD. Zero is relative
## to the centered training response, not a DMSO activity threshold. Positive
## and negative excursions are evaluated separately; abs(score) is never used.
score_sd <- apply(D[metadata$role == "training", ], 2, sd)
if (any(!is.finite(score_sd) | score_sd <= 0)) fail("Invalid training score SD.")
series <- metadata %>% group_by(drug, depmap_id, role) %>% summarise(
  n_doses = n(), n_supported = sum(activity_supported), .groups = "drop") %>%
  filter(n_doses == 3L)
series$id <- paste(series$drug, series$depmap_id, sep = "__")
row_series <- paste(metadata$drug, metadata$depmap_id, sep = "__")
scenarios <- expand.grid(sign_direction = c("positive", "negative"),
  sd_threshold = c(1, 0.5, 0), series_filter = c("all_doses_supported",
    "any_dose_supported", "all_complete"), stringsAsFactors = FALSE)
scenarios$scenario <- with(scenarios, paste(sign_direction, sd_threshold,
  series_filter, sep = "__"))
order_parts <- list(); test_parts <- list()
add_test <- function(x) test_parts[[length(test_parts) + 1L]] <<- as.data.table(x)

## Leave-one-group-out direction must survive removal of EVERY drug and line.
## Undefined comparisons are explicit diagnostic statuses, not model fallbacks.
order_count <- function(first, meta, a, b) {
  aa <- first[, a]; bb <- first[, b]
  forward <- aa < bb & bb < 4L; reverse <- bb < aa & aa < 4L
  n_f <- sum(forward); n_r <- sum(reverse); informative <- forward | reverse
  deletion_min <- function(group) {
    fs <- tapply(as.integer(forward), group, sum)
    rs <- tapply(as.integer(reverse), group, sum)
    den <- n_f + n_r - fs - rs
    if (any(den == 0)) return(NA_real_)
    min((n_f - fs) / den)
  }
  marginal_fraction <- function(group) {
    cells <- split(seq_along(aa), group)
    expected <- vapply(cells, function(ix) {
      ta <- tabulate(aa[ix], nbins = 4L); tb <- tabulate(bb[ix], nbins = 4L)
      c(forward = ta[1] * (tb[2] + tb[3]) + ta[2] * tb[3],
        reverse = tb[1] * (ta[2] + ta[3]) + tb[2] * ta[3]) / length(ix)
    }, numeric(2))
    sums <- rowSums(expected)
    if (sum(sums) == 0) return(NA_real_)
    sums[1] / sum(sums)
  }
  data.frame(A = a, B = b, n_series = nrow(first),
    B_only_fraction_among_B_recruited = if (any(bb < 4L))
      sum(bb < 4L & aa == 4L) / sum(bb < 4L) else NA_real_,
    independent_timing_fraction_within_line = marginal_fraction(meta$depmap_id),
    independent_timing_fraction_within_drug = marginal_fraction(meta$drug),
    n_A_before_B = n_f, n_B_before_A = n_r,
    n_same_dose = sum(aa == bb & aa < 4L),
    n_A_only = sum(aa < 4L & bb == 4L),
    n_B_only = sum(bb < 4L & aa == 4L), n_neither = sum(aa == 4L & bb == 4L),
    n_informative = sum(informative),
    forward_fraction = if (sum(informative)) n_f / sum(informative) else NA_real_,
    n_drugs = n_distinct(meta$drug[informative]),
    n_lines = n_distinct(meta$depmap_id[informative]),
    min_leave_drug_out_fraction = deletion_min(meta$drug),
    min_leave_line_out_fraction = deletion_min(meta$depmap_id))
}
metric_row <- function(rows, pred) {
  tab <- as.data.table(outcomes[rows, ])
  tab[, `:=`(drug = tahoe_drug, observed_function = prism_sensitization,
    predicted_function = pred[rows])]
  res <- function_residual_table(tab, "within_perturbation")
  r <- safe_cor(res$observed_residual, res$predicted_residual)
  if (!is.finite(r)) fail("Undefined model correlation.")
  list(pearson = r, rmse = sqrt(mean((res$observed_residual -
    res$predicted_residual)^2)), n_observations = nrow(res),
    n_drugs = uniqueN(res$drug), n_lines = uniqueN(res$depmap_id))
}
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]; message("Diagnostic scenario: ", sc$scenario)
  orientation <- if (sc$sign_direction == "positive") 1 else -1
  active <- sweep(orientation * D, 2, sc$sd_threshold * score_sd, ">")
  eligible <- switch(sc$series_filter,
    all_doses_supported = series$n_supported == 3L,
    any_dose_supported = series$n_supported >= 1L,
    all_complete = rep(TRUE, nrow(series)))
  ss <- series[eligible, ]
  dose_positions <- lapply(1:3, function(d) match(ss$id,
    row_series[metadata$dose_panel == d]))
  dose_active <- lapply(1:3, function(d) active[which(metadata$dose_panel == d)[
    dose_positions[[d]]], , drop = FALSE])
  if (anyNA(unlist(dose_positions))) fail("Incomplete series entered recruitment analysis.")
  first <- matrix(4L, nrow(ss), length(programs), dimnames = list(NULL, programs))
  for (d in 3:1) first[dose_active[[d]]] <- d
  for (role in c("training", "held_out")) for (p in programs) {
    use <- ss$role == role
    persistent <- !(dose_active[[1]][use, p] & !dose_active[[2]][use, p]) &
      !(dose_active[[2]][use, p] & !dose_active[[3]][use, p])
    add_test(data.frame(record_type = "recruitment_diagnostic", scenario = sc$scenario,
      role = role, program = p, signed_threshold = orientation * sc$sd_threshold * score_sd[p],
      training_score_sd = score_sd[p], n_series = sum(use),
      fraction_dose1 = mean(first[use, p] == 1), fraction_dose2 = mean(first[use, p] == 2),
      fraction_dose3 = mean(first[use, p] == 3), fraction_never = mean(first[use, p] == 4),
      fraction_loses_activity = mean(!persistent)))
  }
  pairs <- expand.grid(A = programs, B = programs, stringsAsFactors = FALSE) %>% filter(A != B)
  counts <- bind_rows(lapply(c("training", "held_out"), function(role) {
    use <- ss$role == role
    bind_rows(lapply(seq_len(nrow(pairs)), function(k)
      order_count(first[use, , drop = FALSE], ss[use, ], pairs$A[k], pairs$B[k]))) %>%
      mutate(role = role)
  }))
  development <- counts %>% filter(role == "training") %>% mutate(
    candidate = n_informative >= 20 & forward_fraction >= 0.70)
  validation <- counts %>% filter(role == "held_out")
  candidates <- development %>% filter(candidate)
  vpos <- match(paste(candidates$A, candidates$B), paste(validation$A, validation$B))
  v <- validation[vpos, ]
  replicated <- !is.na(v$forward_fraction) & v$forward_fraction > 0.5 &
    v$n_drugs >= 2 & v$n_lines >= 2 & candidates$n_drugs >= 2 & candidates$n_lines >= 2 &
    !is.na(v$min_leave_drug_out_fraction) & v$min_leave_drug_out_fraction > 0.5 &
    !is.na(v$min_leave_line_out_fraction) & v$min_leave_line_out_fraction > 0.5 &
    candidates$min_leave_drug_out_fraction > 0.5 &
    candidates$min_leave_line_out_fraction > 0.5
  replicated[is.na(replicated)] <- FALSE
  orders <- candidates[replicated, ]
  key <- paste(counts$A, counts$B)
  counts$development_candidate <- key %in% paste(candidates$A, candidates$B)
  counts$replicated_broad_direction <- key %in% paste(orders$A, orders$B)
  counts$scenario <- sc$scenario
  counts$sign_direction <- sc$sign_direction
  counts$sd_threshold <- sc$sd_threshold
  counts$series_filter <- sc$series_filter
  counts$excess_forward_over_line_marginals <- counts$forward_fraction -
    counts$independent_timing_fraction_within_line
  counts$excess_forward_over_drug_marginals <- counts$forward_fraction -
    counts$independent_timing_fraction_within_drug
  ## These cutoffs summarize nesting sensitivity only; candidate selection
  ## continues to use the requested 20-series/0.70 directional screen.
  counts$B_only_at_most_10pct <- counts$B_only_fraction_among_B_recruited <= 0.10
  counts$B_only_at_most_20pct <- counts$B_only_fraction_among_B_recruited <= 0.20
  order_parts[[i]] <- as.data.table(counts)
  interactions <- matrix(numeric(), nrow(X), 0L)
  n_effect_replicated <- 0L
  for (k in seq_len(nrow(orders))) {
    a <- orders$A[k]; b <- orders$B[k]
    aa <- active[position, a]; bb <- active[position, b]
    interactions <- cbind(interactions, as.numeric(aa & bb))
    colnames(interactions)[ncol(interactions)] <- paste0(a, "_and_", b)
    effects <- list()
    for (role in c("training", "held_out")) {
      use <- outcomes$role == role & outcomes$activity_supported & aa
      tab <- as.data.table(outcomes[use, ])
      tab[, state := as.numeric(bb[use])]
      ## OLS with drug-dose fixed intercepts, implemented by exact within-group
      ## demeaning in the A-active subset. Outcome is existing sensitization;
      ## negate coefficient to express A+B minus A-only log2 viability.
      tab[, `:=`(x = state - mean(state),
        y = prism_sensitization - mean(prism_sensitization)), by = perturbation_id]
      informative <- tab[x != 0]
      estimable <- nrow(informative) > 0L
      beta <- if (estimable) sum(tab$x * tab$y) / sum(tab$x^2) else NA_real_
      effect <- -beta
      status <- if (estimable) "estimated" else "not_estimable_no_within_dose_state_contrast"
      if (!estimable) message("  No conditional contrast: ", a, " -> ", b, " / ", role)
      cluster_se <- NA_real_
      if (estimable && uniqueN(informative$tahoe_drug) > 1L) {
        tab[, score_residual := x * (y - beta * x)]
        cluster_scores <- tab[, .(u = sum(score_residual)), by = tahoe_drug]
        g <- uniqueN(informative$tahoe_drug)
        cluster_se <- sqrt(g / (g - 1) * sum(cluster_scores$u^2)) / sum(tab$x^2)
      }
      effects[[role]] <- effect
      add_test(data.frame(record_type = "conditional_effect", scenario = sc$scenario,
        A = a, B = b, role = role, status = status,
        viability_difference_AB_minus_A = effect,
        drug_cluster_standard_error = cluster_se,
        n_observations = nrow(tab), n_drugs = uniqueN(tab$tahoe_drug),
        n_lines = uniqueN(tab$depmap_id), n_A_only = sum(tab$state == 0),
        n_A_plus_B = sum(tab$state == 1), n_contrast_observations = nrow(informative),
        n_contrast_drugs = uniqueN(informative$tahoe_drug),
        n_contrast_lines = uniqueN(informative$depmap_id),
        n_A_only_drugs = uniqueN(tab[state == 0, tahoe_drug]),
        n_AB_drugs = uniqueN(tab[state == 1, tahoe_drug]),
        n_A_only_lines = uniqueN(tab[state == 0, depmap_id]),
        n_AB_lines = uniqueN(tab[state == 1, depmap_id])))
    }
    n_effect_replicated <- n_effect_replicated +
      as.integer(isTRUE(effects$training * effects$held_out > 0))
  }
  ext_prediction <- base_prediction
  if (ncol(interactions)) {
    if (any(apply(interactions[train, , drop = FALSE], 2, sd) == 0))
      fail("Constant conditional term in development data.")
    design <- cbind(X, interactions)
    fit <- fit_ridge(design[train, ], outcomes$prism_sensitization[train], frozen$lambda)
    ext_prediction <- predict_ridge(fit, design)
    if (any(!is.finite(ext_prediction))) fail("Extended model returned nonfinite predictions.")
  }
  for (cohort in c("activity_supported", "all_matched")) {
    use <- !train & (if (cohort == "activity_supported") outcomes$activity_supported else TRUE)
    bm <- metric_row(use, base_prediction); em <- metric_row(use, ext_prediction)
    add_test(data.frame(record_type = "model_comparison", scenario = sc$scenario,
      cohort = cohort, role = "held_out", current_pearson = bm$pearson,
      extended_pearson = em$pearson, delta_pearson = em$pearson - bm$pearson,
      current_rmse = bm$rmse, extended_rmse = em$rmse,
      n_observations = bm$n_observations, n_drugs = bm$n_drugs, n_lines = bm$n_lines,
      n_development_candidates = nrow(candidates), n_replicated_orders = nrow(orders),
      n_held_orders_B_only_at_most_10pct = sum(counts$role == "held_out" &
        counts$replicated_broad_direction & counts$B_only_at_most_10pct, na.rm = TRUE),
      n_same_direction_effects = n_effect_replicated, frozen_lambda = frozen$lambda,
      status = if (ncol(interactions)) "fitted" else "no_replicated_terms; models_identical"))
    ## Cluster jackknife diagnoses whether the apparent improvement depends on
    ## one held-out line or drug. Models remain fixed; no refitting or selection.
    if (cohort == "activity_supported") for (group in c("depmap_id", "tahoe_drug")) {
      groups <- unique(outcomes[[group]][use])
      delta <- vapply(groups, function(g) {
        keep <- use & outcomes[[group]] != g
        metric_row(keep, ext_prediction)$pearson - metric_row(keep, base_prediction)$pearson
      }, numeric(1))
      add_test(data.frame(record_type = "performance_group_deletion",
        scenario = sc$scenario, cohort = cohort, omitted_group = group,
        omitted_id = groups, delta_pearson = delta))
      add_test(data.frame(record_type = "performance_deletion_diagnostic",
        scenario = sc$scenario, cohort = cohort, omitted_group = group,
        min_delta_pearson = min(delta), max_delta_pearson = max(delta),
        fraction_deletions_improving = mean(delta > 0)))
    }
  }
}
order_table <- rbindlist(order_parts, fill = TRUE)
prism_table <- rbindlist(test_parts, fill = TRUE)
prism_table[record_type == "model_comparison", conclusion := fifelse(
  n_replicated_orders == 0, "no_directional_order_support",
  fifelse(delta_pearson <= 0, "directional_order_only_no_prediction_gain",
    "exploratory_gain_requires_nesting_and_group_robustness_review"))]
for (scenario_id in scenarios$scenario) {
  deletion <- prism_table[scenario == scenario_id &
    record_type == "performance_deletion_diagnostic"]
  broad_gain <- nrow(deletion) == 2L && all(deletion$min_delta_pearson > 0)
  prism_table[scenario == scenario_id & record_type == "model_comparison" &
    cohort == "activity_supported", robust_prediction_gain := delta_pearson > 0 & broad_gain]
}
prism_table[, interpretation := "exploratory sensitivity; activity thresholds not previously validated"]
prism_table[, outcome_definition := "model: within-drug-dose sensitization; effect: A+B minus A-only log2 viability with drug-dose fixed intercepts"]
prism_table[, selection_note := "directions screened in development; holdout molecular replication selects terms; holdout is therefore not untouched validation"]
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
fwrite(order_table, file.path(RESULT, "program_order_tests.tsv"), sep = "\t")
fwrite(prism_table, file.path(RESULT, "program_hierarchy_prism_tests.tsv"), sep = "\t")

## One panel per PDF. Matrix shows the specified 1-SD positive reference
## sensitivity, not whichever scenario produces the highest held-out gain.
source(file.path(WORK_ROOT, "src", "R", "lib_heatmap.R"))
reference <- "positive__1__any_dose_supported"
mm <- matrix(NA_real_, length(programs), length(programs), dimnames = list(programs, programs))
labels <- matrix("", length(programs), length(programs), dimnames = dimnames(mm))
cc <- order_table[scenario == reference & role == "held_out"]
for (k in seq_len(nrow(cc))) {
  mm[cc$A[k], cc$B[k]] <- cc$forward_fraction[k]
  labels[cc$A[k], cc$B[k]] <- sprintf("%.2f\nn=%d%s", cc$forward_fraction[k],
    cc$n_informative[k], ifelse(cc$replicated_broad_direction[k], "*", ""))
}
block_order <- hclust(as.dist(1 - cor(D[metadata$role == "training", ])))$order
mm <- mm[block_order, block_order]; labels <- labels[block_order, block_order]
ht <- draw_heatmap(mm, scale_rows = FALSE, hm_name = "Forward fraction",
  value_range = c(0, 0.5, 1), value_colors = c("#2166AC", "white", "#B2182B"),
  show_row_names = TRUE, show_column_names = TRUE,
  title = "Held-out order: positive > 1 SD; any dose supported\n* broad replicated direction; n = informative series",
  row_title = "A (upstream)", column_title_gp = grid::gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill)
    grid::grid.text(labels[i, j], x, y, gp = grid::gpar(fontsize = 6.5)))
grDevices::cairo_pdf(file.path(RESULT, "program_order_matrix.pdf"), width = 6.2, height = 5.5)
ComplexHeatmap::draw(ht)
grDevices::dev.off()
performance <- prism_table[record_type == "model_comparison" & cohort == "activity_supported"]
performance[, scenario := factor(scenario, levels = rev(scenarios$scenario))]
p <- ggplot(performance, aes(y = scenario)) +
  geom_segment(aes(x = current_pearson, xend = extended_pearson, yend = scenario), color = "grey60") +
  geom_point(aes(x = current_pearson, color = "Current"), size = 1.8) +
  geom_point(aes(x = extended_pearson, color = "Extended"), size = 1.8, shape = 17) +
  scale_y_discrete(labels = function(x) {
    x <- gsub("positive__", "Positive, ", x, fixed = TRUE)
    x <- gsub("negative__", "Negative, ", x, fixed = TRUE)
    x <- gsub("__all_doses_supported", " SD | all doses pass", x, fixed = TRUE)
    x <- gsub("__any_dose_supported", " SD | any dose passes", x, fixed = TRUE)
    gsub("__all_complete", " SD | all complete series", x, fixed = TRUE)
  }) +
  scale_color_manual(values = c(Current = "#333333", Extended = "#0072B2"), name = NULL) +
  labs(x = "Held-out within-drug-dose Pearson r", y = NULL,
    title = "Exploratory signed-threshold and series-filter sensitivity",
    subtitle = "Same frozen ridge penalty; same 646 activity-supported observations") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(),
    legend.position = "bottom", plot.title = element_text(size = 10))
ggsave(file.path(RESULT, "program_hierarchy_prism.pdf"), p,
  width = 7.1, height = 5.6, device = cairo_pdf)
print(performance[, .(scenario, n_development_candidates, n_replicated_orders,
  n_same_direction_effects, current_pearson, extended_pearson, delta_pearson,
  current_rmse, extended_rmse)])
message("Exploratory diagnostics complete; no confirmatory hierarchy claim is made.")
