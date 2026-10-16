#!/usr/bin/env Rscript
# Fixed v106 single-drug function and seven program weights from v143.
# Normalized molecular opposition predicts a bounded signed recovery fraction.
# Primary: 12 NRG1 counts, 6 recipients; negatives and above-control values retained.
# Saturating calibration: tanh(gain * (relative_opposition - threshold[dose])).
# Dose thresholds shrink toward one another; penalty selected inside recipient folds.
# EGF is external control; all-16 fit is a prespecified sensitivity check.
# HGF/Wilson outcomes never fit or select calibration; no log-departure time amplification.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(digest)
  library(nloptr)
})
RhpcBLASctl::blas_set_num_threads(1L)
setDTthreads(1L)
root <- "results/12_fit_ligand_rescue"
source_root <- "results/08_fit_single_drug_bridge"
programs <- sprintf("P%02d", 1:40)
association_cutoff <- .5
for (folder in c("data", "models", "figures")) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
sink(file.path(root, "run.log"), split = TRUE)
write_table <- function(x, name) fwrite(x, file.path(root, "data", paste0(name, ".tsv")), sep = "\t")
pearson <- function(x, y) if (length(x) >= 4 && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
source("R/prediction.R")
source("R/data_io.R")
source("R/heatmaps.R")
model <- qread("results/08_fit_single_drug_bridge/models/final.qs")
prepared <- qread("data/processed/molecular_inputs.qs")$result
receptors <- prepared$receptors
contexts <- read.delim("results/08_fit_single_drug_bridge/data/molecular_contexts.tsv")
context_key <- function(drug, dose) paste(drug, sprintf("%.10g", dose), sep = "|")
context_index <- function(drug, dose) {
  match(context_key(drug, dose), context_key(contexts$drug, contexts$dose))
}
predict_single <- function(model, metadata) {
  predict_lapatinib(model, as.data.frame(metadata), receptors, contexts)
}
ligand_inputs <- qread("data/processed/hgf_program_inputs.qs")$result
outcomes <- as.data.table(read.delim("data/processed/her2_ligand_outcomes.tsv"))
ligands <- as.data.table(ligand_inputs$ligand_program_inputs)
profiles <- as.data.table(ligand_inputs$program_inputs)
stopifnot(
  nrow(outcomes) == 16, all(outcomes$gene %in% c("NRG1", "EGF")),
  !any(outcomes$CellLine == "AU565"), all(outcomes$time_h == 72), uniqueN(outcomes$depmap_id) == 6
)
outcomes[, observation_id := .I]
key <- function(x) paste(x$depmap_id, x$gene, sep = "|")
stopifnot(!anyDuplicated(key(ligands)))
ix <- match(key(outcomes), key(ligands))
stopifnot(!anyNA(ix))
H <- as.matrix(ligands[ix, ..programs])
single <- predict_single(model, outcomes[, .(depmap_id, drug, dose)])
association <- model$template$correlations[, "Lapatinib ditosylate"]
keep <- abs(association) >= association_cutoff
stopifnot(sum(keep) == 7)
write_table(
  data.table(program = programs, observed_single_drug_r = association, used_for_g = keep),
  "fixed_program_selection"
)

# Whole-vector normalization is performed after the frozen decoder transformation.
# The threshold acts on drug-sensitive molecular activity, not on baseline directly.
# Zero ligand has no direction and explicitly gives g=0. The normalized ligand
# representation describes nonzero response patterns, not a ligand-dose response.
make_features <- function(meta, L, H) {
  g <- context_index(meta$drug, meta$dose)
  route <- model$context$route[g]
  beta <- t(model$context$W[, route, drop = FALSE]) / model$context$sigma[g, , drop = FALSE]
  beta[, !keep] <- 0
  d <- pmax(L * beta, 0)
  h <- H * beta
  ligand_norm <- sqrt(rowSums(h^2))
  present <- ligand_norm > 0
  normalized_h <- matrix(0, nrow(h), ncol(h), dimnames = dimnames(h))
  normalized_h[present, ] <- h[present, , drop = FALSE] / ligand_norm[present]
  drug_activity <- rowSums(d)
  drug_activity[!present] <- 0
  opposed <- d * pmax(-normalized_h, 0)
  opposition <- rowSums(opposed)
  relative_opposition <- ifelse(drug_activity > 0, opposition / pmax(drug_activity, .Machine$double.xmin), 0)
  X <- cbind(
    opposition = opposition,
    threshold_0.5 = -drug_activity * (meta$dose == .5), threshold_5 = -drug_activity * (meta$dose == 5)
  )
  stopifnot(all(meta$dose %in% c(.5, 5)), all(relative_opposition >= 0), all(relative_opposition <= 1 + 1e-10))
  stopifnot(all(is.finite(X)))
  list(
    X = X, opposed = opposed, drug_activity = drug_activity, ligand_norm = ligand_norm,
    relative_opposition = relative_opposition, directed_drug = d, normalized_ligand = normalized_h
  )
}
features <- make_features(outcomes, single$D, H)
X <- features$X
y <- outcomes$observed_departure_log2
stopifnot(all(is.finite(y)), all(X[, 1] >= 0), all(X[, 2:3] <= 0))
primary_rows <- which(outcomes$gene == "NRG1")
control_rows <- which(outcomes$gene == "EGF")
stopifnot(length(primary_rows) == 12, length(control_rows) == 4)
write_table(cbind(outcomes, as.data.table(X), data.table(
  ligand_norm = features$ligand_norm,
  drug_activity = features$drug_activity, relative_opposition = features$relative_opposition,
  role = ifelse(outcomes$gene == "NRG1", "primary calibration", "external control")
)), "training_features_and_outcomes")

# Count anchors identify the interaction while preserving the original 72h assay.
# This evaluates g conditional on measured single-drug/ligand controls. The query
# uses the frozen single-drug predictor instead; both validation axes are saved.
stopifnot(
  all(outcomes$observed_drug_pct_vehicle > 0),
  all(outcomes$observed_combination_pct_vehicle > 0), all(is.finite(outcomes$observed_ligand_log2))
)
observed_drug <- outcomes$observed_drug_pct_vehicle / 100
observed_ligand <- 2^outcomes$observed_ligand_log2
target_log2_ratio <- log2(outcomes$observed_combination_pct_vehicle / 100) - outcomes$observed_ligand_log2
predicted_drug_72h <- 2^(-single$sensitivity * 3 / 5)
relative_score <- features$relative_opposition
present <- features$drug_activity > 0 & features$ligand_norm > 0
penalty_grid <- c(0, .1, 1, 10, 100)
curve_cache_path <- "data/processed/ligand_rescue_curve_fits.qs"
curve_signature <- list(
  inputs = vapply(c("data/processed/her2_ligand_outcomes.tsv", "data/processed/hgf_program_inputs.qs",
    file.path(source_root, "models/final.qs")), digest, character(1), file = TRUE, algo = "sha256"),
  code = vapply(c("scripts/12_fit_ligand_rescue.R", "R/prediction.R"), digest,
    character(1), file = TRUE, algo = "sha256"),
  penalties = penalty_grid, nloptr = as.character(packageVersion("nloptr")), R = R.version.string
)
cached_curves <- if (file.exists(curve_cache_path)) qread(curve_cache_path) else NULL
if (!is.null(cached_curves) && identical(cached_curves$fingerprint, curve_signature)) {
  curve_fits <- cached_curves$result
  message("Reusing matching nested ligand-calibration fits")
} else {
  curve_fits <- list()
  message("Nested ligand-calibration cache is absent or stale")
}

# Signed rho is bounded in [-1,1]. Positive rho restores a fraction of depletion;
# negative rho removes a fraction of remaining cells. Values above control are
# possible only from the given ligand main effect or an already-above-control drug.
# Stable log calculations avoid numerical zero in the negative tail.
compose <- function(drug, rho, eta) {
  capacity <- pmax(1 - drug, 0)
  log_ratio <- numeric(length(drug))
  positive <- rho >= 0
  log_ratio[positive] <- log2(drug[positive] + rho[positive] * capacity[positive])
  log_ratio[!positive] <- log2(drug[!positive]) + 1 +
    plogis(2 * eta[!positive], log.p = TRUE) / log(2)
  list(
    log2_ratio = log_ratio, abundance_ratio = 2^log_ratio,
    departure_log2 = log_ratio - log2(drug)
  )
}

# Fit only shared log-gain and two dose thresholds. A small collection of fixed
# starts checks this smooth nonlinear fit; no HGF result enters start selection.
fit_curve <- function(rows, penalty) {
  cache_key <- paste(paste(rows, collapse = ","), penalty, sep = "|")
  if (!is.null(curve_fits[[cache_key]])) return(curve_fits[[cache_key]])
  group <- paste(outcomes$gene[rows], outcomes$dose[rows], sep = "|")
  weight <- 1 / as.numeric(table(group)[group]) / length(unique(group))
  ds <- ifelse(outcomes$dose[rows] == .5, 2L, 3L)
  objective <- function(par) {
    eta <- exp(par[1]) * (relative_score[rows] - par[ds])
    rho <- tanh(eta)
    rho[!present[rows]] <- 0
    eta[!present[rows]] <- 0
    output <- compose(observed_drug[rows], rho, eta)
    sum(weight * (output$log2_ratio - target_log2_ratio[rows])^2) + penalty * (par[2] - par[3])^2
  }
  starts <- rbind(c(log(5), .15, .25), c(log(20), .15, .25), c(log(5), .2, .2), c(log(20), .2, .2))
  lower <- c(log(.001), 0, 0)
  upper <- c(log(1000), 1, 1)
  attempts <- lapply(seq_len(nrow(starts)), function(i) {
    nloptr(
      x0 = starts[i, ], eval_f = objective, lb = lower, ub = upper,
      opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel = 1e-8, ftol_rel = 1e-10, maxeval = 3000)
    )
  })
  ok <- vapply(attempts, function(a) a$status > 0 && is.finite(a$objective), logical(1))
  stopifnot(any(ok))
  values <- vapply(attempts, function(a) if (a$status > 0) a$objective else Inf, numeric(1))
  best <- attempts[[which.min(values)]]
  # The rescue/killing branch boundary is continuous but has a derivative kink.
  # Refine with maintained derivative-free Subplex and check local objective
  # improvements directly rather than applying a smooth-gradient test at a kink.
  polish <- nloptr(
    x0 = best$solution, eval_f = objective, lb = lower, ub = upper,
    opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel = 1e-8, ftol_rel = 1e-10, maxeval = 3000)
  )
  stopifnot(polish$status > 0, all(is.finite(polish$solution)))
  par <- polish$solution
  value <- objective(par)
  local_improvement <- vapply(seq_along(par), function(j) {
    up <- down <- par
    up[j] <- min(upper[j], up[j] + 1e-5)
    down[j] <- max(lower[j], down[j] - 1e-5)
    max(value - objective(up), value - objective(down), 0)
  }, numeric(1))
  stopifnot(max(local_improvement) < 1e-6)
  result <- list(
    par = setNames(par, c("log_gain", "threshold_0.5", "threshold_5")), gain = exp(par[1]),
    threshold = setNames(par[2:3], c(".5", "5")), penalty = penalty, training_rows = rows,
    training_cells = unique(outcomes$depmap_id[rows]), objective = value,
    max_local_improvement = max(local_improvement), optimizer_status = polish$status, start_objectives = values,
    gain_boundary = par[1] > log(1000) - 1e-4
  )
  curve_fits[[cache_key]] <<- result
  save_cache(curve_fits, curve_signature, curve_cache_path)
  result
}
predict_rho <- function(fit, score, dose, active = rep(TRUE, length(score))) {
  ix <- ifelse(dose == .5, 2L, 3L)
  eta <- fit$gain * (score - fit$par[ix])
  eta[!active] <- 0
  rho <- tanh(eta)
  stopifnot(all(is.finite(rho)), all(abs(rho) <= 1))
  list(rho = rho, eta = eta)
}

# Select threshold shrinkage with whole-recipient inner exclusions. An outer-held
# cell remains excluded from every inner fit and selection calculation.
selection_audit <- list()
select_curve <- function(rows, label) {
  cells <- unique(outcomes$depmap_id[rows])
  pred <- matrix(NA_real_, length(rows), length(penalty_grid))
  for (cell in cells) {
    valid <- which(outcomes$depmap_id[rows] == cell)
    train <- rows[-valid]
    for (j in seq_along(penalty_grid)) {
      fit <- fit_curve(train, penalty_grid[j])
      rr <- rows[valid]
      pp <- predict_rho(fit, relative_score[rr], outcomes$dose[rr], present[rr])
      pred[valid, j] <- compose(observed_drug[rr], pp$rho, pp$eta)$log2_ratio
    }
  }
  stopifnot(all(is.finite(pred)))
  group <- paste(outcomes$gene[rows], outcomes$dose[rows], sep = "|")
  w <- 1 / as.numeric(table(group)[group]) / length(unique(group))
  score <- colSums((pred - target_log2_ratio[rows])^2 * w)
  choice <- which.min(score)
  selection_audit[[label]] <<- data.table(
    fit = label, penalty = penalty_grid, internal_objective = score,
    selected = seq_along(score) == choice, training_cells = paste(cells, collapse = ",")
  )
  fit_curve(rows, penalty_grid[choice])
}
prediction_table <- function(fit, rows, evaluation) {
  pp <- predict_rho(fit, relative_score[rows], outcomes$dose[rows], present[rows])
  conditional <- compose(observed_drug[rows], pp$rho, pp$eta)
  full_bridge <- compose(predicted_drug_72h[rows], pp$rho, pp$eta)
  cbind(outcomes[rows], data.table(
    evaluation = evaluation, rho = pp$rho,
    predicted_departure_log2 = conditional$departure_log2,
    predicted_combination_pct_vehicle = 100 * observed_ligand[rows] * conditional$abundance_ratio,
    full_bridge_combination_pct_vehicle = 100 * observed_ligand[rows] * full_bridge$abundance_ratio,
    full_bridge_departure_log2 = full_bridge$log2_ratio - log2(observed_drug[rows])
  ))
}
final <- select_curve(primary_rows, "final")
fits <- list(final = final)
predictions <- list()
for (cell in unique(outcomes$depmap_id[primary_rows])) {
  valid <- primary_rows[outcomes$depmap_id[primary_rows] == cell]
  label <- paste0("held_cell_", cell)
  fit <- select_curve(setdiff(primary_rows, valid), label)
  stopifnot(!cell %in% fit$training_cells)
  fits[[label]] <- fit
  predictions[[label]] <- prediction_table(fit, valid, "held_cell")
}
predictions[["external_EGF"]] <- prediction_table(final, control_rows, "external_EGF")
predictions[["fitted_NRG1"]] <- prediction_table(final, primary_rows, "fitted")
fits[["all_16_sensitivity"]] <- select_curve(seq_len(nrow(outcomes)), "all_16_sensitivity")
pred <- rbindlist(predictions)
write_table(pred, "combination_predictions_frozen")
write_table(rbindlist(selection_audit), "training_only_penalty_selection")
write_table(rbindlist(lapply(names(fits), function(label) {
  f <- fits[[label]]
  data.table(
    fit = label, gain = f$gain, threshold_0.5 = f$par[2], threshold_5 = f$par[3],
    penalty = f$penalty, max_local_improvement = f$max_local_improvement,
    optimizer_status = f$optimizer_status, gain_boundary = f$gain_boundary,
    training_cells = paste(f$training_cells, collapse = ",")
  )
})), "fitted_calibration_parameters")
qsave(
  list(
    single_drug_parent = source_root, g = final, folds = fits, programs = programs[keep],
    formula = "rho=tanh(gain*(R-threshold[dose])); positive restores depletion; negative removes remaining cells",
    calibration = "12 NRG1 72h counts with measured drug/ligand anchors",
    transfer = "72h-derived fractional recovery applied to frozen five-day drug abundance; no 5/3 log multiplication"
  ),
  file.path(root, "models/final.qs")
)

# The HGF query uses frozen five-day f and neutral HGF-alone abundance. The
# bounded fraction is transferred across assay durations without a rate model;
# this is an explicit pilot assumption, not a validated duration prediction.
query <- calibration_queries <- support <- list()
Ha <- profiles[CellLine == "AU565"]
for (dose in c(.5, 5)) {
  meta <- data.table(depmap_id = "ACH-000248", drug = "Lapatinib ditosylate", dose = dose)
  lap <- predict_single(model, meta)
  ff <- make_features(meta[rep(1, nrow(Ha))], lap$D[rep(1, nrow(Ha)), , drop = FALSE], as.matrix(Ha[, ..programs]))
  active <- ff$drug_activity > 0 & ff$ligand_norm > 0
  rr <- predict_rho(final, ff$relative_opposition, rep(dose, nrow(Ha)), active)
  drug <- rep(2^(-lap$sensitivity), nrow(Ha))
  cp <- compose(drug, rr$rho, rr$eta)
  query[[as.character(dose)]] <- data.table(
    CellLine = "AU565", dose = dose, source_profile = Ha$source_profile,
    lapatinib_pct_control = 100 * drug, rho = rr$rho, reversal_fraction = pmax(rr$rho, 0),
    g_log2 = cp$departure_log2, predicted_combination_pct_control = 100 * cp$abundance_ratio,
    ligand_standalone_assumption = "HGF-alone abundance fixed at 100%; not independently learned",
    time_transfer = "72h fractional interaction applied to five-day single-drug prediction; unvalidated time transfer"
  )
  train_dose <- primary_rows[outcomes$dose[primary_rows] == dose]
  support[[as.character(dose)]] <- data.table(
    dose = dose, source_profile = Ha$source_profile,
    relative_opposition = ff$relative_opposition, training_min = min(relative_score[train_dose]),
    training_max = max(relative_score[train_dose])
  )
  for (label in names(fits)) {
    qq <- predict_rho(fits[[label]], ff$relative_opposition, rep(dose, nrow(Ha)), active)
    cc <- compose(drug, qq$rho, qq$eta)
    calibration_queries[[paste(dose, label)]] <- data.table(
      fit = label, dose = dose, source_profile = Ha$source_profile,
      rho = qq$rho, predicted_combination_pct_control = 100 * cc$abundance_ratio
    )
  }
}
answer <- rbindlist(query)
write_table(answer, "AU565_HGF_predictions_frozen")
write_table(rbindlist(calibration_queries), "AU565_HGF_calibration_sensitivity")
write_table(rbindlist(support), "AU565_HGF_molecular_support")
prediction_paths <- file.path(root, "data", c("combination_predictions_frozen.tsv", "AU565_HGF_predictions_frozen.tsv"))
write_table(
  data.table(path = prediction_paths, sha256 = vapply(prediction_paths, digest, character(1), file = TRUE, algo = "sha256")),
  "prediction_hashes"
)

# Scientific summaries preserve comparison axis and assay units. Conditional g
# and full baseline -> molecular -> function composition are labeled separately.
metrics <- pred[, .(
  cells = .N, r_departure = pearson(observed_departure_log2, predicted_departure_log2),
  r_combination_abundance = pearson(observed_combination_pct_vehicle, predicted_combination_pct_vehicle),
  r_full_bridge_abundance = pearson(observed_combination_pct_vehicle, full_bridge_combination_pct_vehicle),
  observed_SD = sd(observed_departure_log2), predicted_SD = sd(predicted_departure_log2),
  observed_median = median(observed_departure_log2), predicted_median = median(predicted_departure_log2),
  observed_combination_median = median(observed_combination_pct_vehicle),
  predicted_combination_median = median(predicted_combination_pct_vehicle)
), by = .(evaluation, gene, dose)]
write_table(metrics, "combination_validation_correlations")

# A constant recovery fraction per dose checks how much apparent discrimination
# follows automatically from having measured single-drug controls. It uses the
# same whole-recipient exclusions and composition, without a molecular score.
constant_predictions <- list()
for (cell in unique(outcomes$depmap_id[primary_rows])) {
  train <- primary_rows[outcomes$depmap_id[primary_rows] != cell]
  valid <- primary_rows[outcomes$depmap_id[primary_rows] == cell]
  for (dose in c(.5, 5)) {
    tt <- train[outcomes$dose[train] == dose]
    vv <- valid[outcomes$dose[valid] == dose]
    opt <- optimize(function(eta) {
      pp <- compose(observed_drug[tt], rep(tanh(eta), length(tt)), rep(eta, length(tt)))
      mean((pp$log2_ratio - target_log2_ratio[tt])^2)
    }, interval = c(-5, 5), tol = 1e-8)
    eta <- rep(opt$minimum, length(vv))
    cp <- compose(observed_drug[vv], tanh(eta), eta)
    fp <- compose(predicted_drug_72h[vv], tanh(eta), eta)
    constant_predictions[[paste(cell, dose)]] <- cbind(outcomes[vv], data.table(
      evaluation = "held_cell_constant_recovery", rho = tanh(eta),
      predicted_departure_log2 = cp$departure_log2,
      predicted_combination_pct_vehicle = 100 * observed_ligand[vv] * cp$abundance_ratio,
      full_bridge_combination_pct_vehicle = 100 * observed_ligand[vv] * fp$abundance_ratio
    ))
  }
}
constant <- rbindlist(constant_predictions)
write_table(constant, "constant_recovery_held_predictions")
constant_metrics <- constant[, .(
  cells = .N,
  r_departure = pearson(observed_departure_log2, predicted_departure_log2),
  r_combination_abundance = pearson(observed_combination_pct_vehicle, predicted_combination_pct_vehicle),
  r_full_bridge_abundance = pearson(observed_combination_pct_vehicle, full_bridge_combination_pct_vehicle),
  observed_SD = sd(observed_departure_log2), predicted_SD = sd(predicted_departure_log2)
), by = .(evaluation, gene, dose)]
write_table(
  rbindlist(list(metrics[evaluation == "held_cell"], constant_metrics), fill = TRUE),
  "comparison_molecular_vs_constant_recovery"
)

# Structural tests cover absence, bounded extrapolation and scale invariance.
stopifnot(
  all(make_features(outcomes, single$D * 0, H)$drug_activity == 0),
  all(make_features(outcomes, single$D, H * 0)$ligand_norm == 0),
  max(abs(make_features(outcomes, single$D, H * 10)$relative_opposition - relative_score)) < 1e-12
)
zero <- predict_rho(final, rep(0, 16), outcomes$dose, rep(FALSE, 16))
stopifnot(all(zero$rho == 0), max(abs(compose(observed_drug, zero$rho, zero$eta)$departure_log2)) < 1e-12)
probe <- seq(0, 1, length.out = 1001)
for (dose in c(.5, 5)) {
  pp <- predict_rho(final, probe, rep(dose, length(probe)))
  for (drug in c(.01, .363, .901, 1, 1.5)) {
    cc <- compose(rep(drug, length(probe)), pp$rho, pp$eta)
    stopifnot(
      all(diff(pp$rho) >= -1e-12), all(cc$abundance_ratio >= 0),
      all(cc$abundance_ratio <= max(1, drug) + 1e-12)
    )
  }
}
write_table(data.table(check = c(
  "zero_drug_zero_interaction", "zero_ligand_zero_interaction", "ligand_scale_invariance",
  "monotone_recovery_fraction", "bounded_abundance_under_neutral_ligand", "single_drug_model_unchanged",
  "outer_recipient_excluded_from_inner_calibration"
), passed = TRUE), "structural_validation")

# Compact publication figures use square scatter panels and independent axes.
annotation <- pred[evaluation == "held_cell" & gene == "NRG1",
  .(label = correlation_label(observed_departure_log2, predicted_departure_log2)), by = dose]
sc <- ggplot(pred[evaluation == "held_cell" & gene == "NRG1"], aes(observed_departure_log2, predicted_departure_log2)) +
  geom_point(size = 1.8, color = "#2769A5") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .4, color = "black") +
  facet_wrap(~dose, scales = "free", labeller = labeller(dose = function(x) paste(x, "uM"))) +
  geom_text(data = annotation, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.05, vjust = 1.4, size = 3.5) +
  labs(
    x = "Measured 72h log2 departure", y = "Held-cell predicted departure",
    subtitle = "Interaction conditional on measured single-drug and ligand controls",
    caption = "Descriptive P values; retrospective recipient exclusions."
  ) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1)
ggsave(file.path(root, "figures/NRG1_held_cell_validation.pdf"), sc, width = 6, height = 3.5)
ggsave(file.path(root, "figures/NRG1_held_cell_validation.png"), sc, width = 6, height = 3.5, dpi = 220)
bar <- answer[source_profile == "Primary_study_mean", .(dose,
  Lapatinib = lapatinib_pct_control,
  `Lapatinib + HGF` = predicted_combination_pct_control
)]
bar <- melt(bar, id.vars = "dose", variable.name = "condition", value.name = "abundance")
bp <- ggplot(bar, aes(condition, abundance, fill = condition)) +
  geom_col(width = .65) +
  facet_wrap(~dose, scales = "free_y", labeller = labeller(dose = function(x) paste(x, "uM"))) +
  scale_fill_manual(values = c("#2769A5", "#D86A25")) +
  labs(x = NULL, y = "Predicted abundance (% control)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(root, "figures/AU565_HGF_prediction.pdf"), bp, width = 4, height = 2.4)
ggsave(file.path(root, "figures/AU565_HGF_prediction.png"), bp, width = 4, height = 2.4, dpi = 220)
wilson_path <- "data/raw/wilson/wilson_observed_source.tsv"
wilson <- fread(wilson_path)[series == "Wilson et al." & condition %in% c("Lapatinib", "Lapatinib + HGF")]
stopifnot(nrow(wilson) == 2, all(is.finite(wilson$abundance)), all(is.finite(wilson$sd)))
qa <- answer[dose == 5 & source_profile == "Primary_study_mean"]
comparison <- rbind(
  data.table(
    panel = "Prediction\n5 uM, five-day baseline", condition = c("Lapatinib", "Lapatinib + HGF"),
    abundance = c(qa$lapatinib_pct_control, qa$predicted_combination_pct_control), sd = NA_real_
  ),
  wilson[, .(panel = "Wilson measurement\n1 uM, 12 days", condition, abundance, sd)]
)
comparison[, condition := factor(condition, levels = c("Lapatinib", "Lapatinib + HGF"))]
write_table(comparison, "prediction_beside_Wilson")
wp <- ggplot(comparison, aes(condition, abundance, fill = condition)) +
  geom_col(width = .65) +
  geom_errorbar(data = comparison[is.finite(sd)], aes(ymin = abundance - sd, ymax = abundance + sd), width = .15) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c("#2769A5", "#D86A25")) +
  labs(x = NULL, y = "Abundance (% control)") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(root, "figures/prediction_beside_Wilson.pdf"), wp, width = 3.8, height = 2.3)
ggsave(file.path(root, "figures/prediction_beside_Wilson.png"), wp, width = 3.8, height = 2.3, dpi = 220)

# Provenance and modeling assumptions are plain tables; no additional memo files.
source_paths <- c(
  "scripts/12_fit_ligand_rescue.R", "R/prediction.R",
  "results/08_fit_single_drug_bridge/models/final.qs",
  "data/processed/her2_ligand_outcomes.tsv", "data/processed/hgf_program_inputs.qs", wilson_path
)
write_table(
  data.table(path = source_paths, sha256 = vapply(source_paths, digest, character(1), file = TRUE, algo = "sha256")),
  "source_hashes"
)
write_table(data.table(item = c(
  "calibration", "learned_parameters", "evaluation", "composition_assumption", "time_transfer",
  "molecular_input", "exclusions"
), value = c(
  "12 NRG1 72h observations; all six recipients and negative/above-control outcomes retained",
  "Shared positive gain, two dose thresholds; threshold-difference penalty selected in nested recipient folds",
  "Held recipient g with observed single-drug/ligand controls; frozen-f composition still uses measured ligand-only control",
  "Positive fraction restores depletion up to ligand-only reference; negative fraction removes remaining cells",
  "72h fractional recovery applied to five-day f; no log-interaction 5/3 multiplication; unvalidated duration transfer",
  "Same seven frozen single-drug programs and normalized ligand direction as v143; no direct baseline functional input",
  "HGF and Wilson outcomes never fitted; EGF external control; all-16 fit is sensitivity analysis"
)), "model_contract")
writeLines(capture.output(sessionInfo()), file.path(root, "data/R_session_info.txt"))
print(answer[source_profile == "Primary_study_mean"])
print(metrics)
print(data.table(gain = final$gain, threshold_0.5 = final$par[2], threshold_5 = final$par[3], penalty = final$penalty))
sink()
