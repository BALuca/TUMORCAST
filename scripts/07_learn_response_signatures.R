#!/usr/bin/env Rscript
# Training-only correlation prototypes -> two fixed scores -> shrunk drug readouts.
# Leave one cell out across all drugs; preserve v94's decoder and inner tuning.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(digest)
})
source("R/heatmaps.R")
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)

root <- "results/07_learn_response_signatures"
previous <- "results/06_prepare_functional_outcomes"
archive_path <- "data/processed/molecular_functional_inputs.rds"
programs <- sprintf("P%02d", 1:40)
seed <- 94012L
lambda_grid <- c(.01, .1, 1, 10)
cache_signature <- list(
  inputs = vapply(c(archive_path, file.path(previous, "data/individual_training_outcomes.tsv"),
    file.path(previous, "data/selected_doses.tsv")), digest, character(1), file = TRUE, algo = "sha256"),
  code = vapply(c("scripts/07_learn_response_signatures.R", "R/heatmaps.R"), digest,
    character(1), file = TRUE, algo = "sha256"),
  seed = seed, lambda_grid = lambda_grid,
  glmnet = as.character(packageVersion("glmnet")), R = R.version.string
)
for (folder in c("data", "models", "figures")) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
sink(file.path(root, "run.log"), split = TRUE)
write_table <- function(x, name) {
  fwrite(x, file.path(root, "data", paste0(name, ".tsv")), sep = "\t")
}
status <- function(message) {
  cat(message, "\n")
  flush.console()
  writeLines(message, file.path(root, "STATUS.txt"))
}
pearson <- function(x, y) {
  if (length(x) < 4 || sd(x) < 1e-10 || sd(y) < 1e-10) {
    return(NA_real_)
  }
  cor(x, y)
}

selection <- fread(file.path(previous, "data/selected_doses.tsv"))
original_folds <- fread(file.path(previous, "data/recipient_folds.tsv"))
outcomes <- fread(file.path(previous, "data/individual_training_outcomes.tsv"))[assay == "secondary"]
drugs <- selection$drug
short <- c(
  "Azacytidine", "Bendamustine", "Bortezomib", "Busulfan", "Cytarabine",
  "Lapatinib", "Topotecan", "Vandetanib", "Vismodegib", "Crizotinib"
)
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
molecular <- merge(as.data.table(archive$prep$m), selection, by = c("drug", "dose"))[
  depmap_id != "ACH-000248"
]
profiles <- merge(molecular, outcomes[, .(depmap_id, drug, dose, sensitivity = -y)],
  by = c("depmap_id", "drug", "dose")
)
loo_cells <- sort(unique(profiles$depmap_id))
folds <- data.table(depmap_id = loo_cells, fold = seq_along(loo_cells))
profiles[, fold := folds$fold[match(depmap_id, folds$depmap_id)]]
stopifnot(
  nrow(profiles) == 286L, length(drugs) == 10L,
  !anyDuplicated(profiles[, .(depmap_id, drug)]), !anyNA(profiles$fold),
  !anyDuplicated(baseline$depmap_id), all(is.finite(as.matrix(profiles[, ..programs])))
)
rm(archive)
write_table(selection, "selected_doses")
write_table(folds, "recipient_folds")
stopifnot(
  length(loo_cells) == 30L,
  setequal(intersect(unique(molecular$depmap_id), unique(outcomes$depmap_id)), loo_cells)
)
status("v95: 286 measured profiles; ten selected contexts; 30 leave-one-cell-out fits.")

# Build two prototypes from training correlations, keeping the program weights fixed.
# Route 1 is always the drug cluster containing lapatinib, solely for stable labels.
learn_templates <- function(train) {
  mu <- sigma <- correlations <- matrix(NA_real_, 10, 40, dimnames = list(drugs, programs))
  reference <- numeric(10)
  counts <- integer(10)
  for (g in seq_along(drugs)) {
    rows <- which(train$drug == drugs[g])
    stopifnot(length(rows) >= 5L)
    D <- as.matrix(train[rows, ..programs])
    mu[g, ] <- colMeans(D)
    sigma[g, ] <- apply(D, 2, sd)
    correlations[g, ] <- as.numeric(cor(D, train$sensitivity[rows]))
    reference[g] <- mean(train$sensitivity[rows])
    counts[g] <- length(rows)
  }
  stopifnot(all(is.finite(correlations)), all(sigma > 1e-10))
  tree <- hclust(dist(correlations), method = "ward.D2")
  cluster <- cutree(tree, 2)
  lap_cluster <- cluster[match("Lapatinib ditosylate", drugs)]
  membership <- ifelse(cluster == lap_cluster, 1L, 2L)
  collapsed <- vapply(1:2, function(k) colMeans(correlations[membership == k, , drop = FALSE]), numeric(40))
  dimnames(collapsed) <- list(programs, c("Signature 1", "Signature 2"))
  norms <- sqrt(colSums(collapsed^2))
  stopifnot(all(norms > 1e-10))
  W <- sweep(collapsed, 2, norms, "/")
  g <- match(train$drug, drugs)
  Z <- (as.matrix(train[, ..programs]) - mu[g, , drop = FALSE]) / sigma[g, , drop = FALSE]
  raw_scores <- Z %*% W
  score_sd <- t(vapply(1:10, function(j) apply(raw_scores[g == j, , drop = FALSE], 2, sd), numeric(2)))
  stopifnot(all(score_sd > 1e-10))
  list(
    mu = mu, sigma = sigma, correlations = t(correlations), collapsed = collapsed, W = W,
    score_sd = score_sd, reference = reference, counts = counts, membership = membership,
    drug_tree = tree, training_cells = unique(train$depmap_id)
  )
}

# Signature score calibration uses training RNA variance, without functional shifts.
scores <- function(template, D, drug) {
  g <- match(drug, drugs)
  stopifnot(!anyNA(g), ncol(D) == 40L)
  Z <- (D - template$mu[g, , drop = FALSE]) / template$sigma[g, , drop = FALSE]
  (Z %*% template$W) / template$score_sd[g, , drop = FALSE]
}

# Global two-coefficient readout + shrunk two-coefficient deviations per drug.
# No intercept: the fitted target is sensitivity relative to its training mean.
fit_readout <- function(template, train, lambda) {
  S <- scores(template, as.matrix(train[, ..programs]), train$drug)
  g <- match(train$drug, drugs)
  target <- train$sensitivity - template$reference[g]
  design <- matrix(0, nrow(train), 22)
  design[, 1:2] <- S
  for (j in 1:10) design[g == j, 2 + (2 * j - 1):(2 * j)] <- S[g == j, , drop = FALSE]
  weights <- 1 / (10 * template$counts[g])
  penalty <- diag(c(rep(.001, 2), rep(lambda / 10, 20)))
  beta <- solve(crossprod(design, design * weights) + penalty, crossprod(design, target * weights))
  coefficients <- t(matrix(beta[-(1:2)], 2, 10))
  coefficients <- sweep(coefficients, 2, beta[1:2], "+")
  list(coefficients = coefficients, shared = as.numeric(beta[1:2]), lambda = lambda)
}
predict_readout <- function(template, readout, D, drug) {
  S <- scores(template, D, drug)
  g <- match(drug, drugs)
  rowSums(S * readout$coefficients[g, , drop = FALSE])
}

# Inner whole-cell folds select only the shrinkage strength. Rebuild prototypes
# in every inner split: inner validation labels never enter its correlation map.
tune_readout <- function(train, job) {
  set.seed(seed + job)
  cells <- sample(sort(unique(train$depmap_id)))
  inner_fold <- setNames(rep(1:3, length.out = length(cells)), cells)
  losses <- list()
  for (k in 1:3) {
    held <- cells[inner_fold == k]
    inner_train <- train[!depmap_id %in% held]
    validation <- train[depmap_id %in% held]
    template <- learn_templates(inner_train)
    stopifnot(!any(held %in% template$training_cells))
    for (lambda in lambda_grid) {
      readout <- fit_readout(template, inner_train, lambda)
      predicted <- predict_readout(template, readout, as.matrix(validation[, ..programs]), validation$drug)
      error <- predicted - (validation$sensitivity - template$reference[match(validation$drug, drugs)])
      loss <- data.table(drug = validation$drug, error)[, .(loss = mean(error^2)), by = drug]
      losses[[length(losses) + 1L]] <- data.table(inner_fold = k, lambda = lambda, loss = mean(loss$loss))
    }
  }
  audit <- rbindlist(losses)
  average <- audit[, .(loss = mean(loss)), by = lambda]
  best <- average[order(loss, -lambda), lambda][1]
  list(lambda = best, audit = audit)
}

# RNA-only baseline -> displacement mapping, fixed independently of viability.
# This secondary replay distinguishes decoder failure from displacement prediction.
fit_rna <- function(train_rna) {
  ids <- sort(unique(train_rna$depmap_id))
  raw <- as.matrix(baseline[match(ids, depmap_id), ..programs])
  bmean <- colMeans(raw)
  bsd <- apply(raw, 2, sd)
  normalized <- sweep(sweep(raw, 2, bmean), 2, bsd, "/")
  basis <- svd(normalized, nu = 0, nv = 8)$v
  bscale <- sqrt(colMeans((normalized %*% basis)^2))
  B <- array(0, c(8, 40, 10))
  mu <- sigma <- matrix(0, 10, 40)
  bm <- matrix(0, 10, 8)
  for (g in 1:10) {
    rows <- which(train_rna$drug == drugs[g])
    D <- as.matrix(train_rna[rows, ..programs])
    mu[g, ] <- colMeans(D)
    sigma[g, ] <- apply(D, 2, sd)
    features <- sweep(normalized[match(train_rna$depmap_id[rows], ids), , drop = FALSE] %*% basis, 2, bscale, "/")
    bm[g, ] <- colMeans(features)
    features <- sweep(features, 2, bm[g, ])
    target <- sweep(sweep(D, 2, mu[g, ]), 2, sigma[g, ], "/")
    fit <- glmnet(features, target,
      family = "mgaussian", alpha = 0, lambda = 1,
      intercept = FALSE, standardize = FALSE
    )
    B[, , g] <- vapply(coef(fit), function(x) as.numeric(x)[-1], numeric(8))
  }
  list(
    bmean = bmean, bsd = bsd, basis = basis, bscale = bscale, B = B,
    mu = mu, sigma = sigma, bm = bm, training_cells = ids
  )
}
predict_rna <- function(model, cells, drug) {
  raw <- as.matrix(baseline[match(cells, depmap_id), ..programs])
  stopifnot(all(is.finite(raw)))
  features <- sweep(
    sweep(sweep(raw, 2, model$bmean), 2, model$bsd, "/") %*% model$basis,
    2, model$bscale, "/"
  )
  g <- match(drug, drugs)
  D <- matrix(0, length(cells), 40)
  for (j in unique(g)) {
    rows <- which(g == j)
    Z <- sweep(features[rows, , drop = FALSE], 2, model$bm[j, ]) %*% model$B[, , j]
    D[rows, ] <- sweep(sweep(Z, 2, model$sigma[j, ], "*"), 2, model$mu[j, ], "+")
  }
  D
}

# Every fold shows its own observed training map and its two collapsed columns.
render_heatmap <- function(template, label) {
  plot_training_signatures(template, short, selection$dose, label)
}

measured_predictions <- baseline_predictions <- coefficient_tables <- membership_tables <- list()
template_tables <- tuning_tables <- audits <- list()
book_path <- file.path(root, "figures/training_correlations_and_collapsed_signatures_per_fold.pdf")
cairo_pdf(book_path, width = 9, height = 9, onefile = TRUE)
book_device <- dev.cur()
final_job <- length(loo_cells) + 1L
for (job in seq_len(final_job)) {
  label <- if (job == final_job) "final" else paste0("loo_", sprintf("%02d", job))
  held <- if (job == final_job) character() else loo_cells[job]
  train <- profiles[!depmap_id %in% held]
  test <- if (job == final_job) profiles else profiles[depmap_id %in% held]
  status(paste(label, ": deriving training-only templates and choosing shrinkage in inner cell folds."))
  model_path <- file.path(root, "models", paste0(label, ".qs"))
  signature <- c(cache_signature, list(held = held, job = job))
  cached <- if (file.exists(model_path)) qread(model_path) else NULL
  reuse_fit <- !is.null(cached) && identical(cached$cache_signature, signature)
  if (reuse_fit) {
    status(paste(label, ": reusing matching fit"))
    template <- cached$template
    tuning <- cached$tuning
    readout <- cached$readout
    rna <- cached$rna
  } else {
    status(paste(label, ": absent or stale fit; rebuilding"))
    template <- learn_templates(train)
    tuning <- tune_readout(train, job)
    readout <- fit_readout(template, train, tuning$lambda)
    rna <- fit_rna(molecular[!depmap_id %in% held])
  }
  stopifnot(!any(held %in% template$training_cells), !any(held %in% rna$training_cells))
  model <- list(template = template, readout = readout, rna = rna, held_cells = held,
    selected = selection, tuning = tuning, cache_signature = signature)
  if (!reuse_fit) {
    qsave(model, paste0(model_path, ".partial"))
    stopifnot(file.rename(paste0(model_path, ".partial"), model_path))
  }
  evaluation <- if (job == final_job) "fitted" else "held_cell"
  predicted <- predict_readout(template, readout, as.matrix(test[, ..programs]), test$drug)
  p <- test[, .(depmap_id, drug, dose, observed_sensitivity = sensitivity)]
  p[, observed_deviation := observed_sensitivity - template$reference[match(drug, drugs)]]
  p[, `:=`(predicted_deviation = predicted, fold = job, evaluation = evaluation)]

  # Individual-drug ridge is a diagnostic comparator on exactly the same cells.
  comparator <- numeric(nrow(test))
  for (g in 1:10) {
    tr <- which(train$drug == drugs[g])
    te <- which(test$drug == drugs[g])
    if (!length(te)) next
    Ztr <- sweep(sweep(as.matrix(train[tr, ..programs]), 2, template$mu[g, ]), 2, template$sigma[g, ], "/")
    Zte <- sweep(sweep(as.matrix(test[te, ..programs]), 2, template$mu[g, ]), 2, template$sigma[g, ], "/")
    fit <- glmnet(Ztr, train$sensitivity[tr] - template$reference[g],
      alpha = 0, lambda = .1,
      intercept = FALSE, standardize = FALSE
    )
    comparator[te] <- as.numeric(predict(fit, Zte, s = .1))
  }
  p[, individual_ridge_deviation := comparator]
  measured_predictions[[label]] <- p
  all_test <- if (job == final_job) outcomes else outcomes[depmap_id %in% held]
  Dhat <- predict_rna(rna, all_test$depmap_id, all_test$drug)
  bp <- all_test[, .(depmap_id, drug, dose, observed_sensitivity = -y)]
  bp[, observed_deviation := observed_sensitivity - template$reference[match(drug, drugs)]]
  bp[, `:=`(
    predicted_deviation = predict_readout(template, readout, Dhat, all_test$drug),
    fold = job, evaluation = evaluation
  )]
  baseline_predictions[[label]] <- bp
  coefficient_tables[[label]] <- data.table(
    fit = label, drug = drugs, dose = selection$dose,
    signature_1_beta = readout$coefficients[, 1], signature_2_beta = readout$coefficients[, 2],
    shared_1 = readout$shared[1], shared_2 = readout$shared[2], lambda = tuning$lambda
  )
  membership_tables[[label]] <- data.table(
    fit = label, drug = drugs, signature = template$membership,
    training_cells = template$counts
  )
  template_tables[[label]] <- data.table(
    fit = label, program = programs,
    collapsed_r_1 = template$collapsed[, 1], collapsed_r_2 = template$collapsed[, 2],
    weight_1 = template$W[, 1], weight_2 = template$W[, 2]
  )
  tuning_tables[[label]] <- copy(tuning$audit)[, fit := label]
  audits[[label]] <- data.table(
    fit = label, training_cells = uniqueN(train$depmap_id),
    held_cell = if (job == final_job) NA_character_ else held,
    held_measured_cells = if (job == final_job) 0L else uniqueN(test$depmap_id),
    held_cell_leakage = length(intersect(held, template$training_cells)), lambda = tuning$lambda
  )
  dev.set(book_device)
  render_heatmap(template, if (job == final_job) label else paste(label, "| held:", held))
  # The PDF contains all 30 folds. Only representative previews are materialized.
  if (job %in% c(1L, 15L, 30L, final_job)) {
    png(file.path(root, "figures", paste0(label, "_training_correlations_and_signatures.png")),
      width = 9, height = 9, units = "in", res = 150, type = "cairo"
    )
    render_heatmap(template, if (job == final_job) label else paste(label, "| held:", held))
    dev.off()
  }
}
dev.set(book_device)
dev.off()
measured <- rbindlist(measured_predictions)
replay <- rbindlist(baseline_predictions)
# Cells outside the RNA/function training set do not require duplicate fits:
# excluding their data changes neither the decoder nor the RNA predictor.
external <- copy(replay[evaluation == "fitted" & !depmap_id %in% loo_cells])
stopifnot(
  !any(external$depmap_id %in% molecular$depmap_id),
  !any(external$depmap_id %in% profiles$depmap_id)
)
external[, evaluation := "held_cell"]
replay <- rbind(replay, external)
stopifnot(
  !anyDuplicated(replay[evaluation == "held_cell", .(depmap_id, drug)]),
  nrow(replay[evaluation == "held_cell"]) == nrow(outcomes),
  nrow(measured[evaluation == "held_cell"]) == nrow(profiles)
)
write_table(measured, "measured_predictions_frozen")
write_table(replay, "baseline_predictions_frozen")
write_table(rbindlist(coefficient_tables), "drug_signature_coefficients")
write_table(rbindlist(membership_tables), "drug_cluster_memberships")
write_table(rbindlist(template_tables), "collapsed_signature_weights")
write_table(rbindlist(tuning_tables), "inner_shrinkage_selection")
write_table(rbindlist(audits), "fold_audit")
prediction_paths <- file.path(root, "data", c("measured_predictions_frozen.tsv", "baseline_predictions_frozen.tsv"))
write_table(data.table(path = prediction_paths, sha256 = vapply(prediction_paths, digest, character(1),
  file = TRUE,
  algo = "sha256"
)), "prediction_hashes")

# Evaluate collected LOO predictions across cells at fixed drug/dose. Centering
# inside a one-cell fold would erase the entire signal and is never performed.
metrics <- measured[, .(
  cells = .N, r = pearson(observed_sensitivity, predicted_deviation),
  r_pooled_deviation = pearson(observed_deviation, predicted_deviation),
  r_individual_ridge = pearson(observed_sensitivity, individual_ridge_deviation)
), by = .(drug, dose, evaluation)]
baseline_metrics <- replay[, .(
  cells = .N, r = pearson(observed_sensitivity, predicted_deviation),
  r_pooled_deviation = pearson(observed_deviation, predicted_deviation)
), by = .(drug, dose, evaluation)]
write_table(metrics, "measured_decoder_correlations")
write_table(baseline_metrics, "baseline_bridge_correlations")
memberships <- rbindlist(membership_tables)[fit != "final"]
stability <- memberships[, .(LOO_fits = .N, fraction_in_lapatinib_cluster = mean(signature == 1L)), by = drug]
write_table(stability, "cluster_membership_stability")

plot_predictions <- function(table, endpoint, stem) {
  data <- copy(table[evaluation == "held_cell"])
  data[, panel := factor(paste0(short[match(drug, drugs)], " | ", dose, " µM"),
    levels = paste0(short, " | ", selection$dose, " µM")
  )]
  stats <- data[, .(label = correlation_label(observed_sensitivity, predicted_deviation)), by = panel]
  plot <- ggplot(data, aes(observed_sensitivity, predicted_deviation)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = .3) +
    geom_point(color = "#0072B2", size = 1.5, alpha = .8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey30", linewidth = .5) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -.05, vjust = 1.15, size = 3.5) +
    facet_wrap(~panel, ncol = 5, scales = "free") +
    theme_bw(base_size = 12) +
    theme(aspect.ratio = 1, legend.position = "bottom", panel.grid.minor = element_blank()) +
    labs(
      title = paste0("v95 | ", endpoint), x = "Measured sensitivity (-log2 PRISM abundance)",
      y = "Predicted sensitivity deviation",
      caption = paste0("Pearson r across held-cell predictions at fixed drug/dose. P values are descriptive.\n",
        "Retrospective comparisons; no centering within one-cell folds.\n",
        "Other PRISM cells use the final model; no drug mean is added to predictions.")
    )
  ggsave(file.path(root, "figures", paste0(stem, ".pdf")), plot, width = 13, height = 6.8, device = cairo_pdf)
  ggsave(file.path(root, "figures", paste0(stem, ".png")), plot, width = 13, height = 6.8, dpi = 150)
}
plot_predictions(measured, "held-cell decoder supplied with measured displacements", "held_measured_decoder")
plot_predictions(replay, "held-cell baseline → RNA-only displacement → fixed signatures", "held_baseline_bridge")
writeLines(
  c(
    "v95 preserves the v94 correlation-collapsed decoder and replaces outer three-fold splits with leave-one-cell-out.",
    "Thirty outer fits hold out each matched cell across every drug and both molecular/function data. AU565 remains excluded.",
    "Training-only signed Pearson maps -> Ward.D2 two drug clusters -> equal-drug mean correlation templates.",
    "Template vectors are L2 normalized and frozen. Scores use measured training drug/program mean and SD, then training drug/score SD.",
    "Only two signed functional coefficients per drug are fitted, shrunk to shared coefficients. No functional intercept or program-wise learned shifts/gains.",
    "Target is within-drug sensitivity deviation, not absolute viability. Drug mean is a reference for observed deviations only, never added to predictions.",
    "Lambda .01/.1/1/10 is selected by inner whole-cell folds, rebuilding templates in each inner training split.",
    "Individual-drug 40-program ridge lambda .1 is a fixed diagnostic comparator, not part of the main predictor.",
    "Secondary baseline replay fits RNA-only 8-axis baseline -> 40 displacement ridge lambda 1, using outer-training RNA only.",
    "Only measured RNA/PRISM matched cells train the functional decoder. Other PRISM cells are secondary evaluation; no ALMANAC or combination refitting here.",
    "Primary statistic: Pearson across collected LOO predictions within each fixed drug/dose; no centering within one-cell folds.",
    "All PRISM cells outside molecular/function training use the final model; their exclusion cannot change fitting. This covers the original 441 lapatinib cells.",
    "Dose selection and previous test labels were already exposed; these are retrospective comparisons, not independent confirmation."
  ),
  file.path(root, "RUN_DETAILS.txt")
)
capture.output(sessionInfo(), file = file.path(root, "SESSION_INFO.txt"))
inputs <- c(
  archive_path, file.path(previous, "data", c(
    "selected_doses.tsv", "recipient_folds.tsv",
    "individual_training_outcomes.tsv"
  )),
  "scripts/07_learn_response_signatures.R"
)
write_table(data.table(path = inputs, sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = "sha256")), "source_hashes")
print(metrics[evaluation == "held_cell"])
print(baseline_metrics[evaluation == "held_cell"])
status("Complete: 30 LOO correlation maps and collapsed signatures, held-cell measured decoder and full PRISM baseline replay.")
sink()
