#!/usr/bin/env Rscript
# Existing two signatures and training-only route assignments are fixed.
# Eight shared coefficients predict molecular-footprint level and positive
# amplitude. No per-drug viability coefficients or outcome-fitted program gains.
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

base_root <- "results/08_fit_single_drug_bridge"
root <- base_root
parent <- "results/07_learn_response_signatures"
source_root <- "results/06_prepare_functional_outcomes/data"
archive_path <- "data/processed/molecular_functional_inputs.rds"
programs <- sprintf("P%02d", 1:40)
signature_amplitude <- .25
calibration_penalty <- .005
assay_penalty <- .01
interaction_penalty <- .05
excluded_drugs <- c("crizotinib", "Vismodegib")
cache_signature <- list(
  inputs = vapply(c(archive_path, list.files(file.path(source_root), full.names = TRUE),
    list.files(file.path(parent, "models"), full.names = TRUE)), digest,
    character(1), file = TRUE, algo = "sha256"),
  code = vapply(c("scripts/08_fit_single_drug_bridge.R", "R/heatmaps.R"), digest,
    character(1), file = TRUE, algo = "sha256"),
  penalties = c(calibration_penalty, assay_penalty, interaction_penalty),
  amplitude = signature_amplitude, R = R.version.string,
  glmnet = as.character(packageVersion("glmnet"))
)
for (folder in c("data", "models", "figures")) dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
sink(file.path(root, "run.log"), split = TRUE)
write_table <- function(x, name) fwrite(x, file.path(root, "data", paste0(name, ".tsv")), sep = "\t")
status <- function(message) {
  cat(format(Sys.time()), message, "\n")
  flush.console()
  writeLines(message, file.path(root, "STATUS.txt"))
}
pearson <- function(x, y) if (length(x) >= 4 && sd(x) > 1e-10 && sd(y) > 1e-10) cor(x, y) else NA_real_
key <- function(cell, drug, dose) paste(cell, drug, sprintf("%.10g", dose), sep = "|")
context_key <- function(drug, dose) paste(drug, sprintf("%.10g", dose), sep = "|")
selection <- fread(file.path(source_root, "selected_doses.tsv"))
single <- fread(file.path(source_root, "individual_training_outcomes.tsv"))
pairs <- fread(file.path(source_root, "combination_training_outcomes.tsv"))
parent_folds <- fread(file.path(parent, "data/recipient_folds.tsv"))
bank_selection <- copy(selection)
drugs <- bank_selection$drug
stopifnot(all(excluded_drugs %in% drugs))
selection <- selection[!drug %in% excluded_drugs]
included_indices <- match(selection$drug, drugs)
single <- single[!drug %in% excluded_drugs]
pairs <- pairs[!drugA %in% excluded_drugs & !drugB %in% excluded_drugs]
short <- c(
  "Azacytidine", "Bendamustine", "Bortezomib", "Busulfan", "Cytarabine",
  "Lapatinib", "Topotecan", "Vandetanib", "Vismodegib", "Crizotinib"
)
archive <- readRDS(archive_path)
baseline <- as.data.table(archive$prep$raw)
receptors <- as.data.table(archive$prep$receptors)
stopifnot(!anyDuplicated(receptors$depmap_id))
molecular_all <- as.data.table(archive$prep$m)[depmap_id != "ACH-000248"]
molecular <- merge(molecular_all, bank_selection, by = c("drug", "dose"))
# Add only lapatinib pairs at supported native molecular doses; no outcome-based selection.
raw_pairs <- as.data.table(archive$pair)[
  depmap_id != "ACH-000248" & drugA %in% selection$drug & drugB %in% selection$drug &
    (drugA == "Lapatinib ditosylate" | drugB == "Lapatinib ditosylate")
]
added_pairs <- coverage <- list()
for (pair_name in unique(raw_pairs$pair)) {
  z <- raw_pairs[pair == pair_name]
  lap_is_A <- z$drugA[1] == "Lapatinib ditosylate"
  partner <- if (lap_is_A) z$drugB[1] else z$drugA[1]
  partner_grid <- sort(unique(if (lap_is_A) z$doseB else z$doseA))
  partner_native <- sort(unique(molecular_all[drug == partner, dose]))
  candidates <- partner_native[partner_native >= min(partner_grid) & partner_native <= max(partner_grid)]
  if (!length(candidates)) {
    coverage[[pair_name]] <- data.table(pair = pair_name, partner = partner, disposition = "no native molecular dose inside cached outcome grid")
    next
  }
  preferred <- selection[drug == partner, dose]
  partner_target <- candidates[which.min(abs(log(candidates / preferred)))]
  targetA <- if (lap_is_A) .5 else partner_target
  targetB <- if (lap_is_A) partner_target else .5
  coverage[[pair_name]] <- data.table(pair = pair_name, partner = partner, disposition = "included", doseA = targetA, doseB = targetB)
  for (cell in unique(z$depmap_id)) {
    zz <- z[depmap_id == cell]
    weights <- function(grid, target) {
      grid <- sort(unique(grid))
      if (target < min(grid) || target > max(grid)) {
        return(NULL)
      }
      if (any(abs(grid - target) < 1e-10)) {
        return(data.table(dose = grid[which.min(abs(grid - target))], weight = 1))
      }
      lo <- max(grid[grid < target])
      hi <- min(grid[grid > target])
      f <- log(target / lo) / log(hi / lo)
      data.table(dose = c(lo, hi), weight = c(1 - f, f))
    }
    wa <- weights(zz$doseA, targetA)
    wb <- weights(zz$doseB, targetB)
    if (is.null(wa) || is.null(wb)) stop("Incomplete dose bracket in added pair")
    corners <- CJ(ia = seq_len(nrow(wa)), ib = seq_len(nrow(wb)))
    corners[, `:=`(doseA = wa$dose[ia], doseB = wb$dose[ib], weight = wa$weight[ia] * wb$weight[ib])]
    corners <- merge(corners, zz[, .(doseA, doseB, y, yA, yB)], by = c("doseA", "doseB"), all.x = TRUE)
    stopifnot(!anyNA(corners), abs(sum(corners$weight) - 1) < 1e-12)
    out <- zz[1, .(depmap_id, pair, drugA, drugB, CELLNAME, PANEL, SCREENER)]
    out[, `:=`(
      doseA = targetA, doseB = targetB, y = sum(corners$weight * corners$y),
      yA = sum(corners$weight * corners$yA), yB = sum(corners$weight * corners$yB)
    )]
    out[, `:=`(I = y - yA - yB, partner = partner, interpolation_corners = nrow(corners))]
    added_pairs[[paste(pair_name, cell)]] <- out
  }
}
added <- rbindlist(added_pairs)
stopifnot(nrow(added) == 192L, uniqueN(added$pair) == 5L)
pairs <- rbind(pairs, added, fill = TRUE)
extra_contexts <- unique(rbind(added[, .(drug = drugA, dose = doseA)], added[, .(drug = drugB, dose = doseB)]))
extra_contexts <- extra_contexts[!context_key(drug, dose) %in% context_key(bank_selection$drug, bank_selection$dose)]
setorder(extra_contexts, drug, dose)
stopifnot(nrow(extra_contexts) == 3L)
contexts <- rbind(bank_selection, extra_contexts)
context_index <- function(drug, dose) match(context_key(drug, dose), context_key(contexts$drug, contexts$dose))
write_table(rbindlist(coverage, fill = TRUE), "added_lapatinib_pair_coverage")
write_table(contexts, "molecular_contexts")
write_table(pairs, "combination_training_outcomes")
prism <- single[assay == "secondary"]
nci <- single[assay != "secondary" & depmap_id != "ACH-000248"]
extra_singles <- unique(rbind(
  added[, .(depmap_id, drug = drugA, dose = doseA, y = yA, assay = paste0("nci_", SCREENER))],
  added[, .(depmap_id, drug = drugB, dose = doseB, y = yB, assay = paste0("nci_", SCREENER))]
))
extra_singles[, `:=`(
  source = "nci_single", control_key = paste(assay, key(depmap_id, drug, dose), sep = "|"),
  group = paste(assay, drug, dose, sep = "|")
)]
nci <- unique(rbind(nci, extra_singles, fill = TRUE), by = c("depmap_id", "drug", "dose", "assay", "y"))
write_table(nci, "nci_single_training_outcomes")
profiles <- merge(molecular, prism[, .(depmap_id, drug, dose, S = -y)], by = c("depmap_id", "drug", "dose"))
stopifnot(
  nrow(selection) == 8L, nrow(profiles) == 226L, nrow(pairs) == 307L,
  !any(prism$drug %in% excluded_drugs),
  !any(c(pairs$drugA, pairs$drugB) %in% excluded_drugs)
)
rm(archive)
stopifnot(
  nrow(profiles) > 0L, nrow(pairs) > 0L,
  !anyDuplicated(profiles[, .(depmap_id, drug)]), !anyDuplicated(baseline$depmap_id),
  all(is.finite(as.matrix(profiles[, ..programs]))), all(is.finite(single$y))
)
loo_cells <- sort(unique(c(profiles$depmap_id, nci$depmap_id, pairs$depmap_id)))
folds <- data.table(depmap_id = loo_cells, fold = seq_along(loo_cells))
write_table(selection, "selected_doses")
write_table(folds, "recipient_folds")

# Lapatinib uses ERBB2/HER2 and MET only, through predicted displacements.
# Other drugs retain the v95 RNA-only map. Nonmatched RNA recipients need
# the same RNA fit with their own cell removed; no functional labels enter it.
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
    fit <- glmnet(features, target, family = "mgaussian", alpha = 0, lambda = 1, intercept = FALSE, standardize = FALSE)
    B[, , g] <- vapply(coef(fit), function(x) as.numeric(x)[-1], numeric(8))
  }
  list(bmean = bmean, bsd = bsd, basis = basis, bscale = bscale, B = B, mu = mu, sigma = sigma, bm = bm, training_cells = ids)
}
fit_lapatinib_receptors <- function(train_rna) {
  train <- train_rna[drug == "Lapatinib ditosylate"]
  X <- as.matrix(receptors[match(train$depmap_id, depmap_id), .(ERBB2, MET)])
  D <- as.matrix(train[, ..programs])
  xm <- colMeans(X)
  xs <- apply(X, 2, sd)
  dm <- colMeans(D)
  ds <- pmax(apply(D, 2, sd), 1e-8)
  X <- sweep(sweep(X, 2, xm), 2, xs, "/")
  Y <- sweep(sweep(D, 2, dm), 2, ds, "/")
  fit <- glmnet(X, Y, family = "mgaussian", alpha = 0, lambda = 1, intercept = FALSE, standardize = FALSE)
  B <- vapply(coef(fit), function(x) as.numeric(x)[-1], numeric(2))
  list(xmean = xm, xsd = xs, mu = dm, sigma = ds, B = B, training_cells = train$depmap_id, lambda = 1)
}
fit_extra_rna <- function(rna, held) {
  maps <- list()
  for (i in seq_len(nrow(extra_contexts))) {
    dr <- extra_contexts$drug[i]
    ds <- extra_contexts$dose[i]
    train <- molecular_all[drug == dr & dose == ds & !depmap_id %in% held]
    stopifnot(nrow(train) >= 10L, !anyDuplicated(train$depmap_id))
    if (dr == "Lapatinib ditosylate") {
      map <- fit_lapatinib_receptors(train)
      map$kind <- "receptors"
    } else {
      raw <- as.matrix(baseline[match(train$depmap_id, depmap_id), ..programs])
      X <- sweep(sweep(sweep(raw, 2, rna$bmean), 2, rna$bsd, "/") %*% rna$basis, 2, rna$bscale, "/")
      xm <- colMeans(X)
      X <- sweep(X, 2, xm)
      D <- as.matrix(train[, ..programs])
      dm <- colMeans(D)
      sdv <- pmax(apply(D, 2, sd), 1e-8)
      Y <- sweep(sweep(D, 2, dm), 2, sdv, "/")
      fit <- glmnet(X, Y, family = "mgaussian", alpha = 0, lambda = 1, intercept = FALSE, standardize = FALSE)
      map <- list(
        kind = "programs", xmean = xm, mu = dm, sigma = sdv,
        B = vapply(coef(fit), function(x) as.numeric(x)[-1], numeric(8)), training_cells = train$depmap_id
      )
    }
    maps[[context_key(dr, ds)]] <- map
  }
  maps
}
predict_rna <- function(model, cells, drug, dose = NULL) {
  if (is.null(dose)) dose <- bank_selection$dose[match(drug, drugs)]
  raw <- as.matrix(baseline[match(cells, depmap_id), ..programs])
  features <- sweep(sweep(sweep(raw, 2, model$bmean), 2, model$bsd, "/") %*% model$basis, 2, model$bscale, "/")
  g <- match(drug, drugs)
  D <- matrix(0, length(cells), 40)
  for (j in unique(g)) {
    rows <- which(g == j)
    Z <- sweep(features[rows, , drop = FALSE], 2, model$bm[j, ]) %*% model$B[, , j]
    D[rows, ] <- sweep(sweep(Z, 2, model$sigma[j, ], "*"), 2, model$mu[j, ], "+")
  }
  rows <- which(drug == "Lapatinib ditosylate")
  if (length(rows)) {
    lap <- model$lap_receptors
    stopifnot(!is.null(lap))
    X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
    Z <- sweep(sweep(X, 2, lap$xmean), 2, lap$xsd, "/") %*% lap$B
    D[rows, ] <- sweep(sweep(Z, 2, lap$sigma, "*"), 2, lap$mu, "+")
  }
  for (i in seq_len(nrow(extra_contexts))) {
    dr <- extra_contexts$drug[i]
    ds <- extra_contexts$dose[i]
    rows <- which(drug == dr & abs(dose - ds) < 1e-10)
    if (!length(rows)) next
    map <- model$extra_rna[[context_key(dr, ds)]]
    stopifnot(!is.null(map))
    if (map$kind == "receptors") {
      X <- as.matrix(receptors[match(cells[rows], depmap_id), .(ERBB2, MET)])
      X <- sweep(sweep(X, 2, map$xmean), 2, map$xsd, "/")
    } else {
      X <- sweep(features[rows, , drop = FALSE], 2, map$xmean)
    }
    D[rows, ] <- sweep(sweep(X %*% map$B, 2, map$sigma, "*"), 2, map$mu, "+")
  }
  stopifnot(all(is.finite(D)), !anyNA(context_index(drug, dose)))
  D
}

# The same three RNA-only footprint coordinates used by v100/v102 supply all
# shared calibration functions. Route labels are frozen training correlations.
make_context <- function(template, rna) {
  fm <- colMeans(template$mu)
  fs <- pmax(apply(template$mu, 2, sd), 1e-8)
  F <- sweep(sweep(template$mu, 2, fm), 2, fs, "/")
  fv <- svd(F, nu = 0, nv = 3)$v
  fscale <- sqrt(colMeans((F %*% fv)^2))
  Q <- cbind(1, sweep(F %*% fv, 2, fscale, "/"))
  mu <- template$mu
  sigma <- template$sigma
  for (i in seq_len(nrow(extra_contexts))) {
    map <- rna$extra_rna[[context_key(extra_contexts$drug[i], extra_contexts$dose[i])]]
    mu <- rbind(mu, map$mu)
    sigma <- rbind(sigma, map$sigma)
  }
  E <- sweep(sweep(mu[11:nrow(mu), , drop = FALSE], 2, fm), 2, fs, "/")
  Q <- rbind(Q, cbind(1, sweep(E %*% fv, 2, fscale, "/")))
  list(
    mu = mu, sigma = sigma, Q = Q, W = template$W * signature_amplitude,
    route = c(template$membership, template$membership[match(extra_contexts$drug, drugs)]),
    fm = fm, fs = fs, fv = fv, fscale = fscale
  )
}
signature_node <- function(D, g, ctx) {
  X <- (D - ctx$mu[g, , drop = FALSE]) / ctx$sigma[g, , drop = FALSE]
  both <- X %*% ctx$W
  route <- ctx$route[g]
  score <- both[cbind(seq_len(nrow(D)), route)]
  # Pair interactions reference zero displacement. Keep each constituent
  # separate and populate only its assigned route; no summed RNA displacement.
  absolute <- (D / ctx$sigma[g, , drop = FALSE]) %*% ctx$W
  activity <- matrix(0, nrow(D), 2)
  activity[cbind(seq_len(nrow(D)), route)] <- absolute[cbind(seq_len(nrow(D)), route)]
  list(score = score, both = both, activity = activity, route = route)
}
pair_features <- function(A, B) cbind(A[, 1] * B[, 1], A[, 1] * B[, 2] + A[, 2] * B[, 1], A[, 2] * B[, 2])
geometry <- function(y, group, control = NULL) {
  g <- match(group, unique(group))
  w <- rep(1, length(y))
  if (!is.null(control)) w <- 1 / as.numeric(table(control)[control])
  totals <- as.numeric(rowsum(w, g, reorder = FALSE))
  w <- w / totals[g] / length(totals)
  gw <- as.numeric(rowsum(w, g, reorder = FALSE))
  ym <- as.numeric(rowsum(y * w, g, reorder = FALSE)) / gw
  list(y = y, w = w, within = max(sum(w * (y - ym[g])^2), 1e-8))
}
loss <- function(pred, geo) {
  error <- pred - geo$y
  list(value = sum(geo$w * error^2) / geo$within, gradient = 2 * geo$w * error / geo$within)
}
unpack <- function(theta) list(level = theta[1:4], log_amplitude = theta[5:8], assay = theta[9:10], Gamma = theta[11:13])
pack <- function(par) unlist(par, use.names = FALSE)
prepare <- function(held, template, rna) {
  train <- profiles[!depmap_id %in% held]
  ns <- nci[!depmap_id %in% held]
  pp <- pairs[!depmap_id %in% held]
  query <- unique(rbind(ns[, .(depmap_id, drug, dose)],
    pp[, .(depmap_id, drug = drugA, dose = doseA)],
    pp[, .(depmap_id, drug = drugB, dose = doseB)]))
  qkey <- key(query$depmap_id, query$drug, query$dose)
  ctx <- make_context(template, rna)
  gobs <- context_index(train$drug, train$dose)
  gquery <- context_index(query$drug, query$dose)
  obs <- signature_node(as.matrix(train[, ..programs]), gobs, ctx)
  native <- signature_node(predict_rna(rna, query$depmap_id, query$drug, query$dose), gquery, ctx)
  is <- match(key(ns$depmap_id, ns$drug, ns$dose), qkey)
  ia <- match(key(pp$depmap_id, pp$drugA, pp$doseA), qkey)
  ib <- match(key(pp$depmap_id, pp$drugB, pp$doseB), qkey)
  products <- pair_features(native$activity[ia, , drop = FALSE], native$activity[ib, , drop = FALSE])
  ctx$interaction_scale <- sqrt(colMeans(products^2))
  ctx$interaction_supported <- ctx$interaction_scale > 1e-10
  # Unrepresented route pairs cannot identify an interaction coefficient.
  # Unit scaling for their zero training column permits an explicitly fixed 0.
  ctx$interaction_scale[!ctx$interaction_supported] <- 1
  stopifnot(any(ctx$interaction_supported), max(abs(tapply(obs$score, train$drug, mean))) < 1e-10)
  list(
    ctx = ctx, train = train, ns = ns, pp = pp, obs = obs, native = native,
    Qobs = ctx$Q[gobs, , drop = FALSE], Qquery = ctx$Q[gquery, , drop = FALSE], is = is, ia = ia, ib = ib,
    products = sweep(products, 2, ctx$interaction_scale, "/"),
    assay = match(ns$assay, c("nci_FF", "nci_FG")), pair_assay = match(pp$SCREENER, c("FF", "FG")),
    geoO = geometry(train$S, train$drug), geoN = geometry(-ns$y, ns$group, ns$control_key),
    geoP = geometry(-pp$y, pp$pair), geoI = geometry(-pp$I, pp$pair),
    training_cells = unique(c(train$depmap_id, ns$depmap_id, pp$depmap_id))
  )
}
fit_decoder <- function(p, template, job) {
  Qactive <- p$ctx$Q[included_indices, , drop = FALSE]
  initial_level <- as.numeric(solve(
    crossprod(Qactive) / nrow(Qactive) + diag(calibration_penalty, 4),
    crossprod(Qactive, template$reference[included_indices]) / nrow(Qactive)
  ))
  initial <- c(initial_level, rep(0, 4), -1, -1, rep(0, 3))
  cache_theta <- cache_value <- NULL
  objective <- function(theta) {
    if (identical(theta, cache_theta)) {
      return(cache_value)
    }
    par <- unpack(theta)
    ao <- as.numeric(exp(p$Qobs %*% par$log_amplitude))
    aq <- as.numeric(exp(p$Qquery %*% par$log_amplitude))
    po <- as.numeric(p$Qobs %*% par$level) + ao * p$obs$score
    native <- as.numeric(p$Qquery %*% par$level) + aq * p$native$score
    scales <- exp(par$assay)
    pn <- native[p$is] * scales[p$assay]
    # Positive amplitude changes interaction magnitude, not its route assignment.
    Qsum <- p$Qquery[p$ia, , drop = FALSE] + p$Qquery[p$ib, , drop = FALSE]
    calibrated_products <- p$products * (aq[p$ia] * aq[p$ib])
    calibrated_RMS <- sqrt(colMeans(calibrated_products^2))
    calibrated_RMS[!p$ctx$interaction_supported] <- 1
    F <- sweep(calibrated_products, 2, calibrated_RMS, "/")
    scale_gradient <- matrix(0, 3, 4)
    for (k in which(p$ctx$interaction_supported)) {
      scale_gradient[k, ] <- colMeans(Qsum * F[, k]^2)
    }
    J <- as.numeric(F %*% par$Gamma)
    pp <- (native[p$ia] + native[p$ib] + J) * scales[p$pair_assay]
    pi <- J * scales[p$pair_assay]
    lo <- loss(po, p$geoO)
    ln <- loss(pn, p$geoN)
    lp <- loss(pp, p$geoP)
    li <- loss(pi, p$geoI)
    go <- lo$gradient
    gn <- .2 * ln$gradient * scales[p$assay]
    gp <- .2 * lp$gradient * scales[p$pair_assay]
    gj <- .2 * (lp$gradient + li$gradient) * scales[p$pair_assay]
    QA <- p$Qquery[p$ia, , drop = FALSE]
    QB <- p$Qquery[p$ib, , drop = FALSE]
    QN <- p$Qquery[p$is, , drop = FALSE]
    level_gradient <- as.numeric(crossprod(p$Qobs, go) + crossprod(QN, gn) + crossprod(QA + QB, gp))
    amplitude_gradient <- as.numeric(crossprod(p$Qobs, go * ao * p$obs$score) +
      crossprod(QN, gn * aq[p$is] * p$native$score[p$is]) +
      crossprod(QA, gp * aq[p$ia] * p$native$score[p$ia]) +
      crossprod(QB, gp * aq[p$ib] * p$native$score[p$ib]) + crossprod(QA + QB, gj * J)) -
      as.numeric(crossprod(scale_gradient, as.numeric(crossprod(F, gj)) * par$Gamma))
    assay_gradient <- vapply(1:2, function(a) {
      .2 * sum(ln$gradient[p$assay == a] * pn[p$assay == a]) +
        .2 * sum(lp$gradient[p$pair_assay == a] * pp[p$pair_assay == a]) +
        .2 * sum(li$gradient[p$pair_assay == a] * pi[p$pair_assay == a])
    }, numeric(1))
    penalty <- calibration_penalty * (sum(par$level^2) + sum(par$log_amplitude^2)) +
      assay_penalty * sum(par$assay^2) + interaction_penalty * sum(par$Gamma^2)
    gradient <- pack(list(
      level = level_gradient + 2 * calibration_penalty * par$level,
      log_amplitude = amplitude_gradient + 2 * calibration_penalty * par$log_amplitude,
      assay = assay_gradient + 2 * assay_penalty * par$assay,
      Gamma = as.numeric(crossprod(F, gj)) + 2 * interaction_penalty * par$Gamma
    ))
    value <- list(value = lo$value + .2 * ln$value + .2 * lp$value + .2 * li$value + penalty, gradient = gradient)
    stopifnot(is.finite(value$value), all(is.finite(value$gradient)))
    cache_theta <<- theta
    cache_value <<- value
    value
  }
  if (job == 1L) {
    probe <- initial + c(.01, -.01, .02, .01, .03, -.02, .01, .02, .01, -.01, .01, -.02, .01)
    analytic <- objective(probe)$gradient
    numeric <- vapply(seq_along(probe), function(i) {
      plus <- minus <- probe
      plus[i] <- plus[i] + 1e-5
      minus[i] <- minus[i] - 1e-5
      (objective(plus)$value - objective(minus)$value) / 2e-5
    }, numeric(1))
    difference <- abs(analytic - numeric) / pmax(1, abs(analytic), abs(numeric))
    write_table(data.table(parameter = 1:13, analytic, numeric, relative_difference = difference), "analytic_gradient_check")
    stopifnot(max(difference) < 1e-5)
  }
  lower <- -c(rep(20, 4), rep(2, 4), rep(4, 2), rep(10, 3))
  upper <- -lower
  unsupported_indices <- 10L + which(!p$ctx$interaction_supported)
  lower[unsupported_indices] <- 0
  upper[unsupported_indices] <- 0
  opt <- optim(pmin(pmax(initial, lower), upper), function(t) objective(t)$value, function(t) objective(t)$gradient,
    method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 1000, factr = 1e8, pgtol = 1e-5)
  )
  final_par <- unpack(opt$par)
  final_aq <- as.numeric(exp(p$Qquery %*% final_par$log_amplitude))
  multiplier <- sqrt(colMeans((p$products * (final_aq[p$ia] * final_aq[p$ib]))^2))
  multiplier[!p$ctx$interaction_supported] <- 1
  p$ctx$interaction_scale <- p$ctx$interaction_scale * multiplier
  list(
    par = final_par, context = p$ctx, training_cells = p$training_cells,
    convergence = opt$convergence, evaluations = opt$counts[1], max_gradient = max(abs(objective(opt$par)$gradient))
  )
}

predict_single <- function(model, meta, measured_D = NULL) {
  dose <- if ("dose" %in% names(meta)) meta$dose else bank_selection$dose[match(meta$drug, drugs)]
  D <- if (is.null(measured_D)) predict_rna(model$rna, meta$depmap_id, meta$drug, dose) else measured_D
  g <- context_index(meta$drug, dose)
  node <- signature_node(D, g, model$context)
  Q <- model$context$Q[g, , drop = FALSE]
  level <- as.numeric(Q %*% model$par$level)
  amplitude <- as.numeric(exp(Q %*% model$par$log_amplitude))
  sensitivity <- level + amplitude * node$score
  assay <- if ("assay" %in% names(meta)) match(meta$assay, c("secondary", "nci_FF", "nci_FG")) else rep(1L, nrow(meta))
  scale <- c(1, exp(model$par$assay))[assay]
  stopifnot(all(amplitude > 0), !anyNA(scale))
  list(
    sensitivity = sensitivity * scale, score = node$score, route = node$route, level = level, amplitude = amplitude,
    activity = node$activity * amplitude, D = D
  )
}
predict_pair <- function(model, meta) {
  A <- predict_single(model, meta[, .(depmap_id, drug = drugA, dose = doseA)])
  B <- predict_single(model, meta[, .(depmap_id, drug = drugB, dose = doseB)])
  products <- sweep(pair_features(A$activity, B$activity), 2, model$context$interaction_scale, "/")
  J <- as.numeric(products %*% model$par$Gamma)
  scale <- exp(model$par$assay)[match(meta$SCREENER, c("FF", "FG"))]
  list(sensitivity = (A$sensitivity + B$sensitivity + J) * scale, departure = J * scale)
}
render_heatmap <- function(template, label) {
  plot_training_signatures(template, short, bank_selection$dose, label, included_indices)
}

measured_predictions <- baseline_predictions <- combo_predictions <- audits <- calibrations <- list()
direction_checks <- coefficient_tables <- interaction_tables <- list()
final_parent <- qread(file.path(parent, "models/final.qs"))
book <- file.path(root, "figures/training_correlations_and_collapsed_signatures_LOO.pdf")
cairo_pdf(book, width = 9, height = 9, onefile = TRUE)
book_device <- dev.cur()
final_job <- length(loo_cells) + 1L
status(paste("v106:", length(loo_cells), "whole-cell exclusions; fixed signatures/routes; shared positive calibration."))
for (job in seq_len(final_job)) {
  held <- if (job == final_job) character() else loo_cells[job]
  label <- if (length(held)) sprintf("loo_%02d", job) else "final"
  parent_index <- match(held, parent_folds$depmap_id)
  if (length(parent_index) && !is.na(parent_index)) {
    ancestor <- qread(file.path(parent, "models", sprintf("loo_%02d.qs", parent_index)))
    template <- ancestor$template
    rna <- ancestor$rna
  } else {
    template <- final_parent$template
    rna <- if (length(held) && held %in% molecular$depmap_id) fit_rna(molecular[depmap_id != held]) else final_parent$rna
  }
  stopifnot(!any(held %in% template$training_cells), !any(held %in% rna$training_cells))
  rna$lap_receptors <- fit_lapatinib_receptors(molecular[!depmap_id %in% held])
  rna$extra_rna <- fit_extra_rna(rna, held)
  stopifnot(!any(held %in% unlist(lapply(rna$extra_rna, `[[`, "training_cells"))))
  stopifnot(!any(held %in% rna$lap_receptors$training_cells))
  prepared <- prepare(held, template, rna)
  model_path <- file.path(root, "models", paste0(label, ".qs"))
  signature <- c(cache_signature, list(held = held, job = job))
  cached <- if (file.exists(model_path)) qread(model_path) else NULL
  reuse_fit <- !is.null(cached) && identical(cached$cache_signature, signature)
  if (reuse_fit) {
    status(paste(label, ": reusing matching fit"))
    fit <- cached[c("par", "context", "training_cells", "convergence", "evaluations", "max_gradient")]
  } else {
    status(paste(label, ": absent or stale fit; rebuilding"))
    fit <- fit_decoder(prepared, template, job)
  }
  model <- c(fit, list(template = template, rna = rna, held_cells = held, selected = selection,
    cache_signature = signature))
  stopifnot(
    !any(held %in% model$training_cells), identical(model$context$W, template$W * signature_amplitude),
    identical(model$context$route[1:10], template$membership)
  )
  if (!reuse_fit) {
    qsave(model, paste0(model_path, ".partial"))
    stopifnot(file.rename(paste0(model_path, ".partial"), model_path))
  }
  # Verify the agreed protection on every training treatment, not just parameter
  # signs: Pearson and centered score/prediction relationships must be preserved.
  training_result <- predict_single(model, prepared$train, as.matrix(prepared$train[, ..programs]))
  training_checks <- data.table(
    drug = prepared$train$drug, observed = prepared$train$S,
    score = training_result$score, predicted = training_result$sensitivity,
    level = training_result$level, amplitude = training_result$amplitude
  )[,
    .(
      score_r = pearson(observed, score), prediction_r = pearson(observed, predicted),
      score_prediction_r = pearson(score, predicted), level = unique(level), amplitude = unique(amplitude)
    ),
    by = drug
  ]
  stopifnot(
    max(abs(training_checks$score_r - training_checks$prediction_r)) < 1e-10,
    max(abs(training_checks$score_prediction_r - 1)) < 1e-10, all(training_checks$amplitude > 0)
  )
  direction_checks[[label]] <- training_checks[, fit := label]
  evaluation <- if (length(held)) "held_cell" else "fitted"
  test <- if (length(held)) profiles[depmap_id == held] else profiles
  if (nrow(test)) {
    result <- predict_single(model, test, as.matrix(test[, ..programs]))
    p <- test[, .(depmap_id, drug, dose, observed_log2 = -S)]
    p[, `:=`(
      predicted_log2 = -result$sensitivity, signature_score = result$score, route = result$route,
      calibration_level = result$level, positive_amplitude = result$amplitude, fold = job, evaluation = evaluation
    )]
    measured_predictions[[label]] <- p
  }
  test <- if (length(held)) prism[depmap_id == held] else prism
  if (nrow(test)) {
    result <- predict_single(model, test)
    p <- test[, .(depmap_id, drug, dose, observed_log2 = y)]
    p[, `:=`(
      predicted_log2 = -result$sensitivity, signature_score = result$score, route = result$route,
      calibration_level = result$level, positive_amplitude = result$amplitude, fold = job, evaluation = evaluation
    )]
    baseline_predictions[[label]] <- p
  }
  test <- if (length(held)) pairs[depmap_id == held] else pairs
  if (nrow(test)) {
    result <- predict_pair(model, test)
    p <- test[, .(depmap_id, pair, SCREENER, drugA, drugB, doseA, doseB, observed_log2 = y, observed_departure_log2 = I)]
    p[, `:=`(predicted_log2 = -result$sensitivity,
      predicted_departure_log2 = -result$departure, fold = job, evaluation = evaluation)]
    combo_predictions[[label]] <- p
  }
  calibration_level <- as.numeric(model$context$Q %*% model$par$level)
  calibration_amplitude <- as.numeric(exp(model$context$Q %*% model$par$log_amplitude))
  calibrations[[label]] <- data.table(
    fit = label, drug = contexts$drug, dose = contexts$dose, route = model$context$route,
    predicted_level = calibration_level, positive_amplitude = calibration_amplitude
  )
  coefficient_tables[[label]] <- data.table(
    fit = label, term = rep(c("intercept", "RNA_PC1", "RNA_PC2", "RNA_PC3"), 2),
    function_name = rep(c("absolute_level", "log_positive_amplitude"), each = 4),
    coefficient = c(model$par$level, model$par$log_amplitude)
  )
  interaction_tables[[label]] <- data.table(
    fit = label, term = c("route1_route1", "cross_route", "route2_route2"),
    training_supported = model$context$interaction_supported,
    product_RMS = model$context$interaction_scale, coefficient = model$par$Gamma
  )
  zeros <- signature_node(matrix(0, nrow(contexts), 40), seq_len(nrow(contexts)), model$context)$activity
  stopifnot(max(abs(zeros)) < 1e-12)
  audits[[label]] <- data.table(
    fit = label, held_cell = if (length(held)) held else NA_character_,
    convergence = fit$convergence, evaluations = fit$evaluations, max_gradient = fit$max_gradient,
    functional_leakage = length(intersect(held, model$training_cells)), RNA_leakage = length(intersect(held, rna$training_cells)),
    minimum_positive_amplitude = min(calibration_amplitude),
    maximum_correlation_change = max(abs(training_checks$score_r - training_checks$prediction_r))
  )
  dev.set(book_device)
  render_heatmap(template, paste(label, if (length(held)) held else "all-data fit"))
  if (job %in% c(1L, 30L, final_job)) {
    png(file.path(root, "figures", paste0(label, "_correlations_and_signatures.png")),
      width = 9, height = 9, units = "in", res = 150, type = "cairo")
    render_heatmap(template, label)
    dev.off()
  }
  write_table(rbindlist(audits), "fit_audit")
  if (job %% 10L == 0L || job == final_job) status(paste("Fits completed:", job, "of", final_job))
}
dev.set(book_device)
dev.off()
measured <- rbindlist(measured_predictions)
replay <- rbindlist(baseline_predictions)
combo <- rbindlist(combo_predictions)
external <- copy(replay[evaluation == "fitted" & !depmap_id %in% loo_cells])
stopifnot(!any(external$depmap_id %in% model$training_cells), !any(external$depmap_id %in% model$rna$training_cells))
external[, evaluation := "held_cell"]
replay <- rbind(replay, external)
stopifnot(
  nrow(measured[evaluation == "held_cell"]) == nrow(profiles),
  nrow(replay[evaluation == "held_cell"]) == nrow(prism), nrow(combo[evaluation == "held_cell"]) == nrow(pairs),
  !anyDuplicated(replay[evaluation == "held_cell", .(depmap_id, drug)])
)
write_table(measured, "measured_predictions_frozen")
write_table(replay, "baseline_predictions_frozen")
write_table(combo, "combination_predictions_frozen")
write_table(rbindlist(direction_checks), "signature_direction_protection_checks")
write_table(rbindlist(calibrations), "molecular_footprint_calibration_values")
write_table(rbindlist(coefficient_tables), "shared_calibration_coefficients")
write_table(rbindlist(interaction_tables), "interaction_feature_scales")

metrics <- measured[, .(
  cells = .N, r = pearson(observed_log2, predicted_log2),
  score_r = pearson(-observed_log2, signature_score), observed_sensitive = sum(observed_log2 <= -1),
  predicted_sensitive = sum(predicted_log2 <= -1), observed_median_pct = median(100 * 2^observed_log2),
  predicted_median_pct = median(100 * 2^predicted_log2)
), by = .(drug, dose, evaluation)]
bridge_metrics <- replay[, .(
  cells = .N, r = pearson(observed_log2, predicted_log2),
  observed_sensitive = sum(observed_log2 <= -1), predicted_sensitive = sum(predicted_log2 <= -1),
  observed_median_pct = median(100 * 2^observed_log2), predicted_median_pct = median(100 * 2^predicted_log2)
), by = .(drug, dose, evaluation)]
combo_metrics <- combo[, .(
  cells = .N, r_total = pearson(observed_log2, predicted_log2),
  r_departure = pearson(observed_departure_log2, predicted_departure_log2),
  observed_departure_SD = sd(observed_departure_log2), predicted_departure_SD = sd(predicted_departure_log2)
),
by = .(pair, SCREENER, evaluation)
]
write_table(metrics, "measured_decoder_correlations")
write_table(bridge_metrics, "baseline_bridge_correlations")
write_table(combo_metrics, "combination_correlations")

plot_single <- function(table, heading, stem) {
  data <- copy(table[evaluation == "held_cell"])
  data[, panel := factor(paste0(short[match(drug, drugs)], " | ", dose, " µM"),
    levels = paste0(short[included_indices], " | ", selection$dose, " µM")
  )]
  stats <- data[, .(label = correlation_label(observed_log2, predicted_log2)), by = panel]
  plot <- ggplot(data, aes(observed_log2, predicted_log2)) +
    geom_point(color = "#0072B2", alpha = .7, size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey30", linewidth = .5) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -.05, vjust = 1.1, size = 3.5) +
    facet_wrap(~panel, ncol = 5, scales = "free") +
    theme_bw(base_size = 12) +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank()) +
    labs(
      title = paste("v106 |", heading), x = "Measured absolute log2 abundance", y = "Predicted absolute log2 abundance",
      caption = paste0("Pearson across held-out cells at fixed drug/dose. P values are descriptive.\n",
        "0 = 100% abundance; -1 = 50%. Retrospective comparisons.\n",
        "Two fixed signatures; eight shared calibration coefficients; positive amplitude.")
    )
  ggsave(file.path(root, "figures", paste0(stem, ".pdf")), plot, width = 13, height = 6.8, device = cairo_pdf)
  ggsave(file.path(root, "figures", paste0(stem, ".png")), plot, width = 13, height = 6.8, dpi = 150)
}
plot_single(measured, "measured-displacement shared decoder", "held_measured_absolute_decoder")
plot_single(replay, "baseline -> displacement -> shared decoder", "held_baseline_absolute_bridge")
for (endpoint in c("total", "departure")) {
  data <- copy(combo[evaluation == "held_cell"])
  data[, observed := if (endpoint == "total") observed_log2 else observed_departure_log2]
  data[, predicted := if (endpoint == "total") predicted_log2 else predicted_departure_log2]
  stats <- data[, .(label = correlation_label(observed, predicted)), by = pair]
  plot <- ggplot(data, aes(observed, predicted)) +
    geom_point(color = "#0072B2", size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey30", linewidth = .5) +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -.05, vjust = 1.1, size = 3.5) +
    facet_wrap(~pair, ncol = 3, scales = "free") +
    theme_bw(base_size = 12) +
    theme(aspect.ratio = 1) +
    labs(
      title = paste("v106 held-cell ALMANAC", endpoint), x = paste("Measured", endpoint, "log2 abundance"),
      y = paste("Predicted", endpoint, "log2 abundance")
    )
  ggsave(file.path(root, "figures", paste0("held_combination_", endpoint, ".pdf")),
    plot, width = 10, height = 3 * ceiling(uniqueN(data$pair) / 3) + .7, device = cairo_pdf)
}

# Show calibration and direction separately on the principal lapatinib panel.
lap_data <- measured[evaluation == "held_cell" & drug == "Lapatinib ditosylate"]
bridge_lap <- replay[evaluation == "held_cell" & drug == "Lapatinib ditosylate"]
lap_plot_data <- rbind(
  lap_data[, .(observed_log2, predicted = signature_score, view = "Held raw signature")],
  lap_data[, .(observed_log2, predicted = -predicted_log2, view = "Held calibrated sensitivity")],
  bridge_lap[, .(observed_log2, predicted = -predicted_log2, view = "Held full baseline bridge")]
)
lap_plot_data[, view := factor(view,
  levels = c("Held raw signature", "Held calibrated sensitivity", "Held full baseline bridge"))]
lap_plot_data[, observed_S := -observed_log2]
lap_stats <- lap_plot_data[, .(label = correlation_label(observed_S, predicted)), by = view]
plot <- ggplot(lap_plot_data, aes(observed_S, predicted)) +
  geom_point(color = "#0072B2", alpha = .7, size = 1.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey30", linewidth = .5) +
  geom_text(data = lap_stats, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.05, vjust = 1.1, size = 3.5) +
  facet_wrap(~view, nrow = 1, scales = "free") +
  theme_bw(base_size = 12) +
  theme(aspect.ratio = 1) +
  labs(
    title = "v106 | Lapatinib 5 µM: fixed signature -> shared positive calibration -> baseline bridge",
    x = "Measured sensitivity (-log2 abundance)", y = "Score / predicted sensitivity"
  )
ggsave(file.path(root, "figures/lapatinib_signature_calibration_bridge.pdf"), plot, width = 10, height = 4.1, device = cairo_pdf)
ggsave(file.path(root, "figures/lapatinib_signature_calibration_bridge.png"), plot, width = 10, height = 4.1, dpi = 150)
write_table(data.table(assay = c("nci_FF", "nci_FG"), shared_endpoint_scale = exp(model$par$assay)), "assay_scale")
writeLines(c(
  "v106: fixed v95 correlation signatures and fold-training route assignments; shared positive affine calibration.",
  sprintf("%d drug/dose contexts; %d measured RNA/PRISM profiles, %d PRISM outcomes, %d ALMANAC pairs.", nrow(selection), nrow(profiles), nrow(prism), nrow(pairs)),
  sprintf("%d whole-cell exclusions across functional sources. AU565 excluded from RNA/functional fits.", length(loo_cells)),
  "Signature score: ((D - training RNA mean) / training RNA SD) times fixed W (v95 W x 0.25), selecting its training-derived route.",
  "Absolute sensitivity: Q beta + exp(Q eta) score. Q = intercept plus three RNA-only PCs of mean drug footprints.",
  "Beta and eta each have four shared coefficients; no drug identity, per-drug functional coefficient, or direct baseline input.",
  "No free gate, per-program outcome-fitted calibration gains/shifts, tanh, or quadratic readout.",
  "All treatment-level training Pearson correlations must agree with their uncalibrated score to 1e-10.",
  "Calibration values vary across whole-cell fits; identical within-fit ranking does not guarantee identical pooled LOO Pearson.",
  "Absolute endpoints fitted jointly: PRISM + .2 NCI single + .2 pair total + .2 pair departure, with original source weighting/geometry.",
  "L2 penalties: beta/eta .005, log NCI assay scales .01, interaction .05; fixed before evaluation.",
  "Two shared NCI assay scales and supported shared route-product interactions are fitted. Unsupported interaction coefficients remain fixed 0.",
  "Added five lapatinib pairs at 0.5 uM, with native partner molecular doses nearest original selected dose inside each cached outcome grid. No extrapolation; 192 added observations, 307 total.",
  "Extra contexts: lapatinib 0.5, azacytidine 0.5, vandetanib 0.5 uM. Original single-drug evaluation doses and frozen signatures/routes preserved.",
  "Interaction products are RMS-normalized AFTER positive single-effect amplitudes within each training fold; amplitude-dependent RMS derivatives included in fitting.",
  "Pair features use separate zero-referenced signature projections, multiplied by positive amplitude. No sum of displacement vectors.",
  "Raw zero displacement is not forced to 100% single-drug abundance in the mean-plus-centered-score model; zero-referenced pair interaction remains exactly 0 for an absent constituent.",
  "Lapatinib RNA predictor: ERBB2/HER2 + MET, standardized in RNA training cells; multivariate ridge lambda 1 predicts 40 standardized displacements, returned to measured units. Receptors never enter the functional decoder.",
  "Other baseline RNA maps unchanged and RNA-only. NCI single/pair training uses predicted RNA inputs.",
  paste("Included drugs:", paste(selection$drug, collapse = "; ")),
  paste("Excluded functional contexts:", paste(excluded_drugs, collapse = "; ")),
  "Fixed historical v95 signature bank and RNA footprint coordinates retained. Exclusion is an ablation of functional fitting, not a new signature-discovery analysis.",
  "HGF route is not identified by training data; both route assumptions are reported without selecting on Wilson labels. No HGF functional labels fitted.",
  "HGF receiver scales RNA-inferred/clamped from original ten footprints as before. HGF queries consistently contrast each constituent with zero displacement. Query lapatinib 5 uM; Wilson reference is 1 uM.",
  "Positive affine functional calibration is not demonstrated plate batch correction.",
  "Previously revealed labels/selected doses make all comparisons retrospective."
), file.path(root, "RUN_DETAILS.txt"))
capture.output(sessionInfo(), file = file.path(root, "SESSION_INFO.txt"))
input_paths <- c(
  archive_path, file.path(source_root, c(
    "selected_doses.tsv", "individual_training_outcomes.tsv",
    "combination_training_outcomes.tsv"
  )),
  file.path(parent, "models/final.qs"), "scripts/08_fit_single_drug_bridge.R"
)
write_table(data.table(path = input_paths, sha256 = vapply(input_paths, digest, character(1), file = TRUE, algo = "sha256")), "source_hashes")
print(metrics[evaluation == "held_cell"])
print(bridge_metrics[evaluation == "held_cell"])
print(combo_metrics[evaluation == "held_cell"])

status(paste("Complete. Solver failures:", sum(rbindlist(audits)$convergence != 0L)))
sink()
