#!/usr/bin/env Rscript
## Stage 38 mockup: can CCLE baseline state predict a frozen RAF-program
## displacement learned in Tahoe, and does that inferred displacement predict
## primary-screen vemurafenib viability in cell lines absent from Tahoe?
suppressPackageStartupMessages({library(data.table);library(matrixStats);library(ggplot2)})

root <- "/work/tahoe_prism_v2"
work <- dirname(root)
out <- file.path(root, "results/38_mockup_baseline_conditioned_braf_transfer")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
unlink(list.files(out, full.names = TRUE, recursive = TRUE), recursive = TRUE,
  force = TRUE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
check <- function(x, message) if (!isTRUE(x)) stop(message, call. = FALSE)
write_table <- function(x, name) fwrite(as.data.table(x), file.path(out, name),
  sep = "\t", na = "NA")
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))

## Reconstruct the exact Stage-36 response scores without changing its output.
read_npy_rows <- function(path, rows_zero_based) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  check(identical(magic, "\x93NUMPY") && version[1L] %in% 1:3,
    "Unsupported NPY file")
  len <- readBin(con, "integer", 1L,
    size = if (version[1L] == 1L) 2L else 4L, signed = FALSE,
    endian = "little")
  header <- rawToChar(readBin(con, "raw", len))
  shape <- as.integer(trimws(strsplit(sub(
    ".*'shape': \\(([^)]*)\\).*", "\\1", header), ",")[[1L]]))
  shape <- shape[!is.na(shape)]
  offset <- seek(con)
  answer <- matrix(NA_real_, length(rows_zero_based), shape[2L])
  for (i in seq_along(rows_zero_based)) {
    seek(con, offset + 4 * shape[2L] * rows_zero_based[i], origin = "start")
    answer[i, ] <- readBin(con, "numeric", n = shape[2L], size = 4L,
      endian = "little")
  }
  answer
}
normalize_pair <- function(treated, reference) {
  treated <- pmax(treated, 0); reference <- pmax(reference, 0)
  tp <- treated / sum(treated); rp <- reference / sum(reference)
  eligible <- tp * 10000 >= .1 & rp * 10000 >= .1
  m <- log2(tp[eligible] / rp[eligible])
  a <- -.5 * log2(tp[eligible] * rp[eligible])
  mq <- quantile(m, c(.3, .7)); aq <- quantile(a, c(.05, .95))
  kept <- m >= mq[1L] & m <= mq[2L] & a >= aq[1L] & a <= aq[2L]
  factor <- 2^median(m[kept])
  log2(1 + tp * 10000 / factor) - log2(1 + rp * 10000)
}

input_art <- file.path(root, "artifacts/26_braf_erbb_bypass")
prepared_art <- file.path(root, "artifacts/27_braf_erbb_hierarchy")
frozen <- readRDS(file.path(root,
  "artifacts/36_source_aware_braf_erbb_programs/frozen_source_aware_programs.rds"))
membership <- as.data.table(frozen$membership)
genes <- as.character(frozen$leaf_order)
signs <- frozen$gene_sign[genes]
programs <- c("P1", "P2", "P3", "P4")
program_drug <- c(P1 = "vemurafenib", P2 = "lapatinib",
  P3 = "cobimetinib", P4 = "vemurafenib")

a <- readRDS(file.path(input_art, "discovery_measurements.rds"))
idx <- match(genes, colnames(a$effects)); check(!anyNA(idx),
  "Frozen gene missing from Tahoe response matrix")
effects <- a$effects[, idx, drop = FALSE]; colnames(effects) <- genes
meta <- a$meta; cells <- a$cells
new_meta <- fread(file.path(prepared_art, "additional_response_metadata.tsv"))
nt <- read_npy_rows(file.path(prepared_art, "additional_treated_mean.npy"),
  new_meta$array_row)
nr <- read_npy_rows(file.path(prepared_art, "additional_reference_mean.npy"),
  new_meta$array_row)
ne <- t(vapply(seq_len(nrow(nt)), function(i) normalize_pair(nt[i, ], nr[i, ]),
  numeric(ncol(nt))))[, idx, drop = FALSE]
colnames(ne) <- genes
new_cells <- unique(new_meta[, .(depmap_id, cell_line, cell_name)])
model_meta <- fread(file.path(work,
  "outputs/dataset_inventory/cache/depmap/Model.csv"))
new_cells <- merge(new_cells, model_meta[, .(depmap_id = ModelID,
  tissue = OncotreeLineage, disease = OncotreePrimaryDisease)], by = "depmap_id")
new_cells[, `:=`(braf_v600_confirmed = TRUE, split = "Development")]
cells <- rbindlist(list(cells, new_cells), fill = TRUE)
cells[, analysis_included := !depmap_id %chin% c("ACH-000311", "ACH-000628")]
new_meta[, `:=`(compound = c("Vemurafenib" = "vemurafenib",
  "Lapatinib ditosylate" = "lapatinib", "Cobimetinib" = "cobimetinib")[drug],
  perturbation_id = paste0(drug, "__", dose, "uM"))]
new_meta[, observation_id := paste(depmap_id, perturbation_id, sep = "__")]
rownames(ne) <- new_meta$observation_id
effects <- rbind(effects, ne)
meta <- rbindlist(list(meta[, .(depmap_id, compound, drug, dose, dose_unit,
  perturbation_id, observation_id, n_plates)], new_meta[, .(depmap_id,
  compound, drug, dose, dose_unit, perturbation_id, observation_id,
  n_plates)]))
meta <- merge(meta, cells[, .(depmap_id, cell_name, tissue, split,
  braf_v600_confirmed, analysis_included)], by = "depmap_id")
setorder(meta, compound, depmap_id, dose)
effects <- effects[meta$observation_id, , drop = FALSE]
prism <- fread(file.path(prepared_art, "primary_prism_responses.tsv"))
prism[passed_str_profiling == FALSE | is.na(passed_str_profiling),
  prism_survival := NA_real_]
meta <- merge(meta, prism[, .(depmap_id, compound, prism_survival)],
  by = c("depmap_id", "compound"), all.x = TRUE)
setorder(meta, compound, depmap_id, dose)
bm <- fread(file.path(prepared_art, "plate_common_response_metadata.tsv"))
bg <- read_npy_rows(file.path(prepared_art, "plate_common_response.npy"),
  bm$array_row)[, idx, drop = FALSE]
pm <- rbindlist(list(fread(file.path(work,
  "cache/tahoe_prism_full_chain/tahoe_pdex_selected/pdex_plate_audit.tsv.gz"))[
    , .(depmap_id, drug, dose, plate)], fread(file.path(prepared_art,
  "additional_plate_audit.tsv"))[, .(depmap_id, drug, dose, plate)]))
for (i in which(meta$analysis_included)) {
  pp <- pm[depmap_id == meta$depmap_id[i] & drug == meta$drug[i] &
    dose == meta$dose[i], unique(plate)]
  jj <- which(bm$depmap_id == meta$depmap_id[i] & bm$plate %in% pp)
  effects[i, ] <- effects[i, ] - colMeans(bg[jj, , drop = FALSE])
}
meta[, dose_weight := c(.25, .5, .25)[match(dose, c(.05, .5, 5))]]
line_drug <- unique(meta[, .(depmap_id, compound, cell_name, tissue, split,
  braf_v600_confirmed, analysis_included, prism_survival)])
setorder(line_drug, compound, depmap_id)
line_drug[, key := paste(depmap_id, compound, sep = "__")]
auc <- matrix(NA_real_, nrow(line_drug), length(genes),
  dimnames = list(line_drug$key, genes))
for (i in seq_len(nrow(line_drug))) {
  ii <- which(meta$depmap_id == line_drug$depmap_id[i] &
    meta$compound == line_drug$compound[i])
  auc[i, ] <- colSums(effects[ii, , drop = FALSE] * meta$dose_weight[ii])
}
program_score <- function(rows, p) rowMeans(sweep(
  auc[rows, membership[program == p, gene], drop = FALSE], 2L,
  signs[membership[program == p, gene]], "*"))
targets <- rbindlist(lapply(programs, function(p) {
  ii <- which(line_drug$compound == program_drug[p])
  cbind(line_drug[ii], program = p, response_score = program_score(ii, p))
}))

## Baseline features from frozen program genes in CCLE 22Q2 expression.
expression_path <- file.path(work,
  "outputs/dataset_inventory/cache/depmap/CCLE_expression_full.csv")
expression_header <- names(fread(expression_path, nrows = 0L))
symbols <- sub(" \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$", "", expression_header,
  perl = TRUE)
select_index <- c(1L, which(symbols %chin% genes))
expression <- fread(expression_path, select = select_index)
setnames(expression, 1L, "depmap_id")
selected_symbols <- symbols[select_index[-1L]]
ex <- as.matrix(expression[, -1L]); storage.mode(ex) <- "double"
if (anyDuplicated(selected_symbols)) {
  ex <- t(rowsum(t(ex), selected_symbols, reorder = FALSE) /
    as.numeric(table(selected_symbols)[unique(selected_symbols)]))
  colnames(ex) <- unique(selected_symbols)
} else colnames(ex) <- selected_symbols
rownames(ex) <- expression$depmap_id
check(all(genes %chin% colnames(ex)), "CCLE expression lacks a frozen gene")
ex <- ex[, genes, drop = FALSE]

mutation_path <- file.path(work,
  "outputs/dataset_inventory/cache/depmap/CCLE_mutations_bool_hotspot.csv")
mh <- names(fread(mutation_path, nrows = 0L))
braf_column <- grep("^BRAF ", mh, value = TRUE)
check(length(braf_column) == 1L, "BRAF hotspot column is ambiguous")
mutation <- fread(mutation_path, select = c(1L, match(braf_column, mh)))
setnames(mutation, c("depmap_id", "braf_hotspot"))

train_ids <- unique(targets[split == "Development" & analysis_included == TRUE,
  depmap_id])
check(length(train_ids) == 28L, "Expected 28 Tahoe development lines")
gene_mean <- colMeans(ex[train_ids, , drop = FALSE])
gene_sd <- apply(ex[train_ids, , drop = FALSE], 2L, sd)
gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- 1
z <- sweep(sweep(ex, 2L, gene_mean, "-"), 2L, gene_sd, "/")
baseline_scores <- rbindlist(lapply(programs, function(p) data.table(
  depmap_id = rownames(z), feature = paste0("baseline_", p),
  value = rowMeans(sweep(z[, membership[program == p, gene], drop = FALSE],
    2L, signs[membership[program == p, gene]], "*")))))
baseline_scores <- dcast(baseline_scores, depmap_id ~ feature,
  value.var = "value")
baseline <- merge(baseline_scores, mutation, by = "depmap_id")

## A second, direct conditioning test in GSE68453.  Baseline and response use
## disjoint vehicle replicates: two vehicle arrays define baseline state and
## the third is subtracted from the PLX4720 mean.  Rotating the held vehicle
## replicate makes the result independent of an arbitrary replicate choice.
v33 <- readRDS(file.path(root,
  "artifacts/33_braf_erbb_programs/verification_inputs.rds"))
m684 <- as.data.table(v33$s684)
hgnc <- fread(file.path(work,
  "outputs/dataset_inventory/cache/cellfie/hgnc_complete_set.txt"))
hgnc <- hgnc[status == "Approved"]
lookup <- rbindlist(lapply(seq_len(nrow(hgnc)), function(i) {
  aliases <- unique(c(hgnc$symbol[i], unlist(strsplit(paste(
    hgnc$prev_symbol[i], hgnc$alias_symbol[i], sep = "|"), "\\|"))))
  data.table(input = aliases[nzchar(aliases) & !is.na(aliases) & aliases != "NA"],
    symbol = hgnc$symbol[i])
}))
lookup <- lookup[, if (uniqueN(symbol) == 1L) .(symbol = symbol[1L]) else
  .(symbol = NA_character_), by = input]
canonical <- function(x) {
  answer <- lookup$symbol[match(x, lookup$input)]
  exact <- x %chin% hgnc$symbol; answer[exact] <- x[exact]; answer
}
official <- canonical(genes)
check(!anyNA(official), "A frozen gene lacks an unambiguous HGNC symbol")
platform <- fread(file.path(input_art, "GSE68453_platform_annotation.tsv"),
  quote = "")
platform[, official_symbol := canonical(Symbol)]
gse_matrix <- suppressWarnings(fread(file.path(input_art,
  "GSE68453_series_matrix.txt.gz"), skip = '"ID_REF"', header = TRUE,
  check.names = FALSE, select = c("ID_REF", m684$gsm)))
gse_matrix <- gse_matrix[ID_REF %chin% platform$ID]
gse_x <- as.matrix(gse_matrix[, m684$gsm, with = FALSE])
storage.mode(gse_x) <- "double"
probe <- data.table(probe = gse_matrix$ID_REF,
  official_symbol = platform$official_symbol[match(gse_matrix$ID_REF,
    platform$ID)], iqr = rowIQRs(gse_x), array_row = seq_len(nrow(gse_x)))
probe <- probe[!is.na(official_symbol)]
setorder(probe, official_symbol, -iqr, probe)
probe <- probe[, head(.SD, 1L), by = official_symbol]
gse_x <- gse_x[probe$array_row[match(official, probe$official_symbol)],,
  drop = FALSE]
rownames(gse_x) <- genes; colnames(gse_x) <- m684$gsm
check(!anyNA(gse_x), "GSE68453 does not measure every frozen gene")
gse_scores <- sapply(programs, function(p) {
  pg <- membership[program == p, gene]
  colMeans(sweep(gse_x[pg,, drop = FALSE], 1L, signs[pg], "*"))
})
gse_long <- melt(data.table(m684[, .(gsm, cell_line, arm, replicate)],
  gse_scores), id.vars = c("gsm", "cell_line", "arm", "replicate"),
  variable.name = "program", value.name = "raw_score")
gse_long[, program := as.character(program)]
gse_lines <- sort(unique(m684$cell_line))
gse_vehicle <- sapply(gse_lines, function(line) rowMeans(gse_x[,
  m684$cell_line == line & m684$arm == "control", drop = FALSE]))
colnames(gse_vehicle) <- gse_lines
gse_gene_center <- rowMeans(gse_vehicle)
gse_gene_sd <- rowSds(gse_vehicle)
gse_gene_sd[!is.finite(gse_gene_sd) | gse_gene_sd == 0] <- 1
gse_vehicle_z <- sweep(sweep(gse_vehicle, 1L, gse_gene_center, "-"),
  1L, gse_gene_sd, "/")
gse_p3_genes <- membership[program == "P3", gene]
gse_baseline_transfer <- data.table(cell_line = gse_lines,
  gene_scaled_P3 = colMeans(sweep(gse_vehicle_z[gse_p3_genes,,
    drop = FALSE], 1L, signs[gse_p3_genes], "*")))
gse_program_p3 <- gse_long[arm == "control" & program == "P3",
  .(program_scaled_P3 = mean(raw_score)), by = cell_line]
gse_baseline_transfer <- merge(gse_baseline_transfer, gse_program_p3,
  by = "cell_line")
bliss <- data.table(cell_line = c("A375", "SKMEL24", "HT144", "DM331",
  "SKMEL28", "DM13", "VMM5A", "VMM18", "DM6", "VMM15", "VMM12",
  "VMM17"), bliss = c(5.7, 29.3, 16.3, 27.3, 49.7, 39.3, 43, 31.7,
  -2.7, -8.7, 3.7, -1))
gse_split <- rbindlist(lapply(1:3, function(held) {
  base <- gse_long[arm == "control" & replicate != held,
    .(baseline_score = mean(raw_score)), by = .(cell_line, program)]
  control <- gse_long[arm == "control" & replicate == held,
    .(held_control = raw_score), by = .(cell_line, program)]
  treated <- gse_long[arm == "plx", .(treated_score = mean(raw_score)),
    by = .(cell_line, program)]
  Reduce(function(x, y) merge(x, y, by = c("cell_line", "program")),
    list(base, control, treated))[, `:=`(held_vehicle_replicate = held,
      response_score = treated_score - held_control)]
}))
gse_split <- merge(gse_split, bliss, by = "cell_line")
loo_linear <- function(x, y) vapply(seq_along(y), function(i) {
  fit <- lm(outcome ~ predictor, data = data.frame(outcome = y[-i],
    predictor = x[-i]))
  predict(fit, newdata = data.frame(predictor = x[i]))
}, numeric(1L))
gse_split[, predicted_response := loo_linear(baseline_score, response_score),
  by = .(program, held_vehicle_replicate)]
gse_stats <- gse_split[, {
  test <- cor.test(baseline_score, response_score)
  .(n_lines = .N,
    baseline_response_pearson = unname(test$estimate),
    baseline_response_p_value = test$p.value,
    loo_response_pearson = safe_cor(predicted_response, response_score),
    loo_response_rmse = rmse(response_score, predicted_response),
    mean_reference_rmse = rmse(response_score, mean(response_score)),
    baseline_bliss_pearson = safe_cor(baseline_score, bliss))
}, by = .(program, held_vehicle_replicate)]
write_table(gse_split, "gse68453_independent_vehicle_splits.tsv")
write_table(gse_stats, "gse68453_independent_vehicle_stats.tsv")
gse_function <- gse_split[, .(baseline_score = mean(baseline_score),
  bliss = bliss[1L]), by = .(cell_line, program)]
gse_function[, predicted_bliss := loo_linear(baseline_score, bliss),
  by = program]
gse_function_stats <- gse_function[, {
  test <- cor.test(baseline_score, bliss)
  .(n_lines = .N, baseline_bliss_pearson = unname(test$estimate),
    association_p_value = test$p.value,
    loo_prediction_pearson = safe_cor(predicted_bliss, bliss),
    loo_prediction_rmse = rmse(bliss, predicted_bliss),
    mean_reference_rmse = rmse(bliss, mean(bliss)))
}, by = program]
write_table(gse_function, "gse68453_baseline_function_predictions.tsv")
write_table(gse_function_stats, "gse68453_baseline_function_stats.tsv")

## Closed-form ridge, with all scaling learned only from the fitting lines.
lambdas <- c(0, 10^seq(-4, 3, length.out = 36L))
fit_ridge <- function(table, outcome, features, lambda) {
  y <- table[[outcome]]
  if (!length(features)) return(list(mean = mean(y), features = features))
  x <- data.matrix(table[, ..features]); storage.mode(x) <- "double"
  center <- colMeans(x)
  scale <- apply(x, 2L, sd); scale[scale == 0 | !is.finite(scale)] <- 1
  xs <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  design <- cbind(1, xs)
  penalty <- diag(c(0, rep(lambda, ncol(xs))))
  beta <- solve(crossprod(design) + penalty, crossprod(design, y))
  list(beta = beta, center = center, scale = scale, features = features)
}
predict_ridge <- function(fit, table) {
  if (!length(fit$features)) return(rep(fit$mean, nrow(table)))
  feature_names <- fit$features
  x <- data.matrix(table[, ..feature_names]); storage.mode(x) <- "double"
  xs <- sweep(sweep(x, 2L, fit$center, "-"), 2L, fit$scale, "/")
  as.numeric(cbind(1, xs) %*% fit$beta)
}
select_ridge <- function(table, outcome, features) {
  if (!length(features)) return(list(lambda = Inf,
    fit = fit_ridge(table, outcome, features, 0), cv = data.table(
      lambda = Inf, cv_rmse = rmse(table[[outcome]], mean(table[[outcome]])))))
  cv <- rbindlist(lapply(lambdas, function(lambda) {
    pred <- vapply(seq_len(nrow(table)), function(i) predict_ridge(
      fit_ridge(table[-i], outcome, features, lambda), table[i]), numeric(1L))
    data.table(lambda = lambda, cv_rmse = rmse(table[[outcome]], pred))
  }))
  best <- cv[order(cv_rmse, -lambda)][1L]
  list(lambda = best$lambda,
    fit = fit_ridge(table, outcome, features, best$lambda), cv = cv)
}
feature_sets <- list(reference = character(), mutation = "braf_hotspot",
  expression = paste0("baseline_P", 1:4),
  expression_mutation = c(paste0("baseline_P", 1:4), "braf_hotspot"))

target_data <- merge(targets, baseline, by = "depmap_id")
target_data <- target_data[split == "Development" & analysis_included == TRUE &
  is.finite(response_score) & is.finite(prism_survival)]
write_table(target_data, "tahoe_baseline_response_training_data.tsv")
vem_rows <- which(line_drug$compound == "vemurafenib")
vem_response <- cbind(line_drug[vem_rows], as.data.table(sapply(programs,
  function(p) program_score(vem_rows, p))))
setnames(vem_response, programs, paste0("response_", programs))
vem_response <- merge(vem_response, baseline, by = "depmap_id")
vem_response <- vem_response[split == "Development" & analysis_included == TRUE &
  is.finite(prism_survival)]
write_table(vem_response, "tahoe_vemurafenib_multivariate_response.tsv")
response_fits <- list(); cv_predictions <- list(); selections <- list()
for (p in programs) for (model_name in names(feature_sets)) {
  d <- target_data[program == p]
  selected <- select_ridge(d, "response_score", feature_sets[[model_name]])
  key <- paste(p, model_name, sep = "__")
  response_fits[[key]] <- selected$fit
  if (!length(feature_sets[[model_name]])) {
    pred <- rep(mean(d$response_score), nrow(d))
  } else {
    pred <- vapply(seq_len(nrow(d)), function(i) {
      fold_fit <- fit_ridge(d[-i], "response_score",
        feature_sets[[model_name]], selected$lambda)
      predict_ridge(fold_fit, d[i])
    }, numeric(1L))
  }
  cv_predictions[[key]] <- d[, .(depmap_id, program, compound,
    braf_hotspot, observed_response = response_score)][,
      `:=`(model = model_name, predicted_response = pred)]
  selections[[key]] <- data.table(program = p, model = model_name,
    n_lines = nrow(d), selected_lambda = selected$lambda,
    training_cv_rmse = rmse(d$response_score, pred),
    training_cv_pearson = safe_cor(d$response_score, pred))
}
cv_predictions <- rbindlist(cv_predictions)
selections <- rbindlist(selections)
write_table(cv_predictions, "tahoe_crossvalidated_response_predictions.tsv")
write_table(selections, "response_model_selection.tsv")

## Freeze one response-to-function line per program using measured Tahoe scores.
decoders <- lapply(programs, function(p) lm(prism_survival ~ response_score,
  data = target_data[program == p]))
names(decoders) <- programs
decoder_table <- rbindlist(lapply(programs, function(p) data.table(program = p,
  term = names(coef(decoders[[p]])), coefficient = unname(coef(decoders[[p]])))))
write_table(decoder_table, "frozen_tahoe_function_decoders.tsv")

## Independent functional test: primary PRISM lines absent from all Tahoe data.
treatment <- fread(file.path(work,
  "data/sciPlex_PRISM/PRISM/primary-screen-replicate-collapsed-treatment-info.csv"))
treatment <- treatment[tolower(name) %chin% unique(unname(program_drug))]
prism_matrix_path <- file.path(work,
  "data/sciPlex_PRISM/PRISM/primary-screen-replicate-collapsed-logfold-change.csv")
ph <- names(fread(prism_matrix_path, nrows = 0L))
py <- fread(prism_matrix_path, select = c(1L,
  match(treatment$column_name, ph)))
setnames(py, c("row_name", as.character(treatment$name)))
cell_info <- fread(file.path(work,
  "data/sciPlex_PRISM/PRISM/primary-screen-cell-line-info.csv"))
py <- merge(py, cell_info[, .(row_name, depmap_id, primary_tissue,
  passed_str_profiling)], by = "row_name")
py <- merge(py, baseline, by = "depmap_id")
tahoe_ids <- unique(cells$depmap_id)
py <- py[passed_str_profiling == TRUE & !depmap_id %chin% tahoe_ids]
external <- rbindlist(lapply(programs, function(p) py[is.finite(get(program_drug[p])),
  .(depmap_id, primary_tissue, braf_hotspot,
    observed_survival = get(program_drug[p]), program = p,
    compound = program_drug[p], baseline_P1, baseline_P2, baseline_P3,
    baseline_P4)]))

external_predictions <- rbindlist(lapply(programs, function(p) {
  d <- external[program == p]
  rbindlist(lapply(names(feature_sets), function(model_name) {
    predicted_response <- predict_ridge(response_fits[[paste(p, model_name,
      sep = "__")]], d)
    predicted_survival <- as.numeric(predict(decoders[[p]],
      newdata = data.frame(response_score = predicted_response)))
    cbind(d, model = model_name, predicted_response, predicted_survival)
  }))
}))
write_table(external_predictions, "external_prism_predictions.tsv.gz")

scope_all <- copy(external_predictions); scope_all[, scope := "all"]
scope_melanoma <- copy(external_predictions[primary_tissue == "skin"])
scope_melanoma[, scope := "melanoma"]
scope_rows <- rbindlist(list(scope_all, scope_melanoma), use.names = TRUE,
  fill = TRUE)
performance <- scope_rows[, .(n_lines = uniqueN(depmap_id),
  pearson = safe_cor(observed_survival, predicted_survival),
  rmse = rmse(observed_survival, predicted_survival),
  observed_sd = sd(observed_survival), predicted_sd = sd(predicted_survival)),
  by = .(scope, program, compound, model)]
performance[, predicted_to_observed_sd := predicted_sd / observed_sd]
write_table(performance, "external_prism_performance.tsv")

external_p1 <- unique(external_predictions[program == "P1",
  .(depmap_id, primary_tissue, braf_hotspot, observed_survival, baseline_P3)])
genotype_rows <- rbindlist(list(
  external_p1[, .(scope = "all", depmap_id, braf_hotspot,
    observed_survival)],
  external_p1[primary_tissue == "skin", .(scope = "melanoma", depmap_id,
    braf_hotspot, observed_survival)]))
genotype <- genotype_rows[, .(n_lines = .N,
  mean_observed_survival = mean(observed_survival),
  sd_observed_survival = sd(observed_survival)), by = .(scope, braf_hotspot)]
genotype_contrast <- genotype_rows[, {
  test <- wilcox.test(observed_survival ~ braf_hotspot, exact = FALSE)
  .(n_lines = .N, n_wild_type = sum(braf_hotspot == 0),
    n_hotspot = sum(braf_hotspot == 1),
    hotspot_minus_wild_type_mean = mean(observed_survival[braf_hotspot == 1]) -
      mean(observed_survival[braf_hotspot == 0]),
    wilcoxon_p_value = test$p.value)
}, by = scope]
write_table(genotype, "braf_genotype_summary.tsv")
write_table(genotype_contrast, "braf_genotype_contrast.tsv")

## Aim 1 pilot: direct held-line prediction from baseline state.  P3 is the
## frozen shared MEK-ERK output program and BRAF hotspot is the genomic feature.
## The model formula is fixed and each line is excluded from coefficient fitting.
loo_model <- function(table, features) vapply(seq_len(nrow(table)), function(i) {
  train <- table[-i]
  if (!length(features)) return(mean(train$observed_survival))
  fit <- lm(reformulate(features, response = "observed_survival"), data = train)
  predict(fit, newdata = table[i])
}, numeric(1L))
ccle_models <- list(mean = character(), mutation = "braf_hotspot",
  expression = "baseline_P3",
  expression_mutation = c("baseline_P3", "braf_hotspot"))
ccle_cv <- rbindlist(lapply(c("all", "melanoma"), function(scope_name) {
  d <- if (scope_name == "all") external_p1 else
    external_p1[primary_tissue == "skin"]
  rbindlist(lapply(names(ccle_models), function(model_name) {
    pred <- loo_model(d, ccle_models[[model_name]])
    d[, .(scope = scope_name, model = model_name, depmap_id, primary_tissue,
      braf_hotspot, baseline_P3, observed_survival,
      predicted_survival = pred)]
  }))
}))
ccle_cv_stats <- ccle_cv[, {
  association <- if (model == "mean") NULL else
    cor.test(predicted_survival, observed_survival)
  .(n_lines = .N,
    pearson = if (is.null(association)) NA_real_ else
      unname(association$estimate),
    correlation_p_value = if (is.null(association)) NA_real_ else
      association$p.value,
    correlation_ci_low = if (is.null(association)) NA_real_ else
      association$conf.int[1L],
    correlation_ci_high = if (is.null(association)) NA_real_ else
      association$conf.int[2L],
    rmse = rmse(observed_survival, predicted_survival),
    prediction_sd = sd(predicted_survival))
}, by = .(scope, model)]
ccle_reference <- ccle_cv_stats[model == "mean", .(scope,
  reference_rmse = rmse)]
ccle_cv_stats <- merge(ccle_cv_stats, ccle_reference, by = "scope")
ccle_cv_stats[, rmse_improvement := reference_rmse - rmse]
write_table(ccle_cv, "ccle_prism_leave_one_line_out_predictions.tsv")
write_table(ccle_cv_stats, "ccle_prism_leave_one_line_out_stats.tsv")
ccle_melanoma <- external_p1[primary_tissue == "skin"]
ccle_fit <- lm(observed_survival ~ baseline_P3 + braf_hotspot,
  data = ccle_melanoma)
ccle_coefficient <- as.data.table(coef(summary(ccle_fit)), keep.rownames = "term")
setnames(ccle_coefficient, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
  c("estimate", "standard_error", "t_value", "p_value"))
ccle_ci <- as.data.table(confint(ccle_fit), keep.rownames = "term")
setnames(ccle_ci, c("2.5 %", "97.5 %"), c("ci_low", "ci_high"))
ccle_coefficient <- merge(ccle_coefficient, ccle_ci, by = "term")
write_table(ccle_coefficient, "ccle_prism_full_model_coefficients.tsv")

## Mechanistic Aim 1 transfer.  Tahoe alone supplies both fitted stages:
## baseline P3 plus BRAF status -> vemurafenib program displacement, followed by
## measured program displacement -> viability.  CCLE/PRISM is never refitted.
condition_programs <- c("P1", "P3", "P4")
condition_fits <- lapply(condition_programs, function(p) lm(reformulate(
  c("baseline_P3", "braf_hotspot"), response = paste0("response_", p)),
  data = vem_response))
names(condition_fits) <- condition_programs
conditioning <- rbindlist(lapply(condition_programs, function(p) {
  outcome <- paste0("response_", p)
  loo <- vapply(seq_len(nrow(vem_response)), function(i) predict(lm(
    reformulate(c("baseline_P3", "braf_hotspot"), response = outcome),
    data = vem_response[-i]), newdata = vem_response[i]), numeric(1L))
  fitted <- predict(condition_fits[[p]], newdata = vem_response)
  data.table(program = p, n_lines = nrow(vem_response),
    fitted_pearson = safe_cor(fitted, vem_response[[outcome]]),
    loo_pearson = safe_cor(loo, vem_response[[outcome]]),
    loo_rmse = rmse(vem_response[[outcome]], loo))
}))
condition_coefficients <- rbindlist(lapply(condition_programs, function(p) {
  x <- as.data.table(coef(summary(condition_fits[[p]])), keep.rownames = "term")
  setnames(x, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "standard_error", "t_value", "p_value"))
  x[, program := p][]
}))
write_table(conditioning, "tahoe_conditioning_model_stats.tsv")
write_table(condition_coefficients, "tahoe_conditioning_model_coefficients.tsv")

decoder_sets <- list(P1 = "response_P1",
  P1_P3 = c("response_P1", "response_P3"),
  P1_P4 = c("response_P1", "response_P4"))
p1_decoder <- lm(prism_survival ~ response_P1, data = vem_response)
transfer_predictions <- rbindlist(lapply(names(decoder_sets), function(name) {
  response_features <- decoder_sets[[name]]
  decoder <- lm(reformulate(response_features, response = "prism_survival"),
    data = vem_response)
  inferred <- as.data.table(lapply(sub("response_", "", response_features),
    function(p) predict(condition_fits[[p]], newdata = external_p1)))
  setnames(inferred, response_features)
  predicted <- predict(decoder, newdata = inferred)
  reference_values <- as.data.frame(as.list(colMeans(
    vem_response[, ..response_features])))
  names(reference_values) <- response_features
  reference <- rep(predict(decoder, newdata = reference_values),
    nrow(external_p1))
  cbind(external_p1, transfer_model = name, inferred,
    predicted_survival = predicted, reference_survival = reference)
}), fill = TRUE)
transfer_scope <- rbindlist(list(
  copy(transfer_predictions)[, scope := "all"],
  copy(transfer_predictions[primary_tissue == "skin"])[, scope := "melanoma"]))
transfer_stats <- transfer_scope[, {
  test <- cor.test(predicted_survival, observed_survival)
  .(n_lines = .N, pearson = unname(test$estimate), p_value = test$p.value,
    ci_low = test$conf.int[1L], ci_high = test$conf.int[2L],
    rmse = rmse(observed_survival, predicted_survival),
    transferred_mean_rmse = rmse(observed_survival, reference_survival),
    predicted_sd = sd(predicted_survival), observed_sd = sd(observed_survival))
}, by = .(scope, transfer_model)]
decoder_coefficients <- rbindlist(lapply(names(decoder_sets), function(name) {
  fit <- lm(reformulate(decoder_sets[[name]], response = "prism_survival"),
    data = vem_response)
  x <- as.data.table(coef(summary(fit)), keep.rownames = "term")
  setnames(x, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "standard_error", "t_value", "p_value"))
  x[, transfer_model := name][]
}))
write_table(transfer_predictions, "tahoe_to_ccle_transfer_predictions.tsv")
write_table(transfer_stats, "tahoe_to_ccle_transfer_stats.tsv")
write_table(decoder_coefficients, "tahoe_viability_decoder_comparison.tsv")

## Apply the same frozen Tahoe chain to the independent GSE68453 functional
## cohort.  Cohort scaling uses vehicle expression only; all lines are BRAF+
## and neither treated expression nor Bliss enters either prediction stage.
gse_transfer <- melt(gse_baseline_transfer, id.vars = "cell_line",
  variable.name = "normalization", value.name = "external_baseline_P3")
gse_transfer <- merge(gse_transfer, bliss, by = "cell_line")
gse_transfer[, baseline_P3 := (external_baseline_P3 -
  mean(external_baseline_P3)) / sd(external_baseline_P3) *
  sd(vem_response$baseline_P3) + mean(vem_response$baseline_P3),
  by = normalization]
gse_transfer[, braf_hotspot := 1L]
gse_transfer[, inferred_P1_displacement := predict(condition_fits$P1,
  newdata = gse_transfer)]
gse_transfer[, predicted_residual_viability := predict(p1_decoder,
  newdata = data.frame(response_P1 = inferred_P1_displacement))]
gse_transfer_stats <- gse_transfer[, {
  test <- cor.test(predicted_residual_viability, bliss)
  .(n_lines = .N, pearson = unname(test$estimate), p_value = test$p.value,
    ci_low = test$conf.int[1L], ci_high = test$conf.int[2L],
    spearman = cor(predicted_residual_viability, bliss, method = "spearman"))
}, by = normalization]
write_table(gse_transfer, "tahoe_to_gse68453_function_predictions.tsv")
write_table(gse_transfer_stats, "tahoe_to_gse68453_function_stats.tsv")

## Use the factorial design to identify the actual bypass phenotype: residual
## P3 during PLX4720 and the extra P3 reduction produced by adding lapatinib.
gse_factorial <- fread(file.path(root,
  "results/37_source_aware_bypass_validation/gse68453_program_factorial.tsv"))
gse_p3_factorial <- gse_factorial[program == "P3",
  .(cell_line, plx_P3_response = plx,
    lapatinib_extra_P3_reduction = -interaction)]
gse_bypass <- merge(gse_transfer, gse_p3_factorial, by = "cell_line")
gse_bypass_stats <- gse_bypass[, {
  plx_test <- cor.test(predicted_residual_viability, plx_P3_response)
  function_test <- cor.test(predicted_residual_viability, bliss)
  factorial_ok <- is.finite(lapatinib_extra_P3_reduction)
  extra_test <- cor.test(predicted_residual_viability[factorial_ok],
    lapatinib_extra_P3_reduction[factorial_ok])
  bridge_test <- cor.test(lapatinib_extra_P3_reduction[factorial_ok],
    bliss[factorial_ok])
  .(n_plx = .N, predicted_to_plx_P3 = unname(plx_test$estimate),
    predicted_to_plx_P3_p = plx_test$p.value,
    n_factorial = sum(factorial_ok),
    predicted_to_extra_P3 = unname(extra_test$estimate),
    predicted_to_extra_P3_p = extra_test$p.value,
    predicted_to_bliss = unname(function_test$estimate),
    predicted_to_bliss_p = function_test$p.value,
    extra_P3_to_bliss = unname(bridge_test$estimate),
    extra_P3_to_bliss_p = bridge_test$p.value)
}, by = normalization]

gse_split_transfer <- gse_split[program == "P3", .(cell_line,
  held_vehicle_replicate, external_baseline_P3 = baseline_score,
  held_control_P3_response = response_score)]
gse_split_transfer[, baseline_P3 := (external_baseline_P3 -
  mean(external_baseline_P3)) / sd(external_baseline_P3) *
  sd(vem_response$baseline_P3) + mean(vem_response$baseline_P3),
  by = held_vehicle_replicate]
gse_split_transfer[, braf_hotspot := 1L]
gse_split_transfer[, inferred_P1_displacement := predict(condition_fits$P1,
  newdata = gse_split_transfer)]
gse_split_transfer[, predicted_residual_viability := predict(p1_decoder,
  newdata = data.frame(response_P1 = inferred_P1_displacement))]
gse_split_bypass_stats <- gse_split_transfer[, .(n_lines = .N,
  predicted_to_held_control_P3 = safe_cor(predicted_residual_viability,
    held_control_P3_response)), by = held_vehicle_replicate]
write_table(gse_bypass, "gse68453_baseline_to_bypass_predictions.tsv")
write_table(gse_bypass_stats, "gse68453_baseline_to_bypass_stats.tsv")
write_table(gse_split_transfer,
  "gse68453_baseline_to_bypass_independent_vehicle_splits.tsv")
write_table(gse_split_bypass_stats,
  "gse68453_baseline_to_bypass_independent_vehicle_stats.tsv")

## Compact grant-oriented figures.
colors <- c(`0` = "#87949C", `1` = "#D88736")
ccle_mel <- ccle_cv[scope == "melanoma" & model == "expression_mutation"]
ccle_mel_stats <- ccle_cv_stats[scope == "melanoma" &
  model == "expression_mutation"]
transfer_mel <- transfer_predictions[primary_tissue == "skin" &
  transfer_model == "P1"]
transfer_mel_stats <- transfer_stats[scope == "melanoma" &
  transfer_model == "P1"]
p <- ggplot(transfer_mel, aes(predicted_survival, observed_survival,
  color = factor(braf_hotspot))) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .8) +
  scale_color_manual(values = colors, name = "BRAF hotspot",
    labels = c("Absent", "Present")) +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1, legend.position = "bottom",
    plot.title = element_text(size = 12.5)) +
  labs(x = "Predicted via inferred P1 displacement",
    y = "Observed primary-PRISM viability score",
    title = "Tahoe conditioning coefficients transfer to CCLE",
    subtitle = sprintf("38 non-Tahoe lines; no refitting, Pearson r = %.2f",
      transfer_mel_stats$pearson))
ggsave(file.path(out, "figures/01_tahoe_to_ccle_p1_transfer.pdf"), p,
  width = 4.7, height = 4.7, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/01_tahoe_to_ccle_p1_transfer.png"), p,
  width = 4.7, height = 4.7, dpi = 240, bg = "white")

mel <- external_p1[primary_tissue == "skin"]
p <- ggplot(mel, aes(factor(braf_hotspot), observed_survival,
  color = factor(braf_hotspot))) + geom_boxplot(width = .48, outlier.shape = NA,
  color = "#263840") + geom_jitter(width = .10, height = 0, size = 2.0) +
  scale_color_manual(values = colors, guide = "none") +
  scale_x_discrete(labels = c("BRAF wild-type", "BRAF hotspot")) +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1,
    plot.title = element_text(size = 12.5)) +
  labs(x = NULL, y = "Observed primary-PRISM viability score",
    title = "BRAF genotype predicts vemurafenib response",
    subtitle = sprintf("%d wild-type and %d hotspot lines; Wilcoxon p = %.3g",
      mel[braf_hotspot == 0, .N], mel[braf_hotspot == 1, .N],
      genotype_contrast[scope == "melanoma", wilcoxon_p_value]))
ggsave(file.path(out, "figures/02_ccle_prism_braf_control.pdf"), p,
  width = 4.6, height = 4.35, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/02_ccle_prism_braf_control.png"), p,
  width = 4.6, height = 4.35, dpi = 240, bg = "white")

p <- ggplot(ccle_mel, aes(predicted_survival, observed_survival,
  color = factor(braf_hotspot))) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .8) +
  scale_color_manual(values = colors, name = "BRAF hotspot",
    labels = c("Absent", "Present")) +
  theme_classic(base_size = 11) +
  theme(aspect.ratio = 1, legend.position = "bottom",
    plot.title = element_text(size = 12.5)) +
  labs(x = "Leave-one-line-out predicted viability score",
    y = "Observed primary-PRISM viability score",
    title = "Direct CCLE fit gives the same held-line result",
    subtitle = sprintf("Baseline P3 + BRAF; LOO r = %.2f, RMSE = %.2f",
      ccle_mel_stats$pearson, ccle_mel_stats$rmse))
ggsave(file.path(out, "figures/03_ccle_prism_direct_sensitivity.pdf"), p,
  width = 4.7, height = 4.7, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/03_ccle_prism_direct_sensitivity.png"), p,
  width = 4.7, height = 4.7, dpi = 240, bg = "white")

p3_split <- gse_split[program == "P3"]
p3_labels <- gse_stats[program == "P3", .(held_vehicle_replicate,
  label = sprintf("r = %.2f", baseline_response_pearson))]
p <- ggplot(p3_split, aes(baseline_score, response_score)) +
  geom_point(size = 2.1, color = "#24745E") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .8) +
  geom_text(data = p3_labels, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.15, vjust = 1.25, size = 3.2) +
  facet_wrap(~held_vehicle_replicate, nrow = 1,
    labeller = labeller(held_vehicle_replicate = function(x)
      paste("Held vehicle", x))) +
  theme_classic(base_size = 10.5) + theme(aspect.ratio = 1,
    plot.title = element_text(size = 12.5)) +
  labs(x = "P3 score in the other two vehicle replicates",
    y = "PLX4720 P3 change vs held vehicle",
    title = "Pretreatment P3 state conditions the measured drug response",
    subtitle = "Twelve BRAF-mutant melanoma lines; disjoint baseline and response controls")
ggsave(file.path(out, "figures/04_gse68453_baseline_conditioned_p3_response.pdf"),
  p, width = 7.4, height = 3.25, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/04_gse68453_baseline_conditioned_p3_response.png"),
  p, width = 7.4, height = 3.25, dpi = 240, bg = "white")

bypass_plot <- gse_bypass[normalization == "program_scaled_P3"]
bypass_stat <- gse_bypass_stats[normalization == "program_scaled_P3"]
bypass_theme <- theme_classic(base_size = 9.5) + theme(aspect.ratio = 1,
  plot.title = element_text(size = 10.5, face = "bold"))
pa <- ggplot(bypass_plot, aes(predicted_residual_viability, plx_P3_response)) +
  geom_point(size = 2.2, color = "#24745E") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) + bypass_theme +
  labs(x = "Baseline-predicted\nbypass propensity",
    y = "P3 response during PLX4720",
    title = sprintf("Residual P3, n=12; r=%.2f",
      bypass_stat$predicted_to_plx_P3))
pb <- ggplot(bypass_plot[is.finite(lapatinib_extra_P3_reduction)],
  aes(predicted_residual_viability, lapatinib_extra_P3_reduction)) +
  geom_point(size = 2.2, color = "#276C9B") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) + bypass_theme +
  labs(x = "Baseline-predicted\nbypass propensity",
    y = "Extra P3 reduction with lapatinib",
    title = sprintf("Factorial RNA, n=5; r=%.2f",
      bypass_stat$predicted_to_extra_P3))
pc <- ggplot(bypass_plot, aes(predicted_residual_viability, bliss)) +
  geom_point(size = 2.2, color = "#5B4775") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .75) + bypass_theme +
  labs(x = "Baseline-predicted\nbypass propensity",
    y = "PLX4720 + lapatinib Bliss benefit",
    title = sprintf("Functional benefit, n=12; r=%.2f",
      bypass_stat$predicted_to_bliss))
p <- ggpubr::ggarrange(pa, pb, pc, nrow = 1, align = "hv")
p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(
  "Pretreatment state predicts molecular and functional bypass",
  face = "bold", size = 13))
ggsave(file.path(out, "figures/05_gse68453_baseline_to_bypass_chain.pdf"),
  p, width = 8.4, height = 3.45, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/05_gse68453_baseline_to_bypass_chain.png"),
  p, width = 8.4, height = 3.45, dpi = 240, bg = "white")

gse_transfer_plot <- copy(gse_transfer)
gse_transfer_plot[, normalization := factor(normalization,
  levels = c("gene_scaled_P3", "program_scaled_P3"),
  labels = c("Gene-wise scaling", "Program-score scaling"))]
gse_transfer_labels <- copy(gse_transfer_stats)
gse_transfer_labels[, normalization := factor(normalization,
  levels = c("gene_scaled_P3", "program_scaled_P3"),
  labels = c("Gene-wise scaling", "Program-score scaling"))]
gse_transfer_labels[, label := sprintf("r = %.2f", pearson)]
p <- ggplot(gse_transfer_plot, aes(predicted_residual_viability, bliss)) +
  geom_point(size = 2.4, color = "#5B4775") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .8) +
  geom_text(data = gse_transfer_labels, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -.15, vjust = 1.25, size = 3.2) +
  facet_wrap(~normalization, nrow = 1) + theme_classic(base_size = 10.5) +
  theme(aspect.ratio = 1, plot.title = element_text(size = 12.5)) +
  labs(x = "Predicted residual viability from Tahoe model",
    y = "Observed PLX4720 + lapatinib Bliss benefit",
    title = "Baseline-conditioned Tahoe model predicts external function",
    subtitle = "Twelve BRAF-mutant melanoma lines; no functional outcome used in prediction")
ggsave(file.path(out, "figures/06_tahoe_to_gse68453_function_transfer.pdf"),
  p, width = 6.8, height = 3.65, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/06_tahoe_to_gse68453_function_transfer.png"),
  p, width = 6.8, height = 3.65, dpi = 240, bg = "white")

p3_function <- gse_function[program == "P3"]
p3_fs <- gse_function_stats[program == "P3"]
p <- ggplot(p3_function, aes(predicted_bliss, bliss)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
    color = "#9AA4A9", linewidth = .6) +
  geom_point(size = 2.5, color = "#5B4775") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
    color = "#263840", linewidth = .8) +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1,
    plot.title = element_text(size = 12.5)) +
  labs(x = "Leave-one-line-out predicted Bliss benefit",
    y = "Observed PLX4720 + lapatinib Bliss benefit",
    title = "Pretreatment P3 predicts function in a held-out line",
    subtitle = sprintf("12 lines, Pearson r = %.2f; RMSE %.1f vs %.1f for mean",
      p3_fs$loo_prediction_pearson, p3_fs$loo_prediction_rmse,
      p3_fs$mean_reference_rmse))
ggsave(file.path(out, "figures/07_gse68453_heldout_function_prediction.pdf"),
  p, width = 4.8, height = 4.45, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "figures/07_gse68453_heldout_function_prediction.png"),
  p, width = 4.8, height = 4.45, dpi = 240, bg = "white")

summary <- data.table(endpoint = c(
  "Tahoe-to-CCLE inferred P1 transfer: melanoma",
  "Tahoe-to-CCLE inferred P1 transfer: all lineages",
  "Direct CCLE/PRISM melanoma P3 plus BRAF sensitivity analysis",
  "CCLE/PRISM melanoma BRAF-only held-line prediction",
  "Non-Tahoe melanoma BRAF hotspot viability contrast",
  "Tahoe-to-GSE68453 functional transfer: gene-wise scaling",
  "Tahoe-to-GSE68453 functional transfer: program-score scaling",
  "GSE68453 P3 baseline-response association: weakest vehicle split",
  "GSE68453 P3 baseline-response association: strongest vehicle split",
  "GSE68453 P3 leave-one-line-out functional prediction",
  "Tahoe P1 baseline-to-response audit"),
  n = c(transfer_mel_stats$n_lines,
    transfer_stats[scope == "all" & transfer_model == "P1", n_lines],
    ccle_mel_stats$n_lines,
    ccle_cv_stats[scope == "melanoma" & model == "mutation", n_lines],
    genotype_contrast[scope == "melanoma", n_lines],
    gse_transfer_stats[normalization == "gene_scaled_P3", n_lines],
    gse_transfer_stats[normalization == "program_scaled_P3", n_lines],
    12L, 12L, p3_fs$n_lines,
    selections[program == "P1" & model == "expression_mutation", n_lines]),
  primary_value = c(transfer_mel_stats$pearson,
    transfer_stats[scope == "all" & transfer_model == "P1", pearson],
    ccle_mel_stats$pearson,
    ccle_cv_stats[scope == "melanoma" & model == "mutation", pearson],
    genotype_contrast[scope == "melanoma", hotspot_minus_wild_type_mean],
    gse_transfer_stats[normalization == "gene_scaled_P3", pearson],
    gse_transfer_stats[normalization == "program_scaled_P3", pearson],
    max(gse_stats[program == "P3",
    baseline_response_pearson]), min(gse_stats[program == "P3",
    baseline_response_pearson]), p3_fs$loo_prediction_pearson,
    selections[program == "P1" & model == "expression_mutation",
      training_cv_pearson]),
  metric = c("External Pearson r", "External Pearson r", "LOO Pearson r",
    "LOO Pearson r",
    "hotspot minus wild-type mean viability", "External Pearson r",
    "External Pearson r", "Pearson r", "Pearson r",
    "LOO Pearson r", "LOO Pearson r"),
  model_rmse = c(transfer_mel_stats$rmse,
    transfer_stats[scope == "all" & transfer_model == "P1", rmse],
    ccle_mel_stats$rmse,
    ccle_cv_stats[scope == "melanoma" & model == "mutation", rmse],
    NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
    p3_fs$loo_prediction_rmse,
    selections[program == "P1" & model == "expression_mutation",
      training_cv_rmse]),
  reference_rmse = c(transfer_mel_stats$transferred_mean_rmse,
    transfer_stats[scope == "all" & transfer_model == "P1",
      transferred_mean_rmse],
    ccle_mel_stats$reference_rmse,
    ccle_cv_stats[scope == "melanoma" & model == "mutation", reference_rmse],
    NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
    p3_fs$mean_reference_rmse,
    selections[program == "P1" & model == "reference", training_cv_rmse]))
write_table(summary, "mockup_summary.tsv")

writeLines(c("# Baseline-conditioned BRAF pilot", "", "## Aim 1: CCLE/PRISM", "",
  sprintf(paste0("The primary two-stage model is learned only in Tahoe: ",
    "baseline P3 score plus BRAF hotspot status predict vemurafenib-induced ",
    "P1 RAF-program displacement; measured P1 displacement predicts viability. ",
    "Freezing both regressions and applying them without refitting to 38 ",
    "non-Tahoe melanoma lines gives r = %.2f (95%% CI %.2f to %.2f, p = %.4g). ",
    "RMSE is %.2f versus %.2f for transfer of the mean Tahoe response."),
    transfer_mel_stats$pearson, transfer_mel_stats$ci_low,
    transfer_mel_stats$ci_high, transfer_mel_stats$p_value,
    transfer_mel_stats$rmse, transfer_mel_stats$transferred_mean_rmse),
  sprintf(paste0("The fitted Tahoe conditioning equation is P1 displacement ",
    "= %.4f %+.4f x BRAF-hotspot %+.4f x baseline-P3. Its in-sample ",
    "molecular correlation is %.2f; leave-one-line-out correlation is %.2f, ",
    "so the conditioning coefficients require confirmation in a larger ",
    "molecular-response cohort."),
    coef(condition_fits$P1)[1L], coef(condition_fits$P1)["braf_hotspot"],
    coef(condition_fits$P1)["baseline_P3"],
    conditioning[program == "P1", fitted_pearson],
    conditioning[program == "P1", loo_pearson]),
  "Adding inferred P3 or P4 displacement to the Tahoe viability decoder does not materially improve external correlation; their added decoder coefficients are not supported. P1 alone is therefore the parsimonious response readout.", "",
  "As a sensitivity analysis fitted within CCLE/PRISM,", 
  sprintf(paste0("A two-feature model using baseline P3 expression and BRAF ",
    "hotspot status predicts the 2.5 uM vemurafenib viability score in 38 non-Tahoe melanoma ",
    "lines under leave-one-line-out testing: r = %.2f (95%% CI %.2f to %.2f), ",
    "RMSE %.2f versus %.2f for the training-mean predictor."),
    ccle_mel_stats$pearson, ccle_mel_stats$correlation_ci_low,
    ccle_mel_stats$correlation_ci_high, ccle_mel_stats$rmse,
    ccle_mel_stats$reference_rmse),
  sprintf(paste0("BRAF status alone gives r = %.2f and RMSE %.2f. Thus the ",
    "pretreatment expression score adds held-line information beyond mutation. ",
    "In the full descriptive fit, the mutation-adjusted P3 coefficient is ",
    "%.2f (p = %.3f)."),
    ccle_cv_stats[scope == "melanoma" & model == "mutation", pearson],
    ccle_cv_stats[scope == "melanoma" & model == "mutation", rmse],
    ccle_coefficient[term == "baseline_P3", estimate],
    ccle_coefficient[term == "baseline_P3", p_value]),
  sprintf(paste0("The panel contains %d BRAF-hotspot and %d wild-type lines. ",
    "Hotspot lines have %.2f lower mean viability (Wilcoxon p = %.3g)."),
    genotype_contrast[scope == "melanoma", n_hotspot],
    genotype_contrast[scope == "melanoma", n_wild_type],
    -genotype_contrast[scope == "melanoma", hotspot_minus_wild_type_mean],
    genotype_contrast[scope == "melanoma", wilcoxon_p_value]), "",
  "## Aim 2: GSE68453", "",
  sprintf(paste0("The same frozen Tahoe conditioning and viability equations ",
    "were applied to vehicle expression in the 12-line GSE68453 functional ",
    "cohort. The baseline prediction correlates with residual P3 response during ",
    "PLX4720 at r = %.2f and with later PLX4720 + lapatinib Bliss benefit at ",
    "r = %.2f. Among the five lines with factorial RNA, it correlates with the ",
    "lapatinib-specific extra P3 reduction at r = %.2f; that measured reduction ",
    "correlates with Bliss benefit at r = %.2f. All lines are BRAF-mutant, so ",
    "their predicted ranking is supplied by pretreatment P3."),
    gse_bypass_stats[normalization == "program_scaled_P3",
      predicted_to_plx_P3],
    gse_bypass_stats[normalization == "program_scaled_P3",
      predicted_to_bliss],
    gse_bypass_stats[normalization == "program_scaled_P3",
      predicted_to_extra_P3],
    gse_bypass_stats[normalization == "program_scaled_P3",
      extra_P3_to_bliss]),
  sprintf(paste0("The conventional vehicle-mean P3 response shares controls ",
    "with baseline. When the three vehicle replicates are rotated so baseline ",
    "uses two and response subtraction uses the held third, the prediction-to-P3 ",
    "correlation remains %.2f to %.2f."),
    min(gse_split_bypass_stats$predicted_to_held_control_P3),
    max(gse_split_bypass_stats$predicted_to_held_control_P3)),
  sprintf(paste0("Cross-platform scaling is the main technical uncertainty: ",
    "gene-wise scaling gives baseline-to-Bliss r = %.2f, while program-score ",
    "scaling gives r = %.2f."),
    gse_transfer_stats[normalization == "gene_scaled_P3", pearson],
    gse_transfer_stats[normalization == "program_scaled_P3", pearson]),
  sprintf(paste0("In 12 BRAF-mutant melanoma lines, the untreated P3 score ",
    "associates with the measured PLX4720 P3 change at r = %.2f to %.2f. ",
    "Each calculation uses two vehicle replicates for baseline and the third ",
    "for response subtraction, so the association does not reuse its control."),
    min(gse_stats[program == "P3", baseline_response_pearson]),
    max(gse_stats[program == "P3", baseline_response_pearson])),
  sprintf(paste0("The mean pretreatment P3 score predicts PLX4720 + lapatinib ",
    "benefit in a held-out line at r = %.2f (RMSE %.1f versus %.1f for the ",
    "training-mean predictor). The full-data association is r = %.2f, p = %.3f."),
    p3_fs$loo_prediction_pearson, p3_fs$loo_prediction_rmse,
    p3_fs$mean_reference_rmse, p3_fs$baseline_bliss_pearson,
    p3_fs$association_p_value),
  "", "## Scope", "",
  "These are retrospective pilots. The Tahoe-to-CCLE coefficients are transferred without refitting, but the specification was refined while inspecting these pilot results; a new locked holdout is needed for a confirmatory claim. The direct CCLE model excludes each tested line from coefficient fitting. The GSE68453 result is small and should be described as external support.", "",
  "BRAF status and baseline expression use cached DepMap/CCLE 22Q2 files. None of the CCLE/PRISM lines overlap Tahoe."),
  file.path(out, "README.md"))

message("Stage 38 mockup complete: ", out)
