suppressPackageStartupMessages(library(data.table))

fail <- function(...) stop(sprintf(...), call. = FALSE)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0)
    return(NA_real_)
  stats::cor(x[ok], y[ok], method = "pearson")
}

weighted_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(x) < 3L) return(NA_real_)
  w <- w / sum(w); mx <- sum(w * x); my <- sum(w * y)
  vx <- sum(w * (x - mx) ^ 2); vy <- sum(w * (y - my) ^ 2)
  if (vx <= 0 || vy <= 0) return(NA_real_)
  sum(w * (x - mx) * (y - my)) / sqrt(vx * vy)
}

balanced_weights <- function(group) {
  counts <- table(group)
  1 / as.numeric(counts[group])
}

weighted_sd <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok] / sum(w[ok]); center <- sum(w * x)
  sqrt(sum(w * (x - center) ^ 2))
}

make_inner_assignment <- function(lines, seed = SEED + OUTER_FOLD,
                                  folds = N_INNER_FOLDS) {
  set.seed(seed)
  shuffled <- sample(sort(unique(lines)), length(unique(lines)), FALSE)
  setNames(rep(seq_len(folds), length.out = length(shuffled)), shuffled)
}

standardize_design <- function(x_train, x_test) {
  center <- colMeans(x_train); scale <- apply(x_train, 2L, stats::sd)
  keep <- is.finite(center) & is.finite(scale) & scale > 0
  if (!any(keep)) fail("Ridge design has no variable feature.")
  list(train = sweep(sweep(x_train[, keep, drop = FALSE], 2L,
    center[keep], "-"), 2L, scale[keep], "/"),
    test = sweep(sweep(x_test[, keep, drop = FALSE], 2L,
      center[keep], "-"), 2L, scale[keep], "/"), center = center,
    scale = scale, keep = keep, feature_names = colnames(x_train))
}

ridge_path <- function(x_train, y_train, x_test, lambda) {
  design <- standardize_design(x_train, x_test)
  y_train <- as.matrix(y_train); y_mean <- colMeans(y_train)
  decomposition <- svd(design$train, nu = min(dim(design$train)),
    nv = min(dim(design$train)))
  uy <- crossprod(decomposition$u, sweep(y_train, 2L, y_mean, "-"))
  cross <- design$test %*% decomposition$v
  result <- matrix(NA_real_, nrow(x_test), length(lambda))
  for (i in seq_along(lambda)) result[, i] <- as.numeric(sweep(
    cross %*% sweep(uy, 1L, decomposition$d /
      (decomposition$d ^ 2 + lambda[i]), "*"), 2L, y_mean, "+"))
  result
}

fit_ridge <- function(x, y, lambda) {
  design <- standardize_design(x, x); y <- as.matrix(y)
  y_mean <- colMeans(y)
  decomposition <- svd(design$train, nu = min(dim(design$train)),
    nv = min(dim(design$train)))
  beta <- decomposition$v %*% sweep(crossprod(decomposition$u,
    sweep(y, 2L, y_mean, "-")), 1L, decomposition$d /
      (decomposition$d ^ 2 + lambda), "*")
  list(lambda = lambda, center = design$center, scale = design$scale,
    keep = design$keep, feature_names = design$feature_names,
    beta_standardized = beta, outcome_mean = y_mean)
}

predict_ridge <- function(fit, x) {
  z <- sweep(sweep(x[, fit$keep, drop = FALSE], 2L,
    fit$center[fit$keep], "-"), 2L, fit$scale[fit$keep], "/")
  as.numeric(sweep(z %*% fit$beta_standardized, 2L,
    fit$outcome_mean, "+"))
}

raw_ridge_coefficients <- function(fit) {
  beta <- setNames(rep(0, length(fit$feature_names)), fit$feature_names)
  beta[fit$feature_names[fit$keep]] <- as.numeric(fit$beta_standardized) /
    fit$scale[fit$keep]
  intercept <- as.numeric(fit$outcome_mean) -
    sum(fit$center[fit$keep] * beta[fit$feature_names[fit$keep]])
  list(intercept = intercept, beta = beta)
}

condition_equal_weights <- function(perturbation_id) {
  counts <- table(perturbation_id)
  weights <- length(perturbation_id) /
    (length(counts) * as.numeric(counts[perturbation_id]))
  if (any(!is.finite(weights)) || any(weights <= 0))
    fail("Condition-equal weights are invalid.")
  weights
}

weighted_standardize <- function(x_train, x_test, weights) {
  w <- weights / sum(weights)
  center <- colSums(x_train * w)
  centered <- sweep(x_train, 2L, center, "-")
  scale <- sqrt(colSums(centered ^ 2 * w))
  keep <- is.finite(center) & is.finite(scale) & scale > 0
  if (!any(keep)) fail("Weighted ridge has no variable feature.")
  list(train = sweep(centered[, keep, drop = FALSE], 2L, scale[keep], "/"),
    test = sweep(sweep(x_test[, keep, drop = FALSE], 2L, center[keep], "-"),
      2L, scale[keep], "/"), center = center, scale = scale, keep = keep,
    feature_names = colnames(x_train))
}

weighted_ridge_path <- function(x_train, y, x_test, weights, lambda) {
  design <- weighted_standardize(x_train, x_test, weights)
  root <- sqrt(weights); y_mean <- sum(weights * y) / sum(weights)
  decomposition <- svd(design$train * root,
    nu = min(dim(design$train)), nv = min(dim(design$train)))
  uy <- as.numeric(crossprod(decomposition$u, (y - y_mean) * root))
  cross <- design$test %*% decomposition$v
  result <- matrix(NA_real_, nrow(x_test), length(lambda))
  for (i in seq_along(lambda)) result[, i] <- y_mean + as.numeric(
    cross %*% (uy * decomposition$d /
      (decomposition$d ^ 2 + lambda[i])))
  result
}

fit_weighted_ridge <- function(x, y, weights, lambda) {
  design <- weighted_standardize(x, x, weights)
  root <- sqrt(weights); y_mean <- sum(weights * y) / sum(weights)
  decomposition <- svd(design$train * root,
    nu = min(dim(design$train)), nv = min(dim(design$train)))
  beta <- decomposition$v %*% (as.numeric(crossprod(decomposition$u,
    (y - y_mean) * root)) * decomposition$d /
      (decomposition$d ^ 2 + lambda))
  list(lambda = lambda, center = design$center, scale = design$scale,
    keep = design$keep, feature_names = design$feature_names,
    beta_standardized = beta, outcome_mean = y_mean)
}

read_npy_float32 <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- rawToChar(readBin(con, "raw", 6L))
  version <- readBin(con, "integer", 2L, size = 1L, signed = FALSE)
  if (!identical(magic, "\x93NUMPY") || !version[1L] %in% c(1L, 2L, 3L))
    fail("Unsupported NPY header: %s", path)
  length_size <- if (version[1L] == 1L) 2L else 4L
  header_length <- readBin(con, "integer", 1L, size = length_size,
    signed = FALSE, endian = "little")
  header <- rawToChar(readBin(con, "raw", header_length))
  fortran <- grepl("'fortran_order': True", header, fixed = TRUE)
  if (!grepl("'descr': '<f4'", header, fixed = TRUE) ||
      !(fortran || grepl("'fortran_order': False", header, fixed = TRUE)))
    fail("Unsupported NPY matrix: %s", path)
  shape_text <- sub(".*'shape': \\(([^)]*)\\).*", "\\1", header)
  shape <- suppressWarnings(as.integer(trimws(strsplit(shape_text, ",")[[1L]])))
  shape <- shape[is.finite(shape)]
  if (length(shape) != 2L) fail("NPY matrix is not two-dimensional: %s", path)
  values <- readBin(con, "numeric", n = prod(shape), size = 4L,
    endian = "little")
  matrix(values, nrow = shape[1L], ncol = shape[2L], byrow = !fortran)
}

read_npy_float32_columns <- function(path, columns_one_based) {
  columns_one_based <- as.integer(columns_one_based)
  if (!length(columns_one_based) || anyNA(columns_one_based) ||
      any(columns_one_based < 1L) || anyDuplicated(columns_one_based))
    fail("Invalid requested NPY columns: %s", path)
  column_file <- tempfile("npy_columns_", fileext = ".txt")
  matrix_file <- tempfile("npy_subset_", fileext = ".npy")
  on.exit(unlink(c(column_file, matrix_file)), add = TRUE)
  writeLines(as.character(columns_one_based - 1L), column_file)
  status <- system2("/usr/bin/python3",
    c(file.path(SCRIPT_DIR, "export_npy_columns.py"), path, column_file,
      matrix_file))
  if (!identical(status, 0L))
    fail("Native NPY column extraction failed with status %d: %s", status,
      path)
  read_npy_float32(matrix_file)
}

function_residual_table <- function(table, metric,
                                    prediction = "predicted_function") {
  table <- copy(table[is.finite(observed_function) &
    is.finite(get(prediction))])
  table[, predicted_value := get(prediction)]
  if (metric == "pooled") {
    table[, `:=`(observed_residual = observed_function,
      predicted_residual = predicted_value)]
  } else if (metric == "within_perturbation") {
    table[, `:=`(observed_residual = observed_function - mean(observed_function),
      predicted_residual = predicted_value - mean(predicted_value)),
      by = perturbation_id]
  } else if (metric == "strict_cell_line_by_perturbation") {
    table[, observed_residual := observed_function -
      ave(observed_function, perturbation_id, FUN = mean) -
      ave(observed_function, depmap_id, FUN = mean) + mean(observed_function)]
    table[, predicted_residual := predicted_value -
      ave(predicted_value, perturbation_id, FUN = mean) -
      ave(predicted_value, depmap_id, FUN = mean) + mean(predicted_value)]
  } else fail("Unknown functional metric: %s", metric)
  table
}

endpoint_metrics <- function(table, observed, predicted, reference = NULL) {
  y <- table[[observed]]; p <- table[[predicted]]
  keep <- is.finite(y) & is.finite(p); table <- table[keep]
  y <- y[keep]; p <- p[keep]; w <- balanced_weights(table$drug)
  denominator <- if (is.null(reference)) sum(y ^ 2) else
    sum((y - table[[reference]]) ^ 2)
  weighted_denominator <- if (is.null(reference)) sum(w * y ^ 2) else
    sum(w * (y - table[[reference]]) ^ 2)
  data.table(n_rows = nrow(table), n_conditions = uniqueN(table$perturbation_id),
    n_drugs = uniqueN(table$drug), n_lines = uniqueN(table$depmap_id),
    pearson = safe_cor(y, p), drug_balanced_pearson = weighted_cor(y, p, w),
    observed_sd = sd(y), predicted_sd = sd(p),
    predicted_to_observed_sd_ratio = sd(p) / sd(y),
    calibration_slope = coef(lm(y ~ p))[2L],
    rmse = sqrt(mean((y - p) ^ 2)), mae = mean(abs(y - p)),
    partial_r2 = 1 - sum((y - p) ^ 2) / denominator,
    drug_balanced_partial_r2 = 1 - sum(w * (y - p) ^ 2) /
      weighted_denominator, correlation_method = "pearson")
}

write_artifact_manifest <- function(folder, producing_script) {
  files <- list.files(folder, recursive = TRUE, full.names = TRUE,
    include.dirs = FALSE)
  files <- files[basename(files) != "artifact_manifest.tsv"]
  dimensions <- function(path) {
    extension <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
    if (extension %in% c("tsv", "csv")) {
      n <- tryCatch(nrow(fread(path, select = 1L, showProgress = FALSE)),
        error = function(e) NA_integer_)
      sprintf("%s rows", format(n, big.mark = ","))
    } else if (extension == "rds") "serialized R object" else
      if (extension == "npy") "NumPy matrix" else "artifact"
  }
  description <- function(path) {
    name <- sub("\\.gz$", "", basename(path))
    name <- tools::file_path_sans_ext(name)
    paste(gsub("[_-]+", " ", name),
      if (grepl("\\.rds$", path)) "serialized R object" else
        if (grepl("\\.tsv(\\.gz)?$", path)) "table" else "artifact")
  }
  manifest <- data.table(relative_path = substring(files, nchar(folder) + 2L),
    description = vapply(files, description, character(1L)),
    dimensions_or_row_count = vapply(files, dimensions, character(1L)),
    file_size_bytes = file.info(files)$size,
    md5 = unname(tools::md5sum(files)), producing_script = producing_script)
  fwrite(manifest, file.path(folder, "artifact_manifest.tsv"), sep = "\t")
}
