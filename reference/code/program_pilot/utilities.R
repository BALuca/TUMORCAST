# Shared I/O, provenance and modeling calculations used by multiple numbered
# stages. Scientific decisions remain visible in the scripts that own them.

# Validate a condition and fail with an informative error.
check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

# Standardize intervention names across source tables; salts are source aliases.
canonical_drug <- function(x) {
  x <- tolower(trimws(as.character(x)))
  aliases <- c(bez235 = "dactolisib", `nvp-bez235` = "dactolisib",
    azd6244 = "selumetinib", `byl-719` = "alpelisib", byl719 = "alpelisib",
    `brd-k54997624` = "alpelisib", `lapatinib ditosylate` = "lapatinib")
  matched <- match(x, names(aliases))
  x[!is.na(matched)] <- aliases[matched[!is.na(matched)]]
  x
}

# Map rounded source doses to nominal micromolar doses, never interpolate.
match_nominal_dose <- function(x) {
  vapply(as.numeric(x), function(value) {
    if (!is.finite(value) || value < 0) return(NA_real_)
    if (value == 0) return(0)
    distance <- abs(log(nominal_doses / value))
    if (min(distance) > dose_log_tolerance) return(NA_real_)
    nominal_doses[which.min(distance)]
  }, numeric(1))
}

# Read the cached two-dimensional NPY float matrices entirely in R. NPY files
# are a data format, not permission to execute the historical Python pipeline.
# Both C and Fortran ordering are explicit; unsupported dtypes fail.
read_npy_matrix <- function(path) {
  check(file.exists(path), paste("Missing matrix:", path))
  connection <- file(path, "rb")
  on.exit(close(connection))
  magic <- readBin(connection, "raw", 6L)
  check(identical(magic, as.raw(c(147, 78, 85, 77, 80, 89))), "Invalid NPY magic")
  version <- readBin(connection, "integer", 2L, size = 1L, signed = FALSE)
  check(version[1] %in% 1:3, "Unsupported NPY version")
  header_size <- readBin(connection, "integer", 1L,
    size = if (version[1] == 1L) 2L else 4L, signed = FALSE, endian = "little")
  header <- rawToChar(readBin(connection, "raw", header_size))
  shape_text <- sub(".*'shape': \\(([^)]*)\\).*", "\\1", header)
  dimensions <- as.integer(trimws(strsplit(shape_text, ",", fixed = TRUE)[[1]]))
  dimensions <- dimensions[!is.na(dimensions)]
  check(length(dimensions) == 2L && all(dimensions > 0), "NPY must be a nonempty matrix")
  check(grepl("'<f[48]'", header), "Only little-endian float32/float64 NPY is supported")
  bytes <- if (grepl("'<f8'", header, fixed = TRUE)) 8L else 4L
  fortran <- grepl("'fortran_order': True", header, fixed = TRUE)
  check(fortran || grepl("'fortran_order': False", header, fixed = TRUE), "Unknown NPY ordering")
  values <- readBin(connection, "numeric", prod(dimensions), size = bytes, endian = "little")
  check(length(values) == prod(dimensions), "Truncated NPY matrix")
  matrix(values, dimensions[1], dimensions[2], byrow = !fortran)
}

# Balance groups by assigning each drug/study equal total weight. Within a group
# each profile has equal weight; return positive weights summing to one.
group_weights <- function(group) {
  check(!anyNA(group) && all(nzchar(group)), "Missing balancing group")
  counts <- table(group)
  weights <- 1 / as.numeric(counts[group])
  weights / sum(weights)
}

# Weighted Spearman correlation of genes (columns) across profiles (rows).
# Missing CytoSig observations are permitted and handled pairwise; missingness
# is reported by preparation. Constant genes fail instead of creating a basis.
weighted_gene_correlation <- function(x, weights) {
  check(nrow(x) == length(weights), "Correlation weights do not align")
  ranks <- matrixStats::colRanks(x, ties.method = "average", preserveShape = TRUE)
  valid <- is.finite(ranks)
  if (all(valid)) {
    centered <- sweep(ranks, 2L, colSums(ranks * weights), "-")
    covariance <- crossprod(centered * sqrt(weights))
    scale <- sqrt(diag(covariance))
    check(all(scale > 0), "A discovery gene has constant ranked responses")
    correlation <- covariance / outer(scale, scale)
  } else {
    ranks[!valid] <- 0
    denominator <- crossprod(valid * sqrt(weights))
    sums <- crossprod(ranks * weights, valid)
    squares <- crossprod(ranks^2 * weights, valid)
    covariance <- crossprod(ranks * sqrt(weights)) - sums * t(sums) / denominator
    variance <- squares - sums^2 / denominator
    correlation <- covariance / sqrt(pmax(variance, 0) * pmax(t(variance), 0))
    check(all(is.finite(correlation)), "Insufficient shared variation in CytoSig gene pairs")
  }
  correlation <- pmax(-1, pmin(1, correlation))
  dimnames(correlation) <- list(colnames(x), colnames(x))
  diag(correlation) <- 1
  correlation
}

# Orthogonal NMF multiplicative updates from the existing D-SPIN implementation:
# A <- A * XS'/(ASS'); S <- S * A'X/(SXA S), implemented with BLAS in R.
# Input is nonnegative genes x signed-correlation features. Orthogonality is a
# soft algorithmic target, not guaranteed exact: residuals are returned and shown.
# Return membership A, feature basis S, and iteration/error diagnostics.
fit_onmf <- function(x, rank, random_seed, iterations) {
  check(all(is.finite(x)) && all(x >= 0), "oNMF input must be finite and nonnegative")
  set.seed(random_seed)
  a <- matrix(runif(nrow(x) * rank), nrow(x), rank)
  s <- abs(t(qr.Q(qr(matrix(runif(ncol(x) * rank), ncol(x), rank)))))
  x_energy <- sum(x^2)
  trace_rows <- list()
  for (iteration in seq_len(iterations)) {
    a <- a * ((x %*% t(s)) / pmax(a %*% tcrossprod(s), 1e-30))
    atx <- crossprod(a, x)
    s <- s * (atx / pmax((s %*% t(atx)) %*% s, 1e-30))
    if (iteration %% 20L == 0L || iteration == iterations) {
      ss <- tcrossprod(s)
      # D-SPIN rescales the fitted product at the end. Report the error after
      # that same optimal scalar correction rather than an arbitrary iterate scale.
      product_energy <- sum(crossprod(a) * ss)
      error <- x_energy - sum(atx * s)^2 / product_energy
      off_diagonal <- ss - diag(diag(ss))
      trace_rows[[length(trace_rows) + 1L]] <- data.table::data.table(
        iteration = iteration, relative_error = sqrt(max(error, 0) / x_energy),
        orthogonality_error = sqrt(sum(off_diagonal^2) / sum(s^2)))
    }
  }
  norms <- sqrt(rowSums(s^2))
  check(all(norms > 0) && all(is.finite(a)), "oNMF produced an empty/nonfinite component")
  s <- s / norms
  a <- sweep(a, 2L, norms, "*")
  multiplier <- sum(crossprod(a, x) * s) / sum(crossprod(a) * tcrossprod(s))
  a <- a * multiplier
  a <- sweep(a, 2L, sqrt(colSums(a^2)), "/")
  list(membership = a, feature_basis = s, trace = data.table::rbindlist(trace_rows))
}

# One-to-one component alignment, shared by repeated starts and refit audits.
# Matrices have genes in rows and programs in columns; return column permutation.
match_programs <- function(reference, candidate) {
  similarity <- cor(reference, candidate)
  check(all(is.finite(similarity)), "Program matching contains undefined correlations")
  as.integer(clue::solve_LSAP((similarity + 1) / 2, maximum = TRUE))
}

# Orient each gene towards the strongest unsigned anchor using signed correlation.
# The resulting named, unit-L2 signed gene basis is frozen for projection.
orient_programs <- function(membership, correlation) {
  anchors <- apply(membership, 2L, which.max)
  signed <- membership
  for (j in seq_len(ncol(membership))) {
    signed[, j] <- membership[, j] * ifelse(correlation[, anchors[j]] < 0, -1, 1)
  }
  signed <- sweep(signed, 2L, sqrt(colSums(signed^2)), "/")
  dimnames(signed) <- dimnames(membership)
  signed
}

# Export a compact plot in both shareable PDF and directly inspectable PNG.
save_plot <- function(plot, name, width, height) {
  ggplot2::ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
    width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
    width = width, height = height, dpi = 160, bg = "white")
}

# Draw a shared-wrapper heatmap, not a custom plotting implementation.
save_heatmap <- function(heatmap, name, width, height, destination = figure_dir, annotation_legends = list()) {
  grDevices::cairo_pdf(file.path(destination, paste0(name, ".pdf")), width, height)
  ComplexHeatmap::draw(heatmap, annotation_legend_list = annotation_legends)
  grDevices::dev.off()
  grDevices::png(file.path(destination, paste0(name, ".png")),
    width = width, height = height, units = "in", res = 160)
  ComplexHeatmap::draw(heatmap, annotation_legend_list = annotation_legends)
  grDevices::dev.off()
}

# Stage provenance is kept alongside results. Record only paths/hashes, no raw
# secrets or external account information. Stage input/output objects stay plain.
write_manifest <- function(stage, inputs, details) {
  check(all(file.exists(inputs)), paste("Missing manifest inputs for", stage))
  files <- data.table::data.table(path = inputs,
    bytes = file.info(inputs)$size,
    sha256 = vapply(inputs, digest::digest, character(1), file = TRUE, algo = "sha256"))
  code <- list.files(code_dir, pattern = "\\.R$", full.names = TRUE)
  code_hashes <- vapply(code, digest::digest, character(1), file = TRUE, algo = "sha256")
  revision <- system2("git", c("-C", work_dir, "rev-parse", "HEAD"), stdout = TRUE)
  jsonlite::write_json(list(stage = stage, timestamp = as.character(Sys.time()),
    seed = seed, inputs = files, code_sha256 = as.list(code_hashes), git_revision = revision,
    details = details, session = capture.output(sessionInfo())),
    file.path(output_dir, paste0(stage, "_manifest.json")), pretty = TRUE, auto_unbox = TRUE)
}
