# Read the original release's two-dimensional NPY matrices without executing Python.
# This reader is retained from src/program_pilot/utilities.R, with explicit checks.
read_npy_matrix <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  magic <- readBin(connection, "raw", 6L)
  stopifnot(identical(magic, as.raw(c(147, 78, 85, 77, 80, 89))))
  version <- readBin(connection, "integer", 2L, size = 1L, signed = FALSE)
  stopifnot(version[1] %in% 1:3)
  header_size <- readBin(connection, "integer", 1L,
    size = if (version[1] == 1L) 2L else 4L, signed = FALSE, endian = "little")
  header <- rawToChar(readBin(connection, "raw", header_size))
  shape_text <- sub(".*'shape': \\(([^)]*)\\).*", "\\1", header)
  dimensions <- as.integer(trimws(strsplit(shape_text, ",", fixed = TRUE)[[1]]))
  dimensions <- dimensions[!is.na(dimensions)]
  stopifnot(length(dimensions) == 2L, all(dimensions > 0), grepl("'<f[48]'", header))
  bytes <- if (grepl("'<f8'", header, fixed = TRUE)) 8L else 4L
  fortran <- grepl("'fortran_order': True", header, fixed = TRUE)
  stopifnot(fortran || grepl("'fortran_order': False", header, fixed = TRUE))
  values <- readBin(connection, "numeric", prod(dimensions), size = bytes, endian = "little")
  stopifnot(length(values) == prod(dimensions))
  matrix(values, dimensions[1], dimensions[2], byrow = !fortran)
}

# Preserve the float32 rounding used by the released arrays and historical fits.
# R continues calculations in double precision after this explicit storage rounding.
round_float32 <- function(x) {
  rounded <- readBin(writeBin(as.numeric(x), raw(), size = 4L), "numeric", length(x), size = 4L)
  attributes(rounded) <- attributes(x)
  rounded
}

# Interrupted writes never replace a completed expensive calculation.
save_cache <- function(result, fingerprint, path) {
  partial <- paste0(path, ".partial")
  qs::qsave(list(fingerprint = fingerprint, result = result), partial, preset = "fast")
  if (!file.rename(partial, path)) stop("Cannot replace cache: ", path)
}
