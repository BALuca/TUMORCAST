# Shared paths for the portable single-drug bridge/ligand rescue export. No server paths or downloads.
config_file <- normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
bundle_root <- dirname(dirname(config_file))
input_root <- file.path(bundle_root, 'inputs')
result_root <- file.path(bundle_root, 'results')
regenerated_root <- file.path(bundle_root, 'regenerated')
programs <- sprintf('P%02d', 1:40)

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Required R packages are missing: ', paste(missing, collapse = ', '),
                            '. See provenance/package_versions.tsv. No packages are installed automatically.')
}
