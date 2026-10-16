#!/usr/bin/env Rscript
# Reproduce frozen HGF transport from corrected source gene responses and recipient RNA.
# These are the same normalization/projection operations as the original preparation.
entry <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
source(file.path(dirname(normalizePath(entry)), 'config.R'))
check_packages('data.table')
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
dir.create(regenerated_root, recursive = TRUE, showWarnings = FALSE)
archive <- readRDS(file.path(input_root, 'molecular_functional_inputs/inputs.rds'))
W <- archive$prep$W
reconstruction <- readRDS(file.path(input_root, 'HGF_gene_sources/HGF_input_reconstruction.rds'))
SX <- reconstruction$corrected_gene_profiles
studies <- colnames(SX)
profiles <- cbind(Primary_study_mean = rowMeans(SX, na.rm = TRUE), SX,
                 vapply(studies, function(st) rowMeans(SX[, studies != st, drop = FALSE], na.rm = TRUE),
                        numeric(nrow(SX))))
colnames(profiles) <- c('Primary_study_mean', studies, paste0('Leave_out_', studies))
profiles[!is.finite(profiles)] <- 0
expression <- fread(file.path(input_root, 'HGF_gene_sources/recipient_expression.tsv.gz'))
saved <- fread(file.path(input_root, 'ligand_response_inputs/data/time_matched_HGF_program_inputs.tsv'))
cells <- unique(saved[, .(CellLine, depmap_id)])
stopifnot(all(cells$depmap_id %in% expression$depmap_id))
outputs <- list()
for (i in seq_len(nrow(cells))) {
  b <- expression[depmap_id == cells$depmap_id[i]]
  stopifnot(nrow(b) == 1)
  genes <- sub(' \\(ENSG[0-9]+(?:\\.[0-9]+)?\\)$', '', names(b)[-1], perl = TRUE)
  raw <- data.table(gene = genes, u = pmax(2^as.numeric(b[, -1]) - 1, 0))[, .(u = mean(u)), by = gene]
  u <- raw$u
  reference <- u / sum(u) * 1e4
  index <- match(raw$gene, rownames(profiles))
  gi <- match(rownames(W), raw$gene)
  stopifnot(!anyNA(gi))
  Z <- matrix(NA_real_, ncol(profiles), ncol(W))
  for (k in seq_len(ncol(profiles))) {
    fc <- profiles[index, k]
    fc[!is.finite(fc)] <- 0
    treated <- u * 2^fc
    treated <- treated / sum(treated) * 1e4
    eligible <- treated >= .1 & reference >= .1
    mv <- log2(treated[eligible] / reference[eligible])
    av <- (log2(treated[eligible]) + log2(reference[eligible])) / 2
    ml <- quantile(mv, c(.3, .7))
    al <- quantile(av, c(.05, .95))
    central <- mv >= ml[1] & mv <= ml[2] & av >= al[1] & av <= al[2]
    factor <- 2^median(mv[central])
    displacement <- log2(1 + treated/factor) - log2(1 + reference)
    Z[k, ] <- as.numeric(displacement[gi] %*% W)
  }
  colnames(Z) <- colnames(W)
  outputs[[i]] <- cbind(cells[rep(i, nrow(Z))], data.table(gene = 'HGF', source_profile = colnames(profiles)),
                         as.data.table(Z))
}
rebuilt <- rbindlist(outputs)
key <- function(x) paste(x$depmap_id, x$source_profile, sep = '|')
stopifnot(!anyDuplicated(key(saved)), !anyDuplicated(key(rebuilt)))
ix <- match(key(saved), key(rebuilt))
stopifnot(!anyNA(ix))
difference <- max(abs(as.matrix(saved[, ..programs]) - as.matrix(rebuilt[ix, ..programs])))
stopifnot(difference < 1e-10)
fwrite(rebuilt, file.path(regenerated_root, 'HGF_program_inputs_rebuilt.tsv'), sep = '\t')
fwrite(data.table(component = 'HGF gene-response transport and projection', rows = nrow(saved),
                  maximum_absolute_difference = difference),
       file.path(regenerated_root, 'HGF_projection_replay_check.tsv'), sep = '\t')
cat('HGF molecular input reconstruction agrees with all frozen program inputs.\n')
