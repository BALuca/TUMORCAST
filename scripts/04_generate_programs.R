#!/usr/bin/env Rscript
# Sensitivity-blind rank screening followed by ten repeated signed oNMF fits.
# Historical method: D-SPIN updates, 500 iterations, ranks 10/20/40/80.

suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(digest)
  library(clue)
  library(cluster)
})
source("R/data_io.R")
source("R/programs.R")

input_path <- "data/processed/gene_relationships.qs"
output_dir <- "results/04_generate_programs"
checkpoint_dir <- "data/processed/program_discovery"
bank_dir <- "data/processed/program_bank"
ranks <- c(10L, 20L, 40L, 80L)
screening_repeats <- 3L
final_repeats <- 10L
iterations <- 500L
workers <- 4L
threads_per_worker <- 4L
# Inputs and factors are checked explicitly. Call BLAS directly rather than
# repeatedly scanning the entire 1.25-GB affinity matrix for missing values.
options(matprod = "blas")
for (path in c(output_dir, checkpoint_dir, bank_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
sink(file.path(output_dir, "run.log"), split = TRUE)
relationships <- qread(input_path)$result
fingerprint <- list(
  input = digest(relationships, algo = "sha256"),
  code = vapply(c("R/programs.R", "R/data_io.R"),
    digest, character(1), file = TRUE, algo = "sha256"),
  iterations = iterations, R = R.version.string,
  packages = vapply(c("bit64", "clue", "qs"), function(x) as.character(packageVersion(x)), character(1))
)
genes <- relationships$genes
joint <- relationships$joint
affinity <- round_float32(cbind(pmax(joint, 0), pmax(-joint, 0)))
stopifnot(nrow(affinity) == length(genes), all(is.finite(affinity)))
# Include the actual fitting matrix so changing its construction invalidates
# restarts even when the source correlations and fitting functions stay fixed.
fingerprint$affinity <- digest(affinity, algo = "sha256")

# Restart files are useful long-job checkpoints. Each carries the complete
# input/code/seed signature and is replaced only after a successful fit.
fit_repeats <- function(rank, repeats) {
  memberships <- vector("list", repeats)
  for (seed in seq_len(repeats)) {
    path <- file.path(checkpoint_dir, sprintf("rank_%02d_seed_%02d.qs", rank, seed))
    signature <- c(fingerprint, list(rank = rank, seed = seed))
    cached <- if (file.exists(path)) qread(path) else NULL
    if (!is.null(cached) && identical(cached$fingerprint, signature)) {
      message("Reusing rank ", rank, ", seed ", seed)
      fitted <- cached$result
    } else {
      message("Fitting rank ", rank, ", seed ", seed, ": absent or stale checkpoint")
      fitted <- fit_program_restart(affinity, rank, seed, iterations)
      save_cache(fitted, signature, path)
    }
    memberships[[seed]] <- fitted$membership
  }
  consensus <- consensus_programs(memberships, genes)
  consensus$signed <- orient_programs(consensus$membership, joint)
  consensus$assignment <- max.col(t(consensus$membership), ties.method = "first")
  consensus
}

# Distinct ranks run independently. Final fits reuse the first three starts.
screen_rank <- function(rank) {
  RhpcBLASctl::blas_set_num_threads(threads_per_worker)
  log_path <- file.path(output_dir, sprintf("rank_%02d.log", rank))
  log_connection <- file(log_path, "wt")
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    close(log_connection)
  })
  consensus <- fit_repeats(rank, screening_repeats)
  a <- t(consensus$membership)
  atx <- crossprod(a, affinity)
  gram <- crossprod(a)
  coefficients <- solve(gram, atx)
  affinity_energy <- sum(affinity^2)
  residual <- max(0, affinity_energy - 2 * sum(coefficients * atx) +
    sum((gram %*% coefficients) * coefficients)) / affinity_energy
  # Internal reconstruction is used for the historical rank elbow only; it is
  # never a functional performance argument. Gene coherence remains diagnostic.
  distance <- as.dist(pmax(1 - abs(joint), 0))
  silhouette <- cluster::silhouette(consensus$assignment, distance)[, "sil_width"]
  coherence <- unlist(lapply(seq_len(rank), function(k) {
    indices <- which(consensus$assignment == k)
    signs <- sign(consensus$signed[k, indices])
    vapply(c("tahoe", "cytosig"), function(source) {
      block <- relationships[[source]][indices, indices] * outer(signs, signs)
      mean(block[upper.tri(block)])
    }, numeric(1))
  }))
  data.frame(rank = rank, relative_reconstruction_error = residual,
    mean_stability = mean(consensus$matched), q10_stability = unname(quantile(consensus$matched, 0.1)),
    minimum_program_size = min(tabulate(consensus$assignment, nbins = rank)),
    negative_silhouette_fraction = mean(silhouette < 0), mean_silhouette = mean(silhouette),
    minimum_source_coherence = min(coherence))
}
screening <- parallel::mclapply(ranks, screen_rank, mc.cores = workers)
if (any(vapply(screening, inherits, logical(1), "try-error"))) stop("Program rank screening failed")
screening <- bind_rows(screening)
gain <- (screening$relative_reconstruction_error[1] - screening$relative_reconstruction_error) /
  max(screening$relative_reconstruction_error[1] - tail(screening$relative_reconstruction_error, 1), 1e-12)
position <- (ranks - min(ranks)) / diff(range(ranks))
screening <- screening %>%
  mutate(elbow_score = gain - position,
    stability_eligible = q10_stability >= 0.5 & minimum_program_size >= 10,
    original_separation_gate = negative_silhouette_fraction <= 0.2,
    original_source_gate = minimum_source_coherence >= 0.05)
eligible <- screening %>% filter(stability_eligible) %>% arrange(desc(elbow_score), rank)
if (!nrow(eligible)) stop("No rank meets the original stability and size criteria")
selected_rank <- eligible$rank[1]
screening$chosen <- screening$rank == selected_rank
write.table(screening, file.path(output_dir, "rank_screening.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
cat("Selected rank", selected_rank, "without functional outcomes.\n")
RhpcBLASctl::blas_set_num_threads(16L)
consensus <- fit_repeats(selected_rank, final_repeats)
programs <- sprintf("P%02d", seq_len(selected_rank))
weights <- bind_rows(lapply(seq_len(selected_rank), function(k) {
  indices <- which(consensus$assignment == k)
  data.frame(program = programs[k], gene = genes[indices], weight = consensus$signed[k, indices],
    absolute_weight = consensus$membership[k, indices], orientation = sign(consensus$signed[k, indices]))
}))
stopifnot(nrow(weights) == length(genes), !anyDuplicated(weights$gene))
write.table(weights, file.path(bank_dir, "frozen_program_weights.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(weights, file.path(output_dir, "program_weights.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(program = programs, mean_repeat_stability = colMeans(consensus$matched),
  minimum_repeat_stability = apply(consensus$matched, 2, min)),
  file.path(output_dir, "repeat_stability.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
save_cache(consensus, c(fingerprint, list(rank = selected_rank, repeats = final_repeats)),
  file.path(checkpoint_dir, "consensus.qs"))

# The reference is evaluated only after independently choosing the rank and fits.
original_path <- "reference/results/program_bank/frozen_program_weights.tsv"
if (file.exists(original_path)) {
  original <- read.delim(original_path)
  index <- match(weights$gene, original$gene)
  comparison <- data.frame(
    genes = nrow(weights), identical_assignments = sum(weights$program == original$program[index]),
    maximum_weight_difference = max(abs(weights$weight - original$weight[index]))
  )
  write.table(comparison, file.path(output_dir, "reference_comparison.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
  print(comparison)
  if (comparison$identical_assignments != length(genes) || comparison$maximum_weight_difference > 1e-5) {
    stop("Regenerated bank differs from the reference; inspect before downstream reproduction")
  }
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
sink()
