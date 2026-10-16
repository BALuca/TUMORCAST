#!/usr/bin/env Rscript
# Repeat signed-correlation oNMF across an explicit rank grid. Inspect all tested
# bases, not only the chosen rank. Excluded development RNA selects complexity.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(ggplot2)
  library(cluster)
})
source("/work/src/program_pilot/config.R")
input_path <- file.path(data_dir, "molecular_inputs.qs")
correlation_path <- file.path(data_dir, "discovery_correlations.qs")
check(all(file.exists(c(input_path, correlation_path))), "Run molecular preparation and correlations first")
input <- qread(input_path)
correlations <- qread(correlation_path)
correlation <- correlations$joint
check(identical(colnames(input$tahoe), rownames(correlation)), "Discovery genes do not align")

# Positive and negative channels retain the signed gene correlation information.
affinity <- cbind(pmax(correlation, 0), pmax(-correlation, 0))
rank_summaries <- list()
validation_rows <- list()
start_rows <- list()
all_membership_rows <- list()
theme_set(theme_classic(base_size = 10))

for (rank in rank_grid) {
  message("Learning rank ", rank, " with ", n_starts, " independent starts")
  starts <- vector("list", n_starts)
  traces <- list()
  for (start in seq_len(n_starts)) {
    fit_path <- file.path(model_dir, sprintf("onmf_rank_%03d_start_%02d.qs", rank, start))
    signature <- digest::digest(list(correlations = digest::digest(correlation), rank = rank,
      start = start, iterations = onmf_iterations, seed = seed), algo = "sha256")
    if (file.exists(fit_path)) {
      fit <- qread(fit_path)
      check(identical(fit$signature, signature), "Incompatible rank checkpoint; choose a new output directory")
    } else {
      fit <- fit_onmf(affinity, rank, seed + start, onmf_iterations)
      fit$signature <- signature
      qsave(fit, fit_path, preset = "fast")
    }
    starts[[start]] <- fit$membership
    traces[[start]] <- cbind(data.table(rank = rank, start = start), fit$trace)
    message("Rank ", rank, "; start ", start, "; relative error ",
      signif(tail(fit$trace$relative_error, 1L), 4))
  }

  # Align to the medoid start using one-to-one component matching. Consensus
  # memberships are soft; dominant assignments are a plotting partition only.
  agreement <- matrix(1, n_starts, n_starts)
  for (a in seq_len(n_starts - 1L)) {
    for (b in (a + 1L):n_starts) {
      permutation <- match_programs(starts[[a]], starts[[b]])
      agreement[a, b] <- agreement[b, a] <- mean(diag(cor(starts[[a]], starts[[b]][, permutation])))
    }
  }
  reference <- starts[[which.max(rowMeans(agreement))]]
  aligned <- lapply(starts, function(x) x[, match_programs(reference, x), drop = FALSE])
  membership <- Reduce("+", aligned) / n_starts
  membership <- sweep(membership, 2L, sqrt(colSums(membership^2)), "/")
  anchor_names <- input$genes[apply(membership, 2L, which.max)]
  order_programs <- order(anchor_names)
  membership <- membership[, order_programs, drop = FALSE]
  dimnames(membership) <- list(input$genes, sprintf("P%02d", seq_len(rank)))
  signed <- orient_programs(membership, correlation)
  assignment <- max.col(membership, ties.method = "first")
  sizes <- tabulate(assignment, nbins = rank)
  matched <- rbindlist(lapply(seq_along(starts), function(start) {
    permutation <- match_programs(membership, starts[[start]])
    data.table(rank = rank, start = start, program = colnames(membership),
      matched_pearson = diag(cor(membership, starts[[start]][, permutation, drop = FALSE])))
  }))
  start_rows[[as.character(rank)]] <- matched

  # Orthogonal projection into the frozen signed span represents excluded
  # whole-line/study responses. No functional labels or held-drug responses enter.
  inverse <- solve(crossprod(signed))
  validation <- list()
  for (source in c("Tahoe", "CytoSig")) {
    x <- if (source == "Tahoe") input$tahoe else input$cytosig
    metadata <- if (source == "Tahoe") input$tahoe_metadata else input$cytosig_metadata
    unit <- if (source == "Tahoe") metadata$depmap_id else metadata$study
    fold <- if (source == "Tahoe") correlations$line_fold[unit] else correlations$study_fold[unit]
    rows <- which(fold == 1L)
    for (group in unique(unit[rows])) {
      xx <- x[rows[unit[rows] == group], , drop = FALSE]
      # Missing ligand genes are not imputed: only fully observed validation
      # profiles enter reconstruction; omitted profiles are counted explicitly.
      complete <- rowSums(!is.finite(xx)) == 0L
      yy <- xx[complete, , drop = FALSE]
      check(nrow(yy) > 0, paste("No complete excluded validation profile in", group))
      scores <- yy %*% signed
      captured <- rowSums((scores %*% inverse) * scores)
      total <- rowSums(yy^2)
      check(all(total > 0), "Zero-energy excluded molecular response")
      validation[[paste(source, group)]] <- data.table(rank = rank, source = source, unit = group,
        profiles = nrow(yy), omitted_incomplete_profiles = sum(!complete),
        relative_squared_error = mean(pmax(total - captured, 0) / total))
    }
  }
  validation <- rbindlist(validation)
  validation_rows[[as.character(rank)]] <- validation
  source_summary <- validation[, .(error = mean(relative_squared_error),
    se = sd(relative_squared_error) / sqrt(.N)), by = source]
  rank_summaries[[as.character(rank)]] <- data.table(rank = rank,
    excluded_response_error = mean(source_summary$error),
    error_se = sqrt(sum(source_summary$se^2)) / 2,
    median_start_agreement = median(matched$matched_pearson),
    q10_start_agreement = as.numeric(quantile(matched$matched_pearson, 0.1)),
    minimum_program_size = min(sizes),
    mean_orthogonality_error = mean(vapply(starts, function(x) {
      gram <- crossprod(x)
      sqrt(sum((gram - diag(diag(gram)))^2) / sum(gram^2))
    }, numeric(1))))
  qsave(list(membership = membership, signed_basis = signed, agreement = agreement,
    traces = rbindlist(traces), assignment = assignment, program_sizes = sizes),
    file.path(model_dir, sprintf("rank_%03d_consensus.qs", rank)), preset = "fast")
  all_membership_rows[[as.character(rank)]] <- data.table(rank = rank, gene = input$genes,
    program = colnames(membership)[assignment], direction = sign(signed[cbind(seq_along(assignment), assignment)]),
    dominant_membership = membership[cbind(seq_along(assignment), assignment)],
    membership_fraction = apply(membership, 1L, max) / rowSums(membership))

  # Persist diagnostics as each rank finishes so long runs are inspectable.
  summary <- rbindlist(rank_summaries)
  fwrite(summary, file.path(data_dir, "rank_diagnostics.tsv"), sep = "\t")
  fwrite(rbindlist(start_rows), file.path(data_dir, "rank_start_agreement.tsv"), sep = "\t")
  fwrite(rbindlist(validation_rows), file.path(data_dir, "rank_excluded_response_errors.tsv"), sep = "\t")
  fwrite(rbindlist(all_membership_rows), file.path(data_dir, "all_rank_gene_assignments.tsv.gz"), sep = "\t")
  p <- ggplot(summary, aes(rank, excluded_response_error)) +
    geom_line() + geom_point() + geom_errorbar(aes(ymin = excluded_response_error - error_se,
      ymax = excluded_response_error + error_se), width = 0) +
    scale_x_log10(breaks = rank_grid) + labs(x = "Candidate program count",
      y = "Excluded response relative squared error", title = "Molecular rank exploration")
  save_plot(p, "03_rank_excluded_response_error", 5.2, 3.7)
}

# Smallest supported rank within one SE of the best excluded RNA representation.
# Do not force a selection if gene-size/stability criteria fail.
summary[, eligible := q10_start_agreement >= minimum_start_agreement &
  minimum_program_size >= minimum_program_genes]
check(any(summary$eligible), "No candidate rank passes molecular size and repeated-start criteria")
best <- summary[eligible == TRUE][which.min(excluded_response_error)]
selected_rank <- min(summary[eligible == TRUE &
  excluded_response_error <= best$excluded_response_error + best$error_se, rank])
summary[, selected := rank == selected_rank]
fwrite(summary, file.path(data_dir, "rank_diagnostics.tsv"), sep = "\t")
bank <- qread(file.path(model_dir, sprintf("rank_%03d_consensus.qs", selected_rank)))
bank$rank <- selected_rank
bank$genes <- input$genes
bank$discovery_scope <- "Discovery cells/studies only; whole-unit molecular fold 1 used for rank selection"
bank$functional_labels_used <- FALSE
bank$rank_selection <- summary
bank$gene_coverage_required <- "Score only measured bank genes; report restricted coverage and all QC before function"
qsave(bank, file.path(model_dir, "frozen_program_bank.qs"), preset = "fast")
write_manifest("03_discover_programs", c(input_path, correlation_path),
  list(rank_grid = rank_grid, starts = n_starts, iterations = onmf_iterations,
    selected_rank = selected_rank, frozen_bank_sha256 = digest::digest(
      file.path(model_dir, "frozen_program_bank.qs"), file = TRUE, algo = "sha256"),
    selection = "minimum eligible rank within one SE of best excluded molecular error",
    functional_labels_used = FALSE))
message("Frozen molecular bank at rank ", selected_rank)
