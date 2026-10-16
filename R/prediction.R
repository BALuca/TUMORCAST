# Transport a named source log2-fold-change profile onto recipient TPM abundance,
# applying the historical paired composition correction before program scoring.
# Missing source genes explicitly receive zero change; the count is reported.
transport_gene_response <- function(baseline, fold_change, weights,
  minimum_cp10k = 0.1, fold_change_trim = c(0.3, 0.7), abundance_trim = c(0.05, 0.95)) {
  abundance <- baseline$abundance
  reference <- abundance / sum(abundance) * 10000
  response <- fold_change[match(baseline$gene, names(fold_change))]
  missing <- sum(!is.finite(response))
  response[!is.finite(response)] <- 0
  treated <- abundance * 2^response
  treated <- treated / sum(treated) * 10000
  eligible <- treated >= minimum_cp10k & reference >= minimum_cp10k
  log_ratio <- log2(treated[eligible] / reference[eligible])
  log_abundance <- (log2(treated[eligible]) + log2(reference[eligible])) / 2
  ratio_bounds <- quantile(log_ratio, fold_change_trim)
  abundance_bounds <- quantile(log_abundance, abundance_trim)
  kept <- log_ratio >= ratio_bounds[1] & log_ratio <= ratio_bounds[2] &
    log_abundance >= abundance_bounds[1] & log_abundance <= abundance_bounds[2]
  stopifnot(any(kept))
  correction <- 2^median(log_ratio[kept])
  displacement <- log2(1 + treated / correction) - log2(1 + reference)
  names(displacement) <- baseline$gene
  index <- match(rownames(weights), baseline$gene)
  stopifnot(!anyNA(index), all(is.finite(displacement)))
  list(programs = as.numeric(displacement[index] %*% weights), genes = displacement,
    correction = correction, normalization_genes = sum(kept), missing_source_genes = missing)
}

# Project recipient ERBB2/MET expression through an already fitted lapatinib RNA
# map. Receptors affect the predicted displacement, not the functional readout.
predict_lapatinib <- function(model, metadata, receptors, contexts) {
  stopifnot(all(metadata$drug == "Lapatinib ditosylate"), all(metadata$dose %in% c(0.5, 5)))
  receptor_index <- match(metadata$depmap_id, receptors$depmap_id)
  stopifnot(!anyNA(receptor_index))
  expression <- as.matrix(receptors[receptor_index, c("ERBB2", "MET"), drop = FALSE])
  displacement <- matrix(NA_real_, nrow(metadata), nrow(model$context$W))
  for (dose in unique(metadata$dose)) {
    rows <- which(metadata$dose == dose)
    fit <- if (dose == 5) model$rna$lap_receptors else
      model$rna$extra_rna[[paste("Lapatinib ditosylate", sprintf("%.10g", dose), sep = "|")]]
    stopifnot(!is.null(fit))
    standardized <- sweep(sweep(expression[rows, , drop = FALSE], 2, fit$xmean), 2, fit$xsd, "/")
    predicted <- standardized %*% fit$B
    displacement[rows, ] <- sweep(sweep(predicted, 2, fit$sigma, "*"), 2, fit$mu, "+")
  }
  colnames(displacement) <- rownames(model$context$W)
  context_key <- paste(contexts$drug, sprintf("%.10g", contexts$dose), sep = "|")
  g <- match(paste(metadata$drug, sprintf("%.10g", metadata$dose), sep = "|"), context_key)
  stopifnot(!anyNA(g), all(is.finite(displacement)))
  context <- model$context
  standardized <- (displacement - context$mu[g, , drop = FALSE]) / context$sigma[g, , drop = FALSE]
  signature_scores <- standardized %*% context$W
  route <- context$route[g]
  score <- signature_scores[cbind(seq_len(nrow(metadata)), route)]
  level <- as.numeric(context$Q[g, , drop = FALSE] %*% model$par$level)
  amplitude <- as.numeric(exp(context$Q[g, , drop = FALSE] %*% model$par$log_amplitude))
  list(D = displacement, sensitivity = level + amplitude * score, score = score,
    route = route, level = level, amplitude = amplitude)
}
