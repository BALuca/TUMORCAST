# Compatibility with the original molecular discovery algorithm. This is an R
# translation of D-SPIN's multiplicative updates, not a different NMF method.
# The old NumPy MT19937 initial draws are reproduced so program identity can be
# compared without confusing new random starts with a scientific change.
legacy_uniforms <- function(n, seed) {
  state <- numeric(624)
  state[1] <- seed
  multiplier <- bit64::as.integer64(1812433253)
  modulus <- bit64::as.integer64(4294967296)
  for (i in 2:624) {
    previous <- state[i - 1L]
    mixed <- previous - previous %% 4 +
      bitwXor(as.integer(previous %% 4), as.integer(floor(previous / 2^30)))
    next_state <- (multiplier * bit64::as.integer64(mixed) + bit64::as.integer64(i - 1L)) %% modulus
    state[i] <- as.double(next_state)
  }
  signed_state <- as.integer(ifelse(state >= 2^31, state - 2^32, state))
  assign(".Random.seed", c(10403L, 624L, signed_state), envir = .GlobalEnv)
  draws <- matrix(runif(2 * n), ncol = 2, byrow = TRUE)
  (floor(draws[, 1] * 2^27) * 2^26 + floor(draws[, 2] * 2^26)) / 2^53
}

# Fit one restart to genes x positive/negative correlation channels.
# Return both factors and normalized gene membership; checkpoints retain seeds.
fit_program_restart <- function(affinity, rank, seed, iterations) {
  n_a <- nrow(affinity) * rank
  initial <- legacy_uniforms(n_a + rank * ncol(affinity), seed)
  a <- matrix(initial[seq_len(n_a)], nrow(affinity), rank, byrow = TRUE)
  s <- matrix(initial[-seq_len(n_a)], rank, ncol(affinity), byrow = TRUE)
  s <- abs(t(svd(t(s), nu = rank, nv = 0)$u))
  rm(initial)
  for (iteration in seq_len(iterations)) {
    a <- a * ((affinity %*% t(s)) / (a %*% tcrossprod(s)))
    a[is.nan(a)] <- 0
    atx <- crossprod(a, affinity)
    s <- s * (atx / ((s %*% t(atx)) %*% s))
    s[is.nan(s)] <- 0
    if (iteration %% 100L == 0L) {
      message("Rank ", rank, ", seed ", seed, ": iteration ", iteration, "/", iterations)
    }
    stopifnot(all(is.finite(a)), all(is.finite(s)))
  }
  component_norms <- sqrt(rowSums(s^2))
  s <- s / component_norms
  a <- sweep(a, 2, component_norms, "*")
  scale <- sum(a * (affinity %*% t(s))) / sum(crossprod(a) * tcrossprod(s))
  a <- a * scale
  membership <- round_float32(t(a))
  # Float32 storage is retained at the same algorithm boundaries as the original.
  norms <- round_float32(sqrt(round_float32(rowSums(round_float32(membership^2)))))
  normalized <- round_float32(membership / norms)
  list(membership = normalized, raw_membership = membership, features = round_float32(s))
}

# Align repeated fits by Pearson correlation and one-to-one Hungarian matching.
# Form medoid-aligned consensus, realign once, then name by the anchor gene.
consensus_programs <- function(memberships, genes) {
  align <- function(left, right) {
    correlation <- cor(t(left), t(right))
    correlation[!is.finite(correlation)] <- -1
    assignment <- as.integer(clue::solve_LSAP(correlation + 1, maximum = TRUE))
    list(assignment = assignment, correlation = correlation)
  }
  average_memberships <- function(matrices) {
    total <- matrices[[1]]
    if (length(matrices) > 1L) {
      for (i in 2:length(matrices)) total <- round_float32(total + matrices[[i]])
    }
    average <- round_float32(total / length(matrices))
    norm <- round_float32(sqrt(round_float32(rowSums(round_float32(average^2)))))
    round_float32(average / norm)
  }
  agreement <- matrix(0, length(memberships), length(memberships))
  for (i in seq_along(memberships)) {
    for (j in seq_along(memberships)) {
      if (j <= i) next
      matched <- align(memberships[[i]], memberships[[j]])
      agreement[i, j] <- mean(matched$correlation[cbind(seq_len(nrow(memberships[[i]])),
        matched$assignment)])
      agreement[j, i] <- agreement[i, j]
    }
  }
  medoid <- which.max(rowMeans(agreement))
  aligned <- lapply(memberships, function(x) x[align(memberships[[medoid]], x)$assignment, ])
  consensus <- average_memberships(aligned)
  aligned <- lapply(memberships, function(x) x[align(consensus, x)$assignment, ])
  consensus <- average_memberships(aligned)
  anchors <- max.col(consensus, ties.method = "first")
  consensus <- consensus[order(genes[anchors]), ]
  matched <- t(vapply(memberships, function(x) {
    result <- align(consensus, x)
    result$correlation[cbind(seq_len(nrow(consensus)), result$assignment)]
  }, numeric(nrow(consensus))))
  list(membership = consensus, matched = matched, medoid = medoid)
}

# Orient each component using correlations to its strongest gene; retain the
# signed unit-length component before assigning each gene to its strongest one.
orient_programs <- function(membership, correlation) {
  signed <- membership
  anchors <- max.col(membership, ties.method = "first")
  for (k in seq_len(nrow(membership))) {
    directions <- ifelse(correlation[, anchors[k]] < 0, -1, 1)
    current <- membership[k, ] * directions
    signed[k, ] <- round_float32(current / sqrt(sum(current^2)))
  }
  signed
}
