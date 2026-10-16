# Shared ridge fitting and observation weights for script 20 evaluation modes.
# Uses the calling workflow settings: lambda_grid, inner_map and programs.
fit_path <- function(x, y, w, intercept) {
  multi <- is.matrix(y)
  center <- if (multi) colMeans(y) else 0
  scale <- if (multi) apply(y, 2L, sd) else 1
  stopifnot(all(is.finite(scale)), all(scale > 1e-12), all(is.finite(x)), all(is.finite(y)))
  target <- if (multi) sweep(sweep(y, 2L, center, '-'), 2L, scale, '/') else y
  fit <- glmnet(x, target, family = if (multi) 'mgaussian' else 'gaussian', alpha = 0,
    lambda = lambda_grid, weights = w, standardize = TRUE, intercept = intercept, thresh = 1e-7)
  list(fit = fit, center = center, scale = scale, multi = multi)
}
pred_path <- function(f, x, s = lambda_grid) {
  p <- predict(f$fit, newx = x, s = s)
  if (!f$multi) return(as.matrix(p))
  for (j in seq_along(f$center)) p[, j, ] <- p[, j, ] * f$scale[j] + f$center[j]
  p
}
fit_ridge <- function(x, y, ids, w = rep(1, nrow(x)), intercept = TRUE) {
  fold <- unname(inner_map[ids])
  stopifnot(!anyNA(fold), length(unique(fold)) >= 3L)
  loss <- numeric(length(lambda_grid))
  mass <- 0
  for (k in sort(unique(fold))) {
    tr <- which(fold != k)
    va <- which(fold == k)
    f <- fit_path(x[tr, , drop = FALSE], if (is.matrix(y)) y[tr, , drop = FALSE] else y[tr], w[tr], intercept)
    p <- pred_path(f, x[va, , drop = FALSE])
    if (is.matrix(y)) {
      for (j in seq_along(lambda_grid)) {
        error <- sweep(p[, , j] - y[va, , drop = FALSE], 2L, f$scale, '/')^2
        loss[j] <- loss[j] + sum(rowMeans(error) * w[va])
      }
    } else loss <- loss + colSums((p - y[va])^2 * w[va])
    mass <- mass + sum(w[va])
  }
  fit <- fit_path(x, y, w, intercept)
  fit$lambda <- lambda_grid[which.min(loss)]
  fit$inner_mse <- loss / mass
  fit$training_ids <- unique(ids)
  fit
}
predict_ridge <- function(f, x) {
  p <- pred_path(f, x, f$lambda)
  if (!f$multi) return(as.numeric(p))
  matrix(p, nrow(x), length(f$center), dimnames = list(NULL, programs))
}
cor_safe <- function(x, y) if (length(x) >= 3L && sd(x) > 0 && sd(y) > 0) cor(x, y) else NA_real_
# Equal total drug mass, equal cell mass within drug, doses retain individual rows.
observation_weights <- function(d, group = 'drug') {
  z <- copy(d[, c(group, 'depmap_id'), with = FALSE])
  z[, n_dose := .N, by = c(group, 'depmap_id')]
  z[, n_cell := uniqueN(depmap_id), by = group]
  w <- 1 / (z$n_dose * z$n_cell)
  w / mean(w)
}
