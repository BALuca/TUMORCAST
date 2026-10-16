#!/usr/bin/env Rscript
# Protect the functional interaction model and learn RNA from cross-fitted errors.
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(glmnet)
  library(ggplot2)
  library(digest)
})
root <- '/work/results/program_pilot/70_protected_molecular_interactions'
previous <- '/work/results/program_pilot/64_nonlinear_response_composition'
links_path <- '/work/results/program_pilot/61_individual_responses_combination_transfer/data/matched_molecular_links.tsv.gz'
basis_path <- file.path(previous, 'models/3h_kernel_bases.qs')
script_path <- '/work/src/program_pilot/70_protected_molecular_interactions.R'
lambda_grid <- 10^seq(3, -4, length.out = 24L)
set.seed(70011L)
setDTthreads(2L)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
options(warn = 1)
for (folder in c('data', 'models', 'figures')) {
  dir.create(file.path(root, folder), recursive = TRUE, showWarnings = FALSE)
}
write_table <- function(x, name) {
  fwrite(x, file.path(root, 'data', paste0(name, '.tsv.gz')), sep = '\t')
}
status <- function(message_text) {
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ' | ', message_text)
  writeLines(message_text, file.path(root, 'STATUS.txt'))
}
options(error = function() {
  status('FAILED; see log')
  traceback(5)
  quit(save = 'no', status = 1L)
})
status('Load existing individual-response kernels; no f refitting')
bases <- qread(basis_path)
obs <- as.data.table(read.delim(gzfile(links_path)))[time_h == 3 & molecular_supported == TRUE]
setorder(obs, anchor, partner, nominal_A_uM, nominal_B_uM, biological_replicate)
stopifnot(nrow(obs) == 496L, uniqueN(obs$condition_id) == 198L,
  !anyDuplicated(obs$observation_id), identical(obs$observation_id, bases$observation_ids),
  all(is.finite(obs$ratio_AB)), all(obs$ratio_AB > 0), all(obs$ratio_A > 0), all(obs$ratio_B > 0))
target <- log2(obs$ratio_AB / (obs$ratio_A * obs$ratio_B))
partners <- sort(unique(obs$partner))
stopifnot(length(partners) == 6L)
for (kind in c('functional', 'RNA')) {
  for (basis in bases[[kind]]) {
    stopifnot(identical(rownames(basis$features), obs$observation_id), all(is.finite(basis$features)))
  }
}
# Kernel construction uses only the individually observed inputs; no combination
# outcomes. All known single profiles are available, including excluded partners.
input_paths <- c(links_path, basis_path, script_path)
write_table(data.table(path = input_paths,
  sha256 = vapply(input_paths, digest, character(1), file = TRUE, algo = 'sha256')), 'input_manifest')

# Preserve64's objective and weighting so the control must reproduce exactly.
balanced_weights <- function(md) {
  group <- interaction(md$anchor, md$partner, md$biological_replicate, drop = TRUE)
  w <- 1 / as.numeric(table(group)[as.character(group)])
  w / sum(w)
}
center_context <- function(x, md) {
  group <- interaction(md$anchor, md$nominal_A_uM, md$biological_replicate, drop = TRUE)
  result <- x
  for (idx in split(seq_len(nrow(md)), group)) {
    if (is.matrix(x)) {
      result[idx, ] <- sweep(x[idx, , drop = FALSE], 2, colMeans(x[idx, , drop = FALSE]), '-')
    } else {
      result[idx] <- x[idx] - mean(x[idx])
    }
  }
  result
}
fit_path <- function(x, y, md) {
  w <- balanced_weights(md)
  augmented_x <- rbind(x, center_context(x, md))
  augmented_y <- c(y, center_context(y, md))
  fit <- glmnet(augmented_x, augmented_y, alpha = 0, lambda = lambda_grid,
    weights = c(w * .5, w * .5), standardize = FALSE, intercept = FALSE,
    thresh = 1e-8, maxit = 100000L)
  stopifnot(fit$jerr == 0L)
  fit
}
predict_path <- function(fit, x) as.matrix(predict(fit, newx = x, s = lambda_grid))
objective <- function(y, prediction, md) {
  w <- balanced_weights(md)
  .5 * colSums((prediction - y)^2 * w) +
    .5 * colSums((center_context(prediction, md) - center_context(y, md))^2 * w)
}

# Reuse fit paths for identical training subsets. Every lookup key encodes the
# complete permitted partner set; no outer held outcome enters a training path.
path_cache <- new.env(parent = emptyenv())
control_cache <- new.env(parent = emptyenv())
residual_cache <- new.env(parent = emptyenv())
subset_key <- function(training_partners) paste(sort(training_partners), collapse = '__')
control_paths <- function(training_partners, width) {
  key <- paste(subset_key(training_partners), width, sep = '::')
  if (exists(key, path_cache, inherits = FALSE)) return(path_cache[[key]])
  idx <- which(obs$partner %chin% training_partners)
  x <- bases$functional[[width]]$features
  fit <- fit_path(x[idx, , drop = FALSE], target[idx], obs[idx])
  result <- list(fit = fit, prediction = predict_path(fit, x))
  path_cache[[key]] <- result
  result
}
control_model <- function(training_partners) {
  key <- subset_key(training_partners)
  if (exists(key, control_cache, inherits = FALSE)) return(control_cache[[key]])
  stopifnot(length(training_partners) >= 3L)
  loss <- array(NA_real_, c(length(training_partners), length(lambda_grid), 3L))
  zero_loss <- numeric(length(training_partners))
  for (i in seq_along(training_partners)) {
    held <- training_partners[i]
    idx <- which(obs$partner == held)
    for (j in 1:3) {
      path <- control_paths(setdiff(training_partners, held), j)
      loss[i, , j] <- objective(target[idx], path$prediction[idx, , drop = FALSE], obs[idx])
    }
    zero_loss[i] <- objective(target[idx], matrix(0, length(idx), 1), obs[idx])
  }
  means <- apply(loss, c(2, 3), mean)
  choice <- which(means == min(means), arr.ind = TRUE)[1, ]
  zero <- mean(zero_loss) <= min(means)
  path <- control_paths(training_partners, choice[2])
  result <- list(prediction = if (zero) rep(0, nrow(obs)) else path$prediction[, choice[1]],
    fit = path$fit, width_index = choice[2], lambda_index = choice[1], zero = zero,
    training_partners = sort(training_partners))
  control_cache[[key]] <- result
  result
}
crossfit_residual <- function(training_partners) {
  key <- subset_key(training_partners)
  if (exists(key, residual_cache, inherits = FALSE)) return(residual_cache[[key]])
  stopifnot(length(training_partners) >= 4L)
  residual <- rep(NA_real_, nrow(obs))
  for (held in training_partners) {
    model <- control_model(setdiff(training_partners, held))
    stopifnot(!held %chin% model$training_partners)
    idx <- which(obs$partner == held)
    residual[idx] <- target[idx] - model$prediction[idx]
  }
  stopifnot(all(is.finite(residual[obs$partner %chin% training_partners])),
    all(is.na(residual[!obs$partner %chin% training_partners])))
  residual_cache[[key]] <- residual
  residual
}

predictions <- selections <- audits <- list()
for (held in partners) {
  status(paste('Exclude', held, '| fit protected functional model and cross-fitted RNA correction'))
  train_partners <- setdiff(partners, held)
  tr <- which(obs$partner %chin% train_partners)
  te <- which(obs$partner == held)
  source_ids <- function(x) unique(unlist(strsplit(x, ';', fixed = TRUE)))
  overlap <- intersect(source_ids(obs$source_sample_ids_AB[tr]), source_ids(obs$source_sample_ids_AB[te]))
  stopifnot(length(overlap) == 0L)
  control <- control_model(train_partners)
  train_residual <- crossfit_residual(train_partners)
  loss <- array(NA_real_, c(5L, length(lambda_grid), 3L))
  zero_loss <- numeric(5L)
  for (i in seq_along(train_partners)) {
    inner_held <- train_partners[i]
    inner_partners <- setdiff(train_partners, inner_held)
    inner_tr <- which(obs$partner %chin% inner_partners)
    inner_te <- which(obs$partner == inner_held)
    inner_residual <- crossfit_residual(inner_partners)
    inner_control <- control_model(inner_partners)
    stopifnot(!held %chin% inner_control$training_partners,
      !inner_held %chin% inner_control$training_partners)
    for (j in 1:3) {
      x <- bases$RNA[[j]]$features
      fit <- fit_path(x[inner_tr, , drop = FALSE], inner_residual[inner_tr], obs[inner_tr])
      candidate <- predict_path(fit, x[inner_te, , drop = FALSE]) + inner_control$prediction[inner_te]
      loss[i, , j] <- objective(target[inner_te], candidate, obs[inner_te])
    }
    zero_loss[i] <- objective(target[inner_te], matrix(inner_control$prediction[inner_te], ncol = 1), obs[inner_te])
  }
  means <- apply(loss, c(2, 3), mean)
  choice <- which(means == min(means), arr.ind = TRUE)[1, ]
  zero <- mean(zero_loss) <= min(means)
  x <- bases$RNA[[choice[2]]]$features
  fit <- fit_path(x[tr, , drop = FALSE], train_residual[tr], obs[tr])
  correction <- if (zero) rep(0, nrow(obs)) else predict_path(fit, x)[, choice[1]]
  model_path <- file.path(root, 'models', paste0('exclude_', held, '.qs'))
  qsave(list(control = control, correction_fit = fit, correction_width = choice[2],
    correction_lambda_index = choice[1], zero_correction = zero,
    training_observation_ids = obs$observation_id[tr], training_partners = train_partners,
    crossfitted_residual = train_residual[tr], basis_path = basis_path), model_path)
  saved <- qread(model_path)
  replay <- if (saved$zero_correction) rep(0, length(te)) else
    predict_path(saved$correction_fit, x[te, , drop = FALSE])[, saved$correction_lambda_index]
  replay_error <- max(abs(replay - correction[te]))
  stopifnot(replay_error < 1e-10, all(is.finite(correction)))
  out <- obs[, .(observation_id, condition_id, anchor, partner, nominal_A_uM,
    nominal_B_uM, biological_replicate)]
  out[, `:=`(fold = held, split = fifelse(partner == held, 'held', 'training'),
    control_g = control$prediction, molecular_correction = correction,
    predicted_g = control$prediction + correction,
    individual_log2 = log2(obs$ratio_A * obs$ratio_B))]
  # Freeze predictions before joining the outer combination outcomes.
  prediction_name <- paste0('prediction_exclude_', held)
  write_table(out, prediction_name)
  prediction_hash <- digest(file.path(root, 'data', paste0(prediction_name, '.tsv.gz')),
    file = TRUE, algo = 'sha256')
  out[, `:=`(observed_g = target, ratio_AB = obs$ratio_AB)]
  predictions[[held]] <- out
  selections[[held]] <- data.table(held_partner = held, zero_correction = zero,
    correction_lambda = if (zero) Inf else lambda_grid[choice[1]],
    correction_width = bases$width_factors[choice[2]], control_zero = control$zero,
    control_lambda = lambda_grid[control$lambda_index],
    control_width = bases$width_factors[control$width_index])
  audits[[held]] <- data.table(held_partner = held, source_overlap = length(overlap),
    replay_error = replay_error, prediction_sha256 = prediction_hash,
    training_pairs = uniqueN(obs[tr, .(anchor, partner)]), held_observations = length(te))
  grid <- data.table(lambda = rep(lambda_grid, 3), width = rep(bases$width_factors, each = length(lambda_grid)),
    internal_objective = as.vector(means), zero_objective = mean(zero_loss))
  write_table(grid, paste0('tuning_exclude_', held))
  write_table(rbindlist(predictions), 'predictions')
  write_table(rbindlist(selections), 'selection')
  write_table(rbindlist(audits), 'validation')
  message('Saved ', held, '; zero molecular correction = ', zero)
}

status('All six folds saved; verify control reproduction and produce diagnostics')
z <- rbindlist(predictions)
stopifnot(nrow(z[split == 'held']) == 496L, !anyDuplicated(z[split == 'held', observation_id]))
old <- as.data.table(read.delim(gzfile(file.path(previous, 'data/predictions.tsv.gz'))))[
  time_h == 3 & fold_type == 'single_partner' & model == 'Nonlinear_control' & split == 'held']
check <- merge(z[split == 'held', .(observation_id, control_g)],
  old[, .(observation_id, predicted_interaction_log2)], by = 'observation_id')
stopifnot(nrow(check) == 496L, max(abs(check$control_g - check$predicted_interaction_log2)) < 1e-8)
safe_cor <- function(x, y) {
  if (length(x) < 3L || sd(x) < 1e-10 || sd(y) < 1e-10) return(NA_real_)
  cor(x, y)
}
means <- z[, .(observed_g = mean(observed_g), control_g = mean(control_g),
  molecular_correction = mean(molecular_correction), predicted_g = mean(predicted_g),
  observed_count = mean(ratio_AB) * 100,
  control_count = mean(2^(individual_log2 + control_g)) * 100,
  predicted_count = mean(2^(individual_log2 + predicted_g)) * 100),
  by = .(split, condition_id, anchor, partner, nominal_A_uM, nominal_B_uM)]
means[, remaining_departure := observed_g - control_g]
summary <- means[, .(conditions = .N, control_interaction_r = safe_cor(observed_g, control_g),
  molecular_interaction_r = safe_cor(observed_g, predicted_g),
  correction_r = safe_cor(remaining_departure, molecular_correction),
  control_count_r = safe_cor(observed_count, control_count),
  molecular_count_r = safe_cor(observed_count, predicted_count)), by = split]
by_partner <- means[split == 'held', .(conditions = .N,
  control_interaction_r = safe_cor(observed_g, control_g),
  molecular_interaction_r = safe_cor(observed_g, predicted_g),
  correction_r = safe_cor(remaining_departure, molecular_correction)), by = partner]
write_table(means, 'figure_condition_predictions')
write_table(summary, 'correlation_summary')
write_table(by_partner, 'held_correlations_by_partner')
by_anchor <- means[split == 'held', .(conditions = .N,
  control_interaction_r = safe_cor(observed_g, control_g),
  molecular_interaction_r = safe_cor(observed_g, predicted_g)), by = anchor]
by_pair <- means[split == 'held', .(conditions = .N,
  observed_mean_departure = mean(observed_g), control_mean_departure = mean(control_g),
  molecular_mean_departure = mean(predicted_g),
  control_interaction_r = safe_cor(observed_g, control_g),
  molecular_interaction_r = safe_cor(observed_g, predicted_g)), by = .(anchor, partner)]
write_table(by_anchor, 'held_correlations_by_anchor')
write_table(by_pair, 'held_correlations_by_pair')
print(summary)
print(by_partner)

theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom'))
save_plot <- function(p, name, width, height) {
  ggsave(file.path(root, 'figures', paste0(name, '.pdf')), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(root, 'figures', paste0(name, '.png')), p, width = width, height = height, dpi = 140)
}
d <- melt(means, id.vars = c('split', 'condition_id', 'anchor', 'partner', 'observed_g'),
  measure.vars = c('control_g', 'predicted_g'), variable.name = 'method', value.name = 'prediction')
d[, method := factor(method, levels = c('control_g', 'predicted_g'),
  labels = c('Function + dose', 'Function + dose + RNA'))]
d[, split := factor(split, levels = c('training', 'held'))]
ann <- d[, .(label = sprintf('Pearson r = %.2f', safe_cor(observed_g, prediction))), by = .(split, method)]
p <- ggplot(d, aes(observed_g, prediction)) + geom_point(aes(color = anchor), alpha = .65, size = 1.4) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.05, vjust = 1.3) +
  facet_wrap(vars(split, method), ncol = 2, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'Does RNA improve prediction of combination departures?',
    subtitle = 'Negative departure: fewer cells than expected from measured individual effects',
    x = 'Observed log2[AB / (A × B)]', y = 'Predicted log2 departure', color = 'Anchor',
    caption = 'Held: all combinations of one partner excluded under both anchors. Known singles allowed.\n198 conditions; biological repeats averaged. Reused outcomes: retrospective development.')
save_plot(p, 'interaction_training_vs_held', 8.5, 8)
held_d <- d[split == 'held']
ann <- held_d[, .(label = sprintf('r = %.2f', safe_cor(observed_g, prediction))), by = .(partner, method)]
p <- ggplot(held_d, aes(observed_g, prediction)) + geom_point(aes(color = anchor), size = 1.2) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .45) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 3) +
  facet_wrap(vars(method, partner), ncol = 6, scales = 'free') + theme(aspect.ratio = 1, strip.text = element_text(size = 8)) +
  labs(title = 'Every excluded partner: actual interaction departures',
    x = 'Observed log2 departure', y = 'Predicted log2 departure', color = 'Anchor')
save_plot(p, 'held_interactions_each_partner', 14, 5.8)
ann <- held_d[, .(label = sprintf('Pearson r = %.2f', safe_cor(observed_g, prediction))), by = .(anchor, method)]
p <- ggplot(held_d, aes(observed_g, prediction)) + geom_point(aes(color = partner), size = 1.5, alpha = .8) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2, size = 3.2) +
  facet_wrap(vars(anchor, method), ncol = 2, scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'Locate the failure: evaluate each anchor separately',
    subtitle = 'All predictions exclude the indicated partner from combination training',
    x = 'Observed log2 departure', y = 'Predicted log2 departure', color = 'Excluded partner')
save_plot(p, 'held_interactions_each_anchor', 8.5, 8)
ann <- means[, .(label = sprintf('Pearson r = %.2f', safe_cor(remaining_departure, molecular_correction))), by = split]
p <- ggplot(means, aes(remaining_departure, molecular_correction)) + geom_point(aes(color = anchor), alpha = .65) +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'black', linewidth = .5) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -.04, vjust = 1.2) +
  facet_wrap(vars(split), scales = 'free') + theme(aspect.ratio = 1) +
  labs(title = 'What information does the molecular correction add?',
    x = 'Observed departure minus function–dose prediction', y = 'RNA correction (log2 units)', color = 'Anchor')
save_plot(p, 'molecular_correction_vs_remaining_departure', 8.5, 4.8)
curves <- melt(means[split == 'held'],
  id.vars = c('anchor', 'partner', 'nominal_A_uM', 'nominal_B_uM'),
  measure.vars = c('observed_g', 'control_g', 'predicted_g'), variable.name = 'source', value.name = 'departure')
curves[, source := factor(source, levels = c('observed_g', 'control_g', 'predicted_g'),
  labels = c('Observed', 'Function + dose', 'Function + dose + RNA'))]
pdf(file.path(root, 'figures/all12_held_interaction_dose_curves.pdf'), width = 9, height = 6)
for (aa in unique(curves$anchor)) {
  for (pp in partners) {
    dd <- curves[anchor == aa & partner == pp]
    p <- ggplot(dd, aes(nominal_A_uM, departure, color = source, group = source)) +
      geom_hline(yintercept = 0, color = 'grey75') + geom_line(linewidth = .6) + geom_point(size = 1.7) +
      scale_x_log10() + facet_wrap(vars(nominal_B_uM), labeller = label_both, ncol = 3) +
      scale_color_manual(values = c('Observed' = 'black', 'Function + dose' = '#009E73', 'Function + dose + RNA' = '#D55E00')) +
      labs(title = paste(aa, '+', pp, '— partner excluded from training'),
        x = paste(aa, 'dose (µM)'), y = 'Log2 combination departure', color = NULL)
    print(p)
  }
}
dev.off()
capture.output(sessionInfo(), file = file.path(root, 'sessionInfo.txt'))
writeLines(system2('git', c('-C', '/work', 'rev-parse', 'HEAD'), stdout = TRUE), file.path(root, 'git_revision.txt'))
report <- c('Version70: protected functional model plus cross-fitted molecular interaction correction.',
  'Original40, BT20 3h RNA / 72h imaged cell counts. Measured individual effects supplied.',
  'No frozen-f refit. No combination RNA. Complete outer-partner exclusion.',
  'Functional comparator reproduces64 to numerical tolerance; all predictions saved before outcome joining.',
  capture.output(print(summary)), capture.output(print(by_partner)),
  capture.output(print(by_anchor)),
  'Correlations describe198 dose conditions from12 pairs in one cell, not198 independent mechanisms.',
  'All outcomes were previously inspected: retrospective model development, not independent validation.')
writeLines(report, file.path(root, 'REPORT.txt'))
status('COMPLETE: all six folds, checks, tables and diagnostic figures saved')
