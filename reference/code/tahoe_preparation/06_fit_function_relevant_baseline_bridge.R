#!/usr/bin/env Rscript

## Fit the accepted global bilinear conditioner from five baseline programs
## and the nine-program reference footprint to the measured molecular target.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
})
root <- normalizePath(Sys.getenv("TAHOE_V2_ROOT", "tahoe_prism_v2"),
  mustWork = TRUE)
source(file.path(root, "config.R"))
source(file.path(root, "scripts", "helpers.R"))
source(file.path(WORK_ROOT, "src", "R", "lib_heatmap.R"))
RESULT <- stage_result("06_baseline_bridge")
ARTIFACT <- stage_artifact("06_baseline_bridge")
dir.create(RESULT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTIFACT, recursive = TRUE, showWarnings = FALSE)

plot_stage06 <- function() {
  theme_stage06 <- function(base_size = 8.5) theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
      plot.title = element_text(face = "bold"))
  square_scatter_theme <- function(base_size = 8.5)
    theme_stage06(base_size) + theme(aspect.ratio = 1)
  save_gg <- function(filename, plot, width, height) ggsave(
    file.path(RESULT, filename), plot, width = width, height = height,
    device = cairo_pdf, bg = "white", limitsize = FALSE)
  save_heatmap_pdf <- function(filename, heatmap, width, height) {
    grDevices::pdf(file.path(RESULT, filename), width = width, height = height,
      useDingbats = FALSE)
    ComplexHeatmap::draw(heatmap)
    grDevices::dev.off()
  }
  fmt <- function(x, digits = 3L) ifelse(is.finite(x),
    formatC(x, digits = digits, format = "f"), "NA")
  source_labels <- c(tahoe_dmso = "Tahoe DMSO", depmap_static = "DepMap")
  source_colors <- c(tahoe_dmso = "#B2182B", depmap_static = "#2166AC")

  selection <- fread(file.path(RESULT, "baseline_bridge_selection.tsv"))
  for (source_name in names(source_labels)) {
    current <- selection[baseline_source == source_name &
      record_type == "aggregate"]
    zero <- current[candidate == "zero_context", mean_validation_mse][1L]
    finite <- current[candidate == "nonzero_ridge"]
    p_selection <- ggplot(finite, aes(log10(lambda),
        mean_validation_mse)) +
      geom_hline(yintercept = zero, linetype = 2, color = "#4D4D4D") +
      geom_line(linewidth = 0.65, color = source_colors[[source_name]]) +
      geom_point(size = 0.9, color = source_colors[[source_name]]) +
      geom_point(data = finite[selected == TRUE], size = 3,
        shape = 21, fill = "gold", color = "black") +
      labs(x = expression(log[10](lambda)),
        y = "Mean condition-balanced validation MSE",
        title = paste(source_labels[[source_name]], "baseline penalty path"),
        subtitle = "Dashed line: exact zero-context model; gold: selected penalty") +
      theme_stage06()
    save_gg(sprintf("06a_%s_baseline_model_selection.pdf", source_name),
      p_selection, 5.8, 4.2)
  }

  coefficients <- fread(file.path(RESULT,
    "baseline_bridge_coefficients.tsv"))
  coefficients[, c("baseline_program", "response_program") := tstrsplit(
    feature, "_x_", fixed = TRUE)]
  baseline_programs <- paste0("B", 1:5)
  response_programs <- c("P1", "P2", "P3", "P4", "P5", "P7", "P8",
    "P9", "P10")
  for (source_name in names(source_labels)) {
    matrix_dt <- dcast(coefficients[baseline_source == source_name],
      baseline_program ~ response_program, value.var = "coefficient", fill = 0)
    matrix_dt <- matrix_dt[match(baseline_programs, baseline_program)]
    coefficient_matrix <- as.matrix(matrix_dt[, ..response_programs])
    rownames(coefficient_matrix) <- matrix_dt$baseline_program
    limit <- max(abs(coefficient_matrix))
    heatmap <- draw_heatmap(coefficient_matrix, scale_rows = FALSE,
      value_range = c(-limit, 0, limit),
      value_colors = c("#2166AC", "white", "#B2182B"),
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = TRUE, show_column_names = TRUE,
      title = paste(source_labels[[source_name]],
        "baseline × footprint coefficients"))
    save_heatmap_pdf(sprintf("06b_%s_conditioner_coefficients.pdf",
      source_name), heatmap, 5.8, 4.4)
  }

  molecular_predictions <- fread(file.path(ARTIFACT,
    "baseline_molecular_predictions.tsv.gz"))
  activity_ids <- fread(file.path(stage_result("04_activity_panel"),
    "condition_activity_census.tsv"))[activity_supported == TRUE,
      perturbation_id]
  for (source_name in names(source_labels)) {
    plot_data <- molecular_predictions[baseline_source == source_name &
      perturbation_id %chin% activity_ids]
    scatter <- ggplot(plot_data, aes(predicted_z, measured_z)) +
      geom_point(alpha = 0.35, size = 0.8,
        color = source_colors[[source_name]]) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
        linewidth = 0.7) +
      labs(title = paste(source_labels[[source_name]],
          "baseline prediction of molecular z"),
        subtitle = sprintf("Pearson r = %s; predicted/observed SD = %s",
          fmt(safe_cor(plot_data$predicted_z, plot_data$measured_z)),
          fmt(sd(plot_data$predicted_z) / sd(plot_data$measured_z))),
        x = "Baseline-predicted z", y = "Measured z") +
      square_scatter_theme()
    save_gg(sprintf("06c_%s_baseline_predicted_vs_measured_score.pdf",
      source_name), scatter, 4.1, 4.1)
  }

  prediction_wide <- dcast(molecular_predictions[
    perturbation_id %chin% activity_ids],
    observation_id + measured_z ~ baseline_source,
    value.var = "predicted_z")
  p_cross <- ggplot(prediction_wide, aes(tahoe_dmso, depmap_static)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_point(alpha = 0.18, size = 0.42, color = "#6A51A3") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
      linewidth = 0.55, color = "black") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
      label = sprintf("r=%s; n=%s", fmt(safe_cor(
        prediction_wide$tahoe_dmso, prediction_wide$depmap_static)),
        format(nrow(prediction_wide), big.mark = ",")), size = 2.5) +
    labs(x = "Tahoe DMSO predicted z", y = "DepMap predicted z",
      title = "Cross-source molecular prediction concordance") +
    square_scatter_theme()
  save_gg("06d_tahoe_depmap_prediction_concordance.pdf", p_cross, 4.6, 4.4)

  performance <- fread(file.path(RESULT,
    "baseline_molecular_performance.tsv"))[
      cohort == "all_activity_supported"]
  performance[, source_label := source_labels[baseline_source]]
  performance_long <- melt(performance,
    id.vars = c("baseline_source", "source_label"),
    measure.vars = c("pearson", "partial_r2",
      "predicted_to_observed_sd_ratio"), variable.name = "metric",
    value.name = "value")
  performance_long[, metric := factor(metric,
    levels = c("pearson", "partial_r2", "predicted_to_observed_sd_ratio"),
    labels = c("Pearson r", "Partial R²", "Predicted / observed SD"))]
  p_performance <- ggplot(performance_long,
    aes(source_label, value, fill = baseline_source)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_col(width = 0.65) +
    geom_text(aes(label = fmt(value)), vjust = -0.25, size = 2.2) +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = source_colors, guide = "none") +
    labs(x = NULL, y = "Performance",
      title = "Baseline prediction of the function-relevant molecular target") +
    theme_stage06()
  save_gg("06e_baseline_molecular_performance.pdf", p_performance, 8.4, 3.8)
}

if ("--plots-only" %in% commandArgs(trailingOnly = TRUE)) {
  plot_stage06()
  message("Stage 06 plots refreshed from existing baseline-bridge artifacts.")
  quit(save = "no", status = 0L)
}

standardize_baseline <- function(raw, training_lines, scoring_lines) {
  center <- colMeans(raw[training_lines, , drop = FALSE])
  scale <- apply(raw[training_lines, , drop = FALSE], 2L, sd)
  if (any(!is.finite(center)) || any(!is.finite(scale)) || any(scale <= 0))
    fail("A baseline program cannot be standardized training-only.")
  scores <- sweep(sweep(raw[scoring_lines, , drop = FALSE], 2L,
    center, "-"), 2L, scale, "/")
  rownames(scores) <- scoring_lines
  list(scores = scores, center = center, scale = scale)
}

make_bilinear_design <- function(row_lines, perturbations, baseline_scores,
                                 footprint, programs) {
  B <- baseline_scores[row_lines, , drop = FALSE]
  Dbar <- footprint[perturbations, programs, drop = FALSE]
  pieces <- lapply(seq_len(ncol(B)), function(index) Dbar * B[, index])
  x <- do.call(cbind, pieces)
  colnames(x) <- unlist(lapply(colnames(B), function(b)
    paste0(b, "_x_", programs)), use.names = FALSE)
  if (any(!is.finite(x))) fail("The bilinear design is nonfinite.")
  x
}

baseline_path <- function(x_train, y, x_test, perturbations) {
  weights <- balanced_weights(perturbations)
  design <- weighted_standardize(x_train, x_test, weights)
  root_weight <- sqrt(weights)
  wx <- design$train * root_weight; wy <- y * root_weight
  decomposition <- svd(wx, nu = min(dim(wx)), nv = min(dim(wx)))
  uy <- as.numeric(crossprod(decomposition$u, wy))
  cross <- design$test %*% decomposition$v
  result <- matrix(0, nrow(x_test), length(BASELINE_LAMBDA_GRID))
  for (i in seq_along(BASELINE_LAMBDA_GRID)) result[, i] <- as.numeric(
    cross %*% (uy * decomposition$d /
      (decomposition$d ^ 2 + BASELINE_LAMBDA_GRID[i])))
  result
}

fit_baseline <- function(x, y, perturbations, lambda) {
  weights <- balanced_weights(perturbations)
  design <- weighted_standardize(x, x, weights)
  root_weight <- sqrt(weights); wx <- design$train * root_weight
  decomposition <- svd(wx, nu = min(dim(wx)), nv = min(dim(wx)))
  beta <- decomposition$v %*% (as.numeric(crossprod(decomposition$u,
    y * root_weight)) * decomposition$d /
      (decomposition$d ^ 2 + lambda))
  list(lambda = lambda, center = design$center, scale = design$scale,
    keep = design$keep, feature_names = design$feature_names,
    beta_standardized = as.numeric(beta))
}

predict_baseline <- function(fit, x) {
  if (is.null(fit)) return(rep(0, nrow(x)))
  z <- sweep(sweep(x[, fit$keep, drop = FALSE], 2L,
    fit$center[fit$keep], "-"), 2L, fit$scale[fit$keep], "/")
  as.numeric(z %*% fit$beta_standardized)
}

baseline_raw_beta <- function(fit) {
  answer <- setNames(rep(0, length(fit$feature_names)), fit$feature_names)
  answer[fit$feature_names[fit$keep]] <- fit$beta_standardized /
    fit$scale[fit$keep]
  answer
}

make_sub_assignment <- function(lines, salt) {
  set.seed(SEED + 100L + salt)
  shuffled <- sample(sort(unique(lines)), length(unique(lines)), FALSE)
  setNames(rep(seq_len(min(5L, length(lines))),
    length.out = length(shuffled)), shuffled)
}

choose_decoder <- function(x, y, lines, salt) {
  assignment <- make_sub_assignment(unique(lines), salt)
  losses <- matrix(NA_real_, max(assignment), length(DECODER_LAMBDA_GRID))
  for (fold in seq_len(max(assignment))) {
    validation_lines <- names(assignment)[assignment == fold]
    train <- !lines %chin% validation_lines
    valid <- lines %chin% validation_lines
    path <- ridge_path(x[train, , drop = FALSE], y[train],
      x[valid, , drop = FALSE], DECODER_LAMBDA_GRID)
    losses[fold, ] <- colMeans(sweep(path, 1L, y[valid], "-") ^ 2)
  }
  index <- which.min(colMeans(losses))
  fit <- fit_ridge(x, y, DECODER_LAMBDA_GRID[index])
  list(fit = fit, raw = raw_ridge_coefficients(fit),
    lambda = DECODER_LAMBDA_GRID[index])
}

make_footprint <- function(D, score_meta, lines, programs) {
  use <- score_meta$depmap_id %chin% lines
  means <- data.table(perturbation_id = score_meta$perturbation_id[use],
    D[use, programs, drop = FALSE])[, lapply(.SD, mean),
    by = perturbation_id, .SDcols = programs]
  result <- as.matrix(means[, ..programs]); rownames(result) <- means$perturbation_id
  result
}

make_measured_z <- function(D, score_meta, lines, decoder, programs) {
  footprint <- make_footprint(D, score_meta, lines, programs)
  position <- match(score_meta$perturbation_id, rownames(footprint))
  z <- as.numeric((D[, programs, drop = FALSE] -
    footprint[position, programs, drop = FALSE]) %*% decoder$beta[programs])
  list(z = z, footprint = footprint)
}

## Frozen response scores and functional rows regenerated in this v2 run.
response_stability <- fread(file.path(stage_result("02_response_programs"),
  "fold_program_stability.tsv"))
programs <- response_stability[outer_fold == OUTER_FOLD & stable == TRUE,
  unique(program)]
programs <- paste0("P", sort(as.integer(sub("P", "", programs))))
score_long <- fread(file.path(stage_artifact("02_response_programs"),
  "fold_program_scores.tsv.gz"))[program %chin% programs]
score_meta <- dcast(score_long, observation_id + depmap_id + cell_line +
  perturbation_id + drug + dose + dose_panel + role ~ program,
  value.var = "score")
D <- as.matrix(score_meta[, ..programs]); rownames(D) <- score_meta$observation_id
function_rows <- fread(file.path(stage_artifact("01_inputs"),
  "function_outcomes.tsv.gz"))
function_rows[, `:=`(observation_id = paste(depmap_id, perturbation_id,
  sep = "__"), drug = tahoe_drug, dose = tahoe_dose_uM)]
function_position <- match(function_rows$observation_id,
  score_meta$observation_id)
chain <- function_rows$chain_observed == 1L & !is.na(function_position)
folds <- fread(file.path(stage_result("01_inputs"), "fold_assignments.tsv"))
training_lines <- folds[outer_fold != OUTER_FOLD, depmap_id]
held_lines <- folds[outer_fold == OUTER_FOLD, depmap_id]
all_lines <- sort(folds$depmap_id)
inner_assignment <- make_inner_assignment(training_lines)
active_ids <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))[activity_supported == TRUE,
  perturbation_id]

decoder_coefficients <- fread(file.path(stage_result(
  "05_program_function_bridge"), "function_decoder_coefficients.tsv"))
outer_decoder <- list(intercept = decoder_coefficients[
  term == "(Intercept)", coefficient],
  beta = setNames(decoder_coefficients[term != "(Intercept)", coefficient],
    decoder_coefficients[term != "(Intercept)", term]))
outer_footprint_dt <- fread(file.path(stage_result(
  "05_program_function_bridge"), "reference_perturbation_footprints.tsv.gz"))
outer_footprint <- as.matrix(outer_footprint_dt[, ..programs])
rownames(outer_footprint) <- outer_footprint_dt$perturbation_id
outer_measured <- fread(file.path(stage_artifact(
  "05_program_function_bridge"), "measured_function_relevant_scores.tsv.gz"))
outer_z <- outer_measured$measured_z[match(score_meta$observation_id,
  outer_measured$observation_id)]

baseline_stability <- fread(file.path(stage_result("03_baseline_programs"),
  "baseline_program_stability.tsv"))
baseline_programs <- baseline_stability[outer_fold == OUTER_FOLD &
  stable == TRUE, unique(program)]
baseline_programs <- paste0("B", sort(as.integer(sub("B", "",
  baseline_programs))))
if (length(baseline_programs) != EXPECTED_BASELINE_PROGRAMS)
  fail("The reproduced baseline bank is not five stable programs.")
baseline_long <- fread(file.path(stage_artifact("03_baseline_programs"),
  "baseline_program_scores.tsv.gz"))
make_baseline_matrix <- function(source, aggregation) {
  x <- dcast(baseline_long[baseline_source == source &
    aggregation_level == aggregation & program %chin% baseline_programs],
    depmap_id ~ program, value.var = "raw_program_score")
  result <- as.matrix(x[, ..baseline_programs]); rownames(result) <- x$depmap_id
  result[all_lines, baseline_programs, drop = FALSE]
}
baseline_raw <- list(tahoe_dmso = make_baseline_matrix("tahoe_dmso",
  "line_mean"), depmap_static = make_baseline_matrix("depmap_static",
  "line_level"))

## One intentionally active outer-training inner fold, with its own decoder,
## reference footprints and baseline scaling.
cv_parts <- list()
for (inner_fold in seq_len(ACTIVE_BASELINE_INNER_FOLDS)) {
  validation_lines <- names(inner_assignment)[inner_assignment == inner_fold]
  learn_lines <- names(inner_assignment)[inner_assignment != inner_fold]
  decoder_rows <- function_rows$depmap_id %chin% learn_lines & chain
  local_decoder <- choose_decoder(
    D[function_position[decoder_rows], , drop = FALSE],
    function_rows$prism_sensitization[decoder_rows],
    function_rows$depmap_id[decoder_rows], inner_fold)
  local_z <- make_measured_z(D, score_meta, learn_lines,
    local_decoder$raw, programs)
  for (source in names(baseline_raw)) {
    standardized <- standardize_baseline(baseline_raw[[source]],
      learn_lines, training_lines)$scores
    train <- score_meta$depmap_id %chin% learn_lines &
      score_meta$perturbation_id %chin% active_ids
    valid <- score_meta$depmap_id %chin% validation_lines &
      score_meta$perturbation_id %chin% active_ids
    x_train <- make_bilinear_design(score_meta$depmap_id[train],
      score_meta$perturbation_id[train], standardized,
      local_z$footprint, programs)
    x_valid <- make_bilinear_design(score_meta$depmap_id[valid],
      score_meta$perturbation_id[valid], standardized,
      local_z$footprint, programs)
    path <- baseline_path(x_train, local_z$z[train], x_valid,
      score_meta$perturbation_id[train])
    predictions <- cbind(0, path)
    candidates <- c(Inf, BASELINE_LAMBDA_GRID)
    cv_parts[[length(cv_parts) + 1L]] <- data.table(
      baseline_source = source, inner_fold = inner_fold,
      candidate = c("zero_context", rep("nonzero_ridge",
        length(BASELINE_LAMBDA_GRID))), lambda = candidates,
      validation_mse = vapply(seq_along(candidates), function(index) {
        w <- balanced_weights(score_meta$perturbation_id[valid])
        sum(w * (local_z$z[valid] - predictions[, index]) ^ 2) / sum(w)
      }, numeric(1L)),
      validation_pearson = vapply(seq_along(candidates), function(index)
        safe_cor(local_z$z[valid], predictions[, index]), numeric(1L)),
      n_training_lines = length(learn_lines),
      n_validation_lines = length(validation_lines),
      n_training_rows = sum(train), n_validation_rows = sum(valid),
      n_features = ncol(x_train), condition_equal_weighting = TRUE,
      outcome = "measured_function_relevant_molecular_score",
      prism_is_baseline_model_outcome = FALSE)
  }
}
cv <- rbindlist(cv_parts)
aggregate <- cv[, .(mean_validation_mse = mean(validation_mse),
  mean_validation_pearson = mean(validation_pearson, na.rm = TRUE),
  n_inner_folds = .N), by = .(baseline_source, candidate, lambda)]
selection <- aggregate[, {
  selected_index <- which.min(mean_validation_mse)
  best_nonzero_index <- which.min(ifelse(candidate == "nonzero_ridge",
    mean_validation_mse, Inf))
  copy(.SD)[, `:=`(selected = seq_len(.N) == selected_index,
    best_nonzero = seq_len(.N) == best_nonzero_index,
    selected_lambda = lambda[selected_index],
    best_nonzero_lambda = lambda[best_nonzero_index],
    zero_context_selected = is.infinite(lambda[selected_index]))]
}, by = baseline_source]
fwrite(rbindlist(list(cv[, `:=`(record_type = "inner_fold")],
  selection[, `:=`(record_type = "aggregate", inner_fold = NA_integer_,
    validation_mse = NA_real_, validation_pearson = NA_real_,
    n_training_lines = NA_integer_, n_validation_lines = NA_integer_,
    n_training_rows = NA_integer_, n_validation_rows = NA_integer_,
    n_features = 45L, condition_equal_weighting = TRUE,
    outcome = "measured_function_relevant_molecular_score",
    prism_is_baseline_model_outcome = FALSE)]), fill = TRUE),
  file.path(RESULT, "baseline_bridge_selection.tsv"), sep = "\t")

## Fit each source on all 26 training lines and predict the seven held lines.
held_score <- score_meta$depmap_id %chin% held_lines
train_score <- score_meta$depmap_id %chin% training_lines &
  score_meta$perturbation_id %chin% active_ids
held_function <- function_rows[outer_fold == OUTER_FOLD]
mean_function <- outer_decoder$intercept + as.numeric(outer_footprint[
  held_function$perturbation_id, programs, drop = FALSE] %*%
  outer_decoder$beta[programs])
measured_function <- fread(file.path(stage_artifact(
  "05_program_function_bridge"), "measured_function_predictions.tsv.gz"))
models <- list(); prediction_parts <- list(); function_parts <- list()
coefficient_parts <- list(); audit_parts <- list()
for (source in names(baseline_raw)) {
  standardized_object <- standardize_baseline(baseline_raw[[source]],
    training_lines, all_lines)
  x_train <- make_bilinear_design(score_meta$depmap_id[train_score],
    score_meta$perturbation_id[train_score], standardized_object$scores,
    outer_footprint, programs)
  selected_lambda <- selection[baseline_source == source & selected == TRUE,
    lambda]
  fit <- if (is.infinite(selected_lambda)) NULL else fit_baseline(x_train,
    outer_z[train_score], score_meta$perturbation_id[train_score],
    selected_lambda)
  models[[source]] <- list(fit = fit, baseline_center =
    standardized_object$center, baseline_scale = standardized_object$scale,
    baseline_programs = baseline_programs, response_programs = programs,
    footprint = outer_footprint)
  x_held <- make_bilinear_design(score_meta$depmap_id[held_score],
    score_meta$perturbation_id[held_score], standardized_object$scores,
    outer_footprint, programs)
  predicted_z <- predict_baseline(fit, x_held)
  prediction_parts[[source]] <- data.table(
    observation_id = score_meta$observation_id[held_score],
    depmap_id = score_meta$depmap_id[held_score],
    cell_line = score_meta$cell_line[held_score],
    perturbation_id = score_meta$perturbation_id[held_score],
    drug = score_meta$drug[held_score], dose = score_meta$dose[held_score],
    baseline_source = source, measured_z = outer_z[held_score],
    predicted_z = predicted_z)
  x_function <- make_bilinear_design(held_function$depmap_id,
    held_function$perturbation_id, standardized_object$scores,
    outer_footprint, programs)
  predicted_function_z <- predict_baseline(fit, x_function)
  function_parts[[source]] <- data.table(
    observation_id = held_function$observation_id,
    depmap_id = held_function$depmap_id, cell_line = held_function$cell_line,
    perturbation_id = held_function$perturbation_id,
    drug = held_function$drug, dose = held_function$dose,
    observed_function = held_function$prism_sensitization,
    evaluation_included = held_function$chain_observed == 1L,
    baseline_source = source,
    measured_z = measured_function$measured_z[match(
      held_function$observation_id, measured_function$observation_id)],
    predicted_z = predicted_function_z,
    mean_prediction = mean_function,
    predicted_function = mean_function + predicted_function_z,
    measured_prediction = measured_function$measured_prediction[match(
      held_function$observation_id, measured_function$observation_id)])
  raw_beta <- if (is.null(fit)) setNames(rep(0, ncol(x_train)),
    colnames(x_train)) else baseline_raw_beta(fit)
  coefficient_parts[[source]] <- data.table(baseline_source = source,
    feature = names(raw_beta), coefficient = as.numeric(raw_beta),
    lambda = selected_lambda, n_training_lines = length(training_lines),
    n_training_rows = sum(train_score), n_conditions = length(active_ids),
    no_intercept = TRUE)
  audit_parts[[source]] <- data.table(baseline_source = source,
    n_baseline_programs = ncol(baseline_raw[[source]]),
    n_response_programs = length(programs), n_features = ncol(x_train),
    training_only_scaling = TRUE, direct_baseline_to_prism = FALSE,
    perturbation_identity_feature = FALSE, per_drug_model = FALSE)
}
molecular_predictions <- rbindlist(prediction_parts)
function_predictions <- rbindlist(function_parts)
fwrite(molecular_predictions, file.path(ARTIFACT,
  "baseline_molecular_predictions.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(function_predictions, file.path(ARTIFACT,
  "baseline_function_predictions.tsv.gz"), sep = "\t", compress = "gzip")
saveRDS(models, file.path(ARTIFACT, "baseline_bridge_models.rds"))
fwrite(rbindlist(coefficient_parts),
  file.path(RESULT, "baseline_bridge_coefficients.tsv"), sep = "\t")
fwrite(rbindlist(audit_parts),
  file.path(RESULT, "baseline_bridge_audit.tsv"), sep = "\t")

activity <- fread(file.path(stage_result("04_activity_panel"),
  "condition_activity_census.tsv"))
molecular_predictions <- merge(molecular_predictions,
  activity[, .(perturbation_id, activity_supported,
    activity_supported_prism_overlap)], by = "perturbation_id",
  all.x = TRUE, sort = FALSE)
performance_parts <- list()
for (source in names(baseline_raw)) {
  all_active <- molecular_predictions[baseline_source == source &
    activity_supported == TRUE]
  function_panel <- function_predictions[baseline_source == source &
    evaluation_included == TRUE & perturbation_id %chin% activity[
      activity_supported_prism_overlap == TRUE, perturbation_id]]
  performance_parts[[length(performance_parts) + 1L]] <- cbind(
    data.table(baseline_source = source, cohort = "all_activity_supported"),
    endpoint_metrics(all_active, "measured_z", "predicted_z"))
  performance_parts[[length(performance_parts) + 1L]] <- cbind(
    data.table(baseline_source = source,
      cohort = "activity_supported_prism"),
    endpoint_metrics(function_panel, "measured_z", "predicted_z"))
}
performance <- rbindlist(performance_parts, fill = TRUE)
fwrite(performance,
  file.path(RESULT, "baseline_molecular_performance.tsv"), sep = "\t")
write_artifact_manifest(ARTIFACT,
  "scripts/06_fit_function_relevant_baseline_bridge.R")
writeLines("complete", file.path(RESULT, "stage06_complete.flag"))
plot_stage06()
message(sprintf("Stage 06 complete: Tahoe r=%.3f; DepMap r=%.3f.",
  performance[baseline_source == "tahoe_dmso" &
    cohort == "all_activity_supported", pearson],
  performance[baseline_source == "depmap_static" &
    cohort == "all_activity_supported", pearson]))
