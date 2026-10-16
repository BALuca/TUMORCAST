# ============================================================================
# Retained wrapper from src/R/lib_heatmap.R, with shared project diagnostics below.
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(matrixStats)
  library(RColorBrewer)
  library(viridisLite)
})

# Established signed-value palette shared by project heatmaps.
hm_value_colors <- c("#006aff", "#0068fe", "#0066fc", "#0063fa", "#0061f8", "#005ff6", "#005df4", "#005af2",
    "#0058f0", "#0056ee", "#0054ec", "#0052ea", "#004fe8", "#004de6", "#004be4", "#0049e2", "#0046e0",
    "#0045de", "#0044dd", "#0042db", "#0041d9", "#0040d8", "#003fd6", "#003ed4", "#003dd3", "#003cd1",
    "#003bcf", "#003ace", "#0039cc", "#0038ca", "#0037c9", "#0036c7", "#0034c5", "#0033c4", "#0131ba",
    "#0c2eae", "#132ca1", "#172995", "#192688", "#1b247c", "#1c2170", "#1c1f65", "#1b1c59", "#1a1a4e",
    "#191743", "#171539", "#15122e", "#130e25", "#100a1b", "#080511", "#000000", "#121205", "#1d1e0a",
    "#26290e", "#30360f", "#3b4210", "#464f0f", "#515c0e", "#5c6a0c", "#677807", "#738601", "#7f9500",
    "#8ba400", "#97b300", "#a4c200", "#b1d200", "#bde100", "#c6ec00", "#c8ee00", "#c9ef00", "#caf000",
    "#cbf100", "#ccf300", "#cef400", "#cff500", "#d0f600", "#d1f800", "#d2f900", "#d4fa00", "#d5fc00",
    "#d6fd00", "#d7fe00", "#d8ff00", "#daff00", "#dbff00", "#ddff00", "#dfff00", "#e0ff00", "#e2ff00",
    "#e4ff00", "#e6ff00", "#e7ff00", "#e9ff00", "#ebff00", "#edff02", "#eeff07", "#f0ff0b", "#f2ff0e",
    "#f4ff12", "#f5ff14", "#f7ff17")

# ------- Discrete palettes (self-contained) ----------------------------------
# draw_heatmap() uses map_discrete()/safe_factor() for annotations. These helpers
# are intentionally defined here so a fresh Rscript session can source this file
# without relying on other libs.

.pub_palettes <- list(
  set_base = c(
    brewer.pal(9, "Set1"),
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Paired")[c(1, 3, 5, 7, 9, 12, 10, 8, 6, 4, 2, 11)],
    brewer.pal(8, "Set3"),
    brewer.pal(8, "Accent")
  ),
  dataset = c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), brewer.pal(8, "Set2")),
  celltype = c(brewer.pal(8, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Dark2")),
  program = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2")),
  binary = c("#00AFBB", "#FC4E07"),
  soft8 = as.vector(hcl.colors(8, "Set2"))
)

.palette_rules <- list(
  Tissue       = "binary",
  Dataset      = "dataset",
  Cluster      = "set_base",
  State        = "set_base",
  CellType     = "celltype",
  `Cell Type`  = "celltype",
  Program      = "program"
)

pal_get <- function(n, name = "set_base") {
  base <- .pub_palettes[[name]]
  if (is.null(base)) base <- .pub_palettes$set_base
  if (n <= length(base)) return(base[seq_len(n)])
  rep <- rep_len(base, n)
  k <- 7L %% n
  if (k == 0) {
    k <- 1L
  }
  rep[c((k + 1):n, 1:k)]
}

map_discrete <- function(x, colname = NULL, fallback = "set_base") {
  x <- safe_factor(x)
  lev <- levels(x)
  pal_name <- if (!is.null(colname) && colname %in% names(.palette_rules)) .palette_rules[[colname]] else fallback
  cols <- pal_get(length(lev), pal_name)
  stats::setNames(cols, lev)
}

safe_factor <- function(x) {
  vals <- as.character(x)
  vals[is.na(vals) | !nzchar(vals)] <- "(missing)"
  factor(vals)
}

# ------- Scaling -------------------------------------------------------------
zscale_rows <- function(mat) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (ncol(mat) < 2) return(mat)
  cm <- rowMeans(mat, na.rm = TRUE)
  cs <- matrixStats::rowSds(mat, na.rm = TRUE)
  cs[cs == 0 | is.na(cs)] <- 1
  out <- sweep(sweep(mat, 1, cm, "-"), 1, cs, "/")
  out[!is.finite(out)] <- 0
  out
}

# ------- Top-gene labels (mean in group minus out) ---------------------------
compute_top_gene_labels <- function(mat, classes_tbl, group_col = "Program", per_cluster = 10) {
  stopifnot(ncol(mat) == nrow(classes_tbl))
  groups <- as.character(classes_tbl[[group_col]])
  lev <- unique(groups)
  keep <- integer(0)
  labs <- character(0)
  for (g in lev) {
    in_g <- groups == g
    mu_in <- rowMeans(mat[, in_g, drop = FALSE], na.rm = TRUE)
    mu_out <- rowMeans(mat[, !in_g, drop = FALSE], na.rm = TRUE)
    sc <- mu_in - mu_out
    ord <- order(sc, decreasing = TRUE)
    sel <- ord[seq_len(min(per_cluster, length(ord)))]
    keep <- c(keep, sel)
    labs <- c(labs, rownames(mat)[sel])
  }
  tibble(which = keep, label = labs) %>%
    dplyr::distinct(.data$label, .keep_all = TRUE) %>%
    dplyr::arrange(.data$which)
}

# ------- Single versatile heatmap wrapper -----------------------------------
# One wrapper for any numeric matrix. It draws heatmaps only (no analysis).
# Always rastered. You can provide **independent** column and row annotations.
# Palettes are routed by annotation column names, and can be overridden.
# Additional ComplexHeatmap::Heatmap args are forwarded via ...
#
# Arguments:
#   mat                 numeric matrix (rows x cols)
#   col_anno_df         data.frame with annotations aligned to **columns** of mat (nrow = ncol(mat))
#   col_anno_cols       character vector: which columns from col_anno_df to show on top
#   row_anno_df         data.frame with annotations aligned to **rows** of mat (nrow = nrow(mat))
#   row_anno_cols       character vector: which columns from row_anno_df to show at left
#   palette_overrides   named list colname -> color vector (applies to both row/col if names match)
#   ann_top_args        named list passed into HeatmapAnnotation(...)
#   ann_row_args        named list passed into rowAnnotation(...)
#   scale_rows          logical; if TRUE z-score rows before plotting
#   value_range         numeric vector of anchors for continuous colors (single knob)
#   value_colors        colors corresponding to value_range (defaults chosen automatically)
#   top_labels          NULL/TRUE/data.frame; if TRUE auto-compute anno_mark labels; if data.frame
#                       must contain (which,label)
#   group_col           column in col_anno_df used when auto-computing top_labels
#   show_row_names/column_names, cluster_rows/columns, title
#   ...                 forwarded to ComplexHeatmap::Heatmap

draw_heatmap <- function(
    mat,
    col_anno_df = NULL,
    col_anno_cols = NULL,
    row_anno_df = NULL,
    row_anno_cols = NULL,
    hm_name = "value",
    palette_overrides = NULL,
    ann_top_args = list(),
    ann_row_args = list(),
    scale_rows = TRUE,
    value_range = NULL,
    value_colors = NULL,
    top_labels = NULL,
    per_cluster = 10,
    group_col = "Program",
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    title = NULL,
    ...) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)

  # ---- Row scaling ----------------------------------------------------------
  mat_plot <- if (isTRUE(scale_rows)) zscale_rows(mat) else mat

  # ---- Colors for the heat tiles (single unified knob) ---------------------
  if (is.null(value_range) || is.null(value_colors)) {
    vals <- mat_plot[is.finite(mat_plot)]
    if (!length(vals)) {
      value_range <- c(0, 1)
      value_colors <- viridisLite::viridis(2)
    } else {
      rng <- quantile(vals, na.rm = TRUE, c(0.025, 0.975))
      # If the matrix is constant (or nearly so), ensure distinct breakpoints.
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || isTRUE(rng[1] == rng[2])) {
        v <- stats::median(vals, na.rm = TRUE)
        delta <- max(1e-6, abs(v) * 1e-6)
        value_range <- c(v - delta, v, v + delta)
        value_colors <- viridisLite::viridis(3)
      } else if (isTRUE(any(mat_plot < 0, na.rm = TRUE) && any(mat_plot > 0, na.rm = TRUE))) {
        # Diverging values: use the established blue-black-yellow palette.
        max_abs <- max(abs(rng))
        value_colors <- hm_value_colors
        value_range <- seq(-max_abs, max_abs, length.out = length(value_colors))
      } else {
        value_range <- c(rng[1], mean(rng), rng[2])
        value_colors <- viridisLite::viridis(3)
      }
    }
  }
  col_fun <- circlize::colorRamp2(value_range, value_colors)

  # ---- Column annotations (top) --------------------------------------------
  ann_top <- NULL
  if (!is.null(col_anno_df) && !is.null(col_anno_cols) && length(col_anno_cols) > 0) {
    stopifnot(ncol(mat) == nrow(col_anno_df))
    top_df <- col_anno_df[, intersect(colnames(col_anno_df), col_anno_cols), drop = FALSE]
    top_col_list <- lapply(colnames(top_df), function(nm) {
      if (!is.null(palette_overrides) && !is.null(palette_overrides[[nm]])) {
        palette_overrides[[nm]]
      } else {
        map_discrete(top_df[[nm]], colname = nm)
      }
    })
    names(top_col_list) <- colnames(top_df)
    ann_top <- do.call(ComplexHeatmap::HeatmapAnnotation, utils::modifyList(list(
      df = data.frame(lapply(top_df, safe_factor), check.names = FALSE),
      col = top_col_list,
      border = TRUE,
      show_annotation_name = TRUE,
      annotation_name_side = "right",
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(fontsize = 9)
    ), ann_top_args))
  }

  # ---- Row annotations (left) ----------------------------------------------
  ann_left <- NULL
  if (!is.null(row_anno_df) && !is.null(row_anno_cols) && length(row_anno_cols) > 0) {
    stopifnot(nrow(mat) == nrow(row_anno_df))
    left_df <- row_anno_df[, intersect(colnames(row_anno_df), row_anno_cols), drop = FALSE]
    left_col_list <- lapply(colnames(left_df), function(nm) {
      if (!is.null(palette_overrides) && !is.null(palette_overrides[[nm]])) {
        palette_overrides[[nm]]
      } else {
        map_discrete(left_df[[nm]], colname = nm)
      }
    })
    names(left_col_list) <- colnames(left_df)
    row_ann_args <- utils::modifyList(list(
      show_annotation_name = TRUE,
      annotation_name_side = "top",
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(fontsize = 9)
    ), ann_row_args)
    ann_left <- do.call(rowAnnotation, c(list(
      df  = data.frame(lapply(left_df, safe_factor), check.names = FALSE),
      col = left_col_list
    ), row_ann_args))
  }

  # ---- Optional row labels via anno_mark (depends on col_anno_df for grouping)
  if (!is.null(top_labels)) {
    if (is.numeric(top_labels)) {
      stopifnot(!is.null(col_anno_df))
      lab_df <- compute_top_gene_labels(
        mat_plot, col_anno_df,
        group_col = group_col, per_cluster = top_labels
      )
    } else if (is.data.frame(top_labels)) {
      stopifnot("gene" %in% names(top_labels), group_col %in% names(top_labels))

      n_per <- if (exists("per_cluster", inherits = TRUE) && is.numeric(per_cluster)) per_cluster else 10

      # 1) take top-N per group *as provided* (no recomputation; keep original order)
      trimmed <- top_labels %>%
        dplyr::group_by(.data[[group_col]]) %>%
        dplyr::slice_head(n = n_per) %>%
        dplyr::ungroup()

      # 2) map by gene name to row indices in mat_plot; drop non-matches
      idx <- match(as.character(trimmed$gene), rownames(mat_plot))
      keep <- !is.na(idx)

      if (!all(keep)) {
        dropped <- unique(as.character(trimmed$gene[!keep]))
        if (length(dropped)) {
          warning(
            sprintf(
              "Dropped %d gene(s) not found in matrix rownames. Example: %s",
              length(dropped),
              paste(utils::head(dropped, 5), collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }

      lab_df <- dplyr::tibble(
        which = idx[keep],
        label = if ("label" %in% names(trimmed)) as.character(trimmed$label[keep]) else rownames(mat_plot)[idx[keep]]
      )
      lab_df[[group_col]] <- as.character(trimmed[[group_col]][keep])
    } else lab_df <- NULL

    if (!is.null(lab_df) && nrow(lab_df) > 0) {
      lab_df <- lab_df[is.finite(lab_df$which) & !is.na(lab_df$which), , drop = FALSE]
      lab_df <- lab_df[order(lab_df$which, lab_df$label), , drop = FALSE]
      if (!nrow(lab_df)) {
        return(ComplexHeatmap::Heatmap(
          mat_plot,
          name = hm_name, top_annotation = ann_top,
          left_annotation = ann_left,
          col = col_fun, show_row_names = show_row_names, show_column_names = show_column_names,
          cluster_rows = cluster_rows, cluster_columns = cluster_columns,
          use_raster = TRUE,
          raster_quality = 3,
          raster_device = "png",
          raster_by_magick = FALSE,
          column_title = title,
          ...
        ))
      }
      cols <- if (!is.null(group_col) && group_col %in% names(lab_df)) {
        map_discrete(lab_df[[group_col]], colname = group_col, fallback = "set_base")
      } else {
        c(default = "#404040")
      }
      col_vec <- if (!is.null(group_col) && group_col %in% names(lab_df)) {
        unname(cols[as.character(lab_df[[group_col]])])
      } else {
        rep(unname(cols[[1L]]), nrow(lab_df))
      }

      # Keep anno_mark labels readable without letting a small label set explode
      # the font size and swamp the heatmap.
      n_lab <- nrow(lab_df)
      avail_pt <- 72 * 2.8
      fs_pt <- min(7, max(3, floor(avail_pt / max(1, n_lab))))
      max_chars <- max(nchar(as.character(lab_df$label)), na.rm = TRUE)
      mark_width_mm <- min(110, max(30, 2.4 * max_chars + 14))

      mark_ann <- rowAnnotation(
        link = anno_mark(
          at          = lab_df$which,
          labels      = lab_df$label,
          side        = "left",
          labels_gp   = grid::gpar(fontsize = fs_pt, lineheight = 1.1, col = "black"),
          link_gp     = grid::gpar(col = col_vec),
          link_width  = grid::unit(6, "mm"),
          padding     = grid::unit(0, "mm"),
          extend      = grid::unit(c(0, 0), "mm")
        ),
        width = grid::unit(mark_width_mm, "mm")
      )
      ann_left <- if (is.null(ann_left)) mark_ann else c(ann_left, mark_ann)
    }
  }


  # ---- Base Heatmap ---------------------------------------------------------
  hm <- ComplexHeatmap::Heatmap(
    mat_plot,
    name = hm_name, top_annotation = ann_top,
    left_annotation = ann_left,
    col = col_fun, show_row_names = show_row_names, show_column_names = show_column_names,
    cluster_rows = cluster_rows, cluster_columns = cluster_columns,
    use_raster = TRUE,
    raster_quality = 3,
    raster_device = "png",
    raster_by_magick = FALSE,
    column_title = title,
    ...
  )

  hm
}

# Reused by scientific scatterplots; comparisons retain their own endpoint and
# grouping axis in the calling script. Axes may use independent ranges.
correlation_label <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || sd(x) == 0 || sd(y) == 0) {
    return(sprintf("r = NA\nP = NA; n = %d", length(x)))
  }
  test <- cor.test(x, y, method = "pearson")
  sprintf("r = %.2f\nP = %s; n = %d", unname(test$estimate),
    format.pval(test$p.value, digits = 2, eps = 0.0001), length(x))
}

# Assign features by their within-versus-outside-class mean difference. The
# returned order is explicit and deterministic; heatmaps must not recluster it.
block_row_order <- function(matrix, classes) {
  classes <- as.character(classes)
  stopifnot(length(classes) == ncol(matrix), !anyNA(classes))
  groups <- sort(unique(classes))
  stopifnot(length(groups) >= 2L)
  delta <- vapply(groups, function(group) {
    rowMeans(matrix[, classes == group, drop = FALSE]) -
      rowMeans(matrix[, classes != group, drop = FALSE])
  }, numeric(nrow(matrix)))
  assignment <- max.col(delta, ties.method = 'first')
  strength <- delta[cbind(seq_len(nrow(matrix)), assignment)]
  list(order = order(assignment, -strength, rownames(matrix)),
    assignment = groups[assignment], delta = strength)
}

# Reused by signature learning and the single-drug bridge. Display changes do
# not alter training templates, route assignments, coefficients or predictions.
plot_training_signatures <- function(template, drug_labels, doses, label,
                                     included = seq_len(ncol(template$correlations))) {
  observed <- template$correlations[, included, drop = FALSE]
  route <- template$membership[included]
  row_blocks <- block_row_order(observed, route)
  within_score <- vapply(seq_len(ncol(observed)), function(j) {
    mean(observed[row_blocks$assignment == as.character(route[j]), j])
  }, numeric(1))
  column_order <- order(route, -within_score, colnames(observed))
  row_order <- row_blocks$order
  observed <- observed[row_order, column_order, drop = FALSE]
  annotation <- data.frame(
    Signature = paste('Signature', route[column_order]),
    Dose = paste(doses[included][column_order], 'uM'),
    row.names = colnames(observed)
  )
  stopifnot(identical(colnames(observed), rownames(annotation)))
  breaks <- seq(-1, 1, length.out = length(hm_value_colors))
  left <- draw_heatmap(observed, scale_rows = FALSE,
    value_range = breaks, value_colors = hm_value_colors, hm_name = 'Training Pearson r',
    col_anno_df = annotation, col_anno_cols = names(annotation),
    show_row_names = TRUE, row_names_side = 'left', row_names_gp = grid::gpar(fontsize = 10),
    show_column_names = TRUE, column_labels = drug_labels[included][column_order],
    column_names_rot = 45, column_names_gp = grid::gpar(fontsize = 11),
    width = grid::unit(105, 'mm'), height = grid::unit(145, 'mm'),
    title = 'Observed training-cell correlations')
  collapsed <- template$collapsed[row_order, , drop = FALSE]
  right <- draw_heatmap(collapsed, scale_rows = FALSE,
    value_range = breaks, value_colors = hm_value_colors, hm_name = 'Collapsed mean r',
    show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 11), width = grid::unit(24, 'mm'),
    height = grid::unit(145, 'mm'), title = 'Fixed signatures', show_heatmap_legend = FALSE)
  drawn <- ComplexHeatmap::draw(left + right, heatmap_legend_side = 'right',
    annotation_legend_side = 'right', merge_legends = TRUE,
    ht_gap = grid::unit(8, 'mm'), column_title = label,
    column_title_gp = grid::gpar(fontsize = 14), padding = grid::unit(c(18, 8, 14, 8), 'mm'))
  stopifnot(identical(as.integer(ComplexHeatmap::row_order(drawn, name = 'Training Pearson r')),
    seq_len(nrow(observed))))
  grid::grid.text('Features ordered by within-versus-outside-signature mean difference. No row z-scoring.\nPositive correlations indicate greater sensitivity (-log2 abundance).',
    x = 0.5, y = grid::unit(4, 'mm'), just = c('centre', 'bottom'),
    gp = grid::gpar(fontsize = 10))
  invisible(drawn)
}
