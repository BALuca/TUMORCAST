#!/usr/bin/env Rscript
# Inspect agreement of the two matched native HGF-control replicate contrasts.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(digest)
})

output_dir <- "results/130_reconstruct_hgf_sources"
cache <- qread("data/processed/hgf_sources.qs")
source_hash <- digest("scripts/130_reconstruct_hgf_sources.R", file = TRUE, algo = "sha256")
if (cache$fingerprint$code_hash != source_hash) stop("HGF cache is stale; rerun its reconstruction first")
current_input_hashes <- vapply(
  names(cache$fingerprint$input_hashes), digest, character(1), file = TRUE, algo = "sha256"
)
if (!identical(current_input_hashes, cache$fingerprint$input_hashes)) {
  stop("HGF source data changed; rerun reconstruction before plotting")
}

# Show the two accepted preprocessing choices. Each point is one mapped gene.
points <- cache$result$native_gene_responses |>
  filter(
    (study == "GSE26732" & preprocessing == "native_RMA_log2") |
      (study == "GSE26735" & preprocessing == "normexp_quantile_log2")
  ) |>
  mutate(platform = if_else(study == "GSE26732", "Affymetrix: native RMA", "Illumina: normexp + quantile"))
stopifnot(all(is.finite(points$replicate1_logFC)), all(is.finite(points$replicate2_logFC)))

# Compute a stat_cor equivalent on exactly the observations in each panel.
platforms <- unique(points$platform)
statistics <- lapply(platforms, function(platform) {
  rows <- points[points$platform == platform, ]
  test <- cor.test(rows$replicate1_logFC, rows$replicate2_logFC, method = "pearson")
  p_label <- if (test$p.value < 2.2e-16) "P < 2.2e-16" else paste0("P = ", format.pval(test$p.value, digits = 2))
  data.frame(
    platform = platform,
    label = paste0("r = ", sprintf("%.2f", unname(test$estimate)), "\n", p_label, "\nn = ", nrow(rows))
  )
}) |>
  bind_rows()

plot <- ggplot(points, aes(replicate1_logFC, replicate2_logFC)) +
  geom_point(size = 0.35, alpha = 0.2, color = "#2568A0") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7, color = "#BF3434") +
  geom_text(
    data = statistics,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.15,
    size = 3.8
  ) +
  facet_wrap(vars(platform), scales = "free", nrow = 1) +
  labs(
    title = "Native HGF response: replicate agreement",
    x = "Replicate 1: HGF - control (log2)",
    y = "Replicate 2: HGF - control (log2)",
    caption = "Each point is a gene. Pearson tests are descriptive across correlated genes."
  ) +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "plain"),
    plot.caption = element_text(size = 10, hjust = 0),
    panel.spacing = grid::unit(7, "mm")
  )

ggsave(file.path(output_dir, "hgf_replicate_agreement.pdf"), plot, width = 7.8, height = 4.4)
ggsave("tmp/hgf_replicate_agreement.png", plot, width = 7.8, height = 4.4, dpi = 150)
message("Saved the native HGF replicate diagnostic with square panels and independent axes.")
