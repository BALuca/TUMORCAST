# Gene-weight and pathway audit for the frozen 40-program bank.
# Run from this directory: Rscript audit_programs.R
# No viability values or prediction performance are used in annotation.

required <- c("readr", "dplyr", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else normalizePath(".")
input <- file.path(root, "input")
out <- file.path(root, "audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_table <- function(name) read_tsv(file.path(input, name), show_col_types = FALSE, progress = FALSE)
weights <- read_table("gene_program_weights.tsv.gz")
programs <- sprintf("P%02d", 1:40)
stopifnot(!anyDuplicated(weights$gene_symbol), all(programs %in% names(weights)))
universe <- weights$gene_symbol
quality <- read_table("program_quality.tsv")

read_gmt <- function(name) {
  lines <- strsplit(readLines(file.path(input, name), warn = FALSE), "\t", fixed = TRUE)
  stopifnot(length(lines) > 0, all(lengths(lines) >= 3))
  setNames(lapply(lines, function(x) intersect(unique(x[-c(1, 2)]), universe)),
           vapply(lines, `[`, character(1), 1))
}
sets <- list(
  Hallmark = read_gmt("h.all.v2026.1.Hs.symbols.gmt"),
  Reactome = read_gmt("c2.cp.reactome.v2026.1.Hs.symbols.gmt")
)

# Perturbation-derived response footprints, with positive/negative target signs kept separate.
env <- new.env(parent = emptyenv())
load(file.path(input, "progeny_model_human_full.rda"), envir = env)
progeny <- as_tibble(env$model_human_full)
sets$PROGENy <- list()
for (top in c(100L, 500L)) {
  net <- progeny |> group_by(pathway) |> slice_min(p.value, n = top, with_ties = FALSE) |> ungroup()
  for (pathway in unique(net$pathway)) for (sign in c("UP", "DOWN")) {
    targets <- net |> filter(.data$pathway == !!pathway, if (sign == "UP") weight > 0 else weight < 0) |> pull(gene)
    sets$PROGENy[[paste("PROGENY", pathway, sign, paste0("TOP", top), sep = "_")]] <- intersect(unique(targets), universe)
  }
}
load(file.path(input, "dorothea_hs.rda"), envir = env)
regulons <- as_tibble(env$dorothea_hs) |> filter(confidence %in% c("A", "B")) |>
  mutate(regulon = paste("DOROTHEA", tf, if_else(mor > 0, "UP", "DOWN"), sep = "_"))
sets$DoRothEA <- lapply(split(regulons$target, regulons$regulon), function(x) intersect(unique(x), universe))

results <- list(); index <- 0L
for (collection in names(sets)) {
  gs <- sets[[collection]][lengths(sets[[collection]]) >= 5 & lengths(sets[[collection]]) <= 2000]
  sizes <- lengths(gs)
  for (program in programs) {
    loading <- setNames(weights[[program]], universe)
    members <- names(loading)[loading > 0]
    ordered <- names(sort(loading[loading > 0], decreasing = TRUE))
    top30 <- head(ordered, 30); top50 <- head(ordered, 50)
    overlap <- vapply(gs, function(x) sum(x %in% members), numeric(1))
    supported_weight <- vapply(gs, function(x) sum(loading[x]), numeric(1))
    maximum_weight <- vapply(gs, function(x) max(loading[x]), numeric(1))
    index <- index + 1L
    results[[index]] <- tibble(
      program = program, collection = collection, pathway = names(gs),
      program_genes = length(members), pathway_genes_in_universe = sizes,
      overlap = overlap, universe_genes = length(universe),
      p = phyper(overlap - 1, sizes, length(universe) - sizes, length(members), lower.tail = FALSE),
      fold_enrichment = overlap / (length(members) * sizes / length(universe)),
      weight_fraction = supported_weight / sum(loading[members]),
      top30_hits = vapply(gs, function(x) sum(top30 %in% x), numeric(1)),
      top50_hits = vapply(gs, function(x) sum(top50 %in% x), numeric(1)),
      top30_weight_fraction = vapply(gs, function(x) sum(loading[intersect(top30, x)]), numeric(1)) / sum(loading[top30]),
      top50_weight_fraction = vapply(gs, function(x) sum(loading[intersect(top50, x)]), numeric(1)) / sum(loading[top50]),
      strongest_supporting_gene_fraction = ifelse(supported_weight > 0, maximum_weight / pmax(supported_weight, .Machine$double.eps), 0),
      supporting_top30_genes = vapply(gs, function(x) paste(intersect(top30, x), collapse = ", "), character(1)),
      supporting_genes = vapply(gs, function(x) paste(intersect(ordered, x), collapse = ", "), character(1))
    )
  }
}
enrichment <- bind_rows(results) |>
  mutate(
    p_top30 = phyper(top30_hits - 1, pathway_genes_in_universe,
                    universe_genes - pathway_genes_in_universe, pmin(30, program_genes), lower.tail = FALSE),
    p_top50 = phyper(top50_hits - 1, pathway_genes_in_universe,
                    universe_genes - pathway_genes_in_universe, pmin(50, program_genes), lower.tail = FALSE),
    p_leave_one_supporting_gene_out = if_else(overlap > 0,
      phyper(overlap - 2, pathway_genes_in_universe, universe_genes - pathway_genes_in_universe,
             program_genes - 1, lower.tail = FALSE), 1)
  ) |>
  group_by(collection) |>
  mutate(q = p.adjust(p, method = "BH"), q_top30 = p.adjust(p_top30, method = "BH"),
         q_top50 = p.adjust(p_top50, method = "BH"),
         q_leave_one_supporting_gene_out = p.adjust(p_leave_one_supporting_gene_out, method = "BH")) |>
  ungroup() |>
  arrange(program, collection, q)
write_tsv(enrichment, file.path(out, "pathway_enrichment_all_programs.tsv"))
chosen <- c("P03", "P16", "P23", "P29", "P31", "P32", "P35", "P36")
write_tsv(enrichment |> filter(program %in% chosen), file.path(out, "pathway_enrichment_selected_programs.tsv"))
signal_pattern <- "MYC|MTOR|PI3K|AKT|ERBB|EGFR|RAF|MEK|ERK|MAPK|RAS_|NFKB|NF_KAPPA|FOXO|EIF|TRANSLATION|RIBOSOM|RESPIRAT|MITOCHONDR|GOLGI|VESICLE|CILI|PROGENY|DOROTHEA"
candidates <- enrichment |>
  filter(program %in% chosen, q < 0.05, overlap >= 3, grepl(signal_pattern, pathway)) |>
  arrange(program, desc(top30_weight_fraction), q)
write_tsv(candidates, file.path(out, "signaling_and_submodule_candidates.tsv"))
stopifnot(all(c("Hallmark", "Reactome", "PROGENy", "DoRothEA") %in% enrichment$collection))
stopifnot(nrow(enrichment |> filter(program == "P31", pathway == "PROGENY_EGFR_DOWN_TOP500", overlap == 27, q < 1e-8)) == 1L)
stopifnot(nrow(enrichment |> filter(program == "P32", pathway == "PROGENY_EGFR_DOWN_TOP500", overlap == 49, q < 1e-15)) == 1L)
for (program in chosen) {
  top <- weights |>
    select(gene_symbol, weight = all_of(program)) |>
    filter(weight > 0) |>
    arrange(desc(weight)) |>
    mutate(rank = row_number(), weight_fraction = weight / sum(weight))
  write_tsv(top, file.path(out, paste0(program, "_ranked_genes.tsv")))
}
write_tsv(quality |> filter(program %in% chosen), file.path(out, "program_quality_selected.tsv"))
for (program in chosen) {
  message("\n", program)
  print(enrichment |> filter(program == !!program, q < 0.05, overlap >= 3) |>
          arrange(desc(top30_weight_fraction), q) |>
          select(pathway, q, overlap, weight_fraction, top30_weight_fraction, supporting_top30_genes) |>
          head(6), width = 160)
}
writeLines(capture.output(sessionInfo()), file.path(out, "audit_session_info.txt"))
