#!/usr/bin/env Rscript
# Diagnostics of the saved43 hierarchy; no gene clustering or regressions refitted.
suppressPackageStartupMessages({library(data.table); library(qs); library(ggplot2); library(cluster); library(patchwork)})
source('/work/src/program_pilot/config.R')
source(file.path(code_dir, 'functional_programs.R'))
setDTthreads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
output_dir <- '/work/results/program_pilot/43_global_mean_drug_relationships_ward_programs'
data_dir <- file.path(output_dir, 'data')
figure_dir <- file.path(output_dir, 'figures')
model_dir <- file.path(output_dir, 'models')
local_minimum_rank <- 2L
rank_grid <- 2:80
bank_path <- file.path(model_dir, 'program_bank.qs')
bank <- qread(bank_path)
current_rank <- bank$diagnostics$selected_rank
aggregate <- qread(bank$aggregate_path)
correlation <- aggregate$clustering_correlation
correlation[aggregate$pair_drug_counts == 0] <- 0
correlation[correlation > 1] <- 1
correlation[correlation < -1] <- -1
distance <- as.dist(sqrt(2 * (1 - correlation)))
check(identical(bank$tree$labels, bank$genes), 'Saved tree genes differ')
cophenetic_distance <- cophenetic(bank$tree)
raw_r <- cor(as.numeric(distance), as.numeric(cophenetic_distance))
fwrite(data.table(genes = length(bank$genes), pairs = length(distance), raw_tree_cophenetic_r = raw_r,
  algorithm = 'One Ward.D2 hierarchy; no consensus resampling'), file.path(data_dir, 'cophenetic_diagnostic.tsv'), sep = '\t')
# Full-matrix coefficient; deterministic subset only for plotting individual pairs.
set.seed(seed)
idx <- sample.int(length(distance), min(20000L, length(distance)))
df <- data.frame(input_distance = as.numeric(distance)[idx], tree_distance = as.numeric(cophenetic_distance)[idx])
p <- ggplot(df, aes(input_distance, tree_distance)) + geom_point(alpha = .08, size = .5, color = '#27818A') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = '#B24C35') +
  theme_classic(base_size = 11) + theme(aspect.ratio = 1) +
  labs(x = 'Distance from averaged gene relationships', y = 'Cophenetic distance in Ward tree',
    title = 'Version43 | fidelity of the single Ward.D2 hierarchy',
    subtitle = sprintf('Pearson r = %.3f across all %s gene pairs', raw_r, format(length(distance), big.mark = ',')),
    caption = '20,000 pairs displayed; coefficient uses every pair.\nUnobserved pairs use neutral correlation0, as in the clustering input.')
save_plot(p, 'cophenetic', 6.5, 6.4)
message('COPHENETIC READY: ', file.path(figure_dir, 'cophenetic.pdf'))
# The useful rank-dependent diagnostic is silhouette, computed for exact cuts
# of the one saved tree. Every cut is checkpointed as soon as it finishes.
nmf_input_signature <- digest::digest(correlation, algo = 'sha256')
rm(aggregate, correlation, cophenetic_distance, df)
invisible(gc())
rank_batches <- split(rank_grid, rep(seq_len(4L), length.out = length(rank_grid)))
parts <- parallel::mclapply(rank_batches, function(grid) {
  fit_ward_cuts(distance, bank$tree, grid, '.', 'pooled', nmf_input_signature, model_dir)
}, mc.cores = 4L)
check(!any(vapply(parts, inherits, logical(1), 'try-error')), 'A Ward diagnostic worker failed')
fits <- unlist(parts, recursive = FALSE)
ranks <- rbindlist(lapply(fits, `[[`, 'diagnostics'))[order(rank)]
ranks <- select_silhouette(ranks, rank_override = current_rank, selection_rule = 'maximum_silhouette')
fwrite(ranks, file.path(data_dir, 'rank_diagnostics.tsv'), sep = '\t')
p <- ggplot(ranks, aes(rank, mean_silhouette)) + geom_line() +
  geom_point(aes(color = minimum_size_qc_pass)) +
  scale_color_manual(values = c(`FALSE` = '#F8766D', `TRUE` = '#00BFC4')) +
  geom_vline(xintercept = current_rank, lty = 2, color = '#B24C35') +
  geom_point(data = ranks[rank == current_rank], shape = 1, size = 4, color = '#B24C35') +
  theme_classic(base_size = 11) +
  labs(x = 'Number of global programs', y = 'Mean gene silhouette', color = 'Every cluster >=12 genes',
    title = 'Version43 | cuts of the averaged-drug Ward.D2 hierarchy',
    subtitle = sprintf('Ranks2–80; current bank uses%d programs; original hierarchy unchanged', current_rank))
save_plot(p, 'silhouette', 8, 4.5)
p <- ggplot(ranks, aes(rank, minimum_size)) + geom_line() + geom_point(size = 1) +
  geom_hline(yintercept = 12, lty = 3) + geom_vline(xintercept = current_rank, lty = 2, color = '#B24C35') +
  theme_classic(base_size = 11) + labs(x = 'Number of global programs', y = 'Smallest cluster: genes',
    title = 'Version43 | minimum program size')
save_plot(p, 'minimum_program_size', 6, 4)
bank$ranks <- ranks
bank$all_partitions <- setNames(lapply(fits, `[[`, 'partition'), vapply(fits, function(f) as.character(f$diagnostics$rank), character(1)))
bank$diagnostics[, raw_tree_cophenetic := raw_r]
check(identical(cutree(bank$tree, current_rank), bank$partition), 'Current program partition changed')
qsave(bank, bank_path, preset = 'high')
message('COMPLETE: cophenetic, silhouette and size diagnostics')
