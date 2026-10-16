# TUMORCAST project instructions

Agreed with the user on 2026-10-16. These instructions consolidate the current
project rules. Newer explicit user instructions take priority. Historical memos
and code are evidence of previous analyses, not additional layers of live rules.

## Project and language

- Recreate the selected molecular bridge and ligand-rescue results from their
  original data. Preserve the accepted figures byte-for-byte in `reference/`.
- Use R for data preparation, modeling, plotting and automation. Prefer base R,
  dplyr, ordinary data frames and named matrices, plus established scientific
  packages when needed. Do not execute Python or use reticulate. An exception
  requires an essential operation that cannot reasonably be done in R, explained
  before use. Existing historical Python code is not an exception.
- Keep this a simple bioinformatics project. Do not introduce data registries,
  workflow frameworks, wrapper stacks, generated pipelines or custom object systems.

## Readable code

- Write linear, numbered scripts in `scripts/`. Each script should be readable
  from top to bottom in the order the analysis is performed.
- One instruction per line. Never separate statements with semicolons, including
  library loads. Use two-space indentation and consistent spaces around operators.
- Break long expressions and function calls across readable lines, usually within
  100-120 characters. Use descriptive names and explicit intermediate objects.
- Put a small configuration block at the top. Run from the project root and use
  relative paths. Keep machine-specific library setup outside analysis scripts.
- Comment every meaningful step: why it is needed, its inputs and outputs, units,
  normalization, exclusions and scientific choices. Comments describe current work.
- Keep reused functions in a few coherent files in `R/`. Explain their arguments
  and returned values. Keep single-use calculations in the owning script.
- Prefer existing maintained packages. Avoid hidden globals, silent fallbacks,
  multiple supported modes without a present need, and thin runner scripts.

## Files, iteration and caching

- `data/raw/` holds original source data grouped by dataset. `data/processed/`
  holds reusable prepared inputs. `results/` holds scientific outputs. `tmp/`
  holds disposable scratch work. `reference/` holds accepted originals and history.
- Explain data origins and execution order briefly in `README.md` and in scripts.
  Do not create a separate metadata or dependency-management system.
- Preserve historical scientific script identifiers where their analysis is being
  reorganized without changing its method. New scientific iterations continue
  after 145 with matching numbered folders directly under `results/`. Preparation
  stages can use their own clearly identified execution order. Plot variants stay
  in the existing analysis folder. Preserve prior scripts and results.
- Save final figures, their necessary source tables, reusable fitted objects and
  diagnostics needed to assess the science. Keep cheap joins and temporary subsets
  in memory. Save an intermediate only for inspection, reuse or substantial time saved.
- Prefer qs for substantial R objects and TSV/TSV.gz for tables. Base R RDS is
  appropriate for small caches. Keep one clear format for each scientific object.
- Cache long-running steps with the relevant input and code fingerprints,
  parameters, seed and package versions stored in the same file as the result.
  Reuse only matching caches. Print whether work is reused or why it is stale.
  Write completed caches to a temporary file, then rename, so interruptions leave
  no apparently valid partial result. Keep this logic local and straightforward.
- Give substantial steps concise progress messages. Keep one useful run log and
  session information per analysis. Check space before creating large outputs.
- Preserve original data and human-authored pilots. Delete obsolete intermediates
  only after their relationship to the retained results has been checked.

## Scientific preservation and diagnostics

- Preserve the selected model's algorithms, training data, response scales,
  program coordinates, dose choices, folds and parameters during reorganization.
  Do not replace a method with a different approximation as a formatting change.
- Keep baseline -> molecular displacement -> function explicit. Preserve the
  selected functional decoder's input restrictions and ligand-composition rule.
- Check identifiers, duplicates, required columns, joins, missingness, matrix
  orientation and sample alignment before modeling. Explain exclusions and fail
  clearly when required data are absent or inconsistent.
- Keep discovery, training and evaluation separate. Do not tune on evaluation
  outcomes. Previously inspected results remain retrospective after rerunning.
- Compare rebuilt inputs, program weights, scores and predictions with accepted
  results before retiring their original dependencies. A frozen-model replay is
  a useful check but does not establish raw-to-result reproduction.
- Lead functional reporting with held-out Pearson correlations and their exact
  endpoint and comparison axis. Distinguish across-drug comparisons within a cell
  from across-cell comparisons at a fixed drug and dose. Use assay-unit effects,
  sensitive/resistant counts and discrimination where appropriate. Do not use
  RMSE/MSE reductions as the scientific performance argument.

## Figures

- Use ggplot2 for ordinary plots and ComplexHeatmap for heatmaps almost
  exclusively. Use the shared heatmap functions in `R/heatmaps.R`; retain their
  established behavior from the original `src/R/lib_heatmap.R` when migrating.
  Do not implement heatmaps with ggplot tiles/raster or unrelated heatmap packages.
- Scatter plotting panels are square, with independently chosen axis limits,
  fitted trend lines and stat_cor() or equivalent annotations showing Pearson r,
  P and n for the plotted observations. Use `theme(aspect.ratio = 1)`, not
  `coord_equal()` or `coord_fixed()`. Label retrospective P values as descriptive.
  Explain scientifically necessary exceptions in the plotting code.
- Heatmaps normally have coherent diagonal blocks. With class labels, assign
  rows by maximum within-versus-outside-class mean difference, order by class
  then decreasing difference, and order columns into matching meaningful blocks.
  Derive stable groups when needed. State why any alternative ordering is needed.
  Ordering must describe the data, not manufacture support for a conclusion.
- Once explicit heatmap order is computed, disable automatic reclustering.
  Gene-gene correlation heatmaps use the same gene order on both axes. Stacked
  heatmaps share one identical column order, normally from the top layer. Check
  annotation identifiers after reordering and report missing samples explicitly.
- Signed or centered heatmaps use the established blue-black-yellow diverging
  palette. Use consistent scales for comparable panels. Explain normalization,
  centering, scaling, clipping, subsets and what zero/white/black means. Do not
  independently rescale predictions to conceal compressed variation.
- Set explicit heatmap body width and height. Rasterize the heatmap body with
  `raster_quality >= 3`. Keep text and annotations readable in the export.
- Show all gene-symbol labels when practical. For large heatmaps, use `anno_mark`
  for roughly 5-10 informative features per block rather than tiny exhaustive labels.
- Row-annotation names sit above their bars; top-annotation names sit on the right.
  Top-annotation legends belong on the right of the figure. Keep stacked layers
  proportionate and reserve space for every label, legend and annotation.
- Use large, explicit text sizes that remain readable at final display size,
  consistent typography and category colors, and compact content-sized canvases.
  Enlarge the figure or reduce panel density when necessary for readability.
- Inspect actual exported figures for clipping, overlap, text size, ordering,
  alignment, trend lines and statistics. Keep the necessary plotting data alongside
  figures. Preserve accepted reference figures unchanged even when updating style.

## Sources of these instructions

User instructions on 2026-10-09 through 2026-10-16; the server coding standard;
the repository's `complexheatmap-diagonal-blocks` skill and shared heatmap code;
and the `r-publication-plots` skill. Earlier conflicting preferences for compressed
code, extensive metadata machinery or changing earlier scientific versions in
place are superseded by the explicit user instructions consolidated above.
