# TUMORCAST

Reconstruct the selected molecular bridge and ligand-rescue analyses from their
source measurements, including gene eligibility, normalization, program
generation, molecular prediction, functional fitting and figures. The three
accepted figure bundles are preserved unchanged in `reference/`.

**Verified on 2026-10-16:** `run_all.R` completed all 17 stages. Every discovery
and modeling stage has been rebuilt from the collected sources, including rank
screening and program generation. All 8,842 original gene assignments were
recovered; the maximum program-weight difference was 1.9e-8.

All 118 accepted figure files are preserved with zero differing bytes, and all
have regenerated counterparts. Functional prediction differences were below
0.00075 log2 abundance; HGF query differences remained below 0.076 percentage
points. See `results/17_verify_reproduction/` for numerical comparisons and figure
coverage, and `results/run_all.log` for the successful complete execution.

## Run the project

Run from `/work/TUMORCAST`. On the current server:

```bash
cd /work/TUMORCAST
export R_LIBS_USER=/work/renv/library/R-4.3/x86_64-pc-linux-gnu
Rscript --vanilla run_all.R
```

`run_all.R` runs the following scripts in fresh R processes. They can also be
run individually in this order. Numbering starts at `01` and describes execution
order; historical analysis numbers remain in reference code and comments.

| Script | Scientific work |
|---|---|
| `01_prepare_expression.R` | Original CCLE expression, duplicate symbols, Ensembl biotypes and the expression-only gene eligibility filter |
| `02_prepare_tahoe.R` | Reconstruct treated/control abundance and sample-local normalization across all discovery profiles |
| `03_build_gene_relationships.R` | Profile centering, drug/study-balanced gene correlations and equal-source Fisher combination |
| `04_generate_programs.R` | Original signed oNMF implemented in R, rank screening, repeated starts, consensus and frozen weights |
| `05_score_molecular_inputs.R` | Complete Tahoe plate pooling, all 4,410 selected molecular profiles, CCLE program scores and receptors |
| `06_prepare_functional_outcomes.R` | Raw PRISM replicate curves, upfront dose selection, raw ALMANAC matching and selected-dose interpolation |
| `07_learn_response_signatures.R` | Recipient-excluded correlation maps, two response signatures, inner shrinkage selection and diagnostics |
| `08_fit_single_drug_bridge.R` | RNA predictors, shared functional calibration, supported drug-pair interactions and whole-recipient exclusions |
| `09_reconstruct_hgf_sources.R` | Native GEO contrasts, source selection, corrected HGF profiles and replicate diagnostics |
| `10_prepare_ligand_assays.R` | Reconstruct all NRG1/EGF assay observations from the original Prism XML |
| `11_transport_hgf_responses.R` | Recipient-specific HGF, NRG1 and EGF gene responses and program projections |
| `12_fit_ligand_rescue.R` | Nested recipient-excluded bounded rescue calibration, external controls, HGF queries and structural checks |
| `13_annotate_programs.R` | Gene weights, pathway enrichment, perturbation footprints and bank quality |
| `14_plot_model_diagnostics.R` | Held-recipient molecular predictions, source relationships and raw drug/HGF opposition |
| `15_plot_polished_panels.R` | Regenerate all 17 polished panels and their four original composite arrangements |
| `16_assemble_figures.R` | Regenerate the additional joined arrangement and molecular-bridge/rescue support figure entirely in R |
| `17_verify_reproduction.R` | Compare rebuilt predictions and figure coverage; verify zero differing bytes in all 118 preserved originals |

Expensive calculations have caches containing input/code fingerprints,
settings, seeds and package versions. Matching caches are reused; stale ones
are identified and rebuilt. Program restarts are checkpointed separately so an
interruption does not discard completed fits. Temporary writes are renamed only
after successful completion. Final verification writes its comparisons to
`results/17_verify_reproduction/` and stops on incomplete reproduction.

## Keep the layout simple

- `scripts/`: readable, sequential scientific stages.
- `R/`: reused data-reading, program-fitting, molecular-prediction and plotting functions.
- `data/raw/`: source measurements, released molecular arrays and fixed annotation resources.
- `data/processed/`: reusable prepared matrices and expensive calculation caches.
- `results/`: matching numbered folders with scientific intermediate results, models, diagnostics and figures.
- `reference/`: accepted bundles, unchanged figures, historical code and comparison targets.
- `tmp/`: disposable scratch work.

Coding and figure instructions are in [AGENTS.md](AGENTS.md). Use R, base data
frames/matrices and dplyr; retain maintained scientific packages where appropriate.
Existing data.table arithmetic is kept where needed to preserve the original
measurement aggregation. Use one instruction per line and two-space indentation.

New scientific iterations get distinct numbered scripts and result folders.
Keep alternate figure views with their existing analysis. Put reusable functions
in the existing small libraries; keep scientific decisions visible in scripts.
Do not add data registries, workflow frameworks or stacks of runner scripts.

The accepted originals in `reference/figures/` are byte-for-byte recoveries.
Freshly generated figures belong in the numbered results. Numerical reproduction
and preservation of original image bytes are checked separately: PDF metadata
and rendering-library versions can change a regenerated file's bytes.

The storage audit and completed cleanup are recorded in
[SPACE_AUDIT.md](SPACE_AUDIT.md). About 58 GiB was recovered; roughly 33 GiB
remained free after collecting the sources and running the reconstruction.

The code repository is [BALuca/TUMORCAST](https://github.com/BALuca/TUMORCAST).
Git contains code, documentation and environment configuration only. Data,
scientific tables, models, generated figures and private Codex state remain local
and are excluded. Commit and push only when explicitly instructed.

## Source measurements and fixed resources

All retained source measurements are inside this folder. Large files were moved
rather than duplicated. Temporary links at their old locations preserve the
historical analyses while the replacement is verified.

| Folder in `data/raw/` | Source retained |
|---|---|
| `depmap/` | Full CCLE expression release, model identities and mutation annotations |
| `tahoe/` | Released selected PDEX mean-response arrays, plate controls/audits, metadata and downloaded source Parquet ranges |
| `cytosig/` | Full released gene-response matrix and original profile metadata |
| `prism/` | Original primary and secondary screens, including secondary replicate wells and treatment metadata |
| `almanac/` | Original November 2017 growth measurements, compound names and field definitions |
| `hgf_geo/` | Original GEO matrices and platform annotations for native HGF contrasts |
| `ligand_assays/` | Original eLife Prism XML measurements |
| `wilson/` | Published figures/XML and the retained measurement transcription |
| `annotations/` | Ensembl v86 and org.Hs.eg.db 3.18 SQLite databases; fixed GMT, PROGENy and DoRothEA resources; original cohort/label choices |
| `lincs/` | Original landmark expression and instance metadata retained with the historical source collection |

The Tahoe source level used here is released pseudobulk expression means and
their controls, not individual-cell FASTQ reprocessing. The Parquet file is a
**sparse partial download**: its 40.6-GB logical size does not mean every row
group is present. Read only validated downloaded groups. The selected arrays
and source groups are complementary; extra topotecan/vinblastine plate profiles
must be included before pooling the downstream molecular panel.

The fixed functional cohort and retained ALMANAC control-arm table contain
identifiers, not copied outcomes. They preserve the historical sample admission
and exact-equality deduplication choices. Every measurement is recalculated from
the raw assay files. Program labels are the accepted human-readable annotations;
their enrichment evidence is regenerated in stage 13.

Wilson's published bar heights and SDs are retained as manual source
transcriptions alongside the original images and article XML. The raw table
contains the three published observations; older model rows from the original
mixed comparison table are archived in `reference/results/`. The pipeline does
not perform automatic image digitization.

Program discovery does not use functional labels or AU565. The historical
nine-context dose selection used already exposed outcomes, and the ligand pilot
was developed retrospectively. Rerunning does not turn those comparisons into
prospective validation. Report held-recipient Pearson correlations across cells
at a fixed drug and dose, with assay units and the precise comparison axis.

## Container and conversation history

The image definition pins R 4.3.2 and the recorded R library in `env/renv.lock`.
One existing development dependency, BiocFileCache 2.99.5, is pinned to its exact
upstream source commit. The build context excludes data and results.

The new Compose setup keeps `/work` on the host and mounts
`/tme/mvp/tme-sim/.codex_state` at `/root/.codex`. Codex local state belongs in
that persistent mount; an image alone does not preserve later session writes.
See the [official Codex state-location documentation](https://learn.chatgpt.com/docs/config-file/config-advanced#config-and-state-locations)
and [Docker bind-mount documentation](https://docs.docker.com/engine/storage/bind-mounts/).

The current history/configuration has been copied to that private host directory,
including consistent SQLite backups. It is excluded from Git. Before the first
container replacement, refresh the copy in the current container:

```bash
Rscript --vanilla tools/preserve_codex_history.R
```

Then, on the host:

```bash
cd /tme/mvp/tme-sim/TUMORCAST
docker compose build analysis
docker compose up -d analysis
```

VS Code can reopen this folder using `.devcontainer/devcontainer.json`. Docker
is unavailable inside the current container, so the new image/mount setup has
not been built or activated here.
