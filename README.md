# TUMORCAST

Reconstruction of the selected molecular bridge and ligand-rescue analyses from
their original data, with the accepted polished figures preserved unchanged.

Use readable, sequential R scripts with base R, dplyr and established scientific
packages. The complete project instructions are in [AGENTS.md](AGENTS.md).

## Folders

- `data/raw/`: original source data grouped by dataset.
- `data/processed/`: reusable prepared inputs and expensive calculation caches.
- `scripts/`: numbered preparation and scientific analysis scripts.
- `R/`: functions reused by the analysis scripts.
- `results/`: numbered scientific results and their necessary diagnostics.
- `reference/`: original result bundles, polished figures and historical code.
- `tmp/`: disposable scratch files.

Code and instructions belong in Git. Large data, results and reference bundles
remain in the project folder and are excluded from Git.

## Selected results and reconstruction order

The accepted scientific results are the original 40-program molecular bank,
recipient-excluded signatures (95), single-drug bridge (106), corrected HGF
sources (130/131), bounded ligand rescue (144), and the supplied polished panels.

The full reconstruction proceeds through original source data, normalization,
program discovery, molecular prediction, functional prediction, ligand rescue,
and figures. Historical scientific identifiers are retained when reorganizing
their method. New scientific methods start after 145. Raw-data preparation uses
its own numbered stages, listed below as they become runnable.

## Run the completed stages

Run from `/work/TUMORCAST`. On the current server, use the existing R library:

```bash
export R_LIBS_USER=/work/renv/library/R-4.3/x86_64-pc-linux-gnu
Rscript --vanilla scripts/00_verify_reference_files.R
Rscript --vanilla scripts/123_prepare_ligand_assays.R
Rscript --vanilla scripts/130_reconstruct_hgf_sources.R
Rscript --vanilla scripts/131_transport_hgf_responses.R
Rscript --vanilla scripts/130_plot_hgf_diagnostics.R
```

The library setting belongs to this server's environment; analysis scripts do
not contain machine-specific library paths. Each scientific stage saves its
useful outputs, numerical comparisons, run log and R session under its numbered
result folder. HGF calculations reuse matching caches and identify stale ones.

- `00`: verifies the three supplied archives and recovers all 118 original PDF,
  PNG and SVG figure files with zero differing bytes. The different B2 Figure 2
  arrangements remain separate.
- `123`: reconstructs the 12 NRG1 and four EGF observations across six recipients
  from the original Prism XML. Replicate pairing, vehicle interpretation and
  within-replicate log-abundance/log-dose interpolation match the original.
- `130`: reconstructs native HGF-control contrasts from GEO matrices and platform
  annotations, selects the original CytoSig source contexts, and rebuilds the
  corrected, precision-weighted HGF gene profiles. Original missingness is retained.
- `131`: projects those profiles onto the original full CCLE expression data for
  14 recipients, using the existing frozen 40-program bank. All 126 molecular
  profiles are checked against the accepted values.
- `130_plot`: displays native replicate agreement with square, independently
  scaled panels, trend lines, descriptive Pearson r/P/n and readable text.

## Source data now collected here

| Current folder | Original location on this system |
|---|---|
| `data/raw/cytosig/` | `/work/data/01_CytoSig_data/` |
| `data/raw/depmap/` | `/work/outputs/dataset_inventory/cache/depmap/` |
| `data/raw/tahoe/` | Selected PDEX data in `/work/cache/tahoe_prism_full_chain/`, source ranges in `/work/almanac_conditionals/diagnosis/`, metadata in `/work/outputs/dataset_inventory/cache/tahoe/` |
| `data/raw/prism/` | `/work/cache/tahoe_prism_full_chain/prism/` and `/work/data/sciPlex_PRISM/PRISM/` |
| `data/raw/almanac/` | `/work/almanac_conditionals/raw/` |
| `data/raw/lincs/` | `/work/v20_pilot_design_2026-10-09/source/` |
| `data/raw/hgf_geo/` | Original source files retained by analysis 130 |
| `data/raw/ligand_assays/` | Original eLife Prism XML retained by analysis 123 |
| `data/raw/wilson/` | Original source figures/XML and the existing measurement transcription |
| `data/raw/annotations/` | Frozen annotation resources supplied with the accepted ggplot2 bundle |

Large source files were moved, not duplicated. Their old locations temporarily
remain symbolic links so historical analyses can still run during reconstruction.
The new scripts read data inside this project. Original code is preserved under
`reference/code/`; scientific comparison targets are under `reference/results/`.

The selected PDEX package contains extracted mean responses and plate controls.
The sparse Parquet source retains only downloaded ranges; it is not a complete
download of the full 40.6 GB release. Preserve its row-group information. These
files must not be mistaken for an already reconstructed full raw-data workflow.

## Still to reconstruct

The bank in `data/processed/program_bank/` is explicitly frozen. Its original
discovery intermediates are preserved in `data/processed/program_bank_original/`;
the required Tahoe normalization records are retained separately. Their R-only
raw-data reconstruction is not yet implemented.

Next are Tahoe normalization and program discovery, baseline/program scoring,
PRISM/ALMANAC preparation, recipient-excluded signatures and the single-drug
bridge, bounded ligand-rescue fitting, and regeneration of the polished figures.
The current stages verify real portions of the chain, not a completed end-to-end
reproduction. The original numerical references never generate the rebuilt HGF
or ligand-assay measurements; they are used only for comparisons.

After the new workflow runs independently, remove temporary legacy links and
checked obsolete intermediates, then finish the reproducible Docker setup and
persistent Codex-history mount.

The preceding Git cleanup removed 149 abandoned temporary pack files (55.4 GiB).
Original data, analysis results and normal Git history were preserved.
