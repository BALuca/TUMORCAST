# Storage audit, 2026-10-16

The mounted filesystem was effectively full at the start. Abandoned Git pack
files accounted for 55.93 GiB. They have been removed after checking Git
connectivity; indexed packs, commits, references and loose objects were retained.
The final temporary pack was also checked against open process descriptors.

Another 1.62 GiB was recovered by removing 3,755 duplicate PNG previews from
seven superseded program-pilot branches. Every corresponding vector PDF was
checked by checksum before and after cleanup and remains unchanged. Scientific
data, fitted objects and scripts were retained. The deletion record is
`reference/cleanup_removed_previews.tsv`.

About 58 GiB was recovered. About 33 GiB is currently available while the
new reconstruction runs. Moving
source measurements into this project avoids a second copy of those datasets.
Old paths temporarily point to the collected files so historical analyses can
still locate their sources.

The largest existing collections, rounded to GiB, are:

| Collection | Approximate allocated space | Interpretation |
|---|---:|---|
| `/work/data` | 100 | Mostly source data from several projects; not disposable |
| `/work/results` | 143 | Historical scientific outputs and repeated analyses |
| `/work/results/parry_core/resource_build` | 34 | Several footprint/program build versions |
| `/work/results/pathway_signatures` | 30 | Many historical analysis versions |
| `/work/cache/linclayer` | 20 | Includes downloaded signature measurements; not all expendable caches |
| `/work/results/program_pilot` | 19 | Selected and abandoned program-pilot branches |
| `/work/TUMORCAST` | 28 | Collected sources, accepted references and the current reconstruction |

Rows overlap; this table is not a sum. Tahoe sparse Parquet files have much
larger logical sizes than their allocated disk usage. Preserve their sparse
layout when copying them.

Automatic packing is disabled locally in the old `/work` repository to prevent
another background pack attempt from filling the filesystem. Its history has
not been rewritten. The new repository excludes raw data, fitted results and
private conversation state from Git.

Further cleanup targets are reproducible full-matrix text exports and
superseded prepared-input caches from abandoned pilot branches. Keep their
scripts, fitted objects, diagnostic summaries and final figures. Unique downloaded
measurements and unrelated project results stay intact. Retire an old dependency
only after checking that the selected reconstruction no longer reads it.

For ongoing work, reuse the collected raw inputs, checkpoint expensive fits,
and keep cheap joins in memory. Save scientific intermediates once in a compact
format. Keep plot variants with their analysis rather than copying entire result
trees. Check available space before new large downloads or matrix construction.
