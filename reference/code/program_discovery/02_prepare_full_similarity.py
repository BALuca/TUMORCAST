#!/usr/bin/env python3
"""Build compound-balanced Tahoe similarity across all eligible shared genes.

No sensitivity ranking admits genes. Discovery retains the original full
perturbation compendium and profile-centering machinery.
"""

from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path("/work/tahoe_prism_v2")
OUT = Path("/work/almanac_conditionals/state_conditioned_bridge/program_first/results")
ART = Path("/work/almanac_conditionals/state_conditioned_bridge/program_first/artifacts")
CACHE = Path("/work/cache/tahoe_prism_full_chain")
TAHOE = CACHE / "tahoe_pdex_selected"


def weighted_rank_correlation(x: np.ndarray, groups: pd.Series) -> np.ndarray:
    """Spearman correlation with equal total weight for every perturbation."""
    ranks = np.empty(x.shape, dtype=np.float64)
    for start in range(0, x.shape[1], 128):
        end = min(start+128, x.shape[1])
        ranks[:, start:end] = rankdata(x[:, start:end], axis=0, method="average")
    counts = groups.value_counts()
    w = groups.map(lambda z: 1.0 / counts[z]).to_numpy(dtype=np.float64)
    w /= w.sum()
    center = np.sum(ranks * w[:, None], axis=0)
    z = ranks - center
    cov = (z * w[:, None]).T @ z
    scale = np.sqrt(np.maximum(np.diag(cov), 1e-30))
    ans = cov / np.outer(scale, scale)
    np.fill_diagonal(ans, 1.0)
    return np.clip(ans, -1.0, 1.0)


def remove_profile_global_displacement(x: np.ndarray) -> np.ndarray:
    """Remove the nonspecific signed amplitude shared by all admitted genes.

    A treatment that moves most admitted genes in the same direction
    creates a broad positive gene--gene background that is unrelated to
    program membership.  Centering each response profile makes similarity ask
    which genes move together *relative to the rest of that same response*.
    It uses no intervention labels and retains every profile.
    """
    return x - np.nanmean(x, axis=1, keepdims=True)


def write_square(path, genes, matrix):
    np.save(ART / "discovery_tahoe_correlation.npy", matrix)


# Use the complete universe frozen by expression, biotype and coverage filters.
genes = pd.read_csv(OUT / "selected_gene_universe.tsv", sep="\t")["gene"].tolist()
gene_table = pd.read_csv(TAHOE / "genes.tsv", sep="\t")
gene_lookup = dict(zip(gene_table["feature"], gene_table["gene_index"]))
gene_index = np.array([gene_lookup[g] for g in genes], dtype=int)

meta = pd.read_csv(
    ROOT / "artifacts/02_response_programs/observed_response_metadata.tsv.gz", sep="\t"
)
meta = meta.loc[meta.depmap_id != "ACH-000248"].copy()
assert not meta.astype(str).apply(lambda s: s.str.contains(r"AU[ _-]?565|ACH-000248", case=False, regex=True)).any().any()
meta = meta.sort_values("selected_response_row").reset_index(drop=True)
norm = pd.read_csv(
    ROOT / "artifacts/02_response_programs/normalization_records.tsv.gz", sep="\t"
)
need = [
    "observation_id",
    "treated_total",
    "reference_total",
    "treated_correction_factor",
    "contributing_plates",
]
meta = meta.merge(norm[need], on="observation_id", how="left", validate="one_to_one")
if meta[need[1:]].isna().any().any():
    raise RuntimeError("Tahoe normalization records did not align")

plate_meta = pd.read_csv(TAHOE / "pdex_plate_reference_metadata.tsv", sep="\t")
plate_lookup = {
    (r.depmap_id, str(r.plate)): int(r.array_row)
    for r in plate_meta.itertuples(index=False)
}
plate_all = np.load(TAHOE / "pdex_plate_reference.npy", mmap_mode="r")
plate_selected = np.asarray(plate_all[:, gene_index], dtype=np.float64)

final = ART / "full_tahoe_candidate_responses.npy"
order_path = OUT / "full_tahoe_gene_order.tsv"
can_reuse = final.exists() and order_path.exists()
if can_reuse:
    old_order = (
        pd.read_csv(order_path, sep="\t").sort_values("matrix_column")["gene"].tolist()
    )
    old_shape = np.load(final, mmap_mode="r").shape
    can_reuse = old_order == genes and old_shape == (len(meta), len(genes))

if not can_reuse:
    partial = ART / "full_tahoe_candidate_responses.npy.partial"
    responses = np.lib.format.open_memmap(
        partial, mode="w+", dtype=np.float32, shape=(len(meta), len(genes))
    )
    for shard_number, (shard_path, block) in enumerate(
        meta.groupby("shard_path", sort=False), 1
    ):
        shard = np.load(CACHE / shard_path, mmap_mode="r")
        row_index = block["shard_row"].to_numpy(dtype=int)
        raw = np.asarray(shard[np.ix_(row_index, gene_index)], dtype=np.float64)
        reference = np.empty_like(raw)
        for j, row in enumerate(block.itertuples(index=False)):
            plates = [p.strip() for p in str(row.contributing_plates).split(",")]
            ref_rows = [plate_lookup[(row.depmap_id, p)] for p in plates]
            reference[j] = plate_selected[ref_rows].mean(axis=0)
        treated = np.maximum(raw + reference, 0.0)
        treated_cp = treated / block["treated_total"].to_numpy()[:, None] * 10000.0
        treated_cp /= block["treated_correction_factor"].to_numpy()[:, None]
        reference_cp = (
            reference / block["reference_total"].to_numpy()[:, None] * 10000.0
        )
        corrected = np.log2(1.0 + treated_cp) - np.log2(1.0 + reference_cp)
        responses[block.index.to_numpy()] = corrected.astype(np.float32)
        print(f"Tahoe shard {shard_number}: {len(block):,} profiles", flush=True)
    responses.flush()
    partial.replace(final)
else:
    print(
        "Reusing complete Tahoe response matrix for unchanged gene universe", flush=True
    )

responses = np.load(final, mmap_mode="r")
tahoe_relative = remove_profile_global_displacement(
    np.asarray(responses, dtype=np.float64)
)
tahoe_cor = weighted_rank_correlation(tahoe_relative, meta["drug"])
write_square(OUT / "discovery_tahoe_correlation.tsv.gz", genes, tahoe_cor)
meta.to_csv(
    OUT / "full_tahoe_profile_metadata.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)
pd.DataFrame({"gene": genes, "matrix_column": np.arange(len(genes))}).to_csv(
    OUT / "full_tahoe_gene_order.tsv", sep="\t", index=False
)
print(
    f"Wrote {len(meta):,} Tahoe profiles, {meta.drug.nunique()} drugs, {len(genes)} genes"
)
