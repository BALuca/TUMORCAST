#!/usr/bin/env python3
"""Prepare independent inputs for the Stage 28 FGF2--fulvestrant bypass study.

Downloads only public processed data and selected processed Tahoe PDEX ranges.
The script never reads GSE153509 expression values while preparing discovery
inputs; it only records the sample design and file checksum.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT.parent
ART = ROOT / "artifacts/28_fgf2_fulvestrant_bypass"
OUT = ROOT / "results/28_fgf2_fulvestrant_bypass"
ART.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

URLS = {
    "gse_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153509/suppl/GSE153509_FGFR.Rescue.seq.counts.QC.csv.gz",
    "gse_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153509/soft/GSE153509_family.soft.gz",
    "shee_supplement": "https://pmc.ncbi.nlm.nih.gov/articles/instance/5839765/bin/JEM_20171818_sm.pdf",
}
FILES = {
    "gse_counts": ART / "GSE153509_FGFR.Rescue.seq.counts.QC.csv.gz",
    "gse_soft": ART / "GSE153509_family.soft.gz",
    "shee_supplement": ART / "JEM_20171818_sm.pdf",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size:
        return
    partial = path.with_suffix(path.suffix + ".partial")
    with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as out:
        shutil.copyfileobj(response, out, 1024 * 1024)
    partial.replace(path)


for key, url in URLS.items():
    download(url, FILES[key])


def soft_samples(path: Path) -> pd.DataFrame:
    rows = []
    current = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    rows.append(current)
                current = {"gsm": line.split(" = ", 1)[1], "characteristics": []}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                current["characteristics"].append(line.split(" = ", 1)[1])
        if current:
            rows.append(current)
    out = pd.DataFrame(rows)
    out["characteristics"] = out.characteristics.map(lambda x: " | ".join(x))
    return out


samples = soft_samples(FILES["gse_soft"])
with gzip.open(FILES["gse_counts"], "rt") as handle:
    count_header = next(csv.reader(handle))[1:]
header = pd.DataFrame({"matrix_column": count_header})
samples = samples.merge(
    header, left_on="title", right_on="matrix_column",
    how="outer", validate="one_to_one", indicator=True
)
if set(samples._merge) != {"both"}:
    raise RuntimeError("SOFT titles and processed count columns do not map one-to-one")
samples = samples.drop(columns="_merge")
samples["construct"] = samples.title.str.split("_", n=1).str[0]
samples["replicate"] = pd.to_numeric(samples.title.str.extract(r"\.\.S(\d+)$")[0])
samples["treatment_raw"] = (
    samples.title.str.replace(r"^[^_]+_", "", regex=True)
    .str.replace(r"\.\.S\d+$", "", regex=True)
)


def arm(raw: str) -> str:
    parts = raw.split(".")
    if "FIIN" in parts:
        i = parts.index("FIIN")
        if i + 1 < len(parts) and parts[i + 1] == "3":
            parts[i:i + 2] = ["FIIN-3"]
    canonical = {
        "DMSO": "Control", "Fulv": "Fulvestrant",
        "Palbo": "Palbociclib", "Tram": "Trametinib",
    }
    return " + ".join(canonical.get(x, x) for x in parts)


samples["arm"] = samples.treatment_raw.map(arm)
samples["parental"] = samples.construct.str.lower().eq("parental")
for token, label in [
    ("Fulv", "fulvestrant"), ("FGF2", "fgf2"),
    ("Palbo", "palbociclib"), ("FIIN", "fiin3"), ("Tram", "trametinib"),
]:
    samples[label] = samples.treatment_raw.str.contains(token, regex=False)
samples.to_csv(OUT / "28a_gse153509_complete_sample_design.tsv", sep="\t", index=False)
samples[samples.parental].to_csv(
    ART / "gse153509_parental_samples.tsv", sep="\t", index=False
)

# Retrieve processed target profiles for every PDEX cancer cell absent from the
# existing local 35-line cache. The PDEX panel contains BT-474 and HS-578T, but
# does not contain MCF-7, T47D, or BT-483; ER-specific interpretation therefore
# remains anchored to BT-474 plus prespecified ER-response biology.
spec = importlib.util.spec_from_file_location(
    "pdex_io", WORK / "scripts/tahoe_prism_full_chain/01_stream_tahoe_full_gene.py"
)
io = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(io)
DRUGS = [
    "Fulvestrant", "palbociclib", "Erdafitinib",
    "Pemigatinib", "Infigratinib", "Futibatinib",
]
DOSES = [0.05, 0.5, 5.0]
TARGETS = {repr([(drug, dose, "uM")]) for drug in DRUGS for dose in DOSES}
source_cache = WORK / "cache/tahoe_prism_full_chain/tahoe_pdex_selected"
existing_cells = pd.read_csv(source_cache / "cells.tsv", sep="\t")
cell_catalog = pd.read_parquet(
    WORK / "outputs/dataset_inventory/cache/tahoe/cell_line_metadata.parquet"
)[["cell_name", "Cell_ID_DepMap", "Cell_ID_Cellosaur", "Organ"]].drop_duplicates()
genes = (
    pd.read_csv(source_cache / "genes.tsv", sep="\t")
    .sort_values("gene_index").feature.tolist()
)
target_outputs = [
    ART / "additional_treated_mean.npy",
    ART / "additional_reference_mean.npy",
    ART / "additional_response_metadata.tsv",
]

if not all(p.exists() for p in target_outputs):
    sparse = ART / ".pdex_ranges.partial"
    fd, parquet = io.initialize_sparse_parquet(io.PDEX_URL, sparse)
    pdex_cells = set()
    for rg in range(parquet.metadata.num_row_groups):
        lo, hi = io.statistic_value(parquet, rg, "cell_line")
        if lo == hi:
            pdex_cells.add(lo)
    missing = pdex_cells - set(existing_cells.cell_line)
    extra = cell_catalog[cell_catalog.Cell_ID_Cellosaur.isin(missing)].copy()
    extra = extra[
        extra.Cell_ID_DepMap.notna()
        & ~extra.Cell_ID_DepMap.isin(["ACH-000311", "ACH-000628"])
    ]
    if extra.Cell_ID_Cellosaur.nunique() != len(missing - {"CVCL_C466"}):
        raise RuntimeError("Could not map every additional PDEX cancer cell")
    CELLS = {
        r.Cell_ID_Cellosaur: (r.Cell_ID_DepMap, r.cell_name)
        for r in extra.drop_duplicates("Cell_ID_Cellosaur").itertuples()
    }
    selected = []
    for rg in range(parquet.metadata.num_row_groups):
        lo, hi = io.statistic_value(parquet, rg, "cell_line")
        if lo != hi or lo not in CELLS:
            continue
        tlo, thi = io.statistic_value(parquet, rg, "target")
        if not any(tlo <= target <= thi for target in TARGETS):
            continue
        start, end = io.parquet_row_group_span(parquet, rg)
        selected.append((rg, lo, start, end))
    print(
        f"Stage 28: {len(selected)} selected Tahoe row groups; "
        f"{sum(e-s+1 for _, _, s, e in selected)/1e6:.1f} MB maximum.",
        flush=True,
    )

    def receive(item):
        rg, cell, start, end = item
        content, _ = io.fetch_small_range(io.PDEX_URL, f"bytes={start}-{end}")
        if len(content) != end - start + 1:
            raise RuntimeError(f"Incomplete Tahoe row group {rg}")
        os.pwrite(fd, content, start)
        return item

    profiles = []
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for index, (rg, cell, start, end) in enumerate(pool.map(receive, selected), 1):
                table = parquet.read_row_group(
                    rg,
                    columns=[
                        "feature", "target", "target_mean", "ref_mean",
                        "target_membership", "ref_membership", "plate", "cell_line",
                    ],
                )
                lengths, targets = io.run_blocks(table["target"])
                offset = 0
                for n, target in zip(lengths, targets):
                    block = table.slice(offset, n)
                    offset += n
                    if target not in TARGETS:
                        continue
                    if n != len(genes):
                        raise RuntimeError("Unexpected Tahoe gene dimension")
                    feature = block["feature"].to_pylist()
                    order = pd.Index(feature).get_indexer(genes)
                    if len(set(feature)) != len(genes) or np.min(order) < 0:
                        raise RuntimeError("Tahoe gene map is not one-to-one")
                    plates = set(block["plate"].to_pylist())
                    cells = set(block["cell_line"].to_pylist())
                    if len(plates) != 1 or cells != {cell}:
                        raise RuntimeError("Ambiguous Tahoe plate or cell identity")
                    drug, dose, unit = io.target_parts(target)
                    profiles.append({
                        "cell_line": cell,
                        "depmap_id": CELLS[cell][0],
                        "cell_name": CELLS[cell][1],
                        "drug": drug,
                        "dose": dose,
                        "dose_unit": unit,
                        "plate": next(iter(plates)),
                        "treated": np.asarray(block["target_mean"].to_numpy())[order],
                        "reference": np.asarray(block["ref_mean"].to_numpy())[order],
                        "target_membership_mean": float(np.mean(block["target_membership"].to_numpy())),
                        "reference_membership_mean": float(np.mean(block["ref_membership"].to_numpy())),
                    })
                if index % 10 == 0 or index == len(selected):
                    print(f"Processed {index}/{len(selected)} row groups.", flush=True)
    finally:
        os.close(fd)
        sparse.unlink(missing_ok=True)

    metadata, treated, reference = [], [], []
    for cell in CELLS:
        for drug in DRUGS:
            for dose in DOSES:
                pp = [
                    x for x in profiles
                    if x["cell_line"] == cell and x["drug"] == drug and x["dose"] == dose
                ]
                if not pp:
                    continue
                metadata.append({
                    "array_row": len(metadata),
                    "cell_line": cell,
                    "depmap_id": CELLS[cell][0],
                    "cell_name": CELLS[cell][1],
                    "drug": drug,
                    "dose": dose,
                    "dose_unit": "uM",
                    "n_plates": len(pp),
                    "plates": ",".join(str(x["plate"]) for x in pp),
                    "target_membership_mean": np.mean([x["target_membership_mean"] for x in pp]),
                    "reference_membership_mean": np.mean([x["reference_membership_mean"] for x in pp]),
                })
                treated.append(np.mean([x["treated"] for x in pp], axis=0))
                reference.append(np.mean([x["reference"] for x in pp], axis=0))
    if not metadata:
        raise RuntimeError("No additional Tahoe target profiles were recovered")
    np.save(target_outputs[0], np.asarray(treated, dtype="<f4"))
    np.save(target_outputs[1], np.asarray(reference, dtype="<f4"))
    pd.DataFrame(metadata).to_csv(target_outputs[2], sep="\t", index=False)
    pd.DataFrame({"gene": genes}).to_csv(
        ART / "additional_gene_order.tsv", sep="\t", index=False
    )
    pd.DataFrame([
        {"cell_line": k, "depmap_id": v[0], "cell_name": v[1]}
        for k, v in CELLS.items()
    ]).to_csv(ART / "additional_cells.tsv", sep="\t", index=False)
else:
    CELLS = {
        r.cell_line: (r.depmap_id, r.cell_name)
        for r in pd.read_csv(ART / "additional_cells.tsv", sep="\t").itertuples()
    }

# Primary PRISM is available for fulvestrant, palbociclib, and erdafitinib.
prism_root = WORK / "data/sciPlex_PRISM/PRISM"
treatments = pd.read_csv(
    prism_root / "primary-screen-replicate-collapsed-treatment-info.csv"
)
treatments = treatments[
    treatments.name.str.lower().isin(["fulvestrant", "palbociclib", "erdafitinib"])
].copy()
if len(treatments) != 3:
    raise RuntimeError("Expected one primary-PRISM column for each principal drug")
values_path = prism_root / "primary-screen-replicate-collapsed-logfold-change.csv"
header = pd.read_csv(values_path, nrows=0).columns
use = [0] + [header.get_loc(x) for x in treatments.column_name]
values = pd.read_csv(values_path, usecols=use).rename(
    columns={header[0]: "row_name"}
)
cell_info = pd.read_csv(prism_root / "primary-screen-cell-line-info.csv")
long = values.melt(
    id_vars="row_name", var_name="column_name", value_name="prism_survival"
)
long = long.merge(
    treatments[["column_name", "name", "dose", "screen_id"]],
    on="column_name", validate="many_to_one",
).merge(
    cell_info[["row_name", "depmap_id", "passed_str_profiling"]],
    on="row_name", validate="many_to_one",
)
long["prism_sensitivity"] = -long.prism_survival
long.to_csv(ART / "primary_prism_responses.tsv", sep="\t", index=False)

manifest = []
for key, path in FILES.items():
    manifest.append({
        "input": key, "path": str(path), "url": URLS[key],
        "bytes": path.stat().st_size, "sha256": sha256(path),
    })
pd.DataFrame(manifest).to_csv(OUT / "28a_download_manifest.tsv", sep="\t", index=False)

summary = {
    "parental_profiles": int(samples.parental.sum()),
    "parental_arms": samples.loc[samples.parental].groupby("arm").size().to_dict(),
    "additional_tahoe_cells": CELLS,
    "tahoe_drugs": DRUGS,
    "tahoe_doses_uM": DOSES,
    "primary_prism_drugs": sorted(treatments.name.str.lower()),
}
(ART / "preparation_manifest.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
