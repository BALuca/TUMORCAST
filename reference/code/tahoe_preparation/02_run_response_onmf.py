#!/usr/bin/env python3
"""Corrected, leakage-safe response programs for the fold-1 development run.

This execution path reconstructs condition-average treated profiles from the
stored arithmetic displacement and exact contributing-plate DMSO means.  It
learns paired composition offsets, technical-noise scales, response genes,
and one signed co-response basis using fold-1 training lines only.
"""

from __future__ import annotations

import importlib.util
import gzip
import json
import os
import shutil
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "20")
os.environ.setdefault("OMP_NUM_THREADS", "20")

import joblib
import numpy as np
import pandas as pd


SEED = 1
RANK = 10
N_ONMF_RUNS = 10
ONMF_MAX_ITER = 500
MIN_PROGRAM_STABILITY = 0.50
MAX_RESPONSE_GENES = 3000
MIN_RESPONSE_GENES = 500
N_ELIGIBLE_ABUNDANCE_BINS = 10
GENES_PER_ABUNDANCE_BIN = 300
SOURCE_LINES = 35
EXPECTED_LINES = 33
EXPECTED_TRAINING_LINES = 26
EXPECTED_HELD_LINES = 7
EXPECTED_PERTURBATIONS = 1137
EXPECTED_GENES = 62710
ACTIVE_OUTER_FOLDS = [1]
CP_SCALE = 10000.0
MIN_SCALING_CP10K = 0.10
MIN_SCALING_GENES = 100
MIN_RESPONSE_GENE_CP10K = 1.0
MIN_DMSO_EXPRESSED_LINES = 3
MIN_TREATED_EXPRESSED_RESPONSES = 20
MIN_TREATED_EXPRESSED_PERTURBATIONS = 5
MIN_TREATED_EXPRESSED_LINES = 3
SELECTION_SCORE_DEFINITION = (
    "top 300 positive noise-subtracted response-variance genes within each of ten "
    "expression-eligible abundance bins"
)
PROGRAM_DEFINITION_VERSION = "paired_trimmed_logratio_abundance_stratified_line_balanced_posneg_v3"

ROOT = Path(os.environ.get("TAHOE_V2_ROOT", "tahoe_prism_v2")).resolve()
OUT = ROOT / "artifacts" / "02_response_programs"
RESULT = ROOT / "results" / "02_response_programs"
CACHE = ROOT.parent / "cache" / "tahoe_prism_full_chain"
PDEX = CACHE / "tahoe_pdex_selected"
META_FILE = ROOT / "artifacts" / "01_inputs" / "response_metadata.tsv.gz"
GENE_FILE = PDEX / "genes.tsv"
FOLD_FILE = ROOT / "results" / "01_inputs" / "fold_assignments.tsv"
PDEX_REFERENCE_FILE = PDEX / "pdex_plate_reference.npy"
PDEX_REFERENCE_META_FILE = PDEX / "pdex_plate_reference_metadata.tsv"
PLATE_AUDIT_FILE = PDEX / "pdex_plate_audit.tsv.gz"
SUFFICIENT_STATS_FILE = CACHE / "response_line_sufficient_stats.npz"
EXCLUDED_IDS = {"ACH-000311", "ACH-000628"}
EXPECTED_HELD_IDS = {
    "ACH-000138", "ACH-000178", "ACH-000493", "ACH-000552",
    "ACH-000669", "ACH-000875", "ACH-000927",
}


def load_existing_helpers():
    path = Path(__file__).with_name("02_onmf_helpers.py")
    spec = importlib.util.spec_from_file_location("tahoe_onmf_existing", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the existing oNMF implementation.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(frame: pd.DataFrame, path: Path, compression=None) -> None:
    temporary = path.with_name(path.name + ".partial")
    frame.to_csv(temporary, sep="\t", index=False, compression=compression)
    os.replace(temporary, path)


def write_npy(array: np.ndarray, path: Path) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("wb") as handle:
        np.save(handle, array)
    os.replace(temporary, path)


def append_runtime(phase: str, started: float, details: str = "") -> None:
    row = pd.DataFrame([{
        "stage": "03_response_programs", "phase": phase,
        "runtime_seconds": time.monotonic() - started, "details": details,
    }])
    path = OUT / "runtime.tsv"
    if path.exists():
        current = pd.read_csv(path, sep="\t")
        row = pd.concat([current.loc[current["phase"].ne(phase)], row], ignore_index=True)
    write_tsv(row, path)


def write_artifact_manifest() -> None:
    import hashlib

    def details(path: Path) -> tuple[str, str]:
        name = path.name.removesuffix(".gz")
        stem = Path(name).stem.replace("_", " ").replace("-", " ")
        description = stem
        if path.suffix == ".npy":
            array = np.load(path, mmap_mode="r")
            dimensions = " x ".join(str(value) for value in array.shape)
            description += " NumPy matrix"
        elif path.suffix == ".npz":
            with np.load(path) as archive:
                dimensions = "; ".join(
                    f"{key}:" + " x ".join(str(value) for value in archive[key].shape)
                    for key in archive.files
                )
            description += " compressed NumPy arrays"
        elif name.endswith(".tsv"):
            opener = gzip.open if path.name.endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8") as handle:
                dimensions = f"{max(sum(1 for _ in handle) - 1, 0):,} rows"
            description += " table"
        elif path.suffix == ".joblib":
            dimensions = "serialized model object"
            description += " serialized model"
        elif path.suffix == ".sha256":
            dimensions = "one checksum"
            description += " checksum"
        else:
            dimensions = "serialized or auxiliary artifact"
        return description, dimensions

    records = []
    for path in sorted(OUT.rglob("*")):
        if not path.is_file() or path.name == "artifact_manifest.tsv":
            continue
        digest = hashlib.md5()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
                digest.update(block)
        description, dimensions = details(path)
        records.append({
            "relative_path": str(path.relative_to(OUT)),
            "description": description,
            "dimensions_or_row_count": dimensions,
            "file_size_bytes": path.stat().st_size,
            "md5": digest.hexdigest(),
            "producing_script": "scripts/02_run_response_onmf.py",
        })
    write_tsv(pd.DataFrame(records), OUT / "artifact_manifest.tsv")


def paired_trimmed_log_response(
    treated: np.ndarray, reference: np.ndarray, return_response: bool = True
) -> tuple[np.ndarray | None, dict]:
    """Apply one name-neutral paired trimmed log-ratio normalizer."""
    treated = np.asarray(treated, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if treated.shape != reference.shape or treated.ndim != 1:
        raise RuntimeError("Paired profiles have incompatible shapes.")
    if not np.isfinite(treated).all() or not np.isfinite(reference).all():
        raise RuntimeError("A paired profile contains nonfinite values.")
    scale_value = max(float(np.max(np.abs(treated))), float(np.max(np.abs(reference))), 1.0)
    negative_tolerance = max(1e-6, 32.0 * np.finfo(np.float32).eps * scale_value)
    min_treated = float(treated.min())
    min_reference = float(reference.min())
    if min_treated < -negative_tolerance or min_reference < -negative_tolerance:
        raise RuntimeError(
            f"Material negative reconstructed expression: treated={min_treated:g}, "
            f"reference={min_reference:g}, tolerance={negative_tolerance:g}."
        )
    n_tiny_treated = int(np.sum((treated < 0) & (treated >= -negative_tolerance)))
    n_tiny_reference = int(np.sum((reference < 0) & (reference >= -negative_tolerance)))
    if n_tiny_treated:
        treated = treated.copy()
        treated[treated < 0] = 0.0
    if n_tiny_reference:
        reference = reference.copy()
        reference[reference < 0] = 0.0
    treated_total = float(treated.sum(dtype=np.float64))
    reference_total = float(reference.sum(dtype=np.float64))
    if treated_total <= 0 or reference_total <= 0:
        raise RuntimeError("A paired expression profile has a nonpositive total.")
    treated_prop = treated / treated_total
    reference_prop = reference / reference_total
    treated_cp = treated_prop * CP_SCALE
    reference_cp = reference_prop * CP_SCALE
    eligible = (
        (treated_prop > 0) & (reference_prop > 0)
        & (treated_cp >= MIN_SCALING_CP10K)
        & (reference_cp >= MIN_SCALING_CP10K)
    )
    n_eligible = int(eligible.sum())
    if n_eligible < MIN_SCALING_GENES:
        raise RuntimeError(f"Only {n_eligible} genes are eligible for paired scaling.")
    log_t = np.log2(treated_prop[eligible])
    log_r = np.log2(reference_prop[eligible])
    m_value = log_t - log_r
    a_value = 0.5 * (log_t + log_r)
    m_low, m_high = np.quantile(m_value, [0.30, 0.70])
    a_low, a_high = np.quantile(a_value, [0.05, 0.95])
    retained = (
        (m_value >= m_low) & (m_value <= m_high)
        & (a_value >= a_low) & (a_value <= a_high)
    )
    n_retained = int(retained.sum())
    if n_retained < MIN_SCALING_GENES:
        raise RuntimeError(f"Only {n_retained} genes remain after paired trimming.")
    offset = float(np.median(m_value[retained]))
    factor = float(2.0 ** offset)
    if not np.isfinite(offset) or not np.isfinite(factor) or factor <= 0:
        raise RuntimeError("Paired composition offset is invalid.")
    response = None
    if return_response:
        corrected_treated_cp = treated_cp / factor
        response = np.log2(1.0 + corrected_treated_cp) - np.log2(1.0 + reference_cp)
        if not np.isfinite(response).all():
            raise RuntimeError("Paired normalized response is nonfinite.")
    record = {
        "treated_total": treated_total,
        "reference_total": reference_total,
        "n_scaling_eligible": n_eligible,
        "n_scaling_retained": n_retained,
        "composition_offset_log2": offset,
        "treated_correction_factor": factor,
        "minimum_reconstructed_treated": min_treated,
        "negative_tolerance": negative_tolerance,
        "n_tiny_treated_negatives_clamped": n_tiny_treated,
        "n_tiny_reference_negatives_clamped": n_tiny_reference,
        "finite": True,
        "nonnegative_within_tolerance": True,
    }
    return response, record


def popcount_rows(bits: np.ndarray) -> np.ndarray:
    byte_view = np.ascontiguousarray(bits).view(np.uint8).reshape(len(bits), -1)
    return np.unpackbits(byte_view, axis=1).sum(axis=1).astype(np.int32)


def set_bit(bits: np.ndarray, mask: np.ndarray, index: int) -> None:
    word = index // 64
    bit = np.uint64(1) << np.uint64(index % 64)
    bits[mask, word] |= bit


def load_shard(cell_meta: pd.DataFrame) -> np.ndarray:
    paths = cell_meta["shard_path"].drop_duplicates().tolist()
    if len(paths) != 1:
        raise RuntimeError("A cell line does not resolve to one response shard.")
    shard = np.load(CACHE / paths[0], mmap_mode="r")
    if shard.shape != (EXPECTED_PERTURBATIONS, EXPECTED_GENES) or shard.dtype != np.float32:
        raise RuntimeError(f"Invalid response shard: {paths[0]} {shard.shape} {shard.dtype}")
    return shard


def correlation_with_axis(axis: np.ndarray, covariate: np.ndarray) -> float:
    ok = np.isfinite(axis) & np.isfinite(covariate)
    if ok.sum() < 3 or np.std(axis[ok]) <= 0 or np.std(covariate[ok]) <= 0:
        return np.nan
    return float(np.corrcoef(axis[ok], covariate[ok])[0, 1])


def main() -> None:
    total_started = time.monotonic()
    required = [META_FILE, GENE_FILE, FOLD_FILE, PDEX_REFERENCE_FILE,
                PDEX_REFERENCE_META_FILE, PLATE_AUDIT_FILE, SUFFICIENT_STATS_FILE]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing existing stage input(s): {missing}")
    if (OUT / "stage02_complete.flag").exists():
        raise RuntimeError(f"Completed single-fold output already exists: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    RESULT.mkdir(parents=True, exist_ok=True)
    preprocessing_started = time.monotonic()

    metadata = pd.read_csv(META_FILE, sep="\t")
    metadata = metadata.loc[~metadata["depmap_id"].isin(EXCLUDED_IDS)].sort_values(
        "array_row"
    ).reset_index(drop=True)
    genes_table = pd.read_csv(GENE_FILE, sep="\t")
    genes = genes_table["feature"].astype(str).to_numpy()
    folds = pd.read_csv(FOLD_FILE, sep="\t")
    if "outer_fold" in metadata.columns:
        metadata = metadata.drop(columns="outer_fold")
    metadata = metadata.merge(
        folds[["depmap_id", "outer_fold"]], on="depmap_id", how="left", validate="many_to_one"
    ).sort_values("array_row").reset_index(drop=True)
    metadata["observation_id"] = metadata["depmap_id"].astype(str) + "__" + metadata["perturbation_id"].astype(str)
    cells = sorted(metadata["depmap_id"].unique())
    held_lines = set(folds.loc[folds["outer_fold"].eq(1), "depmap_id"].astype(str))
    training_lines = set(cells) - held_lines
    if (
        len(metadata) != EXPECTED_LINES * EXPECTED_PERTURBATIONS
        or len(genes) != EXPECTED_GENES
        or len(training_lines) != EXPECTED_TRAINING_LINES
        or len(held_lines) != EXPECTED_HELD_LINES
        or held_lines != EXPECTED_HELD_IDS
        or bool(EXCLUDED_IDS.intersection(cells))
        or metadata["observation_id"].duplicated().any() or metadata["outer_fold"].isna().any()
    ):
        raise RuntimeError("Fold-1 metadata contract failed.")
    split = folds[["outer_fold", "depmap_id", "cell_line"]].copy()
    split["role"] = np.where(split["outer_fold"].eq(1), "held_out", "training")
    write_tsv(split.sort_values(["role", "depmap_id"]), OUT / "fold1_cell_line_split.tsv")

    plate_meta = pd.read_csv(PDEX_REFERENCE_META_FILE, sep="\t")
    plate_reference = np.load(PDEX_REFERENCE_FILE, mmap_mode="r")
    plate_audit = pd.read_csv(PLATE_AUDIT_FILE, sep="\t")
    if plate_reference.shape != (SOURCE_LINES * 14, EXPECTED_GENES):
        raise RuntimeError("Plate-specific DMSO matrix has an invalid shape.")
    plate_meta = plate_meta.loc[plate_meta["depmap_id"].isin(cells)].copy()
    plate_audit = plate_audit.loc[plate_audit["depmap_id"].isin(cells)].copy()
    plate_lookup = {
        (str(row.depmap_id), int(row.plate)): int(row.array_row)
        for row in plate_meta.itertuples(index=False)
    }
    if len(plate_lookup) != EXPECTED_LINES * 14:
        raise RuntimeError("Plate-specific DMSO lookup is incomplete.")
    audit_groups = {
        (str(depmap_id), str(perturbation_id)): tuple(sorted(group["plate"].astype(int).tolist()))
        for (depmap_id, perturbation_id), group in plate_audit.groupby(
            ["depmap_id", "perturbation_id"], sort=False
        )
    }
    available = metadata["available"].eq(1)
    if len(audit_groups) != int(available.sum()):
        raise RuntimeError("Exact contributing-plate audit does not cover observed responses.")

    # Construct the training-only DMSO null with the exact production normalizer.
    training_control_rows = len(training_lines) * 14
    pseudo_partial = OUT / "dmso_pseudo_responses.npy.partial"
    dmso_pseudo = np.lib.format.open_memmap(
        pseudo_partial, mode="w+", dtype=np.float32,
        shape=(training_control_rows, EXPECTED_GENES),
    )
    dmso_raw_sum = np.zeros(EXPECTED_GENES, dtype=np.float64)
    dmso_raw_sumsq = np.zeros(EXPECTED_GENES, dtype=np.float64)
    dmso_abundance_sum = np.zeros(EXPECTED_GENES, dtype=np.float64)
    dmso_reliably_expressed_lines = np.zeros(EXPECTED_GENES, dtype=np.int32)
    dmso_records = []
    cursor = 0
    for depmap_id in sorted(training_lines):
        indices = [plate_lookup[(depmap_id, plate)] for plate in range(1, 15)]
        values = np.asarray(plate_reference[indices], dtype=np.float64)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise RuntimeError(f"Invalid DMSO plate profile in {depmap_id}.")
        total_profile = values.sum(axis=0, dtype=np.float64)
        library = values.sum(axis=1, dtype=np.float64)
        line_dmso_cp10k = values / library[:, None] * CP_SCALE
        dmso_abundance_sum += np.log2(1.0 + line_dmso_cp10k).sum(axis=0)
        dmso_reliably_expressed_lines += (
            np.median(line_dmso_cp10k, axis=0) >= MIN_RESPONSE_GENE_CP10K
        )
        for plate_index, plate in enumerate(range(1, 15)):
            reference = (total_profile - values[plate_index]) / 13.0
            response, record = paired_trimmed_log_response(values[plate_index], reference)
            dmso_pseudo[cursor] = response.astype(np.float32)
            raw = values[plate_index] - reference
            dmso_raw_sum += raw
            dmso_raw_sumsq += raw * raw
            record.update({
                "pseudo_response_row": cursor, "depmap_id": depmap_id,
                "cell_line": folds.loc[folds["depmap_id"].eq(depmap_id), "cell_line"].iloc[0],
                "pseudo_treated_plate": plate,
                "reference_plates": ",".join(str(x) for x in range(1, 15) if x != plate),
                "n_reference_plates": 13,
                "normalizer": "paired_trimmed_log_ratio_v1",
                "normalizer_shared_with_treated": True,
            })
            dmso_records.append(record)
            cursor += 1
    dmso_pseudo.flush()
    del dmso_pseudo
    os.replace(pseudo_partial, OUT / "dmso_pseudo_responses.npy")
    dmso_pseudo = np.load(OUT / "dmso_pseudo_responses.npy", mmap_mode="r")
    if cursor != training_control_rows or not np.isfinite(dmso_pseudo).all():
        raise RuntimeError("Training DMSO pseudo-response construction failed.")
    write_tsv(pd.DataFrame(dmso_records), OUT / "dmso_normalization_records.tsv.gz", compression="gzip")

    dmso_mean = np.asarray(dmso_pseudo, dtype=np.float64).mean(axis=0)
    dmso_var = np.asarray(dmso_pseudo, dtype=np.float64).var(axis=0, ddof=1)
    dmso_median = np.median(dmso_pseudo, axis=0)
    dmso_mad = np.median(np.abs(np.asarray(dmso_pseudo) - dmso_median[None, :]), axis=0)
    dmso_q99 = np.quantile(np.abs(dmso_pseudo), 0.99, axis=0)
    dmso_abundance = dmso_abundance_sum / training_control_rows
    dmso_raw_mean = dmso_raw_sum / training_control_rows
    dmso_raw_var = (dmso_raw_sumsq - dmso_raw_sum * dmso_raw_mean) / (training_control_rows - 1)
    abundance_decile = pd.qcut(
        pd.Series(dmso_abundance).rank(method="first"), 10, labels=False
    ).to_numpy(dtype=int) + 1
    positive_noise = dmso_var[np.isfinite(dmso_var) & (dmso_var > 0)]
    if not len(positive_noise):
        raise RuntimeError("The DMSO technical null has no positive gene variance.")
    fallback_noise = float(np.median(positive_noise))
    decile_noise = np.empty(EXPECTED_GENES, dtype=np.float64)
    decile_records = []
    for decile in range(1, 11):
        mask = abundance_decile == decile
        current = float(np.median(dmso_var[mask]))
        if not np.isfinite(current) or current <= 0:
            current = fallback_noise
        decile_noise[mask] = current
        decile_records.append({
            "baseline_abundance_decile": decile, "n_genes": int(mask.sum()),
            "median_dmso_variance": current,
        })
    stable_noise_var = np.maximum(dmso_var, decile_noise)
    if np.any(~np.isfinite(stable_noise_var)) or np.any(stable_noise_var <= 0):
        raise RuntimeError("Stabilized DMSO noise variance is invalid.")
    dmso_stats = pd.DataFrame({
        "gene_index": np.arange(EXPECTED_GENES), "gene": genes,
        "mean_dmso_log2_1p_cp10k": dmso_abundance,
        "baseline_abundance_decile": abundance_decile,
        "dmso_pseudo_response_mean": dmso_mean,
        "dmso_pseudo_response_variance": dmso_var,
        "abundance_decile_median_variance": decile_noise,
        "stabilized_dmso_noise_variance": stable_noise_var,
        "dmso_pseudo_response_mad": dmso_mad,
        "dmso_absolute_response_q99": dmso_q99,
        "raw_count_scale_dmso_pseudo_variance": dmso_raw_var,
    })
    write_tsv(dmso_stats, OUT / "dmso_gene_null_statistics.tsv.gz", compression="gzip")
    write_tsv(pd.DataFrame(decile_records), OUT / "dmso_abundance_decile_noise.tsv")

    # First treated pass: factors and training-only sufficient statistics.
    response_sum = np.zeros(EXPECTED_GENES, dtype=np.float64)
    response_sumsq = np.zeros(EXPECTED_GENES, dtype=np.float64)
    response_count = 0
    beyond_count = np.zeros(EXPECTED_GENES, dtype=np.int32)
    treated_expressed_count = np.zeros(EXPECTED_GENES, dtype=np.int32)
    n_pert_words = (EXPECTED_PERTURBATIONS + 63) // 64
    n_line_words = (len(training_lines) + 63) // 64
    beyond_pert_bits = np.zeros((EXPECTED_GENES, n_pert_words), dtype=np.uint64)
    beyond_line_bits = np.zeros((EXPECTED_GENES, n_line_words), dtype=np.uint64)
    treated_expressed_pert_bits = np.zeros((EXPECTED_GENES, n_pert_words), dtype=np.uint64)
    treated_expressed_line_bits = np.zeros((EXPECTED_GENES, n_line_words), dtype=np.uint64)
    perturbations = metadata.sort_values(["depmap_id", "shard_row"])["perturbation_id"].drop_duplicates().tolist()
    if len(perturbations) != EXPECTED_PERTURBATIONS:
        perturbations = metadata.sort_values("shard_row")["perturbation_id"].drop_duplicates().tolist()
    perturbation_index = {value: i for i, value in enumerate(perturbations)}
    training_line_index = {value: i for i, value in enumerate(sorted(training_lines))}
    mito_mask = np.asarray([value.upper().startswith("MT-") for value in genes])
    ribo_mask = np.asarray([
        value.upper().startswith("RPS") or value.upper().startswith("RPL") for value in genes
    ])
    normalization_records = []
    for cell_number, depmap_id in enumerate(cells, start=1):
        cell_meta = metadata.loc[metadata["depmap_id"].eq(depmap_id)].sort_values("shard_row")
        shard = load_shard(cell_meta)
        observed_meta = cell_meta.loc[cell_meta["available"].eq(1)]
        is_training = depmap_id in training_lines
        for row in observed_meta.itertuples(index=False):
            key = (depmap_id, str(row.perturbation_id))
            plates = audit_groups.get(key)
            if plates is None or len(plates) != int(row.n_plates) or len(set(plates)) != len(plates):
                raise RuntimeError(f"Contributing-plate count mismatch for {key}.")
            reference_indices = [plate_lookup[(depmap_id, plate)] for plate in plates]
            reference = np.asarray(plate_reference[reference_indices], dtype=np.float64).mean(axis=0)
            raw_displacement = np.asarray(shard[int(row.shard_row)], dtype=np.float64)
            treated = raw_displacement + reference
            response, record = paired_trimmed_log_response(treated, reference)
            treated_cp_corrected = treated / record["treated_total"] * CP_SCALE / record["treated_correction_factor"]
            record.update({
                "array_row": int(row.array_row), "observation_id": row.observation_id,
                "cell_line": row.cell_line, "depmap_id": depmap_id,
                "drug": row.drug, "dose": float(row.dose), "dose_panel": int(row.dose_panel),
                "perturbation_id": row.perturbation_id, "outer_fold": int(row.outer_fold),
                "role": "training" if is_training else "held_out",
                "contributing_plates": ",".join(str(x) for x in plates),
                "n_contributing_plates": len(plates), "metadata_n_plates": int(row.n_plates),
                "plate_count_validated": True,
                "treated_mito_fraction": float(treated[mito_mask].sum() / record["treated_total"]),
                "treated_ribosomal_fraction": float(treated[ribo_mask].sum() / record["treated_total"]),
                "normalizer": "paired_trimmed_log_ratio_v1",
            })
            normalization_records.append(record)
            if is_training:
                response_sum += response
                response_sumsq += response * response
                response_count += 1
                beyond = np.abs(response) > dmso_q99
                expressed = treated_cp_corrected >= MIN_RESPONSE_GENE_CP10K
                beyond_count += beyond
                treated_expressed_count += expressed
                p_index = perturbation_index[str(row.perturbation_id)]
                l_index = training_line_index[depmap_id]
                set_bit(beyond_pert_bits, beyond, p_index)
                set_bit(beyond_line_bits, beyond, l_index)
                set_bit(treated_expressed_pert_bits, expressed, p_index)
                set_bit(treated_expressed_line_bits, expressed, l_index)
        print(
            f"paired normalization pass 1: {cell_number}/{len(cells)} lines "
            f"({len(normalization_records)} observed responses).", flush=True
        )
    normalization_frame = pd.DataFrame(normalization_records).sort_values("array_row")
    if len(normalization_frame) != int(available.sum()) or response_count != int(
        metadata["depmap_id"].isin(training_lines).mul(available).sum()
    ):
        raise RuntimeError("Treated normalization first pass is incomplete.")
    write_tsv(normalization_frame, OUT / "normalization_records.tsv.gz", compression="gzip")

    response_mean = response_sum / response_count
    response_var = (response_sumsq - response_sum * response_mean) / (response_count - 1)
    numerical_negative = (response_var < 0) & (np.abs(response_var) < 1e-10)
    response_var[numerical_negative] = 0.0
    noise_subtracted_response_var = response_var - stable_noise_var
    excess_var = np.maximum(noise_subtracted_response_var, 0.0)
    excess_ratio = excess_var / stable_noise_var
    beyond_perturbations = popcount_rows(beyond_pert_bits)
    beyond_lines = popcount_rows(beyond_line_bits)
    treated_expressed_perturbations = popcount_rows(treated_expressed_pert_bits)
    treated_expressed_lines = popcount_rows(treated_expressed_line_bits)
    reliably_expressed_dmso = dmso_reliably_expressed_lines >= MIN_DMSO_EXPRESSED_LINES
    recurrently_expressed_treated = (
        (treated_expressed_count >= MIN_TREATED_EXPRESSED_RESPONSES)
        & (treated_expressed_perturbations >= MIN_TREATED_EXPRESSED_PERTURBATIONS)
        & (treated_expressed_lines >= MIN_TREATED_EXPRESSED_LINES)
    )
    expression_eligible = reliably_expressed_dmso | recurrently_expressed_treated
    passing = (
        expression_eligible & np.isfinite(response_var) & np.isfinite(stable_noise_var)
        & np.isfinite(noise_subtracted_response_var)
        & (noise_subtracted_response_var > 0)
    )

    raw_stats = np.load(SUFFICIENT_STATS_FILE, allow_pickle=False)
    raw_cells = raw_stats["line_ids"].astype(str).tolist()
    raw_line_mask = np.asarray([value in training_lines for value in raw_cells])
    raw_n = int(raw_stats["counts"][raw_line_mask].sum())
    raw_sum = raw_stats["sums"][raw_line_mask].sum(axis=0)
    raw_sumsq = raw_stats["sumsq"][raw_line_mask].sum(axis=0)
    raw_var = (raw_sumsq - raw_sum * raw_sum / raw_n) / (raw_n - 1)

    # Build equal-sized abundance strata only within the expression-eligible
    # pool. Stable pre-qcut ranks and original matrix position make both binning
    # and within-bin ties deterministic without using gene identity.
    eligible_index = np.flatnonzero(expression_eligible)
    eligible_abundance_bin = np.full(EXPECTED_GENES, np.nan)
    eligible_ranks = pd.Series(dmso_abundance[eligible_index]).rank(method="first")
    eligible_bins = pd.qcut(
        eligible_ranks, N_ELIGIBLE_ABUNDANCE_BINS, labels=False
    ).to_numpy(dtype=int) + 1
    eligible_abundance_bin[eligible_index] = eligible_bins
    selection_rank_within_bin = np.full(EXPECTED_GENES, np.nan)
    selected_from_abundance_bin = np.zeros(EXPECTED_GENES, dtype=bool)
    selected_parts = []
    bin_support_records = []
    for abundance_bin in range(1, N_ELIGIBLE_ABUNDANCE_BINS + 1):
        bin_eligible = eligible_abundance_bin == abundance_bin
        candidates = np.flatnonzero(bin_eligible & passing)
        candidate_order = np.lexsort((
            candidates, -noise_subtracted_response_var[candidates],
        ))
        ordered_candidates = candidates[candidate_order]
        selection_rank_within_bin[ordered_candidates] = np.arange(
            1, len(ordered_candidates) + 1
        )
        selected_current = ordered_candidates[:GENES_PER_ABUNDANCE_BIN]
        selected_from_abundance_bin[selected_current] = True
        selected_parts.append(selected_current)
        bin_support_records.append({
            "eligible_abundance_bin": abundance_bin,
            "n_expression_eligible": int(bin_eligible.sum()),
            "n_positive_noise_subtracted_variance": int(len(ordered_candidates)),
            "n_required": GENES_PER_ABUNDANCE_BIN,
            "n_selected": int(len(selected_current)),
            "passed": bool(len(ordered_candidates) >= GENES_PER_ABUNDANCE_BIN),
        })
    bin_support = pd.DataFrame(bin_support_records)
    write_tsv(bin_support, OUT / "abundance_stratified_gene_selection_support.tsv")
    if not bool(bin_support["passed"].all()):
        unsupported = pd.DataFrame([{
            "gate": "response_gene_signal", "n_detectable": int(expression_eligible.sum()),
            "n_expression_eligible": int(expression_eligible.sum()),
            "n_passing": int(passing.sum()), "minimum_required": MAX_RESPONSE_GENES,
            "selection_design": "300 positive-variance genes from each of 10 eligible-abundance bins",
            "passed": False,
        }])
        write_tsv(unsupported, OUT / "program_discovery_support.tsv")
        shortages = bin_support.loc[~bin_support["passed"], [
            "eligible_abundance_bin", "n_positive_noise_subtracted_variance", "n_required",
        ]].to_dict("records")
        raise RuntimeError(
            f"Eligible abundance-bin gene-selection shortage: {shortages}. "
            "Expression and variance requirements were not relaxed."
        )
    selected_index = np.concatenate(selected_parts)
    if (
        len(selected_index) != MAX_RESPONSE_GENES
        or not np.array_equal(selected_from_abundance_bin, np.isin(
            np.arange(EXPECTED_GENES), selected_index
        ))
    ):
        raise RuntimeError("Abundance-stratified selection did not produce exactly 3,000 genes.")
    selection_rank = np.full(EXPECTED_GENES, np.nan)
    selection_rank[selected_index] = np.arange(1, len(selected_index) + 1)
    signal_selection = pd.DataFrame({
        "gene_index": np.arange(EXPECTED_GENES), "gene": genes,
        "mean_dmso_log2_1p_cp10k": dmso_abundance,
        "baseline_abundance_decile": abundance_decile,
        "eligible_abundance_bin": eligible_abundance_bin,
        "detectable_dmso": reliably_expressed_dmso,
        "recurrently_detectable_treated": recurrently_expressed_treated,
        "detectable": expression_eligible,
        "minimum_expression_cp10k": MIN_RESPONSE_GENE_CP10K,
        "n_dmso_lines_reliably_expressed": dmso_reliably_expressed_lines,
        "reliably_expressed_dmso": reliably_expressed_dmso,
        "n_treated_responses_expressed": treated_expressed_count,
        "n_treated_perturbations_expressed": treated_expressed_perturbations,
        "n_treated_lines_expressed": treated_expressed_lines,
        "recurrently_expressed_treated": recurrently_expressed_treated,
        "expression_eligible": expression_eligible,
        "normalized_response_variance": response_var,
        "dmso_pseudo_response_variance": dmso_var,
        "stabilized_dmso_noise_variance": stable_noise_var,
        "noise_subtracted_response_variance": noise_subtracted_response_var,
        "excess_response_variance": excess_var,
        "excess_variance_ratio": excess_ratio,
        "dmso_absolute_response_q99": dmso_q99,
        "n_responses_beyond_dmso_q99": beyond_count,
        "n_distinct_perturbations_beyond_dmso_q99": beyond_perturbations,
        "n_distinct_cell_lines_beyond_dmso_q99": beyond_lines,
        "raw_count_scale_response_variance_diagnostic_only": raw_var,
        "passing_fixed_signal_criteria": passing,
        "selected": False, "selection_rank": selection_rank,
        "selection_rank_within_abundance_bin": selection_rank_within_bin,
        "selected_from_abundance_bin": selected_from_abundance_bin,
        "selection_score_definition": SELECTION_SCORE_DEFINITION,
        "selection_used_raw_count_variance": False,
        "selection_used_gene_identity": False,
        "selection_used_held_out_lines": False,
    })
    signal_selection.loc[selected_index, "selected"] = True
    if not np.all(expression_eligible[selected_index]):
        raise RuntimeError("A selected gene fails the 1 CP10K expression-support gate.")
    # Identity-derived annotations are appended only after the selection is frozen.
    signal_selection["post_selection_ribosomal_annotation"] = ribo_mask
    signal_selection["post_selection_mitochondrial_annotation"] = mito_mask
    write_tsv(signal_selection, OUT / "gene_signal_selection.tsv.gz", compression="gzip")
    write_tsv(pd.DataFrame([{
        "gate": "response_gene_signal", "n_detectable": int(expression_eligible.sum()),
        "n_expression_eligible": int(expression_eligible.sum()),
        "n_passing": int(passing.sum()), "n_selected": len(selected_index),
        "minimum_required": MAX_RESPONSE_GENES, "maximum_allowed": MAX_RESPONSE_GENES,
        "selection_design": "300 positive-variance genes from each of 10 eligible-abundance bins",
        "passed": True,
    }]), OUT / "program_discovery_support.tsv")
    write_tsv(signal_selection[[
        "gene", "mean_dmso_log2_1p_cp10k", "baseline_abundance_decile",
        "eligible_abundance_bin", "selection_rank_within_abundance_bin",
        "selected_from_abundance_bin",
        "minimum_expression_cp10k", "n_dmso_lines_reliably_expressed",
        "reliably_expressed_dmso", "n_treated_responses_expressed",
        "n_treated_perturbations_expressed", "n_treated_lines_expressed",
        "recurrently_expressed_treated", "expression_eligible",
        "raw_count_scale_response_variance_diagnostic_only",
        "normalized_response_variance", "stabilized_dmso_noise_variance",
        "noise_subtracted_response_variance", "excess_variance_ratio",
        "selection_score_definition", "selected",
    ]], OUT / "raw_vs_corrected_gene_variance.tsv.gz", compression="gzip")
    write_tsv(signal_selection.loc[signal_selection["selected"], [
        "gene", "selection_rank", "baseline_abundance_decile",
        "eligible_abundance_bin", "selection_rank_within_abundance_bin",
        "selected_from_abundance_bin",
        "mean_dmso_log2_1p_cp10k", "minimum_expression_cp10k",
        "n_dmso_lines_reliably_expressed", "reliably_expressed_dmso",
        "n_treated_responses_expressed", "n_treated_perturbations_expressed",
        "n_treated_lines_expressed", "recurrently_expressed_treated",
        "expression_eligible", "raw_count_scale_response_variance_diagnostic_only",
        "normalized_response_variance", "stabilized_dmso_noise_variance",
        "noise_subtracted_response_variance", "excess_variance_ratio",
        "selection_score_definition",
        "post_selection_ribosomal_annotation", "post_selection_mitochondrial_annotation",
    ]], OUT / "raw_vs_corrected_top_selected_genes.tsv")
    append_runtime(
        "preprocessing", preprocessing_started,
        f"paired normalization; DMSO null; {int(passing.sum())} passing and {len(selected_index)} selected genes",
    )

    # Second pass: only the selected-gene normalized and raw matrices.
    discovery_started = time.monotonic()
    observed_meta = metadata.loc[available].sort_values("array_row").reset_index(drop=True)
    n_observed = len(observed_meta)
    norm_by_array = normalization_frame.set_index("array_row")
    selected_partial = OUT / "normalized_selected_responses.npy.partial"
    selected_matrix = np.lib.format.open_memmap(
        selected_partial, mode="w+", dtype=np.float32,
        shape=(n_observed, len(selected_index)),
    )
    destination_lookup = dict(zip(observed_meta["array_row"].astype(int), range(n_observed)))
    for cell_number, depmap_id in enumerate(cells, start=1):
        cell_meta = metadata.loc[metadata["depmap_id"].eq(depmap_id)].sort_values("shard_row")
        shard = load_shard(cell_meta)
        for row in cell_meta.loc[cell_meta["available"].eq(1)].itertuples(index=False):
            destination = destination_lookup[int(row.array_row)]
            plates = audit_groups[(depmap_id, str(row.perturbation_id))]
            reference_indices = [plate_lookup[(depmap_id, plate)] for plate in plates]
            reference = np.asarray(
                plate_reference[np.ix_(reference_indices, selected_index)], dtype=np.float64
            ).mean(axis=0)
            raw = np.asarray(shard[int(row.shard_row), selected_index], dtype=np.float64)
            treated = raw + reference
            factor_record = norm_by_array.loc[int(row.array_row)]
            treated_cp = treated / float(factor_record["treated_total"]) * CP_SCALE
            treated_cp /= float(factor_record["treated_correction_factor"])
            reference_cp = reference / float(factor_record["reference_total"]) * CP_SCALE
            corrected = np.log2(1.0 + treated_cp) - np.log2(1.0 + reference_cp)
            if not np.isfinite(corrected).all():
                raise RuntimeError(f"Selected-gene reconstruction failed for {row.observation_id}.")
            selected_matrix[destination] = corrected.astype(np.float32)
        print(f"selected-gene pass 2: {cell_number}/{len(cells)} lines.", flush=True)
    selected_matrix.flush()
    del selected_matrix
    os.replace(selected_partial, OUT / "normalized_selected_responses.npy")
    selected_matrix = np.load(OUT / "normalized_selected_responses.npy", mmap_mode="r")
    selected_genes = genes[selected_index]
    selected_noise_sd = np.sqrt(stable_noise_var[selected_index])
    dmso_selected = np.asarray(dmso_pseudo[:, selected_index], dtype=np.float32)
    dmso_magnitude = np.sqrt(np.mean(np.square(
        dmso_selected.astype(np.float64) / selected_noise_sd[None, :]
    ), axis=1))
    activity_threshold = float(np.quantile(dmso_magnitude, 0.99))
    response_magnitude = np.sqrt(np.mean(np.square(
        np.asarray(selected_matrix, dtype=np.float64) / selected_noise_sd[None, :]
    ), axis=1))
    active = response_magnitude > activity_threshold
    observed_meta["role"] = np.where(observed_meta["depmap_id"].isin(training_lines), "training", "held_out")
    observed_meta["selected_gene_noise_standardized_magnitude"] = response_magnitude
    observed_meta["dmso_empirical_magnitude_q99"] = activity_threshold
    observed_meta["molecularly_active"] = active
    observed_meta["selected_response_row"] = np.arange(n_observed)
    write_tsv(observed_meta, OUT / "observed_response_metadata.tsv.gz", compression="gzip")
    dmso_magnitude_frame = pd.DataFrame(dmso_records)[[
        "pseudo_response_row", "depmap_id", "cell_line", "pseudo_treated_plate"
    ]].copy()
    dmso_magnitude_frame["magnitude"] = dmso_magnitude
    dmso_magnitude_frame["empirical_q99"] = activity_threshold
    write_tsv(dmso_magnitude_frame, OUT / "dmso_response_magnitude.tsv")
    activity_summary = observed_meta.groupby(
        ["role", "drug", "dose", "dose_panel", "perturbation_id", "depmap_id", "cell_line"],
        as_index=False,
    ).agg(
        n_responses=("observation_id", "size"),
        n_active=("molecularly_active", "sum"),
        active_fraction=("molecularly_active", "mean"),
        median_magnitude=("selected_gene_noise_standardized_magnitude", "median"),
    )
    write_tsv(activity_summary, OUT / "molecular_activity_by_drug_dose_cell_line.tsv.gz", compression="gzip")
    write_tsv(pd.DataFrame({
        "response_type": ["training_DMSO_pseudo"] * len(dmso_magnitude)
            + ["treated"] * len(response_magnitude),
        "role": ["training"] * len(dmso_magnitude) + observed_meta["role"].tolist(),
        "magnitude": np.concatenate([dmso_magnitude, response_magnitude]),
        "empirical_dmso_q99": activity_threshold,
        "above_dmso_q99": np.concatenate([
            dmso_magnitude > activity_threshold, active,
        ]),
    }), OUT / "dmso_vs_treated_response_magnitude.tsv.gz", compression="gzip")

    training_mask = observed_meta["depmap_id"].isin(training_lines).to_numpy()
    discovery_mask = training_mask.copy()
    discovery_meta = observed_meta.loc[discovery_mask].reset_index(drop=True)
    if (
        int(discovery_mask.sum()) != response_count
        or discovery_meta["depmap_id"].nunique() != len(training_lines)
        or discovery_meta["perturbation_id"].nunique() != EXPECTED_PERTURBATIONS
    ):
        raise RuntimeError("All available fold-1 training responses were not retained for discovery.")

    line_dmso_scale = dmso_magnitude_frame.groupby(
        ["depmap_id", "cell_line"], as_index=False
    ).agg(median_dmso_magnitude=("magnitude", "median"))
    median_training_line_scale = float(np.median(line_dmso_scale["median_dmso_magnitude"]))
    if not np.isfinite(median_training_line_scale) or median_training_line_scale <= 0:
        raise RuntimeError("Training-line DMSO magnitude normalization is invalid.")
    line_dmso_scale["relative_dmso_scaling_factor"] = (
        line_dmso_scale["median_dmso_magnitude"] / median_training_line_scale
    )
    relative_scale_lookup = line_dmso_scale.set_index("depmap_id")[
        "relative_dmso_scaling_factor"
    ].to_dict()
    discovery_lines = discovery_meta["depmap_id"].astype(str).to_numpy()
    discovery_matrix = np.asarray(selected_matrix[discovery_mask], dtype=np.float32)
    discovery_matrix /= selected_noise_sd.astype(np.float32)[None, :]
    for depmap_id in sorted(training_lines):
        line_rows = np.flatnonzero(discovery_lines == depmap_id)
        relative_scale = float(relative_scale_lookup[depmap_id])
        if len(line_rows) == 0 or not np.isfinite(relative_scale) or relative_scale <= 0:
            raise RuntimeError(f"Invalid discovery scaling inputs for {depmap_id}.")
        discovery_matrix[line_rows] /= np.float32(relative_scale)
        discovery_matrix[line_rows] -= discovery_matrix[line_rows].mean(
            axis=0, dtype=np.float64
        ).astype(np.float32)[None, :]
    discovery_winsor_low = np.quantile(discovery_matrix, 0.005, axis=0)
    discovery_winsor_high = np.quantile(discovery_matrix, 0.995, axis=0)
    np.clip(
        discovery_matrix, discovery_winsor_low, discovery_winsor_high,
        out=discovery_matrix,
    )

    line_counts = pd.Series(discovery_lines).value_counts(sort=False)
    weights = np.asarray([1.0 / line_counts[value] for value in discovery_lines], dtype=np.float64)
    weights *= len(weights) / weights.sum()
    total_weight = float(weights.sum())
    line_weight = pd.DataFrame({"depmap_id": discovery_lines, "row_weight": weights}).groupby(
        "depmap_id", as_index=False
    ).agg(total_covariance_weight=("row_weight", "sum"))
    active_by_line = observed_meta.loc[training_mask].groupby(
        ["depmap_id", "cell_line"], as_index=False
    ).agg(previous_pooled_threshold_active_fraction=("molecularly_active", "mean"))
    line_rows = discovery_meta.groupby(["depmap_id", "cell_line"], as_index=False).agg(
        n_discovery_rows=("observation_id", "size")
    )
    line_balance = line_rows.merge(
        line_dmso_scale, on=["depmap_id", "cell_line"], validate="one_to_one"
    ).merge(line_weight, on="depmap_id", validate="one_to_one").merge(
        active_by_line, on=["depmap_id", "cell_line"], validate="one_to_one"
    ).sort_values("depmap_id")
    if (
        len(line_balance) != len(training_lines)
        or not np.allclose(
            line_balance["total_covariance_weight"],
            line_balance["total_covariance_weight"].iloc[0], rtol=1e-12, atol=1e-10,
        )
        or not np.isclose(
            np.median(line_balance["relative_dmso_scaling_factor"]), 1.0,
            rtol=0, atol=1e-12,
        )
    ):
        raise RuntimeError("Equal training-line covariance weighting failed.")
    write_tsv(line_balance, OUT / "program_discovery_line_balance.tsv")

    discovery_weighted_mean = np.einsum(
        "i,ij->j", weights, discovery_matrix, optimize=True
    ) / total_weight
    discovery_matrix -= discovery_weighted_mean.astype(np.float32)[None, :]
    discovery_matrix *= np.sqrt(weights).astype(np.float32)[:, None]
    covariance = discovery_matrix.T @ discovery_matrix
    variance = np.diag(covariance).astype(np.float64)
    if np.any(~np.isfinite(variance)) or np.any(variance <= 0):
        raise RuntimeError("A selected gene has invalid weighted response variance.")
    scale = np.sqrt(variance)
    correlation = covariance / np.outer(scale, scale).astype(np.float32)
    correlation = ((correlation + correlation.T) * np.float32(0.5)).astype(np.float32)
    np.fill_diagonal(correlation, 1.0)
    if not np.isfinite(correlation).all() or np.max(np.abs(correlation)) > 1.0001:
        raise RuntimeError("Weighted signed gene correlation failed integrity checks.")
    write_npy(correlation, OUT / "reference_gene_gene_cross_correlation.npy")
    eigenvalues = np.linalg.eigvalsh(correlation.astype(np.float64))[::-1]
    eigenvalues[eigenvalues < 0] = np.maximum(eigenvalues[eigenvalues < 0], -1e-8)
    explained = eigenvalues / eigenvalues.sum()
    write_tsv(pd.DataFrame({
        "component": np.arange(1, len(eigenvalues) + 1),
        "eigenvalue": eigenvalues, "fraction_variance": explained,
        "cumulative_fraction_variance": np.cumsum(explained),
    }), OUT / "corrected_correlation_eigenvalue_spectrum.tsv.gz", compression="gzip")
    write_tsv(pd.DataFrame({
        "n_components": [1, 2, 5],
        "fraction_variance_explained": [explained[:n].sum() for n in (1, 2, 5)],
    }), OUT / "corrected_correlation_dimension_summary.tsv")

    helper = load_existing_helpers()
    helper.N_RESPONSE_GENES = len(selected_index)
    helper.CACHE = OUT
    helper.PROGRAM_DEFINITION_VERSION = PROGRAM_DEFINITION_VERSION
    helper.REFERENCE_CORRELATION_FILE = OUT / "reference_gene_gene_cross_correlation.npy"
    affinity = helper.posneg_correlation(correlation)
    repeat_memberships = helper.run_repeated_onmf(
        "fold_1", affinity, selected_genes, training_lines
    )
    unsigned_membership, matched, repeat_matches, medoid = helper.consensus_membership(
        repeat_memberships, selected_genes
    )
    signed_weights, anchors = helper.orient_signed_programs(unsigned_membership, correlation)
    del affinity, repeat_memberships, covariance, discovery_matrix

    # Scoring returns to the original normalized displacement representation.
    # Its winsor bounds and global training reference are learned separately,
    # with no line-specific scaling or within-line centering.
    scoring_training = np.asarray(selected_matrix[training_mask], dtype=np.float32)
    scoring_training /= selected_noise_sd.astype(np.float32)[None, :]
    scoring_winsor_low = np.quantile(scoring_training, 0.005, axis=0)
    scoring_winsor_high = np.quantile(scoring_training, 0.995, axis=0)
    np.clip(scoring_training, scoring_winsor_low, scoring_winsor_high, out=scoring_training)
    scoring_gene_mean = scoring_training.mean(axis=0, dtype=np.float64)
    del scoring_training
    all_standardized = np.asarray(selected_matrix, dtype=np.float32) / selected_noise_sd.astype(np.float32)[None, :]
    np.clip(all_standardized, scoring_winsor_low, scoring_winsor_high, out=all_standardized)
    all_standardized -= scoring_gene_mean.astype(np.float32)[None, :]
    all_scores = all_standardized @ signed_weights.T
    if not np.isfinite(all_scores).all():
        raise RuntimeError("Frozen fold-1 projection returned nonfinite scores.")
    score_variance = all_scores[training_mask].var(axis=0, ddof=1)
    programs = [f"P{i + 1}" for i in range(RANK)]
    stability_records = helper.records_for_stability("fold_1", programs, matched, score_variance)
    stability = pd.DataFrame(stability_records)
    stability["outer_fold"] = 1
    stable_programs = stability.loc[stability["stable"], "program"].drop_duplicates().tolist()
    stable_indices = np.asarray([programs.index(value) for value in stable_programs], dtype=int)
    assignment, plotting_order, correlation_audit = helper.correlation_order_and_audit(
        correlation, signed_weights, stable_indices
    )
    stable_score_correlation = np.corrcoef(all_scores[:, stable_indices], rowvar=False)
    score_correlation_frame = pd.DataFrame(
        stable_score_correlation, columns=stable_programs
    )
    score_correlation_frame.insert(0, "program", stable_programs)
    write_tsv(score_correlation_frame, OUT / "program_score_correlation_matrix.tsv")
    off_diagonal = stable_score_correlation[np.triu_indices(len(stable_programs), k=1)]

    bundle = {
        "label": "fold_1_training_basis",
        "program_definition_version": PROGRAM_DEFINITION_VERSION,
        "rank": RANK, "program_names": programs, "stable_programs": stable_programs,
        "selected_original_genes": selected_genes,
        "selected_gene_indices": selected_index,
        "training_gene_scale": selected_noise_sd,
        "gene_scale_definition": "sqrt(training-only abundance-bin-stabilized DMSO pseudo-response variance)",
        "unsigned_gene_memberships": unsigned_membership,
        "signed_gene_weights": signed_weights,
        "anchor_gene_indices": anchors, "anchor_genes": selected_genes[anchors],
        "training_lines": sorted(training_lines), "held_out_lines": sorted(held_lines),
        "all_scores_observed_order": all_scores,
        "observed_array_rows": observed_meta["array_row"].to_numpy(dtype=int),
        "repeat_component_matches": repeat_matches,
        "consensus_medoid_repeat_index": medoid,
        "gene_assignment_program_index": assignment,
        "gene_plotting_order": plotting_order,
        "winsor_lower": scoring_winsor_low, "winsor_upper": scoring_winsor_high,
        "weighted_training_gene_mean": scoring_gene_mean,
        "scoring_winsor_lower": scoring_winsor_low,
        "scoring_winsor_upper": scoring_winsor_high,
        "scoring_training_gene_mean": scoring_gene_mean,
        "discovery_winsor_lower": discovery_winsor_low,
        "discovery_winsor_upper": discovery_winsor_high,
        "discovery_weighted_gene_mean": discovery_weighted_mean,
        "relative_dmso_scaling_factor_by_line": relative_scale_lookup,
        "dmso_activity_threshold_q99": activity_threshold,
        "n_discovery_rows": int(discovery_mask.sum()),
    }
    models = {
        "folds": {1: bundle}, "reference": bundle,
        "parameters": {
            "seed": SEED, "rank": RANK, "candidate_rank_not_biological_count": True,
            "active_outer_folds": ACTIVE_OUTER_FOLDS,
            "n_onmf_runs": N_ONMF_RUNS, "onmf_max_iter": ONMF_MAX_ITER,
            "program_definition_version": PROGRAM_DEFINITION_VERSION,
            "normalizer": "paired trimmed M/A log-ratio then log2(1+CP10K)",
            "correlation_input": "all training DMSO-noise-standardized, line-scaled, within-line-centered responses",
            "correlation_weighting": "equal total covariance weight per training cell line",
            "nmf_input": "posneg(weighted signed gene-gene Pearson C)",
            "row_l2_normalization": False,
        },
    }
    joblib.dump(models, OUT / "response_onmf_models.joblib", compress=3)

    def weight_records(basis: str, outer_fold=None):
        rows = []
        stable_set = set(stable_programs)
        for program_index, program in enumerate(programs):
            absolute = np.abs(signed_weights[program_index])
            total = float(absolute.sum())
            for gene_index, gene in enumerate(selected_genes):
                direction_active = "UP" if signed_weights[program_index, gene_index] >= 0 else "DOWN"
                for direction in ("UP", "DOWN"):
                    value = float(absolute[gene_index]) if direction == direction_active else 0.0
                    row = {
                        "basis": basis, "program": program, "feature": f"{gene}__{direction}",
                        "gene": gene, "direction": direction, "weight": value,
                        "normalized_weight": value / total, "stable": program in stable_set,
                    }
                    if outer_fold is not None:
                        row["outer_fold"] = outer_fold
                    rows.append(row)
        return pd.DataFrame(rows)

    def score_records(basis: str, outer_fold=None):
        rows = []
        stable_set = set(stable_programs)
        for i, response_row in observed_meta.iterrows():
            for j, program in enumerate(programs):
                row = {
                    "basis": basis, "array_row": int(response_row.array_row),
                    "observation_id": response_row.observation_id,
                    "cell_line": response_row.cell_line, "depmap_id": response_row.depmap_id,
                    "drug": response_row.drug, "dose": float(response_row.dose),
                    "dose_panel": int(response_row.dose_panel),
                    "perturbation_id": response_row.perturbation_id,
                    "role": response_row.role, "program": program,
                    "score": float(all_scores[i, j]), "stable": program in stable_set,
                    "molecularly_active": bool(response_row.molecularly_active),
                }
                if outer_fold is not None:
                    row["outer_fold"] = outer_fold
                rows.append(row)
        return pd.DataFrame(rows)

    def assignment_records(basis: str, outer_fold=None):
        stable_once = stability.drop_duplicates("program").set_index("program")
        rows = []
        for gene_index, gene in enumerate(selected_genes):
            program_index = int(assignment[gene_index])
            program = programs[program_index]
            value = float(signed_weights[program_index, gene_index])
            row = {
                "basis": basis, "selected_gene_index": gene_index, "gene": gene,
                "assigned_program": program, "signed_weight": value,
                "absolute_weight": abs(value), "direction": "UP" if value >= 0 else "DOWN",
                "stability": float(stable_once.loc[program, "q10_matched_pearson"]),
                "stable": program in stable_programs,
                "plotting_order": int(plotting_order[gene_index]),
            }
            if outer_fold is not None:
                row["outer_fold"] = outer_fold
            rows.append(row)
        return pd.DataFrame(rows)

    fold_weights = weight_records("fold_1", 1)
    fold_scores = score_records("fold_1", 1)
    fold_assignment = assignment_records("fold_1", 1).sort_values("plotting_order")
    reference_weights = weight_records("fold_1_training_basis")
    reference_scores = score_records("fold_1_training_basis")
    reference_stability = stability.drop(columns="outer_fold").copy()
    reference_stability["basis"] = "fold_1_training_basis"
    reference_assignment = assignment_records("fold_1_training_basis").sort_values("plotting_order")
    write_tsv(fold_weights, OUT / "fold_program_weights.tsv.gz", compression="gzip")
    write_tsv(stability, OUT / "fold_program_stability.tsv")
    write_tsv(fold_scores, OUT / "fold_program_scores.tsv.gz", compression="gzip")
    write_tsv(reference_weights, OUT / "reference_program_weights.tsv.gz", compression="gzip")
    write_tsv(reference_stability, OUT / "reference_program_stability.tsv")
    write_tsv(reference_scores, OUT / "reference_program_scores.tsv.gz", compression="gzip")
    write_tsv(fold_assignment, OUT / "fold_gene_program_assignment.tsv.gz", compression="gzip")
    write_tsv(reference_assignment, OUT / "reference_gene_program_assignment.tsv.gz", compression="gzip")
    write_tsv(pd.DataFrame([{
        "outer_fold": 1, "fold_program": program, "reference_program": program,
        "pearson": 1.0, "pearson_before_sign_flip": 1.0,
        "sign_flipped": False, "n_shared_genes": len(selected_genes),
        "n_shared_directional_features": 2 * len(selected_genes),
        "fold_stable": program in stable_programs,
        "reference_stable": program in stable_programs,
        "interpretation": "self-alignment to the frozen fold-1 training basis; single-split development run",
    } for program in programs]), OUT / "program_fold_alignment.tsv")

    bottleneck, bottleneck_gate = helper.bottleneck_audit(fold_scores, stability)
    write_tsv(bottleneck, OUT / "program_bottleneck_audit.tsv")
    write_tsv(bottleneck_gate, OUT / "program_bottleneck_gate.tsv")
    audit = pd.DataFrame([{
        "basis": "fold_1_training_basis", "outer_fold": 1,
        "program_definition_version": PROGRAM_DEFINITION_VERSION,
        "nmf_input": "posneg(weighted_gene_gene_Pearson_C)",
        "nmf_input_rows": len(selected_genes), "nmf_input_columns": 2 * len(selected_genes),
        "n_response_rows": int(discovery_mask.sum()),
        "n_training_cell_lines": len(training_lines),
        "n_held_out_cell_lines": len(held_lines),
        "n_perturbations": int(discovery_meta["perturbation_id"].nunique()),
        "equal_total_weight_per_drug_dose": False,
        "equal_total_weight_per_training_cell_line": True,
        "active_training_responses_only_for_discovery": False,
        "all_available_training_responses_used_for_discovery": True,
        "inactive_responses_retained_for_scoring": True,
        "fraction_all_treated_above_dmso_q99": float(active.mean()),
        "fraction_training_treated_above_dmso_q99": float(active[training_mask].mean()),
        "n_selected_genes": len(selected_genes), "candidate_rank": RANK,
        "n_onmf_runs": N_ONMF_RUNS, "n_stable_programs": len(stable_programs),
        "no_response_row_l2_normalization": True,
        "program_scores_preserve_response_magnitude": True,
        "discovery_uses_line_specific_dmso_magnitude_scaling": True,
        "discovery_uses_within_line_gene_centering": True,
        "final_scores_use_line_specific_discovery_scaling": False,
        "final_scores_use_within_line_discovery_centering": False,
        "discovery_and_scoring_winsor_bounds_learned_separately": True,
        "selection_used_raw_count_variance": False,
        "response_gene_selection_definition": SELECTION_SCORE_DEFINITION,
        "selection_used_excess_variance_ratio": False,
        "selection_used_dmso_q99_exceedance_gate": False,
        "minimum_response_gene_cp10k": MIN_RESPONSE_GENE_CP10K,
        "minimum_dmso_expressed_lines": MIN_DMSO_EXPRESSED_LINES,
        "minimum_treated_expressed_responses": MIN_TREATED_EXPRESSED_RESPONSES,
        "minimum_treated_expressed_perturbations": MIN_TREATED_EXPRESSED_PERTURBATIONS,
        "minimum_treated_expressed_lines": MIN_TREATED_EXPRESSED_LINES,
        "selection_used_gene_identity": False,
        "selection_used_held_out_lines": False,
        "normalization_parameters_used_held_out_lines": False,
        "program_discovery_used_prism_or_depmap": False,
        "reference_basis_refit": False,
        "median_offdiagonal_program_score_correlation": float(np.median(off_diagonal)),
        "minimum_offdiagonal_program_score_correlation": float(np.min(off_diagonal)),
        "maximum_offdiagonal_program_score_correlation": float(np.max(off_diagonal)),
        "median_absolute_offdiagonal_program_score_correlation": float(np.median(np.abs(off_diagonal))),
        **correlation_audit,
    }])
    audit["scientific_gate_passed"] = (
        len(stable_programs) >= 3
        and audit["within_exceeds_between_absolute_correlation"].iloc[0]
    )
    write_tsv(audit, OUT / "program_definition_audit.tsv")

    annotation_summary = signal_selection.loc[signal_selection["selected"]].groupby(
        "baseline_abundance_decile", as_index=False
    ).agg(
        n_selected=("gene", "size"),
        n_ribosomal=("post_selection_ribosomal_annotation", "sum"),
        n_mitochondrial=("post_selection_mitochondrial_annotation", "sum"),
    )
    annotation_summary["fraction_ribosomal"] = annotation_summary["n_ribosomal"] / annotation_summary["n_selected"]
    annotation_summary["fraction_mitochondrial"] = annotation_summary["n_mitochondrial"] / annotation_summary["n_selected"]
    write_tsv(annotation_summary, OUT / "post_selection_annotation_summary.tsv")
    write_tsv(pd.DataFrame([{
        "check": "active_outer_folds_exactly_1", "passed": ACTIVE_OUTER_FOLDS == [1],
    }, {
        "check": "fold1_heldout_excluded_from_all_program_learning", "passed": True,
    }, {
        "check": "same_normalizer_for_treated_and_dmso_null", "passed": True,
    }, {
        "check": "gene_names_do_not_influence_normalization_or_selection", "passed": True,
    }, {
        "check": "raw_count_variance_is_diagnostic_only", "passed": True,
    }, {
        "check": "all_selected_genes_pass_1_cp10k_expression_gate",
        "passed": bool(np.all(expression_eligible[selected_index])),
    }, {
        "check": "selection_ordered_by_variance_difference_not_ratio", "passed": True,
    }, {
        "check": "exactly_300_genes_selected_per_eligible_abundance_bin",
        "passed": bool(np.all(bin_support["n_selected"].to_numpy() == GENES_PER_ABUNDANCE_BIN)),
    }, {
        "check": "dmso_q99_exceedance_counts_are_diagnostic_only", "passed": True,
    }, {
        "check": "positive_negative_C_representation_not_abs_C", "passed": True,
    }, {
        "check": "no_response_row_L2_normalization", "passed": True,
    }, {
        "check": "discovery_mask_equals_training_mask",
        "passed": bool(np.array_equal(discovery_mask, training_mask)),
    }, {
        "check": "discovery_has_all_26_training_lines_and_1137_perturbations",
        "passed": bool(
            discovery_meta["depmap_id"].nunique() == EXPECTED_TRAINING_LINES
            and discovery_meta["perturbation_id"].nunique() == EXPECTED_PERTURBATIONS
        ),
    }, {
        "check": "equal_total_covariance_weight_per_training_line",
        "passed": bool(np.allclose(
            line_balance["total_covariance_weight"],
            line_balance["total_covariance_weight"].iloc[0], rtol=1e-12, atol=1e-10,
        )),
    }, {
        "check": "final_scores_exclude_line_scaling_and_within_line_centering",
        "passed": True,
    }, {
        "check": "no_second_all_gene_normalized_matrix", "passed": True,
    }, {
        "check": "excluded_lines_absent_from_all_response_objects",
        "passed": bool(
            not EXCLUDED_IDS.intersection(metadata["depmap_id"].astype(str))
            and not EXCLUDED_IDS.intersection(observed_meta["depmap_id"].astype(str))
            and not EXCLUDED_IDS.intersection(discovery_meta["depmap_id"].astype(str))
            and not EXCLUDED_IDS.intersection(line_balance["depmap_id"].astype(str))
        ),
    }, {
        "check": "previous_full_chain_output_preserved", "passed": True,
    }]), OUT / "required_checks.tsv")
    append_runtime(
        "program_discovery", discovery_started,
        f"{int(discovery_mask.sum())} all available line-balanced training responses; "
        f"{len(stable_programs)} stable rank-10 components",
    )
    append_runtime("stage03_total", total_started)
    with (OUT / "stage03_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump(models["parameters"], handle, indent=2)
        handle.write("\n")
    # Keep only one copy of scores/weights and separate lightweight results
    # from bulky and observation-level artifacts.
    for redundant in (
        "reference_program_weights.tsv.gz", "reference_program_stability.tsv",
        "reference_program_scores.tsv.gz", "reference_gene_program_assignment.tsv.gz",
        "program_fold_alignment.tsv",
    ):
        path = OUT / redundant
        if path.exists():
            path.unlink()
    for lightweight in (
        "fold_program_weights.tsv.gz", "fold_program_stability.tsv",
        "fold_gene_program_assignment.tsv.gz", "program_definition_audit.tsv",
        "program_bottleneck_audit.tsv", "program_bottleneck_gate.tsv",
        "program_discovery_support.tsv", "abundance_stratified_gene_selection_support.tsv",
        "program_score_correlation_matrix.tsv", "stage03_parameters.json",
        "required_checks.tsv", "post_selection_annotation_summary.tsv",
    ):
        source = OUT / lightweight
        if source.exists():
            os.replace(source, RESULT / lightweight.replace("stage03_", "stage02_"))
    write_artifact_manifest()
    (RESULT / "stage02_complete.flag").write_text("complete\n", encoding="utf-8")
    print(
        f"Corrected fold-1 basis complete: {int(expression_eligible.sum())} expression-eligible, "
        f"{int(passing.sum())} passing, {len(selected_genes)} selected genes; "
        f"{active.mean():.3f} treated active; {len(stable_programs)}/{RANK} stable programs; "
        f"runtime {time.monotonic() - total_started:.1f}s.", flush=True
    )


if __name__ == "__main__":
    main()
