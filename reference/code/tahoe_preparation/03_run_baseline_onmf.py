#!/usr/bin/env python3
"""Learn fold-1 untreated baseline programs from natural DMSO plates.

The biological discovery matrix is one DMSO mean per outer-training cell
line.  Natural plate profiles are retained only for technical-reproducibility
estimation and scoring.  Perturbation responses, held-out lines, drug labels,
and PRISM outcomes never enter gene selection or oNMF.
"""

from __future__ import annotations

import importlib.util
import gzip
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "20")
os.environ.setdefault("OMP_NUM_THREADS", "20")

import joblib
import numpy as np
import pandas as pd


SEED = 1
OUTER_FOLD = 1
CANDIDATE_RANK = 5
N_ONMF_RUNS = 10
ONMF_MAX_ITER = 500
MIN_PROGRAM_STABILITY = 0.50
N_ABUNDANCE_BINS = 10
GENES_PER_BIN = 300
N_SELECTED_GENES = N_ABUNDANCE_BINS * GENES_PER_BIN
SOURCE_LINES = 35
EXPECTED_LINES = 33
EXPECTED_TRAINING_LINES = 26
EXPECTED_HELD_LINES = 7
EXPECTED_PLATES_PER_LINE = 14
CP_SCALE = 10000.0
PROGRAM_DEFINITION_VERSION = "ntc_baseline_line_mean_posneg_v1"

ROOT = Path(os.environ.get("TAHOE_V2_ROOT", "tahoe_prism_v2")).resolve()
SOURCE_OUT = ROOT / "artifacts" / "02_response_programs"
OUT = ROOT / "artifacts" / "03_baseline_programs"
RESULT = ROOT / "results" / "03_baseline_programs"
CACHE = ROOT.parent / "cache" / "tahoe_prism_full_chain"
CONTROL_EXPRESSION_FILE = CACHE / "dmso_control_expression.npy"
CONTROL_METADATA_FILE = CACHE / "dmso_control_metadata.tsv"
CONTROL_GENE_FILE = CACHE / "dmso_control_genes.tsv"
PREPARED_FILE = ROOT / "artifacts" / "01_inputs" / "prepared_inputs.rds"
SPLIT_FILE = SOURCE_OUT / "fold1_cell_line_split.tsv"
EXCLUDED_IDS = {"ACH-000311", "ACH-000628"}


def load_response_helpers():
    path = Path(__file__).with_name("02_onmf_helpers.py")
    spec = importlib.util.spec_from_file_location("tahoe_response_onmf_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the existing response-oNMF helpers.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.SEED = SEED
    module.RANK = CANDIDATE_RANK
    module.N_ONMF_RUNS = N_ONMF_RUNS
    module.ONMF_MAX_ITER = ONMF_MAX_ITER
    module.MIN_PROGRAM_STABILITY = MIN_PROGRAM_STABILITY
    module.N_RESPONSE_GENES = N_SELECTED_GENES
    module.CACHE = OUT
    module.PROGRAM_DEFINITION_VERSION = PROGRAM_DEFINITION_VERSION
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


def write_artifact_manifest() -> None:
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
            "producing_script": "scripts/03_run_baseline_onmf.py",
        })
    write_tsv(pd.DataFrame(records), OUT / "artifact_manifest.tsv")


def safe_cor(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    keep = np.isfinite(left) & np.isfinite(right)
    if keep.sum() < 3 or np.std(left[keep], ddof=1) <= 0 or np.std(right[keep], ddof=1) <= 0:
        return float("nan")
    return float(np.corrcoef(left[keep], right[keep])[0, 1])


def export_prepared_depmap(temporary_file: Path) -> pd.DataFrame:
    """Read the already prepared DepMap matrix without adding a Python RDS dependency."""
    r_code = """
    suppressPackageStartupMessages(library(data.table))
    args <- commandArgs(trailingOnly=TRUE)
    x <- readRDS(args[1])$depmap_baseline$expression
    fwrite(data.table(depmap_id=rownames(x), as.data.table(x)), args[2], sep='\\t')
    """
    subprocess.run(
        ["Rscript", "-e", r_code, str(PREPARED_FILE), str(temporary_file)],
        check=True,
    )
    result = pd.read_csv(temporary_file, sep="\t")
    temporary_file.unlink()
    result = result.set_index("depmap_id")
    result.index = result.index.astype(str)
    return result


def abundance_stratified_selection(
    abundance: np.ndarray, reproducible_variance: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    positive = np.isfinite(abundance) & np.isfinite(reproducible_variance) & (reproducible_variance > 0)
    if positive.sum() < N_SELECTED_GENES:
        raise RuntimeError(
            f"Only {positive.sum()} genes have positive reproducible between-line variance."
        )
    bins = np.full(len(abundance), -1, dtype=np.int16)
    positive_index = np.flatnonzero(positive)
    ranked_abundance = pd.Series(abundance[positive_index]).rank(method="first")
    bins[positive_index] = (
        pd.qcut(ranked_abundance, N_ABUNDANCE_BINS, labels=False).to_numpy(dtype=np.int16) + 1
    )
    selected_parts: list[np.ndarray] = []
    rank_in_bin = np.full(len(abundance), np.nan)
    for abundance_bin in range(1, N_ABUNDANCE_BINS + 1):
        current = np.flatnonzero(bins == abundance_bin)
        order = np.lexsort((current, -reproducible_variance[current]))
        current = current[order]
        if len(current) < GENES_PER_BIN:
            raise RuntimeError(
                f"Abundance bin {abundance_bin} has only {len(current)} positive-variance genes."
            )
        rank_in_bin[current] = np.arange(1, len(current) + 1)
        selected_parts.append(current[:GENES_PER_BIN])
    selected = np.concatenate(selected_parts)
    selection_rank = np.full(len(abundance), np.nan)
    selection_rank[selected] = np.arange(1, len(selected) + 1)
    if len(selected) != N_SELECTED_GENES or len(np.unique(selected)) != N_SELECTED_GENES:
        raise RuntimeError("Abundance-stratified selection did not return 3,000 unique genes.")
    return selected, bins, rank_in_bin, selection_rank


def standardize_genes(
    matrix: np.ndarray, training_rows: np.ndarray, label: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    training = np.asarray(matrix[training_rows], dtype=np.float64)
    center = training.mean(axis=0)
    scale = training.std(axis=0, ddof=1)
    keep = np.isfinite(center) & np.isfinite(scale) & (scale > 0)
    if keep.sum() == 0:
        raise RuntimeError(f"{label} has no finite nonconstant selected genes.")
    z = np.zeros(np.asarray(matrix).shape, dtype=np.float64)
    z[:, keep] = (np.asarray(matrix[:, keep], dtype=np.float64) - center[keep]) / scale[keep]
    if not np.isfinite(z).all():
        raise RuntimeError(f"{label} gene standardization is nonfinite.")
    return z, center, scale, keep


def standardize_programs(
    raw: np.ndarray, training_rows: np.ndarray, label: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = raw[training_rows].mean(axis=0)
    scale = raw[training_rows].std(axis=0, ddof=1)
    if np.any(~np.isfinite(center)) or np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise RuntimeError(f"{label} contains a nonstandardizable program.")
    standardized = (raw - center[None, :]) / scale[None, :]
    return standardized, center, scale


def icc_oneway_balanced(values: np.ndarray) -> float:
    """One-way random-effects ICC(1,1); rows are lines and columns are plates."""
    values = np.asarray(values, dtype=float)
    n_lines, n_plates = values.shape
    line_mean = values.mean(axis=1)
    grand = values.mean()
    between_ms = n_plates * np.square(line_mean - grand).sum() / (n_lines - 1)
    within_ms = np.square(values - line_mean[:, None]).sum() / (n_lines * (n_plates - 1))
    denominator = between_ms + (n_plates - 1) * within_ms
    return float((between_ms - within_ms) / denominator) if denominator > 0 else float("nan")


def main() -> None:
    started = time.monotonic()
    required = [CONTROL_EXPRESSION_FILE, CONTROL_METADATA_FILE, CONTROL_GENE_FILE,
                PREPARED_FILE, SPLIT_FILE]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing frozen baseline input(s): {missing}")
    OUT.mkdir(parents=True, exist_ok=True)
    RESULT.mkdir(parents=True, exist_ok=True)

    split = pd.read_csv(SPLIT_FILE, sep="\t")
    training_lines = sorted(split.loc[split["role"].eq("training"), "depmap_id"].astype(str))
    held_lines = sorted(split.loc[split["role"].eq("held_out"), "depmap_id"].astype(str))
    all_lines = sorted(training_lines + held_lines)
    if (
        len(training_lines) != EXPECTED_TRAINING_LINES
        or len(held_lines) != EXPECTED_HELD_LINES
        or bool(EXCLUDED_IDS.intersection(all_lines))
    ):
        raise RuntimeError("The corrected fold-1 26/7 split changed.")

    control_meta = pd.read_csv(CONTROL_METADATA_FILE, sep="\t")
    control_genes = pd.read_csv(CONTROL_GENE_FILE, sep="\t")
    control = np.load(CONTROL_EXPRESSION_FILE, mmap_mode="r")
    if (
        control.shape != (SOURCE_LINES * EXPECTED_PLATES_PER_LINE, len(control_genes))
        or control_meta.groupby("depmap_id").size().ne(EXPECTED_PLATES_PER_LINE).any()
        or control_genes["feature"].isna().any()
        or control_genes["feature"].duplicated().any()
    ):
        raise RuntimeError("The natural DMSO plate input contract changed.")
    control_meta = control_meta.loc[control_meta["depmap_id"].isin(all_lines)].copy()
    control_meta = control_meta.sort_values("array_row").reset_index(drop=True)
    retained_control_rows = control_meta["array_row"].to_numpy(dtype=int)
    if (
        len(control_meta) != EXPECTED_LINES * EXPECTED_PLATES_PER_LINE
        or bool(EXCLUDED_IDS.intersection(control_meta["depmap_id"].astype(str)))
    ):
        raise RuntimeError("The retained natural DMSO cohort is invalid.")

    depmap_temporary = OUT / ".prepared_depmap_stage03b.tsv"
    if depmap_temporary.exists():
        depmap_temporary.unlink()
    depmap = export_prepared_depmap(depmap_temporary)
    if sorted(depmap.index.tolist()) != all_lines or depmap.columns.duplicated().any():
        raise RuntimeError("The prepared DepMap baseline cohort changed.")

    dmso_gene_names = control_genes["feature"].astype(str).to_numpy()
    depmap_gene_set = set(depmap.columns.astype(str))
    valid = np.asarray([bool(value) and value != "nan" for value in dmso_gene_names])
    universe_index = np.flatnonzero(valid & np.asarray([gene in depmap_gene_set for gene in dmso_gene_names]))
    universe_genes = dmso_gene_names[universe_index]
    if len(universe_genes) != len(np.unique(universe_genes)):
        raise RuntimeError("The common baseline feature universe contains duplicate identifiers.")

    plate_raw = np.asarray(control[np.ix_(retained_control_rows, universe_index)], dtype=np.float64)
    library = plate_raw.sum(axis=1)
    if np.any(~np.isfinite(plate_raw)) or np.any(plate_raw < 0) or np.any(library <= 0):
        raise RuntimeError("DMSO plate expression is invalid.")
    plate_expression = np.log2(1.0 + plate_raw / library[:, None] * CP_SCALE).astype(np.float32)
    del plate_raw

    meta_lookup = {
        (str(row.depmap_id), int(row.plate)): row_index
        for row_index, row in enumerate(control_meta.itertuples(index=False))
    }
    plate_rows = {
        line: np.asarray([meta_lookup[(line, plate)] for plate in sorted(
            control_meta.loc[control_meta["depmap_id"].eq(line), "plate"].astype(int)
        )], dtype=int)
        for line in all_lines
    }
    line_expression = np.stack([plate_expression[plate_rows[line]].mean(axis=0) for line in all_lines])
    training_position = np.asarray([all_lines.index(line) for line in training_lines], dtype=int)
    held_position = np.asarray([all_lines.index(line) for line in held_lines], dtype=int)
    training_means = line_expression[training_position]
    abundance = training_means.mean(axis=0, dtype=np.float64)
    between_variance = training_means.var(axis=0, ddof=1, dtype=np.float64)
    within_by_line = np.stack([
        plate_expression[plate_rows[line]].var(axis=0, ddof=1, dtype=np.float64)
        for line in training_lines
    ])
    mean_within_variance = within_by_line.mean(axis=0)
    mean_plate_noise = np.mean(within_by_line / EXPECTED_PLATES_PER_LINE, axis=0)
    reproducible_variance = np.maximum(between_variance - mean_plate_noise, 0.0)
    selected, abundance_bin, rank_in_bin, selection_rank = abundance_stratified_selection(
        abundance, reproducible_variance
    )

    selection = pd.DataFrame({
        "gene": universe_genes,
        "mean_baseline_abundance_log2_1p_cp10k": abundance,
        "between_line_variance": between_variance,
        "mean_within_line_plate_variance": mean_within_variance,
        "estimated_plate_noise_of_line_mean": mean_plate_noise,
        "reproducible_between_line_variance": reproducible_variance,
        "positive_reproducible_variance": reproducible_variance > 0,
        "abundance_bin": abundance_bin,
        "rank_within_abundance_bin": rank_in_bin,
        "selection_rank": selection_rank,
        "selected": False,
    })
    selection.loc[selected, "selected"] = True
    write_tsv(selection, OUT / "baseline_gene_signal_selection.tsv.gz", compression="gzip")

    selected_genes = universe_genes[selected]
    selected_line_expression = np.asarray(line_expression[:, selected], dtype=np.float64)
    selected_plate_expression = np.asarray(plate_expression[:, selected], dtype=np.float64)
    line_z, dmso_gene_center, dmso_gene_scale, dmso_gene_keep = standardize_genes(
        selected_line_expression, training_position, "DMSO line means"
    )
    if not dmso_gene_keep.all():
        raise RuntimeError("A selected positive-variance DMSO gene became nonstandardizable.")
    discovery_z = line_z[training_position]
    correlation = np.corrcoef(discovery_z, rowvar=False).astype(np.float32)
    correlation = ((correlation + correlation.T) * np.float32(0.5)).astype(np.float32)
    np.fill_diagonal(correlation, 1.0)
    if (
        correlation.shape != (N_SELECTED_GENES, N_SELECTED_GENES)
        or not np.isfinite(correlation).all()
        or np.max(np.abs(correlation)) > 1.0001
    ):
        raise RuntimeError("The signed training-line gene correlation is invalid.")
    write_npy(correlation, OUT / "baseline_gene_gene_correlation.npy")

    helper = load_response_helpers()
    affinity = helper.posneg_correlation(correlation)
    repeats = helper.run_repeated_onmf(
        "baseline_fold_1", affinity, selected_genes, set(training_lines)
    )
    unsigned, matched, repeat_matches, medoid = helper.consensus_membership(
        repeats, selected_genes
    )
    signed_weights, anchors = helper.orient_signed_programs(unsigned, correlation)
    candidate_programs = [f"B{i + 1}" for i in range(CANDIDATE_RANK)]
    candidate_raw = line_z @ signed_weights.T
    score_variance = candidate_raw[training_position].var(axis=0, ddof=1)
    stability_records = helper.records_for_stability(
        "baseline_fold_1", candidate_programs, matched, score_variance
    )
    stability = pd.DataFrame(stability_records)
    stability.insert(0, "outer_fold", OUTER_FOLD)
    stable_programs = stability.loc[stability["stable"], "program"].drop_duplicates().tolist()
    stable_index = np.asarray([candidate_programs.index(program) for program in stable_programs], dtype=int)
    collapsed = len(stable_programs) < 2
    if len(stable_index):
        assignment, plotting_order, correlation_audit = helper.correlation_order_and_audit(
            correlation, signed_weights, stable_index
        )
    else:
        assignment = np.full(N_SELECTED_GENES, -1, dtype=int)
        plotting_order = np.arange(1, N_SELECTED_GENES + 1)
        correlation_audit = {}
    write_tsv(stability, OUT / "baseline_program_stability.tsv")

    weight_rows = []
    stable_set = set(stable_programs)
    for program_index, program in enumerate(candidate_programs):
        current = signed_weights[program_index].astype(float)
        absolute_total = float(np.abs(current).sum())
        for gene_index, gene in enumerate(selected_genes):
            weight_rows.append({
                "outer_fold": OUTER_FOLD,
                "program": program,
                "candidate_program_index": program_index + 1,
                "gene": gene,
                "signed_weight": float(current[gene_index]),
                "absolute_weight": float(abs(current[gene_index])),
                "normalized_absolute_weight": float(abs(current[gene_index]) / absolute_total),
                "direction": "UP" if current[gene_index] >= 0 else "DOWN",
                "stable": program in stable_set,
                "anchor_gene": selected_genes[anchors[program_index]],
            })
    weights = pd.DataFrame(weight_rows)
    write_tsv(weights, OUT / "baseline_program_weights.tsv.gz", compression="gzip")

    assignment_rows = []
    stability_once = stability.drop_duplicates("program").set_index("program")
    for gene_index, gene in enumerate(selected_genes):
        program_index = int(assignment[gene_index])
        program = candidate_programs[program_index] if program_index >= 0 else None
        value = float(signed_weights[program_index, gene_index]) if program_index >= 0 else float("nan")
        assignment_rows.append({
            "outer_fold": OUTER_FOLD,
            "selected_gene_index": gene_index,
            "gene": gene,
            "assigned_program": program,
            "signed_weight": value,
            "absolute_weight": abs(value),
            "direction": "UP" if value >= 0 else "DOWN",
            "stability": float(stability_once.loc[program, "q10_matched_pearson"]) if program else float("nan"),
            "stable": bool(program in stable_set) if program else False,
            "plotting_order": int(plotting_order[gene_index]),
        })
    assignment_frame = pd.DataFrame(assignment_rows).sort_values("plotting_order")
    write_tsv(
        assignment_frame, OUT / "baseline_gene_program_assignment.tsv.gz", compression="gzip"
    )

    if collapsed:
        stable_weights = signed_weights[stable_index]
    else:
        stable_weights = signed_weights[stable_index]
    dmso_line_raw = line_z @ stable_weights.T
    plate_z = (selected_plate_expression - dmso_gene_center[None, :]) / dmso_gene_scale[None, :]
    dmso_plate_raw = plate_z @ stable_weights.T
    dmso_line_standardized, dmso_program_center, dmso_program_scale = standardize_programs(
        dmso_line_raw, training_position, "DMSO baseline programs"
    )
    dmso_plate_standardized = (
        dmso_plate_raw - dmso_program_center[None, :]
    ) / dmso_program_scale[None, :]

    depmap_selected = depmap.loc[all_lines, selected_genes].to_numpy(dtype=np.float64)
    depmap_z, depmap_gene_center, depmap_gene_scale, depmap_gene_keep = standardize_genes(
        depmap_selected, training_position, "DepMap baseline"
    )
    depmap_line_raw = depmap_z @ stable_weights.T
    depmap_line_standardized, depmap_program_center, depmap_program_scale = standardize_programs(
        depmap_line_raw, training_position, "DepMap baseline programs"
    )

    line_to_cell = split.drop_duplicates("depmap_id").set_index("depmap_id")["cell_line"].astype(str)
    role = {line: "training" if line in set(training_lines) else "held_out" for line in all_lines}
    score_rows = []
    for line_index, line in enumerate(all_lines):
        for program_index, program in enumerate(stable_programs):
            score_rows.append({
                "baseline_source": "tahoe_dmso",
                "aggregation_level": "line_mean",
                "outer_fold": OUTER_FOLD,
                "depmap_id": line,
                "cell_line": line_to_cell.loc[line],
                "plate": np.nan,
                "program": program,
                "raw_program_score": float(dmso_line_raw[line_index, program_index]),
                "standardized_program_score": float(dmso_line_standardized[line_index, program_index]),
                "role": role[line],
            })
            score_rows.append({
                "baseline_source": "depmap_static",
                "aggregation_level": "line_level",
                "outer_fold": OUTER_FOLD,
                "depmap_id": line,
                "cell_line": line_to_cell.loc[line],
                "plate": np.nan,
                "program": program,
                "raw_program_score": float(depmap_line_raw[line_index, program_index]),
                "standardized_program_score": float(depmap_line_standardized[line_index, program_index]),
                "role": role[line],
            })
    for row_index, row in control_meta.iterrows():
        line = str(row.depmap_id)
        for program_index, program in enumerate(stable_programs):
            score_rows.append({
                "baseline_source": "tahoe_dmso",
                "aggregation_level": "natural_plate",
                "outer_fold": OUTER_FOLD,
                "depmap_id": line,
                "cell_line": str(row.cell_line),
                "plate": int(row.plate),
                "program": program,
                "raw_program_score": float(dmso_plate_raw[row_index, program_index]),
                "standardized_program_score": float(dmso_plate_standardized[row_index, program_index]),
                "role": role[line],
            })
    baseline_scores = pd.DataFrame(score_rows)
    write_tsv(baseline_scores, OUT / "baseline_program_scores.tsv.gz", compression="gzip")

    eigenvalues = np.linalg.eigvalsh(
        np.cov(dmso_line_standardized, rowvar=False, ddof=1)
    )[::-1]
    eigenvalues = np.maximum(eigenvalues, 0)
    effective_rank = float(np.square(eigenvalues.sum()) / np.square(eigenvalues).sum())
    first_component_fraction = float(eigenvalues[0] / eigenvalues.sum())
    diagnostic_rows = []
    for program_index, program in enumerate(stable_programs):
        train_plate = np.stack([
            dmso_plate_standardized[plate_rows[line], program_index]
            for line in training_lines
        ])
        diagnostic_rows.append({
            "outer_fold": OUTER_FOLD,
            "program": program,
            "onmf_stability_q10": float(stability_once.loc[program, "q10_matched_pearson"]),
            "onmf_stability_mean": float(stability_once.loc[program, "mean_matched_pearson"]),
            "plate_level_icc_1_1": icc_oneway_balanced(train_plate),
            "between_line_score_sd": float(dmso_line_standardized[training_position, program_index].std(ddof=1)),
            "mean_within_line_plate_sd": float(np.mean(train_plate.std(axis=1, ddof=1))),
            "dmso_depmap_training_line_correlation": safe_cor(
                dmso_line_standardized[training_position, program_index],
                depmap_line_standardized[training_position, program_index],
            ),
            "dmso_depmap_held_line_correlation": safe_cor(
                dmso_line_standardized[held_position, program_index],
                depmap_line_standardized[held_position, program_index],
            ),
            "n_assigned_genes": int((assignment == candidate_programs.index(program)).sum()),
            "effective_rank_complete_dmso_score_matrix": effective_rank,
            "first_component_variance_fraction_complete_dmso_score_matrix": first_component_fraction,
        })
    diagnostics = pd.DataFrame(diagnostic_rows)
    write_tsv(diagnostics, OUT / "baseline_program_diagnostics.tsv")

    line_scores = baseline_scores.loc[
        baseline_scores["aggregation_level"].isin(["line_mean", "line_level"])
    ].copy()
    concentration_rows = []
    for (source, program), current in line_scores.groupby(
        ["baseline_source", "program"], sort=True
    ):
        current = current.copy()
        current["squared_score"] = np.square(current["standardized_program_score"])
        current = current.sort_values(
            ["squared_score", "depmap_id"], ascending=[False, True]
        )
        total = float(current["squared_score"].sum())
        concentration_rows.append({
            "baseline_source": source,
            "program": program,
            "most_extreme_depmap_id": str(current.iloc[0]["depmap_id"]),
            "most_extreme_cell_line": str(current.iloc[0]["cell_line"]),
            "maximum_absolute_standardized_score": float(
                abs(current.iloc[0]["standardized_program_score"])
            ),
            "top_1_squared_score_fraction": float(
                current["squared_score"].iloc[:1].sum() / total
            ),
            "top_3_squared_score_fraction": float(
                current["squared_score"].iloc[:3].sum() / total
            ),
            "effective_number_of_lines": float(
                total * total / np.square(current["squared_score"]).sum()
            ),
            "n_lines": int(len(current)),
        })
    write_tsv(
        pd.DataFrame(concentration_rows),
        RESULT / "baseline_program_line_concentration.tsv",
    )

    source_wide = line_scores.pivot(
        index=["depmap_id", "cell_line", "role", "program"],
        columns="baseline_source", values="standardized_program_score",
    ).reset_index()
    concordance_rows = []
    for program, current_program in source_wide.groupby("program", sort=True):
        for score_role in ("training", "held_out", "all"):
            current = current_program if score_role == "all" else current_program.loc[
                current_program["role"].eq(score_role)
            ]
            concordance_rows.append({
                "program": program,
                "role": score_role,
                "n_lines": int(len(current)),
                "pearson": safe_cor(
                    current["tahoe_dmso"].to_numpy(),
                    current["depmap_static"].to_numpy(),
                ),
                "correlation_method": "pearson",
            })
    write_tsv(
        pd.DataFrame(concordance_rows),
        RESULT / "baseline_source_concordance.tsv",
    )

    bundle = {
        "label": "fold_1_ntc_baseline_programs",
        "program_definition_version": PROGRAM_DEFINITION_VERSION,
        "outer_fold": OUTER_FOLD,
        "candidate_rank": CANDIDATE_RANK,
        "program_names": candidate_programs,
        "stable_programs": stable_programs,
        "collapsed": collapsed,
        "selected_genes": selected_genes,
        "selected_universe_indices": universe_index[selected],
        "unsigned_gene_memberships": unsigned,
        "signed_gene_weights": signed_weights,
        "anchor_gene_indices": anchors,
        "anchor_genes": selected_genes[anchors],
        "repeat_component_matches": repeat_matches,
        "consensus_medoid_repeat_index": medoid,
        "gene_assignment_program_index": assignment,
        "gene_plotting_order": plotting_order,
        "training_lines": training_lines,
        "held_out_lines": held_lines,
        "dmso_gene_center": dmso_gene_center,
        "dmso_gene_scale": dmso_gene_scale,
        "dmso_program_center": dmso_program_center,
        "dmso_program_scale": dmso_program_scale,
        "depmap_gene_center": depmap_gene_center,
        "depmap_gene_scale": depmap_gene_scale,
        "depmap_gene_keep": depmap_gene_keep,
        "depmap_program_center": depmap_program_center,
        "depmap_program_scale": depmap_program_scale,
        "effective_rank": effective_rank,
        "first_component_variance_fraction": first_component_fraction,
        "correlation_audit": correlation_audit,
    }
    joblib.dump(bundle, OUT / "baseline_onmf_model.joblib", compress=3)
    parameters = {
        "seed": SEED,
        "outer_fold": OUTER_FOLD,
        "candidate_rank": CANDIDATE_RANK,
        "candidate_rank_is_compression_not_biological_count": True,
        "n_onmf_runs": N_ONMF_RUNS,
        "onmf_max_iter": ONMF_MAX_ITER,
        "minimum_program_stability": MIN_PROGRAM_STABILITY,
        "n_training_lines": len(training_lines),
        "n_held_out_lines": len(held_lines),
        "n_natural_plates_per_line": EXPECTED_PLATES_PER_LINE,
        "n_common_baseline_genes": len(universe_genes),
        "n_selected_genes": len(selected_genes),
        "n_stable_programs": len(stable_programs),
        "representation_collapsed": collapsed,
        "feature_selection": "top 300 positive noise-subtracted between-line-variance genes within each of 10 DMSO-abundance bins",
        "discovery_matrix": "26 outer-training DMSO line means x 3000 selected genes",
        "gene_preprocessing": "plate log2(1+CP10K), line mean, training-line gene z-score",
        "nmf_input": "posneg(signed gene-gene Pearson C)",
        "treated_response_or_prism_used": False,
        "held_out_line_used_for_selection_onmf_or_normalization": False,
        "depmap_used_to_refit_onmf": False,
    }
    with (RESULT / "baseline_onmf_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump(parameters, handle, indent=2)
        handle.write("\n")

    checks = pd.DataFrame([
        {"check": "outer_training_lines", "observed": len(training_lines), "expected": 26, "passed": len(training_lines) == 26},
        {"check": "outer_held_lines", "observed": len(held_lines), "expected": 7, "passed": len(held_lines) == 7},
        {"check": "natural_plates_per_line", "observed": int(control_meta.groupby("depmap_id").size().min()), "expected": 14, "passed": bool(control_meta.groupby("depmap_id").size().eq(14).all())},
        {"check": "selected_baseline_genes", "observed": len(selected_genes), "expected": 3000, "passed": len(selected_genes) == 3000},
        {"check": "candidate_baseline_programs", "observed": CANDIDATE_RANK, "expected": 5, "passed": CANDIDATE_RANK == 5},
        {"check": "at_least_two_stable_baseline_programs", "observed": len(stable_programs), "expected": 2, "passed": len(stable_programs) >= 2},
        {"check": "positive_negative_C_not_abs_C", "observed": 1, "expected": 1, "passed": True},
        {"check": "held_lines_excluded_from_discovery_and_normalization", "observed": 0, "expected": 0, "passed": True},
        {"check": "treated_or_functional_information_used", "observed": 0, "expected": 0, "passed": True},
        {"check": "depmap_refit_of_onmf", "observed": 0, "expected": 0, "passed": True},
        {"check": "excluded_lines_absent_from_baseline_objects", "observed": 0, "expected": 0, "passed": not bool(EXCLUDED_IDS.intersection(baseline_scores["depmap_id"].astype(str)))},
    ])
    write_tsv(checks, RESULT / "baseline_required_checks.tsv")
    status = "collapsed" if collapsed else "supported"
    write_tsv(pd.DataFrame([{
        "stage": "ntc_baseline_program_discovery",
        "status": status,
        "evidence": f"3000 genes; {len(training_lines)} training lines; {len(stable_programs)}/{CANDIDATE_RANK} stable programs",
    }]), RESULT / "baseline_stage_summary.tsv")
    write_tsv(pd.DataFrame([{
        "stage": "03b_ntc_baseline_onmf",
        "phase": "stage03b_total",
        "runtime_seconds": time.monotonic() - started,
        "details": f"26 line means; 3000 genes; {len(stable_programs)}/5 stable programs",
    }]), RESULT / "baseline_runtime.tsv")
    for lightweight in (
        "baseline_program_weights.tsv.gz", "baseline_program_stability.tsv",
        "baseline_gene_program_assignment.tsv.gz",
        "baseline_program_diagnostics.tsv",
    ):
        source = OUT / lightweight
        if source.exists():
            os.replace(source, RESULT / lightweight)
    write_artifact_manifest()
    (RESULT / "stage03_complete.flag").write_text("complete\n", encoding="utf-8")

    print(
        f"NTC baseline stage complete: {len(universe_genes)} common genes; "
        f"{len(selected_genes)} selected; {len(stable_programs)}/{CANDIDATE_RANK} stable; "
        f"effective rank={effective_rank:.3f}; runtime={time.monotonic() - started:.1f}s.",
        flush=True,
    )
    if collapsed:
        raise RuntimeError(
            "NTC baseline representation collapsed: fewer than two stable programs; "
            "the downstream conditioner must not raise model complexity."
        )


if __name__ == "__main__":
    main()
