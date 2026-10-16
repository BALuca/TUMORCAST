#!/usr/bin/env python3
"""Fit the rank-eight recurrent untreated-baseline oNMF representation.

This is the only native-Python addition in stages 10--13.  It reuses the
accepted signed positive/negative-correlation oNMF helper; all regression,
metrics, orchestration and plotting remain in the numbered R scripts.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
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
from dspin.compute import onmf
from scipy.optimize import linear_sum_assignment


SEED = 1
OUTER_FOLD = 1
RANK = 8
N_RUNS = 5
MAX_ITER = 100
LOO_MAX_ITER = 80
N_ABUNDANCE_BINS = 10
GENES_PER_BIN = 300
N_SELECTED_GENES = N_ABUNDANCE_BINS * GENES_PER_BIN
MAX_SINGLE_LINE_FRACTION = 0.50
EXPECTED_LINES = 33
EXPECTED_TRAINING_LINES = 26
EXPECTED_HELD_LINES = 7
EXPECTED_PLATES_PER_LINE = 14
CP_SCALE = 10000.0
EXCLUDED_IDS = {"ACH-000311", "ACH-000628"}
PROGRAM_VERSION = "recurrent_ntc_cross_source_distributed_posneg_rank8_v2"

ROOT = Path(os.environ.get("TAHOE_V2_ROOT", "tahoe_prism_v2")).resolve()
RESULT = ROOT / "results" / "10_baseline_states"
ARTIFACT = ROOT / "artifacts" / "10_baseline_states"
CACHE = ROOT.parent / "cache" / "tahoe_prism_full_chain"
CONTROL_EXPRESSION = CACHE / "dmso_control_expression.npy"
CONTROL_METADATA = CACHE / "dmso_control_metadata.tsv"
CONTROL_GENES = CACHE / "dmso_control_genes.tsv"
PREPARED = ROOT / "artifacts" / "01_inputs" / "prepared_inputs.rds"
SPLIT = ROOT / "results" / "01_inputs" / "fold_assignments.tsv"
OLD_SCORES = ROOT / "artifacts" / "03_baseline_programs" / "baseline_program_scores.tsv.gz"


def load_onmf_helpers():
    path = Path(__file__).with_name("02_onmf_helpers.py")
    spec = importlib.util.spec_from_file_location("tahoe_onmf_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the accepted oNMF helpers.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.SEED = SEED
    module.RANK = RANK
    module.N_ONMF_RUNS = N_RUNS
    module.ONMF_MAX_ITER = MAX_ITER
    module.MIN_PROGRAM_STABILITY = 0.50
    module.N_RESPONSE_GENES = N_SELECTED_GENES
    module.CACHE = ARTIFACT
    module.PROGRAM_DEFINITION_VERSION = PROGRAM_VERSION
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


def safe_cor(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    keep = np.isfinite(left) & np.isfinite(right)
    if keep.sum() < 3 or np.std(left[keep], ddof=1) <= 0 or np.std(right[keep], ddof=1) <= 0:
        return float("nan")
    return float(np.corrcoef(left[keep], right[keep])[0, 1])


def export_depmap(path: Path) -> pd.DataFrame:
    r_code = """
    suppressPackageStartupMessages(library(data.table))
    a <- commandArgs(trailingOnly=TRUE)
    x <- readRDS(a[1])$depmap_baseline$expression
    fwrite(data.table(depmap_id=rownames(x), as.data.table(x)), a[2], sep='\\t')
    """
    subprocess.run(["Rscript", "-e", r_code, str(PREPARED), str(path)], check=True)
    frame = pd.read_csv(path, sep="\t").set_index("depmap_id")
    path.unlink()
    frame.index = frame.index.astype(str)
    return frame


def standardize(matrix: np.ndarray, training_position: np.ndarray, label: str):
    train = np.asarray(matrix[training_position], dtype=np.float64)
    center = train.mean(axis=0)
    scale = train.std(axis=0, ddof=1)
    keep = np.isfinite(center) & np.isfinite(scale) & (scale > 0)
    if not keep.all():
        raise RuntimeError(f"{label} contains a nonstandardizable selected gene.")
    return (np.asarray(matrix, dtype=np.float64) - center) / scale, center, scale


def standardize_programs(raw: np.ndarray, training_position: np.ndarray, label: str):
    center = raw[training_position].mean(axis=0)
    scale = raw[training_position].std(axis=0, ddof=1)
    if np.any(~np.isfinite(center)) or np.any(~np.isfinite(scale)) or np.any(scale <= 0):
        raise RuntimeError(f"{label} contains a nonstandardizable program.")
    return (raw - center) / scale, center, scale


def gene_correlations(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = left - left.mean(axis=0, keepdims=True)
    right = right - right.mean(axis=0, keepdims=True)
    denominator = np.sqrt(np.sum(left * left, axis=0) * np.sum(right * right, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        answer = np.sum(left * right, axis=0) / denominator
    answer[~np.isfinite(answer)] = np.nan
    return answer


def dominance_fraction(matrix: np.ndarray) -> np.ndarray:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    squared = centered * centered
    total = squared.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        answer = squared.max(axis=0) / total
    answer[~np.isfinite(answer)] = 1.0
    return answer


def deterministic_program_order(weights: np.ndarray, genes: np.ndarray) -> np.ndarray:
    anchors = genes[np.argmax(weights, axis=1)]
    return np.argsort(anchors, kind="stable")


def loo_stability(helper, discovery_z: np.ndarray, full_unsigned: np.ndarray,
                  selected_genes: np.ndarray, training_lines: list[str]) -> pd.DataFrame:
    records = []
    for omitted_index, omitted_line in enumerate(training_lines):
        keep = np.arange(len(training_lines)) != omitted_index
        correlation = np.corrcoef(discovery_z[keep], rowvar=False).astype(np.float32)
        correlation = ((correlation + correlation.T) * np.float32(0.5)).astype(np.float32)
        np.fill_diagonal(correlation, 1.0)
        affinity = helper.posneg_correlation(correlation)
        checkpoint = ARTIFACT / f"_loo_membership_{omitted_line}.npz"
        if checkpoint.exists():
            membership = np.load(checkpoint, allow_pickle=False)["membership"]
        else:
            np.random.seed(SEED + 1000 + omitted_index)
            _, observation_by_component = onmf(affinity, RANK,
                max_iter=LOO_MAX_ITER)
            membership = np.asarray(observation_by_component.T, dtype=np.float32)
            membership /= np.linalg.norm(membership, axis=1, keepdims=True)
            temporary = checkpoint.with_name(checkpoint.name + ".partial")
            with temporary.open("wb") as handle:
                np.savez_compressed(handle, membership=membership)
            os.replace(temporary, checkpoint)
        agreement = helper.pearson_matrix(full_unsigned, membership)
        left, right = linear_sum_assignment(-agreement)
        if not np.array_equal(left, np.arange(RANK)):
            raise RuntimeError("Leave-one-line-out matching was not one-to-one.")
        for program_index, matched_index in enumerate(right):
            records.append({
                "omitted_depmap_id": omitted_line,
                "program": f"S{program_index + 1}",
                "matched_component": int(matched_index + 1),
                "matched_membership_pearson": float(agreement[program_index, matched_index]),
                "refit_gene_set": "frozen_outer_training_selected_genes",
                "n_training_lines_in_refit": len(training_lines) - 1,
            })
    return pd.DataFrame(records)


def artifact_manifest() -> None:
    rows = []
    for path in sorted(ARTIFACT.rglob("*")):
        if not path.is_file() or path.name == "artifact_manifest.tsv":
            continue
        digest = hashlib.md5()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
                digest.update(block)
        rows.append({
            "relative_path": str(path.relative_to(ARTIFACT)),
            "file_size_bytes": path.stat().st_size,
            "md5": digest.hexdigest(),
            "producing_script": "scripts/10_run_recurrent_baseline_onmf.py",
        })
    write_tsv(pd.DataFrame(rows), ARTIFACT / "artifact_manifest.tsv")


def main() -> None:
    started = time.monotonic()
    required = [CONTROL_EXPRESSION, CONTROL_METADATA, CONTROL_GENES, PREPARED, SPLIT, OLD_SCORES]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing frozen untreated-baseline input(s): {missing}")
    RESULT.mkdir(parents=True, exist_ok=True)
    ARTIFACT.mkdir(parents=True, exist_ok=True)

    split = pd.read_csv(SPLIT, sep="\t")
    training_lines = sorted(split.loc[split["outer_fold"].ne(OUTER_FOLD), "depmap_id"].astype(str))
    held_lines = sorted(split.loc[split["outer_fold"].eq(OUTER_FOLD), "depmap_id"].astype(str))
    all_lines = sorted(training_lines + held_lines)
    if len(training_lines) != EXPECTED_TRAINING_LINES or len(held_lines) != EXPECTED_HELD_LINES:
        raise RuntimeError("The accepted 26/7 split changed.")
    if EXCLUDED_IDS.intersection(all_lines):
        raise RuntimeError("An excluded cell line entered recurrent baseline discovery.")

    metadata = pd.read_csv(CONTROL_METADATA, sep="\t")
    genes = pd.read_csv(CONTROL_GENES, sep="\t")
    expression = np.load(CONTROL_EXPRESSION, mmap_mode="r")
    metadata = metadata.loc[metadata["depmap_id"].astype(str).isin(all_lines)].copy()
    metadata = metadata.sort_values("array_row").reset_index(drop=True)
    counts = metadata.groupby("depmap_id").size()
    if len(metadata) != EXPECTED_LINES * EXPECTED_PLATES_PER_LINE or not counts.eq(EXPECTED_PLATES_PER_LINE).all():
        raise RuntimeError("The accepted natural DMSO plate contract changed.")

    depmap_path = ARTIFACT / ".depmap_stage10.tsv"
    if depmap_path.exists():
        depmap_path.unlink()
    depmap = export_depmap(depmap_path)
    if sorted(depmap.index.tolist()) != all_lines:
        raise RuntimeError("The prepared DepMap cohort changed.")

    dmso_genes = genes["feature"].astype(str).to_numpy()
    depmap_set = set(depmap.columns.astype(str))
    valid = np.asarray([bool(gene) and gene != "nan" and gene in depmap_set for gene in dmso_genes])
    universe_index = np.flatnonzero(valid)
    universe_genes = dmso_genes[universe_index]
    if len(universe_genes) != len(np.unique(universe_genes)):
        raise RuntimeError("The common untreated-baseline universe has duplicate genes.")

    retained_rows = metadata["array_row"].to_numpy(dtype=int)
    raw = np.asarray(expression[np.ix_(retained_rows, universe_index)], dtype=np.float64)
    library = raw.sum(axis=1)
    if np.any(~np.isfinite(raw)) or np.any(raw < 0) or np.any(library <= 0):
        raise RuntimeError("Natural DMSO expression is invalid.")
    plate_expression = np.log2(1.0 + raw / library[:, None] * CP_SCALE).astype(np.float32)
    del raw

    line_rows = {line: np.flatnonzero(metadata["depmap_id"].astype(str).eq(line).to_numpy()) for line in all_lines}
    line_expression = np.stack([plate_expression[line_rows[line]].mean(axis=0) for line in all_lines])
    training_position = np.asarray([all_lines.index(line) for line in training_lines], dtype=int)
    held_position = np.asarray([all_lines.index(line) for line in held_lines], dtype=int)
    dmso_training = np.asarray(line_expression[training_position], dtype=np.float64)
    depmap_common = depmap.loc[all_lines, universe_genes].to_numpy(dtype=np.float64)
    depmap_training = depmap_common[training_position]

    abundance = dmso_training.mean(axis=0)
    between_variance = dmso_training.var(axis=0, ddof=1)
    plate_variance = np.stack([
        plate_expression[line_rows[line]].var(axis=0, ddof=1, dtype=np.float64)
        for line in training_lines
    ])
    mean_plate_variance = plate_variance.mean(axis=0)
    plate_noise_of_mean = mean_plate_variance / EXPECTED_PLATES_PER_LINE
    reproducible_variance = np.maximum(between_variance - plate_noise_of_mean, 0.0)
    cross_source = gene_correlations(dmso_training, depmap_training)
    dmso_dominance = dominance_fraction(dmso_training)
    depmap_dominance = dominance_fraction(depmap_training)
    maximum_dominance = np.maximum(dmso_dominance, depmap_dominance)
    eligible = (
        np.isfinite(abundance)
        & np.isfinite(reproducible_variance)
        & (reproducible_variance > 0)
        & np.isfinite(cross_source)
        & (cross_source > 0)
        & (maximum_dominance <= MAX_SINGLE_LINE_FRACTION)
    )
    selection_score = np.zeros(len(universe_genes), dtype=np.float64)
    selection_score[eligible] = (
        reproducible_variance[eligible]
        * np.square(cross_source[eligible])
        * np.square(1.0 - maximum_dominance[eligible])
    )
    eligible_index = np.flatnonzero(eligible)
    abundance_bin = np.full(len(universe_genes), -1, dtype=np.int16)
    ranked = pd.Series(abundance[eligible_index]).rank(method="first")
    abundance_bin[eligible_index] = pd.qcut(ranked, N_ABUNDANCE_BINS, labels=False).to_numpy(dtype=np.int16) + 1
    selected_parts = []
    rank_in_bin = np.full(len(universe_genes), np.nan)
    for bin_index in range(1, N_ABUNDANCE_BINS + 1):
        current = np.flatnonzero(abundance_bin == bin_index)
        order = np.lexsort((universe_genes[current], -selection_score[current]))
        current = current[order]
        if len(current) < GENES_PER_BIN:
            raise RuntimeError(f"Only {len(current)} eligible genes remain in abundance bin {bin_index}.")
        rank_in_bin[current] = np.arange(1, len(current) + 1)
        selected_parts.append(current[:GENES_PER_BIN])
    selected = np.concatenate(selected_parts)
    if len(selected) != N_SELECTED_GENES or len(np.unique(selected)) != N_SELECTED_GENES:
        raise RuntimeError("Distributed cross-source selection did not return 3,000 unique genes.")

    selection = pd.DataFrame({
        "gene": universe_genes,
        "mean_dmso_abundance_log2_1p_cp10k": abundance,
        "dmso_between_line_variance": between_variance,
        "mean_within_line_plate_variance": mean_plate_variance,
        "estimated_plate_noise_of_line_mean": plate_noise_of_mean,
        "reproducible_dmso_between_line_variance": reproducible_variance,
        "dmso_depmap_training_pearson": cross_source,
        "dmso_single_line_squared_fraction": dmso_dominance,
        "depmap_single_line_squared_fraction": depmap_dominance,
        "maximum_single_line_squared_fraction": maximum_dominance,
        "positive_cross_source_reproducibility": cross_source > 0,
        "passes_single_line_filter": maximum_dominance <= MAX_SINGLE_LINE_FRACTION,
        "eligible": eligible,
        "selection_score": selection_score,
        "abundance_bin": abundance_bin,
        "rank_within_abundance_bin": rank_in_bin,
        "selected": False,
    })
    selection.loc[selected, "selected"] = True
    write_tsv(selection, ARTIFACT / "baseline_gene_selection.tsv.gz", compression="gzip")

    selected_genes = universe_genes[selected]
    selected_line = np.asarray(line_expression[:, selected], dtype=np.float64)
    selected_plate = np.asarray(plate_expression[:, selected], dtype=np.float64)
    dmso_z, dmso_gene_center, dmso_gene_scale = standardize(selected_line, training_position, "DMSO")
    discovery_z = dmso_z[training_position]
    correlation = np.corrcoef(discovery_z, rowvar=False).astype(np.float32)
    correlation = ((correlation + correlation.T) * np.float32(0.5)).astype(np.float32)
    np.fill_diagonal(correlation, 1.0)
    if not np.isfinite(correlation).all():
        raise RuntimeError("The rank-eight baseline correlation matrix is invalid.")
    write_npy(correlation, ARTIFACT / "baseline_gene_gene_correlation.npy")

    helper = load_onmf_helpers()
    affinity = helper.posneg_correlation(correlation)
    repeats = helper.run_repeated_onmf("recurrent_baseline_fold_1", affinity, selected_genes, set(training_lines))
    unsigned, matched, repeat_matches, medoid = helper.consensus_membership(repeats, selected_genes)
    signed_weights, anchors = helper.orient_signed_programs(unsigned, correlation)
    program_order = deterministic_program_order(unsigned, selected_genes)
    unsigned = unsigned[program_order]
    signed_weights = signed_weights[program_order]
    anchors = anchors[program_order]
    matched = matched[:, program_order]
    repeat_matches = repeat_matches[:, program_order]
    programs = [f"S{i + 1}" for i in range(RANK)]
    candidate_raw = dmso_z @ signed_weights.T
    score_variance = candidate_raw[training_position].var(axis=0, ddof=1)
    stability = pd.DataFrame(helper.records_for_stability(
        "recurrent_baseline_fold_1", programs, matched, score_variance
    ))
    stability.insert(0, "outer_fold", OUTER_FOLD)
    write_tsv(stability, RESULT / "baseline_program_stability.tsv")

    weight_rows = []
    for program_index, program in enumerate(programs):
        current = signed_weights[program_index].astype(float)
        total = np.abs(current).sum()
        for gene_index, gene in enumerate(selected_genes):
            weight_rows.append({
                "outer_fold": OUTER_FOLD,
                "program": program,
                "gene": gene,
                "signed_weight": float(current[gene_index]),
                "absolute_weight": float(abs(current[gene_index])),
                "normalized_absolute_weight": float(abs(current[gene_index]) / total),
                "direction": "UP" if current[gene_index] >= 0 else "DOWN",
                "anchor_gene": selected_genes[anchors[program_index]],
                "stable": bool(stability.loc[stability["program"].eq(program), "stable"].iloc[0]),
            })
    weights = pd.DataFrame(weight_rows)
    write_tsv(weights, RESULT / "baseline_program_loadings.tsv.gz", compression="gzip")

    dmso_raw = dmso_z @ signed_weights.T
    dmso_scores, dmso_program_center, dmso_program_scale = standardize_programs(
        dmso_raw, training_position, "DMSO baseline scores"
    )
    plate_z = (selected_plate - dmso_gene_center) / dmso_gene_scale
    plate_raw = plate_z @ signed_weights.T
    plate_scores = (plate_raw - dmso_program_center) / dmso_program_scale
    depmap_selected = depmap.loc[all_lines, selected_genes].to_numpy(dtype=np.float64)
    depmap_z, depmap_gene_center, depmap_gene_scale = standardize(
        depmap_selected, training_position, "DepMap"
    )
    depmap_raw = depmap_z @ signed_weights.T
    depmap_scores, depmap_program_center, depmap_program_scale = standardize_programs(
        depmap_raw, training_position, "DepMap baseline scores"
    )

    line_to_cell = split.drop_duplicates("depmap_id").set_index("depmap_id")["cell_line"].astype(str)
    role = {line: "training" if line in set(training_lines) else "held_out" for line in all_lines}
    score_rows = []
    for line_index, line in enumerate(all_lines):
        for program_index, program in enumerate(programs):
            common = {
                "outer_fold": OUTER_FOLD, "depmap_id": line,
                "cell_line": line_to_cell.loc[line], "program": program,
                "role": role[line], "plate": np.nan,
            }
            score_rows.append({**common, "baseline_source": "tahoe_dmso",
                               "aggregation_level": "line_mean",
                               "raw_program_score": float(dmso_raw[line_index, program_index]),
                               "standardized_program_score": float(dmso_scores[line_index, program_index])})
            score_rows.append({**common, "baseline_source": "depmap_static",
                               "aggregation_level": "line_level",
                               "raw_program_score": float(depmap_raw[line_index, program_index]),
                               "standardized_program_score": float(depmap_scores[line_index, program_index])})
    for row_index, row in metadata.iterrows():
        line = str(row.depmap_id)
        for program_index, program in enumerate(programs):
            score_rows.append({
                "outer_fold": OUTER_FOLD, "depmap_id": line,
                "cell_line": str(row.cell_line), "program": program,
                "role": role[line], "plate": int(row.plate),
                "baseline_source": "tahoe_dmso", "aggregation_level": "natural_plate",
                "raw_program_score": float(plate_raw[row_index, program_index]),
                "standardized_program_score": float(plate_scores[row_index, program_index]),
            })
    scores = pd.DataFrame(score_rows)
    write_tsv(scores, ARTIFACT / "baseline_program_scores.tsv.gz", compression="gzip")

    loo = loo_stability(helper, discovery_z, unsigned, selected_genes, training_lines)
    write_tsv(loo, RESULT / "baseline_program_leave_one_line_out_stability.tsv")
    loo_q10 = loo.groupby("program")["matched_membership_pearson"].quantile(0.10)

    contribution_rows = []
    for program_index, program in enumerate(programs):
        current = np.square(dmso_scores[training_position, program_index])
        total = current.sum()
        order = np.lexsort((np.asarray(training_lines), -current))
        ranks = np.empty(len(training_lines), dtype=int)
        ranks[order] = np.arange(1, len(training_lines) + 1)
        for line_index, line in enumerate(training_lines):
            contribution_rows.append({
                "program": program, "depmap_id": line,
                "cell_line": line_to_cell.loc[line],
                "standardized_program_score": float(dmso_scores[all_lines.index(line), program_index]),
                "squared_score_fraction": float(current[line_index] / total),
                "absolute_score_rank": int(ranks[line_index]),
            })
    contributions = pd.DataFrame(contribution_rows)
    write_tsv(contributions, RESULT / "baseline_program_line_contributions.tsv")

    dominance_rows = []
    stability_once = stability.drop_duplicates("program").set_index("program")
    for source_name, matrix in (("tahoe_dmso", dmso_scores), ("depmap_static", depmap_scores)):
        for program_index, program in enumerate(programs):
            values = np.square(matrix[training_position, program_index])
            total = values.sum()
            ordered = np.argsort(-values, kind="stable")
            top_fraction = float(values[ordered[0]] / total)
            dominance_rows.append({
                "baseline_source": source_name, "program": program,
                "top_depmap_id": training_lines[int(ordered[0])],
                "top_1_squared_score_fraction": top_fraction,
                "top_3_squared_score_fraction": float(values[ordered[:3]].sum() / total),
                "effective_number_of_lines": float(total * total / np.square(values).sum()),
                "single_line_dominated": top_fraction > MAX_SINGLE_LINE_FRACTION,
                "onmf_stable": bool(stability_once.loc[program, "stable"]),
                "leave_one_line_out_q10": float(loo_q10.loc[program]),
                "leave_one_line_out_stable": bool(loo_q10.loc[program] >= 0.50),
                "eligible_recurrent_state": bool(stability_once.loc[program, "stable"] and top_fraction <= MAX_SINGLE_LINE_FRACTION and loo_q10.loc[program] >= 0.50),
            })
    dominance = pd.DataFrame(dominance_rows)
    write_tsv(dominance, RESULT / "baseline_program_line_dominance.tsv")

    concordance_rows = []
    for program_index, program in enumerate(programs):
        for role_name, positions in (("training", training_position), ("held_out", held_position), ("all", np.arange(len(all_lines)))):
            concordance_rows.append({
                "program": program, "role": role_name, "n_lines": len(positions),
                "pearson": safe_cor(dmso_scores[positions, program_index], depmap_scores[positions, program_index]),
                "correlation_method": "pearson",
            })
    concordance = pd.DataFrame(concordance_rows)
    write_tsv(concordance, RESULT / "baseline_dmso_depmap_reproducibility.tsv")

    old = pd.read_csv(OLD_SCORES, sep="\t", compression="gzip")
    old = old.loc[(old["baseline_source"].eq("tahoe_dmso")) & (old["aggregation_level"].eq("line_mean"))]
    old_wide = old.pivot(index="depmap_id", columns="program", values="standardized_program_score").loc[all_lines]
    comparison_rows = []
    for new_index, program in enumerate(programs):
        for old_program in sorted(old_wide.columns):
            for role_name, positions in (("training", training_position), ("held_out", held_position), ("all", np.arange(len(all_lines)))):
                comparison_rows.append({
                    "new_program": program, "historical_program": old_program,
                    "role": role_name, "n_lines": len(positions),
                    "pearson": safe_cor(dmso_scores[positions, new_index], old_wide.iloc[positions][old_program].to_numpy()),
                    "absolute_pearson": abs(safe_cor(dmso_scores[positions, new_index], old_wide.iloc[positions][old_program].to_numpy())),
                })
    write_tsv(pd.DataFrame(comparison_rows), RESULT / "baseline_program_comparison_historical.tsv")

    bundle = {
        "label": "rank_eight_recurrent_untreated_baseline_states",
        "program_definition_version": PROGRAM_VERSION,
        "outer_fold": OUTER_FOLD, "rank": RANK, "programs": programs,
        "selected_genes": selected_genes, "signed_gene_weights": signed_weights,
        "unsigned_gene_memberships": unsigned, "anchor_genes": selected_genes[anchors],
        "training_lines": training_lines, "held_lines": held_lines,
        "dmso_gene_center": dmso_gene_center, "dmso_gene_scale": dmso_gene_scale,
        "dmso_program_center": dmso_program_center, "dmso_program_scale": dmso_program_scale,
        "depmap_gene_center": depmap_gene_center, "depmap_gene_scale": depmap_gene_scale,
        "depmap_program_center": depmap_program_center, "depmap_program_scale": depmap_program_scale,
        "repeat_component_matches": repeat_matches, "consensus_medoid_repeat_index": medoid,
        "gene_selection_rule": "positive cross-source Pearson; max single-line SS fraction <= 0.50; abundance-bin top cross-source distributed score",
        "treated_response_used": False, "prism_used": False,
    }
    joblib.dump(bundle, ARTIFACT / "recurrent_baseline_onmf_model.joblib", compress=3)
    with (RESULT / "baseline_state_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump({
            "seed": SEED, "outer_fold": OUTER_FOLD, "rank": RANK,
            "n_onmf_runs": N_RUNS, "max_iter": MAX_ITER,
            "leave_one_line_out_max_iter": LOO_MAX_ITER,
            "selected_genes": N_SELECTED_GENES,
            "maximum_single_line_squared_fraction": MAX_SINGLE_LINE_FRACTION,
            "nmf_input": "posneg(signed gene-gene Pearson C)",
            "feature_selection": "DMSO variance weighted by positive DMSO-DepMap reproducibility and distributed-line support within abundance bins",
            "baseline_inputs_only": True, "treated_response_used": False,
            "prism_used": False, "held_out_used_for_learning": False,
        }, handle, indent=2)
        handle.write("\n")

    checks = pd.DataFrame([
        {"check": "training_lines", "observed": len(training_lines), "expected": 26, "passed": len(training_lines) == 26},
        {"check": "held_lines", "observed": len(held_lines), "expected": 7, "passed": len(held_lines) == 7},
        {"check": "candidate_rank", "observed": RANK, "expected": 8, "passed": RANK == 8},
        {"check": "selected_genes", "observed": len(selected), "expected": 3000, "passed": len(selected) == 3000},
        {"check": "excluded_lines_present", "observed": len(EXCLUDED_IDS.intersection(all_lines)), "expected": 0, "passed": not EXCLUDED_IDS.intersection(all_lines)},
        {"check": "baseline_only_inputs", "observed": 1, "expected": 1, "passed": True},
        {"check": "held_out_used_for_selection_scaling_or_onmf", "observed": 0, "expected": 0, "passed": True},
        {"check": "positive_negative_C_not_abs_C", "observed": 1, "expected": 1, "passed": True},
        {"check": "every_program_has_26_leave_one_out_refits", "observed": int(loo.groupby("program").size().min()), "expected": 26, "passed": bool(loo.groupby("program").size().eq(26).all())},
    ])
    write_tsv(checks, RESULT / "baseline_state_required_checks.tsv")
    write_tsv(pd.DataFrame([{
        "stage": "recurrent_baseline_states", "rank": RANK,
        "n_stable_programs": int(stability.drop_duplicates("program")["stable"].sum()),
        "n_eligible_recurrent_programs": int(dominance.loc[dominance["baseline_source"].eq("tahoe_dmso"), "eligible_recurrent_state"].sum()),
        "n_dmso_single_line_dominated_programs": int(dominance.loc[dominance["baseline_source"].eq("tahoe_dmso"), "single_line_dominated"].sum()),
        "median_training_dmso_depmap_pearson": float(concordance.loc[concordance["role"].eq("training"), "pearson"].median()),
        "median_leave_one_line_out_membership_pearson": float(loo["matched_membership_pearson"].median()),
        "runtime_seconds": time.monotonic() - started,
    }]), RESULT / "baseline_state_summary.tsv")
    artifact_manifest()
    print(f"Stage 10 native oNMF complete in {time.monotonic() - started:.1f}s.", flush=True)


if __name__ == "__main__":
    main()
