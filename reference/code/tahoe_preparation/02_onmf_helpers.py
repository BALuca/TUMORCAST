#!/usr/bin/env python3
"""Small shared signed-oNMF helpers retained from the accepted pipeline."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "20")
os.environ.setdefault("OMP_NUM_THREADS", "20")

import numpy as np
import pandas as pd
from dspin.compute import onmf
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform

SEED = 1
RANK = 10
N_ONMF_RUNS = 10
ONMF_MAX_ITER = 500
MIN_PROGRAM_STABILITY = 0.50
N_RESPONSE_GENES = 3000
PROGRAM_DEFINITION_VERSION = "gene_gene_co_response_posneg_v1"
CACHE = Path(".")


def pearson_matrix(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    left_centered = left - left.mean(axis=1, keepdims=True)
    right_centered = right - right.mean(axis=1, keepdims=True)
    left_norm = np.sqrt(np.sum(left_centered**2, axis=1, keepdims=True))
    right_norm = np.sqrt(np.sum(right_centered**2, axis=1, keepdims=True))
    denominator = left_norm @ right_norm.T
    with np.errstate(divide="ignore", invalid="ignore"):
        answer = (left_centered @ right_centered.T) / denominator
    answer[~np.isfinite(answer)] = -1.0
    return answer


def posneg_correlation(correlation: np.ndarray) -> np.ndarray:
    """Retain signed C as separate nonnegative positive/negative channels."""
    affinity = np.concatenate(
        [np.maximum(correlation, 0.0), np.maximum(-correlation, 0.0)], axis=1
    ).astype(np.float32, copy=False)
    if (
        affinity.shape != (N_RESPONSE_GENES, 2 * N_RESPONSE_GENES)
        or not np.isfinite(affinity).all()
        or np.any(affinity < 0)
    ):
        raise RuntimeError("Positive/negative correlation lift is invalid.")
    return affinity


def onmf_signature(label: str, training_lines: set[str], genes: np.ndarray) -> str:
    signature = hashlib.sha256()
    signature.update(PROGRAM_DEFINITION_VERSION.encode("utf-8"))
    signature.update(b"\0")
    signature.update(label.encode("utf-8"))
    signature.update(b"\0")
    signature.update("\n".join(sorted(training_lines)).encode("utf-8"))
    signature.update(b"\0")
    signature.update("\n".join(genes.tolist()).encode("utf-8"))
    signature.update(f"\0rank={RANK};runs={N_ONMF_RUNS};iter={ONMF_MAX_ITER}".encode())
    return signature.hexdigest()


def run_repeated_onmf(
    label: str,
    affinity: np.ndarray,
    selected_genes: np.ndarray,
    training_lines: set[str],
) -> list[np.ndarray]:
    run_dir = CACHE / f"_onmf_{PROGRAM_DEFINITION_VERSION}_{label}"
    run_dir.mkdir(parents=False, exist_ok=True)
    signature_file = run_dir / "input_signature.sha256"
    expected_signature = onmf_signature(label, training_lines, selected_genes)
    if signature_file.exists():
        if signature_file.read_text(encoding="utf-8").strip() != expected_signature:
            raise RuntimeError(f"Refusing incompatible checkpoint directory: {run_dir}")
    else:
        if any(run_dir.iterdir()):
            raise RuntimeError(f"Unsigned oNMF checkpoint directory is not empty: {run_dir}")
        signature_file.write_text(expected_signature + "\n", encoding="utf-8")
    memberships = []
    for seed in range(SEED, SEED + N_ONMF_RUNS):
        checkpoint = run_dir / f"gene_membership_seed_{seed}.npz"
        if checkpoint.exists():
            saved = np.load(checkpoint, allow_pickle=False)
            gene_membership = saved["gene_membership"]
            feature_membership = saved["feature_membership"]
            if (
                gene_membership.shape != (RANK, N_RESPONSE_GENES)
                or feature_membership.shape != (RANK, 2 * N_RESPONSE_GENES)
                or not np.isfinite(gene_membership).all()
                or not np.isfinite(feature_membership).all()
            ):
                raise RuntimeError(f"Invalid oNMF checkpoint: {checkpoint}")
            print(f"{label}: reusing validated co-response oNMF seed {seed}.", flush=True)
        else:
            print(
                f"{label}: co-response pos/neg oNMF seed {seed} "
                f"({seed - SEED + 1}/{N_ONMF_RUNS}).",
                flush=True,
            )
            np.random.seed(seed)
            feature_by_component, observation_by_component = onmf(
                affinity, RANK, max_iter=ONMF_MAX_ITER
            )
            gene_membership = np.asarray(observation_by_component.T, dtype=np.float32)
            feature_membership = np.asarray(feature_by_component.T, dtype=np.float32)
            if (
                gene_membership.shape != (RANK, N_RESPONSE_GENES)
                or feature_membership.shape != (RANK, 2 * N_RESPONSE_GENES)
                or not np.isfinite(gene_membership).all()
                or not np.isfinite(feature_membership).all()
                or np.any(gene_membership < 0)
                or np.any(feature_membership < 0)
            ):
                raise RuntimeError(f"{label}: oNMF seed {seed} returned invalid factors.")
            temporary = checkpoint.with_name(checkpoint.name + ".partial")
            with temporary.open("wb") as handle:
                np.savez_compressed(
                    handle,
                    gene_membership=gene_membership,
                    feature_membership=feature_membership,
                )
            os.replace(temporary, checkpoint)
        membership_norm = np.linalg.norm(gene_membership, axis=1, keepdims=True)
        if np.any(membership_norm <= 0):
            raise RuntimeError(f"{label}: oNMF seed {seed} has an empty component.")
        memberships.append(gene_membership / membership_norm)
    return memberships


def consensus_membership(
    repeat_memberships: list[np.ndarray], selected_genes: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    n_runs = len(repeat_memberships)
    agreement = np.zeros((n_runs, n_runs), dtype=float)
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            correlation = pearson_matrix(repeat_memberships[i], repeat_memberships[j])
            left, right = linear_sum_assignment(-correlation)
            if not np.array_equal(left, np.arange(RANK)):
                raise RuntimeError("Repeat-to-repeat oNMF alignment was not one-to-one.")
            agreement[i, j] = agreement[j, i] = float(correlation[left, right].mean())
    medoid = int(np.argmax(agreement.mean(axis=1)))
    reference = repeat_memberships[medoid]
    aligned = []
    for membership in repeat_memberships:
        correlation = pearson_matrix(reference, membership)
        left, right = linear_sum_assignment(-correlation)
        if not np.array_equal(left, np.arange(RANK)):
            raise RuntimeError("Medoid oNMF alignment was not one-to-one.")
        aligned.append(membership[right])
    consensus = np.mean(np.stack(aligned), axis=0)
    consensus /= np.linalg.norm(consensus, axis=1, keepdims=True)
    realigned = []
    for membership in repeat_memberships:
        correlation = pearson_matrix(consensus, membership)
        left, right = linear_sum_assignment(-correlation)
        if not np.array_equal(left, np.arange(RANK)):
            raise RuntimeError("Consensus oNMF alignment was not one-to-one.")
        realigned.append(membership[right])
    consensus = np.mean(np.stack(realigned), axis=0)
    consensus /= np.linalg.norm(consensus, axis=1, keepdims=True)
    anchor = np.argmax(consensus, axis=1)
    consensus = consensus[np.argsort(selected_genes[anchor], kind="stable")]
    matched = np.empty((n_runs, RANK), dtype=float)
    final_matches = []
    for run_index, membership in enumerate(repeat_memberships):
        correlation = pearson_matrix(consensus, membership)
        left, right = linear_sum_assignment(-correlation)
        if not np.array_equal(left, np.arange(RANK)):
            raise RuntimeError("Final oNMF repeat alignment was not one-to-one.")
        matched[run_index] = correlation[left, right]
        final_matches.append(right)
    return consensus, matched, np.asarray(final_matches), medoid


def orient_signed_programs(
    unsigned_membership: np.ndarray, correlation: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    anchors = np.argmax(unsigned_membership, axis=1)
    signed = np.empty_like(unsigned_membership, dtype=np.float32)
    for program_index, anchor in enumerate(anchors):
        directions = np.where(correlation[:, anchor] < 0, -1.0, 1.0)
        current = unsigned_membership[program_index] * directions
        norm = float(np.linalg.norm(current))
        if not np.isfinite(norm) or norm <= 0:
            raise RuntimeError("A signed gene program has invalid norm.")
        current /= norm
        if current[anchor] < 0:
            current *= -1
        signed[program_index] = current.astype(np.float32)
    if not np.allclose(np.linalg.norm(signed, axis=1), 1.0, atol=2e-6):
        raise RuntimeError("Signed gene programs are not unit L2 normalized.")
    return signed, anchors


def records_for_stability(
    label: str, programs: list[str], matched: np.ndarray, score_variance: np.ndarray
) -> list[dict]:
    records = []
    for j, program in enumerate(programs):
        values = matched[:, j]
        mean_value = float(np.mean(values))
        q10_value = float(np.quantile(values, 0.1))
        minimum_value = float(np.min(values))
        stable = bool(
            q10_value >= MIN_PROGRAM_STABILITY
            and np.isfinite(score_variance[j])
            and score_variance[j] > 0
        )
        for repeat, seed in enumerate(range(SEED, SEED + N_ONMF_RUNS), start=1):
            records.append({
                "basis": label,
                "program": program,
                "repeat": repeat,
                "seed": seed,
                "matched_pearson": float(values[repeat - 1]),
                "n_runs": N_ONMF_RUNS,
                "mean_matched_pearson": mean_value,
                "q10_matched_pearson": q10_value,
                "minimum_matched_pearson": minimum_value,
                "training_score_variance": float(score_variance[j]),
                "stability_threshold": MIN_PROGRAM_STABILITY,
                "stable": stable,
            })
    return records


def correlation_order_and_audit(
    correlation: np.ndarray,
    signed_weights: np.ndarray,
    stable_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict]:
    if len(stable_indices) == 0:
        raise RuntimeError("No stable program is available for gene assignment.")
    assigned = stable_indices[np.argmax(np.abs(signed_weights[stable_indices]), axis=0)]
    plotting_order = np.empty(N_RESPONSE_GENES, dtype=int)
    cursor = 1
    for program_index in stable_indices:
        members = np.flatnonzero(assigned == program_index)
        if len(members) > 1:
            distance = 1.0 - np.abs(correlation[np.ix_(members, members)]).astype(np.float64)
            distance = np.clip((distance + distance.T) * 0.5, 0.0, None)
            np.fill_diagonal(distance, 0.0)
            tree = linkage(squareform(distance, checks=False), method="average", optimal_ordering=True)
            members = members[leaves_list(tree)]
        plotting_order[members] = np.arange(cursor, cursor + len(members))
        cursor += len(members)
    if cursor != N_RESPONSE_GENES + 1 or np.any(plotting_order <= 0):
        raise RuntimeError("Program-block gene ordering is incomplete.")
    upper_i, upper_j = np.triu_indices(N_RESPONSE_GENES, k=1)
    same = assigned[upper_i] == assigned[upper_j]
    values = correlation[upper_i, upper_j].astype(np.float64)
    if not np.any(same) or not np.any(~same):
        raise RuntimeError("Within/between-program correlation audit is undefined.")
    audit = {
        "mean_within_program_signed_correlation": float(values[same].mean()),
        "mean_between_program_signed_correlation": float(values[~same].mean()),
        "mean_within_program_absolute_correlation": float(np.abs(values[same]).mean()),
        "mean_between_program_absolute_correlation": float(np.abs(values[~same]).mean()),
        "median_within_program_absolute_correlation": float(np.median(np.abs(values[same]))),
        "median_between_program_absolute_correlation": float(np.median(np.abs(values[~same]))),
    }
    audit["within_exceeds_between_absolute_correlation"] = bool(
        audit["mean_within_program_absolute_correlation"]
        > audit["mean_between_program_absolute_correlation"]
    )
    return assigned, plotting_order, audit


def bottleneck_audit(
    fold_scores: pd.DataFrame, fold_stability: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    stable = fold_stability.drop_duplicates(["outer_fold", "program"])
    stable = stable.loc[stable["stable"], ["outer_fold", "program"]]
    records = []
    for row in stable.itertuples(index=False):
        current = fold_scores.loc[
            (fold_scores["outer_fold"] == row.outer_fold)
            & (fold_scores["program"] == row.program)
            & (fold_scores["role"] == "training")
        ]
        by_drug = current.groupby("drug", sort=True)["score"].agg(
            mean_score="mean", score_energy=lambda x: float(np.square(x).sum())
        )
        strength = np.abs(by_drug["mean_score"])
        threshold = float(strength.mean() + strength.std())
        energy = by_drug["score_energy"].to_numpy(dtype=float)
        energy_fraction = energy / energy.sum()
        sorted_fraction = np.sort(energy_fraction)[::-1]
        effective_drugs = float(1.0 / np.square(energy_fraction).sum())
        n_strong = int((strength >= threshold).sum())
        barcode_like = bool(
            n_strong < 3 or sorted_fraction[0] > 0.50 or effective_drugs < 5
        )
        records.append({
            "outer_fold": int(row.outer_fold),
            "program": row.program,
            "n_training_drugs": int(len(by_drug)),
            "strong_activation_rule": "absolute drug mean score >= across-drug absolute mean + 1 SD",
            "strong_activation_threshold": threshold,
            "n_strongly_activating_drugs": n_strong,
            "top_drug_energy_fraction": float(sorted_fraction[0]),
            "top_five_drug_energy_fraction": float(sorted_fraction[:5].sum()),
            "energy_hhi": float(np.square(energy_fraction).sum()),
            "effective_drug_count": effective_drugs,
            "barcode_like": barcode_like,
        })
    audit = pd.DataFrame(records)
    barcode_fraction = float(audit["barcode_like"].mean()) if len(audit) else np.nan
    median_effective = float(audit["effective_drug_count"].median()) if len(audit) else np.nan
    passed = bool(len(audit) and barcode_fraction <= 0.20 and median_effective >= 10)
    gate = pd.DataFrame([{
        "gate": "response_program_bottleneck",
        "n_stable_fold_programs": len(audit),
        "barcode_like_fraction": barcode_fraction,
        "maximum_allowed_barcode_like_fraction": 0.20,
        "median_effective_drug_count": median_effective,
        "minimum_required_median_effective_drug_count": 10.0,
        "passed": passed,
        "rule": "barcode-like if <3 strong drugs, top-drug energy >0.50, or effective drug count <5",
    }])
    return audit, gate
