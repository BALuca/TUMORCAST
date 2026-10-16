#!/usr/bin/env python3
"""Independent structural checks for the stage-30 BRAF--ERBB experiment."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pymupdf
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "30_braf_erbb_programs"
checks: list[tuple[str, bool, str]] = []


def record(name: str, passed: bool, detail: str) -> None:
    checks.append((name, bool(passed), detail))


membership = pd.read_csv(OUT / "frozen_program_genes.tsv", sep="\t")
record("gene_count", len(membership) == 60, f"observed={len(membership)}")
sizes = membership.groupby("program").size()
record("program_structure", len(sizes) == 6 and sizes.min() >= 6,
       f"sizes={','.join(map(str, sorted(sizes.tolist())))}")
role_programs = membership[["program", "role"]].drop_duplicates().groupby("role").size()
record("resolved_y", set(role_programs.index) == {"B", "E", "R"} and (role_programs == 2).all(),
       role_programs.to_dict().__repr__())
generic = re.compile(r"^(RPS|RPL|MRPS|MRPL|MT-|HIST|EIF|EEF|LOC|LINC|MIR|SNOR|SNAR|SCARNA|RNU|RNA5|RNA18|RNA28)")
bad = [g for g in membership.gene if generic.search(str(g))]
record("generic_genes_removed", not bad, f"flagged={bad}")

branch = pd.read_csv(OUT / "branch_prism_response_associations.tsv", sep="\t")
record("branch_prism_association", len(branch) == 3 and (branch.nominal_p < 0.05).all(),
       "; ".join(f"{r.role}: rho={r.rho:.3f}, p={r.nominal_p:.3g}" for r in branch.itertuples()))
coh = pd.read_csv(OUT / "program_context_coherence.tsv", sep="\t")
two = coh[coh.source.isin(["Tahoe_balanced", "CytoSig"])]
record("cross_source_coherence", len(two) == 12 and (two.signed_mean_rho > 0).all(),
       f"range={two.signed_mean_rho.min():.3f}..{two.signed_mean_rho.max():.3f}")

with gzip.open(OUT / "joint_gene_correlations.tsv.gz", "rt") as handle:
    joint = pd.read_csv(handle, sep="\t")
record("square_joint_matrix", joint.shape == (60, 61), f"table_shape={joint.shape}")
coverage = pd.read_csv(OUT / "external_gene_coverage.tsv", sep="\t")
coverage_summary = coverage.groupby("dataset").measured.mean()
record("external_gene_coverage", (coverage_summary >= 0.95).all(), coverage_summary.to_dict().__repr__())

main_pdfs = sorted((OUT / "figures" / "main").glob("*.pdf"))
main_numbers = [int(p.name[:2]) for p in main_pdfs]
record("main_figure_sequence", main_numbers == list(range(1, 11)),
       f"observed={main_numbers}")
all_pdfs = sorted((OUT / "figures").rglob("*.pdf"))
one_page = True
pdf_details = []
for path in all_pdfs:
    with pymupdf.open(path) as document:
        pages = document.page_count
    one_page &= pages == 1
    pdf_details.append(f"{path.name}:{pages}")
record("single_figure_per_pdf", one_page, ", ".join(pdf_details))

for number in (2, 3, 4):
    path = main_pdfs[number - 1]
    with pymupdf.open(path) as document:
        rect = document[0].rect
        width, height = rect.width, rect.height
    record(f"square_correlation_figure_{number:02d}", abs(width - height) < 1.0,
           f"points={width:.1f}x{height:.1f}")

record("output_hygiene", not list(OUT.glob("*.pdf")) and not (OUT / "grant_figures").exists(),
       "figures are confined to figures/main and figures/supplementary")
record("viewer_and_story", (OUT / "index.html").exists() and (OUT / "STORY.md").exists(),
       "index.html and STORY.md")

result = pd.DataFrame(checks, columns=["check", "passed", "detail"])
result.to_csv(OUT / "independent_verification.tsv", sep="\t", index=False)
failed = result.loc[~result.passed]
if not failed.empty:
    raise SystemExit("Verification failed:\n" + failed.to_string(index=False))
print(f"All {len(result)} independent checks passed")
