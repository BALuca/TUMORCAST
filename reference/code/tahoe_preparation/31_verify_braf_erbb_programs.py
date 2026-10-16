#!/usr/bin/env python3
"""Independent structural and result checks for the stage-31 BRAF--ERBB experiment."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pymupdf
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "31_braf_erbb_programs"
checks: list[tuple[str, bool, str]] = []


def record(name: str, passed: bool, detail: str) -> None:
    checks.append((name, bool(passed), detail))


membership = pd.read_csv(OUT / "frozen_program_genes.tsv", sep="\t")
record("gene_count", len(membership) == 88, f"observed={len(membership)}")
sizes = membership.groupby("program").size()
record("program_structure", len(sizes) == 6 and sizes.min() >= 10,
       f"sizes={','.join(map(str, sorted(sizes.tolist())))}")
role_programs = membership[["program", "role"]].drop_duplicates().groupby("role").size()
record("resolved_y", set(role_programs.index) == {"B", "E", "R"} and (role_programs == 2).all(),
       role_programs.to_dict().__repr__())
generic = re.compile(r"^(RPS|RPL|MRPS|MRPL|MT-|HIST|EIF|EEF|LOC|LINC|MIR|SNOR|SNAR|SCARNA|RNU|RNA5|RNA18|RNA28)")
bad = [g for g in membership.gene if generic.search(str(g))]
record("generic_genes_removed", not bad, f"flagged={bad}")

branch = pd.read_csv(OUT / "branch_prism_response_associations.tsv", sep="\t")
record("branch_prism_association", len(branch) == 3 and (branch.pearson.abs() >= 0.4).all() and (branch.pearson_p < 0.05).all(),
       "; ".join(f"{r.role}: Pearson={r.pearson:.3f}, p={r.pearson_p:.3g}" for r in branch.itertuples()))
held = pd.read_csv(OUT / "root_heldout_model_comparison.tsv", sep="\t").iloc[0]
record("heldout_prism_association", abs(held.pearson_controlling_drug) >= 0.4,
       f"Pearson={held.pearson_controlling_drug:.3f}")
resistance = pd.read_csv(OUT / "gse68453_plx_resistance_test.tsv", sep="\t").iloc[0]
adaptive = pd.read_csv(OUT / "gse68453_adaptive_output_benefit_test.tsv", sep="\t").iloc[0]
functional = pd.read_csv(OUT / "gse68453_bypass_test.tsv", sep="\t")
record("functional_chain", abs(resistance.pearson) >= 0.4 and abs(adaptive.pearson) >= 0.4 and
       len(functional) == 2 and (functional.pearson.abs() >= 0.4).all(),
       f"resistance={resistance.pearson:.3f}, 12-line benefit={adaptive.pearson:.3f}, " +
       ", ".join(f"{r.role} interaction={r.pearson:.3f}" for r in functional.itertuples()))

phospho = pd.read_csv(OUT / "gerosa_panorama_prediction_tests.tsv", sep="\t")
erk = phospho[phospho.analyte.str.startswith("ERK")]
sos = phospho[phospho.analyte.str.startswith("SOS1")]
record("phosphoproteomic_bridge", len(erk) == 2 and len(sos) == 1 and
       (erk.one_uM_vem_vs_DMSO < 0.5).all() and
       (erk.EGF_rescue_over_one_uM_vem > 2).all() and
       (sos.one_uM_vem_vs_DMSO < 0.5).all(),
       "; ".join(f"{r.analyte}: vem/vehicle={r.one_uM_vem_vs_DMSO:.3f}, "
                 f"EGF/vem={r.EGF_rescue_over_one_uM_vem:.3f}" for r in phospho.itertuples()))

claim_chain = pd.read_csv(OUT / "grant_claim_evidence_chain.tsv", sep="\t")
record("grant_claim_chain", len(claim_chain) == 6 and claim_chain.step.tolist() == list(range(1, 7)) and
       claim_chain.result.notna().all() and claim_chain.primary_link.notna().all(),
       "six claims connect cross-source recovery, Y hierarchy, molecular, phospho, and functional evidence")

locked_bliss = {"A375": 5.7, "SKMEL24": 29.3, "HT144": 16.3, "DM331": 27.3,
                "SKMEL28": 49.7, "DM13": 39.3, "VMM5A": 43.0, "VMM18": 31.7,
                "DM6": -2.7, "VMM15": -8.7, "VMM12": 3.7, "VMM17": -1.0}
joined = pd.read_csv(OUT / "gse68453_expression_bliss_join.tsv", sep="\t")
observed_bliss = joined[["cell_line", "bliss"]].drop_duplicates().set_index("cell_line")["bliss"].to_dict()
record("bliss_identity_mapping", observed_bliss == locked_bliss,
       "digitized outcomes are joined by explicit cell-line identifier")
coh = pd.read_csv(OUT / "program_context_coherence.tsv", sep="\t")
two = coh[coh.source.isin(["Tahoe_balanced", "CytoSig"])]
record("cross_source_coherence", len(two) == 12 and (two.signed_mean_rho > 0).all(),
       f"range={two.signed_mean_rho.min():.3f}..{two.signed_mean_rho.max():.3f}")

with gzip.open(OUT / "joint_gene_correlations.tsv.gz", "rt") as handle:
    joint = pd.read_csv(handle, sep="\t")
record("square_joint_matrix", joint.shape == (88, 89), f"table_shape={joint.shape}")
coverage = pd.read_csv(OUT / "external_gene_coverage.tsv", sep="\t")
coverage_summary = coverage.groupby("dataset").measured.mean()
record("external_gene_coverage", (coverage_summary >= 0.95).all(), coverage_summary.to_dict().__repr__())

main_pdfs = sorted((OUT / "figures" / "main").glob("*.pdf"))
main_numbers = [int(p.name[:2]) for p in main_pdfs]
record("main_figure_sequence", main_numbers == list(range(1, 12)),
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

for number in (3, 4, 5):
    path = main_pdfs[number - 1]
    with pymupdf.open(path) as document:
        rect = document[0].rect
        width, height = rect.width, rect.height
    record(f"square_correlation_figure_{number:02d}", abs(width - height) < 1.0,
           f"points={width:.1f}x{height:.1f}")

record("output_hygiene", not list(OUT.glob("*.pdf")) and not (OUT / "grant_figures").exists(),
       "figures are confined to figures/main and figures/supplementary")
record("viewer_and_story", (OUT / "index.html").exists() and (OUT / "STORY.md").exists() and
       (OUT / "PHOSPHOPROTEOMICS_AND_FUNCTION.md").exists(),
       "index.html, STORY.md, and PHOSPHOPROTEOMICS_AND_FUNCTION.md")

result = pd.DataFrame(checks, columns=["check", "passed", "detail"])
result.to_csv(OUT / "independent_verification.tsv", sep="\t", index=False)
failed = result.loc[~result.passed]
if not failed.empty:
    raise SystemExit("Verification failed:\n" + failed.to_string(index=False))
print(f"All {len(result)} independent checks passed")
