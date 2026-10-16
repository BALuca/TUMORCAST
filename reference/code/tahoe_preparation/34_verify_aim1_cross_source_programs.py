#!/usr/bin/env python3
from pathlib import Path
import re
import pandas as pd
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "34_aim1_cross_source_programs"
checks = []
def check(name, ok, detail): checks.append((name, bool(ok), detail))

m = pd.read_csv(OUT / "frozen_program_genes.tsv", sep="\t")
sizes = m.groupby("program").size().to_dict()
check("three_sizeable_programs", sizes == {"P1": 20, "P2": 20, "P3": 14}, str(sizes))
bad = [g for g in m.gene if re.match(r"^(RPS|RPL|MRPS|MRPL|MT-|HIST|EIF|EEF|LOC|LINC|MIR|SNOR|RNU)", str(g))]
check("no_generic_genes", not bad, str(bad))

audit = pd.read_csv(OUT / "correlation_and_heatmap_sample_audit.tsv", sep="\t").set_index("source")
check("all_tahoe_profiles", audit.loc["Tahoe", "profiles"] == 252, str(audit.loc["Tahoe"].to_dict()))
check("all_cytosig_profiles", audit.loc["CytoSig", "profiles"] == 2056, str(audit.loc["CytoSig"].to_dict()))
method = pd.read_csv(OUT / "correlation_method_audit.tsv", sep="\t")
check("individual_profile_correlations", len(method) == 5 and not method.treatment_medians_used.astype(bool).any() and
      method.loc[method.source.eq("CytoSig"), "individual_profiles"].iloc[0] == 2056,
      "; ".join(f"{r.source}: {r.correlation}" for r in method.itertuples()))

coh = pd.read_csv(OUT / "program_source_coherence.tsv", sep="\t")
check("cross_source_coherence", len(coh) == 6 and (coh.mean_signed_rho >= .10).all(),
      f"range={coh.mean_signed_rho.min():.3f}..{coh.mean_signed_rho.max():.3f}")
train = pd.read_csv(OUT / "matched_program_training_viability.tsv", sep="\t")
check("matched_training_function", len(train) == 3 and (train.pearson >= .4).all() and (train.p < .05).all(),
      "; ".join(f"{r.program}/{r.drug}: r={r.pearson:.3f}" for r in train.itertuples()))
held = pd.read_csv(OUT / "heldout_prism_test.tsv", sep="\t").iloc[0]
check("heldout_function", held.pearson >= .4 and held.p < .05 and held.program_rmse < held.baseline_rmse,
      f"r={held.pearson:.3f}; RMSE={held.baseline_rmse:.3f}->{held.program_rmse:.3f}")

pdfs = sorted((OUT / "figures" / "main").glob("*.pdf"))
check("main_figure_sequence", [int(p.name[:2]) for p in pdfs] == list(range(1, 12)), str([p.name for p in pdfs]))
single = all(pymupdf.open(p).page_count == 1 for p in pdfs)
check("one_figure_per_pdf", single, f"n={len(pdfs)}")
for n in (4, 5, 6):
    p = [x for x in pdfs if int(x.name[:2]) == n][0]
    with pymupdf.open(p) as d: box = d[0].rect
    check(f"square_correlation_{n:02d}", abs(box.width-box.height) < 1, f"{box.width:.1f}x{box.height:.1f}")
check("viewer", (OUT / "index.html").exists() and (OUT / "README.md").exists(), "viewer and summary present")

result = pd.DataFrame(checks, columns=["check", "passed", "detail"])
result.to_csv(OUT / "independent_verification.tsv", sep="\t", index=False)
failed = result[~result.passed]
if len(failed): raise SystemExit("Verification failed:\n" + failed.to_string(index=False))
print(f"All {len(result)} independent checks passed")
