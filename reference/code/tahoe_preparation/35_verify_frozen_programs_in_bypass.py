#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

root = Path(__file__).resolve().parents[1]
out = root / "results/35_frozen_program_bypass_validation"
aim1 = pd.read_csv(root / "results/34_aim1_cross_source_programs/frozen_program_genes.tsv", sep="\t")
aim2 = pd.read_csv(out / "frozen_program_genes_exact_copy.tsv", sep="\t")
assert aim1.equals(aim2), "Aim 2 changed the frozen Aim 1 program definition"
coverage = pd.read_csv(out / "external_gene_coverage.tsv", sep="\t")
assert coverage.query("dataset == 'GSE127988'").n_measured.tolist() == [20, 19, 14]
auc = pd.read_csv(out / "gse127988_egf_auc_and_cobimetinib_effect.tsv", sep="\t")
for p in ["P2", "P3"]:
    z = auc.query("program == @p")
    assert (z.vem_auc > 0).all() and (z.cobimetinib_change < 0).all()
assert (auc.query("program == 'P1'").vem_auc < 0).all()
fs = pd.read_csv(out / "gse68453_functional_prediction_stats.tsv", sep="\t")
assert (fs.pearson >= 0.4).all()
v = pd.read_csv(out / "verification.tsv", sep="\t")
assert v.passed.astype(bool).all()
assert len(list((out / "figures/main").glob("*.pdf"))) == 6
print("All Aim-2 frozen-transfer checks passed")
