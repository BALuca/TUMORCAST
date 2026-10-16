#!/usr/bin/env python3
"""Independent structural and numerical checks for Stage 36."""
from pathlib import Path
import csv, math, re

ROOT = Path('/work/tahoe_prism_v2')
OUT = ROOT / 'results/36_source_aware_braf_erbb_programs'

def rows(name):
    with (OUT / name).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

checks = []
def add(name, passed, detail=''):
    checks.append((name, bool(passed), detail))

genes = rows('frozen_program_genes.tsv')
counts = {p: sum(r['program'] == p for r in genes) for p in ('P1','P2','P3','P4')}
add('58 unique frozen genes', len(genes) == 58 and len({r['gene'] for r in genes}) == 58, str(counts))
add('Programs are sizeable', min(counts.values()) >= 10, str(counts))
add('RAF genes retained', {'ARAF','RAF1','HRAS'} <= {r['gene'] for r in genes})
add('No ribosomal genes', not any(re.match(r'^RP[SL]\d', r['gene']) for r in genes))

coh = rows('program_source_coherence.tsv')
required = [float(r['mean_signed_rho']) for r in coh if r['required_for_admission'] == 'TRUE']
add('Required coherence at least 0.10', min(required) >= 0.10, f'min={min(required):.3f}')
raf_cyto = [r for r in coh if r['program']=='P1' and r['source']=='CytoSig'][0]
add('RAF CytoSig evidence is descriptive', raf_cyto['required_for_admission']=='FALSE')

assoc = rows('matched_program_training_viability.tsv')
add('Every matched training Pearson at least 0.40', min(float(r['pearson']) for r in assoc) >= 0.40)
held = rows('heldout_prism_test.tsv')[0]
add('Heldout Pearson at least 0.40', float(held['pearson']) >= 0.40, held['pearson'])
add('Heldout RMSE improves', float(held['program_rmse']) < float(held['baseline_rmse']))

tahoe = rows('06_tahoe_all_profiles_gene_heatmap_sample_order.tsv')
cyto = rows('07_cytosig_all_profiles_gene_heatmap_sample_order.tsv')
add('All individual profiles plotted', len(tahoe)==252 and len(cyto)==2056, f'{len(tahoe)}/{len(cyto)}')
ordered = True
for p in {r['dominant_program'] for r in cyto}:
    values = [float(r['dominant_activity']) for r in cyto if r['dominant_program']==p]
    ordered &= all(a >= b - 1e-12 for a,b in zip(values, values[1:]))
add('CytoSig enrichment decreases within each program block', ordered)

pdfs = sorted((OUT/'figures/main').glob('*.pdf'))
pngs = sorted((OUT/'figures/main').glob('*.png'))
add('Nine main figures with PDF and PNG', len(pdfs)==9 and {p.stem for p in pdfs}=={p.stem for p in pngs})
add('Every figure is nonempty', all(p.stat().st_size > 5000 for p in pdfs+pngs))

with (OUT/'independent_verification.tsv').open('w', newline='') as f:
    w=csv.writer(f, delimiter='\t'); w.writerow(['check','passed','detail']); w.writerows((n,str(p).upper(),d) for n,p,d in checks)
failed=[n for n,p,_ in checks if not p]
if failed: raise SystemExit('FAILED: '+', '.join(failed))
print(f'Stage 36 independent verification passed ({len(checks)} checks)')

