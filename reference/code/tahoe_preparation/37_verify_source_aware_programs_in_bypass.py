#!/usr/bin/env python3
"""Independent structural and numerical checks for Stage 37."""
from pathlib import Path
import csv

ROOT=Path('/work/tahoe_prism_v2')
OUT=ROOT/'results/37_source_aware_bypass_validation'
SRC=ROOT/'results/36_source_aware_braf_erbb_programs'
def rows(path):
    with path.open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
checks=[]
def add(name,passed,detail=''): checks.append((name,bool(passed),detail))

source=rows(SRC/'frozen_program_genes.tsv'); copied=rows(OUT/'frozen_program_genes_exact_copy.tsv')
fields=('gene','program','orientation','plotting_order')
add('Exact Stage-36 program transfer', [[r[k] for k in fields] for r in source]==[[r[k] for k in fields] for r in copied])
cov=rows(OUT/'external_gene_coverage.tsv')
g127=[r for r in cov if r['dataset']=='GSE127988']; g684=[r for r in cov if r['dataset']=='GSE68453']
add('GSE127988 absences explicit', sum(int(r['n_frozen'])-int(r['n_measured']) for r in g127)==2 and {r['missing_genes'] for r in g127 if r['missing_genes']}=={'KIT','EREG'})
add('GSE68453 measures all genes', sum(int(r['n_measured']) for r in g684)==58)

auc=rows(OUT/'gse127988_egf_response_auc.tsv')
for p in ('P2','P3'):
    rr=[r for r in auc if r['program']==p]
    add(f'{p} activated by EGF on both days', len(rr)==2 and all(float(r['vemurafenib_auc'])>0 for r in rr))
    add(f'Cobimetinib attenuates {p} on both days', all(float(r['cobimetinib_change'])<0 for r in rr))

stats=rows(OUT/'gse68453_functional_prediction_stats.tsv'); p3=[r for r in stats if r['program']=='P3'][0]
add('P3 functional Pearson at least 0.40', float(p3['pearson'])>=0.40, p3['pearson'])
phos=rows(OUT/'gerosa_phosphoproteomic_tests.tsv'); erk=[r for r in phos if r['analyte'].startswith('ERK')]
add('EGF restores ERK at least threefold', len(erk)>=2 and all(float(r['EGF_rescue_over_one_uM_vem'])>=3 for r in erk))
pdfs=sorted((OUT/'figures/main').glob('*.pdf'));pngs=sorted((OUT/'figures/main').glob('*.png'))
add('Seven main figures with PDF and PNG',len(pdfs)==7 and {p.stem for p in pdfs}=={p.stem for p in pngs})
add('Single Aim 2 design schema present',(OUT/'figures/main/02_aim2_experimental_design.pdf').exists())
add('Collapsed EGF program heatmap present',(OUT/'figures/main/04_gse127988_collapsed_program_activity.pdf').exists())
collapsed=rows(OUT/'gse127988_collapsed_program_response.tsv')
add('Collapsed heatmap retains P1 context but excludes P4',len(collapsed)==12 and {r['program'] for r in collapsed}=={'P1','P2','P3'})
add('Every figure is nonempty',all(p.stat().st_size>5000 for p in pdfs+pngs))
with (OUT/'independent_verification.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['check','passed','detail']);w.writerows((n,str(p).upper(),d) for n,p,d in checks)
failed=[n for n,p,_ in checks if not p]
if failed:raise SystemExit('FAILED: '+', '.join(failed))
print(f'Stage 37 independent verification passed ({len(checks)} checks)')
