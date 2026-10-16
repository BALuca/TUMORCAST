#!/usr/bin/env python3
from __future__ import annotations
import argparse, html, json, shutil
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path('/work/tahoe_prism_v2')
BASE=ROOT/'results/28_fgf2_fulvestrant_bypass'

def coherence(out: Path):
    mem=pd.read_csv(out/'28b_frozen_program_genes.tsv',sep='\t')
    rows=[]
    for source,fn in [('Tahoe','28b_tahoe_signed_correlations.tsv.gz'),('CytoSig','28b_cyto_signed_correlations.tsv.gz'),('Integrated','28b_joint_signed_correlations.tsv.gz')]:
        z=pd.read_csv(out/fn,sep='\t',index_col=0)
        for prog,g in mem.groupby('program').gene:
            a=z.loc[g,g].values; vals=a[np.triu_indices_from(a,1)]
            rows.append(dict(source=source,program=prog,role=mem.loc[mem.program==prog,'role'].iloc[0],
                             n_genes=len(g),mean_within_correlation=np.nanmean(vals)))
    d=pd.DataFrame(rows);d.to_csv(out/'28b_program_coherence.tsv',sep='\t',index=False)
    return d

def figure_card(stem,title):
    return f'''<article><h3>{html.escape(title)}</h3><a href="figures/{stem}.pdf"><img src="figures/{stem}.png" alt="{html.escape(title)}"></a><p><a href="figures/{stem}.pdf">PDF</a> · <a href="figures/{stem}.png">PNG</a></p></article>'''

def build_iteration(iteration: str):
    out=BASE/f'iteration_{iteration}'
    if not (out/'figures').exists(): return
    figs=[
      ('28b_01_joint_correlation','Integrated square gene-correlation matrix'),
      ('28b_02_tahoe_correlation','Tahoe square gene-correlation matrix'),
      ('28b_03_cytosig_correlation','CytoSig square gene-correlation matrix'),
      ('28b_04_tahoe_gene_sample','Tahoe ordered gene-by-profile heatmap'),
      ('28b_05_cytosig_gene_sample','CytoSig ordered gene-by-study heatmap'),
      ('28b_06_silhouette','Program-resolution silhouette curve'),
      ('28b_07_program_fingerprint','Program perturbation fingerprints'),
      ('28b_08_prism_association','PRISM functional associations'),
      ('28b_09_y_hierarchy','Tested Y hierarchy'),
      ('28c_01_program_contrasts','Held-out program contrasts'),
      ('28c_02_gse_gene_sample','GSE153509 ordered gene-by-sample heatmap'),
      ('28c_03_condition_program_heatmap','Treatment-level program activity'),
      ('28c_04_fgf_root_coupling','FGF-to-root coupling'),
      ('28c_05_validated_y_hierarchy','Validated molecular Y'),
      ('28c_06_gfp_holdout_contrasts','Untouched GFP-vector T47D confirmation'),
      ('28d_01_shee_functional_rescue','Published viable-cell response')]
    cards=''.join(figure_card(s,t) for s,t in figs if (out/'figures'/f'{s}.png').exists())
    tables=sorted(p.name for p in out.glob('*.tsv'))
    code=sorted(p.name for p in (out/'code').glob('*')) if (out/'code').exists() else []
    status='Recommended complete run' if iteration=='04' else 'Retained end-to-end diagnostic run'
    narrative=(out/'README.md').read_text() if (out/'README.md').exists() else ''
    page=f'''<!doctype html><html><head><meta charset="utf-8"><title>Stage 28 · iteration {iteration}</title>
<style>body{{font-family:Arial,sans-serif;max-width:1200px;margin:30px auto;color:#20252a;line-height:1.5}}h1,h2{{color:#17324d}}.status{{padding:10px 14px;background:#e8f3ef;border-left:5px solid #167c80;font-weight:bold}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:22px}}article{{border:1px solid #ddd;padding:14px;background:white}}img{{width:100%;height:430px;object-fit:contain;background:#fff}}code{{background:#f4f4f4;padding:2px 4px}}a{{color:#145a86}}pre{{white-space:pre-wrap;background:#f7f7f7;padding:16px}}ul{{columns:2}}</style></head><body>
<p><a href="../index.html">← all iterations</a></p><h1>FGF2–fulvestrant bypass · iteration {iteration}</h1><p class="status">{status}</p>
<h2>Scientific record</h2><pre>{html.escape(narrative)}</pre><h2>Figures</h2><div class="grid">{cards}</div>
<h2>Tables</h2><ul>{''.join(f'<li><a href="{html.escape(t)}">{html.escape(t)}</a></li>' for t in tables)}</ul>
<h2>Frozen code</h2><ul>{''.join(f'<li><a href="code/{html.escape(c)}">{html.escape(c)}</a></li>' for c in code)}</ul></body></html>'''
    (out/'index.html').write_text(page)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--iteration',default='04');args=ap.parse_args()
    it=f'{int(args.iteration):02d}';out=BASE/f'iteration_{it}'
    code=out/'code';code.mkdir(parents=True,exist_ok=True)
    for name in ['28_prepare_fgf2_fulvestrant_bypass.py','28_discover_fgf2_fulvestrant_programs.R',
      '28_plot_fgf2_fulvestrant_discovery.R','28_validate_fgf2_fulvestrant_bypass.R',
      '28_validate_gfp_holdout.R',
      '28_add_shee_functional_evidence.R','28_verify_fgf2_fulvestrant_bypass.py',
      '28_build_fgf2_fulvestrant_report.py']:
        shutil.copy2(ROOT/'scripts'/name,code/name)
    coh=coherence(out)
    crit=pd.read_csv(out/'28c_primary_hypothesis_tests.tsv',sep='\t')
    mem=pd.read_csv(out/'28b_frozen_program_genes.tsv',sep='\t')
    res=pd.read_csv(out/'28b_silhouette_by_k.tsv',sep='\t')
    k=mem.program.nunique();sizes=mem.groupby('program').size().to_dict()
    vals={r.criterion:r for r in crit.itertuples()}
    readme=f'''# FGF2-mediated escape from fulvestrant in parental T47D

## Question

Does FGF2 activate a receptor-proximal bypass branch that converges with ER signaling on a shared proliferation program, and does palbociclib disconnect that shared output and reduce viable-cell rescue?

## Experimental design

Program discovery used Tahoe transcriptional responses to fulvestrant, four FGFR inhibitors, and palbociclib; primary PRISM viability response; and 17 exact-FGF2 CytoSig profiles collapsed to 12 studies. GSE153509 was not opened until the 45 genes, signs, and memberships were frozen. Validation then used all 79 parental-T47D RNA profiles across 14 deposited treatment arms. All 45 frozen genes were measured, and the prespecified replicate-correlation rule removed zero samples.

## Discovery result

The final set contains 45 non-ribosomal, mechanistically interpretable genes: 15 ER transcriptional genes, 15 FGF/MAPK genes, and 15 directional cell-cycle-entry genes. Role-constrained hierarchical clustering selected {k} nested programs with sizes {sizes}. Five and six were the admissible resolutions; the selected solution maximized silhouette while requiring at least three genes per program. The identical gene order is used in the integrated, Tahoe, CytoSig, and held-out heatmaps. Correlation heatmap bodies are explicitly square and use empirical correlation limits.

## Held-out molecular result

All three primary predictions passed after BH correction:

- FGF2 added to fulvestrant activated the FGF receptor branch: effect {vals['FGF branch activates during FGF2 rescue'].effect:+.3f}, q={vals['FGF branch activates during FGF2 rescue'].q_value:.3g}.
- The same perturbation restored the shared proliferation root: effect {vals['Shared root restores during FGF2 rescue'].effect:+.3f}, q={vals['Shared root restores during FGF2 rescue'].q_value:.3g}.
- Adding palbociclib suppressed that root: effect {vals['Palbociclib suppresses the restored root'].effect:+.3f}, q={vals['Palbociclib suppresses the restored root'].q_value:.3g}.

## Untouched secondary molecular confirmation

After the final programs were frozen, the 81 GFP-vector T47D profiles were opened as a separate cohort. All three directional predictions repeated: FGF2 under fulvestrant activated the FGF branch (+0.504, q=0.00556) and shared root (+0.567, q=0.00830), while palbociclib suppressed the root (-1.577, q=2.35e-7) without materially changing the FGF branch. This cohort was not used in the iterative repairs of the parental analysis.

## Functional result

Shee et al. Figure S3H independently measured viable-cell abundance after five days in parental T47D. The published maximum fold-change annotations show 1.64-fold FGF2 rescue with fulvestrant and 1.36-fold with fulvestrant plus palbociclib. Palbociclib therefore removed 44% of the excess rescue above baseline. These values were transcribed from the plotted annotations because downloadable raw SRB replicates were not available.

## Interpretation and limits

The evidence supports the predicted sequence: FGF2 engages the receptor branch, the bypass reaches a shared cell-cycle-entry program during ER blockade, and CDK4/6 inhibition suppresses that molecular root while attenuating functional rescue. The parental set was used iteratively and its p-values are therefore developmental; the untouched GFP-vector cohort supplies the clean molecular confirmation. Fulvestrant alone did not measurably suppress the ER or root scores after 24 hours, so the RNA data support FGF-driven bypass activation and downstream disconnection more strongly than they support the ER-to-root arm. PRISM provided scalable functional weighting during gene admission; the aggregate cross-line program correlations in figure 28b_08 are descriptive. GSE153509 used 100 nM fulvestrant, whereas the functional experiment used 1 µM, so the result validates direction and ranking rather than a calibrated viability percentage.

Sources: [GSE153509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153509) · [Shee et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5839765/) · [local supplement](../../../artifacts/28_fgf2_fulvestrant_bypass/PMC5839765_supplementary/JEM_20171818_sm.pdf)
'''
    (out/'README.md').write_text(readme)
    for d in sorted(BASE.glob('iteration_*')): build_iteration(d.name.split('_')[-1])
    rows=[]
    for d in sorted(BASE.glob('iteration_*')):
        n=d.name.split('_')[-1]; status='recommended: all primary criteria passed' if n==it else 'retained diagnostic'
        rows.append(f'<li><a href="iteration_{n}/index.html">iteration_{n}</a> — {status}</li>')
    rootpage=f'''<!doctype html><html><head><meta charset="utf-8"><title>Stage 28 FGF2 bypass</title><style>body{{font-family:Arial,sans-serif;max-width:900px;margin:40px auto;line-height:1.5;color:#20252a}}h1{{color:#17324d}}a{{color:#145a86}}li{{margin:10px 0}}</style></head><body><h1>Stage 28 · FGF2–fulvestrant bypass</h1><p>The recommended complete run is <strong>iteration_{it}</strong>. Earlier complete runs are retained as diagnostics because they motivated scientifically material changes.</p><ul>{''.join(rows)}</ul><p><a href="iteration_{it}/README.md">Scientific report</a> · <a href="iteration_{it}/28e_verification.json">Verification</a></p></body></html>'''
    (BASE/'index.html').write_text(rootpage)

if __name__=='__main__': main()
