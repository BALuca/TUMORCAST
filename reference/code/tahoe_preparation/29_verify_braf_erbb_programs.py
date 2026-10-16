#!/usr/bin/env python3
"""Independently reconstruct Stage 29 transformations, scores and key statistics.

Uses NumPy/Pandas/SciPy, including original processed Tahoe arrays and author
RNA matrix. It does not fit alternative predictors or modify analysis outputs.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import re
from urllib.parse import unquote
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parents[1];WORK=ROOT.parent
OUT=ROOT/'results/29_braf_erbb_programs';ART=ROOT/'artifacts/29_braf_erbb_programs';PREP=ROOT/'artifacts/27_braf_erbb_hierarchy'
CACHE=WORK/'cache/tahoe_prism_full_chain';P=CACHE/'tahoe_pdex_selected'
checks=[]
def record(name,condition,detail=''):
    checks.append(dict(check=name,passed=bool(condition),detail=str(detail)))
def equal(name,a,b,tol=1e-9):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    ok=a.shape==b.shape and np.allclose(a,b,rtol=tol,atol=tol,equal_nan=True)
    error=np.nanmax(abs(a-b)) if a.shape==b.shape and np.isfinite(a-b).any() else None
    record(name,ok,f'max absolute difference: {error}')
def table(name,folder=OUT):return pd.read_csv(folder/name,sep='\t')
def matrix(name):return table('verify_'+name+'.tsv.gz',ART).set_index('id' if name in {'effects','raw_effects','auc','background'} else 'gene')

m=table('frozen_program_genes.tsv').sort_values('plotting_order');genes=m.gene.tolist()
meta=table('tahoe_primary_prism_observations.tsv');ld=table('dose_integrated_molecular_metadata.tsv')
cells=table('experiment_cell_lines.tsv');excluded={'ACH-000311','ACH-000628'}
record('Unique genes and six nontrivial subprograms',len(m)==120 and m.gene.nunique()==120 and m.program.nunique()==6 and m.groupby('program').size().min()>=14)
record('Frozen display order is contiguous',m.plotting_order.tolist()==list(range(1,len(m)+1)))
record('Outliers absent from every analysis cohort',not set(cells.loc[cells.analysis_included,'depmap_id'])&excluded)
record('Whole-line development/holdout separation',not set(cells.loc[cells.split=='Development','depmap_id'])&set(cells.loc[cells.split=='Holdout','depmap_id']))
record('Three independent V600 melanoma contexts',sum(cells.analysis_included&cells.melanoma_v600)==3)
for p,g in m.groupby('program'):
    equal(f'{p}: equal signed frozen weights',g.weight,g.orientation/len(g))
fp=table('raw_perturbation_coordinates.tsv.gz').set_index('gene')
expected_sign=np.where(m.role.eq('B'),-np.sign(fp.loc[genes,'vemurafenib']),np.sign(fp.loc[genes,'EGF']))
equal('All gene signs derived from discovery perturbations',m.orientation,expected_sign)
record('Nominal holdout contains seven lines',sum(cells.split.eq('Holdout'))==7)
record('Primary PRISM remains a single assay dose',set(meta.functional_dose_uM.dropna())=={2.5})
equal('Primary PRISM sensitivity is negative survival',meta.prism_sensitivity,-meta.prism_survival)

raw=matrix('raw_effects');effect=matrix('effects');auc=matrix('auc');bg=matrix('background')
bgm=table('plate_common_response_metadata.tsv',PREP);applied=table('plate_calibration_applied.tsv')
background_profiles=table('plate_background_observations.tsv',PREP)
record('Calibration excludes all three target compounds',not set(background_profiles.drug)&{'Vemurafenib','Lapatinib ditosylate','Cobimetinib'})
record('Calibration is available for all 35 retained lines',bgm.depmap_id.nunique()==35 and len(bgm)==210)
for r in applied.itertuples():
    plates=[int(p) for p in str(r.plates).split(',')]
    ix=bgm.index[(bgm.depmap_id==r.depmap_id)&bgm.plate.isin(plates)]
    equal('Plate subtraction: '+r.observation_id,effect.loc[r.observation_id],raw.loc[r.observation_id]-bg.iloc[ix].mean())
for r in ld.itertuples():
    rows=meta[(meta.depmap_id==r.depmap_id)&(meta.compound==r.compound)].sort_values('dose')
    record('All doses: '+r.key,rows.dose.tolist()==[.05,.5,5.])
    equal('Trapezoidal AUC: '+r.key,auc.loc[r.key],np.array([.25,.5,.25])@effect.loc[rows.observation_id].to_numpy())

# Reconstruct selected original and added molecular profiles from source means.
def normalize(t,r):
    t=np.maximum(np.asarray(t,dtype=float),0);r=np.maximum(np.asarray(r,dtype=float),0)
    tp=t/t.sum();rp=r/r.sum();ok=(tp*1e4>=.1)&(rp*1e4>=.1)
    mm=np.log2(tp[ok]/rp[ok]);aa=.5*np.log2(tp[ok]*rp[ok]);mq=np.quantile(mm,[.3,.7]);aq=np.quantile(aa,[.05,.95])
    mask=(mm>=mq[0])&(mm<=mq[1])&(aa>=aq[0])&(aa<=aq[1]);f=2**np.median(mm[mask])
    return np.log2(1+tp*1e4/f)-np.log2(1+rp*1e4)
source_meta=table('response_effect_metadata.tsv',CACHE)
pa=table('pdex_plate_audit.tsv.gz',P);rm=table('pdex_plate_reference_metadata.tsv',P)
ref=np.load(P/'pdex_plate_reference.npy',mmap_mode='r')
fullgenes=table('genes.tsv',P).sort_values('gene_index').feature
ix=pd.Index(fullgenes).get_indexer(raw.columns)
chosen_ids=['ACH-000750',cells.loc[cells.split.eq('Holdout'),'depmap_id'].iloc[0]]
for depmap in chosen_ids:
    sm=source_meta[source_meta.depmap_id==depmap];shard=np.load(CACHE/sm.shard_path.iloc[0],mmap_mode='r')
    for r in meta[meta.depmap_id==depmap].itertuples():
        sr=sm[sm.perturbation_id==r.perturbation_id].iloc[0]
        pp=pa[(pa.depmap_id==depmap)&(pa.perturbation_id==r.perturbation_id)].plate
        corrected=[]
        for plate in pp:
            rr=ref[rm[(rm.depmap_id==depmap)&(rm.plate==plate)].array_row.item()].astype(float)
            corrected.append(normalize(shard[int(sr.shard_row)].astype(float)+rr,rr)[ix])
        equal('Original processed-array reconstruction: '+r.observation_id,raw.loc[r.observation_id],np.mean(corrected,axis=0),tol=2e-6)
nt=np.load(PREP/'additional_treated_mean.npy');nr=np.load(PREP/'additional_reference_mean.npy');nm=table('additional_response_metadata.tsv',PREP)
for r in nm.itertuples():
    obs=meta[(meta.depmap_id==r.depmap_id)&(meta.drug==r.drug)&(meta.dose==r.dose)].observation_id.item()
    equal('Additional processed-array reconstruction: '+obs,raw.loc[obs],normalize(nt[r.array_row],nr[r.array_row])[ix],tol=2e-6)

# Independently reconstruct real-sample integrated signed correlation matrix.
ct=[]
for drug in ['vemurafenib','lapatinib','cobimetinib']:
    rows=meta[meta.compound.eq(drug)&meta.split.eq('Development')&meta.analysis_included&(meta.melanoma_v600 if drug=='vemurafenib' else True)]
    ct.append(effect.loc[rows.observation_id].corr(method='spearman'))
cyto=matrix('cyto_x');cc=cyto.T.corr(method='spearman')
joint=(sum(ct)/3+cc)/2
sign=m.set_index('gene').orientation.loc[joint.index]
joint=joint*np.outer(sign,sign)
saved=table('joint_gene_correlations.tsv.gz').set_index('gene')
equal('Joint source correlations reconstructed in Python',saved,joint.loc[saved.index,saved.columns],tol=1e-10)

# Held-out predictions independently refit from development rows only.
def scores(x,role):
    mm=m[m.role==role].set_index('gene');return x.loc[:,mm.index].mul(mm.orientation,axis=1).mean(axis=1)
root_scores=scores(auc,'R')
ld['R']=root_scores.loc[ld.key].to_numpy()
train=ld[ld.split.eq('Development')&ld.analysis_included&ld.eligible]
test=ld[ld.split.eq('Holdout')&ld.analysis_included&ld.eligible].copy()
drugs=sorted(train.compound.unique())
def design(df,root=False):
    cols=[np.ones(len(df))]+[df.compound.eq(d).to_numpy().astype(float) for d in drugs[1:]]
    if root:cols.append(df.R.to_numpy())
    return np.column_stack(cols)
y=train.prism_survival.to_numpy()
beta=np.linalg.lstsq(design(train),y,rcond=None)[0];beta_r=np.linalg.lstsq(design(train,True),y,rcond=None)[0]
pred=table('root_heldout_predictions.tsv').set_index('key').loc[test.key]
equal('Held-out baseline predictions from development-only fit',pred.predicted_baseline,design(test)@beta)
equal('Held-out root predictions from development-only fit',pred.predicted_root,design(test,True)@beta_r)
record('Effective functional holdout is six lines and eighteen observations',test.depmap_id.nunique()==6 and len(test)==18)
rx=rankdata(test.R);ry=rankdata(test.prism_survival);d=design(test)
rx-=d@np.linalg.lstsq(d,rx,rcond=None)[0];ry-=d@np.linalg.lstsq(d,ry,rcond=None)[0]
equal('Partial rank correlation adjusts both variables for drug',table('root_heldout_model_comparison.tsv').partial_spearman_controlling_drug,[np.corrcoef(rx,ry)[0,1]])

# Reconstruct GSE127988 from authors' numeric matrix and all individual wells.
rpm_path=ROOT/'artifacts/26_braf_erbb_bypass/gse127988/GSE127988_RPM_Gerosa_et_al.csv.gz'
with gzip.open(rpm_path,'rt') as handle:headers=[next(csv.reader(handle)) for _ in range(6)]
samples=table('gse127988_complete_sample_audit.tsv')
colix=[]
for r in samples.itertuples():
    candidates=[j for j in range(2,len(headers[0])) if headers[3][j]=='Individual' and headers[4][j]==str(r.sample_id)]
    record('Unique individual RNA well '+r.gsm,len(candidates)==1)
    j=candidates[0];colix.append(j)
    record('Author RNA condition '+r.gsm,[float(headers[k][j]) for k in [0,1,2,5]]==[r.vem_uM,r.cobi_uM,r.egf_h,r.biological_set])
rpm=pd.read_csv(rpm_path,skiprows=7,header=None,usecols=[0,1]+colix)
# pandas sorts integer usecols; explicitly restore the SOFT/GSM order.
rna_map=table('gse127988_gene_mapping.tsv.gz')
values=rpm.loc[:,colix].copy();values.columns=samples.gsm.tolist();values['gene']=rna_map.official_symbol.to_numpy()
rna=np.log2(1+values.dropna(subset=['gene']).groupby('gene').sum())
fm=table('external_gene_coverage.tsv');fm=fm[fm.dataset.eq('GSE127988')].set_index('gene').official_symbol
rx=rna.reindex(fm.loc[genes]).copy();rx.index=genes
full=samples[samples.biological_set.isin([3,4])];absolute=rx.loc[:,full.gsm].copy()
for day in [3,4]:
    controls=full.loc[full.biological_set.eq(day)&full.arm.eq('DMSO'),'gsm'];record(f'Four same-day untreated controls on day {day}',len(controls)==4)
    ids=full.loc[full.biological_set.eq(day),'gsm'];absolute.loc[:,ids]=absolute.loc[:,ids].sub(rx.loc[:,controls].mean(axis=1),axis=0)
actual=matrix('absolute127');equal('All absolute RNA displacements from author matrix',actual,absolute.loc[actual.index,actual.columns],tol=1e-8)
for role in ['B','E','R']:
    mm=m[m.role.eq(role)].set_index('gene');score=absolute.loc[mm.index].mul(mm.orientation,axis=0).mean(axis=0)
    saved=table('gse127988_hierarchy_vs_untreated.tsv');saved=saved[saved.role.eq(role)]
    equal('Independent RNA hierarchy score: '+role,saved.score,score.loc[saved.gsm],tol=1e-8)

bp=table('gse68453_bypass_prediction_vs_bliss.tsv')
expected=np.minimum((bp.E-bp.E.mean())/bp.E.std(ddof=1),(bp.R-bp.R.mean())/bp.R.std(ddof=1))
equal('Frozen AND predictor reconstructed',bp.BypassIndex,expected)
equal('All twelve locked efficacy values unchanged',bp.set_index('cell_line').loc[['A375','SKMEL24','HT144','DM331','SKMEL28','DM13','VMM5A','VMM18','DM6','VMM15','VMM12','VMM17'],'bliss'],[5.7,29.3,16.3,27.3,49.7,39.3,43,31.7,-2.7,-8.7,3.7,-1])
equal('Combination prediction statistic reconstructed',table('gse68453_bypass_test.tsv').spearman,[spearmanr(bp.BypassIndex,bp.bliss).statistic])

# Viewer/source accessibility and prohibited plotting patterns.
code=(ROOT/'scripts/29_refine_braf_erbb_programs.R').read_text()
record('No forbidden ggplot heatmaps or fixed-unit scatter axes',not re.search(r'geom_(?:tile|raster)\s*\(|coord_(?:fixed|equal)\s*\(',code))
html=(OUT/'index.html').read_text()
live_html=re.sub(r'<pre>.*?</pre>','',html,flags=re.S)
for link in set(re.findall(r'(?:href|src)="([^"]+)"',live_html)):
    if re.match(r'https?://|#',link) or link=='independent_verification.tsv':continue
    record('Viewer link: '+link,(OUT/unquote(link)).exists())
record('Viewer embeds the exact analysis source',''.join(code.split()) in ''.join(html.replace('&gt;','>').replace('&lt;','<').replace('&amp;','&').split()))
record('Every scientific PDF has a nonempty rendered PNG',all(p.with_suffix('.png').exists() and p.with_suffix('.png').stat().st_size>1000 for p in OUT.glob('29*.pdf')))
result=pd.DataFrame(checks);result.to_csv(OUT/'independent_verification.tsv',sep='\t',index=False)
files=[ROOT/'scripts'/name for name in ['27_prepare_braf_context.py','27_calibrate_plate_background.py','29_refine_braf_erbb_programs.R','29_verify_braf_erbb_programs.py']]+[ART/'frozen_programs.rds']
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]).to_csv(OUT/'source_and_freeze_sha256.tsv',sep='\t',index=False)
print(f'{result.passed.sum()}/{len(result)} independent checks passed.')
if not result.passed.all():
    print(result[~result.passed].to_string(index=False));raise SystemExit(1)
