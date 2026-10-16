"""Verify final ten-repeat stability and that frozen weights match checkpoints."""
from pathlib import Path
import hashlib, importlib.util, json
import numpy as np
import pandas as pd
BASE=Path('/work/almanac_conditionals/state_conditioned_bridge/program_first');OUT=BASE/'results';ART=BASE/'artifacts'
manifest=json.loads((OUT/'program_bank_freeze.json').read_text())
for name,expected in manifest.items():
    if name.endswith('.tsv'): assert hashlib.sha256((OUT/name).read_bytes()).hexdigest()==expected,name
sp=importlib.util.spec_from_file_location('helper','/work/tahoe_prism_v2/scripts/02_onmf_helpers.py');helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
helper.RANK=rank=manifest['rank'];genes=pd.read_csv(OUT/'selected_gene_universe.tsv',sep='\t').gene.to_numpy()
reps=[]
for seed in range(1,11):
    f=ART/f'_onmf_all_eligible_sensitivity_blind_v1_full_universe_rank_{rank}_runs_10'/f'gene_membership_seed_{seed}.npz'
    u=np.load(f)['gene_membership'];assert u.shape==(rank,len(genes));reps.append(u/np.linalg.norm(u,axis=1,keepdims=True))
u,matched,_,_=helper.consensus_membership(reps,genes)
signed,_=helper.orient_signed_programs(u,np.load(ART/'discovery_joint_correlation.npy',mmap_mode='r'))
a=u.argmax(axis=0);w=pd.read_csv(OUT/'frozen_program_weights.tsv',sep='\t');error=0
records=[]
for k in range(rank):
    program=f'P{k+1:02d}';ix=np.flatnonzero(a==k);actual=w[w.program==program].set_index('gene').loc[genes[ix]]
    error=max(error,float(np.max(np.abs(actual.weight.to_numpy()-signed[k,ix]))))
    q10=float(np.quantile(matched[:,k],.1));records.append(dict(program=program,ten_repeat_mean=float(matched[:,k].mean()),ten_repeat_q10=q10,ten_repeat_minimum=float(matched[:,k].min()),stable=q10>=.5,genes=len(ix)))
assert error<1e-10,error
pd.DataFrame(records).to_csv(OUT/'final_program_stability.tsv',sep='\t',index=False)
audit={'rank':rank,'repeats':10,'global_q10_stability':float(np.quantile(matched,.1)),'stable_programs':sum(x['stable'] for x in records),'minimum_program_size':min(x['genes'] for x in records),'negative_silhouette_fraction':float((w.silhouette<0).mean()),'frozen_weight_reconstruction_max_error':error,'frozen_file_hashes_match':True}
(OUT/'final_bank_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
