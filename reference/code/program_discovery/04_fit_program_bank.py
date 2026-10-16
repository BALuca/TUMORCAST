"""Sensitivity-blind rank screening, followed by ten-repeat signed oNMF."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','16')
os.environ.setdefault('OMP_NUM_THREADS','16')
from pathlib import Path
import importlib.util, json, hashlib, time, shutil, sys, subprocess
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_samples
BASE=Path('/work/almanac_conditionals/state_conditioned_bridge/program_first'); OUT=BASE/'results'; ART=BASE/'artifacts'
sp=importlib.util.spec_from_file_location('helper','/work/tahoe_prism_v2/scripts/02_onmf_helpers.py');helper=importlib.util.module_from_spec(sp);sp.loader.exec_module(helper)
genes=pd.read_csv(OUT/'selected_gene_universe.tsv',sep='\t').gene.to_numpy()

def digest_file(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for block in iter(lambda:f.read(16*1024*1024),b''):h.update(block)
 return h.hexdigest()
source_files=[OUT/'selected_gene_universe.tsv',ART/'discovery_tahoe_correlation.npy',ART/'discovery_cytosig_correlation.npy']
source_hashes={f.name:digest_file(f) for f in source_files}
proof=OUT/'discovery_input_hashes.json'
if proof.exists():
 assert json.loads(proof.read_text())==source_hashes, 'Discovery inputs changed; refuse old factorization checkpoints.'
else: proof.write_text(json.dumps(source_hashes,indent=2))

tahoe=np.load(ART/'discovery_tahoe_correlation.npy',mmap_mode='r');cyto=np.load(ART/'discovery_cytosig_correlation.npy',mmap_mode='r');assert tahoe.shape==cyto.shape==(len(genes),len(genes))
joint_path=ART/'discovery_joint_correlation.npy'
if not joint_path.exists():
 joint=np.tanh((np.arctanh(np.clip(tahoe,-.999,.999))+np.arctanh(np.clip(cyto,-.999,.999)))/2);np.fill_diagonal(joint,1);np.save(joint_path,joint)
joint=np.load(joint_path,mmap_mode='r');assert np.isfinite(joint).all()
helper.CACHE=ART;helper.N_RESPONSE_GENES=len(genes);helper.ONMF_MAX_ITER=500;helper.PROGRAM_DEFINITION_VERSION='all_eligible_sensitivity_blind_v1'
# Same signed affinity and original dspin algorithm; double precision avoids repeated dtype conversion inside dspin.
affinity=helper.posneg_correlation(joint).astype(np.float64);distance=np.clip(1-np.abs(joint),0,None);np.fill_diagonal(distance,0)
protocol={'ranks':[10,20,40,80],'screening_repeats':3,'final_repeats':10,'iterations':500,'rank_rule':'Elbow of relative reconstruction improvement versus rank, among q10 stability >=0.50 and minimum program size >=10. Original separation/source gates reported separately. No sensitivity or Wilson outcomes used.','program_selection':'After freeze: Pearson association of weighted measured dose-summary displacement with sensitivity; BH q<0.05 across all bank programs and both focal drugs.'}
(OUT/'discovery_protocol.json').write_text(json.dumps(protocol,indent=2))
def fit(rank,nruns):
 helper.RANK=rank;helper.N_ONMF_RUNS=nruns
 # Original checkpoint signature includes repeat count, so separate stage paths avoid reuse ambiguities.
 if nruns==10:
  screen=ART/f'_onmf_{helper.PROGRAM_DEFINITION_VERSION}_full_universe_rank_{rank}_runs_3'
  final=ART/f'_onmf_{helper.PROGRAM_DEFINITION_VERSION}_full_universe_rank_{rank}_runs_10'
  if not final.exists():
   final.mkdir()
   signature=helper.onmf_signature(f'full_universe_rank_{rank}_runs_10',{'all_eligible_gene_compendia_AU565_excluded'},genes)
   (final/'input_signature.sha256').write_text(signature+'\n')
   for seed in range(1,4): shutil.copy2(screen/f'gene_membership_seed_{seed}.npz',final/f'gene_membership_seed_{seed}.npz')
 reps=helper.run_repeated_onmf(f'full_universe_rank_{rank}_runs_{nruns}',affinity,genes,{'all_eligible_gene_compendia_AU565_excluded'})
 consensus,matched,_,_=helper.consensus_membership(reps,genes);signed,_=helper.orient_signed_programs(consensus,joint);assigned=consensus.argmax(axis=0)
 return consensus,signed,assigned,matched
screen_rank=int(sys.argv[2]) if len(sys.argv)>2 and sys.argv[1]=='--screen-rank' else None
if screen_rank is None:
 jobs=[]
 for rank in protocol['ranks']:
  if (OUT/f'rank_screening_rank_{rank}.json').exists(): continue
  log=open(ART/f'rank_{rank}_screen.log','w')
  env=dict(os.environ,OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
  jobs.append((rank,subprocess.Popen([sys.executable,__file__,'--screen-rank',str(rank)],stdout=log,stderr=subprocess.STDOUT,env=env),log))
 for rank,job,log in jobs:
  code=job.wait();log.close()
  if code: raise RuntimeError(f'Rank {rank} screening failed: {code}')
rows=[];solutions={}
for rank in ([screen_rank] if screen_rank is not None else protocol['ranks']):
 if screen_rank is None:
  rows.append(json.loads((OUT/f'rank_screening_rank_{rank}.json').read_text()))
  continue
 started=time.monotonic();u,s,a,matched=fit(rank,3);solutions[rank]=(u,s,a,matched)
 # Least-squares reconstruction of the unchanged signed affinity from the consensus membership.
 # This measures how much of the discovery matrix this consensus spans; it never uses functional labels.
 A=u.T;AtX=A.T@affinity;gram=A.T@A;coef=np.linalg.solve(gram,AtX);norm2=np.sum(affinity*affinity);residual=max(0.,norm2-2*np.sum(coef*AtX)+np.sum((gram@coef)*coef))/norm2
 sil=silhouette_samples(distance,a,metric='precomputed')
 source=[]
 for k in range(rank):
  ix=np.flatnonzero(a==k);sign=np.sign(s[k,ix]);src=[]
  for C in (tahoe,cyto):
   block=C[np.ix_(ix,ix)]*np.outer(sign,sign);src.append(float(block[np.triu_indices(len(ix),1)].mean()) if len(ix)>1 else np.nan)
  source.append(src)
 row={'rank':rank,'relative_reconstruction_error':residual,'mean_stability':float(matched.mean()),'q10_stability':float(np.quantile(matched,.1)),'minimum_program_size':int(np.bincount(a,minlength=rank).min()),'negative_silhouette_fraction':float(np.mean(sil<0)),'mean_silhouette':float(sil.mean()),'minimum_source_coherence':float(np.nanmin(source)),'seconds':time.monotonic()-started}
 rows.append(row);(OUT/f'rank_screening_rank_{rank}.json').write_text(json.dumps(row,indent=2));print('RANK_DIAGNOSTIC '+json.dumps(row),flush=True)
if screen_rank is not None: sys.exit(0)
df=pd.DataFrame(rows);ranks=df['rank'].to_numpy();err=df.relative_reconstruction_error.to_numpy();gain=(err[0]-err)/max(err[0]-err[-1],1e-12);position=(ranks-ranks[0])/(ranks[-1]-ranks[0]);df['elbow_score']=gain-position;df['stability_eligible']=(df.q10_stability>=.5)&(df.minimum_program_size>=10);df['original_separation_gate']=df.negative_silhouette_fraction<=.2;df['original_source_gate']=df.minimum_source_coherence>=.05
eligible=df[df.stability_eligible];assert len(eligible), 'No rank meets stability/size criteria; do not force a bank.'
chosen=int(eligible.sort_values(['elbow_score','rank'],ascending=[False,True]).iloc[0]['rank']);df['chosen']=df['rank']==chosen;df.to_csv(OUT/'rank_screening.tsv',sep='\t',index=False)
print(f'RANK_FROZEN {chosen}; starting ten-repeat final bank',flush=True)
u,s,a,matched=fit(chosen,10);sil=silhouette_samples(distance,a,metric='precomputed');records=[];coherence=[]
for k in range(chosen):
 ix=np.flatnonzero(a==k);sign=np.sign(s[k,ix]);name=f'P{k+1:02d}'
 for j in ix:records.append({'program':name,'gene':genes[j],'weight':float(s[k,j]),'absolute_weight':float(u[k,j]),'orientation':int(np.sign(s[k,j])),'silhouette':float(sil[j])})
 vals=[]
 for C in (tahoe,cyto):
  block=C[np.ix_(ix,ix)]*np.outer(sign,sign);vals.append(float(block[np.triu_indices(len(ix),1)].mean()) if len(ix)>1 else np.nan)
 coherence.append({'program':name,'genes':len(ix),'mean_repeat_stability':float(matched[:,k].mean()),'minimum_repeat_stability':float(matched[:,k].min()),'Tahoe_coherence':vals[0],'CytoSig_coherence':vals[1]})
pd.DataFrame(records).to_csv(OUT/'frozen_program_weights.tsv',sep='\t',index=False);pd.DataFrame(coherence).to_csv(OUT/'program_quality.tsv',sep='\t',index=False);(OUT/'selected_rank.txt').write_text(str(chosen)+'\n')
manifest={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [OUT/'frozen_program_weights.tsv',OUT/'rank_screening.tsv',OUT/'selected_gene_universe.tsv']};manifest['sensitivity_used']=False;manifest['AU565_used']=False;manifest['rank']=chosen;(OUT/'program_bank_freeze.json').write_text(json.dumps(manifest,indent=2));print('BANK_FROZEN '+json.dumps(manifest),flush=True)
