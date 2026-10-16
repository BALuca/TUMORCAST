#!/usr/bin/env python3
"""Estimate a common response from non-target compounds on matching Tahoe plates.

Each retained cell line is calibrated independently from its own profiles. This separates a
compound-specific displacement from the response shared by other perturbations;
it cannot distinguish technical plate bias from generic treatment biology.
Original LOX-IMVI measurements are local. The two added melanomas use range
reads of processed PDEX, never raw single-cell data.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import os
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];WORK=ROOT.parent
ART=ROOT/'artifacts/27_braf_erbb_hierarchy';CACHE=WORK/'cache/tahoe_prism_full_chain'
P=CACHE/'tahoe_pdex_selected'
DRUGS={'Vemurafenib','Lapatinib ditosylate','Cobimetinib'}
PLATES={1,2,3,7,8,9}
CELLS={'CVCL_1666':('ACH-000348','RPMI-7951'),'CVCL_1097':('ACH-000580','C32')}
genes=pd.read_csv(P/'genes.tsv',sep='\t').sort_values('gene_index').feature.tolist()
paths=[ART/'plate_common_response.npy',ART/'plate_common_response_metadata.tsv',ART/'plate_background_observations.tsv']
cached=all(p.exists() for p in paths)
existing_meta=pd.read_csv(paths[1],sep='\t') if cached else pd.DataFrame()
existing_array=np.load(paths[0]) if cached else None
existing_audit=pd.read_csv(paths[2],sep='\t').to_dict('records') if cached else []
local_cells=pd.read_csv(P/'cells.tsv',sep='\t')
local_cells=local_cells[~local_cells.depmap_id.isin({'ACH-000311','ACH-000628'})]
expected_ids=set(local_cells.depmap_id)|{v[0] for v in CELLS.values()}
if cached and set(existing_meta.depmap_id)==expected_ids and len(existing_meta)==6*len(expected_ids):
    print('All retained lines have plate calibration. Nothing downloaded.')
    raise SystemExit(0)
previous_manifest=json.loads((ART/'plate_background_manifest.json').read_text()) if (ART/'plate_background_manifest.json').exists() else {}
existing_keys=set(zip(existing_meta.depmap_id,existing_meta.plate)) if cached else set()


def effect(t,r):
    t=np.maximum(np.asarray(t,dtype=float),0);r=np.maximum(np.asarray(r,dtype=float),0)
    tp=t/t.sum();rp=r/r.sum();ok=(tp*1e4>=.1)&(rp*1e4>=.1)
    mm=np.log2(tp[ok]/rp[ok]);aa=.5*np.log2(tp[ok]*rp[ok]);mq=np.quantile(mm,[.3,.7]);aq=np.quantile(aa,[.05,.95])
    keep=(mm>=mq[0])&(mm<=mq[1])&(aa>=aq[0])&(aa<=aq[1]);factor=2**np.median(mm[keep])
    return (np.log2(1+tp*1e4/factor)-np.log2(1+rp*1e4)).astype(np.float32)

blocks={};audit=[]
def receive_profile(depmap,name,plate,drug,dose,t,r):
    if plate not in PLATES or drug in DRUGS or drug.upper() in {'DMSO','CONTROL'}:return
    key=(depmap,name,plate);blocks.setdefault(key,[]).append(effect(t,r))
    audit.append(dict(depmap_id=depmap,cell_name=name,plate=plate,drug=drug,dose=dose))

# Local background profile reconstruction is exactly treatment-minus-reference
# plus its matched plate reference, using all genes for composition scaling.
meta=pd.read_csv(CACHE/'response_effect_metadata.tsv',sep='\t')
meta=meta[meta.depmap_id.isin(local_cells.depmap_id)&(meta.available==1)&~meta.drug.isin(DRUGS)]
pa=pd.read_csv(P/'pdex_plate_audit.tsv.gz',sep='\t');pa=pa[pa.depmap_id.isin(local_cells.depmap_id)]
rm=pd.read_csv(P/'pdex_plate_reference_metadata.tsv',sep='\t')
refs=np.load(P/'pdex_plate_reference.npy',mmap_mode='r')
for cell in local_cells.itertuples():
    cm=meta[meta.depmap_id==cell.depmap_id];cp=pa[pa.depmap_id==cell.depmap_id]
    raw=np.load(CACHE/cm.shard_path.iloc[0],mmap_mode='r')
    for row in cm.itertuples():
        pp=cp[cp.perturbation_id==row.perturbation_id].plate.tolist()
        if len(pp)!=1 or pp[0] not in PLATES or (row.depmap_id,pp[0]) in existing_keys:continue
        ri=rm[(rm.depmap_id==row.depmap_id)&(rm.plate==pp[0])].array_row.item()
        ref=np.asarray(refs[ri],dtype=float)
        receive_profile(row.depmap_id,cell.cell_name,pp[0],row.drug,row.dose,raw[row.shard_row]+ref,ref)
    print(f'Local background: {cell.cell_name}; {len(audit)} new profiles.',flush=True)

selected=[]
if not all((v[0],p) in existing_keys for v in CELLS.values() for p in PLATES):
    spec=importlib.util.spec_from_file_location('pdex_io',WORK/'scripts/tahoe_prism_full_chain/01_stream_tahoe_full_gene.py')
    io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
    sparse=ART/'.plate_background_ranges.partial';fd,pq=io.initialize_sparse_parquet(io.PDEX_URL,sparse)
    selected=[]
    for rg in range(pq.metadata.num_row_groups):
        lo,hi=io.statistic_value(pq,rg,'cell_line')
        if lo!=hi or lo not in CELLS:continue
        plo,phi=io.statistic_value(pq,rg,'plate')
        if plo is not None and phi is not None and not any(float(plo)<=p<=float(phi) for p in PLATES):continue
        start,end=io.parquet_row_group_span(pq,rg);selected.append((rg,lo,start,end))
    print(f'Reading {len(selected)} processed row groups, {sum(e-s+1 for _,_,s,e in selected)/1e6:.1f} MB.',flush=True)
    def fetch(item):
        rg,cell,start,end=item
        data,_=io.fetch_small_range(io.PDEX_URL,f'bytes={start}-{end}')
        if len(data)!=end-start+1:raise RuntimeError('Incomplete processed range')
        os.pwrite(fd,data,start);return item
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for i,(rg,cell,start,end) in enumerate(pool.map(fetch,selected),1):
                tab=pq.read_row_group(rg,columns=['feature','target','target_mean','ref_mean','plate','cell_line'])
                lengths,targets=io.run_blocks(tab['target']);offset=0
                for n,target in zip(lengths,targets):
                    b=tab.slice(offset,n);offset+=n
                    if n!=len(genes):raise RuntimeError('Unexpected gene count')
                    plates=set(b['plate'].to_pylist())
                    if len(plates)!=1:raise RuntimeError('Ambiguous plate')
                    plate=int(next(iter(plates)))
                    drug,dose,unit=io.target_parts(target)
                    if plate not in PLATES or drug in DRUGS:continue
                    features=b['feature'].to_pylist();order=pd.Index(features).get_indexer(genes)
                    if len(set(features))!=len(genes) or min(order)<0:raise RuntimeError('Gene mapping mismatch')
                    t=b['target_mean'].to_numpy()[order];r=b['ref_mean'].to_numpy()[order]
                    if not np.isfinite(t).all() or not np.isfinite(r).all():raise RuntimeError('Invalid processed means')
                    receive_profile(*CELLS[cell],plate,drug,dose,t,r)
                if i%25==0 or i==len(selected):print(f'{i}/{len(selected)} row groups; {len(audit)} background profiles.',flush=True)
    finally:
        os.close(fd);sparse.unlink(missing_ok=True)
rows=existing_meta.to_dict('records') if cached else [];mat=list(existing_array) if cached else []
for (depmap,name,plate),xx in sorted(blocks.items()):
    if len(xx)<20:raise RuntimeError(f'Insufficient background profiles: {depmap}, plate {plate}, n={len(xx)}')
    rows.append(dict(array_row=len(rows),depmap_id=depmap,cell_name=name,plate=plate,n_profiles=len(xx)))
    mat.append(np.median(np.stack(xx),axis=0))
if len(rows)!=6*len(expected_ids):raise RuntimeError('Incomplete retained-line by six-plate calibration')
np.save(paths[0],np.stack(mat).astype(np.float32));pd.DataFrame(rows).to_csv(paths[1],sep='\t',index=False)
pd.DataFrame(existing_audit+audit).to_csv(paths[2],sep='\t',index=False)
(ART/'plate_background_manifest.json').write_text(json.dumps(dict(source='https://huggingface.co/datasets/tahoebio/tahoe-de-rhaister/resolve/main/pdex/all_plates_pdex.parquet',range_bytes_this_run=sum(e-s+1 for _,_,s,e in selected),prior_download=previous_manifest,background_profiles=len(existing_audit+audit),excluded_compounds=sorted(DRUGS),interpretation='Median non-target-compound response; technical and generic biological components are not separately identifiable'),indent=2))
print('Saved plate common-response calibration.',flush=True)
