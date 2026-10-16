#!/usr/bin/env python3
"""Retrieve processed Tahoe profiles for two primary-PRISM-matched melanomas.

No single-cell reads are downloaded/reprocessed. Parquet footer statistics prune
cells and drug/dose row groups before HTTP range reads. Existing streamer code
is used only for byte-range and parquet I/O, never for cohort/analysis choices.
"""
from pathlib import Path
import importlib.util
import json
import os
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT.parent
ART = ROOT / 'artifacts/27_braf_erbb_hierarchy'
ART.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location('pdex_io', WORK / 'scripts/tahoe_prism_full_chain/01_stream_tahoe_full_gene.py')
io = importlib.util.module_from_spec(spec)
spec.loader.exec_module(io)
CELLS = {'CVCL_1666': ('ACH-000348', 'RPMI-7951'), 'CVCL_1097': ('ACH-000580', 'C32')}
DRUGS = ['Vemurafenib', 'Lapatinib ditosylate', 'Cobimetinib']
DOSES = [.05, .5, 5.]
requested = [repr([(d, q, 'uM')]) for d in DRUGS for q in DOSES]
source_cache = WORK / 'cache/tahoe_prism_full_chain/tahoe_pdex_selected'
genes = pd.read_csv(source_cache / 'genes.tsv', sep='\t').sort_values('gene_index').feature.tolist()
outputs = [ART / f for f in ['additional_treated_mean.npy', 'additional_reference_mean.npy', 'additional_response_metadata.tsv', 'additional_plate_audit.tsv']]
if all(p.exists() for p in outputs):
    print('Additional processed inputs already available. Nothing downloaded.')
else:
    sparse = ART / '.pdex_ranges.partial'
    fd, parquet = io.initialize_sparse_parquet(io.PDEX_URL, sparse)
    selected = []
    for rg in range(parquet.metadata.num_row_groups):
        lo, hi = io.statistic_value(parquet, rg, 'cell_line')
        if lo != hi or lo not in CELLS:
            continue
        tlo, thi = io.statistic_value(parquet, rg, 'target')
        if not any(tlo <= t <= thi for t in requested):
            continue
        start, end = io.parquet_row_group_span(parquet, rg)
        selected.append((rg, lo, start, end))
    print(f'{len(selected)} row groups; at most {sum(e-s+1 for _,_,s,e in selected)/1e6:.1f} MB of processed data.', flush=True)
    audit, profiles = [], []
    def receive(item):
        rg, cell, start, end = item
        content, content_range = io.fetch_small_range(io.PDEX_URL, f'bytes={start}-{end}')
        if len(content) != end-start+1:
            raise RuntimeError(f'Incomplete row group {rg}')
        os.pwrite(fd, content, start)
        return item
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for index, (rg, cell, start, end) in enumerate(pool.map(receive, selected), 1):
                table = parquet.read_row_group(rg, columns=['feature','target','target_mean','ref_mean','target_membership','ref_membership','plate','cell_line'])
                lengths, targets = io.run_blocks(table['target'])
                if any(n != len(genes) for n in lengths):
                    raise RuntimeError('Unexpected gene-block dimensions')
                offset = 0
                for n, target in zip(lengths, targets):
                    block = table.slice(offset, n); offset += n
                    drug, dose, unit = io.target_parts(target)
                    if drug not in DRUGS or dose not in DOSES:
                        continue
                    feature = block['feature'].to_pylist()
                    if set(feature) != set(genes) or len(set(feature)) != len(genes):
                        raise RuntimeError('Gene universe mismatch')
                    order = pd.Index(feature).get_indexer(genes)
                    plates = set(block['plate'].to_pylist())
                    cells = set(block['cell_line'].to_pylist())
                    if len(plates) != 1 or cells != {cell}:
                        raise RuntimeError('Ambiguous plate or cell identity')
                    plate = next(iter(plates))
                    treated = np.asarray(block['target_mean'].to_numpy())[order].astype(np.float64)
                    reference = np.asarray(block['ref_mean'].to_numpy())[order].astype(np.float64)
                    if not np.isfinite(treated).all() or not np.isfinite(reference).all():
                        raise RuntimeError('Non-finite processed abundance')
                    profiles.append((cell, drug, dose, plate, treated, reference))
                    audit.append(dict(cell_line=cell,depmap_id=CELLS[cell][0],cell_name=CELLS[cell][1],drug=drug,dose=dose,
                                      plate=plate,target_membership_mean=float(np.mean(block['target_membership'].to_numpy())),
                                      reference_membership_mean=float(np.mean(block['ref_membership'].to_numpy())),row_group=rg))
                if index % 10 == 0 or index == len(selected):
                    print(f'Processed {index}/{len(selected)} row groups; {len(profiles)} requested plate profiles.', flush=True)
        metadata, treated_rows, reference_rows = [], [], []
        for cell in CELLS:
            for drug in DRUGS:
                for dose in DOSES:
                    pp = [p for p in profiles if p[:3] == (cell,drug,dose)]
                    if not pp or len({p[3] for p in pp}) != len(pp):
                        raise RuntimeError(f'Missing/duplicate plate profiles: {cell} {drug} {dose}')
                    metadata.append(dict(array_row=len(metadata),cell_line=cell,depmap_id=CELLS[cell][0],cell_name=CELLS[cell][1],
                                         drug=drug,dose=dose,dose_unit='uM',n_plates=len(pp),plates=','.join(str(p[3]) for p in pp)))
                    treated_rows.append(np.mean([p[4] for p in pp],axis=0))
                    reference_rows.append(np.mean([p[5] for p in pp],axis=0))
        np.save(outputs[0], np.asarray(treated_rows,dtype='<f4'))
        np.save(outputs[1], np.asarray(reference_rows,dtype='<f4'))
        pd.DataFrame(metadata).to_csv(outputs[2],sep='\t',index=False)
        pd.DataFrame(audit).to_csv(outputs[3],sep='\t',index=False)
        pd.DataFrame({'gene': genes}).to_csv(ART/'additional_gene_order.tsv',sep='\t',index=False)
        (ART/'processed_download_manifest.json').write_text(json.dumps(dict(source=io.PDEX_URL,
            kind='Official processed treatment/reference means; no raw single-cell reprocessing',
            cells=CELLS,row_groups=len(selected),download_bytes=sum(e-s+1 for _,_,s,e in selected),
            profiles=len(profiles),gene_count=len(genes)),indent=2))
    finally:
        os.close(fd)
        sparse.unlink(missing_ok=True)

# All functional responses come from ONE primary-screen assay at 2.5 uM.
prism = WORK/'data/sciPlex_PRISM/PRISM'
treatments = pd.read_csv(prism/'primary-screen-replicate-collapsed-treatment-info.csv')
treatments = treatments[treatments.name.str.lower().isin(['vemurafenib','lapatinib','cobimetinib'])].copy()
assert len(treatments) == 3 and (treatments.dose == 2.5).all()
values = pd.read_csv(prism/'primary-screen-replicate-collapsed-logfold-change.csv',usecols=[0]+list(pd.read_csv(prism/'primary-screen-replicate-collapsed-logfold-change.csv',nrows=0).columns.get_indexer(treatments.column_name)))
values = values.rename(columns={values.columns[0]:'row_name'})
cell_info = pd.read_csv(prism/'primary-screen-cell-line-info.csv')
assert not cell_info.row_name.duplicated().any()
original = pd.read_csv(source_cache/'cells.tsv',sep='\t')
wanted = set(original.depmap_id) | {v[0] for v in CELLS.values()}
cell_info = cell_info[cell_info.depmap_id.isin(wanted)].copy()
long = values.melt(id_vars='row_name',var_name='column_name',value_name='prism_survival')
long = long.merge(treatments[['column_name','name','dose','screen_id']],on='column_name',validate='many_to_one')
long = long.merge(cell_info[['row_name','depmap_id','passed_str_profiling']],on='row_name',validate='many_to_one')
long = long.rename(columns={'name':'compound','dose':'functional_dose_uM'})
long['prism_sensitivity'] = -long.prism_survival
long.to_csv(ART/'primary_prism_responses.tsv',sep='\t',index=False)
print('Primary PRISM coverage:',long.groupby('compound').prism_survival.count().to_dict(),flush=True)
