#!/usr/bin/env python3
"""Retrieve Wilson Figure 1a and export its categorical cell-line x ligand matrix.

This is a transcription of the published figure, not gene-expression data.
Every category is independently checked against the original colored cell.
"""
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from PIL import Image

root = Path('/work/tahoe_prism_v2')
source = root / 'artifacts/21_program_drug_hierarchy/wilson_source'
result = root / 'results/21_program_drug_hierarchy'
source.mkdir(parents=True, exist_ok=True)
result.mkdir(parents=True, exist_ok=True)
urls = {
    'supplement.pdf': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11249/MediaObjects/41586_2012_BFnature11249_MOESM52_ESM.pdf',
    'figure1.ppt': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature11249/MediaObjects/41586_2012_BFnature11249_MOESM48_ESM.ppt',
    'figure1.jpg': 'https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnature11249/MediaObjects/41586_2012_Article_BFnature11249_Fig1_HTML.jpg',
    'pubmed_geo_links.xml': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gds&id=22763448',
}
manifest = []
for name, url in urls.items():
    path = source / name
    if not path.exists():
        with urlopen(Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=30) as response:
            path.write_bytes(response.read())
    data = path.read_bytes()
    manifest.append({'file': str(path), 'url': url, 'bytes': len(data),
                     'sha256': hashlib.sha256(data).hexdigest()})
assert (source / 'supplement.pdf').read_bytes().startswith(b'%PDF')
assert (source / 'figure1.ppt').read_bytes().startswith(bytes.fromhex('d0cf11e0a1b11ae1'))
(source / 'download_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

# Original Figure 1a row order, labels, six ligand calls and image row centers.
# HGF, FGF, IGF, PDGF, NRG1, EGF. N=NR, P=partial rescue, R=complete rescue.
groups = [
 ('HER2 amplified', 'lapatinib', [
  ('KPL4',       'PPNNRN', 127),
  ('EFM192A',    'NNNNRN', 154),
  ('SKBR3',      'NNNNRR', 180),
  ('BT474',      'NNNNRP', 207),
  ('AU565',      'PPNNRP', 234),
  ('HCC1419',    'NNNNRN', 260),
  ('UACC893',    'NPNNRR', 286),
  ('HCC1954',    'PNNNNP', 313),
  ('HCC202',     'NPNNRN', 340)]),
 ('NRG1 autocrine', 'lapatinib', [
  ('CHL-1',      'RRPNNN', 393),
  ('SAT',        'RRNNNR', 420),
  ('FaDu',       'RRNNNP', 447),
  ('PCI-6A',     'RRNNNR', 474)]),
 ('EGFR mutation', 'erlotinib', [
  ('PC9',        'RRPNNN', 527),
  ('HCC827',     'RPNNNN', 554),
  ('HCC4006',    'RRNNPN', 581),
  ('KHM-3S',     'RRPNNN', 607)]),
 ('EGF-like ligand autocrine', 'erlotinib', [
  ('CALU-3',     'PNNNNN', 660),
  ('H358',       'NNNNNN', 687),
  ('H1648',      'NNNNNN', 714),
  ('H1666',      'RRNNRN', 741)]),
 ('MET amplified', 'crizotinib', [
  ('KATO II',    'NRNNRR', 794),
  ('GTL16',      'NPNNRR', 821),
  ('EBC-1',      'NNNNPP', 848)]),
 ('PDGFR amplified', 'sunitinib', [
  ('NCI-H1703',  'RRNNNR', 901),
  ('A-204',      'NRNNNN', 928)]),
 ('FGFR amplified', 'PD173074', [
  ('MFM223',     'NPNNRP', 981),
  ('Kato III',   'NNNNNN', 1008),
  ('RT112',      'RNNNRR', 1035),
  ('RT4',        'NNNNRR', 1061)]),
 ('ALK translocation/mutation', 'TAE684', [
  ('KELLY',      'NPNNNR', 1115),
  ('SH-SY5Y',    'NPNNNN', 1141),
  ('NCI-H3122',  'RNNNRR', 1168),
  ('NCI-H2228',  'NNNNNN', 1195)]),
 ('BRAF(V600E) mutation', 'PLX4032', [
  ('M14',        'NPNNPN', 1248),
  ('NAE',        'RNNNRP', 1275),
  ('A-375',      'PRNNPN', 1301),
  ('SK-MEL-28',  'PNNNRN', 1328),
  ('Malme-3M',   'PRNNNP', 1355),
  ('COLO-201',   'PNNNNN', 1382),
  ('COLO-818',   'PNNNNN', 1408)]),
]
ligands = ['HGF', 'FGF', 'IGF', 'PDGF', 'NRG1', 'EGF']
im = np.array(Image.open(source / 'figure1.jpg').convert('RGB'))
assert im.shape == (1428, 946, 3), 'Image changed: recheck transcription coordinates.'
bounds = [(289, 335), (336, 378), (379, 415), (416, 473), (474, 531), (532, 575)]
wide, long, checks = [], [], []
for group, drug, records in groups:
    for cell, calls, y in records:
        assert len(calls) == 6 and set(calls) <= set('NPR')
        row = {'figure_row': len(wide) + 1, 'cell_line': cell, 'kinase_dependency': group,
               'kinase_inhibitor': drug}
        for ligand, call, (left, right) in zip(ligands, calls, bounds):
            rgb = np.median(im[y-8:y+9, left+4:right-4], axis=(0, 1))
            r, g, b = rgb
            observed = 'R' if g > r + 40 and g > b + 30 else (
                'P' if b > r + 40 and b > g + 20 else ('N' if min(rgb) > 220 else '?'))
            assert observed == call, (cell, ligand, call, observed, rgb)
            category = 'NR' if call == 'N' else call
            row[ligand] = category
            long.append({**{k: row[k] for k in ['figure_row', 'cell_line', 'kinase_dependency', 'kinase_inhibitor']},
                         'ligand': ligand, 'rescue_category': category,
                         'assay': 'categorical viability rescue under inhibitor plus ligand',
                         'source': 'Wilson et al. 2012, Figure 1a',
                         'doi': '10.1038/nature11249',
                         'extraction': 'manual figure transcription; independent pixel-color check'})
            checks.append({'cell_line': cell, 'ligand': ligand, 'category': category,
                           'image_y': y, 'image_x_left': left, 'image_x_right': right,
                           'median_r': r, 'median_g': g, 'median_b': b, 'pixel_verified': True})
        wide.append(row)
assert len(wide) == 41 and len(long) == 246
assert len({row['cell_line'] for row in wide}) == 41

def write(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

write(result / 'wilson_ligand_rescue_matrix.tsv', wide)
write(result / 'wilson_ligand_rescue_long.tsv', long)
write(source / 'figure1a_transcription_checks.tsv', checks)
note = '''# Wilson 2012 downloaded data

Source: https://www.nature.com/articles/nature11249 (PMID 22763448).
The publisher's Figure 1 image, Figure 1 PowerPoint and 34-page supplement were downloaded. URLs and SHA256 digests are saved in download_manifest.json.

wilson_ligand_rescue_matrix.tsv is a transcription of Figure 1a: 41 cell lines by six RTK ligands, with the relevant kinase inhibitor for each line. NR=no rescue; P=partial rescue; R=complete rescue. The published cells were read manually and all 246 categorical calls were independently verified against their image background colors. Original row, cell-line and ligand labels are retained. These categories have no continuous measurement units.

This is a functional rescue matrix, not a genes-by-samples expression matrix. The study reports viability, immunoblots and associated assays. The supplementary material contains Figures 1-13, compound synthesis and a list of the 446 screened secreted factors. No perturbation transcriptome matrix or expression accession was identified in the article or supplement. The saved NCBI PubMed-to-GEO link response contains no linked dataset.

These data can test functional rescue predictions after a molecular model is specified; they cannot directly measure the eight gene-program scores. No program scores, expression values or combined drug/ligand transcriptomes are inferred from the rescue categories. No Tahoe/CytoSig cell-line substitute is used as Wilson validation.

Run: /work/.venv/bin/python /work/tahoe_prism_v2/scripts/download_wilson_matrix.py
'''
(source / 'README.md').write_text(note)
(result / 'wilson_data_availability.md').write_text(note)
print('Downloaded original Wilson sources and exported the 41 x 6 functional rescue matrix.')
print('Verified every one of the 246 categories against the original figure; no gene-expression matrix identified.')
