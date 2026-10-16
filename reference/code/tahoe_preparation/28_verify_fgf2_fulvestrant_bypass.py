#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/work/tahoe_prism_v2')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--iteration', default='04')
    args = ap.parse_args()
    iteration = f"{int(args.iteration):02d}"
    out = ROOT/'results/28_fgf2_fulvestrant_bypass'/f'iteration_{iteration}'
    checks = []

    def check(name, ok, detail=''):
        checks.append({'check': name, 'passed': bool(ok), 'detail': str(detail)})
        if not ok:
            raise AssertionError(f'{name}: {detail}')

    mem = pd.read_csv(out/'28b_frozen_program_genes.tsv', sep='\t')
    check('gene_count_30_to_60', 30 <= len(mem) <= 60, len(mem))
    bad = mem.gene.str.match(r'^(RPS(?!6K)\d|RPL\d|MRP[SL]\d|MT-|HIST|IG[HKL]|HB[ABDEGQZ]\d*$)')
    check('no_generic_or_ribosomal_genes', not bad.any(), ','.join(mem.loc[bad, 'gene']))
    check('three_biological_roles', set(mem.role) == {'ER', 'FGF', 'Root'}, sorted(mem.role.unique()))
    check('five_or_six_programs', mem.program.nunique() in (5, 6), mem.program.nunique())
    sizes = mem.groupby('program').size()
    check('program_minimum_size', sizes.min() >= 3, sizes.to_dict())
    genes = mem.gene.tolist()
    for fn in ('28b_tahoe_signed_correlations.tsv.gz', '28b_cyto_signed_correlations.tsv.gz',
               '28b_joint_signed_correlations.tsv.gz'):
        matrix = pd.read_csv(out/fn, sep='\t', compression='gzip', index_col=0)
        check(f'{fn}_square', matrix.shape == (len(genes), len(genes)), matrix.shape)
        check(f'{fn}_fixed_order', matrix.index.tolist() == genes and matrix.columns.tolist() == genes,
              'order mismatch')
        check(f'{fn}_symmetric', np.allclose(matrix.values, matrix.values.T, atol=1e-8, equal_nan=True),
              np.nanmax(abs(matrix.values-matrix.values.T)))
    coverage = pd.read_csv(out/'28c_gse153509_gene_coverage.tsv', sep='\t')
    measured = coverage.measured.astype(str).str.upper().eq('TRUE')
    check('all_frozen_genes_measured_in_gse153509', measured.all(), f'{measured.sum()}/{len(measured)}')
    criteria = pd.read_csv(out/'28c_primary_hypothesis_tests.tsv', sep='\t')
    passed = criteria.passed.astype(str).str.upper().eq('TRUE')
    check('all_primary_molecular_criteria_pass', passed.all(), criteria[['criterion', 'passed']].to_dict('records'))
    gfp = pd.read_csv(out/'28c_gfp_holdout_primary_criteria.tsv', sep='\t')
    gfp_passed = gfp.passed.astype(str).str.upper().eq('TRUE')
    check('untouched_gfp_t47d_criteria_pass', gfp_passed.all(), gfp[['contrast','role','passed']].to_dict('records'))
    manifest = dict(line.split('\t', 1) for line in (out/'28c_validation_manifest.tsv').read_text().splitlines())
    complete = int(manifest['n_extreme_outliers_removed']) == 0 and int(manifest['n_parental_samples']) == 79
    check('complete_parental_design_used', complete, manifest)
    functional = pd.read_csv(out/'28d_shee_functional_rescue.tsv', sep='\t')
    folds = functional.max_viable_cell_fold_change.tolist()
    check('published_functional_direction', folds[0] > folds[1] > 1, folds)
    figures = out/'figures'
    required = [
      '28b_01_joint_correlation','28b_02_tahoe_correlation','28b_03_cytosig_correlation',
      '28b_04_tahoe_gene_sample','28b_05_cytosig_gene_sample','28b_06_silhouette',
      '28b_07_program_fingerprint','28b_08_prism_association','28b_09_y_hierarchy',
      '28c_01_program_contrasts','28c_02_gse_gene_sample','28c_03_condition_program_heatmap',
      '28c_04_fgf_root_coupling','28c_05_validated_y_hierarchy','28d_01_shee_functional_rescue']
    required.insert(-1, '28c_06_gfp_holdout_contrasts')
    for stem in required:
        for ext in ('pdf', 'png'):
            path = figures/f'{stem}.{ext}'
            check(f'figure_{stem}_{ext}', path.exists() and path.stat().st_size > 5000,
                  path.stat().st_size if path.exists() else 'missing')
    code = '\n'.join((ROOT/'scripts'/name).read_text() for name in [
      '28_plot_fgf2_fulvestrant_discovery.R','28_validate_fgf2_fulvestrant_bypass.R'])
    check('no_ggplot_heatmap_implementation', not re.search(r'geom_(tile|raster)\s*\(', code),
          'forbidden heatmap geom')
    check('no_fixed_scatter_coordinates', not re.search(r'coord_(fixed|equal)\s*\(', code),
          'forbidden fixed coordinate')
    check('scatter_requirements', all(token in code for token in ['aspect.ratio=1','geom_smooth','stat_cor']),
          'missing required scatter layer')
    plot_code = (ROOT/'scripts/28_plot_fgf2_fulvestrant_discovery.R').read_text()
    dims = "width=unit(90,'mm'),height=unit(90,'mm')"
    check('square_correlation_heatmap_body', dims in plot_code, 'body dimensions not explicit/equal')
    report = {'iteration': iteration, 'passed': all(x['passed'] for x in checks), 'checks': checks}
    (out/'28e_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'iteration': iteration, 'passed': report['passed'], 'n_checks': len(checks)}, indent=2))

if __name__ == '__main__':
    main()
