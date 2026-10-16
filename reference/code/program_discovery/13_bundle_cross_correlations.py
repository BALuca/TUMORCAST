from pathlib import Path
import fitz
base=Path('/work/almanac_conditionals/state_conditioned_bridge/program_first/figures/cross_correlations')
names=['01_joint_gene_correlations','02_tahoe_gene_correlations','03_cytosig_gene_correlations']
merged=fitz.open()
for name in names:
    with fitz.open(base/f'{name}.pdf') as doc:
        assert len(doc)==1
        merged.insert_pdf(doc)
        if name!=names[0]:doc[0].get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False).save(base/f'{name}.png')
assert len(merged)==3
merged.save(base/'cross_correlations_all_sources.pdf',garbage=4,deflate=True)
merged.close()
print('Bundled three checked source heatmaps.')
