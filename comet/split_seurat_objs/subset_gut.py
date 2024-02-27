import pandas as pd
import numpy as np
import os
from scipy.io import mmread

mtx_dir = '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/raw'
mtx_path = os.path.join(mtx_dir, 'matrix.mtx.gz')
mtx = mmread(mtx_path)
print('matrix read', flush=True)

bc_meta = pd.read_csv('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/meta/bc_meta.tsv', sep='\t', index_col=0)
adult_bcs = list(pd.read_csv('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/meta/adult_bcs.txt').values)

idx = []
for i, bc in enumerate(bc_meta.index.values):
    if i % 10000 == 0: print(i, flush=True)
    if bc in adult_bcs:
        idx.append(i)
print('indexing done', flush=True)

mtx = mtx.todense()
print('mtx densed', flush=True)

mtx = mtx[:, idx]
mmwrite('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/raw/adults/matrix.mtx', mtx)
print('mtx saved', flush=True)

