import pandas as pd
import numpy as np
import os, csv
from scipy.io import mmread
from scipy.io import mmwrite

bc_meta_path = '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/meta/bc_meta.tsv'
bc_meta = pd.read_csv(bc_meta_path, sep='\t', index_col=0)
bc_meta.loc[:, 'idx'] = np.arange(bc_meta.shape[0])

adult_bcs_path = '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/meta/adult_bcs.txt'
with open(adult_bcs_path, 'rt') as f:
    reader = csv.reader(f)
    adult_bcs = [x[0] for x in reader]

idx = list(bc_meta.loc[adult_bcs, :].loc[:, 'idx'].values)
print('meta read and idx created', flush=True)

mtx_dir = '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/raw'
mtx_path = os.path.join(mtx_dir, 'matrix.mtx.gz')
mtx = mmread(mtx_path)
print('matrix read', flush=True)

mtx = mtx.todense()
print('mtx densed', flush=True)

mtx = mtx[:, idx]
mmwrite('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/gut/healthy/sanger/raw/adults/matrix.mtx', mtx)
print('mtx saved', flush=True)

