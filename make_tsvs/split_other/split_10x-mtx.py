from mpi4py import MPI
import pandas as pd
import numpy as np
import os, re
 

meta_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/meta/annot_humanAll.csv'
meta = pd.read_csv(meta_path, sep=',', header=0, index_col=None)
meta = meta.loc[meta.loc[:, 'typeSample'] == 'scRnaSeq', :]
valid_bcs = meta.loc[:, 'cell'].values
celltypes = meta.loc[:, 'annot']
celltypes.index = valid_bcs


mtx_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/raw/countTable_human'
mtx_path = os.path.join(mtx_dir, 'matrix.mtx')

genes_path = os.path.join(mtx_dir, 'features.tsv')
genes = pd.read_csv(genes_path, header=None, index_col=None, dtype=str)
genes = np.squeeze(genes.values)
num_genes = len(genes)

barcodes_path = os.path.join(mtx_dir, 'barcodes.tsv')
barcodes = pd.read_csv(barcodes_path, header=None, index_col=None).values[:, 0]


mtx_out_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/comet_mtx'


uniq_celltypes = np.unique(celltypes)
celltype_bc_sets = {}
for celltype in uniq_celltypes:
    celltype_bc_sets[celltype] = set()


with open(mtx_path, 'r') as mtx:
    # lines have the format 'gene_idx cell_idx num_umis'
    for i, line in enumerate(mtx):
        if i < 2: continue
        gene_idx, cell_idx, num_umis = line.strip().split(' ')
        bc = barcodes[int(cell_idx) - 1]
        if bc in celltypes.index:
            celltype = celltypes.loc[bc]
            celltype_bc_sets[celltype].add(bc)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, celltype in enumerate(uniq_celltypes):
    if i != rank: continue
    celltype_name = celltype.replace(' ', '_').replace('/', '-')
    print(f'{celltype_name} started...', flush=True)
    data = np.zeros((num_genes, len(celltype_bc_sets[celltype])), dtype=int)
    columns = list(celltype_bc_sets[celltype])
    comet_mtx = pd.DataFrame(data, index=genes, columns=columns)
    with open(mtx_path, 'r') as mtx:
        # lines have the format 'gene_idx cell_idx num_umis'
        for i, line in enumerate(mtx):
            if i < 2: continue
            gene_idx, cell_idx, num_umis = line.strip().split(' ')
            bc = barcodes[int(cell_idx) - 1]
            if bc in comet_mtx.columns:
                col_idx = comet_mtx.columns.get_loc(bc)
                comet_mtx.iloc[int(gene_idx) - 1, col_idx] = num_umis
    comet_mtx_out_path = os.path.join(mtx_out_dir, f'liver_{celltype_name}_comet-mtx.tsv')
    comet_mtx.to_csv(comet_mtx_out_path, sep='\t')
    print(f'{celltype_name} done', flush=True)
