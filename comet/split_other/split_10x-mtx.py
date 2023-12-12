import pandas as pd
import numpy as np
import os
 

meta_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/meta/annot_humanAll.csv'
meta = pd.read_csv(meta_path, sep=',', header=0, index_col=None)
meta = meta.loc[meta.loc[:, 'typeSample'] == 'scRnaSeq', :]
valid_bcs = meta.loc[:, 'cell'].values
celltypes = meta.loc[:, 'annot']
celltypes.index = valid_bcs


mtx_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/raw/countTable_human'
mtx_path = os.path.join(mtx_dir, 'matrix.mtx')

genes_path = os.path.join(mtx_dir, 'features.tsv')
genes = pd.read_csv(genes_path, header=None, index_col=None)
num_genes = len(genes)

barcodes_path = os.path.join(mtx_dir, 'barcodes.tsv')
barcodes = pd.read_csv(barcodes_path, header=None, index_col=None).values[:, 0]


mtx_out_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/comet_mtx'


uniq_celltypes = np.unique(celltypes)
celltype_bc_sets = {}
for celltype in uniq_celltypes:
    celltype_bc_sets[celltype] = set()


with open(mtx_path) as mtx:
    # lines have the format 'gene_idx cell_idx num_umis'
    for i, line in enumerate(mtx):
        if i < 2: continue
        gene_idx, cell_idx, num_umis = line.strip().split(' ')
        bc = barcodes[int(cell_idx) - 1]
        if bc in celltypes.index:
            celltype = celltypes.loc[bc]
            celltype_bc_sets[celltype].add(bc)


for i, celltype in enumerate(uniq_celltypes):
    print(celltype, flush=True)
    data = np.zeros((num_genes, len(celltype_bc_sets[celltype])), dtype=int)
    columns = [''] * celltype_bc_sets[celltype]
    comet_mtx = pd.DataFrame(data, index=genes, columns=columns)
    colnames = np.full(len(celltype_bc_sets[celltype]), '', dtype=object)
    with open(mtx_path) as mtx:
        # lines have the format 'gene_idx cell_idx num_umis'
        for i, line in enumerate(mtx):
            if i < 2: continue
            gene_idx, cell_idx, num_umis = line.strip().split(' ')
            bc = barcodes[int(cell_idx)]
            if bc in celltypes.index and celltypes.loc[bc] == celltype:
                comet_mtx.iloc[int(gene_idx) - 1, int(cell_idx) - 1] = num_umis
                colnames[int(cell_idx) - 1] = bc
    comet_mtx.columns = colnames
    comet_mtx_out_path = os.path.join(mtx_out_dir, f'liver_{celltype}_comet-mtx.tsv')
    comet_mtx.to_csv(comet_mtx_out_path, sep='\t')
