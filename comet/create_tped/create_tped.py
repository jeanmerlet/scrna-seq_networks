import subprocess
import argparse
import re
import os
import numpy as np 
import pandas as pd

data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_mtx/heart_immune-dc-macrophage_comet-mtx.tsv'

def make_thresholds(count_mtx):
    # make a matrix uniform quantile transform-produced thresholds
    # each gene has a unique set of thresholds
    thresholds = np.empty([count_mtx.shape[0],20])
    for gene in range(count_mtx.shape[0]):
        # remove zero values
        g_nonzero = count_mtx[gene, :][count_mtx[gene, :] != 0]
        # make thresholds using only nonzero values
        idx = np.linspace(0, len(g_nonzero)-1, 21).astype(int)
        g_sorted = np.sort(g_nonzero)
        thresholds[gene, :] =  g_sorted[idx[0:20]]
    return thresholds.T


def make_celltype_tped(in_file):
    # load in .tsv of one celltype's gene expression values
    gene_by_cell = pd.read_csv(in_file, delimiter ='\t', header = 0)

    # remove rows where with genes with low variance
    gene_by_cell = gene_by_cell[gene_by_cell.var(axis=1)>0.000000001]
    
    # grab the gene names to make labels later
    gene_IDs = gene_by_cell.index.tolist()
    gene_by_cell = gene_by_cell.to_numpy()
    
    # construct the output file's path
    file_dir, file_name = os.path.split(in_file)
    celltype = file_name.split("_")[1]
    out_file_name = re.sub(r'tsv$', 'tped', file_name)
    out_file_dir = re.sub(r'comet_mtx', 'tped/'+celltype, file_dir)
    out_path = os.path.join(out_file_dir, '0.000000001_var_'+out_file_name)
     
    # get thresholds
    t = make_thresholds(gene_by_cell)

    # create output directory if it doesn't already exist
    try:
        os.mkdir(out_file_dir)
    except OSError:
        pass
    if os.path.exists(out_path):
        os.remove(out_path)

    # binarize gene expression values and write them into tped file
    at = np.array(['A','T'])
    for gene in range(gene_by_cell.shape[0]):
        expr = gene_by_cell[gene, :]
        g_binary = np.flip(((expr[:, None] >= t[:,gene][None, :]).astype(int)), axis=1)
        g_flat = g_binary.reshape(-1)
        line = list(at[g_flat])
        label = gene_IDs[gene]
        line = ['0', label, '0', '0'] + line 
        line[-1] = line[-1] + '\n'
        line = '\t'.join(line)
        with open(out_path, 'a') as f:
            f.writelines(str(line))

make_celltype_tped(data_dir)

#in_file_list = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if '-mtx.tsv' in f and 'other' not in f]
#in_file_list.sort()

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()
'''
for i, in_file_path in enumerate(in_file_list):
    # distribute across ranks
    if i % size != rank: continue
    head, file_name = os.path.split(in_file_path)
    print(f'making tped of {file_name} ({rank}/{size})')
    make_celltype_tped(in_file_path)
'''
