from mpi4py import MPI
import numpy as np 
import pandas as pd
import subprocess
import argparse
import re, os


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'


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
    # construct the output file's path
    file_dir, file_name = os.path.split(in_file)
    out_file_name = re.sub(r'tsv$', 'tped', file_name)
    if 'brain' in file_dir:
        celltype = file_name.split("_")[2]
    else:
        celltype = file_name.split("_")[1]
    out_tped_dir = re.sub(r'comet_mtx', 'tped', file_dir)
    out_file_dir = os.path.join(out_tped_dir, celltype)
    out_path = os.path.join(out_file_dir, '0.000000001_var_' + out_file_name)
    # create output directories if they don't already exist
    try:
        os.mkdir(out_tped_dir)
    except OSError:
        pass
    try:
        os.mkdir(out_file_dir)
    except OSError:
        pass
    # delete the out filepath if it already exists
    if os.path.exists(out_path):
        os.remove(out_path)

    # load in .tsv of one celltype's gene expression values
    gene_by_cell = pd.read_csv(in_file, delimiter='\t', header=0)

    # remove low gene variance rows
    gene_by_cell = gene_by_cell[ gene_by_cell.var(axis=1) > 0.000000001 ]
    
    # grab the gene names to make labels later
    gene_IDs = gene_by_cell.index.tolist()
    gene_by_cell = gene_by_cell.to_numpy()
     
    # get thresholds
    t = make_thresholds(gene_by_cell)

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


in_file_list = []
for r, d, f in os.walk(data_dir):
    for comet_mtx in f:
        if 'comet-mtx.tsv' in comet_mtx:
            if 'heart' not in comet_mtx:
                in_file_list.append(os.path.join(r, comet_mtx))
in_file_list.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # distribute filepaths to ranks
    # this assumes the number of files equals the number of ranks
    if i != rank: continue
    head, file_name = os.path.split(in_file_path)
    print(f'making tped from {file_name} ({rank}/{size})', flush=True)
    make_celltype_tped(in_file_path)
