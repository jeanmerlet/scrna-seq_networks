from mpi4py import MPI
import numpy as np 
import pandas as pd
import subprocess
import argparse
import re, os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-r', '--run', action='store_true')
args = parser.parse_args()


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
    dir_up, _ = os.path.split(file_dir)
    tped_dir = os.path.join(dir_up, 'tped')
    celltype = '_'.join(file_name.split('_')[:-1])
    celltype_dir = os.path.join(tped_dir, celltype)
    out_file_name = re.sub(r'tsv$', 'tped', file_name)
    out_path = os.path.join(celltype_dir, '0.000000001_var_' + out_file_name)
    # create output directories if they don't already exist
    os.makedirs(tped_dir, exist_ok=True)
    os.makedirs(celltype_dir, exist_ok=True)
    # delete the out filepath if it already exists
    if os.path.exists(out_path):
        os.remove(out_path)

    # load in .tsv of one celltype's gene expression values
    gene_by_cell = pd.read_csv(in_file, delimiter='\t', header=0, index_col=0)

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


def check_tped_created(path):
    file_dir, file_name = os.path.split(path)
    dir_up, _ = os.path.split(file_dir)
    celltype = '_'.join(file_name.split('_')[:-1])
    tped_dir = os.path.join(dir_up, 'tped', celltype)
    out_file_name = re.sub(r'tsv$', 'tped', file_name)
    out_path = os.path.join(tped_dir, '0.000000001_var_' + out_file_name)
    if os.path.isfile(out_path):
        return True
    return False


in_file_list = []
for r, d, f in os.walk(args.data_dir):
    for comet_mtx in f:
        if 'comet-mtx.tsv' in comet_mtx:
            in_file_list.append(os.path.join(r, comet_mtx))
in_file_list.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


if args.run:
    new_file_list = []
    for f in in_file_list:
        if not check_tped_created(f):
            new_file_list.append(f)
    in_file_list = new_file_list
    for i, in_file_path in enumerate(in_file_list):
        # distribute filepaths to ranks
        # this assumes the number of files equals the number of ranks
        if i != rank: continue
        head, file_name = os.path.split(in_file_path)
        print(f'making tped from {file_name} ({rank}/{size})', flush=True)
        make_celltype_tped(in_file_path)
        print(f'Finished making tped from {file_name} ({rank}/{size})', flush=True)
else:
    print(f'total comet mtxs: {len(in_file_list)}')
    count = 0
    for path in in_file_list:
        if not check_tped_created(path):
            count += 1
    print(f'total tpeds to create: {count}')
