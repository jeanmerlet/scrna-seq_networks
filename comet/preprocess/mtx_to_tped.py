from mpi4py import MPI
import numpy as np 
import pandas as pd
import subprocess
import argparse
import re, os


def make_thresholds(count_mtx, max_bins=20):
    # apply a uniform quantile transform to each gene
    # to generate num_bins thresholds per gene
    # and return those thresholds in a matrix
    num_genes = count_mtx.shape[0]
    nz_genesums = np.sum(count_mtx != 0, axis=1)
    # restrict bin size to smallest number of nonzero values in a gene
    num_bins = min(max_bins, min(nz_genesums))
    thresholds = np.empty([num_genes, num_bins])
    for i in range(num_genes):
        gene = count_mtx[i, :]
        gene_nz = np.sort(gene[gene != 0])
        idx = np.ceil(np.linspace(0, len(gene_nz)-1, num_bins, endpoint=False))
        # prevent IndexError from floats
        idx = idx.astype(int)
        assert(len(idx) == len(np.unique(idx)))
        thresholds[i, :] = gene_nz[idx]
    return thresholds


def make_tped(mtx_path):
    # construct the output file's path
    file_dir, file_name = os.path.split(mtx_path)
    dir_up, _ = os.path.split(file_dir)
    tped_dir = os.path.join(dir_up, 'tped')
    celltype = '_'.join(file_name.split('_')[:-1])
    celltype_dir = os.path.join(tped_dir, celltype)
    out_file_name = re.sub(r'tsv$', 'tped', file_name)
    out_path = os.path.join(celltype_dir, out_file_name)
    # create output directories if they don't already exist
    os.makedirs(tped_dir, exist_ok=True)
    os.makedirs(celltype_dir, exist_ok=True)
    # delete the tped if it already exists
    if os.path.exists(out_path):
        os.remove(out_path)

    # read mtx
    gene_by_cell = pd.read_csv(mtx_path, delimiter='\t', header=0, index_col=0)
    # remove genes with very low variance
    gene_by_cell = gene_by_cell[gene_by_cell.var(axis=1) > 0.0000000001]
    # remove genes expressed in fewer than 5% of cells
    pct_cell_exp = np.sum(gene_by_cell != 0, axis=1) / gene_by_cell.shape[1]
    gene_by_cell = gene_by_cell.loc[pct_cell_exp.to_numpy() >= 0.05, :]
    # grab the gene names to make labels later
    gene_ids = gene_by_cell.index.tolist()
    gene_by_cell = gene_by_cell.to_numpy()
    # get thresholds
    t = make_thresholds(gene_by_cell)
    t = t.T
    # binarize gene expression values and write them into tped file
    at = np.array(['A','T'])
    for gene in range(gene_by_cell.shape[0]):
        expr = gene_by_cell[gene, :]
        g_binary = np.flip(((expr[:, None] >= t[:,gene][None, :]).astype(int)), axis=1)
        g_flat = g_binary.reshape(-1)
        line = list(at[g_flat])
        label = gene_ids[gene]
        line = ['0', label, '0', '0'] + line 
        line[-1] = line[-1] + '\n'
        line = '\t'.join(line)
        with open(out_path, 'a') as f:
            f.writelines(str(line))


def check_tped_created(mtx_path):
    file_dir, file_name = os.path.split(mtx_path)
    dir_up, _ = os.path.split(file_dir)
    tped_dir = os.path.join(dir_up, 'tped')
    celltype = '_'.join(file_name.split('_')[:-1])
    celltype_dir = os.path.join(tped_dir, celltype)
    out_file_name = re.sub(r'tsv$', 'tped', file_name)
    out_path = os.path.join(celltype_dir, out_file_name)
    if os.path.isfile(out_path):
        return True
    return False


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-r', '--run', action='store_true')
args = parser.parse_args()


mtx_paths = []
for r, d, f in os.walk(args.data_dir):
    for comet_mtx in f:
        if 'comet-mtx.tsv' in comet_mtx:
            mtx_paths.append(os.path.join(r, comet_mtx))
mtx_paths.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


if args.run:
    unprocessed_mtx_paths = []
    for path in mtx_paths:
        if not check_tped_created(path):
            unprocessed_mtx_paths.append(path)
    for i, path in enumerate(unprocessed_mtx_paths):
        # distribute filepaths to ranks
        # this assumes the number of files equals the number of ranks
        if i != rank: continue
        head, file_name = os.path.split(path)
        print(f'Creating tped from {file_name} ({rank}/{size})...', flush=True)
        make_tped(path)
        print(f'Finished making tped from {file_name} ({rank}/{size})', flush=True)
else:
    print(f'total comet mtxs: {len(mtx_paths)}')
    count = 0
    for path in mtx_paths:
        if not check_tped_created(path):
            print(path)
            count += 1
    print(f'total tpeds to create: {count}')
