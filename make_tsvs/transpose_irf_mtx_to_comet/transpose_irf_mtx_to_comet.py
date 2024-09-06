from mpi4py import MPI
import pandas as pd
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
args = parser.parse_args()


def transpose(mtx_path):
    root, name = os.path.split(mtx_path)
    name = name.replace('irf-loop-mtx.tsv', 'comet-mtx.tsv')
    root, _ = os.path.split(root)
    out_dir = os.path.join(root, 'comet_mtx')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, name)
    mtx = pd.read_csv(mtx_path, sep='\t', index_col=0, header=0)
    mtx = mtx.T
    mtx.to_csv(out_path, sep='\t', index=True, header=True)


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'
mtx_paths = []
for r, d, f in os.walk(data_dir):
    for tsv in f:
        if 'irf-loop-mtx.tsv' in tsv:
            if 'Archive' not in r:
                mtx_paths.append(os.path.join(r, tsv))

mtx_paths.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


for i, mtx_path in enumerate(mtx_paths):
    if i != rank: continue
    transpose(mtx_path)
