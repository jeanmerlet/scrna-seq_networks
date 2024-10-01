import numpy as np
import pandas as pd
import subprocess
import argparse
import re, os

# this script combines the tsvs then max normalizes


def get_max(path):
    weights = pd.read_csv(path, sep='\t', usecols=[2], index_col=None,
                          header=None).values
    return np.max(weights)


def max_norm(in_path, out_path):
    max_weight = get_max(in_path)
    with open(out_path, 'wt') as out_file:
        with open(in_path, 'rt') as in_file:
            for line in in_file:
                id1, id2, weight = line.strip().split('\t')
                weight = str(round(float(weight) / max_weight, 6))
                out_line = '\t'.join([id1, id2, weight]) + '\n'
                out_file.write(out_line)


def combine(tsv_dir, out_path):
    command = f'cat {tsv_dir}/out*.txt > {out_path}'
    subprocess.run(command, shell=True)


def check_combined(tsv_dir):
    filepaths = [os.path.join(tsv_dir, f) for f in os.listdir(tsv_dir)]
    for path in filepaths:
        if re.search('combined.tsv$', path):
            return True
    return False


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-r', '--run', action='store_true')
args = parser.parse_args()


tsv_dirs = []
for r, d, f in os.walk(args.data_dir):
    for tsv_dir in d:
        if 'comet_out' in r:
            tsv_dirs.append(os.path.join(r, tsv_dir))
tsv_dirs.sort()


if not args.run:
    counter = 0
    for d in tsv_dirs:
        if check_combined(d):
            counter += 0
    print(f'total num tsv dirs: {len(tsv_dirs)}')
    print(f'num tsvs to combine: {counter}')
else:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    for i, tsv_dir in enumerate(tsv_dirs):
        # match ranks
        if i != rank: continue
        head, dir_name = os.path.split(tsv_dir)
        print(f'postprocessing {dir_name} ({rank}/{size})', flush=True)
        out_path = os.path.join(tsv_dir, 'combined_all.tsv')
        out_norm_path = os.path.join(tsv_dir, 'combined_all_max-norm.tsv')
        combine(tsv_dir, out_path)
        max_norm(out_path, out_norm_path)
        print(f'finished {dir_name} ({rank}/{size})', flush=True)
