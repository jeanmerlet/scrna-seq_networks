from mpi4py import MPI
import subprocess
import argparse
import re, os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-r', '--run', action='store_true')
args = parser.parse_args()


def combine(tsv_dir):
    command = f'cat {tsv_dir}/HHLL_out*.tsv > {tsv_dir}/combined_all.tsv'
    subprocess.run(command, shell=True)


tsv_dirs = []
for r, d, f in os.walk(args.data_dir):
    for tsv_dir in d:
        if 'comet_out' in r:
            if 'archive' not in r:
                tsv_dirs.append(os.path.join(r, tsv_dir))
tsv_dirs.sort()


if not args.run:
    for d in tsv_dirs:
        print(d)
    print(len(tsv_dirs))
    raise SystemExit()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
for i, path in enumerate(tsv_dirs):
    # match ranks
    if i != rank: continue
    head, dir_name = os.path.split(path)
    print(f'postprocessing {dir_name} ({rank}/{size})', flush=True)
    combine(path)
    print(f'finished {dir_name} ({rank}/{size})', flush=True)
