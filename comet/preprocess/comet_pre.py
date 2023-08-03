import subprocess
import os, re
from mpi4py import MPI

comet_tools_dir = '/ccs/proj/syb111/sw/frontier/comet/src/install_single_release_frontier/bin/'

def run_comet_pre(path):
    prefix = path[:-5]
    subprocess.run([os.path.join(comet_tools_dir, 'line_labels.sh'), path, prefix + '_line_labels.txt'])
    subprocess.run([os.path.join(comet_tools_dir, 'line_indices'), path, prefix + '_line_indices.bin'])
    out_file_path = prefix + '_allele_labels.txt'
    out_file = open(out_file_path, 'w')
    subprocess.run(['sed', '-e', 's/.*/AT/', prefix + '_line_labels.txt'], stdout=out_file)
    out_file.close()
    subprocess.run([os.path.join(comet_tools_dir, 'preprocess'), 'duo', path, prefix + '.bin'])
    subprocess.run([os.path.join(comet_tools_dir, 'preprocess_validate'), 'duo', prefix + '.bin', path, prefix + '_allele_labels.txt'])

root_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/tped'
tped_paths = []
for r, d, f in os.walk(root_dir):
    for tped in f:
        if '.tped' in tped:
            tped_paths.append(os.path.join(r, tped))

tped_paths.sort()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, path in enumerate(tped_paths):
    # distribute across ranks
    if i % size != rank: continue
    head, file_name = os.path.split(path)
    print(f'preprocessing {file_name} ({rank}/{size})')
    run_comet_pre(path)


