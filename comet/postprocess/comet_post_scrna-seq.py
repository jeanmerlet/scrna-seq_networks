import subprocess
import argparse
import re
import os
from mpi4py import MPI

# data directory
data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_out'

# postprocessing tools binary dir
comet_tools_dir = '/ccs/proj/syb111/sw/frontier/comet/src/genomics_gpu/tools'
num_way = '2'

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--start_idx')
parser.add_argument('-e', '--end_idx')
args = parser.parse_args()

def comet_postprocess(in_file):
    celltype = os.path.basename(os.path.dirname(in_file))
    tped_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(in_file))),'tped/'+celltype)
    postprocess = os.path.join(comet_tools_dir, 'postprocess')
    allele_labels = os.path.join(tped_dir, 'heart_'+celltype+'_comet-mtx_allele_labels.txt')
    line_labels =  os.path.join(tped_dir, 'heart_'+celltype+'_comet-mtx_line_labels.txt')
    out_file = re.sub(r'bin$', 'txt', in_file)
    subprocess.run([postprocess,
                    num_way,
                    allele_labels,
                    line_labels,
                    in_file,
                    out_file])
in_file = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_out/adipocyte/out_0.bin'
comet_postprocess(in_file)

"""
in_file_list = []
for r, d, f in os.walk(data_dir):
  for bin_file in f:
    if '.bin' in bin_file:
      in_file_list.append(os.path.join(r, bin_file))

in_file_list.sort()
in_file_list = in_file_list[int(args.start_idx):int(args.end_idx)]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # distribute across ranks

    print("rank: " + rank)
    print("size: " + size)
    if i % size != rank: continue
    file_path, file_name = os.path.split(in_file_path)
    print(f'postprocessing {file_name} ({rank}/{size})')
    comet_postprocess(in_file_path)
"""
