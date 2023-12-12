from mpi4py import MPI
import subprocess
import argparse
import re, os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-t', '--tped_dir')
parser.add_argument('-c', '--comet_tools_dir')
parser.add_argument('-w', '--num_way')
args = parser.parse_args()


def comet_postprocess(in_file, comet_tools_dir, tped_prefix, num_way):
    postprocess = os.path.join(comet_tools_dir, 'postprocess')
    allele_labels = tped_prefix + '_allele_labels.txt'
    line_labels = tped_prefix + '_line_labels.txt'
    out_file = re.sub(r'bin$', 'txt', in_file)
    subprocess.run([postprocess,
                    num_way,
                    allele_labels,
                    line_labels,
                    in_file,
                    out_file])


in_file_list = []
for r, d, f in os.walk(args.data_dir):
    for out_file in f:
        if '.bin' in out_file:
            in_file_list.append(os.path.join(r, out_file))
in_file_list.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, in_file_path in enumerate(in_file_list):
    # distribute across ranks
    if i % size != rank: continue
    head, file_name = os.path.split(in_file_path)
    print(f'postprocessing {file_name} ({rank}/{size})', flush=True)
    _, celltype = os.path.split(head)
    tped_prefix = os.path.join(args.tped_dir, celltype, 'heart_' + celltype + '_comet-mtx')
    comet_postprocess(in_file_path, args.comet_tools_dir, tped_prefix, args.num_way)
