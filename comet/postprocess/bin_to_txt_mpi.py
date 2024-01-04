import subprocess
import argparse
import re, os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-t', '--tped_dir')
parser.add_argument('-c', '--comet_tools_dir')
parser.add_argument('-w', '--num_way')
parser.add_argument('-r', '--run')
args = parser.parse_args()


def comet_postprocess(in_file, comet_tools_dir, tped_prefix, num_way):
    postprocess = os.path.join(comet_tools_dir, 'postprocess')
    allele_labels = tped_prefix + '_allele_labels.txt'
    line_labels = tped_prefix + '_line_labels.txt'
    out_file = re.sub(r'bin$', 'txt', in_file)
    subprocess.run([postprocess,
                    'duo',
                    num_way,
                    allele_labels,
                    line_labels,
                    in_file,
                    out_file])


def file_is_postprocessed(path):
    out_path = path[:-3] + 'txt'
    if os.path.isfile(out_path):
        return True
    return False


comet_out_paths = []
for r, d, f in os.walk(args.data_dir):
    for out_file in f:
        if re.match('^out.*bin$', out_file):
            comet_out_paths.append(os.path.join(r, out_file))
comet_out_paths.sort()


if args.run:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    for i, path in enumerate(comet_out_paths):
        # distribute across ranks
        if i % size != rank: continue
        head, file_name = os.path.split(path)
        print(f'postprocessing {file_name} ({rank}/{size})', flush=True)
        _, celltype = os.path.split(head)
        tped_prefix = os.path.join(args.tped_dir, celltype, 'heart_' + celltype + '_comet-mtx')
        comet_postprocess(path, args.comet_tools_dir, tped_prefix, args.num_way)
        print(f'finished {file_name} ({rank}/{size})', flush=True)
else:
    print(f'total number of comet out files: {len(comet_out_paths)}')
    count = 0
    for path in comet_out_paths:
        if not file_is_postprocessed(path):
            count += 1
        else:
            print(path)
    print(f'total files to postprocess: {count}')
