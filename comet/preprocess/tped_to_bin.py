import subprocess
import argparse
import os


def tped_to_bin(path, validate=False):
    # drop '.tped' from filename to create prefix
    prefix = path[:-5]
    subprocess.run([os.path.join(args.comet_tools_dir, 'line_labels.sh'), path, prefix + '_line_labels.txt'])
    subprocess.run([os.path.join(args.comet_tools_dir, 'line_indices'), path, prefix + '_line_indices.bin'])
    out_file_path = prefix + '_allele_labels.txt'
    out_file = open(out_file_path, 'w')
    subprocess.run(['sed', '-e', 's/.*/AT/', prefix + '_line_labels.txt'], stdout=out_file)
    out_file.close()
    subprocess.run([os.path.join(args.comet_tools_dir, 'preprocess'), 'duo', path, prefix + '.bin'])
    # validation can be slow
    if validate:
        subprocess.run([os.path.join(args.comet_tools_dir, 'preprocess_validate'), 'duo', prefix + '.bin', path, prefix + '_allele_labels.txt'])


def tped_to_bin_done(path):
    if os.path.isfile(path[:-4] + 'bin'):
        return True
    return False


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-c', '--comet_tools_dir')
parser.add_argument('-r', '--run', action='store_true')
args = parser.parse_args()


tped_paths = []
for r, d, f in os.walk(args.data_dir):
    for tped in f:
        if '.tped' in tped:
            tped_paths.append(os.path.join(r, tped))

tped_paths.sort()


if args.run:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    paths = []
    for f in tped_paths:
        if not tped_to_bin_done(f):
            paths.append(f)
    tped_paths = paths
    for i, tped_path in enumerate(tped_paths):
        # distribute filepaths to ranks
        # this assumes the number of files is equal to the number of ranks
        if i != rank: continue
        head, filename = os.path.split(tped_path)
        print(f'preprocessing tped {filename} ({rank}/{size})', flush=True)
        tped_to_bin(tped_path)
        print(f'Finished preprocessing tped {filename} ({rank}/{size})', flush=True)
else:
    print(f'number of tpeds: {len(tped_paths)}')
    count = 0
    for f in tped_paths:
        if not tped_to_bin_done(f):
            #print(f)
            count += 1
    print(f'total tpeds to binarize: {count}')
