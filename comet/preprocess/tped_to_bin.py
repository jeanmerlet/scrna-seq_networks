from mpi4py import MPI
import subprocess
import os


comet_tools_dir = '/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/genomics_gpu/tools'


def tped_to_bin(path, validate=False):
    # drop '.tped from filename to create prefix
    prefix = path[:-5]
    subprocess.run([os.path.join(comet_tools_dir, 'line_labels.sh'), path, prefix + '_line_labels.txt'])
    subprocess.run([os.path.join(comet_tools_dir, 'line_indices'), path, prefix + '_line_indices.bin'])
    out_file_path = prefix + '_allele_labels.txt'
    out_file = open(out_file_path, 'w')
    subprocess.run(['sed', '-e', 's/.*/AT/', prefix + '_line_labels.txt'], stdout=out_file)
    out_file.close()
    subprocess.run([os.path.join(comet_tools_dir, 'preprocess'), 'duo', path, prefix + '.bin'])
    if validate:
        subprocess.run([os.path.join(comet_tools_dir, 'preprocess_validate'), 'duo', prefix + '.bin', path, prefix + '_allele_labels.txt'])


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'
tped_paths = []
for r, d, f in os.walk(data_dir):
    for tped in f:
        if '.tped' in tped:
            if 'heart' not in tped:
                tped_paths.append(os.path.join(r, tped))

tped_paths.sort()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, tped_path in enumerate(tped_paths):
    # distribute filepaths to ranks
    # this assumes the number of files is gve to the number of ranks
    if i != rank: continue
    head, filename = os.path.split(tped_path)
    print(f'preprocessing tped {filename} ({rank}/{size})', flush=True)
    tped_to_bin(tped_path)
