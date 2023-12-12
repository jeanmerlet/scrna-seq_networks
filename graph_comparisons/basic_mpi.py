from mpi4py import MPI
import numpy as np
import os, subprocess
import network_metrics as nm


# count number of unique nodes
# TODO: check max gene length across all genes
def count_uniq_nodes(path, num_lines):
    genes = np.zeros(2*num_lines, dtype="S30")
    with open(path) as in_file:
        for i, line in enumerate(in_file):
            gene_a, gene_b, _ = line.strip().split('\t')
            genes[i] = gene_a
            genes[i+1] = gene_b
    num_uniq_nodes = len(np.unique(genes))
    return num_uniq_nodes


# chop up file at <path> into chunks of <chunk_size> number of lines
def get_chunks(path, num_chunks=4):
    num_lines_process = subprocess.run(['wc', '-l', path], capture_output=True, text=True)
    num_lines = int(num_lines_process.stdout.split(' ')[0])
    num_nodes = count_uniq_nodes(path, num_lines)
    chunk_size = num_nodes // num_chunks
    chunks = [list(range(i*chunk_size, (i+1)*chunk_size)) for i in range(num_chunks)]
    chunks[-1] = list(range(chunks[-1][0], num_nodes))
    return chunks
    

if __name__ == '__main__':
    # options
    in_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/networks/heart_adipocyte_irf-loop-mtx_no_correlated_data_top_0.01_pct_correlates_readded.tsv'
    out_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/networks'
    chunks = get_chunks(in_path)
    # mpi init
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # distribute across ranks
    for i, indices in enumerate(chunks):
        if i % size != rank: continue
        print(f'i: {i}, rank: {rank}, size: {size}')
        out_path = os.path.join(out_dir, f'heart_adipocyte_netx_irf-loop_{indices[0]}-{indices[-1]}.tsv')
        nm.get_metrics(in_path, indices, out_path)
