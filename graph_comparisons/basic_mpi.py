from mpi4py import MPI
import os, subprocess
import network_metrics as nm


# chop up file at <path> into chunks of <chunk_size> number of lines
def get_chunks(path, chunk_size):
    num_lines_process = subprocess.run(['wc', '-l', path], capture_output=True, text=True)
    num_lines = int(num_lines_process.stdout.split(' ')[0])
    num_chunks = (num_lines // chunk_size) + 1
    chunks = [list(range(i*1000, (i+1)*1000)) for i in range(num_chunks)]
    chunks[-1] = list(range(chunks[-1][0], num_lines))
    return chunks

    

if __name__ == '__main__':
    # options
    chunk_size = 50000
    in_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/networks/heart_adipocyte_irf-loop-mtx_no_correlated_data_top_0.01_pct_correlates_readded.tsv'
    out_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/networks'
    chunks = get_chunks(in_path, chunk_size)


    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    for i, indices in enumerate(chunks):
        # distribute across ranks
        if i % size != rank: continue
        print(f'i: {i}, rank: {rank}, size: {size}')
        out_path = os.path.join(out_dir, f'heart_adipocyte_netx_irf-loop_{indices[0]}-{indices[-1]}.tsv')
        nm.get_metrics(in_path, indices, out_path)
