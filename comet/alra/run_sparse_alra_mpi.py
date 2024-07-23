from mpi4py import MPI
import subprocess


def run_alra_on_gene_set(mtx):
    # Rscript sparse_alra passing it the gene_indices




comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, pair in enumerate(paired_fastq_list):
    if i % size != rank: continue

