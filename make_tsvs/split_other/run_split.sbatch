#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N 17
#SBATCH -t 2:00:00
#SBATCH -J ct_split
#SBATCH -o logs/ct_split.%J.out
#SBATCH -e logs/ct_split.%J.err

module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main

#module load python

mpi_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/split_other/split_10x-mtx.py"

srun -n 17 -c 1 --ntasks-per-node=1 python $mpi_script
#srun -n 17 -c 8 python $mpi_script
