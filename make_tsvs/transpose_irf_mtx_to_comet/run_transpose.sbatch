#!/bin/bash

#SBATCH -A SYB111
#SBATCH -t 30:00 
#SBATCH -N 175
#SBATCH -J irf-mtx_to_comet-mtx
#SBATCH -o ./logs/irf-mtx_to_comet-mtx.%J.out
#SBATCH -e ./logs/irf-mtx_to_comet-mtx.%J.err
#SBATCH -p batch


module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main


script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/transpose_irf_mtx_to_comet/transpose_irf_mtx_to_comet.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"


srun -n 699 -c 8 --ntasks-per-node=4 python $script -d $data_dir
