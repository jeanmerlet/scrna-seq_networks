#!/bin/bash

#SBATCH -A SYB111
#SBATCH -t 1:00:00 
#SBATCH -N 1
#SBATCH -J mtx_to_tped
#SBATCH -o ./logs/comet_mtx_to_tped.%J.out
#SBATCH -e ./logs/comet_mtx_to_tped.%J.err
#SBATCH -p batch


module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main


py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/create_tped/mtx_to_tped.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"


srun -n 2 -c 8 --ntasks-per-node=2 python $py_script -d $data_dir --run
