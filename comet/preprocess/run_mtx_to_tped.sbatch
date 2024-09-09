#!/bin/bash

#SBATCH -A SYB114
#SBATCH -t 1:00:00 
#SBATCH -N 137
#SBATCH -J mtx_to_tped
#SBATCH -o ./logs/comet_mtx_to_tped.%J.out
#SBATCH -e ./logs/comet_mtx_to_tped.%J.err
#SBATCH -p batch


source activate frontier_mpi4py


py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/preprocess/mtx_to_tped.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"


srun -n 274 -c 8 --ntasks-per-node=2 python $py_script -d $data_dir --run
