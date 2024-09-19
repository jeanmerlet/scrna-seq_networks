#!/bin/bash

#SBATCH -A SYB114
#SBATCH -t 01:00:00 
#SBATCH -N 137
#SBATCH -J tped_to_bin
#SBATCH -o ./logs/tped_to_bin.%J.out
#SBATCH -e ./logs/tped_to_bin.%J.err


source activate frontier_mpi4py


py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/preprocess/tped_to_bin.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"
comet_tools_dir="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/genomics_gpu/tools"


srun -n 274 --ntasks-per-node=2 python $py_script -d $data_dir -c $comet_tools_dir --run
