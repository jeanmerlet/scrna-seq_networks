#!/bin/bash 
  
#SBATCH -A SYB114
#SBATCH -t 00:30:00
#SBATCH -J combine_tsvs
#SBATCH -o ./logs/combine_tsvs.%J.out
#SBATCH -e ./logs/combine_tsvs.%J.err
#SBATCH -N 274

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_mpi4py

py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/postprocess/combine_tsvs.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"

srun -n 274 -c 8 --ntasks-per-node=1 python $py_script -d $data_dir --run
