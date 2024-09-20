#!/bin/bash 
  
#SBATCH -A SYB114
#SBATCH -t 1:00:00
#SBATCH -J comet_to_edge
#SBATCH -o ./logs/comet_to_edge.%J.out
#SBATCH -e ./logs/comet_to_edge.%J.err
#SBATCH -N 274

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_mpi4py

py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/postprocess/group_HHLL_jean.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"

#srun -n 1864 -c 7 --ntasks-per-node=8 python $py_script -d $data_dir --run
