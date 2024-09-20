#!/bin/bash 
  
#SBATCH -A SYB114
#SBATCH -t 02:00:00
#SBATCH -J postprocess
#SBATCH -o ./logs/bin_to_txt_all.%J.out
#SBATCH -e ./logs/bin_to_txt_all.%J.err
#SBATCH -N 274

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_mpi4py

py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/postprocess/bin_to_txt_mpi.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"
comet_tools_dir="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/genomics_gpu/tools"
num_way="2"

srun -n 2192 -c 7 --ntasks-per-node=8 python $py_script -d $data_dir \
     -c $comet_tools_dir -w $num_way --run
