#!/bin/bash

#SBATCH -A SYB111
#SBATCH -t 00:30:00 
#SBATCH -N 117
#SBATCH -J tped_to_bin
#SBATCH -o ./logs/tped_to_bin.%J.out
#SBATCH -e ./logs/tped_to_bin.%J.err
#SBATCH -p batch


module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main


py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/preprocess/tped_to_bin.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"
comet_tools_dir="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/genomics_gpu/tools"


srun -n 234 -c 8 --ntasks-per-node=2 python $py_script -d $data_dir -c $comet_tools_dir --run
