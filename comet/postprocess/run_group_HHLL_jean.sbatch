#!/bin/bash 
  
#SBATCH -A SYB111
#SBATCH -t 1:00:00
#SBATCH -J comet_to_edge
#SBATCH -o ./logs/comet_to_edge.%J.out
#SBATCH -e ./logs/comet_to_edge.%J.err
#SBATCH -N 233


module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main


py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/postprocess/group_HHLL_jean.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"

srun -n 1864 -c 7 --ntasks-per-node=8 python $py_script -d $data_dir --run
