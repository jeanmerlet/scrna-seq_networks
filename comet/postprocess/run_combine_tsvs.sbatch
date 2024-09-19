#!/bin/bash 
  
#SBATCH -A SYB111
#SBATCH -t 30:00
#SBATCH -J combine_tsvs
#SBATCH -o ./logs/combine_tsvs.%J.out
#SBATCH -e ./logs/combine_tsvs.%J.err
#SBATCH -N 260


module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main


py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/postprocess/combine_tsvs.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"

srun -n 260 -c 8 --ntasks-per-node=1 python $py_script -d $data_dir --run
