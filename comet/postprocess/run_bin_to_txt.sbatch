#!/bin/bash 
  
#SBATCH -A SYB111
#SBATCH -t 02:00:00
#SBATCH -J postprocess
#SBATCH -o ./logs/bin_to_txt_all.%J.out
#SBATCH -e ./logs/bin_to_txt_all.%J.err
#SBATCH -N 16

module load cray-mpich
unset MPICH_OFI_NIC_POLICY
source activate frontier_main

py_script="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/postprocess/bin_to_txt_mpi.py"
data_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human"
comet_tools_dir="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/genomics_gpu/tools"
num_way="2"

#srun -n 1864 -c 7 --ntasks-per-node=8 python $py_script -d $data_dir \
srun -n 16 -c 7 --ntasks-per-node=1 python $py_script -d $data_dir \
     -c $comet_tools_dir -w $num_way --run
