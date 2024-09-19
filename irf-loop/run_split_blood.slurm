#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N 64
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH -J mpi_blood
#SBATCH --ntasks-per-node=1
#SBATCH -o logs/mpi_blood.%J.out
#SBATCH -e logs/mpi_blood.%J.err

source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /ccs/home/atown/Apps/Andes/anaconda3/envs/deseq2_andes

script_path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/split_by_celltype_human-blood.R"

mpiexec -np 64 Rscript $script_path
