#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J split_blood
#SBATCH -N 8
#SBATCH \-\-mem=0
#SBATCH -t 6:00:00
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/logs/blood_split.%J.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/logs/blood_split.%J.err

echo "loading conda env"

source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /ccs/home/atown/Apps/Andes/anaconda3/envs/deseq2_andes 

echo "running R script"

mpiexec --mca btl self,tcp -np 8 --map-by ppr:1:node Rscript split_by_celltype_human-blood.R 

echo "R script complete.. beep bop beep"

