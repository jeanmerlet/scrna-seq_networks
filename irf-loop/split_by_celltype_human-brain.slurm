#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J split_brain
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH \-\-mem=0
#SBATCH -t 48:00:00
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/logs/brain_split.%J.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/logs/brain_split.%J.err

echo "loading conda env"

source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /ccs/home/atown/Apps/Andes/anaconda3/envs/deseq2_andes 

echo "running R script"

Rscript split_by_celltype_human-brain.R 

echo "R script complete.. beep bop beep"

