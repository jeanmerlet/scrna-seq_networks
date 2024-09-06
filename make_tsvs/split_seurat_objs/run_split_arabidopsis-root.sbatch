#!/bin/bash

#SBATCH -A SYB112
#SBATCH -N 1
#SBATCH -t 8:00:00
#SBATCH -J ct_split_ara-root
#SBATCH -o logs/ct_split_ara-root.%J.out
#SBATCH -e logs/ct_split_ara-root.%J.err

source activate /ccs/home/atown/Apps/Andes/anaconda3/envs/deseq2_andes

path="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/split_by_celltype_arabidopsis-root.R"

srun -n 1 Rscript $path
