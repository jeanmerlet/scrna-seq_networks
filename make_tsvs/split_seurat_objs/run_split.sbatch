#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -J ct_split
#SBATCH -o logs/ct_split.%J.out
#SBATCH -e logs/ct_split.%J.err

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/split_by_celltype_liver.R"

srun -n 1 Rscript $path
