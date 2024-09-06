#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N 1
#SBATCH -t 2:00:00
#SBATCH -J ct_split_heart
#SBATCH -o logs/ct_split_heart.%J.out
#SBATCH -e logs/ct_split_heart.%J.err

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/split_by_celltype_heart-sanger.R"

srun -n 1 Rscript $path
