#!/bin/bash

#SBATCH -A SYB112
#SBATCH -N 1
#SBATCH -t 8:00:00
#SBATCH -J ct_split_kidney
#SBATCH -o logs/ct_split_kidney.%J.out
#SBATCH -e logs/ct_split_kidney.%J.err
#SBATCH -p gpu

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/split_by_celltype_kidney-kpmp.R"

srun -n 1 Rscript $path
