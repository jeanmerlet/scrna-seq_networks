#!/bin/bash

#SBATCH -A SYB112
#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH -J ct_split_nuc-acc
#SBATCH -o logs/ct_split_nuc-acc.%J.out
#SBATCH -e logs/ct_split_nuc-acc.%J.err

source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

path="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/split_by_celltype_nucleus-acc.R"

srun -n 1 Rscript $path
