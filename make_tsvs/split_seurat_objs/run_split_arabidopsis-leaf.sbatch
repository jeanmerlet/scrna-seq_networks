#!/bin/bash

#SBATCH -A SYB111
#SBATCH -N 1
#SBATCH -t 8:00:00
#SBATCH -J ct_split_ara-leaf
#SBATCH -o logs/ct_split_ara-leaf.%J.out
#SBATCH -e logs/ct_split_ara-leaf.%J.err

#source activate /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/seurat
source activate /ccs/home/atown/Apps/Andes/anaconda3/envs/deseq2_andes

path="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/split_by_celltype_arabidopsis-leaf.R"

srun -n 1 Rscript $path
