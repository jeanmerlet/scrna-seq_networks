#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J tissue_ct_split_gtex
#SBATCH -N 8
#SBATCH \-\-mem=0
#SBATCH -t 48:00:00
#SBATCH -o /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/logs/tissue_ct_split_gtex.%J.out
#SBATCH -e /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/logs/tissue_ct_split_gtex.%J.err

source activate /gpfs/alpine/syb105/proj-shared/Personal/atown/Libraries/Andes/Anaconda3/envs/deseq2_andes

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/breast/healthy/gtex/seurat/breast_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/esophagus_mucosa/healthy/gtex/seurat/esophagus_mucosa_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/esophagus_muscularis/healthy/gtex/seurat/esophagus_muscularis_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/heart/healthy/gtex/seurat/heart_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/lung/healthy/gtex/seurat/lung_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/prostate/healthy/gtex/seurat/prostate_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/skeletal_muscle/healthy/gtex/seurat/skeletal_muscle_imputed.rds" &

srun -n1 -N1 -c16 --cpu-bind=threads --exclusive Rscript tissue_cell_type_split.R -d "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/skin/healthy/gtex/seurat/skin_imputed.rds" &
wait
