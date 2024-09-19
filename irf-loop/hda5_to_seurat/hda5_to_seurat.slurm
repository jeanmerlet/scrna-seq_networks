#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J hda5_to_seurat
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH \-\-mem=0
#SBATCH -t 48:00:00
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/logs/hda5_to_seurat.%J.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/logs/hda5_to_seurat.%J.err

echo "loading conda env"

source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /lustre/orion/syb111/proj-shared/Tools/andes/anaconda3/envs/python  

echo "running R script"

python hda5_to_seurat.py /lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/expression_matrices/WHB-10Xv3/20240330/WHB-10Xv3-Neurons-raw.h5ad /lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/expression_matrices/WHB-10Xv3/20240330/neurons_mat 

echo "python script complete.. beep bop beep"

