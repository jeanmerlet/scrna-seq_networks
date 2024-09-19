#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J cerebral-cortex_vascular-smooth-muscle-cell_preprocess
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/preprocess/human/brain/logs/cerebral-cortex_vascular-smooth-muscle-cell_preprocess.%j.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/preprocess/human/brain/logs/cerebral-cortex_vascular-smooth-muscle-cell_preprocess.%j.err
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --mem=0


unset SLURM_EXPORT_ENV


# Filepath to the preprocessing data
INPUT_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/irf-loop_mtx/unimputed/cerebral-cortex_vascular-smooth-muscle-cell_irf-loop-mtx.tsv"
PROJECT_DIRECTORY_PATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/networks/unimputed/cerebral-cortex_vascular-smooth-muscle-cell"
OUTPUT_FILE_PATH="${PROJECT_DIRECTORY_PATH}/preprocessed.tsv"


source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_loop_env_andes


cd $PROJECT_DIRECTORY_PATH


SECONDS=0
srun python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/preprocess.py --infile $INPUT_FILEPATH --has_indices --remove_low_variance --remove_high_corr --corr_thresh 0.95 --outfile $OUTPUT_FILE_PATH
echo $SECONDS elapsed
