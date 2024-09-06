#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J leaf_hydathode_preprocess
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/preprocess/arabidopsis/leaf/logs/leaf_hydathode_preprocess.%j.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/preprocess/arabidopsis/leaf/logs/leaf_hydathode_preprocess.%j.err
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --mem=0


unset SLURM_EXPORT_ENV


# Filepath to the preprocessing data
INPUT_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/leaf/healthy/chory/irf-loop_mtx/leaf_hydathode_irf-loop-mtx.tsv"
PROJECT_DIRECTORY_PATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/leaf/healthy/chory/networks/leaf_hydathode"
OUTPUT_FILE_PATH="${PROJECT_DIRECTORY_PATH}/preprocessed.tsv"


source /lustre/orion/syb111/world-shared/frontier_hack/source_env


cd $PROJECT_DIRECTORY_PATH


SECONDS=0
srun python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/preprocess.py --infile $INPUT_FILEPATH --has_indices --remove_low_variance --remove_high_corr --corr_thresh 0.95 --outfile $OUTPUT_FILE_PATH
echo $SECONDS elapsed
