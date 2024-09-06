#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J leaf_unknown_postprocess
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/postprocess/arabidopsis/leaf/logs/leaf_unknown_postprocess.%j.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/postprocess/arabidopsis/leaf/logs/leaf_unknown_postprocess.%j.err
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 2
#SBATCH --mem=0


unset SLURM_EXPORT_ENV


# Filepath to the processed data
INPUT_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/leaf/healthy/chory/networks/leaf_unknown/processed.tsv"
PROJECT_DIRECTORY_PATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/leaf/healthy/chory/networks/leaf_unknown"
REPMAP_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/leaf/healthy/chory/irf-loop_mtx/leaf_unknown_irf-loop-mtx_ge0.01variance_nonrep_to_rep_map.tsv"
iRF_CODE="python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/postprocess.py"
threshold=0.01


source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /lustre/orion/syb111/proj-shared/Tools/irf_hp/andes_mpi


cd $PROJECT_DIRECTORY_PATH


srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --verbose &
wait
