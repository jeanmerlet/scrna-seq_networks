#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J blood_b-cells-switched-memory_postprocess
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/postprocess/human/blood/logs/blood_b-cells-switched-memory_postprocess.%j.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/postprocess/human/blood/logs/blood_b-cells-switched-memory_postprocess.%j.err
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 8
#SBATCH --mem=0


unset SLURM_EXPORT_ENV


# Filepath to the processed data
INPUT_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/networks/unimputed/blood_b-cells-switched-memory/processed.tsv"
PROJECT_DIRECTORY_PATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/networks/unimputed/blood_b-cells-switched-memory"
REPMAP_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/irf-loop_mtx/unimputed/blood_b-cells-switched-memory_irf-loop-mtx_ge0.01variance_nonrep_to_rep_map.tsv"
iRF_CODE="python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/postprocess.py"


source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh
conda activate /lustre/orion/syb111/proj-shared/Tools/irf_hp/andes_mpi


cd $PROJECT_DIRECTORY_PATH


srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.005 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.005 --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.005 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.005 --rep_map_path $REPMAP_FILEPATH --verbose &


srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.05 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.05 --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.05 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.05 --rep_map_path $REPMAP_FILEPATH --verbose &


srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.1 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.1 --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.1 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.1 --rep_map_path $REPMAP_FILEPATH --verbose &


srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.2 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.2 --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.2 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.2 --rep_map_path $REPMAP_FILEPATH --verbose &


srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.25 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.25 --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.25 --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &
srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold 0.25 --rep_map_path $REPMAP_FILEPATH --verbose &
wait
