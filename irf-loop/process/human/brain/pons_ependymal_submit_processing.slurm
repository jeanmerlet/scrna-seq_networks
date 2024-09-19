#!/bin/bash
#SBATCH -A SYB114
#SBATCH -J pons_ependymal_process
#SBATCH -o /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/process/human/brain/logs/pons_ependymal_process.%j.out
#SBATCH -e /lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/process/human/brain/logs/pons_ependymal_process.%j.err
#SBATCH -t 2:00:00
#SBATCH -p batch
#SBATCH -N 2000
#SBATCH --threads-per-core=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=28


unset SLURM_EXPORT_ENV


module load amd-mixed/5.6.0
module load boost
module load ums/default
module load cray-mpich/8.1.28
export MPI4PY_RC_RECV_MPROBE=0
export OMP_NUM_THREADS=2
export MPICH_OFI_NIC_POLICY="ROUND-ROBIN"


# Filepath to the preprocessed data
INPUT_FILEPATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/networks/unimputed/pons_ependymal/preprocessed.tsv"
PROJECT_DIRECTORY_PATH="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/networks/unimputed/pons_ependymal"
OUTPUT_FILE_PATH="${PROJECT_DIRECTORY_PATH}/processed.tsv"


source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh
source activate /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_loop_env


cd $PROJECT_DIRECTORY_PATH


SECONDS=0
srun -N 2000 --ntasks-per-node=2 --cpus-per-task=28 --threads-per-core=2 python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/process.py --infile $INPUT_FILEPATH --outfile $OUTPUT_FILE_PATH --n_estimators 100 --max_depth 75 --header_row_idx 0 --n_jobs 28
echo $SECONDS elapsed


