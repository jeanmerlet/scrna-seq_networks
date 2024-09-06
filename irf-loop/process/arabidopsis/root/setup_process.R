
code_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/process/arabidopsis/root/"
log_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/process/arabidopsis/root/logs/"
if(!dir.exists(log_dir)) {
    dir.create(log_dir)
}
tissue_folder <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/root/healthy/ohler/"
tissue_network_folder <- paste0(tissue_folder,"networks/")
if(!dir.exists(tissue_network_folder)) {
    dir.create(tissue_network_folder)
}
irf_matrices <- list.files(tissue_network_folder)
#irf_matrices <- irf_matrices[grepl("_irf-loop-mtx.tsv",irf_matrices)]

lapply(1:length(irf_matrices),function(cell_type) {
    irf_matrix_cell_type <- paste0(tissue_network_folder,irf_matrices[cell_type],"/preprocessed.tsv")
    cell_type <- irf_matrices[cell_type]
    print(paste0("cell type: ",cell_type))
    cell_type_dir <- paste0(tissue_network_folder,cell_type,"/")
    if(!dir.exists(cell_type_dir)) {
        dir.create(cell_type_dir)
    }
    script <- c(
	"#!/bin/bash",
	"#SBATCH -A SYB111",
	paste0("#SBATCH -J ",cell_type,"_process"),
	paste0("#SBATCH -o ",log_dir,cell_type,"_process.%j.out"),
	paste0("#SBATCH -e ",log_dir,cell_type,"_process.%j.err"),
	"#SBATCH -t 12:00:00",
	"#SBATCH -p batch",
	"#SBATCH -N 200",
	"#SBATCH --threads-per-core=2",
	"\n",
	"unset SLURM_EXPORT_ENV",
	"\n",
	"module load amd-mixed/5.6.0",
	"module load boost",
	"module load ums/default",
	"module load cray-mpich/8.1.28",
	"export MPI4PY_RC_RECV_MPROBE=0",
	"export OMP_NUM_THREADS=28",
	'export MPICH_OFI_NIC_POLICY="ROUND-ROBIN"',
	"\n",
	"# Filepath to the preprocessed data",
	paste0('INPUT_FILEPATH="',irf_matrix_cell_type,'"'),
	paste0('PROJECT_DIRECTORY_PATH="',substr(cell_type_dir,1,nchar(cell_type_dir)-1),'"'),
	'OUTPUT_FILE_PATH="${PROJECT_DIRECTORY_PATH}/processed.tsv"',
	"\n",
	"source /lustre/orion/syb111/proj-shared/Tools/frontier/load_anaconda.sh",
	"source activate /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_loop_env",
	"\n",
	"cd $PROJECT_DIRECTORY_PATH",
	"\n",
	"SECONDS=0",
	"srun -N 200 --ntasks-per-node=4 --cpus-per-task=14 --threads-per-core=2 python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/process.py --infile $INPUT_FILEPATH --outfile $OUTPUT_FILE_PATH --n_estimators 100 --max_depth 75 --header_row_idx 0 --n_jobs 28",
	"echo $SECONDS elapsed",
	"\n"
    )
    write.table(
        script,
	file = paste0(code_dir,cell_type,"_submit_processing.slurm"),
	sep = "\n",
	row.names = FALSE,
	col.names = FALSE,
	quote = FALSE
    )
})
