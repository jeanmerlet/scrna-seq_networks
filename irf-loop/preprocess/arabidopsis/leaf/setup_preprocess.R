
code_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/preprocess/arabidopsis/leaf/"
log_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/preprocess/arabidopsis/leaf/logs/"
if(!dir.exists(log_dir)) {
    dir.create(log_dir)
}
tissue_folder <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/arabidopsis/leaf/healthy/chory/"
tissue_network_folder <- paste0(tissue_folder,"networks/")
if(!dir.exists(tissue_network_folder)) {
    dir.create(tissue_network_folder)
}
irf_matrices <- list.files(paste0(tissue_folder,"irf-loop_mtx/"))
irf_matrices <- irf_matrices[grepl("_irf-loop-mtx.tsv",irf_matrices)]

lapply(1:length(irf_matrices),function(cell_type) {
    irf_matrix_cell_type <- paste0(tissue_folder,"irf-loop_mtx/",irf_matrices[cell_type])
    cell_type <- gsub("_irf-loop-mtx.tsv","",irf_matrices[cell_type])
    print(paste0("cell type: ",cell_type))
    cell_type_dir <- paste0(tissue_network_folder,cell_type,"/")
    if(!dir.exists(cell_type_dir)) {
        dir.create(cell_type_dir)
    }
    script <- c(
	"#!/bin/bash",
	"#SBATCH -A SYB111",
	paste0("#SBATCH -J ",cell_type,"_preprocess"),
	paste0("#SBATCH -o ",log_dir,cell_type,"_preprocess.%j.out"),
	paste0("#SBATCH -e ",log_dir,cell_type,"_preprocess.%j.err"),
	"#SBATCH -t 2:00:00",
	"#SBATCH -p batch",
	"#SBATCH -N 1",
	"#SBATCH --mem=0",
	"\n",
	"unset SLURM_EXPORT_ENV",
	"\n",
	"# Filepath to the preprocessing data",
	paste0('INPUT_FILEPATH="',irf_matrix_cell_type,'"'),
	paste0('PROJECT_DIRECTORY_PATH="',substr(cell_type_dir,1,nchar(cell_type_dir)-1),'"'),
	'OUTPUT_FILE_PATH="${PROJECT_DIRECTORY_PATH}/preprocessed.tsv"',
	"\n",
	"source /lustre/orion/syb111/world-shared/frontier_hack/source_env",
	"\n",
	"cd $PROJECT_DIRECTORY_PATH",
	"\n",
	"SECONDS=0",
	"srun python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/preprocess.py --infile $INPUT_FILEPATH --has_indices --remove_low_variance --remove_high_corr --corr_thresh 0.95 --outfile $OUTPUT_FILE_PATH",
	"echo $SECONDS elapsed"
    )
    write.table(
        script,
	file = paste0(code_dir,cell_type,"_submit_preprocessing.slurm"),
	sep = "\n",
	row.names = FALSE,
	col.names = FALSE,
	quote = FALSE
    )
})
