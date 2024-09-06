
code_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/postprocess/human/blood/"
log_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/postprocess/human/blood/logs/"
if(!dir.exists(log_dir)) {
    dir.create(log_dir)
}
tissue_folder <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/"
tissue_network_folder <- paste0(tissue_folder,"networks/")
if(!dir.exists(tissue_network_folder)) {
    dir.create(tissue_network_folder)
}
cell_types <- list.files(tissue_network_folder)
rep_maps <- list.files(paste0(tissue_folder,"irf-loop_mtx/"))
rep_maps <- rep_maps[grepl("variance_nonrep_to_rep_map.tsv",rep_maps)]

lapply(1:length(cell_types),function(cell_type_index) {
    cell_type <- cell_types[cell_type_index]
    print(paste0("cell type: ",cell_type))
    cell_type_dir <- paste0(tissue_network_folder,cell_type,"/") 
    if(!dir.exists(cell_type_dir)) {
        dir.create(cell_type_dir)
    }
    rep_map <- rep_maps[grepl(cell_type,rep_maps)]
    script <- c(
	"#!/bin/bash",
	"#SBATCH -A SYB111",
	paste0("#SBATCH -J ",cell_type,"_postprocess"),
	paste0("#SBATCH -o ",log_dir,cell_type,"_postprocess.%j.out"),
	paste0("#SBATCH -e ",log_dir,cell_type,"_postprocess.%j.err"),
	"#SBATCH -t 2:00:00",
	"#SBATCH -p batch",
	"#SBATCH -N 2",
	"#SBATCH --mem=0",
	"\n",
	"unset SLURM_EXPORT_ENV",
	"\n",
	"# Filepath to the processed data",
	paste0('INPUT_FILEPATH="',cell_type_dir,'processed.tsv"'),
	paste0('PROJECT_DIRECTORY_PATH="',substr(cell_type_dir,1,nchar(cell_type_dir)-1),'"'),
	paste0('REPMAP_FILEPATH="',tissue_folder,"irf-loop_mtx/",rep_map,'"'),
	paste0('iRF_CODE="python /lustre/orion/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/postprocess.py"'),
	"threshold=0.01",
	"\n",
	"source /lustre/orion/syb111/proj-shared/Tools/andes/load_anaconda.sh",
        "conda activate /lustre/orion/syb111/proj-shared/Tools/irf_hp/andes_mpi",
	"\n",
	"cd $PROJECT_DIRECTORY_PATH",
	"\n",
	"srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --make_undirected --verbose &",
	"srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --make_undirected --verbose &",
	"srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --weigh_edges_by_acc --verbose &",
	"srun -N1 -n1 -c21 $iRF_CODE --infile $INPUT_FILEPATH --threshold $threshold --rep_map_path $REPMAP_FILEPATH --verbose &",
	"wait"
    )
    write.table(
        script,
	file = paste0(code_dir,cell_type,"_submit_postprocessing.slurm"),
	sep = "\n",
	row.names = FALSE,
	col.names = FALSE,
	quote = FALSE
    )
})
