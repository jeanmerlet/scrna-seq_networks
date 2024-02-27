gtex_tissues <- c("protstate")

lapply(1:length(gtex_tissues),function(x){
	print(paste0("tissue: ",gtex_tissues[x]))
	tissue_folder <- paste0("/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/",gtex_tissues[x],"/healthy/gtex/")
	tissue_network_folder <- paste0(tissue_folder,"networks/")
	if(!dir.exists(tissue_network_folder)) {
		dir.create(tissue_network_folder)
	}
	irf_matrices <- list.files(paste0(tissue_folder,"irf-loop_mtx/"))
	irf_rep_non_rep <- irf_matrices[grepl("mtx_ge0.01variance_nonrep_to_rep_map.tsv",irf_matrices)]
	irf_matrices <- irf_matrices[grepl("_irf-loop-mtx.tsv",irf_matrices)]
	lapply(1:length(irf_matrices),function(y){
		irf_matrix_cell_type <- paste0(tissue_folder,"irf-loop_mtx/",irf_matrices[y])
		irf_rep_non_rep_cell_type <- paste0(tissue_folder,"irf-loop_mtx/",irf_rep_non_rep[y])
		cell_type <- gsub("_irf-loop-mtx.tsv","",irf_matrices[y])
		print(paste0("cell type: ",cell_type))
		cell_type_dir <- paste0(tissue_network_folder,cell_type,"/")
		submits_dir <- paste0(cell_type_dir,"submits/")
		if(!dir.exists(cell_type_dir)) {
			dir.create(cell_type_dir)
		}
		if(!dir.exists(submits_dir)) {
			dir.create(submits_dir)
		}
		script <- c(
			# all of this below needs to be changed
			"#!/bin/bash",
			"#BSUB -P SYB112",
			paste0("#BSUB -J ",cell_type,"_postprocessing"),
			paste0("#BSUB -o ",cell_type,"_postprocessing.out"),
			paste0("#BSUB -e ",cell_type,"_postprocessing.err"),
			"#BSUB -W 2:00",
			"#BSUB -nnodes 4",
			"\n",
			"source /ccs/home/atown/Scripts/loadCondaOnSummit.sh",
			"conda activate /gpfs/alpine2/syb111/proj-shared/environments/irf_loop_env",
			"\n",
			"threshold=0.01",
			paste0('infile="',cell_type_dir,"processed_data.tsv",'"'),
			paste0('repmap_path="',irf_rep_non_rep_cell_type,'"'),
			paste0('PROJECT_DIRECTORY_PATH="',substr(cell_type_dir,1,nchar(cell_type_dir)-1),'"'),
			'iRF_CODE="/gpfs/alpine2/syb111/proj-shared/Tools/irf_hp/irf_network_gen/src/postprocess.py"',
			"\n",
			"cd $PROJECT_DIRECTORY_PATH",
			"\n",
			"jsrun -n1 -a1 -c21 -g3 -dpacked python $iRF_CODE --infile $infile --delim , --threshold $threshold --rep_map_path $repmap_path --weigh_edges_by_acc --make_undirected --verbose &",
			"\n",
			"jsrun -n1 -a1 -c21 -g3 -dpacked python $iRF_CODE --infile $infile --delim , --threshold $threshold --rep_map_path $repmap_path --make_undirected --verbose &",
			"\n",
			"jsrun -n1 -a1 -c21 -g3 -dpacked python $iRF_CODE --infile $infile --delim , --threshold $threshold --rep_map_path $repmap_path --weigh_edges_by_acc --verbose &",
			"\n",
			"jsrun -n1 -a1 -c21 -g3 -dpacked python $iRF_CODE --infile $infile --delim , --threshold $threshold --rep_map_path $repmap_path --verbose &",
			"\n",
			"wait"
		)
		write.table(script,file = paste0(submits_dir,"submit_postprocessing.lsf"), sep = "\n",row.names = FALSE, col.names = FALSE,quote = FALSE)
	})
})
