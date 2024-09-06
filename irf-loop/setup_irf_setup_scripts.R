runs <- list.files("/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc")
runs <- runs[runs != "README.txt"]
matrices <- list.files("/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/ct_no_correlated_data/")

script <- c(
    "#!/bin/bash",
    "#SBATCH -A SYB111",
    "#SBATCH -J setup_irf",
    "#SBATCH -N 1",
    "#SBATCH -t 48:00:00",
    "#SBATCH -o /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/logs/setup_irf.%J.out",
    "#SBATCH -e /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/logs/setup_irf.%J.err",
    "\n",
    "source activate python_andes",
    "\n",
	do.call("c",lapply(1:length(runs),function(x){
		paste0("cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/",runs[x]
			,"\n\n",
			"python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/ct_no_correlated_data/",matrices[x]," --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111","\n"
		)
	}))
)

write.table(script,file = "/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/setup-irf-setup-scripts.slurm",sep = "\n",row.names = FALSE,col.names = FALSE,quote = FALSE)
