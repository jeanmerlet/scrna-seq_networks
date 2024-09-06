library(Seurat)
library(SeuratDisk)
library(stringr)
library(Matrix)
library(dplyr)


obj_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/raw/synapse_download/'
cell_types <- list.files(obj_path)
cell_types <- cell_types[cell_types != "cache"]
cell_type_names <- gsub("_","-",cell_types)
tissue <- "blood"
out_irf_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/irf-loop_mtx/"


# read in borad cell type objects and split further into sub-types 
# note: there are too many rows per broad cell type to run irf-loop
# note: the individual broad cell type objects are alread normalized
# note: too many cells to run imputation across all matrices


objects <- lapply(1:length(cell_types),function(x){
    print(paste0("working on: ",cell_types[x]))
    if(cell_types[x] %in% c("cd4_t_cells","conventional_cd8_t_cells")) {
        meta <- paste0(obj_path,cell_types[x],"/",gsub("_t_cells","",cell_types[x]),"_metadata.csv")
        cell_type_path <- paste0(obj_path,cell_types[x],"/",gsub("_t_cells","",cell_types[x]),"_rna.rds")
    } else if (cell_types[x] == "cd4_t_helper_memory_cells") {
        meta <- paste0(obj_path,cell_types[x],"/",gsub("t_helper_memory_cells","helper_memory",cell_types[x]),"_metadata.csv")
        cell_type_path <- paste0(obj_path,cell_types[x],"/",gsub("t_helper_memory_cells","helper_memory",cell_types[x]),"_rna.rds")
    } else {
	meta <- paste0(obj_path,cell_types[x],"/",cell_types[x],"_metadata.csv")
        cell_type_path <- paste0(obj_path,cell_types[x],"/",cell_types[x],"_rna.rds")
    }
    obj <- readRDS(cell_type_path)
    meta <- read.csv(meta,row.names = 'X')
    #if("Cluster_names" %in% names(meta)) {
    #    meta <- meta %>% group_by(Cluster_names) %>% summarize(total = n())
    #    print(meta)
    #}
    obj <- CreateSeuratObject(counts = obj,meta.data = meta)
    if(cell_types[x] %in% c("mait_cells","progenitor_cells")) {
	y_file <- cell_type_names[x]
        out_irf_path <- paste0(out_irf_dir,tissue,'_',y_file,'_irf-loop-mtx.tsv')
        write.table(t(as.data.frame(as.matrix(GetAssayData(object = obj)))),file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
    } else {
	subs <- unique(meta$Cluster_names)
        lapply(1:length(subs),function(y){
            y_file <- paste0(cell_type_names[x],"-",tolower(gsub("/","-",gsub("_","-",gsub(" ","-",subs[y])))))
            obj_sub <- subset(obj,subset = Cluster_names == subs[y])
            out_irf_path <- paste0(out_irf_dir,tissue,'_',y_file,'_irf-loop-mtx.tsv') 
            tryCatch({
	        write.table(as.data.frame(as.matrix(t(GetAssayData(object = obj_sub)))),file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
	    },error = function(e) return(NULL))
        })
    }
})


