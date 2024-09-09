library(Seurat)
library(SeuratDisk)
library(stringr)
library(Matrix)
library(dplyr)
library(zellkonverter)


obj_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/expression_matrices/WHB-10Xv3/20240330/'
obj <- list.files(obj_path)
# change below for the neuron matrix after non-neuron is completed
obj <- obj[2]

out_irf_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/irf-loop_mtx/"

cluster_anno <- read.csv("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/metadata/WHB-taxonomy/20240330/cluster_annotation_term.csv",header = TRUE)
#clust_to_clust <- read.csv("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/metadata/WHB-taxonomy/20240330/cluster_to_cluster_annotation_membership.csv",header = TRUE)
cluster <- read.csv("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/metadata/WHB-taxonomy/20240330/cluster.csv",header = TRUE)

cell_meta <- read.csv("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/metadata/WHB-10Xv3/20240330/cell_metadata.csv",header = TRUE)
#regions <- read.csv("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/metadata/WHB-10Xv3/20240330/region_of_interest_structure_map.csv",header = TRUE)
#anatom <- read.csv("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas/metadata/WHB-10Xv3/20240330/anatomical_division_structure_map.csv",header = TRUE)

# summarize cell count by region
#cell_meta %>% group_by(anatomical_division_label) %>% dplyr::summarize(total = n()) %>% data.frame()

# summarize cell count by cell type
#cell_meta %>% group_by(description) %>% dplyr::summarize(total = n()) %>% data.frame()

# summarize cell count by region/cell type
#cell_meta %>% group_by(anatomical_division_label,description) %>% dplyr::summarize(total = n()) %>% data.frame()

# join cell meta to cluster information
cell_meta$row_order <- 1:nrow(cell_meta)
cell_meta <- merge(cell_meta,cluster,by = "cluster_alias",all.x = TRUE) 
cell_meta <- merge(cell_meta,cluster_anno,by = "label",all.x = TRUE)
cell_meta <- cell_meta[order(cell_meta$row_order),]
cell_meta$description <- gsub(" \\(subcluster [0-9]+\\)","",cell_meta$description)

obj <- paste0(obj_path,obj)
obj <- Read10X(obj)
obj <- CreateSeuratObject(obj)

#cell_meta <- cell_meta[cell_meta$cell_label %in% colnames(obj),]
cell_meta <- cell_meta[cell_meta$feature_matrix_label == "WHB-10Xv3-Nonneurons",]

if(all(colnames(obj) == cell_meta$cell_label)) {
    obj_matrix <- GetAssayData(object = obj)
}

for(region in unique(cell_meta$anatomical_division_label)) {
    for(cell_type in unique(cell_meta$description)) {
	#skip <- FALSE
	print(paste0("working on region: ",region," / cell type: ",cell_type))
	subset <- cell_meta$anatomical_division_label == region & cell_meta$description == cell_type
        tryCatch({
            matrix_subset <- t(as.data.frame(as.matrix(obj_matrix[,subset])))
	    region_name <- gsub(" ","-",tolower(region))
	    cell_type_name <- gsub("\\(|\\)","",gsub(" ","-",tolower(cell_type)))
	    out_path <- paste0(out_irf_dir,region_name,"_",cell_type_name,"_irf-loop-mtx.tsv")
	    if(nrow(matrix_subset) >= 60) {
	        write.table(matrix_subset,file = out_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
	    } else {
		write.table(matrix_subset,file = out_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
                print(paste0("WARNING not enough data for region: ",region_name," / cell type: ",cell_type_name," with ",nrow(matrix_subset)," cells"))
	    }},
            error = function(e) { 
		#skip <<- TRUE
	        print(paste0("WARNING matrix too large to save with ",sum(subset)," cells")) 
	    }
	)
	#if(skip) { next }
    }
}

