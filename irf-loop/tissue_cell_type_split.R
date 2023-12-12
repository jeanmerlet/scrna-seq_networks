
library(Seurat)
library(stringr)
library(dplyr)
library(optparse)

option_list <- list(
	make_option(
		c("-d","--dirin"), 
		type = "character",
		default = NULL, 
        help = "directory for tissue", 
        metavar = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# folder_names <- list.files(directory)
# broad cell type in gtex non-broad tissue specific cell type

split_by_celltype <- function(obj,out_prefix) {
    split_obj <- SplitObject(obj,split.by = 'Granular.cell.type')
    # split_obj <- SplitObject(obj,split.by = 'final_cell_type')
    # ct <- names(split_obj)
    # ct <- gsub('\\(', '', gsub('\\)', '', gsub(' ', '-', tolower(ct))))
    # ct <- paste0("/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/",ct)
    # lapply(1:length(ct),function(x){
    #     if (!dir.exists(ct[x])) {
    #         print(paste0("creating directory for: ",ct[x]))
    #         dir.create(ct[x])
    #     }
    # })
    for (i in seq_along(split_obj)) {
        obj <- split_obj[[i]]
        celltype <- gsub('\\(|,|[.]', '', gsub('\\)', '', gsub(' ', '-', tolower(names(split_obj)[i]))))
	celltype <- gsub("/","-",celltype)
        print(paste0("Working on cell type: ",celltype))
        out_path <- paste0(out_prefix, '_', celltype, '_irf-loop-mtx.tsv')
        count_matrix <- as.data.frame(GetAssayData(object = obj))
        count_matrix <- t(count_matrix)
        vars <- data.frame(do.call("rbind",lapply(1:ncol(count_matrix),function(z){
                c(colnames(count_matrix)[z],var(count_matrix[,z]))
        })))
        names(vars) <- c("Gene","Var")
        vars$Var <- as.numeric(vars$Var)
        vars <- vars[order(vars$Var,decreasing = TRUE),]
        vars <- vars[vars$Var > 0.0000000001,]
        vars_keep <- vars$Gene
        count_matrix <- count_matrix[,which(colnames(count_matrix) %in% vars_keep)]
        names(count_matrix) <- c("",names(count_matrix))
        write.table(count_matrix,out_path,sep = "\t",col.names = NA,row.names = TRUE,quote = FALSE)
    }
}

obj_path <- opt$dirin
out_dir <- paste0(dirname(dirname(obj_path)),'/irf-loop_mtx')
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}
obj <- readRDS(obj_path)
filename <- basename(obj_path)
# tissue <- str_match(filename,'([a-zA-Z-]+)_imputed.rds')[1, 2]
tissue <- gsub("_imputed.rds","",filename)
print(paste0('Splitting ',tissue,'...'))
# print("adjusting metadata beforehand...")

# adjusting the metadata
# fct <- read.table("/gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/meta/cell_meta_imputed_final_labels.tsv",sep = "\t",header = TRUE)
# obj@meta.data <- fct
# obj@meta.data$final_cell_type_broad <- gsub("BA8_|BA9_|BA10_","",obj@meta.data$final_cell_type)

# run the script
out_prefix <- paste0(out_dir,'/',tissue)
split_by_celltype(obj,out_prefix)

# meta <- obj@meta.data
# meta$row_order <- 1:nrow(meta)
# rownames <- rownames(meta)
# meta$seq_id <- rownames
# meta_new <- merge(
#     meta,
#     fct[,c("seq_id","spec_cell_types")],
#     by = "seq_id",
#     all.x = TRUE
# )
# meta_new <- meta_new[,which(colnames(meta_new) != "seq_id")]
# meta_new <- meta_new[order(meta_new$row_order,decreasing = FALSE),]
# rownames(meta_new) <- rownames
# meta_new <- meta_new[,which(colnames(meta_new) != "row_order")]
# all(colnames(obj) == rownames(meta_new))
# obj@meta.data <- meta_new
# cts <- do.call("c",lapply(1:nrow(obj@meta.data),function(x){
#     current_row <- obj@meta.data[x,]
#     if(current_row$sub_clust == "Oligo_mat") {
#         paste0(current_row$Brain.Regions.,"_",current_row$spec_cell_types,"_",current_row$sub_clust)
#     } else {
#         if(current_row$sub_clust == "PN_dev") {
#             paste0(current_row$Brain.Regions.,"_",current_row$major_clust)
#             } else {
#                 paste0(current_row$Brain.Regions.,"_",current_row$spec_cell_types,"_",current_row$major_clust)        
#             }
#     }
# }))
# cts <- gsub("Glial/Non-Neuron","glial-non-neuron",cts)
# cts <- gsub("Excitatory Principal Neurons","excitatory-principal-neurons",cts)
# cts <- gsub("MGE Inhibitory Interneurons","mge-inhibitory-interneurons",cts)
# cts <- gsub("CGE Inhibitory Interneurons","cge-inhibitory-interneurons",cts)
# obj <- AddMetaData(obj,cts,col.name = "final_cell_type")
