
library(stringr)
library(Matrix)
library(dplyr)
library(purrr)
library(Seurat)


pbdMPI::init()


source("/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/irf-loop/impute.R")


tissue <- "blood"
obj_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/raw/synapse_download/'
cell_types <- list.files(obj_path)
cell_types <- cell_types[!(cell_types %in% c("cache","extra"))]
cell_type_names <- gsub("_","-",cell_types)


out_irf_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/irf-loop_mtx/"
out_comet_dir <- "/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/blood/healthy/broad_inst/comet_mtx/"


# read in all objects and combine into one large object
if (pbdMPI::comm.rank() == 0) {
    print("reading in matrices and joining")
}  
    mtx <- lapply(1:length(cell_types),function(x){
    #mtx <- lapply(1:2,function(x){
        print(paste0("loading: ",cell_types[x]))
        if(cell_types[x] %in% c("cd4_t_cells","conventional_cd8_t_cells")) {
            meta <- paste0(obj_path,cell_types[x],"/",gsub("_t_cells","",cell_types[x]),"_metadata.csv")
            cell_type_path <- paste0(obj_path,cell_types[x],"/",gsub("_t_cells","",cell_types[x]),"_rna.rds")
        } else {
            meta <- paste0(obj_path,cell_types[x],"/",cell_types[x],"_metadata.csv")
            cell_type_path <- paste0(obj_path,cell_types[x],"/",cell_types[x],"_rna.rds")
        }
        obj <- readRDS(cell_type_path)
        meta <- read.csv(meta,row.names = 'X')
        meta$cell_type <- cell_type_names[x]
	obj <- CreateSeuratObject(counts = obj,meta.data = meta)
        obj <- as.spam.dgCMatrix(t(GetAssay(obj,assay = "RNA")@data))
        # for Seurat v5
        #obj <- as.spam.dgCMatrix(t(obj[["RNA"]]$counts))
        #obj <- as.spam.dgCMatrix((obj[["RNA"]]$counts))
        return(list(
	    obj,
	    meta
        ))
    })
    meta_all <- bind_rows(map(mtx,2))
    mtx_ <- do.call("rbind.spam",map(mtx,1))
#} else {
#    meta_all <- NULL
#    mtx <- NULL
#}


# extract components of spam matrix then rebuild on all other ranks
#if (pbdMPI::comm.rank() == 0) {
#    values <- mtx@entries
#    row_indices <- mtx@rowpointers
#    col_indices <- mtx@colindices
#    dimensions <- dim(mtx)
#} else {
#    values <- NULL
#    row_indices <- NULL
#    col_indices <- NULL
#    dimensions <- NULL
#}


# broadcast the data
#meta_all <- pbdMPI::bcast(meta_all, rank.source = 0)
#mtx <- pbdMPI::bcast(mtx, rank.source = 0)
#values <- pbdMPI::bcast(values, rank.source = 0)
#row_indices <- pbdMPI::bcast(row_indices, rank.source = 0)
#col_indices <- pbdMPI::bcast(col_indices, rank.source = 0)
#dimensions <- pbdMPI::bcast(dimensions, rank.source = 0)


#if (pbdMPI::comm.rank() != 0) {
#    mtx <- spam(matrix(0, nrow = dimensions[1], ncol = dimensions[2]))
#    mtx@entries <- values
#    mtx@rowpointers <- row_indices
#    mtx@colindices <- col_indices
#}


# run imputation
if (pbdMPI::comm.rank() == 0) {
    num_cells <- dim(mtx)[1]
    num_genes <- dim(mtx)[2]
    sprintf('Matrix has %s cells and %s genes', num_cells, num_genes)
    k <- choose_k(mtx)
} else {
    num_cells <- NULL
    num_genes <- NULL
    k <- NULL
}
num_cells <- pbdMPI::bcast(num_cells, rank.source = 0)
num_genes <- pbdMPI::bcast(num_genes, rank.source = 0)
k <- pbdMPI::bcast(k, rank.source = 0)
if (pbdMPI::comm.rank() == 0) {
    print('Starting rsvd...')
    start <- Sys.time()
    rsvd_out <- irlba(mtx, k)
    u <- as.spam(rsvd_out$u)
    d <- as.spam(rsvd_out$d)
    v <- as.spam(rsvd_out$v)
    rm(rsvd_out)
} else {
    u <- NULL
    d <- NULL
    v <- NULL
}
u <- pbdMPI::bcast(u, rank.source = 0)
d <- pbdMPI::bcast(d, rank.source = 0)
v <- pbdMPI::bcast(v, rank.source = 0)
u_mtx <- u[,1:k]
d_mtx <- as.spam(diag(d[1:k]))
ud_mtx <- u_mtx %*% d_mtx
v_mtx <- t(v[,1:k])
all_indices <- 1:num_genes
chunk_size <- num_genes %/% (num_nodes - 1)
idx_sets <- split(all_indices, ceiling(seq_along(all_indices)/chunk_size))
if (pbdMPI::comm.rank() == 0) {
    print('Imputing...')
    start <- Sys.time()
}
job_ids <- seq(1, length(idx_sets))
gene_colsets <- pbdMPI::task.pull(job_ids,mpi_split_colsets,idx_sets,mtx,ud_mtx,v_mtx)
if (pbdMPI::comm.rank() == 0) {
    print('Imputation done')
    print(Sys.time() - start)
    print('Cbinding...')
    start <- Sys.time()
    mtx <- do.call('cbind.spam',gene_colsets)
    print('Cbinding done')
    print(Sys.time() - start)
    print(dim(mtx))
}


# split objects again with meta and save cell types
print("splitting by cell type and saving")
if( pbdMPI::comm.rank() == 0) {
    lapply(unique(meta_all$cell_type),function(cell_type){
        meta_ct <- meta_all[meta_all$cell_type == cell_type,]
        if(!all(is.na(meta_ct$Cluster_names))) {
            lapply(unique(meta_ct$Cluster_names),function(sub_type){
		print(paste0("saving cell type: ",cell_type,"/",sub_type))
                idx <- meta_all$cell_type == cell_type & meta_all$Cluster_names == sub_type
                file <- paste0(cell_type,"-",tolower(gsub("/","-",gsub("_","-",gsub(" ","-",sub_type)))))
                mtx_sub <- mtx[idx,]
                out_irf_path <- paste0(out_irf_dir,tissue,'_',file,'_irf-loop-mtx.tsv')
	        out_comet_path <- paste0(out_comet_dir,tissue,'_',file,'_comet-mtx.tsv')
                tryCatch({
                    write.table(
			as.data.frame(as.matrix(mtx_sub)),
			file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE
		    )
	            write.table(
			as.data.frame(t(as.matrix(mtx_sub))),
			file = out_comet_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE
		    )},error = function(e) print("matrix to large to save")
		)
            })
        } else {
	    print(paste0("saving cell type: ",cell_type))
            idx <- meta_all$cell_type == cell_type
            file <- cell_type
	    mtx_sub <- mtx[idx,]
	    out_irf_path <- paste0(out_irf_dir,tissue,'_',file,'_irf-loop-mtx.tsv')
	    out_comet_path <- paste0(out_comet_dir,tissue,'_',file,'_comet-mtx.tsv')
            tryCatch({
                write.table(
                    as.data.frame(as.matrix(mtx_sub)),
                    file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE
                )   
                write.table(
                    as.data.frame(t(as.matrix(mtx_sub))),
                    file = out_comet_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE
                )},error = function(e) print("matrix to large to save")
            )
	}
    })
}


pbdMPI::finalize()


#objects <- lapply(1:length(cell_types),function(x){
#    print(paste0("working on: ",cell_types[x]))
#    if(cell_types[x] %in% c("cd4_t_cells","conventional_cd8_t_cells")) {
#        meta <- paste0(obj_path,cell_types[x],"/",gsub("_t_cells","",cell_types[x]),"_metadata.csv")
#        cell_type_path <- paste0(obj_path,cell_types[x],"/",gsub("_t_cells","",cell_types[x]),"_rna.rds")
#    } else if (cell_types[x] == "cd4_t_helper_memory_cells") {
#        meta <- paste0(obj_path,cell_types[x],"/",gsub("t_helper_memory_cells","helper_memory",cell_types[x]),"_metadata.csv")
#        cell_type_path <- paste0(obj_path,cell_types[x],"/",gsub("t_helper_memory_cells","helper_memory",cell_types[x]),"_rna.rds")
#    } else {
#	meta <- paste0(obj_path,cell_types[x],"/",cell_types[x],"_metadata.csv")
#        cell_type_path <- paste0(obj_path,cell_types[x],"/",cell_types[x],"_rna.rds")
#    }
#    obj <- readRDS(cell_type_path)
#    meta <- read.csv(meta,row.names = 'X')
#    if("Cluster_names" %in% names(meta)) {
#        meta_summary <- meta %>% group_by(Cluster_names) %>% summarize(total = n())
#        print(meta_summary)
#    }
#    obj <- CreateSeuratObject(counts = obj,meta.data = meta)
#    if(cell_types[x] %in% c("mait_cells","progenitor_cells")) {
#	y_file <- cell_type_names[x]
#        out_irf_path <- paste0(out_irf_dir,tissue,'_',y_file,'_irf-loop-mtx.tsv')
#        write.table(t(as.data.frame(as.matrix(GetAssayData(object = obj)))),file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
#	# add comet matrix to output as well
#	# comet_mtx/_comet-mtx.tsv
#    } else {
#	subs <- unique(meta$Cluster_names)
#        lapply(1:length(subs),function(y){
             y_file <- paste0(cell_type_names[x],"-",tolower(gsub("/","-",gsub("_","-",gsub(" ","-",subs[y])))))
             obj_sub <- subset(obj,subset = Cluster_names == subs[y])
             out_irf_path <- paste0(out_irf_dir,tissue,'_',y_file,'_20k_irf-loop-mtx.tsv')
             obj_sub <- t(GetAssayData(object = obj_sub))
	     set.seed(123)
             obj_sub <- obj_sub[sample(1:dim(obj_sub)[1],20000,replace = FALSE),]
	     write.table(as.data.frame(as.matrix(obj_sub)),file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
#            tryCatch({
#	        write.table(as.data.frame(as.matrix(t(GetAssayData(object = obj_sub)))),file = out_irf_path,col.names = NA,row.names = TRUE,sep = '\t',quote = FALSE)
#	    },error = function(e) return(NULL))
#        })
#    }
#})


