library(Seurat)
library(stringr)
library(readr)
library(DropletUtils)
library(Matrix)
library(R.utils)


# /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat

obj_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/kidney/healthy/kpmp/raw/all_kidney.rds'
celltype_colname <- 'cell_type'
tissue <- 'kidney'


apply_imputation <- function(obj) {
    print('starting alra')
    source('/lustre/orion/syb111/proj-shared/Personal/jmerlet/projects/sc-flow/scripts/preprocess/alra.R')
    obj <- NormalizeData(obj, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
    alra_result <- alra(t(as.matrix(GetAssayData(object=obj, slot = "data"))))
    obj_alra <- t(alra_result[[3]])
    colnames(obj_alra) <- rownames(obj@meta.data)
    obj_alra <- Matrix(obj_alra,sparse = T)
    obj <- SetAssayData(object = obj,slot = "data",new.data = obj_alra)
    #out_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/seurat/all_liver_imputed.rds'
    #print('saving rds')
    #saveRDS(obj, out_path)
    print('alra done')
    return(obj)
}

split_by_celltype <- function(tissue,obj, out_comet_dir, out_irf_dir, celltype_colname) {
    # subset to healthy
    #obj <- subset(obj, subset = disease == 'normal')
    print('starting split')
    split_obj <- SplitObject(obj, split.by=celltype_colname)
    low_cellcount_objs <- list()
    for (i in seq_along(split_obj)) {
        obj <- split_obj[[i]]
        num_cells <- length(Cells(obj))
        celltype <- gsub('\\(', '', gsub('\\)', '', gsub('/', '-', gsub(' ', '-', tolower(names(split_obj)[i])))))
        if (num_cells > 100) {
            out_comet_path <- paste0(out_comet_dir, tissue,'_', celltype, '_comet-mtx.tsv')
            out_irf_path <- paste0(out_irf_dir, tissue,'_', celltype, '_irf-loop-mtx.tsv')
	    # matrices already have rownames and column names!
	    print(out_comet_path)
	    print(out_irf_path)
            write.table(GetAssayData(object=obj), file=out_comet_path,col.names = NA,row.names = TRUE, sep='\t', quote=FALSE)
            write.table(t(GetAssayData(object=obj)), file=out_irf_path,col.names = NA,row.names = TRUE, sep='\t', quote=FALSE)
        } else {
            low_cellcount_objs <- append(low_cellcount_objs, obj)
        }
    }
    if (length(low_cellcount_objs) > 0) {
        obj <- merge(x=low_cellcount_objs[[1]], y=low_cellcount_objs)
        out_comet_path <- paste0(out_comet_dir, tissue,'_other_comet-mtx.tsv')
        out_irf_path <- paste0(out_irf_dir, tissue,'_other_irf-loop-mtx.tsv')
        write.table(GetAssayData(object=obj), file=out_comet_path,col.names = NA,row.names = TRUE, sep='\t', quote=FALSE)
        write.table(t(GetAssayData(object=obj)), file=out_irf_path,col.names = NA,row.names = TRUE, sep='\t', quote=FALSE)
    }
    print('split done')
}


obj <- readRDS(obj_path)
# added forward slash to end of paths (due to paste!)
out_comet_dir <- paste0(dirname(dirname(obj_path)), '/comet_mtx/')
out_irf_dir <- paste0(dirname(dirname(obj_path)), '/irf-loop_mtx/')
if (!dir.exists(out_comet_dir)) {
    dir.create(out_comet_dir)
}
if (!dir.exists(out_irf_dir)) {
    dir.create(out_irf_dir)
}

obj <- apply_imputation(obj)
split_by_celltype(tissue, obj, out_comet_dir, out_irf_dir, celltype_colname)
