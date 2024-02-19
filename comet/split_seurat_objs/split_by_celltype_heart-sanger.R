library(Seurat)
library(stringr)
library(readr)
library(DropletUtils)
library(Matrix)
library(R.utils)

source('/gpfs/alpine2/syb112/proj-shared/Projects/sc_pens/scripts/comet/split_seurat_objs/alra_loop.R')

# general steps:
# read10x
# createseuratobj
# addmetadata
# subset if needed

mtx_dir <- '/gpfs/alpine2/syb112/proj-shared/Projects/sc_pens/data/human/heart/sanger/raw'
meta_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/sc_pens/data/human/heart/sanger/meta/bc_meta.tsv'
out_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/sc_pens/data/human/heart/sanger/seurat/seurat_heart_imputed.rds'
celltype_colname <- 'cell_type'
tissue <- 'heart'


tenx_to_seurat <- function(obj) {
    metadata <- read.table(meta_path, sep='\t', header=TRUE)
    print(str(metadata))
    obj <- CreateSeuratObject(obj)
    obj <- AddMetaData(obj, metadata)
    print(str(obj))
    return(obj)
}


apply_imputation <- function(obj) {
    print('starting alra')
    obj <- NormalizeData(obj, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
    alra_result <- alra(t(as.matrix(GetAssayData(object=obj, slot='data'))), k=0)
    #obj_alra <- t(alra_result[[3]])
    obj_alra <- t(alra_result)
    colnames(obj_alra) <- rownames(obj@meta.data)
    obj_alra <- Matrix(obj_alra, sparse=T)
    obj <- SetAssayData(object=obj, slot='data', new.data=obj_alra)
    print('saving rds')
    saveRDS(obj, out_path)
    print('alra done')
    return(obj)
}


split_by_celltype <- function(tissue,obj, out_comet_dir, out_irf_dir, celltype_colname) {
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
}


obj <- Read10X(mtx_dir)
print('10x read in')
# added forward slash to end of paths (due to paste!)
out_comet_dir <- paste0(dirname(mtx_dir), '/comet_mtx/')
out_irf_dir <- paste0(dirname(mtx_dir), '/irf-loop_mtx/')
if (!dir.exists(out_comet_dir)) {
    dir.create(out_comet_dir)
}
if (!dir.exists(out_irf_dir)) {
    dir.create(out_irf_dir)
}

obj <- tenx_to_seurat(obj)
print('metadata added')
obj <- apply_imputation(obj)
print('imputed')
split_by_celltype(tissue, obj, out_comet_dir, out_irf_dir, celltype_colname)
print('split done')
