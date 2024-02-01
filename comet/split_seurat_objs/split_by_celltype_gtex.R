library(Seurat)
library(stringr)

# /lustre/orion/syb111/proj-shared/Personal/jmerlet/envs/conda/frontier/frontier_seurat


data_dir <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'
obj_paths <- list.files(data_dir, pattern='*imputed.rds', full.names=TRUE, recursive=TRUE)


split_by_celltype <- function(obj, out_prefix) {
    split_obj <- SplitObject(obj, split.by='Granular.cell.type')
    low_cellcount_objs <- list()
    for (i in seq_along(split_obj)) {
        obj <- split_obj[[i]]
        num_cells <- length(Cells(obj))
        celltype <- gsub('\\(', '', gsub('\\)', '', gsub('/', '-', gsub(' ', '-', tolower(names(split_obj)[i])))))
        if (num_cells > 100) {
            out_path <- paste0(out_prefix, '_', celltype, '_comet-mtx.tsv')
            write.table(GetAssayData(object=obj), file=out_path, row.names=rownames(obj), col.names=colnames(obj), sep='\t')
        } else {
            low_cellcount_objs <- append(low_cellcount_objs, obj)
        }
    }
    if (length(low_cellcount_objs) > 0) {
        obj <- merge(x=low_cellcount_objs[[1]], y=low_cellcount_objs)
        out_path <- paste0(out_prefix, '_other_comet-mtx.tsv')
        write.table(GetAssayData(object=obj), file=out_path, row.names=rownames(obj), col.names=colnames(obj), sep='\t')
    }
}


for (obj_path in obj_paths) {
    if (!grepl('gtex', obj_path, fixed=TRUE)) { next }
    out_dir <- paste0(dirname(dirname(obj_path)), '/comet_mtx')
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    obj <- readRDS(obj_path)
    filename <- basename(obj_path)
    tissue <- str_match(filename, '([a-zA-Z-]+)_imputed.rds')[1, 2]
    print(paste0('Splitting ', tissue, '...'))
    out_prefix <- paste0(out_dir, '/', tissue)
    split_by_celltype(obj, out_prefix)
}
