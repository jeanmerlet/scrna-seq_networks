library(Seurat)
library(stringr)
library(readr)
library(DropletUtils)
library(Matrix)
library(R.utils)

source('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/split_seurat_objs/alra_loop.R')

# general steps:
# read10x
# createseuratobj
# addmetadata
# subset if needed

mtx_dir <- '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/heart/healthy/sanger/raw'
meta_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/heart/healthy/sanger/meta/bc_meta.tsv'
out_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/heart/healthy/sanger/seurat/heart_imputed.tsv'
celltype_colname <- 'cell_type'
tissue <- 'heart'


tenx_to_seurat <- function(obj) {
    metadata <- read.table(meta_path, sep='\t', header=TRUE)
    obj <- CreateSeuratObject(obj)
    obj <- AddMetaData(obj, metadata)
    return(obj)
}


randomized.svd <- function(A, k, method, q, mkl.seed=-1) {
    out <- setNames(vector('list', 3), c('u', 'd', 'v'))
    if (method == 'rsvd') {
        library(rsvd)
        out <- rsvd(A, k, q=q)
    } else if (method == 'rsvd-mkl') {
        library(fastRPCA)
        fastPCAOut <- fastPCA(A, k=k, its=q, l=(k+10), seed=mkl.seed)
        out$u <- fastPCAOut$U
        out$v <- fastPCAOut$V
        out$d <- diag(fastPCAOut$S)
    } else if (method == 'irlba') {
        # not implemented
        library(irlba)
    }
    mat_rank_k <- out$u[,1:k] %*% diag(out$d[1:k]) %*% t(out$v[,1:k])
    return(mat_rank_k)
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


scale_genes <- function(mat, mat_rank_k, quantile.prob) {
    cat(sprintf('Find the %f quantile of each gene\n', quantile.prob))
    mat_rank_k_mins <- abs(apply(mat_rank_k, 2, FUN=function(x) quantile(x, quantile.prob)))
    cat('Sweep\n')
    mat_rank_k_cor <- replace(mat_rank_k, mat_rank_k <= mat_rank_k_mins[col(mat_rank_k)], 0)
    sd_nonzero <- function(x) sd(x[!x == 0])
    sigma_1 <- apply(mat_rank_k_cor, 2, sd_nonzero)
    sigma_2 <- apply(mat, 2, sd_nonzero)
    mu_1 <- colSums(mat_rank_k_cor)/colSums(!!mat_rank_k_cor)
    mu_2 <- colSums(mat)/colSums(!!mat)
    toscale <- !is.na(sigma_1) & !is.na(sigma_2) & !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)
    sigma_1_2 <- sigma_2/sigma_1
    toadd  <- -1*mu_1*sigma_2/sigma_1 + mu_2
    mat_rank_k_cor_sc <- mat_rank_k_cor
    cat(sprintf('Scaling all except for %d columns (genes)\n', sum(!toscale)))
    idx <- 0
    for (value in toscale) {
        idx <- idx + 1
        if (idx %% 1000 == 0) { print(idx) }
        if (!value) { next }
        gene_vector <- mat_rank_k_cor[,idx] * sigma_1_2[idx]
        gene_vector <- gene_vector + toadd[idx]
        gene_vector[mat_rank_k_cor[,idx] == 0] = 0
        mat_rank_k_cor_sc[,idx] <- gene_vector
    }
    mat_rank_k_cor_sc[mat_rank_k_cor==0] = 0
    lt0 <- mat_rank_k_cor_sc < 0
    mat_rank_k_cor_sc[lt0] <- 0
    cat(sprintf('%.2f%% of the values became negative in the scaling process and were set to zero\n', 100*sum(lt0)/(nrow(mat)*ncol(mat))))
    return(mat_rank_k_cor_sc)
}


run_alra <- function(obj, k, method, q, quantile.prob) {
    mat <- t(as.matrix(GetAssayData(object=obj, slot='data')))
    mat_rank_k <- randomized.svd(mat, k, method, q)
    mat_rank_k_cor_sc <- scale_genes(mat, mat_rank_k, quantile.prob)
    remove(mat_rank_k)
    # set originally nonzero values that are now 0 from thresholding
    # back to their pre-imputation nonzero values
    originally_nonzero <- mat > 0
    mat_rank_k_cor_sc[originally_nonzero & mat_rank_k_cor_sc == 0] <- mat[originally_nonzero & mat_rank_k_cor_sc == 0]
    colnames(mat_rank_k_cor_sc) <- colnames(mat)
    #original_nz <- sum(mat > 0)/(nrow(mat)*ncol(mat))
    #completed_nz <- sum(mat_rank_k_cor_sc > 0)/(nrow(mat)*ncol(mat))
    #cat(sprintf('The matrix went from %.2f%% nonzero to %.2f%% nonzero\n', 100*original_nz, 100*completed_nz))
    remove(mat)
    mat_rank_k_cor_sc <- t(mat_rank_k_cor_sc)
    colnames(mat_rank_k_cor_sc) <- rownames(obj@meta.data)
    #print('saving imputed matrix to tsv')
    write.table(mat_rank_k_cor_sc, file=out_path, sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
    #saveRDS(mat_rank_k_cor_sc, out_path)
    # cannot dgcmatrix on matrices with more than ~180k cells x ~35k genes
    #mat_rank_k_cor_sc <- Matrix(mat_rank_k_cor_sc, sparse=T)
    #obj <- SetAssayData(object=obj, slot='data', new.data=mat_rank_k_cor_sc)
    return(mat_rank_k_cor_sc)
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
obj <- NormalizeData(obj, normalization.method='LogNormalize', scale.factor=10000, verbose=FALSE)
print('seurat obj normalized')
mat <- run_alra(obj, k=69, method='rsvd', q=10, quantile.prob=0.001)
print('alra done')
#print('saving rds')
##saveRDS(obj, out_path)
#print('splitting by cell type')
#split_by_celltype(tissue, obj, out_comet_dir, out_irf_dir, celltype_colname)
print('done')
