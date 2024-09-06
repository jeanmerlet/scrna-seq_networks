library(Seurat)
library(spam)
library(spam64)
library(irlba)
library(doParallel)


seurat <- TRUE
if (!seurat) {
    # 10x test matrix directory
    mtx_dir <- '/lustre/orion/syb111/proj-shared/Personal/jmerlet/projects/atopic_dermatitis/ko/data/count-matrices/SRR14253412_Solo.out/Gene/filtered'
    mtx <- Read10X(mtx_dir)
    #mtx <- as.spam(as.matrix(mtx))
    mtx <- as.matrix(mtx)
    # expecting cells x genes
    mtx <- t(mtx)
} else {
    # test a large matrix
    print('loading seurat object...')
    obj <- readRDS('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/arabidopsis/root/healthy/ohler/seurat/root_atlas_seu4.rds')
    print('converting seurat object...')
    DefaultAssay(obj) <- 'RNA'
    #mtx <- as.spam(t(as.matrix(GetAssayData(object=obj, slot='data'))))
    mtx <- t(as.matrix(GetAssayData(object=obj, slot='data')))
    rm(obj)
}


normalize <- function(mtx) {
    print('normalizing...')
    total_umis_per_cell <- rowSums(mtx)
    mtx <- sweep(mtx, 1, total_umis_per_cell, '/')
    mtx <- mtx * 1E4
    mtx <- log(mtx + 1)
    mtx <- as.spam(mtx)
    print('normalized!')
    return(mtx)
}


choose_k <- function(mtx, k=100, noise_start=80, thresh=6) {
    print('choosing k...')
    noise_svals <- noise_start:k
    rsvd_out <- irlba(mtx, nv=k)
    diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]
    mu <- mean(diffs[noise_svals-1])
    sigma <- sd(diffs[noise_svals-1])
    num_of_sds <- (diffs-mu)/sigma
    k <- max(which(num_of_sds > thresh))
    print(paste0('k=', k, ' chosen!'))
    return(k)
}


impute_geneset <- function(orig_colset, gene_colset, quantile.prob=0.001) {
    # loop over each gene column in gene_colset
    gene_colset_zero <- c()
    for (i in 1:dim(gene_colset)[2]) {
        orig_col <- orig_colset[,i]
        gene_col <- gene_colset[,i]
        #orig_col <- as.spam(orig_colset[,i])
        #gene_col <- as.spam(gene_colset[,i])
        
        gene_quants <- abs(quantile(gene_col, quantile.prob))
        gene_col_zero <- replace(gene_col, gene_col <= gene_quants, 0)
        
        sigma_1 <- sd(gene_col_zero[which(gene_col_zero != 0)])
        sigma_2 <- sd(orig_col[which(orig_col != 0)])
        mu_1 <- mean(gene_col_zero[which(gene_col_zero > 0)])
        mu_2 <- mean(orig_col[which(orig_col != 0)])

        sigma_1_2 <- sigma_2/sigma_1
        mu_1_2 <- -1*mu_1*sigma_2/sigma_1 + mu_2

        # determine whether or not imputation is needed for this gene
        # and assign to boolean do_imputation
        do_imputation <- !is.na(sigma_1) & !is.na(sigma_2) &
                         !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)

        if (do_imputation) {
            # do imputation
            zeros <- which(gene_col_zero == 0)
            gene_col_zero <- gene_col_zero * sigma_1_2
            gene_col_zero <- gene_col_zero + mu_1_2
            gene_col_zero[zeros] = 0
        }
        negative <- gene_col_zero < 0
        gene_col_zero[negative] <- 0
        gene_col_zero[orig_col > 0 & gene_col_zero == 0] <- gene_col[orig_col > 0 & gene_col_zero == 0]
        gene_colset_zero <- c(gene_colset_zero, as.spam(gene_col_zero))
    }
    gene_colset_zero <- do.call('cbind.spam', gene_colset_zero)
    return(gene_colset_zero)
}


# main function
sparse_alra <- function(mtx, normalize=FALSE, choose_k=TRUE) {
    num_cells <- dim(mtx)[1]
    num_genes <- dim(mtx)[2]
    sprintf('imputing matrix with %s cells and %s genes', num_cells, num_genes)
    ## need the count mtx to be cells x genes (rows x cols) ##
    # normalize
    if (normalize) { mtx <- normalize(mtx) }
    # pick k heuristically
    if (choose_k) {
        k <- choose_k(mtx)
    } else {
        k <- 50
    }
    # randomized svd
    start <- Sys.time()
    print('rsvd...')
    rsvd_out <- irlba(mtx, k)
    #u <- rsvd_out$u
    #d <- rsvd_out$d
    #v <- rsvd_out$v
    u <- as.spam(rsvd_out$u)
    d <- as.spam(rsvd_out$d)
    v <- as.spam(rsvd_out$v)
    rm(rsvd_out)
    # note: dim(u[,1:k] %*% diag(d[1:k]) %*% t(v[,1:k])) == dim(mtx)
    # note for linear algebra retards (us): matrix mult needs
    # num columns 1st matrix == num rows 2nd matrix
    # u 3197 x k
    # d k x k
    # v k x 57186
    
    # test code below
    #i <- 1
    #orig_col <- mtx[,i]
    #gene_col <- u[,1:k] %*% diag(d[1:k]) %*% as.matrix(t(v[, 1:k])[,i])
    #gene_col <- impute_gene(orig_col, gene_col)
    # test code above

    # original code to compare to
    #test_result <- u[,1:k] %*% diag(d[1:k]) %*% t(v[,1:k])

    # doparallel attempt

    # linear algrebra
    u_mtx <- u[,1:k]
    d_mtx <- as.spam(diag(d[1:k]))
    ud_mtx <- u_mtx %*% d_mtx
    v_mtx <- t(v[,1:k])

    print('rsvd done!')
    print(Sys.time() - start)

    num_cores <- detectCores(logical=FALSE) - 1
    print(num_cores)
    cluster <- makeCluster(num_cores)
    print(cluster)
    registerDoParallel(cluster)
    num_genes <- dim(mtx)[2]
    chunk_size <- num_genes %/% (num_cores - 1)
    all_indices <- 1:num_genes
    idx_sets <- split(all_indices, ceiling(seq_along(all_indices)/chunk_size))

    print('Imputing...')
    start <- Sys.time()
    imputed_mtx <- foreach(i=(1:length(idx_sets)), .export=c('impute_geneset'), .packages=c('spam', 'spam64', 'irlba'), .combine='cbind.spam') %dopar% {
        idx_set <- idx_sets[[i]]
        orig_colset <- mtx[,idx_set]
        gene_colset <- ud_mtx %*% v_mtx[,idx_set]
        gene_colset <- impute_geneset(orig_colset, gene_colset)
        return(gene_colset)
    }
    print('Imputation done!')
    print(Sys.time() - start)
    stopCluster(cluster)
    return(imputed_mtx)
}

sparse_result <- sparse_alra(mtx, normalize=TRUE, choose_k=FALSE)
