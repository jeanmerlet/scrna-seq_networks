suppressMessages({
    library(pbdMPI)
    library(Seurat)
    library(spam)
    library(spam64)
    library(irlba)
    library(doParallel)
})


# num_nodes should be number of nodes you allocate minus 1
#TODO: find Sys() call to detect num of nodes
num_nodes <- 7
num_cpus_per_node <- 2


#TODO: fix normalize to work on spam matrices (just use Seurat normalize)
#TODO: break up and modularize sparse_alra, impute_geneset
#TODO: add comments / modify comments
#TODO: address this seurat 5.0 warning message: 'The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0. Please use the `layer` argument instead.'
#TODO: create package skeleton


check_rank_zero <- function(...) {
    if (pbdMPI::comm.rank() == 0) {
        for (arg in list(...)) {
            eval(arg)
        }
    }
}


seurat <- FALSE
if (!seurat) {
    # 10x test matrix directory
    mtx_dir <- '/lustre/orion/syb111/proj-shared/Personal/jmerlet/projects/atopic_dermatitis/ko/data/count-matrices/SRR14253412_Solo.out/Gene/filtered'
    mtx <- Read10X(mtx_dir)
    if (normalize) { mtx <- NormalizeData(mtx) }
    mtx <- t(as.spam.dgCMatrix(mtx))
} else {
    # test a large matrix
    check_rank_zero(print('Loading seurat object...'), start <- Sys.time())
    obj <- readRDS('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/arabidopsis/root/healthy/ohler/seurat/root_atlas_seu4.rds')
    check_rank_zero(print('Seurat object loaded'), print(Sys.time() - start),
                    print('Converting seurat object...'), start <- Sys.time())
    DefaultAssay(obj) <- 'RNA'
    #mtx <- as.spam(t(as.matrix(GetAssayData(object=obj, slot='data'))))
    mtx <- t(as.matrix(GetAssayData(object=obj, slot='data')))
    rm(obj)
    check_rank_zero(print('Seurat object converted'), print(Sys.time() - start))
}


choose_k <- function(mtx, k, noise_start=80, max_sds=6) {
    check_rank_zero(print('Choosing k...'), start <- Sys.time())
    noise_svals <- noise_start:k
    rsvd_out <- irlba(mtx, nv=k)
    diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]
    mu <- mean(diffs[noise_svals-1])
    sigma <- sd(diffs[noise_svals-1])
    num_of_sds <- (diffs-mu)/sigma
    k <- max(which(num_of_sds > max_sds))
    check_rank_zero(print(paste0('k=', k, ' chosen')), print(Sys.time() - start))
    return(k)
}


impute_geneset <- function(orig_colset, gene_colset, quantile.prob=0.001) {
    # loop over each gene column in gene_colset
    gene_colset_zero <- c()
    for (i in 1:dim(gene_colset)[2]) {
        orig_col <- orig_colset[,i]
        gene_col <- gene_colset[,i]
        
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


impute_geneset_foreach <- function(orig_colset, gene_colset, quantile.prob=0.001) {
    #num_cores <- detectCores(logical=FALSE) - 1
    #if (pbdMPI::comm.rank() == 1) { print(num_cores) }
    #cluster <- makeCluster(num_cores)
    cluster <- makeCluster(num_cpus_per_node)
    if (pbdMPI::comm.rank() == 1) { print(cluster) }
    registerDoParallel(cluster)
    # loop over each gene column in gene_colset
    #gene_colset_zero <- c()
    gene_colset_zero <- foreach(i=(1:dim(gene_colset)[2]), .packages=c('spam', 'spam64'), .combine='cbind.spam') %dopar% {
    #for (i in 1:dim(gene_colset)[2]) {
        orig_col <- orig_colset[,i]
        gene_col <- gene_colset[,i]
        
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
        #gene_colset_zero <- c(gene_colset_zero, as.spam(gene_col_zero))

        return(gene_col_zero)
    }
    #gene_colset_zero <- do.call('cbind.spam', gene_colset_zero)
    return(gene_colset_zero)
    stopCluster(cluster)
}


mpi_split_colsets <- function(job_id, idx_sets, mtx, ud_mtx, v_mtx) {
    idx_set <- idx_sets[[job_id]]
    orig_colset <- mtx[,idx_set]
    gene_colset <- ud_mtx %*% v_mtx[,idx_set]
    #gene_colset <- impute_geneset(orig_colset, gene_colset)
    gene_colset <- impute_geneset_foreach(orig_colset, gene_colset)
    return(gene_colset)
}


# main function
sparse_alra <- function(mtx, normalize=FALSE, choose_k=TRUE) {
    num_cells <- dim(mtx)[1]
    num_genes <- dim(mtx)[2]
    check_rank_zero(sprintf('Matrix has %s cells and %s genes', num_cells, num_genes))
    ## need the count mtx to be cells x genes (rows x cols) ##
    # normalize
    if (normalize) { mtx <- normalize(mtx) }
    # pick k heuristically
    if (pbdMPI::comm.rank() == 0) {
        if (choose_k) {
            k <- choose_k(mtx)
        } else {
            k <- 50
        }
    } else {
        k <- NULL
    }
    k <- pbdMPI::bcast(k, rank.source = 0)
    # randomized svd
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
    #u <- rsvd_out$u
    #d <- rsvd_out$d
    #v <- rsvd_out$v
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

    if (pbdMPI::comm.rank() == 0) {
        print('Rsvd done')
        print(Sys.time() - start)
    }

    num_genes <- dim(mtx)[2]
    all_indices <- 1:num_genes
    chunk_size <- num_genes %/% (num_nodes - 1)
    idx_sets <- split(all_indices, ceiling(seq_along(all_indices)/chunk_size))

    if (pbdMPI::comm.rank() == 0) {
        print('Imputing...')
        start <- Sys.time()
    }
    job_ids <- seq(1, length(idx_sets))
    gene_colsets <- pbdMPI::task.pull(job_ids, mpi_split_colsets,
                                      idx_sets, mtx, ud_mtx, v_mtx)
    if (pbdMPI::comm.rank() == 0) {
        print('Imputation done')
        print(Sys.time() - start)
        print('Cbinding...')
        start <- Sys.time()
        imputed_mtx <- do.call('cbind.spam', gene_colsets)
        print('Cbinding done')
        print(Sys.time() - start)
        print(str(imputed_mtx))
        print(class(imputed_mtx))
        print(dim(imputed_mtx))
        # change this to save?
        return(imputed_mtx)
    }
}


sparse_result <- sparse_alra(mtx, normalize=TRUE, choose_k=FALSE)
#check_rank_zero(print(dim(sparse_result)))

finalize()
