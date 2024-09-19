suppressMessages({
    #library(pbdMPI)
    library(Seurat)
    library(SeuratDisk)
    library(spam)
    library(spam64)
    library(irlba)
    library(doParallel)
})

num_nodes <- 64
num_cpus_per_node <- 16

check_rank_zero <- function(...) {
    if (pbdMPI::comm.rank() == 0) {
        for (arg in list(...)) {
            eval(arg)
        }
    }
}

choose_k <- function(mtx, k=100, noise_start=80, max_sds=6) {
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

