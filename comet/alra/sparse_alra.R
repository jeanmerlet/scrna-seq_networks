library(Seurat)
library(spam)
library(spam64)
library(irlba)
library(doParallel)


# 10x test matrix directory
mtx_dir <- '/lustre/orion/syb111/proj-shared/Personal/jmerlet/projects/atopic_dermatitis/ko/data/count-matrices/SRR14253412_Solo.out/Gene/filtered'
mtx <- Read10X(mtx_dir)
mtx <- as.spam(as.matrix(mtx))
# expecting cells x genes
mtx <- t(mtx)
# eventually remove this code ^


normalize <- function(mtx) {
    print('normalizing...')
    total_umis_per_cell <- rowSums(mtx)
    mtx <- sweep(mtx, 1, total_umis_per_cell, '/')
    mtx <- mtx * 1E4
    mtx <- log(mtx + 1)
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


get_genes_to_scale_index <- function(mtx, rank_k_mtx, quantile.prob=0.001) {
    # COMMENTED OUT HAS BEEN CONVERTED TO PER-GENE
    #cat(sprintf("Find the %f quantile of each gene\n", quantile.prob))
    #rank_k_mins <- abs(apply(rank_k_mtx, 2, FUN=function(x) quantile(x, quantile.prob)))
    #print("Sweep")
    #rank_k_mtx_cor <- replace(rank_k_mtx, rank_k_mtx <= rank_k_mins[col(rank_k_mtx)], 0)
    #sd_nonzero <- function(x) sd(x[!x == 0])
    #sigma_1 <- apply(rank_k_mtx_cor, 2, sd_nonzero)
    #sigma_2 <- apply(mtx, 2, sd_nonzero)
    mu_1 <- colSums(rank_k_mtx_cor)/colSums(!!rank_k_mtx_cor)
    mu_2 <- colSums(mtx)/colSums(!!mtx)
    to_scale <- !is.na(sigma_1) & !is.na(sigma_2) & !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)
    sigma_1_2 <- sigma_2/sigma_1
    to_add  <- -1*mu_1*sigma_2/sigma_1 + mu_2
    return(list(rank_k_mtx_cor, to_scale, sigma_1_2, to_add))
}


scale_genes <- function(rank_k_mtx_cor, to_scale, sigma_1_2, to_add) {
    rank_k_mtx_tmp <- rank_k_mtx_cor[,toscale]
    rank_k_mtx_tmp <- sweep(rank_k_mtx_tmp,2, sigma_1_2[to_scale],FUN = "*")
    rank_k_mtx_tmp <- sweep(rank_k_mtx_tmp,2, toadd[to_scale],FUN = "+")
    rank_k_mtx_cor_sc <- rank_k_mtx_cor
    rank_k_mtx_cor_sc[,toscale] <- rank_k_mtx_tmp
    rank_k_mtx_cor_sc[rank_k_mtx_cor==0] = 0
    lt0 <- rank_k_mtx_cor_sc < 0
    rank_k_mtx_cor_sc[lt0] <- 0
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
        #sigma_2 <- sd(mtx[,idx][which(mtx[,idx] != 0)])
        sigma_2 <- sd(orig_col[which(orig_col != 0)])
        mu_1 <- mean(gene_col_zero[which(gene_col_zero > 0)])
        #mu_2 <- mean(mtx[,idx][which(mtx[,idx] != 0)])
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
        
        orig_values <- gene_col[gene_col > 0 & gene_col_zero == 0]
        gene_col_zero[gene_col > 0 & gene_col_zero == 0] <- orig_values
        gene_colset_zero <- c(gene_colset_zero, gene_col_zero)
    }
    gene_colset_zero <- cbind.spam(gene_colset_zero)
    return(gene_colset_zero)
}


# main function
alra <- function(mtx, normalize=FALSE, choose_k=TRUE) {
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
    rsvd_out <- irlba(mtx, k)
    u <- rsvd_out$u
    d <- rsvd_out$d
    v <- rsvd_out$v
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
    d_mtx <- diag(d[1:k])
    ud_mtx <- u_mtx %*% d_mtx
    v_mtx <- t(v[,1:k])

    num_cores <- detectCores(logical=FALSE) - 1
    cluster <- makeCluster(num_cores)
    registerDoParallel(cluster)
    num_genes <- dim(mtx)[2]
    chunk_size <- num_genes %/% (num_cores - 1)
    all_indices <- 1:num_genes
    idx_sets <- split(all_indices, ceiling(seq_along(all_indices)/chunk_size))

    print('starting per gene-col imputation')
    start <- Sys.time()
    imputed_mtx <- foreach(i=(1:length(idx_sets)), .export=c('impute_geneset'), .packages=c('spam', 'spam64', 'irlba'), .combine='cbind.spam') %dopar% {
        idx_set <- idx_sets[[i]]
        orig_colset <- mtx[,idx_set]
        gene_colset <- as.spam(ud_mtx %*% as.matrix(v_mtx[,idx_set]))
        gene_colset <- impute_geneset(orig_colset, gene_colset)
        #gene_col <- impute_gene(orig_col, gene_col)
        return(gene_colset)
    }
    print(Sys.time() - start)
    stopCluster(cluster)
    return(imputed_mtx)
}

result <- alra(mtx, normalize=TRUE, choose_k=FALSE)
