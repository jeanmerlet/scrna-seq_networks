library(Seurat)
library(spam)
library(spam64)
library(irlba)
library(doParallel)
library(ggplot2)
library(rsvd)


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


# main function
compare_alras <- function(mtx, normalize=TRUE, k=50) {
    num_cells <- dim(mtx)[1]
    num_genes <- dim(mtx)[2]
    sprintf('imputing matrix with %s cells and %s genes', num_cells, num_genes)
    ## need the count mtx to be cells x genes (rows x cols) ##
    # normalize
    if (normalize) { mtx <- normalize(mtx) }
    # randomized svd
    before <- FALSE
    if (before) {
        print('irlbaing...')
        irlba_out <- irlba(mtx, nv=k)
        u <- irlba_out$u
        d <- irlba_out$d
        v <- irlba_out$v
        irlba_result <- u[,1:k] %*% diag(d[1:k]) %*% t(v[,1:k])

        print('rsvding...')
        rsvd_out <- rsvd(mtx, k=k, q=10)
        u <- rsvd_out$u
        d <- rsvd_out$d
        v <- rsvd_out$v
        rsvd_result <- u[,1:k] %*% diag(d[1:k]) %*% t(v[,1:k])

        diff_before <- abs(irlba_result - rsvd_result)
        # histogram before
        print('histogramming...')
        out_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/alra/irlba-rsvd_histogram_before.png'
        # flatten mtx
        flat_mtx <- data.frame('diffs'=as.vector(diff_before))
        hist <- ggplot(data=flat_mtx, aes(x=diffs)) + geom_histogram(bins=100)
        ggsave(out_path, plot=hist)
    }

    # heatmap before
    #print('heatmapping...')
    #out_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/alra/irlba-rsvd_heatmap_before.png'
    #png(file=out_path)
    #pal <- colorRampPalette(c("red", "yellow"), space = "rgb")
    #max_diff <- max(diff_before)
    #print(levelplot(diff_before, main="", xlab="", ylab="", col.regions=pal(4), cuts=3, at=seq(0,max_diff,max_diff/2)))
    #dev.off()



    # after imputation
    source('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/alra/alra.R')
    source('/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/alra/sparse_alra.R')
    alra_mtx <- alra(as.matrix(mtx), k=50, q=10)
    sparse_mtx <- sparse_alra(mtx, normalize=FALSE, choose_k=FALSE)
    print(dim(alra_mtx))
    print(dim(sparse_mtx))
    diff_after <- alra_mtx[[3]] - as.matrix(sparse_mtx)

    print('histogramming...')
    out_path <- '/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/code/comet/alra/irlba-rsvd_histogram_after.png'
    # flatten mtx
    flat_mtx <- data.frame('diffs'=as.vector(diff_after))
    hist <- ggplot(data=flat_mtx, aes(x=diffs)) + geom_histogram(bins=100)
    ggsave(out_path, plot=hist)
    
    
    return()
}

compare_alras(mtx)
