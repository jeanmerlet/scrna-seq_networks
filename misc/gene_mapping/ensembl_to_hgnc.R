library(jsonlite)

genes_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/kidney/healthy/kpmp/meta/kidney_ensembl_genes.txt'
increment <- 1000


get_geneID_table <- function(genes) {
  
  genes <- paste0("\"",genes,"\"",collapse = ",")
  request <- paste0("curl -H 'Content-Type: text/json' -d '{\"Symbols\":[",genes,"]}' https://toppgene.cchmc.org/API/lookup")
  api_json <- system(request,intern = TRUE)
  result <- fromJSON(api_json)$Genes
  return(result)
  
}


genes <- read.table(genes_path, row.names=NULL, header=FALSE, sep='\n')
genes <- genes$V1
base_out_path <- '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/kidney/healthy/kpmp/meta/'

split_genes <- split(genes, ceiling(seq_along(genes)/increment))
for (i in 1:length(split_genes)) {
    if (i != 2) { next }
    print(i)
    current_split <- split_genes[[i]]
    print(length(current_split))
    print(length(unique(current_split)))
    out_path <- paste0(base_out_path, i, '_hgnc_genes.tsv')
    result <- get_geneID_table(current_split)
    #print(result)
    write.table(result, out_path, sep='\t', row.names=FALSE, col.names=TRUE)
}

