library(pbdMPI)


#print(Sys.info()["nodename"])


print_node <- function(job_ids) {
    Sys.sleep(1)
    print(Sys.info()["nodename"])
}


job_ids <- 1:7
dummy_list <- task.pull(job_ids, print_node)
