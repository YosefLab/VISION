#' Utility methods


#' Helper utility to group list items into batches
#'
#' This is used inside the batchSigEval function
#'
#' @param items list of items to group into batches
#' @param per_batch approximate target for size of batches
#' @param n_workers number of batches is made divisible by n_workers
#' @return list of list of items
batchify <- function(items, per_batch, n_workers = 1) {
    n_iterations <- round(length(items) / (per_batch * n_workers))

    if (n_iterations == 0){
        n_batches <- ceiling(length(items) / per_batch)
    } else {
        n_batches <- n_iterations * n_workers
    }
    per_batch <- ceiling(length(items) / n_batches)

    out <- lapply(seq(n_batches), function(i) {
        start_i <- (i - 1) * per_batch + 1
        end_i <- i * per_batch
        if (end_i > length(items)){
            end_i  <- length(items)
        }
        return(items[start_i:end_i])
    })

    return(out)
}

#' Gets number of processes available for multicore work
#'
#' @return integer number of cores
getWorkerCount <- function() {

    backends <- BiocParallel::registered()

    if ("MulticoreParam" %in% backends) {
        return(backends$MulticoreParam$workers)
    } else {
        return(1)
    }

}
