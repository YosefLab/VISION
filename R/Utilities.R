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

#' Check's the version of the FastProject object and displays error if necessary
#'
#' @param object FastProject object
#' @return NULL
versionCheck <- function(object) {

    templateStr <- paste0(
        "This FastProject object was created with an older version of the library.",
        "  To view, either install an older version (commit #COMMIT) from GitHub (www.github.com/YosefLab/FastProjectR) or",
        " recreate the object and re-run analyze"
    )

    if (!.hasSlot(object, "version")) {
        msg <- gsub("#COMMIT", "0cd5268", templateStr)
        stop(msg, call. = FALSE)
    }

    if(object@version < 1.0) {
        msg <- gsub("#COMMIT", "0cd5268", templateStr)
        stop(msg, call. = FALSE)
    }

    # Add new commit hashes here as version increases and breaks backward compatibility

    return()
}
