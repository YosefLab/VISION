## Functions that are used to select genes

## Functions here use a variety of criteria to reduce the number
## of genes to a more manageable size - ideally extracting the
## genes that are more biologically informative.

#' Applies filters to the inputted expression data (may remove rows)
#'
#' @param data ExpressionData object
#' @param threshold minimum number of samples gene must be detected in to pass
#' @param filterInput list of filters to compute
#' @return The updated ExpressionData object
applyFilters <- function(data, threshold, filterInput) {
    expr <- getExprData(data)

    # If the filterInput is a list of items, assume it's a list
    # of genes and just use that to filter the matrix
    if(length(filterInput) > 1) {
        f_expr = expr[filterInput, ];
        data@fanoFilter = f_expr
        return(data)
    }

    for (filter in filterInput) {
        if (tolower(filter) == "novar") {
            f_expr <- filterGenesNovar(expr)
            data@fanoFilter <- f_expr
        } else if (tolower(filter) == "threshold") {
            f_expr <- filterGenesThreshold(expr, threshold)
            data@fanoFilter <- f_expr
        } else if (tolower(filter) == "fano") {
            t_expr <- filterGenesThreshold(expr, threshold)
            f_expr <- filterGenesFano(t_expr)
            data@fanoFilter <- f_expr
        } else {
            stop("Filter not recognized")
        }
    }

    return(data)

}

#' Eliminate genes whose sample variance is equal to 0 (may remove rows);
#' run when --nofilter option
#' is selected
#' @importFrom stats var
#' @param data expression matrix
#' @return filtered expression matrix
filterGenesNovar <- function(data) {
    message("Applying no variance filter...", appendLF=FALSE)
    fdata <- data[apply(data, 1, var) != 0,]
    message(paste(nrow(fdata), "Genes Retained"))
    return(fdata)
}

#' Filter genes whose values sum to less than some threshold value (may remove rows)
#'
#' @param data (data.frame) expression matrix
#' @param threshold (int) threshold value to filter by
#' @return filtered expression matrix
filterGenesThreshold <- function(data, threshold) {
    message("Applying threshold filter...", appendLF=FALSE)
    fdata <- data[rowSums(data > 0) > threshold,]
    message(paste(nrow(fdata), "Genes Retained"))
    return(fdata)
}

#' Applies the Fano filter to the input data (may remove rows)
#' @importFrom stats median sd
#' @param data NUM_GENES x NUM_SAMPLES expression matrix
#' @param num_mad number of median absolute deviations
#' @param plot whether or not to generate a diagnostic plot
#' @return NUM_GENES_PASSED_FANO_FILTER x NUM_SAMPLES filtered expression matrix
filterGenesFano <- function(data, num_mad=2, plot=FALSE) {

    if(plot){
        if(!requireNamespace("ggplot2", quietly = TRUE)){
            stop("ggplot2 is required to plot filter results.")
        }
    }

    message("Applying fano filter...", appendLF=FALSE)
    sub_data <- data
    # if too many samples, subsample for fano filter
    if (ncol(data) > 50000) {
    sub_data <- data[,sample(ncol(data), 50000)]
    }

    mu <- Matrix::rowMeans(sub_data)
    fano <- rowVars(sub_data) / mu


    aa <- order(mu)
    mu_sort <- mu[aa]
    fano_sort <- fano[aa]


    N_QUANTS <- 30
    m <- floor(length(mu_sort) / N_QUANTS)

    gene_passes <- rep(FALSE, nrow(sub_data))

    genePassList <- lapply(0:N_QUANTS, function(i) {
        if (i == N_QUANTS-1) {
            rr <- seq(i*m+1, length(mu_sort))
        } else {
            rr <- seq(i*m+1, (i+1)*m)
        }

        mu_quant <- mu_sort[rr]
        mu_quant[mu_quant == 0] <- 1
        fano_quant <- fano_sort[rr]
        mad_quant <- median(abs(fano_quant - median(fano_quant)))
        gene_passes_quant <- (fano_quant > (median(fano_quant)
                                            + num_mad * mad_quant))
        gene_passes_i <- which(gene_passes_quant != 0) + (i*m)
        return(gene_passes_i)

    })

    gene_passes[unlist(genePassList)] <- TRUE
    gene_passes <- gene_passes[order(aa)]
    fdata <- data[gene_passes,]
    message(paste(nrow(fdata), "Genes Retained"))

    if(plot) {
        g <- ggplot() + aes(x=mu, y=fano, color=gene_passes) +
            geom_point(size=.5, alpha=.5) +
            scale_x_log10() +
            ylim(0, quantile(fano, .99)) +
            scale_color_manual(values=c("FALSE"="black", "TRUE"="darkcyan"))
        print(g)
    }
    return(fdata)

}

rowVars <- function(x) {
    return (rowSums(( x - Matrix::rowMeans(x))^2) / (ncol(x) - 1))
}
