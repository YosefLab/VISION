## Functions that are used to select genes

## Functions here use a variety of criteria to reduce the number
## of genes to a more manageable size - ideally extracting the
## genes that are more biologically informative.

#' Applies filters to the inputted expression data (may remove rows)
#'
#' @param expr a numeric expression matrix (genes x cells)
#' @param filterInput list of filters to compute
#' @param threshold minimum number of samples gene must be detected in to pass
#' @param num_mad number of median absolute deviations to use in fano filter
#' @return character vector of gene names passing filter
applyFilters <- function(expr, filterInput, threshold, num_mad) {

    for (filter in filterInput) {
        if (tolower(filter) == "novar") {
            gene_passes <- filterGenesNovar(expr)
        } else if (tolower(filter) == "threshold") {
            gene_passes <- filterGenesThreshold(expr, threshold)
        } else if (tolower(filter) == "fano") {
            gene_passes <- filterGenesThreshold(expr, threshold)
            t_expr <- expr[gene_passes, ]
            gene_passes <- filterGenesFano(t_expr, num_mad)
        } else {
            stop("Filter not recognized")
        }
    }

    return(gene_passes)

}

#' Eliminate genes whose sample variance is equal to 0 (may remove rows);
#' run when --nofilter option
#' is selected
#' @importFrom stats var
#' @param data expression matrix
#' @return character vector of gene names passing filter
filterGenesNovar <- function(data) {
    message("Applying no variance filter...", appendLF = FALSE)
    genes_passing <- rownames(data)[rowVarsSp(data) != 0]
    message(paste(length(genes_passing), "Genes Retained"))
    return(genes_passing)
}

#' Filter genes whose values sum to less than some threshold value (may remove rows)
#'
#' @param data (data.frame) expression matrix
#' @param threshold (int) threshold value to filter by
#' @return character vector of gene names passing filter
filterGenesThreshold <- function(data, threshold) {
    message(
        sprintf("    Applying Threshold filter...removing genes detected in less than %i cells", threshold)
    )

    if ( is(data, "sparseMatrix") ){
        valid_rows <- rowSums(data > 0) >= threshold
    } else {
        valid_rows <- matrixStats::rowCounts(data > 0) >= threshold
    }

    genes_passing <- rownames(data)[valid_rows]

    message(
        sprintf("      Genes Retained: %i", length(genes_passing))
    )

    return(genes_passing)
}

#' Applies the Fano filter to the input data (may remove rows)
#' @importFrom stats median sd
#' @param data NUM_GENES x NUM_SAMPLES expression matrix
#' @param num_mad number of median absolute deviations
#' @param plot whether or not to generate a diagnostic plot
#' @return character vector of gene names passing filter
filterGenesFano <- function(data, num_mad=2, plot=FALSE) {

    if(plot){
        if(!requireNamespace("ggplot2", quietly = TRUE)){
            stop("ggplot2 is required to plot filter results.")
        }
    }

    message(
        sprintf("    Applying Fano filter...removing genes with Fano < %.1f MAD in each of 30 bins", num_mad)
    )

    sub_data <- data
    # if too many samples, subsample for fano filter
    if (ncol(data) > 50000) {
        # This is done using boolean indexing so it works on
        # sparse matrices

        selected_columns <- logical(ncol(data))
        selected_columns[sample(ncol(data), 50000)] <- TRUE
        sub_data <- data[,selected_columns]
    }

    mu <- Matrix::rowMeans(sub_data)
    fano <- rowVarsSp(sub_data) / mu


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

    gene_pass_names <- names(unlist(genePassList))

    message(
        sprintf("      Genes Retained: %i", length(gene_pass_names))
    )

    if(plot) {
        g <- ggplot2::ggplot() + ggplot2::aes(x=mu, y=fano, color=gene_passes) +
            ggplot2::geom_point(size=.5, alpha=.5) +
            ggplot2::scale_x_log10() +
            ggplot2::ylim(0, quantile(fano, .99)) +
            ggplot2::scale_color_manual(values=c("FALSE"="black", "TRUE"="darkcyan"))
        print(g)
    }
    return(gene_pass_names)

}
