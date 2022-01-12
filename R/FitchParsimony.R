#' Comptue fitch parsimony for each node in a tree
#' 
#' @param tree a rooted tree of class `phylo`
#' @param metaData a named list mapping each leaf of the tree to a category
#' @return a named list mapping each node in the tree to a normalized parsimony score
computeFitchHartiganParsimonyPerNode <- function(tree, metaData) {

    parsimonyScores <- list()
    for (node in depthFirstTraverse(tree)) {
        if (!is_tip(tree, node)) {
            score <- computeFitchHartiganParsimony(tree, metaData, source=node)
            parsimonyScores[[as.character(node)]] <- score
        }
    }
    return(parsimonyScores)

}

#' Computes the fitch parsimony of a tree with respect to categorical meta data
#' 
#' @param tree a rooted tree of class `phylo`
#' @param metaData a named list mapping each leaf of the tree to a category
#' @param source node to start analysis
#' @return a normalized parsimony score
computeFitchHartiganParsimony <- function(tree, metaData, source=NULL) {

    if (is.numeric(metaData)) {
        stop("Meta data must be character or factor data.")
    }

    possibleLabels <- bottomUpFitchHartigan(tree, metaData, source=source)
    assignments <- topDownFitchHartigan(tree, possibleLabels, source=source)
    parsimony <- scoreParsimony(tree, assignments, source=source)

    nNode <- length(depthFirstTraverse(tree, source=source))
    return(parsimony / (nNode - 1))
}

#' Bottom up Fitch-Hartigan
#' 
#' Infers a set of possible labels for each internal node of a tree according
#' to the Fitch-Hartigan algorithm.
#' 
#' @param tree a rooted tree of class `phylo`
#' @param metaData a named list mapping each leaf to a category
#' @return a named list mapping each node in the tree to a set of possible labels
bottomUpFitchHartigan <- function(tree, metaData, source=NULL) {

    possibleLabels <- list()
    for (node in depthFirstTraverse(tree, source=source)) {
        
        if (is_tip(tree, node)) {
            cell.name <- tree$tip.label[[node]]
            possibleLabels[as.character(node)] <- list(metaData[[cell.name]])
        } else {
            children <- as.character(get_children(tree, node))
            labels <- unlist(lapply(children, function(x) possibleLabels[[x]]))
            frequencies <- table(labels)
            optimal.labels <- list(names(which(frequencies == max(frequencies))))
            possibleLabels[as.character(node)] <- optimal.labels
        }
    }
    return(possibleLabels)

}

#' Top Down Fitch-Hartigan
#' 
#' Performs an iteration of the Fitch-Hartigan top-down algorithm, providing a
#' maximum-parsimony assignmenet for each node in the tree.
#' 
#' @param tree a rooted tree of class `phylo`
#' @param possibleLabels possible labels for each node, a named list
#' @return a named list mapping nodes to a maximum parsimony assignment
topDownFitchHartigan <- function(tree, possibleLabels, source=NULL) {

    assignments <- list()
    if (is.null(source)) {
        root <- find_root(tree)
    } else {
        root <- source
    }
    
    for (node in depthFirstTraverse(tree, postorder=F, source=root)) {
        if (node == root) {
            assignments[as.character(node)] <- sample(possibleLabels[[as.character(node)]], 1)
        } else {
            parent <- get_parent(tree, node)

            parent.state <- assignments[[as.character(parent)]]
            node.states <- possibleLabels[[as.character(node)]]

            if (parent.state %in% node.states) {
                assignments[[as.character(node)]] <- parent.state
            } else {
                assignments[[as.character(node)]] <- sample(node.states, 1)
            }
        }
    }
    return(assignments)
}

#' Scores parsimony from a tree with assignments
#' 
#' Counts the number of transitions on a tree with assignments.
#' 
#' @param tree a rooted tree of class `phylo`
#' @param assignments a named list mapping node names to an assignment
#' @return a parsimony score
scoreParsimony <- function(tree, assignments, source=NULL) {

    parsimony <- 0
    if (is.null(source)) {
        root <- find_root(tree)
    } else {
        root <- source
    }

    for (node in depthFirstTraverse(tree, postorder=F, source=root)) {
        if (node == root) {
            next
        } else {
            parent <- get_parent(tree, node)

            parent.state <- assignments[[as.character(parent)]]
            node.state <- assignments[[as.character(node)]]

            if (parent.state != node.state) {
                parsimony <- parsimony + 1
            }
        }
    }
    return(parsimony)
}
