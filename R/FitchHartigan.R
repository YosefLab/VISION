#' Compute parsimony with Fitch-Hartigan algorithm.
#' 
#' @param tree a rooted tree of the class `phylo`
#' @param tipMeta a list of labels for each tip of the tree
#' @return parsimony score
computeFitchHartiganParsimony <- function(tree, tipMeta) {

    if (!(is.character(tipMeta[,1]) | is.factor(tipMeta[,1]))) {
        stop('tipMeta must be categorical data')
    }
    ancestral_states <- bottomUpFitchHartigan(tree, tipMeta)

    assignments <- topDownFitchHartigan(tree, ancestral_states)

    parsimony <- countTransitions(tree, assignments)

    return(parsimony)
}

#' Runs bottom-up Fitch-Hartigan algorithm
#' 
#' @param tree a rooted tree of the class `phylo` with tip meta data
#' @param tipMeta a data frame of labels for each tip of the tree
#' @return a list mapping nodes to possible ancestral states
bottomUpFitchHartigan <- function(tree, tipMeta) {

    ancestral_states <- list() 

    postorder_traversal <- depth_first_traverse_nodes(tree)
    for (node in postorder_traversal) {
        if (node <= length(tree$tip.label)) {
            tip.name <- tree$tip.label[[node]]
            ancestral_states[as.character(node)] <- list(as.character(tipMeta[tip.name,]))
        } else {
            children <- get_children(tree, node)
            if (length(children) == 1) {
                ancestral_states[as.character(node)] <- ancestral_states[as.character(children[[1]])]
            } else{
                states_of_children <- unlist(lapply(children, function(x) ancestral_states[as.character(x)]))
                frequency_table <- table(states_of_children)
                labels <- names(frequency_table)[as.numeric(which(frequency_table == max(frequency_table)))]
                ancestral_states[[as.character(node)]] <- labels
            }
        }
    }
    return(ancestral_states)
}

#' Runs top-down Fitch-Hartigan algorithm
#' 
#' @importFrom phytools getParent
#' @param tree a rooted tree of the class `phylo`
#' @param ancestral_states a list mapping each node to a set of possible ancestral states
#' @return a list mapping each node to an assignment that satisfies the maximum parsimony criterion
topDownFitchHartigan <- function(tree, ancestral_states) {

    assignments <- list()
    preorder_traversal <- depth_first_traverse_nodes(tree, postorder=F)

    for (node in preorder_traversal) {
        
        possible_states <- ancestral_states[[as.character(node)]]
        if (find_root(tree) != node) {
            parent <- getParent(tree, node)
            parent_state <- assignments[[as.character(parent)]]
            intersection <- intersect(parent_state, possible_states)
            if (length(intersection) > 0) {
                possible_states <- c(parent_state)
            }
        }
        assignments[as.character(node)] <- sample(possible_states, 1)   
    }
    return(assignments)
}

#' Counts transitions from Fitch-Hartigan assignments
#' 
#' @param tree a rooted tree of the class `phylo`
#' @param assignments Assignments for each node in the tree
#' @return parsimony score
countTransitions <- function(tree, assignments) {
    
    score <- 0
    preorder <- depth_first_traverse_nodes(tree, postorder=F)

    for (node in preorder) {
        

        if (find_root(tree) != node) {
            parent <- getParent(tree, node)
            parent_state <- assignments[[as.character(parent)]]
            node_state <- assignments[[as.character(node)]]
            score <- score + as.integer(!(parent_state == node_state))
        }
    }
    return(score)
}