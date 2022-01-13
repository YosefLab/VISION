context("Single cell plasticity")

test_that("Fitch bottom up algorithm works on binary small input", {

    nwk <- "( ((1, 2)7, 3)8, ((4,5)9, 6)10)11;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "B", "4" = "C", "5" = "A", "6" = "C")

    possibleLabels <- bottomUpFitchHartigan(tree, metaData)
    expected_labels <-  list(
            "1" = c("A"), "2" = c("A"), "3" = c("B"), "4" = c("C"), "5" = c("A"), "6" = c("C"),
            "7" = c("A"), "8" = c("A", "B"), "10" = c("C"), "9" = c("A", "C"), "11" = c("A", "B", "C")
        )

    # check leaves are assigned correctly
    for (i in 1:6) {
        observed_state = possibleLabels[[as.character(i)]]
        expected_state = expected_labels[[tree$tip.label[[i]]]]
        expect_equal(observed_state, expected_state)
    }

    # test internal nodes are assigned correctly
    for (i in 7:11) {
        observed_state = possibleLabels[[as.character(i)]]
        expected_state = expected_labels[[tree$node.label[[i - 6]]]]
        expect_equal(observed_state, expected_state)
    }

})


test_that("Fitch top down algorithm works on binary small input", {

    nwk <- "( ((1, 2)7, 3)8, ((4,5)9, 6)10)11;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "B", "4" = "B", "5" = "B", "6" = "B")

    possibleLabels <- bottomUpFitchHartigan(tree, metaData)
    assignments <- topDownFitchHartigan(tree, possibleLabels)

    # check leaves are assigned correctly
    for (i in 1:6) {
        observed_state = assignments[[as.character(i)]]
        expected_state = metaData[[tree$tip.label[[i]]]]
        expect_equal(observed_state, expected_state)
    }

    # check internal nodes are accurately assigned
    expected_assignments = list(
        "11" = "B", "10" = "B", "9" = "B", "8" = "B", "7" = "A"
    )
    for (i in 7:11) {
        observed_state = assignments[[as.character(i)]]
        expected_state = expected_assignments[[tree$node.label[[i - 6]]]]
        expect_equal(observed_state, expected_state)
    }

})

test_that("Fitch-Hartigan parsimony works on small example", {

    nwk <- "( ((1, 2)7, 3)8, ((4,5)9, 6)10)11;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "B", "4" = "C", "5" = "A", "6" = "C")

    possibleLabels <- bottomUpFitchHartigan(tree, metaData)
    assignments <- topDownFitchHartigan(tree, possibleLabels)
    score <- scoreParsimony(tree, assignments)

    expect_equal(score, 3)

    # test normalized score
    rate <- computeNormalizedFitchHartiganParsimony(tree, metaData)
    expect_equal(rate, 0.3)

})

test_that("Fitch-Hartigan rates per node", {

    nwk <- "( ((1, 2)7, 3)8, ((4,5)9, 6)10)11;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "B", "4" = "C", "5" = "A", "6" = "C")

    score_per_node <- computeFitchHartiganParsimonyPerNode(tree, metaData)
    
    expected_scores = list("7" = 0, "8" = 0.25, "9" = 0.5, "10" = 0.25, "11" = 0.3)
    for (i in 7:11) {
        observed_score = score_per_node[[as.character(i)]]
        expected_score = expected_scores[[tree$node.label[[i - 6]]]]
        expect_equal(observed_score, expected_score)
    }

})

test_that("Single cell plasticities on an basic tree", {

    nwk <- "( ((1, 2)7, 3)8, ((4,5)9, 6)10)11;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "B", "4" = "C", "5" = "A", "6" = "C")

    score_per_node <- computeFitchHartiganParsimonyPerNode(tree, metaData)
    leaf_plasticities <- computeSingleCellFitchScores(tree, score_per_node)

    expected_scores = list("1" = 0.183, "2" = 0.183, "3" = 0.275, "4" = 0.35, "5" = 0.35, "6" = 0.275)
    for (leafName in tree$tip.label) {
        observed_score <- leaf_plasticities[[leafName]]
        expected_score <- expected_scores[[leafName]]
        expect_equal(observed_score, expected_score, tolerance=0.3)
    }

})

test_that("Parsimony scores work on a general tree", {

    nwk <- "(((1, 2)11, 3, (4, 5)13)12, (6, (7, (8, 9, 10)14)15)16)17;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "C", "4" = "D", "5" = "D", "6" = "B", "7" = "B", "8" = "C", "9" = "C", "10" = "B")
    possibleLabels <- bottomUpFitchHartigan(tree, metaData)
    assignments <- topDownFitchHartigan(tree, possibleLabels)
    parsimony <- scoreParsimony(tree, assignments)

    expect_equal(parsimony, 5)

    rate <- computeNormalizedFitchHartiganParsimony(tree, metaData)
    expect_equal(rate, 5/16)

})

test_that("Single cell plasticities work on a general tree", {

    nwk <- "(((1, 2)11, 3, (4, 5)13)12, (6, (7, (8, 9, 10)14)15)16)17;"
    tree <- read.tree(text=nwk)

    metaData <- list("1" = "A", "2" = "A", "3" = "C", "4" = "D", "5" = "D", "6" = "B", "7" = "B", "8" = "C", "9" = "C", "10" = "B")

    node.scores <- computeFitchHartiganParsimonyPerNode(tree, metaData)
    expected.scores <- list("11" = 0, "12" = 2/7, "13" = 0, "14" = 1/3, "15" = 0.4, "16" = 2/7, "17" = 5/16)
    for (i in 11:17) {
        observed_score <- node.scores[[as.character(i)]]
        internal_node_name <- tree$node.label[[i - length(tree$tip.label)]]
        expected_score <- expected.scores[[internal_node_name]]
        expect_equal(observed_score, expected_score)
    }

    plasticities <- computeSingleCellFitchScores(tree, node.scores)
    expected_plasticities <- list(
        "1" = 0.199, "2" = 0.199, "3" = 0.299,
        "4" = 0.199, "5" = 0.199, "6" = 0.299, "7" = 0.333, "8" = 0.332,
        "9" = 0.332, "10" = 0.332
    )
    for (leafName in tree$tip.label) {
        expected_score <- expected_plasticities[[leafName]]
        observed_score <- plasticities[[leafName]]
        expect_equal(observed_score, expected_score, tolerance=3)
    }
})