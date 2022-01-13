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
    rate <- computeFitchHartiganParsimony(tree, metaData)
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