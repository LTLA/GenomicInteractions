# This tests the constructor, getters and setters for the GenomicInteractions class.
# library(testthat); library(GenomicInteractions); source("setup.R"); source("test-class.R")

set.seed(199999)

####################################
####################################

test_that("GenomicInteractions() constructor works correctly with GRanges", {
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(15)
    X <- sample(reg1, 50, replace=TRUE)
    Y <- sample(reg2, 50, replace=TRUE)

    gi <- GenomicInteractions(X, Y)
    expect_identical(featureSets(gi)[[1]], featureSets(gi)[[2]])
    expect_identical(X, partnerFeatures(gi, 1))
    expect_identical(Y, partnerFeatures(gi, 2))

    # Works when we don't ask for common regions.
    gi2 <- GenomicInteractions(X, Y, common=FALSE)
    expect_identical(featureSets(gi2)[[1]], unique(X))
    expect_identical(featureSets(gi2)[[2]], unique(Y))
    expect_identical(X, partnerFeatures(gi2, 1))
    expect_identical(Y, partnerFeatures(gi2, 2))

    # Works on empty inputs.
    expect_identical(length(GenomicInteractions(X[0], Y[0])), 0L)
    expect_identical(length(GenomicInteractions(X[0], Y[0], common=FALSE)), 0L)
})

test_that("GenomicInteractions() constructor works correctly with integer inputs", {
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(8)
    X <- sample(length(reg1), 50, replace=TRUE)
    Y <- sample(length(reg2), 50, replace=TRUE)

    gi <- GenomicInteractions(X, Y, reg1)
    expect_identical(featureSets(gi)[[1]], reg1)
    expect_identical(featureSets(gi)[[2]], reg1)
    expect_identical(reg1[X], partnerFeatures(gi, 1))
    expect_identical(reg1[Y], partnerFeatures(gi, 2))

    # Works with GRangesList inputs.
    gi2 <- GenomicInteractions(X, Y, List(reg1, reg2))
    expect_identical(featureSets(gi2)[[1]], reg1)
    expect_identical(featureSets(gi2)[[2]], reg2)
    expect_identical(reg1[X], partnerFeatures(gi2, 1))
    expect_identical(reg2[Y], partnerFeatures(gi2, 2))

    # Works on empty inputs.
    expect_identical(length(GenomicInteractions(X[0], Y[0], reg1)), 0L)
    expect_identical(length(GenomicInteractions(X[0], Y[0], List(reg1, reg2))), 0L)
})

test_that("validity check triggers for invalid object", {
    reg <- IRanges(1:10, 1:10)
    out <- IndexedRelations(list(1, 1), List(reg, reg))
    expect_error(as(out, "GenomicInteractions"), "not valid")

    greg <- spawn_regions(10)
    out <- IndexedRelations(list(1, 1, 1), List(greg, greg, greg))
    expect_error(validObject(as(out, "GenomicInteractions")), "must contain pairwise")
})

####################################
####################################

test_that("anchor getters work correctly", {
    # Not using spawn_gi(), because we want to compare to 'X' and 'Y'.
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(15)
    X <- sample(reg1, 50, replace=TRUE)
    Y <- sample(reg2, 50, replace=TRUE)
    gi <- GenomicInteractions(X, Y)

    p <- Pairs(X, Y)
    mcols(p) <- NULL
    expect_identical(anchors(gi), p)
    expect_identical(anchors(gi, type=1), X)
    expect_identical(anchors(gi, type=2), Y)
    expect_identical(anchors(gi, type="first"), X)
    expect_identical(anchors(gi, type="second"), Y)

    expect_identical(anchors(gi, type=1), first(gi))
    expect_identical(anchors(gi, type=2), second(gi))

    expect_identical(anchors(gi, id=TRUE), partners(gi))
    expect_identical(anchors(gi, id=TRUE, type=1), partners(gi)[,1])
    expect_identical(anchors(gi, id=TRUE, type=2), partners(gi)[,2])
    expect_identical(anchors(gi, id=TRUE, type="first"), partners(gi)[,1])
    expect_identical(anchors(gi, id=TRUE, type="second"), partners(gi)[,2])

    expect_warning(out <- anchors(gi, type="both"), "deprecated")
    expect_identical(out, anchors(gi, type=NULL))
})

test_that("anchor setters work correctly with regions", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(reg1, N, replace=TRUE)
    Y <- sample(reg2, N, replace=TRUE)

    gi <- spawn_gi(N)
    anchors(gi) <- Pairs(X, Y)
    expect_identical(first(gi), X)
    expect_identical(second(gi), Y)

    # Replacing one-by-one.
    gi <- spawn_gi(N)
    gi2 <- gi

    anchors(gi, 1) <- X
    expect_identical(first(gi), X)
    expect_identical(second(gi), second(gi2))

    anchors(gi, 2) <- Y
    expect_identical(first(gi), X)
    expect_identical(second(gi), Y)

    # Using the convenience accessors.
    gi <- spawn_gi(N)
    first(gi) <- X
    expect_identical(first(gi), X)
    second(gi) <- Y
    expect_identical(second(gi), Y)

    # Checking deprecation.
    gi <- spawn_gi(N)
    expect_warning(anchors(gi, type="both") <- Pairs(X, Y), "deprecated")
    expect_identical(anchors(gi, 1), X)
    expect_identical(anchors(gi, 2), Y)
})

test_that("anchor setters work correctly with indices", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(length(reg1), N, replace=TRUE)
    Y <- sample(length(reg2), N, replace=TRUE)
    original <- GenomicInteractions(X, Y, List(reg1, reg2))

    # Replacing all at once.
    gi <- original
    X2 <- sample(length(reg1), N, replace=TRUE)
    Y2 <- sample(length(reg2), N, replace=TRUE)

    anchors(gi, id=TRUE) <- DataFrame(X2, Y2)
    expect_identical(first(gi), reg1[X2])
    expect_identical(second(gi), reg2[Y2])

    # Replacing one by one.
    gi <- original
    X2 <- sample(length(reg1), N, replace=TRUE)
    Y2 <- sample(length(reg2), N, replace=TRUE)

    anchors(gi, id=TRUE, type=1) <- X2
    expect_identical(first(gi), reg1[X2])
    expect_identical(second(gi), reg2[Y])

    anchors(gi, id=TRUE, type=2) <- Y2
    expect_identical(first(gi), reg1[X2])
    expect_identical(second(gi), reg2[Y2])
})

####################################
####################################

test_that("region getters work correctly", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(length(reg1), N, replace=TRUE)
    Y <- sample(length(reg2), N, replace=TRUE)
    gi <- GenomicInteractions(X, Y, List(reg1, reg2))

    expect_identical(regions(gi, type=NULL), List(first=reg1, second=reg2))
    expect_identical(regions(gi, type=1), reg1)
    expect_identical(regions(gi, type=2), reg2)

    expect_warning(old <- regions(gi), "deprecated")
    expect_identical(old, reg1)
})

test_that("region setters work correctly", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(length(reg1), N, replace=TRUE)
    Y <- sample(length(reg2), N, replace=TRUE)
    original <- GenomicInteractions(X, Y, List(reg1, reg2))

    # Replacing all elements at once.
    areg1 <- spawn_regions(10)
    areg2 <- spawn_regions(5)
    gi <- original
    regions(gi, type=NULL) <- List(areg1, areg2)

    expect_identical(regions(gi, type=1), areg1)
    expect_identical(first(gi), areg1[X])
    expect_identical(regions(gi, type=2), areg2)
    expect_identical(second(gi), areg2[Y])

    # Going one by one.
    gi <- original
    regions(gi, type=1) <- areg1
    expect_identical(regions(gi, type=1), areg1)
    expect_identical(regions(gi, type=2), reg2)

    regions(gi, type=2) <- areg2
    expect_identical(regions(gi, type=1), areg1)
    expect_identical(regions(gi, type=2), areg2)
    
    # Dep warning triggered.
    gi <- original
    expect_warning(regions(gi) <- areg1, "deprecated")
    expect_identical(regions(gi, type=1), areg1)
    expect_identical(regions(gi, type=2), areg1)
})
