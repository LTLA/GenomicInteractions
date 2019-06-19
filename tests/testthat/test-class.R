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
    expect_identical(X, unfactor(first(gi)))
    expect_identical(Y, unfactor(second(gi)))
    expect_identical(levels(first(gi)), levels(second(gi)))

    # Works when we don't ask for common regions.
    gi2 <- GenomicInteractions(X, Y, common=FALSE)
    expect_identical(X, unfactor(first(gi2)))
    expect_identical(Y, unfactor(second(gi2)))
    expect_identical(levels(first(gi2)), unique(X))
    expect_identical(levels(second(gi2)), unique(Y))

    # Works on empty inputs.
    expect_identical(length(GenomicInteractions(X[0], Y[0])), 0L)
    expect_identical(length(GenomicInteractions(X[0], Y[0], common=FALSE)), 0L)
})

test_that("GenomicInteractions() constructor works correctly with integer inputs", {
    reg <- spawn_regions(10)
    X <- sample(length(reg), 50, replace=TRUE)
    Y <- sample(length(reg), 50, replace=TRUE)

    gi <- GenomicInteractions(X, Y, reg)
    ref <- GenomicInteractions(reg[X], reg[Y])
    expect_identical(unfactor(first(gi)), unfactor(first(ref)))
    expect_identical(unfactor(second(gi)), unfactor(second(ref)))

    # Works on empty inputs.
    expect_identical(length(GenomicInteractions(X[0], Y[0], reg)), 0L)
})

test_that("GenomicInteractions() constructor deals correctly with mcols", {
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(15)
    X <- sample(reg1, 50, replace=TRUE)
    Y <- sample(reg2, 50, replace=TRUE)
    mcols(X)$blah <- 2
    mcols(Y)$blah <- 3

    gi <- GenomicInteractions(X, Y)
    expect_true(all(mcols(gi)$anchor1.blah==2))
    expect_true(all(mcols(gi)$anchor2.blah==3))
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
    o <- anchors(gi)
    expect_s4_class(o, "Pairs")
    expect_identical(unfactor(first(o)), X)
    expect_identical(unfactor(second(o)), Y)

    expect_identical(unfactor(anchors(gi, type=1)), X)
    expect_identical(unfactor(anchors(gi, type=2)), Y)
    expect_identical(unfactor(anchors(gi, type="first")), X)
    expect_identical(unfactor(anchors(gi, type="second")), Y)

    expect_identical(anchors(gi, type=1), first(gi))
    expect_identical(anchors(gi, type=2), second(gi))
})

test_that("anchor setters work correctly with regions", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(reg1, N, replace=TRUE)
    Y <- sample(reg2, N, replace=TRUE)

    gi <- spawn_gi(N)
    anchors(gi) <- Pairs(X, Y)
    expect_identical(unfactor(first(gi)), X)
    expect_identical(unfactor(second(gi)), Y)

    # Replacing one-by-one.
    gi <- spawn_gi(N)
    gi2 <- gi

    anchors(gi, 1) <- X
    expect_identical(unfactor(first(gi)), X)
    expect_identical(unfactor(second(gi)), unfactor(second(gi2)))

    anchors(gi, 2) <- Y
    expect_identical(unfactor(first(gi)), X)
    expect_identical(unfactor(second(gi)), Y)

    # Using the convenience accessors.
    gi <- spawn_gi(N)
    first(gi) <- X
    expect_identical(unfactor(first(gi)), X)
    second(gi) <- Y
    expect_identical(unfactor(second(gi)), Y)
})

test_that("anchor setters work correctly with indices", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(length(reg1), N, replace=TRUE)
    Y <- sample(length(reg2), N, replace=TRUE)
    original <- GenomicInteractions(Factor(index=X, levels=reg1), Factor(index=Y, levels=reg2))

    # Replacing all at once.
    gi <- original
    X2 <- sample(length(reg1), N, replace=TRUE)
    Y2 <- sample(length(reg2), N, replace=TRUE)

    anchors(gi, id=TRUE) <- DataFrame(X2, Y2)
    expect_identical(unfactor(first(gi)), reg1[X2])
    expect_identical(unfactor(second(gi)), reg2[Y2])

    # Replacing one by one.
    gi <- original
    X2 <- sample(length(reg1), N, replace=TRUE)
    Y2 <- sample(length(reg2), N, replace=TRUE)

    anchors(gi, id=TRUE, type=1) <- X2
    expect_identical(unfactor(first(gi)), reg1[X2])
    expect_identical(unfactor(second(gi)), reg2[Y])

    anchors(gi, id=TRUE, type=2) <- Y2
    expect_identical(unfactor(first(gi)), reg1[X2])
    expect_identical(unfactor(second(gi)), reg2[Y2])
})

####################################
####################################

test_that("region getters work correctly", {
    N <- 20
    reg1 <- spawn_regions(10)
    reg2 <- spawn_regions(5)
    X <- sample(length(reg1), N, replace=TRUE)
    Y <- sample(length(reg2), N, replace=TRUE)
    gi <- GenomicInteractions(Factor(index=X, levels=reg1), Factor(index=Y, levels=reg2))

    expect_identical(sort(regions(gi, type=1)), sort(reg1))
    expect_identical(sort(regions(gi, type=2)), sort(reg2))
    expect_identical(regions(gi, type="both"), 
        SimpleList(first=regions(gi,1), second=regions(gi,2)))

    expect_warning(old <- regions(gi), "deprecated")
    expect_identical(old, regions(gi, type=1))
})

test_that("region setters work correctly", {
    N <- 20
    reg1 <- unique(spawn_regions(10))
    reg2 <- unique(spawn_regions(5))
    X <- sample(length(reg1), N, replace=TRUE)
    Y <- sample(length(reg2), N, replace=TRUE)
    original <- GenomicInteractions(Factor(index=X, levels=reg1), Factor(index=Y, levels=reg2))

    # Replacing all elements at once.
    areg1 <- unique(spawn_regions(10))
    areg2 <- unique(spawn_regions(5))
    gi <- original
    regions(gi, type="both") <- SimpleList(areg1, areg2)

    expect_identical(regions(gi, type=1), areg1)
    expect_identical(unfactor(first(gi)), areg1[X])
    expect_identical(regions(gi, type=2), areg2)
    expect_identical(unfactor(second(gi)), areg2[Y])

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
