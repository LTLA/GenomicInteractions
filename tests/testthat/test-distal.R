# Testing the findDistalAnchors function.
# library(GenomicInteractions); library(testthat); source("setup.R"); source("test-distal.R")

set.seed(200)

test_that("identification of distal anchors works correctly", {
    x <- spawn_gi()
    combined <- unlist(regions(x, type="both"))

    for (i in seq_along(combined)) {
        distals <- findDistalAnchors(x, combined[i])
        olap1 <- overlapsAny(first(x), combined[i])
        olap2 <- overlapsAny(second(x), combined[i])

        o <- c(which(olap2), which(olap1))
        ref <- c(first(x)[olap2], second(x)[olap1])
        expect_identical(distals, unfactor(ref[order(o)]))
    }
})

test_that("first anchors go before second anchors when local=TRUE", {
    x <- GenomicInteractions(GRanges("chrA:10-12"), GRanges("chrA:1-5"))
    out <- findDistalAnchors(x, GRanges("chrA:1-20"))
    expect_identical(out, unfactor(c(first(x), second(x))))
})

test_that("identification of distal anchors works correctly when local=FALSE", {
    x <- spawn_gi()
    combined <- unlist(regions(x, type="both"))

    for (i in seq_along(combined)) {
        distals <- findDistalAnchors(x, combined[i], local=FALSE)
        olap1 <- overlapsAny(first(x), combined[i])
        olap2 <- overlapsAny(second(x), combined[i])

        keep1 <- olap1 & !olap2
        keep2 <- olap2 & !olap1
        o <- c(which(keep2), which(keep1))
        ref <- c(first(x)[keep2], second(x)[keep1])
        expect_identical(distals, unfactor(ref[order(o)]))
    }
})

test_that("findDistalAnchors passes arguments along correctly", {
    x <- spawn_gi()

    thing <- GRanges("chrA:1-20")
    distal1 <- findDistalAnchors(x, thing, maxgap=100)
    olap1 <- overlapsAny(first(x), thing, maxgap=100)
    olap2 <- overlapsAny(second(x), thing, maxgap=100)

    o <- c(which(olap2), which(olap1))
    ref <- c(first(x)[olap2], second(x)[olap1])
    expect_identical(distal1, unfactor(ref[order(o)]))

    distal2 <- findDistalAnchors(x, thing)
    expect_false(identical(distal1, distal2))
})

test_that("findDistalAnchors handles metadata correctly", {
    x <- spawn_gi()
    thing <- GRanges("chrA:1-20")

    # Doesn't like it when metadata fields are not consistent.
    mcols(levels(first(x)))$thing <- 10
    expect_warning(out <- findDistalAnchors(x, thing), "mismatching")
    expect_identical(ncol(mcols(out)), 0L)

    # Becomes a little happier when metadata fields are present in both regions.
    mcols(levels(second(x)))$thing <- 20
    expect_warning(out <- findDistalAnchors(x, thing), NA)

    F1 <- unfactor(first(x), ignore.mcols=TRUE)
    F2 <- unfactor(second(x), ignore.mcols=TRUE)
    olap1 <- overlapsAny(F1, thing)
    olap2 <- overlapsAny(F2, thing)
    o <- c(which(olap2), which(olap1))
    expect_identical(out$thing, c(F1$thing[olap2], F2$thing[olap1])[order(o)])

    # Includes metadata from 'x' itself.
    mcols(x)$stuff <- seq_along(x)
    out <- findDistalAnchors(x, thing)
    expect_identical(out$stuff, sort(o))
})

test_that("linearization behaves correctly with empty inputs", {
    x <- spawn_gi()
    expect_identical(length(findDistalAnchors(x[0,], GRanges())), 0L)
    expect_identical(length(suppressWarnings(findDistalAnchors(x, GRanges("chrC", IRanges(1,1))))), 0L)
})
