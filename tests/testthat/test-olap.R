# Testing the various overlap methods for GenomicInteractions objects.
# library(GenomicInteractions); library(testthat); source("setup.R"); source("test-olap.R")

expect_identical_hits <- function(x, y) {
    expect_identical(sort(x), sort(y))
}

set.seed(300)

########################################
########################################

CHECK_HITS_LEFT <- function(x, y, ...) {
    expect_identical_hits(findOverlaps(x, y, ..., use.region="first"), olap1 <- findOverlaps(first(x), y, ...))
    expect_identical_hits(findOverlaps(x, y, ..., use.region="second"), olap2 <- findOverlaps(second(x), y, ...))
    expect_identical_hits(ref <- findOverlaps(x, y, ...), union(olap1, olap2))
    expect_false(is.unsorted(queryHits(ref)))
}

CHECK_HITS_RIGHT <- function(x, y, ...) {
    expect_identical_hits(findOverlaps(x, y, ..., use.region="first"), olap1 <- findOverlaps(x, first(y), ...))
    expect_identical_hits(findOverlaps(x, y, ..., use.region="second"), olap2 <- findOverlaps(x, second(y), ...))
    expect_identical_hits(ref <- findOverlaps(x, y, ...), union(olap1, olap2))
    expect_false(is.unsorted(queryHits(ref)))
}

CHECK_SELECT <- function(x, y, ..., select) {
    out <- findOverlaps(x, y, ..., select=select)
    expect_identical(length(out), length(x))
    expect_type(out, "integer")
}

test_that("findOverlaps works with 1D overlaps (LEFT)", {
    x <- spawn_gi()
    y <- spawn_regions()

    # Basic checks.
    CHECK_HITS_LEFT(x, y)
    CHECK_HITS_LEFT(x, y, maxgap=20)
    CHECK_HITS_LEFT(x, y, minoverlap=10)
    CHECK_HITS_LEFT(x, y, type="within")

    # Respects 'select'.
    CHECK_SELECT(x, y, select="first")
    CHECK_SELECT(x, y, select="last")
    CHECK_SELECT(x, y, select="arbitrary")
})

test_that("findOverlaps works with 1D overlaps (RIGHT)", {
    x <- spawn_regions()
    y <- spawn_gi()

    # Basic checks.
    CHECK_HITS_RIGHT(x, y)
    CHECK_HITS_RIGHT(x, y, maxgap=20)
    CHECK_HITS_RIGHT(x, y, minoverlap=10)
    CHECK_HITS_RIGHT(x, y, type="within")

    # Respects 'select'.
    CHECK_SELECT(x, y, select="first")
    CHECK_SELECT(x, y, select="last")
    CHECK_SELECT(x, y, select="arbitrary")
})

test_that("1D overlaps respect strandedness", {
    x <- spawn_gi()
    y <- spawn_regions()

    lref <- findOverlaps(x, y)
    rref <- findOverlaps(y, x)

    strand(regions(x, 1)) <- sample(c("+", "-"), length(regions(x, 1)), replace=TRUE)
    strand(regions(x, 2)) <- sample(c("+", "-"), length(regions(x, 2)), replace=TRUE)

    # From the left:
    alt <- findOverlaps(x, y)
    expect_identical(lref, alt)
    alt2 <- findOverlaps(x, y, ignore.strand=FALSE)
    expect_false(identical(lref, alt2))

    CHECK_HITS_LEFT(x, y, ignore.strand=FALSE)

    # From the right:
    alt <- findOverlaps(y, x)
    expect_identical(rref, alt)
    alt2 <- findOverlaps(y, x, ignore.strand=FALSE)
    expect_false(identical(rref, alt2))

    CHECK_HITS_RIGHT(y, x, ignore.strand=FALSE)
})

test_that("1D overlaps behave with empty inputs", {
    x <- spawn_gi()
    y <- spawn_regions()
    expect_identical(length(findOverlaps(x[0], y)), 0L)
    expect_identical(length(findOverlaps(x, y[0])), 0L)
    expect_identical(length(findOverlaps(x[0], y[0])), 0L)
})

########################################
########################################

CHECK_2D_HITS <- function(x, y, ...) {
    olap11 <- findOverlaps(first(x), first(y), ...)
    olap22 <- findOverlaps(second(x), second(y), ...)
    olapA <- BiocGenerics::intersect(olap11, olap22)

    olap12 <- findOverlaps(first(x), second(y), ...)
    olap21 <- findOverlaps(second(x), first(y), ...)
    olapB <- BiocGenerics::intersect(olap12, olap21)

    expect_identical_hits(findOverlaps(x, y, ...), union(olapA, olapB))
    expect_identical_hits(findOverlaps(x, y, ..., use.region="match"), olapA)
    expect_identical_hits(findOverlaps(x, y, ..., use.region="reverse"), olapB)
}

test_that("findOverlaps works with 2D overlaps", {
    x <- spawn_gi()
    y <- spawn_gi()

    CHECK_2D_HITS(x, y)
    CHECK_2D_HITS(x, y, type="within")
    CHECK_2D_HITS(x, y, maxgap=20)
    CHECK_2D_HITS(x, y, minoverlap=10)

    # Respects 'select'.
    CHECK_SELECT(x, y, select="first")
    CHECK_SELECT(x, y, select="last")
    CHECK_SELECT(x, y, select="arbitrary")
})

test_that("findOverlaps works with 2D overlaps", {
    x <- spawn_gi()
    y <- spawn_gi()
    ref <- findOverlaps(x, y)

    strand(regions(x, 1)) <- sample(c("+", "-"), length(regions(x, 1)), replace=TRUE)
    strand(regions(x, 2)) <- sample(c("+", "-"), length(regions(x, 2)), replace=TRUE)
    strand(regions(y, 1)) <- sample(c("+", "-"), length(regions(y, 1)), replace=TRUE)
    strand(regions(y, 2)) <- sample(c("+", "-"), length(regions(y, 2)), replace=TRUE)

    alt <- findOverlaps(x, y)
    expect_identical(ref, alt)
    alt2 <- findOverlaps(x, y, ignore.strand=FALSE)
    expect_false(identical(ref, alt2))

    CHECK_2D_HITS(x, y, ignore.strand=FALSE)
    CHECK_2D_HITS(x, y, ignore.strand=FALSE, type="within")
    CHECK_2D_HITS(x, y, ignore.strand=FALSE, maxgap=20)
    CHECK_2D_HITS(x, y, ignore.strand=FALSE, minoverlap=10)
})

test_that("2D overlaps behave with empty inputs", {
    x <- spawn_gi()
    y <- spawn_gi()
    expect_identical(length(findOverlaps(x[0], y)), 0L)
    expect_identical(length(findOverlaps(x, y[0])), 0L)
    expect_identical(length(findOverlaps(x[0], y[0])), 0L)
})
