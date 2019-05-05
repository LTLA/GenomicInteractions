# Testing the bounding box behaviour.
# library(GenomicInteractions); library(testthat); source("setup.R"); source("test-box.R")

set.seed(10000)

MANUAL <- function(x, f, reflect=TRUE) {
    if (reflect) x <- swapAnchors(x)
    ref1 <- unlist(range(S4Vectors::split(first(x), f)))
    ref2 <- unlist(range(S4Vectors::split(second(x), f)))
    ref <- GenomicInteractions(unname(ref1), unname(ref2), common=reflect)
    names(ref) <- names(ref1)
    ref
}

test_that("boundingBox works in default settings", {
    for (i in 1:3) { 
        x <- spawn_gi()
        f0 <- paste0(seqnames(first(x)), seqnames(second(x)))
        f <- paste0(f0, sample(i, length(x), replace=TRUE))
        expect_identical(MANUAL(x, f), boundingBox(x, f))

        # Ensure that inter-chromosomals are reflected correctly.
        f0[!intrachr(x)] <- "inter"
        f <- paste0(f0, sample(i, length(x), replace=TRUE))
        expect_identical(MANUAL(x, f), boundingBox(x, f))
    }
})

test_that("boundingBox works when reflection is disabled", {
    for (i in 1:3) { 
        x <- spawn_gi()
        f <- paste0(seqnames(first(x)), seqnames(second(x)), sample(i, length(x), replace=TRUE))
        ref <- MANUAL(x, f, reflect=FALSE)
        expect_identical(ref, boundingBox(x, f, reflect=FALSE))

        # Checking that it actually does give different results.
        expect_false(all(ref==boundingBox(x, f)))
    }
})

test_that("boundingBox works without a factor", {
    for (left in c("chrA", "chrB")) {
        for (right in c("chrA", "chrB")) {
            x <- spawn_gi()
            fchr <- seqnames(first(x))
            schr <- seqnames(second(x))

            only <- as.logical((fchr==left & schr==right) | (fchr==right & schr==left))
            x <- x[only]

            y <- swapAnchors(x)
            ref1 <- range(first(y))
            ref2 <- range(second(y))
            ref <- GenomicInteractions(unname(ref1), unname(ref2))
            names(ref) <- 1
            expect_identical(boundingBox(x), ref)
        }
    }
})

test_that("disabling reflection breaks with mixed chromosomes", {
    x <- GenomicInteractions(
        GRanges(c("chrA:1-1", "chrB:2-2")),
        GRanges(c("chrB:1-1", "chrA:2-2"))
    )

    expect_error(boundingBox(x, reflect=FALSE), "multiple chromosomes")

    out <- boundingBox(x)
    expect_identical(out, MANUAL(x, f=c(1,1)))
})

test_that("boundingBox reports correct seqinfo", {
    x <- GenomicInteractions(
        GRanges(c("chrA:1-1", "chrA:2-2")),
        GRanges(c("chrA:1-1", "chrB:2-2")),
        common=FALSE
    )

    out <- boundingBox(x, c(1,2), reflect=FALSE)
    expect_identical(seqlevels(first(out)), "chrA")
    expect_identical(seqlevels(second(out)), c("chrA", "chrB"))
    expect_false(identical(regions(out, 1), regions(out, 2)))

    out <- boundingBox(x, c(1,2), reflect=TRUE)
    expect_identical(seqlevels(first(out)), c("chrA", "chrB"))
    expect_identical(seqlevels(second(out)), c("chrA", "chrB"))
    expect_identical(regions(out, 1), regions(out, 2))
})

test_that("boundingBox breaks with silly inputs", {
    x <- spawn_gi()
    expect_error(boundingBox(x), "multiple chromosomes for group '1'")
    f <- rep("whee", length(x))
    f[-1] <- "YAY"
    expect_error(boundingBox(x,f), "multiple chromosomes for group 'YAY'")
    
    ref <- GenomicInteractions(first(x)[0], second(x)[0])
    names(ref) <- character(0)
    expect_as_if(boundingBox(x[0]), ref)
})
