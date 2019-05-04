# Tests the linkOverlaps method
# library(GenomicInteractions); library(testthat); source("setup.R"); source("test-link.R")

set.seed(9000)

MANUAL <- function(obj, regions1, regions2, ..., ignore.strand=TRUE, use.region=c("any", "match")) {
    olap11 <- findOverlaps(first(obj), regions1, ..., ignore.strand=ignore.strand)
    olap22 <- findOverlaps(second(obj), regions2, ..., ignore.strand=ignore.strand)
    combo <- base::merge(olap11, olap22, by.x=1, by.y=1)

    use.region <- match.arg(use.region)
    if (use.region=="any") { 
        olap12 <- findOverlaps(first(obj), regions2, ..., ignore.strand=ignore.strand)
        olap21 <- findOverlaps(second(obj), regions1, ..., ignore.strand=ignore.strand)
        combo2 <- base::merge(olap21, olap12, by.x=1, by.y=1)
        combo <- rbind(combo, combo2)
    }

    combo <- DataFrame(combo)
    colnames(combo) <- c("query", "subject1", "subject2")
    rownames(combo) <- NULL

    # Obtaining a sorted and de-duplicated DataFrame.
    is.dup <- duplicated(sprintf("%s.%s.%s", combo$query, combo$subject1, combo$subject2))
    combo <- combo[!is.dup,]
    o <- order(combo$query, combo$subject1, combo$subject2)
    combo[o,]
}

test_that("linkOverlaps works correctly in the basic case", {
    x <- spawn_gi()
    gene.regions <- spawn_regions(10)
    enh.regions <- spawn_regions(8)

    expect_identical(MANUAL(x, gene.regions, enh.regions), 
        linkOverlaps(x, gene.regions, enh.regions))

    expect_identical(MANUAL(x, gene.regions, enh.regions, maxgap=10), 
        linkOverlaps(x, gene.regions, enh.regions, maxgap=10))

    expect_identical(MANUAL(x, gene.regions, enh.regions, minoverlap=10), 
        linkOverlaps(x, gene.regions, enh.regions, minoverlap=10))

    expect_identical(MANUAL(x, gene.regions, enh.regions, type="within"), 
        linkOverlaps(x, gene.regions, enh.regions, type="within"))
})

test_that("linkOverlaps works with matched hits", {
    x <- spawn_gi()
    gene.regions <- spawn_regions(15)
    enh.regions <- spawn_regions(10)
    expect_identical(MANUAL(x, gene.regions, enh.regions, use.region="match"),
        linkOverlaps(x, gene.regions, enh.regions, use.region="match"))

    # Checking consistency with use.region="any"
    ref <- linkOverlaps(x, gene.regions, enh.regions)
    left <- linkOverlaps(x, gene.regions, enh.regions, use.region="match")
    right <- linkOverlaps(x, enh.regions, gene.regions, use.region="match")

    right <- right[,c(1,3,2)]
    colnames(right) <- colnames(left)
    expect_identical(ref, unique(sort(rbind(left, right))))
})

test_that("linkOverlaps works correctly with self hits", {
    x <- spawn_gi()

    gene.regions <- spawn_regions(10)
    ref <- MANUAL(x, gene.regions, gene.regions)
    ref <- ref[ref$subject1 < ref$subject2,]
    expect_identical(linkOverlaps(x, gene.regions), ref)

    # Correctly avoid duplicates.
    out <- linkOverlaps(x, resize(gene.regions, 100))
    expect_false(any(duplicated(out)))
})

test_that("linkOverlaps works correctly with stranded information", {
    # By default, we just ignore the strand.
    x <- spawn_gi()
    gene.regions <- spawn_regions(15)
    enh.regions <- spawn_regions(10)

    ref <- linkOverlaps(x, gene.regions, enh.regions)
    strand(regions(x, 1)) <- sample(c("+", "-"), length(regions(x, 1)), replace=TRUE)
    strand(regions(x, 2)) <- sample(c("+", "-"), length(regions(x, 2)), replace=TRUE)
    expect_identical(ref, linkOverlaps(x, gene.regions, enh.regions))

    # Including the strand.
    alt <- linkOverlaps(x, gene.regions, enh.regions, ignore.strand=FALSE)
    expect_false(identical(alt, ref))
    expect_identical(alt, MANUAL(x, gene.regions, enh.regions, ignore.strand=FALSE))
})

test_that("linking overlaps behaves with empty inputs", {
    obj <- spawn_gi()
    gene.regions <- spawn_regions(10)

    expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), 
        linkOverlaps(obj, gene.regions[0], gene.regions))

    expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), 
        linkOverlaps(obj, gene.regions, gene.regions[0]))

    expect_identical(DataFrame(query=integer(0), subject1=integer(0), subject2=integer(0)), 
        linkOverlaps(obj, gene.regions[0], gene.regions[0]))
})
