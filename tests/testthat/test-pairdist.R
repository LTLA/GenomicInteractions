# Tests the distance calculation methods for an InteractionSet.
# library(GenomicInteractions); library(testthat); source("setup.R"); source("test-pairdist.R")

set.seed(100)

test_that("Paired distances are correct for GI objects", {
    x <- spawn_gi()
    a1 <- first(x)
    a2 <- second(x)
    
    is.intra <- intrachr(x)
    expect_identical(is.intra, as.logical(seqnames(a1)==seqnames(a2)))
    expect_identical(is.intra, pairdist(x, type="intra"))

    # Don't use 'mid', it does its own rounding.
    expect_identical(pairdist(x), 
        ifelse(is.intra, abs(start(a1)+end(a1)-start(a2)-end(a2))/2L, NA_real_)) 

    expect_identical(pairdist(x,type="gap"), 
        ifelse(is.intra, pmax(start(a1), start(a2)) - pmin(end(a1), end(a2)) -1L, as.integer(NA))) 

    expect_identical(pairdist(x,type="span"), 
        ifelse(is.intra, pmax(end(a1), end(a2)) - pmin(start(a1), start(a2)) +1L, as.integer(NA))) 
})

test_that("Diagonal extraction works as expected", {
    x <- spawn_gi()
    expect_error(pairdist(x, type="diag"), "same region set")

    y <- GenomicInteractions(unfactor(first(x)), unfactor(second(x)))
    ax1 <- anchors(y, type=1, id=TRUE)
    ax2 <- anchors(y, type=2, id=TRUE)
    expect_identical(old <- pairdist(y, type="diag"), ifelse(intrachr(y), abs(ax1-ax2), as.integer(NA)))

    # Handles unordered output in both regions.
    regions <- regions(y, 1)
    r <- sample(length(regions))

    converter <- integer(length(r))
    converter[r] <- seq_along(r)

    y <- GenomicInteractions(converter[ax1], converter[ax2], regions[r])
    expect_identical(pairdist(y, type="diag"), old) 

    # Handles unordered output in one region.
    y <- GenomicInteractions(
        Factor(index=ax1, levels=regions), 
        Factor(index=converter[ax2], levels=regions[r])
    )
    expect_identical(pairdist(y, type="diag"), old) 
}) 

test_that("pairdist behaves with empty inputs", {    
    x <- spawn_gi()
    expect_identical(pairdist(x[0,]), numeric(0))
    expect_identical(pairdist(x[0,], type="intra"), logical(0))

    is.intra <- intrachr(x)
    expect_identical(pairdist(x[!is.intra,]), rep(NA_real_, sum(!is.intra)))
})
