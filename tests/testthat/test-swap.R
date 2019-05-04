# This tests the swapAnchors() function.
# library(testthat); library(GenomicInteractions); source("setup.R"); source("test-swap.R")

test_that("default anchor swapping works", {
    x <- spawn_gi()
    new.x <- swapAnchors(x)
    expect_true(all(first(new.x) <= second(new.x)))
    expect_identical(partnerNames(x), partnerNames(new.x))

    left <- first(x)
    right <- second(x)
    expect_identical(first(new.x), pmin(left, right))
    expect_identical(second(new.x), pmax(left, right))
})

test_that("reversed anchor swapping works", {
    x <- spawn_gi()
    new.x <- swapAnchors(x, mode="reverse")
    expect_true(all(first(new.x) >= second(new.x)))
    expect_identical(partnerNames(x), partnerNames(new.x))

    left <- first(x)
    right <- second(x)
    expect_identical(first(new.x), pmax(left, right))
    expect_identical(second(new.x), pmin(left, right))
})

test_that("simple anchor swapping works", {
    x <- spawn_gi()
    new.x <- swapAnchors(x, mode="all")
    expect_identical(first(new.x), second(x))
    expect_identical(second(new.x), first(x))
    expect_identical(partnerNames(x), partnerNames(new.x))
})

test_that("anchor swapping works with empty inputs", {
    x <- spawn_gi()
    expect_as_if(x[0], swapAnchors(x[0]))
    expect_as_if(x[0], swapAnchors(x[0], mode="reverse"))
    expect_as_if(x[0], swapAnchors(x[0], mode="all"))
})

