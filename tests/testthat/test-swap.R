# This tests the swapAnchors() function.
# library(testthat); library(GenomicInteractions); source("setup.R"); source("test-swap.R")

test_that("default anchor swapping works", {
    x <- spawn_gi()
    new.x <- swapAnchors(x)
    expect_true(all(first(new.x) <= second(new.x)))

    left <- unfactor(first(x))
    right <- unfactor(second(x))
    expect_identical(unfactor(first(new.x)), pmin(left, right))
    expect_identical(unfactor(second(new.x)), pmax(left, right))
})

test_that("reversed anchor swapping works", {
    x <- spawn_gi()
    new.x <- swapAnchors(x, mode="reverse")
    expect_true(all(first(new.x) >= second(new.x)))

    left <- unfactor(first(x))
    right <- unfactor(second(x))
    expect_identical(unfactor(first(new.x)), pmax(left, right))
    expect_identical(unfactor(second(new.x)), pmin(left, right))
})

test_that("simple anchor swapping works", {
    x <- spawn_gi()
    new.x <- swapAnchors(x, mode="all")
    expect_identical(unfactor(first(new.x)), unfactor(second(x)))
    expect_identical(unfactor(second(new.x)), unfactor(first(x)))
})

test_that("anchor swapping works with empty inputs", {
    x <- spawn_gi()
    expect_as_if(x[0], swapAnchors(x[0]))
    expect_as_if(x[0], swapAnchors(x[0], mode="reverse"))
    expect_as_if(x[0], swapAnchors(x[0], mode="all"))
})

