# library(testthat); library(GenomicInteractions)

spawn_gi <- function(N1=30, N2=15, Np=50) {
    all.starts <- round(runif(N1, 1, 100))
    all.ends <- all.starts + round(runif(N1, 5, 20))
    regions1 <- GRanges(sample(c("chrA", "chrB"), N1, replace=TRUE), IRanges(all.starts, all.ends),
        seqlengths=c(chrA=1000, chrB=2000))

    all.starts <- round(runif(N2, 1, 100))
    all.ends <- all.starts + round(runif(N2, 5, 20))
    regions2 <- GRanges(sample(c("chrA", "chrB"), N2, replace=TRUE), IRanges(all.starts, all.ends),
        seqlengths=c(chrA=1000, chrB=2000))

    Np <- 20
    all.anchor1 <- sample(N1, Np, replace=TRUE)
    all.anchor2 <- sample(N2, Np, replace=TRUE)
    GenomicInteractions(all.anchor1, all.anchor2, list(regions1, regions2))
}

expect_as_if <- function(x, y) {
    expect_identical(first(x), first(y))
    expect_identical(second(x), second(y))
    expect_identical(mcols(x), mcols(y))
    expect_identical(metadata(x), metadata(y))
    expect_identical(names(x), names(y))
}
