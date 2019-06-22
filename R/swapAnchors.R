#' Swap anchor regions 
#' 
#' Swap anchor regions within an interaction to ensure that the first anchor is not greater than the second.
#'
#' @param x A \linkS4class{GenomicInteractions} object.
#' @param mode String specifying the type of swapping to perform.
#'
#' @details
#' \code{mode="order"} will ensure that the first anchor region is not ordered above the second anchor in the output object.
#'
#' \code{mode="reverse"} will ensure that the first anchor region is not ordered below the second anchor in the output.
#'
#' \code{mode="all"} will simply swap the first and second anchor regions without any consideration of ordering.
#'
#' @return A GenomicInteractions object with some or all of first and second anchor regions swapped.
#'
#' @author Aaron Lun
#' 
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr2", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(5, 50, 3, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' swapAnchors(test)
#' swapAnchors(test, mode="reverse")
#' swapAnchors(test, mode="all")
#'
#' @export
#' @name swapAnchors
#' @aliases swapAnchors swapAnchors,GenomicInteractions-method
#' @importFrom utils head tail
#' @importFrom S4Vectors first second
setMethod("swapAnchors", "GenomicInteractions", function(x, mode=c("order", "reverse", "all")) {
    mode <- match.arg(mode)
    y <- initialize(x, first=second(x), second=first(x))

    if (mode=="order" || mode=="reverse") {
        combined <- c(x, y)
        x <- head(combined, length(x))
        y <- tail(combined, length(y))

        if (mode=="order") {
            replace <- x > y
        } else {
            replace <- x < y
        }
        x[replace] <- y[replace]

    } else {
        x <- y
    }

    x 
})
