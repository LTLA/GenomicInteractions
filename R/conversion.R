#' Converting to/from pairs
#'
#' Convert a \linkS4class{GenomicInteractions} object to a \linkS4class{Pairs} object, or vice versa.
#'
#' @param x A Pairs object with two \linkS4class{GenomicRanges} entries.
#' 
#' @return An object of the converted type.
#' 
#' @author Aaron Lun
#'
#' @details
#' An automated coercion from a Pairs object to a GenomicInteractions object is not possible,
#' as the former may not contain two GenomicRanges objects (in which case no standard conversion exists).
#' Hence the need for the more wordy function to be explicit about the expected inputs.
#'
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#'
#' p <- as(test, "Pairs")
#' p
#'
#' makeGInteractionsFromGRangesPairs(p)
#'
#' @name conversion
#' @aliases coerce,GenomicInteractions,Pairs-method
NULL

#' @export
#' @rdname conversion
#' @importClassesFrom S4Vectors Pairs
#' @importFrom S4Vectors Pairs first second mcols
setAs("GenomicInteractions", "Pairs", function(from) {
     out <- Pairs(first(from), second(from), names=names(from))
     mcols(out) <- mcols(from)
     out
})

#' @export
#' @rdname conversion
makeGInteractionsFromGRangesPairs <- function(x) {
    if (!is(x, "Pairs")) { 
        stop("'x' must be a Pairs object")
    }
    if (!is(first(x), "GRanges") || !is(second(x), "GRanges")) {
        stop("both paired elements must be GRanges")
    }
    out <- GenomicInteractions(anchor1=first(x), anchor2=second(x))
    mcols(out) <- mcols(x)
    names(out) <- names(x)
    out
}

