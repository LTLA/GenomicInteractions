#' Converting to/from pairs
#'
#' Convert a \linkS4class{GenomicInteractions} object to a \linkS4class{Pairs} object, or vice versa.
#'
#' @param x A Pairs object with two \linkS4class{GRanges} and/or \linkS4class{GRangesFactor} entries.
#' 
#' @return An object of the converted type.
#' 
#' @author Aaron Lun
#'
#' @details
#' An automated coercion from a Pairs object to a GenomicInteractions object is not possible,
#' as the former may not contain two \linkS4class{GRangesFactor} objects (in which case no standard conversion exists).
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
#' makeGInteractionsFromGRangesPairs(p)
#'
#' # Also works with two GRanges.
#' alt <- Pairs(unfactor(first(p)), unfactor(second(p)))
#' makeGInteractionsFromGRangesPairs(alt)
#'
#' @name conversion
#' @aliases coerce,GenomicInteractions,Pairs-method
NULL

#' @export
#' @name conversion
#' @importFrom S4Vectors first second mcols<-
#' @importClassesFrom S4Vectors Pairs
makeGInteractionsFromGRangesPairs <- function(x) {
    if (!is(x, "Pairs")) { 
        stop("'x' must be a Pairs object")
    }
    f1 <- .convert_to_factor(first(x), "first")
    f2 <- .convert_to_factor(second(x), "second")

    out <- GenomicInteractions(anchor1=f1, anchor2=f2)
    mcols(out) <- mcols(x)
    names(out) <- names(x)
    out
}

#' @importFrom methods is
#' @importFrom S4Vectors Factor
#' @importClassesFrom GenomicRanges GRangesFactor GRanges
.convert_to_factor <- function(x, msg) {
    if (!is(x, "GRangesFactor")) {
        if (!is(x, "GRanges")) {
            stop(sprintf("%s element must be a GRanges or GRangesFactor", msg))
        }
        x <- Factor(x)
    }
    x
}
