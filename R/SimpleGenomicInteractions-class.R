#' @title The SimpleGenomicInteractions class
#'
#' @description
#' The SimpleGenomicInteractions class is a subclass of the \linkS4class{GenomicInteractions} class,
#' where all indices are guaranteed to point to a single set of regions.
#'
#' This simplifies preservation of backwards compatibility for many methods with the old GenomicInteractions class,
#' especially if they perform operations on the \code{regions}.
#' 
#' @section Constructors:
#' \code{SimpleGenomicInteractions(anchor1, anchor2, regions, ...)} will create a SimpleGenomicInteractions object, given:
#' \itemize{
#' \item A \linkS4class{GenomicRanges} object in each of \code{anchor1} and \code{anchor2}, and missing \code{regions}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, 
#' and a \linkS4class{SimpleGenomicRangesList} of length 1 in \code{regions}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, and a GenomicRanges in \code{regions}.
#' }
#'
#' @section Getters and setters:
#' The only getter that differs from those for \linkS4class{GenomicInteractions} is the \code{regions} method.
#' This has the form of \code{regions(x, as.list=FALSE)} where \code{x} is a SimpleGenomicInteractions method.
#' It returns the \linkS4class{GenomicRanges} of the set of regions directly, 
#' unless \code{as.list=TRUE}, in which case a \linkS4class{GenomicRangesList} is returned 
#' (consistent with the GenomicInteractions method).
#'
#' Similarly, the setter has the form \code{regions(x, as.list=FALSE) <- value}.
#' If \code{as.list=FALSE}, \code{value} should be a GenomicRanges object of the same length as \code{regions(x)}.
#' If \code{as.list=TRUE}, \code{value} should be a GenomicRangesList object of the same lengths as \code{regions(x)}.
#' 
#' @author Aaron Lun
#' @name SimpleGenomicInteractions
#' @docType class
#' @examples
#'
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- SimpleGenomicInteractions(anchor1, anchor2)
#' test
#'
#' regions(test)
#' regions(test) <- resize(regions(test), 1000)
#' test
#'
#' @aliases SimpleGenomicInteractions
#' regions,SimpleGenomicInteractions-method
#' regions<-SimpleGenomicInteractions-method
NULL

#' @export
#' @rdname SimpleGenomicInteractions
#' @importFrom BiocGenerics match unique
#' @importFrom IndexedRelations IndexedRelations
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importFrom S4Vectors mcols<- DataFrame
SimpleGenomicInteractions <- function(anchor1, anchor2, regions, ...) {
    if (missing(regions)) {
        regions <- unique(c(anchor1, anchor2))
        anchor1 <- match(anchor1, regions)
        anchor2 <- match(anchor2, regions)
    } 
    if (is(regions, "GenomicRanges")) {
        regions <- List(regions)
    }
    out <- IndexedRelations(list(anchor1, anchor2), regions, mapping=c(1L, 1L))

    meta <- list(...)
    if (length(meta)) {
        mcols(out) <- do.call(DataFrame, meta)
    }
    new("SimpleGenomicInteractions", out)
}

#' @importFrom IndexedRelations featureSets
setValidity2("SimpleGenomicInteractions", function(object) {
    if (length(featureSets(object))!=1L) {
        return("only one set of regions should be present")
    }
    TRUE
})

#' @export
#' @rdname SimpleGenomicInteractions
setMethod("regions", "SimpleGenomicInteractions", function(x, as.list=FALSE) {
    reg <- callNextMethod()
    if (!as.list) { 
        reg[[1]]
    } else {
        reg
    } 
})

#' @export
#' @rdname SimpleGenomicInteractions
#' @importFrom S4Vectors List
setReplaceMethod("regions", "SimpleGenomicInteractions", function(x, as.list=FALSE, ..., value) {
    if (!as.list) { 
        value <- List(value)
    } 
    callNextMethod()
})

