#' Find distal anchor regions
#'
#' Given a reference region, find the distal anchor regions that interact with that region.
#'
#' @param x A \linkS4class{GenomicInteractions} object.
#' @param ref A \linkS4class{GenomicRanges} object of length 1.
#' @param local A logical scalar specifying whether interactions should be retained if both anchors overlap \code{ref}.
#' @param ... Further arguments to pass to \code{\link{findOverlaps}}.
#'
#' @details
#' This function will identify all interactions in \code{x} that have one-dimensional overlaps with \code{ref},
#' and it will return the non-overlapping anchor region for each such interaction.
#' Any per-interaction metadata in \code{x} is carried over to the output object,
#' along with per-feature metadata in \code{regions(x)}.
#' 
#' The returned genomic intervals represent the distal interactors to \code{ref} in \code{x}.
#' One can view this as taking a \dQuote{cross-section} of the interaction space at the specified \code{ref},
#' thus projecting all the interactions involving \code{ref} onto the linear genome.
#' This is analogous to what is done in 4C experiments.
#' 
#' @return A \linkS4class{GenomicRanges} object containing all distal anchor regions.
#' 
#' @author Aaron Lun
#' @examples
#' anchor1 <- GRanges("chr1", IRanges(sample(1000, 100), width=5))
#' anchor2 <- GRanges("chr1", IRanges(sample(1000, 100), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' mcols(test)$original <- seq_along(test)
#' test
#' 
#' findDistalAnchors(test, GRanges("chr1:1-100"))
#' findDistalAnchors(test, GRanges("chr1:1-100"), local=FALSE)
#' findDistalAnchors(test, GRanges("chr1:1-100"), type="within")
#' 
#' @export
#' @name findDistalAnchors
#' @aliases findDistalAnchors findDistalAnchors,GenomicInteractions-method
#' @importFrom IndexedRelations partners featureSets
#' @importFrom IRanges findOverlaps overlapsAny
#' @importFrom S4Vectors mcols mcols<-
setMethod("findDistalAnchors", "GenomicInteractions", function(x, ref, local=TRUE, ...) {
    keep.first <- overlapsAny(x, ref, use.region="second", ...) # yes, the first/second flip is intended.
    keep.second <- overlapsAny(x, ref, use.region="first", ...)

    if (!local) {
        tmp <- keep.first
        keep.first <- keep.first & !keep.second
        keep.second <- keep.second & !tmp
    }

    all.first <- partners(x)[keep.first,1]
    all.first <- featureSets(x)[[1]][all.first]
    all.second <- partners(x)[keep.second,2]
    all.second <- featureSets(x)[[2]][all.second]

    if (!identical(colnames(mcols(all.second)), colnames(mcols(all.first)))) {
        warning("removing mismatching 'mcols' between feature sets")
        mcols(all.second) <- mcols(all.first) <- NULL
    }

    collected <- c(all.first, all.second)
    keep <- c(which(keep.first), which(keep.second))

    o <- order(keep)
    collected <- collected[o]
    keep <- keep[o]

    if (!is.null(mcols(x))) {
        mcols(collected) <- cbind(mcols(collected), mcols(x)[keep,,drop=FALSE])
    }
    collected
})
