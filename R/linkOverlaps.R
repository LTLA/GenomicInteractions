#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
.linkOverlaps <- function(query, subject1, subject2, ...) {
    olap1 <- findOverlaps(query, subject1, ..., use.region="first")
    olap2 <- findOverlaps(query, subject2, ..., use.region="second")
    tab1 <- data.frame(query=queryHits(olap1), subject1=subjectHits(olap1))
    tab2 <- data.frame(query=queryHits(olap2), subject2=subjectHits(olap2))
    output <- merge(tab1, tab2, by="query")
    as(output, "DataFrame")
}

#' Link overlapping regions
#' 
#' Identify interactions that link two sets of regions by having anchor regions overlapping one entry in each set.
#' 
#' @param query A \linkS4class{GenomicInteractions} object.
#' @param subject1 A \linkS4class{Vector} object supporting one-dimensional \code{\link{findOverlaps}} with \code{query}.
#' @param subject2 A \linkS4class{Vector} object supporting one-dimensional \code{\link{findOverlaps}} with \code{query}.
#' This can be missing.
#' @param ... Additional arguments to be passed to \code{\link{findOverlaps}}.
#' @param use.region String specifying which anchor regions of \code{query} should be matched to which \code{subject}s.
#' Ignored when \code{subject2} is missing.
#' 
#' @details
#' This function identifies all interactions in \code{query} where one anchor overlaps an entry in \code{subject1} and the other anchor overlaps an entry in \code{subject2}.
#' It is designed to be used to identify regions that are linked by interactions in \code{query}.
#' For example, one might specify genes as \code{subject1} and enhancers as \code{subject2}, to identify all gene-enhancer contacts present in \code{query}.
#' This is useful when the exact pairings between \code{subject1} and \code{subject2} are undefined.
#' 
#' The function returns a DataFrame specifying the index of the interaction in \code{query}; 
#' the index of the overlapped region in \code{subject1};
#' and the index of the overlapped region in \code{subject2}.
#' If multiple regions in \code{subject1} and/or \code{subject2} overlap the anchor regions of a particular interaction,
#' all combinations of two overlapping regions (one from each \code{subject*} set) are reported for that interaction.
#' 
#' If \code{subject2} is not specified, links within \code{subject1} are identified instead, i.e., \code{subject2} is set to \code{subject1}.
#' In such cases, the returned DataFrame is such that the first subject index is always less than the second subject index, to avoid redundant permutations.
#' Any trivial interaction of a subject with itself is discarded from the output.
#' 
#' If \code{use.region="any"}, there is no constraint on the overlaps between anchor regions and entries in \code{subject1} and \code{subject2}.
#' If \code{use.region="match"}, linking overlaps are only considered between the first anchor region and \code{subject1},
#' and the second anchor region and \code{subject2}.
#' (Obviously, users can simply swap the arguments in \code{subject1} and \code{subject2} to achieve the reverse effect.)
#' All settings of \code{use.region} are ignored when \code{subject2} is missing, 
#' as both anchors must overlap \code{subject1}.
#' 
#' @return
#' A DataFrame of integer indices indicating which elements of \code{query} link which elements of \code{subject1} and \code{subject2}.
#' 
#' @seealso
#' \code{\link{findOverlaps}} for GenomicInteractions instances.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' all.genes <- GRanges("chr1", IRanges(0:9*20, 1:10*20))
#' all.enhancers <- GRanges("chr2", IRanges(0:9*20, 1:10*20))
#' 
#' out <- linkOverlaps(test, all.genes, all.enhancers)
#' out
#' 
#' out <- linkOverlaps(test, all.genes)
#' out
#' 
#' @export
#' @aliases linkOverlaps
#' @name linkOverlaps
#' @importFrom BiocGenerics unique rbind
setMethod("linkOverlaps", c("GenomicInteractions", "Vector", "Vector"), 
    function(query, subject1, subject2, ..., use.region=c("any", "match"))
{
    out <- .linkOverlaps(query, subject1, subject2, ...)

    use.region <- match.arg(use.region)
    if (use.region=="any") {
        flipped <- .linkOverlaps(query, subject2, subject1, ...)

        tmp <- flipped$subject1
        flipped$subject1 <- flipped$subject2
        flipped$subject2 <- tmp
        out <- rbind(out, flipped)
        out <- unique(out)
    }

    out
})

#' @export
#' @rdname linkOverlaps
setMethod("linkOverlaps", c("GenomicInteractions", "Vector", "missing"), 
    function(query, subject1, subject2, ..., use.region=NULL)
{
    out <- .linkOverlaps(query, subject1, subject1, ...)
    out[out$subject1 < out$subject2,]
})
