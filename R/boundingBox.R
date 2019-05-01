#' Get the bounding box
#' 
#' Compute a minimum bounding box for groups of pairwise genomic interactions.
#' 
#' @param x A \linkS4class{GenomicInteractions} object.
#' @param f A factor or vector of length equal to that of \code{x}, 
#' indicating the group to which each pairwise interaction belongs.
#' @param reflect A logical scalar indicating whether intra-chromosomal interactions can be flipped around the diagonal to achieve a smaller bounding box.
#' 
#' @return A \linkS4class{GenomicInteractions} object containing the coordinates of each bounding box.
#' 
#' @details
#' For any group of pairwise interactions, the minimum bounding box is the smallest rectangle in the interaction space that contains all interactions in the group.
#' Each side of the box has coordinates spanning the most extreme anchor regions on the corresponding chromosome.
#' This is often useful for summarizing clusters of interactions.
#' 
#' Grouping of interactions is specified using \code{f}, where interactions in \code{x} with the same level of \code{f} are considered to be in the same group.
#' If \code{f} is not specified, all interactions in \code{x} are assumed to be in a single group (named as ``1'').
#' An error will be raised if a group spans multiple chromosomes for either the first or second anchor regions.
#' 
#' The function returns a \linkS4class{GenomicInteractions} object containing the coordinates of the bounding boxes for all groups.
#' Each interaction represents a bounding box for a group, where the anchor regions represent the sides of the box.
#' Entries are named according to the levels of \code{f}, in order to specify which bounding box corresponds to which group.
#' 
#' \code{reflect=TRUE} will ensure that all interactions lie on one side of the diagonal of the interaction space.
#' (Specifically, the first anchor will always start before the second anchor, see \code{\link{swapAnchors}} for details.)
#' This has a number of implications:
#' \itemize{
#' \item For intra-chromosomal groups, this generally yields a smaller bounding box, i.e., with a smaller maximum dimension.
#' \item For inter-chromosomal groups, setting \code{reflect=TRUE} means that the function will work properly if the chromosomes are swapped between the first and second anchors for different interactions.
#' Otherwise, an error will be raised.
#' \item If \code{reflect=TRUE}, the output GenomicInteractions object will have the same universe of regions for both anchors, equivalent to \code{common=TRUE} in the \code{\link{GenomicInteractions}} constructor.
#' Otherwise, it will have two separate GenomicRanges for the first and second anchor regions.
#' }
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' gi <- GenomicInteractions(anchor1, anchor2)
#' 
#' # Making up a sensible grouping, e.g., chromosome pairs.
#' f <- c(1,1,1,2)
#' boundingBox(gi, f)
#' 
#' # Fails for multiple chromosomes
#' try(out <- boundingBox(gi))
#'
#' @export
#' @name boundingBox
#' @aliases boundingBox
#' @importFrom GenomeInfoDb seqnames seqinfo
#' @importFrom BiocGenerics start end
#' @importFrom IndexedRelations partners
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors merge
setMethod("boundingBox", "GenomicInteractions", function(x, f, reflect=TRUE) {
    if (missing(f)) { 
        f <- rep(1L, length(x))
    }
    if (length(f)!=length(x)) { 
        stop("length of 'f' must be equal to number of interactions")
    } 

    # Get the run lengths and values in a slightly convoluted way to handle factors (as rle() does not).
    by.f <- split(seq_along(f), f)
    o <- unlist(by.f)
    f.runs <- lengths(by.f)
    f.values <- names(by.f)

    reg1 <- featureSets(x)[[1]]
    reg2 <- featureSets(x)[[2]]

    out <- bounding_box(f.runs, f.values,
        partners(x)[,1][o], as.character(seqnames(reg1)), start(reg1), end(reg1),
        partners(x)[,2][o], as.character(seqnames(reg2)), start(reg2), end(reg2),
        reflect)
    bound1 <- out[[1]]
    bound2 <- out[[2]]

    if (reflect) {
        si1 <- si2 <- merge(seqinfo(reg1), seqinfo(reg2))
    } else {
        si1 <- seqinfo(reg1)
        si2 <- seqinfo(reg2)
    }

    gr1 <- GRanges(bound1[[1]], IRanges(bound1[[2]], bound1[[3]]), seqinfo=si1)
    gr2 <- GRanges(bound2[[1]], IRanges(bound2[[2]], bound2[[3]]), seqinfo=si2)
    out <- GenomicInteractions(gr1, gr2, common=reflect) 
    names(out) <- f.values
    out
})

# Note on reflect=TRUE:
#
# It is fairly straightforward to show that the bounding box of the _start_
# positions (i.e., the starts of the two anchor regions) minimizes its maximum
# dimension when all starts lie on one side of the diagonal.
#
# 1. Pick one point and fix it to one side of the diagonal. This simply 
#    avoids the need to consider both sides of the diagonal.
# 2. The pairwise minimum bounding box containing the fixed point and 
#    any other point is minimized (in terms of its larger dimension)
#    when both points are on the same side of the diagonal.
# 3. The larger dimension of the bounding box of all points is equal
#    to the maximum of the larger dimensions of the pairwise boxes.
# 4. Keeping all points on one side of the diagonal (same as that 
#    of the fixed point) minimizes the pairwise boxes, and thus 
#    minimizes the bounding box of call points.
#
# Unfortunately, this does not generalize to the intervals defined
# by each pairwise interaction. A simple counterexample:
# 
#   ir1 <- GenomicInteractions(GRanges("A:1-1"), GRanges("A:10-10"))
#   ir2 <- GenomicInteractions(GRanges("A:9-1000"), GRanges("A:10-10"))
#   boundingBox(c(ir1, ir2))
#   boundingBox(c(ir1, rearrangePartners(ir2, 2:1)), reflect=FALSE)
# 
# You can see the bounding box for the second instance has a smaller
# maximum dimension (991 vs 1000).

