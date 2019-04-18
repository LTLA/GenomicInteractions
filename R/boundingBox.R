#' Get the bounding box
#' 
#' Compute a minimum bounding box for groups of pairwise genomic interactions.
#' 
#' @param A \linkS4class{GenomicInteractions} object.
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
#' For intra-chromosomal groups, \code{reflect=TRUE} will ensure that all interactions lie on one side of the diagonal of the intra-chromosomal interaction space.
#' (Specifically, the first anchor will always start before the second anchor.)
#' This yields a minimum bounding box with the smallest possible larger dimension,
#' which will only increase in size if interactions are placed on the other side of the diagonal.
#' 
#' % Bit of a pain to prove, but basically, if you flip a point to be above the diagonal, the Chebyshev distance to a point below the diagonal will always increase.
#' % This means that you must increase the size of one of your sides of your bounding box.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' example(GenomicInteractions, echo=FALSE)
#' 
#' # Making up a sensible grouping, e.g., chromosome pairs.
#' gi <- sort(gi)
#' all.chrs <- as.character(seqnames(regions(gi)))
#' f <- paste0(all.chrs[anchors(gi, type="first", id=TRUE)], ".",
#'             all.chrs[anchors(gi, type="second", id=TRUE)])
#' 
#' boundingBox(gi, f)
#' 
#' # Fails for multiple chromosomes
#' try(out <- boundingBox(gi))
#' in.A <- f=="chrA.chrA"
#' out <- boundingBox(gi[in.A])
#'
#' @export
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom IndexedRelations partners mapping
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
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

    reg1 <- regions(x)[[mapping(x)[1]]]
    reg2 <- regions(x)[[mapping(x)[2]]]

    out <- bounding_box(f.runs, f.values, 
        partners(x)[,1][o], as.character(seqnames(reg1)), start(reg1), end(reg1),
        partners(x)[,2][o], as.character(seqnames(reg2)), start(reg2), end(reg2))
    bound1 <- out[[1]]
    bound2 <- out[[2]]

    gr1 <- GRanges(bound1[[1]], IRanges(bound1[[2]], bound1[[3]]), seqinfo=seqinfo(x)) 
    gr2 <- GRanges(bound2[[1]], IRanges(bound2[[2]], bound2[[3]]), seqinfo=seqinfo(x))
    out <- GenomicInteractions(list(gr1, gr2))
    names(out) <- f.values
    out
})
