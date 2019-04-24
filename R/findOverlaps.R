#' Find overlaps between GenomicInteractions
#'
#' Find one-dimensional overlaps between GenomicInteractions and other genomic intervals,
#' or two-dimensional overlaps between two GenomicInteractions objects.
#'
#' @param query A \linkS4class{GenomicInteractions} instance or any \linkS4class{Vector}-like object that can be overlapped with a \linkS4class{GenomicRanges} via \code{\link{findOverlaps}}.
#' @param subject Same as \code{query}.
#' @param maxgap Integer scalar specifying the maximum gap between two genomic intervals to consider them as overlapping, see \code{?\link{findOverlaps}} for details.
#' @param minoverlap Integer scalar specifying the minimum overlap between two genomic intervals to consider them as overlapping, see \code{?\link{findOverlaps}} for details.
#' @param type String specifying the type of overlap, see \code{?\link{findOverlaps}} for details.
#' @param select String specifying what kind of output to return, see \code{?\link{findOverlaps}} for details.
#' @param ignore.strand Logical scalar indicating whether strand information should be ignored.
#' If \code{FALSE}, stranded intervals must be on the same strand to be considered as overlapping.
#' @param ... Further arguments to pass to \code{\link{findOverlaps}}.
#' @param use.region String specifying which regions should be used to perform the overlap, see below.
#'
#' @return
#' If \code{type="any"}, a \linkS4class{Hits} object is returned specifying the overlaps between \code{query} and \code{subject}.
#'
#' Otherwise, an integer vector is returned of length equal to \code{length(query)}, containing the selected index of \code{subject} that is overlapped by each entry of \code{query} (or \code{NA}, if there are no overlaps).
#'
#' @section One-dimensional overlaps:
#' If only one of \code{query} or \code{subject} is an \linkS4class{IndexedRelations} object,
#' a one-dimensional overlap is performed. 
#' This involves identifying overlaps between the individual anchor regions without consideration of the pairing of anchor regions in each interaction.
#' One-dimensional overlaps are useful for identifying any kind of overlap between the interactions' anchor regions and the query/subject of interest.
#' 
#' Let's say that \code{query} is the IndexedRelations object, in which case:
#' \itemize{
#' \item If \code{use.region="any"}, an element of \code{subject} is considered to overlap an element of \code{query} if the former overlaps either of the anchor regions of the latter.
#' \item If \code{use.region="first"}, an element of \code{subject} is considered to overlap an element of \code{query} if the former overlaps the first anchor region of the latter.
#' \item If \code{use.region="second"}, an element of \code{subject} is considered to overlap an element of \code{query} if the former overlaps the second anchor region of the latter.
#' }
#' The same principles apply when \code{subject} is the IndexedRelations object.
#' 
#' Overlaps between genomic regions are defined based on the various parameters passed to \code{\link{findOverlaps}},
#' e.g., \code{maxgap}, \code{minoverlap}.
#' If \code{query} is the GenomicInteractions object, its anchor regions will also be the \code{query} in the \code{\link{findOverlaps}} call.
#' Conversely, if \code{subject} is the IndexedRelations, the anchor regions will be the \code{subject}.
#' This has implications for overlap settings that are not symmetric, e.g., \code{type="within"}.
#' 
#' By default, all overlaps with GenomicInteractions objects are performed with \code{ignore.strand=TRUE}.
#' This reflects the fact that interactions generally occur between unstranded genomic loci rather than stranded transcripts.
#'
#' @section Two-dimensional overlaps:
#' If both \code{query} and \code{subject} are \linkS4class{IndexedRelations} objects,
#' a two-dimensional overlap can be performed.
#' An interation in \code{query} only overlaps an interaction in \code{subject} if both of the \code{query}'s anchor regions overlaps both of the \code{subject}'s anchor regions.
#' This is useful for identifying interactions that span the same part of the two-dimensional interaction space.
#'
#' The exact nature of the overlap can be controlled with \code{use.region}:
#' \itemize{
#' \item If \code{use.region="match"}, the first query anchor must overlap the first subject anchor, and
#' while the second query anchor must overlap the second subject anchor.
#' \item If \code{use.region="reverse"}, the first query anchor must overlap the second subject anchor, and
#' while the second query anchor must overlap the first subject anchor.
#' \item If \code{use.region="any"}, the first query anchor can overlap either one subject anchor, 
#' while the second query anchor must overlap the other subject anchor.
#' }
#'
#' Again, overlaps between genomic regions are defined based on the various parameters passed to \code{\link{findOverlaps}},
#' e.g., \code{maxgap}, \code{minoverlap}.
#' The query's anchor regions will be used as the \code{query} in the \code{\link{findOverlaps}} call,
#' which ensures that asymmetric modes like \code{type="within"} are correctly handled.
#'
#' @section More one-dimensional overlaps:
#' Even if both \code{query} and \code{subject} are \linkS4class{IndexedRelations} objects,
#' a one-dimensional overlap can still be performed.
#' This will identify overlaps between the anchor regions of one object with the anchor regions of the other object,
#' without consideration of the pairing between anchors.
#' It gives the same results as calling \code{findOverlaps} on the \code{partnerFeatures} of one of the objects,
#' but is more efficient as it does not explicitly create the GenomicRanges object for a given set of anchor regions.
#' 
#' The exact nature of the overlap can be controlled with \code{use.region}.
#' This has the form of \code{"QMODE-SMODE"} where \code{QMODE} and \code{SMODE} can be any of \code{"any"}, \code{"first"} or \code{"second"}, and determines the anchor regions that are used from each object.
#' For example, \code{"any-any"} will define overlaps between two interactions if any of the anchor regions in the query overlap any of the anchor regions in the subject.
#' Setting \code{"first-second"} will only consider overlaps between the first anchor region in the query and the second anchor region in the subject.
#'
#' @examples
#' ######################
#' ## One-dimensional ###
#' ######################
#'
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' findOverlaps(test, GRanges("chr1:25-100"))
#' findOverlaps(test, GRanges("chr1:25-100"), use.region="first")
#' findOverlaps(test, GRanges("chr1:25-100"), use.region="second")
#' 
#' findOverlaps(GRanges("chr1:25-100"), test)
#' findOverlaps(GRanges("chr1:25-100"), test, use.region="first")
#' findOverlaps(GRanges("chr1:25-100"), test, use.region="second")
#'
#' ######################
#' ## Two-dimensional ###
#' ######################
#' 
#' alt <- GenomicInteractions(
#'    GRanges("chr1:1-20"), GRanges("chr1:100-150")
#' )
#' findOverlaps(test, alt)
#' findOverlaps(test, alt, use.region="match")
#' findOverlaps(test, alt, use.region="reverse")
#'
#' alt2 <- GenomicInteractions(
#'    GRanges("chr2:1-60"), GRanges("chr1:1-30")
#' )
#' suppressWarnings(findOverlaps(alt2, test))
#' suppressWarnings(findOverlaps(alt2, test, use.region="match"))
#' findOverlaps(alt2, test, use.region="reverse")
#' 
#' ###########################
#' ## More one-dimensional ###
#' ###########################
#'
#' findOverlaps(alt, test, use.region="any-any")
#' findOverlaps(alt, test, use.region="first-any")
#' findOverlaps(alt, test, use.region="any-second")
#' findOverlaps(alt, test, use.region="first-first")
#' findOverlaps(alt, test, use.region="second-second")
#' 
#' @author Aaron Lun
#' @export
#' @name findOverlaps
#' @aliases findOverlaps
#' @importMethodsFrom IRanges findOverlaps
setMethod("findOverlaps", c("GenomicInteractions", "ANY"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ignore.strand=TRUE, ..., 
        use.region=c("any", "first", "second", "both"))
{
    .find_single_overlap_wrapper(query, subject, FUN=.find_single_overlap_left,
        use.region=match.arg(use.region), select=match.arg(select), 
        maxgap=maxgap, minoverlap=minoverlap, type=match.arg(type), 
        ignore.strand=ignore.strand, ...)
})

#' @export
#' @rdname findOverlaps
#' @importMethodsFrom IRanges findOverlaps
setMethod("findOverlaps", c("ANY", "GenomicInteractions"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ignore.strand=TRUE, ..., 
        use.region=c("any", "first", "second", "both"))
{
    .find_single_overlap_wrapper(query, subject, FUN=.find_single_overlap_right,
        use.region=match.arg(use.region), select=match.arg(select), 
        maxgap=maxgap, minoverlap=minoverlap, type=match.arg(type),
        ignore.strand=ignore.strand, ...)
})

#' @export
#' @rdname findOverlaps
#' @importMethodsFrom IRanges findOverlaps
setMethod("findOverlaps", c("GenomicInteractions", "GenomicInteractions"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ignore.strand=TRUE, ..., use.region="any")
{
    use.region <- match.arg(use.region, c(.options_2d, .options_1.5d))
    select <- match.arg(select)

    if (use.region %in% .options_2d) {
        if (use.region=="any") {
            do.same <- do.reverse <- TRUE
        } else if (use.region=="match") {
            do.same <- TRUE
            do.reverse <- FALSE
        } else if (use.region=="reverse") {
            do.same <- FALSE 
            do.reverse <- TRUE
        }

        hits <- .find_double_overlap(query, subject, do.same=do.same, do.reverse=do.reverse,
            maxgap=maxgap, minoverlap=minoverlap, type=match.arg(type), 
            ignore.strand=ignore.strand, ..., find.arbitrary=(select=="arbitrary"))

    } else {
        split.arg <- strsplit(use.region, "-")[[1]]
        q.opt <- split.arg[[1]]
        s.opt <- split.arg[[2]]

        if (q.opt=="any") {
            q.left <- q.right <- TRUE
        } else if (q.opt=="first") {
            q.left <- TRUE
            q.right <- FALSE
        } else if (q.opt=="second") {
            q.left <- FALSE 
            q.right <- TRUE
        }
        
        if (s.opt=="any") {
            s.left <- s.right <- TRUE
        } else if (s.opt=="first") {
            s.left <- TRUE
            s.right <- FALSE
        } else if (s.opt=="second") {
            s.left <- FALSE
            s.right <- TRUE
        }

        hits <- .find_single_overlap_IR(query, subject, 
            query.left=q.left, query.right=q.right, subject.left=s.left, subject.right=s.right,
            maxgap=maxgap, minoverlap=minoverlap, type=match.arg(type), 
            ignore.strand=ignore.strand, ...)
    }

    selectHits(hits, select=match.arg(select))
})
