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
#' If \code{query} is the IndexedRelations object, the anchor regions will also be the \code{query} in the \code{\link{findOverlaps}} call.
#' Conversly, if \code{subject} is the IndexedRelations, the anchor regions will be the \code{subject}.
#' This has implications for overlap settings that are not symmetric, e.g., \code{type="within"}.
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
#' @export
#' @rdname findOverlaps
#' @aliases findOverlaps
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors selectHits Hits
#' @importFrom BiocGenerics unique intersect union sort
setMethod("findOverlaps", c("GenomicInteractions", "ANY"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ..., use.region=c("any", "first", "second", "both"))
{
    use.region <- match.arg(use.region)
    if (use.region=="both") {
        .Deprecated(msg="'use.region=\"both\" is deprecated.\nUse 'use.region=\"any\"' instead.")
        use.region <- "any"
    }
    arg.pack <- list(maxgap = maxgap, minoverlap = minoverlap, type = match.arg(type), ...)

    if (use.region=="first" || use.region=="any") {
        out <- do.call(.find_single_overlap_left, c(list(query, 1L, subject), arg.pack))
        hits1 <- Hits(out$query, out$subject, length(query), length(subject))
    } 
    
    if (use.region=="second" || use.region=="any") {
        out <- do.call(.find_single_overlap_left, c(list(query, 2L, subject), arg.pack))
        hits2 <- Hits(out$query, out$subject, length(query), length(subject))
    }

    if (use.region=="any") { 
        hits <- union(hits1, hits2)
    } else if (use.region=="first") {
        hits <- hits1
    } else if (use.region=="second") {
        hits <- hits2
    }

    hits <- sort(hits)
    selectHits(hits, select = match.arg(select))
})

#' @export
#' @rdname findOverlaps
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors selectHits Hits
#' @importFrom BiocGenerics unique intersect union sort
setMethod("findOverlaps", c("ANY", "GenomicInteractions"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ..., use.region=c("any", "first", "second", "both"))
{
    use.region <- match.arg(use.region)
    if (use.region=="both") {
        .Deprecated(msg="'use.region=\"both\" is deprecated.\nUse 'use.region=\"any\"' instead.")
        use.region <- "any"
    }
    arg.pack <- list(maxgap = maxgap, minoverlap = minoverlap, type = match.arg(type), ...)

    if (use.region=="first" || use.region=="any") {
        out <- do.call(.find_single_overlap_right, c(list(query, 1L, subject), arg.pack))
        hits1 <- Hits(out$query, out$subject, length(query), length(subject))
    } 
    
    if (use.region=="second" || use.region=="any") {
        out <- do.call(.find_single_overlap_right, c(list(query, 2L, subject), arg.pack))
        hits2 <- Hits(out$query, out$subject, length(query), length(subject))
    }

    if (use.region=="any") { 
        hits <- union(hits1, hits2)
    } else if (use.region=="first") {
        hits <- hits1
    } else if (use.region=="second") {
        hits <- hits2
    }

    hits <- sort(hits)
    selectHits(hits, select = match.arg(select))
})


