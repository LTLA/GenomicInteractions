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
#' If \code{select="all"}, a \linkS4class{Hits} object is returned specifying the overlaps between \code{query} and \code{subject}.
#'
#' Otherwise, an integer vector is returned of length equal to \code{length(query)}, containing the selected index of \code{subject} that is overlapped by each entry of \code{query} (or \code{NA}, if there are no overlaps).
#'
#' @section One-dimensional overlaps:
#' If only one of \code{query} or \code{subject} is an \linkS4class{GenomicInteractions} object,
#' a one-dimensional overlap is performed. 
#' This involves identifying overlaps between the individual anchor regions without consideration of the pairing of anchor regions in each interaction.
#' One-dimensional overlaps are useful for identifying any kind of overlap between the interactions' anchor regions and the query/subject of interest.
#' 
#' Let's say that \code{query} is the GenomoicInteractions object, in which case:
#' \itemize{
#' \item If \code{use.region="any"}, an element of \code{subject} is considered to overlap an element of \code{query} if the former overlaps either of the anchor regions of the latter.
#' \item If \code{use.region="first"}, an element of \code{subject} is considered to overlap an element of \code{query} if the former overlaps the first anchor region of the latter.
#' \item If \code{use.region="second"}, an element of \code{subject} is considered to overlap an element of \code{query} if the former overlaps the second anchor region of the latter.
#' }
#' The same principles apply when \code{subject} is the GenomicInteractions object.
#' 
#' Overlaps between genomic regions are defined based on the various parameters passed to \code{\link{findOverlaps}},
#' e.g., \code{maxgap}, \code{minoverlap}.
#' If \code{query} is the GenomicInteractions object, its anchor regions will also be the \code{query} in the \code{\link{findOverlaps}} call.
#' Conversely, if \code{subject} is the GenomicInteractions, the anchor regions will be the \code{subject}.
#' This has implications for overlap settings that are not symmetric, e.g., \code{type="within"}.
#' 
#' By default, all overlaps with GenomicInteractions objects are performed with \code{ignore.strand=TRUE}.
#' This reflects the fact that interactions generally occur between unstranded genomic loci rather than stranded transcripts.
#'
#' @section Two-dimensional overlaps:
#' If both \code{query} and \code{subject} are \linkS4class{GenomicInteractions} objects,
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
#' @author Aaron Lun
#' @name findOverlaps
#' @aliases findOverlaps findOverlaps,GenomicInteractions,GenomicRanges-method
#' findOverlaps,GenomicRanges,Genomicinteractions-method
#' findOverlaps,GenomicInteractions,GenomicInteractions-method
NULL

#' @export
#' @rdname findOverlaps
#' @importClassesFrom S4Vectors SortedByQueryHits 
#' @importFrom S4Vectors first second selectHits
#' @importFrom BiocGenerics union intersect
#' @importMethodsFrom IRanges findOverlaps
setMethod("findOverlaps", c("GenomicInteractions", "GenomicRanges"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ignore.strand=TRUE, ..., 
        use.region=c("any", "first", "second", "both"))
{
    use.region <- match.arg(use.region)
    select <- match.arg(select)
    type <- match.arg(type)
    FUN <- function(Query, Select="all") {
        findOverlaps(Query, subject, select=Select, maxgap=maxgap, minoverlap=minoverlap,
            type=type, ignore.strand=ignore.strand, ...)
    }

    if (use.region=="first") {
        hits <- FUN(first(query), select)
    } else if (use.region=="second") {
        hits <- FUN(second(query), select)
    } else {
        hits1 <- FUN(first(query))
        hits2 <- FUN(second(query))
        if (use.region=="any") {
            hits <- union(hits1, hits2)
            hits <- as(hits, "SortedByQueryHits")
        } else {
            hits <- intersect(hits1, hits2)
        }
        hits <- selectHits(hits, select=select)
    }
    hits
})

#' @export
#' @rdname findOverlaps
#' @importClassesFrom S4Vectors SortedByQueryHits 
#' @importFrom S4Vectors first second selectHits
#' @importFrom BiocGenerics union intersect
#' @importMethodsFrom IRanges findOverlaps
setMethod("findOverlaps", c("GenomicRanges", "GenomicInteractions"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ignore.strand=TRUE, ..., 
        use.region=c("any", "first", "second", "both"))
{
    use.region <- match.arg(use.region)
    select <- match.arg(select)
    type <- match.arg(type)
    FUN <- function(Subject, Select="all") {
        findOverlaps(query, Subject, select=Select, maxgap=maxgap, minoverlap=minoverlap,
            type=type, ignore.strand=ignore.strand, ...)
    }

    if (use.region=="first") {
        hits <- FUN(first(subject), select)
    } else if (use.region=="second") {
        hits <- FUN(second(subject), select)
    } else {
        hits1 <- FUN(first(subject))
        hits2 <- FUN(second(subject))
        if (use.region=="any") {
            hits <- union(hits1, hits2)
            hits <- as(hits, "SortedByQueryHits")
        } else {
            hits <- intersect(hits1, hits2)
        }
        hits <- selectHits(hits, select=select)
    }
    hits
})

#' @export
#' @rdname findOverlaps
#' @importClassesFrom S4Vectors SortedByQueryHits 
#' @importFrom S4Vectors first second selectHits
#' @importFrom BiocGenerics union intersect
#' @importMethodsFrom IRanges findOverlaps
setMethod("findOverlaps", c("GenomicInteractions", "GenomicInteractions"), 
    function (query, subject, maxgap = -1L, minoverlap = 0L, type = c("any",
        "start", "end", "within", "equal"), select = c("all", "first",
        "last", "arbitrary"), ignore.strand=TRUE, ..., 
        use.region=c("any", "match", "reverse"))
{
    use.region <- match.arg(use.region)
    type <- match.arg(type)
    FUN <- function(Query, Subject) {
        findOverlaps(Query, Subject, maxgap=maxgap, minoverlap=minoverlap,
            type=type, ignore.strand=ignore.strand, ...)
    }

    if (use.region=="any" || use.region=="match") {
        hits.match.1 <- FUN(first(query), first(subject))
        hits.match.2 <- FUN(second(query), second(subject))
        hits.match <- intersect(hits.match.1, hits.match.2)
    }
    if (use.region=="any" || use.region=="reverse") {
        hits.rev.1 <- FUN(first(query), second(subject))
        hits.rev.2 <- FUN(second(query), first(subject))
        hits.rev <- intersect(hits.rev.1, hits.rev.2)
    }

    if (use.region=="any") {
        hits <- union(hits.match, hits.rev)
        hits <- as(hits, "SortedByQueryHits")
    } else if (use.region=="match") {
        hits <- hits.match
    } else {
        hits <- hits.rev
    }

    selectHits(hits, select=match.arg(select))
})
