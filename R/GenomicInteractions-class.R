#' @title The GenomicInteractions class
#'
#' @description
#' The GenomicInteractions class is dedicated to representing interactions between pairs of genomic intervals.
#' Each genomic interval is referred to as an \dQuote{anchor}, and each interaction is defined between a pair of anchors.
#' We will refer to the two anchors in a pair as the \dQuote{first} and \dQuote{second} anchors (arbitrarily defined).
#'
#' All first anchors in a GenomicInteractions instance are represented by one \linkS4class{GRangesFactor} object.
#' (The same applies for all second anchors.)
#' This improves efficiency by avoiding storage and processing of redundant copies of the same genomic coordinates and metadata;
#' especially for interaction data, where one interval may be involved in multiple interactions.
#' 
#' It is worth clarifying the terminology at this point.
#' The \dQuote{anchors} refer to the genomic intervals themselves;
#' the \dQuote{anchor type} refers to the anchor ordering, i.e., first or second;
#' the \dQuote{anchor levels} are the universe of possible anchors across all interactions for a given anchor type;
#' and the \dQuote{anchor indices} are integer indices to the anchor levels,
#' such that subsetting the anchor levels by the anchor indices yields the anchors.
#'
#' @section Constructors:
#' \code{GenomicInteractions(anchor1, anchor2, regions, ..., metadata=list())}
#' will create a GenomicInteractions object, given:
#' \enumerate{
#' \item A \linkS4class{GRangesFactor} object in each of \code{anchor1} and \code{anchor2}, and missing \code{regions}.
#' \item A \linkS4class{GRanges} object in each of \code{anchor1} and \code{anchor2}, and missing \code{regions}.
#' \item A \linkS4class{GRanges} object in each of \code{anchor1} and \code{anchor2}, 
#' and a GRanges in \code{regions} that contains a superset of intervals in \code{anchor1} and \code{anchor2}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, and a GRanges in \code{regions}.
#' Anchor regions are defined by using integers to index \code{regions}.
#' }
#'
#' Any arguments in \code{...} should be vector-like and are added to the element-wise metadata of the output object.
#' For options 1-3, any metadata in \code{anchor1} or \code{anchor2} are moved to the element-wise metadata,
#' with the anchor of origin prepended to the name.
#'
#' Option 2 has an additional \code{common=TRUE} argument,
#' where the anchor levels for both the first and second anchors are defined as the union of the regions used in both anchors.
#' This is provided for backwards compatibility with the old GenomicInteractions behaviour,
#' and provides identical \code{regions} regardless of \code{type} (see below).
#'
#' The \code{metadata} argument is expected to be a list-like object,
#' to be used as the \code{\link{metadata}} slot of the output object.
#'
#' @section Getters:
#' In the following code snippets, \code{x} is a GenomicInteractions object:
#' \describe{
#' \item{\code{anchors(x, type="both", id=FALSE)}:}{Returns the anchor regions.
#' If \code{type="both"}, it returns a \linkS4class{Pairs} of length 2,
#' where the first and second entries contain a GRanges of the first and second anchors, respectively.
#'
#' A GRanges of only the first anchor regions can be returned directly with \code{type=1} or \code{type="first"}.
#' The same applies for the second anchors with \code{type=2} or \code{type="second"}.
#' 
#' \code{id=TRUE} will return the integer indices pointing to the anchor regions in \code{regions(x)}.
#' If \code{type="both"}, this will be a \linkS4class{DataFrame} of integer vectors, with one column per partner.
#' Otherwise, the appropriate integer vector is directly returned for the partner specified by \code{type}.
#' }
#' \item{\code{regions(x, type=NA)}:}{Returns the universe of all regions corresponding to anchor region \code{type}.
#' Setting \code{type=1} or \code{"first"} returns a GRanges of all first anchors,
#' while setting \code{type=2} or \code{"second"} does the same for all second anchors.
#' Setting \code{type="both"} returns a List of length 2 containing the regions used in the first and second anchors.
#' 
#' \code{type=NA} behaves like \code{type=1} with a deprecation warning.
#' This is intended to encourage users to actually specify a value of \code{type} in calls to \code{regions(x)},
#' as this option was not available in previous versions.
#'
#' Note that \code{type=1} and \code{type=2} may not yield the same result!
#' Anchor indices are only directly comparable if the universe of regions is the same for both anchors.
#' }
#' \item{\code{first(x)}:}{A synonym for \code{anchors(x, 1)}.}
#' \item{\code{second(x)}:}{A synonym for \code{anchors(x, 2)}.}
#' }
#' All getter methods applicable to \linkS4class{Vector} objects can also be used, 
#' e.g., \code{\link{mcols}}, \code{\link{metadata}}. 
#' 
#' @section Setters:
#' In the following code snippets, \code{x} is a GenomicInteractions object:
#' \describe{
#' \item{\code{anchors(x, type="both", id=FALSE) <- value}:}{Replaces the anchor regions specified by \code{type} with \code{value}.
#' If \code{type="both"}, \code{value} should be a Pairs of two GRanges,
#' where the first and second entries are to be used as the replacement first and second anchor regions, respectively.
#' If \code{type=1} or \code{"first"}, \code{value} should be a GRanges to replace the first anchors;
#' same for \code{type=2} or \code{"second"} for the second anchors.
#'
#' \code{id=TRUE} will replace the integer indices pointing to the anchor regions in \code{regions(x)}.
#' If \code{type="both"}, \code{value} should be a \linkS4class{DataFrame} of integer vectors, with one column per partner.
#' Otherwise, it should be an integer vector to replace the indices for the partner specified by \code{type}.
#' }
#' \item{\code{regions(x, type=NA) <- value}:}{Replaces the anchor levels in \code{x}.
#' If \code{type=1} or \code{"first"}, \code{value} should be a GRanges to replace the first anchor levels;
#' same for \code{type=2} or \code{"second"} for the second anchor levels.
#' If \code{type="both"}, \code{value} should be a List of length 2 to replace both levels simultaneously.
#' 
#' \code{type=NA} is provided for backwards compatibility, 
#' and will replace the region sets for \emph{both} anchors with the GRanges \code{value}.
#' This will trigger a deprecation warning - 
#' users wanting this behaviour should use \code{type="both"} with the replacement value set to \code{List(value, value)}.
#' }
#' \item{\code{first(x) <- value}:}{A synonym for \code{anchors(x, 1) <- value}.}
#' \item{\code{second(x) <- value}:}{A synonym for \code{anchors(x, 2) <- value}.}
#' }
#' Again, all setter methods applicable to Vector objects can also be used here, 
#' e.g., \code{\link{mcols<-}}, \code{\link{metadata<-}}.
#'
#' @author Aaron Lun
#' @seealso
#' \linkS4class{IndexedRelations}, from which this class is derived.
#'
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor1 <- sample(anchor1, 20, replace=TRUE)
#' 
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' anchor2 <- sample(anchor2, 20, replace=TRUE)
#' 
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' # Alternative methods of construction:
#' combined <- sort(unique(c(anchor1, anchor2)))
#' m1 <- match(anchor1, combined)
#' m2 <- match(anchor1, combined)
#' test2 <- GenomicInteractions(m1, m2, combined)
#' test2
#'
#' # Getting
#' anchors(test, 1)
#' anchors(test, 2) 
#' regions(test, 1)
#' regions(test, 2)
#' 
#' # Setting
#' anchors(test, 1) <- rev(anchors(test, 1))
#' test
#'
#' first.regions <- regions(test, type=1)
#' regions(test, type=1) <- resize(first.regions, 1000)
#' test
#'
#' @rdname GenomicInteractions
#' @name GenomicInteractions
#' @aliases GenomicInteractions
#' GenomicInteractions,GRangesFactor,GRangesFactor,missing-method
#' GenomicInteractions,GRanges,GRanges,missing-method
#' GenomicInteractions,GRanges,GRanges,GRanges-method
#' GenomicInteractions,integer,integer,GRanges-method
#' GenomicInteractions-class
#' anchors anchors,GenomicInteractions-method anchors<- anchors<-,GenomicInteractions-method
#' regions regions,GenomicInteractions-method regions<- regions<-,GenomicInteractions-method
#' first,GenomicInteractions-method second,GenomicInteractions-method
#' first<-,GenomicInteractions-method second<-,GenomicInteractions-method
NULL

#' @export
#' @importFrom S4Vectors DataFrame mcols mcols<- metadata<-
#' @importFrom BiocGenerics match
setMethod("GenomicInteractions", c("GRangesFactor", "GRangesFactor", "missing"), 
    function(anchor1, anchor2, regions, ..., metadata=list()) 
{
    out <- new("GenomicInteractions", first=anchor1, second=anchor2)

    mcol1 <- mcols(anchor1)
    colnames(mcol1) <- sprintf("anchor1.%s", colnames(mcol1))
    mcol2 <- mcols(anchor2)
    colnames(mcol2) <- sprintf("anchor2.%s", colnames(mcol2))

    meta <- c(list(...), as.list(mcol1), as.list(mcol2))
    if (length(meta)) {
        mcols(out) <- do.call(DataFrame, meta)
    }

    metadata(out) <- as.list(metadata)
    out 
})

#' @export
#' @importFrom GenomicRanges GRangesFactor
setMethod("GenomicInteractions", c("integer", "integer", "GRanges"), 
    function(anchor1, anchor2, regions, ..., metadata=list()) 
{
    anchor1 <- GRangesFactor(index=anchor1, levels=regions)
    anchor2 <- GRangesFactor(index=anchor2, levels=regions)
    GenomicInteractions(anchor1, anchor2, ..., metadata=metadata)
})

#' @export
#' @importFrom GenomicRanges GRangesFactor
setMethod("GenomicInteractions", c("GRanges", "GRanges", "missing"),
    function(anchor1, anchor2, regions, ..., metadata=list(), common=TRUE) 
{
    if (common) {
        regions <- unique(sort(c(anchor1, anchor2)))
        anchor1 <- GRangesFactor(anchor1, levels=regions)
        anchor2 <- GRangesFactor(anchor2, levels=regions)
    } else {
        anchor1 <- GRangesFactor(anchor1)
        anchor2 <- GRangesFactor(anchor2)
    }
    GenomicInteractions(anchor1, anchor2)
})

#' @export
#' @importFrom GenomicRanges GRangesFactor
setMethod("GenomicInteractions", c("GRanges", "GRanges", "GRanges"),
    function(anchor1, anchor2, regions, ..., metadata=list()) 
{
    anchor1 <- GRangesFactor(anchor1, regions)
    anchor2 <- GRangesFactor(anchor2, regions)
    GenomicInteractions(anchor1, anchor2)
})

#' @importFrom S4Vectors parallelSlotNames
setMethod("parallelSlotNames", "GenomicInteractions", function(x) {
    c("first", "second", callNextMethod())    
})

###############################################################################
###############################################################################
# Basic getters and setters.

.convert_type <- function(type) {
    if (is.character(type)) {
        return(match.arg(type, c("first", "second", "both")))
    } else if (is.numeric(type)) {
        if (type==1) {
            return("first")
        } else if (type==2) {
            return("second")
        }
    }
    stop("unknown 'type'")
}

#' @export
#' @importFrom S4Vectors DataFrame Pairs first second 
setMethod("anchors", "GenomicInteractions", function(x, type="both", id=FALSE) {
    type <- .convert_type(type)
    if (type=="both") {
        a1 <- first(x)
        a2 <- second(x)
        if (id) {
            return(DataFrame(first=as.integer(a1), second=as.integer(a2)))
        } else {
            return(Pairs(first=a1, second=a2))
        }
    } else if (type=="first") {
        a1 <- first(x)
        if (id) {
            return(as.integer(a1))
        } else {
            return(a1)
        }
    } else {
        a2 <- second(x)
        if (id) {
            return(as.integer(a2))
        } else {
            return(a2)
        }
    }
})

#' @export
#' @importFrom GenomicRanges GRangesFactor
#' @importFrom S4Vectors first second first<- second<- levels
setMethod("anchors<-", "GenomicInteractions", function(x, type="both", id=FALSE, ..., value) {
    type <- .convert_type(type)
    if (type=="both") {
        if (id) {
            first(x) <- GRangesFactor(levels=levels(first(x)), value[[1]])
            second(x) <- GRangesFactor(levels=levels(second(x)), value[[2]])
        } else {
            first(x) <- first(value)
            second(x) <- second(value)
        }
    } else if (type=="first") {
        if (id) {
            first(x) <- GRangesFactor(levels=levels(first(x)), value)
        } else {
            first(x) <- value
        }
    } else {
        if (id) {
            second(x) <- GRangesFactor(levels=levels(second(x)), value)
        } else {
            second(x) <- value
        }
    }
    x
})

#' @export
#' @importFrom S4Vectors levels first second SimpleList
setMethod("regions", "GenomicInteractions", function(x, type=NA) {
    if (is.na(type)) {
        type <- 1L
        .Deprecated(msg="'regions(..., type=NA)' is deprecated.\nSee '?regions' for alternatives.")
    }
    type <- .convert_type(type)
    if (type=="first") {
        levels(first(x))
    } else if (type=="second") {
        levels(second(x))
    } else {
        SimpleList(first=levels(first(x)), second=levels(second(x)))
    }
})

#' @export
#' @importFrom S4Vectors first<- second<-
setMethod("regions<-", "GenomicInteractions", function(x, type=NA, ..., value) {
    if (is.na(type)) {
        type <- 1L
        .Deprecated(msg="'regions(..., type=NA)' is deprecated.\nSee '?regions' for alternatives.")
    }
    type <- .convert_type(type)
    if (type=="first") {
        levels(first(x)) <- value
    } else if (type=="second") {
        levels(second(x)) <- value
    } else {
        levels(first(x)) <- value[[1]]
        levels(second(x)) <- value[[2]]
    }
    x
})

#' @export
#' @importFrom S4Vectors first
setMethod("first", "GenomicInteractions", function(x) x@first)

#' @export
#' @importFrom S4Vectors first<-
#' @importClassesFrom GenomicRanges GRangesFactor
setReplaceMethod("first", "GenomicInteractions", function(x, ..., value) {
    x@first <- as(value, "GRangesFactor")
    x
})

#' @export
#' @importFrom S4Vectors second
setMethod("second", "GenomicInteractions", function(x) x@second)

#' @export
#' @importFrom S4Vectors second<-
#' @importClassesFrom GenomicRanges GRangesFactor
setReplaceMethod("second", "GenomicInteractions", function(x, ..., value) {
    x@second <- as(value, "GRangesFactor")
    x
})
