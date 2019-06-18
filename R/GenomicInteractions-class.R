#' @title The GenomicInteractions class
#'
#' @description
#' The GenomicInteractions class is a subclass of the \linkS4class{IndexedRelations} class,
#' dedicated to representing interactions between pairs of genomic intervals.
#'
#' @section Overview:
#' The main specializations from the base \linkS4class{IndexedRelations} class are:
#' \itemize{
#' \item The \code{featureSets} slot must be a \linkS4class{SimpleGenomicRangesList} object.
#' \item There must be exactly two partners for each relationship.
#' }
#' 
#' In short, each individual feature set is now a \linkS4class{GenomicRanges}.
#' Partners are referred to as \dQuote{anchor regions}, so called as they anchor the interactions between genomic loci.
#'
#' @section Constructors:
#' \code{GenomicInteractions(anchor1, anchor2, regions, ..., metadata=list(), common=TRUE)}
#' will create a GenomicInteractions object, given:
#' \enumerate{
#' \item A \linkS4class{GenomicRanges} object in each of \code{anchor1} and \code{anchor2}, and missing \code{regions}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, 
#' and a \linkS4class{SimpleGenomicRangesList} of length 2 in \code{regions}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, and a GenomicRanges in \code{regions}.
#' }
#'
#' Any arguments in \code{...} should be vector-like and are added to the element-wise metadata of the output object.
#' For option 1, any metadata in \code{anchor1} or \code{anchor2} are moved to the element-wise metadata,
#' with the anchor of origin prepended to the name.
#'
#' If \code{common=TRUE} for option 1, the region set for each anchor will be the union of the regions used in both anchors.
#' This is provided for backwards compatibility with the old GenomicInteractions behaviour,
#' and provides identical \code{regions} regardless of \code{type} (see below).
#'
#' The \code{metadata} argument is expected to be a list-like object,
#' to be used as the \code{\link{metadata}} slot of the output object.
#'
#' @section Getters:
#' In the following code snippets, \code{x} is a GenomicInteractions object:
#' \describe{
#' \item{\code{anchors(x, type=NULL, id=FALSE)}:}{Returns the anchor regions.
#' If \code{type=NULL}, it returns a \linkS4class{Pairs} of length 2,
#' where the first and second entries contain a GenomicRanges of the first and second anchor regions, respectively.
#' A GenomicRanges of only the first or second anchor regions can be returned directly with \code{type=1} or \code{type=2},
#' respectively.
#' 
#' It is also possible to use strings as \code{type} to refer to specific anchors.
#' By default, the \code{GenomicInteractions} constructor will name the anchors as \code{"first"} and \code{"second"},
#' so these can be used instead of \code{type=1} and \code{type=2} if the names have not been changed since.
#' If the names have been changed (e.g., with \code{\link{partnerNames<-}}, \code{type} can be used to refer to those custom names.
#' 
#' For backwards compatibility, \code{type="both"} has the same behavior as \code{type=NULL}.
#' This will be deprecated.
#' 
#' \code{id=TRUE} will return the integer indices pointing to the anchor regions in \code{regions(x)}.
#' If \code{type=NULL}, this will be a \linkS4class{DataFrame} of integer vectors, with one column per partner.
#' Otherwise, the appropriate integer vector is directly returned for the partner specified by \code{type}.
#' }
#' \item{\code{regions(x, type=NA)}:}{Returns the universe of all regions corresponding to anchor region \code{type}.
#' \code{type=1} or \code{2} returns a GenomicRanges object of all regions used in the first or second anchors, respectively.
#' \code{type=NULL} returns a GenomicRangesList of length 2, containing the regions used in the first and second anchors.
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
#' All getter methods applicable to \linkS4class{IndexedRelations} objects can also be used, 
#' e.g., \code{\link{partners}}, \code{\link{featureSets}}. 
#' 
#' @section Setters:
#' In the following code snippets, \code{x} is a GenomicInteractions object:
#' \describe{
#' \item{\code{anchors(x, type=NULL, id=FALSE) <- value}:}{Replaces the anchor regions specified by \code{type} with \code{value}.
#' If \code{type=NULL}, \code{value} should be a Pairs of two GenomicRanges,
#' where the first and second entries are to be used as the replacement first and second anchor regions, respectively.
#' If \code{type=1} or \code{type=2}, \code{value} should be a GenomicRanges to replace the first or second anchors, respectively.
#'
#' \code{type} can also be character strings such as \code{"first"} or \code{"second"}, depending on \code{partnerNames(x)}.
#' See above for related comments on \code{anchors}.
#' Again, \code{type="both"} can be used in place of \code{type=NULL}, though this is deprecated behavior. 
#' 
#' \code{id=TRUE} will replace the integer indices pointing to the anchor regions in \code{regions(x)}.
#' If \code{type=NULL}, \code{value} should be a \linkS4class{DataFrame} of integer vectors, with one column per partner.
#' Otherwise, it should be an integer vector to replace the indices for the partner specified by \code{type}.
#' }
#' \item{\code{regions(x, type=NA) <- value}:}{Replaces the universe of anchor regions in \code{x}.
#' If \code{type=1} or \code{2}, \code{value} should be a GenomicRanges object
#' which will replace the regions used for the first or second anchors, respectively.
#' If \code{type=NULL}, \code{value} should be a GenomicRangesList of length 2 to replace both region sets simultaneously.
#' 
#' \code{type=NA} is provided for backwards compatibility, 
#' and will replace the region sets for \emph{both} anchors with the GenomicRanges \code{value}.
#' This will trigger a deprecation warning - 
#' users wanting this behaviour should use \code{type=NULL} with the replacement value set to \code{List(value, value)}.
#' }
#' \item{\code{first(x) <- value}:}{A synonym for \code{anchors(x, 1) <- value}.}
#' \item{\code{second(x) <- value}:}{A synonym for \code{anchors(x, 2) <- value}.}
#' }
#' Again, all setter methods applicable to IndexedRelations objects can also be used, 
#' e.g., \code{\link{partners<-}}, \code{\link{featureSets<-}}.
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
#' GenomicInteractions-class
#' anchors anchors,GenomicInteractions-method anchors<- anchors<-,GenomicInteractions-method
#' regions regions,GenomicInteractions-method regions<- regions<-,GenomicInteractions-method
#' first,GenomicInteractions-method second,GenomicInteractions-method
#' first<-,GenomicInteractions-method second<-,GenomicInteractions-method
NULL

#' @export
#' @importFrom S4Vectors DataFrame mcols<- metadata<-
#' @importFrom BiocGenerics match
GenomicInteractions <- function(anchor1, anchor2, regions, ..., metadata=list(), common=TRUE) { 
    if (missing(regions) && common) {
        regions <- unique(sort(c(anchor1, anchor2)))
        anchor1 <- match(anchor1, regions)
        anchor2 <- match(anchor2, regions)
    }

    grf1 <- .create_GRF(anchor1, regions, "1")
    grf2 <- .create_GRF(anchor2, regions, "2")
    out <- new("GenomicInteractions", first=grf1$factor, second=grf2$factor)

    meta <- c(list(...), grf1$mcol, grf2$mcol)
    if (length(meta)) {
        mcols(out) <- do.call(DataFrame, meta)
    }
    metadata(out) <- as.list(metadata)
    out 
}

#' @importFrom GenomicRanges GRangesFactor
#' @importFrom S4Vectors mcols mcols<-
.create_GRF <- function(anchor, regions, n) {
    if (is(anchor, "GenomicRanges")) {
        mcol <- mcols(anchor)
        mcols(anchor) <- NULL
        colnames(mcol) <- sprintf("anchor%s.%s", n, colnames(mcol1))
        
        if (missing(regions)) {
            out <- GRangesFactor(anchor)
        } else {
            out <- GRangesFactor(anchor, regions)
        }
    } else if (is.numeric(anchor)) {
        out <- GRangesFactor(levels=regions, index=anchor)
        mcol <- NULL
    }
    list(factor=out, mcol=mcol)
}

#' @importFrom S4Vectors parallelSlotNames
setMethod("parallelSlotNames", "GenomicInteractions", function(x) {
    c("first", "second", callNextMethod())    
})

###############################################################################
###############################################################################
# Basic getters and setters.

.convert_type <- function(type) {
    match.arg(type, c("first", "second", "both"))
}

#' @export
#' @importFrom S4Vectors List first second
setMethod("anchors", "GenomicInteractions", function(x, type=NULL, id=FALSE) {
    type <- .convert_type(type)
    if (type=="both") {
        a1 <- first(x)
        a2 <- second(x)
        if (id) {
            return(List(first=as.integer(a1), second=as.integer(a2)))
        } else {
            return(List(first=a1, second=a2))
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
#' @importFrom S4Vectors first second
setMethod("anchors<-", "GenomicInteractions", function(x, type=NULL, id=FALSE, ..., value) {
    type <- .convert_type(type)
    if (type=="both") {
        if (id) {
            first(x) <- GRangesFactor(levels=levels(first(x)), values[[1]])
            second(x) <- GRangesFactor(levels=levels(second(x)), values[[2]])
        } else {
            first(x) <- value[[1]]
            second(x) <- value[[2]]
        }
    } else if (type=="first") {
        if (id) {
            first(x) <- GRangesFactor(levels=levels(first(x)), values)
        } else {
            first(x) <- values
        }
    } else {
        if (id) {
            second(x) <- GRangesFactor(levels=levels(second(x)), values)
        } else {
            second(x) <- values
        }
    }
    x
})

#' @export
#' @importFrom IndexedRelations featureSets
setMethod("regions", "GenomicInteractions", function(x, type=NA) {
    if (is.null(type)) {
        featureSets(x)
    } else {
        if (is.na(type)) {
            type <- 1L
            .Deprecated(msg="'regions(..., type=NA)' is deprecated.\nSee '?regions' for alternatives.")
        }
        featureSets(x)[[type]]
    }
})

#' @export
#' @importFrom IndexedRelations featureSets<-
setMethod("regions<-", "GenomicInteractions", function(x, type=NA, ..., value) {
    if (is.null(type)) {
        featureSets(x) <- value
    } else {
        if (is.na(type)) {
            .Deprecated(msg="'regions(..., type=NA)<-' is deprecated.\nSee '?regions' for alternatives.")
            featureSets(x) <- List(value, value)
        } else {
            featureSets(x)[[type]] <- value
        }
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
