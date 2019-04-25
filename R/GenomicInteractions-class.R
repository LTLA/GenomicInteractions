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
#' \code{GenomicInteractions(anchor1, anchor2, regions, ...)} will create a GenomicInteractions object, given:
#' \itemize{
#' \item A \linkS4class{GenomicRanges} object in each of \code{anchor1} and \code{anchor2}, and missing \code{regions}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, 
#' and a \linkS4class{SimpleGenomicRangesList} of length 2 in \code{regions}.
#' \item An integer vector in each of \code{anchor1} and \code{anchor2}, and a GenomicRanges in \code{regions}.
#' }
#'
#' Note that only the last option will yield a GenomicInteractions object with a single GenomicRanges as the feature set.
#' This is useful in situations where the first and second anchors have the same origin, e.g., genomic bins.
#' For the other options, the output object will contain a List of two GenomicRanges (one per anchor region).
#' This is useful in situations where the anchors are distinct, e.g., promoters versus enhancers.
#'
#' Any arguments in \code{...} are added to the element-wise metadata of the output object.
#'
#' @section Getters:
#' In the following code snippets, \code{x} is a GenomicInteractions object:
#' \describe{
#' \item{\code{anchors(x, type)}:}{Returns the anchor regions specified by \code{type} as a GenomicRanges object.
#' \code{type} should be an integer scalar specifying the first (\code{1}) or second anchor (\code{2}).}
#' \item{\code{regions(x, as.list=NULL)}:}{If \code{as.list=FALSE}, returns the feature sets of all regions as a SimpleGenomicRangesList object.
#' The returned object can be of length 1 or 2, depending on how \code{x} was constructed.
#'
#' If \code{as.list=NULL} or \code{TRUE}, only the first feature set is returned.
#' This is the short-term default for backwards compatibility,
#' but is deprecated as \code{x} may contain multiple features - in which case, returning only the first is rarely useful.}
#' \item{\code{first(x)}:}{A synonym for \code{anchors(x, 1)}.}
#' \item{\code{second(x)}:}{A synonym for \code{anchors(x, 2)}.}
#' }
#' All getter methods applicable to IndexedRelations objects can also be used, e.g., \code{\link{partners}}, \code{\link{featureSets}}, \code{\link{mapping}}.
#' 
#' @section Setters:
#' In the following code snippets, \code{x} is a GenomicInteractions object:
#' \describe{
#' \item{\code{anchors(x, type) <- value}:}{Replaces the anchor regions specified by \code{type} with the entries of \code{value}, a GenomicRanges object.
#' \code{type} should be an integer scalar specifying the first (\code{1}) or second anchor (\code{2}).}
#' \item{\code{regions(x, as.list=NULL) <- value}:}{If \code{as.list=TRUE}, replaces the feature sets of all regions as with \code{value}, a SimpleGenomicRangesList object.
#' \code{lengths(value)} should be of the same as \code{lengths(regions(x))} to ensure validity of the modified object.
#' 
#' If \code{as.list=NULL} or \code{TRUE}, \code{value} should be a GenomicRanges object that is used to replace the first feature set only.
#' This is the short-term default for backwards compatibility,
#' but is deprecated as \code{x} may contain multiple region sets - in which case, replacing the first alone is rarely useful.}
#' \item{\code{first(x) <- value}:}{A synonym for \code{anchors(x, 1) <- value}.}
#' \item{\code{second(x) <- value}:}{A synonym for \code{anchors(x, 2) <- value}.}
#' }
#' Again, all setter methods applicable to IndexedRelations objects can also be used, e.g., \code{\link{partners<-}}, \code{\link{featureSets<-}}.
#'
#' @author Aaron Lun
#' @seealso
#' \linkS4class{IndexedRelations}, from which this class is derived.
#'
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' # Getting
#' anchors(test, 1)
#' anchors(test, 2) 
#' regions(test, as.list=TRUE)
#' 
#' # Setting
#' anchors(test, 1) <- rev(anchors(test, 1))
#' test
#'
#' # Set 'as.list=TRUE' in getters/setters to avoid deprecation warnings.
#' first.regions <- regions(test, as.list=TRUE)[[1]]
#' regions(test, as.list=TRUE)[[1]] <- resize(first.regions, 1000)
#' test
#'
#' @rdname GenomicInteractions
#' @name GenomicInteractions
#' @aliases GenomicInteractions
#' anchors anchors,GenomicInteractions-method anchors<- anchors<-,GenomicInteractions-method
#' regions regions,GenomicInteractions-method regions<- regions<-,GenomicInteractions-method
#' first,GenomicInteractions-method second,GenomicInteractions-method
#' first<-,GenomicInteractions-method second<-,GenomicInteractions-method
NULL

#' @export
#' @importFrom S4Vectors List mcols<- DataFrame
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importFrom IndexedRelations IndexedRelations
GenomicInteractions <- function(anchor1, anchor2, regions, ...) { 
    if (missing(regions)) {
        out <- IndexedRelations(list(anchor1, anchor2))
    } else if (is(regions, "GenomicRanges")) {
        out <- IndexedRelations(list(anchor1, anchor2), List(regions), mapping=c(1L, 1L))
    } else {
        out <- IndexedRelations(list(anchor1, anchor2), regions)
    }

    meta <- list(...)
    if (length(meta)) {
        mcols(out) <- do.call(DataFrame, meta)
    }
    new("GenomicInteractions", out)
}

#' @importFrom S4Vectors setValidity2
setValidity2("GenomicInteractions", function(object) {
    msg <- character(0)

    if (npartners(object)!=2L) {
        stop("GenomicInteractions objects must contain pairwise interactions")
    }

    if (length(msg)) {
        return(msg)
    }
    TRUE
})

###############################################################################
###############################################################################
# Basic getters and setters.

#' @export
#' @importFrom IndexedRelations partnerFeatures
setMethod("anchors", "GenomicInteractions", function(x, type) partnerFeatures(x, type))

#' @export
#' @importFrom IndexedRelations partnerFeatures<-
setMethod("anchors<-", "GenomicInteractions", function(x, type, ..., value) {
    partnerFeatures(x, type) <- value
    x
})

#' @export
#' @importFrom IndexedRelations featureSets
setMethod("regions", "GenomicInteractions", function(x, as.list=NULL) {
    if (is.null(as.list) || !as.list) {
        if (is.null(as.list)) {
            .Deprecated(msg="'as.list=NULL' is deprecated.\nSee '?regions' for details.")
        } else {
            .Deprecated(msg="'as.list=FALSE' is deprecated.\nSee '?regions' for details.")
        }
        featureSets(x)[[1]]    
    } else {
        featureSets(x)
    }
})

#' @export
#' @importFrom IndexedRelations featureSets<-
setMethod("regions<-", "GenomicInteractions", function(x, as.list=NULL, ..., value) {
    if (is.null(as.list) || !as.list) {
        if (is.null(as.list)) {
            .Deprecated(msg="'as.list=NULL' is deprecated.\nSee '?regions' for details.")
        } else {
            .Deprecated(msg="'as.list=FALSE' is deprecated.\nSee '?regions' for details.")
        }
        featureSets(x)[[1]] <- value
    } else {
        featureSets(x) <- value
    }
    x
})

#' @export
#' @importFrom S4Vectors first
setMethod("first", "GenomicInteractions", function(x) anchors(x, 1))

#' @export
#' @importFrom S4Vectors first<-
setReplaceMethod("first", "GenomicInteractions", function(x, ..., value) {
    anchors(x, 1) <- value
    x
})

#' @export
#' @importFrom S4Vectors second
setMethod("second", "GenomicInteractions", function(x) anchors(x, 2))

#' @export
#' @importFrom S4Vectors second<-
setReplaceMethod("second", "GenomicInteractions", function(x, ..., value) {
    anchors(x, 2) <- value
    x
})
