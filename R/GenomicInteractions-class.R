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
#' \code{GenomicInteractions(anchor1, anchor2, regions, ..., metadata=list())}
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
#' For backwards compatibility, it is also possible to set \code{type="first"}, \code{"second"} or \code{"both"}.
#' This has the same effect as \code{type=1}, {2} and \code{NULL}, respectively.
#' This behaviour is deprecated as it is possible to set alternative names for each anchor region,
#' and to retrieve them by using those names as \code{type} - see \code{\link{partnerNames}} for more details.
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
#' \code{type} can also be \code{"first"}, \code{"second"} or \code{"both"}, 
#' though this is deprecated behaviour - see above for related comments on \code{anchors}.
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
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' # Alternative methods of construction:
#' combined <- c(anchor1, anchor2)
#' m1 <- match(anchor1, combined)
#' m2 <- match(anchor1, combined)
#' test2 <- GenomicInteractions(m1, m2, combined)
#' test2
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
#' @importFrom S4Vectors List mcols<- DataFrame metadata<-
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom IndexedRelations IndexedRelations
#' @importFrom BiocGenerics match
GenomicInteractions <- function(anchor1, anchor2, regions, ..., metadata=list()) { 
    meta <- list(...)
    if (is(anchor1, "GenomicRanges")) {
        mcol1 <- mcols(anchor1)
        mcols(anchor1) <- NULL
        colnames(mcol1) <- sprintf("anchor1.%s", colnames(mcol1))
        meta <- c(meta, lapply(mcol1, I))
    }
    if (is(anchor2, "GenomicRanges")) {
        mcol2 <- mcols(anchor2)
        mcols(anchor2) <- NULL
        colnames(mcol2) <- sprintf("anchor2.%s", colnames(mcol2))
        meta <- c(meta, lapply(mcol2, I))
    }

    x <- list(first=anchor1, second=anchor2)
    if (missing(regions)) {
        out <- IndexedRelations(x)
    } else if (is(regions, "GenomicRanges")) {
        out <- IndexedRelations(x, List(regions, regions))
    } else {
        out <- IndexedRelations(x, regions)
    }

    if (length(meta)) {
        mcols(out) <- do.call(DataFrame, meta)
    }
    metadata(out) <- as.list(metadata)
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

.convert_type <- function(type) {
    if (is.character(type)) {
        type <- switch(type, first=1, second=2, both=NULL, type)
        if (!is.character(type)) {
            .Deprecated("'type=\"first\"', etc. is deprecated.\nUse integer 'type' instead.")
        }
    }
    type
}

#' @export
#' @importFrom IndexedRelations partnerFeatures
setMethod("anchors", "GenomicInteractions", function(x, type=NULL, id=FALSE) {
    type <- .convert_type(type)
    if (id) {
        out <- partners(x)
        if (!is.null(type)) {
            out <- out[,type]
        }
    } else {
        if (is.null(type)) {
            out <- as(x, "Pairs")
        } else {
            out <- partnerFeatures(x, type)
        }
    }
    out
})

#' @export
#' @importFrom IndexedRelations partnerFeatures<-
setMethod("anchors<-", "GenomicInteractions", function(x, type=NULL, id=FALSE, ..., value) {
    type <- .convert_type(type)
    if (id) {
        if (is.null(type)) {
            partners(x) <- value
        } else {
            partners(x)[,type] <- value
        }
    } else {
        if (is.null(type)) {
            partnerFeatures(x, 1) <- value[[1]]
            partnerFeatures(x, 2) <- value[[2]]
        } else {
            partnerFeatures(x, type) <- value
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
            warning("'type=NA' is deprecated.\nSee '?regions' for alternatives.")
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
            warning("'type=NA' is deprecated.\nSee '?regions' for alternatives.")
            featureSets(x)[[type]] <- List(value, value)
        } else {
            featureSets(x)[[type]] <- value
        }
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
