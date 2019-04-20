#' Get the linear distance for each interaction
#'
#' Compute the distance between interacting regions on the linear genome, for each pairwise interaction contained in a \linkS4class{GenomicInteractions} object.
#'
#' @param x A \linkS4class{GenomicInteractions} object.
#' @param type String pecifying the type of distance to compute.
#' Can take values of \code{"mid"}, \code{"gap"}, \code{"span"}, \code{"diag"} or \code{"intra"}.
#'
#' @param
#' An integer vector of base-pair distances for \code{type="gap"} or \code{"span"}.
#' 
#' An integer vector of \dQuote{unit-based} distances for \code{type="diag"}.
#' 
#' A numeric vector of base-pair distances for \code{type="mid"}.
#'
#' A logical vector for \code{type="intra"} or from \code{intrachr}.
#'
#' @details
#' For each interaction in \code{x}, the \code{pairdist} method computes the distance between the two interacting regions.
#' An integer vector is returned, with values computed according to the specified value of \code{type}:
#' \describe{
#' \item{\code{"mid"}:}{The distance between the midpoints of the two regions (rounded down to the nearest integer) is returned.}
#' \item{\code{"gap"}:}{The length of the gap between the closest points of the two regions is computed - 
#' negative lengths are returned for overlapping regions, indicating the length of the overlap.}
#' \item{\code{"span"}:}{The distance between the furthermost points of the two regions is computed.}
#' \item{\code{"diag"}:}{The difference between the anchor indices is returned.
#' This corresponds to a diagonal on the interaction space when bins are used in \code{regions(x)}.
#' Note that this only makes sense when both anchors refer to the same set of regions.}
#' }
#' 
#' Interchromosomal interactions are marked with \code{NA}.
#' Alternatively, if \code{type="intra"}, a logical vector is returned indicating whether the interaction occurs between two regions on the same chromosome.
#' \code{intrachr(x)} is an alias for \code{pairdist(x, type="intra")}.
#' 
#' @examples
#' anchor1 <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor2 <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor1, anchor2)
#' test
#'
#' pairdist(test)
#' pairdist(test, type="gap")
#' pairdist(test, type="span")
#' intrachr(test)
#' 
#' # Common regions.
#' regions <- GRanges("chr1", IRanges(sort(sample(1000, 10)), width=5))
#' i1 <- sample(length(regions), 10)
#' i2 <- sample(length(regions), 10)
#' test.common <- GenomicInteractions(i1, i2, regions)
#' pairdist(test.common, type="diag")
#' 
#' @author
#' Aaron Lun
#' 
#' @export
#' @name pairdist
#' @importFrom IndexedRelations partners mapping featureSets
#' @importFrom BiocGenerics start end
#' @importFrom GenomeInfoDb seqnames
setMethod("pairdist", "GenomicInteractions", function(x, type="mid") {
    type <- match.arg(type, c("mid", "gap", "span", "diag", "intra"))

    # Getting locations.
    all.starts <- all.ends <- all.chr <- vector("list", 2L)
    for (i in seq_len(2L)) {
        chosen <- partners(x)[,i]
        ftype <- mapping(x)[i]
        regions <- featureSets(x)[[ftype]]
        all.starts[[i]] <- start(regions)[chosen]
        all.ends[[i]] <- end(regions)[chosen]
        all.chr[[i]] <- as.character(seqnames(regions)[chosen])
    }

    # Protection when all inter's.
    is.same <- all.chr[[1]]==all.chr[[2]]
    if (type=="intra") { 
        return(is.same)
    }

    output <- rep(NA_integer_, length(x))
    if (!any(is.same)) { 
        return(output) 
    }

    s1 <- all.starts[[1]][is.same]
    s2 <- all.starts[[2]][is.same]
    e1 <- all.ends[[1]][is.same]
    e2 <- all.ends[[2]][is.same]

    if (type=="gap") {
        output[is.same] <- pmax(s1, s2) - pmin(e1, e2) - 1L
    } else if (type=="span") {
        output[is.same] <- pmax(e1, e2) - pmin(s1, s2) + 1L
    } else if (type=="mid") {
        output[is.same] <- abs(s1 + e1 - s2 - e2)/2
    } else if (type=="diag") {
        if (length(featureSets(x))!=1L) {
            stop("'type=\"diag\"' only makes sense for objects with a single region set")
        }
        subpartners <- partners(x)[is.same,,drop=FALSE]
        output[is.same] <- subpartners[,1] - subpartners[,2]
    }
    output
})

#' @export
#' @rdname pairdist
intrachr <- function(x) pairdist(x, type="intra")