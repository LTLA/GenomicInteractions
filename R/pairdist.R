#' Get the linear distance for each interaction
#'
#' Compute the distance between interacting regions on the linear genome, for each pairwise interaction contained in a \linkS4class{GenomicInteractions} object.
#'
#' @param x A \linkS4class{GenomicInteractions} object.
#' @param type String pecifying the type of distance to compute.
#' Can take values of \code{"mid"}, \code{"gap"}, \code{"span"}, \code{"diag"} or \code{"intra"}.
#'
#' @return
#' An integer vector of base-pair distances for \code{type="gap"} or \code{"span"}.
#' 
#' An integer vector of the number of intervening regions for \code{type="diag"}.
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
#' \item{\code{"diag"}:}{The number of regions in \code{regions(x)} between the two anchors is returned.
#' This requires that both anchors have the same universe of regions (an error is raised otherwise).}
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
#' @importFrom BiocGenerics start end
#' @importFrom GenomeInfoDb seqnames
#' @aliases pairdist pairdist,GenomicInteractions-method
setMethod("pairdist", "GenomicInteractions", function(x, type="mid") {
    type <- match.arg(type, c("mid", "gap", "span", "diag", "intra"))

    # Getting locations.
    a1 <- first(x)
    a2 <- second(x)
    f1 <- levels(a1)
    f2 <- levels(a2)
    i1 <- as.integer(a1)
    i2 <- as.integer(a2)

    # Protection when all inter's.
    c1 <- as.character(seqnames(f1))[i1]
    c2 <- as.character(seqnames(f2))[i2]
    is.same <- c1==c2
    if (type=="intra") { 
        return(is.same)
    }

    if (type=="mid") {
        default <- NA_real_
    } else {
        default <- NA_integer_
    }
    output <- rep(default, length(x))
    if (!any(is.same)) { 
        return(output) 
    }

    # Special behaviour with diagonal type.
    if (type=="diag") {
        o1 <- order(f1)
        o2 <- order(f2)
        if (length(f1)!=length(f2) || !all(f1[o1]==f2[o2])) {
            stop("'type=\"diag\"' only makes sense when anchors have the same region set")
        }

        re1 <- re2 <- integer(length(f1))
        re1[o1] <- seq_along(re1)
        re2[o2] <- seq_along(re2)

        output[is.same] <- abs(re1[i1[is.same]] - re2[i2[is.same]])
        return(output)
    }

    s1 <- start(f1)[i1[is.same]]
    s2 <- start(f2)[i2[is.same]]
    e1 <- end(f1)[i1[is.same]]
    e2 <- end(f2)[i2[is.same]]

    if (type=="gap") {
        output[is.same] <- pmax(s1, s2) - pmin(e1, e2) - 1L
    } else if (type=="span") {
        output[is.same] <- pmax(e1, e2) - pmin(s1, s2) + 1L
    } else if (type=="mid") {
        output[is.same] <- abs(s1 + e1 - s2 - e2)/2
    }
    output
})

#' @export
#' @rdname pairdist
intrachr <- function(x) pairdist(x, type="intra")
