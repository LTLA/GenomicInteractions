#' @importFrom S4Vectors queryHits subjectHits
expand_single_overlap <- function(hits, query_idx, subject_idx) {
    oq <- order(query_idx)
    os <- order(subject_idx)
    out <- expand_hits(queryHits(hits), subjectHits(hits), query_idx[oq], subject_idx[os])
    list(query=oq[out[[1]]], subject=os[out[[2]]])
}

#' @importFrom IndexedRelations featureSets mapping partners
.get_used_regions <- function(x, i) 
# Only search for overlaps to the used subset of regions.
# Avoids wasting time if only a few regions are used.
{
    freg <- featureSets(x)[[mapping(x)[i]]]
    p <- partners(x)[,i]

    to.use <- logical(length(freg))
    to.use[p] <- TRUE
    if (!all(to.use)) {
        freg <- freg[to.use]
        p <- match(p, which(to.use))
    }

    list(region=freg, index=p)
}

#' @importFrom IRanges findOverlaps
.find_single_overlap_left <- function(query, i, subject, ...) {
    used <- .get_used_regions(query, i)
    olap <- findOverlaps(used$region, subject, ...)
    expand_single_overlap(olap, used$index, seq_along(subject))
}

#' @importFrom IRanges findOverlaps
.find_single_overlap_right <- function(query, i, subject, ...) {
    used <- .get_used_regions(subject, i)
    olap <- findOverlaps(query, used$region, ...)
    expand_single_overlap(olap, seq_along(query), used$index)
}
