#' @importFrom S4Vectors queryHits subjectHits
expand_single_overlap <- function(hits, query_idx, subject_idx) {
    oq <- order(query_idx)
    os <- order(subject_idx)
    out <- expand_hits(queryHits(hits), subjectHits(hits), query_idx[oq], subject_idx[os])
    list(query=oq[out[[1]]], subject=os[out[[2]]])
}

#' @importFrom IRanges findOverlaps
#' @importFrom IndexedRelations featureSets mapping partners
.find_single_overlap_left <- function(query, i, subject, ...) {
    freg <- featureSets(query)[[mapping(query)[i]]]

    # Only search for overlaps to the used subset of regions.
    to.use <- logical(length(freg))
    p <- partners(query)[,i]
    to.use[p] <- TRUE
    if (!all(to.use)) {
        freg <- freg[to.use]
        p <- match(p, which(to.use))
    }

    olap <- findOverlaps(freg, subject, ...)
    expand_single_overlap(olap, p, seq_along(subject))
}
