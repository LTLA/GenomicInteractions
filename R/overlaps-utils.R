#######################################
### One-dimensional overlap methods ###
#######################################

#' @importFrom IndexedRelations featureSets mapping partners
.get_used_regions <- function(x, i) 
# Only search for overlaps to the used subset of regions.
# Avoids wasting time if only a few regions are used.
{
    freg <- featureSetByPartner(x, i)
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
#' @importFrom S4Vectors queryHits subjectHits Hits
.find_single_overlap_left <- function(query, i, subject, ...) {
    used <- .get_used_regions(query, i)
    hits <- findOverlaps(used$region, subject, ...)

    query_idx <- used$index
    oq <- order(query_idx)
    out <- expand_1D_hits(queryHits(hits), subjectHits(hits), query_idx[oq], seq_along(subject))
    Hits(oq[out[[1]]], out[[2]], length(query), length(subject))
}

#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits Hits
.find_single_overlap_right <- function(query, i, subject, ...) {
    used <- .get_used_regions(subject, i)
    hits <- findOverlaps(query, used$region, ...)

    subject_idx <- used$index
    os <- order(subject_idx)
    out <- expand_1D_hits(queryHits(hits), subjectHits(hits), seq_along(query), subject_idx[os])
    Hits(out[[1]], os[out[[2]]], length(query), length(subject))
}

#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors selectHits Hits
#' @importFrom BiocGenerics union sort
.find_single_overlap_wrapper <- function(query, subject, FUN, use.region, select, ...) {
    if (use.region=="both") {
        .Deprecated(msg="'use.region=\"both\" is deprecated.\nUse 'use.region=\"any\"' instead.")
        use.region <- "any"
    }
    if (use.region=="first" || use.region=="any") {
        hits1 <- FUN(query, 1L, subject, ...)
    } 
    if (use.region=="second" || use.region=="any") {
        hits2 <- FUN(query, 2L, subject, ...)
    }

    if (use.region=="any") { 
        hits <- union(hits1, hits2)
    } else if (use.region=="first") {
        hits <- hits1
    } else if (use.region=="second") {
        hits <- hits2
    }

    hits <- sort(hits)
    selectHits(hits, select = select)
}

#######################################
### Two-dimensional overlap methods ###
#######################################

#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits Hits
.find_double_overlap <- function(query, subject, ..., do.same=TRUE, do.reverse=TRUE) {
    qused_left <- .get_used_regions(query, 1)
    qused_right <- .get_used_regions(query, 2)
    qleft <- qused_left$index
    qright <- qused_right$index

    sused_left <- .get_used_regions(subject, 1)
    sused_right <- .get_used_regions(subject, 2)
    sleft <- sused_left$index
    sright <- sused_right$index

    # Ordering the indices.
    ol <- order(sleft)
    or <- order(sright)
    sleft <- sleft[ol]
    sright <- sright[or]

    if (do.same) {
        left_hits <- findOverlaps(qused_left$region, sused_left$region, ...)
        right_hits <- findOverlaps(qused_right$region, sused_right$region, ...)

        out.same <- collate_2D_hits(qleft, qright,
            queryHits(left_hits), subjectHits(left_hits),
            queryHits(right_hits), subjectHits(right_hits),
            sleft, ol, sright, or)
        hits.same <- Hits(out.same[[1]], out.same[[2]], length(query), length(subject))
    }

    if (do.reverse){ 
        left_rhits <- findOverlaps(qused_right$region, sused_left$region, ...)
        right_rhits <- findOverlaps(qused_left$region, sused_right$region, ...)

        out.rev <- collate_2D_hits(qright, qleft,
            queryHits(left_rhits), subjectHits(left_rhits),
            queryHits(right_rhits), subjectHits(right_rhits),
            sleft, ol, sright, or)
        hits.rev <- Hits(out.rev[[1]], out.rev[[2]], length(query), length(subject))
    }

    if (do.same && do.reverse){
        hits <- union(hits.same, hits.rev)
        sort(hits)
    } else if (do.same) {
        hits.same
    } else {
        hits.rev
    }
}

