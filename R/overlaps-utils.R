#######################################
### One-dimensional overlap methods ###
#######################################

#' @importFrom IndexedRelations featureSets partners
.get_used_regions <- function(x, i) 
# Only search for overlaps to the used subset of regions.
# Avoids wasting time if only a few regions are used.
{
    freg <- featureSets(x)[[i]]
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
    olap <- findOverlaps(used$region, subject, ...)

    query_idx <- used$index
    oq <- order(query_idx)
    subject_idx <- os <- seq_along(subject)

    out <- expand_1D_hits(queryHits(olap), subjectHits(olap), query_idx[oq], oq, subject_idx, os)
    Hits(out[[1]], out[[2]], length(query), length(subject), sort.by.query=TRUE)
}

#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits Hits
.find_single_overlap_right <- function(query, i, subject, ...) {
    used <- .get_used_regions(subject, i)
    olap <- findOverlaps(query, used$region, ...)

    query_idx <- oq <- seq_along(query)
    subject_idx <- used$index
    os <- order(subject_idx)

    out <- expand_1D_hits(queryHits(olap), subjectHits(olap), query_idx, oq, subject_idx[os], os)
    Hits(out[[1]], out[[2]], length(query), length(subject), sort.by.query=TRUE)
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

    selectHits(hits, select = select)
}

#######################################
### Two-dimensional overlap methods ###
#######################################

.options_2d <- c("any", "match", "reverse")

#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits Hits
.find_double_overlap <- function(query, subject, ..., do.same=TRUE, do.reverse=TRUE, find.arbitrary=FALSE) {
    q_used_left <- .get_used_regions(query, 1)
    q_used_right <- .get_used_regions(query, 2)
    q_reg_left <- q_used_left$region
    q_ind_left <- q_used_left$index
    q_reg_right <- q_used_right$region
    q_ind_right <- q_used_right$index

    s_used_left <- .get_used_regions(subject, 1)
    s_used_right <- .get_used_regions(subject, 2)
    s_reg_left <- s_used_left$region
    s_ind_left <- s_used_left$index
    s_reg_right <- s_used_right$region
    s_ind_right <- s_used_right$index

    ol <- order(s_ind_left)
    or <- order(s_ind_right)
    s_ind_left <- s_ind_left[ol]
    s_ind_right <- s_ind_right[or]

    if (do.same) {
        left_olap <- findOverlaps(q_reg_left, s_reg_left, ...)
        right_olap <- findOverlaps(q_reg_right, s_reg_right, ...)

        out.same <- collate_2D_hits(q_ind_left, q_ind_right,
            queryHits(left_olap), subjectHits(left_olap),
            queryHits(right_olap), subjectHits(right_olap),
            s_ind_left, ol, s_ind_right, or, 
            find.arbitrary)
        hits.same <- Hits(out.same[[1]], out.same[[2]], length(query), length(subject), sort.by.query=TRUE)
    }

    if (do.reverse){ 
        left_olap <- findOverlaps(q_reg_right, s_reg_left, ...)
        right_olap <- findOverlaps(q_reg_left, s_reg_right, ...)

        out.rev <- collate_2D_hits(q_ind_right, q_ind_left,
            queryHits(left_olap), subjectHits(left_olap),
            queryHits(right_olap), subjectHits(right_olap),
            s_ind_left, ol, s_ind_right, or,
            find.arbitrary)
        hits.rev <- Hits(out.rev[[1]], out.rev[[2]], length(query), length(subject), sort.by.query=TRUE)
    }

    if (do.same && do.reverse){
        union(hits.same, hits.rev)
    } else if (do.same) {
        hits.same
    } else {
        hits.rev
    }
}

#######################################
### 1.5-dimensional overlap methods ###
#######################################

.options_1.5d <- c("any-any", 
    "first-any", "second-any", 
    "any-first", "any-second",
    "first-first", "first-second", 
    "second-first", "second-second")

.process_anchors <- function(IR, i) {
    used <- .get_used_regions(IR, i)
    i <- used$index
    o <- order(i)
    used$index <- i[o]
    used$order <- o
    used
}

#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits Hits
#' @importFrom BiocGenerics union sort
.find_single_overlap_IR <- function(query, subject, ..., query.left=TRUE,
    query.right=TRUE, subject.left=TRUE, subject.right=TRUE)
{
    q_info <- list()
    if (query.left) {
        q_info$left <- .process_anchors(query, 1)
    }
    if (query.right) {
        q_info$right <- .process_anchors(query, 2)
    }

    s_info <- list()
    if (subject.left) {
        s_info$left <- .process_anchors(subject, 1)
    }
    if (subject.right) {
        s_info$right <- .process_anchors(subject, 2)
    }

    collected_hits <- list()
    for (q in seq_along(q_info)) {
        cur_q_region <- q_info[[q]]$region
        cur_q_index <- q_info[[q]]$index
        cur_q_order <- q_info[[q]]$order

        for (s in seq_along(s_info)) {
            cur_s_region <- s_info[[s]]$region
            cur_s_index <- s_info[[s]]$index
            cur_s_order <- s_info[[s]]$order

            olap <- findOverlaps(cur_q_region, cur_s_region, ...)
            out <- expand_1D_hits(queryHits(olap), subjectHits(olap), cur_q_index, cur_q_order, cur_s_index, cur_s_order)
            hits <- Hits(out[[1]], out[[2]], length(query), length(subject), sort.by.query=TRUE)
            collected_hits <- append(collected_hits, list(hits))
        }
    }

    Reduce(union, collected_hits)
}
