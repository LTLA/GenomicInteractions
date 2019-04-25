#' Update old GenomicInteractions 
#' 
#' Update an old GenomicInteractions instance (usually saved to file in serialized form)
#' into a backwards-compatible \linkS4class{GenomicInteractions} object.
#'
#' @param object Any instance of an old GenomicInteractions class, 
#' generated prior to version 2.0.0 of the \pkg{GenomicInteractions} package.
#' @param ... Further arguments, not used.
#' @param verbose Logical scalar indicating whether to display messages during conversion.
#'
#' @return A \linkS4class{GenomicInteractions} object.
#' 
#' @details
#' We return a \linkS4class{GenomicInteractions} instance rather than a \linkS4class{GenomicInteractions} instance,
#' simply to maximize backwards-compatibility with code that calls \code{\link{regions}} 
#' and expects a \linkS4class{GenomicRanges} in return.
#' 
#' @author Aaron Lun,
#' based on code by Malcolm Perry and Liz Ing-Simmons.
#' @examples
#' \dontrun{
#' updateObject(readRDS('some_old_objects.rds'))
#' }
#'
#' @export
#' @rdname updateObject
#' @importFrom S4Vectors DataFrame mcols<- metadata<-
#' @importFrom BiocGenerics updateObject getObjectSlots unique sort
setMethod("updateObject", "GenomicInteractions", function(object, ..., verbose = FALSE){
    if (verbose) message("updating GenomicInteractions object")
    
    all.slot.names <- names(getObjectSlots(object))
    if ("anchor_one" %in% all.slot.names){
        anchor1 <- object@anchor_one
        anchor2 <- object@anchor_two
        mcols(anchor1) <- NULL
        mcols(anchor2) <- NULL

        em <- DataFrame(counts = object@counts, object@elementMetadata)
        meta <- object@metadata

        combined <- unique(sort(c(anchor1, anchor2)))
        object <- GenomicInteractions(match(anchor1, combined), match(anchor2, combined), combined)
        mcols(object) <- em
        metadata(object) <- meta
    } else if ("anchor1" %in% all.slot.names) {
        anchor1 <- object@anchor1
        anchor2 <- object@anchor2
        regions <- object@regions

        em <- object@elementMetadata
        meta <- object@metadata
        
        object <- GenomicInteractions(anchor1, anchor2, regions)
        mcols(object) <- em
        metadata(object) <- meta
    }

    object
})
