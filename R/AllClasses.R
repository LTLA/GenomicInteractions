#' @export
#' @import methods
#' @importClassesFrom IndexedRelations IndexedRelations
setClass("GenomicInteractions", contains="IndexedRelations", slots=c(featureSets="SimpleGenomicRangesList"))

#' @export
setClass("SimpleGenomicInteractions", contains="GenomicInteractions")
