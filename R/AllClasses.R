#' @export
#' @import methods
#' @importClassesFrom S4Vectors Vector 
#' @importClassesFrom GenomicRanges GRangesFactor
setClass("GenomicInteractions", contains="Vector", slots=c(first="GRangesFactor", second="GRangesFactor"))
