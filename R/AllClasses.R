#' @export
#' @import methods
#' @importClassesFrom S4Vectors Pairs 
#' @importClassesFrom GenomicRanges GRangesFactor
setClass("GenomicInteractions", contains="Pairs", slots=c(first="GRangesFactor", second="GRangesFactor"))
