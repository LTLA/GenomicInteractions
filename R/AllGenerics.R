#' @export
setGeneric("anchors", function(x, ...) standardGeneric("anchors"))

#' @export
setGeneric("anchors<-", function(x, ..., value) standardGeneric("anchors<-"))

#' @export
setGeneric("regions", function(x, ...) standardGeneric("regions"))

#' @export
setGeneric("regions<-", function(x, ..., value) standardGeneric("regions<-"))

#' @export
setGeneric("pairdist", function(x, ...) standardGeneric("pairdist"))

#' @export
setGeneric("boundingBox", function(x, ...) standardGeneric("boundingBox"))

#' @export
setGeneric("linkOverlaps", function(query, subject1, subject2, ...) standardGeneric("linkOverlaps"))

#' @export
setGeneric("swapAnchors", function(x, ...) standardGeneric("swapAnchors"))

#' @export
setGeneric("findDistalAnchors", function(x, ...) standardGeneric("findDistalAnchors"))
