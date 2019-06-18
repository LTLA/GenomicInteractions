# The show method involves enough submethods that I've just
# put it into a separate file; this is documented in GenomicInteractions-class.R

#' @importFrom S4Vectors mcols showAsCell
.makeNakedMatFromGenomicInteractions <- function(x) {
    ans <- cbind(first=showAsCell(first(x)), second=showAsCell(second(x)))

    x_mcols <- mcols(x, use.names = FALSE)
    if (!is.null(x_mcols) && ncol(x_mcols) > 0L) {
        tmp <- do.call(data.frame, c(lapply(x_mcols, showAsCell), list(check.names = FALSE)))
        ans <- cbind(ans, `|` = rep.int("|", length(x)), as.matrix(tmp))
    }

    ans
}

#' @importFrom S4Vectors mcols
showGenomicInteractions <- function(x, margin = "", print.classinfo = FALSE) {
    x_class <- class(x)
    x_len <- length(x)
    x_mcols <- mcols(x, use.names = FALSE)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    cat(x_class, " object with ", 
        x_len, " relation", ifelse(x_len == 1L, "", "s"), " and ", 
        x_nmc, " metadata column", ifelse(x_nmc == 1L, "", "s"), 
        ":\n", sep = "")
    
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromGenomicInteractions)

    if (print.classinfo) {
        .COL2CLASS <- c(first=class(first(x)), second=class(second(x)))
        classinfo <- S4Vectors:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L) {
        rownames(out) <- paste0(margin, rownames(out))
    }

    print(out, quote = FALSE, right = TRUE, max = length(out))
}

#' @export
#' @importFrom methods show
setMethod("show", "GenomicInteractions", function(object) {
    showGenomicInteractions(object, margin="  ", print.classinfo=TRUE)
})
