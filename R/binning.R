#' @importFrom stats quantile

# get bin boundaries for bin(x, binmode="equalN", ...)
.breaksEqualN <- function(x, nElements, minAbsX = NULL) {
    if (!is.null(minAbsX)) {
        n1 <- round(sum(x < -minAbsX) / nElements) * nElements
        n2 <- round(sum(x >  minAbsX) / nElements) * nElements
        x1 <- sort(x, decreasing = FALSE)[1:n1]
        x2 <- sort(x, decreasing = TRUE)[1:n2]
        bin.breaks <- c(quantile(x1, seq(0, 1, length.out = n1 / nElements + 1)),
                        quantile(x2, seq(0, 1, length.out = n2 / nElements + 1)))
    } else {
        bin.breaks <- quantile(x, seq(0, 1, length.out = round(length(x) / nElements) + 1L))
    }
    bin.breaks
}

# get bin boundaries for bin(x, binmode="equalWidth", ...)
.breaksEqualWidth <- function(x, nBins, minAbsX) {
    seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = nBins + 1)
}

#' @title Bin elements of \code{x}.
#'
#' @description \code{bin} groups elements of \code{x} into bins with either
#'     a constant number of elments per bin, a constant bin width or according
#'     to user-provided bin boundaries.
#'
#' @param x A numerical vector with the values used for binning.
#' @param binmode The algorithm to be used for binning. Possible values are:
#'     "equalN" (default), "equalWidth" or "breaks" (see Details).
#' @param nElements The number of elements per bin (only for \code{binmode="equalN"}).
#'     The width of bins is adjusted accordingly.
#' @param nBins The number of bins (only for \code{binmode="equalWidth"}).
#'     The number of elments per bin will be variable.
#' @param minAbsX The minimal absolute value in \code{x} for exlements to be binned
#'     using the \code{binmode="equalN"} (ignored for other values of \code{binmode}).
#'     Elements with \code{x} values in \code{[-minAbsX,minAbsX]} will be collected in a single bin.
#' @param breaks Numerical vector with bin boundaries (only for \code{binmode="breaks"}).
#'     \code{breaks} has to be ordered and strictly increasing, and has to be of
#'     length (number of bins) + 1.
#'
#' @details Elements are binned according to the values in \code{x} depending on
#'     \code{binmode}:
#'     \describe{
#'         \item{equalN}{Items are grouped into a variable number of bins with
#'             \code{nElements} elements each. If \code{minAbsX} is not \code{NULL},
#'             elements with \code{x}-values in \code{[-minAbsX,minAbsX]} will
#'             first collected in a single bin before binning the remaining elements.
#'             The boundaries of this single bin may be slightly adjusted in order
#'             to respect the \code{nElements} elements per bin.}
#'         \item{equalWidth}{Items are group into \code{nBins} bins with a variable
#'             number of elements each.}
#'         \item{breaks}{Items are grouped into bins using \code{cut(x, breaks, include.lowest = TRUE)}}
#'     }
#'
#' @seealso \code{\link{cut}} which is used internally.
#'
#' @return A factor of the same length as \code{x}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' summary(bin(x, "equalN", nElements=20))
#' summary(bin(x, "equalN", nElements=20, minAbsX=0.5))
#' summary(bin(x, "equalWidth", nBins=5))
#' summary(bin(x, "breaks", breaks=c(-10,-1,0,1,10)))
#'
#' @export
bin <- function(x, binmode = c("equalN", "equalWidth", "breaks"),
                nElements = NULL, nBins = NULL, minAbsX = NULL, breaks = NULL) {
    stopifnot(is.numeric(x) && length(x) > 1)
    binmode <- match.arg(binmode)
    stopifnot(is.null(nElements) || (is.numeric(nElements) && length(nElements) == 1))
    stopifnot(is.null(nBins) || (is.numeric(nBins) && length(nBins) == 1))
    breaks <- switch(binmode,
                     "equalN" = .breaksEqualN(x, nElements, minAbsX),
                     "equalWidth" = .breaksEqualWidth(x, nBins, minAbsX),
                     "breaks" = breaks)
    cut(x, breaks = breaks, include.lowest = TRUE)
}

