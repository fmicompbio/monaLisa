#' @importFrom stats quantile

# get bin boundaries for bin(x, binmode="equalN", ...)
.breaksEqualN <- function(x, nElements, minAbsX = NULL) {
    if (!is.null(minAbsX)) {
        n1 <- round(sum(x < -minAbsX, na.rm = TRUE) / nElements) * nElements
        n2 <- round(sum(x >  minAbsX, na.rm = TRUE) / nElements) * nElements
        x1 <- sort(x, decreasing = FALSE, na.last = NA)[seq_len(n1)]
        x2 <- sort(x, decreasing = TRUE,  na.last = NA)[seq_len(n2 + 1)]
        bin.breaks <- c(if (n1 > 0) quantile(x1, seq(0, 1, length.out = n1 / nElements + 1)) else min(x, na.rm = TRUE),
                        if (n2 > 0) quantile(x2, seq(0, 1, length.out = n2 / nElements + 1)) else max(x, na.rm = TRUE))
        attr(bin.breaks, "bin0") <- ceiling(n1 / nElements + 1)
    } else {
        bin.breaks <- quantile(x, seq(0, 1, length.out = round(length(x) / nElements) + 1L))
        attr(bin.breaks, "bin0") <- NA
    }
    bin.breaks
}

# get bin boundaries for bin(x, binmode="equalWidth", ...)
.breaksEqualWidth <- function(x, nBins, minAbsX) {
    if (!is.null(minAbsX)) {
        e1 <- sum(x < -minAbsX, na.rm = TRUE)
        e2 <- sum(x >  minAbsX, na.rm = TRUE)
        b1 <- round((nBins - 1) * e1 / (e1 + e2))
        b2 <- nBins - 1 - b1
        bw12 <- c((-minAbsX - min(x, na.rm = TRUE)) / b1,
                  (max(x, na.rm = TRUE) - minAbsX) / b2)
        bw <- max(bw12[is.finite(bw12)])
        almostOne <- (length(x) - 1) / length(x)
        bin.breaks <- c(if (b1 > 0) rev(seq(-minAbsX, min(x, na.rm = TRUE) - almostOne * bw, by = -bw)) else min(x, na.rm = TRUE),
                        if (b2 > 0) seq(minAbsX, max(x, na.rm = TRUE) + almostOne * bw, by = bw) else max(x, na.rm = TRUE))
        attr(bin.breaks, "bin0") <- b1 + 1
    } else {
        bin.breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = nBins + 1)
        attr(bin.breaks, "bin0") <- NA
    }
    bin.breaks
}

#' @title Bin elements of \code{x}.
#'
#' @description \code{bin} groups elements of \code{x} into bins with either a
#'   constant number of elments per bin, a constant bin width or according to
#'   user-provided bin boundaries.
#'
#' @param x A numerical vector with the values used for binning.
#' @param binmode The algorithm to be used for binning. Possible values are:
#'   "equalN" (default), "equalWidth" or "breaks" (see Details).
#' @param nElements The number of elements per bin (only for
#'   \code{binmode="equalN"}). The width of bins is adjusted accordingly.
#' @param nBins The number of bins (only for \code{binmode="equalWidth"}). The
#'   number of elments per bin will be variable.
#' @param minAbsX The minimal absolute value in \code{x} for elements to be
#'   binned using the \code{binmode="equalN"} or \code{binmode="equalWidth"}
#'   (ignored for other values of \code{binmode}). Elements with \code{x} values
#'   in \code{[-minAbsX,minAbsX]} will be collected in a single bin.
#' @param breaks Numerical vector with bin boundaries (only for
#'   \code{binmode="breaks"}). \code{breaks} has to be ordered and strictly
#'   increasing, and has to be of length (number of bins) + 1.
#' @param ... further arguments to be passed to \code{cut(x, breaks,
#'   include.lowest = TRUE, ...)}, such as \code{labels=FALSE}.
#'
#' @details Elements are binned according to the values in \code{x} depending on
#'   \code{binmode}: \describe{ \item{equalN}{Items are grouped into a variable
#'   number of bins with \code{nElements} elements each. If \code{minAbsX} is
#'   not \code{NULL}, elements with \code{x}-values in \code{[-minAbsX,minAbsX]}
#'   will first be collected in a single bin before binning the remaining
#'   elements. The boundaries of this single bin may be slightly adjusted in
#'   order to respect the \code{nElements} elements in the other bins.}
#'   \item{equalWidth}{Items are group into \code{nBins} bins with a variable
#'   number of elements each.} \item{breaks}{Items are grouped into bins using
#'   \code{cut(x, breaks, include.lowest = TRUE)}} }
#'
#' @seealso \code{\link{cut}} which is used internally.
#'
#' @return The return value from \code{cut(x, ...)}, typically a factor of the
#'   same length as \code{x}. Binning mode, bin boundaries and the "neutral" bin
#'   are available from \code{attr(..., "binmode")}, \code{attr(..., "breaks")}
#'   and \code{attr(..., "bin0")}. For \code{binmode = "breaks"}, the latter
#'   will be \code{NA}.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' summary(bin(x, "equalN", nElements=10))
#' summary(bin(x, "equalN", nElements=10, minAbsX=0.5))
#' summary(bin(x, "equalWidth", nBins=5))
#' summary(bin(x, "breaks", breaks=c(-10,-1,0,1,10)))
#'
#' @export
bin <- function(x, binmode = c("equalN", "equalWidth", "breaks"),
                nElements = round(length(x)/5), nBins = NULL, minAbsX = NULL, breaks = NULL,
                ...) {
    stopifnot(is.numeric(x) && length(x) > 1)
    binmode <- match.arg(binmode)
    stopifnot(is.null(nElements) || (is.numeric(nElements) && length(nElements) == 1))
    stopifnot(is.null(nBins) || (is.numeric(nBins) && length(nBins) == 1))
    stopifnot(is.null(minAbsX) || (is.numeric(minAbsX) && length(minAbsX) == 1 && minAbsX > 0))
    breaks <- switch(binmode,
                     "equalN" = .breaksEqualN(x, nElements, minAbsX),
                     "equalWidth" = .breaksEqualWidth(x, nBins, minAbsX),
                     "breaks" = breaks)
    res <- cut(x, breaks = breaks, include.lowest = TRUE, ...)
    attr(res, "binmode") <- binmode
    attr(res, "breaks") <- unname(breaks)
    attr(res, "bin0") <- if (binmode == "breaks") NA else attr(breaks, "bin0")
    res
}

#' Get and set the zero bin manually
#' 
#' @param bins Factor, typically the return value of \code{\link[monaLisa]{bin}}.
#' @param zeroBin Numeric or character scalar indicating the level to use as 
#'   the zero bin.
#' 
#' @examples 
#' set.seed(1)
#' x <- rnorm(100)
#' bins <- bin(x, "equalN", nElements = 10, minAbsX = 0.5)
#' getZeroBin(bins)
#' bins <- setZeroBin(bins, 2)
#' 
#' @return 
#' For \code{getZeroBin}, the index of the level representing the zero bin. 
#' For \code{setZeroBin}, a modified factor with the zero bin set to the 
#' provided value.
#' 
#' @name getSetZeroBin
NULL

#' @export
#' @rdname getSetZeroBin
getZeroBin <- function(bins) {
    attr(bins, "bin0")
}

#' @export 
#' @rdname getSetZeroBin
setZeroBin <- function(bins, zeroBin) {
    .assertVector(x = bins, type = "factor")
    if (is.numeric(zeroBin)) {
        .assertScalar(x = zeroBin, type = "numeric",
                      validValues = seq_len(nlevels(bins)))
        attr(bins, "bin0") <- as.integer(zeroBin)
    } else if (is.character(zeroBin)) {
        .assertScalar(x = zeroBin, type = "character",
                      validValues = levels(bins))
        attr(bins, "bin0") <- match(zeroBin, levels(bins))
    } else {
        stop("'zeroBin' must be of type 'character' or 'numeric'")
    }
    bins
}
