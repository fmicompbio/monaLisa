# Coordinates for an A-polygon (internal)
letterA <- function (x.pos, y.pos, ht, wt) {
    x <- 0.1 * c(0, 4, 6, 10, 8, 6.8, 3.2, 2, 0, 3.6, 5, 6.4, 3.6)
    y <- 0.1 * c(0, 10, 10, 0, 0, 3, 3, 0, 0, 4, 7.5, 4, 4)
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    id <- c(rep(1, 9), rep(2, 4))
    list(x = x, y = y, id = id, fill = c("#4DAF4A", "white"))
}

# Coordinates for an C-polygon (internal)
letterC <- function (x.pos, y.pos, ht, wt) {
    angle1 <- seq(0.3 + pi/2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    id <- rep(1, length(x))
    list(x = x, y = y, id = id, fill = "#377EB8")
}

# Coordinates for an G-polygon (internal)
letterG <- function (x.pos, y.pos, ht, wt) {
    angle1 <- seq(0.3 + pi/2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    h1 <- max(y.l1)
    r1 <- max(x.l1)
    h1 <- 0.4
    x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
    y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
    id <- c(rep(1, length(x)), rep(2, length(x.add)))
    x <- c(rev(x), x.add)
    y <- c(rev(y), y.add)
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    list(x = x, y = y, id = id, fill = c("#FFA500", "#FFA500"))
}

# Coordinates for an T-polygon (internal)
letterT <- function (x.pos, y.pos, ht, wt) {
    x <- 0.1 * c(0, 10, 10, 6, 6, 4, 4, 0)
    y <- 0.1 * c(10, 10, 9, 9, 0, 0, 9, 9)
    x <- x.pos + wt * x
    y <- y.pos + ht * y
    id <- rep(1, 8)
    list(x = x, y = y, id = id, fill = "#E41A1C")
}

# Add coordinates for a new base polygon to the coordinaes in 'letters' (internal)
addLetter <- function (letters, which = c("A", "C", "G", "T"), x.pos, y.pos, ht, wt) {
    which <- match.arg(which)
    letter <- switch(which,
                     "A" = letterA(x.pos, y.pos, ht, wt),
                     "C" = letterC(x.pos, y.pos, ht, wt),
                     "G" = letterG(x.pos, y.pos, ht, wt),
                     "T" = letterT(x.pos, y.pos, ht, wt))
    letters$x <- c(letters$x, letter$x)
    letters$y <- c(letters$y, letter$y)
    lastID <- ifelse(is.null(letters$id), 0, max(letters$id))
    letters$id <- c(letters$id, lastID + letter$id)
    letters$fill <- c(letters$fill, letter$fill)
    letters
}

# Calculate the information content for each position in a PFMatrix (internal)
pfm2ic <- function (pfm) {
    npos <- ncol(pfm)
    ic <- numeric(length = npos)
    for (i in 1:npos)
        ic[i] <- 2 + sum(sapply(pfm[, i], function(x) if (x > 0) (x * log2(x)) else 0))
    ic
}

#' @title Create a simple sequence logo grob.
#'
#' @description Create a simple sequence logo grob (grid-graphics object) for a
#'     transcription factor from a position frequency matrix. The logo drawing
#'     code is a simplified version from \code{\link[seqLogo]{seqLogo}} and for
#'     example can be used to embedd sequence logos within other plots.
#'
#' @param x A \code{\link[TFBSTools]{PFMatrix}} object
#'
#' @return A polygon grob.
#'
#' @example
#' if (require(JASPAR2018)) {
#'     pfm1 <- getMatrixByID(JASPAR2018, "MA0139")
#'     pfm2 <- getMatrixByID(JASPAR2018, "MA0531")
#'
#'     g1 <- seqLogoGrob(pfm1)
#'     g2 <- seqLogoGrob(pfm2)
#'
#'     gridExtra::grid.arrange(g1, g2)
#' }
#' @seealso \code{\link[seqLogo]{seqLogo}} for the original, more flexible version
#'     of this function.
#'
#' @export
seqLogoGrob <- function(x) {
    stopifnot(is(x, "PFMatrix"))

    xm <- TFBSTools::Matrix(x)
    xm <- sweep(xm, MARGIN = 2, colSums(xm), "/")
    xm[is.nan(xm)] <- 0.25
    chars <- c("A", "C", "G", "T")
    letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
    npos <- ncol(xm)
    facs <- pfm2ic(xm)
    wt <- 1
    x.pos <- 0
    for (j in 1:npos) {
        hts <- 0.95 * xm[, j] * facs[j]
        letterOrder <- order(hts)
        y.pos <- 0
        for (i in 1:4) {
            letter <- chars[letterOrder[i]]
            ht <- hts[letterOrder[i]]
            if (ht > 0)
                letters <- addLetter(letters, letter, x.pos, y.pos, ht, wt)
            y.pos <- y.pos + ht + 0.01
        }
        x.pos <- x.pos + wt
    }

    x <- unit(letters$x / max(letters$x), "npc")
    y <- unit(letters$y / max(letters$y), "npc")
    polygonGrob(x = x, y = y, id = letters$id, gp = gpar(fill = letters$fill, col = "transparent"))
}

# draw a seqLogoGrob (fill the current graphics device, internal)
drawSeqLogoGrob <- function(g, newpage = TRUE) {
    if (newpage)
        grid.newpage()
    pushViewport(plotViewport(c(0, 0, 0, 0)))
    pushViewport(dataViewport(0:1, 0:1, name = "vp1"))
    grid.draw(g)
    popViewport(2)
}

