#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis hist lines par plot rect rug segments barplot 
#'   matplot abline legend text
#' @importFrom stats density dist hclust
#' @importFrom S4Vectors isEmpty
NULL


#' @title Get colors by bin.
#'
#' @description Get colors for elements according to their bin.
#'     Colors are assigned to bins forming a gradient from \code{col1}
#'     to \code{col2} in the order of \code{levels{b}}. \code{col0} is assigned
#'     to the neutral bin (attribute \code{""}) if available.
#'
#' @param b A factor that groups elements into bins (typically the output of
#'     \code{\link{bin}}).
#' @param col1 First color.
#' @param col2 Second color.
#' @param col0 Neutral color.
#'
#' @seealso \code{\link{bin}}.
#'
#' @return A character vector with colors for the elements in \code{b}.
#' 
#' @examples 
#' set.seed(1)
#' x <- rnorm(100)
#' b <- bin(x, "equalN", nElements = 10)
#' cols <- getColsByBin(b)
#' 
#' @export
getColsByBin <- function(b,
                         col1 = c("#003C30", "#01665E", "#35978F", 
                                  "#80CDC1", "#C7EAE5"),
                         col2 = c("#F6E8C3", "#DFC27D", "#BF812D", 
                                  "#8C510A", "#543005"),
                         col0 = "#F5F5F5") {
    if (!is.factor(b)) {
        b <- factor(b, levels = unique(b))
        b <- setZeroBin(b, NA)
    }

    if (!is.null(getZeroBin(b)) && !is.na(getZeroBin(b))) {
        bin0 <- getZeroBin(b)
        cols <- c(colorRampPalette(col1)(bin0 - 1L),
                  "#AAAAAA33",
                  colorRampPalette(col2)(nlevels(b) - bin0))
    } else {
        nh <- round(nlevels(b) / 2)
        cols <- c(colorRampPalette(col1)(nh),
                  colorRampPalette(col2)(nlevels(b) - nh))
    }

    res <- cols[b]
    names(cols) <- levels(b)
    attr(res, "cols") <- cols
    return(res)
}


#' @title Histogram of binned elements.
#'
#' @description Plot a histogram of binned elements with binning information.
#'
#' @param x A numerical vector with the values used for binning.
#' @param b A factor that groups elements of \code{x} into bins (typically 
#'     the output of \code{\link{bin}}).
#' @param breaks Controls the histogram breaks (passed to \code{hist(...)}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legend If not \code{NULL}, draw a legend with binning information 
#'     (will be passed to \code{legend(x=legend)} to control legend position).
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}, \code{\link[graphics]{hist}}
#'
#' @return Invisibly the return value of \code{hist(...)} that generated the 
#'     plot.
#' 
#' @examples 
#' set.seed(1)
#' x <- rnorm(100)
#' b <- bin(x, "equalN", nElements = 10)
#' plotBinHist(x, b)
#' 
#' @export
plotBinHist <- function(x, b, breaks = 10 * nlevels(b),
                        xlab = deparse(substitute(x, env = as.environment(-1))),
                        ylab = "Frequency",
                        main = "", legend = "topright", legend.cex = 1.0, ...) {
    .assertVector(x = b, type = "factor", len = length(x))
    stopifnot("breaks" %in% names(attributes(b)))
    .assertScalar(x = legend.cex, type = "numeric", rngExcl = c(0, Inf))
    cols <- getColsByBin(b, ...)
    binbreaks <- attr(b, "breaks")
    bincols <- attr(cols, "cols")
    h <- hist(x, breaks = breaks, plot = FALSE)
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    ret <- hist(x, breaks = breaks, 
                col = bincols[findInterval(h$mids, binbreaks, 
                                           all.inside = TRUE)],
                xlab = xlab, ylab = ylab, main = main)
    pusr <- par('usr')
    segments(x0 = pusr[c(1,1)], y0 = pusr[c(4,3)],
             x1 = pusr[c(1,2)], y1 = pusr[c(3,3)])
    rug(binbreaks, col = "black")
    if (!is.null(legend) && legend[1] != FALSE)
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)),
               fill = bincols, bty = "n", cex = legend.cex)
    invisible(ret)
}


#' @title Density plot of binned elements.
#'
#' @description Plot the density of binned elements with binning information.
#'
#' @param x A numerical vector with the values used for binning.
#' @param b A factor that groups elements of \code{x} into bins (typically the 
#'     output of \code{\link{bin}}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legend If not \code{NULL}, draw a legend with binning information
#'     (will be passed to \code{legend(x=legend)} to control legend position).
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}
#'
#' @return Invisibly the return value of \code{density(x)} that generated the 
#'     plot.
#'
#' @examples 
#' set.seed(1)
#' x <- rnorm(100)
#' b <- bin(x, "equalN", nElements = 10)
#' plotBinDensity(x, b)
#' 
#' @export
plotBinDensity <- function(x, b,
                           xlab = deparse(substitute(x, 
                                                     env = as.environment(-1))),
                           ylab = "Density",
                           main = "", legend = "topright", 
                           legend.cex = 1.0, ...) {
    .assertVector(x = b, type = "factor", len = length(x))
    stopifnot("breaks" %in% names(attributes(b)))
    .assertScalar(x = legend.cex, type = "numeric", rngExcl = c(0, Inf))
    cols <- getColsByBin(b, ...)
    binbreaks <- attr(b, "breaks")
    bincols <- attr(cols, "cols")
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    ret <- density(x)
    plot(ret$x, ret$y, type = "l", col = "black", xlab = xlab, ylab = ylab, 
         main = main, axes = FALSE)
    axis(1)
    axis(2)
    pusr <- par('usr')
    segments(x0 = pusr[c(1,1)], y0 = pusr[c(4,3)],
             x1 = pusr[c(1,2)], y1 = pusr[c(3,3)])
    rug(binbreaks, col = "black")
    dx <- diff(ret$x[seq_len(2)]) / 2
    rect(xleft = ret$x - dx, ybottom = 0, xright = ret$x + dx, ytop = ret$y,
         col = bincols[findInterval(ret$x, binbreaks, all.inside = TRUE)], 
         border = NA)
    lines(ret$x, ret$y)

    if (!is.null(legend) && legend[1] != FALSE)
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)),
               fill = bincols, bty = "n", cex = legend.cex)
    invisible(ret)
}


#' @title Scatter plot (xy-plot) of binned elements.
#'
#' @description Plot a scatter (xy-plot) of binned elements with binning 
#'     information.
#'
#' @param x A numerical vector with x values.
#' @param y A numerical vector with y values (the values used for binning).
#' @param b A factor that groups elements of \code{x,y} into bins (typically 
#'     the output of \code{\link{bin}(y)}).
#' @param cols A color vector (will be computed based on \code{b} by default 
#'     using \code{\link{getColsByBin}(b)}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legend If not \code{NULL}, draw a legend with binning information 
#'     (will be passed to \code{legend(x=legend)} to control legend position).
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
#' @param ... Further arguments passed to \code{plot(x, y, ...)}.
#'
#' @seealso \code{\link{bin}}, \code{\link{getColsByBin}}
#'
#' @return \code{TRUE} (invisibly).
#'
#' @examples 
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' b <- bin(y, "equalN", nElements = 10)
#' plotBinScatter(x, y, b)
#' 
#' @export
plotBinScatter <- function(x, y, b,
                           cols = getColsByBin(b),
                           xlab = deparse(substitute(x, 
                                                     env = as.environment(-1))),
                           ylab = deparse(substitute(y, 
                                                     env = as.environment(-1))),
                           main = "", legend = "topright", 
                           legend.cex = 1.0, ...) {
    .assertVector(x = y, len = length(x))
    .assertVector(x = b, len = length(x))
    .assertScalar(x = legend.cex, type = "numeric", rngExcl = c(0, Inf))
    if (length(cols) == 1L)
        cols <- rep(cols, length(x))
    stopifnot(length(x) == length(cols))
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    plot(x, y, pch = 16, cex = 0.6, col = cols,
         xlab = xlab, ylab = ylab, main = main, axes = FALSE, ...)
    axis(1)
    axis(2)
    pusr <- par('usr')
    segments(x0 = pusr[c(1,1)], y0 = pusr[c(4,3)],
             x1 = pusr[c(1,2)], y1 = pusr[c(3,3)])
    if (!is.null(legend) && legend[1] != FALSE) {
        stopifnot("cols" %in% names(attributes(cols)))
        bincols <- attr(cols, "cols")
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)),
               fill = bincols, bty = "n", cex = legend.cex)
    }
    invisible(TRUE)
}


#' @title Heatmap of motif enrichments.
#'
#' @description Plot motif enrichments (e.g. significance or magnitude) as a 
#'     heatmap.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with 
#'     numerical matrices (motifs-by-bins) in its \code{assays()}, typically 
#'     the return value of \code{\link{calcBinnedMotifEnrR}} or 
#'     \code{\link{calcBinnedMotifEnrHomer}}.
#' @param which.plots Selects which heatmaps to plot (one or several from 
#'     \code{"negLog10P"}, \code{"negLog10Padj"}, \code{"pearsonResid"} and 
#'     \code{"log2enr"}).
#' @param width The width (in inches) of each individual heatmap, without 
#'     legend.
#' @param col.enr Colors used for enrichment heatmap ("pearsonResid" and 
#'     "log2enr").
#' @param col.sig Colors used for significance hetmaps ("negLog10P" and 
#'     "negLog10Padj").
#' @param col.gc Colors used for motif GC content (for 
#'     \code{show_motif_GC = TRUE}).
#' @param maxEnr Cap color mapping at enrichment = \code{maxEnr}
#'     (default: 99.5th percentile).
#' @param maxSig Cap color mapping at -log10 P value or -log10 FDR = 
#'     \code{maxSig} (default: 99.5th percentile).
#' @param highlight A logical vector indicating motifs to be highlighted.
#' @param cluster If \code{TRUE}, the order of transcription factors will be 
#'     determined by hierarchical clustering of the \code{"pearsonResid"} 
#'     component. Alternatively, an \code{hclust}-object can be supplied which 
#'     will determine the motif ordering.
#'     No reordering is done for \code{cluster = FALSE}.
#' @param show_dendrogram If \code{cluster != FALSE}, controls whether to show
#'     a row dendrogram for the clustering of motifs. Ignored for 
#'     \code{cluster = FALSE}.
#' @param show_motif_GC If \code{TRUE}, show a column with the percent G+C of 
#'     the motif as part of the heatmap.
#' @param show_seqlogo If \code{TRUE}, show a sequence logo next to each motif 
#'     label. This will likely only make sense for a heatmap with a low number 
#'     of motifs.
#' @param width.seqlogo The width (in inches) for the longest sequence logo 
#'     (shorter logos are drawn to scale).
#' @param use_raster \code{TRUE} or \code{FALSE} (default). Passed to 
#'     \code{use_raster} of \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param na_col "white" (default). Passed to \code{na_col} of 
#'     \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param doPlot If \code{TRUE} (default), plot the generated heatmap(s)
#'     using \code{Reduce(ComplexHeatmap::add_heatmap, heatmapList)}. If
#'     \code{FALSE}, just return the list of heatmap(s) (\code{heatmapList}) in
#'     example before), allowing to modify them further before plotting.
#' @param ... Further arguments passed to \code{\link[ComplexHeatmap]{Heatmap}}
#'     when creating the main heatmaps selected by \code{which.plots}. 
#'
#' @details The heatmaps are created using the \pkg{ComplexHeatmap} package
#'     and plotted side-by-side.
#'
#'     Each heatmap will be \code{width} inches wide, so the total plot needs a
#'     graphics device with a width of at least 
#'     \code{length(which.plots) * width} plus the space used for motif names 
#'     and legend. The height will be auto-adjusted to the graphics device.
#'
#' @seealso \code{\link{bin}}, \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @references Gu, Z. Complex heatmaps reveal patterns and correlations in 
#'     multidimensional genomic data. Bioinformatics 2016.
#'
#' @return A list of \code{ComplexHeatmap::Heatmap} objects.
#'
#' @examples 
#' se <- readRDS(system.file("extdata", 
#'                           "results.binned_motif_enrichment_LMRs.rds", 
#'                           package = "monaLisa"))
#' i <- which(SummarizedExperiment::assay(se, "negLog10Padj")[, 8] > 4)
#' plotMotifHeatmaps(se[i, ], which.plots = "pearsonResid",
#'                   width = 2, show_seqlogo = TRUE)
#' 
#' @importFrom methods is show
#' @importFrom stats hclust dist quantile
#' @importFrom TFBSTools Matrix
#' @importFrom grDevices colorRampPalette
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assayNames assay rowData
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap add_heatmap
#' @importFrom grid unit
#' @importFrom circlize colorRamp2
#'
#' @export
plotMotifHeatmaps <- function(x,
                              which.plots = c("negLog10P", "pearsonResid", 
                                              "negLog10Padj", "log2enr"),
                              width = 4,
                              col.enr = c("#053061", "#2166AC", "#4393C3",
                                          "#92C5DE", "#D1E5F0", "#F7F7F7",
                                          "#FDDBC7", "#F4A582", "#D6604D",
                                          "#B2182B", "#67001F"),
                              col.sig = c("#F0F0F0", "#D9D9D9", "#BDBDBD",
                                          "#969696", "#737373", "#525252",
                                          "#252525", "#000000"),
                              col.gc = c("#F7FCF5", "#E5F5E0", "#C7E9C0",
                                         "#A1D99B", "#74C476", "#41AB5D",
                                         "#238B45", "#006D2C", "#00441B"),
                              maxEnr = NULL,
                              maxSig = NULL,
                              highlight = NULL,
                              cluster = FALSE,
                              show_dendrogram = FALSE,
                              show_motif_GC = FALSE,
                              show_seqlogo = FALSE,
                              width.seqlogo = 1.5,
                              use_raster = FALSE,
                              na_col = "white", 
                              doPlot = TRUE,
                              ...) {
    stopifnot(exprs = {
        is(x, "SummarizedExperiment")
        all(which.plots %in% assayNames(x))
        "bins" %in% names(metadata(x))
        (!show_motif_GC || "motif.percentGC" %in% colnames(rowData(x)))
    })
    b <- metadata(x)$bins
    .assertScalar(x = width, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = show_dendrogram, type = "logical")
    .assertScalar(x = show_motif_GC, type = "logical")
    .assertScalar(x = show_seqlogo, type = "logical")
    .assertScalar(x = width.seqlogo, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = use_raster, type = "logical")
    .assertScalar(x = na_col, type = "character")
    .assertScalar(x = doPlot, type = "logical")
    stopifnot(exprs = {
        ncol(x) == nlevels(b)
        all(which.plots %in% c("negLog10P", "negLog10Padj", 
                               "pearsonResid", "log2enr"))
        is.null(highlight) || (is.logical(highlight) && 
                                 length(highlight) == nrow(x))
    })
    bincols <- attr(getColsByBin(b), "cols")
    if (identical(cluster, TRUE)) {
        clAssayName <- "pearsonResid"
        clAssay <- assay(x, clAssayName)
        allNA <- rowSums(is.na(clAssay)) == ncol(clAssay)
        if (any(allNA)) {
            warning("removing motifs without finite values in '",
                    clAssayName, "': ",
                    paste(rownames(clAssay)[allNA], collapse = ", "))
            x <- x[!allNA, ]
            clAssay <- clAssay[!allNA, ]
        }
        clres <- hclust(dist(clAssay))
    } else if (identical(cluster, FALSE)) {
        clres <- FALSE
    } else if (is(cluster, "hclust")) {
        clres <- cluster
    } else {
        stop("'cluster' must be either TRUE, FALSE or an hclust-object.")
    }
    hmBin <- HeatmapAnnotation(df = data.frame(bin = colnames(x)), name = "bin",
                               col = list(bin = bincols),
                               show_annotation_name = FALSE,
                               which = "column", width = unit(width,"inch"),
                               annotation_height = unit(width / 16, "inch"),
                               show_legend = FALSE)
    tmp <- matrix(if (!is.null(highlight)) {
        as.character(highlight) 
    } else {
        rep(NA, nrow(x))
    },
    ncol = 1, dimnames = list(unname(rowData(x)$motif.name), NULL))
    hmSeqlogo <- NULL
    if (show_seqlogo) {
        pfms <- rowData(x)$motif.pfm
        maxwidth <- max(vapply(TFBSTools::Matrix(pfms), ncol, 0L))
        grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
        hmSeqlogo <- HeatmapAnnotation(
            logo = annoSeqlogo(grobL = grobL, which = "row",
                               space = unit(0.5, "mm"),
                               width = unit(width.seqlogo, "inch")),
            show_legend = FALSE, show_annotation_name = FALSE, which = "row")
    }
    hmMotifs <- Heatmap(
        matrix = tmp, name = "names",
        width = unit(if (!is.null(highlight)) .2 else 0, "inch"),
        na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
        cluster_rows = clres, show_row_dend = show_dendrogram,
        cluster_columns = FALSE, show_row_names = TRUE,
        row_names_side = "left", show_column_names = FALSE,
        show_heatmap_legend = FALSE, left_annotation = hmSeqlogo
    )

    assayNameMap1 <- c(negLog10P = "P value",
                       negLog10Padj = "adj. P value",
                       pearsonResid = "Pearson residual",
                       log2enr = "log2 enrichment")
    assayNameMap2 <- c(negLog10P = "P value (-log10)",
                       negLog10Padj = "adj. P value (-log10)",
                       pearsonResid = "Pearson residual (o-e)/sqrt(e)",
                       log2enr = "enrichment (log2)")
    L <- list(labels = hmMotifs)
    if (show_motif_GC) {
        tmp <- as.matrix(rowData(x)[, "motif.percentGC", drop = FALSE])
        hmPercentGC <- Heatmap(
            matrix = tmp, name = "Percent G+C",
            width = unit(0.2, "inch"), na_col = NA,
            col = colorRamp2(breaks = c(0, seq(20, 80, length.out = 254), 100),
                             colors = colorRampPalette(col.gc)(256)),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(color_bar = "continuous"),
            use_raster = use_raster
        )
        L <- c(L, list("percentGC" = hmPercentGC))
    }
    ret <- c(L, lapply(which.plots, function(w) {
        dat <- assay(x, w)
        if ((w == "pearsonResid") | (w == "log2enr")) {
            rng <- c(-1, 1) * if (is.null(maxEnr)) {
                quantile(abs(dat), .995, na.rm = TRUE) 
            } else {
                maxEnr
            }
            cols <- col.enr
        } else {
            rng <- c(0, 
                     if (is.null(maxSig)) {
                         quantile(dat, .995, na.rm = TRUE) 
                     } else {
                         maxSig
                     })
            cols <- col.sig
        }
        Heatmap(
            matrix = dat,
            name = assayNameMap1[w],
            width = unit(width,"inch"),
            column_title = assayNameMap2[w],
            col = colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
                             colors = colorRampPalette(cols)(256)),
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            # column_names_side = "bottom", 
            # column_names_max_height = unit(1.5,"inch"),
            top_annotation = hmBin, show_heatmap_legend = TRUE,
            heatmap_legend_param = list(color_bar = "continuous"),
            use_raster = use_raster,
            na_col = na_col, 
            ...
        )
    }))
    names(ret)[seq(length(ret) - length(which.plots) + 1L, length(ret))] <- 
        which.plots
    if (doPlot) {
        show(Reduce(ComplexHeatmap::add_heatmap, ret))
    }
    invisible(ret)
}


#' @title Plot Stability Paths
#'
#' @description Plot the stability paths of each variable (predictor), 
#'   showing the selection probability as a function of the regularization step.
#'
#' @param se the \code{SummarizedExperiment} object resulting from stability 
#'   selection, by running \code{\link[monaLisa]{randLassoStabSel}}.
#' @param selProbMin A numerical scalar in [0,1]. Predictors with a selection
#'   probability greater than \code{selProbMin} are shown as colored lines. The
#'   color is defined by the \code{col} argument.
#' @param col color of the selected predictors.
#' @param lwd line width (default = 1).
#' @param lty line type (default = 1).
#' @param ylim limits for y-axis (default = c(0,1.1)).
#' @param ... additional parameters to pass on to \code{matplot}.
#'
#' @return \code{TRUE} (invisibly).
#' 
#' @examples 
#' ## create data set
#' Y <- rnorm(n = 500, mean = 2, sd = 1)
#' X <- matrix(data = NA, nrow = length(Y), ncol = 50)
#' for (i in seq_len(ncol(X))) {
#'   X[ ,i] <- runif(n = 500, min = 0, max = 3)
#' }
#' s_cols <- sample(x = seq_len(ncol(X)), size = 10, 
#'   replace = FALSE)
#' for (i in seq_along(s_cols)) {
#'   X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
#' }
#'   
#' ## reproducible randLassoStabSel() with 1 core
#' set.seed(123)
#' ss <- randLassoStabSel(x = X, y = Y)
#' plotStabilityPaths(ss)
#'
#' @seealso \code{\link[stabs]{stabsel}} and \code{\link[graphics]{matplot}}
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom graphics matplot
#'
#' @export
plotStabilityPaths <- function(se,
                               selProbMin = metadata(se)$stabsel.params.cutoff,
                               col = "cadetblue", 
                               lwd = 1, lty = 1, ylim = c(0, 1.1), ...) {
    # checks
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment")
    }

    # set plot parameters
    mat <- as.matrix(colData(se))
    mat <- t(mat[, grep(pattern = "^regStep", x = colnames(mat))])
    cols <- rep("black", ncol(mat))
    sel <- se$selProb > selProbMin
    cols[sel] <- col
  
    # plot stability paths
    graphics::matplot(mat, col = cols, type = "l", lty = lty,
            ylab = "Selection Probability", xlab = "Regularization Step",
            ylim = ylim, lwd = lwd, ...)
    abline(h = selProbMin, lty = 5, col = "red", lwd = lwd)
    legend("topleft", legend = c("not selected", "selected", "selProbMin"),
           col = c("black", col, "red"), lty = c(1, 1, 5), bty = "n", lwd = lwd)
    
    # return TRUE
    invisible(TRUE)
}


#' @title Plot selection probabilities of predictors
#'
#' @description This function plots the selection probabilities of predictors
#'   (for example the selected motifs), optionally multiplied with either +1 or
#'   -1 to give a sense of both the strength and the directionality of the
#'   associated effects. The directionality is estimated from the sign of the
#'   correlation coefficient between each predictor and the response vector.
#'
#' @param se The \code{SummarizedExperiment} object with the results from
#'   stability selection (typically returned by \code{\link{randLassoStabSel}}).
#' @param directional A logical scalar. If \code{TRUE}, selection probabilities
#'   are plotted with the sign of the marginal correlation between a predictor
#'   and the response.
#' @param selProbMin A numerical scalar in [0,1]. Predictors with a selection
#'   probability greater than \code{selProbMin} are shown as colored bars. The
#'   color is defined by \code{col[1]}. By default, \code{selProbMin} is
#'   extracted from the parameters stored in \code{se}.
#' @param selProbMinPlot A numerical scalar in [0,1] less than 
#'   \code{selProbMin}.
#'   Predictors with a selection probability greater than \code{selProbMinPlot}
#'   but less than \code{selProbMin} are shown as bars with color \code{col[2]}.
#'   \code{selProbMinPlot} is useful to include additional predictors in the 
#'   plot that were not selected according to \code{selProbMin} but may be 
#'   close to that cutoff. Setting \code{selProbMinPlot = 0} will create a plot 
#'   including all predictors.
#' @param showSelProbMin A logical scalar. If \code{TRUE}, the value of
#'   \code{selProbMin} is shown by a horizontal dashed line of color 
#'   \code{col[3]}.
#' @param col A color vector giving the three colors used for predictors with
#'   selection probability greater than \code{selProbMin}, additional predictors
#'   with selection probability greater than \code{selProbMinPlot}, and the
#'   selection probability cutoff line.
#' @param method A character scalar with the correlation method to use in the
#'   calculation of predictor-response marginal correlations. One of "pearson",
#'   "kendall" or "spearman" (see \code{\link[stats]{cor}}).
#' @param ylimext A numeric scalar defining how much the y axis limits should be
#'   expanded beyond the plotted probabilities to allow for space for the
#'   bar labels.
#' @param legend the position of the legend in the bar plot (will
#'     be passed to \code{legend(x=legend)} to control legend position).
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
#' @param ... additional parameters passed to \code{\link[graphics]{barplot}}. 
#'
#' @details This function creates a bar plot using the 
#'   \code{\link[graphics]{barplot}} function.
#'   Each bar corresponds to a predictor (motif) and the colors correspond to 
#'   whether or not it was selected. The y-axis shows the selection 
#'   probabilities (\code{directional=FALSE}) or selection probabilities with 
#'   the sign of the marginal correlation to the response 
#'   (\code{directional=TRUE}). 
#'
#' @return a \code{matrix} with one column, containing the coordinates of the 
#'   bar midpoints, or \code{NULL} if no bar plot is drawn. 
#'
#' @examples 
#' ## create data set
#' Y <- rnorm(n = 500, mean = 2, sd = 1)
#' X <- matrix(data = NA, nrow = length(Y), ncol = 50)
#' for (i in seq_len(ncol(X))) {
#'   X[ ,i] <- runif(n = 500, min = 0, max = 3)
#' }
#' s_cols <- sample(x = seq_len(ncol(X)), size = 10, 
#'   replace = FALSE)
#' for (i in seq_along(s_cols)) {
#'   X[ ,s_cols[i]] <- X[ ,s_cols[i]] + Y
#' }
#'   
#' ## reproducible randLassoStabSel() with 1 core
#' set.seed(123)
#' ss <- randLassoStabSel(x = X, y = Y)
#' plotSelectionProb(ss)
#'
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom S4Vectors metadata
#' @importFrom stats cor
#' @importFrom graphics barplot abline legend text axis
#'
#' @export
plotSelectionProb <- function(se,
                              directional = TRUE,
                              selProbMin = metadata(se)$stabsel.params.cutoff, 
                              selProbMinPlot = 0.4,
                              showSelProbMin = TRUE,
                              col = c("cadetblue", "grey", "red"),
                              method = c("pearson", "kendall", "spearman"),
                              ylimext = 0.25,
                              legend = "topright", 
                              legend.cex = 1.0, 
                              ...) {
    
    # checks
    .assertScalar(x = directional, type = "logical")
    .assertScalar(x = selProbMin, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = selProbMinPlot, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = showSelProbMin, type = "logical")
    .assertScalar(x = legend, type = "character")
    .assertScalar(x = legend.cex, type = "numeric", rngExcl = c(0, Inf))
    stopifnot(exprs = {
        is(se, "SummarizedExperiment")
        selProbMin >= selProbMinPlot
    })
    .assertVector(x = col, len = 3L)
    method <- match.arg(method)
    .assertScalar(x = ylimext, type = "numeric", rngIncl = c(0, Inf))

    # selection probabilities * sign(correlation to y)
    probs <- se$selProb
    cols <- ifelse(probs > selProbMin, col[1], col[2])
    if (directional) {
        corcoef <- as.vector(cor(x = SummarizedExperiment::rowData(se)$y,
                                 y = SummarizedExperiment::assay(se, "x"),
                                 method = method))
        probs <- probs * sign(corcoef)
    }

    # kept and ordered
    keep <- which(abs(probs) >= selProbMinPlot)
    keep <- keep[order(probs[keep], decreasing = TRUE)]
    cols <- cols[keep]
    predNames <- colnames(se)[keep]
    probs <- probs[keep]
    up <- probs > 0

    # plot
    if (any(keep)) {
        ret <- graphics::barplot(probs, col = cols, border = NA,
                                        ylab = ifelse(
                                            directional, 
                                            "Directional selection probability",
                                            "Selection probability"
                                        ),
                                        names.arg = NA, axes = FALSE,
                                        ylim = c(min(probs) - ylimext,
                                                 max(probs) + ylimext),
                                        ...)
        ys <- pretty(x = c(0, probs))
        graphics::axis(side = 2, at = ys)
        if (showSelProbMin) {
            hval <- if (directional) c(-1, 1) * selProbMin else selProbMin
            graphics::abline(h = hval, lty = 5, col = col[3])
        }
        graphics::legend(x = legend, bty = "n", fill = col[seq_len(2)], 
                         border = NA, legend = c("selected", "not selected"), 
                         cex = legend.cex)
        if (any(up)) {
            graphics::text(x = ret[up], y = probs[up] + par("cxy")[2] / 3,
                           labels = predNames[up], col = cols[up],
                           xpd = TRUE, srt = 90, adj = c(0, 0.5))
        }
        if (any(!up)) {
            graphics::text(x = ret[!up], y = probs[!up] - par("cxy")[2] / 3,
                           labels = predNames[!up], col = cols[!up],
                           xpd = TRUE, srt = 90, adj = c(1, 0.5))
        }
    } else{
        ret <- NULL
    }
    invisible(ret)
}
