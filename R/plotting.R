#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis hist lines par plot rect rug segments barplot matplot abline legend
#' @importFrom stats density dist hclust
NULL

#' @title Get colors by bin.
#'
#' @description Get colors for elements according to their bin.
#'
#' @param b A factor that groups elements into bins (typically the output of
#'     \code{\link{bin}}).
#' @param col1 First color.
#' @param col2 Second color.
#' @param col0 Neutral color.
#'
#' @description Colors are assigned to bins forming a gradient from \code{col1}
#'     to \code{col2} in the order of \code{levels{b}}. \code{col0} is assigned
#'     to the neutral bin (attribute \code{""}) if available.
#'
#' @seealso \code{\link{bin}}.
#'
#' @return A character vector with colors for the elements in \code{b}.
#'
#' @export
getColsByBin <- function(b,
                         col1 = c("#3F007D","#54278F","#6A51A3","#807DBA","#9E9AC8","#BCBDDC","#DADAEB"),
                         col2 = c("#FDD0A2","#FDAE6B","#FD8D3C","#F16913","#D94801","#A63603","#7F2704"),
                         col0 = "AAAAAA33") {
    if (!is.factor(b))
        b <- factor(b, levels=unique(b))

    if (!is.na(attr(b, "bin0"))) {
        bin0 <- attr(b, "bin0")
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
#' @param b A factor that groups elements of \code{x} into bins (typically the output of
#'     \code{\link{bin}}).
#' @param breaks Controls the histogram breaks (passed to \code{hist(...)}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legend If not \code{NULL}, draw a legend with binning information (will
#'     be passed to \code{legend(x=legend)} to control legend position).
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}, \code{\link[graphics]{hist}}
#'
#' @return Invisibly the return value of \code{hist(...)} that generated the plot.
#'
#' @export
plotBinHist <- function(x, b, breaks = 10 * nlevels(b),
                        xlab = deparse(substitute(x)), ylab = "Frequency",
                        main = "", legend = "topright", ...) {
    stopifnot(length(x) == length(b))
    if (!is.factor(b))
        b <- factor(b, levels=unique(b))
    cols <- getColsByBin(b, ...)
    binbreaks <- attr(b, "breaks")
    bincols <- attr(cols, "cols")
    h <- hist(x, breaks = breaks, plot = FALSE)
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    ret <- hist(x, breaks = breaks, col = bincols[findInterval(h$mids, binbreaks, all.inside = TRUE)],
                xlab = xlab, ylab = ylab, main = main)
    pusr <- par('usr'); segments(x0=pusr[c(1,1)], y0=pusr[c(4,3)], x1=pusr[c(1,2)], y1=pusr[c(3,3)])
    rug(binbreaks, col="black")
    if (!is.null(legend) && legend[1] != FALSE)
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)), fill = bincols, bty = "n")
    invisible(ret)
}

#' @title Denstity plot of binned elements.
#'
#' @description Plot the density of binned elements with binning information.
#'
#' @param x A numerical vector with the values used for binning.
#' @param b A factor that groups elements of \code{x} into bins (typically the output of
#'     \code{\link{bin}}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legend If not \code{NULL}, draw a legend with binning information (will
#'     be passed to \code{legend(x=legend)} to control legend position).
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}
#'
#' @return Invisibly the return value of \code{density(x)} that generated the plot.
#'
#' @export
plotBinDensity <- function(x, b,
                           xlab = deparse(substitute(x)), ylab = "Density",
                           main = "", legend = "topright", ...) {
    stopifnot(length(x) == length(b))
    if (!is.factor(b))
        b <- factor(b, levels=unique(b))
    cols <- getColsByBin(b, ...)
    binbreaks <- attr(b, "breaks")
    bincols <- attr(cols, "cols")
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    ret <- density(x)
    plot(ret$x, ret$y, type = "l", col = "black", xlab = xlab, ylab = ylab, main = main, axes = FALSE)
    axis(1)
    axis(2)
    pusr <- par('usr')
    segments(x0=pusr[c(1,1)], y0=pusr[c(4,3)], x1=pusr[c(1,2)], y1=pusr[c(3,3)])
    rug(binbreaks, col="black")
    dx <- diff(ret$x[1:2]) / 2
    rect(xleft = ret$x - dx, ybottom = 0, xright = ret$x + dx, ytop = ret$y,
         col = bincols[findInterval(ret$x, binbreaks, all.inside = TRUE)], border = NA)
    lines(ret$x, ret$y)

    if (!is.null(legend) && legend[1] != FALSE)
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)), fill = bincols, bty = "n")
    invisible(ret)
}

#' @title Scatter plot (xy-plot) of binned elements.
#'
#' @description Plot a scatter (xy-plot) of binned elements with binning information.
#'
#' @param x A numerical vector with x values.
#' @param y A numerical vector with y values (the values used for binning).
#' @param b A factor that groups elements of \code{x,y} into bins (typically the output
#'     of \code{\link{bin}(y)}).
#' @param cols A color vector (will be computed based on \code{b} by default using
#'     \code{\link{getColsByBin}(b)}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legend If not \code{NULL}, draw a legend with binning information (will
#'     be passed to \code{legend(x=legend)} to control legend position).
#' @param ... Further arguments passed to \code{plot(x, y, ...)}.
#'
#' @seealso \code{\link{bin}}, \code{\link{getColsByBin}}
#'
#' @return Invisibly the return value of \code{plot(x, y, ...)} that generated the plot.
#'
#' @export
plotBinScatter <- function(x, y, b,
                           cols = getColsByBin(b),
                           xlab = deparse(substitute(x)),
                           ylab = deparse(substitute(y)),
                           main = "", legend = "topright", ...) {
    stopifnot(length(x) == length(y))
    stopifnot(length(x) == length(b))
    if (length(cols) == 1L)
        cols <- rep(cols, length(x))
    stopifnot(length(x) == length(cols))
    if (!is.factor(b))
        b <- factor(b, levels=unique(b))
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    ret <- plot(x, y, pch = 16, cex = 0.6, col = cols,
                xlab = xlab, ylab = ylab, main = main, axes = FALSE, ...)
    axis(1)
    axis(2)
    pusr <- par('usr')
    segments(x0=pusr[c(1,1)], y0=pusr[c(4,3)], x1=pusr[c(1,2)], y1=pusr[c(3,3)])
    if (!is.null(legend) && legend[1] != FALSE) {
        bincols <- attr(cols, "cols")
        if (is.null(bincols))
            stop("cannot create legend automatically (missing attributes from 'b')")
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)), fill = bincols, bty = "n")
    }
    invisible(ret)
}

#' @title Heatmap of motif enrichments.
#'
#' @description Plot motif enrichments (e.g. significance or magnitude) as a heatmap.
#'
#' @param x A list with numerical matrices (motifs-by-bins), typically the return value
#'     of \code{\link{runHomer}} or \code{\link{parseHomerOutput}}.
#' @param b A factor that groups elements of \code{x,y} into bins (typically the output
#'     of \code{\link{bin}()}).
#' @param which.plots Selects which heatmaps to plot (one or several from \code{"p"}, \code{"FDR"},
#'     \code{"enr"} and \code{"log2enr"}).
#' @param width The width (in inches) of each individual heatmap, without legend.
#' @param col.enr Colors used for enrichment heatmap.
#' @param col.sig Colors used for significance hetmaps (P values and FDR).
#' @param maxEnr Cap color mapping at enrichment = \code{maxEnr} (default: 99.5th percentile).
#' @param maxSig Cap color mapping at -log10 P value or -log10 FDR = \code{maxSig}
#'     (default: 99.5th percentile).
#' @param highlight A logical vector indicating motifs to be highlighted.
#' @param cluster If \code{TRUE}, the order of transcription factors will be determined by
#'     hierarchical clustering of the \code{"enr"} component.
#'
#' @details The heatmaps are plotted side-by-side and are created internally using
#'     the \pkg{ComplexHeatmap} package.
#'
#'     Each heatmap will be \code{width} inches wide, so the total plot needs a
#'     graphics device with a width of at least \code{length(which.plots) * width}
#'     plus the space used for motif names and legend. The height will be auto-adjusted to
#'     the graphics device.
#'
#' @seealso \code{\link{bin}}, \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @references Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
#'     genomic data. Bioinformatics 2016.
#'
#' @return A list of \code{ComplexHeatmap::Heatmap} objects.
#'
#' @export
plotMotifHeatmaps <- function(x, b, which.plots = c("p", "enr", "FDR", "log2enr"), width = 4,
														 col.enr = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
														 						"#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"),
														 col.sig = c("#FFF5EB","#FEE6CE","#FDD0A2","#FDAE6B","#FD8D3C",
														 						"#F16913","#D94801","#A63603","#7F2704"),
														 maxEnr = NULL, maxSig = NULL, highlight = NULL, cluster = FALSE) {
	stopifnot(requireNamespace("ComplexHeatmap"))
	stopifnot(requireNamespace("circlize"))
	stopifnot(requireNamespace("grid"))
	stopifnot(is.list(x))
	stopifnot(all(sapply(x, ncol) == nlevels(b)))
	stopifnot(all(sapply(x, function(xx) all(dim(xx) == dim(x[[1]])))))
	if (!is.factor(b))
		b <- factor(b, levels=unique(b))
	stopifnot(all(which.plots %in% c("p", "FDR", "enr", "log2enr")))
	stopifnot(all(which.plots %in% names(x)))
	stopifnot(is.null(highlight) || (is.logical(highlight) && length(highlight) == nrow(x[[1]])))
	bincols <- attr(getColsByBin(b), "cols")
	if (cluster) {
	    clres <- stats::hclust(stats::dist(x[["enr"]]))
	    o <- clres$order
	} else {
	    o <- seq.int(nrow(x[[1]]))
	}
	hmBin <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(bin = colnames(x[[1]])), name="bin",
											   col = list(bin = bincols),
											   which = "column", width = grid::unit(width / 16,"inch"),
											   show_legend=FALSE)
	tmp <- matrix(if (!is.null(highlight)) as.character(highlight) else rep(NA, nrow(x[[1]])),
								ncol = 1, dimnames = list(rownames(x[[1]]), NULL))
	hmMotifs <- ComplexHeatmap::Heatmap(matrix = tmp, name = "names",
										width = grid::unit(if (!is.null(highlight)) .2 else 0, "inch"),
										na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
										cluster_rows = FALSE, cluster_columns = FALSE,
										show_row_names = TRUE, row_names_side = "left",
										show_column_names = FALSE, show_heatmap_legend = FALSE)

	ret <- c(list(labels = hmMotifs), lapply(which.plots, function(w) {
		dat <- x[[w]]
		if ((w == "enr") | (w == "log2enr")) {
			rng <- c(-1, 1) * if (is.null(maxEnr)) quantile(abs(dat), .995) else maxEnr
			cols <- col.enr
		} else {
			rng <- c(0, if (is.null(maxSig)) quantile(dat, .995) else maxSig)
			cols <- col.sig
		}
		hm <- ComplexHeatmap::Heatmap(matrix = dat, name = c(p="P value", FDR="FDR",
																												 enr="enrichment", log2enr="log2 enrichment")[w],
																	width = grid::unit(width,"inch"),
																	column_title = c(p = "P value (-log10)", FDR = "FDR (-log10)",
																									 enr = "enrichment (o-e)/sqrt(e)", log2enr="enrichment (log2)")[w],
																	col = circlize::colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
																														 colors = colorRampPalette(cols)(256)),
																	cluster_rows = FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE,
																	##column_names_side = "bottom", column_names_max_height = grid::unit(1.5,"inch"),
																	top_annotation = hmBin, top_annotation_height = grid::unit(width / 16, "inch"),
																	show_heatmap_legend = TRUE, heatmap_legend_param = list(color_bar="continuous"),
																	use_raster = TRUE)
		hm
	}))
	names(ret)[-1] <- which.plots
	show(Reduce(ComplexHeatmap::add_heatmap, ret))
	invisible(ret)
}


#'@title Plot Stability Paths
#'
#'@description Plot the stability paths of each variable (predictor), showing the selection probability
#'as a function of the lambda parameter.
#'
#'@param stabs_object the \code{stabs} object resulting from stability selection.
#'@param cols color vector for the varaiables from the predictor matrix. By default, it's set to NULL
#'and the function colors the variables by whether or not they were selected.
#'@param ... additional parameters to pass on to \code{matplot}
#'
#'@return plot of stability paths.
#'@export
plot_stabilityPaths <- function(stabs_object, cols=NULL, lwd = 1, lty=1, ylim=c(0,1.1), ...) {
  
  # ... checks
  if(!class(stabs_object)=="stabsel"){stop("stabs_object must be of class 'stabsel', the resulting object from running stability selection with the `stabs` package")}
  
  # set plot parameters
  mat <- t(stabs_object$phat)
  cols <- rep("black", ncol(mat))
  names(cols) <- rep("Not Selected", length(cols))
  cols[stabs_object$selected] <- "Steelblue"
  names(cols)[stabs_object$selected] <- "Selected"
  
  # plot stability paths
  matplot(mat, col=cols, type="l", lty=lty, ylab = "Selection Probability", xlab = "Regularization Step", ylim=ylim, lwd=lwd, ...)
  abline(h = stabs_object$cutoff, lty = 5, col = "red", lwd=lwd)
  legend("topleft", legend = c(unique(names(cols)), "cutoff"), col = c(unique(cols), "red"), lty=c(1,1,5), bty = "n", lwd=lwd)
  
}


#'@title Barplot Selection Probabilities
#'
#'@description Do a bar plot of the selection probabilities in order from greatest to least.
#'
#'@param stabs_object the \code{stabs} object resulting from stability selection.
#'@param ylim the limits for the y-axis.
#'@param onlySelected logical (default=TRUE) indicating if only selected predictors' selection probabilities 
#'whould be plotted.
#'@param las (2 by default) plot labels vertically or horizontally.
#'@param ... additional parameters for the \code{barplot} function.
#'
#'@return barplot of selection probabilities.
#'@export
barplot_selectionProbability <- function(stabs_object, ylim=c(0,1), onlySelected=TRUE, las=2, ...) {
  
  # ... checks
  if(!class(stabs_object)=="stabsel"){stop("stabs_object must be of class 'stabsel', the resulting object from running stability selection with the `stabs` package")}
  
  phat <- t(stabs_object$phat)
  TF_prob <- phat[nrow(phat), ]
  
  if(onlySelected) {
    TF_prob <- TF_prob[stabs$selected]
    
  }
  
  # check if empty
  if(isEmpty(TF_prob)){stop("The input for the barplot is empty")}
  
  # order 
  TF_prob <- TF_prob[order(TF_prob, decreasing = TRUE)]
  
  # plot
  barplot(TF_prob, ylim=ylim, ylab = "Selection Probability", las=las, ...)
  abline(h=stabs_object$cutoff, lty = 5, col = "red")
  legend("topright", legend = "cutoff", lty = 5, col = "red", bty = "n")
  
}


















