#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis hist lines par plot rect rug segments barplot matplot abline legend text
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
#' @export
getColsByBin <- function(b,
                         col1 = c("#3F007D","#54278F","#6A51A3","#807DBA","#9E9AC8","#BCBDDC","#DADAEB"),
                         col2 = c("#FDD0A2","#FDAE6B","#FD8D3C","#F16913","#D94801","#A63603","#7F2704"),
                         col0 = "AAAAAA33") {
    if (!is.factor(b)) {
        b <- factor(b, levels=unique(b))
        attr(b, "bin0") <- NA
    }

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
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}, \code{\link[graphics]{hist}}
#'
#' @return Invisibly the return value of \code{hist(...)} that generated the plot.
#'
#' @export
plotBinHist <- function(x, b, breaks = 10 * nlevels(b),
                        xlab = deparse(substitute(x)), ylab = "Frequency",
                        main = "", legend = "topright", legend.cex = 1.0, ...) {
    stopifnot(length(x) == length(b))
    stopifnot(exprs = { is.factor(b); "breaks" %in% names(attributes(b)) })
    stopifnot(exprs = { is.numeric(legend.cex); length(legend.cex) == 1 })
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
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)),
               fill = bincols, bty = "n", cex = legend.cex)
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
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}
#'
#' @return Invisibly the return value of \code{density(x)} that generated the plot.
#'
#' @export
plotBinDensity <- function(x, b,
                           xlab = deparse(substitute(x)), ylab = "Density",
                           main = "", legend = "topright", legend.cex = 1.0, ...) {
    stopifnot(length(x) == length(b))
    stopifnot(exprs = { is.factor(b); "breaks" %in% names(attributes(b)) })
    stopifnot(exprs = { is.numeric(legend.cex); length(legend.cex) == 1 })
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
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)),
               fill = bincols, bty = "n", cex = legend.cex)
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
#' @param legend.cex A scalar that controls the text size in the legend relative
#'     to the current \code{par("cex")} (see \code{\link{legend}}).
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
                           main = "", legend = "topright", legend.cex = 1.0, ...) {
    stopifnot(length(x) == length(y))
    stopifnot(length(x) == length(b))
    stopifnot(exprs = { is.numeric(legend.cex); length(legend.cex) == 1 })
    if (length(cols) == 1L)
        cols <- rep(cols, length(x))
    stopifnot(length(x) == length(cols))
    par(mar = c(5, 4, 4 - if (main == "") 3 else 0, 2) + 0.1, cex = 1.25)
    ret <- plot(x, y, pch = 16, cex = 0.6, col = cols,
                xlab = xlab, ylab = ylab, main = main, axes = FALSE, ...)
    axis(1)
    axis(2)
    pusr <- par('usr')
    segments(x0=pusr[c(1,1)], y0=pusr[c(4,3)], x1=pusr[c(1,2)], y1=pusr[c(3,3)])
    if (!is.null(legend) && legend[1] != FALSE) {
        stopifnot("cols" %in% names(attributes(cols)))
        bincols <- attr(cols, "cols")
        legend(x = legend, legend = sprintf("%s : %d", levels(b), table(b)),
               fill = bincols, bty = "n", cex = legend.cex)
    }
    invisible(ret)
}

#' @title Heatmap of motif enrichments.
#'
#' @description Plot motif enrichments (e.g. significance or magnitude) as a heatmap.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with numerical matrices
#'     (motifs-by-bins) in its \code{assays()}, typically the return value
#'     of \code{\link{runHomer}}.
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
#'     hierarchical clustering of the \code{"enr"} component. Alternatively, an
#'     \code{hclust}-object can be supplied which will determine the motif ordering.
#'     No reordering is done for \code{cluster = FALSE}.
#' @param show_dendrogram If \code{cluster != FALSE}, controls whether to show
#'     a row dendrogram for the clustering of motifs. Ignored for \code{cluster = FALSE}.
#' @param show_motif_GC If \code{TRUE}, show a column with the percent G+C of the motif
#'     as part of the heatmap.
#' @param show_seqlogo If \code{TRUE}, show a sequence logo next to each motif label.
#'     This will likely only make sense for a heatmap with a low number of motifs.
#' @param width.seqlogo The width (in inches) for the longest sequence logo (shorter
#'     logos are drawn to scale).
#' @param use_raster \code{TRUE} or \code{FALSE} (default). Passed to \code{use_raster}
#'     of \code{\link[ComplexHeatmap]{Heatmap}}.
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
#' @importFrom methods is
#' @importFrom grDevices colorRampPalette
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assayNames assay rowData
#'
#' @export
plotMotifHeatmaps <- function(x, which.plots = c("p", "enr", "FDR", "log2enr"), width = 4,
                              col.enr = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
                                          "#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"),
                              col.sig = c("#FFF5EB","#FEE6CE","#FDD0A2","#FDAE6B","#FD8D3C",
                                          "#F16913","#D94801","#A63603","#7F2704"),
                              maxEnr = NULL, maxSig = NULL, highlight = NULL, cluster = FALSE,
                              show_dendrogram = FALSE, show_motif_GC = FALSE,
                              show_seqlogo = FALSE, width.seqlogo = 1.5,
                              use_raster = FALSE) {
	stopifnot(exprs = {
	    is(x, "SummarizedExperiment")
	    all(assayNames(x) == c("p", "FDR", "enr", "log2enr"))
	    "bins" %in% names(metadata(x))
	    (!show_motif_GC || "motif.percentGC" %in% colnames(rowData(x)))
	})
	b <- metadata(x)$bins
	stopifnot(exprs = {
	    ncol(x) == nlevels(b)
	    all(which.plots %in% c("p", "FDR", "enr", "log2enr"))
	    is.numeric(width)
	    length(width) == 1
	    width > 0
	    is.null(highlight) || (is.logical(highlight) && length(highlight) == nrow(x))
	    is.logical(show_dendrogram)
	    length(show_dendrogram) == 1L
	    is.logical(show_motif_GC)
	    length(show_motif_GC) == 1L
	    is.logical(show_seqlogo)
	    length(show_seqlogo) == 1L
	    is.numeric(width.seqlogo)
	    length(width.seqlogo) == 1
	    width.seqlogo > 0
	    is.logical(use_raster)
	    length(use_raster) == 1
	})
	bincols <- attr(getColsByBin(b), "cols")
	if (is.logical(cluster) && length(cluster) == 1 && cluster[1] == TRUE) {
	    clres <- stats::hclust(stats::dist(SummarizedExperiment::assay(x, "enr")))
	} else if (is.logical(cluster) && length(cluster) == 1 && cluster[1] == FALSE) {
	    clres <- FALSE
	} else if (is(cluster, "hclust")) {
	    clres <- cluster
	} else {
	    stop("'cluster' must be either TRUE, FALSE or an hclust-object.")
	}
	hmBin <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(bin = colnames(x)), name="bin",
											   col = list(bin = bincols),
											   show_annotation_name = FALSE,
											   which = "column", width = grid::unit(width,"inch"),
											   annotation_height = grid::unit(width / 16, "inch"),
											   show_legend = FALSE)
	tmp <- matrix(if (!is.null(highlight)) as.character(highlight) else rep(NA, nrow(x)),
								ncol = 1, dimnames = list(rownames(x), NULL))
	hmSeqlogo <- NULL
	if (show_seqlogo) {
	    pfms <- rowData(x)$motif.pfm
	    maxwidth <- max(sapply(TFBSTools::Matrix(pfms), ncol))
	    grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
	    hmSeqlogo <- ComplexHeatmap::HeatmapAnnotation(
	        logo = anno_seqlogo(grobL = grobL, which = "row",
	                            space = grid::unit(0.5, "mm"), width = grid::unit(width.seqlogo, "inch")),
	        show_legend = FALSE, show_annotation_name = FALSE, which = "row")
	}
	hmMotifs <- ComplexHeatmap::Heatmap(matrix = tmp, name = "names",
										width = grid::unit(if (!is.null(highlight)) .2 else 0, "inch"),
										na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
										cluster_rows = clres, show_row_dend = show_dendrogram, cluster_columns = FALSE,
										show_row_names = TRUE, row_names_side = "left",
										show_column_names = FALSE, show_heatmap_legend = FALSE,
										left_annotation = hmSeqlogo)

	assayNameMap1 <- c(p="P value", FDR="FDR", enr="enrichment", log2enr="log2 enrichment")
	assayNameMap2 <- c(p = "P value (-log10)", FDR = "FDR (-log10)",
	                   enr = "enrichment (o-e)/sqrt(e)", log2enr="enrichment (log2)")
	L <- list(labels = hmMotifs)
	if (show_motif_GC) {
	    tmp <- as.matrix(SummarizedExperiment::rowData(x)[, "motif.percentGC", drop = FALSE])
	    gccols <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
	    hmPercentGC <- ComplexHeatmap::Heatmap(matrix = tmp, name = "Percent G+C",
	                                           width = grid::unit(0.2, "inch"), na_col = NA,
	                                           col = circlize::colorRamp2(breaks = c(0, seq(20, 80, length.out = 254), 100),
	                                                                      colors = colorRampPalette(gccols)(256)),
	                                           cluster_rows = FALSE, cluster_columns = FALSE,
	                                           show_row_names = FALSE, show_column_names = FALSE,
	                                           show_heatmap_legend = TRUE, heatmap_legend_param = list(color_bar="continuous"),
	                                           use_raster = use_raster)
	    L <- c(L, list("percentGC" = hmPercentGC))
	}
	ret <- c(L, lapply(which.plots, function(w) {
		dat <- SummarizedExperiment::assay(x, w)
		if ((w == "enr") | (w == "log2enr")) {
			rng <- c(-1, 1) * if (is.null(maxEnr)) quantile(abs(dat), .995) else maxEnr
			cols <- col.enr
		} else {
			rng <- c(0, if (is.null(maxSig)) quantile(dat, .995) else maxSig)
			cols <- col.sig
		}
		hm <- ComplexHeatmap::Heatmap(matrix = dat,
		                              name = assayNameMap1[w],
		                              width = grid::unit(width,"inch"),
		                              column_title = assayNameMap2[w],
		                              col = circlize::colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
		                                                         colors = colorRampPalette(cols)(256)),
		                              cluster_rows = FALSE, cluster_columns=FALSE, show_row_names=FALSE, show_column_names=FALSE,
		                              ##column_names_side = "bottom", column_names_max_height = grid::unit(1.5,"inch"),
		                              top_annotation = hmBin,
		                              show_heatmap_legend = TRUE, heatmap_legend_param = list(color_bar="continuous"),
		                              use_raster = use_raster)
		hm
	}))
	names(ret)[seq(length(ret) - length(which.plots) + 1L, length(ret))] <- which.plots
	show(Reduce(ComplexHeatmap::add_heatmap, ret))
	invisible(ret)
}


#'@title Plot Stability Paths
#'
#'@description Plot the stability paths of each variable (predictor), showing the selection probability
#'as a function of the regularization step.
#'
#'@param stabs_object the \code{stabs} object resulting from stability selection.
#'@param cols color vector for the varaiables from the predictor matrix. By default, it's set to NULL
#'and the function colors the variables by whether or not they were selected.
#'@param lwd line width (default 1).
#'@param lty line type (default 1).
#'@param ylim limits for y-axis (default c(0,1.1)).
#'@param ... additional parameters to pass on to \code{matplot}.
#'
#'@return plot of stability paths.
#'
#'@seealso \code{\link[stabs]{stabsel}} and \code{\link[graphics]{matplot}}
#'
#'@export
plotStabilityPaths <- function(stabs_object, cols=NULL, lwd = 1, lty=1, ylim=c(0,1.1), ...) {

  # ... checks
  if (!base::inherits(stabs_object, what="stabsel")){stop("stabs_object must be of class 'stabsel', the resulting object from running stability selection with the `stabs` package")}

  # set plot parameters
  mat <- t(stabs_object$phat)
  if(is.null(cols)){
    cols <- rep("black", ncol(mat))
    names(cols) <- rep("Not Selected", length(cols))
    cols[stabs_object$selected] <- "cadetblue"
    names(cols)[stabs_object$selected] <- "Selected"
  }

  # plot stability paths
  matplot(mat, col = cols, type = "l", lty = lty, ylab = "Selection Probability", xlab = "Regularization Step", ylim = ylim, lwd = lwd, ...)
  abline(h = stabs_object$cutoff, lty = 5, col = "red", lwd = lwd)
  legend("topleft", legend = c(unique(names(cols)), "cutoff"), col = c(unique(cols), "red"), lty = c(1, 1, 5), bty = "n", lwd = lwd)
  invisible(TRUE)

}


#'@title Barplot Selection Probabilities
#'
#'@description Create a bar plot of the selection probabilities in descending order.
#'
#'@param stabs_object the \code{stabs} object resulting from stability selection.
#'@param ylim the limits for the y-axis.
#'@param onlySelected logical (default=TRUE) indicating if only selected predictors' selection probabilities
#'    should be plotted.
#'@param sel_color color of the selected predictors.
#'@param las (2 by default) plot labels vertically or horizontally.
#'@param ... additional parameters for the \code{barplot} function.
#'
#'@seealso \code{\link[graphics]{barplot}}
#'
#'@return barplot of selection probabilities.
#'@export
plotSelectionProb <- function(stabs_object, ylim = c(0,1.1), onlySelected = TRUE, sel_color="cadetblue", las = 2, ...) {

  # ... checks
  if (!base::inherits(stabs_object, what="stabsel")) {stop("stabs_object must be of class 'stabsel', the resulting object from running stability selection with the `stabs` package")}

  phat <- t(stabs_object$phat)
  TF_prob <- phat[nrow(phat), ]
  cols <- rep("grey", length(TF_prob))
  cols[stabs_object$selected] <- sel_color

  if (onlySelected) {
    TF_prob <- TF_prob[stabs_object$selected]
    cols <- cols[stabs_object$selected]
  }

  # check if empty
  if (S4Vectors::isEmpty(TF_prob)) {stop("The input for the barplot is empty")}

  # order
  TF_prob <- TF_prob[order(TF_prob, decreasing = TRUE)]

  # plot
  graphics::barplot(TF_prob, ylim = ylim, ylab = "Selection Probability", las = las, col = cols, border = NA, ...)
  abline(h = stabs_object$cutoff, lty = 5, col = "red")
  legend("topright", legend = "cutoff", lty = 5, col = "red", bty = "n")
  invisible(TRUE)

}


#'@title Plot Directionality of Predictor Effect
#'
#'@description This function plots the selectiong probabilities of the chosen predictors (for example the selected motifs)
#'and assigns a + or - sign to these probabilities to give a sense of directionality of the effect. The assumption is that 
#'the response vector on which stability selection was performed is a measure of fold-change. The correlation (pearson by default)
#'of each predictor to the response vector is calculated. The selection probabilities of the chosen predictors multiplied by the 
#'sign of the correlation is plotted to indicate the directionality.
#'
#'@param stabs_obj the \code{stabs} object resulting from stability selection.
#'@param response the response vector that was used for the stability selection (like the log-fold change of a measure of interest).
#'@param predictor_matrix the predictor matrix that was used for the stability selection (like the number of predicted TFBS of all motifs across the regions of interest).
#'@param sel_color the color for the selected predictors from stability selection.
#'@param min_sel_prob predictors with a selection probability greater than or equal to this are included in the plot.
#'@param cor_method the correlation method to be used.
#'@param ... additional parameters for the \code{barplot} function.
#'
#'@seealso \code{\link[graphics]{barplot}}
#'
#'@return a barplot indicating the directionality of the motifs with respect to the correlation to the response vector.
#'
#'@export
plotMotifDirectionality <- function(stabs_obj = NULL, response = NULL, predictor_matrix = NULL, sel_color="cadetblue", min_sel_prob=0.4, cor_method="pearson", ...) {
    
    # checks
    # ... NULL checks
    stopifnot(!is.null(stabs_obj))
    stopifnot(!is.null(response))
    stopifnot(!is.null(predictor_matrix))
    # ... class checks
    stopifnot(class(stabs_obj)=="stabsel")
    stopifnot(class(response)=="numeric")
    stopifnot(any(class(predictor_matrix)=="matrix"))
    # ... compatibility checks
    stopifnot(length(response)==nrow(predictor_matrix))
    if(!is.null(colnames(predictor_matrix))&!is.null(names(response))){
        stopifnot(all(rownames(stabs_obj$phat)==colnames(predictor_matrix)))
    }
    
    # correlation 
    cor <- as.vector(stats::cor(x = response, y = predictor_matrix, method = cor_method))
    cols <- rep("grey", ncol(predictor_matrix))
    cols[stabs_obj$selected] <- sel_color
    if(!is.null(colnames(predictor_matrix))) {
        tf_names <- colnames(predictor_matrix)
    } else {
        tf_names <- paste0("pred", 1:ncol(predictor_matrix))
    }
    # probabilities with directionality
    probs <- stabs_obj$phat[, ncol(stabs_obj$phat)]
    probs <- probs*sign(cor)
    
    # kept and ordered
    keep <- stabs_obj$phat[,ncol(stabs_obj$phat)]>=min_sel_prob
    cor <- cor[keep]
    cols <- cols[keep]
    tf_names <- tf_names[keep]
    probs <- probs[keep]
    o <- order(probs, decreasing = TRUE)
    cor <- cor[o]
    cols <- cols[o]
    tf_names <- tf_names[o]
    probs <- probs[o]
    up <- probs>0
    
    # plot
    bar <- graphics::barplot(probs, col = cols, border = NA, ylab = "Sel Prob * sign(cor to response)", names.arg = NA, 
                             ylim = c(min(0, range(probs)[1]-abs(0.3*range(probs)[1])), max(1, range(probs[2]+abs(0.3*range(probs)[2])))), ...)
    legend("topright", bty = "n", lty = 1, legend = c("selected", "not selected"), col = c(sel_color, "grey"))
    
    if(!(sum(up)==0)&!isEmpty(probs[up])){
      graphics::text(x = bar[up], y = probs[up], labels = tf_names[up], col = cols[up], xpd = TRUE, srt=90, adj = 0)
    }
    if(!(sum(!up)==0)&!isEmpty(probs[!up])){
      graphics::text(x = bar[!up], y = probs[!up], labels = tf_names[!up], col = cols[!up], xpd = TRUE, srt=90, adj = 1)
    }
    
    invisible(TRUE)
    
}

















