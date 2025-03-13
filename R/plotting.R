#' @importFrom grDevices colorRampPalette
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
#' @importFrom grDevices col2rgb
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
#' @param breaks A \code{numeric} scalar controlling the histogram breaks
#'     (passed to \code{geom_hist(..., bins = breaks)}).
#' @param xlab,ylab,main \code{character} scalars that set the x-axis label,
#'     y-axis label and the main title. Use \code{""} to suppress the label.
#' @param legendPosition A \code{character} scalar.
#'     If not \code{"none"}, draw a legend with binning information. The value
#'     is used to control the legend position and will be passed to
#'     \code{theme(legend.position = legendPosition)}.
#' @param legend Deprecated (ignored). Please use \code{legendPosition} to
#'     control the drawing and position of the legend.
#' @param legend.cex Deprecated (ignored). You can use
#'     \code{\link[ggplot2]{theme}} to set legend and other graphical
#'     parameters.
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}, \code{\link[ggplot2]{geom_histogram}}
#'
#' @return The generated histogram as a \code{ggplot} object.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' b <- bin(x, "equalN", nElements = 10)
#' plotBinHist(x, b)
#'
#' @importFrom ggplot2 ggplot aes geom_histogram element_blank theme_classic theme
#' @importFrom cli cli_warn
#'
#' @export
plotBinHist <- function(x, b,
                        breaks = 6 * nlevels(b),
                        xlab = deparse(substitute(x, env = as.environment(-1))),
                        ylab = "Frequency",
                        main = "",
                        legendPosition = "right",
                        legend = NULL,
                        legend.cex = NULL, ...) {
    .assertVector(x = x, type = "numeric")
    .assertVector(x = b, type = "factor", len = length(x))
    .assertScalar(x = breaks, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = xlab, type = "character")
    .assertScalar(x = ylab, type = "character")
    .assertScalar(x = main, type = "character")
    if (!is.null(legend)) {
        cli_warn(c("{.arg legend} is deprecated and ignored.",
                   "i" = "You can use {.arg legendPosition} to control legend drawing and position"))
    }
    if (!is.null(legend.cex)) {
        cli_warn(c("{.arg legend.cex} is deprecated and ignored.",
                   "i" = "You can use {.fn theme} to control legend and other graphical parameters"))
    }

    cols <- getColsByBin(b, ...)
    bincols <- attr(cols, "cols")

    # add number of elements to bin names
    bn <- unclass(table(b))
    levels(b) <- names(bincols) <- paste0(levels(b), ": ", bn)

    p <- ggplot(data = data.frame(x = x, b = b), mapping = aes(x = x, fill = b)) +
        geom_histogram(bins = breaks, colour = "gray20") +
        scale_fill_manual(values = bincols) +
        labs(x = ifelse(xlab != "", xlab, element_blank()),
             y = ifelse(ylab != "", ylab, element_blank()),
             main = ifelse(main != "", main, element_blank()),
             fill = "Bins") +
        theme_classic() +
        theme(legend.position = legendPosition)

    return(p)
}


#' @title Density plot of binned elements.
#'
#' @description Plot the density of binned elements with binning information.
#'
#' @param x A numerical vector with the values used for binning.
#' @param b A factor that groups elements of \code{x} into bins (typically the
#'     output of \code{\link{bin}}).
#' @param xlab,ylab,main \code{character} scalars that set the x-axis label,
#'     y-axis label and the main title. Use \code{""} to suppress the label.
#' @param legendPosition A \code{character} scalar.
#'     If not \code{"none"}, draw a legend with binning information. The value
#'     is used to control the legend position and will be passed to
#'     \code{theme(legend.position = legendPosition)}.
#' @param legend Deprecated (ignored). Please use \code{legendPosition} to
#'     control the drawing and position of the legend.
#' @param legend.cex Deprecated (ignored). You can use
#'     \code{\link[ggplot2]{theme}} to set legend and other graphical
#'     parameters.
#' @param ... Further arguments passed to \code{\link{getColsByBin}}.
#'
#' @seealso \code{\link{getColsByBin}}, \code{\link[ggplot2]{geom_density}}
#'
#' @return The generated density plot as a \code{ggplot} object.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' b <- bin(x, "equalN", nElements = 10)
#' plotBinDensity(x, b)
#'
#' @importFrom stats density
#' @importFrom ggplot2 ggplot aes geom_area geom_line geom_rug element_blank
#'     theme_classic theme
#' @importFrom cli cli_warn
#' @importFrom rlang .data
#'
#' @export
plotBinDensity <- function(x, b,
                           xlab = deparse(
                               substitute(x, env = as.environment(-1))),
                           ylab = "Density",
                           main = "",
                           legendPosition = "right",
                           legend = NULL,
                           legend.cex = NULL,
                           ...) {
    .assertVector(x = x, type = "numeric")
    .assertVector(x = b, type = "factor", len = length(x))
    stopifnot("breaks" %in% names(attributes(b)))
    .assertScalar(x = xlab, type = "character")
    .assertScalar(x = ylab, type = "character")
    .assertScalar(x = main, type = "character")
    if (!is.null(legend)) {
        cli_warn(c("{.arg legend} is deprecated and ignored.",
                   "i" = "You can use {.arg legendPosition} to control legend drawing and position"))
    }
    if (!is.null(legend.cex)) {
        cli_warn(c("{.arg legend.cex} is deprecated and ignored.",
                   "i" = "You can use {.fn theme} to control legend and other graphical parameters"))
    }

    cols <- getColsByBin(b, ...)
    binbreaks <- attr(b, "breaks")
    bincols <- attr(cols, "cols")

    # add number of elements to bin names
    bn <- unclass(table(b))
    levels(b) <- names(bincols) <- paste0(levels(b), ": ", bn)

    # calculate density and assign bins
    dens <- data.frame(density(x = x)[c("x", "y")])
    binbreaks[c(1, length(binbreaks))] <- range(dens$x)
    dens$b <- factor(x = levels(b)[cut(dens$x, breaks = binbreaks,
                                       include.lowest = TRUE, labels = FALSE)],
                     levels = levels(b))

    p <- ggplot(dens) +
        geom_area(data = dens, mapping = aes(.data$x, .data$y,
                                             colour = .data$b, fill = .data$b),
                  outline.type = "full") +
        geom_line(data = dens, mapping = aes(.data$x, .data$y), colour = "gray20") +
        geom_rug(data = data.frame(x = binbreaks),
                 mapping = aes(.data$x), colour = "gray20") +
        scale_colour_manual(values = bincols) +
        scale_fill_manual(values = bincols) +
        labs(x = ifelse(xlab != "", xlab, element_blank()),
             y = ifelse(ylab != "", ylab, element_blank()),
             main = ifelse(main != "", main, element_blank()),
             colour = "Bins",
             fill = "Bins") +
        theme_classic() +
        theme(legend.position = legendPosition)

    return(p)
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
#' @param cols \code{NULL} or a color vector defining the colors of points.
#'     If \code{NULL}, the colors will be computed based on \code{b} using
#'     \code{\link{getColsByBin}(b)}).
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param main Main title.
#' @param legendPosition A \code{character} scalar.
#'     If not \code{"none"}, draw a legend with binning information. The value
#'     is used to control the legend position and will be passed to
#'     \code{theme(legend.position = legendPosition)}.
#' @param legend Deprecated (ignored). Please use \code{legendPosition} to
#'     control the drawing and position of the legend.
#' @param legend.cex Deprecated (ignored). You can use
#'     \code{\link[ggplot2]{theme}} to set legend and other graphical
#'     parameters.
#' @param ... Further arguments passed to \code{\link{getColsByBin}} (only used
#'     if \code{cols} is \code{NULL}).
#'
#' @seealso \code{\link{bin}}, \code{\link{getColsByBin}},
#'     \code{\link[ggplot2]{geom_point}}
#'
#' @return The generated scatter plot as a \code{ggplot} object.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' b <- bin(y, "equalN", nElements = 10)
#' plotBinScatter(x, y, b)
#'
#' @importFrom ggplot2 ggplot aes geom_point element_blank theme_classic theme
#' @importFrom cli cli_warn
#'
#' @export
plotBinScatter <- function(x, y, b,
                           cols = NULL,
                           xlab = deparse(
                               substitute(x, env = as.environment(-1))),
                           ylab = deparse(
                               substitute(y, env = as.environment(-1))),
                           main = "",
                           legendPosition = "right",
                           legend = NULL,
                           legend.cex = NULL,
                           ...) {
    .assertVector(x = x, type = "numeric")
    .assertVector(x = y, type = "numeric", len = length(x))
    .assertVector(x = b, len = length(x))
    bincols <- NULL
    if (is.null(cols)) {
        cols <- getColsByBin(b, ...)
        bincols <- attr(cols, "cols")
    } else if (length(cols) == 1L) {
        cols <- rep(cols, length(x))
    }
    .assertVector(x = cols, len = length(x))
    .assertScalar(x = xlab, type = "character")
    .assertScalar(x = ylab, type = "character")
    .assertScalar(x = main, type = "character")
    if (is.null(bincols) && !identical(legendPosition, "none")) {
        cli_warn(paste0(
            "Setting {.arg legendPosition} to 'none' - ",
            "cannot use custom colors ({.arg cols}) with bin-based legend"))
    }
    if (!is.null(legend)) {
        cli_warn(c("{.arg legend} is deprecated and ignored.",
                   "i" = "You can use {.arg legendPosition} to control legend drawing and position"))
    }
    if (!is.null(legend.cex)) {
        cli_warn(c("{.arg legend.cex} is deprecated and ignored.",
                   "i" = "You can use {.fn theme} to control legend and other graphical parameters"))
    }

    # add number of elements to bin names
    bn <- unclass(table(b))
    levels(b) <- paste0(levels(b), ": ", bn)

    p <- ggplot(data = data.frame(x = x, y = y, b = b, cols = cols),
                mapping = aes(x, y)) +
        labs(x = ifelse(xlab != "", xlab, element_blank()),
             y = ifelse(ylab != "", ylab, element_blank()),
             main = ifelse(main != "", main, element_blank()),
             colour = "Bins",
             fill = "Bins") +
        theme_classic() +
        theme(legend.position = legendPosition)

    if (is.null(bincols)) {
        p <- p + geom_point(colour = cols)
    } else {
        names(bincols) <- levels(b)
        p <- p + geom_point(aes(colour = b)) +
            scale_colour_manual(values = bincols)
    }

    return(p)
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
#' @param show_bin_legend If \code{TRUE}, show a legend for the bin labels.
#'     If FALSE (default), the bin legend will be hidden.
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
#'     when creating the main heatmaps selected by \code{which.plots}. For
#'     example, the following will set the font size of the motif names:
#'     \code{plotMotifHeatmaps(..., row_names_gp = gpar(fontsize = 12))}
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
#' @importFrom cli cli_abort
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
                              show_bin_legend = FALSE,
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
    .assertScalar(x = show_bin_legend, type = "logical")
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
        cli_abort("{.arg cluster} must be either {.code TRUE}, {.code FALSE} or an {.cls hclust}-object.")
    }
    hmBin <- HeatmapAnnotation(
        df = data.frame(bin = factor(colnames(x),
                                     levels = colnames(x))),
        name = "bin",
        col = list(bin = bincols),
        show_annotation_name = FALSE,
        which = "column", width = unit(width,"inch"),
        annotation_height = unit(width / 16, "inch"),
        show_legend = show_bin_legend)
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
        show_heatmap_legend = FALSE, left_annotation = hmSeqlogo,
        ...
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
            use_raster = use_raster,
            ...
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
#' @param selColor color for the selected predictors which have a selection
#'    probability greater than \code{selProbMin}.
#' @param notSelColor color for the rest of the (un-selected) predictors.
#' @param selProbCutoffColor color for the line depicting the selection 
#'    probability cutoff.
#' @param linewidth line width.
#' @param alpha line transparency of the stability paths.
#' @param ylim limits for y-axis.
#' @param labelPaths if TRUE, the predictor labels will be shown at the end of 
#'    the stability paths. The predictor labels given in \code{labels} will
#'    be shown. If unspecified, the labels corresponding to the selected
#'    predictors will be added.
#' @param labels the predictors which should be labelled. If \code{NULL}, the 
#'    selected predictors greater than \code{metadata(se)$stabsel.params.cutoff}  
#'    will be shown.
#'
#' @return a \code{ggplot2} object.
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
#' @seealso \code{\link[stabs]{stabsel}}
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual geom_segment labs
#'     scale_y_continuous scale_x_continuous guides guide_legend 
#'     theme_classic
#' @importFrom tidyr pivot_longer starts_with
#' @importFrom cli cli_abort
#'
#' @export
plotStabilityPaths <- function(se,
                               selProbMin = metadata(se)$stabsel.params.cutoff,
                               selColor = "cadetblue", 
                               notSelColor = "grey",
                               selProbCutoffColor = "firebrick",
                               linewidth = 0.5, 
                               alpha = 1, 
                               ylim = c(0, 1), 
                               labelPaths = FALSE, 
                               labels = NULL) {
    # checks
    if (!is(se, "SummarizedExperiment")) {
        cli_abort("{.arg se} must be a {.cls SummarizedExperiment}")
    }
    .assertScalar(x = selProbMin, type = "numeric", rngIncl = c(0, 1))
    .assertColor(x = selColor, len = 1)
    .assertColor(x = notSelColor, len = 1)
    .assertColor(x = selProbCutoffColor, len = 1)
    .assertScalar(x = linewidth, type = "numeric")
    .assertScalar(x = alpha, type = "numeric", rngIncl = c(0, 1))
    .assertVector(x = ylim, type = "numeric", len = 2)
    .assertVector(x = colnames(se), type = "character", len = ncol(se))
    .assertVector(x = rownames(se), type = "character", len = nrow(se))
    .assertVector(x = colnames(colData(se)), type = "character")
    if (!any(grepl(pattern = "regStep", x = colnames(colData(se))))) {
        cli_abort("the columns in {.code colData(se)} containing the selection 
                probabilities for each regularization step {.code i} must have
                column names starting with {.emph regStep}. See 
                {.fn randLassoStabSel} for more details.")
    }
    if (is.null(rownames(colData(se)))) {
        cli_abort("{.code rownames(colData(se))} must not be empty.")
    }

    # prepare dataframe to plot
    df <- as.data.frame(colData(se))
    df$predictor <- rownames(df)
    df$selected <- df$selProb > selProbMin
    df <- tidyr::pivot_longer(df,
                              cols = tidyr::starts_with(match = "regStep"),
                              names_to = "regStep",
                              values_to = "selectionProbability")
    df$regStep <- as.integer(gsub(pattern = "regStep",
                                  replacement = "",
                                  x = df$regStep))
    df$yintercept <- selProbMin
    df$cutoff <- paste0("selProbMin = ", selProbMin)

    # plot stability paths (returns ggplot object)
    gg <- ggplot(data = df, 
           mapping = aes(x = .data$regStep, y = .data$selectionProbability)) +
        geom_line(mapping = aes(group = .data$predictor, color = .data$selected),
                  linewidth = linewidth, 
                  alpha = alpha) +
        scale_color_manual(values = c("TRUE" = selColor, 
                                      "FALSE" = notSelColor)) +
        geom_segment(mapping = aes(x = min(.data$regStep), xend = max(.data$regStep) + 1, 
                                 y = .data$yintercept, yend = .data$yintercept, 
                                 linetype = .data$cutoff),
                   color = selProbCutoffColor, 
                   linewidth = linewidth) +
        labs(x = "Regularization Step",
             y = "Selection Probability") +
        scale_y_continuous(expand = c(0, 0), limits = ylim) + 
        scale_x_continuous(expand = c(0, 0)) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme_classic()
    if(labelPaths){
      .assertPackagesAvailable("ggrepel")
      
      if(is.null(labels)){
        labels <- unique(df$predictor[df$selected])
      }
      else {
        # check that the labels are correct and exist
        .assertVector(x = labels, type = "character", 
                      validValues = unique(df$predictor))
      }
      gg <- gg + 
        ggrepel::geom_text_repel(data = df[(df$predictor %in% labels) & 
                                             (df$regStep == max(df$regStep)), ], 
                                 aes(label = .data$predictor, 
                                     color = .data$selected), 
                                 nudge_x = 8, na.rm = TRUE, 
                                 size = 3, direction = "y", 
                                 hjust = 0, segment.linetype = "dotted", 
                                 segment.size = 0.7, segment.curvature = -0.1, 
                                 segment.angle = 20, box.padding = 0.4, 
                                 segment.alpha = 0.5, vjust = 0, 
                                 show.legend = FALSE) 
      
    }
    gg
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
#'   probability greater than \code{selProbMin} are considered selected and 
#'   colored by the input from \code{selColor}. By default, \code{selProbMin} is
#'   extracted from the parameters stored in \code{se}.
#' @param selProbMinPlot A numerical scalar in [0,1] less than 
#'   \code{selProbMin}. Predictors with a selection probability greater 
#'   than \code{selProbMinPlot} but less than \code{selProbMin} are shown as 
#'   bars with color defined in \code{notSelColor}. \code{selProbMinPlot} is 
#'   useful in order to include additional predictors in the barplot, that were 
#'   not selected according to \code{selProbMin} but may be close to that cutoff
#'   or are simply nice to visualize alongside the selected predictors. Setting 
#'   \code{selProbMinPlot = 0} will include all predictors.
#' @param showSelProbMin A logical scalar. If \code{TRUE}, the value of
#'   \code{selProbMin} is shown by a horizontal line with the color defined 
#'   by \code{selProbCutoffColor}.
#' @param selColor Color for the selected predictors which have a selection
#'    probability greater than \code{selProbMin}.
#' @param notSelColor Color for the rest of the (unselected) predictors which 
#'    will be show in the barplot.
#' @param selProbCutoffColor Color for the line depicting the selection 
#'    probability cutoff.
#' @param method A character scalar with the correlation method to use in the
#'   calculation of predictor-response marginal correlations. One of "pearson",
#'   "kendall" or "spearman" (see \code{\link[stats]{cor}}).
#' @param ylimext A numeric scalar defining how much the y axis limits should be
#'   expanded beyond the plotted probabilities to allow for space for the
#'   bar labels. This value can be increased if the predictor names above the 
#'   bars are too long and not showing in the plot.
#'
#' @details This function creates a bar plot with \code{ggplot}.
#'   Each bar corresponds to a predictor (motif) and the colors correspond to 
#'   whether or not it was selected. The y-axis shows the selection 
#'   probabilities (\code{directional=FALSE}) or selection probabilities with 
#'   the sign of the marginal correlation to the response 
#'   (\code{directional=TRUE}). 
#'
#' @return a \code{ggplot2} object.
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
#' @importFrom stats cor reorder
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual geom_hline labs 
#'     geom_text ylim theme theme_classic
#'
#' @export
plotSelectionProb <- function(se,
                              directional = TRUE,
                              selProbMin = metadata(se)$stabsel.params.cutoff,
                              selProbMinPlot = 0.4,
                              showSelProbMin = TRUE,
                              selColor = "cadetblue", 
                              notSelColor = "grey",
                              selProbCutoffColor = "firebrick", 
                              method = c("pearson", "kendall", "spearman"),
                              ylimext = 0.2) {

    # checks
    .assertScalar(x = directional, type = "logical")
    .assertScalar(x = selProbMin, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = selProbMinPlot, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = showSelProbMin, type = "logical")
    stopifnot(exprs = {
        is(se, "SummarizedExperiment")
        selProbMin >= selProbMinPlot
    })
    .assertColor(x = selColor, len = 1)
    .assertColor(x = notSelColor, len = 1)
    .assertColor(x = selProbCutoffColor, len = 1)
    method <- match.arg(method)
    .assertScalar(x = ylimext, type = "numeric", rngIncl = c(0, Inf))
    
    # prepare dataframe to plot
    df <- data.frame(
      corcoef = as.vector(cor(x = SummarizedExperiment::rowData(se)$y,
                              y = SummarizedExperiment::assay(se, "x"),
                              method = method)), 
      probs = se$selProb,
      predNames = colnames(se), 
      selected = se$selProb > selProbMin
    )
    df$yintercept <- selProbMin
    df$cutoff <- paste0("selProbMin = ", selProbMin)

    # ... multiply prob by the sign of corcoef if directional=TRUE
    if(directional){
      df$probs <- df$probs * sign(df$corcoef)
    }
    
    # ... only keep abs(probs) >= selProbMinPlot
    df <- df[abs(df$probs) >= selProbMinPlot, ]
    ymin <- min(min(df$probs), 0) - ylimext*max(df$probs)
    ymax <- max(df$probs) + ylimext*max(df$probs)
    
    # plot (returns ggplot object)
    gg <- ggplot(data = df, 
                 mapping = aes(x = stats::reorder(.data$predNames, 
                                                  -.data$probs), 
                               y = .data$probs)) + 
      geom_bar(mapping = aes(fill = .data$selected), stat = "identity") + 
      scale_fill_manual(values = c("TRUE" = selColor, "FALSE" = notSelColor)) + 
      labs(x = element_blank(), 
           y = ifelse(
             directional, 
             "Directional selection probability",
             "Selection probability"
           )) + 
      ylim(c(ymin, ymax)) + 
      geom_text(mapping = aes(label = .data$predNames, 
                              hjust = ifelse(.data$probs >= 0, 0, 1)),  
                angle = 90) + 
      theme_classic() + 
      theme(axis.text.x = element_blank(), 
            axis.line.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.ticks.length.y  = unit(0.2, "cm"))
    if(showSelProbMin){
      gg <- gg + 
        geom_hline(mapping = aes(yintercept = .data$yintercept, 
                                 linetype = .data$cutoff), 
                   color = selProbCutoffColor, linewidth = 1)
      if(directional & min(df$probs) <= -selProbMin){
        gg <- gg + 
          geom_hline(yintercept = -selProbMin, 
                     color = selProbCutoffColor, 
                     linewidth = 1)
      }

    }
    
    gg
    
}
