#' Plot diagnostics of binned sequences
#'
#' Plot various diagnostics of binned sequences. Three plot types are
#' available:
#' \describe{
#' \item{\code{length}}{ plots the distribution of sequence lengths within
#'   each bin.}
#' \item{\code{GCfrac}}{ plots the distribution of GC fractions within each
#'   bin.}
#' \item{\code{dinucfreq}}{ plots a heatmap of the relative frequency of each
#'   dinucleotide, averaged across the sequences within each bin. The values
#'   are centered for each dinucleotide to better highlight differences
#'   between the bins. The average relative frequency of each dinucleotide
#'   (across the bins) is indicated as well.}
#' }
#' @param seqs \code{\link[Biostrings]{DNAStringSet}} object with sequences.
#' @param bins Factor of the same length and order as \code{seqs}, indicating
#'   the bin for each sequence. Typically the return value of \code{bin}.
#' @param aspect The diagnostic to plot. Should be one of \code{"length"},
#'   \code{"GCfrac"} and \code{"dinucfreq"}, to plot the distribution of
#'   sequence lengths, the distribution of GC fractions and the average
#'   relative dinucleotide frequencies across the bins.
#' @param draw_quantiles For aspect=\code{"length"} or \code{"GCfrac"},
#'   draw vertical lines at the given quantiles of the density estimate.
#'   If \code{NULL}, no quantile lines will be drawn.
#' @param ... Additional argument passed to \code{getColsByBin}.
#'
#' @export
#'
#' @return For aspect=\code{"length"} or \code{"GCfrac"}, returns a
#'   \code{\link[ggplot2]{ggplot}} object. For aspect=\code{"dinucfreq"},
#'   returns (invisibly) a \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @examples
#' seqs <- Biostrings::DNAStringSet(
#'   vapply(1:250, function(i) paste(sample(x = c("A", "C", "G", "T"),
#'                                          size = round(stats::rnorm(1, 20, 5)),
#'                                          replace = TRUE), collapse = ""), "")
#' )
#' bins <- factor(rep(c("a", "b", "c", "d", "e"), each = 50))
#' plotBinDiagnostics(seqs, bins, aspect = "length")
#' plotBinDiagnostics(seqs, bins, aspect = "GCfrac", draw_quantiles = NULL)
#' plotBinDiagnostics(seqs, bins, aspect = "dinucfreq")
#'
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_manual
#'   scale_colour_manual labs element_blank theme_classic
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom BiocGenerics width
#' @importFrom cli cli_abort
#' @importFrom rlang .data
plotBinDiagnostics <- function(seqs, bins,
                               aspect = c("length", "GCfrac", "dinucfreq"),
                               draw_quantiles = c(0.25, 0.5, 0.75),
                               ...) {
    .assertVector(x = seqs, type = "DNAStringSet")
    .assertVector(x = bins, type = "factor")
    if (length(seqs) != length(bins)) {
        cli_abort("{.arg seqs} and {.arg bins} must have equal length and order")
    }
    aspect <- match.arg(aspect)
    .assertVector(x = draw_quantiles, type = "numeric", rngIncl = c(0, 1),
                  allowNULL = TRUE)

    binCols <- getColsByBin(bins, ...)

    if (aspect %in% c("length", "GCfrac")) {
        pd <- data.frame(bin = bins)
        if (identical(aspect, "length")) {
            pd$xvalue <- width(seqs)
            xlab <- "Length (bp)"
        } else if (identical(aspect, "GCfrac")) {
            onf <- oligonucleotideFrequency(seqs, width = 1, as.prob = TRUE)
            pd$xvalue <- onf[, "G"] + onf[, "C"]
            xlab <- "GC fraction"
        }
        p <- ggplot(data = pd,
                    mapping = aes(.data$xvalue, .data$bin, fill = .data$bin)) +
            geom_violin(draw_quantiles = draw_quantiles, show.legend = FALSE) +
            scale_fill_manual(values = attr(binCols, "cols")) +
            labs(x = xlab, y = element_blank()) +
            theme_classic()
        return(p)
    } else if (aspect == "dinucfreq") {
        dnf <- Biostrings::oligonucleotideFrequency(seqs, width = 2,
                                                    as.prob = TRUE)
        dnfsplit <- lapply(split.data.frame(dnf, bins), colMeans)
        dnfmat <- t(do.call(rbind, dnfsplit))
        cols <- circlize::colorRamp2(breaks = range(rowMeans(dnfmat)),
                                     colors = c("white", "grey50"))
        hm <- ComplexHeatmap::Heatmap(
            t(scale(t(dnfmat), center = TRUE, scale = FALSE)),
            right_annotation = ComplexHeatmap::rowAnnotation(
                OverallFreq = rowMeans(dnfmat),
                col = list(OverallFreq = cols),
                show_annotation_name = FALSE,
                annotation_label = list(
                    OverallFreq = "Average\nrelative\nfrequency\nacross bins")),
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            name = "Mean\nrelative\nfrequency\nin bin\n(centered)"
        )
        show(hm)
        return(invisible(hm))
    }
}
