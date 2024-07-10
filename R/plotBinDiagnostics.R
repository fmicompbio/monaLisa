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
#' @param seqs DNAStringSet object with sequences.
#' @param bins factor of the same length and order as seqs, indicating the bin
#'   for each sequence. Typically the return value of \code{bin}.
#' @param aspect The diagnostic to plot. Should be one of \code{"length"},
#'   \code{"GCfrac"} and \code{"dinucfreq"}, to plot the distribution of
#'   sequence lengths, the distribution of GC fractions and the average
#'   relative dinucleotide frequencies across the bins.
#' @param ... Additional argument passed to \code{getColsByBin}.
#'
#' @export
#'
#' @return For aspect=\code{"length"} or \code{"GCfrac"}, returns (invisibly)
#'   the output of \code{vioplot()}, which generates the plot. For
#'   aspect=\code{"dinucfreq"}, returns (invisibly) the ComplexHeatmap
#'   object.
#'
#' @examples
#' seqs <- Biostrings::DNAStringSet(
#'   vapply(1:100, function(i) paste(sample(c("A", "C", "G", "T"), 10,
#'                                          replace = TRUE), collapse = ""), "")
#' )
#' bins <- factor(rep(1:2, each = 50))
#' plotBinDiagnostics(seqs, bins, aspect = "GCfrac")
#' plotBinDiagnostics(seqs, bins, aspect = "dinucfreq")
#'
#' @importFrom vioplot vioplot
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom Biostrings oligonucleotideFrequency
plotBinDiagnostics <- function(seqs, bins,
                               aspect = c("length", "GCfrac", "dinucfreq"),
                               ...) {
    .assertVector(x = seqs, type = "DNAStringSet")
    .assertVector(x = bins, type = "factor")
    if (length(seqs) != length(bins)) {
        stop("'seqs' and 'bins' must be of equal length and in the same order")
    }

    aspect <- match.arg(aspect)
    binCols <- getColsByBin(bins, ...)
    if (aspect == "length") {
        vp <- vioplot::vioplot(split(width(seqs), bins), ylab = "",
                               col = attr(binCols, "cols"),
                               xlab = "Length", axes = FALSE, las = 2,
                               horizontal = TRUE, cex = 0.5, cex.axis = 0.75,
                               cex.names = 0.75, par(mar = c(4, 6, 2, 2) + 0.1))
        return(invisible(vp))
    } else if (aspect == "GCfrac") {
        onf <- Biostrings::oligonucleotideFrequency(seqs, width = 1,
                                                    as.prob = TRUE)
        gcfrac <- onf[, "G"] + onf[, "C"]
        vp <- vioplot::vioplot(split(gcfrac, bins), ylab = "",
                               col = attr(binCols, "cols"),
                               xlab = "GC fraction", axes = FALSE, las = 2,
                               horizontal = TRUE, cex = 0.5, cex.axis = 0.75,
                               cex.names = 0.75, par(mar = c(4, 6, 2, 2) + 0.1))
        return(invisible(vp))
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
