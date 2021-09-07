#' @title Binned motif enrichment analysis.
#' 
#' @description Calculate motif enrichment analysis on bins of sequences.
#'   Motif occurrences in each bin are contrasted to the motif occurrences
#'   in a background set of sequences (by default the sequences in all other
#'   bins).
#' 
#' @param seqs The sequences to analyze. Either a \code{DNAStringSet} (for
#'   \code{method = "R"}) or a \code{GRanges} object defining the genomic
#'   regions from which to extract the sequences (for \code{method = "Homer"}).
#' @param bins Factor of the same length and order as \code{seqs}, indicating
#'   the bin for each sequence or genomic region. Typically the return value of
#'   \code{\link[monaLisa]{bin}}.
#' @param motifs The motifs for which to calculate enrichments. Either a
#'   \code{PWMatrixList} (for \code{method = "R"}) or a file with motifs in
#'   Homer format (for \code{method = "Homer"}).
#' @param method A \code{character} scalar specifying the backend to use for
#'   enrichment calculations. One of \code{"R"} (default) or \code{"Homer"}.
#' @param pseudocount.log2enr A numerical scalar with the pseudocount to add to
#'   foreground and background counts when calculating log2 motif enrichments
#' @param pseudofreq.pearsonResid A numerical scalar with the pseudo-frequency
#'   to add to background frequencies when calculating Pearson residuals.
#'   The value needs to be in [0,1] and corresponds to the minimal expected
#'   fraction of background sequences that contain at least one motif hit.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#' @param BPPARAM Specifies the number of CPU cores to use for parallel
#'   processing. Either a \code{\link[BiocParallel]{SerialParam}} or a
#'   \code{\link[BiocParallel]{MulticoreParam}} object instance, or for
#'   \code{method = "Homer"} a numeric scalar.
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
#' @param ... Additional arguments for passed to
#'   \code{\link{calcBinnedMotifEnrR}} or
#'   \code{\link{calcBinnedMotifEnrHomer}}.
#'
#' @details This function is a wrapper for
#'   \code{\link{calcBinnedMotifEnrR}} or
#'   \code{\link{calcBinnedMotifEnrHomer}}. For additional supported arguments
#'   (via \code{...}) please see the help pages of these functions.
#'   
#' @return A \code{SummarizedExperiment} object with motifs in rows and bins
#'   in columns, containing six assays: \itemize{
#'   \item{negLog10P}{: -log10 P values}
#'   \item{negLog10Padj}{: -log10 adjusted P values}
#'   \item{pearsonResid}{: motif enrichments as Pearson residuals}
#'   \item{log2enr}{: motif enrichments as log2 ratios}
#'   \item{sumForegroundWgtWithHits}{: Sum of foreground sequence weights
#'     in a bin that have motif hits}
#'   \item{sumBackgroundWgtWithHits}{: Sum of background sequence weights
#'     in a bin that have motif hits}
#' }
#' 
#' @examples 
#' seqs <- Biostrings::DNAStringSet(c("GTCAGTCGATC", "CAGTCTAGCTG",
#'                                    "CGATCGTCAGT", "AGCTGCAGTCT"))
#' bins <- factor(rep(1:2, each = 2))
#' m <- rbind(A = c(2, 0, 0),
#'            C = c(1, 1, 0),
#'            G = c(0, 2, 0),
#'            T = c(0, 0, 3))
#' pwms <- TFBSTools::PWMatrixList(
#'     TFBSTools::PWMatrix(ID = "m1", profileMatrix = m),
#'     TFBSTools::PWMatrix(ID = "m2", profileMatrix = m[, 3:1])
#' )
#' calcBinnedMotifEnr(seqs = seqs, bins = bins, motifs = pwms,
#'                    method = "R", min.score = 3)

#'
#' @seealso \code{\link{bin}} for binning of sequences;
#'   \code{\link{calcBinnedMotifEnrR}} that is used for
#'   \code{method = "R"} and \code{\link{calcBinnedMotifEnrHomer}} used
#'   for \code{method = "Homer"}.
#' 
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam bpnworkers
#'
#' @export
calcBinnedMotifEnr <- function(seqs,
                               bins,
                               motifs,
                               method =  c("R", "Homer"),
                               pseudocount.log2enr = 8,
                               pseudofreq.pearsonResid = 0.001,
                               p.adjust.method = "BH",
                               BPPARAM = SerialParam(),
                               verbose = FALSE,
                               ...) {
    method <- match.arg(method)
    if (identical(method, "R")) {

        if (verbose)
            message("using R backend for binned motif enrichment analysis")

        se <- calcBinnedMotifEnrR(seqs = seqs,
                                  bins = bins,
                                  pwmL = motifs,
                                  pseudocount.log2enr = pseudocount.log2enr,
                                  pseudofreq.pearsonResid = pseudofreq.pearsonResid,
                                  p.adjust.method = p.adjust.method,
                                  BPPARAM = BPPARAM,
                                  verbose = verbose,
                                  ...)

    } else if (identical(method, "Homer")) {
        
        if (verbose)
            message("using Homer backend for binned motif enrichment analysis")
        
        if (is(BPPARAM, "SerialParam") || is(BPPARAM, "MulticoreParam")) {
            ncpu <- bpnworkers(BPPARAM)
        } else if (is.numeric(BPPARAM) && length(BPPARAM) == 1L && BPPARAM > 0) {
            ncpu <- as.integer(BPPARAM)
        } else {
            stop("For method = 'Homer', BPPARAM must be either a numeric scalar, ",
                 "a SerialParam() or a MulticoreParam() object.")
        }

        se <- calcBinnedMotifEnrHomer(gr = seqs,
                                      b = bins,
                                      motifFile = motifs,
                                      pseudocount.log2enr = pseudocount.log2enr,
                                      pseudofreq.pearsonResid = pseudofreq.pearsonResid,
                                      p.adjust.method = p.adjust.method,
                                      Ncpu = ncpu,
                                      verbose = verbose,
                                      ...)
    }
    
    return(se)
}
