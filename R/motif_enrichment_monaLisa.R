#' @title Calculate motif enrichment
#'
#' @description Given motif counts, foreground/background labels and
#'   weights for a set of sequences, calculate the enrichment of each motif
#'   in foreground compared to background. This function is called by
#'   \code{calcBinnedMotifEnrR()} for each bin.
#'
#'   The default type of test is \code{"fisher"}, which is also what
#'   \code{Homer} uses if "-h" is specified for a hypergeometric test.
#'   Alternatively, a binomial test can be used by \code{test = "binomial"}
#'   (what \code{Homer} does by default). Using Fisher's exact test has 
#'   the advantage that special cases such as zero background counts are 
#'   handled without ad-hoc adjustments to the frequencies.
#'
#'   For \code{test = "fisher"}, \code{fisher.test} is used with
#'   \code{alternative = "greater"}, making it a one-sided test for enrichment,
#'   as is the case with the binomial test.
#'
#' @param motifHitMatrix matrix with 0 and 1 entries for absence or presence of
#'   motif hits in each sequence.
#' @param df a \code{DataFrame} with sequence information as returned by
#'   \code{.iterativeNormForKmers()}.
#' @param test type of motif enrichment test to perform.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @return a \code{data.frame} containing the motifs as rows and the columns:
#'   \itemize{
#'     \item{motifName}{: the motif name}
#'     \item{logP}{: the log p-value for enrichment (natural logarithm).
#'        If \code{test="binomial"} (default), this log p-value is identical to
#'        the one returned by Homer.}
#'     \item{sumForegroundWgtWithHits}{: the sum of the weights of the 
#'        foreground sequences that have at least one instance of a specific 
#'        motif (ZOOPS mode).}
#'     \item{sumBackgroundWgtWithHits}{: the sum of the weights of the 
#'        background sequences that have at least one instance of a specific 
#'        motif (ZOOPS mode).}
#'     \item{totalWgtForeground}{: the total sum of weights of foreground
#'        sequences.}
#'     \item{totalWgtBackground}{: the total sum of weights of background
#'        sequences.}
#'   }
#'
#' @importFrom stats pbinom fisher.test
#' 
#' @keywords internal
.calcMotifEnrichment <- function(motifHitMatrix,
                                 df,
                                 test = c("fisher", "binomial"),
                                 verbose = FALSE){

    # checks
    if (!is.matrix(motifHitMatrix)) {
        stop("'motifHitMatrix' has to be a matrix")
    }
    if (nrow(motifHitMatrix) != nrow(df)) {
        stop("'motifHitMatrix' and 'df' must have the same number of rows")
    }
    if (!is.null(rownames(motifHitMatrix)) && !is.null(rownames(df)) &&
        !identical(rownames(motifHitMatrix), rownames(df))) {
        stop("'motifHitMatrix' and 'df' must have identical rownames")
    }
    .checkDfValidity(df)
    test <- match.arg(test)
    .assertScalar(x = verbose, type = "logical")

    totalWgtForeground <- sum(df$seqWgt[df$isForeground])
    totalWgtBackground <- sum(df$seqWgt[!df$isForeground])

    motifHitMatrixWeighted <- motifHitMatrix * df$seqWgt
    TFmatchedSeqCountForeground <- 
      colSums(motifHitMatrixWeighted[df$isForeground, , drop = FALSE])
    TFmatchedSeqCountBackground <- 
      colSums(motifHitMatrixWeighted[!df$isForeground, , drop = FALSE])

    # calculate motif enrichment
    if (identical(test, "binomial")) {
        logP <- .binomEnrichmentTest(
            matchCountBg = TFmatchedSeqCountBackground,
            totalWeightBg = totalWgtBackground,
            matchCountFg = TFmatchedSeqCountForeground,
            totalWeightFg = totalWgtForeground,
            verbose = verbose
        )
    } else if (identical(test, "fisher")) {
        logP <- .fisherEnrichmentTest(
            matchCountBg = TFmatchedSeqCountBackground, 
            totalWeightBg = totalWgtBackground,
            matchCountFg = TFmatchedSeqCountForeground,
            totalWeightFg = totalWgtForeground,
            verbose = verbose
        )
    }

    return(data.frame(motifName = names(logP),
                      logP = logP,
                      sumForegroundWgtWithHits = TFmatchedSeqCountForeground,
                      sumBackgroundWgtWithHits = TFmatchedSeqCountBackground,
                      totalWgtForeground = totalWgtForeground,
                      totalWgtBackground = totalWgtBackground))
}


#' @title Binned Motif Enrichment Analysis with \code{monaLisa}
#'
#' @description This function performs a motif enrichment analysis on bins of
#'   sequences. For each bin, the sequences in all other bins are used as
#'   background.
#'
#' @param seqs \code{\link[Biostrings]{DNAStringSet}} object with sequences to 
#'   test
#' @param bins factor of the same length and order as \code{seqs}, indicating
#'   the bin for each sequence. Typically the return value of
#'   \code{\link[monaLisa]{bin}}. For \code{background = "genome"}, \code{bins}
#'   can be omitted.
#' @param pwmL PWMatrixList with motifs for which to calculate enrichments.
#' @param background A \code{character} scalar specifying the background 
#'   sequences to use. One of \code{"otherBins"} (default), \code{"allBins"}, 
#'   \code{"zeroBin"} or \code{"genome"} (see "Details").
#' @param test A \code{character} scalar specifying the type of enrichment test
#'   to perform. One of \code{"fisher"} (default) or \code{"binomial"}. The
#'   enrichment test is one-sided (enriched in foreground).
#' @param maxFracN A numeric scalar with the maximal fraction of N bases allowed
#'   in a sequence (defaults to 0.7). Sequences with higher fractions are
#'   excluded from the analysis.
#' @param maxKmerSize the maximum k-mer size to consider, when adjusting
#'   background sequence weights for k-mer composition compared to the
#'   foreground sequences. The default value (3) will correct for mono-, di-
#'   and tri-mer composition.
#' @param min.score the minimal score for motif hits, used in
#'   \code{\link[monaLisa]{findMotifHits}}.
#' @param matchMethod the method used to scan for motif hits, passed to the
#'   \code{method} parameter in \code{\link[monaLisa]{findMotifHits}}.
#' @param GCbreaks The breaks between GC bins. The default value is based on
#'   the hard-coded bins used in Homer.
#' @param pseudocount.log2enr A numerical scalar with the pseudocount to add to
#'   foreground and background counts when calculating log2 motif enrichments
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#' @param genome A \code{BSgenome} or \code{DNAStringSet} object with the
#'   genome sequence. Only used for \code{background = "genome"} for extracting
#'   background sequences.
#' @param genome.regions An optional \code{\link[GenomicRanges]{GRanges}} object
#'   defining the intervals in \code{genome} from which background sequences are
#'   sampled for \code{background = "genome"}. If \code{NULL}, background
#'   sequences are sampled randomly from \code{genome}.
#' @param genome.oversample A \code{numeric} scalar of at least 1.0 defining how
#'   many background sequences will be sampled per foreground sequence for
#'   \code{background = "genome"}. Larger values will take longer but improve
#'   the sequence composition similarity between foreground and background
#'   (see \code{"Details"}).
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation.
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
#' @param ... Additional arguments for  \code{\link[monaLisa]{findMotifHits}}.
#'
#' @details This function implements a binned motif enrichment analysis. In each
#'   enrichment analysis, the sequences in a specific bin are used as foreground
#'   sequences to test for motif enrichments comparing to background sequences
#'   (defined by \code{background}, see below). The logic follows the
#'   \code{findMotifsGenome.pl} tool from \code{Homer} version 4.11, with
#'   \code{-size given -nomotif -mknown} and additionally \code{-h} if using 
#'   \code{test = "fisher"}, and gives very similar results. As in the 
#'   \code{Homer} tool, sequences are weighted to correct for GC and k-mer 
#'   composition differences between fore- and background sets.
#'   
#'   The background sequences are defined according to the value of the
#'   \code{background} argument:
#'   \itemize{
#'     \item{otherBins}{: sequences from all other bins (excluding the current 
#'       bin)}
#'     \item{allBins}{: sequences from all bins (including the current bin)}
#'     \item{zeroBin}{: sequences from the "zero bin", defined by the
#'       \code{maxAbsX} argument of \code{\link[monaLisa]{bin}}. If \code{bins}
#'       does not define a "zero bin", for example because it was created by
#'       \code{bin(..., maxAbsX = NULL)}, selecting this background definition
#'       will abort with an error.}
#'     \item{genome}{: sequences randomly sampled from the genome (or the
#'       intervals defined in \code{genome.regions} if given). For each
#'       foreground sequence, \code{genome.oversample} background sequences
#'       of the same size are sampled (on average). From these, one per
#'       foreground sequence is selected trying to match the G+C composition.
#'       In order to make the sampling deterministic, a seed number needs to be
#'       provided to the \code{RNGseed} parameter in 
#'       \code{\link[BiocParallel]{SerialParam}}
#'       or \code{\link[BiocParallel]{MulticoreParam}} when creating the 
#'       \code{BiocParallelParam} instance in \code{BPPARAM}.}
#'   }
#'
#'   Motif hits are predicted using \code{\link[monaLisa]{findMotifHits}} and
#'   multiple hits per sequence are counted as just one hit (ZOOPS mode). For
#'   each motif, the weights of sequences that have a hit are summed separately
#'   for foreground (\code{sumForegroundWgtWithHits}) and background
#'   (\code{sumBackgroundWgtWithHits}). The total foreground 
#'   (\code{totalWgtForeground}) and background (\code{totalWgtBackground})
#'   sum of sequence weights is also calculated. If a motif has zero 
#'   \code{sumForegroundWgtWithHits} and \code{sumBackgroundWgtWithHits},
#'   then any values (p-values and enrichment) that are calculated using 
#'   these two numbers are set to NA.
#'
#'   Two statistical tests for the calculation of enrichment log p-value are
#'   available: \code{test = "fisher"} (default) to perform Fisher's exact
#'   tests, or \code{test = "binomial"} to perform binomial tests
#'   (default in \code{Homer}), using:
#'   \itemize{
#'     \item{fisher}{: \code{fisher.test(x = tab, alternative =
#'       "greater")}, where \code{tab} is the contingency table with the summed
#'       weights of sequences in foreground or background sets (rows), and with
#'       or without a hit for a particular motif (columns).}
#'     \item{binomial}{: \code{pbinom(q = sumForegroundWgtWithHits - 1, size =
#'       totalWgtForeground, 
#'       prob = sumBackgroundWgtWithHits / totalWgtBackground,
#'       lower.tail = FALSE, log.p = TRUE)}}
#'   }
#'   
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   with motifs in rows and bins in columns, containing seven assays: \itemize{
#'   \item{negLog10P}{: -log10 P values}
#'   \item{negLog10Padj}{: -log10 adjusted P values}
#'   \item{pearsonResid}{: motif enrichments as Pearson residuals}
#'   \item{expForegroundWgtWithHits}{: expected number of foreground 
#'     sequences with motif hits}
#'   \item{log2enr}{: motif enrichments as log2 ratios}
#'   \item{sumForegroundWgtWithHits}{: Sum of foreground sequence weights
#'     in a bin that have motif hits}
#'   \item{sumBackgroundWgtWithHits}{: Sum of background sequence weights
#'     in a bin that have motif hits}
#' }
#' The \code{rowData} of the object contains annotations (name, PFMs, PWMs 
#' and GC fraction) for the motifs, while the \code{colData} slot contains 
#' summary information about the bins. 
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
#' calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
#'                     min.score = 3)
#'
#' @importFrom TFBSTools ID name
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqnames seqlevels
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom stats p.adjust
#'
#' @export
calcBinnedMotifEnrR <- function(seqs,
                                bins = NULL,
                                pwmL = NULL,
                                background = c("otherBins", "allBins", 
                                               "zeroBin", "genome"),
                                test = c("fisher", "binomial"),
                                maxFracN = 0.7,
                                maxKmerSize = 3L,
                                min.score = 10,
                                matchMethod = "matchPWM",
                                GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                             0.45, 0.5, 0.6, 0.7, 0.8),
                                pseudocount.log2enr = 8,
                                p.adjust.method = "BH",
                                genome = NULL,
                                genome.regions = NULL,
                                genome.oversample = 2,
                                BPPARAM = SerialParam(),
                                verbose = FALSE,
                                ...) {

    # checks
    .assertVector(x = seqs, type = "DNAStringSet")
    if (is.null(bins) && identical(background, "genome")) {
        bins <- factor(rep(1, length(seqs)))
    }
    .assertVector(x = bins, type = "factor")
    if (length(seqs) != length(bins)) {
        stop("'seqs' and 'bins' must be of equal length and in the same order")
    }
    .assertVector(x = pwmL, type = "PWMatrixList")
    background <- match.arg(background)
    .assertScalar(x = pseudocount.log2enr, type = "numeric", 
                  rngIncl = c(0, Inf))
    .assertScalar(x = p.adjust.method, type = "character", 
                  validValues = stats::p.adjust.methods)
    if (identical(background, "zeroBin") &&
        (is.null(getZeroBin(bins)) || is.na(getZeroBin(bins)))) {
        stop("For background = 'zeroBin', 'bins' has to define a zero bin ",
             "(see 'maxAbsX' arugment of 'bin' function).")
    }
    if (identical(background, "genome")) {
        if (is.null(genome) || !(is(genome, "DNAStringSet") || 
                                 is(genome, "BSgenome"))) {
            stop("For background = 'genome', 'genome' must be either a ",
                 "DNAStringSet or a BSgenome object.")
        }
        if (!is.null(genome.regions)) {
            if (!is(genome.regions, "GRanges")) {
                stop("For background = 'genome', 'genome.regions' must be ",
                     "either NULL or a GRanges object.")
            }
            if (!all(seqlevels(genome.regions) %in% names(genome))) {
                stop("'genome.regions' contains seqlevels not contained in ", 
                     "'genome'")
            }
        }
        .assertScalar(x = genome.oversample, type = "numeric", 
                      rngIncl = c(1, Inf))
    }
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")
    if (is.null(names(seqs))) {
        names(seqs) <- paste0("s", seq_along(seqs))
    }
    .checkIfSeqsAreEqualLength(x = seqs)


    # filter sequences
    if (verbose) {
        message("Filtering sequences ...")
    }
    keep <- .filterSeqs(seqs, maxFracN = maxFracN, verbose = verbose)
    battr <- attributes(bins) # rescue attributes dropped by subsetting
    bin0 <- getZeroBin(bins)
    bins <- bins[keep]
    attr(bins, "binmode") <- battr$binmode
    attr(bins, "breaks") <- battr$breaks
    if (!is.null(bin0)) {
        bins <- setZeroBin(bins, bin0)
    }
    seqs <- seqs[keep]
    
    # stop if all sequences were filtered out
    if (sum(keep) == 0) {
      stop("No sequence passed the filtering step. Cannot proceed with the ", 
           "enrichment analysis ...")
    }
    
    # scan sequences with motif
    if (verbose) {
        message("Scanning sequences for motif hits...")
    }
    hits <- findMotifHits(query = pwmL, subject = seqs, min.score = min.score,
                          method = matchMethod, BPPARAM = BPPARAM, ...)
    if (isEmpty(hits)) {
        stop("No motif hits found in any of the sequences - aborting.")
    }


    # create motif hit matrix
    if (verbose) {
        message("Create motif hit matrix...")
    }
    hitmatrix <- unclass(table(
        factor(seqnames(hits), levels = seqlevels(hits)),
        factor(hits$pwmid, levels = TFBSTools::ID(pwmL))
    ))
    # zoops (zero or one per sequence) mode
    hitmatrix[hitmatrix > 0] <- 1

    # iterate over bins
    enrichL <- bplapply(structure(seq.int(nlevels(bins)), names = levels(bins)),
                        function(i) {

        if (verbose)
            message("starting analysis of bin ", levels(bins)[i])
        verbose1 <- verbose && bpnworkers(BPPARAM) == 1L

        # define background set and create sequence info data frame
        if (verbose1) {
            message("Defining background sequence set (", background, ")...")
        }
        df <- .defineBackground(sqs = seqs,
                                bns = bins,
                                bg = background,
                                currbn = i,
                                gnm = genome,
                                gnm.regions = genome.regions,
                                gnm.oversample = genome.oversample,
                                maxFracN = maxFracN)

        # for background = "genome", scan sampled background sequences for 
        # motifs
        if (identical(background, "genome")) {
            if (verbose1) {
              message("Scanning genomic background sequences for motif hits...")
            }
            hits.genome <- findMotifHits(
                query = pwmL, subject = df$seqs[!df$isForeground],
                min.score = min.score, method = matchMethod,
                BPPARAM = BPPARAM, ...)
            if (isEmpty(hits.genome)) {
              stop("No motif hits found in any of the genomic background ", 
                   "sequences - aborting.")
            }

            # create motif hit matrix
            hitmatrix.genome <- unclass(
                table(factor(seqnames(hits.genome),
                             levels = seqlevels(hits.genome)),
                      factor(hits.genome$pwmid,
                             levels = TFBSTools::ID(pwmL)))
            )
            # zoops (zero or one per sequence) mode
            hitmatrix.genome[hitmatrix.genome > 0] <- 1
            
            hitmatrix2 <- rbind(hitmatrix, hitmatrix.genome)

        } else {
            hitmatrix2 <- hitmatrix
        }
        
        # calculate initial background sequence weights based on G+C composition
        if (verbose1) {
            message("Correcting for GC differences to the background ", 
                    "sequences...")
        }
        df <- .calculateGCweight(df = df,
                                 GCbreaks = GCbreaks,
                                 verbose = verbose1)
        
        # if df is empty, then all seqs were filtered out in the GC weight 
        # calculation step
        if (nrow(df) == 0) {
          stop("No sequences remained after the GC weight calculation ", 
               "step in bin ", levels(bins)[i], 
               " due to no GC bin containing both fore- and background ", 
               "sequences. Cannot proceed with the enrichment analysis ...")
        }

        # update background sequence weights based on k-mer composition
        if (verbose1) {
            message("Correcting for k-mer differences between fore- and ", 
                    "background sequences...")
        }
        df <- .iterativeNormForKmers(df = df,
                                     maxKmerSize = maxKmerSize,
                                     verbose = verbose1)

        # calculate motif enrichments
        if (verbose1) {
            message("Calculating motif enrichment...")
        }
        enrich1 <- .calcMotifEnrichment(
            motifHitMatrix = hitmatrix2[rownames(df), , drop = FALSE],
            df = df, test = test, verbose = verbose1
        )
        return(enrich1)

    }, BPPARAM = BPPARAM)


    # summarize results as SE
    # ... -log10 of P value
    P <- do.call(cbind, lapply(enrichL, function(enrich1) {
        logpVals <- -enrich1[, "logP"] / log(10)
        names(logpVals) <- enrich1[, "motifName"]
        logpVals
    }))

    # ... -log10 of adjusted P value
    padj <- matrix(-log10(p.adjust(as.vector(10**(-P)), 
                                   method = p.adjust.method)), nrow = nrow(P))
    dimnames(padj) <- dimnames(P)
    padj[which(padj == Inf, arr.ind = TRUE)] <- max(padj[is.finite(padj)])

    # ... Pearson residuals
    enrTF <- do.call(cbind, lapply(enrichL, function(enrich1) {
        enr <- .calcPearsonResiduals(
            matchCountBg = enrich1[, "sumBackgroundWgtWithHits"], 
            totalWeightBg = unique(enrich1[, "totalWgtBackground"]),
            matchCountFg = enrich1[, "sumForegroundWgtWithHits"],
            totalWeightFg = unique(enrich1[, "totalWgtForeground"])
        )
        names(enr) <- enrich1[, "motifName"]
        enr
    }))

    # expected foreground weights
    expFG <- do.call(cbind, lapply(enrichL, function(enrich1) {
        expfg <- .calcExpFg(
            matchCountBg = enrich1[, "sumBackgroundWgtWithHits"], 
            totalWeightBg = unique(enrich1[, "totalWgtBackground"]),
            matchCountFg = enrich1[, "sumForegroundWgtWithHits"],
            totalWeightFg = unique(enrich1[, "totalWgtForeground"])
        )
        names(expfg) <- enrich1[, "motifName"]
        expfg
    }))
    
    # log2 enrichments
    log2enr <- do.call(cbind, lapply(enrichL, function(enrich1) {
        l2e <- .calcLog2Enr(
            matchCountBg = enrich1[, "sumBackgroundWgtWithHits"],
            totalWeightBg = unique(enrich1[, "totalWgtBackground"]),
            matchCountFg = enrich1[, "sumForegroundWgtWithHits"],
            totalWeightFg = unique(enrich1[, "totalWgtForeground"]), 
            pseudocount = pseudocount.log2enr
        )
        names(l2e) <- enrich1[, "motifName"]
        l2e
    }))

    # return SummarizedExperiment
    pfmL <- do.call(TFBSTools::PFMatrixList, lapply(pwmL, function(x) {
        m <- TFBSTools::Matrix(x)
        m <- 2^m * 0.25 # assuming uniform background and no pseudocounts
        PFMatrix(ID = TFBSTools::ID(x),
                 name = TFBSTools::name(x),
                 profileMatrix = m * 100)
    }))
    percentGC <- unlist(lapply(pfmL, function(x) {
        x <- TFBSTools::Matrix(x)
        100 * sum(x[c("C","G"), ]) / sum(x)
    }), use.names = FALSE)
    rdat <- DataFrame(motif.id = TFBSTools::ID(pwmL),
                      motif.name = TFBSTools::name(pwmL),
                      motif.pfm = pfmL,
                      motif.pwm = pwmL,
                      motif.percentGC = percentGC)
    if (is.null(attr(bins, "breaks"))) {
        binL <- binH <- rep(NA, nlevels(bins))
    } else {
        binL <- attr(bins, "breaks")[-(nlevels(bins) + 1)]
        binH <- attr(bins, "breaks")[-1]
    }
    cdat <- DataFrame(
        bin.names = levels(bins),
        bin.lower = binL,
        bin.upper = binH,
        bin.nochange = seq.int(nlevels(bins)) %in% getZeroBin(bins),
        totalWgtForeground = do.call(c, lapply(enrichL, function(x) {
            x$totalWgtForeground[1]
        })), 
        totalWgtBackground = do.call(c, lapply(enrichL, function(x) {
            x$totalWgtBackground[1]
        })))
    mdat <- list(bins = bins,
                 bins.binmode = attr(bins, "binmode"),
                 bins.breaks = as.vector(attr(bins, "breaks")),
                 bins.bin0 = getZeroBin(bins),
                 param = list(method = "R",
                              background = background,
                              test = test,
                              maxFracN = maxFracN,
                              maxKmerSize = maxKmerSize,
                              min.score = min.score,
                              matchMethod = matchMethod,
                              pseudocount.log2enr = pseudocount.log2enr,
                              p.adj.method = p.adjust.method,
                              genome.class = class(genome),
                              genome.regions = genome.regions,
                              genome.oversample = genome.oversample,
                              BPPARAM.class = class(BPPARAM),
                              BPPARAM.bpnworkers = bpnworkers(BPPARAM),
                              verbose = verbose))
    assaySumForegroundWgtWithHits <- 
        do.call(cbind, lapply(enrichL, function(x){x$sumForegroundWgtWithHits}))
    assaySumBackgroundWgtWithHits <- 
        do.call(cbind, lapply(enrichL, function(x){x$sumBackgroundWgtWithHits}))
    
    # ... set motifs with zero fore- and background sums to NA
    assayFgBgSum <- assaySumForegroundWgtWithHits + 
        assaySumBackgroundWgtWithHits
    set_NA <- assayFgBgSum == 0
    P[set_NA] <- NA
    padj[set_NA] <- NA
    enrTF[set_NA] <- NA
    expFG[set_NA] <- NA
    log2enr[set_NA] <- NA

    se <- SummarizedExperiment(
      assays = list(negLog10P = P, 
                    negLog10Padj = padj, 
                    pearsonResid = enrTF,
                    expForegroundWgtWithHits = expFG,
                    log2enr = log2enr, 
                    sumForegroundWgtWithHits = assaySumForegroundWgtWithHits, 
                    sumBackgroundWgtWithHits = assaySumBackgroundWgtWithHits),
      rowData = rdat, colData = cdat, metadata = mdat
    )

    return(se)
}
