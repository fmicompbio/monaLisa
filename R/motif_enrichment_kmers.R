#' @title Create matrix from consensus sequence
#'
#' @description Given a nucleotide sequence of A,C,G,T letter corresponding
#'   to a motif's consensus string, construct a positional frequency matrix.
#'   This matrix can for example be used as the \code{profileMatrix} argument
#'   in the constructor for a \code{TFBSTools::PFMatrix} object.
#'
#' @return A positional frequency matrix.
#' 
#' @param x Character scalar with the motif the consensus sequence.
#' @param n Integer scalar giving the columns sums in the constructed
#'   matrix (number of observed bases at each position).
#' 
#' @keywords internal
.cons2matrix <- function(x, n = 100L) {
    .assertScalar(x = x, type = "character")
    .assertScalar(x = n, type = "integer", rngIncl = c(1L, Inf))
    m <- matrix(0L, nrow = 4, ncol = nchar(x),
                dimnames = list(c("A","C","G","T"), NULL))
    xx <- strsplit(x, "", fixed = TRUE)[[1]]
    ok <- which(xx %in% rownames(m))
    m[cbind(match(xx[ok], rownames(m)), ok)] <- n
    m
}


#' @title Calculate observed and expected k-mer frequencies
#'
#' @description Given a set of sequences, calculate observed and expected k-mer
#'   frequencies. Expected frequencies are based on a Markov model of order
#'   \code{MMorder}.
#'
#' @param seqs Set of sequences, either a \code{character} vector or a
#'   \code{\link{DNAStringSet}}.
#' @param kmerLen A \code{numeric} scalar giving the k-mer length.
#' @param MMorder A \code{numeric} scalar giving the order of the Markov model
#'   used to calculate the expected frequencies.
#' @param pseudocount A \code{numeric} scalar - will be added to the observed
#'   counts for each k-mer to avoid zero values.
#' @param zoops A \code{logical} scalar. If \code{TRUE} (the default), only one
#'   or zero occurences of a k-mer are considered per sequence.
#' @param strata A \code{factor} or a \code{numeric} scalar defining the strata
#'   of sequences. A separate Markov model and expected k-mer frequencies are
#'   estimated for the set of sequences in each stratum (level in a \code{strata}
#'   factor). If \code{strata} is a scalar value, it will be interpreted as the
#'   number of strata to split the sequences into according to their CpG
#'   observed-over-expected counts using \code{kmeans(CpGoe, centers = strata)}.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#' @param includeRevComp A \code{logical} scalar. If \code{TRUE} (default),
#'   count k-mer occurrences in both \code{seqs} and their reverse-complement,
#'   by concatenating \code{seqs} and their reverse-complemented versions
#'   before the counting. This is useful if motifs can be expected to occur
#'   on any strand (e.g. DNA sequences of ChIP-seq peaks). If motifs are only
#'   expected on the forward strand (e.g. RNA sequences of CLIP-seq peaks),
#'   \code{includeRevComp = FALSE} should be used. Note that if \code{strata}
#'   is a vector of the same length as \code{seqs}, each reverse-complemented
#'   sequence will be assigned to the same stratum as the forward sequence.
#'
#' @return A \code{list} with observed and expected k-mer frequencies (\code{freq.obs}
#'   and \code{freq.exp}, respectively), and enrichment statistics for each k-mer.
#'
#' @examples 
#' res <- getKmerFreq(seqs = c("AAAAATT", "AAATTTT"), kmerLen = 3)
#' names(res)
#' head(res$freq.obs)
#' head(res$freq.exp)
#' 
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency reverseComplement
#' @importFrom XVector subseq
#' @importFrom stats ppois kmeans
#'
#' @export
getKmerFreq <- function(seqs,
                        kmerLen = 5,
                        MMorder = 1,
                        pseudocount = 1,
                        zoops = TRUE,
                        strata = rep(1L, length(seqs)),
                        p.adjust.method = "BH",
                        includeRevComp = TRUE) {

    ##comments

    ## When plotting the log2 enrichments or something similar against observed frequencies,
    ## there will be a correlation as several kmers have the same expected frequency and thus their
    ## enrichment and observed frequencies will naturally correlate. It makes thus more sense to plot
    ## the log2 enrichment against expected frequencies

    ## Generally, higher MMorder will take care of CpG bias and similar effects. Zoops will take care of
    ## higher-order repeats.

    ## pre-flight checks
    if (is.character(seqs))
        seqs <- DNAStringSet(seqs)
    .assertScalar(x = kmerLen, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = MMorder, type = "numeric", rngExcl = c(0, kmerLen - 1L))
    .assertScalar(x = pseudocount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = zoops, type = "logical")
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)
    .assertScalar(x = includeRevComp, type = "logical")
    stopifnot(exprs = {
        is(seqs, "DNAStringSet")
        round(kmerLen, 0L) == kmerLen
        round(MMorder, 0L) == MMorder
        length(strata) == length(seqs) || (is.numeric(strata) && length(strata) == 1L)
    })
    
    ## include reverse-complement sequences?
    if (identical(includeRevComp, TRUE)) {
        seqsrc <- Biostrings::reverseComplement(seqs)
        names(seqsrc) <- paste0(names(seqs), "_rc")
        
        if (length(strata) == length(seqs)) {
            # recycle strata (same stratum for forward and rev-comp sequence)
            strata <- rep(strata, 2L)
        }
        seqs <- c(seqs, seqsrc)
    }
    
    ## split sequences into strata
    CpGoe <- NA
    if (identical(length(strata), 1L)) {
        n1 <- oligonucleotideFrequency(seqs, width = 1L)
        n2 <- oligonucleotideFrequency(seqs, width = 2L)
        CpGoe <- n2[, "CG"] / (n1[,"C"] / rowSums(n1) * n1[,"G"])
        strata <- kmeans(x = CpGoe, centers = strata)$cluster
    }

    ## for each sequence stratum, calculate...
    res.strata <- lapply(split(seqs, strata), function(seqs.stratum) {
        ## ... observed k-mer frequencies
        kmerFreqRaw.stratum <- oligonucleotideFrequency(seqs.stratum, width = kmerLen)
        if (zoops) {
            kmerFreq.stratum <- colSums(kmerFreqRaw.stratum > 0)
            #make new sequences that are simply the kmers detected
            seqs.stratum.zoops <- DNAStringSet(rep(names(kmerFreq.stratum), kmerFreq.stratum))
            lp_long  <- colSums(oligonucleotideFrequency(seqs.stratum.zoops, width = MMorder + 1L)) + pseudocount
            lp_short <- colSums(oligonucleotideFrequency(seqs.stratum.zoops, width = MMorder)     ) + pseudocount
        } else {
            kmerFreq.stratum <- colSums(kmerFreqRaw.stratum)
            lp_long  <- colSums(oligonucleotideFrequency(seqs.stratum, width = MMorder + 1L)) + pseudocount
            lp_short <- colSums(oligonucleotideFrequency(seqs.stratum, width = MMorder))      + pseudocount
        }
        lp_long  <- log2(lp_long / sum(lp_long))
        lp_short <- log2(lp_short / sum(lp_short))

        ## ... expected k-mer frequencies (log2-probabilities with a pseudocount)
        n <- nchar(names(kmerFreq.stratum)[1]) - MMorder
        log2pMM <- vapply(names(kmerFreq.stratum), function(current.kmer) {
            ii_long <- substr(rep(current.kmer, n),
                              start = seq_len(n), stop = seq_len(n) + MMorder)
            ii_short <- substr(rep(current.kmer, n - 1L),
                               start = seq(2, n), stop = seq_len(n - 1L) + MMorder)
            sum(lp_long[ii_long]) - sum(lp_short[ii_short])
        }, numeric(1))
        kmerFreqMM.stratum <- (2 ** log2pMM) * sum(kmerFreq.stratum)
        cbind(obs = kmerFreq.stratum, exp = kmerFreqMM.stratum)
    })

    ## sum frequencies over strata
    tmp <- Reduce("+", res.strata)
    kmerFreq   <- tmp[, "obs"]
    kmerFreqMM <- tmp[, "exp"]

    ## calculate enrichment statistics
    ## ... log2 (obs/exp)
    lenr <- log2((kmerFreq + pseudocount) / (kmerFreqMM + pseudocount))
    ## ... z value (Pearson residuals)
    z <- (kmerFreq - kmerFreqMM) / sqrt(kmerFreqMM)
    ## ... use square-root as variance-stablizing function
    sDelta <- sqrt(kmerFreq) - sqrt(kmerFreqMM)
    ## ... P value
    p <- ppois(q = kmerFreq, lambda = kmerFreqMM, lower.tail = FALSE)
    padj <- p.adjust(p, method = p.adjust.method)

    ## return results
    list(freq.obs = kmerFreq, freq.exp = kmerFreqMM,
         log2enr = lenr, sqrtDelta = sDelta, z = z, p = p, padj = padj,
         strata = strata, freq.strata = res.strata, CpGoe = CpGoe)
}


#' @title Calculate k-mer enrichment
#'
#' @description Given sequences, foreground/background labels and
#'   weights, calculate the enrichment of each k-mer
#'   in foreground compared to background. This function is called by
#'   \code{calcBinnedKmerEnr()} for each bin if \code{background != "model"}.
#'
#'   The default type of test is \code{"fisher"}.
#'   Alternatively, a binomial test can be used by \code{test = "binomial"}.
#'   Using Fisher's exact test has the advantage that special cases such as
#'   zero background counts are handled without ad-hoc adjustments to the
#'   k-mer frequencies.
#'
#'   For \code{test = "fisher"}, \code{fisher.test} is used with
#'   \code{alternative = "greater"}, making it a one-sided test for enrichment,
#'   as is the case with the binomial test.
#'
#' @param k Numeric scalar giving the length of k-mers to analyze.
#' @param df a \code{DataFrame} with sequence information as returned by
#'   \code{.iterativeNormForKmers()}.
#' @param test type of motif enrichment test to perform.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#' 
#' @details The function works in ZOOPS mode, which means only one
#'   or zero occurrences of a k-mer are considered per sequence. This is helpful
#'   to reduce the impact of simple sequence repeats occurring in few sequences.
#'
#' @return A \code{data.frame} containing the motifs as rows and the columns:
#'   \itemize{
#'     \item{motifName}{: the motif name}
#'     \item{logP}{: the log p-value for enrichment (natural logarithm).
#'        If \code{test="binomial"} (default), this log p-value is identical to
#'        the one returned by Homer.}
#'     \item{sumForegroundWgtWithHits}{: the weighted number of k-mer hits in
#'        foreground sequences.}
#'     \item{sumBackgroundWgtWithHits}{: the weighted number of k-mer hits in
#'        background sequences.}
#'     \item{totalWgtForeground}{: the total sum of weights of foreground
#'        sequences.}
#'     \item{totalWgtBackground}{: the total sum of weights of background
#'        sequences.}
#'   }
#'
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom stats pbinom fisher.test
#' 
#' @keywords internal
.calcKmerEnrichment <- function(k,
                                df,
                                test = c("fisher", "binomial"),
                                verbose = FALSE){
    
    # checks
    .assertScalar(x = k, type = "numeric", rngIncl = c(1, Inf))
    .checkDfValidity(df)
    test <- match.arg(test)
    .assertScalar(x = verbose, type = "logical")
    
    totalWgtForeground <- sum(df$seqWgt[df$isForeground])
    totalWgtBackground <- sum(df$seqWgt[!df$isForeground])
    
    kmerHitMatrix <- oligonucleotideFrequency(x = df$seqs,
                                               width = k,
                                               as.prob = FALSE,
                                               with.labels = TRUE)
    
    kmerHitMatrix[kmerHitMatrix > 0] <- 1 # ZOOPS
    kmerHitMatrixWeighted <- kmerHitMatrix * df$seqWgt
    KmatchedSeqCountForeground <- colSums(kmerHitMatrixWeighted[df$isForeground, ])
    KmatchedSeqCountBackground <- colSums(kmerHitMatrixWeighted[!df$isForeground, ])
    
    # calculate k-mer enrichment
    if (identical(test, "binomial")) {
        logP <- .binomEnrichmentTest(matchCountBg = KmatchedSeqCountBackground,
                                     totalWeightBg = totalWgtBackground,
                                     matchCountFg = KmatchedSeqCountForeground,
                                     totalWeightFg = totalWgtForeground,
                                     verbose = verbose)
    } else if (identical(test, "fisher")) {
        logP <- .fisherEnrichmentTest(matchCountBg = KmatchedSeqCountBackground, 
                                      totalWeightBg = totalWgtBackground,
                                      matchCountFg = KmatchedSeqCountForeground,
                                      totalWeightFg = totalWgtForeground,
                                      verbose = verbose)
    }
    
    return(data.frame(motifName = names(logP),
                      logP = logP,
                      sumForegroundWgtWithHits = KmatchedSeqCountForeground,
                      sumBackgroundWgtWithHits = KmatchedSeqCountBackground,
                      totalWgtForeground = totalWgtForeground,
                      totalWgtBackground = totalWgtBackground))
}


#' @title Calculate k-mer enrichment in bins of sequences.
#'
#' @description Given a set of sequences and corresponding bins, identify
#'   enriched k-mers (n-grams) in each bin. The sequences can be given either
#'   directly or as genomic coordinates.
#'
#' @param seqs \code{\link[Biostrings]{DNAStringSet}} object with sequences to test
#' @param bins factor of the same length and order as \code{seqs}, indicating
#'   the bin for each sequence. Typically the return value of
#'   \code{\link[monaLisa]{bin}}. For \code{background = "genome"} or
#'   \code{background = "model"}, \code{bins} can be omitted.
#' @param kmerLen A \code{numeric} scalar giving the k-mer length.
#' @param background A \code{character} scalar specifying the background sequences
#'   to use. One of \code{"otherBins"} (default), \code{"allBins"}, \code{"zeroBin"},
#'   \code{"genome"} or \code{"model"} (see "Details").
#' @param MMorder A \code{numeric} scalar giving the order of the Markov model
#'   used to calculate the expected frequencies for \code{background = "model"}.
#' @param test A \code{character} scalar specifying the type of enrichment test
#'   to perform. One of \code{"fisher"} (default) or \code{"binomial"}. The
#'   enrichment test is one-sided (enriched in foreground).
#' @param includeRevComp A \code{logical} scalar. If \code{TRUE} (default),
#'   count k-mer occurrences in both \code{seqs} and their reverse-complement,
#'   by concatenating \code{seqs} and their reverse-complemented versions
#'   before the counting. This is useful if motifs can be expected to occur
#'   on any strand (e.g. DNA sequences of ChIP-seq peaks). If motifs are only
#'   expected on the forward strand (e.g. RNA sequences of CLIP-seq peaks),
#'   \code{includeRevComp = FALSE} should be used. Note that \code{bins}
#'   will be recycled for the reverse complemented sequences, which means that
#'   each reverse-complemented sequence will be assigned to the same bib as the
#'   corresponding forward sequence.
#' @param maxFracN A numeric scalar with the maximal fraction of N bases allowed
#'   in a sequence (defaults to 0.7). Sequences with higher fractions are
#'   excluded from the analysis.
#' @param maxKmerSize the maximum k-mer size to consider, when adjusting
#'   background sequence weights for k-mer composition compared to the
#'   foreground sequences. The default value (3) will correct for mono-, di-
#'   and tri-mer composition.
#' @param GCbreaks The breaks between GC bins. The default value is based on
#'   the hard-coded bins used in Homer.
#' @param pseudocount.kmers A \code{numeric} scalar - will be added to the
#'   observed and expected counts for each k-mer to avoid zero values.
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
#' @param verbose A \code{logical} scalar. If \code{TRUE}, report on progress.
#'
#' @details This function implements a binned k-mer enrichment analysis. In each
#'   enrichment analysis, the sequences in a specific bin are used as foreground
#'   sequences to test for k-mer enrichments comparing to background sequences
#'   (defined by \code{background}, see below), similarly as in done for motifs
#'   in \code{\link{calcBinnedMotifEnrR}}. Sequences are weighted to correct for
#'   GC and shorter k-mer composition differences between fore- and background
#'   sets.
#'   
#'   The background sequences are defined according to the value of the
#'   \code{background} argument:
#'   \itemize{
#'     \item{otherBins}{: sequences from all other bins (excluding the current bin)}
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
#'       In order to make the sampling deterministic, the random number
#'       generator has to be seeded using \code{set.seed} before calling 
#'       this function.}
#'     \item{model}{: a Markov model of the order \code{MMorder} is estimated
#'       from the foreground sequences and used to estimate expected k-mer
#'       frequencies. K-mer enrichments are then calculated comparing observed
#'       to these expected frequencies.}
#'   }
#'
#'   For each k-mer, the weights of sequences is multiplied with the number
#'   of k-mer occurrences in each sequence and summed, separately for foreground
#'   (\code{sumForegroundWgtWithHits}) and background
#'   (\code{sumBackgroundWgtWithHits}) sequences. The function works in ZOOPS
#'   (Zero-Or-One-Per-Sequence) mode, so at most one occurrence per
#'   sequence is counted, which helps reduce the impact of sequence repeats.
#'   The total foreground (\code{totalWgtForeground}) and background
#'   (\code{totalWgtBackground}) sum of sequence weights is also calculated. If
#'   a k-mer has zero \code{sumForegroundWgtWithHits} and
#'   \code{sumBackgroundWgtWithHits}, then any values (p-values and enrichment)
#'   that are calculated using these two numbers are set to NA.
#'
#'   Two statistical tests for the calculation of enrichment log p-value are
#'   available: \code{test = "fisher"} (default) to perform Fisher's exact
#'   tests, or \code{test = "binomial"} to perform binomial tests, using:
#'   \itemize{
#'     \item{fisher}{: \code{fisher.test(x = tab, alternative =
#'       "greater")}, where \code{tab} is the contingency table with the summed
#'       weights of sequences in foreground or background sets (rows), and with
#'       or without a occurrences of a particular k-mer (columns).}
#'     \item{binomial}{: \code{pbinom(q = sumForegroundWgtWithHits - 1, size =
#'       totalWgtForeground, prob = sumBackgroundWgtWithHits / totalWgtBackground,
#'       lower.tail = FALSE, log.p = TRUE)}}
#'   }
#'   
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'   with motifs in rows and bins in columns, containing seven assays: \itemize{
#'   \item{negLog10P}{: -log10 P values}
#'   \item{negLog10Padj}{: -log10 adjusted P values}
#'   \item{pearsonResid}{: k-mer enrichments as Pearson residuals}
#'   \item{expForegroundWgtWithHits}{: expected number of foreground 
#'     sequences with motif hits}
#'   \item{log2enr}{: k-mer enrichments as log2 ratios}
#'   \item{sumForegroundWgtWithHits}{: Sum of foreground sequence weights
#'     in a bin that have k-mer occurrences}
#'   \item{sumBackgroundWgtWithHits}{: Sum of background sequence weights
#'     in a bin that have k-mer occurrences}
#' }
#' #' The \code{rowData} of the object contains annotations (name, PFMs, PWMs 
#' and GC fraction) for the k-mers, while the \code{colData} slot contains 
#' summary information about the bins. 
#' 
#' @examples 
#' seqs <- Biostrings::DNAStringSet(c("GCATGCATGC", "CATGCGCATG"))
#' bins <- factor(1:2)
#' calcBinnedKmerEnr(seqs = seqs, bins = bins, kmerLen = 3)
#' 
#' @seealso \code{\link{getKmerFreq}} used to calculate k-mer enrichments;
#'   \code{\link[BSgenome]{getSeq,BSgenome-method}} which is used to extract
#'   sequences from \code{genomepkg} if \code{x} is a \code{GRanges} object;
#'   \code{\link[BiocParallel]{bplapply}} that is used for parallelization;
#'   \code{\link{bin}} for binning of regions
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom Biostrings reverseComplement
#' @importFrom S4Vectors DataFrame split
#' @importFrom TFBSTools PFMatrix PFMatrixList ID name toPWM
#' @importFrom stats ppois p.adjust
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#'
#' @export
calcBinnedKmerEnr <- function(seqs,
                              bins = NULL,
                              kmerLen = 5,
                              background = c("otherBins", "allBins", "zeroBin", "genome", "model"),
                              MMorder = 1,
                              test = c("fisher", "binomial"),
                              includeRevComp = TRUE,
                              maxFracN = 0.7,
                              maxKmerSize = 3L,
                              GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                           0.45, 0.5, 0.6, 0.7, 0.8),
                              pseudocount.kmers = 1,
                              pseudocount.log2enr = 8,
                              p.adjust.method = "BH",
                              genome = NULL,
                              genome.regions = NULL,
                              genome.oversample = 2,
                              BPPARAM = SerialParam(),
                              verbose = FALSE) {
    ## pre-flight checks
    .assertVector(x = seqs, type = "DNAStringSet")
    if (is.null(bins) && (background %in% c("genome", "model"))) {
        bins <- factor(rep(1, length(seqs)))
    }
    .assertVector(x = bins, type = "factor")
    if (length(seqs) != length(bins)) {
        stop("'seqs' and 'bins' must be of equal length and in the same order")
    }
    .assertScalar(x = kmerLen, type = "numeric", rngIncl = c(1, Inf))
    background <- match.arg(background)
    test <- match.arg(test)
    .assertScalar(x = includeRevComp, type = "logical")
    .assertScalar(x = pseudocount.kmers, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pseudocount.log2enr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)
    if (identical(background, "zeroBin") &&
        (is.null(getZeroBin(bins)) || is.na(getZeroBin(bins)))) {
        stop("For background = 'zeroBin', 'bins' has to define a zero bin (see ",
             "'maxAbsX' arugment of 'bin' function).")
    }
    if (identical(background, "genome")) {
        if (is.null(genome) || !(is(genome, "DNAStringSet") || is(genome, "BSgenome"))) {
            stop("For background = 'genome', 'genome' must be either a ",
                 "DNAStringSet or a BSgenome object.")
        }
        if (!is.null(genome.regions)) {
            if (!is(genome.regions, "GRanges")) {
                stop("For background = 'genome', 'genome.regions' must be either ",
                     "NULL or a GRanges object.")
            }
            if (!all(seqlevels(genome.regions) %in% names(genome))) {
                stop("'genome.regions' contains seqlevels not contained in 'genome'")
            }
        }
        .assertScalar(x = genome.oversample, type = "numeric", rngIncl = c(1, Inf))
    }
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")
    if (is.null(names(seqs))) {
        names(seqs) <- paste0("s", seq_along(seqs))
    }
    .checkIfSeqsAreEqualLength(x = seqs)

    
    ## include reverse-complement sequences?
    if (identical(includeRevComp, TRUE)) {
        seqsrc <- Biostrings::reverseComplement(seqs)
        names(seqsrc) <- paste0(names(seqs), "_rc")
        
        battr <- attributes(bins) # rescue attributes dropped by c()
        bin0 <- getZeroBin(bins)
        bins <- c(bins, bins)
        attr(bins, "binmode") <- battr$binmode
        attr(bins, "breaks") <- battr$breaks
        if (!is.null(bin0)) {
            bins <- setZeroBin(bins, bin0)
        }
        seqs <- c(seqs, seqsrc)
    }


    ## filter sequences
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
        stop("No sequence passed the filtering step. Cannot proceed with the enrichment analysis ...")
    }
    
    
    # iterate over bins
    enrichL <- bplapply(structure(seq.int(nlevels(bins)), names = levels(bins)),
                        function(i) {

        if (verbose)
            message("starting analysis of bin ", levels(bins)[i])
        verbose1 <- verbose && bpnworkers(BPPARAM) == 1L

        if (identical(background, "model")) {
            # no need to calculate sequence weights for background = "model"
            is.fg <- as.numeric(bins) == i
            Nfg <- sum(is.fg)
            res1 <- getKmerFreq(seqs = seqs[is.fg],
                                kmerLen = kmerLen,
                                MMorder = MMorder,
                                pseudocount = pseudocount.kmers,
                                zoops = TRUE,
                                includeRevComp = FALSE)

            if (identical(test, "binomial")) {
                logP <- .binomEnrichmentTest(matchCountBg = res1$freq.exp, 
                                             totalWeightBg = Nfg,
                                             matchCountFg = res1$freq.obs, 
                                             totalWeightFg = Nfg, 
                                             verbose = FALSE)
            } else if (identical(test, "fisher")) {
                logP <- .fisherEnrichmentTest(matchCountBg = res1$freq.exp,
                                              totalWeightBg = Nfg,
                                              matchCountFg = res1$freq.obs,
                                              totalWeightFg = Nfg,
                                              verbose = FALSE)
            }

            return(data.frame(motifName = names(logP),
                              logP = logP,
                              sumForegroundWgtWithHits = res1$freq.obs,
                              sumBackgroundWgtWithHits = res1$freq.exp,
                              totalWgtForeground = Nfg,
                              totalWgtBackground = Nfg))
            
        } else {
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

            # calculate initial background sequence weights based on G+C composition
            if (verbose1) {
                message("Correcting for GC differences to the background sequences...")
            }
            df <- .calculateGCweight(df = df,
                                     GCbreaks = GCbreaks,
                                     verbose = verbose1)

            # if df is empty, then all seqs were filtered out in the GC weight calculation step
            if (nrow(df) == 0) {
                stop("No sequences remained after the GC weight calculation step in bin ", levels(bins)[i], 
                     " due to no GC bin containing both fore- and background sequences. ", 
                     "Cannot proceed with the enrichment analysis ...")
            }

            # update background sequence weights based on k-mer composition
            if (verbose1) {
                message("Correcting for k-mer differences between fore- and background sequences...")
            }
            df <- .iterativeNormForKmers(df = df,
                                         maxKmerSize = maxKmerSize,
                                         verbose = verbose1)

            # calculate motif enrichments
            if (verbose1) {
                message("Calculating ", kmerLen, "-mer enrichment...")
            }
            enrich1 <- .calcKmerEnrichment(k = kmerLen,
                                           df = df,
                                           test = test,
                                           verbose = verbose1)
            return(enrich1)
        }
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
    kmers <- enrichL[[1]]$motifName
    pfmL <- do.call(TFBSTools::PFMatrixList, lapply(kmers, function(kmer) {
        TFBSTools::PFMatrix(ID = kmer, name = kmer,
                            profileMatrix = .cons2matrix(kmer))
    }))
    pwmL <- TFBSTools::toPWM(pfmL)
    percentGC <- unlist(lapply(pfmL, function(x) {
        m <- TFBSTools::Matrix(x)
        100 * sum(m[c("C","G"), ]) / sum(m)
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
    cdat <- DataFrame(bin.names = levels(bins),
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
                 param = list(kmerLen = kmerLen,
                              background = background,
                              MMorder = MMorder,
                              test = test,
                              includeRevComp = includeRevComp,
                              maxFracN = maxFracN,
                              maxKmerSize = maxKmerSize,
                              pseudocount.kmers = pseudocount.kmers,
                              pseudocount.log2enr = pseudocount.log2enr,
                              zoops = TRUE,
                              p.adj.method = p.adjust.method,
                              genome.class = class(genome),
                              genome.regions = genome.regions,
                              genome.oversample = genome.oversample,
                              BPPARAM.class = class(BPPARAM),
                              BPPARAM.bpnworkers = bpnworkers(BPPARAM),
                              verbose = verbose))
    assaySumForegroundWgtWithHits <- do.call(cbind, lapply(enrichL, function(x){x$sumForegroundWgtWithHits}))
    assaySumBackgroundWgtWithHits <- do.call(cbind, lapply(enrichL, function(x){x$sumBackgroundWgtWithHits}))
    
    # ... set motifs with zero fore- and background sums to NA
    assayFgBgSum <- assaySumForegroundWgtWithHits + assaySumBackgroundWgtWithHits
    set_NA <- assayFgBgSum == 0
    P[set_NA] <- NA
    padj[set_NA] <- NA
    enrTF[set_NA] <- NA
    log2enr[set_NA] <- NA
    expFG[set_NA] <- NA
    
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
    rownames(se) <- kmers
    
    return(se)
}
