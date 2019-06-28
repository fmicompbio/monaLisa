# given a string of A,C,G,T letters, construct a profileMatrix (internal)
.cons2matrix <- function(x, n = 100L) {
    stopifnot(exprs = {
        is.character(x)
        length(x) == 1L
    })
    m <- matrix(0L, nrow = 4, ncol = nchar(x), dimnames = list(c("A","C","G","T"), NULL))
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
#' @param pseudoCount A \code{numeric} scalar - will be added to the observed
#'   counts for each k-mer to avoid zero values.
#'
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @importFrom tidyr %>%
#' @importFrom XVector subseq
#' @importFrom stats ppois
#'
#' @export
getKmerFreq <- function(seqs, kmerLen = 4, MMorder = 2, pseudoCount = 0.5) {
    ## pre-flight checks
    if (is.character(seqs))
        seqs <- DNAStringSet(seqs)
    stopifnot(exprs = {
        is(seqs, "DNAStringSet")
        is.numeric(kmerLen)
        length(kmerLen) == 1L
        round(kmerLen, 0L) == kmerLen
        is.numeric(MMorder)
        length(MMorder) == 1L
        round(MMorder, 0L) == MMorder
        MMorder > 1
        MMorder < kmerLen
    })

    ## observed k-mer frequencies
    kmerFreq <- oligonucleotideFrequency(seqs, width = kmerLen) %>% colSums

    ## expected k-mer frequencies (log2-probabilities with a pseudocount)
    lp_long  <- log2(oligonucleotideFrequency(seqs, width = MMorder) %>%
                     { colSums(.) + pseudoCount } %>%
                     { . / sum(.) })
    lp_short <- log2(oligonucleotideFrequency(seqs, width = MMorder - 1L) %>%
                     { colSums(.) + pseudoCount } %>%
                     { . / sum(.) })
    log2pMM <- sapply(names(kmerFreq), function(current.kmer) {
        n <- nchar(current.kmer) - MMorder + 1L
        ii_long <- substr(rep(current.kmer, n),
                          start = 1:n, stop = 0:(n - 1L) + MMorder)
        ii_short <- substr(rep(current.kmer, n - 1L),
                           start = 2:n, stop = 1:(n - 1L) + MMorder - 1L)
        sum(lp_long[ii_long]) - sum(lp_short[ii_short])
    })
    kmerFreqMM <- (2 ** log2pMM) * sum(kmerFreq)

    ## calculate enrichment statistics
    ## ... log2 (obs/exp)
    lenr <- log2((kmerFreq + pseudoCount) / (kmerFreqMM + pseudoCount))
    ## ... z value (Pearson residuals)
    z <- (kmerFreq - kmerFreqMM) / sqrt(kmerFreqMM)
    ## ... P value
    p <- ppois(q = kmerFreq, lambda = kmerFreqMM, lower.tail = FALSE)
    padj <- p.adjust(p, method = "fdr")

    ## return results
    data.frame(freq.obs=kmerFreq, freq.exp=kmerFreqMM,
               log2enr = lenr, z = z, p = p, FDR = padj)
}


#' @title Run a k-mer enrichment analysis.
#'
#' @description Given a set of sequences and corresponding bins, identify
#'   enriched k-mers (n-grams) in each bin. The sequences can be given either
#'   directly or as genomic coordinates.
#'
#' @param x A \code{character} vector, \code{\link[Biostrings]{DNAStringSet}} or
#'   a \code{\link[GenomicRanges]{GRanges}} object with the sequences to analyze.
#' @param b A vector of the same length as \code{x} that groups its elements
#'   into bins (typically a factor, such as the one returned by
#'   \code{\link{bin}}).
#' @param genomepkg Only used if \code{x} is a \code{GRanges} object: A
#'   \code{character} scalar with the name of a \code{BSgenome} package from
#'   which to extract the sequences.
#' @param kmerLen A \code{numeric} scalar giving the k-mer length.
#' @param MMorder A \code{numeric} scalar giving the order of the Markov model
#'   used to calculate the expected frequencies.
#' @param pseudoCount A \code{numeric} scalar - will be added to the observed
#'   counts for each k-mer to avoid zero values.
#' @param Ncpu Number of parallel threads to use.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @seealso \code{\link{getKmerFreq}} used to calculate k-mer enrichments;
#'   \code{\link[BSgenome]{getSeq,BSgenome-method}} which is used to extract
#'   sequences from \code{genomepkg} if \code{x} is a \code{GRanges} object;
#'   \code{\link[parallel]{mclapply}} that is used for parallelization;
#'   \code{\link{bin}} for binning of regions
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} \code{y}
#'   with \describe{
#'     \item{assays(y)}{containing the four components \code{p}, \code{FDR},
#'        \code{enr} and \code{log2enr}), each a k-mer (rows) by bin (columns)
#'        matrix with raw -log10 P values, -log10 false discovery rates and
#'        k-mer enrichments as Pearson residuals (\code{enr}) and as log2 ratios
#'        (\code{log2enr}).}
#'     \item{rowData(x)}{containing information about the k-mers.}
#'     \item{colData(x)}{containing information about the bins.}
#'     \item{metaData(x)}{containing meta data on the object (e.g. parameter values).}
#'   }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame split
#' @importFrom parallel mclapply
#' @importFrom BSgenome getSeq
#' @importFrom TFBSTools PFMatrix PFMatrixList
#'
#' @export
kmerEnrichments <- function(x, b, genomepkg = NULL, kmerLen = 4, MMorder = 2,
                            pseudoCount = 1, Ncpu = 2L, verbose = TRUE) {
    ## pre-flight checks
    stopifnot(exprs = { is.logical(verbose); length(verbose) == 1L })
    if (is.character(x)) {
        if (!all(grepl("^[ACGTNacgtn]+$", x))) {
            stop("'x' must contain only A, C, G, T or N letters")
        }
        if (verbose) {
            message("converting 'x' to 'DNAStringSet'")
        }
        x <- DNAStringSet(x)

    } else if (is(x, "GRanges"))  {
        if (is.null(genomepkg) || !is.character(genomepkg) ||
            length(genomepkg) != 1L || !require(genomepkg, character.only = TRUE)) {
            stop("'genomepkg' must be a character scalar with the name",
                 "of an installed BSgenome package")
        }
        if (verbose) {
            message("extracting sequences for regions in 'x' from '", genomepkg, "'")
        }
        x <- getSeq(get(genomepkg), x)

    } else if (!is(x, "DNAStringSet")) {
        stop("'x' needs to be either a 'character', 'GRanges' or 'DNAStringSet'")
    }
    if (!is.factor(b)) {
        if (verbose) {
            message("converting 'b' to a 'factor'")
        }
        b <- factor(b, levels=unique(b))
    }
    stopifnot(exprs = {
        length(x) == length(b)
        is.numeric(Ncpu)
        length(Ncpu) == 1
        Ncpu > 0
    })

    ## identify enriched k-mers in each bin
    if (verbose) {
        message("searching for enriched ", kmerLen, "-mers in ", nlevels(b),
                " bins using ", Ncpu, if (Ncpu > 1) " cores..." else " core...",
                appendLF = FALSE)
    }
    resL <- parallel::mclapply(split(x, b)[levels(b)], getKmerFreq, kmerLen = kmerLen,
                               MMorder = MMorder, pseudoCount = pseudoCount, mc.cores = Ncpu)
    if (verbose) {
        message("done")
    }

    ## ... create SummarizedExperiment
    brks <- attr(b, "breaks")
    if (is.null(brks)) {
        brks <- rep(NA, nlevels(b) + 1L)
    }
    cdat <- S4Vectors::DataFrame(bin.names = levels(b),
                                 bin.lower = brks[-(nlevels(b)+1)],
                                 bin.upper = brks[-1],
                                 bin.nochange = seq.int(nlevels(b)) %in% attr(b, "bin0"))
    kmers <- rownames(resL[[1]])
    pfms <- do.call(TFBSTools::PFMatrixList, lapply(kmers, function(kmer) {
        TFBSTools::PFMatrix(ID = kmer, name = kmer,
                            profileMatrix = .cons2matrix(kmer))
    }))
    percentGC <- unlist(lapply(pfms, function(x) {
        m <- TFBSTools::Matrix(x)
        100 * sum(m[c("C","G"), ]) / sum(m)
    }), use.names = FALSE)
    rdat <- S4Vectors::DataFrame(motif.name = rownames(resL[[1]]),
                                 motif.pfm = pfms,
                                 motif.percentGC = percentGC)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(p = -log10(do.call(cbind, lapply(resL, "[[", "p"))),
                      FDR = -log10(do.call(cbind, lapply(resL, "[[", "FDR"))),
                      enr = do.call(cbind, lapply(resL, "[[", "z")),
                      log2enr = do.call(cbind, lapply(resL, "[[", "log2enr"))),
        colData = cdat, rowData = rdat,
        metadata = list(sequences = x,
                        bins = b,
                        bins.binmode = attr(b, "binmode"),
                        bins.breaks = as.vector(attr(b, "breaks")),
                        bins.bin0 = attr(b, "bin0"),
                        param.genomepkg = genomepkg,
                        param.kmerLen = kmerLen,
                        param.MMorder = MMorder,
                        param.pseudoCount = pseudoCount,
                        param.Ncpu = Ncpu,
                        motif.distances = NULL)
    )
    rownames(se) <- rownames(resL[[1]])
    return(se)
}
