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
#' @param zoops A \code{logical} scalar. If \code{TRUE} (the default), only one
#'   or zero occurences of a k-mer are considered per sequence.
#'
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @importFrom XVector subseq
#' @importFrom stats ppois
#'
#' @export
getKmerFreq <- function(seqs, kmerLen = 5, MMorder = 1, pseudoCount = 1, zoops = TRUE) {
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
        MMorder > 0
        MMorder < kmerLen - 1L
        is.logical(zoops)
        length(zoops) == 1L
    })

    ## observed k-mer frequencies
    kmerFreqRaw <- oligonucleotideFrequency(seqs, width = kmerLen)
    if (zoops) {
        kmerFreq <- colSums(kmerFreqRaw > 0)
        #lp_long  <- colSums(oligonucleotideFrequency(seqs, width = MMorder + 1L) > 0) + pseudoCount
        #lp_short <- colSums(oligonucleotideFrequency(seqs, width = MMorder)      > 0) + pseudoCount
    } else {
        kmerFreq <- colSums(kmerFreqRaw)
        #lp_long  <- colSums(oligonucleotideFrequency(seqs, width = MMorder + 1L)) + pseudoCount
        #lp_short <- colSums(oligonucleotideFrequency(seqs, width = MMorder))      + pseudoCount
    }
    lp_long  <- colSums(oligonucleotideFrequency(seqs, width = MMorder + 1L)) + pseudoCount
    lp_short <- colSums(oligonucleotideFrequency(seqs, width = MMorder))      + pseudoCount
    lp_long  <- log2(lp_long / sum(lp_long))
    lp_short <- log2(lp_short / sum(lp_short))

    ## expected k-mer frequencies (log2-probabilities with a pseudocount)
    n <- nchar(names(kmerFreq)[1]) - MMorder
    log2pMM <- sapply(names(kmerFreq), function(current.kmer) {
        ii_long <- substr(rep(current.kmer, n),
                          start = 1:n, stop = 1:n + MMorder)
        ii_short <- substr(rep(current.kmer, n - 1L),
                           start = 2:n, stop = 1:(n - 1L) + MMorder)
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
kmerEnrichments <- function(x, b, genomepkg = NULL, kmerLen = 5, MMorder = 1,
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

#' @title Transform k-mer enrichments to motif enrichments.
#'
#' @description Using a set of know motifs and the result of a k-mer enrichment
#'   analysis (k-mers by bins), transform the results into motif enrichments
#'   (motifs by bins).
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} with the
#'   results of a k-mer enrichment analysis (typically generated by
#'   \code{\link{kmerEnrichments}}).
#' @param m Either a \code{\link[TFBSTools]{PFMatrixList}}, or a character
#'   scalar with a file containing motifs in HOMER format (loaded into a
#'   \code{\link[TFBSTools]{PFMatrixList}} by
#'   \code{\link{homerToPFMatrixList}}).
#' @param Ncpu The number of CPU cores to use when calculating similarities
#'   between motifs and k-mers. This uses \code{\link[parallel]{mclapply}}.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @seealso \code{\link{kmerEnrichments}} for performing a k-mer enrichment
#'   analysis, \code{\link[parallel]{mclapply}} for how parallelization is done.
#'
#' @importFrom SummarizedExperiment assayNames assay SummarizedExperiment
#'   colData metadata
#' @importFrom S4Vectors DataFrame
#'
#' @export
convertKmersToMotifs <- function(x, m, Ncpu = 1L, verbose = TRUE) {
    ## pre-flight checks
    stopifnot(is.logical(verbose) && length(verbose) == 1L)
    if (is.character(m) && length(m) == 1L && file.exists(m)) {
        if (verbose) {
            message("reading motifs from ", basename(m))
        }
        m <- homerToPFMatrixList(m)
    }
    stopifnot(exprs = {
        is(x, "SummarizedExperiment")
        "param.kmerLen" %in% names(metadata(x))
        nrow(x) == 4^metadata(x)$param.kmerLen
        all(c("enr","log2enr") %in% assayNames(x))
        is(m, "PFMatrixList")
        is.numeric(Ncpu)
        length(Ncpu) == 1
        Ncpu > 0
    })

    ## calculate motif-by-kmer matrix (probabilities)
    m.k <- motifKmerSimilarity(m, kmerLen = metadata(x)$param.kmerLen,
                               Ncpu = Ncpu, verbose = verbose)

    ## calculate motif-by-bin = motif-by-kmer %*% kmer-by-bin
    m.b.enr <- m.k %*% assay(x, "enr")
    m.b.log2enr <- m.k %*% assay(x, "log2enr")

    ## create and return new SummarizedExperiment object
    percentGC <- unlist(lapply(m, function(x) {
        m <- TFBSTools::Matrix(x)
        100 * sum(m[c("C","G"), ]) / sum(m)
    }), use.names = FALSE)
    rdat <- DataFrame(motif.name = name(m),
                      motif.pfm = m,
                      motif.percentGC = percentGC)
    mdat <- metadata(x)
    mdat[["param.Ncpu"]] <- Ncpu
    mzero <- matrix(0, nrow = nrow(m.b.enr), ncol = ncol(m.b.enr),
                    dimnames = dimnames(m.b.enr))
    se <- SummarizedExperiment(
        assays = list(p = mzero,
                      FDR = mzero,
                      enr = m.b.enr,
                      log2enr = m.b.log2enr),
        colData = colData(x), rowData = rdat,
        metadata = mdat
    )
    return(se)
}
