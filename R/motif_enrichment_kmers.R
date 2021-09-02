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
#' @param strata A \code{factor} or a \code{numeric} scalar defining the strata
#'   of sequences. A separate Markov model and expected k-mer frequencies are
#'   estimated for the set of sequences in each stratum (level in a \code{strata}
#'   factor). If \code{strata} is a scalar value, it will be interpreted as the
#'   number of strata to split the sequences into according to their CpG
#'   observed-over-expected counts using \code{kmeans(CpGoe, centers = strata)}.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#'
#' @return A \code{list} with observed and expected k-mer frequencies (\code{freq.obs}
#'   and \code{freq.exp}, respectively), and enrichment statistics for each k-mer.
#'
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @importFrom XVector subseq
#' @importFrom stats ppois kmeans
#'
#' @export
getKmerFreq <- function(seqs, kmerLen = 5, MMorder = 1, pseudoCount = 1, zoops = TRUE,
                        strata = rep(1L, length(seqs)), p.adjust.method = "BH") {

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
    stopifnot(exprs = {
        is(seqs, "DNAStringSet")
        round(kmerLen, 0L) == kmerLen
        round(MMorder, 0L) == MMorder
        length(strata) == length(seqs) || (is.numeric(strata) && length(strata) == 1L)
    })
    .assertScalar(x = kmerLen, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = MMorder, type = "numeric", rngExcl = c(0, kmerLen - 1L))
    .assertScalar(x = pseudoCount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = zoops, type = "logical")
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)
    
    ## split sequences into strata
    CpGoe <- NA
    if (length(strata) == 1L) {
        n1 <- oligonucleotideFrequency(seqs, width = 1L)
        n2 <- oligonucleotideFrequency(seqs, width = 2L)
        CpGoe <- n2[, "CG"] / (n1[,"C"] / rowSums(n1) * n1[,"G"])
        strata <- kmeans(x = CpGoe, centers = strata)$cluster
    }

    ## for each sequence stratum, calcluate...
    res.strata <- lapply(split(seqs, strata), function(seqs.stratum) {
        ## ... observed k-mer frequencies
        kmerFreqRaw.stratum <- oligonucleotideFrequency(seqs.stratum, width = kmerLen)
        if (zoops) {
            kmerFreq.stratum <- colSums(kmerFreqRaw.stratum > 0)
            #make new sequences that are simply the kmers detected
            seqs.stratum.zoops <- DNAStringSet(rep(names(kmerFreq.stratum), kmerFreq.stratum))
            lp_long  <- colSums(oligonucleotideFrequency(seqs.stratum.zoops, width = MMorder + 1L)) + pseudoCount
            lp_short <- colSums(oligonucleotideFrequency(seqs.stratum.zoops, width = MMorder)     ) + pseudoCount
        } else {
            kmerFreq.stratum <- colSums(kmerFreqRaw.stratum)
            lp_long  <- colSums(oligonucleotideFrequency(seqs.stratum, width = MMorder + 1L)) + pseudoCount
            lp_short <- colSums(oligonucleotideFrequency(seqs.stratum, width = MMorder))      + pseudoCount
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
    lenr <- log2((kmerFreq + pseudoCount) / (kmerFreqMM + pseudoCount))
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



#' @title Cluster k-mers
#'
#' @description Given enriched k-mers (typically the result of
#'   \code{\link{getKmerFreq}}), cluster k-mers into groups of similar k-mers
#'   that are likely to originate from the same motif.
#'
#' @param x Enriched k-mers, either a \code{character} vector or the return
#'   value of \code{\link{getKmerFreq}}, from which enriched k-mers can be
#'   extracted.
#' @param method A \code{character} scalar specifying the clustering method to use.
#'   Currently, either \code{"cooccurrence"} (co-occurrence of k-mers in a
#'   set of sequences, the default) or \code{"similarity"} (using k-mer
#'   similarities, see Details below).
#' @param allowReverseComplement A \code{logical} scalar. If \code{TRUE}, also
#'   include the reverse complement of enriched k-mers in the analysis. For
#'   \code{method = "cooccurrence"}, this will extract also reverse complement
#'   k-mers from the co-occurrence count matrix for graph estimation. For
#'   \code{method = "similarity"}, this will use the reverse complement of a
#'   k-mer if it has a higher toal similarity to all other k-mers than the
#'   original k-mer. In that case, the resulting k-mer is reported with the
#'   suffix \code{_rc} in the results, and the reverse-complement k-mer sequence
#'   is used to calculate pairwise similarities.
#' @param nKmers A \code{numeric} scalar. If \code{x} is a \code{list} as
#'   returned by \code{\link{getKmerFreq}}, the number of top enriched k-mers to
#'   use ordered by padj. If \code{NULL} (the default), it will use \code{nKmers
#'   <- max(10, sum(x$padj < 0.05))}. \code{nKmers} is ignored if \code{x} is a
#'   character vector.
#' @param maxShift (Only used for \code{method = "similarity"}.) A \code{numeric}
#'   scalar with the maximal number of shifts to perform when calculating k-mer
#'   similarities. If \code{NULL} (the default), it will use \code{maxShift <- kmerLen - 2}.
#' @param minSim (Only used for \code{method = "similarity"}.) A \code{numeric}
#'   scalar with the minimal k-mer similarity to link two k-mer nodes in the
#'   graph by an edge. If \code{NULL}, it will use \code{minSim <- kmerLen - maxShift + 1}.
#' @param seqs (Only used for \code{method = "cooccurrence"}.) A \code{DNAStringSet}
#'   with the sequences, in which k-mer co-occurrences will be counted using
#'   \code{\link{countKmerPairs}}.
#' @param zoops (Only used for \code{method = "cooccurrence"}.) A \code{logical}
#'   scalar passed to \code{\link{countKmerPairs}}.
#' @param n (Only used for \code{method = "cooccurrence"}.) An integer scalar
#'   defining the maximum downstream distance of second k-mers, relative to the
#'   start position of the first k-mer (default: \code{n = NULL}).
#'
#' @details The clustering is performed depending on \code{method}.
#'   For \code{method = "cooccurrence"},
#'   \code{countKmerPairs(x = seqs, k = k, n = 1, zoops = zoops)}
#'   will be used to first get a pairwise co-occurrence count matrix.
#'   Rows and columns corresponding to enriched k-mers from \code{x} (and their
#'   reverse complements, if \code{allowReverseComplement = TRUE}) will
#'   be extracted, and the resulting adjacency matrix will be converted to a
#'   graph, in which clusters will be identified as communities.
#'   For \code{method = "similarity"}, all pairwise k-mer distances for all
#'   possible shifts (defined by \code{maxShift}) are calculated, defined as
#'   the Hamming distance of the overlapping substring plus the number of shifts.
#'   For each k-mer pair, the minimal distance over shifts is retained. If
#'   \code{allowReverseComplement == TRUE}, this procedure is repeated to compare
#'   each k-mer to all reverse-complemented k-mers, and replace it with the
#'   reverse-complemented version if this yields a lower sum of pairwise distances.
#'   The resulting distance matrix is then converted into a similarity matrix by
#'   subtracting it from the length of the k-mers. Elements less than \code{minSim}
#'   are set to zero, and the matrix is used as an adjacency matrix to construct
#'   a k-mer graph. The clusters are identified as the connected components in
#'   this graph.
#'
#' @return A named \code{numeric} vector, with k-mer as names and values
#'   indicating the k-mer cluster memberships. In \code{attr(y, "graph")} (where
#'   \code{y} is the return value), you can get the k-mer graph on which the
#'   clustering is based.
#'
#' @seealso \code{\link{getKmerFreq}} for finding enriched k-mers in a set of
#'   sequences; \code{\link{countKmerPairs}} for counting k-mer co-occurrences
#'   in a set of sequences.
#' 
#' @importFrom stringdist stringdistmatrix
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain clusters communities
#' @importFrom Biostrings reverseComplement DNAStringSet
#'
#' @export
clusterKmers <- function(x, method = c("cooccurrence", "similarity"),
                         allowReverseComplement = FALSE,
                         nKmers = NULL, maxShift = NULL, minSim = NULL,
                         seqs = NULL, zoops = TRUE, n = 1L) {
    ## pre-flight checks
    method <- match.arg(method)
    if (is(x, "list")) {
        stopifnot(exprs = {
            "padj" %in% names(x)
            !is.null(names(x$padj))
        })
        if (is.null(nKmers)) {
            nKmers <- max(10, sum(x$padj < 0.05))
        }
        x <- names(x$padj)[order(x$padj)[seq.int(nKmers)]]
    }
    stopifnot(exprs = {
        is(x, "character")
        all(nchar(x) == nchar(x[1]))
        all(grepl("^[ACGT]+$", x))
    })
    .assertScalar(x = allowReverseComplement, type = "logical")
    kmerLen <- nchar(x[1])
    if (method == "cooccurrence") {
        stopifnot(is(seqs, "DNAStringSet"))
        .assertScalar(x = zoops, type = "logical")
        .assertScalar(x = n, type = "numeric", rngIncl = c(1, Inf))
        if (n >= kmerLen)
            warning("Using n > k (n = ", n, ", k = ", kmerLen,
                    ") will count co-occurrences of non-overlapping k-mers.")
        message("clustering ", length(x), " ", kmerLen,
                "-mers (co-occurrence-based in ", length(seqs), " sequences)")
    } else if (method == "similarity") {
        if (is.null(maxShift)) {
            maxShift <- kmerLen - 2L
        }
        .assertScalar(x = maxShift, type = "numeric", rngIncl = c(0, kmerLen - 1))
        if (is.null(minSim)) {
            minSim <- kmerLen - maxShift + 1
        }
        .assertScalar(x = minSim, type = "numeric", rngIncl = c(1, kmerLen - 1))
        message("clustering ", length(x), " ", kmerLen,
                "-mers (similarity-based, maxShift: ", maxShift, ", minSim: ", minSim, ")")
    }
    
    if (method == "cooccurrence") {
        ## calculate pairwise co-occurrence matrix
        co <- countKmerPairs(x = seqs, k = kmerLen, n = n, zoops = zoops)
        ## select enriched k-mers
        if (allowReverseComplement) {
            xrc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
            xsel <- unique(c(x, xrc))
        } else {
            xsel <- x
        }
        co <- co[xsel, xsel]
        ## construct graph
        G <- igraph::graph_from_adjacency_matrix(co, mode = "undirected", weighted = TRUE)
        ## find communities
        comm <- igraph::cluster_louvain(G)
        res <- igraph::membership(comm)
        attr(res, "graph") <- G
        
    } else if (method == "similarity") {
        ## calculate pairwise k-mer distances
        ## distance metric: "hamming + shifts"
        ##   k-mers are allowed to be shifted relative to one another
        ##   a shift of one position is identical to one mismatch position
        ## ... first take the minimum hamming distance over all possible shifts
        d <- Reduce(pmin, lapply(seq(from = 0, to = min(kmerLen - 1, maxShift)),
                                 function(s) {
                                     sufx <- substr(x, start = s + 1, stop = kmerLen)
                                     prex <- substr(x, start = 1, stop = kmerLen - s)
                                     stringdist::stringdistmatrix(sufx, prex, method = "hamming") + s
                                 }))
        ## ... then take the minimal distance between left and right shifts
        d[upper.tri(d)] <- pmin(d[upper.tri(d)], t(d)[upper.tri(d)])
        d[lower.tri(d)] <- t(d)[lower.tri(d)]
        colnames(d) <- x
        ## ... repeat by comparing x to reverseComplement(x) and keep the minium distance
        if (allowReverseComplement) {
            message("also considering reverse complement k-mers")
            xrc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
            drc <- Reduce(pmin, lapply(seq(from = 0, to = min(kmerLen - 1, maxShift)),
                                       function(s) {
                                           sufx <- substr(x, start = s + 1, stop = kmerLen)
                                           prex <- substr(xrc, start = 1, stop = kmerLen - s)
                                           stringdist::stringdistmatrix(sufx, prex, method = "hamming") + s
                                       }))
            drc[upper.tri(drc)] <- pmin(drc[upper.tri(drc)], t(drc)[upper.tri(drc)])
            drc[lower.tri(drc)] <- t(drc)[lower.tri(drc)]
            ## ... for each k-mer, decide wether to use the forward or reverse-complement
            useRC <- which(colSums(drc) < colSums(unname(d)))
            if (length(useRC) > 0) {
                colnames(d)[useRC] <- paste0(x[useRC], "_rc")
                d[, useRC] <- drc[, useRC]
                d[useRC, ] <- drc[useRC, ]
            }
        }
        
        ## construct and prune k-mer graph (adjacency matrix)
        ## ... convert the distance into a similarity (adjacency) matrix
        a <- kmerLen - d
        ## ... prune and construct the graph
        a[a < minSim] <- 0
        diag(a) <- 0
        G <- igraph::graph_from_adjacency_matrix(a, mode = "undirected", weighted = TRUE)
        
        ## find communities
        #comm <- igraph::cluster_louvain(G)
        comm <- igraph::clusters(G)
        res <- igraph::membership(comm)
        attr(res, "graph") <- G
    }
    
    ## return results
    return(res)
}


#' @title Calculate k-mer enrichment in bins of sequences.
#'
#' @description Given a set of sequences and corresponding bins, identify
#'   enriched k-mers (n-grams) in each bin. The sequences can be given either
#'   directly or as genomic coordinates.
#'
#' @param x A \code{character} vector, \code{\link[Biostrings]{DNAStringSet}} or
#'   a \code{\link[GenomicRanges]{GRanges}} object with the sequences to
#'   analyze.
#' @param b A vector of the same length as \code{x} that groups its elements
#'   into bins (typically a factor, such as the one returned by
#'   \code{\link{bin}}).
#' @param genomepkg Only used if \code{x} is a \code{GRanges} object: A
#'   \code{character} scalar with the name of a \code{BSgenome} package from
#'   which to extract the sequences.
#' @param kmerLen A \code{numeric} scalar giving the k-mer length.
#' @param background A \code{character} scalar. If \code{"other"} (the default),
#'   the enrichments in a bin are calculated compared to the sequences in all
#'   other bins. If \code{"model"}, the enrichments are relative to the expected
#'   frequencies obtained from a Markov model of order \code{MMorder}.
#' @param MMorder A \code{numeric} scalar giving the order of the Markov model
#'   used to calculate the expected frequencies for \code{background = "model"}.
#' @param zoops A \code{logical} scalar. If \code{TRUE} (the default), only one
#'   or zero occurences of a k-mer are considered per sequence. This is helpful
#'   to reduce the impact of simple sequence repeats occurring in few sequences.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#' @param pseudoCount A \code{numeric} scalar - will be added to the observed
#'   counts for each k-mer to avoid zero values.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation.
#' @param verbose A \code{logical} scalar. If \code{TRUE}, report on progress.
#'
#' @example 
#' seqs <- DNAStringSet(c("GCATGCATGC", "CTAGCTAGCTG"))
#' bins <- factor(1:2)
#' calcBinnedKmerEnr(x = seqs, b = bins, kmerLen = 3)
#' 
#' @seealso \code{\link{getKmerFreq}} used to calculate k-mer enrichments;
#'   \code{\link[BSgenome]{getSeq,BSgenome-method}} which is used to extract
#'   sequences from \code{genomepkg} if \code{x} is a \code{GRanges} object;
#'   \code{\link[BiocParallel]{bplapply}} that is used for parallelization;
#'   \code{\link{bin}} for binning of regions
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with
#'   k-mers in rows and bins in columns, containing four assays: \itemize{
#'   \item{negLog10P}{: -log10 P values}
#'   \item{negLog10Padj}{: -log10 adjusted P values}
#'   \item{pearsonResid}{: motif enrichments as Pearson residuals}
#'   \item{log2enr}{: motif enrichments as log2 ratios}
#' }.
#' rowData(x) contains information about the k-mers, colData(x) contains
#' information about the bins, and metaData(x) contains meta data on the object
#' and how it was generated (e.g. parameter values).
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame split
#' @importFrom BSgenome getSeq
#' @importFrom TFBSTools PFMatrix PFMatrixList
#' @importFrom stats ppois p.adjust
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#'
#' @export
calcBinnedKmerEnr <- function(x,
                              b,
                              genomepkg = NULL,
                              kmerLen = 5,
                              background = c("other", "model"),
                              MMorder = 1,
                              zoops = TRUE,
                              p.adjust.method = "BH",
                              pseudoCount = 1,
                              BPPARAM = SerialParam(),
                              verbose = FALSE) {
    ## pre-flight checks
    background <- match.arg(background)
    .assertScalar(x = verbose, type = "logical")
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
        b <- factor(b, levels = unique(b))
    }
    .assertVector(x = b, len = length(x))
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)
    .assertScalar(x = pseudoCount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = zoops, type = "logical")
    stopifnot(is(BPPARAM, "BiocParallelParam"))
    
    ## identify enriched k-mers in each bin
    if (verbose) {
        tmpmsg <- paste0("searching for enriched ", kmerLen, "-mers in ", nlevels(b),
                         " bins using ", bpnworkers(BPPARAM),
                         if (bpnworkers(BPPARAM) > 1) " cores" else " core",
                         " (background: ", c("other" = "sequences in other bins",
                                             "model" = paste0("Markov model of order ", MMorder))[background],
                         ")...")
        message(tmpmsg, appendLF = FALSE)
    }
    resL <- bplapply(split(x, b)[levels(b)], getKmerFreq, kmerLen = kmerLen,
                     MMorder = MMorder, pseudoCount = pseudoCount,
                     zoops = zoops, BPPARAM = BPPARAM)
    if (verbose) {
        message("done")
    }
    
    ## calculate motif enrichments
    if (background == "other") {
        f.obs <- do.call(cbind, lapply(resL, "[[", "freq.obs"))
        f.exp <- do.call(cbind, lapply(seq_along(resL), function(i) {
            rowSums(do.call(cbind, lapply(resL[-i], "[[", "freq.obs")))
        }))
        colnames(f.exp) <- colnames(f.obs)
        n.obs <- colSums(f.obs)
        n.exp <- colSums(f.exp)
        f.obs <- t(t(f.obs) / n.obs * pmin(n.obs, n.exp))
        f.exp <- t(t(f.exp) / n.exp * pmin(n.obs, n.exp))
        lenr <- log2((f.obs + pseudoCount) / (f.exp + pseudoCount))
        z <- (f.obs - f.exp) / sqrt(f.exp)
        p <- ppois(q = f.obs, lambda = f.exp, lower.tail = FALSE)
        padj <- matrix(p.adjust(p, method = p.adjust.method), nrow = nrow(p), dimnames = dimnames(p))
        assayL <- list(negLog10P = -log10(p), negLog10Padj = -log10(padj),
                       pearsonResid = z, log2enr = lenr)
        
    } else if (background == "model") {
        p <- do.call(cbind, lapply(resL, "[[", "p"))
        padj <- matrix(p.adjust(p, method = p.adjust.method), nrow = nrow(p), dimnames = dimnames(p))
        assayL <- list(negLog10P = -log10(p), negLog10Padj = -log10(padj),
                       pearsonResid = do.call(cbind, lapply(resL, "[[", "z")),
                       log2enr = do.call(cbind, lapply(resL, "[[", "log2enr")))
    }
    
    ## create SummarizedExperiment
    brks <- attr(b, "breaks")
    if (is.null(brks)) {
        brks <- rep(NA, nlevels(b) + 1L)
    }
    cdat <- S4Vectors::DataFrame(bin.names = levels(b),
                                 bin.lower = brks[-(nlevels(b) + 1)],
                                 bin.upper = brks[-1],
                                 bin.nochange = seq.int(nlevels(b)) %in% getZeroBin(b))
    kmers <- names(resL[[1]][[1]])
    pfms <- do.call(TFBSTools::PFMatrixList, lapply(kmers, function(kmer) {
        TFBSTools::PFMatrix(ID = kmer, name = kmer,
                            profileMatrix = .cons2matrix(kmer))
    }))
    percentGC <- unlist(lapply(pfms, function(x) {
        m <- TFBSTools::Matrix(x)
        100 * sum(m[c("C","G"), ]) / sum(m)
    }), use.names = FALSE)
    rdat <- S4Vectors::DataFrame(motif.name = kmers,
                                 motif.pfm = pfms,
                                 motif.percentGC = percentGC)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assayL,
        colData = cdat, rowData = rdat,
        metadata = list(bins = b,
                        bins.binmode = attr(b, "binmode"),
                        bins.breaks = as.vector(attr(b, "breaks")),
                        bins.bin0 = getZeroBin(b),
                        param = list(genomepkg = genomepkg,
                                     kmerLen = kmerLen,
                                     background = background,
                                     MMorder = MMorder,
                                     p.adjust.method = p.adjust.method,
                                     pseudoCount = pseudoCount,
                                     zoops = zoops,
                                     Ncpu = bpnworkers(BPPARAM)),
                        motif.distances = NULL)
    )
    rownames(se) <- kmers
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
#'   \code{\link{calcBinnedKmerEnr}}).
#' @param m Either a \code{\link[TFBSTools]{PFMatrixList}}, or a character
#'   scalar with a file containing motifs in HOMER format (loaded into a
#'   \code{\link[TFBSTools]{PFMatrixList}} by
#'   \code{\link{homerToPFMatrixList}}).
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'     with \code{length(m)} rows (motifs) and \code{ncol(x)} columns (bins).
#'
#' @seealso \code{\link{calcBinnedKmerEnr}} for performing a k-mer enrichment
#'   analysis, \code{\link[BiocParallel]{bplapply}} used for parallelization.
#'
#' @importFrom SummarizedExperiment assayNames assay SummarizedExperiment
#'   colData
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#'
#' @export
convertKmersToMotifs <- function(x, m, BPPARAM = SerialParam(), verbose = FALSE) {
    ## pre-flight checks
    .assertScalar(x = verbose, type = "logical")
    if (is.character(m) && length(m) == 1L && file.exists(m)) {
        if (verbose) {
            message("reading motifs from ", basename(m))
        }
        m <- homerToPFMatrixList(m)
    }
    stopifnot(exprs = {
        is(x, "SummarizedExperiment")
        "param" %in% names(metadata(x)) && "kmerLen" %in% names(metadata(x)$param)
        nrow(x) == 4^metadata(x)$param.kmerLen
        all(c("pearsonResid","log2enr") %in% assayNames(x))
        is(m, "PFMatrixList")
        is(BPPARAM, "BiocParallelParam")
    })
    
    ## calculate motif-by-kmer matrix (probabilities)
    m.k <- motifKmerSimilarity(m, kmerLen = metadata(x)$param$kmerLen,
                               BPPARAM = BPPARAM, verbose = verbose)
    
    ## calculate motif-by-bin = motif-by-kmer %*% kmer-by-bin
    m.b.enr <- m.k %*% assay(x, "pearsonResid")
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
    mdat[["param.Ncpu"]] <- bpnworkers(BPPARAM)
    mzero <- matrix(0, nrow = nrow(m.b.enr), ncol = ncol(m.b.enr),
                    dimnames = dimnames(m.b.enr))
    se <- SummarizedExperiment(
        assays = list(negLog10P = mzero,
                      negLog10Padj = mzero,
                      pearsonResid = m.b.enr,
                      log2enr = m.b.log2enr),
        colData = colData(x), rowData = rdat,
        metadata = mdat
    )
    return(se)
}


#' @title Match a set of k-mers to input sequences and determine frequencies of overlapping matches.
#'
#' @description Using a set of k-mers, search input sequences for matches and retrieve
#'      the frequencies of sequences overlapping with 1 or more overlapping k-mers.
#'
#' @param seqs Set of sequences, either a \code{character} vector or a
#'    \code{\link{DNAStringSet}}.
#'   
#' @param x A \code{character} vector of enriched k-mers.
#'
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation.
#'
#' @return A named integer vector with the number of observed occurrences of
#'   sequences overlapping enriched k-mers. The sequences are given as the
#'   names, and the elements are sorted by decreasing frequency.
#'
#' @seealso \code{\link{calcBinnedKmerEnr}} for performing a k-mer enrichment
#'     analysis, \code{\link[BiocParallel]{bplapply}} used for parallelization. 
#'
#' @importFrom Biostrings matchPDict reverseComplement DNAStringSet
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom IRanges reduce
#' 
#' @export
extractOverlappingKmerFrequecies <- function(seqs, x, BPPARAM = SerialParam()) {
    ## pre-flight checks
    x <- toupper(x)
    stopifnot(exprs = {
        is(seqs, "DNAStringSet")
        is(x, "character")
        all(grepl("^[ACGT]+$", x))
        is(BPPARAM, "BiocParallelParam")
    })
    #also include reverse complements
    # enriched.kmers <- unique(c(x, as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))))
    # tmp <- unlist(bplapply(seq_along(seqs), function(i) {
    #     tmp.range <- IRanges::reduce(do.call(c, lapply(enriched.kmers, function(kmer) {
    #         Biostrings::vmatchPattern(kmer, seqs[i])[[1]]
    #     })))
    #     tmp.range
    #     res <- NULL
    #     if (length(tmp.range) > 0) {
    #         res <- as.character(do.call(c, lapply(seq_along(tmp.range), function(ii) {
    #             subseq(seqs[i], start = start(tmp.range)[ii], end = end(tmp.range)[ii])
    #         })))
    #     }
    #     res
    # }, BPPARAM = BPPARAM))
    xx <- Biostrings::DNAStringSet(x)
    xxrc <- unique(c(xx, Biostrings::reverseComplement(xx)))
    xxdict <- Biostrings::PDict(xxrc)
    tmp <- unlist(bplapply(seq_along(seqs), function(i) {
        hitranges <- IRanges::reduce(unlist(matchPDict(pdict = xxdict, subject = seqs[[i]])))
        if (length(hitranges)) {
            substring(as.character(seqs[[i]]), first = start(hitranges), last = end(hitranges))
        } else {
            character(0)
        }
    }, BPPARAM = BPPARAM))
    ttmp <- table(tmp)
    extended.seqs <- structure(as.vector(ttmp), names = names(ttmp))
    extended.seqs[order(extended.seqs, decreasing = TRUE)]
}





