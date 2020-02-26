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
                        strata = rep(1L, length(seqs))) {

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
        length(strata) == length(seqs) || (is.numeric(strata) && length(strata) == 1L)
    })

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
        log2pMM <- sapply(names(kmerFreq.stratum), function(current.kmer) {
            ii_long <- substr(rep(current.kmer, n),
                              start = 1:n, stop = 1:n + MMorder)
            ii_short <- substr(rep(current.kmer, n - 1L),
                               start = 2:n, stop = 1:(n - 1L) + MMorder)
            sum(lp_long[ii_long]) - sum(lp_short[ii_short])
        })
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
    padj <- p.adjust(p, method = "fdr")

    ## return results
    list(freq.obs = kmerFreq, freq.exp = kmerFreqMM,
         log2enr = lenr, sqrtDelta=sDelta, z = z, p = p, FDR = padj,
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
#'   use ordered by FDR. If \code{NULL} (the default), it will use \code{nKmers
#'   <- max(10, sum(x$FDR < 0.05))}. \code{nKmers} is ignored if \code{x} is a
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
            "FDR" %in% names(x)
            !is.null(names(x$FDR))
        })
        if (is.null(nKmers)) {
            nKmers <- max(10, sum(x$FDR < 0.05))
        }
        x <- names(x$FDR)[order(x$FDR)[seq.int(nKmers)]]
    }
    stopifnot(exprs = {
        is(x, "character")
        all(nchar(x) == nchar(x[1]))
        all(grepl("^[ACGT]+$", x))
        is(allowReverseComplement, "logical")
        length(allowReverseComplement) == 1L
    })
    kmerLen <- nchar(x[1])
    if (method == "cooccurrence") {
        stopifnot(exprs = {
            is(seqs, "DNAStringSet")
            is(zoops, "logical")
            length(zoops) == 1L
            is(n, "numeric")
            length(n) == 1L
        })
        if (n >= kmerLen)
            warning("Using n > k (n = ", n, ", k = ", kmerLen,
                    ") will count co-occurrences of non-overlapping k-mers.")
        message("clustering ", length(x), " ", kmerLen,
                "-mers (co-occurrence-based in ", length(seqs), " sequences)")
    } else if (method == "similarity") {
        if (is.null(maxShift)) {
            maxShift <- kmerLen - 2L
        }
        stopifnot(exprs = {
            is(maxShift, "numeric")
            length(maxShift) == 1L
            maxShift >= 0 && maxShift < kmerLen
        })
        if (is.null(minSim)) {
            minSim <- kmerLen - maxShift + 1
        }
        stopifnot(exprs = {
            is(minSim, "numeric")
            length(minSim) == 1L
            minSim >= 1 && minSim <= (kmerLen - 1L)
        })
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


#' @title Run a k-mer enrichment analysis.
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
#' @param pseudoCount A \code{numeric} scalar - will be added to the observed
#'   counts for each k-mer to avoid zero values.
#' @param Ncpu Number of parallel threads to use.
#' @param verbose A \code{logical} scalar. If \code{TRUE}, report on progress.
#'
#' @seealso \code{\link{getKmerFreq}} used to calculate k-mer enrichments;
#'   \code{\link[BSgenome]{getSeq,BSgenome-method}} which is used to extract
#'   sequences from \code{genomepkg} if \code{x} is a \code{GRanges} object;
#'   \code{\link[parallel]{mclapply}} that is used for parallelization;
#'   \code{\link{bin}} for binning of regions
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} \code{y}
#'   with \describe{ \item{assays(y)}{containing the four components \code{p},
#'   \code{FDR}, \code{enr} and \code{log2enr}), each a k-mer (rows) by bin
#'   (columns) matrix with raw -log10 P values, -log10 false discovery rates and
#'   k-mer enrichments as Pearson residuals (\code{enr}) and as log2 ratios
#'   (\code{log2enr}).} \item{rowData(x)}{containing information about the
#'   k-mers.} \item{colData(x)}{containing information about the bins.}
#'   \item{metaData(x)}{containing meta data on the object (e.g. parameter
#'   values).} }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame split
#' @importFrom parallel mclapply
#' @importFrom BSgenome getSeq
#' @importFrom TFBSTools PFMatrix PFMatrixList
#' @importFrom stats ppois p.adjust
#'
#' @export
kmerEnrichments <- function(x, b, genomepkg = NULL, kmerLen = 5,
                            background = c("other", "model"), MMorder = 1,
                            zoops = TRUE, pseudoCount = 1,
                            Ncpu = 2L, verbose = TRUE) {
    ## pre-flight checks
    background <- match.arg(background)
    stopifnot(exprs = {
        is.logical(verbose)
        length(verbose) == 1L
    })
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
        is.numeric(pseudoCount)
        length(pseudoCount) == 1L
        is.logical(zoops)
        length(zoops) == 1L
        is.numeric(Ncpu)
        length(Ncpu) == 1L
        Ncpu > 0
    })

    ## identify enriched k-mers in each bin
    if (verbose) {
        message("searching for enriched ", kmerLen, "-mers in ", nlevels(b),
                " bins using ", Ncpu, if (Ncpu > 1) " cores" else " core",
                " (background: ", c("other" = "sequences in other bins",
                                    "model" = paste0("Markov model of order ", MMorder))[background],
                ")...", appendLF = FALSE)
    }
    resL <- parallel::mclapply(split(x, b)[levels(b)], getKmerFreq, kmerLen = kmerLen,
                               MMorder = MMorder, pseudoCount = pseudoCount,
                               zoops = zoops, mc.cores = Ncpu)
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
        padj <- matrix(p.adjust(p, method="fdr"), nrow = nrow(p), dimnames = dimnames(p))
        assayL <- list(p = -log10(p), FDR = -log10(padj), enr = z, log2enr = lenr)

    } else if (background == "model") {
        p <- do.call(cbind, lapply(resL, "[[", "p"))
        padj <- matrix(p.adjust(p, method="fdr"), nrow = nrow(p), dimnames = dimnames(p))
        assayL <- list(p = -log10(p), FDR = -log10(padj),
                       enr = do.call(cbind, lapply(resL, "[[", "z")),
                       log2enr = do.call(cbind, lapply(resL, "[[", "log2enr")))
    }
    
    ## create SummarizedExperiment
    brks <- attr(b, "breaks")
    if (is.null(brks)) {
        brks <- rep(NA, nlevels(b) + 1L)
    }
    cdat <- S4Vectors::DataFrame(bin.names = levels(b),
                                 bin.lower = brks[-(nlevels(b)+1)],
                                 bin.upper = brks[-1],
                                 bin.nochange = seq.int(nlevels(b)) %in% attr(b, "bin0"))
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
        metadata = list(sequences = x,
                        bins = b,
                        bins.binmode = attr(b, "binmode"),
                        bins.breaks = as.vector(attr(b, "breaks")),
                        bins.bin0 = attr(b, "bin0"),
                        param.genomepkg = genomepkg,
                        param.kmerLen = kmerLen,
                        param.background = background,
                        param.MMorder = MMorder,
                        param.pseudoCount = pseudoCount,
                        param.zoops = zoops,
                        param.Ncpu = Ncpu,
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
#'   colData
#' @importFrom S4Vectors DataFrame metadata
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


#' @title Match a set of kmers to input sequences and determine frequencies of overlapping matches.
#'
#' @description Using a set of kmers, search input sequences for matches and retrieve
#'    the frequencies of sequences overlapping with 1 or more overlapping kmers.
#'
#'@param seqs Set of sequences, either a \code{character} vector or a
#'   \code{\link{DNAStringSet}}.
#'   
#' @param x A \code{character} vector of enriched kmers.
#'
#' @param Ncpu The number of CPU cores to used. \code{\link[parallel]{mclapply}} is used.
#'
#' @seealso \code{\link{kmerEnrichments}} for performing a k-mer enrichment
#'   analysis, \code{\link[parallel]{mclapply}} for how parallelization is done. 
#'
#' @export
extractOverlappingKmerFrequecies <- function(seqs, x, Ncpu = 1L) {
    ## pre-flight checks
    stopifnot(exprs = {
        class(seqs) == "DNAStringSet"
        class(x) == "character"
    })
    #also include reverse complements
    enriched.kmers <- unique(c(x, as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))))
    tmp <- unlist(parallel::mclapply(seq_along(seq), function(i){
        tmp.range <- IRanges::reduce(do.call(c, lapply(enriched.kmers, function(kmer){
            vmatchPattern(kmer, seq[i])[[1]]
        })))
        tmp.range
        res <- NULL
        if(length(tmp.range) > 0){
            res <- as.character(do.call(c, lapply(1:length(tmp.range), function(ii){
                subseq(seq[i], start=start(tmp.range)[ii], end=end(tmp.range)[ii])
            })))
        }
        res
    }, mc.cores=Ncpu))
    extended.seqs <- table(tmp)
    extended.seqs[order(extended.seqs, decreasing=TRUE)]
}





