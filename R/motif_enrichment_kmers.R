# dummy import to suppress R CMD check NOTE:
# "Namespace in Imports field not imported from: 'Rcpp'"
#' @importFrom Rcpp sourceCpp


#' @title Create matrix from consensus sequence
#'
#' @description Given a nucleotide sequence of A,C,G,T letter corresponding
#'   to a motif's consensus string, construct a positional frequency matrix.
#'   This matrix can for example be used as the \code{profileMatrix} argument
#'   in the constructor for a \code{TFBSTools::PFMatrix} object.
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
#' @param includeRevComp A \code{logcial} scalar. If \code{TRUE} (default),
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
        log2pMM <- sapply(names(kmerFreq.stratum), function(current.kmer) {
            ii_long <- substr(rep(current.kmer, n),
                              start = seq_len(n), stop = seq_len(n) + MMorder)
            ii_short <- substr(rep(current.kmer, n - 1L),
                               start = seq(2, n), stop = seq_len(n - 1L) + MMorder)
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
#' @param includeRevComp A \code{logical} scalar. If \code{TRUE}, also
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
#'   reverse complements, if \code{includeRevComp = TRUE}) will
#'   be extracted, and the resulting adjacency matrix will be converted to a
#'   graph, in which clusters will be identified as communities.
#'   For \code{method = "similarity"}, all pairwise k-mer distances for all
#'   possible shifts (defined by \code{maxShift}) are calculated, defined as
#'   the Hamming distance of the overlapping substring plus the number of shifts.
#'   For each k-mer pair, the minimal distance over shifts is retained. If
#'   \code{includeRevComp == TRUE}, this procedure is repeated to compare
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
clusterKmers <- function(x,
                         method = c("cooccurrence", "similarity"),
                         includeRevComp = FALSE,
                         nKmers = NULL,
                         maxShift = NULL,
                         minSim = NULL,
                         seqs = NULL,
                         zoops = TRUE,
                         n = 1L) {
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
    .assertScalar(x = includeRevComp, type = "logical")
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
        if (identical(includeRevComp, TRUE)) {
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
        if (identical(includeRevComp, TRUE)) {
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
#' @param zoops A \code{logical} scalar. If \code{TRUE} (the default), only one
#'   or zero occurrences of a k-mer are considered per sequence. This is helpful
#'   to reduce the impact of simple sequence repeats occurring in few sequences.
#' @param test type of motif enrichment test to perform.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @return a \code{data.frame} containing the motifs as rows and the columns:
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
                                zoops = TRUE,
                                test = c("fisher", "binomial"),
                                verbose = FALSE){
    
    # checks
    .assertScalar(x = k, type = "numeric", rngIncl = c(1, Inf))
    .checkDfValidity(df)
    .assertScalar(x = zoops, type = "logical")
    test <- match.arg(test)
    .assertScalar(x = verbose, type = "logical")
    
    totalWgtForeground <- sum(df$seqWgt[df$isForeground])
    totalWgtBackground <- sum(df$seqWgt[!df$isForeground])
    
    kmerHitMatrix <- oligonucleotideFrequency(x = df$seqs,
                                               width = k,
                                               as.prob = FALSE,
                                               with.labels = TRUE)
    
    if (zoops) {
        kmerHitMatrix[kmerHitMatrix > 0] <- 1
    }

    kmerHitMatrixWeighted <- kmerHitMatrix * df$seqWgt
    KmatchedSeqCountForeground <- colSums(kmerHitMatrixWeighted[df$isForeground, ])
    KmatchedSeqCountBackground <- colSums(kmerHitMatrixWeighted[!df$isForeground, ])
    
    # calculate k-mer enrichment
    if (identical(test, "binomial")) {
        
        if (verbose) {
            message("using binomial test to calculate ",
                    "log(p-values) for k-mer enrichments")
        }
        
        prob <- KmatchedSeqCountBackground / totalWgtBackground
        minProb <- 1 / totalWgtBackground
        maxProb <- (totalWgtBackground - 1) / totalWgtBackground
        if (any(i <- (prob < minProb))) {
            # warning("some background k-mer match probabilities are below ",
            #         "minProb (for example when there were zero hits) ",
            #         "and will be given a value of minProb=1/totalWgtBackground")
            prob[i] <- minProb
        }
        if (any(i <- (prob > maxProb))) {
            # warning("some k-mer match probabilities a above",
            #         "maxProb (for example when all sequences had hits) ",
            #         "and will be givena value of ",
            #         "maxProb=(totalWgtBackground-1)/totalWgtBackground")
            prob[i] <- maxProb
        }
        
        logP <- pbinom(q = KmatchedSeqCountForeground - 1,
                       size = totalWgtForeground,
                       prob = prob, lower.tail = FALSE, log.p = TRUE)
        
    } else if (identical(test, "fisher")) {
        
        if (verbose) {
            message("using fisher's exact test (one-sided) to calculate ",
                    "log(p-values) for motif enrichments")
        }
        
        # contingency table per motif for fisher's exact test (rounded to integer):
        #              withHit  noHit
        #   foreground    x       y
        #   background    z       w
        #
        logP <- log(vapply(structure(seq_along(KmatchedSeqCountForeground),
                                     names = names(KmatchedSeqCountForeground)),
                           function(i) {
                               ctab <- rbind(c(KmatchedSeqCountForeground[i],
                                               totalWgtForeground - KmatchedSeqCountForeground[i]),
                                             c(KmatchedSeqCountBackground[i],
                                               totalWgtBackground - KmatchedSeqCountBackground[i]))
                               ctab <- round(ctab)
                               fisher.test(x = ctab, alternative = "greater")$p.value
                           }, FUN.VALUE = numeric(1)))
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
#' @param pseudocount.pearsonResid A numerical scalar with the pseudocount to add
#'   to foreground and background frequencies when calculating expected counts
#'   and Pearson residuals.
#' @param zoops A \code{logical} scalar. If \code{TRUE} (the default), only one
#'   or zero occurrences of a k-mer are considered per sequence. This is helpful
#'   to reduce the impact of simple sequence repeats occurring in few sequences.
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
#' @param genome.seed An \code{integer} scalar used to seed the random number
#'   generator before sampling regions.
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
#'       generator is seeded using \code{genome.seed}.}
#'     \item{model}{: a Markov model of the order \code{MMorder} is estimated
#'       from the foreground sequences and used to estimate expected k-mer
#'       frequencies. K-mer enrichments are then calculated comparing observed
#'       to these expected frequencies.}
#'   }
#'
#'   For each k-mer, the weights of sequences is multiplied with the number
#'   of k-mer occurrences in each sequence and summed, separately for foreground
#'   (\code{sumForegroundWgtWithHits}) and background
#'   (\code{sumBackgroundWgtWithHits}) sequences. For \code{zoops = TRUE}
#'   (Zero-Or-One-Per-Sequence mode, default), at most one occurrence per
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
#'   with motifs in rows and bins in columns, containing six assays: \itemize{
#'   \item{negLog10P}{: -log10 P values}
#'   \item{negLog10Padj}{: -log10 adjusted P values}
#'   \item{pearsonResid}{: k-mer enrichments as Pearson residuals}
#'   \item{log2enr}{: k-mer enrichments as log2 ratios}
#'   \item{sumForegroundWgtWithHits}{: Sum of foreground sequence weights
#'     in a bin that have k-mer occurrences}
#'   \item{sumBackgroundWgtWithHits}{: Sum of background sequence weights
#'     in a bin that have k-mer occurrences}
#' }
#'
#' @seealso \code{\link{getKmerFreq}} used to calculate k-mer enrichments;
#'   \code{\link[BiocParallel]{bplapply}} that is used for parallelization;
#'   \code{\link{bin}} for binning of regions
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
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
                              maxFracN = 0.7,
                              maxKmerSize = 3L,
                              GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                           0.45, 0.5, 0.6, 0.7, 0.8),
                              pseudocount.kmers = 1,
                              pseudocount.log2enr = 8,
                              pseudocount.pearsonResid = 0.001,
                              zoops = TRUE,
                              p.adjust.method = "BH",
                              genome = NULL,
                              genome.regions = NULL,
                              genome.oversample = 2,
                              genome.seed = 42L,
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
    .assertScalar(x = pseudocount.kmers, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pseudocount.log2enr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pseudocount.pearsonResid, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = zoops, type = "logical")
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)
    if (identical(background, "zeroBin") &&
        (!"bin0" %in% names(attributes(bins)) || is.na(attr(bins, "bin0")))) {
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
        .assertScalar(x = genome.seed, type = "integer")
    }
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")
    if (is.null(names(seqs))) {
        names(seqs) <- paste0("s", seq_along(seqs))
    }
    
    
    ## filter sequences
    if (verbose) {
        message("Filtering sequences ...")
    }
    keep <- .filterSeqs(seqs, maxFracN = maxFracN, verbose = verbose)
    battr <- attributes(bins) # rescue attributes dropped by subsetting
    bins <- bins[keep]
    attr(bins, "binmode") <- battr$binmode
    attr(bins, "breaks") <- battr$breaks
    attr(bins, "bin0") <- battr$bin0
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
                                zoops = zoops)

            if (identical(test, "binomial")) {
                prob <- res1$freq.exp / Nfg
                minProb <- 1 / Nfg
                maxProb <- (Nfg - 1) / Nfg
                if (any(i <- (prob < minProb))) {
                    prob[i] <- minProb
                }
                if (any(i <- (prob > maxProb))) {
                    prob[i] <- maxProb
                }
                
                logP <- pbinom(q = res1$freq.obs - 1,
                               size = Nfg,
                               prob = prob, lower.tail = FALSE, log.p = TRUE)
                
            } else if (identical(test, "fisher")) {
                # contingency table per motif for fisher's exact test (rounded to integer):
                #              withHit  noHit
                #   foreground    x       y
                #   background    z       w
                #
                logP <- log(vapply(structure(seq_along(res1$freq.obs),
                                             names = names(res1$freq.obs)),
                                   function(i) {
                                       ctab <- rbind(c(res1$freq.obs[i],
                                                       Nfg - res1$freq.obs[i]),
                                                     c(res1$freq.exp[i],
                                                       Nfg - res1$freq.exp[i]))
                                       ctab <- round(ctab)
                                       fisher.test(x = ctab, alternative = "greater")$p.value
                                   }, FUN.VALUE = numeric(1)))
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
                                    gnm.seed = genome.seed + i,
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
                stop(paste0("No sequences remained after the GC weight calculation step in bin ", levels(bins)[i], 
                            " due to no GC bin containing both fore- and background sequences. ", 
                            "Cannot proceed with the ernichment analysis ..."))
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
                                           zoops = zoops,
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
        fracForeground <- enrich1[, "sumForegroundWgtWithHits"] / enrich1[, "totalWgtForeground"]
        fracBackground <- enrich1[, "sumBackgroundWgtWithHits"] / enrich1[, "totalWgtBackground"]
        obsTF <- enrich1[, "sumForegroundWgtWithHits"]
        expTF <- obsTF / (fracForeground + pseudocount.pearsonResid) * (fracBackground + pseudocount.pearsonResid)
        enr <- (obsTF - expTF) / sqrt(expTF)
        enr[ is.na(enr) ] <- 0
        names(enr) <- enrich1[, "motifName"]
        enr
    }))
    
    # log2 enrichments
    log2enr <- do.call(cbind, lapply(enrichL, function(enrich1) {
        D <- enrich1[, c("sumForegroundWgtWithHits", "sumBackgroundWgtWithHits")]
        nTot <- unlist(enrich1[1, c("totalWgtForeground", "totalWgtBackground")])
        D.norm <- t(t(D) / nTot * min(nTot))
        DL <- log2(D.norm + pseudocount.log2enr)
        log2enr <- DL[, 1] - DL[, 2]
        log2enr
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
                      bin.nochange = seq.int(nlevels(bins)) %in% attr(bins, "bin0"),
                      totalWgtForeground = do.call(c, lapply(enrichL, function(x){x$totalWgtForeground[1]})), 
                      totalWgtBackground = do.call(c, lapply(enrichL, function(x){x$totalWgtBackground[1]})))
    mdat <- list(bins = bins,
                 bins.binmode = attr(bins, "binmode"),
                 bins.breaks = as.vector(attr(bins, "breaks")),
                 bins.bin0 = attr(bins, "bin0"),
                 param = list(kmerLen = kmerLen,
                              background = background,
                              MMorder = MMorder,
                              test = test,
                              maxFracN = maxFracN,
                              maxKmerSize = maxKmerSize,
                              pseudocount.kmers = pseudocount.kmers,
                              pseudocount.log2enr = pseudocount.log2enr,
                              pseudocount.pearsonResid = pseudocount.pearsonResid,
                              zoops = zoops,
                              p.adj.method = p.adjust.method,
                              genome.class = class(genome),
                              genome.regions = genome.regions,
                              genome.oversample = genome.oversample,
                              genome.seed = genome.seed,
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
    
    se <- SummarizedExperiment(assays = list(negLog10P = P, 
                                             negLog10Padj = padj, 
                                             pearsonResid = enrTF,
                                             log2enr = log2enr, 
                                             sumForegroundWgtWithHits = assaySumForegroundWgtWithHits, 
                                             sumBackgroundWgtWithHits = assaySumBackgroundWgtWithHits),
                               rowData = rdat, colData = cdat, metadata = mdat)
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


#' @title Match a set of k-mers to input sequences and determine frequencies of
#'     overlapping matches.
#'
#' @description Using a set of k-mers, search input sequences for matches and
#'     retrieve the frequencies of sequences overlapping with any of the k-mers.
#'     Overlapping k-mer matches will result in extracted sequences that are
#'     longer than the input k-mers.
#'
#' @param seqs \code{\link{DNAStringSet}} with sequences to search.
#' @param x A \code{character} vector of k-mers to search in \code{seqs}.
#' @param includeRevComp A \code{logcial} scalar. If \code{TRUE} (default),
#'     scan both \code{seqs} and their reverse-complement for matches to \code{x}.
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
#' @importFrom Biostrings DNAStringSet reverseComplement vmatchPattern
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BSgenome getSeq
#' 
#' @export
extractOverlappingKmerFrequencies <- function(seqs,
                                              x,
                                              includeRevComp = TRUE,
                                              BPPARAM = SerialParam()) {
    ## pre-flight checks
    x <- toupper(x)
    stopifnot(exprs = {
        is(seqs, "DNAStringSet")
        is(x, "character")
        all(grepl("^[ACGT]+$", x))
        is(BPPARAM, "BiocParallelParam")
    })
    if (is.null(names(seqs)))
        names(seqs) <- paste0("s", seq_along(seqs))
    if (identical(includeRevComp, TRUE)) {
        seqsrc <- Biostrings::reverseComplement(seqs)
        names(seqsrc) <- paste0(names(seqs), "_rc")
        seqs <- c(seqs, seqsrc)
    }
    gr <- GenomicRanges::reduce(do.call(c, bplapply(seq_along(x), function(i) {
        hits <- Biostrings::vmatchPattern(pattern = x[i], subject = seqs)
        GenomicRanges::GRanges(seqnames = rep(names(seqs), lengths(hits)),
                               ranges = unlist(hits),
                               seqlengths = GenomeInfoDb::seqlengths(seqs))
    }, BPPARAM = BPPARAM)))
    tmp <- BSgenome::getSeq(seqs, gr)
    ttmp <- table(as.character(tmp))
    extended.seqs <- structure(as.vector(ttmp), names = names(ttmp))
    extended.seqs[order(extended.seqs, decreasing = TRUE)]
}


#' Build a directed graph from enriched k-mers
#' 
#' Build a directed graph from k-mers enriched in a set of sequences. The 
#' graph is based on a de Bruijn graph, in which edges are weighted based on 
#' the co-occurrence of neighboring k-mers in the provided input sequences. 
#' Putative motifs can be found by only retaining edges with weight above a 
#' user-defined threshold (see \code{getMotifsFromDirGraph}). 
#' 
#' @param seqs A \code{\link{DNAStringSet}} with sequences within which the 
#'   kmers in \code{x} are enriched. 
#' @param x A \code{character} vector of k-mers enriched in the sequences 
#'   in \code{seqs}.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'   instance determining the parallel back-end to be used during evaluation.
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' g <- buildDirGraphFromKmers(seqs, x)
#' hist(E(g)$weights)
#' getMotifsFromDirGraph(g, 10)
#' }
#' 
#' @importFrom igraph make_de_bruijn_graph vertex_attr induced_subgraph as_ids
#'   E
#' @importFrom Biostrings DNAStringSet mkAllStrings
#' @importFrom rlang .data
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter %>%
#' @importFrom tidyr gather unite
#' 
buildDirGraphFromKmers <- function(seqs, x, BPPARAM = SerialParam()) {
    x <- toupper(x)
    stopifnot(exprs = {
        is(seqs, "DNAStringSet")
        is(x, "character")
        all(grepl("^[ACGT]+$", x))
        length(unique(nchar(x))) == 1
    })
    
    # Build de Bruijn graph
    bases <- c("A", "C", "G", "T")
    ## m = number of letters to choose from
    ## n = k-mer length (the length of the sequences in x)
    k <- nchar(x[1])
    g <- igraph::make_de_bruijn_graph(m = length(bases), n = k)
    ## Set vertex names
    kmn <- Biostrings::mkAllStrings(alphabet = bases, width = k)
    igraph::vertex_attr(g, "name") <- kmn
    
    # Map k-mers back to the sequences
    # This is useful in case a central k-mer does not pass the 
    # enrichment threshold
    overlapkmers <- extractOverlappingKmerFrequencies(seqs, x,
                                                      BPPARAM = BPPARAM)
    
    # Get subgraph induced by kmers in the 'overlapping k-mers' only
    cooccs <- Biostrings::oligonucleotideFrequency(
        Biostrings::DNAStringSet(names(overlapkmers)), width = k
    )
    cooccs <- cooccs[, colSums(cooccs) > 0]
    gi <- igraph::induced_subgraph(
        g, vids = match(colnames(cooccs), igraph::vertex_attr(g, "name"))
    )
    
    # Assign a weight to each edge, representing the number of times the 
    # pair of k-mers appear just after each other in the input sequences
    kmpairs <- countKmerPairs(
        x = Biostrings::DNAStringSet(x = rep(names(overlapkmers),
                                             overlapkmers)), 
        k = k, n = 1, zoops = TRUE
    )
    
    ## Add weights to the edges based on the number of co-occurrences
    igraph::E(gi)$weight <- 0
    kmp <- as.data.frame(kmpairs) %>% 
        tibble::rownames_to_column("km1") %>% 
        tidyr::gather(key = "km2", value = "weight", -.data$km1) %>% 
        tidyr::unite(.data$km1, .data$km2, col = "node", sep = "|") %>%
        dplyr::filter(.data$node %in% igraph::as_ids(igraph::E(gi)))
    igraph::E(gi)$weight[match(kmp$node, igraph::as_ids(igraph::E(gi)))] <-
        kmp$weight

    gi
}


#' Filter directed graph
#' 
#' Filter directed graph, by only retaining edges with a 
#' weight exceeding a given threshold.
#' 
#' @param g A directed k-mer graph (e.g., from \code{buildDirGraphFromKmers})
#' @param edge_weight_thr An edge weight threshold. Edges with weight below 
#'   this threshold will be removed. 
#'   
#' @export
#' 
#' @importFrom igraph subgraph.edges V
#' 
filterDirGraph <- function(g, edge_weight_thr) {
    ## Plot graph, and subset to only the most abundant edges
    gisub <- igraph::subgraph.edges(
        g, eids = which(igraph::E(g)$weight > edge_weight_thr)
    )
    gisub
}

#' Infer putative motifs from directed graph
#' 
#' Infer putative motifs from a directed graph
#' 
#' @param seqs A \code{DNAStringSet} with sequences
#' @param g A directed k-mer graph, e.g. from \code{buildDirGraphFromKmers}
#'   or \code{filterDirGraph}
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'   instance determining the parallel back-end to be used during evaluation.
#' 
#' @export
#' 
#' @importFrom igraph vertex_attr as_ids
#' @importFrom BiocParallel SerialParam
#' @importFrom Biostrings DNAStringSet extractAt
#' @importFrom IRanges IRanges
#' 
getMotifsFromDirGraph <- function(seqs, g, BPPARAM = SerialParam()) {
    ## Get only the interesting parts of the input sequences
    overlapkmers <- extractOverlappingKmerFrequencies(
        seqs, igraph::vertex_attr(g, "name"),
        BPPARAM = BPPARAM
    )
    y <- Biostrings::DNAStringSet(x = names(overlapkmers))
    ## Extract all k-mers in order from each input sequence
    kmerlen <- nchar(igraph::vertex_attr(g, "name")[1])
    irl <- as(lapply(seq_along(y), function(i) {
        IRanges::IRanges(start = seq_len(width(y)[i] - (kmerlen - 1)),
                         width = kmerlen)
    }), "IRangesList")
    kmers <- Biostrings::extractAt(y, at = irl)
    ## Find matches in graph
    matches_to_graph <- lapply(kmers, function(kl) {
        igraph::vertex_attr(g, "name")[match(kl, igraph::vertex_attr(g, "name"))]
    })
    ## Insert NA between any k-mers that are not connected in the graph
    matches_to_graph <- lapply(matches_to_graph, function(x) {
        tmp <- sapply(seq_along(x), function(i) paste(x[i], x[i + 1], sep = "|"))
        tmp <- tmp %in% igraph::as_ids(igraph::E(g))
        w <- which(!tmp) + 0.5
        idx <- c(seq_along(x), w)
        vals <- c(x, rep(NA, length(w)))
        vals[order(idx)]
    })
    # matches_to_graph
    ## Create motifs (graph paths)
    helpfun <- function(v) {
        idx1 <- c(1, which(is.na(v)) + 1)
        idx2 <- setdiff(seq_along(v), idx1)
        v[idx2] <- substr(v[idx2], nchar(v[idx2]), nchar(v[idx2]))
        v[is.na(v)] <- "X"
        strsplit(paste(v, collapse = ""), "X")
    }
    motifs <- unlist(lapply(matches_to_graph, helpfun))
    # motifs
    ## Remove any motif that is a substring of another
    unique(motifs[!sapply(motifs, 
                          function(m) any(grepl(m, motifs[motifs != m], 
                                                fixed = TRUE)))])
}
