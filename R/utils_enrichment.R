.binomEnrichmentTest <- function(matchCountBg, totalWeightBg, matchCountFg,
                                 totalWeightFg, verbose) {
    .assertVector(matchCountBg, type = "numeric")
    .assertScalar(totalWeightBg, type = "numeric")
    .assertVector(matchCountFg, type = "numeric")
    .assertScalar(totalWeightFg, type = "numeric")
    .assertScalar(verbose, type = "logical")
    
    if (verbose) {
        message("using binomial test to calculate ",
                "log(p-values) for enrichments")
    }
    
    prob <- matchCountBg / totalWeightBg
    minProb <- 1 / totalWeightBg
    maxProb <- (totalWeightBg - 1) / totalWeightBg
    if (any(i <- (prob < minProb))) {
        # warning("some background match probabilities are below ",
        #         "minProb (for example when there were zero hits) ",
        #         "and will be given a value of minProb=1/totalWeightBg")
        prob[i] <- minProb
    }
    if (any(i <- (prob > maxProb))) {
        # warning("some match probabilities a above",
        #         "maxProb (for example when all sequences had hits) ",
        #         "and will be given a value of ",
        #         "maxProb=(totalWeightBg-1)/totalWeightBg")
        prob[i] <- maxProb
    }
    
    pbinom(q = matchCountFg - 1,
           size = totalWeightFg,
           prob = prob, lower.tail = FALSE, log.p = TRUE)
}

.fisherEnrichmentTest <- function(matchCountBg, totalWeightBg, matchCountFg,
                                  totalWeightFg, verbose) {
    .assertVector(matchCountBg, type = "numeric")
    .assertScalar(totalWeightBg, type = "numeric")
    .assertVector(matchCountFg, type = "numeric")
    .assertScalar(totalWeightFg, type = "numeric")
    .assertScalar(verbose, type = "logical")
    
    if (verbose) {
        message("using Fisher's exact test (one-sided) to calculate ",
                "log(p-values) for enrichments")
    }
    
    # contingency table per sequence for Fisher's exact test (rounded to integer):
    #              withHit  noHit
    #   foreground    x       y
    #   background    z       w
    #
    logP <- log(vapply(structure(seq_along(matchCountFg),
                                 names = names(matchCountFg)),
                       function(i) {
                           ctab <- rbind(c(matchCountFg[i],
                                           totalWeightFg - matchCountFg[i]),
                                         c(matchCountBg[i],
                                           totalWeightBg - matchCountBg[i]))
                           ctab <- round(ctab)
                           fisher.test(x = ctab, alternative = "greater")$p.value
                       }, FUN.VALUE = numeric(1)))
}

.calcPearsonResiduals <- function(matchCountBg, totalWeightBg, matchCountFg,
                                  totalWeightFg) {
    .assertVector(matchCountBg, type = "numeric")
    .assertScalar(totalWeightBg, type = "numeric")
    .assertVector(matchCountFg, type = "numeric")
    .assertScalar(totalWeightFg, type = "numeric")

    obsTF <- matchCountFg
    expTF <- totalWeightFg * (matchCountFg + matchCountBg) /
        (totalWeightFg + totalWeightBg)
    N <- totalWeightFg + totalWeightBg
    enr <- (obsTF - expTF) / sqrt(expTF * (1 - totalWeightFg / N) *
                                      (1 - (matchCountFg +
                                                matchCountBg) / N))
    enr[is.na(enr)] <- 0
    enr
}

.calcExpFg <- function(matchCountBg, totalWeightBg, matchCountFg, 
                       totalWeightFg) {
    .assertVector(matchCountBg, type = "numeric")
    .assertScalar(totalWeightBg, type = "numeric")
    .assertVector(matchCountFg, type = "numeric")
    .assertScalar(totalWeightFg, type = "numeric")

    totalWeightFg * (matchCountFg + matchCountBg) / 
        (totalWeightFg + totalWeightBg)
}

.calcLog2Enr <- function(matchCountBg, totalWeightBg, matchCountFg, 
                         totalWeightFg, pseudocount) {
    .assertVector(matchCountBg, type = "numeric")
    .assertScalar(totalWeightBg, type = "numeric")
    .assertVector(matchCountFg, type = "numeric")
    .assertScalar(totalWeightFg, type = "numeric")

    minTot <- min(totalWeightFg, totalWeightBg)
    normFg <- log2(matchCountFg/totalWeightFg * minTot + pseudocount) 
    normBg <- log2(matchCountBg/totalWeightBg * minTot + pseudocount)
    log2enr <- normFg - normBg
    log2enr
    
    # D <- enrich1[, c("sumForegroundWgtWithHits", "sumBackgroundWgtWithHits")]
    # nTot <- unlist(enrich1[1, c("totalWgtForeground", "totalWgtBackground")])
    # D.norm <- t(t(D) / nTot * min(nTot))
    # DL <- log2(D.norm + pseudocount.log2enr)
    # log2enr <- DL[, 1] - DL[, 2]
    # log2enr
}

#' @title Check if seqinfo DataFrame is valid
#'
#' @description Check if the DataFrame with sequence information is valid,
#'   i.e. is of the correct object type (DataFrame) and has all expected
#'   columns and attributes.
#'
#' @param df Input object to be checked. It should have an attribute \code{err}
#'   and columns:
#'   \itemize{
#'     \item{\code{seqs}}{: a \code{DNAStringSet} object.}
#'     \item{\code{isForeground}}{ that indicates if a sequence is in the
#'                                 foreground group.}
#'     \item{\code{GCfrac}}{: the fraction of G+C bases per sequence.}
#'     \item{\code{GCbin}}{: the GC bin for each sequence.}
#'     \item{\code{GCwgt}}{: the sequence weight to adjust for GC
#'       differences between foreground and background sequences.}
#'     \item{\code{seqWgt}}{: the sequence weight to adjust for k-mer
#'       differences between foreground and background sequences.}
#'   }
#'
#' @return \code{TRUE} (invisibly) if \code{df} is valid, otherwise it
#'   raises an exception using \code{stop()}
#' 
#' @keywords internal
.checkDfValidity <- function(df) {
    expected_cols <- c("seqs", "isForeground",
                       "GCfrac", "GCbin", "GCwgt", "seqWgt")
    expected_types <- c("DNAStringSet", "logical",
                        "numeric", "numeric", "numeric", "numeric")
    expected_attrs <- c("err")
    
    if (!is(df, "DataFrame")) {
        stop("'df' should be a DataFrame, but it is a ", class(df))
        
    } else if (!all(expected_cols %in% colnames(df))) {
        stop("'df' has to have columns: ",
             paste(expected_cols, collapse = ", "))
        
    } else if (!all(unlist(lapply(seq_along(expected_cols), function(i) {
        is(df[, expected_cols[i]], expected_types[i])
    })))) {
        stop("Not all columns in 'df' have the expected types:\n ",
             paste(paste0(expected_cols, ": '", expected_types, "'"),
                   collapse = "\n "))
        
    } else if (!all(expected_attrs %in% names(attributes(df)))) {
        stop("'df' has to have attributes: ",
             paste(expected_attrs, collapse = ", "))
    }
    
    return(invisible(TRUE))
}

#' @title Filter Sequences
#'
#' @description Filter sequences that are unlikely to be useful for motif
#'   enrichment analysis. The current defaults are based on HOMER (version 4.11).
#'
#' @param seqs a \code{DNAStringSet} object.
#' @param maxFracN A numeric scalar with the maximal fraction of N bases allowed
#'   in a sequence (defaults to 0.7).
#' @param minLength The minimum sequence length (default from Homer).
#'   Sequences shorter than this will be filtered out.
#' @param maxLength The maximum sequence length (default from Homer).
#'   Sequences bigger than this will be filtered out.
#' @param verbose A logical scalar. If \code{TRUE}, report on filtering.
#'
#' @details The filtering logic is based on \code{removePoorSeq.pl} from Homer.
#'
#' @return a logical vector of the same length as \code{seqs} with \code{TRUE}
#'   indicated to keep the sequence and \code{FALSE} to filter it out.
#'
#' @importFrom Biostrings alphabetFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
#' 
#' @keywords internal
.filterSeqs <- function(seqs, maxFracN = 0.7, minLength = 5L,
                        maxLength = 100000L, verbose = FALSE) {
    
    if (!is(seqs, "DNAStringSet")) {
        stop("'seqs' must be a DNAStringSet object.")
    }
    .assertScalar(x = maxFracN,  type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minLength, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = maxLength, type = "numeric", rngIncl = c(max(minLength, 0), Inf))
    .assertScalar(x = verbose,   type = "logical")
    
    # fraction of N bases per sequence
    observedFracN <- alphabetFrequency(seqs, as.prob = TRUE)[, "N"]
    
    f1 <- width(seqs) > 0 & observedFracN > maxFracN
    if (sum(f1) > 0 && verbose) {
        message("  ", sum(f1), " of ", length(seqs),
                " sequences (", round(100 * sum(f1) / length(seqs), 1), "%)",
                " have too many N bases")
    }
    
    f2 <- (width(seqs) < minLength) | (width(seqs) > maxLength)
    if (sum(f2) > 0 && verbose) {
        message("  ", sum(f2), " of ", length(seqs),
                " sequences (", round(100 * sum(f2) / length(seqs), 1), "%)",
                " are too short or too long")
    }
    
    res <- !(f1 | f2)
    if (verbose) {
        message("  in total filtering out ", sum(!res), " of ", length(seqs),
                " sequences (", round(100 * sum(!res) / length(seqs), 1), "%)")
    }
    
    return(res)
}

#' @title Define background sequence set for a single motif enrichment calculation
#' 
#' @description Define the background set for the motif enrichment calculation
#'   in a single bin, depending on the background mode and given foreground
#'   sequences.
#' 
#' @param sqs,bns,bg The \code{seqs}, \code{bins} and \code{background} arguments
#'   from \code{calcBinnedMotifEnrR}.
#' @param currbn An \code{integer} scalar with the current bin defining the
#'   foreground sequences.
#' @param gnm,gnm.regions,gnm.oversample The \code{genome}, \code{genome.regions}
#'   and \code{genome.oversample} arguments from \code{calcBinnedMotifEnrR}.
#' @param maxFracN The \code{maxFracN} argument from \code{calcBinnedMotifEnrR}.
#' @param GCbreaks The breaks between GC bins. The default value is based on
#'   the hard-coded bins used in Homer.
#' 
#' @return a \code{DataFrame} with sequences represented by rows and columns
#'   \code{seqs}, \code{isForeground}, \code{GCfrac}, \code{GCbin}, \code{GCwgt}
#'   and \code{seqWgt}. Only the first three are already filled in.
#' 
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
#' @importFrom BSgenome getSeq
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom BiocGenerics width
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom stats quantile
#' 
#' @keywords internal
.defineBackground <- function(sqs,
                              bns,
                              bg,
                              currbn,
                              gnm,
                              gnm.regions,
                              gnm.oversample,
                              maxFracN = 0.7,
                              GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                           0.45, 0.5, 0.6, 0.7, 0.8)) {
    
    # calculate G+C frequencies of input sequences
    fmono <- oligonucleotideFrequency(sqs, width = 1, as.prob = TRUE)
    gcf <- fmono[, "G"] + fmono[, "C"]
    
    inCurrBin <- as.integer(bns) == currbn
    
    if (identical(bg, "genome")) {
        
        # (over-)sample random sequences from the genome
        n <- round(gnm.oversample * sum(inCurrBin))
        w <- width(sqs)[inCurrBin]
        mw <- mean(w)
        if (is.null(gnm.regions)) {
            slen <- seqlengths(gnm)
            gnm.regions <- GRanges(seqnames = names(slen),
                                   ranges = IRanges(start = 1, end = slen))
        }
        gnm.regions.width <- width(gnm.regions)
        gnmsqs <- DNAStringSet()
        
        iter <- 0
        while (length(gnmsqs) < n && iter < 10) {
            # sample coordinates
            n1 <- n - length(gnmsqs)
            ridx <- sort(sample(x = length(gnm.regions), size = n1,
                                replace = TRUE,
                                prob = pmax(0, gnm.regions.width - mean(w) + 1)),
                         decreasing = FALSE)
            rw <- sample(x = w, size = n1, replace = TRUE)
            rst <- start(gnm.regions)[ridx] +
                unlist(lapply(pmax(gnm.regions.width[ridx] - rw, 0),
                              function(w1) sample(x = w1, size = 1)))
            reg <- GRanges(seqnames = seqnames(gnm.regions)[ridx],
                           ranges = IRanges(start = rst,
                                            end = pmin(rst + rw - 1,
                                                       end(gnm.regions)[ridx])),
                           seqlengths = seqlengths(gnm))
            
            # extract sequences
            s <- getSeq(gnm, reg)
            
            # filter sequences
            keep <- .filterSeqs(seqs = s, maxFracN = maxFracN)
            
            # add to already sampled sequences
            gnmsqs <- c(gnmsqs, s[keep])
            
            iter <- iter + 1
        }
        names(gnmsqs) <- paste0("g", seq_along(gnmsqs))
        
        # calculate G+C composition of sampled sequences
        fmono.gnm <- oligonucleotideFrequency(gnmsqs, width = 1, as.prob = TRUE)
        gcf.gnm <- fmono.gnm[, "G"] + fmono.gnm[, "C"]
        
        # select a set that matches sqs[inCurrBin]
        sqs.gcbin <- findInterval(x = gcf[inCurrBin], vec = GCbreaks, all.inside = TRUE)
        gnm.gcbin <- findInterval(x = gcf.gnm, vec = GCbreaks, all.inside = TRUE)
        sqs.tab <- tabulate(sqs.gcbin, nbins = length(GCbreaks) - 1) / length(sqs.gcbin)
        gnm.tab <- tabulate(gnm.gcbin, nbins = length(GCbreaks) - 1) / length(gnm.gcbin)
        sel <- sample(x = length(gnmsqs), size = sum(inCurrBin), replace = FALSE,
                      prob = ((sqs.tab + 0.5) / (gnm.tab + 0.5))[gnm.gcbin])
        
        gnm.tab.sel <- tabulate(gnm.gcbin[sel], nbins = length(GCbreaks) - 1) / length(sel)
        gcdist <- sqrt(sum((sqs.tab - gnm.tab.sel)^2))
        if (gcdist > (3 / (length(GCbreaks) - 1))) {
            warning("The background sequences sampled from 'genome' do not match ",
                    "well\nthe G+C content of 'seqs' (normdist = ", round(gcdist, 3),
                    ").\n\n",
                    "Consider increasing 'genome.oversample', and/or focus the ",
                    "sampling on a subset\nwith similar sequence composition as ",
                    "'seqs' using 'genome.regions'.\n\n",
                    "It may also be useful to split 'seqs' into subsets with ",
                    "homogenous G+C\n(e.g. with/without CpG islands) and use ",
                    "an appropriate 'genome' for them\n(e.g. all CpG islands and ",
                    "the rest of the genome, respectively).")
        }
        
        isFg <- rep(c(TRUE, FALSE), each = sum(inCurrBin))
        sqs <- c(sqs[inCurrBin], gnmsqs[sel])
        gcf <- c(gcf[inCurrBin], gcf.gnm[sel])
        
    } else {
        
        if (identical(bg, "otherBins")) {
            isFg <- inCurrBin
            
        } else if (identical(bg, "allBins")) {
            isFg <- rep(c(TRUE, FALSE),
                        c(sum(inCurrBin), length(sqs)))
            sqs <- c(sqs[inCurrBin], sqs)
            gcf <- c(gcf[inCurrBin], gcf)
            
        } else if (identical(bg, "zeroBin")) {
            inZeroBin <- as.integer(bns) == getZeroBin(bns)
            isFg <- rep(c(TRUE, FALSE),
                        c(sum(inCurrBin), sum(inZeroBin)))
            sqs <- c(sqs[inCurrBin], sqs[inZeroBin])
            gcf <- c(gcf[inCurrBin], gcf[inZeroBin])
        }
    } 
    
    # create sequence DataFrame
    df <- DataFrame(seqs = sqs,
                    isForeground = isFg,
                    GCfrac = gcf,
                    GCbin = rep(NA_integer_, length(sqs)),
                    GCwgt = rep(NA_real_, length(sqs)),
                    seqWgt = rep(NA_real_, length(sqs)))
    attr(df, "err") <- NA
    
    # return DataFrame  
    return(df)
}

#' @title Get background sequence weights for GC bins
#'
#' @description The logic is based on Homer (version 4.11). All sequences
#'   binned depending on GC content (\code{GCbreaks}). The numbers of
#'   foreground and background sequences in each bin are counted, and weights
#'   for background sequences in bin i are defined as:
#'   weight_i = (number_fg_seqs_i / number_bg_seqs_i) * (number_bg_seqs_total /
#'   number_fg_seqs_total)
#'
#' @param df a \code{DataFrame} with sequence information.
#' @param GCbreaks The breaks between GC bins. The default value is based on
#'   the hard-coded bins used in Homer.
#' @param verbose A logical scalar. If \code{TRUE}, report on GC weight
#'   calculation.
#'
#' @return a \code{DataFrame} of the same dimensions as the input \code{df},
#'   with the columns \code{GCfrac}, \code{GCbin} and \code{GCwgt}
#'   filled in with the sequence GC content, assigned GC bins and weights to
#'   correct differences in GC distributions between foreground and background
#'   sequences.
#'
#' @importFrom BiocGenerics width
#' @importFrom stats median
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
#' 
#' @keywords internal
.calculateGCweight <- function(df,
                               GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                            0.45, 0.5, 0.6, 0.7, 0.8),
                               verbose = FALSE) {
    
    .checkDfValidity(df)
    GCbreaks <- sort(GCbreaks, decreasing = FALSE)
    .assertVector(x = GCbreaks, type = "numeric", rngIncl = c(0, 1))
    if (length(GCbreaks) < 2) {
        stop("'GCbreaks' must be of length 2 or greater")
    }
    .assertScalar(x = verbose,   type = "logical")
    
    # calculate G+C frequencies
    if (any(is.na(df$GCfrac))) {
        fmono <- oligonucleotideFrequency(df$seqs, width = 1, as.prob = TRUE)
        df$GCfrac <- fmono[, "G"] + fmono[, "C"]
    }
    
    # assign sequences to a GC bin
    df$GCbin <- findInterval(x = df$GCfrac, vec = GCbreaks, all.inside = TRUE)
    
    # keep bins that have at least 1 foreground and 1 background sequence
    # and filter out sequences from unused bins
    used_bins <- sort(intersect(df$GCbin[df$isForeground],
                                df$GCbin[!df$isForeground]))
    keep <- df$GCbin %in% used_bins
    if (verbose) {
        message("  ", length(used_bins), " of ", length(GCbreaks) - 1,
                " GC-bins used (have both fore- and background sequences)\n",
                "  ", sum(!keep), " of ", nrow(df), " sequences (",
                round(100 * sum(!keep) / nrow(df), 1),
                "%) filtered out from unused GC-bins.")
    }
    df <- df[keep, ]
    
    # total number of foreground and background sequences
    total_fg <- sum(df$isForeground)
    total_bg <- sum(!df$isForeground)
    
    # calculate GC weight per bin
    n_fg_b <- table(df$GCbin[df$isForeground])
    n_bg_b <- table(df$GCbin[!df$isForeground])
    weight_per_bin <- (n_fg_b / n_bg_b) * (total_bg / total_fg)
    
    # assign calculated GC weight to each background sequence
    # (foreground sequences get a weight of 1)
    df$GCwgt <- ifelse(df$isForeground,
                       1.0,
                       weight_per_bin[as.character(df$GCbin)])
    
    # return df
    return(df)
}

#' @title Adjust for k-mer composition (single iteration)
#'
#' @description Adjust background sequence weights for differences in k-mer
#'   composition compared to the foreground sequences. This function
#'   implements a single iteration, and is called iteratively by
#'   \code{.iterativeNormForKmers} to get to the final set of adjusted
#'   weights, which will be the result of adjusting for GC and k-mer
#'   composition. The logic is based on Homer's
#'   \code{normalizeSequenceIteration()} function found in \code{Motif2.cpp}.
#'
#' @param kmerFreq a \code{list} with of matrices. The matrix at index \code{i}
#'   in the list contains the probability of k-mers of length \code{i}, for each
#'   k-mer (columns) and sequence (rows).
#' @param goodKmers a \code{list} of \code{numeric} vectors; the element at index
#'   \code{i} contains the number of good (non-N-containing) k-mers of length
#'   \code{i} for each sequence.
#' @param kmerRC a \code{list} of character vectors; the element at index
#'   \code{i} contains the reverse complement sequences of all k-mers of length
#'   \code{i}.
#' @param seqWgt a \code{numeric} vector with starting sequence weights
#'   at the beginning of the iteration.
#' @param isForeground logical vector of the same length as \code{seqs}.
#'   \code{TRUE} indicates that the sequence is from the foreground,
#'   \code{FALSE} that it is a background sequence.
#' @param minSeqWgt Numeric scalar greater than zero giving the
#'   minimal weight of a sequence. The default value (0.001) is based on
#'   \code{Homer} (HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).
#' @param maxSeqWgt Numeric scalar greater than zero giving the
#'   maximal weight of a sequence. The default value (1000) is based on
#'   \code{HOMER} (1 / HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).
#'
#' @return a named \code{list} with elements \code{seqWgt} (updated
#'   weights) and \code{err} (error measuring difference of foreground
#'   and weighted background sequence compositions).
#' 
#' @keywords internal
.normForKmers <- function(kmerFreq,
                          goodKmers,
                          kmerRC,
                          seqWgt,
                          isForeground,
                          minSeqWgt = 0.001,
                          maxSeqWgt = 1000) {
    
    # set starting error
    err <- 0
    
    # iterate over the k-mer sizes
    for (k in seq_along(kmerFreq)) {
        
        # sum k-mer weights for fore- and background sequences
        kmerWgtPerSeq <- kmerFreq[[k]] * seqWgt
        kmerWgtSumForeground <- colSums(kmerWgtPerSeq[isForeground, ])
        kmerWgtSumBackground <- colSums(kmerWgtPerSeq[!isForeground, ])
        
        totalWgtForeground <- sum(kmerWgtSumForeground)
        totalWgtBackground <- sum(kmerWgtSumBackground)
        
        # average weight of a k-mer and its reverse complement
        # (cap at the bottom at 0.5 / total)
        kmerWgtSumForeground <- pmax((kmerWgtSumForeground + kmerWgtSumForeground[kmerRC[[k]]]) / 2,
                                     0.5 / totalWgtForeground)
        kmerWgtSumBackground <- pmax((kmerWgtSumBackground + kmerWgtSumBackground[kmerRC[[k]]]) / 2,
                                     0.5 / totalWgtBackground)
        
        
        # Calculate kmerNormFactors
        kmerNormFactors <- (kmerWgtSumForeground / kmerWgtSumBackground) *
            (totalWgtBackground / totalWgtForeground)
        
        # update error
        err <- err + mean((kmerNormFactors - 1)^2)
        
        # calculate new weights for background sequences
        # ... sum the kmerNormFactors of all k-mers per background sequence
        newSeqWgtBackground <- rowSums(sweep(x = kmerFreq[[k]][!isForeground, ],
                                             MARGIN = 2, STATS = kmerNormFactors, FUN = "*"))
        
        # ... new weight for each background sequence
        seqWgtBackground <- seqWgt[!isForeground]
        newSeqWgtBackground <- newSeqWgtBackground * seqWgtBackground
        
        # ... apply minimum and maximum weights
        newSeqWgtBackground[newSeqWgtBackground < minSeqWgt] <- minSeqWgt
        newSeqWgtBackground[newSeqWgtBackground > maxSeqWgt] <- maxSeqWgt
        
        # ... calculate correction
        penalityBackground <- (ifelse(newSeqWgtBackground > 1,
                                      newSeqWgtBackground,
                                      1 / newSeqWgtBackground))^2
        diffBackground <-  newSeqWgtBackground - seqWgtBackground
        
        toCorrect <- penalityBackground > 1 &
            ((diffBackground > 0 & newSeqWgtBackground > 1) |
                 (diffBackground < 0 & newSeqWgtBackground < 1))
        newSeqWgtBackground[toCorrect] <- seqWgtBackground[toCorrect] +
            diffBackground[toCorrect] / penalityBackground[toCorrect]
        
        # update seqWgt with newSeqWgtBackground
        seqWgt[!isForeground] <- newSeqWgtBackground
    }
    
    # return new seqWgt and error
    list(seqWgt = seqWgt, err = err)
}

#' @title Adjust for k-mer composition (multiple iterations)
#'
#' @description Here we run `.normForKmers` multiple times to converge to
#'   the final weights that will be used to correct the background
#'   sequences for k-mer composition differences compared to the foreground. We
#'   closely follow \code{HOMER}'s \code{normalizeSequence()} function found in
#'   \code{Motif2.cpp}. Note that \code{HOMER} runs the
#'   \code{normalizeSequence()} one last time after going through all iterations
#'   or reaching a low error, which we do not do here.
#'
#' @param df a \code{DataFrame} with sequence information as returned by
#'   \code{.calculateGCweight}.
#' @param maxKmerSize Integer scalar giving the maximum k-mer size to
#'   consider. The default is set to 3 (like in \code{HOMER}), meaning that
#'   k-mers of size 1, 2 and 3 are considered.
#' @param minSeqWgt Numeric scalar greater than zero giving the
#'   minimal weight of a sequence. The default value (0.001) was also used by
#'   \code{HOMER} (HOMER_MINIMUM_SEQ_WEIGHT  constant in Motif2.h).
#' @param maxIter An integer scalar giving the maximum number if
#'   times to run \code{.normForKmers}. the default is set to 160 (as in
#'   \code{HOMER}).
#' @param verbose A logical scalar. If \code{TRUE}, report on k-mer composition
#'   adjustment.
#'
#' @return a DataFrame containing: \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} containing the sequence GC content, GC bins they were
#'   assigned to, the weight to correct for GC differences between foreGround
#'   and background sequences, the weight to adjust for kmer composition, and
#'   the the error term} \item{sequenceNucleotides}{: a \code{DNAStringSet}
#'   object containing the raw sequences} }
#'
#' @importFrom Biostrings oligonucleotideFrequency reverseComplement
#'   DNAStringSet
#' @importFrom S4Vectors DataFrame
#' 
#' @keywords internal
.iterativeNormForKmers <- function(df,
                                   maxKmerSize = 3L,
                                   minSeqWgt = 0.001,
                                   maxIter = 160L,
                                   verbose = FALSE) {
    
    .checkDfValidity(df)
    .assertScalar(x = maxKmerSize, type = "integer", rngIncl = c(1, Inf))
    .assertScalar(x = minSeqWgt, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = maxIter, type = "integer", rngExcl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")
    
    # initialize seqWgt using GCwgt
    lastErr <- Inf
    curWgt <- df$GCwgt
    
    # pre-calculate k-mer frequencies used in iterations
    kmerFreq <- goodKmers <- kmerRC <- list()
    for (k in seq_len(maxKmerSize)) {
        kmerFreq[[k]] <- oligonucleotideFrequency(x = df$seqs, width = k)
        goodKmers[[k]] <- rowSums(kmerFreq[[k]])
        kmerFreq[[k]] <- kmerFreq[[k]] / goodKmers[[k]]
        kmers <- DNAStringSet(colnames(kmerFreq[[k]]))
        kmerRC[[k]] <- as.character(reverseComplement(x = kmers))
    }
    
    # run .normForKmers() up to maxIter times or
    # stop when new error is bigger than the error from the previous iteration
    if (verbose) {
        message("  starting iterative adjustment for k-mer composition (up to ",
                maxIter, " iterations)")
    }
    
    res <- list()
    for (i in seq_len(maxIter)) {
        
        res <- .normForKmers(kmerFreq = kmerFreq,
                             goodKmers = goodKmers,
                             kmerRC = kmerRC,
                             seqWgt = curWgt,
                             isForeground = df$isForeground,
                             minSeqWgt = minSeqWgt)
        
        if (res$err >= lastErr) {
            if (verbose) {
                tmpmsg <- paste0("    detected increasing error - stopping after ", i, " iterations")
                message(tmpmsg)
            }
            break
        } else {
            if (verbose && (i %% 40 == 0)) {
                message("    ", i, " of ", maxIter, " iterations done")
            }
            curWgt <- res$seqWgt
            lastErr <- res$err
        }
    }
    if (verbose) {
        message("    iterations finished")
    }
    
    # return final weights
    df$seqWgt <- curWgt
    attr(df, "err") <- lastErr
    return(df)
}

