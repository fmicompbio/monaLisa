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
#' @param gnm.seed An \code{integer} scalar to seed the random number generator.
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
                              gnm.seed,
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
        if (is.null(gnm.regions)) {
            slen <- seqlengths(gnm)
            gnm.regions <- GRanges(seqnames = names(slen),
                                   ranges = IRanges(start = 1, end = slen))
        }
        gnm.regions.width <- width(gnm.regions)
        set.seed(gnm.seed)
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
            ridxL <- split(seq_along(ridx), ridx)[as.character(unique(ridx))]
            reg <- GRanges(seqnames = seqnames(gnm.regions)[ridx],
                           ranges = do.call(c,
                              unname(lapply(ridxL, function(i) {
                  rst <- sample(x = gnm.regions.width[ridx[i[1]]] - mean(w) + 1,
                                size = length(i), replace = FALSE)
                  IRanges(start = rst,
                          end = pmin(rst + rw[i] - 1, gnm.regions.width[ridx[i[1]]]))
            }))), seqlengths = seqlengths(gnm))

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
            inZeroBin <- as.integer(bns) == attr(bns, "bin0")
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
#' @param normalizeByLength A logical scalar. If \code{TRUE}, the weight calculated
#'   for each background sequence in a specific GC bin, to account for GC differences 
#'   between foreground and background, is multiplied by the length of the background
#'   sequence divided by the median length of foreground sequences in that GC bin.
#' @param verbose A logical scalar. If \code{TRUE}, report on GC weight
#'   calculation.
#'
#' @return a \code{DataFrame} of the same dimensions as the input \code{df},
#'   with the columns \code{GCfrac}, \code{GCbin} and \code{GCwgt}
#'   filled in with the sequence GC content, assigned GC bins and weights to
#'   correct differences in GC distributions between foreground and background
#'   sequences.
#'
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
#' 
#' @keywords internal
.calculateGCweight <- function(df,
                               GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                            0.45, 0.5, 0.6, 0.7, 0.8),
                               normalizeByLength = TRUE, 
                               verbose = FALSE) {

    .checkDfValidity(df)
    GCbreaks <- sort(GCbreaks, decreasing = FALSE)
    .assertVector(x = GCbreaks, type = "numeric", rngIncl = c(0, 1))
    if (length(GCbreaks) < 2) {
        stop("'GCbreaks' must be of length 2 or greater")
    }
    .assertScalar(x = normalizeByLength, type = "logical")
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
    
    # adjust for sequence lengths (per GC bin)
    if (normalizeByLength) {
        # ... get sequence lengths
        seq_lengths <- width(df$seqs)
        # ... get median length of foreground sequences per GC bin
        median_FG_length_per_GCbin <- vapply(
            X = as.character(names(weight_per_bin)), 
            FUN = function(x){
                median(seq_lengths[df$isForeground & as.character(df$GCbin) == x])
            }, 
            FUN.VALUE = 0)
        median_FG_length_per_GCbin_vector <- 
            ifelse(df$isForeground,
                   NA,
                   median_FG_length_per_GCbin[as.character(df$GCbin)])
        # ... correct background sequence weights
        w_bg <- !df$isForeground
        df$GCwgt[w_bg] <- df$GCwgt[w_bg] * seq_lengths[w_bg] / median_FG_length_per_GCbin_vector[w_bg]
    }
    
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
          message("    detected increasing error - stopping after ", i, " iterations")
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
#'   the advantage that special cases such as zero background counts are handled 
#'   without ad-hoc adjustments to the frequencies.
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
#'     \item{sumForegroundWgtWithHits}{: the sum of the weights of the foreground
#'        sequences that have at least one instance of a specific motif
#'        (ZOOPS mode).}
#'     \item{sumBackgroundWgtWithHits}{: the sum of the weights of the background
#'        sequences that have at least one instance of a specific motif
#'        (ZOOPS mode).}
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
    method <- match.arg(test)
    .assertScalar(x = verbose, type = "logical")

    totalWgtForeground <- sum(df$seqWgt[df$isForeground])
    totalWgtBackground <- sum(df$seqWgt[!df$isForeground])

    motifHitMatrixWeighted <- motifHitMatrix * df$seqWgt
    TFmatchedSeqCountForeground <- colSums(motifHitMatrixWeighted[df$isForeground, ])
    TFmatchedSeqCountBackground <- colSums(motifHitMatrixWeighted[!df$isForeground, ])

    # calculate motif enrichment
    if (method == "binomial") {

      if (verbose) {
          message("using binomial test to calculate ",
                  "log(p-values) for motif enrichments")
      }

      prob <- TFmatchedSeqCountBackground / totalWgtBackground
      minProb <- 1 / totalWgtBackground
      maxProb <- (totalWgtBackground - 1) / totalWgtBackground
      if (any(i <- (prob < minProb))) {
          # warning("some background TF match probabilities are below ",
          #         "minProb (for example when there were zero hits) ",
          #         "and will be given a value of minProb=1/totalWgtBackground")
          prob[i] <- minProb
      }
      if (any(i <- (prob > maxProb))) {
          # warning("some TF match probabilities a above",
          #         "maxProb (for example when all sequences had hits) ",
          #         "and will be givena value of ",
          #         "maxProb=(totalWgtBackground-1)/totalWgtBackground")
          prob[i] <- maxProb
      }

      logP <- pbinom(q = TFmatchedSeqCountForeground - 1,
                     size = totalWgtForeground,
                     prob = prob, lower.tail = FALSE, log.p = TRUE)

    } else if (method == "fisher") {

        if (verbose) {
            message("using fisher's exact test (one-sided) to calculate ",
                    "log(p-values) for motif enrichments")
        }

        # contingency table per motif for fisher's exact test (rounded to integer):
        #              withHit  noHit
        #   foreground    x       y
        #   background    z       w
        #
        logP <- log(vapply(structure(seq_along(TFmatchedSeqCountForeground),
                                     names = names(TFmatchedSeqCountForeground)),
                               function(i) {
            ctab <- rbind(c(TFmatchedSeqCountForeground[i],
                            totalWgtForeground - TFmatchedSeqCountForeground[i]),
                          c(TFmatchedSeqCountBackground[i],
                            totalWgtBackground - TFmatchedSeqCountBackground[i]))
            ctab <- round(ctab)
            fisher.test(x = ctab, alternative = "greater")$p.value
        }, FUN.VALUE = numeric(1)))
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
#' @param seqs DNAStringSet object with sequences to test
#' @param bins factor of the same length and order as \code{seqs}, indicating
#'   the bin for each sequence. Typically the return value of
#'   \code{\link[monaLisa]{bin}}. For \code{background = "genome"}, \code{bins}
#'   can be omitted.
#' @param pwmL PWMatrixList with motifs for which to calculate enrichments.
#' @param background A \code{character} scalar specifying the background sequences
#'   to use. One of \code{"otherBins"} (default), \code{"allBins"}, \code{"zeroBin"}
#'   or \code{"genome"} (see "Details").
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
#' @param normalizeByLength A logical scalar. If \code{TRUE}, account for 
#'   sequence length differences between foreground and background sequences
#'   (see Details).
#' @param pseudocount.log2enr A numerical scalar with the pseudocount to add to
#'   foreground and background counts when calculating log2 motif enrichments
#' @param pseudocount.pearsonResid A numerical scalar with the pseudocount to add
#'   to foreground and background frequencies when calculating expected counts
#'   and Pearson residuals.
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
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
#' @param ... Additional arguments for  \code{\link[monaLisa]{findMotifHits}}.
#'
#' @details This function implements a binned motif enrichment analysis. In each
#'   enrichment analysis, the sequences in a specific bin are used as foreground
#'   sequences to test for motif enrichment comparing to background sequences
#'   (defined by \code{background}, see below). The logic follows the
#'   \code{findMotifsGenome.pl} tool from \code{Homer} version 4.11, with
#'   \code{-size given -nomotif -mknown} and additionally \code{-h} if using 
#'   \code{test = "fisher"}, and gives very similar results when 
#'   \code{normalizeByLength = FALSE}. As in the \code{Homer} tool, sequences 
#'   are weighted to correct for GC and k-mer composition differences between 
#'   fore- and background sets. With \code{normalizeByLength = TRUE}, the sequence
#'   weights are additionally corrected for length differences.
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
#'   }
#'
#'   Motif hits are predicted using \code{\link[monaLisa]{findMotifHits}} and
#'   multiple hits per sequence are counted as just one hit (ZOOPS mode). For
#'   each motif, the weights of sequences that have a hit are summed separately
#'   for foreground (\code{sumForegroundWgtWithHits}) and background
#'   (\code{sumBackgroundWgtWithHits}). The total foreground 
#'   (\code{totalWgtForeground}) and background (\code{totalWgtBackground})
#'   sum of sequence weights is also calculated. If a motif has zero 
#'   (\code{sumForegroundWgtWithHits}) and (\code{sumBackgroundWgtWithHits}), 
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
#'       totalWgtForeground, prob = sumBackgroundWgtWithHits / totalWgtBackground,
#'       lower.tail = FALSE, log.p = TRUE)}}
#'   }
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
                                background = c("otherBins", "allBins", "zeroBin", "genome"),
                                test = c("fisher", "binomial"),
                                maxFracN = 0.7,
                                maxKmerSize = 3L,
                                min.score = 10,
                                matchMethod = "matchPWM",
                                GCbreaks = c(0.2, 0.25, 0.3, 0.35, 0.4,
                                             0.45, 0.5, 0.6, 0.7, 0.8),
                                normalizeByLength = TRUE, 
                                pseudocount.log2enr = 8,
                                pseudocount.pearsonResid = 0.001,
                                p.adjust.method = "BH",
                                genome = NULL,
                                genome.regions = NULL,
                                genome.oversample = 2,
                                genome.seed = 42L,
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
    .assertScalar(x = pseudocount.log2enr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pseudocount.pearsonResid, type = "numeric", rngIncl = c(0, Inf))
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

    
    # filter sequences
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
    hitmatrix <- unclass(table(factor(seqnames(hits), levels = seqlevels(hits)),
                               factor(hits$pwmid, levels = TFBSTools::ID(pwmL))))
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
                                gnm.seed = genome.seed + i,
                                maxFracN = maxFracN)

        # for background = "genome", scan sampled background sequences for motifs
        if (identical(background, "genome")) {
            if (verbose1) {
              message("Scanning genomic background sequences for motif hits...")
            }
            hits.genome <- findMotifHits(query = pwmL, subject = df$seqs[!df$isForeground],
                                         min.score = min.score, method = matchMethod,
                                         BPPARAM = BPPARAM, ...)
            if (isEmpty(hits.genome)) {
              stop("No motif hits found in any of the genomic background sequences - aborting.")
            }

            # create motif hit matrix
            hitmatrix.genome <- unclass(table(factor(seqnames(hits.genome),
                                                     levels = seqlevels(hits.genome)),
                                        factor(hits.genome$pwmid,
                                               levels = TFBSTools::ID(pwmL))))
            # zoops (zero or one per sequence) mode
            hitmatrix.genome[hitmatrix.genome > 0] <- 1
            
            hitmatrix2 <- rbind(hitmatrix, hitmatrix.genome)

        } else {
            hitmatrix2 <- hitmatrix
        }
        
        # calculate initial background sequence weights based on G+C composition
        if (verbose1) {
            message("Correcting for GC differences to the background sequences...")
        }
        df <- .calculateGCweight(df = df,
                                 GCbreaks = GCbreaks,
                                 normalizeByLength = normalizeByLength, 
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
            message("Calculating motif enrichment...")
        }
        enrich1 <- .calcMotifEnrichment(motifHitMatrix = hitmatrix2[rownames(df), ],
                                        df = df, test = test, verbose = verbose1)
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
    padj <- matrix(-log10(p.adjust(as.vector(10**(-P)), method = p.adjust.method)), nrow = nrow(P))
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
                 param = list(method = "R",
                              background = background,
                              test = test,
                              maxFracN = maxFracN,
                              maxKmerSize = maxKmerSize,
                              min.score = min.score,
                              matchMethod = matchMethod,
                              pseudocount.log2enr = pseudocount.log2enr,
                              pseudocount.pearsonResid = pseudocount.pearsonResid,
                              p.adj.method = p.adjust.method,
                              genome.class = class(genome),
                              genome.regions = genome.regions,
                              genome.oversample = 2,
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

    return(se)
}
