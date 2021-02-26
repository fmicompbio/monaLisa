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
#'     \item{\code{gc_frac}}{: the fraction of G+C bases per sequence.}
#'     \item{\code{gc_bin}}{: the GC bin for each sequence.}
#'     \item{\code{GCwgt}}{: the sequence weight to adjust for GC
#'       differences between foreground and background sequences.}
#'     \item{\code{seqWgt}}{: the sequence weight to adjust for k-mer
#'       differences between foreground and background sequences.}
#'   }
#'
#' @return \code{TRUE} (invisibly) if \code{df} is valid, otherwise it
#'   raises an exception using \code{stop()}
.checkDfValidity <- function(df) {
    expected_cols <- c("seqs", "isForeground",
                       "gc_frac", "gc_bin", "GCwgt", "seqWgt")
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

    f1 <- observedFracN > maxFracN
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
        message("  in total filtering out ", sum(res), " of ", length(seqs),
                " sequences (", round(100 * sum(res) / length(seqs), 1), "%)")
    }

    return(res)
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
#'   with the columns \code{gc_frac}, \code{gc_bin} and \code{GCwgt}
#'   filled in with the sequence GC content, assigned GC bins and weights to
#'   correct differences in GC distributions between foreground and background
#'   sequences.
#'
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
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
  
    # calculate GC fraction for each sequence
    fmono <- oligonucleotideFrequency(df$seqs, width = 1, as.prob = TRUE)
    df$gc_frac <- fmono[, "G"] + fmono[, "C"]

    # assign sequences to a GC bin
    df$gc_bin <- findInterval(x = df$gc_frac, vec = GCbreaks, all.inside = TRUE)

    # keep bins that have at least 1 foreground and 1 background sequence
    # and filter out sequences from unused bins
    used_bins <- sort(intersect(df$gc_bin[df$isForeground],
                                df$gc_bin[!df$isForeground]))
    keep <- df$gc_bin %in% used_bins
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
    n_fg_b <- table(df$gc_bin[df$isForeground])
    n_bg_b <- table(df$gc_bin[!df$isForeground])
    weight_per_bin <- (n_fg_b / n_bg_b) * (total_bg / total_fg)

    # assign calculated GC weight to each background sequence
    # (foreground sequences get a weight of 1)
    df$GCwgt <- ifelse(df$isForeground,
                       1.0,
                       weight_per_bin[as.character(df$gc_bin)])

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
#'   \code{get_binned_motif_enrichment()} for each bin.
#'
#'   The defaults type of test is \code{"binomial"}, which is also what
#'   \code{Homer} uses by default. Alternatively, Fisher's exact test can be
#'   used by \code{test = "fisher"}, which has the advantage that special cases
#'   such as zero background counts are handled without ad-hoc adjustments to
#'   the frequencies.
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
.calcMotifEnrichment <- function(motifHitMatrix,
                                 df,
                                 test = c("binomial", "fisher"),
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
        #          TF_hit  not_TF_hit
        #   is_fg     x         y
        #   is_bg     z         w
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
#'   \code{\link[monaLisa]{bin}}.
#' @param pwmL PWMatrixList with motifs for which to calculate enrichments.
#' @param test A \code{character} scalar specifying the type of enrichment test
#'   to perform. One of \code{"binomial"} (default) or \code{"fisher"}. The
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
#' @param Ncpu Number of CPUs to use for parts of the analysis that can run
#'   in parallel (e.g. the scanning for motif hits).
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
#' @param ... Additional arguments for  \code{\link[monaLisa]{findMotifHits}}.
#'
#' @details This function implements a binned motif enrichment analysis. In each
#'   enrichment analysis, the sequences in a specific bin are used as foreground
#'   sequences to test for motif enrichment using sequences in the remaining
#'   bins as background. The logic follows the \code{findMotifsGenome.pl} tool
#'   from \code{Homer} version 4.11, with \code{-size given -nomotif -mknown}
#'   and gives very similar results. As in the \code{Homer} tool, sequences are
#'   weighted to correct for GC and k-mer composition differences between fore-
#'   and background sets.
#'
#'   Motif hits are predicted using \code{\link[monaLisa]{findMotifHits}} and
#'   multiple hits per sequence are counted as just one hit (ZOOPS mode). For
#'   each motif, the weights of sequences that have a hit are summed separately
#'   for foreground (\code{sumForegroundWgtWithHits}) and background
#'   (\code{sumBackgroundWgtWithHits}). The total foreground 
#'   (\code{totalWgtForeground}) and background (\code{totalWgtBackground})
#'   sum of sequence weights is also calculated.
#'
#'   To statistical tests for the calculation of enrichment log p-value are
#'   available: \code{test = "binomial"} (default) to perform binomial test
#'   like \code{Homer}, or \code{test = "fisher"} to perform Fisher's exact
#'   tests, using:
#'   \itemize{
#'     \item{binomial}{: \code{pbinom(q = sumForegroundWgtWithHits - 1, size =
#'       totalWgtForeground, prob = sumBackgroundWgtWithHits / totalWgtBackground,
#'       lower.tail = FALSE, log.p = TRUE)}}
#'     \item{fisher}{: \code{fisher.test(x = tab, alternative =
#'       "greater")}, where \code{tab} is the contingency table with the summed
#'       weights of sequences in foreground or background sets (rows), and with
#'       or without a hit for a particular motif (columns).}
#'   }
#'   
#'   TODO:  - Use Ncpu to also parallelize across list of DFs (enrichment per bin)
#'          - parallellize the other functions using Ncpu.
#'          - change functions to do one-time calculations when using in binned mode?
#'          - filter seqs: make one-time thing.
#'          - think if we want to parametrize some of the pseudo-counts used (leaning towards not)
#'          - change the names of the assays in the SE or return the p-values and FDR (not -log10)
#'          - have a look again at using motif IDs and not names, but keeping both information
#'          - make test="fisher" the default in .calcMotifEnrichment
#'          - better assay names? P and fdr are -log10(...), yes or transform into p-values
#'          - add assay with fraction of sequences containing hits for a given motif
#'          - add motifs without hits to assays (as NA values)

#'
#' @return A \code{SummarizedExperiment} object where the rows are the motifs
#'   and the columns are bins. The four assays are: \itemize{
#'   \item{p}{: -log10 P values}
#'   \item{FDR}{: -log10 false discovery rates}
#'   \item{enr}{: motif enrichments as Pearson residuals}
#'   \item{log2enr}{: motif enrichments as log2 ratios}
#' }
#'
#' @importFrom TFBSTools ID name
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqnames seqlevels
#'
#' @export
get_binned_motif_enrichment <- function(seqs,
                                        bins,
                                        pwmL,
                                        test = c("binomial", "fisher"),
                                        maxFracN = 0.7,
                                        maxKmerSize = 3L,
                                        min.score = 10,
                                        matchMethod = "matchPWM.concat",
                                        Ncpu = 1L,
                                        verbose = FALSE,
                                        ...) {

    # checks
    if (!is(seqs, "DNAStringSet")) {
        stop("class of 'seqs' must be DNAStringSet")
    }
    if (!is(bins, "factor")) {
        stop("'bins' must be of class 'factor'")
    }
    if (length(seqs) != length(bins)) {
        stop("'seqs' and 'bins' must be of equal length and in the same order")
    }
    if (!is(pwmL, "PWMatrixList")) {
        stop("'pwmL' must be of class 'PWMatrixList'")
    }
    .assertScalar(x = Ncpu, type = "integer", rngIncl = c(0,Inf))
    .assertScalar(x = verbose, type = "logical")
    if (is.null(names(seqs))) {
        names(seqs) <- paste0("s", seq_along(seqs))
    }
    if (is.null(names(pwmL))) {
        stop("names(pwmL) is NULL, please name the PWMs, preferably ",
             "with their unique ID.")
    }

    # filter sequences
    if (verbose) {
        message("filtering out bad sequences ...")
    }
    keep <- .filterSeqs(seqs, maxFracN = maxFracN, verbose = verbose)
    seqs <- seqs[keep]
    
    # scan sequences with motif
    if (verbose) {
        message("Scanning sequences for motif hits...")
    }
    hits <- findMotifHits(query = pwmL, subject = seqs, min.score = min.score,
                          method = matchMethod, Ncpu = Ncpu, ...)
    if (isEmpty(hits)) {
        stop("No motif hits found in any of the sequences - aborting.")
    }
    mat <- as.matrix(as.data.frame.matrix(table(seqnames(hits), as.character(hits$pwmid))))
    # ... add missing rows (sequences that had no hit)
    missing_row_names <- names(seqs)[!names(seqs) %in% rownames(mat)]
    missing_mat <- matrix(data = 0, nrow = length(missing_row_names), ncol = ncol(mat))
    rownames(missing_mat) <- missing_row_names
    colnames(missing_mat) <- colnames(mat)
    complete_hit_mat <- rbind(mat, missing_mat)
    complete_hit_mat <- complete_hit_mat[names(seqs), , drop = FALSE]
    # ... ZOOPS mode
    w <- complete_hit_mat > 1
    complete_hit_mat[w] <- 1

    # create list of DataFrames, one for each bin
    bin_levels <- levels(bins)
    DF_list <- list()
    for (i in seq_along(bin_levels)) {
        isForeground <- logical(length = length(seqs))
        isForeground[bins == bin_levels[i]] <- TRUE
        df <- DataFrame(seqs = seqs,
                        isForeground = isForeground,
                        gc_frac = NA_real_,
                        gc_bin = NA_integer_,
                        GCwgt = NA_real_,
                        seqWgt = NA_real_)
        attr(df, "err") <- NA

        DF_list[[i]] <- df
        names(DF_list)[i] <- bin_levels[i]
    }

    # calculate weight to adjust for GC differences between foreground and background per bin
    # ... in this step, sequences may be filtered out (if a GC bin contains one sequence only, that sequence is filtered out).
    # ... Since the kept seqs may differ per bin, we select the kept seqs at the enrichment per bin.
    if (verbose) {
        message("Correcting for GC differences to the background sequences per bin ...")
    }
    DF_list <- lapply(DF_list, function(x) {.calculateGCweight(df = x, verbose = verbose)})

    # update weight to in addition adjust for kmer composition differences between foreground and background per bin
    if (verbose) {
        message("Correcting for kmer differences to the background sequences per bin ...")
    }
    DF_list <- lapply(DF_list, function(x) {.iterativeNormForKmers(df = x, maxKmerSize = maxKmerSize, verbose = verbose)})

    # get log(p-values) for motif enrichment
    if (verbose) {
        message("Calculating motif enrichment per bin ...")
    }
    enrich_list <- lapply(DF_list, function(x) {
        .calcMotifEnrichment(motifHitMatrix = complete_hit_mat[rownames(x), , drop = FALSE],
                             df = x, test = test, verbose = verbose)
    })

    # summarize results to return as SE

    # - log10 of p-value
    P <- do.call(cbind, lapply(enrich_list, function(bin_res) {
        D <- bin_res[, "logP"]
        logpVals <- -D / log(10)
        names(logpVals) <- bin_res[, "motifName"]
        logpVals
    }))

    # - log10 of FDR
    tmp <-  as.vector(10**(-P))
    fdr <- matrix(-log10(p.adjust(tmp, method = "BH")), nrow = nrow(P)) # add as parameter? the method of choice for multiple testing correction?
    dimnames(fdr) <- dimnames(P)
    fdr[which(fdr == Inf, arr.ind = TRUE)] <- max(fdr[is.finite(fdr)])

    # enrTF (Pearson residuals)
    enrTF <- do.call(cbind, lapply(enrich_list, function(bin_res) {
        D <- as.matrix(
            data.frame(
              sumForegroundWgtWithHits = bin_res[, "sumForegroundWgtWithHits"],
              frac_fg_seq_with_motif = (bin_res[, "sumForegroundWgtWithHits"] / bin_res[, "totalWgtForeground"]),
              frac_bg_seq_with_motif = (bin_res[, "sumBackgroundWgtWithHits"] / bin_res[, "totalWgtBackground"]))
        )
        rownames(D) <-  bin_res[, "motifName"]
        obsTF <- D[, "sumForegroundWgtWithHits"]
        expTF <- D[, "sumForegroundWgtWithHits"] / (D[, "frac_fg_seq_with_motif"] + 0.001) * (D[, "frac_bg_seq_with_motif"] + 0.001) # keep pseudo count fixed or make parameter? exp = numb_fg_total*bg_frac
        enr <- (obsTF - expTF) / sqrt(expTF)
        enr[ is.na(enr) ] <- 0
        enr
    }))

    # log2enr
    log2enr <- do.call(cbind, lapply(enrich_list, function(bin_res) {
        D <- bin_res[, c("sumForegroundWgtWithHits", "sumBackgroundWgtWithHits")]
        nTot <- c(bin_res[1, "totalWgtForeground"], bin_res[1, "totalWgtBackground"]) # Do I round the bg total to the nearest integer?
        D.norm <- t(min(nTot) * t(D) / nTot) # scale to smaller number (usually number of target sequences) #
        DL <- log2(D.norm + 8) # keep pseudo count fixed? --> yes
        log2enr <- DL[, 1] - DL[, 2]
        names(log2enr) <- bin_res[, "motifName"]
        log2enr
    }))

    # return SummarizedExperiment
    se <- SummarizedExperiment(assays = list(p = P, FDR = fdr, enr = enrTF,
                                             log2enr = log2enr))
    percentGC <- unlist(lapply(pwmL, function(x) {
      m <- TFBSTools::Matrix(x)
      m <- 2^m * 0.25
      100 * sum(m[c("C","G"), ]) / sum(m)
    }), use.names = FALSE)
    rowData(se) <- S4Vectors::DataFrame(motif.id = TFBSTools::ID(pwmL),
                                        motif.name = TFBSTools::name(pwmL),
                                        motif.pwm = pwmL,
                                        motif.percentGC = percentGC)
    metadata(se) <- list(params = list(test = test, maxFracN = maxFracN,
                                       maxKmerSize = maxKmerSize,
                                       min.score = min.score,
                                       matchMethod = matchMethod,
                                       Ncpu = Ncpu, verbose = verbose))

    return(se)
}
