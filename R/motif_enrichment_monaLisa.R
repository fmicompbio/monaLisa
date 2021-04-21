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
        message("  in total filtering out ", sum(!res), " of ", length(seqs),
                " sequences (", round(100 * sum(!res) / length(seqs), 1), "%)")
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
    df$GCfrac <- fmono[, "G"] + fmono[, "C"]

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
#' 
#' @keywords internal
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
#' @param bins factor of the 
#' same length and order as \code{seqs}, indicating
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
#' @param pseudocount.log2enr A numerical scalar with the pseudocount to add to
#'   foreground and background counts when calculating log2 motif enrichments
#' @param pseudocount.pearsonResid A numerical scalar with the pseudocount to add
#'   to foreground and background frequencies when calculating expected counts
#'   and Pearson residuals.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation.
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
                                bins,
                                pwmL,
                                test = c("binomial", "fisher"),
                                maxFracN = 0.7,
                                maxKmerSize = 3L,
                                min.score = 10,
                                matchMethod = "matchPWM",
                                pseudocount.log2enr = 8,
                                pseudocount.pearsonResid = 0.001,
                                p.adjust.method = "BH",
                                BPPARAM = SerialParam(),
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
    .assertScalar(x = pseudocount.log2enr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pseudocount.pearsonResid, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)
    if (!is(BPPARAM, "BiocParallelParam")) {
        stop("'BPPARAM' must be of class 'BiocParallelParam'")
    }
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

        # create sequence info data frame
        df <- DataFrame(seqs = seqs,
                        isForeground = (as.integer(bins) == i),
                        GCfrac = rep(NA_real_, length(seqs)),
                        GCbin = rep(NA_integer_, length(seqs)),
                        GCwgt = rep(NA_real_, length(seqs)),
                        seqWgt = rep(NA_real_, length(seqs)))
        attr(df, "err") <- NA

        # calculate initial background sequence weights based on G+C composition
        if (verbose1) {
            message("Correcting for GC differences to the background sequences...")
        }
        df <- .calculateGCweight(df = df,
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
            message("Calculating motif enrichment...")
        }
        enrich1 <- .calcMotifEnrichment(motifHitMatrix = hitmatrix[rownames(df), , drop = FALSE],
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
    mdat <- list(sequences = seqs,
                 bins = bins,
                 bins.binmode = attr(bins, "binmode"),
                 bins.breaks = as.vector(attr(bins, "breaks")),
                 bins.bin0 = attr(bins, "bin0"),
                 param.test = test,
                 param.maxFracN = maxFracN,
                 param.maxKmerSize = maxKmerSize,
                 param.min.score = min.score,
                 param.matchMethod = matchMethod,
                 param.pseudocount.log2enr = pseudocount.log2enr,
                 param.pseudocount.pearsonResid = pseudocount.pearsonResid,
                 param.p.adj.method = p.adjust.method,
                 param.BPPARAM.class = class(BPPARAM),
                 param.BPPARAM.bpnworkers = bpnworkers(BPPARAM),
                 param.verbose = verbose)
    assaySumForegroundWgtWithHits <- do.call(cbind, lapply(enrichL, function(x){x$sumForegroundWgtWithHits}))
    assaySumBackgroundWgtWithHits <- do.call(cbind, lapply(enrichL, function(x){x$sumBackgroundWgtWithHits}))
    se <- SummarizedExperiment(assays = list(negLog10P = P, 
                                             negLog10Padj = padj, 
                                             pearsonResid = enrTF,
                                             log2enr = log2enr, 
                                             sumForegroundWgtWithHits = assaySumForegroundWgtWithHits, 
                                             sumBackgroundWgtWithHits = assaySumBackgroundWgtWithHits),
                               rowData = rdat, colData = cdat, metadata = mdat)

    return(se)
}
