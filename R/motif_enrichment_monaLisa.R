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
#'     \item{\code{is_foreground}}{ that indicates if a sequence is in the
#'                                 foreground group.}
#'     \item{\code{gc_frac}}{: the fraction of G+C bases per sequence.}
#'     \item{\code{gc_bin}}{: the GC bin for each sequence.}
#'     \item{\code{gc_weight}}{: the sequence weight to adjust for GC
#'       differences between foreground and background sequences.}
#'     \item{\code{kmer_weight}}{: the sequence weight to adjust for k-mer
#'       differences between foreground and background sequences.}
#'   }
#'
#' @return \code{TRUE} (invisibly) if \code{df} is valid, otherwise it
#'   raises an exception using \code{stop()}
.checkDfValidity <- function(df) {
    expected_cols <- c("seqs", "is_foreground",
                       "gc_frac", "gc_bin", "gc_weight", "kmer_weight")
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
#' @param df a \code{DataFrame} with sequence information.
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
#' @return the filtered \code{df}.
#'
#' @importFrom Biostrings alphabetFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
.filterSeqs <- function(df, maxFracN = 0.7, minLength = 5L,
                        maxLength = 100000L, verbose = FALSE) {

    .checkDfValidity(df)
    if (!is.numeric(maxFracN) || length(maxFracN) != 1L || maxFracN < 0 || maxFracN > 1) {
        stop("'maxFracN' has to be a numerical scalar with a value in [0,1]")
    }
    if (!is.numeric(minLength) || length(minLength) != 1L || minLength < 0) {
        stop("'minLength' has to be a non-negative numerical scalar")
    }
    if (!is.numeric(maxLength) || length(maxLength) != 1L || maxLength < max(minLength, 0)) {
        stop("'maxLength' has to be a non-negative numerical scalar",
             " greater than 'minLength' (", minLength, ")")
    }
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("'verbose' has to be either TRUE or FALSE")
    }

    # fraction of N bases per sequence
    observedFracN <- alphabetFrequency(df$seqs, as.prob = TRUE)[, "N"]

    # remove sequences with observedFracN > maxFracN
    w <- which(observedFracN > maxFracN)
    if (length(w) > 0) {
        if (verbose) {
            message("  filtering out ", length(w), " of ", nrow(df),
                    " sequences (", round(100 * length(w) / nrow(df), 1), "%)",
                    " with too many N bases")
        }
        df <- df[-w, ]
    }

    # remove too long/short sequences
    w <- which(width(df$seqs) < minLength | width(df$seqs) > maxLength)
    if (length(w) > 0) {
        if (verbose) {
            message("  filtering out ", length(w), " of ", nrow(df),
                    " sequences (", round(100 * length(w) / nrow(df), 1), "%)",
                    " that are too short or too long")
        }
        df <- df[-w, ]
    }

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
#'   with the columns \code{gc_frac}, \code{gc_bin} and \code{gc_weight}
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
    if (!is.numeric(GCbreaks) || length(GCbreaks) < 2 ||
        any(diff(GCbreaks) <= 0) || any(GCbreaks < 0 || GCbreaks > 1)) {
        stop("'GCbreaks' have to be a strictly increasing numerical vector of ",
             "at least length two and with values in (0, 1).")
    }
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("'verbose' has to be either TRUE or FALSE")
    }
    
    # calculate GC fraction for each sequence
    fmono <- oligonucleotideFrequency(df$seqs, width = 1, as.prob = TRUE)
    df$gc_frac <- fmono[, "G"] + fmono[, "C"]
  
    # assign sequences to a GC bin
    df$gc_bin <- findInterval(x = df$gc_frac, vec = GCbreaks, all.inside = TRUE)
  
    # keep bins that have at least 1 foreground and 1 background sequence
    # and filter out sequences from unused bins
    used_bins <- sort(intersect(df$gc_bin[df$is_foreground],
                                df$gc_bin[!df$is_foreground]))
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
    total_fg <- sum(df$is_foreground)
    total_bg <- sum(!df$is_foreground)
  
    # calculate GC weight per bin
    n_fg_b <- table(df$gc_bin[df$is_foreground])
    n_bg_b <- table(df$gc_bin[!df$is_foreground])
    weight_per_bin <- (n_fg_b / n_bg_b) * (total_bg / total_fg)
  
    # assign calculated GC weight to each background sequence
    # (foreground sequences get a weight of 1)
    df$gc_weight <- ifelse(df$is_foreground,
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
#' @param kmer_freq a \code{list} with of matrices. The matrix at index \code{i}
#'   in the list contains the counts of k-mers of length \code{i}, for each
#'   k-mer (columns) and sequence (rows).
#' @param g_oligos a \code{list} of \code{numeric} vectors; the element at index
#'   \code{i} contains the number of good (non-N-containing) k-mers of length
#'   \code{i} for each sequence.
#' @param kmer_seq_rc a \code{list} of character vectors; the element at index
#'   \code{i} contains the reverse complement sequences of all k-mers of length
#'   \code{i}.
#' @param kmer_weight a \code{numeric} vector with starting sequence weights
#'   at the beginning of the iteration.
#' @param is_foreground logical vector of the same length as \code{seqs}.
#'   \code{TRUE} indicates that the sequence is from the foreground,
#'   \code{FALSE} that it is a background sequence.
#' @param minimum_seq_weight Numeric scalar greater than zero giving the
#'   minimal weight of a sequence. The default value (0.001) is based on
#'   \code{Homer} (HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).
#' @param maximum_seq_weight Numeric scalar greater than zero giving the
#'   maximal weight of a sequence. The default value (1000) is based on
#'   \code{HOMER} (1 / HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).
#'
#' @return a named \code{list} with elements \code{kmer_weight} (updated
#'   weights) and \code{err} (error measuring difference of foreground
#'   and weighted background sequence compositions).
.normForKmers <- function(kmer_freq,
                          g_oligos,
                          kmer_seq_rc,
                          kmer_weight,
                          is_foreground,
                          minimum_seq_weight = 0.001,
                          maximum_seq_weight = 1000) {

    # set starting error
    err <- 0
  
    # iterate over the k-mer sizes
    for (k in seq_along(kmer_freq)) {
  
        # divide the current weight of each sequence by its g_oligos
        div_weight <- kmer_weight / g_oligos[[k]]
    
        # for each sequence multiply the frequency of each oligo by its div_weight
        # (distributing the weights of each sequence to its oligos)
        oligo_weights_per_seq <- kmer_freq[[k]] * div_weight
    
        # sum weights per oligo over foreground (foreground_levels) and background
        # (background_levels) sequences
        foreground_levels <- colSums(oligo_weights_per_seq[is_foreground, ])
        background_levels <- colSums(oligo_weights_per_seq[!is_foreground, ])
    
        # sum weights in foreground_levels and background_levels
        total_foreground <- sum(foreground_levels)
        total_background <- sum(background_levels)
    
        # min values as used by Homer
        min_foreground_levels <- 0.5 / total_foreground
        min_background_levels <- 0.5 / total_background
    
        # Average the weight of a k-mer with its reverse complement
        f_level <- (foreground_levels + foreground_levels[kmer_seq_rc[[k]]]) / 2
        b_level <- (background_levels + background_levels[kmer_seq_rc[[k]]]) / 2
    
        # check if less than set minimum
        f_level[f_level < min_foreground_levels] <- min_foreground_levels
        b_level[b_level < min_background_levels] <- min_background_levels
    
        # Calculate normFactor (to be used to correct background sequences)
        norm_factors <- (f_level / b_level) * (total_background / total_foreground)
    
        # update error
        err <- err + sum((norm_factors - 1)^2 / length(foreground_levels))
    
        # calculate new weights for background sequences
    
        # ... sum the norm_factors of all k-mers per background sequence
        bg_new_score <- rowSums(sweep(x = kmer_freq[[k]][!is_foreground, ],
                                      MARGIN = 2, STATS = norm_factors, FUN = "*"))
    
        # ... Homer check: if number of good oligos is > 0.5
        bg_g_oligos <- g_oligos[[k]][!is_foreground]
        g <- bg_g_oligos > 0.5
        bg_new_score[g] <- bg_new_score[g] / bg_g_oligos[g]
    
        # ... new weight for each background sequence
        # ... ... newWeight = newScore*currentWeight
        bg_cur_weight <- kmer_weight[!is_foreground]
        bg_new_weight <- bg_new_score * bg_cur_weight
    
        # ... Homer's minimum and maximum weight cutoffs
        bg_new_weight[bg_new_weight < minimum_seq_weight] <- minimum_seq_weight
        bg_new_weight[bg_new_weight > maximum_seq_weight] <- maximum_seq_weight
    
        # ... penalty (as used by Homer)
        bg_penalty <- bg_new_weight
        g <- bg_penalty < 1
        bg_penalty[g] <- 1 / bg_penalty[g]
        bg_penalty <- bg_penalty^2
    
        # ... delta (as used by Homer)
        bg_delta <-  bg_new_weight - bg_cur_weight
    
        # ... bg_new_weight1 (as used by Homer)
        bg_new_weight1 <- bg_cur_weight + bg_delta
        g <- bg_penalty > 1 & ((bg_delta > 0 & bg_new_weight > 1) |
                                 (bg_delta < 0 & bg_new_weight < 1))
        bg_new_weight1[g] <- bg_cur_weight[g] + bg_delta[g] / bg_penalty[g]
    
        # ... update kmer_weight with bg_new_weight1
        kmer_weight[!is_foreground] <- bg_new_weight1
    }
  
    # return new kmer_weight and error
    list(kmer_weight = kmer_weight, err = err)
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
#' @param max_kmer_size Integer scalar giving the maximum k-mer size to
#'   consider. The default is set to 3 (like in \code{HOMER}), meaning that
#'   k-mers of size 1, 2 and 3 are considered.
#' @param minimum_seq_weight Numeric scalar greater than zero giving the
#'   minimal weight of a sequence. The default value (0.001) was also used by
#'   \code{HOMER} (HOMER_MINIMUM_SEQ_WEIGHT  constant in Motif2.h).
#' @param max_autonorm_iters An integer scalar giving the maximum number if
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
                                   max_kmer_size = 3L,
                                   minimum_seq_weight = 0.001,
                                   max_autonorm_iters = 160L,
                                   verbose = FALSE) {

  .checkDfValidity(df)
  if (!is.integer(max_kmer_size) ||
      length(max_kmer_size) != 1L ||
      max_kmer_size < 1) {
    stop("'max_kmer_size' must be an integer scalar greater than zero")
  }
  if (!is.numeric(minimum_seq_weight) ||
      length(minimum_seq_weight) != 1L ||
      minimum_seq_weight <= 0) {
    stop("'minimum_seq_weight' must be a numeric scalar greater than zero")
  }
  if (!is.integer(max_autonorm_iters) ||
      length(max_autonorm_iters) != 1L ||
      max_autonorm_iters < 1) {
    stop("'max_autonorm_iters' must be an integer scalar greater than zero")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' has to be either TRUE or FALSE")
  }
  
  # initialize kmer_weight using gc_weight
  last_error <- Inf
  cur_weight <- df$gc_weight

  # pre-calculate k-mer frequencies used in iterations
  kmer_freq <- g_oligos <- kmer_seq_rc <- list()
  for (k in seq_len(max_kmer_size)) {
    kmer_freq[[k]] <- oligonucleotideFrequency(x = df$seqs, width = k)
    g_oligos[[k]] <- rowSums(kmer_freq[[k]])
    kmers <- DNAStringSet(colnames(kmer_freq[[k]]))
    kmer_seq_rc[[k]] <- as.character(reverseComplement(x = kmers))
  }

  # run .normForKmers() up to max_autonorm_iters times or
  # stop when new error is bigger than the error from the previous iteration
  if (verbose) {
    message("  starting iterative adjustment for k-mer composition (up to ",
            max_autonorm_iters, " iterations)")
  }

  res <- list()
  for (i in seq_len(max_autonorm_iters)) {

    # run .normForKmers
    res <- .normForKmers(kmer_freq = kmer_freq,
                         g_oligos = g_oligos,
                         kmer_seq_rc = kmer_seq_rc,
                         kmer_weight = cur_weight,
                         is_foreground = df$is_foreground,
                         minimum_seq_weight = minimum_seq_weight)

    # if current error is bigger than the last one, stop
    if (res$err >= last_error) {
      if (verbose) {
        message("    detected increasing error - stopping after ", i, " iterations")
      }
      break
    } else {
      if (verbose && (i %% 40 == 0)) {
        message("    ", i, " of ", max_autonorm_iters, " iterations done")
      }
      cur_weight <- res$kmer_weight
      last_error <- res$err
    }
  }
  if (verbose) {
    message("    iterations finished")
  }

  # return final weights
  df$kmer_weight <- cur_weight
  attr(df, "err") <- last_error
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
#' @param motif_matrix matrix with 0 and 1 entries for absence or presence of 
#'   motif hits in each sequence.
#' @param df a \code{DataFrame} with sequence information as returned by
#'   \code{.iterativeNormForKmers()}.
#' @param test type of motif enrichment test to perform. 
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#' 
#' @return a \code{data.frame} containing the motifs as rows and the columns:
#'   \itemize{
#'     \item{motif_name}{: the motif name}
#'     \item{log_p_value}{: the log p-value for enrichment (natural logarithm).
#'        If \code{test="binomial"} (default), this log p-value is identical to
#'        the one returned by Homer.}
#'     \item{fg_weight_sum}{: the sum of the weights of the foreground sequences
#'        that have at least one instance of a specific motif (ZOOPS mode).}
#'     \item{bg_weight_sum}{: the sum of the weights of the background sequences
#'        that have at least one instance of a specific motif (ZOOPS mode).}
#'     \item{fg_weight_sum_total}{: the total sum of weights of foreground
#'        sequences.}
#'     \item{bg_weight_sum_total}{: the total sum of weights of background
#'        sequences.}
#'   }
#' 
#' @importFrom stats pbinom fisher.test
.calcMotifEnrichment <- function(motif_matrix, 
                                 df, 
                                 test = c("binomial", "fisher"), 
                                 verbose = FALSE){
  
    # checks
    if (!is.matrix(motif_matrix)) {
        stop("'motif_matrix' has to be a matrix")
    }
    if (!nrow(motif_matrix) == nrow(df)) {
        stop("'motif_matrix' and 'df' must have the same number of rows")
    }
    if (!is.null(rownames(motif_matrix)) && !is.null(rownames(df)) &&
        !identical(rownames(motif_matrix), rownames(df))) {
        stop("'motif_matrix' and 'df' must have identical rownames")
    }
    .checkDfValidity(df)
    method <- match.arg(test)
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("'verbose' has to be either TRUE or FALSE")
    }
  
    # total sum of sequence weights for fg and bg
    total_wgt_fg <- sum(df$kmer_weight[df$is_foreground])
    total_wgt_bg <- sum(df$kmer_weight[!df$is_foreground])
    
    # sum of sequence weights for fg and bg per TF (for seqs with hits)
    tf_wgt_fg <- apply(X = motif_matrix, MARGIN = 2,
                       FUN = function(x) {
                          sum((df$kmer_weight * x)[df$is_foreground])
                       })
    
    tf_wgt_bg <- apply(X = motif_matrix, 
                           MARGIN = 2, 
                           FUN = function(x){
                             sum((df$kmer_weight * x)[!df$is_foreground])
                           })
  
    # calculate motif enrichment
    if (method == "binomial") {
    
      if (verbose) {
          message("using binomial test to calculate ",
                  "log(p-values) for motif enrichments")
      }
     
      # follow Homer in terms of setting upper and lower limits for the 
      # prob value in pbinom (see logbinomial function in statistics.cpp file
      # and scoreEnrichmentBinomial function in Motif2.cpp file)
      
      # limits
      lower_prob_limit <- 1 / total_wgt_bg
      upper_prob_limit <- (total_wgt_bg - 1) / total_wgt_bg
      
      # prob
      prob <- tf_wgt_bg / total_wgt_bg
      
      if (any(i <- (prob < lower_prob_limit))) {
          warning("some probabilities (TF_bgSum/total_bgSum) have a ",
                  "value less than lower_prob_limit (example when TF_bgSum=0) ",
                  "and will be given a value of lower_prob_limit=1/total_bgSum")
          prob[i] <- lower_prob_limit
      }
      if (any(i <- (prob > upper_prob_limit))) {
          warning("some probabilities (TF_bgSum/total_bgSum) have a ",
                  "value greater than upper_prob_limit and will be given ",
                  "a value of upper_prob_limit=(total_bgSum-1)/total_bgSum")
          prob[i] <- upper_prob_limit
      }
      
      enr_logp <- pbinom(q = tf_wgt_fg - 1,
                         size = total_wgt_fg,
                         prob = prob,
                         lower.tail = FALSE, log.p = TRUE)
      
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
        enr_logp <- log(vapply(structure(seq_along(tf_wgt_fg),
                                         names = names(tf_wgt_fg)),
                               function(i) {
            cont_table <- rbind(c(tf_wgt_fg[i], total_wgt_fg - tf_wgt_fg[i]),
                                c(tf_wgt_bg[i], total_wgt_bg - tf_wgt_bg[i]))
            cont_table <- round(cont_table)
            fisher.test(x = cont_table, alternative = "greater")$p.value
        }, FUN.VALUE = numeric(1)))
      
    }
    
    return(data.frame(motif_name = names(enr_logp), 
                      log_p_value = enr_logp, 
                      fg_weight_sum = tf_wgt_fg, 
                      bg_weight_sum = tf_wgt_bg, 
                      fg_weight_sum_total = total_wgt_fg, 
                      bg_weight_sum_total = total_wgt_bg))
}


#' @title Do Binned Motif Enrichment Analysis with \code{monaLisa} a la \code{Homer}
#'
#' @description This function does a motif enrichment analysis on the set of sequences belonging
#'   to the same bin, using all the sequences in the rest of the bins as background. 
#'   In each enrichment analysis (per bin), \code{get_binned_motif_enrichment} uses other 
#'   functions within `monaLisa`, which do a motif enrichment analysis as implemented in `Homer` 
#'   (version 4.11). See Details for more.
#' 
#' @param seqs DNAStringSet object with sequences to test
#' @param bins factor vector of the same length and order as \code{seqs}, indicating the bin
#'   each sequence belongs to. See the \code{\link[monaLisa]{bin}} function. 
#' @param pwmL PWMatrixList object with PWMs of motifs.
#' @param enrichment_test A \code{character} scalar specifying the type of
#'   enrichment test to perform. One of \code{"binomial"} (default) or
#'   \code{"fisher"}. The enrichment test is one-sided (enriched in foreground).
#' @param frac_N_allowed A numeric scalar with the maximal fraction of N bases allowed 
#'   in a sequence (defaults to 0.7). Sequences with a fraction of N bases greater than
#'   this will be excluded from the analysis.
#' @param max_kmer_size the maximum kmer size to consider, when adjusting background weights
#'   for kmer composition compared to the foreground sequences (default 1-mer, 2-mers and 3-mers).
#' @param min.score the \code{min.score} parameter in \code{\link[monaLisa]{findMotifHits}}. It
#'   is the minimum score for counting a match (default set to 10).
#' @param match_method the \code{method} parameter in \code{\link[monaLisa]{findMotifHits}}. It specifies
#'   the method used for motif searching.
#' @param Ncpu Number of CPUs to use (default set to 1). This can enter the \code{findMotifHits}
#'   function when searching for motif matches. 
#' @param verbose A logical scalar. If \code{TRUE}, report progress.
#' @param ... Additional arguments for  \code{\link[monaLisa]{findMotifHits}}.
#'
#' @details This function implements a binned motif enrichment analysis. In each enrichment
#'   analysis, the sequences in a specific bin are used as foreground sequences to test for 
#'   motif enrichment using all other sequences (in the remaining bins) as background sequences.
#'   The function implements the \code{findMotifsGenome.pl} function from \code{Homer} version 4.11, 
#'   with \code{-size given -nomotif -mknown}, using the given list of known PWMs to test for 
#'   their enrichment. These \code{Homer} functions have been re-iplemented in R within \code{monaLisa} and 
#'   can reproduce the same output as \code{Homer} in this mode of use. Namely, weights that correct for
#'   GC, and k-mer composition differences between foreground and background are calculated for each
#'   sequence. Motif hits are scanned across the specified genome using \code{findMotifHits}, resulting
#'   in a sequence (row) by motif (column) matrix in ZOOPS mode, with a 1 entry if at least 1 hit is found, 
#'   and a 0 entry otherwise. For each motif, the weights of sequences that have a hit are summed separately
#'   for foreground (\code{fg_weight_sum}) and background (\code{bg_weight_sum}). The total foreground 
#'   (\code{fg_weight_sum_total}) and background (\code{bg_weight_sum_total}) sum is also calculated.
#'   
#'   The enrichment log p-value calculation can be done in two modes: Binomial (default) or Fisher's exact test: \itemize{
#'     \item{binomial}{: \code{pbinom(q = fg_weight_sum - 1, size = fg_weight_sum_total, prob = bg_weight_sum / bg_weight_sum_total, lower.tail = FALSE, log.p = TRUE)}}
#'     \item{fisher}{: \code{fisher.test(x = cont_table, alternative = "greater")}, where cont_table is the contingency table with 
#'        rows specifying foreground or background, and columns indicating hit or no hit for a particular motif. The entries in this table
#'        are the sum of the weights of sequences meeting the row and column specifications.}
#'   }
#'   
#'   TODO: - Use Ncpu to also parallelize across list of DFs (enrichment per bin)
#'                - parallellize the other functions using Ncpu..
#'                - change functions to do one-time calculations when using in binned mode?
#'                - filter seqs: make one-time thing.
#'                - think if we want to parametrize some of the pseudo-counts used (leaning towards not)
#'                - change the names of the assays in the SE or return the p-values and FDR (not -log10)
#'                - have a look again at using motif IDs and not names, but keeping both information
#'
#' @return A \code{SummarizedExperiment} object where the rows are the motifs and the columns are bins. The
#'   four assays are: \itemize{
#'   \item{p}{: -log10 P values}
#'   \item{FDR}{: -log10 false discovery rates}
#'   \item{enr}{: motif enrichments as Pearson residuals}
#'   \item{log2enr}{: motif enrichments as log2 ratios}
#' }
#'  
#' @importFrom TFBSTools ID name
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#'
#' @export
get_binned_motif_enrichment <- function(seqs, 
                                        bins, 
                                        pwmL, 
                                        enrichment_test = c("binomial", "fisher"), 
                                        frac_N_allowed = 0.7, 
                                        max_kmer_size = 3L, 
                                        min.score = 10, 
                                        match_method = "matchPWM.concat",
                                        Ncpu = 1L, 
                                        verbose = FALSE, 
                                        ...) {
  
    ### Questions: - do `seqs` have to have names?
    
    # checks
    # ... correct classes
    if (!is(seqs, "DNAStringSet")) {
        stop("class of 'seqs' must be DNAStringSet")
    }
    if (!is(bins, "factor")) {
        stop("'bins' must be of class 'factor'")
    }
    if (!is(pwmL, "PWMatrixList")) {
        stop("'pwmL' must be of class 'PWMatrixList'")
    }
    if (!is(frac_N_allowed, "numeric")) {
        stop("'frac_N_allowed' must be of class 'numeric'")
    }
    if (!is(max_kmer_size, "integer")) {
        stop("'max_kmer_size' must be of class 'integer'")
    }
    if (!is(Ncpu, "integer")) {
        stop("'Ncpu' must be of class 'integer'")
    }
    if (!is(verbose, "logical")) {
      stop("'verbose' must be of class 'logical'")
    }
    # ... supply/check for names
    if (is.null(names(seqs))) {
        if (verbose) {
            message("names(seqs) is empty, naming the sequences ...")
        }
        names(seqs) <- paste0("seq_", seq_along(seqs))
    }
    if (is.null(names(pwmL))) {
        stop("names(pwmL) is NULL, please name the PWMs, preferably with their unique ID.")
    }
    # ... number of sequences and bins match
    if (length(seqs) != length(bins)) {
        stop("'seqs' and 'bins' must be of equal length and in the same order")
    }
    
    # create data.frame of motif names and symbols
    TF_df <- data.frame(motif_ID = TFBSTools::ID(pwmL), 
                        motif_symbol = TFBSTools::name(pwmL))
    
    # create list of DataFrames, one for each bin
    bin_levels <- levels(bins)
    DF_list <- list()
    for (i in 1:length(bin_levels)) {
        is_foreground <- logical(length = length(seqs))
        is_foreground[bins == bin_levels[i]] <- TRUE
        df <- DataFrame(seqs = seqs,
                        is_foreground = is_foreground,
                        gc_frac = NA_real_,
                        gc_bin = NA_integer_,
                        gc_weight = NA_real_,
                        kmer_weight = NA_real_)
        attr(df, "err") <- NA
        
        DF_list[[i]] <- df
        names(DF_list)[i] <- bin_levels[i]
    }
    
    # filter 'bad' sequences per bin (can be done only once for all seqs--> need to change filtering function)
    if (verbose) {
        message("filtering out bad sequences ...")
    }
    DF_list <- lapply(DF_list, function(x){.filterSeqs(df = x, maxFracN = frac_N_allowed, verbose = verbose)})
    
    # calculate weight to adjust for GC differences between foreground and background per bin 
    # ... in this step, sequences may be filtered out (if a GC bin contains one sequence only, that sequence is filtered out).
    # ... Since the kept seqs may differ per bin, we select the kept seqs at the enrichment per bin.
    if (verbose) {
        message("Correcting for GC differences to the background sequences per bin ...")
    }
    DF_list <- lapply(DF_list, function(x){.calculateGCweight(df = x, verbose = verbose)})
    
    # update weight to in addition adjust for kmer composition differences between foreground and background per bin
    if (verbose) {
        message("Correcting for kmer differences to the background sequences per bin ...")
    }
    DF_list <- lapply(DF_list, function(x){.iterativeNormForKmers(df = x, max_kmer_size = max_kmer_size, verbose = verbose)})
    
    # get motif hits matrix in ZOOPS mode for all seqs
    if (verbose) {
        message("Scanning sequences for motif hits...")
    }
    hits <- findMotifHits(query = pwmL, subject = seqs, min.score = min.score,
                          method = match_method, Ncpu = Ncpu, ...)
    if (isEmpty(hits)) {
        stop("motif hits matrix is empty")
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
    
    # get log(p-values) for motif enrichment
    if (verbose) {
        message("Calculating motif enrichment per bin ...")
    }
    enrich_list <- lapply(DF_list, function(x) {
        .calcMotifEnrichment(motif_matrix = complete_hit_mat[rownames(x), , drop = FALSE],
                             df = x, test = enrichment_test, verbose = verbose)
    })
    
    # summarize results to return as SE
  
    # -log10(p-value)
    P <- do.call(cbind, lapply(enrich_list, function(bin_res) {
        D <- bin_res[, "log_p_value"]
        logpVals <- -D / log(10)
        names(logpVals) <- bin_res[, "motif_name"]
        logpVals
        # logpVals[order(names(logpVals))]
    }))
    
    # -log10(FDR)
    tmp <-  as.vector(10**(-P))
    fdr <- matrix(-log10(p.adjust(tmp, method = "BH")), nrow = nrow(P)) # add as parameter? the method of choice for multiple testing correction?
    dimnames(fdr) <- dimnames(P)
    fdr[which(fdr == Inf, arr.ind = TRUE)] <- max(fdr[is.finite(fdr)])
    
    # enrTF (Pearson residuals)
    enrTF <- do.call(cbind, lapply(enrich_list, function(bin_res) {
        D <- as.matrix(
            data.frame(
              fg_weight_sum = bin_res[, "fg_weight_sum"], 
              frac_fg_seq_with_motif = (bin_res[, "fg_weight_sum"] / bin_res[, "fg_weight_sum_total"]), 
              frac_bg_seq_with_motif = (bin_res[, "bg_weight_sum"] / bin_res[, "bg_weight_sum_total"]))
        )
        rownames(D) <-  bin_res[, "motif_name"]
        obsTF <- D[, "fg_weight_sum"]
        expTF <- D[, "fg_weight_sum"] / (D[, "frac_fg_seq_with_motif"] + 0.001) * (D[, "frac_bg_seq_with_motif"] + 0.001) # keep pseudo count fixed or make parameter? exp = numb_fg_total*bg_frac 
        enr <- (obsTF - expTF) / sqrt(expTF)
        enr[ is.na(enr) ] <- 0
        enr
    }))
    
    # log2enr
    log2enr <- do.call(cbind, lapply(enrich_list, function(bin_res) {
        D <- bin_res[, c("fg_weight_sum", "bg_weight_sum")]
        nTot <- c(bin_res[1, "fg_weight_sum_total"], bin_res[1, "bg_weight_sum_total"]) # Do I round the bg total to the nearest integer?
        D.norm <- t(min(nTot)*t(D)/nTot) # scale to smaller number (usually number of target sequences) # 
        DL <- log2(D.norm + 8) # keep pseudo count fixed? --> yes
        log2enr <- DL[, 1] - DL[, 2]
        names(log2enr) <- bin_res[, "motif_name"]
        log2enr
    }))
                       
    # return SummarizedExperiment
    se <- SummarizedExperiment(assays = list(p = P, FDR = fdr, enr = enrTF,
                                             log2enr = log2enr)) # better names? P and fdr are -log10(p-value), yes or transform into p-values
    m <- match(rownames(se), rownames(TF_df))
    rowData(se) <- TF_df[m, ]
    
    # ?order by p-value or enrichment in the future?
    return(se)
}

