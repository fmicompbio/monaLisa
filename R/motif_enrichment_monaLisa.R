#' @title Check input object
#'
#' @description Check if the input object is valid, i.e. is of the correct
#'   object type (DataFrame) and has all expected columns and attributes.
#'
#' @param x Input object to be checked.
#'
#' @return \code{TRUE} if \code{x} is valid, \code{FALSE} otherwise
is_valid_df <- function(x) {
  expected_cols <- c("seqs", "is_foreground",
                     "gc_frac", "gc_bin", "gc_weight",
                     "kmer_weight")
  expected_attrs <- c("err")
  if (!is(df, "DataFrame")) {
    message("'df' should be a DataFrame, but it is a ", class(df))
    return(FALSE)
  } else if (!all(expected_cols %in% colnames(df))) {
    message("'df' has to have columns: ",
            paste(expected_cols, collapse = ", "))
    return(FALSE)
  } else if (!all(expected_attrs %in% names(attributes(x)))) {
    message("'df' has to have attributes: ",
            paste(expected_attrs, collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}


#' @title Filter Bad Sequences
#'
#' @description We filter sequences similarly to how HOMER does it. Namely,
#'   sequences with more than 70% (default) N are removed.
#'
#' @param df a \code{DataFrame} with an attribute \code{err} and columns
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
#' @param frac A numeric scalar with the maximal fraction of N bases allowed in
#'   a sequence (defaults to 0.7).
#'
#' @return the filtered \code{df}.
#'
#' @importFrom Biostrings alphabetFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
#'
#' @export
filter_seqs <- function(df, frac = 0.7) {

  # checks
  if (!is_valid_df(df)) {
    stop("'df' is not a valid input object")
  }
  if (!is.numeric(frac) || length(frac) != 1L || frac < 0 || frac > 1) {
    stop("'frac' has to be a numerical scalar with a value in [0,1]")
  }

  # fraction of N bases per sequence
  frac_N <- alphabetFrequency(df$seqs)[, "N"] / nrow(df)

  # remove sequences with frac_N > frac
  w <- which(frac_N > frac)
  if (length(w) > 0) {
    df <- df[-w, ]
  }

  # return filtered df
  df
}


#' @title Get GC weights for background regions
#'
#' @description In this implementation we still follow HOMER. But this can be
#'   re-implemented in a more continuous fashion. Each sequence is put in a
#'   specific GC bin depending on its GC content. Then, for each GC bin, the
#'   number of foreGround and background sequences in that bin is calculated.
#'   Weights are calculated for the background sequences in bin i as follows:
#'   weight_i = (number_fg_seqs_i / number_bg_seqs_i) * (number_bg_seqs_total /
#'   number_fg_seqs_total)
#'
#' @param df a \code{DataFrame} with an attribute \code{err} and columns
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
#' @return a \code{DataFrame} of the same dimensions as the input \code{df},
#'   with the columns \code{gc_frac}, \code{gc_bin} and \code{gc_weight}
#'   filled in with the sequence GC content, GC bins they were
#'   assigned to, and the weight to correct for GC differences between
#'   foreground and background sequences, respecitvely.
#'
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#' @importFrom S4Vectors DataFrame
#'
#' @export
calculate_GC_weight <- function(df) {

  # checks
  if (!is_valid_df(df)) {
    stop("'df' is not a valid input object")
  }

  # calculate GC fraction for each sequence
  f_mono <- oligonucleotideFrequency(df$seqs, width = 1, as.prob = TRUE)
  df$gc_frac <- f_mono[, "G"] + f_mono[, "C"]

  # HOMER's GC breaks/bins
  gc_breaks <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8)

  # assign each sequence to a GC bin (for foreGround and background)
  df$gc_bin <- findInterval(x = df$gc_frac, vec = gc_breaks)

  # keep bins that have at least 1 foreground and 1 background sequence
  used_bins <- sort(intersect(df$gc_bin[df$is_foreground],
                              df$gc_bin[!df$is_foreground]))

  # keep sequences belonging to these bins
  df <- df[df$gc_bin %in% used_bins, ]

  # total number of foreGround and background sequences
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

  # return result
  df
}



#' @title Adjust for k-mer composition (single iteration)
#'
#' @description Here we correct the background sequence weights, adjusting for
#'   k-mer composition compared to the foreground sequences. This function
#'   implements a single iteration, and is called iteratively by
#'   \code{iterate_norm_for_kmer_comp} to get to the final set of adjusted
#'   weights, which will be the result of adjusting for GC and k-mer
#'   composition. We closely follow HOMER's \code{normalizeSequenceIteration()}
#'   function found in \code{Motif2.cpp}.
#'
#' @param df a \code{DataFrame} as returned by \code{calculate_GC_weight}.
#' @param max_kmer_size Integer scalar giving the maximum k-mer size to
#'   consider. The default is set to 3 (like in \code{HOMER}), meaning that
#'   k-mers of size 1, 2 and 3 are considered.
#' @param minimum_seq_weight Numeric scalar greater than zero giving the
#'   minimal weight of a sequence. The default value (0.001) was also used by
#'   \code{HOMER} (HOMER_MINIMUM_SEQ_WEIGHT  constant in Motif2.h).
#'
#' @return a \code{DataFrame} of the same dimensions as the input \code{df},
#'   with the weight to adjust for k-mer composition and the current error
#'   stored in the column \code{kmer_weight} and the attribute \code{err}.
#'
#' @importFrom Biostrings oligonucleotideFrequency reverseComplement
#'   DNAStringSet
#' @importFrom S4Vectors DataFrame
#'
#' @export
norm_for_kmer_comp <- function(df,
                               max_kmer_size = 3L,
                               minimum_seq_weight = 0.001) {

  # checks
  # df is checked by iterate_norm_for_kmer_comp()
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

  # set starting error
  error <- 0

  # set starting weight per sequence
  cur_weight <- df$kmer_weight

  # iterate over the kmer sizes: 1 till max_kmer_size
  for (cur_len in seq_len(max_kmer_size)) {

    # frequency of each k-mer of size cur_len per sequence
    kmer_freq <- oligonucleotideFrequency(x = df$seqs, width = cur_len)

    # number of good (non-N containing) oligos per sequence
    g_oligos <- rowSums(kmer_freq)

    # divide the current weight of each sequence by its g_oligos
    div_weight <- cur_weight / g_oligos

    # For each sequence multiply the frequency of each oligo by its div_weight.
    # This is the same as summing the weights for each oligo in a sequence.
    oligo_weights_per_seq <- sweep(x = kmer_freq, MARGIN = 1,
                                   STATS = div_weight, FUN = "*")

    # sum weights per oligo over foreground (foreground_levels) and background
    # (background_levels) sequences
    foreground_levels <- colSums(oligo_weights_per_seq[df$is_foreground, ])
    background_levels <- colSums(oligo_weights_per_seq[!df$is_foreground, ])
    kmers <- DNAStringSet(names(foreground_levels))

    # sum weights in foreground_levels and background_levels
    total_foreground <- sum(foreground_levels)
    total_background <- sum(background_levels)

    # min values given by HOMER
    min_foreground_levels <- 0.5 / total_foreground
    min_background_levels <- 0.5 / total_background

    # Average the weight of a kmer with its reverse complement
    rev_kmers_fg <- as.character(reverseComplement(x = kmers))
    rev_kmers_bg <- as.character(reverseComplement(x = kmers))
    f_level <- (foreground_levels + foreground_levels[rev_kmers_fg]) / 2
    b_level <- (background_levels + background_levels[rev_kmers_bg]) / 2

    # check if less than set minimum
    f_level[f_level < min_foreground_levels] <- min_foreground_levels
    b_level[b_level < min_background_levels] <- min_background_levels

    # Calculate normFactor (to be used to correct background sequences)
    norm_factors <- (f_level / b_level) * (total_background / total_foreground)

    # update error
    error <- error + sum((norm_factors - 1)^2 / length(foreground_levels))

    # calculate new weights for background sequences

    # ... sum the norm_factors of all k-mers per background sequence
    bg_new_score <- rowSums(sweep(x = kmer_freq[!df$is_foreground,
                                                names(norm_factors)],
                                  MARGIN = 2, STATS = norm_factors, FUN = "*"))

    # ... HOMER check: if number of good oligos is > 0.5
    bg_g_oligos <- g_oligos[!df$is_foreground]
    g <- bg_g_oligos > 0.5
    bg_new_score[g] <- bg_new_score[g] / bg_g_oligos[g]

    # ... new weight for each background sequence
    # ... ... newWeight = newScore*currentWeight
    bg_cur_weight <- cur_weight[!df$is_foreground]
    bg_new_weight <- bg_new_score * bg_cur_weight

    # ... HOMERs minimum weight cutoffs
    g <- bg_new_weight < minimum_seq_weight # lower bound
    bg_new_weight[g] <- minimum_seq_weight
    g <- bg_new_weight > 1 / minimum_seq_weight # upper bound
    bg_new_weight[g] <- 1 / minimum_seq_weight

    # ... penalty (still following HOMER)
    bg_penalty <- bg_new_weight
    g <- bg_penalty < 1
    bg_penalty[g] <- 1 / bg_penalty[g]
    bg_penalty <- bg_penalty^2

    # ... delta (still following HOMER)
    bg_delta <-  bg_new_weight - bg_cur_weight

    # ... bg_new_weight1 (still following HOMER)
    bg_new_weight1 <- bg_cur_weight + bg_delta
    g <- bg_penalty > 1 & ((bg_delta > 0 & bg_new_weight > 1) |
                             (bg_delta < 0 & bg_new_weight < 1))
    bg_new_weight1[g] <- bg_cur_weight[g] + bg_delta[g] / bg_penalty[g]

    # ... update cur_weight with bg_new_weight1
    cur_weight[!df$is_foreground] <- bg_new_weight1

  }

  # update kmer_weight and err in df
  df$kmer_weight <- cur_weight
  attr(df, "err") <- error

  # return new df, containing new k-mer weights and the error
  df
}




#' @title Adjust for k-mer composition (multiple iterations)
#'
#' @description Here we run `norm_for_kmer_comp` multiple times to converge to
#'   the final weights that will be used to correct the background
#'   sequences for k-mer composition differences compared to the foreground. We
#'   closely follow \code{HOMER}'s \code{normalizeSequence()} function found in
#'   \cod{Motif2.cpp}. Note that \code{HOMER} runs the
#'   \code{normalizeSequence()} one last time after going through all iterations
#'   or reaching a low error, which we do not do here.
#'
#' @param df a \code{DataFrame} as returned by \code{calculate_GC_weight}
#' @param max_autonorm_iters An integer scalar giving the maximum number if
#'   times to run \code{norm_for_kmer_comp}. the default is set to 160 (as in
#'   \code{HOMER}).
#'
#' @return a list containing: \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} containing the sequence GC content, GC bins they were
#'   assigned to, the weight to correct for GC differences between foreGround
#'   and background sequences, the weight to adjust for kmer composition, and
#'   the the error term} \item{sequenceNucleotides}{: a \code{DNAStringSet}
#'   object containing the raw sequences} }
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors DataFrame
#'
#' @export
iterate_norm_for_kmer_comp <- function(df,
                                       max_autonorm_iters = 160L) {

  # check
  if (!is_valid_df(df)) {
    stop("'df' is not a valid input object")
  }
  if (!is.integer(max_autonorm_iters) ||
      length(max_autonorm_iters) != 1L ||
      max_autonorm_iters < 1) {
    stop("'max_autonorm_iters' must be an integer scalar greater than zero")
  }

  # set starting df_final
  df_final <- df
  last_error <- Inf

  # initialize kmer_weight using gc_weight
  df_final$kmer_weight <- df_final$gc_weight

  # run norm_for_kmer_comp() up to max_autonorm_iters times or
  # stop when new error is bigger than the error from the previous iteration
  for (i in seq_len(max_autonorm_iters)) {

    # run norm_for_kmer_comp
    df_final <- norm_for_kmer_comp(df_final)
    cur_error <- attr(df_final, "err")

    # if current error is bigger than the last one, stop
    if (cur_error >= last_error) {
      break
    } else {
      last_error <- cur_error
    }
  }

  # return final weights
  df_final
}



#' @title Do Motif Enrichment Analysis with `monaLisa` (for simple fg vs bg
#'   situation, this wil lbe changed for bins later)
#'
#' @param seqs DNAStringSet object with sequences to analyze
#' @param is_forground logical vector of the same length as \code{seqs}.
#'   \code{TRUE} indicates that the sequence corresponds to the foreground set,
#'   \code{FALSE} indicates a background set sequence.
#'
#'
#'
#'
#' @export
run_monaLisa <- function(seqs, is_foreground) {

  # checks
  if (!is(seqs, "DNAStringSet")) {
    stop("class of 'seqs' must be DNAStringSet")
  }
  if (!is.logical(is_foreground)) {
    stop("'is_foreground' must be a logical vector with TRUE to indicate",
         "foreground sequences and FALSE to indicate background sequences")
  }
  if (length(seqs) != length(is_foreground)) {
    stop("'seqs' and 'is_foreground' must be of equal length")
  }
  if (is.null(names(seqs))) {
    names(seqs) <- paste0("seq_", seq_along(seqs))
  }

  # create list of the inputs
  df <- DataFrame(seqs = seqs,
                  is_foreground = is_foreground,
                  gc_frac = NA,
                  gc_bin = NA,
                  gc_weight = NA,
                  kmer_weight = NA)
  attr(df, "err") <- NA

  # filter 'bad' sequences
  df <- filter_seqs(df)

  # calculate weight to adjust for GC differences between foreground and
  # background
  df <- calculate_GC_weight(df)

  # calculate weight to in addition adjust for kmer composition differences
  df <- iterate_norm_for_kmer_comp(df)

}
