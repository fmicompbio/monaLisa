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
filterSeqs <- function(df, maxFracN = 0.7, minLength = 5L,
                        maxLength = 100000L, verbose = TRUE) {

    .checkDfValidity(df)
    if (!is.numeric(maxFracN) || length(maxFracN) != 1L || maxFracN < 0 || maxFracN > 1) {
        stop("'maxFracN' has to be a numerical scalar with a value in [0,1]")
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
#'   binned depending on GC content (hard-coded bin breaks). The numbers of
#'   foreground and background sequences in each bin are counted, and weights
#'   for background sequences in bin i are defined as:
#'   weight_i = (number_fg_seqs_i / number_bg_seqs_i) * (number_bg_seqs_total /
#'   number_fg_seqs_total)
#'
#' @param df a \code{DataFrame} with sequence information.
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
#'
#' @export
.calculateGCweight <- function(df, verbose = TRUE) {

  .checkDfValidity(df)
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' has to be either TRUE or FALSE")
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
  if (verbose) {
    message(sprintf("  %d of %d GC-bins used (have both foreground and background sequences)",
                    length(used_bins), length(gc_breaks) - 1))
  }
  
  # keep sequences belonging to these bins
  if (verbose) {
    message(sprintf("  %d of %d (%.1f%%) sequences filtered out from unused GC-bins",
                    sum(!df$gc_bin %in% used_bins), nrow(df),
                    100 * sum(!df$gc_bin %in% used_bins) / nrow(df)))
  }
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
  if (verbose) {
    message("  GC-weight calculation finished")
  }

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
#' @param kmer_freq a \code{list} with of matrices. The matrix at index \code{i}
#'   contains the counts of k-mers of length \code{i} (columns) for each
#'   sequence (rows).
#' @param g_oligos a \code{list} of \code{numeric} vectors; the element at index
#'   \code{i} contains the number of good (non-N-containing) k-mers of length
#'   \code{i} for each sequence.
#' @param kmer_seq_rc a \code{list} of character vectors; the element at index
#'   \code{i} contains the reverse complement sequences of all k-mers of length
#'   \code{i}.
#' @param kmer_weight a \code{numeric} vector with sequence weights adjusting
#'   for k-mer composition at the beginning of the iteration.
#' @param is_foreground logical vector of the same length as \code{seqs}.
#'   \code{TRUE} indicates that the sequence corresponds to the foreground set,
#'   \code{FALSE} indicates a background set sequence.
#' @param minimum_seq_weight Numeric scalar greater than zero giving the
#'   minimal weight of a sequence. The default value (0.001) was also used by
#'   \code{HOMER} (HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).
#' @param maximum_seq_weight Numeric scalar greater than zero giving the
#'   maximal weight of a sequence. The default value (1000) was also used by
#'   \code{HOMER} (1 / HOMER_MINIMUM_SEQ_WEIGHT constant in Motif2.h).
#'
#' @return a \code{DataFrame} of the same dimensions as the input \code{df},
#'   with the weight to adjust for k-mer composition and the current error
#'   stored in the column \code{kmer_weight} and the attribute \code{err}.
#'
#' @export
norm_for_kmer_comp <- function(kmer_freq,
                               g_oligos,
                               kmer_seq_rc,
                               kmer_weight,
                               is_foreground,
                               minimum_seq_weight = 0.001,
                               maximum_seq_weight = 1000) {

  # checks
  # arguments are all checked by iterate_norm_for_kmer_comp()

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

    # min values given by HOMER
    min_foreground_levels <- 0.5 / total_foreground
    min_background_levels <- 0.5 / total_background

    # Average the weight of a kmer with its reverse complement
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

    # ... HOMER check: if number of good oligos is > 0.5
    bg_g_oligos <- g_oligos[[k]][!is_foreground]
    g <- bg_g_oligos > 0.5
    bg_new_score[g] <- bg_new_score[g] / bg_g_oligos[g]

    # ... new weight for each background sequence
    # ... ... newWeight = newScore*currentWeight
    bg_cur_weight <- kmer_weight[!is_foreground]
    bg_new_weight <- bg_new_score * bg_cur_weight

    # ... HOMERs minimum and maximum weight cutoffs
    bg_new_weight[bg_new_weight < minimum_seq_weight] <- minimum_seq_weight
    bg_new_weight[bg_new_weight > maximum_seq_weight] <- maximum_seq_weight

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

    # ... update kmer_weight with bg_new_weight1
    kmer_weight[!is_foreground] <- bg_new_weight1

  }

  # return new kmer_weight and error
  list(kmer_weight = kmer_weight, err = err)
}




#' @title Adjust for k-mer composition (multiple iterations)
#'
#' @description Here we run `norm_for_kmer_comp` multiple times to converge to
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
#'   times to run \code{norm_for_kmer_comp}. the default is set to 160 (as in
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
#' @export
iterate_norm_for_kmer_comp <- function(df,
                                       max_kmer_size = 3L,
                                       minimum_seq_weight = 0.001,
                                       max_autonorm_iters = 160L,
                                       verbose = TRUE) {

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

  # run norm_for_kmer_comp() up to max_autonorm_iters times or
  # stop when new error is bigger than the error from the previous iteration
  if (verbose) {
    message("  starting iterative adjustment for k-mer composition (up to ",
            max_autonorm_iters, " iterations)")
  }

  res <- list()
  for (i in seq_len(max_autonorm_iters)) {

    # run norm_for_kmer_comp
    res <- norm_for_kmer_comp(kmer_freq = kmer_freq,
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

#' @title Get Motif Hits
#'
#' @param df a \code{DataFrame} with sequence information.
#' @param pwmL PWMatrixList object with PWMs of motifs
#' @param homerfile \code{character} scalar with path to homer2 binary, passed on
#'     to \code{findMotifHits}.
#' @param method the \code{method} parameter in the \code{findMotifHits} function
#' @param min.score the \code{min.score} parameter in the \code{findMotifHits} function
#' @param Ncpu number of CPUs to use
#' @param verbose A logical scalar. If \code{TRUE}, describe motif matrix being created
#'
#' @return a matrix where the rows are the sequences and the columns the motifs. This matrix
#'   consists of 0 and 1 entries for a motif hit or not, respectively.
#'   
#' @details TODO: should the returned matrix be a sparse matrix? and should we run this once for all
#'   sequences in all bins (when doing motif enrichment across bins)?
#'      
#' @importFrom S4Vectors DataFrame
#' 
#' @export
get_motif_hits_in_ZOOPS_mode <- function(df, 
                                         pwmL, 
                                         method=c("homer2"), 
                                         homerfile = findHomer("homer2"), 
                                         min.score = 10L, 
                                         Ncpu = 1L, 
                                         verbose = TRUE){
  
  .checkDfValidity(df)
  if (!is(pwmL, "PWMatrixList")) {
    stop("pwmL must be of class PWMatrixList")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' has to be either TRUE or FALSE")
  }
  
  # verbose=TRUE
  if (verbose) {
    message("creating matrix of motif hits in ZOOPS mode (sequences as rows and 
    motifs as columns), where 1 means at least one motif is present and 0 means no hit")
  }
  
  # find motifs
  hits <- monaLisa::findMotifHits(pwmL, df$seqs, method = "homer2", homerfile = homerfile, min.score = min.score)
  hits_matrix <- as.matrix(as.data.frame.matrix(table(seqnames(hits), as.character(hits$pwmname))))
  
  # match rows in df (add missing rows with no TF hits)
  missing <- rownames(df)[!rownames(df)%in%rownames(hits_matrix)]
  if(!isEmpty(missing)) {
    missing_mat <- matrix(data = 0, nrow = length(missing), ncol = ncol(hits_matrix))
    rownames(missing_mat) <- missing
    colnames(missing_mat) <- colnames(hits_matrix)
    # add to hits_matrix and order as in df
    hits_matrix <- rbind(hits_matrix, missing_mat)
    hits_matrix <- hits_matrix[rownames(df), ]
  }
  
  # ZOOPS mode: 1/0 matrix for hit or no hit
  w <- hits_matrix > 1
  if(!isEmpty(w)) {
    hits_matrix[w] <- 1
  }
  
  # return
  hits_matrix
  
}


#' @title Do Motif Enrichment
#'
#' @param motif_matrix matrix with 0 and 1 entries for absence or presence of 
#'   motif hits per sequence.
#' @param df a \code{DataFrame} with sequence information as returned by
#'   \code{iterate_norm_for_kmer_comp()}.
#' @param test type of test to do for the motif enrichment. By default it is the binomial, which 
#'   is what \code{Homer} uses by default. Fisher's exact test is another 
#'   option which allows for testing enrichment, without having to account for special cases of zero 
#'   background counts for a motif. \code{fisher.test} is used with \code{alternative="greater"}, making
#'   it a one-sided test for enrichment, as is the case with the binomial test. 
#' @param verbose A logical scalar. If \code{TRUE}, report motif enrichment test.
#' 
#' @return a \code{data.frame} containing the motifs as rows and the following columns: \itemize{
#'   \item{motif_name}{: the motif name}
#'   \item{log_p_value}{: the log p-value for enrichment. If test=binomial (default), this log p-value is identical to the one
#'      that Homer returns.}
#'   \item{fg_weight_sum}{: the sum of the weights of the foreGround sequences that have at least one instance of a specific motif (ZOOPS mode).}
#'   \item{bg_weight_sum}{: the sum of the weights of the backGround sequences that have at least one instance of a specific motif (ZOOPS mode).}
#'   \item{fg_weight_sum_total}{: the total sum of the weights of all foreGround sequences.}
#'   \item{bg_weight_sum_total}{: the total sum of the weights of all backGround sequences.}
#' }
#' 
#' @importFrom stats pbinom fisher.test
#' 
#' @export 
get_motif_enrichment <- function(motif_matrix=NULL, 
                                 df=NULL, 
                                 test=c("binomial", "fishers_exact"), 
                                 verbose = TRUE){
  
  # checks
  if (!nrow(motif_matrix) == nrow(df)) {
    stop("'motif_matrix' and 'df' must have the same number of rows")
  }
  if (!all(rownames(motif_matrix) == rownames(df))) {
    stop("'motif_matrix' and 'df' must have matching names (same order)")
  }
  if (is.null(motif_matrix) | is.null(df)) {
    stop("'motif_matrix' and 'df' have to be provided")
  }
  if (!is.matrix(motif_matrix)) {
    stop("'motif_matrix' has to be a matrix")
  }
  .checkDfValidity(df)
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' has to be either TRUE or FALSE")
  }
  method <- match.arg(test)

  # calculate sum of sequence weights for fg and bg per TF (for seqs with hits)
  total_foreground <- sum(df$kmer_weight[df$is_foreground])
  total_background <- sum(df$kmer_weight[!df$is_foreground])
  
  tf_foreground <- apply(X = motif_matrix, 
                         MARGIN = 2, 
                         FUN = function(x){
                           sum(df$kmer_weight[df$is_foreground & x==1])
                         })
  
  tf_background <- apply(X = motif_matrix, 
                         MARGIN = 2, 
                         FUN = function(x){
                           sum(df$kmer_weight[!df$is_foreground & x==1])
                         })

  # calculate motif enrichment
  if(method=="binomial") {
  
    # follow Homer in terms of setting upper and lower limits for the 
    # prob value in pbinom (see logbinomial function in statistics.cpp file
    # and scoreEnrichmentBinomial function in Motif2.cpp file)
    
    # verbose=TRUE
    if (verbose) {
      message("using the binomial test to calculate log(p-values) for motif enrichments")
    }
   
    # limits
    lower_prob_limit <- 1/total_background
    upper_prob_limit <- (total_background-1)/total_background
    
    # prob
    prob <- tf_background / total_background
    
    # print warnings if necessary
    if(sum(prob < lower_prob_limit) > 0){
      warning("some probabilities (TF_bgSum/total_bgSum) have a 
      value less than lower_prob_limit (example when TF_bgSum=0) 
      and will be given a value of lower_prob_limit=1/total_bgSum")
    }
    if(sum(prob > upper_prob_limit) > 0){
      warning("some probabilities (TF_bgSum/total_bgSum) have a 
      value greater than upper_prob_limit and will be given 
      a value of upper_prob_limit=(total_bgSum-1)/total_bgSum")
    }
    
    # update prob if necessary
    prob[prob < lower_prob_limit] <- lower_prob_limit
    prob[prob > upper_prob_limit] <- upper_prob_limit
    
    
    # enrichment
    enrichment_log_p_value <- stats::pbinom(q = tf_foreground - 1, size = total_foreground, prob = prob, lower.tail = FALSE, log.p = TRUE)
    
  }
  if(method=="fishers_exact") {
    
    # verbose=TRUE
    if (verbose) {
      message("using fisher's exact test (two-sided) to calculate log(p-values) for motif enrichments")
    }
    
    # index
    ind <- 1:length(tf_foreground)
    names(ind) <- names(tf_foreground)
    
    # contingency table per motif for fisher's exact test (x, y, z and w are rounded to the nearest integer):
    #          TF_hit  not_TF_hit
    #   is_fg     x         y
    #   is_bg     z         w
    #
    enrichment_log_p_value <- log(vapply(ind, function(i){
      
      # contingency table
      cont_table <-  matrix(data = c(tf_foreground[i], (total_foreground-tf_foreground[i]), 
                                     tf_background[i], (total_background-tf_background[i])), 
                            ncol = 2, 
                            byrow = TRUE)
      
      # round to integer for fisher's exact test
      cont_table <- round(cont_table)
      
      # fisher's exact test: get p-value
      stats::fisher.test(x = cont_table, alternative = "greater")$p.value}, 
      FUN.VALUE = 0.02)
    )
    
  }
  
  # # sort
  # o <- order(enrichment_log_p_value)
  o <- 1:length(enrichment_log_p_value)

  # return sorted log-pvalues and what was used to calculate them per motif
  data.frame(motif_name=names(enrichment_log_p_value)[o], 
             log_p_value=enrichment_log_p_value[o], 
             fg_weight_sum=tf_foreground[o], 
             bg_weight_sum=tf_background[o], 
             fg_weight_sum_total=rep(total_foreground, length(o)), 
             bg_weight_sum_total=rep(total_background, length(o)))
  
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
#' @param genome a BSgenome object indicating the genome to which the \code{seqs} 
#'   belong to, and which will be scanned for motif matches.
#' @param enrichment_test A \code{character} scalar specifying the type of
#'   enrichment test to perform. One of \code{"binomial"} (default) or
#'   \code{"fisher_exact"}. The enrichment test is one-sided (enriched in foreground).
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
#'     \item{fisher_exact}{: \code{fisher.test(x = cont_table, alternative = "greater")}, where cont_table is the contingency table with 
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
get_binned_motif_enrichment <- function(seqs = NULL, 
                                        bins = NULL, 
                                        pwmL = NULL, 
                                        genome = NULL, 
                                        enrichment_test = c("binomial", "fishers_exact"), 
                                        frac_N_allowed = 0.7, 
                                        max_kmer_size = 3L, 
                                        min.score = 10, 
                                        match_method = "matchPWM.concat",
                                        Ncpu = 1L, 
                                        verbose = TRUE, 
                                        ...) {
  # checks
  # ... missing arguments
  if(is.null(seqs)){
    stop("'seqs' argument must be supplied")
  }
  if(is.null(bins)){
    stop("'bins' argument must be supplied")
  }
  if(is.null(pwmL)){
    stop("'pwmL' argument must be supplied")
  }
  if(is.null(genome)){
    stop("'genome' argument must be supplied")
  }
  # ... correct classes
  if (!is(seqs, "DNAStringSet")) {
    stop("class of 'seqs' must be DNAStringSet")
  }
  if(!is(bins, "factor")){
    stop("'bins' must be of class 'factor'")
  }
  if(!is(pwmL, "PWMatrixList")){
    stop("'pwmL' must be of class 'PWMatrixList'")
  }
  if(!is(genome, "BSgenome")){
    stop("'genome' must be of class 'BSgenome'")
  }
  if(!is(frac_N_allowed, "numeric")){
    stop("'frac_N_allowed' must be of class 'numeric'")
  }
  if(!is(max_kmer_size, "integer")){
    stop("'max_kmer_size' must be of class 'integer'")
  }
  if(!is(Ncpu, "integer")){
    stop("'Ncpu' must be of class 'integer'")
  }
  if(!is(verbose, "logical")){
    stop("'verbose' must be of class 'logical'")
  }
  # ... supply/check for names
  if (is.null(names(seqs))) {
    if(verbose){
      message("names(seqs) is empty, naming the sequences ...")
    }
    names(seqs) <- paste0("seq_", seq_along(seqs))
  }
  if(is.null(names(pwmL))){
    stop("names(pwmL) is NULL, please name the PWMs, preferably with their unique ID.")
  }
  # ... number of sequences and bins match
  if (length(seqs) != length(bins)) {
    stop("'seqs' and 'bins' must be of equal length and in the same order")
  }
  
  # create data.frame of motif names and symbols
  TF_df <- data.frame(motif_ID=TFBSTools::ID(pwmL), 
                      motif_symbol=TFBSTools::name(pwmL))
  
  # create list of DataFrames, one for each bin
  bin_levels <- levels(bins)
  DF_list <- list()
  for(i in 1:length(bin_levels)){
    is_foreground <- logical(length = length(seqs))
    is_foreground[bins==bin_levels[i]] <- TRUE
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
  DF_list <- lapply(DF_list, function(x){monaLisa::filterSeqs(df = x, maxFracN = frac_N_allowed, verbose = verbose)})
  
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
  DF_list <- lapply(DF_list, function(x){monaLisa::iterate_norm_for_kmer_comp(df = x, max_kmer_size = max_kmer_size, verbose = verbose)})
  
  # get motif hits matrix in ZOOPS mode for all seqs
  if(verbose){
    message("Finding motif hits across seqs ...")
  }
  hits <- findMotifHits(query = pwmL, subject = seqs, min.score = min.score, method = match_method, Ncpu = Ncpu, genome = genome, ...)
  if(isEmpty(hits)){
    stop("motif hits matrix is empty")
  }
  mat <- as.matrix(as.data.frame.matrix(table(seqnames(hits), as.character(hits$pwmid))))
  # ... add missing rows (sequences that had no hit)
  missing_row_names <- names(seqs)[!names(seqs)%in%rownames(mat)]
  missing_mat <- matrix(data = 0, nrow = length(missing_row_names), ncol = ncol(mat))
  rownames(missing_mat) <- missing_row_names
  colnames(missing_mat) <- colnames(mat)
  complete_hit_mat <- rbind(mat, missing_mat)
  complete_hit_mat <- complete_hit_mat[names(seqs), , drop=FALSE]
  # ... ZOOPS mode
  w <- complete_hit_mat > 1
  complete_hit_mat[w] <- 1
  
  # get log(p-values) for motif enrichment
  if(verbose){
    message("Calculating motif enrichment per bin ...")
  }
  enrich_list <- lapply(DF_list, function(x){get_motif_enrichment(motif_matrix = complete_hit_mat[rownames(x), , drop=FALSE], df = x, test = enrichment_test, verbose = verbose)})
  
  # summarize results to return as SE

  # -log10(p-value)
  P <- do.call(cbind, lapply(enrich_list, function(bin_res){
    D <- bin_res[, "log_p_value"]
    logpVals <- -D/log(10)
    names(logpVals) <- bin_res[, "motif_name"]
    logpVals
    # logpVals[order(names(logpVals))]
  })
  )
  
  # -log10(FDR)
  tmp <-  as.vector(10**(-P))
  fdr <- matrix(-log10(p.adjust(tmp, method="BH")), nrow=nrow(P)) # add as parameter? the method of choice for multiple testing correction?
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
  })
  )
  
  # log2enr
  log2enr <- do.call(cbind, lapply(enrich_list, function(bin_res) {
    D <- bin_res[, c("fg_weight_sum", "bg_weight_sum")]
    nTot <- c(bin_res[1, "fg_weight_sum_total"], bin_res[1, "bg_weight_sum_total"]) # Do I round the bg total to the nearest integer?
    D.norm <- t(min(nTot)*t(D)/nTot) # scale to smaller number (usually number of target sequences) # 
    DL <- log2(D.norm + 8) # keep pseudo count fixed? --> yes
    log2enr <- DL[, 1] - DL[, 2]
    names(log2enr) <- bin_res[, "motif_name"]
    log2enr
  })
  )
                     
  # return SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(p=P, FDR=fdr, enr=enrTF, log2enr=log2enr)) # better names? P and fdr are -log10(p-value), yes or transform into p-values
  m <- match(rownames(se), rownames(TF_df))
  SummarizedExperiment::rowData(se) <- TF_df[m, ]
  
  # ?order by p-value or enrichment in the future?
  se
  
}

