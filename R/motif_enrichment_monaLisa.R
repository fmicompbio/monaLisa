#' @title Filter Bad Sequences
#'
#' @description We filter sequences similarly to how HOMER does it. Namely,
#'   sequences with more than 70% (default) N are removed.
#'
#' @param inputList a `list` containing \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} with nrows equal to the number of sequences, with a column
#'   called `foreGround` that is 1 if a sequence is in the foreGround group, and
#'   0 otherwise, and a column containing the name of the sequnce.}
#'   \item{sequenceNucleotides}{: a \code{DNAStringSet} object containing the
#'   raw sequences} }
#'
#' @param frac the fraction of Ns allowed in a sequence (default set to 0.7)
#'
#' @return the filtered inputList.
#'
#' @importFrom Biostrings alphabetFrequency DNAStringSet
#'
#' @export
filter_seqs <- function(inputList=NULL, frac=0.7) {

  # checks
  if (is.null(inputList)){stop("'inputList' in NULL")}

  # fraction of Ns per sequence
  frac_N <- Biostrings::alphabetFrequency(inputList$sequenceNucleotides)[, "N"] / lengths(inputList$sequenceNucleotides)

  # remove sequences with fraction > frac
  w <- which(frac_N>frac)
  if (!isEmpty(w)){
    inputList$sequenceNucleotides <- inputList$sequenceNucleotides[-w]
    inputList$sequenceWeights <- inputList$sequenceWeights[-w, ]
  }

  # return filtered inputList
  inputList

}


#' @title Get GC weights for background regions
#'
#' @description In this implementation we still follow HOMER. But this can be
#'   re-implemented in a more continuous fashion. Each sequence is put in a
#'   specific GC bin depending on its GC content. Then, for each GC bin, the
#'   number of foreGround and backGround sequences in that bin is calculated.
#'   Weights are calculated for the backGround sequences in bin i as follows:
#'   weight_i =
#'   (number_fg_seqs_i/number_bg_seqs_i)*(number_bg_seqs_total/number_fg_seqs_total)
#'
#'
#' @param inputList a `list` containing \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} with nrows equal to the number of sequences, with a column
#'   called `foreGround` that is 1 if a sequence is in the foreGround group, and
#'   0 otherwise, and a column containing the name of the sequnce.}
#'   \item{sequenceNucleotides}{: a \code{DNAStringSet} object containing the
#'   raw sequences} }
#'
#' @return a list containing: \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} containing the sequence GC content, GC bins they were
#'   assigned to, and the weight to correct for GC differences between
#'   foreGround and backGround sequences} \item{sequenceNucleotides}{: a
#'   \code{DNAStringSet} object containing the raw sequences} }
#'
#' @importFrom Biostrings oligonucleotideFrequency DNAStringSet
#'
#' @export
get_GC_weight <- function(inputList=NULL) {

  # checks
  if (is.null(inputList)){stop("'inputList' is NULL")}

  # calculate GC fraction for each sequence
  f_mono <- Biostrings::oligonucleotideFrequency(inputList$sequenceNucleotides, width=1, as.prob=TRUE)
  gc_frac <- f_mono[,"G"]+f_mono[,"C"]

  # HOMER's GC breaks/bins
  gc_breaks <- c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8)

  # assign each sequence to a GC bin (for foreGround and backGround)
  gc_bin <- findInterval(x = gc_frac, vec = gc_breaks)

  # keep bins that have at least 1 foreGround and 1 backGround sequence
  bins <- unique(gc_bin)
  keep <- bins%in%unique(gc_bin[inputList$sequenceWeights$foreGround==1]) &
    bins%in%unique(gc_bin[inputList$sequenceWeights$foreGround==0])
  bins <- bins[keep]

  # keep sequences belonging to these bins
  keep_seq <- gc_bin%in%bins
  gc_bin <- gc_bin[keep_seq]
  inputList$sequenceNucleotides <- inputList$sequenceNucleotides[keep_seq]
  inputList$sequenceWeights <- inputList$sequenceWeights[keep_seq, ]
  gc_frac <- gc_frac[keep_seq]

  # total number of foreGround and backGround sequences
  total_fg <- sum(inputList$sequenceWeights$foreGround==1)
  total_bg <- sum(inputList$sequenceWeights$foreGround==0)

  # calculate GC weight per bin
  weight_per_bin <- sapply(bins, function(b){
    n_fg_b <- sum(gc_bin[inputList$sequenceWeights$foreGround==1]%in%b) # number of fg seqs in b
    n_bg_b <- sum(gc_bin[inputList$sequenceWeights$foreGround==0]%in%b) # number of bg seqs in b
    (n_fg_b/n_bg_b)* (total_bg/total_fg)
  })

  # assign calculated GC weight to each backGround sequence (foreGround get a weight of 1)
  df <- inputList$sequenceWeights
  df$gc_fraction <- gc_frac
  df$gc_bin <- gc_bin
  df$gc_weight <- rep(1, nrow(df))

  for (i in 1:length(bins)){
    b <- bins[i]
    w <- weight_per_bin[i]
    df$gc_weight[gc_bin%in%b & df$foreGround==0] <- w
  }

  # update list
  inputList$sequenceWeights <- df

  # return list
  inputList

}



#' @title Adjust for Kmer Composition
#'
#' @description Here we correct the backGround sequence weights, adjusting for
#'   kmer composition compared to the foreGroud sequences. This function
#'   implements one run, and it will be iteratively called multiple times by
#'   another function to get to the final set of adjusted weights. These weights
#'   will the result of adjusting for GC and kmer composition. We closely follow
#'   `HOMER`'s `normalizeSequenceIteration()` function found in `Motif2.cpp`.
#'
#' @param inputList a `list`.
#' @param maxKmerSize the maximum kmer size to consider. The default is set to 3
#'   (like in `HOMER`), meaning that kmers of size 1,2 and 3 are considered.
#'
#' @return a list containing: \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} containing the sequence GC content, GC bins they were
#'   assigned to, the weight to correct for GC differences between foreGround
#'   and backGround sequences, the weight to adjust for kmer composition, and
#'   the the error term} \item{sequenceNucleotides}{: a \code{DNAStringSet}
#'   object containing the raw sequences} }
#'
#' @importFrom Biostrings oligonucleotideFrequency reverseComplement
#'   DNAStringSet
#'
#' @export
norm_for_kmer_comp <- function(inputList=NULL, maxKmerSize=3){

  # checks


  # given cutoff values from HOMER
  HOMER_MINIMUM_SEQ_WEIGHT <- 0.001 # value from 'Motif2.h' file

  # set starting error
  error <- 0

  # set starting weight per sequence
  cur_weight <- inputList$sequenceWeights$KmerAdjWeight 

  # iterate over the kmer sizes: 1 till maxKmerSize
  for (curLen in 1:maxKmerSize){

    # frequency of each kmer of size curLen per sequence
    kmer_freq <- Biostrings::oligonucleotideFrequency(x = inputList$sequenceNucleotides, width = curLen)

    # number of good oligos per sequence
    g_oligos <- rowSums(kmer_freq)

    # divide the current weight of each sequence by its g_oligos
    div_weight <- cur_weight/g_oligos

    # for each sequence multiply the frequency of each oligo by its div_weight.
    # This is the same as summing the weights for each oligo in a sequence.
    oligo_wights_per_seq <- sweep(x = kmer_freq, MARGIN = 1, STATS = div_weight, FUN = "*")

    # sum all weights per oligo for foreGround (target_levels) and backGround (background_levels)
    target_levels <- colSums(oligo_wights_per_seq[inputList$sequenceWeights$foreGround==1, ])
    background_levels <- colSums(oligo_wights_per_seq[inputList$sequenceWeights$foreGround==0, ])

    # sum weights in target_levels and background_levels
    total_target <- sum(target_levels)
    total_background <- sum(background_levels)

    # min values given by HOMER
    min_target_levels <- 0.5/total_target
    min_background_levels <- 0.5/total_background

    # Average the weight of a kmer with its reverse complement
    rev_kmers_fg <- as.character(Biostrings::reverseComplement(x = Biostrings::DNAStringSet(names(target_levels))))
    rev_kmers_bg <- as.character(Biostrings::reverseComplement(x = Biostrings::DNAStringSet(names(background_levels))))
    t_level <- (target_levels + target_levels[rev_kmers_fg])/2
    b_level <- (background_levels + background_levels[rev_kmers_bg])/2

    # check if less than set minimum
    t_level[t_level<min_target_levels] <- min_target_levels
    b_level[b_level<min_target_levels] <- min_background_levels

    # Calculate normFactor (to be used to correct backGround sequences)
    norm_factors <- (t_level/b_level)* (total_background/total_target)

    # update error
    error <- error + sum((norm_factors-1)^2/length(target_levels))

    # calculate new weights for background sequences

    # ... sum the norm_factors of all kmers per backGround sequence
    bg_new_score <- rowSums(sweep(x = kmer_freq[inputList$sequenceWeights$foreGround==0, names(norm_factors)], MARGIN = 2, STATS = norm_factors, FUN = "*"))

    # ... HOMER check: if number of good oligos is > 0.5
    bg_g_oligos <- g_oligos[inputList$sequenceWeights$foreGround==0]
    g <- bg_g_oligos > 0.5
    bg_new_score[g] <- bg_new_score[g]/bg_g_oligos[g]

    # ... new weight for each background sequence
    # ... ... newWeight = newScore*currentWeight
    bg_cur_weight <- cur_weight[inputList$sequenceWeights$foreGround==0]
    bg_new_weight <- bg_new_score*bg_cur_weight

    # ... HOMERs minimum weight cutoffs
    g <- bg_new_weight < HOMER_MINIMUM_SEQ_WEIGHT # lower bound
    bg_new_weight[g] <- HOMER_MINIMUM_SEQ_WEIGHT
    g <- bg_new_weight > 1/HOMER_MINIMUM_SEQ_WEIGHT # upper bound
    bg_new_weight[g] <- 1/HOMER_MINIMUM_SEQ_WEIGHT

    # ... penalty (still following HOMER)
    bg_penalty <- bg_new_weight
    g <- bg_penalty < 1
    bg_penalty[g] <- 1/bg_penalty[g]
    bg_penalty <- bg_penalty^2

    # ... delta (still following HOMER)
    bg_delta <-  bg_new_weight - bg_cur_weight

    # ... newWeight1 (still following HOMER)
    bg_new_weight1 <- bg_cur_weight + bg_delta
    g <- bg_penalty > 1 & ((bg_delta>0 & bg_new_weight>1) | (bg_delta<0 & bg_new_weight < 1))
    bg_new_weight1[g] <- bg_cur_weight[g] + bg_delta[g]/bg_penalty[g]

    # ... update cur_weight with newWeight1
    cur_weight[inputList$sequenceWeights$foreGround==0] <- bg_new_weight1

  }

  # update KmerAdjWeight in inputList$sequenceWeights
  inputList$sequenceWeights$KmerAdjWeight <- cur_weight
  inputList$error <- error

  # return new inputList, containing new kmer weights, and the error
  inputList

}




#' @title Get Final Sequence Weights
#'
#' @description Here we run `norm_for_kmer_comp` multiple times to get
#'   the final set of weights that will be used to correct the background
#'   sequences for kmer composition differences compared to the foreGround. We
#'   closely follow `HOMER`'s `normalizeSequence()` function found in
#'   `Motif2.cpp`. Note that `HOMER` runs the `normalizeSequence()` one last
#'   time after going through all iterations or reaching a low error, which we
#'   do not do here.
#'
#' @param inputList a `list` resulting from running `get_GC_weight`
#' @param max_autonorm_iters the maximum number if times to run
#'   `norm_for_kmer_comp`. the default is set to 160 (as in `HOMER`).
#' @param last_error teh starting value for the last error. Once
#'   `norm_for_kmer_comp` is run, the aim is to stop when teh current
#'   error is bigger than the last error.
#'
#' @return a list containing: \itemize{ \item{sequenceWeights}{: a
#'   \code{dataframe} containing the sequence GC content, GC bins they were
#'   assigned to, the weight to correct for GC differences between foreGround
#'   and backGround sequences, the weight to adjust for kmer composition, and
#'   the the error term} \item{sequenceNucleotides}{: a \code{DNAStringSet}
#'   object containing the raw sequences} }
#'
#' @importFrom Biostrings DNAStringSet
#'
#' @export
iterate_norm_for_kmer_comp <- function(inputList=NULL, max_autonorm_iters=160, last_error=1e100){

  # check

  # set starting l_final
  l_final <- inputList

  # set current KmerAdjWeight equal to the gc_weight. This will change later
  l_final$sequenceWeights$KmerAdjWeight <- l_final$sequenceWeights$gc_weight

  # run norm_for_kmer_comp() max_autonorm_iters times or stop when error
  # is bigger than previous run
  for (i in 1:max_autonorm_iters){

    # run norm_for_kmer_comp
    l_final <- norm_for_kmer_comp(inputList = l_final)
    cur_error <- l_final$error

    # if current error is bigger than the last one, stop
    if (cur_error >= last_error){
      break
    } else {
      last_error <- cur_error
    }
  }

  # return final weights
  l_final

}



#' @title Do Motif Enrichment Analysis with `monaLisa` (for simple fg vs bg
#'   situation, this wil lbe changed for bins later)
#'
#'
#'
#'
#'
#'
#' @export
run_monaLisa <- function(seqs=NULL, foreGround=NULL){

  # checks
  if (class(seqs)!="DNAStringSet"){stop("class of 'seqs' must be DNAStringSet")}
  if (class(foreGround)!="numeric"){stop("'foreGround' must be a numeric vector with 1 to indicate foreGround sequences and 0 to indicate backGround sequences")}
  if (length(unique(foreGround))!=2){stop("make sure that the 'foreGround' vector only contains 1s and 0s")}
  if (!all(foreGround%in%c(1,0))){stop("make sure that the 'foreGround' vector only contains 1s and 0s")}
  if (length(seqs)!=length(foreGround)){stop("'seqs' and 'foreGround' must be of equal length")}
  if (is.null(names(seqs))){
    nm <- paste0("seq_", 1:length(seqs))
  } else {
    nm <- names(seqs)
  }

  # create list of the inputs
  df <- data.frame(seqName=nm, foreGround=foreGround)
  l <- list(sequenceWeights=df, sequenceNucleotides=seqs)

  # filter 'bad' sequences
  l <- filter_seqs(inputList = l) ## CHANGE FUNCTION ABOVE

  # calculate weight to adjust for GC differences between foreGround and
  # backGround
  l <- get_GC_weight(inputList = l)

  # calculate weight to in addition adjust for kmer composition differences
  l <- iterate_norm_for_kmer_comp(inputList = l)

}
