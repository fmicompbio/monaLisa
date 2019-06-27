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
#'
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @importFrom tidyr %>%
#' @importFrom XVector subseq
#'
#' @export
getKmerFreq <- function(seqs, kmerLen = 4, MMorder = 2, pseudoCount = 1) {
    if (is.character(seqs))
        seqs <- DNAStringSet(seqs)
    stopifnot(is(seqs, "DNAStringSet"))
    stopifnot(exprs = {
        is.numeric(kmerLen)
        length(kmerLen) == 1L
        round(kmerLen, 0L) == kmerLen
    })
    stopifnot(exprs = {
        is.numeric(MMorder)
        length(MMorder) == 1L
        round(MMorder, 0L) == MMorder
        MMorder > 1
        MMorder < kmerLen
    })

    #kmer frequencies
    kmerFreq <- oligonucleotideFrequency(seqs, width = kmerLen) %>% colSums

    #calculate probabilities for bg model (with a pseudocount)
    p_long  <- oligonucleotideFrequency(seqs, width = MMorder)      %>% {colSums(. + pseudoCount)} %>% {./sum(.)}
    p_short <- oligonucleotideFrequency(seqs, width = MMorder - 1L) %>% {colSums(. + pseudoCount)} %>% {./sum(.)}

    log2pMM <- sapply(names(kmerFreq), function(current.kmer){
        n <- nchar(current.kmer) - MMorder + 1L
        ii_long <- substr(rep(current.kmer, n), start = 1:n, stop = 0:(n - 1L) + MMorder)
        ii_short <- substr(rep(current.kmer, n - 1L), start = 2:n, stop = 1:(n - 1L) + MMorder - 1L)
        sum(log2(p_long[ii_long])) - sum(log2(p_short[ii_short]))
    })
    kmerFreqMM <- (2**log2pMM)*sum(kmerFreq)

    data.frame(kmerFreq=kmerFreq, kmerFreqMM=kmerFreqMM)
}
