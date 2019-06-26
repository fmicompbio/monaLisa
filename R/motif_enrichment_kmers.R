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
#'   counts for each k-mer to avoid zero counts.
#'
#' @importFrom Biostrings DNAStringSet oligonucleotideFrequency
#' @importFrom tidyr %>%
#'
#' @export
getKmerFreq <- function(seqs, kmerLen = 4, MMorder = 1, pseudoCount = 1) {
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
        MMorder < kmerLen - 1L
    })

    #kmer frequencies
    kmerFreq <- oligonucleotideFrequency(seqs, width = kmerLen) %>% colSums

    #calculate probabilities for bg model (with a pseudocount of 1)
    p_long  <- oligonucleotideFrequency(seqs, width = MMorder + 1) %>% {colSums(. + pseudoCount)} %>% {./sum(.)}
    p_short <- oligonucleotideFrequency(seqs, width = MMorder)     %>% {colSums(. + pseudoCount)} %>% {./sum(.)}

    log2pMM <- sapply(names(kmerFreq), function(current.kmer){
        ii_long <- sapply(1:(nchar(current.kmer) - MMorder), function(i) {
            subseq(current.kmer, start = i, width = MMorder + 1)
        })
        ii_short <- sapply(2:(nchar(current.kmer) - MMorder), function(i) {
            subseq(current.kmer, start = i, width = MMorder)
        })
        sum(log2(p_long[ii_long])) - sum(log2(p_short[ii_short]))
    })
    kmerFreqMM <- (2**log2pMM)*sum(kmerFreq)

    data.frame(kmerFreq=kmerFreq, kmerFreqMM=kmerFreqMM)
}
