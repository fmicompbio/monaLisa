getKmerFreq <- function(seq, kmerLen=4, MMorder=1, pseudoCount=0) {
    #seq is DNAStringSet
    #only works if Markov order at least twice smaller than kmerLen
    require(Biostrings)
    require(tidyverse)
    #kmer frequencies
    kmerFreq <- oligonucleotideFrequency(seq, width = kmerLen) %>% colSums

    #calculate probabilities for bg model (with a pseudocount of 1)
    p_long <- oligonucleotideFrequency(seq, width=MMorder+1) %>% {colSums(. + pseudoCount)} %>% {./sum(.)}
    p_short <- oligonucleotideFrequency(seq, width=MMorder) %>% {colSums(. + pseudoCount)} %>% {./sum(.)}

    log2pMM <- sapply(names(kmerFreq), function(current.kmer){
        ii_long <- sapply(1:(nchar(current.kmer)-(MMorder)), function(i){
            subseq(current.kmer, start = i, width = MMorder + 1)
        })
        ii_short <- sapply(2:(nchar(current.kmer)-MMorder), function(i){
            subseq(current.kmer, start = i, width = MMorder)
        })
        sum(log2(p_long[ii_long])) - sum(log2(p_short[ii_short]))
    })
    kmerFreqMM <- (2**log2pMM)*sum(kmerFreq)
    data.frame(kmerFreq=kmerFreq, kmerFreqMM=kmerFreqMM)
}
