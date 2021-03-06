# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Count k-mer pairs in sequences
#'
#' Over a set of sequences, for each k-mer a, count occurences of all k-mers b
#' within a defined distance downstream of the start of a.
#'
#' @param x DNAStringSet with sequences.
#' @param k An integer scalar with the length of sequences.
#' @param n An integer scalar defining the maximum downstream distance of
#'     second k-mers, relative to the start position of the first k-mer.
#' @param zoops A logical scalar. If TRUE, count each observed k-mer pair
#'     only once per sequence.
#'
#' @examples
#' countKmerPairs(Biostrings::DNAStringSet(c("AACCGGTT")), k = 2, n = 1)
#'
#' @return A numeric matrix with observed k-mer pairs counts.
#' @export
countKmerPairs <- function(x, k = 6L, n = 5L, zoops = FALSE) {
    .Call(`_monaLisa_countKmerPairs`, x, k, n, zoops)
}

