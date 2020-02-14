#' monaLisa - MOtif aNAlysis with Lisa.
#'
#' \pkg{monaLisa} is a collection of tools that simplify motif enrichment
#' analyses in genomic regions of interest.
#'
#' She makes use of her father Homer (http://homer.ucsd.edu/homer/index.html)
#' and other algorithms to search for motif hits and look for enriched motifs in
#' sets of genomic regions, compared to all other regions.
#'
#' Known motifs can for example be obtained from a collection of transcription
#' factor binding site specificities, such as \pkg{JASPAR2018}.
"_PACKAGE"

## usethis namespace: start
#' @useDynLib monaLisa, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import methods
#' @import IRanges
#' @import XVector
#' @import S4Vectors
#' @import Biostrings
## usethis namespace: end
NULL
