# Calculate similarities between motifs and k-mers.

For each motif, calculate it's similarity to all k-mers of length
`kmerLen`, defined as the maximal probability of observing the k-mer
given the base frequencies of the motif (the maximum is taken over for
all possible ungapped alignments between motif and k-mer). If necessary
matrices are padded on the sides with background base frequencies
(assuming all bases to have a frequency of 0.25).

## Usage

``` r
motifKmerSimilarity(
  x,
  kmerLen = 5,
  kmers = NULL,
  includeRevComp = FALSE,
  BPPARAM = SerialParam(),
  verbose = FALSE
)
```

## Arguments

- x:

  Either a
  [`PFMatrixList`](https://rdrr.io/pkg/TFBSTools/man/XMatrixList-class.html),
  or a character scalar with a file containing motifs in HOMER format
  (used directly `method = "HOMER"`, loaded into a
  [`PFMatrixList`](https://rdrr.io/pkg/TFBSTools/man/XMatrixList-class.html)
  by
  [`homerToPFMatrixList`](https://fmicompbio.github.io/monaLisa/reference/homerToPFMatrixList.md)
  for `method = "R"`).

- kmerLen:

  A `numeric` scalar giving the k-mer length.

- kmers:

  Either a character vector of k-mers for which to calculate the
  similarity to each motif, or `NULL`, in which case all k-mers of
  length `kmerLen` are used.

- includeRevComp:

  A `logical` scalar. If set to `TRUE`, each k-mer as well as its
  reverse complement is compared to each motif, and the larger of the
  two similarities is returned.

- BPPARAM:

  An optional
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance determining the parallel back-end to be used during
  evaluation.

- verbose:

  A logical scalar. If `TRUE`, report on progress.

## Value

A matrix of probabilties for each motif - k-mer pair.

## See also

[`bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html) used
for parallelization.

## Examples

``` r
m <- rbind(A = c(12,  0,  0),
           C = c( 3,  2,  0),
           G = c( 0, 14,  0),
           T = c( 0,  0, 15))
pfms <- TFBSTools::PFMatrixList(
    TFBSTools::PFMatrix(name = "m1", profileMatrix = m),
    TFBSTools::PFMatrix(name = "m2", profileMatrix = m[, 3:1])
)
motifKmerSimilarity(pfms, kmerLen = 3)[, c("AGT", "TGA")]
#>       AGT    TGA
#> m1 0.7000 0.0625
#> m2 0.0625 0.7000
```
