# Calculate similarities between pairs of motifs.

For each pair of motifs, calculate the similarity defined as the maximal
Pearson's correlation coefficient between base frequencies over all
possible shifts (relative positions of the two matrices with at least
one overlapping position). If necessary matrices are padded on the sides
with background base frequencies (assuming all bases to have a frequency
of 0.25) to enable comparison of all positions in both matrices.

## Usage

``` r
motifSimilarity(
  x,
  y = NULL,
  method = c("R", "HOMER"),
  homerfile = findHomer("compareMotifs.pl"),
  homerOutfile = NULL,
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

- y:

  Either a
  [`PFMatrixList`](https://rdrr.io/pkg/TFBSTools/man/XMatrixList-class.html)
  or `NULL` (default). If `y = NULL`, then similarities will be
  calucalted for all pairs of motifs within `x`. Otherwise, `method`
  must be `"R"` and similarities will be calculated between any motif
  from `x` to any motif from `y`.

- method:

  A character scalar specifying the method for similarity calculations.
  Either `"R"` (pure R implementation) or `"HOMER"` (will call the
  `compareMotifs.pl` script from HOMER). Results are identical (apart
  from rounding errors), and the R implementation is usually faster and
  can be parallelized (`BPPARAM` argument).

- homerfile:

  Path to the HOMER script `compareMotifs.pl` (only used for
  `method = "HOMER"`.

- homerOutfile:

  A character scalar giving the file to save the similarity scores (only
  for `metho = "HOMER"`). If `NULL`, scores will be stored into a
  temporary file.

- BPPARAM:

  An optional
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance determining the parallel back-end to be used during
  evaluation (only used for `method = "R"`).

- verbose:

  A logical scalar. If `TRUE`, report on progress.

## Value

A matrix of Pearson's correlation coefficients for each pair of motifs.

## See also

[`bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html) used
for parallelization for `method = "R"`, documentation of HOMER's
`compareMotifs.pl` for details on `method = "HOMER"`.

## Examples

``` r
m <- rbind(A = c(12,  0,  0),
           C = c( 3,  2,  0),
           G = c( 0, 14,  0),
           T = c( 0,  0, 15))
pfms <- TFBSTools::PFMatrixList(
    TFBSTools::PFMatrix(name = "m1", profileMatrix = m),
    TFBSTools::PFMatrix(name = "m2", profileMatrix = m + 10),
    TFBSTools::PFMatrix(name = "m3", profileMatrix = m[, 3:1])
)
motifSimilarity(pfms)
#>           m1        m2        m3
#> m1 1.0000000 0.9997644 0.4382761
#> m2 0.9997644 1.0000000 0.4317897
#> m3 0.4382761 0.4317897 1.0000000
```
