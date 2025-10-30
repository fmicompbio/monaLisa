# Calculate motif enrichment

Given motif counts, foreground/background labels and weights for a set
of sequences, calculate the enrichment of each motif in foreground
compared to background. This function is called by
[`calcBinnedMotifEnrR()`](https://fmicompbio.github.io/monaLisa/reference/calcBinnedMotifEnrR.md)
for each bin.

The default type of test is `"fisher"`, which is also what `Homer` uses
if "-h" is specified for a hypergeometric test. Alternatively, a
binomial test can be used by `test = "binomial"` (what `Homer` does by
default). Using Fisher's exact test has the advantage that special cases
such as zero background counts are handled without ad-hoc adjustments to
the frequencies.

For `test = "fisher"`, `fisher.test` is used with
`alternative = "greater"`, making it a one-sided test for enrichment, as
is the case with the binomial test.

## Usage

``` r
.calcMotifEnrichment(
  motifHitMatrix,
  df,
  test = c("fisher", "binomial"),
  verbose = FALSE
)
```

## Arguments

- motifHitMatrix:

  Matrix with 0 and 1 entries for absence or presence of motif hits in
  each sequence.

- df:

  A `DataFrame` with sequence information as returned by
  [`.iterativeNormForKmers()`](https://fmicompbio.github.io/monaLisa/reference/dot-iterativeNormForKmers.md).

- test:

  Type of motif enrichment test to perform.

- verbose:

  A logical scalar. If `TRUE`, report on progress.

## Value

A `data.frame` containing the motifs as rows and the columns:

- motifName:

  : the motif name

- logP:

  : the log p-value for enrichment (natural logarithm). If
  `test="binomial"` (default), this log p-value is identical to the one
  returned by Homer.

- sumForegroundWgtWithHits:

  : the sum of the weights of the foreground sequences that have at
  least one instance of a specific motif (ZOOPS mode).

- sumBackgroundWgtWithHits:

  : the sum of the weights of the background sequences that have at
  least one instance of a specific motif (ZOOPS mode).

- totalWgtForeground:

  : the total sum of weights of foreground sequences.

- totalWgtBackground:

  : the total sum of weights of background sequences.
